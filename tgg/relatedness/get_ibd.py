#!/usr/bin/env python

import argparse
import hail as hl
import io
import logging

from collections import defaultdict
from os.path import dirname
from gnomad.sample_qc.pipeline import filter_rows_for_qc
from gnomad.utils.file_utils import file_exists
from gnomad.utils.reference_genome import get_reference_genome
from gnomad_qc.v2.resources.sample_qc import qc_mt_path
from gnomad_qc.v3.resources.sample_qc import qc
from typing import Dict, List, Tuple

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("check pedigree")
logger.setLevel(logging.INFO)

logger.info("Setting hail flag to avoid array index out of bounds error...")
# Setting this flag isn't generally recommended, but is needed (since at least Hail version 0.2.75) to avoid an array index out of bounds error until changes are made in future versions of Hail
# TODO: reassess if this flag is still needed for future versions of Hail
hl._set_flags(no_whole_stage_codegen="1")


def subset_samples(
    input_mt: hl.MatrixTable,
    pedigree: hl.Table,
    sex_ht: hl.Table,
    output_dir: str,
    output_name: str,
) -> Tuple[hl.MatrixTable, hl.Table, list, list]:
    """
    Filter the MatrixTable and sex Table to only samples in the pedigree.

    :param input_mt: MatrixTable
    :param pedigree: Pedigree file from seqr loaded as a Hail Table
    :param sex_ht: Table of inferred sexes for each sample
    :param output_dir: Path to directory to output results
    :param output_name: Output prefix to use for results
    :return: MatrixTable and sex ht subsetted to the samples given in the pedigree, list of samples in the pedigree, list of samples in the VCF
    """
    # Get sample names to subset from the pedigree
    samples_to_subset = hl.set(pedigree.Individual_ID.collect())
    # Subset mt and ht to samples in the pedigree
    mt_subset = input_mt.filter_cols(samples_to_subset.contains(input_mt["s"]))
    sex_ht = sex_ht.filter(samples_to_subset.contains(sex_ht["s"]))

    # Filter to variants that have at least one alt call after the subsetting
    mt_subset = mt_subset.filter_rows(hl.agg.any(mt_subset.GT.is_non_ref()))

    # Check that the samples in the pedigree are present in the VCF subset and output samples that are missing
    out_missing_samples = hl.hadoop_open(
        f"{output_dir}/{output_name}_missing_samples_in_subset.txt", "w"
    )
    expected_samples = pedigree.Individual_ID.collect()
    vcf_samples = mt_subset.s.collect()

    missings = set(expected_samples) - set(vcf_samples)
    for i in missings:
        out_missing_samples.write(i + "\n")
    out_missing_samples.close()

    return (mt_subset, sex_ht, expected_samples, vcf_samples)


def remap_samples(
    original_mt_path: str,
    input_mt: hl.MatrixTable,
    pedigree: hl.Table,
    inferred_sex: str,
) -> Tuple[hl.MatrixTable, hl.Table]:
    """
    Rename `s` col in the MatrixTable and inferred sex ht.

    :param original_mt_path: Path to original MatrixTable location
    :param input_mt: MatrixTable 
    :param pedigree: Pedigree file from seqr loaded as a Hail Table
    :param inferred_sex: Path to text file of inferred sexes
    :return: mt and sex ht with sample names remapped
    """
    base_path = "/".join(dirname(original_mt_path).split("/")[:-1]) + ("/base/projects")
    project_list = list(set(pedigree.Project_GUID.collect()))

    # Get the list of hts containing sample remapping information for each project
    remap_hts = []

    sex_ht = hl.import_table(inferred_sex)

    for i in project_list:
        remap = f"{base_path}/{i}/{i}_remap.tsv"
        if hl.hadoop_is_file(remap):
            remap_ht = hl.import_table(remap)
            remap_ht = remap_ht.key_by("s", "seqr_id")
            remap_hts.append(remap_ht)

    logger.info("Found %d projects that need to be remapped.", len(remap_hts))
    
    if len(remap_hts) > 0:
        ht = remap_hts[0]
        for next_ht in remap_hts[1:]:
            ht = ht.join(next_ht, how="outer")

        # If a sample has a non-missing value for seqr_id, rename it to the sample name for the mt and sex ht
        ht = ht.key_by("s")
        input_mt = input_mt.annotate_cols(seqr_id=ht[input_mt.s].seqr_id)
        input_mt = input_mt.key_cols_by(
            s=hl.if_else(hl.is_missing(input_mt.seqr_id), input_mt.s, input_mt.seqr_id)
        )

        sex_ht = sex_ht.annotate(seqr_id=ht[sex_ht.s].seqr_id).key_by("s")
        sex_ht = sex_ht.key_by(
            s=hl.if_else(hl.is_missing(sex_ht.seqr_id), sex_ht.s, sex_ht.seqr_id)
        )
    else:
        sex_ht = sex_ht.key_by("s")

    return input_mt, sex_ht


def ld_prune(input_mt: hl.MatrixTable, build: str, gnomad_ld: bool) -> hl.MatrixTable:
    """
    LD prune the MatrixTable.

    :param input_mt: MatrixTable
    :param build: Build for the input MatrixTable
    :param gnomad_ld: Whether or not to use LD data from gnomAD dataset for the pruning step
    :return: ld-pruned MatrixTable
    """
    if gnomad_ld == False:
        mm_pruned = hl.ld_prune(input_mt.GT, r2=0.1)
        input_mt = input_mt.filter_rows(hl.is_defined(mm_pruned[input_mt.row_key]))
    else:
        # Borrow from gnomAD ld pruning
        if build == "GRCh37":
            pruned_mt = hl.read_matrix_table(qc_mt_path("joint", ld_pruned=True))

        elif build == "GRCh38":
            pruned_mt = hl.read_matrix_table(qc.path)

        input_mt = input_mt.filter_rows(
            hl.is_defined(pruned_mt.index_rows(input_mt.row_key))
        )

    return input_mt


def add_project_and_family_annotations(
    ht: hl.Table, seqr_projects: dict, family_ids: dict
) -> hl.Table:

    """
    Add seqr project and family ID annotations to the kinship table.

    :param ht: Hail Table of kinship values
    :param seqr_projects: Dictionary of seqr projects for each sample
    :param family_ids: Dictionary of family ids for each sample
    :return: Table with seqr project and family id annotations added
    """
    # Add annotation for seqr projects of sample i and sample j
    hl_seqr_projects = hl.literal(seqr_projects)
    ht = ht.annotate(
        seqr_proj_i=hl_seqr_projects.get(ht.i), seqr_proj_j=hl_seqr_projects.get(ht.j),
    )

    # Add annotation for family ids of sample i and sample j
    hl_family_ids = hl.literal(family_ids)
    ht = ht.annotate(
        fam_id_i=hl_family_ids.get(ht.i), fam_id_j=hl_family_ids.get(ht.j),
    )

    return ht


def filter_kin_ht(
    ht: hl.Table,
    out_summary: io.TextIOWrapper,
    first_degree_pi_hat: float = 0.40,
    grandparent_pi_hat: float = 0.20,
    grandparent_ibd1: float = 0.25,
    grandparent_ibd2: float = 0.15,
) -> hl.Table:

    """
    Filter the kinship table to relationships of grandparents and above.

    :param ht: hl.Table
    :param out_summary: Summary file with a summary statistics and notes
    :param first_degree_pi_hat: Minimum pi_hat threshold to use to filter the kinship table to first degree relatives
    :param grandparent_pi_hat: Minimum pi_hat threshold to use to filter the kinship table to grandparents
    :param grandparent_ibd1: Minimum IBD1 threshold to use to filter the kinship table to grandparents
    :param grandparent_ibd2: Maximum IBD2 threshold to use to filter the kinship table to grandparents
    :return: Table containing only relationships of grandparents and above
    """
    # Filter to anything above the relationship of a grandparent
    ht = ht.filter(
        (ht.pi_hat > first_degree_pi_hat)
        | (
            (ht.pi_hat > grandparent_pi_hat)
            & (ht.ibd1 > grandparent_ibd1)
            & (ht.ibd2 < grandparent_ibd2)
        )
    )
    ht = ht.annotate(pair=hl.sorted([ht.i, ht.j]))

    out_summary.write(
        f"NOTE: kinship table was filtered to:\n(kin > {first_degree_pi_hat}) or kin > {grandparent_pi_hat} and IBD1 > {grandparent_ibd1} and IBD2 > {grandparent_ibd2})\n"
    )
    out_summary.write(f"relationships not meeting this critera were not evaluated\n\n")

    return ht


def check_sex(sex_ht: hl.Table, output_dir: str, output_name: str,) -> None:

    """
    Compare inferred to given sex and output file with column added for discrepancies.

    Output directory and name here are used to locate the functioning pedigree with given sexes.

    :param sex_ht: Table of inferred sexes for each sample
    :param output_dir: Path to directory to output results
    :param output_name: Output prefix to use for results
    :return: None
    """
    # Read in functioning pedigree with given sexes
    ped_ht = hl.import_table(f"{output_dir}/{output_name}_functioning_pedigree.ped")
    ped_ht = ped_ht.key_by(s=ped_ht.Individual_ID).select("Sex")

    ped_ht = ped_ht.annotate(
        given_sex=hl.case()
        .when(ped_ht.Sex == "M", "XY")
        .when(ped_ht.Sex == "F", "XX")
        .default(ped_ht.Sex)
    ).drop("Sex")

    sex_ht = sex_ht.join(ped_ht, how="outer")
    sex_ht = sex_ht.annotate(discrepant_sex=sex_ht.sex != sex_ht.given_sex)
    sex_ht.export(f"{output_dir}/{output_name}_sex_check.txt")


def write_functional_pedigree(
    input_pedigree: str, vcf_samples: list, output_dir: str, output_name: str
) -> Tuple[dict, dict, dict]:

    """
    Write a functional pedigree (pedigree with samples not in the VCF removed) and create dictionary of seqr projects, family IDs, and given sex.

    :param input_pedigree: Pedigree
    :param vcf_samples: Dictionary of samples found in the VCF
    :return: Dictionary of project IDs for each sample
    """
    seqr_projects = defaultdict(str)
    family_ids = defaultdict(str)
    given_sex = defaultdict(str)

    out_new_ped = hl.hadoop_open(
        f"{output_dir}/{output_name}_functioning_pedigree.ped", "w"
    )
    out_new_ped.write("Family_ID\tIndividual_ID\tPaternal_ID\tMaternal_ID\tSex\n")

    with hl.hadoop_open(input_pedigree, "r") as infile:
        next(infile)
        for line in infile:
            line = line.rstrip("\n")
            items = line.split("\t")
            (
                Project_GUID,
                Family_ID,
                Individual_ID,
                Paternal_ID,
                Maternal_ID,
                Sex,
            ) = items[0:6]

            if Individual_ID not in vcf_samples:
                Individual_ID = "."
            if Paternal_ID not in vcf_samples:
                Paternal_ID = "."
            if Maternal_ID not in vcf_samples:
                Maternal_ID = "."

            # Only output line from pedigree if the proband is not missing
            if Individual_ID != ".":
                seqr_projects[Individual_ID] = Project_GUID
                family_ids[Individual_ID] = Family_ID
                given_sex[Individual_ID] = Sex

                out_new_ped.write(
                    f"{Family_ID}\t{Individual_ID}\t{Paternal_ID}\t{Maternal_ID}\t{Sex}\n"
                )

    out_new_ped.close()

    return seqr_projects, family_ids, given_sex


def main(args):
    output_dir = args.output_dir
    output_name = args.output_name
    inferred_sex = args.inferred_sex
    mt_path = args.mt_path
    input_pedigree = args.input_pedigree

    gnomad_ld = args.gnomad_ld
    run_ibd = args.run_ibd
    first_degree_pi_hat = args.first_degree_pi_hat
    grandparent_pi_hat = args.grandparent_pi_hat
    grandparent_ibd1 = args.grandparent_ibd1
    grandparent_ibd2 = args.grandparent_ibd2
    filter_kinship_ht = args.filter_kinship_ht

    logger.info("Reading in inputs...")
    mt = hl.read_matrix_table(mt_path)
    pedigree = hl.import_table(input_pedigree, impute=True)

    # Infer build of the MatrixTable
    build = get_reference_genome(mt.locus).name

    logger.info("Filtering to biallelic SNVs on autosomes and performing LD pruning...")
    mt = filter_rows_for_qc(
        mt, min_af=0.001, min_callrate=0.99, apply_hard_filters=False
    )
    mt = ld_prune(mt, build, gnomad_ld)
    out_mt = f"{output_dir}/{output_name}_processed_mt.mt"

    logger.info("Remapping sample names...")
    mt, sex_ht = remap_samples(mt_path, mt, pedigree, inferred_sex)

    mt = mt.checkpoint(out_mt, overwrite=True)

    if run_ibd:
        logger.info("Running identity by descent...")
        ibd_results_ht = hl.identity_by_descent(mt, maf=mt.af, min=0.10, max=1.0)
        ibd_results_ht = ibd_results_ht.annotate(
            ibd0=ibd_results_ht.ibd.Z0,
            ibd1=ibd_results_ht.ibd.Z1,
            ibd2=ibd_results_ht.ibd.Z2,
            pi_hat=ibd_results_ht.ibd.PI_HAT,
        ).drop("ibs0", "ibs1", "ibs2", "ibd")
        out_ht = f"{output_dir}/{output_name}_ibd_kinship.tsv"
        ibd_results_ht.export(out_ht)

    else:
        logger.warn("Skipping IBD - using previous calculations...")
        if not file_exists(f"{output_dir}/{output_name}_ibd_kinship.tsv"):
            logger.warning(
                "IBD calculation was skipped but no file with previous calculations was found...",
                sample,
            )

    logger.info("Reading in kinship ht...")
    kin_ht = hl.import_table(f"{output_dir}/{output_name}_ibd_kinship.tsv", impute=True)

    # Subset MatrixTable and sex ht to the samples in the pedigree
    mt_subset, sex_ht, expected_samples, vcf_samples = subset_samples(
        mt, pedigree, sex_ht, output_dir, output_name
    )

    # Subset Table to the samples in the pedigree
    subset = hl.set(expected_samples)
    kin_ht = kin_ht.filter(subset.contains(kin_ht.i) | subset.contains(kin_ht.j))

    # Key the Table
    kin_ht = kin_ht.key_by("i", "j")

    # Setup output file
    out_summary = hl.hadoop_open(
        f"{output_dir}/{output_name}_ped_check_summary.txt", "w"
    )

    if filter_kinship_ht:
        logger.info(
            "Filtering kinship table to remove unrelated individuals from analysis..."
        )
        kin_ht = filter_kin_ht(kin_ht, out_summary, first_degree_pi_hat, grandparent_pi_hat, grandparent_ibd1, grandparent_ibd2)

    # Output basic stats
    out_summary.write(
        "Number individuals in pedigree: " + str(len(expected_samples)) + "\n"
    )
    out_summary.write(
        "Number individuals in subset from the VCF: " + str(len(vcf_samples)) + "\n"
    )
    out_summary.write(
        "Number of relationships in the kinship table: " + str(kin_ht.count()) + "\n\n"
    )
    out_summary.close()

    seqr_projects, family_ids, given_sex = write_functional_pedigree(
        input_pedigree, vcf_samples, output_dir, output_name
    )

    # Compare inferred and given sex
    check_sex(sex_ht, output_dir, output_name)

    kin_ht = add_project_and_family_annotations(kin_ht, seqr_projects, family_ids)

    logger.info("Writing kinship ht per project...")
    # Output original ht per project
    for project in set(seqr_projects.values()):
        full_ht = kin_ht.filter(
            (kin_ht.seqr_proj_i == project) | (kin_ht.seqr_proj_j == project)
        )
        full_ht.drop("seqr_proj_i", "seqr_proj_j").export(
            f"{output_dir}/{project}/{output_name}_{project}_annotated_kin.txt"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script checks a given pedigree against IBD kinship output"
    )
    parser.add_argument(
        "-d", "--output-dir", help="Path to directory to output results"
    )
    parser.add_argument("-n", "--output-name", help="Output prefix to use for results")
    parser.add_argument(
        "-s", "--inferred-sex", help="Path to text file of inferred sex output"
    )
    parser.add_argument("-m", "--mt-path", help="Path to input MatrixTable")
    parser.add_argument("-p", "--input-pedigree", help="Path to input pedigree")
    parser.add_argument(
        "--gnomad-ld", help="Use variants from gnomAD for ld prune", action="store_true"
    )
    parser.add_argument("--run-ibd", help="Run IBD", action="store_true")
    parser.add_argument(
        "--filter-kinship-ht",
        help="Filter the kinship table to relationships of grandparents and above",
        action="store_true",
    )
    parser.add_argument(
        "--first-degree_pi_hat",
        help="Minimum pi_hat threshold to use to filter the kinship table to first degree relatives",
        type=float,
        default=0.40,
    )
    parser.add_argument(
        "--grandparent-pi-hat",
        help="Minimum pi_hat threshold to use to filter the kinship table to grandparents",
        type=float,
        default=0.20,
    )
    parser.add_argument(
        "--grandparent-ibd1",
        help="Minimum IBD1 threshold to use to filter the kinship table to grandparents",
        type=float,
        default=0.25,
    )
    parser.add_argument(
        "--grandparent-ibd2",
        help="Maximum IBD2 threshold to use to filter the kinship table to grandparents",
        type=float,
        default=0.15,
    )

    args = parser.parse_args()

    main(args)
