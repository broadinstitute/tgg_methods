#!/usr/bin/env python

import argparse
import hail as hl
import io
import logging

from collections import defaultdict
from os.path import dirname
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.reference_genome import get_reference_genome
from gnomad_qc.v2.resources.sample_qc import qc_mt_path
from typing import Dict, List, Tuple

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("check pedigree")
logger.setLevel(logging.INFO)

logger.info("Setting hail flag to avoid array index out of bounds error error...")
# Setting this flag isn't generally recommended, but is needed (since at least Hail version 0.2.75) to avoid an array index out of bounds error until changes are made in future versions of Hail
# TODO: reassess if this flag is still needed for future versions of Hail
hl._set_flags(no_whole_stage_codegen="1")


def check_subset(
    input_mt: hl.MatrixTable, pedigree: hl.Pedigree, output_dir: str, output_name: str
) -> Tuple[hl.MatrixTable, list, list]:
    """
    Filter the MatrixTable to biallelics SNVs on the autosomes.

    :param input_mt: MatrixTable
    :param pedigree: Pedigree file
    :return: MatrixTable subsetted to the samples given in the pedigree, list of samples in the pedigree, list of samples in the VCF
    """
    # Get sample names to subset from the pedigree
    samples_to_subset = hl.literal(pedigree.Individual_ID.collect())
    mt_subset = input_mt.filter_cols(samples_to_subset.contains(input_mt["s"]))

    # Filter to variants that have at least one alt call after the subsetting
    mt_subset = mt_subset.filter_rows(hl.agg.any(mt_subset.GT.is_non_ref()))

    # Check that the samples in the pedigree are present in the VCF subset
    out_missing_samples = hl.hadoop_open(
        f"{output_dir}/{output_name}_missing_samples_in_subset.txt", "w"
    )
    expected_samples = pedigree.Individual_ID.collect()
    vcf_samples = mt_subset.s.collect()

    missings = set(expected_samples) - set(vcf_samples)
    for i in missings:
        out_missing_samples.write(i + "\n")
    out_missing_samples.close()

    return (mt_subset, expected_samples, vcf_samples)


def remap_samples(
    original_mt_path: str,
    input_mt: hl.MatrixTable,
    pedigree: hl.Pedigree,
    kin_ht: hl.Table,
    sample_sexes: dict,
) -> Tuple[hl.MatrixTable, hl.Table]:
    """
    Rename `s` col in the MatrixTable and `i` and `j` cols in the kinship ht based on seqr remap files.

    :param original_mt_path: Path to original MatrixTable location
    :param input_mt: MatrixTable 
    :param pedigree: Pedigree file
    :param kin_ht: Kinship table
    :param sample_sexes: Dictionary of sample sexes
    :return: mt, kinship ht, and sex dictionary with sample names remapped
    """
    base_path = "/".join(dirname(original_mt_path).split("/")[:-1]) + ("/base/projects")
    project_list = list(set(pedigree.Project_GUID.collect()))

    # Get the list of hts containing sample remapping information for each project
    remap_hts = []
    for i in project_list:
        remap = f"{base_path}/{i}/{i}_remap.tsv"
        if hl.hadoop_is_file(remap):
            remap_ht = hl.import_table(remap)
            remap_ht = remap_ht.key_by("s", "seqr_id")
            remap_hts.append(remap_ht)

    if len(remap_hts) > 0:
        ht = remap_hts[0]
        for next_ht in remap_hts[1:]:
            ht = ht.join(next_ht, how="outer")

        # Create dictionary for sample name and it's remapped sample name
        rename_dict = dict(hl.tuple([ht.s, ht.seqr_id]).collect())

        # If a sample has a non-missing value for seqr_id, rename it to the sample name
        ht = ht.key_by("s")
        input_mt = input_mt.annotate_cols(seqr_id=ht[input_mt.col_key].seqr_id)
        input_mt = input_mt.key_cols_by(
            s=hl.if_else(hl.is_missing(input_mt.seqr_id), input_mt.s, input_mt.seqr_id)
        )
        kin_ht = kin_ht.annotate(seqr_id_i=ht[kin_ht.i].seqr_id)
        kin_ht = kin_ht.annotate(
            i=hl.if_else(hl.is_missing(kin_ht.seqr_id_i), kin_ht.i, kin_ht.seqr_id_i)
        )
        kin_ht = kin_ht.annotate(seqr_id_j=ht[kin_ht.j].seqr_id)
        kin_ht = kin_ht.annotate(
            j=hl.if_else(hl.is_missing(kin_ht.seqr_id_j), kin_ht.j, kin_ht.seqr_id_j)
        )
        kin_ht = kin_ht.drop("seqr_id_j", "seqr_id_i")

        # Remap samples in the sex dictionary
        sexes_renamed = {}
        for sample, sex in sample_sexes.items():
            if sample in rename_dict:
                new_sample_id = rename_dict[sample]
                sexes_renamed[new_sample_id] = sex
            else:
                sexes_renamed[sample] = sex

    return (input_mt, kin_ht, sexes_renamed)


def filter_to_biallelic_auto_snvs(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter the MatrixTable to biallelics SNVs on the autosomes.

    :param input_mt: MatrixTable
    :return: mt containing only SNV, biallelic variants, high callrate variants, and no rare variants
    """
    # Filter to biallelic SNVs
    input_mt = input_mt.filter_rows(
        (hl.len(input_mt.alleles) == 2)
        & hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])
    )
    input_mt = hl.variant_qc(input_mt)
    # Only take variants with a high callrate
    input_mt = input_mt.filter_rows(input_mt.variant_qc.call_rate > 0.99)

    # Recalculate AF since samples were removed
    input_mt = input_mt.annotate_rows(
        AC=hl.agg.sum(input_mt.GT.n_alt_alleles()),
        AN=2 * hl.agg.count_where(hl.is_defined(input_mt.GT)),
    )
    input_mt = input_mt.annotate_rows(AF=input_mt.AC / input_mt.AN)

    # Remove variants with low AF
    input_mt = input_mt.filter_rows(input_mt.AF > 0.001)

    input_mt = filter_to_autosomes(input_mt)

    return input_mt


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
        pruned_mt = hl.read_matrix_table(qc_mt_path("joint", ld_pruned=True))

        if build == "GRCh38":
            pruned_mt_path = (
                "gs://seqr-kml/pedigree_inference_resources/pruned_mt_lift.mt"
            )
            pruned_mt = hl.read_matrix_table(pruned_mt_path)

        input_mt = input_mt.filter_rows(
            hl.is_defined(pruned_mt.index_rows(input_mt.row_key))
        )

    return input_mt


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


def write_functional_pedigree(
    input_pedigree: str, vcf_samples: list, output_dir: str, output_name: str
) -> dict:

    """
    Write a functional pedigree (pedigree with samples not in the VCF removed) and create dictionary of sample as key and seqr project as value.

    :param input_pedigree: Pedigree
    :param vcf_samples: Dictionary of samples found in the VCF
    :return: Dictionary of project IDs for each sample
    """
    seqr_projects = defaultdict(str)

    out_new_ped = hl.hadoop_open(
        f"{output_dir}/{output_name}_functioning_pedigree.ped", "w"
    )

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

            # Remove line from pedigree if the proband is missing
            if Individual_ID != ".":
                seqr_projects[Individual_ID] = Project_GUID
                out_new_ped.write(
                    "{Family_ID}\t{Individual_ID}\t{Paternal_ID}\t{Maternal_ID}\t{Sex}\n".format(
                        **locals()
                    )
                )

    out_new_ped.close()

    return seqr_projects


def main(args):

    output_dir = args.output_dir
    output_name = args.output_name
    sex_check = args.sex_check
    input_mt = args.input_mt
    input_pedigree = args.input_pedigree

    gnomad_ld = args.gnomad_ld
    run_ibd = args.run_ibd
    first_degree_pi_hat = args.first_degree_pi_hat
    grandparent_pi_hat = args.grandparent_pi_hat
    grandparent_ibd1 = args.grandparent_ibd1
    grandparent_ibd2 = args.grandparent_ibd2

    logger.info("Reading in inputs...")
    mt = hl.read_matrix_table(input_mt)
    pedigree = hl.import_table(input_pedigree, impute=True)

    # Infer build of the mt
    build = get_reference_genome(mt.locus).name

    logger.info("Filtering to biallelic SNVs on autosomes and performing LD pruning...")
    mt = filter_to_biallelic_auto_snvs(mt)
    mt = ld_prune(mt, build, gnomad_ld)
    out_mt = "{output_dir}/{output_name}_processed_mt.mt".format(**locals())

    if run_ibd:
        logger.info("Running identity by descent...")
        ibd_results_ht = hl.identity_by_descent(mt, maf=mt.AF, min=0.10, max=1.0)
        ibd_results_ht = ibd_results_ht.annotate(
            ibd0=ibd_results_ht.ibd.Z0,
            ibd1=ibd_results_ht.ibd.Z1,
            ibd2=ibd_results_ht.ibd.Z2,
            pi_hat=ibd_results_ht.ibd.PI_HAT,
        ).drop("ibs0", "ibs1", "ibs2", "ibd")
        out_ht = f"{output_dir}/{output_name}_ibd_kinship.tsv"
        ibd_results_ht.export(out_ht)

        mt.write(out_mt, overwrite=True)

    else:
        logger.warn("Skipping IBD - using previous calculations...")

    logger.info("Reading in mt...")
    mt = hl.read_matrix_table(
        "{output_dir}/{output_name}_processed_mt.mt".format(**locals())
    )
    mt_path = input_mt

    logger.info("Reading in kinship ht...")
    kin_ht = hl.import_table(
        "{output_dir}/{output_name}_ibd_kinship.tsv".format(**locals()), impute=True
    )

    logger.info(
        "Reading sex check into dictionary..."
    )  # code here is based on old sex check, need new sex check that works on dense data
    sample_sexes = {}
    with hl.hadoop_open(sex_check, "r") as infile:
        next(infile)
        for line in infile:
            line = line.rstrip()
            items = line.split("\t")
            sample = items[0]
            sex = items[6]
            if sex == "female":
                sample_sexes[sample] = True
            elif sex == "male":
                sample_sexes[sample] = False

    logger.info("Remapping sample names...")
    mt, kin_ht, sample_sexes = remap_samples(
        mt_path, mt, pedigree, kin_ht, sample_sexes
    )

    # Subset mt to the samples in the pedigree
    mt_subset, expected_samples, vcf_samples = check_subset(
        mt, pedigree, output_dir, output_name
    )
    num_ped = len(expected_samples)
    num_vcf = len(vcf_samples)

    # Subset ht to the samples in the pedigree
    subset = hl.literal(expected_samples)
    kin_ht = kin_ht.filter(subset.contains(kin_ht.i) | subset.contains(kin_ht.j))

    # Key the ht
    kin_ht = kin_ht.key_by("i", "j")
    # Save the original kin_ht as it will be filtered downstream
    original_kin_ht = kin_ht

    # Setup output file
    out_summary = hl.hadoop_open(
        f"{output_dir}/{output_name}_ped_check_summary.txt", "w"
    )

    logger.info(
        "Filtering kinship table to remove unrelated individuals from analysis..."
    )
    kin_ht = filter_kin_ht(kin_ht, out_summary)

    # Output basic stats
    out_summary.write("Number individuals in pedigree: " + str(num_ped) + "\n")
    out_summary.write(
        "Number individuals in subset from the VCF: " + str(num_vcf) + "\n"
    )
    out_summary.write(
        "Number of relationships in the kinship table: " + str(kin_ht.count()) + "\n\n"
    )
    out_summary.close()

    seqr_projects = write_functional_pedigree(
        input_pedigree, vcf_samples, output_dir, output_name
    )

    # count how many of each discrepancy or confirmed
    hl_seqr_projects = hl.literal(seqr_projects)

    # add annotation for seqr projects of sample i and sample j
    original_kin_ht = original_kin_ht.annotate(
        seqr_proj_i=hl_seqr_projects.get(original_kin_ht.i),
        seqr_proj_j=hl_seqr_projects.get(original_kin_ht.j),
    )
    logger.info("Writing kinship ht per project...")
    # Output original ht per project
    for project in set(seqr_projects.values()):
        project_ht = original_kin_ht.filter(
            (original_kin_ht.seqr_proj_i == project)
            | (original_kin_ht.seqr_proj_j == project)
        )
        project_ht.export(
            "{output_dir}/{project}/{output_name}_{project}_annotated_kin_TEST.txt".format(
                **locals()
            )
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script checks a given pedigree against IBD kinship output"
    )
    parser.add_argument(
        "-d", "--output-dir", help="Path to directory to output results"
    )
    parser.add_argument("-n", "--output-name", help="Output prefix to use for results")
    parser.add_argument("-s", "--sex-check", help="Path to sex check output")
    parser.add_argument("-i", "--input-mt", help="Path to input MatrixTable")
    parser.add_argument("-p", "--input-pedigree", help="Path to input pedigree")
    parser.add_argument(
        "--gnomad-ld", help="Use variants from gnomAD for ld prune", action="store_true"
    )
    parser.add_argument("--run-ibd", help="Run IBD", action="store_true")
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
