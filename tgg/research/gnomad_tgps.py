"""
This script checks the number of common (AF > 1%) variants per gnomAD sample.

Eimear Kenny (one of the co-chairs of the Human Pangenome Reference Consortium; HPRC)
requested that we compute the number of common variants in v3 samples that aren't
present in the pangenome reference samples. They want to use these counts (grouped by
genetic ancestry group) to assess the new pangenome reference.

Additional context from Eimear: Karen Miga (cc’d) and I are co-chairing the
population sampling and representation (PSR) working group in the NHGRI-funded [Human
Pangenome Reference Consortium]

The goal of the consortium is to update the current human genome reference by
producing a pangenome of high quality long read sequencing of 350 diverse humans
which will become the new pangenome reference. The PSR group is responsible to
selecting samples to be included in the pangenome, which should comprehensively cover
most common variants, defined as variants >1%, in human populations globally. We are
halfway through (phase 1) of the project, and most of the sampling has been selected
from the 1 000 Genome Cell lines (N=174 samples). We have developed a simple variant
counting metric that helps us evaluate anticipated performance of the pangenome using
short read data in phase one samples and additional out-of-sample genomes from the
1KG project. Here is a link to a brief slide deck that describes this: [
https://www.dropbox.com/s/fo5nrxerxwukw7x/20220908.HPRC-PopulationSampling.Gnomad.pdf
?dl=0]

We would love to additionally be able to assess the performance in gnomAD – not only
is this an independently ascertained dataset, but it would also give us a look ahead
to how aligned the pangenome will be to the top resource used by clinical testing
labs. Ideally, we would like to run the analysis on individual-level data in gnomAD,
and I’m reaching out to see if you would be willing to collaborate with us on this
analysis. The script is a very simple variant count algorithm, so hopefully very
straightforward to run. We are currently writing up this work, and of course,
anybody who participates in this analysis would be included as a co-author.
"""


import argparse
import logging

import hail as hl
import hailtop.fs as hfs
from gnomad.resources.resource_utils import DataException
from gnomad_qc.v3.resources.basics import get_gnomad_v3_vds
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.release import release_sites

NUM_OF_GNOMAD_SAMPLES = 76156
"""
The number of samples included in gnomAD v3.
See: https://gnomad.broadinstitute.org/help/what-populations-are-represented-in-the-
gnomad-data
"""

NUM_OF_PARTITIONS = 10000
"""
The number of partitions to read checkpointed files.
"""
HGDP_SUBSET_MT = (
    "gs://gcp-public-data--gnomad/release/3.1.2/mt/genomes/"
    "gnomad.genomes.v3.1.2.hgdp_1kg_subset_sparse.mt"
)
OUTPUT_FILENAME = "pangenome_assessment_stats.txt"
PANGENOME_CHECKPOINT_FILEPATH = "gs://gnomad-tmp-4day/pangenome/pangenome_checkpoint.ht"
NOT_IN_PANGENOME_CHECKPOINT_FILEPATH = (
    "gs://gnomad-tmp-4day/pangenome/not_in_pangenome_checkpoint.ht"
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Pangenome_assessment")
logger.setLevel(logging.INFO)


def get_pangenome_ids(
    id_path: str,
) -> hl.expr.SetExpression:
    """
    Import Hail set of pangenome sample IDs from specified file path.

    :param id_path: Path to text file (stored in Google bucket) containing
    pangenome sample IDs.
    :return: SetExpression of pangenome sample IDs.
    """
    ids = []
    with hfs.open(id_path, "r") as file:
        for line in file:
            sid = line.split(" ")[1].strip()
            ids.append(sid)
    return hl.set(ids)


def write_output(
    var_stats: hl.expr.StructExpression,
    file_dir: str,
    filename: str,
) -> None:
    """
    Write the output Struct containing computed stats to a .txt file.

    :param var_stats: Data containing computed variant stats.
    :param file_dir: Path to a Google bucket for output file.
    :param filename: Output filename.
    :return: None; function writes output to the output path.
    """
    with hfs.open(f"{file_dir}/{filename}", "w") as f:
        print(var_stats, file=f)


def compute_stats(
    pangenome_ids_path: str,
    output_dir: str,
) -> None:
    """
    Count the number of common (AF > 1%) variants per gnomAD sample.

    Function produces counts only for gnomAD samples that are not part of the pangenome.

    :param pangenome_ids_path: Path to text file (stored in Google bucket)
    containing pangenome sample IDs.
    :param output_dir: Path to a Google bucket for output file.
    :return: None; function writes a text file to output bucket.
    """
    # Get public gnomAD release sites HT.
    gnomad_freq_ht = release_sites().ht()

    # Get gnomAD v3 metadata and filter to released samples.
    meta_ht = meta.ht()
    meta_ht = meta_ht.filter(meta_ht.release == True)

    # Read gnomAD v3 raw data and split multi-allelics.
    mt = get_gnomad_v3_vds().variant_data
    mt = hl.experimental.sparse_split_multi(mt)

    # Use the filtered metadata to filter gnomAD v3 data.
    mt = mt.filter_cols(hl.is_defined(meta_ht[mt.col_key]))
    num_of_v3_samples = mt.count_cols()
    logger.info("The number of samples in gnomAD v3 is: %i", num_of_v3_samples)
    if num_of_v3_samples != NUM_OF_GNOMAD_SAMPLES:
        raise DataException(
            f"Found {num_of_v3_samples} samples in VDS but was expecting"
            f" {NUM_OF_GNOMAD_SAMPLES} samples; please double check!"
        )

    # Get pangenome sample IDs.
    pangenome_ids = get_pangenome_ids(pangenome_ids_path)
    logger.info(
        "The total number of samples in pangenome is %i", hl.eval(hl.len(pangenome_ids))
    )

    # Filter mt to contain gnomAD samples not in pangenome.
    mt = mt.filter_cols(~pangenome_ids.contains(mt.s))

    # Annotate mt with allele frequencies and filter out rare SNVs (AF lower than 1%).
    # Also filter out sites where no samples had a variant call.
    mt = mt.annotate_rows(freq=gnomad_freq_ht[mt.row_key].freq)
    mt = mt.filter_rows(mt.freq.AF[0] > 0.01)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    # Read the HGDP+1KG subset matrix table and split multi-allelics.
    pg_mt = hl.read_matrix_table(HGDP_SUBSET_MT)
    pg_mt = hl.experimental.sparse_split_multi(pg_mt)
    pg_mt = pg_mt.filter_rows(hl.len(pg_mt.alleles) == 1)

    # Filter out sites where no samples had a variant call.
    pg_mt = pg_mt.filter_cols(pangenome_ids.contains(pg_mt.s))
    pg_mt = pg_mt.filter_rows(hl.agg.any(pg_mt.GT.is_non_ref()))

    # Log the numbers of gnomAD samples not in and in the pangenome reference.
    num_in_samples = pg_mt.count_cols()
    num_out_of_samples = mt.count_cols()
    logger.info(
        "The number of samples both in pangenome and gnomAD v3 is: %i", num_in_samples
    )
    logger.info("The number of samples only in gnomAD v3 is: %i", num_out_of_samples)

    # Checkpoint the rows from the pangenome MT if no checkpoint file exists.
    pg_ht = pg_mt.rows()
    if not hfs.is_file(PANGENOME_CHECKPOINT_FILEPATH):
        pg_ht.checkpoint(PANGENOME_CHECKPOINT_FILEPATH)

    # Annotate the gnomAD matrix table to count the number of minor alleles carried
    # by each participant.
    mt = mt.annotate_cols(total_var_per_sample=hl.agg.count_where(mt.GT.is_non_ref()))

    # Filter out SNVs that are in the reference group.
    pg_ht = hl.read_table(
        PANGENOME_CHECKPOINT_FILEPATH, _n_partitions=NUM_OF_PARTITIONS
    )
    logger.info(
        "The number of variants for the pangenome data set is %i", pg_ht.count()
    )
    mt_filtered = mt.anti_join_rows(pg_ht)

    # Annotate the SNV filtered gnomAD dataset to count the number of minor alleles
    # carried by each participant.
    mt_filtered = mt_filtered.annotate_cols(
        n_var_per_sample=hl.agg.count_where(mt_filtered.GT.is_non_ref())
    )

    # Keep only the columns and annotate with population data, both continental and
    # subcontinental ancestry.
    ht = mt_filtered.cols()
    ht = ht.annotate(population=meta_ht[ht.key].population_inference.pop)

    # Checkpoint the ht to perform the computation if no checkpoint file exists.
    if not hfs.is_file(NOT_IN_PANGENOME_CHECKPOINT_FILEPATH):
        ht.checkpoint(NOT_IN_PANGENOME_CHECKPOINT_FILEPATH)

    # Read the computed ht.
    ht = hl.read_table(
        NOT_IN_PANGENOME_CHECKPOINT_FILEPATH, _n_partitions=NUM_OF_PARTITIONS
    )

    # Annotate with fraction of variants carried by each sample that is not included
    # in the pangenome.
    ht = ht.annotate(
        fraction_not_in_pangenome=ht.n_var_per_sample / ht.total_var_per_sample
    )

    # Compute aggregation stats.
    res = ht.aggregate(
        hl.struct(
            mean_var_per_sample_excluding_pangenome_var=hl.agg.group_by(
                ht.population, hl.agg.mean(ht.n_var_per_sample)
            ),
            mean_var_per_sample=hl.agg.group_by(
                ht.population, hl.agg.mean(ht.tn_var_per_sample)
            ),
            mean_fraction_of_var_not_in_pangenome=hl.agg.group_by(
                ht.population, hl.agg.mean(ht.fraction_not_in_pangenome)
            ),
            stats_var_per_sample_excluding_pangenome_var=hl.agg.group_by(
                ht.population, hl.agg.stats(ht.n_var_per_sample)
            ),
            stats_var_per_sample=hl.agg.group_by(
                ht.population, hl.agg.stats(ht.tn_var_per_sample)
            ),
            stats_fraction_of_var_not_in_pangenome=hl.agg.group_by(
                ht.population, hl.agg.stats(ht.fraction_not_in_pangenome)
            ),
        )
    )

    # Write to output.
    write_output(res, output_dir, OUTPUT_FILENAME)


def main(args):
    """
    Calculate the number of variants per gnomAD sample for gnomAD samples not present
    in pangenome.
    """
    logging_path = f"{args.tmp_dir}/logs"
    try:
        # Initialize Hail.
        hl.init(
            log="/pangenome_assessment_counts.log",
            tmp_dir=args.tmp_dir,
            default_reference="GRCh38",
            backend="spark",
        )
        compute_stats(args.pangenome_ids, args.out_dir)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path)


if __name__ == "__main__":
    # Parse input arguments.
    parser = argparse.ArgumentParser(
        "This script calculates the number of common variants per gnomAD sample "
        "excluding pangenom samples."
    )
    parser.add_argument(
        "--tmp-dir",
        help="Temporary directory for the project.",
        default="gs://gnomad-tmp-4day/pangenome",
    )
    parser.add_argument(
        "--pangenome-ids",
        help="Text file containing pangenome sample IDs.",
        default="gs://jialan-tmp-7day/sample_ids_phase2_n55_all.txt",
    )
    parser.add_argument(
        "--out-dir",
        help="Output directory.",
        default="gs://gnomad-tmp-4day/pangenome",
    )
    args = parser.parse_args()

    # Perform calculations and write outputs to the output directory.
    main(args)
