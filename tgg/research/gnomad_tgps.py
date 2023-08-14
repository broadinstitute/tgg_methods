"""This script checks the number of common (AF > 1%) variants per gnomAD sample. The
gnomAD samples used are not present in the pangenome reference samples.

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
anybody who participates in this analysis would be included as a co-author."""

import hail as hl
from gnomad_qc.v3.resources.basics import get_gnomad_v3_vds
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.meta import meta
import hailtop.fs as hfs
import argparse
import logging

NUM_OF_GNOMAD_SAMPLES = 76156
OUTPUT_FILENAME = "pangemoe_assessment_stats.txt"

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Pangenome_assessment")
logger.setLevel(logging.INFO)


def get_pangenome_ids(
        id_path: str,
) -> hl.expr.expressions.typed_expressions.SetExpression:
    """
    Generate a Hail set of sample IDs in the pangenome.

    :param str id_path: Cloud path to read in the .txt file containing pangenome sample
    IDs.
    :return: A set of sample IDs as hl.expr.expressions.typed_expressions.SetExpression.
    """
    ids = []
    with hfs.open(id_path, 'r') as file:
        for line in file:
            sid = line.split(' ')[1].strip()
            ids.append(sid)
    ids = hl.set(ids)
    return ids


def write_output(
        dat: hl.utils.struct.Struct,
        file_dir: str,
        filename: str,
) -> None:
    """
    Write the output Struct containing computed stats to a .txt file.

    :param hl.utils.struct.Struct dat: Data containing computed stats.
    :param str file_dir: Output directory.
    :param str filename: Output filename.
    :return: None; function writes output to the output path.
    """
    with hfs.open(f'{file_dir}/{filename}', 'w') as f:
        print(dat, file=f)


def compute_stats(
        pangenome_ids_path: str,
        output_dir: str,
) -> None:
    """
    Compute the stats of variants per gnomAD sample for gnomAD samples not in the
    pangenome.

    :param str pangenome_ids_path: Cloud path for the text file containing pangenome
    sample IDs.
    :param str output_dir: Cloud path to output the result file.
    :return: None; function writes .txt file to output directory.
    """
    # Get variant site frequency table.
    gnomad_freq_ht = release_sites().ht()

    # Get gnomAD v3 data in split form.
    mt = get_gnomad_v3_vds().variant_data
    mt = hl.experimental.sparse_split_multi(mt)

    # Get gnomAD v3 metadata and filter to released samples.
    meta_ht = meta.ht()
    meta_ht = meta_ht.filter(meta_ht.release == True)

    # Use the filtered metadata to filter gnomAD v3 data.
    mt = mt.filter_cols(hl.is_defined(meta_ht[mt.col_key]))
    num_of_v3_samples = mt.count_cols()
    logger.info("The number of samples in gnomAD v3 is: ", num_of_v3_samples)

    # Get pangenome sample IDs.
    pangenome_ids = get_pangenome_ids(pangenome_ids_path)
    logger.info("The total number of samples in pangenome is ", hl.eval(hl.len(
        pangenome_ids)))

    # Filter mt to two new matrix tables: one containing gnomAD samples not in
    # pangenome, and one containing samples in pangenome.
    mt = mt.annotate_cols(s_in_pangenome=pangenome_ids.contains(mt.s))
    pg_mt = mt.filter_cols(mt.s_in_pangenome)
    mt = mt.filter_cols(~mt.s_in_pangenome)

    # Log the numbers of gnomAD samples not in and in the pangenome reference.
    num_in_samples = pg_mt.count_cols()
    num_out_of_samples = mt.count_cols()
    logger.info("The number of samples both in pangenome and gnomAD v3 is: ",
                num_in_samples)
    logger.info("The number of samples only in gnomAD v3 is: ", num_out_of_samples)

    # Annotate the gnomAD matrix table with allele frequencies and
    # filter out SNVs with minor allele frequencies lower than 1% and
    # filter out samples with no variant call.
    mt = mt.annotate_rows(freq=gnomad_freq_ht[mt.row_key].freq)
    mt = mt.filter_rows(mt.freq.AF[1] > 0.01)  # ##why AF[1]?
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    # Filter out samples with no variant call in the pangenome matrix table.
    pg_mt = pg_mt.filter_rows(hl.agg.any(pg_mt.GT.is_non_ref()))

    # Checkpoint the rows from the pangenome MT.
    # Issue: the job is still quite large with 115376 partitions.
    # logger.info("The number of variants for the pangenome data set is ",
    # pg_mt.count_rows())

    # Annotate the gnomAD matrix table to count the number of minor alleles carried
    # by each participant.
    mt = mt.annotate_cols(total_var_per_sample=hl.agg.count_where(mt.GT.is_non_ref()))

    # Filter out SNVs that are in the reference group.
    mt_filtered = mt.anti_join_rows(pg_mt.rows())

    # Annotate the SNV filtred gnomAD dataset to count the number of minor alleles
    # carried by each participant.
    mt_filtered = mt_filtered.annotate_cols(
        n_var_per_sample=hl.agg.count_where(mt_filtered.GT.is_non_ref())
    )

    # Keep only the columns and annotate with population data, both continental and
    # subcontinental ancestry.
    ht = mt_filtered.cols()
    ht = ht.annotate(population=meta_ht[ht.key].population_inference.pop)

    # Annotate with fraction of variants carried by each sample that is not included
    # in the pangenome.
    ht = ht.annotate(fraction=ht.n_var_per_sample / ht.total_var_per_sample)

    # Compute aggregation stats.
    res = ht.aggregate(
            hl.struct(
                ght=hl.agg.group_by(ht.population,
                                    hl.agg.mean(ht.n_var_per_sample)),
                ghtt=hl.agg.group_by(ht.population,
                                     hl.agg.mean(ht.tn_var_per_sample)),
                ghtf=hl.agg.group_by(ht.population,
                                     hl.agg.mean(ht.fraction)),
                ght_stats=hl.agg.group_by(ht.population,
                                          hl.agg.stats(ht.n_var_per_sample)),
                ghtt_stats=hl.agg.group_by(ht.population,
                                           hl.agg.stats(ht.fraction)),
                ghtf_stats=hl.agg.group_by(ht.population,
                                           hl.agg.stats(ht.tn_var_per_sample))
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
        hl.init(log="/pangenome_assessment_counts.log",
                tmp_dir=args.tmp_dir,
                default_reference="GRCh38",
                backend="spark")
        compute_stats(args.pangenome_ids, args.out_dir)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path)


if __name__ == '__main__':
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
