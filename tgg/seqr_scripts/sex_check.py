import argparse
import logging
from os.path import dirname
import matplotlib.pyplot as plt

import hail as hl

from gnomad.utils.reference_genome import get_reference_genome


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("sex_check")
logger.setLevel(logging.INFO)


def get_chr_cov(
    mt: hl.MatrixTable,
    build: str,
    chr_name: str,
    af_field: str = "AF",
    call_rate_threshold: float = 0.25,
    af_threshold: float = 0.01,
) -> hl.Table:
    """
    Calculate mean chromosome coverage.

    :param mt: MatrixTable containing samples with chrY variants
    :param build: Reference used, either GRCh37 or GRCh38
    :param chr_name: Chosen chromosome. Must be either autosome (number only) or sex chromosome (X, Y)
    :param af_field: Name of field containing allele frequency information. Default is "AF"
    :param call_rate_threshold: Minimum call rate threshold. Default is 0.25
    :param af_threshold: Minimum allele frequency threshold. Default is 0.01
    :return: Table annotated with mean coverage of specified chromosome
    """

    logger.warning(
        "This function expects chrY to be at index 23 and chrX to be at index 22 in the reference genome contigs list!"
    )

    if chr_name == "Y":
        chr_place = 23
    elif chr_name == "X":
        chr_place = 22
    else:
        try:
            # Chromosome index in '.contigs' list should be one less than the chromosome number
            chr_place = int(chr_name) - 1
        except ValueError:
            logger.error("chr_name cannot be converted to an integer")
            return -99

    chr_name = hl.get_reference(build).contigs[chr_place]

    logger.info(
        "Filtering to chromosome (and filtering to non-par regions if chromosome is X or Y)..."
    )
    sex_mt = hl.filter_intervals(
        mt,
        [hl.parse_locus_interval(chr_name, reference_genome=build)],
    )

    if chr_place == 22:
        sex_mt  = sex_mt.filter_rows(sex_mt.locus.in_x_nonpar())
    if chr_place == 23:
        sex_mt  = sex_mt.filter_rows(sex_mt.locus.in_y_nonpar())    

    # Filter to common SNVs above defined callrate (should only have one index in the array because the MT only contains biallelic variants)
    sex_mt = sex_mt.filter_rows(sex_mt[af_field] > af_threshold)

    # TODO: Make callrate filtering optional before adding code to gnomad_methods
    sex_mt = hl.variant_qc(sex_mt)
    sex_mt = sex_mt.filter_rows(sex_mt.variant_qc.call_rate > call_rate_threshold)

    logger.info("Returning mean coverage on chromosome %s...", chr_name)
    sex_mt = sex_mt.annotate_cols(**{f"{chr_name}_mean_dp": hl.agg.mean(sex_mt.DP)})
    sex_ht = sex_mt.cols()
    return sex_ht


def run_hails_impute_sex(
    mt: hl.MatrixTable,
    build: str,
    outdir: str,
    callset_name: str,
    xy_fstat_threshold: float = 0.75,
    xx_fstat_threshold: float = 0.5,
    aaf_threshold: float = 0.05,
) -> hl.Table:
    """
    Impute sex, annotate MatrixTable with results, and output a histogram of fstat values.

    :param MatrixTable mt: MatrixTable containing samples to be ascertained for sex
    :param build: Reference used, either GRCh37 or GRCh38
    :param callset_name: Basename for callset and output results
    :param outdir: Directory to output results
    :param xy_fstat_threshold: F-stat threshold above which a sample will be called XY. Default is 0.75
    :param xx_fstat_threshold: F-stat threshold below which a sample will be called XX. Default is 0.5
    :param aaf_threshold: Alternate allele frequency threshold for `hl.impute_sex`. Default is 0.05
    :return: Table with imputed sex annotations
    """

    logger.warning(
        "User needs to decide on fstat thresholds for XY/XX by looking at fstat plots!"
    )

    # Filter to the X chromosome and impute sex
    mt = hl.filter_intervals(
        mt,
        [
            hl.parse_locus_interval(
                get_reference_genome(mt.locus).x_contigs[0], reference_genome=build
            )
        ],
    )
    sex_ht = hl.impute_sex(
        mt.GT,
        aaf_threshold=aaf_threshold,
        male_threshold=xy_fstat_threshold,
        female_threshold=xx_fstat_threshold,
    )
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    sex_ht = mt.cols()

    # Plot histogram of fstat values
    df = sex_ht.to_pandas()
    plt.clf()
    plt.hist(df["f_stat"])
    plt.xlabel("Fstat")
    plt.ylabel("Frequency")
    plt.axvline(xy_fstat_threshold, color="blue", linestyle="dashed", linewidth=1)
    plt.axvline(xx_fstat_threshold, color="red", linestyle="dashed", linewidth=1)

    outplot = f"{outdir}/fstat_{callset_name}.png"
    with hl.hadoop_open(outplot, "wb") as out:
        plt.savefig(out)

    return sex_ht


def call_sex(
    callset: str,
    use_y_cov: bool = False,
    y_cov_threshold: float = 0.1,
    normalization_contig: str = "20",
    xy_fstat_threshold: float = 0.75,
    xx_fstat_threshold: float = 0.5,
    aaf_threshold: float = 0.05,
    call_rate_threshold: float = 0.25,
) -> hl.Table:
    """
    Call sex for the samples in a given callset and export results file to the callset directory.

    :param str callset: String of full MatrixTable path for the callset
    :param use_y_cov: Set to True to calculate and use chrY coverage for sex inference. Default is False
    :param y_cov_threshold: Coverage on chrY above which supports male call. Default is 0.1
    :param normalization_contig: Chosen chromosome for calculating normalized coverage. Default is "20"
    :param xy_fstat_threshold: F-stat threshold above which a sample will be called XY. Default is 0.75
    :param xx_fstat_threshold: F-stat threshold below which a sample will be called XX. Default is 0.5
    :param aaf_threshold: Alternate allele frequency threshold for `hl.impute_sex`. Default is 0.05
    :param call_rate_threshold: Minimum required call rate. Default is 0.25
    :return: Table with sex annotations
    """

    # Read in matrix table and define output directory
    # TODO: Generalize before moving into gnomad_methods
    outdir = dirname(callset)
    mt_name = callset.split("/")[-1].strip("\.mt")
    logger.info("Reading matrix table for callset: %s", callset)
    logger.info("Using chromosome Y coverage? %s", use_y_cov)

    mt = hl.read_matrix_table(callset)

    # Filter to SNVs and biallelics
    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
    )

    # Filter to pass variants only (empty set)
    # TODO: Make this an optional argument before moving to gnomad_methods
    mt = mt.filter_rows(mt.filters.length() == 0, keep=True)

    # Infer build:
    build = get_reference_genome(mt.locus).name
    logger.info("Build inferred as %s", build)

    logger.info("Inferring sex...")
    # TODO: Change "female" and "male" to "XX" and "XY"
    if use_y_cov:
        sex_ht = get_chr_cov(mt, "GRCh38", normalization_contig)
        mt = mt.annotate_cols(**sex_ht[mt.col_key])
        sex_ht = get_chr_cov(mt, "GRCh38", "Y")
        mt = mt.annotate_cols(**sex_ht[mt.col_key])
        mt = mt.annotate_cols(
            normalized_y_coverage=hl.or_missing(
                mt[f"chr{normalization_contig}_mean_dp"] > 0,
                mt.chrY_mean_dp / mt[f"chr{normalization_contig}_mean_dp"],
            )
        )
        sex_ht = run_hails_impute_sex(
            mt,
            build,
            outdir,
            mt_name,
            xy_fstat_threshold,
            xx_fstat_threshold,
            aaf_threshold,
        )
        sex_ht = sex_ht.annotate(
            ambiguous_sex=hl.is_missing(sex_ht.is_female),
            sex_aneuploidy=(sex_ht.is_female)
            & hl.is_defined(sex_ht.normalized_y_coverage)
            & (sex_ht.normalized_y_coverage > y_cov_threshold)
            | (~sex_ht.is_female)
            & hl.is_defined(sex_ht.normalized_y_coverage)
            & (sex_ht.normalized_y_coverage < y_cov_threshold),
        )

        sex_expr = (
            hl.case()
            .when(sex_ht.ambiguous_sex, "ambiguous_sex")
            .when(sex_ht.sex_aneuploidy, "sex_aneuploidy")
            .when(sex_ht.is_female, "female")
            .default("male")
        )

        sex_ht = sex_ht.annotate(sex=sex_expr)
        sex_ht = sex_ht.select(
            sex_ht.is_female,
            sex_ht.f_stat,
            sex_ht.n_called,
            sex_ht.expected_homs,
            sex_ht.observed_homs,
            sex_ht.sex,
            sex_ht.chrY_mean_dp,
            sex_ht.chr20_mean_dp,
            sex_ht.normalized_y_coverage,
        )

    else:
        sex_ht = run_hails_impute_sex(
            mt,
            build,
            outdir,
            mt_name,
            xy_fstat_threshold,
            xx_fstat_threshold,
            aaf_threshold,
        )
        sex_ht = sex_ht.annotate(ambiguous_sex=hl.is_missing(sex_ht.is_female))
        sex_expr = hl.if_else(
            sex_ht.ambiguous_sex,
            "ambiguous_sex",
            hl.if_else(sex_ht.is_female, "female", "male"),
        )
        sex_ht = sex_ht.annotate(sex=sex_expr)
        sex_ht = sex_ht.select(
            sex_ht.is_female,
            sex_ht.f_stat,
            sex_ht.n_called,
            sex_ht.expected_homs,
            sex_ht.observed_homs,
            sex_ht.sex,
        )

    outfile = f"{outdir}/sex_{mt_name}.txt"
    sex_ht.export(outfile)
    return sex_ht


def main(args):
    """
    Call `call_sex` to infer sex of samples in input MatrixTable.

    :param args: User's command line inputs
    """

    call_sex(**vars(args))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script infers the sex of samples"
    )
    parser.add_argument(
        "-i", "--callset", required=True, help="Path to Callset MatrixTable"
    )
    parser.add_argument(
        "-u",
        "--use-y-cov",
        help="Whether to use chromosome Y coverage when inferring sex. Note that Y coverage is required to infer sex aneuploidies",
        action="store_true",
    )
    parser.add_argument(
        "-y",
        "--y-cov-threshold",
        help="Y coverage threshold used to infer sex aneuploidies (XY samples below and XX samples above this threshold will be inferred as having aneuploidies)",
        default=0.1,
    )
    parser.add_argument(
        "-m",
        "--xy-fstat-threshold",
        help="F-stat threshold above which a sample will be called XY. Default is 0.75",
        default=0.75,
    )
    parser.add_argument(
        "-f",
        "--xx-fstat-threshold",
        help="F-stat threshold below which a sample will be called XX. Default is 0.5",
        default=0.50,
    )
    parser.add_argument(
        "-a",
        "--aaf-threshold",
        help="Alternate allele frequency threshold for `hl.impute_sex`. Default is 0.05",
        default=0.05,
    )
    parser.add_argument(
        "-c",
        "--call-rate-threshold",
        help="Minimum variant call rate threshold. Default is 0.25",
        default=0.25,
    )
    # NOTE: this is an integer here because get_chr_cov expects chromosomes to be specified using their numbers
    parser.add_argument(
        "-n",
        "--normalization-contig",
        help="Autosome to use to normalize sex chromosome coverage. Default is chromosome 20",
        default="20",
    )

    args = parser.parse_args()
    main(args)
