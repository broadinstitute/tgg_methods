import argparse
import logging
import hail as hl


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("subset")
logger.setLevel(logging.INFO)


# Adapted from hail-elasticsearch-pipelines
def subset_samples_and_variants(
    mt: hl.MatrixTable, sample_path: str, header: bool = True, table_key: str = "s"
) -> hl.MatrixTable:
    """
    Subset the MatrixTable to the provided list of samples and their variants
    
    :param MatrixTable mt: Input MatrixTable
    :param str sample_path: Path to a file with list of samples
    :param bool header: Whether file with samples has a header. Default is True
    :param str table_key: Key to sample Table. Default is "s"
    :return: MatrixTable subsetted to list of samples and their variants
    :rtype: hl.Matrixtable
    """
    sample_ht = hl.import_table(sample_path, no_header=not header, key=table_key)
    sample_count = sample_ht.count()
    anti_join_ht = sample_ht.anti_join(mt.cols())
    anti_join_ht_count = anti_join_ht.count()
    full_count = mt.count_cols()

    if anti_join_ht_count != 0:
        missing_samples = anti_join_ht.s.collect()
        raise MatrixTableSampleSetError(
            f"Only {sample_count - anti_join_ht_count} out of {sample_count} "
            "subsetting-table IDs matched IDs in the variant callset.\n"
            f'IDs that aren"t in the callset: {missing_samples}\n'
            f"All callset sample IDs:{mt.s.collect()}",
            missing_samples,
        )

    mt = mt.semi_join_cols(sample_ht)
    mt = mt.filter_rows((hl.agg.any(mt.GT.is_non_ref())) > 0)

    logger.info(
        f"Finished subsetting samples. Kept {mt.count_cols()} "
        f"out of {full_count} samples in MT"
    )
    return mt


def main(args):

    hl.init(default_reference="GRCh38", log="/subset.log")

    if args.vcf_path:
        logger.info("Importing VCF...")
        # Note: always assumes file is bgzipped
        mt = hl.import_vcf(
            args.vcf_path, force_bgz=True, reference_genome=args.import_build
        )

    else:
        logger.info("Reading in MT...")
        mt = hl.read_matrix_table(args.mt_path)

    logger.info("Subsetting to specified samples and their variants...")
    mt = subset_samples_and_variants(mt, args.samp, args.header, args.table_key)

    logger.info("Exporting VCF...")
    if "bgz" not in args.vcf_out:
        logger.warning(
            "Path to output VCF does not contain '.bgz'; export might be really slow"
        )
    hl.export_vcf(mt, args.vcf_out, parallel=args.parallel)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        "This script subsets a MatrixTable and exports a VCF"
    )
    parser.add_argument("--vcf_path", help="Path to input VCF")
    parser.add_argument(
        "--import_build",
        help="Reference build to use when importing VCF",
        default="GRCh38",
    )
    parser.add_argument("-m", "--mt_path", help="Path to input MatrixTable")
    parser.add_argument(
        "-s", "--samp", help="Path to file with list of samples to subset"
    )
    parser.add_argument(
        "--header", help="Whether sample list file has a header", action="store_false"
    )
    parser.add_argument(
        "--table_key", help="Field used to key sample Table", default="s"
    )
    parser.add_argument("-v", "--vcf_out", help="Path to output VCF")
    parser.add_argument(
        "-p",
        "--parallel",
        help="Export sharded VCF (option to pass to hail for parallel output). Default is None for non-parallel output",
        choices=("separate_header", "header_per_shard"),
    )
    args = parser.parse_args()

    if not (args.vcf_path or args.mt_path):
        parser.error(
            "Need to specify at least one input (one of --vcf_path or --mt_path)"
        )

    main(args)
