import argparse
import logging

import hail as hl

from gnomad.utils.filtering import subset_samples_and_variants


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("subset")
logger.setLevel(logging.INFO)


def remap_sample_ids(mt: hl.MatrixTable, remap_path: str):
    """
    Maps the MatrixTable's sample ID field ('s') to the sample ID used within seqr ('seqr_id').

    If the sample does not have a mapping in the remap file, their ID becomes their seqr ID.

    :param hl.MatrixTable mt: Input MatrixTable.
    :param str remap_path: Path to a file with two columnsL 's' and 'seqr_id'.
    :return: MatrixTable with VCF sample IDs mapped to seqr IDs and keyed with seqr ID.
    :rtype: hl.MatrixTable
    """
    remap_ht = hl.import_table(remap_path, key="s")
    missing_samples = remap_ht.anti_join(mt.cols()).collect()
    remap_count = remap_ht.count()

    if len(missing_samples) != 0:
        logger.error(
            f"Only {remap_ht.semi_join(mt.cols()).count()} out of {remap_count} "
            "remap IDs matched IDs in the variant callset.\n"
            f"IDs that aren't in the callset: {missing_samples}\n"
            f"All callset sample IDs:{mt.s.collect()}",
            missing_samples,
        )

    mt = mt.annotate_cols(**remap_ht[mt.s])
    remap_expr = hl.cond(hl.is_missing(mt.seqr_id), mt.s, mt.seqr_id)
    mt = mt.annotate_cols(seqr_id=remap_expr, vcf_id=mt.s)
    mt = mt.key_cols_by(s=mt.seqr_id)
    logger.info(f"Remapped {remap_count} sample ids...")
    return mt


def main(args):
    """
    Subset joint-called VCF or MT to desired samples.

    Used when CMG/MGRC collaborators request VCFs of their data.
    """
    hl.init(log="/subset.log", default_reference="GRCh38")

    if args.vcf_path:
        logger.info("Importing VCF...")
        logger.warning("Assuming VCF is bgzipped!")
        mt = hl.import_vcf(
            args.vcf_path, force_bgz=True, reference_genome=args.import_build
        )

    else:
        logger.info("Reading in MT...")
        mt = hl.read_matrix_table(args.mt_path)

    logger.info(f"Input MT counts: {mt.count()}")
    mt.describe()

    if args.mapping:
        logger.info("Mapping VCF IDs to seqr IDs...")
        logger.warning("Assuming mapping file the field names s and seqr_id!")
        mt = remap_sample_ids(mt, args.mapping)

    logger.info("Subsetting to specified samples and their variants...")
    mt = subset_samples_and_variants(
        mt, args.sample_list, header=args.no_header, table_key=args.table_key
    )
    logger.info(f"MT counts after subsetting: {mt.count()}")

    logger.info("Exporting VCF...")
    if "bgz" not in args.vcf_out:
        logger.warning(
            "Path to output VCF does not contain '.bgz'; export might be really slow!"
        )
    hl.export_vcf(mt, args.vcf_out, parallel=args.parallel)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        "This script subsets a MatrixTable and exports a VCF"
    )
    parser.add_argument("--vcf-path", help="Path to input VCF")
    parser.add_argument(
        "--import-build",
        help="Reference build to use when importing VCF",
        default="GRCh38",
    )
    parser.add_argument("-m", "--mt-path", help="Path to input MatrixTable")
    parser.add_argument(
        "-s", "--sample-list", help="Path to file with list of sample IDs to subset"
    )
    parser.add_argument("--mapping", help="Path to file with VCF to seqr ID mapping")
    parser.add_argument(
        "--no-header",
        help="Whether sample list file is missing a header. Specify only if False",
        action="store_false",
    )
    parser.add_argument(
        "--table-key", help="Field used to key sample Table", default="s"
    )
    parser.add_argument("-v", "--vcf-out", help="Path to output VCF")
    parser.add_argument(
        "-p",
        "--parallel",
        help="Export sharded VCF (option to pass to hail for parallel output). Default is None for non-parallel output",
        choices=("separate_header", "header_per_shard"),
    )
    args = parser.parse_args()

    if not (args.vcf_path or args.mt_path):
        parser.error(
            "Need to specify at least one input (one of --vcf-path or --mt-path)"
        )

    main(args)
