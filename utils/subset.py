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


def remap_sample_ids(mt, remap_path):
    """
    Remap the MatrixTable's sample ID, 's', field to the sample ID used within seqr, 'seqr_id'
    If the sample 's' does not have a 'seqr_id' in the remap file, 's' becomes 'seqr_id'
    :param mt: MatrixTable from VCF
    :param remap_path: Path to a file with two columns 's' and 'seqr_id'
    :return: MatrixTable remapped and keyed to use seqr_id
    """
    remap_ht = hl.import_table(remap_path, key='s')
    missing_samples = remap_ht.anti_join(mt.cols()).collect()
    remap_count = remap_ht.count()

    if len(missing_samples) != 0:
        logger.error(
            f'Only {remap_ht.semi_join(mt.cols()).count()} out of {remap_count} '
            'remap IDs matched IDs in the variant callset.\n'
            f'IDs that aren\'t in the callset: {missing_samples}\n'
            f'All callset sample IDs:{mt.s.collect()}', missing_samples
        )

    mt = mt.annotate_cols(**remap_ht[mt.s])
    remap_expr = hl.cond(hl.is_missing(mt.seqr_id), mt.s, mt.seqr_id)
    mt = mt.annotate_cols(seqr_id=remap_expr, vcf_id=mt.s)
    mt = mt.key_cols_by(s=mt.seqr_id)
    logger.info(f'Remapped {remap_count} sample ids...')
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

    logger.info(f"Input MT counts: {mt.count()}")

    if args.mapping:
        logger.info("Mapping sample IDs")
        logger.warning("Assuming mapping file has header!")
        mt = remap_sample_ids(mt, args.mapping)       
 
    logger.info("Subsetting to specified samples and their variants...")
    mt = subset_samples_and_variants(mt, args.samp, header=args.header, table_key=args.table_key)
    logger.info(f"Subset MT counts: {mt.count()}")

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
    parser.add_argument("--mapping", help="Path to file with sample ID mapping")
    parser.add_argument(
        "--header", help="Whether sample list file has a header. Specify only if False", action="store_false"
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
