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


def main(args):

    hl.init(default_reference="GRCh38", log="/import_vcf.log", tmp_dir='hdfs:///import_vcf.tmp/')

    logger.info("Importing filters info...")
    filters_ht = hl.import_vcf(
        args.part_two_path,
        reference_genome='GRCh38', 
        force_bgz=True, 
        find_replace=('nul','.')
    ).rows()
    logger.info("Importing VCF...")
    # Note: always assumes file is bgzipped
    mt = hl.import_vcf(
        args.part_one_path, 
        force_bgz=True, 
        reference_genome=args.import_build,
        find_replace=('nul','.'),
    )
    mt = mt.annotate_rows(filters=filters_ht[mt.row_key].filters)
    logger.info(f'MT count: {mt.count()}')
    hl.export_vcf(mt, args.vcf_out, parallel='header_per_shard')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        "This script subsets a MatrixTable and exports a VCF"
    )
    parser.add_argument("--part_two_path", help="Path to part one outputs from callset (contain filter information)")
    parser.add_argument("--part_one_path", help="Path to input VCFs")
    parser.add_argument(
        "--import_build",
        help="Reference build to use when importing VCF",
        default="GRCh38",
    )
    parser.add_argument("--vcf_out", help="Output path for VCF")
    args = parser.parse_args()

    main(args)
