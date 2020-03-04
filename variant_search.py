import logging
import argparse

import hail as hl
import hail.expr.aggregators as agg
from gnomad_hail import *
from gnomad_hail.utils import *

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('seqr_sample_qc')
logger.setLevel(logging.INFO)


def filter_to_locus(mt: hl.MatrixTable, vcf_path: str,  project: str, locus: str, build: str):
    """
    This module retrieves the specified variant locus.
    :param hl.MatrixTable mt: MatrixTable to filter
    :param str vcf_path: path to the vcf for writing purposes
    :param str project: Project name for writing purpose
    :param str locus: locus to filter rows
    :param str build: 37 or 38
    :return:
    """
    mt = hl.filter_intervals(r_mt, [hl.parse_locus_interval(locus, 'GRCh38')])
    mt = mt.annotate_cols(variant_present=hl.agg.count_where(mt.GT.is_non_ref()))
    mt = mt.filter_cols(mt.variant_present > 0)
    variants, samples = mt.count()
    logger.info(f'{samples} samples have this variant')
    mt.s.show()
    mt.GT.show()


def main(args):
    """
    This module subsets a vcf to the supplied list of samples and retrieves the specified variant locus. The list of
    samples containing a variant at the supplied locus is written out as a tsv. Either step can be skipped using the
    appropriate argument
    """
    mt_path = args.mt_path

    mt = hl.read_matrix_table(f'{mt_path}')  # TODO to avoid writing out vcf again, read this back in if present'''
    filter_to_locus(mt, vcf_path, project, args.variant_locus, build)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--mt-path', help='Path to mt')
    parser.add_argument('--variant-locus', help='Locus of the variant, chrom:position, i.e. 19:47259533')#TODO Replace with file which will be read in and parsed in the code to avoid needs more flags for alleles
    parser.add_argument('--slack-channel', help='Slack channel to post results and notifications to.')


    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
