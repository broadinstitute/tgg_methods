import argparse
import hail as hl
import logging
import typing
from gnomad_hail.resources.basics import *


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('get_gnomAD_popmax')
logger.setLevel(logging.INFO)


def import_var(seqr: str) -> hl.Table:
    '''
    Reads in tsv of variants downloaded from seqr into a hail Table.

    :param str seqr: Path to  variants tsv 
    :return: Table of variants
    :rtype: hl.Table
    '''
    ht = hl.import_table(seqr, impute=True)

    # add 'chr' in front of chromosome; hail will not recognize a locus as a valid b38 locus unless the chromosome is prefixed with 'chr'
    ht = ht.transmute(chrom=hl.format('chr%s', ht.chrom))
    
    # create locus and alleles (need these two fields to succesfully join with gnomAD data)
    ht = ht.transmute(
                locus=hl.parse_locus(hl.format('%s:%s', ht.chrom, ht.pos)),
                alleles=[ht.ref, ht.alt]
                 )
    ht = ht.key_by('locus', 'alleles')
    ht.describe()
    return ht


def join_tables(ht: hl.Table, exomes: bool) -> hl.Table:
    '''
    Joins seqr variant table to gnomAD table. NOTE code was written assuming most recent gnomAD release is v3

    :param Table ht: Table with variants downloaded from seqr
    :param bool exomes: Whether to join with gnomAD exomes table or genomes table
    :return: seqr variants Table joined with gnomAD table
    :rtype: hl.Table
    '''
    if exomes:
        # read in exomes table
        gnomad_ht = hl.read_table(get_gnomad_liftover_data_path('exomes', version='2.1.1'))
        gnomad_ht = gnomad_ht.select('freq', 'popmax')
        gnomad_ht = gnomad_ht.select_globals()
        gnomad_ht = gnomad_ht.transmute(
                        gnomad_exomes_AC=gnomad_ht.freq[0].AC,
                        gnomad_exomes_AN=gnomad_ht.freq[0].AN,
                        gnomad_exomes_popmax_AF=gnomad_ht.popmax[0].AF,
                        gnomad_exomes_popmax_pop=gnomad_ht.popmax[0].pop)
        gnomad_ht.describe()
    else:
        # read in genomes table
        gnomad_ht = hl.read_table('gs://gnomad-public/release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht')
        gnomad_ht = gnomad_ht.select('freq')
        gnomad_ht = gnomad_ht.transmute(
                        gnomad_genomes_AC=gnomad_ht.freq[0].AC,
                        gnomad_genomes_AN=gnomad_ht.freq[0].AN)
        gnomad_ht = gnomad_ht.select_globals()
        gnomad_ht.describe()

    ht = ht.annotate(**gnomad_ht[ht.key])
    ht.describe()
    return ht


def add_global_af(ht: hl.Table, temp: str) -> hl.Table:
    '''
    Adds gnomAD global AF annotation to Table

    :param Table ht: Input Table
    :param str temp: Path to temp bucket (to store intermediary files)
    :return: Table with gnomAD global AF annotation
    :rtype: Table
    '''
    # checkpoint table after completing both gnomAD exomes and gnomAD genomes join
    temp_path = f'{temp}/join.ht'
    ht = ht.checkpoint(temp_path)

    # set gnomAD ACs and ANs to 0 if they are missing after the join 
    ht = ht.transmute(
                gnomad_exomes_AC=hl.if_else(hl.is_defined(ht.gnomad_exomes_AC), ht.gnomad_exomes_AC, 0),
                gnomad_genomes_AC=hl.if_else(hl.is_defined(ht.gnomad_genomes_AC), ht.gnomad_genomes_AC, 0),
                gnomad_exomes_AN=hl.if_else(hl.is_defined(ht.gnomad_exomes_AN), ht.gnomad_exomes_AN, 0),
                gnomad_genomes_AN=hl.if_else(hl.is_defined(ht.gnomad_genomes_AN), ht.gnomad_genomes_AN, 0),
                )

    ht = ht.annotate(
    gnomad_global_AF=(hl.if_else(
                      ((ht.gnomad_exomes_AN == 0) & (ht.gnomad_genomes_AN == 0)),
                        0.0,
                        hl.float((ht.gnomad_exomes_AC + ht.gnomad_genomes_AC) / 
                         (ht.gnomad_exomes_AN + ht.gnomad_genomes_AN))
                        )
                     )
    )
    ht.describe()
    return ht


def main(args):

    hl.init(default_reference='GRCh38', log='/get_gnomad_popmax.log')

    logger.info('Importing seqr variants')
    ht = import_var(args.seqr)

    logger.info('Joining seqr table to gnomAD exomes table')
    ht = join_tables(ht, True)

    logger.info('Joining seqr table to gnomAD genomes table')
    ht = join_tables(ht, False)

    logger.info('Adding gnomAD global AF annotation')
    ht = add_global_af(ht, args.temp)

    logger.info('Exporting table to tsv')
    ht.export(args.out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--seqr', help='Path to variants tsv downloaded from seqr', required=True)
    parser.add_argument('-t', '--temp', help='Path to temp bucket', required=True)
    parser.add_argument('-o', '--out', help='Path to output tsv', required=True)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
