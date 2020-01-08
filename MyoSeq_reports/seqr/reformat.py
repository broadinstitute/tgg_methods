#!/usr/bin/env python


import argparse
import logging
import os
import typing
from utils import get_samples


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('bigquery')
logger.setLevel(logging.INFO)



def fix_blanks(tsv, out) -> None:
    '''
    Fixes blanks in tsv (for import into hail)

    :param str tsv: Path to file downloadedfrom seqr
    :param str out: Output path
    :return: None
    :rtype: None
    '''
    with open(tsv) as t, open(out, 'w') as o:

        # get header of seqr export file
        header = t.readline()
        o.write(header)
        header = header.strip().split('\t')

        for line in t:
            line = line.strip().split('\t')

            # skip malformed lines
            if len(line) != len(header):
                logger.warning (f'Skipping the following line (length differs from header length): \n{line}')

            for i in range(len(line)):
                if line[i] == '':
                    line[i] = '.'

            o.write('\t'.join(line) + '\n')


def main(args):

    logger.info('Fixing input file')
    fix_blanks(args.tsv, args.out)


if __name__ == '__main__':
    
    # Get args from command line
    parser = argparse.ArgumentParser(description='Fixes blanks in seqr downloads for import into hail')
    parser.add_argument('-t', '--tsv', help='TSV downloaded from seqr', required=True)
    parser.add_argument('-o', '--out', help='output TSV', required=True)
    args = parser.parse_args()

    main(args)
