#!/usr/bin/env python


import argparse
import logging
import os
from utils import get_samples


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('bigquery')
logger.setLevel(logging.INFO)


def parse_var(var, samples):
    """
    Gets chrom, pos, ref, alt, sample ID from seqr variant export

    :param str var: Name of file exported from seqr
    :param list samples: List of requested sample IDs
    :return: Dictionary; key: variant, value: list of samples with variant
    :rtype: dict
    """
    variants = {}
    with open(var) as v:

        # get header of seqr export file
        header = v.readline().strip().split('\t')
        # this is to double check that seqr download format has not changed
        logger.info('seqr format check: {}, {}'.format(header[0:6], header[19:23]))
        
        for line in v:
            line = line.strip().split('\t')
            
            # extract necessary information for variant
            # chrom, pos, ref, alt, gene, worst_consequence, rsid, hgvs, hgvps, clinvar, gold stars
            variant_info = line[0:6] + line[19:23]
            variant = '{}-{}-{}-{}'.format(variant_info[0], variant_info[1], variant_info[2], variant_info[3])

            # every field at index 25 and after should have samle level information
            sample_info = line[27:]

            for sample in sample_info:
                if ':' in sample:
                    sample = sample.split(':')
                    if sample[0] in samples:
                        if variant not in variants:
                            variants[variant] = []
                        variants[variant].append(sample[0])
    return variants


def write_bigquery_tsv(variants, out):
    """
    Writes seqr variants to output TSV (to be uploaded into bigquery)

    :param dict variants: Dictionary of variants and samples that have each variant
    :param str out: Name of output TSV
    :return: None
    :rtype: None    
    """
    # open file to write to output
    with open(out, 'w') as o:

        for variant in variants:
            v = variant.split('-')
            samples = str(variants[variant][0])
            if len(variants[variant]) > 1:
                samples = ','.join(variants[variant])
            line = '{}\t{}\t{}\t{}\t{}\t{}\n'.format(v[0], v[1], v[2], v[3], variant, samples)
            o.write(line)


def main(args):

    logger.info('Getting individual IDs from sample list')
    samples = get_samples(args.samp)

    logger.info('Getting variants from seqr')
    variants = parse_var(args.var, samples)
    
    logger.info('Writing variants to output')
    write_bigquery_tsv(variants, args.out)


if __name__ == '__main__':
    
    # Get args from command line
    parser = argparse.ArgumentParser(description='Reformats seqr variant export for biqquery AF lookup')
    parser.add_argument('-v', '--var', help='input (seqr variants)', required=True)
    parser.add_argument('-s', '--samp', help='input sample list', required=True)
    parser.add_argument('-o', '--out', help='output TSV', required=True)
    args = parser.parse_args()

    main(args)
