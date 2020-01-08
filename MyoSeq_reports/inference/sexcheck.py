#!/usr/bin/env python

import argparse
import logging
from utils import get_samples, check_missing_samples

logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('sexcheck')
logger.setLevel(logging.INFO)


def ped_sex(ped, samples):
    """
    Gets reported sex for all samples from pedigree file

    :param str ped: Name of pedigree file
    :param list samples: List of requested sample IDs
    :return: Dictionary; key: sample, value: reported sex
    :rtype: dict
    """
    reported = {}
    with open(ped) as p:
        for line in p:
            line = line.strip().split('\t')
            fid, sid, did, mid, sex, aff = line[0:6]
            if sid in samples:
                reported[sid] = sex

    # check if any samples are missing and print missing samples to stdout
    check_missing_samples(samples, reported.keys(), 'ped')

    return reported


def compare(infer, samples, reported, out):
    """
    Gets inferred sex for all samples, compares to reported sex, and writes to output

    :param str infer: Name of file with inferred sex
    :param list samples: List of requested sample IDs
    :param dict reported: Dictionary of sample IDs and reported sex
    :param str out: Name of output TSV
    :return: None
    :rtype: None
    """
    # get inferred sex for all samples of interest
    # s is_female   f_stat  n_called    expected_homs   observed_homs   sex y_cov   twenty_cov  normalized_y_coverage
    with open(infer) as i, open(out, 'w') as o:
        found = []
        for line in i:
            s, is_female, f_stat, n_called, e_homs, o_homs, sex, y_cov, twenty_cov, norm_y_cov  = line.strip().split('\t')

            if s in samples:
                found.append(s)
                inferred = sex.capitalize()

                # check for conflicts between reported and inferred sex and write to output
                if inferred == reported[s]:
                    o.write('{}\t{}\t{}\t{}\t{}\n'.format(s, reported[s], inferred, f_stat, 'CONCORD'))
                else:
                    o.write('{}\t{}\t{}\t{}\t{}\n'.format(s, reported[s], inferred, f_stat, 'CONFLICT'))

    # check if any samples are missing from inferred file
    check_missing_samples(samples, found, 'inferred sex')
                

def main(args):

    logger.info('Getting samples from file')
    samples = get_samples(args.samp)

    logger.info('Getting reported sex for each sample')
    reported = ped_sex(args.ped, samples)

    logger.info('Comparing inferred and reported sex and writing to output')
    compare(args.infer, samples, reported, args.out)


if __name__ == '__main__':
    
    # Get args from command line
    parser = argparse.ArgumentParser(description='Checks inferred sex against reported sex for MYOSEQ reports')
    parser.add_argument('-p', '--ped', help='ped file', required=True)
    parser.add_argument('-s', '--samp', help='input list of samples', required=True)
    parser.add_argument('-i', '--infer', help='input list of inferred sex', required=True)
    parser.add_argument('-o', '--out', help='output TSV', required=True)
    args = parser.parse_args()

    main(args)
