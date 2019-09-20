#!/usr/bin/env python

import argparse
import logging
from subprocess import Popen, PIPE
import sys


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('make_myoseq_report')
logger.setLevel(logging.INFO)


# define global line break for pdflatex
# the line break is just two (\\); using the extra two to escape the tow needed 
LINE_BREAK = '\\\\'


def awk(sample: str, fname: str) -> str:
    '''
    Quickly looks up one line from a file using awk

    :param str sample: Sample ID to find
    :param str fname: Name of file to search
    :return: Line from file
    :rtype: str
    '''
    cmd = '''awk $1 == "{}"' {}'''.format(sample, fname)
    line, err = Popen([cmd], stdout=PIPE, stderr=PIPE, shell=True).communicate()
    return line


def get_patient_details(proband: str, dirname: str) -> str:
    '''
    Formats proband ID, sex, and ancestry for first page of report

    :param str proband: Proband ID 
    :param str dirname: Name of top level directory with all of sample information for reports
    :return: Formatted string to print out to .tex file
    :rtype: str
    '''

    proband_id = proband.replace('_', '')
    patient_details = '{\\large \\textbf{Proband ID:} ' + f'{proband_id}' + '}' + f'\n{LINE_BREAK}\n'

    sex_check = f'{dirname}/inference/MYOSEQ_sex.tsv'
    sex_line = awk(proband, sex_check)
    sample,reported,inferred,fstat,status = sex_line.replace('_', '').strip().split('\t')
    patient_details += '{\\large \\textbf{Inferred Sex (Reported):} ' + f'{inferred} ({reported})' + '}' + f'\n{LINE_BREAK}\n'

    ancestry_check = f'{dirname}/inference/MYOSEQ_pop.tsv'
    ancestry_line = awk(proband, ancestry_check)
    sample,ancestry = ancestry_line.strip().split('\t')
    patient_details += '{\\large \\textbf{Inferred Ancestry (Reported):} ' + f'{ancestry} (European)' + '}' + f'\n{LINE_BREAK} {LINE_BREAK} {LINE_BREAK} \n' 

    return patient_details


def main(args):

    dirname = args.dirname
    proband = args.proband

    # create string to be written to output tex file
    out_string = ''
    logger.info('Preparing top of first page of reports (ID, sex, ancestry)')
    out_string += get_patient_details(proband, dirname)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Creates patient section of tex file for pdflatex')
    parser.add_argument('-c', '--cnv', help='Samples have CNV candidates', action='store_true', default=False)
    parser.add_argument('-s', '--sma', help='Samples are SMA carriers', action='store_true', default=False)
    parser.add_argument('-d', '--dirname', help='Top level directory that contains all data for reports', required=True)
    parser.add_argument('-p', '--proband', help='Proband ID', required=True)
    parser.add_argument('-o', '--out', help='Output directory', required=True)
    args = parser.parse_args()

    main(args)
