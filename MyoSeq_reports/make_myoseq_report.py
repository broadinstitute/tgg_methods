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


def cat(fname: str, outname: str) -> None:
    '''
    cats one file into output file.

    :param fname: File name with information to be cat
    :param str outname: Output file name
    :return: None
    :rtype: None
    '''
    cmd = '''cat {} >> {}'''.format(fname, outname)
    Popen([cmd], shell=True).communicate()


def append_out(info: str, outname: str) -> None:
    '''
    Opens output file and appends information

    :param str info: Information to be added into file
    :param str outname: Output file name
    :return: None
    :rtype: None
    '''
    with open(outname, 'a') as o:
        o.write(info)


def get_patient_details(proband: str, dirname: str, outname: str) -> None:
    '''
    Formats proband ID, sex, and ancestry for first page of report

    :param str proband: Proband ID 
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Name of output file
    :return: None
    :rtype: None
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

    append_out(patient_details, outname)


def get_report_variants(proband: str, dirname: str, unsolved: bool, outname: str, resources: str) -> None:
    '''
    Formats any candidate variants (tagged 'REPORT' in seqr) for first page of report

    :param str proband: Proband ID
    :param str dirname: Name of top level directory with all sample information for reports
    :param bool unsolved: Whether patient is unsolved
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :return: None
    :rtype: None
    '''
    if unsolved:
        cat(f'{resources}/myoseq_template_no_candidate_variants_notes.tex', outname)
    else:
        report_file = f'{dirname}/seqr/variants/{proband}.flagged.txt' 

    append_out(report_str, outname)


def main(args):

    proband = args.proband
    dirname = args.dirname
    resources = args.resources
    unsolved = args.unsolved
    outname = f'{args.out}/{proband}.tex'

    logger.info('Setting up the report header')
    cat(f'{resources}/myoseq_template_header.tex', outname)

    logger.info('Preparing top of first page of reports (ID, sex, ancestry)')
    get_patient_details(proband, dirname, outname)

    logger.info('Finishing first page of reports (REPORT genes box)')
    get_report_variants(proband, dirname, unsolved, outname, resources)

    logger.info('Preparing appendix of all variants in gene list')
    cat(f'{resources}/myoseq_template_appendix_gene_list.tex', outname)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Creates patient section of tex file for pdflatex')
    parser.add_argument('-c', '--cnv', help='Samples have CNV candidates', action='store_true')
    parser.add_argument('-s', '--sma', help='Samples are SMA carriers', action='store_true')
    parser.add_argument('-u', '--unsolved', help='Sample is unsolved (no candidates)', action='store_false')
    parser.add_argument('-d', '--dirname', help='Top level directory that contains all data for reports', required=True)
    parser.add_argument('-r', '--resources', help='Directory containing tex files', required=True)
    parser.add_argument('-p', '--proband', help='Proband ID', required=True)
    parser.add_argument('-o', '--out', help='Output directory', required=True)
    args = parser.parse_args()

    main(args)
