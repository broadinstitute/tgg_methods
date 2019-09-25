#!/usr/bin/env python

import argparse
import logging
import pandas as pd
from subprocess import Popen, PIPE
import sys


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('make_myoseq_report')
logger.setLevel(logging.INFO)


# define globals for pdflatex
# the line break is just two (\\); using the extra two to escape the tow needed 
LINE_BREAK = '\\\\'
END_BOX = '\\end{tabular}\n\\end{small}\n'


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


def get_patient_details(proband: str, dirname: str, outname: str) -> str:
    '''
    Formats proband ID, sex, and ancestry for first page of report

    :param str proband: Proband ID 
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Name of output file
    :return: Sample sex
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

    append_out(patient_details, outname)
    return inferred


def report_cnvs(proband, dirname, outname, resources):
    '''
    Formats any candidate CNVs for first page of report

    :param str proband: Proband ID
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :return: None
    :rtype: None
    '''
    cat(f'{resources}/myoseq_template_report_cnv_table_header.tex', outname)
    
    # add empty line to top of CNV box for formatting
    cnv_info = f'{} & {} & {} & {} & {} {LINE_BREAK}\n'

    logger.info('Extracting information from CNV file')
    cnv_file = f'{dirname}/summary/{proband}_CNV.tsv'
    with open(cnv_file) as c:
        # previous format:
        # CAPN3	Loss (CN=1)	WES	15:42676429-42686791	10.362
        # new format (from gCNV output):
        # chr   start   end CN  gene
        for line in c:
            chrom, start, end, cn, gene = line.strip().split('\t')
            cnv = f'{chrom}:{start}-{end}'
            # size is reported in kb 
            size = float(start - end) / 1000
            if int(cn) <= 1:
                cnv_type = f'Loss (CN={cn})'
            else:
                cnv_type = f'Gain (CN={cn})'
            cnv_info += f'{gene} & ' + '\\textbf{\\color{red}' + f'{cnv_type}' + '} & ' + f'WES & {cnv} & {size} {LINE_BREAK} \\hline \n'
    cnv_info += f'{END_BOX}'
    append_out(cnv_info, outname)

    logger.info('Closing out CNV section of first page of reports')
    cat(f'{resources}/myoseq_template_candidate_cnv_notes.tex')
            

def report_sma(proband, dirname, outname, resources):
    '''
    Formats SMA information for first page of report

    :param str proband: Proband ID
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :return: None
    :rtype: None
    '''
    cat(f'{resources}/myoseq_template_report_sma_table_header.tex', outname)

    logger.info('Creating SMA box')    
    sma_info = f'{} & {} & {} {LINE_BREAK}\n'
    sma_info += 'SMN1 & \\textbf{\\color{red}Loss (CN=0)} & WES ' + f'{LINE_BREAK}' + ' \\hline\n'
    sma_info += f'{END_BOX}'
    append_out(sma_info, outname)

    logger.info('Closing out SMA section of first page of reports')
    cat(f'{resources}/myoseq_template_candidate_cnv_notes.tex')


def get_report_variants(proband: str, dirname: str, unsolved: bool, outname: str, resources: str, sex: str) -> None:
    '''
    Formats any candidate SNVs/indels (tagged 'REPORT' in seqr) for first page of report

    :param str proband: Proband ID
    :param str dirname: Name of top level directory with all sample information for reports
    :param bool unsolved: Whether patient is unsolved
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :param str sex: Sample's inferred sex
    :return: None
    :rtype: None
    '''
    if unsolved:
        logger.info('Adding None to REPORT genes box')
        cat(f'{resources}/myoseq_template_no_candidate_variants_notes.tex', outname)
    else:
        # set up pdflatex formatting for special characters
        gt = '\\textgreater'
        underscore = '\\textunderscore'

        logger.info('Starting REPORT genes box')
        cat(f'{resources}/myoseq_template_report_variants_table_header.tex', outname)
        report_file = f'{dirname}/seqr/variants/{proband}.flagged.txt'
        report_info = pandas.read_csv(report_file, sep='\t').to_dict('records')

        # break variant into component chr:pos, ref, alt and get rest of variant information
        #'variant': 'chr1:236883444 T>A',
        chr_pos = report_info[0]['variant'].split(' ')[0]
        ref, alt = report_info[0]['variant'].split(' ')[1].split('>')
        genotype = report_info[0]['genotype']
        transcript, hgvsc = report_info[0]['hgvsc'].replace('>', f'{gt}').replace('_', f'{underscore}').split(':')
        hgvsp = report_info[0]['hgvsp'].split(':')[1].replace('>', f'{gt}').replace('_', f'{underscore}')
        function = report_info[0]['functional_class']
        global_af = report_info[0]['gnomad_global_af']
        popmax_af = report_info[0]['gnomad_pop_max_af']
        popmax_pop = report_info[0]['gnomad_pop_max_population']
        stars = report_info[0]['number_of_stars']

        # truncate long indels
        indel_fmt = '\\textsuperscript{+}'
        if len(ref) > 10:
            ref = ref[:10] + f'{indel_fmt}'
        if len(alt) > 10:
            alt = alt[:10] + f'{indel_fmt}'
        if len(hgvsc) > 30:
            hgvsc = hgsvsc[:30] + f'{indel_fmt}'

        # finish formatting variant
        variant = f'{chr_pos} {ref} {gt} {alt}'
        variant_info = '\\bf{' + f"{report_info[0]['gene_name']}" + '} & ' + f'{variant} & {round(global_af, 2)} &'

        # correct hom to hem for males on X
        if 'X' in chr_pos and genotype == 'hom' and sex == 'Male':
            genotype = 'hem'

        # parse ClinVar status and determine color
        # Monkol split on "/" to split up entries that are too long (pathogenic/likely pathogenic)
        clinsig = report_info[0]['clinvar_clinsig']).split('/')
        revstat = report_info[0]['clinvar_clnrevstat']
        if clinsig == ['.']:
            variant_info += f' NA & {comment_foot} {LINE_BREAK}\n'
        else:
            if clinsig.contains('pathogenic'):
                color = 'red'
            elif clinsig.contains('conflicting'):
                color = 'red'
            elif clinsig.contains('uncertain'):
                color = 'Orange'
            elif clinsig.contains('benign'):
                color = 'Green'
            # other clinvar statuses default to black text
            else:
                color = ''

            if color != '': 
                for i in len(clinsig):
                    variant_info += '\\textcolor{' + f'{color}' + '}{\\textbf{' + f'{clinsig[i]}' + '}} & ' + f'{comment_foot} {LINE_BREAK}\n'
            else:
                for i in len(clinsig):
                    variant_info += '{\\textbf{' + f'{clinsig[i]}' + '}} & ' + f'{comment_foot} {LINE_BREAK}\n'

        # color variant's hgvsc red if it is a splice variant; otherwise display in black
        if function == 'splice_donor_variant' or function == 'slice_acceptor_variant':
            variant_info += f'{transcript} & ' + '\\textcolor{red}{' + f'{hgvsc}' + '} & ' + f'{round(popmax_af, 2)} & ' + f'{revstat} & {LINE_BREAK}\n'
        else:
            variant_info += f'{transcript} & {hgvsc} & {round(popmax_af, 2)} & {revstat} & {LINE_BREAK}\n'
        
        # color variant's hgvsp  according to functional impact
        if len(hgvsp) > 0:
            color = ''
            if function == 'synonymous_variant':
                color = 'Green'
            elif function == 'frameshift_variant' or function == 'stop_gained' or function == 'splice_donor_variant' or function == 'slice_acceptor_variant':
                color = 'red'
            elif function == 'missense_variant' or function == 'inframe_deletion' or function == 'inframe_insertion':
                color = 'Orange'
            else:
                color = ''

            if color != '':
                hgvsp =  '\\textcolor{' + f'{color}' + '}{' + f'{hgvps}' + '}'
        variant_info += '\\textit{' + f'{genotype}' + '} & ' + f'{hgvps} & & {stars} & {LINE_BREAK}\n'
        variant_info += f'& {report_info[0]["rsid"]} & & & {LINE_BREAK}\n\\hline\n'
        variant_info += f'{END_BOX}'      

        # format notes from seqr
        # create footer for notes made in seqr
        # Monkol only accounts for up to 6 notes (a-g), so I've done the same
        alpha_foot = ['a', 'b', 'c', 'd', 'e', 'f', 'g']

        comments = report_info[0]['comments']
        if comments != '.':
            comments = comments.split('|')
            # 	kchao@broadinstitute.org: https://www.ncbi.nlm.nih.gov/pubmed/28554942|kchao@broadinstitute.org: https://www.ncbi.nlm.nih.gov/pubmed/28544275|kchao@broadinstitute.org: not seen in gnomAD
            # truncating number of comments at 6, per Monkol's code
            for i in range(7):
                comment = comments[i].split(': ')[1].replace('http', '\\\\\nhttp')
                variant_info += '\\textsuperscript{' + f'{alpha_foot[i]}' + '} ' + f'{comment}\n{LINE_BREAK}\n' 

        append_out(variant_info, outname)
        logger.info('Finishing off REPORT genes box with notes .tex')
        cat(f'{resources}/myoseq_template_candidate_variants_notes.tex', outname)


def add_cnv_plot():
    pass


def add_sma_plot():
    pass


def main(args):

    proband = args.proband
    dirname = args.dirname
    resources = args.resources
    unsolved = args.unsolved
    outname = f'{args.out}/{proband}.tex'

    logger.info('Setting up the report header')
    cat(f'{resources}/myoseq_template_header.tex', outname)

    logger.info('Preparing top of first page of reports (ID, sex, ancestry)')
    sex = get_patient_details(proband, dirname, outname)

    logger.info('Finishing first page of reports (REPORT genes box)')
    get_report_variants(proband, dirname, unsolved, outname, resources, sex)

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
