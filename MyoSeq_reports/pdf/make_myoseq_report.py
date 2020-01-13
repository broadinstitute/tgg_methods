#!/usr/bin/env python

import argparse
import logging
import pandas as pd
from decimal import Decimal
from subprocess import Popen, PIPE
import sys


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('make_myoseq_report')
logger.setLevel(logging.INFO)


# define globals for pdflatex
# the line break is just two (\\); using the extra two to escape the tow needed 
LINE_BREAK = '\\\\'
END_BOX = '\\end{tabular}\n\\end{small}\n'
NEW_PAGE = '\\newpage\n'
GT = '\\textgreater '
UNDERSCORE = ' \\textunderscore '


def awk(sample: str, fname: str) -> str:
    '''
    Quickly looks up one line from a file using awk

    :param str sample: Sample ID to find
    :param str fname: Name of file to search
    :return: Line from file
    :rtype: str
    '''
    cmd = '''awk '$1 == "{}"' {}'''.format(sample, fname)
    line, err = Popen([cmd], stdout=PIPE, stderr=PIPE, shell=True).communicate()
    return line.decode('UTF-8')


def cat(fname: str, outname: str) -> None:
    '''
    cats one file into output file.

    :param fname: File name with information to be cat
    :param str outname: Output file name
    :return: None
    '''
    cmd = '''cat {} >> {}'''.format(fname, outname)
    Popen([cmd], shell=True).communicate()


def append_out(info: str, outname: str) -> None:
    '''
    Opens output file and appends information

    :param str info: Information to be added into file
    :param str outname: Output file name
    :return: None
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
    patient_details += '{\\large \\textbf{Inferred Ancestry (Reported):} ' + f'{ancestry} (European)' + '}' + f'\n{LINE_BREAK} {LINE_BREAK} {LINE_BREAK}\n' 
    
    append_out(patient_details, outname)
    return inferred


def report_cnvs(proband: str, dirname: str, outname: str, resources: str) -> tuple:
    '''
    Formats any candidate CNVs for first page of report

    :param str proband: Proband ID
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :return: Gene name of CNV, copy number
    :rtype: str
    '''
    cat(f'{resources}/myoseq_template_report_cnv_table_header.tex', outname)
    
    # add empty line to top of CNV box for formatting
    cnv_info = '{} & {} & {} & {} & {}' + f' {LINE_BREAK}\n'

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
            size = float(start) - float(end) / 1000
            if int(cn) <= 1:
                cnv_type = f'Loss (CN={cn})'
            else:
                cnv_type = f'Gain (CN={cn})'
            cnv_info += f'{gene} & ' + '\\textbf{\\color{red}' + f'{cnv_type}' + '} & ' + f'WES & {cnv} & {size} {LINE_BREAK} \\hline \n'
    cnv_info += f'{END_BOX}'
    append_out(cnv_info, outname)

    logger.info('Closing out CNV section of first page of reports')
    cat(f'{resources}/myoseq_template_candidate_cnv_notes.tex', outname)
    return (gene, int(cn))
            

def report_sma(proband: str, dirname: str, outname: str, resources: str) -> None:
    '''
    Formats SMA information for first page of report

    :param str proband: Proband ID
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :return: None
    '''
    cat(f'{resources}/myoseq_template_report_sma_table_header.tex', outname)

    logger.info('Creating SMA box')    
    sma_info = '{} & {} & {} '+ f'{LINE_BREAK}\n'
    sma_info += 'SMN1 & \\textbf{\\color{red}Loss (CN=0)} & WES ' + f'{LINE_BREAK}' + ' \\hline\n'
    sma_info += f'{END_BOX}'
    append_out(sma_info, outname)

    logger.info('Closing out SMA section of first page of reports')
    cat(f'{resources}/myoseq_template_candidate_cnv_notes.tex')


def get_variants(line: dict, report: bool, comment_index: int) -> str:
    '''
    Formats variant information for one row in variant table

    :param dict line: Line with variant information from patient.flagged.txt (REPORT variants) or patient.genes.txt (all rare variants called)
    :param bool report: Whether this variant is going on the first page of the reports (was tagged REPORT in seqr)
    :param int comment_index: Which index to use for comments
    :return: Formatted variant information for .tex file
    :rtype: str
    '''

    # break variant into component chr:pos, ref, alt and get rest of variant information
    #'variant': 'chr1:236883444 T>A',
    chr_pos = line['variant'].split(' ')[0]
    ref, alt = line['variant'].split(' ')[1].split('>')
    genotype = line['genotype']
    transcript, hgvsc = line['hgvs_c'].replace('>', f'{GT}').replace('_', f'{UNDERSCORE}').split(':')
    try:
        hgvsp = line['hgvs_p'].split(':')[1].replace('>', f'{GT}').replace('_', f'{UNDERSCORE}')
    except:
        hgvsp = ''
    function = line['functional_class']
    global_af = '%.2E' % Decimal(line['gnomad_global_af'])
    popmax_af = '%.2E' % Decimal(line['gnomad_pop_max_af'])
    popmax_pop = line['gnomad_pop_max_population']
    stars = line['number_of_stars']
    rsid = str(line['rsid'])
    if rsid == 'nan':
        rsid = ''

    # format notes from seqr
    # NOTE: assumes you won't have more comments than letters in the alphabet
    comments = line['comments']
    if comments != '.' and report:
        comment_foot = chr(97+comment_index)
        comment_info = ''
        comments = comments.split('|')
        # 	kchao@broadinstitute.org: https://www.ncbi.nlm.nih.gov/pubmed/28554942|kchao@broadinstitute.org: https://www.ncbi.nlm.nih.gov/pubmed/28544275|kchao@broadinstitute.org: not seen in gnomAD
        for c in comments:
            c = c.split(': ')[1].replace('http', '\\\\\nhttp')
            comment_info += c + ' ' 
    else:
        comment_foot = ''
        comment_info = ''

    # truncate long indels
    indel_fmt = '\\textsuperscript{+}'
    if len(ref) > 10:
        ref = ref[:10] + f'{indel_fmt}'
    if len(alt) > 10:
        alt = alt[:10] + f'{indel_fmt}'
    if len(hgvsc) > 30:
        hgvsc = hgvsc[:30] + f'{indel_fmt}'

    # format variant into chr:pos ref>alt
    variant = f'{chr_pos} {ref} {GT} {alt}'
    variant_info = '\\bf{' + f"{line['gene_name']}" + '} & ' + f'{variant} & {global_af} &'

    # correct hom to hem for males on X
    if 'X' in chr_pos and genotype == 'hom' and sex == 'Male':
        genotype = 'hem'

    # parse ClinVar status and determine color
    # Monkol split on "/" to split up entries that are too long (pathogenic/likely pathogenic)
    clinsig = line['clinvar_clinsig'].split('/')
    revstat = line['clinvar_clnrevstat']
    if clinsig == ['.']:
        variant_info += f' NA & {comment_foot} {LINE_BREAK}\n'
    else:
        color = ''
        for i in clinsig:
            if 'pathogenic' in i or 'conflicting' in i:
                color = 'red'
            elif 'uncertain' in i:
                color = 'Orange'
            elif 'benign' in i:
                color = 'Green'
        
        if color != '':
            for i in range(len(clinsig)):
                variant_info += '\\textcolor{' + f'{color}' + '}{\\textbf{' + f'{clinsig[i]}' + '}} & ' + f'{comment_foot} {LINE_BREAK}\n'
        else:
            for i in range(len(clinsig)):
                variant_info += '{\\textbf{' + f'{clinsig[i]}' + '}} & ' + f'{comment_foot} {LINE_BREAK}\n'

    # color variant's hgvsc red if it is a splice variant; otherwise display in black
    if function == 'splice_donor_variant' or function == 'slice_acceptor_variant':
        variant_info += f'{transcript} & ' + '\\textcolor{red}{' + f'{hgvsc}' + '} & ' + f'{popmax_af} & ' + f'{revstat} & {LINE_BREAK}\n'
    else:
        variant_info += f'{transcript} & {hgvsc} & {popmax_af} & {revstat} & {LINE_BREAK}\n'

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
            hgvsp =  '\\textcolor{' + f'{color}' + '}{' + f'{hgvsp}' + '}'
    variant_info += '\\textit{' + f'{genotype}' + '} & ' + f'{hgvsp} & & {stars} & {LINE_BREAK}\n'
    variant_info += f'& {rsid} & & & {LINE_BREAK}\n\\hline\n'

    if comments != '.':
        variant_info += '\\textsuperscript{' + f'{comment_foot}' + '} ' + f'{comments}\n{LINE_BREAK}\n' 

    return variant_info


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
    '''
    if unsolved:
        logger.info('Adding None to REPORT genes box')
        cat(f'{resources}/myoseq_template_no_candidate_variants_notes.tex', outname)
    else:
        logger.info('Starting REPORT genes box')
        cat(f'{resources}/myoseq_template_report_variants_table_header.tex', outname)
        report_file = f'{dirname}/seqr/variants/{proband}.flagged.txt'
        report_info = pd.read_csv(report_file, sep='\t').to_dict('records')
        variant_info = ''
        for i in range(len(report_info)):
            variant_info += get_variants(report_info[i], True, i)
        variant_info += f'{END_BOX}'      
        append_out(variant_info, outname)

        logger.info('Finishing off REPORT genes box with notes .tex')
        cat(f'{resources}/myoseq_template_candidate_variants_notes.tex', outname)


def add_cnv_plot(proband: str, gene: str, cn: int, dirname: str, resources: str, outname: str) -> None:
    '''
    Adds CNV plot to second page of report

    :param str proband: Proband ID
    :param str gene: Gene of interest
    :param int cn: Copy number of CNV
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :return: None
    '''
    cnv_figure = f'{dirname}/images/cnv/{proband}'
    cnv_info = '\\includegraphics[width=14cm, height=14cm]{' + f'{cnv_figure}' + '} {LINE_BREAK}\n'
    cnv_type = 'Deletion'
    if cn > 1:
        cnv_type = 'Duplication' 
    cnv_info += '\\textbf{Figure 1. ' + f'{cnv_type} of exons in {gene} detected by gCNV.' + '}' + f'{LINE_BREAK}\n'
    append_out(cnv_info, outname)
    cat(f'{resources}/myoseq_template_cnv_figure_caption.tex', outname)
    append_out(NEW_PAGE, outname)


def add_sma_plot(proband: str, dirname: str, outname: str) -> None:
    '''
    Adds SMA plot to second page of report

    :param str proband: Proband ID
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :return: None
    '''
    sma_info = '{\\Large \\textbf{{SMA Carrier Detection}}} ' + f'{LINE_BREAK} {LINE_BREAK} \n'
    sma_info += 'Detected SMN1 Copy Number deletion (CN=0) in ' + f'{proband} {LINE_BREAK} {LINE_BREAK} {LINE_BREAK}\n'
    sma_figure = f'{dirname}/images/sma/{proband}_carrier_probabilities_plot.pdf'
    sma_info += '\\includegraphics[width=14cm, height=8cm]{' + f'{sma_figure} {LINE_BREAK}\n'
    sma_info += '\\textbf{Figure 1. SMA copy number detection showing all samples in MYOSEQ project}' + f'{LINE_BREAK}\n'
    append_out(sma_info, outname)
    cat(f'{resources}/myoseq_template_sma_figure_caption.tex', outname)
    append_out(NEW_PAGE, outname)   
 

def get_all_variants(proband: str, dirname: str, outname: str, resources: str) -> None:
    '''
    Formats any candidate SNVs/indels (tagged 'REPORT' in seqr) for first page of report

    :param str proband: Proband ID
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :return: None
    '''
    logger.info('Starting appendix (all rare variants) box')
    header_text = '{\\Large \\textbf{\\underline{Appendix - Candidate Gene List}}}\n' + f'{LINE_BREAK} {LINE_BREAK}\n'
    append_out(header_text, outname)

    cat(f'{resources}/myoseq_template_all_variants_table_header.tex', outname)
    myoseq_file = f'{dirname}/seqr/variants/{proband}.genes.txt'
    myoseq_info = pd.read_csv(myoseq_file, sep='\t').to_dict('records')
    variant_info = ''
    for i in range(len(myoseq_info)):
        variant_info += get_variants(myoseq_info[i], False, i)
    variant_info += '\\end{longtable}\n\\end{small}\n'      
    append_out(variant_info, outname)

    logger.info('Finishing off appendix with notes .tex')
    cat(f'{resources}/myoseq_template_additional_variants_notes.tex', outname)


def get_coverage_table(proband: str, dirname: str, outname: str, resources: str) -> None:
    '''
    Formats sample coverage across MyoSeq gene list into table for appendix

    :param str proband: Proband ID
    :param str dirname: Name of top level directory with all sample information for reports
    :param str outname: Output file name
    :param str resources: Directory with tex files
    :return: None
    '''
    logger.info('Starting gene coverage table')
    cov_header = '{\\large \\textbf{\\underline{Appendix}}} ' + f'{LINE_BREAK} {LINE_BREAK} {LINE_BREAK}\n'
    cov_header += '{\\large \\textbf{\\textit{Gene Coverage}}}\n'
    append_out(cov_header, outname)
    cat(f'{resources}/myoseq_template_coverage_table_header.tex', outname)

    cov_file = f'{dirname}/per_sample/{proband}_coverage.tsv'
    cov_info = ''
    with open(cov_file) as c:
        # Gene	Cohort_Mean	Sample_Mean	Cohort_Callable	Sample_Callable	Cohort_Uncallable	Sample_Uncallable
        header = c.readline()
        for line in c:
            line = line.strip().split('\t')
            row = ' & '.join(line)

            # color line if sample has uncallable sites
            if int(line[-1]) > 0:
                cov_info += '\\rowcolor{Lavender} ' 
            cov_info += f'{row} {LINE_BREAK}\n'
    cov_info += '\\end{longtable}\n'
    append_out(cov_info, outname)
    cat(f'{resources}/myoseq_template_coverage_notes.tex', outname)


def main(args):

    proband = args.proband
    dirname = args.dirname
    resources = args.resources
    unsolved = args.unsolved
    cnv = args.cnv
    sma = args.sma
    outname = f'{args.out}/{proband}.tex'

    logger.info('Setting up the report header')
    cat(f'{resources}/myoseq_template_header.tex', outname)

    logger.info('Preparing top of first page of reports (ID, sex, ancestry)')
    sex = get_patient_details(proband, dirname, outname)

    if cnv:
        logger.info('Adding candidate CNVs to first page of reports')
        cnv = report_cnvs(proband, dirname, outname, resources) 
        gene = cnv[0]
        cn = cnv[1]
    if sma:
        logger.info('Adding SMA carrier status to first page of reports')
        report_sma(proband, dirname, outname, resources)

    logger.info('Finishing first page of reports (REPORT genes box)')
    get_report_variants(proband, dirname, unsolved, outname, resources, sex)

    if cnv:
        logger.info('Adding CNV plot to second page of reports')
        add_cnv_plot(proband, gene, cn, dirname, resources, outname) 
    if sma:
        logger.info('Adding SMA plot to second page of reports')
        add_sma_plot(proband, dirname, outname) 

    logger.info('Adding MyoSeq gene list to reports (starting appendix)')
    cat(f'{resources}/myoseq_template_appendix_gene_list.tex', outname)

    logger.info('Creating table for all rare variants called in MyoSeq gene list')
    get_all_variants(proband, dirname, outname, resources)

    logger.info('Adding gene coverage table to appendix')
    get_coverage_table(proband, dirname, outname, resources)
 
    logger.info('Adding methods section to reports')
    cat(f'{resources}/myoseq_template_general_methodology_cnv_sma_nogenelist.tex', outname)

    logger.info('Finishing report')
    append_out('\\end{document}\n', outname)

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Creates patient section of tex file for pdflatex')
    parser.add_argument('-c', '--cnv', help='Samples have CNV candidates', action='store_true')
    parser.add_argument('-s', '--sma', help='Samples are SMA carriers', action='store_true')
    parser.add_argument('-u', '--unsolved', help='Sample is unsolved (no candidates)', action='store_true')
    parser.add_argument('-d', '--dirname', help='Top level directory that contains all data for reports', required=True)
    parser.add_argument('-r', '--resources', help='Directory containing tex files', required=True)
    parser.add_argument('-p', '--proband', help='Proband ID', required=True)
    parser.add_argument('-o', '--out', help='Output directory', required=True)
    args = parser.parse_args()

    main(args)
