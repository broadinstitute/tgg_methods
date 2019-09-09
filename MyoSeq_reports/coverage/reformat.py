#!/usr/bin/env python

import argparse
import json
import logging
import os
import sys
from utils import get_samples


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('reformat')
logger.setLevel(logging.INFO)


def get_gnomAD_af(jfile, variants, popmax):
    """
    Gets gnomAD popmax AF OR global AC/AN from bigquery JSON output

    :param str jfile: Name of JSON
    :param dict variants: Dictionary; key: variants, value: AC, AN, popmax AF + population
    :param bool popmax: Whether to pull popmax AF
    :return: Dictionary of variants; key: variants, value: AC, AN, popmax AF + pop
    :rtype: dict
    """
     # read in json file as dictionary; format:
    '''
    {
        'chrom': '4',
        'pos': '3494833',
        'ref': 'A',
        'alt': 'AGCCT',
        'global': {
          'adj': 'true',
          'pop': 'all',
          'sex': 'all',
          'ac': '123',
          'an': '193594',
          'af': '0.00063535026911991073',
          'hom': '1'
        },
        'popmax': {
          'pop': 'nfe',
          'ac': '87',
          'af': '0.0010328616202868268',
          'hom': '1'
        }
    }
    '''
    with open(jfile) as j:
        jdict = json.loads(j.read())
        
    # cycle through json dictionary to extract variant information + AC/AN 
    for i in xrange(len(jdict)):
        variant = '{}:{} {}>{}'.format(jdict[i]['chrom'], jdict[i]['pos'], jdict[i]['ref'], jdict[i]['alt'])
        
        # add variant if it isn't already in dict
        if variant not in variants:
            variants[variant] = {}
            variants[variant]['ac'] = 0
            variants[variant]['an'] = 0

        # check which mode the function is being run (popmax == True or popmax == False)
        # add popmax to dict if popmax is True
        if popmax:
            variants[variant]['popmax'] = jdict[i]['popmax']['af']
            variants[variant]['pop'] = jdict[i]['popmax']['pop']

        # add ac/an to variant dict if popmax is False
        else:
            variants[variant]['ac'] += int(jdict[i]['global']['ac'])
            variants[variant]['an'] += int(jdict[i]['global']['an'])

    return variants


def get_fields(seqr, all_variants, samples):
    """
    Reads in seqr export file and stores each sample's variants in dict

    :param str seqr: Name of seqr export TSV
    :param dict all_variants: Dictionary of variants and their frequencies (AC, AN, popmax AF + pop)
    :param list samples: List of requested sample IDs

    :return: Dictionary of samples (key) and their variants with frequencies (value)
    :rtype: dict
    """

    sample_var = {}
    with open(seqr) as s:

        # get header of seqr export file
        header = s.readline().strip().split('\t')

        # cycle through seqr export file and pull relevant fields
        for line in s:
            line = line.strip().split('\t')
            chrom, pos, ref, alt, gene, worst_consequence = line[0:6]
            rsid, hgvsc, hgvsp, clinvar_clinsig, clinvar_gold_stars = line[19:24]

            # set clinvar values to '.' if they don't exist
            clinvar_clinsig = clinvar_clinsig.replace('_', ' ').lower()
            if len(clinvar_clinsig) == 0:
                clinvar_clinsig = '.'
            if len(clinvar_gold_stars) == 0:
                clinvar_gold_stars = '.'

            # get variant
            variant = '{}:{} {}>{}'.format(chrom, pos, ref, alt)

            # get gene name from column - check for gene name bug
            if ',' in gene:
                gene = gene.replace("'","")

            # get popmax and global_af from files; if key doesn't exist, then use '0.0' (not seen in gnomAD)
            try:
                global_af = str(float(all_variants[variant]['ac']) / all_variants[variant]['an'])
                popmax = all_variants[variant]['popmax']
                pop = all_variants[variant]['pop']
            except KeyError:
                global_af = '0.0'
                popmax = '0.0'
                pop = 'N/A'

            # columns needed: gene_name genotype    variant functional_class    hgvs_c  hgvs_p  rsid    exac_global_af  exac_pop_max_af exac_pop_max_population clinvar_clinsig clinvar_clnrevstat  number_of_stars clinvar_url comments
            # get gene name from column - check for gene name bug
            if ',' in gene:
                gene = gene.replace("'","")
    
            # get samples with variant and genotypes for all samples
            sample_info = line[25:]
            for sample in sample_info:
                sample = sample.split(':')

                # skip if sample isn't in list of requested samples
                if sample[0] not in samples:
                    continue
                if sample[0] not in sample_var:
                    sample_var[sample[0]] = []

                # get sample genotype
                try:     
                    geno = int(sample[1])
                except:
                    logger.warn('{} has a non-int genotype'.format(sample))
                    continue
                
                # only keep het or homvar genotypes
                if geno == 0:
                    continue
                elif geno == 1:
                    gt = 'het'
                else:
                    if geno != 2:
                        logger.warn('Genotype field is {}; not 0, 1, or 2'.format(geno))
                        gt = 'het'
                    else:
                        gt = 'hom'

                # columns needed: gene_name genotype    variant functional_class    hgvs_c  hgvs_p  rsid    exac_global_af  exac_pop_max_af exac_pop_max_population clinvar_clinsig clinvar_clnrevstat  number_of_stars clinvar_url comments
                # add formatted line to sample dict
                sample_var[sample].append([gene, gt, variant, worst_consequence, hgvsc, hgvsp, rsid, global_af, popmax, pop, clinvar_clinsig, '.', clinvar_gold_stars, '.'])

    return sample_var
        

def write_out(sample_var, out):
    """
    Writes out files for each sample for MYOSEQ reports

    :param dict sample_var: Dictionary of samples (key) and their variants (with frequencies; value)
    :param str out: Name of output file
    :return: None
    :rtype: None
    """

    logger.info('Number of samples found: {}'.format(sample_var))
    for sample in sample_var:
        outfile = '{}/report_for_MYOSEQ_v20_{}.genes.txt'.format(out, sample)

        # open a file for each sample in specified output directory
        with open(outfile, 'w') as o:

            # write header to output
            o.write('gene_name\tgenotype\tvariant\tfunctional_class\thgvs_c\thgvs_p\trsid\tgnomad_global_af\tgnomad_pop_max_af\tgnomad_pop_max_population\tclinvar_clinsig\tclinvar_clnrevstat\tnumber_of_stars\tclinvar_url\tcomments\n')

            # write variants to output
            for line in sample_var[sample]:
                o.write('\t'.join(line) + '\n')


def main(args):

    variants = {}
    
    logging.info('Getting individual IDs from ped file')
    samples = get_samples(args.slist)

    logging.info('Getting global AC and AN from gnomAD exomes')
    wes_variants = get_gnomAD_af(args.wes, variants, False)

    logging.info('Getting global AC and AN from gnomAD genomes')
    wes_wgs_variants = get_gnomAD_af(args.wgs, wes_variants, False)

    logging.info('Getting gnomAD popmax AF')
    all_variants = get_gnomAD_af(args.popmax, wes_wgs_variants, True)

    logging.info('Getting seqr export information')
    sample_var = get_fields(args.seqr, all_variants, samples)
    
    logging.info('Writing output files')
    write_out(sample_var, args.out)


if __name__ == '__main__':
    
    # Get args from command line
    parser = argparse.ArgumentParser(description='Reformats seqr variant export (after bigquery) to format for Monkol\'s scripts')
    parser.add_argument('-s', '--seqr', help='input file of variants exported from seqr', required=True)
    parser.add_argument('-e', '--wes', help='json file output from bigquery of gnomAD exome frequencies', required=True)
    parser.add_argument('-g', '--wgs', help='json file output from bigquery of gnomAD genome frequencies', required=True)
    parser.add_argument('-p', '--popmax', help='json file output from bigquery of gnomAD popmax', required=True)
    parser.add_argument('-l', '--slist', help='input sample list', required=True)
    parser.add_argument('-o', '--out', help='output directory', required=True)
    args = parser.parse_args()

    main(args)
