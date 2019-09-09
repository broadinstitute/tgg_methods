#!/usr/bin/env python

import argparse
import logging
from utils import get_samples,check_missing_samples 


logging.basicConfig(format='%(asctime)s (%(name)s %(lineno)s): %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger('ancestry')
logger.setLevel(logging.INFO)


def parse_pca(pca, samples, out):
    """
    Get inferred ancestry from seqr sample QC file and write to output TSV

    :param str pca: Name of seqr sample QC file containing ancestry PCA results
    :param list samples: List of sample IDs
    :param str out: Name of output TSV
    :return: None
    :rtype: None
    """
    # create dict to map pca ancestry to ancestry text for report
    # s PROJECT DATA_TYPE   RESEARCH_PROJECT_NUMBER PRODUCT INDIVIDUAL_NAME ANALYSIS_END    PCT_CONTAMINATION   AL_PCT_CHIMERAS WGS_MEAN_COVERAGE   HS_MEAN_TARGET_COVERAGE HS_PCT_TARGET_BASES_20X MYRANK  seqr_id filtered_callrate   filter_flags    qc_platform qc_popprob_afr  prob_amr    prob_asj    prob_eas    prob_est_b1 prob_est_b2 prob_fin    prob_nfe    prob_oth    prob_sas    pop_PC1 pop_PC2 pop_PC3 pop_PC4 pop_PC5 pop_PC6 data_type   known_pop   source  fail_n_snp  fail_r_ti_tv    fail_r_insertion_deletion   fail_n_insertion    fail_n_deletion fail_r_het_hom_var  pop_platform_filters    sample_qc.dp_stats.mean sample_qc.dp_stats.stdev    sample_qc.dp_stats.min  sample_qc.dp_stats.max  sample_qc.gq_stats.mean sample_qc.gq_stats.stdev    sample_qc.gq_stats.minsample_qc.gq_stats.max    sample_qc.call_rate sample_qc.n_called  sample_qc.n_not_called  sample_qc.n_filtered    sample_qc.n_hom_ref sample_qc.n_het sample_qc.n_hom_var sample_qc.n_non_ref sample_qc.n_singleton   sample_qc.n_snp sample_qc.n_insertion   sample_qc.n_deletion    sample_qc.n_transition  sample_qc.n_transversion    sample_qc.n_star    sample_qc.r_ti_tv   sample_qc.r_het_hom_var sample_qc.r_insertion_deletion  sample_qc.f_inbreeding.f_stat   sample_qc.f_inbreeding.n_called sample_qc.f_inbreeding.expected_homs    sample_qc.f_inbreeding.observed_homs
    popmap = {'amr': 'Admixed American', 'asj': 'Ashkenazi Jewish', 'eas': 'East Asian', 'fin': 'Finnish', 'nfe': 'European', 'oth': 'Other', 'sas': 'Southeast Asian', 'afr': 'African'}

    with open(pca) as p, open(out, 'w') as o:

        # pull indices for relevant columns from header
        header = p.readline().strip().split('\t')
        sidx = header.index('s')
        pidx = header.index('qc_pop')
        
        # get inferred ancestry for desired samples and write to output
        found = []
        for line in p:
            line = line.strip().split('\t')
            s = line[sidx]
            pop = line[pidx]
            if s in samples:
                found.append(s)
                o.write('{}\t{}\n'.format(s, popmap[pop]))

    # check if any samples are missing and print missing samples to stdout
    check_missing_samples(samples, found, 'pca')
    

def main(args):

    logger.info('Getting samples from input file')
    samples = get_samples(args.samp)

    logger.info('Getting inferred ancestry from PCA and writing to output')
    parse_pca(args.pca, samples, args.out)


if __name__ == '__main__':
    
    # Get args from command line
    parser = argparse.ArgumentParser(description='Pulls ancestry from PCA file for MYOSEQ reports')
    parser.add_argument('-p', '--pca', help='input TSV (PCA)', required=True)
    parser.add_argument('-s', '--samp', help='input list of samples', required=True)
    parser.add_argument('-o', '--out', help='output TSV', required=True)
    args = parser.parse_args()

    main(args)
