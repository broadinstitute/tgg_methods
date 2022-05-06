import logging
import hail as hl
import hail.expr.aggregators as agg

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from os.path import dirname

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

from gnomad.utils.reference_genome import get_reference_genome

def get_y_cov(mt: hl.MatrixTable, build: str, call_rate_threshold: float=0.25) -> hl.Table:
    """
    Calculates mean chromosome Y coverage, mean chromosome 20 coverage, and normalized chromosome Y coverage for each sample
    :param MatrixTable mt: MatrixTable containing samples with chrY variants
    :param str build: Reference used, either GRCh37 or GRCh38
    :param float call_rate_threshold: Minimum required call rate
    :return: Table with coverage annotations
    :rtype: Table
    """
    #could utilize hail's reference genome chr names
    if build == "GRCh38":
        chr20 = "chr20"
        chrY = "chrY"
    else:
        chr20 = "20"
        chrY = "Y"

    logging.info("Filtering to chromosome 20 and non-par regions on chromosome Y to calculate normalized Y coverage...")
    sex_mt = hl.filter_intervals(mt, [hl.parse_locus_interval(chr20, reference_genome=build), hl.parse_locus_interval(chrY, reference_genome=build)])
    sex_mt = sex_mt.filter_rows((sex_mt.locus.contig == chr20) | (sex_mt.locus.in_y_nonpar()), keep = True)


    # filter to common snvs above defined callrate (because only biallelic, should now only have one index in the array)
    sex_mt = sex_mt.filter_rows(sex_mt.AF > 0.01)

    #optional for gnomad methods
    sex_mt = hl.variant_qc(sex_mt) 
    sex_mt = sex_mt.filter_rows(sex_mt.variant_qc.call_rate > call_rate_threshold)

    logging.info("Calculating mean coverage on chromosome 20 and Y...")
    sex_mt = sex_mt.annotate_cols(chrY_mean_cov=hl.agg.filter(sex_mt.locus.contig == chrY, hl.agg.mean(sex_mt.DP)),
                                  chr20_mean_cov=hl.agg.filter(sex_mt.locus.contig == chr20, hl.agg.mean(sex_mt.DP)))
    
    #hail.cond is outdated
    sex_mt = sex_mt.annotate_cols(normalized_y_coverage=hl.cond(sex_mt.chr20_mean_cov > 0, sex_mt.chrY_mean_cov/sex_mt.chr20_mean_cov, -99))
    sex_ht = sex_mt.cols()

    return(sex_ht)

def get_x_cov(mt: hl.MatrixTable, build: str, call_rate_threshold: float=0.25) -> hl.Table:
    """
    Calculates mean chromosome X coverage,and normalized chromosome X coverage for each sample
    :param MatrixTable mt: MatrixTable containing samples with chrY variants
    :param str build: Reference used, either GRCh37 or GRCh38
    :param float call_rate_threshold: Minimum required call rate
    :return: Table with coverage annotations
    :rtype: Table
    """
    if build == "GRCh38":
        chr20 = "chr20"
        chrX = "chrX"
    else:
        chr20 = "20"
        chrX = "X"

    logging.info("Filtering to chromosome 20 and non-par regions on chromosome X to calculate normalized Y coverage...")
    sex_mt = hl.filter_intervals(mt, [hl.parse_locus_interval(chr20, reference_genome=build), hl.parse_locus_interval(chrX, reference_genome=build)])
    sex_mt = sex_mt.filter_rows((sex_mt.locus.contig == chr20) | (sex_mt.locus.in_x_nonpar()), keep = True)


    # filter to common snvs above defined callrate (because only biallelic, should now only have one index in the array)
    sex_mt = sex_mt.filter_rows(sex_mt.AF > 0.01)

    #optional for gnomad methods
    sex_mt = hl.variant_qc(sex_mt) 
    sex_mt = sex_mt.filter_rows(sex_mt.variant_qc.call_rate > call_rate_threshold)

    logging.info("Calculating mean coverage on chromosome 20 and X...")
    sex_mt = sex_mt.annotate_cols(chrX_mean_cov=hl.agg.filter(sex_mt.locus.contig == chrX, hl.agg.mean(sex_mt.DP)),
                                  chr20_mean_cov=hl.agg.filter(sex_mt.locus.contig == chr20, hl.agg.mean(sex_mt.DP)))
    
    #hail.cond is outdated
    sex_mt = sex_mt.annotate_cols(normalized_x_coverage=hl.cond(sex_mt.chr20_mean_cov > 0, sex_mt.chrX_mean_cov/sex_mt.chr20_mean_cov, -99))
    sex_ht = sex_mt.cols()

    return(sex_ht)
    

def run_hails_impute_sex(mt: hl.MatrixTable, 
    build: str ,
    outdir: str, 
    callset_name: str,
    male_fstat_threshold: float=0.75, 
    female_fstat_threshold: float=0.5, 
    aaf_threshold: float=0.05) -> hl.Table:
    """
    Imputes sex and annotates MatrixTable with results, outputs a histogram of fstat values
    :param MatrixTable mt: MatrixTable containing samples to be ascertained for sex
    :param str build: Reference used, either GRCh37 or GRCh38
    :param str callset_name: Basename for callset and output results
    :param str outdir: Directory to output results
    :param str male_fstat_threshold: Fstat threshold above which a sample will be called male
    :param str female_fstat_threshold: Fstat threshold below which a sample will be called female
    :param str aaf_threshold: Minimum alternate allele frequency required 
    :return: Table with imputed sex annotations
    :rtype: Table
    """
    
    logging.warn("User needs to decide on fstat thresholds for male/female by looking at fstat plots")

    
    # filter to the X chromosome and impute sex
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(get_reference_genome(mt.locus).x_contigs[0], reference_genome=build)])   
    sex_ht = hl.impute_sex(mt.GT, aaf_threshold=aaf_threshold, male_threshold=male_fstat_threshold, female_threshold=female_fstat_threshold)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    sex_ht = mt.cols()
    
    # plot histogram of fstat values
    df = sex_ht.to_pandas()
    plt.clf()
    plt.hist(df["f_stat"])
    plt.xlabel("Fstat")
    plt.ylabel("Frequency")
    plt.axvline(male_fstat_threshold, color='blue', linestyle='dashed', linewidth=1)
    plt.axvline(female_fstat_threshold, color='red', linestyle='dashed', linewidth=1)

    outplot = outdir + "/fstat_{callset_name}.png".format(**locals())
    with hl.hadoop_open(outplot, 'wb') as out:
        plt.savefig(out)   
        
    return sex_ht  

    
def call_sex(callset: str, 
             use_y_cov: bool=False,
             y_cov_threshold: float=0.1,
             male_fstat_threshold: float=0.75, 
             female_fstat_threshold: float=0.5, 
             aaf_threshold: float=0.05,
             call_rate_threshold: float=0.25) -> hl.Table:
    
    """
    Calls sex for the samples in a given callset and exports results file to the callset directory
    :param str callset: String of full matrix table path for the callset
    :param bool use_y_cov: Set to True to calculate and use chrY coverage for sex inference
    :param float y_cov_threshold: Coverage on chrY above which supports male call
    :param str male_fstat_threshold: Fstat threshold above which a sample will be called male
    :param str female_fstat_threshold: Fstat threshold below which a sample will be called female
    :param str aaf_threshold: Minimum alternate allele frequency required 
    :param float call_rate_threshold: Minimum required call rate
    :return: Table with sex annotations
    :rtype: Table
    """

    # read in matrix table and define output directory
    # will need to generalize more for gnomad methods
    outdir = dirname(callset)
    mt_name = callset.split("/")[-1].strip("\.mt")
    logging.info("Reading matrix table for callset: {callset}".format(**locals()))
    logging.info("Using chromosome Y coverage? {use_y_cov}".format(**locals()))

    mt = hl.read_matrix_table(callset)
    
    # filter to SNVs and biallelics
    # filter to biallelics might be a function
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]))
    
    # filter to pass variants only (empty set)
    # gnomad methods - will need to make this an optional argument
    mt = mt.filter_rows(mt.filters.length() == 0, keep=True)
    
    #infer build:
    build = get_reference_genome(mt.locus).name
    logging.info("Build inferred as {build}".format(**locals()))
    
    logging.info("Inferring sex...")
    #for production change female and male to XX and XY
    if use_y_cov:
        cov_ht = get_y_cov(mt, build, call_rate_threshold)
        mt = mt.annotate_cols(**cov_ht[mt.col_key])
        sex_ht = run_hails_impute_sex(mt, build, outdir, mt_name, male_fstat_threshold, female_fstat_threshold, aaf_threshold)  
        sex_ht = sex_ht.annotate(ambiguous_sex=hl.is_missing(sex_ht.is_female),
                                 sex_aneuploidy=(sex_ht.is_female) & hl.is_defined(sex_ht.normalized_y_coverage) & (sex_ht.normalized_y_coverage > y_cov_threshold) |
                                 (~sex_ht.is_female) & hl.is_defined(sex_ht.normalized_y_coverage) & (sex_ht.normalized_y_coverage < y_cov_threshold))

        sex_expr = (hl.case()
            .when(sex_ht.ambiguous_sex, "ambiguous_sex")
            .when(sex_ht.sex_aneuploidy, "sex_aneuploidy")
            .when(sex_ht.is_female, "female")
            .default("male"))

        sex_ht = sex_ht.annotate(sex=sex_expr)
        sex_ht = sex_ht.select(sex_ht.is_female, sex_ht.f_stat, sex_ht.n_called, sex_ht.expected_homs, sex_ht.observed_homs, sex_ht.sex, sex_ht.chrY_mean_cov, sex_ht.chr20_mean_cov, sex_ht.normalized_y_coverage)
        
    else:
        sex_ht = run_hails_impute_sex(mt, build, outdir, mt_name, male_fstat_threshold, female_fstat_threshold, aaf_threshold)
        sex_ht = sex_ht.annotate(ambiguous_sex=hl.is_missing(sex_ht.is_female))
        sex_expr = hl.cond(sex_ht.ambiguous_sex, "ambiguous_sex", hl.cond(sex_ht.is_female, "female", "male"))
        sex_ht = sex_ht.annotate(sex=sex_expr)
        sex_ht = sex_ht.select(sex_ht.is_female, sex_ht.f_stat, sex_ht.n_called, sex_ht.expected_homs, sex_ht.observed_homs, sex_ht.sex)

    outfile = outdir + f"/sex_{mt_name}.txt"
    sex_ht.export(outfile)
    return(sex_ht)

def main(args):
    call_sex(**vars(args))

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='This script infers the sex of samples')
    parser.add_argument('-i', '--callset', required=True, help='Path to Callset MatrixTable')
    parser.add_argument('-u', '--use-y-cov', help='bool for whether or not to use chrY coverage', action='store_true')
    parser.add_argument('-y', '--y-cov-threshold', help='chrY coverage threshold', default=0.1)
    parser.add_argument('-m', '--male-fstat-threshold', help='male_fstat_threshold for hails impute_sex', default=0.75)
    parser.add_argument('-f', '--female-fstat-threshold', help='female_fstat_threshold for hails impute_sex', default=0.50)
    parser.add_argument('-a', '--aaf-threshold', help='aaf_threshold for hails impute_sex', default=0.05)
    parser.add_argument('-c', '--call-rate-threshold', help='call_rate_threshold required to use chrY variant', default=0.25)
    
    args = parser.parse_args()
    main(args)
    