import hail as hl
import hail.expr.aggregators as agg
hl.init()

import numpy as np
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter, defaultdict
from os.path import basename, splitext, isfile, dirname
from math import log, isnan
from pprint import pprint


def call_sex(callset: str) -> hl.MatrixTable:
    """
    Calls sex for the samples in a given callsetand exports sex.txt file to the callset directory
    :param str callset: string of full matrix table path for the callset
    :return: None
    """
    
    def get_y_cov(mt: hl.MatrixTable) -> hl.Table:
        """
        Calculates mean Y coverage, mean 20 coverage, and normalized Y coverage for each sample
        :param MatrixTable mt: matrix table containing exome samples
        :return: Hail table with coverage annotations
        :rtype: Table
        """
        # filter to chromosome 20 and non-par Y regions
        snv_mt = mt.filter_rows((mt.locus.contig == "20") | (mt.locus.in_y_nonpar()),keep = True)

        # pull out the number of samples 
        num_samples = snv_mt.count_cols()

        # filter to pass variants only (empty set)
        snv_mt = snv_mt.filter_rows(snv_mt.filters.length() == 0, keep=True)
    
        # filter to SNVs and biallelics
        snv_mt = snv_mt.filter_rows((hl.len(snv_mt.alleles) == 2) & hl.is_snp(snv_mt.alleles[0], snv_mt.alleles[1]))
    
        # filter to common snvs(since only biallelic, should now only have one index in the array)
        snv_mt = snv_mt.filter_rows(snv_mt.info.AF[0] > 0.05)
    
        # filter to sites where there was a call in at least x% of the samples
        snv_mt = snv_mt.filter_rows(snv_mt.info.AN > (.5*num_samples))

        # calculate coverage on Y and 20
        snv_mt = snv_mt.annotate_cols(
        y_cov = hl.agg.filter(snv_mt.locus.contig == "Y",hl.agg.mean(snv_mt.DP)))
    
        snv_mt = snv_mt.annotate_cols(
        twenty_cov = hl.agg.filter(snv_mt.locus.contig == "20",hl.agg.mean(snv_mt.DP)))
    
        snv_mt = snv_mt.annotate_cols(normalized_y_coverage = snv_mt.y_cov/snv_mt.twenty_cov)
        snv_ht = snv_mt.cols()

        return(snv_ht)
    
    def annotate_sex(mt: hl.MatrixTable, build: str ,outdir: str, male_threshold: float = 0.75, female_threshold: float = 0.5, aaf_threshold = 0.05) -> hl.Table:
        """
        Imputes sex and annotates matrix table with results, outputs a histogram of fstat values
        NOTE: decide on thresholds for male/female by looking at fstat plots
        :param MatrixTable mt: matrix table containing samples to be ascertained for sex
        :param str build: reference used, either GRCh37 or GRCh38
        :param str outdir: directory to output results
        :param str male_threshold: threshold above which a sample will be called male
        :param str female_threshold: threshold below which a sample will be called female
        :param str aaf_threshold: minimum alternate allele frequency 

        :return: HailTable with imputed sex annotations stashed in column annotation 'sex'
        :rtype: HailTable
        """
    
        # filter to SNVs and biallelics
        mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]))
    
        # filter to the X chromosome and impute sex
        if build == "GRCh37":
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval('X', reference_genome="GRCh37")])
            sex_ht = hl.impute_sex(mt.GT, aaf_threshold = aaf_threshold, male_threshold = male_threshold, female_threshold = female_threshold)

            mt = mt.annotate_cols(**sex_ht[mt.col_key])
            sex_ht = mt.cols()

            sex_ht = sex_ht.annotate(ambiguous_sex = hl.is_missing(sex_ht.is_female),
            sex_aneuploidy = (sex_ht.is_female) & hl.is_defined(sex_ht.normalized_y_coverage) & (sex_ht.normalized_y_coverage > 0.1) |
                                (~sex_ht.is_female) & hl.is_defined(sex_ht.normalized_y_coverage) & (sex_ht.normalized_y_coverage < 0.1))

            sex_expr = (hl.case()
                                .when(sex_ht.ambiguous_sex, "ambiguous_sex")
                                .when(sex_ht.sex_aneuploidy, "sex_aneuploidy")
                                .when(sex_ht.is_female, "female")
                                .default("male"))

            sex_ht = sex_ht.annotate(sex = sex_expr)        
        
        else:
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval('chrX', reference_genome="GRCh38")])
            sex_ht = hl.impute_sex(mt.GT, aaf_threshold = aaf_threshold, male_threshold = male_threshold, female_threshold = female_threshold)
        
            mt = mt.annotate_cols(**sex_ht[mt.col_key])
            sex_ht = mt.cols()
        
            sex_ht = sex_ht.annotate(ambiguous_sex = hl.is_missing(sex_ht.is_female))
            sex_expr = hl.cond(sex_ht.ambiguous_sex, "ambiguous_sex",  hl.cond(sex_ht.is_female, "female", "male"))
            sex_ht = sex_ht.annotate(sex = sex_expr)
        

        # plot histogram of fstat values
        df = sex_ht.to_pandas()
        plt.clf()
        plt.hist(df["f_stat"])
        plt.xlabel("Fstat")
        plt.ylabel("Frequency")
        plt.axvline(male_threshold, color='blue', linestyle='dashed', linewidth=1)
        plt.axvline(female_threshold, color='red', linestyle='dashed', linewidth=1)

        outplot = outdir + "/fstat.png"
        with hl.hadoop_open(outplot, 'wb') as out:
            plt.savefig(out)

        return sex_ht   
    
    
    # read in matrix table and define output directory
    outdir = dirname(callset)
    mt = hl.read_matrix_table(callset)
    print(callset) 

    # grab y coverage and evaluate sex for exomes
    if re.search("_WES_", callset):
        # get y cov
        cov_ht = get_y_cov(mt) # change to mt, return mt
        # annotated the matrix table with the coverage results and impute sex
        mt = mt.annotate_cols(**cov_ht[mt.col_key])        
        sex_ht = annotate_sex(mt, "GRCh37", outdir)
        sex_ht = sex_ht.select(sex_ht.is_female,sex_ht.f_stat,sex_ht.n_called,sex_ht.expected_homs,sex_ht.observed_homs,sex_ht.sex, sex_ht.y_cov,sex_ht.twenty_cov,sex_ht.normalized_y_coverage)

        outfile = outdir + "/sex.txt"
        sex_ht.export(outfile)
        
    # evaluate sex for genomes   
    elif re.search("_WGS_", callset):
        # impute sex
        sex_ht = annotate_sex(mt, "GRCh38", outdir)
        sex_ht = sex_ht.select(sex_ht.is_female,sex_ht.f_stat,sex_ht.n_called,sex_ht.expected_homs,sex_ht.observed_homs,sex_ht.sex)
        
        # output sex results without taking into account y coverage (because genomes do not have y calls)
        outfile = outdir + "/sex.txt"
        sex_ht.export(outfile)



