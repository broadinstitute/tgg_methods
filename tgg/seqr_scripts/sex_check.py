import argparse
import logging
import matplotlib.pyplot as plt
import pandas as pd
from os.path import dirname

import hail as hl
import hail.expr.aggregators as agg

from gnomad.utils.reference_genome import get_reference_genome

from gnomad.utils.reference_genome import get_reference_genome


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("sex_check")
logger.setLevel(logging.INFO) 


def get_chr_cov(mt: hl.MatrixTable, build: str, call_rate_threshold: float=0.25, chr_name: str) -> hl.Table:
    """
    Calculate mean chromosome coverage.
    
    :param mt: MatrixTable containing samples with chrY variants
    :param build: Reference used, either GRCh37 or GRCh38
    :param call_rate_threshold: Minimum required call rate
    :param chr_name: chosen chromosome, either 1, 2, 3, ..., X, or Y
    :return: Table with coverage annotations
    """

    logger.warning("This function expects the chrY to be at index 23 and chrX to be at index 22")
    
    if chr_name == "Y"
        filter_nonpar_expr = sex_mt.locus.in_y_nonpar()
        chr_place = 23
    elif chr_name = "X"
        filter_nonpar_expr = sex_mt.locus.in_x_nonpar()
        chr_place = 22
    else:
        try:
            # Chromsome index in `.contigs` list should be one less than the chromsome number
            chr_place = int(chr_name) - 1
        except:
            logger.warning("chr_name cannot be converted to an integer")

    chr_name = hl.get_reference(build).contigs[chr_place]

    logging.info(f"Filtering to non-par regions on {chr_name} to calculate normalized {chr_name} coverage...")
    sex_mt = hl.filter_intervals(mt, [hl.parse_locus_interval(chr_name, reference_genome=build)])
    if chr_place == 22 or chr_place == 23 :
        sex_mt = sex_mt.filter_rows((filter_nonpar_expr), keep = True)


    #filter to common snvs above defined callrate (because only biallelic, should now only have one index in the array)
    sex_mt = sex_mt.filter_rows(sex_mt.AF > 0.01)

    #optional for gnomad methods
    sex_mt = hl.variant_qc(sex_mt) 
    sex_mt = sex_mt.filter_rows(sex_mt.variant_qc.call_rate > call_rate_threshold)

    logging.info(f"Calculating mean coverage on {chr_name}...")
    sex_mt = sex_mt.annotate_cols(**{f"{chr_name}_mean_dp": mt.aggregate(hl.agg.mean(mt.DP))})
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
    :param mt: MatrixTable containing samples to be ascertained for sex
    :param build: Reference used, either GRCh37 or GRCh38
    :param callset_name: Basename for callset and output results
    :param outdir: Directory to output results
    :param male_fstat_threshold: Fstat threshold above which a sample will be called male
    :param female_fstat_threshold: Fstat threshold below which a sample will be called female
    :param aaf_threshold: Minimum alternate allele frequency required 
    :return: Table with imputed sex annotations
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

    outplot = f"{outdir}/fstat_{callset_name}.png"
    with hl.hadoop_open(outplot, 'wb') as out:
        plt.savefig(out)   
        
    return sex_ht  

    
def call_sex(callset: str, 
             use_y_cov: bool=False,
             y_cov_threshold: float=0.1,
             normalization_contig: str = "chr20",
             male_fstat_threshold: float=0.75, 
             female_fstat_threshold: float=0.5, 
             aaf_threshold: float=0.05,
             call_rate_threshold: float=0.25) -> hl.Table:
    
    """
    Calls sex for the samples in a given callset and exports results file to the callset directory
    :param callset: String of full matrix table path for the callset
    :param use_y_cov: Set to True to calculate and use chrY coverage for sex inference
    :param y_cov_threshold: Coverage on chrY above which supports male call
    :param male_fstat_threshold: Fstat threshold above which a sample will be called male
    :param female_fstat_threshold: Fstat threshold below which a sample will be called female
    :param aaf_threshold: Minimum alternate allele frequency required 
    :param call_rate_threshold: Minimum required call rate
    :return: Table with sex annotations
    """

    # read in matrix table and define output directory
    # will need to generalize more for gnomad methods
    outdir = dirname(callset)
    mt_name = callset.split("/")[-1].strip("\.mt")
    logging.info(f"Reading matrix table for callset: {callset}")
    logging.info(f"Using chromosome Y coverage? {use_y_cov}")

    mt = hl.read_matrix_table(callset)
    
    # filter to SNVs and biallelics
    # filter to biallelics might be a function
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]))
    
    # filter to pass variants only (empty set)
    # gnomad methods - will need to make this an optional argument
    mt = mt.filter_rows(mt.filters.length() == 0, keep=True)
    
    # Infer build:
    build = get_reference_genome(mt.locus).name
    logging.info(f"Build inferred as {build}")
    
    logging.info("Inferring sex...")
    # TODO: Change "female" and "male" to "XX" and "XY"
    if use_y_cov:
        mt = mt.annotate_cols(
            **{
                f"{normalization_contig}_mean_dp": get_chr_cov(mt, build, call_rate_threshold, normalization_contig),
                "chrY_mean_dp": get_chr_cov(mt, build, call_rate_threshold, chrY_name),
            }
        )
        mt = mt.annotate_cols(chry_normalized_cov=hl.or_missing(mt[f"{normalization_contig}_mean_dp"] > 0, chry_mean_dp, mt[f"{normalization_contig}_mean_dp"]/chrY_mean_dp,))
        mt = mt.annotate_cols(**mt[mt.col_key])
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
        sex_expr = hl.if_else(sex_ht.ambiguous_sex, "ambiguous_sex", hl.if_else(sex_ht.is_female, "female", "male"))
        sex_ht = sex_ht.annotate(sex=sex_expr)
        sex_ht = sex_ht.select(sex_ht.is_female, sex_ht.f_stat, sex_ht.n_called, sex_ht.expected_homs, sex_ht.observed_homs, sex_ht.sex)

    outfile = f"{outdir}/sex_{mt_name}.txt"
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
    parser.add_argument('-n', '--chr-name', requried=True, help='Input chromosome name for get_chr_cov')

    args = parser.parse_args()
    main(args)
