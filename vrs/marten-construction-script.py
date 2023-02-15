"""
Script that constructs a ~100k Variant Test VCF for GA4GH-VRS Project.

Contact: Daniel Marten (marten@broadinstitute.org)
"""

import argparse
import logging
import random
from random import sample

import hail as hl
from gnomad.resources.grch38.gnomad import public_release
from gnomad_qc.v3.resources.basics import gnomad_v3_genotypes_vds

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_vrs_test_dataset")
logger.setLevel(logging.INFO)


def main(args):

    # Set randomness for Python package "random"
    random.seed(args.py_rand_seed)
    # Initialize Hail with a seed passed as an argument
    hl.init(
        default_reference="GRCh38",
        global_seed=args.hail_rand_seed,
        tmp_dir=args.tmp_dir,
        quiet=True,
        spark_conf={
            "spark.hadoop.fs.gs.requester.pays.mode": "AUTO",
            "spark.hadoop.fs.gs.requester.pays.project.id": f"{args.google_cloud_project}",
        },
    )

    # Import the gnomAD v3.1.2 full variant Table (approx 759,000,000 variants)
    # NOTE: No need to repartition the whole Table, as it was recently repartitioned to 10,000
    ht_whole = public_release("genomes").ht()

    read_if_exists = not args.overwrite

    """
    Previously Calculated Frequency Statistics

    Frequency of the following variant types in gnomAD v3:
    - indels
    - long reference alleles (>5bp; long reference)
    - long alternate alleles (>5bp; long variant)
    - X chromosome variants 
    - Y chromosome variants

    Calculated with the following code:
    whole_count = ht_whole.count()
    indel_freq = ht_whole.aggregate(hl.agg.count_where(hl.is_indel(ht_whole.alleles[0], ht_whole.alleles[1]))) / whole_count
    ref_freq = ht_whole.aggregate(hl.agg.count_where(hl.len(ht_whole.alleles[0])>5)) / whole_count
    var_freq = ht_whole.aggregate(hl.agg.count_where(hl.len(ht_whole.alleles[1])>5)) / whole_count
    x_freq = ht_whole.aggregate(hl.agg.count_where(ht_whole.locus.contig == "chrX")) / whole_count
    y_freq = ht_whole.aggregate(hl.agg.count_where(ht_whole.locus.contig == "chrY")) / whole_count

    logger.info("Indel Frequency: %f , Long Ref Frequency: %f , Long Var Frequency: %f , X Frequency: %f , Y Frequency: %f", 
                indel_freq, ref_freq, var_freq, x_freq, y_freq)

    # NOTE: this can take a few minutes to run, depending on cluster size and internet connection, so it was ran and numbers included ahead of time down below
    # Feel free to uncomment and run to check the values here for yourself if you would like 
    """

    # Frequencies of the variant types in gnomAD v3.1.2
    # This will need to be updated with different releases or versions
    indel_freq = 0.14525709561723196
    ref_freq = 0.02122722096390106
    var_freq = 0.0381956465303165
    x_freq = 0.0398928560027597
    y_freq = 0.0015374206699161612
    random_freq = 1.00
    whole_count = 759302267
    # e.g.: For 1000 variants, 145 will be indels
    # for 1000 variants, 40 will be on chromosome X and only 1.5 will be on chromosome Y

    # Downsample the matrix table if the user doesn't have the capacity to annotate and go through 759million variants
    if args.downsample <= 1.00:
        ht_whole = ht_whole.sample(args.downsample)
        whole_count *= args.downsample
    else:
        ValueError(
            "Downsample value provided is greater than 1.00, which is not allowed! Float must be below 1.00"
        )

    # Construct Tables for each of the targeted variant types desired for the test subset:
    # Default: 50k random, 10k extra indels, 10k additional long references, 10k  additional long variants, 10k on sex chromosomes (9k on X, 1k on Y)

    # Create freq dictionary with variant type as the key and pre-calculated frequencies as the value
    freq_dict = {
        "indel": indel_freq,
        "long_ref": ref_freq,
        "long_var": var_freq,
        "on_x": x_freq,
        "on_y": y_freq,
        "random": random_freq,
    }

    # Create desired counts dictionary with variant type as the key and desired number of that variant type in the test subset as the value
    # NOTE: would like to incorporate multiallelic and min-rep into this dictionary (and above)
    # BUT since they are generated very differently, it's difficult to work out.
    # UPDATE 02/02: added args <3
    desired_counts = {
        "random": args.n_random,
        "indel": args.n_indel,
        "long_ref": args.n_long_ref,
        "long_var": args.n_long_var,
        "on_x": args.n_xchr,
        "on_y": args.n_ychr,
    }

    logger.info(f"Counts as: {desired_counts}")

    # Annotate the Table with the status of allele we are looking at
    # Random is set to 'True' for all variants because we would like to sample all variants randomly
    ht_whole = ht_whole.annotate(
        indel=hl.is_indel(ht_whole.alleles[0], ht_whole.alleles[1]),
        long_ref=hl.len(ht_whole.alleles[0]) > 5,
        long_var=hl.len(ht_whole.alleles[1]) > 5,
        on_x=(ht_whole.locus.contig == "chrX"),
        on_y=(ht_whole.locus.contig == "chrY"),
        random=True,
    )
    logger.info("Annotation finished")

    # Create Table with Reference, Indel, Long Reference, Long Variant, and Sex Chromosome variants
    hts = []
    for variant_type in freq_dict:
        variant_freq = freq_dict[variant_type]
        desired_count = desired_counts[variant_type]
        sampling_freq = desired_count / (whole_count * variant_freq)

        ht_var = ht_whole.filter(ht_whole[variant_type])
        ht_var = ht_var.sample(sampling_freq)
        ht_var = ht_var.checkpoint(
            f"{args.tmp_dir}/{variant_type}_temp.ht",
            _read_if_exists=read_if_exists,
            overwrite=args.overwrite,
        )
        hts.append(ht_var)
        logger.info(f"Finished creating table for %s", variant_type)

    ht_union = hts[0].union(*hts[1:])

    logger.info(
        f"Total number of variants (including random, indel, long ref, long var, sex chr) in combined Table as: %d",
        ht_union.count(),
    )
    ht_union = ht_union.checkpoint(
        f"{args.tmp_dir}/union_path_temp.ht",
        _read_if_exists=read_if_exists,
        overwrite=args.overwrite,
    )

    # Construct a Table of only multiallelic sites

    # For 3 random partitions, only take the highly multiallelic sites (>5 variants per loci)
    # NOTE: Construct Table only by partitions, and each partition has ~3200 multiallelic variants
    # 3 partitions to yield 9600 multiallelic variants (default, closest to 10,000)

    ht_filt = ht_whole._filter_partitions(
        sample(range(ht_whole.n_partitions()), args.n_multiallelic_partitions)
    )
    group_ht = ht_filt.group_by(ht_filt.locus).aggregate(n=hl.agg.count())
    group_ht = group_ht.filter(group_ht.n > 5)
    n_multiallelic = group_ht.count()

    # Collect all variants at highly multiallelic loci
    # Do collect() and filter() since it is fairly easy just for one field (locus) , instead of an inner join
    multi_list = group_ht.locus.collect(_localize=False)
    ht_multiallelic = ht_filt.filter(multi_list.contains(ht_filt.locus))
    logger.info(
        f"There are %d very multiallelic sites with %d total variants that will be added.",
        n_multiallelic,
        ht_multiallelic.count(),
    )

    # Checkpoint for final multiallelic variants
    ht_multiallelic = ht_multiallelic.checkpoint(
        f"{args.tmp_dir}/multiallelic_path_temp.ht",
        _read_if_exists=read_if_exists,
        overwrite=args.overwrite,
    )

    # Assemble a table of variants altered by minimum representation
    # NOTE: Construct Table only by partition, assume each partition has ~600 minrep variants
    ht_minrep = gnomad_v3_genotypes_vds.vds().variant_data.rows()
    ht_minrep = ht_minrep.annotate(original_alleles=ht_minrep.alleles)
    ht_minrep = hl.split_multi(ht_minrep)
    ht_minrep = ht_minrep.filter(
        ht_minrep.alleles[1] != ht_minrep.original_alleles[ht_minrep.a_index]
    )

    ht_minrep = ht_minrep._filter_partitions(
        sample(range(ht_minrep.n_partitions()), args.n_minrep_partitions)
    )

    # Perform an inner join, which will join on both locus and allele
    ht_minrep = ht_minrep.select()
    ht_minrep = ht_minrep.join(ht_whole, how="inner")

    # Checkpoint the final minrep variants
    ht_minrep = ht_minrep.checkpoint(
        f"{args.tmp_dir}/minrep-checkpoint-v2-checking.ht",
        _read_if_exists=read_if_exists,
        overwrite=args.overwrite,
    )

    # Counts and further assembly
    logger.info(f"There are %d variants in the minrep set", ht_minrep.count())
    ht_final = ht_union.union(ht_multiallelic)
    ht_final = ht_final.union(ht_minrep)

    # Check for duplicates and de-duplication
    org_count = ht_final.count()
    ht_final = ht_final.distinct()
    duplicates_removed = org_count - ht_final.count()
    logger.info(f"Number of duplicates removed: %d", duplicates_removed)

    # Run naive_coalesce() or repartition() to reduce VCF output times
    # QUESTION: if the default value is set to 100, would I even need this 'if' statement? 
    logger.info(f"args.naive_coalesce as: %i",args.naive_coalesce)
    if args.naive_coalesce:
        ht_final = ht_final.naive_coalesce(args.naive_coalesce)

    # Export final Table and VCF (both zipped and unzipped) and append header info for missing FILTER descriptions
    # As noted in the arg parser section, gnomAD includes AC0 and AS_VQSR filters, which are absent in VCF Tools
    ht_final = ht_final.checkpoint(
        f"{args.final_dir}/annotation-testing-vcf-export.ht",
        overwrite=args.overwrite,
        _read_if_exists=False,
    )
    logger.info(f"Total number of variants in VCF: %d", ht_final.count())

    # Export VCF
    final_vcf_path = f"{args.final_dir}/annotation-testing-vcf-export.vcf"
    if args.export_bgz:
        logger.info("Will zip VCF with .bgz extension")
        final_vcf_path += ".bgz"
    else:
        logger.info("Will export without .bgz compression")

    hl.export_vcf(
        ht_final,
        final_vcf_path,
        append_to_header=args.header_fix_path,
    )
    logger.info(f"Finished exporting VCF as: {final_vcf_path}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--google-cloud-project",
        help="Google Cloud Project for billing purposes. Default is broad-mpg-gnomad.",
        default="broad-mpg-gnomad",
    )
    parser.add_argument(
        "--tmp-dir",
        help="Path for storing checkpoints. Default is gs://gnomad-tmp-4day/vrs , which is temp directory set for Hail.",
        default="gs://gnomad-tmp-4day/vrs",
    )
    parser.add_argument(
        "--final-dir",
        help="Final directory for outputting vcf, vcf.bgz, and mt into.",
        default="gs://gnomad-tmp-4day/vrs",
    )
    parser.add_argument(
        "--py-rand-seed",
        help="Random seed for python. Default is 505.",
        default=505,
        type=int,
    )
    parser.add_argument(
        "--hail-rand-seed",
        help="Random seed for hail. Default is 5.",
        default=5,
        type=int,
    )
    parser.add_argument(
        "--n-random",
        help="Number of random variants. Default is 50,000.",
        type=int,
        default=50000,
    )
    parser.add_argument(
        "--n-indel",
        help="Number of additional spiked-in indels. Default is 10,000.",
        type=int,
        default=10000,
    )
    parser.add_argument(
        "--n-long-ref",
        help="Number of long reference allele variants (>5 bp). Default is 10,000.",
        type=int,
        default=10000,
    )
    parser.add_argument(
        "--n-long-var",
        help="Number of long alternate allele variants (>5 bp). Default is 10,000.",
        type=int,
        default=10000,
    )
    parser.add_argument(
        "--n-xchr",
        help="Number of variants on chromosome X. Default is 9,000.",
        type=int,
        default=9000,
    )
    parser.add_argument(
        "--n-ychr",
        help="Number of variants on chromosome Y. Default is 1,000.",
        type=int,
        default=1000,
    )
    parser.add_argument(
        "--n-multiallelic-partitions",
        help="Number of gnomAD v3.1.2 partitions to search for multiallelic variants. Default is 3 partitions, approx 3200 variants/partition.",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--n-minrep-partitions",
        help="Number of gnomAD v3.1.2 partitions to search for variants altered by min-rep. Default is 5 partitions, approx 600 variants/partition.",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--overwrite",
        help="Pass to --overwrite when checkpointing, overwrites when checkpointing, highly recommended.",
        action="store_true",
    )
    parser.add_argument(
        "--header-fix-path",
        help="Path for lines to append to header: gnomAD filters AC0 and AS_VSQR are not defaults in VCF Headers and need to be added.",
        default="gs://gnomad-marten/outputs-and-finals-01-20-23/marten_filter_header_0120.txt",
    )
    parser.add_argument(
        "--downsample",
        help="Proportion to which to downsample the original whole Hail table at the start. Default is 1.00 (no downsampling)",
        type=float,
        default=1.00,
    )
    parser.add_argument(
        "--naive-coalesce",
        help="Integer to pass to naive_coalesce before exporting VCF.",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--export-bgz",
        help="Block gzip output VCF (file will end with .vcf.bgz). If not set, will output uncompressed VCF.",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
