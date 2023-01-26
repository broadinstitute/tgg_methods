
"""
Script that constructs a ~100k Variant Test VCF for GA4GH-VRS Project

Contact: Daniel Marten (marten@broadinstitute.org)
"""
# NOTE: Ran isort then black at EOD Thursday 01/26

import argparse
import logging
import random
from random import sample

# Import and set up important packages/modules
import hail as hl
from gnomad.resources.grch38.gnomad import public_release
from gnomad.resources.resource_utils import (
    MatrixTableResource,
    VariantDatasetResource,
    VersionedMatrixTableResource,
    VersionedVariantDatasetResource,
)
from gnomad_qc.v3.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v3_vds,
    gnomad_v3_genotypes_vds,
    qc_temp_prefix,
)
from gnomad_qc.v3.resources.constants import CURRENT_RELEASE, CURRENT_VERSION
from gnomad_qc.v3.resources.release import (
    append_to_vcf_header_path,
    hgdp_tgp_subset,
    release_header_path,
    release_sites,
    release_vcf_path,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_vrs_test_dataset")
logger.setLevel(logging.INFO)

TMP_DIR = "gs://gnomad-tmp-4day/vrs"
"""
Bucket to store temporary data.
"""


def main(args):

    # Set randomness for Python package "random"
    random.seed(args.py_rand_seed)
    # Initialize Hail with a seed passed as an argument
    hl.init(
        default_reference="GRCh38",
        global_seed=args.hail_rand_seed,
        tmp_dir=TMP_DIR,
        quiet=True,
        spark_conf={
            "spark.hadoop.fs.gs.requester.pays.mode": "AUTO",
            "spark.hadoop.fs.gs.requester.pays.project.id": "broad-mpg-gnomad",
        },
    )

    # Import the gnomAD v3.1.2 full variant Table (approx 759,000,000 variants)
    # NOTE: No need to repartition the whole Table, as it was recently repartitioned to 10,000
    ht_whole = public_release("genomes").ht()

    """
    PREVIOUSLY CALCULATED FREQUENCY STATISTICS

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
    # This may need to be changed with future releases or versions
    indel_freq = 0.14525709561723196
    ref_freq = 0.02122722096390106
    var_freq = 0.0381956465303165
    x_freq = 0.0398928560027597
    y_freq = 0.0015374206699161612
    random_freq = 1.00
    whole_count = 759302267
    # e.g.: For 1000 variants, 145 will be indels
    # for 1000 variants, 40 will be on chromosome X and only 1.5 will be on chromosome Y

    # Downsample the matrix table if the User doesn't have the capacity to annotate and go through 759million variants
    if args.downsample is float and args.downsample < 1.00:
        ht_whole = ht_whole.sample(args.downsample)
        whole_count *= args.downsample

    # Construct Tables for each of the targeted variant types desired for the test subset:
    # DEFAULT: 50k random, 10k extra indels, 10k additional long references, 10k  additional long variants, 10k on sex chromosomes (9k on X, 1k on Y)

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
    # NOTE: would like to incorporate Multiallelic and Min-Rep into this dictionary (and above)
    # BUT since they are generated very differently, it's difficult to work out.
    desired_counts = {
        "random": 50000,
        "indel": 10000,
        "long_ref": 10000,
        "long_var": 10000,
        "on_x": 9000,
        "on_y": 1000,
    }

    # Reset the desired counts if the user chooses to change them - Multiallelic stored at [-2], Min-Rep at [-1]
    counts_list = []
    if args.custom:
        counts_list = []
        if args.custom:
            argument_list = args.counts.split(",")
            for argo in argument_list:
                counts_list.append(argo)

        for obj in enumerate(desired_counts.keys()):
            desired_counts[obj[1]] = int(counts_list[obj[0]])

    logger.info(f"Counts as: {desired_counts}")

    # Annotate the Table with the status of Allele we are looking at
    ht_whole = ht_whole.annotate(
        indel=hl.is_indel(ht_whole.alleles[0], ht_whole.alleles[1]),
        long_ref=hl.len(ht_whole.alleles[0]) > 5,
        long_var=hl.len(ht_whole.alleles[1]) > 5,
        on_x=(ht_whole.locus.contig == "chrX"),
        on_y=(ht_whole.locus.contig == "chrY"),
        random=(True),
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
            _read_if_exists=args.read_if_exists,
            overwrite=args.overwrite,
        )
        hts.append(ht_var)
        logger.info(f"Finished creating table for %s", variant_type)

    ht_union = hts[0].union(*hts[1:])

    logger.info(
        f"Total size of random, indel, long ref, long var, sex chr combined Table as: %i",
        ht_union.count(),
    )
    ht_union = ht_union.checkpoint(
        f"{args.tmp_dir}/union_path_temp.ht",
        _read_if_exists=args.read_if_exists,
        overwrite=args.overwrite,
    )

    # Construct a Table of only multiallelic sites

    # For 3 random partitions, only take the highly multiallelic sites (>5 variants per loci)
    # NOTE: Construct Table only by partitions, and each partition has ~3200 multiallelic variants
    # 3 partitions to yield 9600 multiallelic variants (default, closest to 10,000)
    multiallelic_partitions = 3

    if args.custom:
        multiallelic_partitions = int(counts_list[-2]) // 3200
        if multiallelic_partitions < 1 and counts_list[-2] != 0:
            multiallelic_partitions = 1
        logger.info(
            f"Number of multiallelic partions to search: %i", multiallelic_partitions
        )

    ht_filt = ht_whole._filter_partitions(
        sample(range(ht_whole.n_partitions()), multiallelic_partitions)
    )
    group_ht = ht_filt.group_by(ht_filt.locus).aggregate(n=hl.agg.count())
    group_ht = group_ht.filter(group_ht.n > 5)
    n_multiallelic = group_ht.count()

    # Collect all variants at highly multiallelic loci
    # Do collect() and filter() since it is fairly easy just for one field (locus) , instead of an inner join
    multi_list = group_ht.locus.collect(_localize=False)
    ht_multiallelic = ht_filt.filter(multi_list.contains(ht_filt.locus))
    logger.info(
        f"There are %i very multiallelic sites with %i total variants that will be added",
        n_multiallelic,
        ht_multiallelic.count(),
    )

    # Checkpoint for final multiallelic variants
    ht_multiallelic = ht_multiallelic.checkpoint(
        f"{args.tmp_dir}/multiallelic_path_temp.ht",
        _read_if_exists=args.read_if_exists,
        overwrite=args.overwrite,
    )

    # Assemble a table of variants altered by Minimum Representation
    # NOTE: Construct Table only by partition, assume each partition has ~600 minrep variants
    ht_minrep = gnomad_v3_genotypes_vds.vds().variant_data.rows()
    ht_minrep = ht_minrep.annotate(original_alleles=ht_minrep.alleles)
    ht_minrep = hl.split_multi(ht_minrep)
    ht_minrep = ht_minrep.filter(
        ht_minrep.alleles[1] != ht_minrep.original_alleles[ht_minrep.a_index]
    )

    # 5 partitions to yield 3000 min-rep variants, the default
    minrep_partitions = 5
    if args.custom:
        minrep_partitions = int(counts_list[-1]) // 600
        if minrep_partitions < 1 and counts_list[-1] != 0:
            minrep_partitions = 1
    logger.info(f"Partitions to search for min-rep: %i", minrep_partitions)

    ht_minrep = ht_minrep._filter_partitions(
        sample(range(ht_minrep.n_partitions()), minrep_partitions)
    )
    ht_minrep = ht_minrep.select()
    ht_minrep = ht_minrep.join(ht_whole, how="inner")
    # an inner join here joins on Locus and Allele
    # doing the same using a collect() and filter() step would be awkward, since it involves two separate fields

    # checkpoint for final minrep variants
    ht_minrep = ht_minrep.checkpoint(
        f"{args.tmp_dir}/minrep-checkpoint-v2-checking.ht",
        _read_if_exists=args.read_if_exists,
        overwrite=args.overwrite,
    )

    # Counts and further assembly
    logger.info(f"There are %i variants in the minrep set", ht_minrep.count())
    ht_final = ht_union.union(ht_multiallelic)
    ht_final = ht_final.union(ht_minrep)

    # Check for duplicates and de-duplication
    org_count = ht_final.count()
    ht_final = ht_final.distinct()
    duplicates_removed = org_count - ht_final.count()
    logger.info(f"Number of duplicates removed: %i", duplicates_removed)

    # Export final Table and VCF (both zipped and unzipped) and append header info for missing FILTER descriptions
    logger.info(f"Before exporting, there are %i total variants", ht_final.count())
    ht_final = ht_final.checkpoint(
        f"{args.final_path}/testing-vcf-export.ht",
        overwrite=args.overwrite,
        _read_if_exists=args.read_if_exists,
    )
    logger.info("VCF as Hail Table written to final path")
    if args.export_vcf:
        hl.export_vcf(
            ht_final,
            f"{args.final_path}/testing-vcf-export.vcf",
            append_to_header=args.header_fix_path,
        )
        logger.info("VCF as Unzipped VCF written to final path")
    if args.export_bgz:
        hl.export_vcf(
            ht_final,
            f"{args.final_path}/testing-vcf-export.vcf.bgz",
            append_to_header=args.header_fix_path,
        )
        logger.info("VCF as Zipped VCF.BGZ written to final path")

    logger.info("Script and output is finished!")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tmp-dir",
        help="Path for storing checkpoints, default to gs://gnomad-tmp-4day/vrs",
        default=TMP_DIR,
    )
    parser.add_argument(
        "--py_rand_seed", help="Random seed for python, def: 505", default=505, type=int
    )
    parser.add_argument(
        "--hail_rand_seed", help="Random seed for hail, def: 5", default=5, type=int
    )
    parser.add_argument(
        "--final_path",
        help="Final path for outputting vcf, vcf.bgz, and mt",
        default="gs://gnomad-tmp-4day/vrs",
    )
    # NOTE: Lots of room for improvement for Custom Counts handling.
    # -> I thought about separat args for rand_count, indel_count, long_ref_count, ... but thought it would be too wordy
    parser.add_argument(
        "--custom",
        help="Set if you would like to include custom counts",
        action="store_true",
    )
    parser.add_argument(
        "--counts",
        help="COMMA SEPARATED LIST for: random, indel, long ref, long var, chrX, chrY, multiallelic sites, and min-rep altered variants",
        type=str,
    )
    parser.add_argument(
        "--read_if_exists",
        help="Pass to --read_if_exists when checkpointing, not recommended",
        action="store_true",
    )  # defaults to False, recommended not to include!
    parser.add_argument(
        "--overwrite",
        help="Pass to --overwrite when checkpointing, overwrites when checkpointing, highly recommended",
        action="store_true",
    )  # defaults to False, highly recommended to include!
    parser.add_argument(
        "--header_fix_path",
        help="Path for files to append to header",
        default="gs://gnomad-marten/outputs-and-finals-01-20-23/marten_filter_header_0120.txt",
    )
    parser.add_argument(
        "--downsample",
        help="Downsample the original whole Hail table at the start",
        type=float,
        default=1.00,
    )
    parser.add_argument(
        "--export_vcf",
        help="Export resulting VCF as a VCF with file ending .vcf",
        action="store_true",
    )
    parser.add_argument(
        "--export_bgz",
        help="Export resulting VCF as BGZ-zipped with file ending .vcf.bgz",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
