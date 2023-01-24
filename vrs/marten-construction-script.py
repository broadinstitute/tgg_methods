"""
Script that constructs a ~100k Variant Test VCF for GA4GH-VRS Project

Contact: Daniel Marten (marten@broadinstitute.org)
"""

# Import and set up important packages/modules
import hail as hl
import random
from random import sample
from gnomad.resources.grch38.gnomad import public_release

from gnomad.resources.resource_utils import (
    MatrixTableResource,
    VariantDatasetResource,
    VersionedMatrixTableResource,
    VersionedVariantDatasetResource,
)

# Import logging for logs (to replace print statements)
import logging

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_vrs_test_dataset")
logger.setLevel(logging.INFO)

# Set randomness for python package random
random.seed(505)
# Set randomness for Hail Table sampling
hl.init(global_seed=5)

# Import the gnomAD v3.1.2 full variant Table (approx 759,000,000 variants)
# Note: No need to repartition the whole Table, as it was recently repartitioned to 10,000
ht_whole = public_release("genomes").ht()

# Set a temporary path and variables to make writing checkpoints easier (and easier to edit for future users)
tmp_path = 'gs://gnomad-tmp-4day/marten-temp-01-20-pr-v2/'
read_if_var = True
overwrite_var = True

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

logger.info("Indel Frequency: %f , Ref Frequency: %f , Var Frequency: %f , X Frequency: %f , Y Frequency: %f", 
            indel_freq, ref_freq, var_freq, x_freq, y_freq)

print(f'Frequencies as: indels {indel_freq} and long ref {ref_freq} and long var {var_freq} and x freq {x_freq} and y freq {y_freq}')

"""
indel_freq = 0.14525709561723196
ref_freq = 0.02122722096390106
var_freq = 0.0381956465303165 
x_freq = 0.0398928560027597
y_freq = 0.0015374206699161612
whole_count = 759302267
# e.g.: For 1000 variants, 145 will be indels
# for 1000 variants, 40 will be on chromosome X and only 1.5 will be on chromosome Y

# Construct Tables for each of the targeted variant types desired for the test subset:
# 50k random, 10k extra indels, 10k additional long references, 10k  additional long variants, 10k on sex chromosomes (9k on X, 1k on Y)

# Create freq dictionary with variant type as the key and pre-calculated frequencies as the value
freq_dict = {"indel": indel_freq,
    "long_ref": ref_freq,
    "long_var": var_freq,
    "on_x":  x_freq,
    "on_y": y_freq}

# Create desired counts dictionary with variant type as the key and desired number of that variant type in the test subset as the value
desired_counts = {"indel": 10000,
    "long_ref": 10000,
    "long_var": 10000,
    "on_x":  9000,
    "on_y": 1000}

# Note: It is important to define everything after this (like the random samples) since if you dont, fields dont match
ht_whole = ht_whole.annotate(indel = hl.is_indel(ht_whole.alleles[0], ht_whole.alleles[1]),
   long_ref = hl.len(ht_whole.alleles[0]) > 5,
   long_var = hl.len(ht_whole.alleles[1]) > 5,
   on_x = (ht_whole.locus.contig == "chrX"),
   on_y = (ht_whole.locus.contig == "chrY")
   )

# Create Table with 50,000 random variants
rand_variants = 50000
rand_sample_freq = rand_variants/whole_count
ht_union = ht_whole.sample(rand_sample_freq)
ht_union = ht_union.checkpoint(f'{tmp_path}ht_rand_temp.ht', _read_if_exists=read_if_var, overwrite = overwrite_var)

# Append Table with Indel, Long Reference, Long Variant, and Sex Chromosome variants
for variant_type in freq_dict:
    variant_freq = freq_dict[variant_type]
    desired_count = desired_counts[variant_type]
    sampling_freq = desired_count / (whole_count * variant_freq)
    
    ht_var = ht_whole.filter(ht_whole[variant_type])
    ht_var = ht_var.sample(sampling_freq)
    ht_var = ht_var.checkpoint(f'{tmp_path}{variant_type}_temp.ht', _read_if_exists=read_if_var, overwrite=overwrite_var)
    ht_union = ht_union.union(ht_var)

logger.info(f'Total size of random, indel, long ref, long var, sex chr combined Table as: %i', ht_union.count())
ht_union = ht_union.checkpoint(f'{tmp_path}union_path_temp.ht', _read_if_exists=read_if_var, overwrite=overwrite_var)

# Construct a table of only multiallelic sites

# For 3 random partitions, only takes the highly multiallelic sites (>5 variants per loci)
ht_filt = ht_whole._filter_partitions(sample(range(ht_whole.n_partitions()), 3))
group_ht = ht_filt.group_by(ht_filt.locus).aggregate(n=hl.agg.count())
group_ht = group_ht.filter(group_ht.n > 5) 
multiAlleleCount = group_ht.count()
 
# Collect all variants at highly multiallelic loci    
# Do collect() and filter() since it is fairly easy just for one field (locus) , instead of an innter join
multiList = group_ht.locus.collect(_localize = False)
ht_multiallelic = ht_filt.filter(multiList.contains(ht_filt.locus))
logger.info(f'There are %i very multiallelic sites with %i total variants that will be added', multiAlleleCount , ht_multiallelic.count())

# Checkpoint for final multiallelic variants
ht_multiallelic = ht_multiallelic.checkpoint(f'{tmp_path}multiallelic_path_temp.ht', _read_if_exists=read_if_var, overwrite = overwrite_var)

# Construct table of variants altered by min-rep

vds = VersionedVariantDatasetResource(
    '3.1',
    {"3.1": VariantDatasetResource("gs://gnomad/v3.1/raw/gnomad_v3.1.vds")},
)

ht_minrep = vds.vds().variant_data.rows()
ht_minrep = ht_minrep.annotate(original_alleles=ht_minrep.alleles)
ht_minrep = hl.split_multi(ht_minrep)
ht_minrep = ht_minrep.filter(ht_minrep.alleles[1] != ht_minrep.original_alleles[ht_minrep.a_index])
ht_minrep = ht_minrep._filter_partitions(sample(range(ht_minrep.n_partitions()), 5))
ht_minrep = ht_minrep.select()
ht_minrep = ht_minrep.join(ht_whole, how = "inner")
# an inner join here joins on Locus and Allele
# doing the same using a collect() and filter() step would be awkward, since it involves two separate fields
# doing the following 
# ^ results in some few thousand variants, all being samples changed by min-rep

# checkpoint for final minrep variants
ht_minrep = ht_minrep.checkpoint(f'{tmp_path}minrep-checkpoint.ht', _read_if_exists=read_if_var, overwrite = overwrite_var)

# Counts and further assembly
minrepcount = ht_minrep.count()
logger.info(f'There are %i variants in the minrep set', minrepcount)
ht_prefinal = ht_union.union(ht_multiallelic)
ht_final = ht_prefinal.union(ht_minrep)

# Check for duplicates and de-duplication
org_count = ht_final.count()
ht_final = ht_final.distinct()
duplicates_removed = org_count - ht_final.count()
logger.info(f'Number of duplicates removed: %i' , duplicates_removed)

# Export final table and VCF and append header info for missing FILTER descriptions
logger.info(f'Before exporting, there are %i total variants',ht_final.count())
final_path = "gs://gnomad-marten/outputs-and-finals-01-20-23/marten-vcf-01-20-23"
final_path_hail = final_path + ".ht"
final_path_vcf = final_path + ".vcf"
# final_path_vcf_zipped
# ht_final = ht_final.checkpoint(final_path_hail,overwrite = True,append_to_header="gs://gnomad-marten/headerFixv3.txt")
hl.export_vcf(ht_final,final_path_vcf,append_to_header="gs://gnomad-marten/outputs-and-finals-01-20-23/marten_filter_header_0120.txt")
hl.export_vcf(ht_final,f'{final_path_vcf}.bgz',append_to_header="gs://gnomad-marten/outputs-and-finals-01-20-23/marten_filter_header_0120.txt")
ht_final.write(final_path_hail)
