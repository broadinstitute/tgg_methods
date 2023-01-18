#!/usr/bin/env python
# coding: utf-8

# Daniel Marten 01-18-2022
# Construct a ~100k Variant Test Database for GA4GH-VRS Annotation Script
# ran in broad-mpg-gnomad/gnomad_marten as marten-assembly-v3.ipynb
# major changes from last version: used Black to format, added temporary and final path variables, 
# added variables for number of variants for each type, de-duplicated, and fixed the header!

import hail as hl
import random
from random import sample
from gnomad.resources.grch38.gnomad import public_release
import logging

# Setting randomness
random.seed(505)
hl.init(global_seed=5)

# Import the gnomAD v3.1.2 full variant table (approx 759,000,000 variants)
# note: no need to repartition the whole table, as it was recently repartitioned to 10,000

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

print(f'Frequencies as: indels {indel_freq} and long ref {ref_freq} and long var {var_freq} and x freq {x_freq} and y freq {y_freq}')

Quite time-intensive, so ran previously to save computation headache

Feel free to run yourself - at your own peril ! 
"""

indel_freq = 0.14525709561723196
ref_freq = 0.02122722096390106
var_freq = 0.0381956465303165
x_freq = 0.0398928560027597
y_freq = 0.0015374206699161612
whole_count = 759302267

# e.g.: for 1000 variants, 145 will be indels
# for 1000 variants, 40 will be X chromosomes and only 1.5 will be Y chromosome

# Set a temporary path to make writing checkpoints easier (and more easy to edit for future users)
tmp_path = "gs://gnomad-tmp-4day/marten-temp-01-18-v3/"

# Construct tables of samples from this
# Makeup: 50k random, 10k extra indels, 10k more long references, 10k long variants, 10k sex chromosome
# --> (9k X, 1k Y)

# 50,000 random variants
rand_variants = 50000
randSampleFreq = rand_variants / whole_count
ht_rand = ht_whole.sample(randSampleFreq)
# ht_rand = ht_rand.naive_coalesce(100)
# TO ADD: ht_50k = ht_50k.repartition(100)
rand_path = tmp_path + "ht_rand_temp.ht"
ht_rand = ht_rand.checkpoint(rand_path, _read_if_exists=True, overwrite=True)

# 10,000 indels
ht_indel = ht_whole.filter(hl.is_indel(ht_whole.alleles[0], ht_whole.alleles[1]))
# ht_whole.filter(ht_whole.vep.variant_class != "SNV")
indel_variants = 10000
indel_sample_freq = indel_variants / (whole_count * indel_freq)
ht_indel = ht_indel.sample(indel_sample_freq)
# ht_indel = ht_indel.naive_coalesce(100)
indel_path = tmp_path + "ht_indel_temp.ht"
ht_indel = ht_indel.checkpoint(indel_path, _read_if_exists=True, overwrite=True)

# 10,000 variants with long reference alleles
ht_longRef = ht_whole.filter(hl.len(ht_whole.alleles[0]) > 5)
longRef_variants = 10000
longRef_sample_freq = longRef_variants / (whole_count * ref_freq)
ht_longRef = ht_longRef.sample(longRef_sample_freq)
# ht_longRef = ht_longRef.naive_coalesce(100)
longRef_path = tmp_path + "ht_longRef_temp.ht"
ht_longRef = ht_longRef.checkpoint(longRef_path, _read_if_exists=True, overwrite=True)

# 10,000 variants with long variant alleles
ht_longVar = ht_whole.filter(hl.len(ht_whole.alleles[1]) > 5)
longVar_variants = 10000
longVar_sample_freq = longVar_variants / (whole_count * var_freq)
ht_longVar = ht_longVar.sample(longVar_sample_freq)
# ht_longVar = ht_longVar.naive_coalesce(100) # was changed to 50, took much longer than the above code
longVar_path = tmp_path + "ht_longVar_temp.ht"
ht_longVar = ht_longVar.checkpoint(longVar_path, _read_if_exists=True, overwrite=True)

# 9,000 X-Chromosome Variants
ht_x = ht_whole.filter(ht_whole.locus.contig == "chrX")
x_variants = 9000
x_sample_freq = x_variants / (whole_count * x_freq)
ht_x = ht_x.sample(x_sample_freq)  # .repartition(100)
# ht_x = ht_x.naive_coalesce(100)
# TO ADD: ht_9k_x = ht_9k_x.repartition(100)
xchr_path = tmp_path + "ht_x_temp.ht"
ht_x = ht_x.checkpoint(xchr_path, _read_if_exists=True, overwrite=True)

# 1,000 Y-Chromosome Variants
ht_y = ht_whole.filter(ht_whole.locus.contig == "chrY")
y_variants = 1000
y_sample_freq = y_variants / (whole_count * y_freq)
ht_y = ht_y.sample(y_sample_freq)
# ht_y = ht_y.naive_coalesce(10)
ychr_path = tmp_path + "ht_y_temp.ht"
ht_y = ht_y.checkpoint(ychr_path, _read_if_exists=True, overwrite=True)

# Combine separate tables into one # 01-17: after running, change ht_50k to ht_rand
tableParts = [ht_rand, ht_indel, ht_longRef, ht_longVar, ht_x, ht_y]
ht_union = tableParts[0]
for item_to_union in tableParts[1:]:
    ht_union = ht_union.union(item_to_union)
    # final table of all but the multiallelic sites
    # at this point, approx ~90,000 variants

print(f"Size of random, indel, long ref, long var, sex chr as: {ht_union.count()}")

ht_union_path = tmp_path + "union_path_temp.ht"
ht_union = ht_union.checkpoint(ht_union_path, _read_if_exists=True, overwrite=True)


# Construct a table of only multiallelic sites

# for 3 random partitions, only takes the highly multiallelic sites (>5 variants per loci)
ht_filt = ht_whole._filter_partitions(sample(range(ht_whole.n_partitions()), 3))
group_ht = ht_filt.group_by(ht_filt.locus).aggregate(n=hl.agg.count())
group_ht = group_ht.filter(group_ht.n > 5)

multiAlleleCount = group_ht.count()
print(
    f"There are {multiAlleleCount} very multiallelic SITES in the table above"
)  # 585 sites, plenty of variables

# collects all variants at highly multiallelic loci
# doing collect() and filter() since it is fairly easy just for one field (locus)
multiList = group_ht.locus.collect(_localize=False)
ht_multiallelic = ht_filt.filter(multiList.contains(ht_filt.locus))
# returns 9602 variants on one trial

# checkpoint for final multiallelic variants
multiallelic_path = tmp_path + "multiallelic.ht"
ht_multiallelic = ht_multiallelic.checkpoint(
    multiallelic_path, _read_if_exists=True, overwrite=True
)  # , overwrite=overwrite, _read_if_exists=read_if_exists)

# Construct table of variants altered by min-rep
# code in this section provided by Julia Goodrich (best boss this side of the Mississippi)
vds = hl.vds.read_vds("gs://gnomad/v3.1/raw/gnomad_v3.1.vds")
ht_minrep = vds.variant_data.rows()
# ht_vds = hl.vds.read_vds(vds_path).variant_data.rows() # never told path, revert back to original
ht_minrep = ht_minrep.annotate(original_alleles=ht_minrep.alleles)
ht_minrep = hl.split_multi(ht_minrep)
ht_minrep = ht_minrep.filter(
    ht_minrep.alleles[1] != ht_minrep.original_alleles[ht_minrep.a_index]
)
ht_minrep = ht_minrep._filter_partitions(sample(range(ht_minrep.n_partitions()), 5))
ht_minrep = ht_minrep.select()
ht_minrep = ht_minrep.join(ht_whole, how="inner")
# an inner join here joins on Locus and Allele
# doing the same using a collect() and filter() step would be awkward, since it involves two separate fields
# doing the following
# ^ results in some few thousand variants, all being samples changed by min-rep

# checkpoint for final minrep variants
minrep_path_final = tmp_path + "minrep-checkpoint.ht"
ht_minrep = ht_minrep.checkpoint(
    minrep_path_final, _read_if_exists=True, overwrite=True
)

# Counts and further assembly
minrepcount = ht_minrep.count()
print(f"There are {minrepcount} variants in the minrep set")

ht_prefinal = ht_union.union(ht_multiallelic)
ht_final = ht_prefinal.union(ht_minrep)

# Check for duplicates and de-duplication

# Simplified de-duplication
final_count = ht_final.count()
ht_final_distinct = ht_final.distinct()
duplicates_removed = final_count - ht_final_distinct.count()
print(f"Number of duplicates removed: {duplicates_removed}")
ht_final = ht_final_distinct

"""
# de-duplication with visuals (visuals for first tenth of dataset, running on full dataset takes prohibitively long)
ht_final_select = ht_final.select().head(10000)
ht_final_annotate = ht_final_select.annotate(variant=hl.format('%s-%s-%s-%s', ht_final_select.locus.contig, ht_final_select.locus.position, ht_final_select.alleles[0], ht_final_select.alleles[1]))
var_counter = ht_final_annotate.aggregate(hl.agg.counter(ht_final_annotate.variant))

duplicate_count = 0
for var in var_counter:
    if var_counter[var] > 1:
        print(f'{var} present in HT more than once')
        duplicate_count += 1  
        
# Show duplication, filter, show same site
ht_final.filter((ht_final.locus.contig == "chr1") & (ht_final.locus.position == 246310002)).show()
ht_final_distinct_visual = ht_final.distinct()
ht_final_distinct_visual.filter((ht_final_distinct_visual.locus.contig == "chr1") & (ht_final_distinct_visual.locus.position == 246310002)).show()
duplicates_removed = final_count - ht_final_distinct_visual.count()
print(f'Number of duplicates removed: {duplicates_removed}')
"""

# Final Export
final_path = "gs://gnomad-marten/marten-dataset-01-18-23"
final_path_hail = final_path + ".ht"
final_path_vcf = final_path + ".vcf.bgz"
# ht_final = ht_final.checkpoint(final_path_hail,overwrite = True,append_to_header="gs://gnomad-marten/headerFixv3.txt")
hl.export_vcf(
    ht_final,
    final_path_vcf,
    append_to_header="gs://gnomad-marten/marten_filter_header_0118.txt",
)
