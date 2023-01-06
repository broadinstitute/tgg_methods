#!/usr/bin/env python
# coding: utf-8

# Daniel Marten 01-04-2022
# Constructing a ~100k Variant Test Database for GA4GH-VRS Annotation Script
# NOTE: adapted from a Jupyter Notebook, so excuse any formatting quirks! 
# NOTE: please make sure to change the paths for checkpointing!!!
# Edits on: 01-06-2022

import hail as hl
import random
from random import sample
from gnomad.resources.grch38.gnomad import public_release

# Setting randomness
random.seed(505)
hl.init(global_seed=5)

# Import the gnomAD v3.1.2 full variant table (approx 759,000,000 variants)

ht_whole = public_release("genomes").ht()
ht_whole = ht_whole.naive_coalesce(10_000)

# PREVIOUSLY CALCULATED FREQUENCY STATISTICS: 
# Frequency of the following in the general v3 database, based off of frequencies in a randomly sampled
# database of 1mil variants
# for: indels, long references (>5bp), long variants (>5bp), x-chromosome variants, y-chromosome variants
# if you would like, uncomment and run the below code: 

# # Generating Frequencies: 
# ht_1mil = ht_whole.sample(0.0013175) # into 1million samples == 1_000_000/759_000_000 (whole count)
# ht_1mil = ht_1mil.checkpoint('gs://gnomad-tmp-4day/marten-temp-01-06/sample_1mil.ht',overwrite=True)
# count_1mil = ht_1mil.count()
# indel_freq = ht_1mil.filter(hl.is_indel(ht_1mil.alleles[0], ht_1mil.alleles[1])).count() / count_1mil # = 0.14531311994611185
# ref_freq = ht_1mil.filter(hl.len(ht_1mil.alleles[0])>5).count() / count_1mil # = 0.02129941204747994
# var_freq = ht_1mil.filter(hl.len(ht_1mil.alleles[1])>5).count() / count_1mil # = 0.03805950683146261
# x_freq = ht_1mil.filter(ht_1mil.locus.contig == "chrX").count() / count_1mil # = 0.04005632671091385
# y_freq = ht_1mil.filter(ht_1mil.locus.contig == "chrY").count() / count_1mil # = 0.0015570797658583776
# print(f'indel freq as: {indel_freq} \n\tlong ref and var freq as: {ref_freq} , {var_freq} \n\tx and y as: {x_freq} , {y_freq}')

indel_freq = 0.1453
ref_freq = 0.02130 
var_freq = 0.03806 
x_freq = 0.04006 
y_freq = 0.001557 
whole_count = 759_000_000

# e.g.: for 1000 variants, 146 will be indels
# for 1000 variants, 40.06 will be X chromosomes and only 1.557 will be Y chromosome

# Construct tables of samples from this 
# Makeup: 50k random, 10k extra indels, 10k more long references, 10k long variants, 10k sex chromosome
# --> (9k X, 1k Y)

# 50,000 random samples
randSampleFreq = 50000/whole_count
ht_50k = ht_whole.sample(randSampleFreq).repartition(100)
# TO ADD: ht_50k = ht_50k.repartition(100)
rand_path = 'gs://gnomad-tmp-4day/marten-temp-01-06/ht_50k_temp.ht'
ht_50k = ht_50k.checkpoint(rand_path,overwrite = True)

# 10,000 indels
ht_10k_indel = ht_whole.filter(hl.is_indel(ht_whole.alleles[0], ht_whole.alleles[1]))
# ht_whole.filter(ht_whole.vep.variant_class != "SNV")
indel_sample_freq = 10000/(whole_count*indel_freq) # due to repeat samples, looking for 1.5x what we were expecting
ht_10k_indel = ht_10k_indel.sample(indel_sample_freq).repartition(100)
# TO ADD: ht_10k_indel = ht_10k_indel.repartition(100)
indel_path = 'gs://gnomad-tmp-4day/marten-temp-01-06/ht_10k_indel_temp.ht'
ht_10k_indel = ht_10k_indel.checkpoint(indel_path,overwrite = True)

# 10,000 variants with long reference alleles
ht_10k_longRef = ht_whole.filter(hl.len(ht_whole.alleles[0])>5)
longRef_sample_freq = 10000/(whole_count*ref_freq)
ht_10k_longRef = ht_10k_longRef.sample(longRef_sample_freq).repartition(10)
# TO ADD: ht_10k_longRef = ht_10k_longRef.repartition(100)
longRef_path = 'gs://gnomad-tmp-4day/marten-temp-01-06/ht_10k_longRef_temp.ht'
ht_10k_longRef = ht_10k_longRef.checkpoint(longRef_path,overwrite = True)

# 10,000 variants with long variant alleles
ht_10k_longVar = ht_whole.filter(hl.len(ht_whole.alleles[1])>5)
longVar_sample_freq = 10000/(whole_count*var_freq)
ht_10k_longVar = ht_10k_longVar.sample(longVar_sample_freq).repartition(100)
# TO ADD: ht_10k_longVar = ht_10k_longVar.repartition(100)
longVar_path = 'gs://gnomad-tmp-4day/marten-temp-01-06/ht_10k_longVar_temp.ht'
ht_10k_longVar = ht_10k_longVar.checkpoint(longVar_path,overwrite = True)

# 9,000 X-Chromosome Variants
ht_9k_x = ht_whole.filter(ht_whole.locus.contig == "chrX")
x_sample_freq = 9000/(whole_count*x_freq)
ht_9k_x = ht_9k_x.sample(x_sample_freq).repartition(100)
# TO ADD: ht_9k_x = ht_9k_x.repartition(100)
xchr_path = 'gs://gnomad-tmp-4day/marten-temp-01-06/ht_9k_x_temp.ht'
ht_9k_x = ht_9k_x.checkpoint(xchr_path,overwrite = True)

# 1,000 Y-Chromosome Variants 
ht_1k_y = ht_whole.filter(ht_whole.locus.contig == "chrY")
y_sample_freq = 1000/(whole_count*y_freq)
ht_1k_y = ht_1k_y.sample(y_sample_freq)
ht_1k_y = ht_1k_y.repartition(100)
ychr_path = 'gs://gnomad-tmp-4day/marten-temp-01-06/ht_1k_y_temp.ht'
ht_1k_y = ht_1k_y.checkpoint(ychr_path,overwrite = True)

# Combine separate tables into one 

tableParts = [ht_50k, ht_10k_indel, ht_10k_longRef, ht_10k_longVar, ht_9k_x, ht_1k_y]
ht_union = tableParts[0]
for item_to_union in tableParts[1:]:
    ht_union = ht_union.union(item_to_union) 
    # final table of all but the multiallelic sites
    # at this point, approx ~90,000 variants

# cannot figure how to make LOGGER work ? 
print(f'Size of ~50k random, ~10k indel, ~10k long ref, ~10k long var, ~10k sex chr as: {ht_union.count()}')

# Construct a table of only multiallelic sites

# for 3 random partitions, only takes the highly multiallelic sites (>5 variants per loci)
ht_filt = ht_whole._filter_partitions(sample(range(ht_whole.n_partitions()), 3))
group_ht = ht_filt.group_by(ht_filt.locus).aggregate(n=hl.agg.count())
group_ht = group_ht.filter(group_ht.n > 5)
 
multiAlleleCount = group_ht.count()
print(f'There are {multiAlleleCount} very multiallelic SITES in the table above') # 585 sites, plenty of variables
 
# collects all variants at highly multiallelic loci    
multiList = group_ht.locus.collect(_localize = False)
ht_multiallelic = ht_filt.filter(multiList.contains(ht_filt.locus))
# returns 9602 variants on one trial

# checkpoint for final multiallelic variants
multiallelic_path = "gs://gnomad-tmp-4day/marten-temp-01-06/multiallelic.ht"
ht_multiallelic = ht_multiallelic.checkpoint(multiallelic_path, overwrite = True)# , overwrite=overwrite, _read_if_exists=read_if_exists)
# ^^^ TEMPORARY PATH

# Construct table of variants altered by min-rep

# code in this section provided by Julia Goodrich (best boss this side of the Mississippi)
vds = hl.vds.read_vds("gs://gnomad/v3.1/raw/gnomad_v3.1.vds")
ht_minrep = vds.variant_data.rows()
ht_minrep = ht_minrep.annotate(original_alleles=ht_minrep.alleles)
ht_minrep = hl.split_multi(ht_minrep)
ht_minrep = ht_minrep.filter(ht_minrep.alleles[1] != ht_minrep.original_alleles[ht_minrep.a_index])

ht_minrep_filt = ht_minrep._filter_partitions(sample(range(ht_minrep.n_partitions()), 5))
ht_minrep_filt = ht_minrep_filt.select()
ht_minrep_final = ht_minrep_filt.join(ht_whole, how = "inner")
# ^ results in some few thousand variants, all being samples changed by min-rep

# checkpoint for final minrep variants
minrep_path_final = "gs://gnomad-tmp-4day/marten-temp-01-06/minrep-checkpoint.ht"
ht_minrep_final = ht_minrep_final.checkpoint(minrep_path_final, overwrite = True)

minrepcount = ht_minrep_final.count()
print(f'There are {minrepcount} variants in the minrep set')

# Final output and construction steps

ht_prefinal = ht_union.union(ht_multiallelic)
ht_final = ht_prefinal.union(ht_minrep_final)
final_path = "gs://gnomad-tmp/marten-temp-01-06/final_output.ht"
ht_final = ht_final.checkpoint(final_path,overwrite = True)
# this also writes the table to where you want it, so you can set this as the final output!!!

ht_final.show(5)
ht_final.describe()
ht_final.count()

## spiked-in Extra ERRORS Coming Soon !!! 
