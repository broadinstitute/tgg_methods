#!/usr/bin/env python
# coding: utf-8

# Daniel Marten 01-04-2022
# Constructing a ~100k Variant Test Database for GA4GH-VRS Annotation Script
# NOTE: adapted from a Jupyter Notebook, so excuse any formatting quirks! 

import hail as hl
from random import sample

# importing the gnomAD v3.1.2 full variant table, approx 750,000,000 variants 

ht_whole = hl.read_table('gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht',_n_partitions=10000)


# PREVIOUSLY CALCULATED FREQUENCY STATISTICS: 
# frequency of the following in the general v3 database , based off of frequencties in random 10k and 100k databases

indel_freq = 0.1460 # based on the amount of non-SNV variants in a 10k and 100k dataset
ref_freq = 0.3704 # fraction of Long References (>3 bp) in a 10k and 100k dataset
var_freq = 0.4650 # fraction of Long Variants (>3 bp) in a 10k and 100k dataset
x_freq = 0.04004 # ^^ for X-Chromosome Variants
y_freq = 0.001531 # ^^ for Y-Chromosome Variants

# Constructing tables of samples from this 
# Makeup: 50k random, 10k extra indels, 10k more long references, 10k long variants, 10k sex chromosome
# --> (9k X, 1k Y)
# + 10k (actual results are 20k) multiallelic variants 

finalSampleFreq = 50000/759000000 # the sampling rate to get 50000 variants of the 759,000,000 read in
final_50k_rand = ht_whole.sample(finalSampleFreq)

final_10k_indel = ht_whole.filter(ht_whole.vep.variant_class != "SNV") # count ~= 111,800,700
indel_sample_freq = 15000/(759000000*indel_freq) # due to repeat samples, looking for 1.5x what we were expecting
final_10k_indel = final_10k_indel.sample(indel_sample_freq)

final_10k_longRef = ht_whole.filter(hl.len(ht_whole.alleles[0])>3)
longRef_sample_freq = 15000/(759000000*ref_freq)
final_10k_longRef = final_10k_longRef.sample(longRef_sample_freq)

final_10k_longVar = ht_whole.filter(hl.len(ht_whole.alleles[1])>3)
longVar_sample_freq = 15000/(759000000*var_freq)
final_10k_longVar = final_10k_longVar.sample(longVar_sample_freq)

final_9k_x = ht_whole.filter(ht_whole.locus.contig == "chrX")
x_sample_freq = 9500/(759000000*x_freq)
final_9k_x = final_9k_x.sample(x_sample_freq)

final_1k_y = ht_whole.filter(ht_whole.locus.contig == "chrY")
y_sample_freq = 1500/(759000000*y_freq)
final_1k_y = final_1k_y.sample(y_sample_freq)

# combining these tables into one final Hail Table to work with 

tableParts = [final_50k_rand, final_10k_indel, final_10k_longRef, final_10k_longVar, final_9k_x, final_1k_y]
final_struct = tableParts[0]
for item_to_union in tableParts[1:]:
    final_struct = final_struct.union(item_to_union) # final table of all but the multiallelic sites

# TEMPORARY PATH BELOW PLEASE USE YOUR OWN PATH 
# final_struct = final_struct.checkpoint("gs://seqr-scratch-temp/marten-temp-01-04/non-multiallelic-checkpoint.ht")
fs = final_struct.count()
print(f'Samples before monoallelic spike-in: {fs}')
# 01-04 @ 5:08 pm: only ~79,000 spiked-in when I wanted 90,000 , probably lots of overlap 
# between random and the extra ones 

# constructing a table of only multiallelic sites

ht1 = ht_whole._filter_partitions(sample(range(ht_whole.n_partitions()), 4)) # reads in 4 partitions (~15k var)
ht1_gb = ht1.group_by(ht1.locus).aggregate(n=hl.agg.count()) # gb for group-by
ht1_filt = ht1_gb.filter(ht1_gb.n > 3) # only the sites for multiallelic sites
# -- 
multiAlleleCount = ht1_filt.count()
print(f'There are {multiAlleleCount} very multiallelic SITES in the table above') # 2292 sites, plenty of variables
# -- 
multiList = ht1_filt.locus.collect(_localize = False)
ht_multiallelic = ht1.filter(multiList.contains(ht1.locus))
# ^ would be awfully slow and expensive IF we weren't working with only 4 partiions

# checkpoint for easier working PLEASE USE YOUR OWN PATH HERE 
# ht_multiallelic = ht_multiallelic.checkpoint("gs://seqr-scratch-temp/marten-temp-01-04/multiallelic-checkpoint-v2.ht")# , overwrite=overwrite, _read_if_exists=read_if_exists)
# ^^^ TEMPORARY PATH

# multiallalic looking/summary

ht_multiallelic.head(n=10).show()
ht_multiallelic.tail(n=10).show()
ht_multiallelic.count()

# outputting! also doing this allows us to count the total size without it melting my computer 

marten_final_output = final_struct.union(ht_multiallelic)

# TEMPORARY PATH PLEASE USE YOUR OWN PATHS 
# marten_final_output = marten_final_output.checkpoint("gs://seqr-scratch-temp/marten-temp-01-04/final-100k-01-04-23.ht")

marten_final_output.head(n=5).show()
marten_final_output.describe()
marten_final_output.count() # 102746 final variants! perfect! fun! nice! 

## spiked-in Extra ERRORS Coming Soon !!! 