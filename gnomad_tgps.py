import hail as hl
from gnomad_qc.v3.resources import basics
import hailtop.fs as hfs

# initialize Hail
hl.init(log='test.log', default_reference='GRCh38', backend='spark')

# get the site table to extract AC and AF frequencies
gnomad_freq_path = "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht"
gnomad_freq_ht = hl.read_table(gnomad_freq_path)

# get gnomad V3 data and split multi-allelic variants
vds = basics.get_gnomad_v3_vds()
vds_mt = vds.variant_data
vds_split = hl.experimental.sparse_split_multi(vds_mt)

# get the metadata of gnomad v3
meta_path = "gs://gnomad/metadata/genomes_v3.1/gnomad_v3.1_sample_qc_metadata.ht"
meta_dat = hl.read_table(meta_path)

# filter out unreleased samples
release_ht = meta_dat.filter(meta_dat.release == True)
vds_split_filterd = vds_split.filter_cols(hl.is_defined(release_ht[vds_split.col_key]))

# load samples in the pangenome
mt1kg = vds_split_filterd
table_1kg = meta_dat

p1_ids = []
id_path = "gs://jialan-tmp-7day/sample_ids_phase2_n55_all.txt"  # read the sample IDs for the reference set
with hfs.open(id_path, 'r') as p1_file:
    for line in p1_file:
        sid = line.strip()
        p1_ids.append(f'{sid}_{sid}')
p1_ids = hl.set(p1_ids)

# Divide the table with genetic data into sample (reference) and out-of-sample populations
mt1kg_oos = mt1kg.filter_cols(p1_ids.contains(mt1kg.s), keep=False)
mt1kg_s = mt1kg.filter_cols(p1_ids.contains(mt1kg.s), keep=True)

# annotate the original matrix tables with AC and AF frequencies
mt1kg_oos = mt1kg_oos.annotate_rows(freq=gnomad_freq_ht[mt1kg_oos.row_key].freq)
mt1kg_s = mt1kg_s.annotate_rows(freq=gnomad_freq_ht[mt1kg_s.row_key].freq)

# In the out-of-sample (OOS) dataset, filter out SNVs with minor allele frequencies lower than 1%
# mt1kg_oos = mt1kg_oos.filter_rows(mt1kg_oos.variant_qc.AF[1] > 0.01)
mt1kg_oos = mt1kg_oos.filter_rows(mt1kg_oos.freq.AF[1] > 0.01)

# In the sample population, filter ouft SNVs with minor allele COUNT lower than 1 -> keep only those SNVs that have a minor allele carrier in the reference group
# mt1kg_s = mt1kg_s.filter_rows(mt1kg_s.variant_qc.AC[1] > 0)
mt1kg_s = mt1kg_s.filter_rows(mt1kg_s.freq.AC[1] > 0)

# Annotate the OOS dataset to count the number of minor alleles carried by each participant
mt1kg_oos = mt1kg_oos.annotate_cols(tn_var_per_sample=hl.agg.count_where(mt1kg_oos.GT.is_non_ref()))
# Create a new matrix with only SNVs that are not included in the reference group
mt1kg_oos_filtered = mt1kg_oos.anti_join_rows(mt1kg_s.rows())
# Annotate the new dataset to count the number of minor alleles carried by each participant
mt1kg_oos_filtered = mt1kg_oos_filtered.annotate_cols(
    n_var_per_sample=hl.agg.count_where(mt1kg_oos_filtered.GT.is_non_ref()))

# Keep only the columns (sample meta data)
mt1kg_oos_filtered_ht = mt1kg_oos_filtered.cols()

# annotate with population data, both continental and subcontinental ancestry
# genetic region?
mt1kg_oos_filtered_ht = mt1kg_oos_filtered_ht.annotate(
    population=table_1kg[mt1kg_oos_filtered_ht.s].population_inference.pop)

ht = mt1kg_oos_filtered_ht
# annotate with fraction of variants carried by each sample that is not included in the reference
ht = ht.annotate(fraction=ht.n_var_per_sample / ht.tn_var_per_sample)

# compute aggretation stats
res = ht.aggregate(hl.struct(ght=hl.agg.group_by(ht.population, hl.agg.mean(ht.n_var_per_sample)),
                             ghtt=hl.agg.group_by(ht.population, hl.agg.mean(ht.tn_var_per_sample)),
                             ghtf=hl.agg.group_by(ht.population, hl.agg.mean(ht.fraction)),
                             ght_stats=hl.agg.group_by(ht.population, hl.agg.stats(ht.n_var_per_sample)),
                             ghtt_stats=hl.agg.group_by(ht.population, hl.agg.stats(ht.fraction)),
                             ghtf_stats=hl.agg.group_by(ht.population, hl.agg.stats(ht.tn_var_per_sample))))

# output directory
tmp_dir = "gs://jialan-tmp-7day/gnomad"

# write to output
with hfs.open(f'{tmp_dir}/stats.txt', 'w') as f:
    print(res, file=f)