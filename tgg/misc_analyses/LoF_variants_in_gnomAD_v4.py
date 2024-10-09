"""
Script to query gnomAD site data for LoF variants within specified genomic regions.

This script takes the gnomAD sites data (v4 exomes, v4 genomes) and filters to
specified genomic regions based on a bed file, followed by filtering to LoF variants 
on any transcript of a gene within those regions.

The script generates a table of annotated variants with the following columns:
    - locus
    - alleles
    - gene_symbol
    - gene_id
    - AF
    - AC
    - AN
    - filters
    - polyphen
    - SIFT
    - cadd
    - revel
    - loftee
    - lof_flags
    - lof_filters
    - lof_info
    - transcript_consequence
    - transcript_id

Outputs the annotated variants table as tab-separated values file (.tsv.gz)

Usage:
    - Provide the file path to the bed file for specified genomic regions
    - Provide the output paths for the TSV files

"""

import hail as hl
hl.init(default_reference = 'GRCh38', log = 'hail.log')

from typing import Dict, List, Tuple, Union
import argparse

# (1) generate results for gnomAD exomes

# load gnomAD v4.1 sites data from exomes

gnomADe_data = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht')
gnomADe_data = gnomADe_data.select_globals()

# load bed file for genomics regions of interest

genelist_hi_intervals = hl.import_bed('gs://fc-2fe0585f-14c6-49a9-9b69-3cb56e680356/clinvar_variants_in_gnomAD/genelist_hi_231130_SG.bed', reference_genome='GRCh38')
gnomADe_hi = hl.filter_intervals(gnomADe_data,genelist_hi_intervals['interval'].collect())

gnomADe_hi = gnomADe_hi.explode(gnomADe_hi.vep.transcript_consequences)

gnomADe_hi = gnomADe_hi.filter(
    (
        (gnomADe_hi.vep.transcript_consequences.consequence_terms.contains('stop_gained')) |
        (gnomADe_hi.vep.transcript_consequences.consequence_terms.contains('frameshift_variant')) |
        (gnomADe_hi.vep.transcript_consequences.consequence_terms.contains('splice_acceptor_variant')) |
        (gnomADe_hi.vep.transcript_consequences.consequence_terms.contains('splice_donor_variant'))
    )
)

# keep required annotation for lof analysis (both HC and LC, all transcripts in genelist_hi)

gnomADe_hi = gnomADe_hi.select(
    gene_symbol = gnomADe_hi.vep.transcript_consequences.gene_symbol,
    gene_id = gnomADe_hi.vep.transcript_consequences.gene_id,
    AF = gnomADe_hi.freq[0:11].AF,
    AC = gnomADe_hi.freq[0:11].AC,
    AN = gnomADe_hi.freq[0:11].AN,
    filters = gnomADe_hi.filters,
    polphen = gnomADe_hi.in_silico_predictors.polyphen_max,
    SIFT = gnomADe_hi.in_silico_predictors.sift_max,
    cadd = gnomADe_hi.in_silico_predictors.cadd.phred,
    revel = gnomADe_hi.in_silico_predictors.revel_max,
    loftee = gnomADe_hi.vep.transcript_consequences.lof,
    lof_flags = gnomADe_hi.vep.transcript_consequences.lof_flags,
    lof_filters = gnomADe_hi.vep.transcript_consequences.lof_filter,
    lof_info = gnomADe_hi.vep.transcript_consequences.lof_info,
    transcript_consequence = gnomADe_hi.vep.transcript_consequences.consequence_terms,
    transcript_id = gnomADe_hi.vep.transcript_consequences.transcript_id
)

# export lof variants in gnomad

gnomADe_hi.export('gs://fc-2fe0585f-14c6-49a9-9b69-3cb56e680356/clinvar_variants_in_gnomAD/genelist_hi_lof_variants_in_gnomAD_v4.1_exomes_12012023.tsv.bgz')

# (2) generate results for gnomAD genomes

# load gnomAD v4.1 sites data from genomes

gnomADg_data = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht')
gnomADg_data = gnomADg_data.select_globals()

# load bed file for genomics regions of interest

genelist_hi_intervals = hl.import_bed('gs://fc-2fe0585f-14c6-49a9-9b69-3cb56e680356/clinvar_variants_in_gnomAD/genelist_hi_231130_SG.bed', reference_genome='GRCh38')
gnomADg_hi = hl.filter_intervals(gnomADg_data,genelist_hi_intervals['interval'].collect())

gnomADg_hi = gnomADg_hi.explode(gnomADg_hi.vep.transcript_consequences)

gnomADg_hi = gnomADg_hi.filter(
    (
        (gnomADg_hi.vep.transcript_consequences.consequence_terms.contains('stop_gained')) |
        (gnomADg_hi.vep.transcript_consequences.consequence_terms.contains('frameshift_variant')) |
        (gnomADg_hi.vep.transcript_consequences.consequence_terms.contains('splice_acceptor_variant')) |
        (gnomADg_hi.vep.transcript_consequences.consequence_terms.contains('splice_donor_variant'))
    )
)

# keep required annotation for lof analysis (both HC and LC, all transcripts in genelist_hi)

gnomADg_hi = gnomADg_hi.select(
    gene_symbol = gnomADg_hi.vep.transcript_consequences.gene_symbol,
    gene_id = gnomADg_hi.vep.transcript_consequences.gene_id,
    AF = gnomADg_hi.freq[0:11].AF,
    AC = gnomADg_hi.freq[0:11].AC,
    AN = gnomADg_hi.freq[0:11].AN,
    filters = gnomADg_hi.filters,
    polphen = gnomADg_hi.in_silico_predictors.polyphen_max,
    SIFT = gnomADg_hi.in_silico_predictors.sift_max,
    cadd = gnomADg_hi.in_silico_predictors.cadd.phred,
    revel = gnomADg_hi.in_silico_predictors.revel_max,
    loftee = gnomADg_hi.vep.transcript_consequences.lof,
    lof_flags = gnomADg_hi.vep.transcript_consequences.lof_flags,
    lof_filters = gnomADg_hi.vep.transcript_consequences.lof_filter,
    lof_info = gnomADg_hi.vep.transcript_consequences.lof_info,
    transcript_consequence = gnomADg_hi.vep.transcript_consequences.consequence_terms,
    transcript_id = gnomADg_hi.vep.transcript_consequences.transcript_id
)

# export lof variants in gnomad

gnomADg_hi.export('gs://fc-2fe0585f-14c6-49a9-9b69-3cb56e680356/clinvar_variants_in_gnomAD/genelist_hi_lof_variants_in_gnomAD_v4.1_genomes_12012023.tsv.bgz')
