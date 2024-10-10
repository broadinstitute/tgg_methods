"""
Script to query gnomAD site data for variants reported in ClinVar.

This script takes the gnomAD sites data (v4 exomes, v4 genomes) and filters to
variants reported in ClinVar. Additional annotation is added on ClinVar clinical significance, 
gnomAD pext score, and whether the site is outside the capture/calling region of Broad/UKB exomes.

The script generates a table of annotated variants with the following columns:
    - locus
    - alleles
    - AF
    - AC
    - AN
    - filters
    - polyphen
    - SIFT
    - cadd
    - revel
    - most_severe_consequence
    - gene_symbol
    - gene_id
    - clinvar_clinical_significance
    - pext_score
    - outside_broad_capture_region
    - outside_ukb_capture_region
    - outside_broad_calling_region
    - outside_ukb_calling_region

Outputs the annotated variants table as tab-separated values file (.tsv.gz)

Usage:
    - Provide the file path to the ClinVar variant table
        - available from: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
    - Provide the output paths for the TSV files

"""

import hail as hl
hl.init(default_reference = 'GRCh38', log = 'hail.log')

from typing import Dict, List, Tuple, Union
import argparse

# load clinvar data
# available from: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

clinvar_data = hl.import_table('gs://fc-2fe0585f-14c6-49a9-9b69-3cb56e680356/clinvar_variants_in_gnomAD/clinvar_variant_summary_12012023.txt', types={'PositionVCF': hl.tint32})

# annotate locus and alleles

clinvar_data = clinvar_data.annotate(
    locus = hl.locus(clinvar_data.Chromosome, clinvar_data.PositionVCF, reference_genome='GRCh38'),
    alleles = [clinvar_data.ReferenceAlleleVCF, clinvar_data.AlternateAlleleVCF]
)

# key by GRCh37 locus in order to annotate pext score

clinvar_data = clinvar_data.key_by('locus_GRCh37')

# load gnomAD pext scores

pext_data = hl.read_table('gs://gcp-public-data--gnomad/papers/2019-tx-annotation/gnomad_browser/all.baselevel.021620.ht')

pext_data = pext_data.select(
    pext_data.mean_proportion
)
pext_data = pext_data.annotate(
    locus_GRCh37 = hl.str(pext_data.locus)
)

pext_data = pext_data.key_by('locus_GRCh37')

pext_data.describe()

# annotate pext scores

clinvar_data = clinvar_data.annotate(
    pext_score = pext_data[clinvar_data.locus_GRCh37].mean_proportion
)
clinvar_data.describe()

# key by GRCh38 locus and alleles for matching with gnomAD allele frequency data

clinvar_data = clinvar_data.key_by('locus', 'alleles')

# (1) generate results for gnomAD exomes

# load gnomAD v4.1 sites data from exomes

gnomADe_data = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.sites.ht')
gnomADe_data = gnomADe_data.select_globals()

# keep required annotation for clinvar analysis

gnomADe_clinvar = gnomADe_data.select(
    AF = gnomADe_data.freq[0:11].AF,
    AC = gnomADe_data.freq[0:11].AC,
    AN = gnomADe_data.freq[0:11].AN,
    AF_non_ukb = gnomADe_data.freq[167].AF,
    AC_non_ukb = gnomADe_data.freq[167].AC,
    AN_non_ukb = gnomADe_data.freq[167].AN,
    filters = gnomADe_data.filters,
    polphen = gnomADe_data.in_silico_predictors.polyphen_max,
    SIFT = gnomADe_data.in_silico_predictors.sift_max,
    cadd = gnomADe_data.in_silico_predictors.cadd.phred,
    revel = gnomADe_data.in_silico_predictors.revel_max,
    most_severe_consequence = gnomADe_data.vep.most_severe_consequence
)

# annotate with clinvar gene symbol, clinical significance, and pext score

gnomADe_clinvar = gnomADe_clinvar.annotate(
    gene_symbol = clinvar_data[gnomADe_clinvar.key].GeneSymbol,
    gene_id = clinvar_data[gnomADe_clinvar.key].GeneID,
    clinvar_clinical_significance = clinvar_data[gnomADe_clinvar.key].ClinicalSignificance,
    pext_score = clinvar_data[gnomADe_clinvar.key].pext_score,
)

# filter to clinvar variants only

gnomADe_clinvar = gnomADe_clinvar.filter(
    hl.is_defined(gnomADe_clinvar.clinvar_clinical_significance)
)

# annotate with whether the site is captured and called by Broad and UKB exomes

gnomAD_sites = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.allele_number_all_sites.ht')

gnomADe_clinvar = gnomADe_clinvar.key_by('locus')

gnomADe_clinvar = gnomADe_clinvar.annotate(
    outside_broad_capture_region = gnomAD_sites[gnomADe_clinvar.key].outside_broad_capture_region,
    outside_ukb_capture_region = gnomAD_sites[gnomADe_clinvar.key].outside_ukb_capture_region,
    outside_broad_calling_region = gnomAD_sites[gnomADe_clinvar.key].outside_broad_calling_region,
    outside_ukb_calling_region = gnomAD_sites[gnomADe_clinvar.key].outside_ukb_calling_region
)

# export clinvar analysis

gnomADe_clinvar.export('gs://fc-2fe0585f-14c6-49a9-9b69-3cb56e680356/clinvar_variants_in_gnomAD/clinvar_variants_in_gnomAD_v4.1_exomes_12012023.tsv.bgz')

# (2) generate results for gnomAD genomes

# load gnomAD v4.1 sites data from genomes

gnomADg_data = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht')
gnomADg_data = gnomADg_data.select_globals()

# keep required annotation for clinvar analysis

gnomADg_clinvar = gnomADg_data.select(
    AF = gnomADg_data.freq[0:11].AF,
    AC = gnomADg_data.freq[0:11].AC,
    AN = gnomADg_data.freq[0:11].AN,
    filters = gnomADg_data.filters,
    polphen = gnomADg_data.in_silico_predictors.polyphen_max,
    SIFT = gnomADg_data.in_silico_predictors.sift_max,
    cadd = gnomADg_data.in_silico_predictors.cadd.phred,
    revel = gnomADg_data.in_silico_predictors.revel_max,
    most_severe_consequence = gnomADg_data.vep.most_severe_consequence
)

# annotate with clinvar gene symbol, clinical significance, and pext score

gnomADg_clinvar = gnomADg_clinvar.annotate(
    gene_symbol = clinvar_data[gnomADg_clinvar.key].GeneSymbol,
    gene_id = clinvar_data[gnomADg_clinvar.key].GeneID,
    clinvar_clinical_significance = clinvar_data[gnomADg_clinvar.key].ClinicalSignificance,
    pext_score = clinvar_data[gnomADg_clinvar.key].pext_score,
)

# filter to clinvar variants only

gnomADg_clinvar = gnomADg_clinvar.filter(
    hl.is_defined(gnomADg_clinvar.clinvar_clinical_significance)
)

# annotate with whether the site is captured and called by Broad and UKB exomes

gnomAD_sites = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/exomes/gnomad.exomes.v4.1.allele_number_all_sites.ht')

gnomADg_clinvar = gnomADg_clinvar.key_by('locus')

gnomADg_clinvar = gnomADg_clinvar.annotate(
    outside_broad_capture_region = gnomAD_sites[gnomADg_clinvar.key].outside_broad_capture_region,
    outside_ukb_capture_region = gnomAD_sites[gnomADg_clinvar.key].outside_ukb_capture_region,
    outside_broad_calling_region = gnomAD_sites[gnomADg_clinvar.key].outside_broad_calling_region,
    outside_ukb_calling_region = gnomAD_sites[gnomADg_clinvar.key].outside_ukb_calling_region
)

# export clinvar analysis

gnomADg_clinvar.export('gs://fc-2fe0585f-14c6-49a9-9b69-3cb56e680356/clinvar_variants_in_gnomAD/clinvar_variants_in_gnomAD_v4.1_genomes_12012023.tsv.bgz')