ClinVar_variants_in_gnomAD_v4.py

Description: This script takes the gnomAD sites data (v4 exomes, v4 genomes) and filters to variants reported in ClinVar. 
Additional annotation is added on ClinVar clinical significance, gnomAD pext score, and whether the site is outside the capture/calling region of Broad/UKB exomes.
Associated study: Gudmundsson et al., medRxiv 2024: https://www.biorxiv.org/content/10.1101/2024.06.12.593113v1.full

LoF_variants_in_gnomAD_v4.py

Description: This script takes the gnomAD sites data (v4 exomes, v4 genomes) and filters to specified genomic regions based on a bed file, 
followed by filtering to LoF variants on any transcript of a gene within those regions.
Associated study: Gudmundsson et al., medRxiv 2024: https://www.biorxiv.org/content/10.1101/2024.06.12.593113v1.full
