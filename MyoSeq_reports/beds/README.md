
# How to generate BED files for MyoSeq gene list
All BED files currently exist in this repo under `resources/beds/`. This README describes how the files were made, and how to update the files with future Gencode releases.

## Part 1: Download gene list
Download the MyoSeq gene list from [_seqr_](https://seqr.broadinstitute.org/dashboard). Extract the second column (gene names).
`cat myoseq_gene_list_genes.tsv | sed 's/"//g' | cut -f2 | grep -v "Symbol" > myoseq_genes.txt`

## Part 2: Download from Gencode and parse Gencode download

 1. Download **Basic gene annotation** gtf from [Gencode]([https://www.gencodegenes.org/human/]).
 2. Grab all lines with **CDS** from gtf: `awk '{if($3=="CDS"){print $0}}' gencode.v32.basic.annotation.gtf > gencode.v32.basic.annotation.CDS.gtf`
 3. Extract only the lines with relevant genes:
    `python get_myoseq_gene_boundaries.py -l myoseq_genes.txt -g gencode.v32.basic.annotation.genes.gtf -o ../resources/beds/`
