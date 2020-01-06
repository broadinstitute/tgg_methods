
# How to generate BED files for MyoSeq gene list

## Part 1: Download gene list
Download the MyoSeq gene list from [_seqr_](https://seqr.broadinstitute.org/dashboard). Extract the second column (gene names).
`cat myoseq_gene_list_genes.tsv | sed 's/"//g' | cut -f2 | grep -v "Symbol" > myoseq_genes.txt`

## Part 1: Download from Gencode

 1. Download **Basic gene annotation** gtf from [Gencode]([https://www.gencodegenes.org/human/]).
 2. Grab all lines with **gene** from gtf: `awk '{if($3=="gene"){print $0}}' gencode.v32.basic.annotation.gtf > gencode.v32.basic.annotation.genes.gtf`
 3. Extract only the lines with relevant genes:
