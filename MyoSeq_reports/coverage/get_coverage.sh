#!/bin/bash

usage() {
cat << EOF

    This program calculates coverage across the MyoSeq gene list in a set of samples..

    Assumptions:
        BEDfiles exist for all MyoSeq genes.

    Inputs:
        -d      Directory containing BEDfile for all MyoSeq genes
        -c      List of crams for the input samples (one cram per line)
        -f      Reference FASTA
        -o      Output directory for bgzipped coverage TSVs

    Outputs:
        Coverage tsvs for each gene, where each column is a sample (in the same order as in crams list)
EOF
}

# check number of arguments
if [[ $# -lt 8 ]]; then
    usage
    exit 1
fi

# parse args
while getopts "d:c:f:o:h" opt; do
    case $opt in
        d)
            dir=$OPTARG
        ;;
        c)
            crams=$OPTARG
        ;;
        f)
            fasta=$OPTARG
        ;;
        o)
            out=$OPTARG
        ;;
        h)
            usage
            exit 0
        ;;
        \?)
            usage
            exit 0
        ;;
    esac
done

# index each cram if its index doesn't exist
while read line; do
    crai="${line}.crai"
    if [[ ! -s $crai ]]; then
        samtools index -c $line > $crai
    fi
done < ${crams}

# calculate coverage for each gene
for file in ${dir}/*bed; do

    gene=$(basename $file | cut -d '.' -f1)
    len=$(awk 'END {print NR}' $file)
    if (( $len > 1 )); then
        chrom=$(head -1 $file | cut -f1)
        first=$(cat $file | cut -f2 | sort -n | head -1)
        last=$(cat $file | cut -f3 | sort -nr | head -1)
    else
        chrom=$(cat $file | cut -f1)
        first=$(cat $file | cut -f2)
        last=$(cat $file | cut -f3)
    fi

    region="${chrom}:${first}-${last}"
    outfile="${out}/${gene}.tsv.bgz"

    echo "samtools depth -b ${file} -r ${region} -q 10 -Q 20 -a -f ${crams} --reference ${fasta} | bgzip > ${outfile}"
    samtools depth -b ${file} -r ${region} -q 10 -Q 20 -a -f ${crams} --reference ${fasta} | bgzip > ${outfile}
done
