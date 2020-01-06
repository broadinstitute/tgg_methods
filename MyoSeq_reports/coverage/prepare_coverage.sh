#!/bin/bash

usage() {
cat << EOF

    This program prepares a batch of MyoSeq samples for coverage calculations across the MyoSeq gene list.

    Assumptions:
        BEDfiles exist for all MyoSeq list genes.

    Inputs:
        -d      Directory containing BEDfile for all MyoSeq list genes
        -c      List of crams for the input samples (one cram per line)
        -f      Reference FASTA
        -b      Output batch file for dsub
        -o      Output bucket for coverage files

    Outputs:
        Batch file for dsub.
EOF
}

# check number of arguments
if [[ $# -lt 10 ]]; then
    usage
    exit 1
fi

# parse args
while getopts "d:c:f:b:o:h" opt; do
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
        b)
            batchfile=$OPTARG
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

# write header of file
echo -e "--input BAMLIST\t--env REGION\t--input FASTA\t--output OUT" > ${batchfile}

# prepare files for dsub
for file in ${dir}/*bed; do

    gene=$(basename $file | cut -d '.' -f1)
    len=$(awk 'END {print NR}' $file)
    if (( $len > 1 )); then
        chrom=$(head -1 $file | cut -f1)
        first=$(head -1 $file | cut -f2)
        last=$(tail -1 $file | cut -f2)
    else
        chrom=$(cat $file | cut -f1)
        first=$(cat $file | cut -f2)
        last=$(cat $file | cut -f3)
    fi

    region="${chrom}:${first}-${last}"
    outfile="${out}/${gene}.tsv.bgz"
    echo -e "${crams}\t${region}\t${fasta}\t${outfile}" >> ${batchfile}
done 
