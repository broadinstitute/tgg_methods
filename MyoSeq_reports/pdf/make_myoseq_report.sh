#!/bin/bash

usage() {
cat << EOF

    This script called make_myoseq_report.py to create .tex files for all requested MyoSeq patients.
    Use make_pdfs.pl to create pdfs from the .tex files.

    Inputs:
        -l      Sample list
        -c      Create reports with candidate CNVs
        -s      Create reports with SMA carrier status
        -u      Samples have no candidate variants
        -d      Top level directory that contains all data for reports
        -r      Directory containing resource .tex files
        -p      Path to make_myoseq_report.py
        -o      Output directory

    Outputs:
        .tex files for all samples
EOF
}

# check number of arguments
if [[ $# -lt 10 ]]; then
    usage
    exit 1
fi

# parse args
while getopts "l:d:r:o:p:ucsh" opt; do
    case $opt in
        l)
            list=$OPTARG
        ;;
        c)
            cnv=true
        ;;
        s)
            sma=true
        ;;
        u)
            unsolved=true
        ;;
        d)
            dir=$OPTARG
        ;;
        r)
            resources=$OPTARG
        ;;
        o)
            out=$OPTARG
        ;;
        p)
            script_path=$OPTARG
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

while read line; do 
    if [ $cnv ]; then
        if [ $unsolved ]; then
            python $script_path/make_myoseq_report.py -p ${line} -d ${dir} -r ${resources} -o ${out} -c -u
        else
            python $script_path/make_myoseq_report.py -p ${line} -d ${dir} -r ${resources} -o ${out} -c
        fi
    elif [ $sma ]; then
         python $script_path/make_myoseq_report.py -p ${line} -d ${dir} -r ${resources} -o ${out} -s
    elif [ $unsolved ]; then
        python $script_path/make_myoseq_report.py -p ${line} -d ${dir} -r ${resources} -o ${out} -u
    else
        python $script_path/make_myoseq_report.py -p ${line} -d ${dir} -r ${resources} -o ${out}
    fi
done < ${list} 
