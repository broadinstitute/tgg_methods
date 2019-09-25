#!/bin/bash

usage() {
cat << EOF

    This script creates the actual reports (pdf files) for requested MyoSeq samples.

    Assumes pdflatex is installed.

    Inputs:
        -d      Directory with all .tex files

    Outputs:
        pdf files for all samples in same directory as .tex files
EOF
}

# check number of arguments
if [[ $# -lt 2 ]]; then
    usage
    exit 1
fi

# parse args
while getopts "d:o:h" opt; do
    case $opt in
        d)
            dir=$OPTARG
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

for file in ${dir}/*tex; do
    pdflatex ${file}
done
