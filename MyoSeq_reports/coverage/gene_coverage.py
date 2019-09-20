#!/usr/bin/env python

import gzip
from os.path import basename
import os
import sys

# output should be named gene.coverage.txt
tsv = os.environ['TSV']
out = os.environ['OUT']

# get gene name
gene = basename(tsv).split('.')[0]

# open file for processing
with gzip.open(tsv) as t, open(out, 'w') as o:

    # write header
    #o.write('gene\tchrom\tpos\tmean\tmedian\t1\t5\t10\t15\t20\t25\t30\t50\t100\n')
    for line in t:
        line = line.strip().split('\t')
        chrom = line[0]
        pos = line[1]
        n_samples = len(line)
        total_cov = 0

        mean_cov = float(total_cov) / (n_samples - 2)
        sorted_cov = sorted(map(int, line[2:n_samples]))
        n_samples -= 2
        if n_samples % 2 == 1:
            median_cov = sorted_cov[n_samples/2]
        else:
            median_cov = (sorted_cov[(n_samples/2) - 1] + sorted_cov[(n_samples/2)]) / 2

        o.write('{}\t{}\t{}\t{}\t{}\n'.format(gene, chrom, pos, mean_cov, median_cov))
