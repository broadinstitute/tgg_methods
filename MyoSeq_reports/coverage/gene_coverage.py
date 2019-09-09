#!/usr/bin/env python

import gzip
from os.path import basename
import os
import sys

tsv = os.environ['TSV']
out = os.environ['OUT']

# get gene name
gene = basename(tsv).split('.')[0]

# set coverage thresholds
thresholds = [1, 5, 10, 15, 20, 25, 30, 50, 100]
with gzip.open(tsv) as t, open(out, 'w') as o:

    # write header
    #o.write('gene\tchrom\tpos\tmean\tmedian\t1\t5\t10\t15\t20\t25\t30\t50\t100\n')
    for line in t:
        line = line.strip().split('\t')
        chrom = line[0]
        pos = line[1]
        n_samples = len(line)
        total_cov = 0
        threshold_counts = [0, 0, 0, 0, 0, 0, 0, 0, 0]

        # check if sample's coverage exceeds thresholds; if so, increment threshold_counts
        for i in xrange(2, n_samples):
            samp_cov = int(line[i])
            total_cov += samp_cov
            for j in xrange(len(thresholds)):
                if samp_cov >= thresholds[j]:
                    threshold_counts[j] += 1

        # get mean and median cov
        mean_cov = float(total_cov) / (n_samples - 2)
        sorted_cov = sorted(map(int, line[2:n_samples]))
        n_samples -= 2
        if n_samples % 2 == 1:
            median_cov = sorted_cov[n_samples/2]
        else:
            median_cov = (sorted_cov[(n_samples/2) - 1] + sorted_cov[(n_samples/2)]) / 2
        for i in xrange(len(threshold_counts)):
            cov = threshold_counts[i]
            threshold_counts[i] = round(float(cov) / len(sorted_cov), 2)

        o.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene, chrom, pos, mean_cov, median_cov, '\t'.join(threshold_counts))
