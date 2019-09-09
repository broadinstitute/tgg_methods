import argparse
import hail as hl
import logging
from gnomad_hail import *


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("summarize_coverage")
logger.setLevel(logging.INFO)


def main(args):

    num_partitions = args.n_partitions
    gene_list = args.gene_list
    hl.init(default_reference='GRCh38', log='/summarize_coverage.log')
    hl._set_flags(newaggs=None)

    bucket = args.bucket

    if args.read_coverage_files:
        with hl.hadoop_open(gene_list) as g:
            for line in g:
                gene = line.strip()
                fname = f'{bucket}/*.tsv'
                logger.info('Loading: {}'.format(fname))

                # guessing 50 partitions is good here
                ht = hl.import_table(fname, impute=True, no_header=True)
                ht = ht.transmute(gene=ht.f0, chrom=ht.f1, pos=ht.f2, mean=ht.f3, median=ht.f4)
                ht = ht.key_by('gene')
                ht = ht.repartition(num_partitions).write(f'{bucket}/coverage.ht', args.overwrite)

    if args.process_coverage:

        # GOAL: Gene	Size(kb)	Mean_Coverage	Median_Coverage	Callable_bases(%)	Uncallable_bases
        logger.info('Reading in coverage ht for processing')
        ht = hl.read_table(f'{bucket}/coverage.ht')

        # Monkol defined callable bases as bases with mean cov > 6
        ht_result = ht.group_by(ht.gene).aggregate(
                                            mean=hl.agg.mean(ht.mean),
                                            median=hl.median(hl.agg.collect(ht.median)),
                                            size_kb=(hl.agg.count() / 1000),
                                            pct_callable_bases=((hl.agg.count_where(ht.mean >= 6) / size_kb) * 100),
                                            pct_uncallable_bases=((hl.agg.count_where(ht.mean < 6) / size kb)* 100))

        # export grouped table
        ht_result.export(f'{bucket}/myoseq_gene_coverage_summary.tsv.bgz')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--read_coverage_files', help='Read raw coverage files', action='store_true')
    parser.add_argument('--process_coverage', help='Process raw, merged coverage ht', action='store_true')
    parser.add_argument('--n_partitions', help='Number of partitions for output', type=int, required=True)
    parser.add_argument('--gene_list', help='(cloud) Location of MyoSeq gene list', required=True)
    parser.add_argument('--bucket', help='Bucket (for input and output)', required=True)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
