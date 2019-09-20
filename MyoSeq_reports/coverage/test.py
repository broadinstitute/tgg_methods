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

    bucket = args.bucket

    myoseq_genes = {}
    with hl.hadoop_open(gene_list) as g:
        gene,size = line.strip().split('\t')
        myoseq_genes[gene] = size

    if args.import_coverage_files:

        logger.info('Importing coverage files')
        sample_names = []
        with hl.hadoop_open(sample_list) as s:
            for line in s:
                sample_names.append(line.strip())

        for g in myoseq_genes:
            size = myoseq_genes[g]
            fname = f'{bucket}/{g}.tsv.bgz'
            logger.info('Loading: {}'.format(fname))
            mt = hl.import_matrix_table(fname, no_header=True, min_partitions=num_partitions, row_fields={'f0': hl.tstr, 'f1': hl.tint})
            mt = mt.transmute_rows(locus=hl.locus(mt.f0, mt.f1)).annotate_rows(gene=g, size_kb=size)

            # key columns by sample
            # default col_key is called col_id, and default entry name is x
            mt = mt.key_cols_by(s=hl.literal(sample_data)[mt.col_id]).drop('col_id')
            mt.transmute_entries(coverage=mt.x).write(f'{bucket}/{g}.mt', args.overwrite)


    if args.merge_coverage:

        logger.info('Reading in coverage mts for merging')
        mt = hl.read_matrix_table(f'{bucket}/{myoseq_genes[0]}.mt')
        for i in range(1, len(myoseq_genes)):
            next_mt = hl.read_matrix_table(f'{bucket}/{myoseq_genes[i]}.mt')

            # need to merge mt rows -- each mt has same samples with coverage across different gene
            # NOTE: hail docs warn union_rows might cause shuffle - best to use workers
            mt = mt.union_rows(next_mt)

        mt = mt.naive_coalesce(num_partitions)
        mt.write(f'{bucket}/myoseq_gene_coverage.mt', args.overwrite)


    if args.process_coverage:
        # GOAL: need to get mean/median coverage across gene and each sample's mean/median cov 
        # group by rows for gene coverage
        # Gene	Size(kb)	Mean_Coverage	Median_Coverage	Callable_bases(%)	Uncallable_bases
        # Monkol defined callable bases as bases with mean cov > 6
        logger.info('Reading in coverage mt for processing')
        mt = hl.read_matrix_table(f'{bucket}/myoseq_gene_coverage.mt')
        logger.info('Getting coverage summaries per gene')
        gene_coverage = mt.group_rows_by(mt.gene)
                                                .aggregate(
                                                            mean=hl.agg.mean(mt.coverage),
                                                            median=hl.median(hl.agg.collect(mt.coverage)),
                                                            pct_callable_bases=((hl.agg.count_where(mt.mean >= 6) / mt.size_kb) * 100),
                                                            uncallable_base=hl.agg.count_where(mt.mean < 6)).result()
        gene_coverage.export(f'{bucket}/myoseq_gene_coverage_summary.tsv.bgz') 

        # Gene	Cohort_Mean	Sample_Mean	Cohort_Callable	Sample_Callable	Cohort_Uncallable	Sample_Uncallable
        logger.info('Getting per sample coverage summary')
        sample_coverage = mt.group_rows_by(mt.gene)
                                                .aggregate_entries(
                                                                mean=hl.agg.mean(mt.coverage),
                                                                median=hl.median(hl.agg.collect(mt.coverage)),
                                                                pct_callable_bases=((hl.agg.count_where(mt.mean >= 6) / mt.size_kb) * 100),
                                                                uncallable_base=hl.agg.count_where(mt.mean < 6)).result()
        sample_coverage.export(f'{bucket}/myoseq_gene_coverage_sumamry.tsv.bgz')
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--import_coverage_files', help='Import raw coverage files', action='store_true')
    parser.add_argument('--merge_coverage', help='Merge coverage mts', action='store_true')
    parser.add_argument('--process_coverage', help='Process raw, merged coverage mt', action='store_true')
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
