from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *
import hail as hl
import hdbscan
from resources.resources_seqr_qc import *
import logging
import argparse

logger = logging.getLogger("seqr_sample_qc")
logger.setLevel(logging.INFO)


def assign_platform_pcs(plat_pc_table: hl.Table, out_filepath: str, cluster_size: int, num_pcs: int = 9) -> hl.Table:
    """
    Function assumes that platform_pc_table contains columns named 'combined_sample',
     'gross_platform' (known labels), 'callratePC<n>'
    :param Table plat_pc_table: Table containing samples and callrate PCs
    :param str out_filepath: filepath for tsv containing samples, callrate PCs, and imputed platform labels
    :param int cluster_size: minimum cluster size for HDBSCAN
    :param int num_pcs: number of callrate PCs to use in platform imputation
    :return: Table containing samples, callrate PCs, and imputed platform labels
    :rtype: Table
    """
    # Read and format data for clustering
    data = plat_pc_table.to_pandas()
    cols = ['PC' + str(i + 1) for i in range(num_pcs)]
    callrate_data = data[cols].as_matrix()
    logger.info('Assigning platforms to {} exome samples in MT...'.format(len(callrate_data)))

    clusterer = hdbscan.HDBSCAN(min_cluster_size=cluster_size)
    cluster_labels = clusterer.fit_predict(callrate_data)
    n_clusters = len(set(cluster_labels)) - (
                -1 in cluster_labels)  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info('Found {} unique platforms during platform imputation...'.format(n_clusters))

    data['qc_platform'] = cluster_labels
    with hl.hadoop_open(out_filepath, 'w') as out:
        data.to_csv(out, sep="\t", index=False)
    new_data = hl.import_table(out_filepath, impute=True, types={'qc_platform': hl.tstr}).key_by('data_type', 's')
    return new_data


def main(args):
    hl.init(log='/platform_pca.log')

    data_type = "exomes"
    external = True if args.external else False
    version = args.callset_version
    cluster_size = args.min_cluster_size
    test = True if args.test else False

    if not args.skip_prepare_data_for_platform_pca:
        logger.info('Preparing data for platform PCA...')
        mt = hl.read_matrix_table(mt_path(data_type, external, version, test))
        mt = filter_to_autosomes(mt)
        intervals = hl.import_locus_intervals(evaluation_intervals_path)
        mt = mt.annotate_rows(interval=intervals[mt.locus].target)
        mt = mt.filter_rows(hl.is_defined(mt.interval) & (hl.len(mt.alleles) == 2))
        mt = mt.select_entries(GT=hl.or_missing(hl.is_defined(mt.GT), hl.struct()))
        callrate_mt = mt.group_rows_by(mt.interval).aggregate(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
        callrate_mt.write(exome_callrate_mt_path(external, version, test), args.overwrite)

    if not args.skip_run_platform_pca:
        logger.info('Running platform pca')
        callrate_mt = hl.read_matrix_table(exome_callrate_mt_path(external, version, test))
        callrate_mt = callrate_mt.annotate_entries(callrate=hl.int(callrate_mt.callrate > 0.25))
        # Center until Hail's PCA does it for you
        callrate_mt = callrate_mt.annotate_rows(mean_callrate=hl.agg.mean(callrate_mt.callrate))
        callrate_mt = callrate_mt.annotate_entries(callrate=callrate_mt.callrate - callrate_mt.mean_callrate)
        eigenvalues, scores_ht, _ = hl.pca(callrate_mt.callrate, compute_loadings=False)
        logger.info('Eigenvalues: {}'.format(eigenvalues))
        # external test [2036450.4046471054, 1706942.3950260752, 1152586.6865198347, 550812.9905120205, 461544.4847840831, 241506.66432095322, 163937.8513164449, 133651.98686930467, 106823.68126452898, 78806.9446078671]
        # internal test [397688.8211215944, 157481.4486448636, 41901.02898184424, 17979.67513516618, 16910.630213155262, 16278.087057597506, 15442.625207568679, 14567.677537599762, 14203.75719936242, 12410.252339536075]
        scores_ht.write(exome_callrate_scores_ht_path(external, version, test))

    logger.info('Annotating with platform PCs and known platform annotations...')
    scores_ht = hl.read_table(exome_callrate_scores_ht_path).annotate_globals(exome_source="external")
    d = {'PC%s' % (i + 1): scores_ht.scores[i] for i in range(10)}
    scores_ht = scores_ht.annotate(**d)
    scores_ht = scores_ht.drop('scores')
    platform_pcs = assign_platform_pcs(scores_ht, qc_temp_data(data_type, external, test, version) +
                                       'assigned_platform_pcs.txt.bgz', cluster_size)
    platform_pcs.write(platform_labels_path(data_type, external, test, version))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--build', help='Reference build, 37 or 38', choices=["37", "38"], required=True) #TODO:Don't need this just yet but will
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--internal', help='Internal sample data in callset', action='store_true')
    parser.add_argument('--external', help='External sample data in callset', action='store_true')
    parser.add_argument('--min_cluster_size', help='Minimum cluster size for pca labeling', defaul=40)
    parser.add_argument('--test', help='To run a test of the pipeline using test files and directories',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--skip_prepare_data_for_platform_pca',
                        help='Skip prepping data for platform imputation (assuming already done)', action='store_true')
    parser.add_argument('--skip_run_platform_pca', help='Skip platform PCA (assuming already done)',
                        action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
