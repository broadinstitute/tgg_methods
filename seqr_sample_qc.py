# TODO main is too long - break out pcas into separate functions
# TODO add in kristen's metric script and confirm output format


import hail as hl
import hail.expr.aggregators as agg
from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *
from gnomad_hail.utils import *
import hdbscan
import logging
import argparse
import pandas as pd
import pickle

from resources_seqr_qc import *
#from call_sex import *
#from retrieve_ped_from_seqr import *

import getpass


logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("seqr_sample_qc")
logger.setLevel(logging.INFO)


def apply_filter_flags_expr(mt: hl.MatrixTable, data_type: str) -> hl.expr.SetExpression:
    """
    Annotates table with flags for elevated contamination and chimera as well as low coverage and call rate
    :param Table mt: input MatrixTable
    :param str data_type: 'WES' or 'WGS' for selecting coverage threshold
    :return: Set of sequencing metric flags
    :rtype: SetExpression
    """
    flags = {
        'callrate_less_than_0.85': mt.filtered_callrate < 0.85,
        'contamination_greater_than_0.05': mt.freemix > 0.05,  # TODO revisit current thresholds and rename to kristen's output
        'chimera_greater_than_0.05': mt.pct_chimeras > 0.05  # TODO(cont) almost everything under 0.05 for contam and chimera
    }
    if data_type == 'WES':
        flags.update({
            'coverage_less_than_85_at_20x': mt.pct_target_bases_20x < 0.85
        })
    else:
        flags.update({
            'coverage_less_than_30x': mt.mean_target_coverage < 30
        })

    return hl.set(hl.filter(lambda x: hl.is_defined(x),
                            [hl.or_missing(filter_expr, name) for name, filter_expr in flags.items()]))


def liftover_to_37(vcf_mt: hl.MatrixTable)->hl.MatrixTable:
    """
    Liftover input mt to GRCh37 and return for gnomAD population assignment
    :param vcf_mt: MatrixTable on GRCh38
    :return: MatrixTable with GRCh37 locus
    :rtype MatrixTable:
    """
    logger.info('Lifting a 38 genome over to 37')
    chain_file = "gs://hail-common/references/grch38_to_grch37.over.chain.gz"
    grch38 = hl.get_reference('GRCh38')
    grch38.add_liftover(chain_file, 'GRCh37')
    vcf_mt = vcf_mt.annotate_rows(liftover_locus=hl.liftover(vcf_mt.locus, 'GRCh37'))
    vcf_mt = vcf_mt.filter_rows(hl.is_defined(vcf_mt.liftover_locus))
    vcf_mt = vcf_mt.rename({'locus': 'locus_grch38', 'liftover_locus': 'locus'})
    vcf_mt = vcf_mt.key_rows_by('locus', 'alleles')
    return vcf_mt


def assign_platform_pcs(plat_pc_table: hl.Table, out_filepath: str, cluster_size: int, num_pcs: int = 9) -> hl.Table:
    """
    Function assumes that platform_pc_table contains columns named 'combined_sample',
     'gross_platform' (known labels), 'callratePC<n>'
    :param Table plat_pc_table: Table containing samples and callrate PCs
    :param str out_filepath: filepath for tsv containing samples, callrate PCs, and imputed platform labels
    :param int cluster_size: minimum cluster size for HDBSCAN - 40 performed best after testing 25,40,50,100
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
    n_clusters = len(set(cluster_labels)) - (-1 in cluster_labels)  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info('Found {} unique platforms during platform imputation...'.format(n_clusters))

    data['qc_platform'] = cluster_labels
    with hl.hadoop_open(out_filepath, 'w') as out:
        data.to_csv(out, sep="\t", index=False)
    new_data = hl.import_table(out_filepath, impute=True, types={'qc_platform': hl.tstr}).key_by('s')  # was seqr_id
    return new_data


def prep_meta(ht: hl.Table)->hl.Table:
    """
    Preps gnomAD exome and genome metadata for population PC by selecting relevant fields only
    :param Table ht: Hail table from gnomAD's metadata column annotations
    :return: Hail Table ready for joining
    :rtype: Table
    """
    ht = ht.key_by('s')
    ht = ht.annotate(source='gnomAD')
    ht = ht.select('qc_pop', 'pop_platform_filters', 'PC1',
                   'PC2', 'PC3', 'PC4', 'PC5', 'PC6',
                   'PC7', 'PC8', 'PC9', 'PC10', 'data_type', 'source')
    ht = ht.rename({'qc_pop': 'gnomad_pop'})
    ht = ht.filter(hl.is_defined(ht.PC1))
    return ht


def assign_population_pcs(
        pop_pc_pd: pd.DataFrame,
        pc_cols: List[str],
        known_col: str = 'gnomad_pop',
        fit: RandomForestClassifier = None,
        seed: int = 42,
        prop_train: float = 0.8,
        n_estimators: int = 100,
        min_prob: float = 0.9,
        output_col: str = 'pop',
        missing_label: str = 'oth') -> Tuple[pd.DataFrame, RandomForestClassifier]:
    """
    This function uses a random forest model to assign population labels based on the results of PCA.
    Default values for model and assignment parameters are those used in gnomAD.
    Note that if PCs come from Hail, they will be stored in a single column named `scores` by default.

    :param pd.DateFrame pop_pc_pd: Pandas dataframe containing population PCs as well as a column with population labels
    :param List[str] pc_cols: Columns storing the PCs to use
    :param str known_col: Column storing the known population labels
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :param int seed: Random seed
    :param float prop_train: Proportion of known data used for training
    :param int n_estimators: Number of trees to use in the RF model
    :param float min_prob: Minimum probability of belonging to a given population for the population to be set (otherwise set to `None`)
    :param str output_col: Output column storing the assigned population
    :param str missing_label: Label for samples for which the assignment probability is smaller than `min_prob`
    :return: Dataframe containing sample IDs and imputed population labels, trained random forest model
    :rtype: Tuple[pd.DataFrame, RandomForestClassifier]
    """

    train_data = pop_pc_pd.loc[~pop_pc_pd[known_col].isnull()]

    n = len(train_data)

    # Split training data into subsets for fitting and evaluating
    if not fit:
        random.seed(seed)
        train_subsample_ridx = random.sample(list(range(0, n)), int(n * prop_train))
        train_fit = train_data.iloc[train_subsample_ridx]
        fit_samples = [x for x in train_fit['s']]
        evaluate_fit = train_data.loc[~train_data['s'].isin(fit_samples)]

        # Train RF
        training_set_known_labels = train_fit[known_col].values
        training_set_pcs = train_fit[pc_cols].values
        evaluation_set_pcs = evaluate_fit[pc_cols].values

        pop_clf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)
        pop_clf.fit(training_set_pcs, training_set_known_labels)
        print('Random forest feature importances are as follows: {}'.format(pop_clf.feature_importances_))

        # Evaluate RF
        predictions = pop_clf.predict(evaluation_set_pcs)
        error_rate = 1 - sum(evaluate_fit[known_col] == predictions) / float(len(predictions))
        print('Estimated error rate for RF model is {}'.format(error_rate))
    else:
        pop_clf = fit

    # Classify data
    pop_pc_pd[output_col] = pop_clf.predict(pop_pc_pd[pc_cols].values)
    probs = pop_clf.predict_proba(pop_pc_pd[pc_cols].values)
    probs = pd.DataFrame(probs, columns=[f'prob_{p}' for p in pop_clf.classes_])
    pop_pc_pd = pd.concat([pop_pc_pd, probs], axis=1)
    probs['max'] = probs.max(axis=1)
    pop_pc_pd.loc[probs['max'] < min_prob, output_col] = missing_label
    return pop_pc_pd, pop_clf


def run_assign_population_pcs(pop_pc_table: hl.Table, outfile: str, picklefile: str, pcs: List[int],
                              fit: RandomForestClassifier = None,
                              seed: int = 42) -> Tuple[hl.Table, RandomForestClassifier]:
    """
    :param Table pop_pc_table: Table containing population PCs ('PC<n>') as well as a column 'known_pop' with pop labels
    :param str outfile: filepath to tsv with input samples and imputed population labels
    :param str picklefile: filepath to which the pickled random forest model is written
    :param List of int pcs: 1-based list of PCs to train the model on
    :param RandomForestClassifier fit: fit from a previously trained RF model
    :param int seed: Random seed
    :return: Table containing sample IDs and imputed population labels, trained random forest model
    :rtype: Table, RandomForestClassifier
    """
    data = pop_pc_table.to_pandas()
    pc_cols = ['PC{}'.format(pc) for pc in pcs]
    new_data, pop_clf = assign_population_pcs(data, pc_cols, fit=fit, seed=seed)

    if not fit:
        # Pickle RF
        with hl.hadoop_open(picklefile, 'wb') as out:
            pickle.dump(pop_clf, out)

    with hl.hadoop_open(outfile, 'w') as out:
        new_data.to_csv(out, sep="\t", na_rep="NA", index=False)
    return hl.import_table(outfile, impute=True).key_by('data_type', 's'), pop_clf


def calculate_metrics(mt: hl.MatrixTable, qc_metrics: List[str], data_type: str, is_external: bool, version: int,
                      test: bool) -> hl.Table:
    """
    Compute median, MAD, and upper and lower thresholds for each metric used in population- and
     platform-specific outlier filtering
    :param MatrixTable mt: MT containing relevant sample QC metric annotations
    :param List[str] qc_metrics: list of metrics for which to compute the critical values for filtering outliers
    :param str data_type: 'WES' or 'WGS'
    :param int version: Version of callset
    :param bool is_external: Whether or not the callset contains external data
    :param bool test: Is this a test run?
    :return: Table grouped by pop and platform, with upper and lower threshold values computed for each sample QC metric
    :rtype: Table
    """
    cols = ['qc_pop', 'qc_platform']
    key = 'seqr_id'
    colnames = cols + ['sample_qc.{}'.format(metric) for metric in qc_metrics]
    new_colnames = cols + ['{}'.format(metric) for metric in qc_metrics]
    ht = mt.cols().flatten().key_by(key)
    ht = ht.select(*colnames).rename(dict(zip(colnames, new_colnames)))

    key_expr = ['{0}_median'.format(metric) for metric in qc_metrics] + ['{0}_mad'.format(metric)
                                                                         for metric in qc_metrics]
    value_expr = [hl.median(hl.agg.collect(ht['{}'.format(metric)])) for metric in qc_metrics] \
                 + [1.4826 * hl.median(hl.abs(hl.agg.collect(ht['{0}'.format(metric)])
                                              - hl.median(hl.agg.collect(ht['{}'.format(metric)]))
                                              )) for metric in qc_metrics]
    ht = ht.group_by(ht.qc_pop, ht.qc_platform).aggregate(**dict(zip(key_expr, value_expr)))

    upper_key_expr = ['{0}_upper'.format(metric) for metric in qc_metrics]
    upper_value_expr = [ht['{0}_median'.format(metric)] + 4 * ht['{0}_mad'.format(metric)] if metric != 'callrate'
                        else 1 for metric in qc_metrics]

    lower_key_expr = ['{0}_lower'.format(metric) for metric in qc_metrics]
    lower_value_expr = [ht['{0}_median'.format(metric)] - 4 * ht['{0}_mad'.format(metric)] if metric != 'callrate'
                        else 0.99 for metric in qc_metrics]

    annotation = dict(zip(upper_key_expr, upper_value_expr))
    annotation.update(dict(zip(lower_key_expr, lower_value_expr)))
    ht = ht.annotate(**annotation)
    ht = ht.annotate(idx=ht.qc_pop + "_" + hl.str(ht.qc_platform))
    ht = ht.filter(hl.is_defined(ht.idx))
    ht.export(metric_MAD_data_path(data_type, is_external, version, test))
    return ht


def apply_pop_platform_filters(ht: hl.Table, mt: hl.MatrixTable, qc_metrics: List[str]) -> hl.MatrixTable:
    """
    Flag samples that fall outside pop- and platform-specific QC metric thresholds
    :param Table ht: Table with upper- and lower-threshold values for each pop- and platform-specific cohort of samples
    :param MatrixTable mt: MT containing samples to which pop- and platform-specific outlier filtering is applied
    :param List['str'] qc_metrics: list of metrics on which to flag outlier samples
    :return: MatrixTable with the pop_platform_filter annotation summarizing metrics on which a sample is an outlier
    :rtype: MatrixTable
    """
    thresh_colnames = ['idx'] + ['{}_upper'.format(metric) for metric in qc_metrics] + [
        '{}_lower'.format(metric) for metric in qc_metrics]
    ht = ht.select(*thresh_colnames).key_by('idx')
    ht.write('/tmp.h', True)  # temporary write/read for hail assert bug... Do we need this?
    ht = hl.read_table('/tmp.h')
    mt = mt.annotate_cols(_thresholds=ht[mt.idx])
    fail_exprs = {
        'fail_{}'.format(metric):
            (mt.sample_qc['{}'.format(metric)] >= mt._thresholds['{}_upper'.format(metric)]) |
            (mt.sample_qc['{}'.format(metric)] <= mt._thresholds['{}_lower'.format(metric)]) for metric in qc_metrics}
    mt = mt.annotate_cols(**fail_exprs).drop('_thresholds')
    pop_platform_filters = make_pop_filters_expr(mt, qc_metrics)
    return mt.annotate_cols(pop_platform_filters=pop_platform_filters).drop('idx')


def make_pop_filters_expr(mt: hl.MatrixTable, qc_metrics: List[str]) -> hl.expr.SetExpression:
    """
    :param MatrixTable mt: MT to which the pop-/platform-specific filtering should apply
    :param List[str] qc_metrics: list of sample QC metrics for which to collect pop-/platform-specific filtering status
    :return: SetExpression value used to create the pop_platform_filters sample annotation
    :rtype: SetExpression
    """
    return hl.set(hl.filter(lambda x: hl.is_defined(x),
                            [hl.or_missing(mt['fail_{}'.format(metric)],
                                           '{}'.format(metric)) for metric in qc_metrics]))


def flag_high_callrate(mt: hl.MatrixTable, data_type: str, is_external: bool, version: int, test: bool, overwrite: bool):
    """
    This function filters the a copy of the MatrixTable to bi-allelic, high-callrate SNPs and then
    annotates the original MatrixTable with the callrate
    :param MatrixTable mt: Callset VCF converted MT
    :param str data_type: WGS or WES for write path
    :param bool is_external: Boolean for write path
    :param int version: Integer representing the version of the callset for write path
    :param bool test: Boolean on whether testing pipeline for write path
    :param bool overwrite: Boolean to overwrite the previous written data if present
    """
    logger.info("Annotating bi-allelic, high-callrate, common SNPs for sample QC...")
    bi_mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
                           & (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
                           (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
    bi_mt = bi_mt.annotate_cols(filtered_callrate=hl.agg.fraction(hl.is_defined(bi_mt.GT)))
    mt = mt.annotate_cols(filtered_callrate=bi_mt.cols()[mt.col_key].filtered_callrate)
    mt.naive_coalesce(5000).write(qc_mt_path(data_type, is_external, version, test), overwrite=overwrite)


def annotate_all_sample_metadata(qc_mt: hl.MatrixTable, data_type: str, is_external: bool, version: int, test: bool) -> hl.MatrixTable:
    """
    Annotate MatrixTable with all current metadata: sample sequencing metrics, sample ID mapping, project ID mapping.
    Flag all samples with sequencing metric outliers
    :param MatrixTable qc_mt: VCF converted to a MatrixTable
    :param str data_type: WGS or WES for write path and flagging metrics
    :param bool is_external: Boolean for write path
    :param int version: Int for write path
    :param bool test: Boolean on whether testing pipeline for write path
    :return: MatrixTable with metric, mapping, and flag annotations
    :rtype: MatrixTable
    """
    logger.info("Importing and annotating with metadata...")
    meta_ht = hl.import_table(metadata_path(data_type, is_external, version, test), impute=True).key_by('sample')
    qc_mt = qc_mt.annotate_cols(**meta_ht[qc_mt.s])

    logger.info("Importing and annotating seqr ID names...")
    remap_ht = hl.import_table(remap_path(data_type, is_external, version, test), impute=True).key_by('vcf_id')
    qc_mt = qc_mt.annotate_cols(**remap_ht[qc_mt.s])
    remap_expr = hl.cond(hl.is_missing(qc_mt.seqr_id), qc_mt.s, qc_mt.seqr_id)
    qc_mt = qc_mt.annotate_cols(seqr_id=remap_expr).key_cols_by('seqr_id')

    logger.info("Mapping C project to seqr project guid")  # TODO Create c_project to project_guid Mapping
    #project_ht = hl.import_table(project_map_path(data_type, is_external, version, test), impute=True).key_by(
    #    'c_project')
    #qc_mt = qc_mt.annotate_cols(**project_ht[qc_mt.project])

    logger.info("Annotating with filter flags...")
    qc_mt = qc_mt.annotate_cols(filter_flags=apply_filter_flags_expr(qc_mt, data_type))
    return qc_mt


def run_platform_pca(qc_mt: hl.MatrixTable, data_type: str, is_external: bool, version: int, cluster_size: int, overwrite: bool, test: bool) -> hl.MatrixTable:
    """
    Prepares a callrate MatrixTable and then runs assign_platform_pcs. Annotates original QC MatrixTable with assigned platform and PCs.
    :param MatrixTable qc_mt: QC MatrixTable
    :param str data_type: WGS or WES for write path
    :param bool is_external: Boolean for write path
    :param int version: Int for write path
    :param int cluster_size: Minimum cluster size for HDBSCAN clustering of platforms
    :param bool overwrite: Boolean for write path
    :param bool test: Boolean on whether testing pipeline for write path
    :return: MatrixTable annotated with assigned platform and PCs
    :rtype: MatrixTable
    """
    mt = hl.read_matrix_table(mt_path(data_type, is_external, version, test))
    mt = filter_to_autosomes(mt)
    intervals = hl.import_locus_intervals(evaluation_intervals_path)
    mt = mt.annotate_rows(interval=intervals[mt.locus].target)
    mt = mt.filter_rows(hl.is_defined(mt.interval) & (hl.len(mt.alleles) == 2))
    mt = mt.select_entries(GT=hl.or_missing(hl.is_defined(mt.GT), hl.struct()))
    callrate_mt = mt.group_rows_by(mt.interval).aggregate(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
    callrate_mt.write(exome_callrate_mt_path(is_external, version, test), overwrite=overwrite)
    callrate_mt = hl.read_matrix_table(exome_callrate_mt_path(is_external, version, test))
    callrate_mt = callrate_mt.annotate_entries(callrate=hl.int(callrate_mt.callrate > 0.25))

    # Center until Hail's PCA does it for you
    callrate_mt = callrate_mt.annotate_rows(mean_callrate=hl.agg.mean(callrate_mt.callrate))
    callrate_mt = callrate_mt.annotate_entries(callrate=callrate_mt.callrate - callrate_mt.mean_callrate)
    eigenvalues, scores_ht, _ = hl.pca(callrate_mt.callrate, compute_loadings=False)
    logger.info('Eigenvalues: {}'.format(eigenvalues))
    # external test [2036450.4046471054, 1706942.3950260752, 1152586.6865198347, 550812.9905120205, 461544.4847840831, 241506.66432095322, 163937.8513164449, 133651.98686930467, 106823.68126452898, 78806.9446078671]
    # internal test [397688.8211215944, 157481.4486448636, 41901.02898184424, 17979.67513516618, 16910.630213155262, 16278.087057597506, 15442.625207568679, 14567.677537599762, 14203.75719936242, 12410.252339536075]
    scores_ht.write(exome_platform_callrate_scores_ht_path(is_external, version, test), overwrite=overwrite)

    logger.info('Annotating with platform PCs and known platform annotations...')
    scores_ht = hl.read_table(exome_platform_callrate_scores_ht_path(is_external, version, test)).annotate_globals(exome_source="external")
    d = {'PC%s' % (i + 1): scores_ht.scores[i] for i in range(10)}
    scores_ht = scores_ht.annotate(**d)
    scores_ht = scores_ht.drop('scores')
    platform_pcs = assign_platform_pcs(scores_ht, qc_temp_data(data_type, is_external, test, version)
                                       + 'assigned_platform_pcs.txt.bgz', cluster_size)
    platform_pcs = platform_pcs.rename({'PC1': 'plat_PC1', 'PC2': 'plat_PC2', 'PC3': 'plat_PC3', 'PC4': 'plat_PC4',
                                        'PC5': 'plat_PC5', 'PC6': 'plat_PC6', 'PC7': 'plat_PC7', 'PC8': 'plat_PC8',
                                        'PC9': 'plat_PC9', 'PC10': 'plat_PC10'})
    platform_pcs.write(platform_labels_path(data_type, is_external, test, version))
    qc_mt = qc_mt.annotate_cols(**platform_pcs[qc_mt.col_key])
    qc_mt.write(qc_mt_path(data_type, is_external, version, test), overwrite=overwrite)
    qc_mt = hl.read_matrix_table(qc_mt_path(data_type, is_external, version, test))
    return qc_mt


def run_population_pca(qc_mt: hl.MatrixTable, data_type: str, is_external: bool, version: int, overwrite: bool, test: bool)-> hl.MatrixTable:
    """
    Preps data for population PCA by joining gnomAD metadata tables with new samples' scores from pc_project and selects
    columns for run_assign_population_pcs
    :param MatrixTable qc_mt: QC MatrixTable
    :param str data_type: WGS or WES for write path
    :param bool is_external: Boolean for write path
    :param int version: Int for write path
    :param bool overwrite: Boolean for write path
    :param bool test: Boolean on whether testing pipeline for write path
    :return: MatrixTable annotated with assigned gnomAD population and PCs
    :rtype: MatrixTable

    """
    qc_mt = qc_mt.select_entries('GT')
    meta_exomes = hl.read_table(metadata_exomes_ht_path(version=CURRENT_EXOME_META))
    meta_genomes = hl.read_table(metadata_genomes_ht_path(version=CURRENT_GENOME_META))
    loadings = hl.read_table(ancestry_pca_loadings_ht_path())
    meta_exomes = prep_meta(meta_exomes)
    meta_genomes = prep_meta(meta_genomes)
    meta = meta_exomes.union(meta_genomes)

    scores = pc_project(qc_mt, loadings)
    d = {'PC%s' % (i + 1): scores.scores[i] for i in range(10)}
    scores = scores.annotate(**d)
    scores = scores.drop('scores')
    scores = scores.annotate(data_type="exomes", gnomad_pop=hl.null(hl.tstr), source="RDG").rename({'seqr_id': 's'})
    meta = meta.select('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'data_type', 'gnomad_pop',
                       'source')
    meta = meta.union(scores)
    pcs = list(range(1, 7))

    logger.info('Assigning gnomAD populations...')
    joint_pca_ht_pop, joint_pca_fit_pop = run_assign_population_pcs(meta,
                                                                    population_assignments_ht_path(data_type,
                                                                                                   is_external,
                                                                                                   version, test),
                                                                    population_RF_fit_path(data_type, is_external,
                                                                                           version, test), pcs=pcs)
    joint_pca_ht_pop.export(population_assignments_tsv_path(data_type, is_external, version, test))
    joint_pca_ht_pop = joint_pca_ht_pop.rename({'PC1': 'pop_PC1', 'PC2': 'pop_PC2', 'PC3': 'pop_PC3', 'PC4': 'pop_PC4',
                                                'PC5': 'pop_PC5', 'PC6': 'pop_PC6', 'PC7': 'pop_PC7', 'PC8': 'pop_PC8',
                                                'PC9': 'pop_PC9', 'PC10': 'pop_PC10', 'gnomad_pop': 'qc_pop'}).key_by('s')
    qc_mt = qc_mt.annotate_cols(**joint_pca_ht_pop[qc_mt.seqr_id])
    qc_mt.write(qc_mt_temp_path(data_type, is_external, version, test), overwrite)
    qc_mt = hl.read_matrix_table(qc_mt_temp_path(data_type, is_external, version, test))
    return qc_mt


def run_hail_sample_qc(qc_mt: hl.MatrixTable, data_type: str, is_external: bool, version: int, test: bool)-> hl.MatrixTable:
    """
    Runs hails built in sample qc function on the MatrixTable. Splits the MatrixTable in order to calculate inbreeding
    coefficient and annotates the result back onto original MatrixTable. Applies flags by population and platform groups.
    :param MatrixTable qc_mt: QC MatrixTable
    :param str data_type: WGS or WES for write path
    :param bool is_external: Boolean for write path
    :param int version: Int for write path
    :param bool test: Boolean on whether testing pipeline for write path
    :return: MatrixTable annotated with hails sample qc metrics as well as pop and platform outliers
    :rtype: MatrixTable
    """
    split_mt = hl.split_multi_hts(qc_mt)
    split_mt = hl.sample_qc(split_mt)
    qc_mt = hl.sample_qc(qc_mt)
    split_mt = split_mt.annotate_cols(
        sample_qc=split_mt.sample_qc.annotate(f_inbreeding=hl.agg.inbreeding(split_mt.GT, split_mt.info.AF[0])))
    qc_mt = qc_mt.annotate_cols(f_inbreeding=split_mt.cols()[qc_mt.col_key].sample_qc.f_inbreeding)
    qc_mt = qc_mt.annotate_cols(idx=qc_mt.qc_pop + "_" + hl.str(qc_mt.qc_platform))
    qc_metrics = ['n_snp', 'r_ti_tv', 'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']
    if data_type == "WGS":
        qc_metrics = qc_metrics + ['call_rate']
    metrics_ht = calculate_metrics(qc_mt, qc_metrics, data_type, is_external, version, test)
    logger.info('Flagging samples failing pop/platform-specific sample qc thresholds...')
    qc_mt = apply_pop_platform_filters(metrics_ht, qc_mt, qc_metrics)
    checkpoint_pass = qc_mt.aggregate_cols(hl.agg.count_where(hl.len(qc_mt.pop_platform_filters) == 0))
    logger.info('{0} samples found passing pop/platform-specific filtering'.format(checkpoint_pass))
    checkpoint_fail = qc_mt.aggregate_cols(hl.agg.count_where(hl.len(qc_mt.pop_platform_filters) != 0))
    logger.info('{0} samples found failing pop/platform-specific filtering'.format(checkpoint_fail))
    return qc_mt


def main(args):

    hl.init(tmp_dir='hdfs:///pc_relate.tmp/')
    logger.info("Beginning seqr sample QC pipeline...")

    data_type = args.sample_type
    build = args.build
    is_external = True if args.data_source == "external" else False
    version = args.callset_version
    test = True if args.test else False
    male_f_threshold = args.male_f_threshold
    female_f_threshold = args.female_f_threshold
    cluster_size = args.plat_min_cluster_size
    overwrite = args.overwrite

    logger.info("Converting vcf to MatrixTable and flagging for high callrate...")
    if not args.skip_write_qc_mt:
        vcf = callset_vcf_path(data_type, is_external, version, test)
        hl.import_vcf(vcf, force_bgz=True, reference_genome=f'GRCh{build}',
                      min_partitions=4).write(mt_path(data_type, is_external, version, test), overwrite=True)
        mt = hl.read_matrix_table(mt_path(data_type, is_external, version, test))
        if build == '38':
            mt = liftover_to_37(mt)
            mt.write(qc_mt_temp_path(data_type, is_external, version, test), overwrite)
            mt = hl.read_matrix_table(qc_mt_temp_path(data_type, is_external, version, test), overwrite)
        flag_high_callrate(mt, data_type, is_external, version, test, overwrite)

    qc_mt = hl.read_matrix_table(qc_mt_path(data_type, is_external, version, test))
    qc_mt = annotate_all_sample_metadata(qc_mt, data_type, is_external, version, test)

    # Need to uncomment and add in Kristen's modules once they are ready
    #logger.info("Computing and annotating sex...")
    #call_sex(mt_path(data_type, is_external, version))
    #sex_ht = hl.import_table(sex_check_path(data_type, is_external, version)).key_by('s')
    #qc_mt = qc_mt.annotate_cols(**sex_ht[qc_mt.s])

    #logger.info("Annotating with pedigree concordance...")
    #ped_ht = hl.import_table(ped_path(data_type, is_external, version, test)).key_by('Individual_ID') # TODO: Make launch script for pipeline where ped retreive runs locally and first and result passed

    #qc_mt = qc_mt.annotate_cols(**ped_ht[qc_mt.seqr_id])
    #qc_mt = qc_mt.rename({'Sex': 'ped_sex'})
    #qc_mt = qc_mt.drop('Notes')
    #concordance_expr = hl.cond(qc_mt.ped_sex == qc_mt.sex, True, False)
    #qc_mt = qc_mt.annotate_cols(sex_concordance=concordance_expr)

    qc_mt.write(qc_mt_temp_path(data_type, is_external, version, test), overwrite)
    qc_mt = hl.read_matrix_table(qc_mt_temp_path(data_type, is_external, version, test))

    logger.info('Assign platform or product')
    if is_external and not args.skip_platform_pca:
        logger.info('Running platform pca...')
        qc_mt = run_platform_pca(qc_mt, data_type, is_external, version, cluster_size, overwrite, test)
    else:
        logger.info('Assigning platform from product in metadata...')
        qc_mt = qc_mt.annotate_cols(qc_platform=qc_mt.PRODUCT)
        qc_mt.write(qc_mt_path(data_type, is_external, version, test), overwrite)
        qc_mt = hl.read_matrix_table(qc_mt_path(data_type, is_external, version, test))

    if not args.skip_pop_pca:
        logger.info('Projecting gnomAD population PCs...')
        qc_mt = run_population_pca(qc_mt, data_type, is_external, version, overwrite, test)

    if not args.skip_calculate_sample_metrics:
        logger.info('Running hail\'s sample qc...')
        qc_mt = run_hail_sample_qc(qc_mt, data_type, is_external, version, test)

    qc_mt.write(qc_mt_path(data_type, is_external, version, test), overwrite)
    qc_mt.cols().write(qc_ht_path(data_type, is_external, version, test))  # TODO Change to .tsv?


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("--sample-type", help="sample type (WES or WGS)", choices=["WES", "WGS"], required=True)
    parser.add_argument('-b', '--build', help='Reference build, 37 or 38', choices=["37", "38"], required=True)
    parser.add_argument('-v', '--callset-version', help='Version of callset vcf', type=int, required=True)
    parser.add_argument('--data-source', help="Data source type (internal or external)", choices=["internal", "external"], required=True)
    parser.add_argument('--test', help='To run a test of the pipeline using test files and directories',
                        action='store_true')
    parser.add_argument('--male-f-threshold', help='Male f threshold for computing sex', type=float, default=0.8)
    parser.add_argument('--female-f-threshold', help='Female f threshold for computing sex', type=float, default=0.5)
    parser.add_argument('--plat-min-cluster-size', help='Minimum cluster size for platform pca labeling', default=40)
    parser.add_argument('--skip-write-qc-mt', help='Skip writing out qc mt', action='store_true')
    parser.add_argument('--skip-platform-pca', help='Skip platform PCA (assuming already done)', action='store_true')
    parser.add_argument('--skip-pop-pca', help='Skip calculating population PCs on unrelated samples',
                        action='store_true')
    parser.add_argument('--skip-calculate-sample-metrics', help='Skip calculating sample metrics (sample_qc)',
                        action='store_true')
    parser.add_argument('--ped', help='Pedigree bucket path')
    parser.add_argument('--project-list', help='List of seqr projects that are in the callset')
    parser.add_argument('--output', help='Override default output directory with given directory.')
    parser.add_argument('--slack-channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help="Overwrite previous paths", action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
