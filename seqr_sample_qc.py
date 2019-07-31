import logging
import argparse
import pickle

import hail as hl
import hail.expr.aggregators as agg

from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *
from gnomad_hail.utils import *
from gnomad_hail.utils.sample_qc import *
from resources_seqr_qc import *

logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("seqr_sample_qc")
logger.setLevel(logging.INFO)


def validate_mt(mt: hl.MatrixTable, build: int, data_type: str, threshold=0.3):
    """
    Validate the mt by checking against a list of common coding and non-coding variants given its
    genome version. This validates genome_version, variants, and the reported sample type.

    :param mt: mt to validate
    :param build: reference build
    :param data_type: WGS or WES
    :param threshold: Threshold percentage of variants matching validation tables
    :return: True or Exception
    """
    def sample_type_stats(mt, build, threshold):
        """
        Calculate stats for sample type by checking against a list of common coding and non-coding variants.
        If the match for each respective type is over the threshold, we return a match.

        :param mt: Matrix Table to check
        :param build: reference build
        :return: a dict of coding/non-coding to dict with 'matched_count', 'total_count' and 'match' boolean.
        :rtype Set
        """
        stats = {}
        types_to_ht_path = {
            'noncoding': val_noncoding_ht_path(build),
            'coding': val_coding_ht_path(build)
        }
        for variant_type, ht_path in types_to_ht_path.items():
            ht = hl.read_table(ht_path)
            stats[variant_type] = ht_stats = {
                'matched_count': mt.semi_join_rows(ht).count_rows(),
                'total_count': ht.count(),
            }

            ht_stats['match'] = (ht_stats['matched_count']/ht_stats['total_count']) >= threshold
        return stats


    sample_type_stats = sample_type_stats(mt, build, threshold)

    for name, stat in sample_type_stats.items():
        logger.info('Table contains %i out of %i common %s variants.' %
                    (stat['matched_count'], stat['total_count'], name))

    has_coding = sample_type_stats['coding']['match']
    has_noncoding = sample_type_stats['noncoding']['match']

    if not has_coding and not has_noncoding:
        raise VCFDataTypeError(
            f'Genome version validation error: dataset specified as GRCh{build} but doesn\'t contain '
            f'the expected number of common GRChf{build} variants'
        )
    elif has_noncoding and not has_coding:
        raise VCFDataTypeError(
            'Sample type validation error: Dataset contains noncoding variants but is missing common coding '
            f'variants for GRCh{build}. Please verify that the dataset contains coding variants.'
        )
    elif has_coding and not has_noncoding:
        if data_type != 'WES':
            raise VCFDataTypeError(
                f'Sample type validation error: dataset sample-type is specified as {data_type} but appears to be '
                'WES because it contains many common coding variants and lacks many non-coding variants'
            )
    elif has_noncoding and has_coding:
        if data_type != 'WGS':
            raise VCFDataTypeError(
                f'Sample type validation error: dataset sample-type is specified as {data_type} but appears to be '
                'WGS because it contains many common non-coding variants'
            )
    return True


def apply_filter_flags_expr(mt: hl.MatrixTable, data_type: str) -> hl.expr.SetExpression:
    """
    Annotates table with flags for elevated contamination and chimera as well as low coverage and call rate
    :param Table mt: input MatrixTable
    :param str data_type: 'WES' or 'WGS' for selecting coverage threshold
    :return: Set of sequencing metric flags
    :rtype: SetExpression
    """
    flags = {
        'callrate': mt.filtered_callrate < 0.85,
        'contamination': mt.PCT_CONTAMINATION > 5,  # TODO revisit current thresholds and rename once have to kristen's script output
        'chimera': mt.AL_PCT_CHIMERAS > 5
    }
    if data_type == 'WES':
        flags.update({
            'coverage': mt.HS_PCT_TARGET_BASES_20X < 85
        })
    else:
        flags.update({
            'coverage': mt.WGS_MEAN_COVERAGE < 30
        })

    return hl.set(hl.filter(lambda x: hl.is_defined(x),
                            [hl.or_missing(filter_expr, name) for name, filter_expr in flags.items()]))


def liftover_to_37(vcf_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Liftover input mt from GRCh38 to GRCh37 and return for gnomAD population assignment
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


def prep_meta(ht: hl.Table) -> hl.Table:
    """
    Preps gnomAD exome and genome metadata for population PC by selecting relevant fields only
    :param Table ht: Hail table from gnomAD's metadata column annotations
    :return: Hail Table ready for joining
    :rtype: Table
    """
    ht = ht.key_by('s')
    ht = ht.annotate(source='gnomAD')
    ht = ht.select('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'data_type', 'qc_pop', 'source')
    ht = ht.filter(hl.is_defined(ht.PC1)).rename({'qc_pop': 'known_pop'})
    return ht


def get_all_sample_metadata(mt: hl.MatrixTable, build: int, data_type: str, data_source: str, version: int) -> hl.Table:
    """
    Annotate MatrixTable with all current metadata: sample sequencing metrics, sample ID mapping,
    and callrate for bi-allelic, high-callrate common SNPs.
    :param MatrixTable mt: VCF converted to a MatrixTable
    :param int build: build for write, 37 or 38
    :param str data_type: WGS or WES for write path and flagging metrics
    :param str data_source: internal or external for write path
    :param int version: Int for write path
    :return: Table with seq metrics and mapping
    :rtype: Table
    """
    logger.info("Importing and annotating with sequencing metrics...")
    meta_ht = hl.import_table(seq_metrics_path(build, data_type, data_source, version), impute=True).key_by('SAMPLE')

    logger.info("Importing and annotating seqr ID names...")
    remap_ht = hl.import_table(remap_path(build, data_type, data_source, version), impute=True).key_by('vcf_id')
    meta_ht = meta_ht.annotate(**remap_ht[meta_ht.key])
    meta_ht = meta_ht.annotate(seqr_id=hl.cond(hl.is_missing(meta_ht.seqr_id), meta_ht.SAMPLE, meta_ht.seqr_id))

    logger.info("Annotating bi-allelic, high-callrate, common SNPs for sample QC...")
    mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
                        & (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
                        (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
    callrate_ht = mt.select_cols(filtered_callrate=hl.agg.fraction(hl.is_defined(mt.GT))).cols()
    meta_ht = meta_ht.annotate(**callrate_ht[meta_ht.key])
    return meta_ht


def run_platform_imputation(mt: hl.MatrixTable, plat_min_cluster_size: int, plat_min_sample_size: int) -> hl.Table:
    """
    Run PCA using sample callrate across gnomAD's evaluation interval and create Hail Table
    with platform PCs and assigned platform.
    :param MatrixTable mt: QC MatrixTable
    :param plat_min_cluster_size: min cluster size for HBDscan clustering
    :param plat_min_sample_size: min sample size for HBDscan clustering
    :return: Table with platform PCs and assigned platform
    :rtype: Table
    """
    intervals = hl.import_locus_intervals(evaluation_intervals_path)
    callrate_mt = compute_callrate_mt(mt, intervals)
    eigenvalues, scores_ht, _ = run_platform_pca(callrate_mt)
    plat_ht = assign_platform_from_pcs(scores_ht, hdbscan_min_cluster_size=plat_min_cluster_size,
                                       hdbscan_min_samples=plat_min_sample_size)
    d = {f'plat_PC{i+1}': scores_ht.scores[i] for i in list(range(0, 6))}
    scores_ht = scores_ht.annotate(**d).drop('scores')
    plat_ht = plat_ht.annotate(**scores_ht[plat_ht.key])
    return plat_ht


def run_population_pca(mt: hl.MatrixTable,
                       build: int,
                       data_type: str,
                       data_source: str,
                       version: int,
                       is_test: bool,
                       fit: Any = None,
                       ) -> hl.Table:
    """
    Preps data for population PCA by joining gnomAD metadata tables with new samples' scores from pc_project and selects
    columns for run_assign_population_pcs
    :param MatrixTable mt: QC MatrixTable
    :param int build: 37 or 38 for write path
    :param str data_type: WGS or WES for write path
    :param str data_source: internal or external for write path
    :param int version: Int for write path
    :param bool is_test: Boolean on whether testing pipeline for write path
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :return: Table annotated with assigned gnomAD population and PCs
    :rtype: Table

    """
    meta_exomes = hl.read_table(metadata_exomes_ht_path(version=CURRENT_EXOME_META))
    meta_genomes = hl.read_table(metadata_genomes_ht_path(version=CURRENT_GENOME_META))
    loadings = hl.read_table(ancestry_pca_loadings_ht_path())

    meta_exomes = prep_meta(meta_exomes)
    meta_genomes = prep_meta(meta_genomes)
    meta = meta_exomes.union(meta_genomes)

    mt = mt.select_entries('GT')
    scores = pc_project(mt, loadings)
    d = {f'PC{i+1}': scores.scores[i] for i in range(6)}
    scores = scores.annotate(**d).drop('scores')
    scores = scores.annotate(data_type=data_type, known_pop=hl.null(hl.tstr), source="RDG").key_by('s')
    meta = meta.select('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'data_type', 'known_pop', 'source')
    meta = meta.union(scores)

    logger.info('Assigning gnomAD populations...')
    data = meta.to_pandas()
    pc_cols = ['PC{}'.format(pc) for pc in range(1, 7)]
    pop_pca_ht, pop_clf = assign_population_pcs(data, pc_cols=pc_cols, output_col='qc_pop', fit=fit)
    pop_pca_ht = hl.Table.from_pandas(pop_pca_ht)

    if not fit:
        # Pickle RF
        with hl.hadoop_open(pop_RF_fit_path(build, data_type, data_source, version, is_test), 'wb') as out:
            pickle.dump(pop_clf, out)

    pop_pca_ht = pop_pca_ht.filter(pop_pca_ht.source == "RDG").drop('known_pop', 'source', 'data_type').key_by('s')
    pop_pca_ht = pop_pca_ht.annotate(**scores[pop_pca_ht.s])
    pop_pca_ht = pop_pca_ht.rename({f'PC{pc}': f'pop_PC{pc}' for pc in range(1,7)})
    logger.info("Found the following sample count after population assignment: {}".format(
        ", ".join(f'{pop}: {count}' for pop, count in pop_pca_ht.aggregate(hl.agg.counter(pop_pca_ht.qc_pop)).items())))
    return pop_pca_ht


def run_hail_sample_qc(mt: hl.MatrixTable, data_type: str) -> hl.MatrixTable:
    """
    Runs Hail's built in sample qc function on the MatrixTable. Splits the MatrixTable in order to calculate inbreeding
    coefficient and annotates the result back onto original MatrixTable. Applies flags by population and platform groups.
    :param MatrixTable mt: QC MatrixTable
    :param str data_type: WGS or WES for write path
    :return: MatrixTable annotated with hails sample qc metrics as well as pop and platform outliers
    :rtype: MatrixTable
    """
    mt = hl.split_multi_hts(mt)
    mt = hl.sample_qc(mt)
    mt = mt.annotate_cols(sample_qc=mt.sample_qc.annotate(f_inbreeding=hl.agg.inbreeding(mt.GT, mt.info.AF[0])))
    mt = mt.annotate_cols(idx=mt.qc_pop + "_" + hl.str(mt.qc_platform))

    qc_metrics = ['n_snp', 'r_ti_tv', 'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']
    if data_type == "WGS":
        qc_metrics = qc_metrics + ['call_rate']

    metric_ht = compute_stratified_metrics_filter(mt.cols(), qc_metrics, ['qc_pop', 'qc_platform'])
    checkpoint_pass = metric_ht.aggregate(hl.agg.count_where(hl.len(metric_ht.pop_platform_filters) == 0))
    logger.info('{0} samples found passing pop/platform-specific filtering'.format(checkpoint_pass))
    checkpoint_fail = metric_ht.aggregate(hl.agg.count_where(hl.len(metric_ht.pop_platform_filters) != 0))
    logger.info('{0} samples found failing pop/platform-specific filtering'.format(checkpoint_fail))
    metric_ht = metric_ht.annotate(sample_qc=mt.cols()[metric_ht.key].sample_qc)
    return metric_ht


def main(args):

    hl.init(tmp_dir='hdfs:///pc_relate.tmp/')
    logger.info("Beginning seqr sample QC pipeline...")

    data_type = args.data_type
    build = args.build
    data_source = args.data_source
    version = args.callset_version
    is_test = args.is_test
    overwrite = args.overwrite

    logger.info("Importing callset...")
    if not args.skip_write_mt:
        logger.info("Converting vcf to MatrixTable...")
        vcf = callset_vcf_path(build, data_type, data_source, version, is_test)
        hl.import_vcf(vcf, force_bgz=True, reference_genome=f'GRCh{build}',
                      min_partitions=4).write(mt_path(build, data_type, data_source, version, is_test), overwrite=True)
    mt = hl.read_matrix_table(mt_path(build, data_type, data_source, version, is_test))

    if not args.skip_validate_mt:
        logger.info("Validating data type...")
        validate_mt(mt, build, data_type)

    if is_test:
        logger.info('Creating test mt...')
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval('20', reference_genome=f'GRCh{build}')]).persist()

    if build == '38':
        mt = liftover_to_37(mt)  # Needed for gnomAD PC project
        mt = mt.checkpoint(liftover_mt_path(build, data_type, data_source, version, is_test), overwrite)
        logger.info('Liftover complete')

    logger.info("Annotating with sequencing metrics and filtered callrate...")
    meta_ht = get_all_sample_metadata(mt, build, data_type, data_source, version)
    mt = mt.annotate_cols(**meta_ht[mt.col_key])

    logger.info("Annotating with filter flags...")
    mt = mt.annotate_cols(filter_flags=apply_filter_flags_expr(mt, data_type))

    logger.info('Assign platform or product')
    if data_type == 'WES' and data_source == 'External':
        logger.info('Running platform imputation...')
        plat_ht = run_platform_imputation(mt, args.plat_min_cluster_size, args.plat_min_sample_size)
        mt = mt.annotate_cols(**plat_ht[mt.col_key])
    else:
        logger.info('Assigning platform from product in metadata...')
        mt = mt.annotate_cols(qc_platform=mt.PRODUCT)

        missing_metrics = mt.filter_cols(hl.is_defined(mt.PRODUCT), keep=False)  # Write out samples missing seq metrics
        missing_metrics.cols().select().export(missing_metrics_path(build, data_type, data_source, version))

        mt = mt.filter_cols(hl.is_defined(mt.PRODUCT))  # Remove control samples and blacklisted samples
        mt = mt.filter_rows(hl.agg.count_where(mt.GT.is_non_ref()) > 0)

    # Checkpoint because mt.to_pandas in run_population_pca nukes dataproc nodes otherwise
    mt = mt.checkpoint(mt_temp_path(build, data_type, data_source, version, is_test), overwrite)

    logger.info('Projecting gnomAD population PCs...')
    pop_ht = run_population_pca(mt, build, data_type, data_source, version, is_test, args.pop_pca_fit)
    mt = mt.annotate_cols(**pop_ht[mt.col_key])

    logger.info('Running Hail\'s sample qc...')
    hail_metric_ht = run_hail_sample_qc(mt, data_type)
    mt = mt.annotate_cols(**hail_metric_ht[mt.col_key])

    logger.info('Exporting sample QC tables...')
    ht = mt.cols()
    ht = ht.checkpoint(sample_qc_ht_path(build, data_type, data_source, version, is_test), overwrite)
    ht.flatten().export(ht_to_tsv_path(build, data_type, data_source, version, is_test))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--data-type', help="Sequencing data type (WES or WGS)", choices=["WES", "WGS"], required=True)
    parser.add_argument('-b', '--build', help='Reference build, 37 or 38', choices=["37", "38"], required=True)
    parser.add_argument('-v', '--callset-version', help='Version of callset vcf', type=int, required=True)
    parser.add_argument('--data-source', help="Data source (Internal or External)", choices=["Internal", "External"], required=True)
    parser.add_argument('--is-test', help='To run a test of the pipeline using test files and directories',
                        action='store_true')
    parser.add_argument('--plat-min-cluster-size', help='Minimum cluster size for platform pca labeling', default=40)
    parser.add_argument('--plat-min-sample-size', help='Minimum sample size for platform pca labeling', default=40)
    parser.add_argument('--skip-write-mt', help='Skip writing out qc mt', action='store_true')
    parser.add_argument('--skip-validate-mt', help='Skip validating the mt against common coding and noncoding variants', action='store_true')
    parser.add_argument('--pop-pca-fit', help='Random Forest classifier for population assignment')
    parser.add_argument('--project-list', help='List of seqr projects that are in the callset')
    parser.add_argument('--slack-channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help="Overwrite previous paths", action='store_true')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
