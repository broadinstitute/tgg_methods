from gnomad_hail import *
from resources.resources_seqr_qc import *
import hail as hl
import logging
import sys
import argparse


logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("seqr_sample_qc")
logger.setLevel(logging.INFO)


def apply_filter_flags_expr(ht: hl.Table, data_type: str) -> hl.expr.SetExpression:
    """
    Annotates table with flags for elevated contamination and chimera as well as low coverage and call rate
    :param: Table ht: input MT
    :param: str data_type: 'WES' or 'WGS'
    :return: output MT
    :rtype: SetExpression
    """
    filters = {
        'callrate_less_than_0.85': ht.callrate < 0.85,
        'contamination_greater_than_0.05': ht.contamination > 0.05,  # TODO determine thresholds for contamination and chimera
        'chimera_greater_than_0.05': ht.chimera > 0.05
    }
    if data_type == 'WES':
        filters.update({
            'coverage_less_than_85_at_20x': ht.bases_at_20x < 0.85
        })
    else:
        filters.update({
            'coverage_less_than_30x': ht.mean_target_coverage < 30
        })

    return hl.set(hl.filter(lambda x: hl.is_defined(x),
                            [hl.or_missing(filter_expr, name) for name, filter_expr in filters.items()]))


def compute_sex(mt: hl.MatrixTable, data_type: str, is_external: bool, version: int, test: bool,
                male_threshold: float, female_threshold: float) -> hl.MatrixTable:
    """
    Imputes sex, exports data, and annotates mt with this data # TODO evaluate F cut offs for females and males - run internal and external
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param string data_type: WES or WGS
    :param boolean is_external: Whether data is external
    :param integer version: Callset version number
    :param string male_threshold: determined from r plots
    :param string female_threshold: determined from r plots
    :param boolean test: boolean for test data
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable
    """

    # filter to the X chromosome
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval('X')])

    # create a sex_ht using hail’s impute sex function and export the resulting table
    sex_ht = hl.impute_sex(mt.GT, aaf_threshold=0.05, female_threshold=female_threshold, male_threshold=male_threshold)
    sex_ht.export(sex_check_path(data_type, is_external, version, test))

    # select three columns from the sex_ht and key by sample
    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])

    if test:
        y_coverage_ht = hl.import_table('gs://seqr-datasets/methods_dev/test_data/GRCh37/topf_y_coverage.txt',
                                        impute=True).key_by('sample_id')
        mt = mt.annotate_cols(**y_coverage_ht[mt.col_key])

    return mt


def annotate_sex(mt: hl.MatrixTable):
    """
    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable
    """
    mt = mt.annotate_cols(ambiguous_sex=((mt.f_stat >= 0.5) & (hl.is_defined(mt.normalized_y_coverage) &
                                                               (mt.normalized_y_coverage <= 0.1))) |
                                        (hl.is_missing(mt.f_stat)) |
                                        ((mt.f_stat >= 0.4) & (mt.f_stat <= 0.6) &
                                         (hl.is_defined(mt.normalized_y_coverage) &
                                          (mt.normalized_y_coverage > 0.1))),
                          sex_aneuploidy=(mt.f_stat < 0.4) & hl.is_defined(mt.normalized_y_coverage) &
                                         (mt.normalized_y_coverage > 0.1))
    sex_expr = hl.cond(mt.ambiguous_sex, "ambiguous_sex",
                       hl.cond(mt.sex_aneuploidy, "sex_aneuploidy", hl.cond(mt.is_female, "female", "male")))
    mt = mt.annotate_cols(sex=sex_expr)
    return mt


def main(args):
    hl.init()

    data_type = "WGS" if args.genomes else "WES"
    build = args.build
    is_external = True if args.is_external else False
    version = args.callset_version
    test = True if args.test else False
    male_threshold = args.male_threshold
    female_threshold = args.female_threshold

    if not args.skip_write_qc_mt:
        vcf = callset_vcf_path(data_type, build, is_external, test)
        hl.import_vcf(vcf, force_bgz=True, reference_genome=f'GRCh{build}',
                      min_partitions=4).write(mt_path(data_type, is_external, version, test), args.overwrite)
        mt = hl.read_matrix_table(mt_path(data_type, is_external, version, test))
        logger.info("Filtering to bi-allelic, high-callrate, common SNPs for sample QC...")
        mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
                            & (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
                            (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
        mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT))).naive_coalesce(5000).write(
            qc_mt_path(data_type, is_external, version, test), args.overwrite)

    qc_mt = hl.read_matrix_table(qc_mt_path(data_type, is_external, version, test))

    # optional flag for path that will override default logic
    logger.info("Importing and annotating with metadata...")
    meta_ht = hl.import_table(metadata_path(data_type, is_external, version, test),
                              impute=True).key_by('sample')
    qc_mt = qc_mt.annotate_cols(**meta_ht[qc_mt.s])

    logger.info("Importing and annotating seqr ID names...")
    remap_ht = hl.import_table(remap_path(data_type, is_external, version, test),
                               impute=True).key_by('sample_id')
    qc_mt = qc_mt.annotate_cols(**remap_ht[qc_mt.s])
    remap_expr = hl.cond(hl.is_missing(qc_mt.seqr_id), qc_mt.s, qc_mt.seqr_id)
    qc_mt = qc_mt.annotate_cols(seqr_id=remap_expr)

    logger.info("Annotating with filter flags...")
    qc_mt = qc_mt.annotate_cols(filter_flags=apply_filter_flags_expr(qc_mt, data_type))

    logger.info("Computing and annotating sex...")
    qc_mt = compute_sex(qc_mt, data_type, is_external, version, test, male_threshold, female_threshold)
    qc_mt = annotate_sex(qc_mt)  # TODO need to add args for variable thresholds

    logger.info("Annotating with pedigree concordance...")
    ped = ped_path()  # Will we make a massive ped file?
    ped_ht = hl.import_table(ped, impute=True).key_by('individual')
    qc_mt = qc_mt.annotate_cols(**ped_ht[qc_mt.s])
    gender_expr = hl.cond(qc_mt.gender == 1, "male", hl.cond(qc_mt.gender == 2, "female", "unknown"))
    qc_mt = qc_mt.annotate_cols(ped_sex=gender_expr)
    concordance_expr = hl.cond(qc_mt.ped_sex == qc_mt.sex, True,
                               False)  # gender is what is given in pedigree and sex is what was imputed
    qc_mt = qc_mt.annotate_cols(sex_concordance=concordance_expr)

    qc_mt.export(get_filepath(data_type, is_external, version, test) + ".annotations.txt.bgz")
    qc_mt.write(qc_mt_path(data_type, is_external, version, test))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--exomes',
                        help='Input MatrixTable contains exomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('-g', '--genomes',
                        help='Input MatrixTable contains genomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--skip_write_qc_mt',
                        help='Skip writing pre-calculated MatrixTable containing high-quality variants',
                        action='store_true')

    parser.add_argument('-b', '--build', help='Reference build, 37 or 38', choices=["37", "38"], required=True)
    parser.add_argument('-o', '--output', help='Override default output directory with given directory.')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('-v', '--callset_version', help='Version of callset vcf', type=int)
    parser.add_argument('--internal', help='Internal sample data in callset', action='store_true')
    parser.add_argument('--external', help='External sample data in callset', action='store_true')
    parser.add_argument('--test', help='To run a test of the pipeline using test files and directories',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite previous QC for callset version', actions='store_true')
    parser.add_argument('--male_threshold', help='Male f threshold for computing sex', type=float, default=0.8)
    parser.add_argument('--female_threshold', help='Female f threshold for computing sex', type=float, default=0.5)  # TODO Add in thresholds for annotate_sex

    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
