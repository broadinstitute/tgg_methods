class VCFDataTypeError(Exception):
    def __init__(self, message):
        super().__init__(message)


def get_filepath(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Returns string with file path to callset version
    :param data_type: String of WES or WGS
    :param data_source: Data is External or Internal
    :param version: int version of callset
    :param is_test: Used if running pipeline on designated test data
    :return:
    """
    if is_test:
        return f'gs://seqr-datasets/v02/GRCh{build}/RDG_{data_type}_Broad_{data_source}/v{version}/sample_qc/test/' # TODO revisit where test should be - want to use same resource folder to avoid copying over agian
    else:
        return f'gs://seqr-datasets/v02/GRCh{build}/RDG_{data_type}_Broad_{data_source}/v{version}/sample_qc/'


def callset_vcf_path(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Returns callset vcf path. Can be internal or external, exomes or genomes.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test)
    if is_test:
        filepath = filepath[:-5]
    filepath = filepath[:-10] + f'RDG_{data_type}_Broad_{data_source}.vcf.bgz'
    return filepath


def mt_path(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Returns MatrixTable for seqr sample qc purposes: can be exomes or genomes, internal or external data.
    Generated from vcf.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test)
    return filepath + f'resources/callset.mt'


def seq_metrics_path(build: str, data_type: str, data_source: str, version: int) -> str:
    """
    Path to metadata file associated with samples in the callset for seqr sample QC.
    Set is_test to false because resource files do not change during testing.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test=False)
    return filepath + 'resources/callset_seq_metrics.txt'


def remap_path(build: str, data_type: str, data_source: str, version: int) -> str:
    """
    Path to remap ID file which needs to be stored in the callset's bucket prior to launching pipeline
    Set is_test to false because resource files do not change during testing.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test=False)
    return filepath + 'resources/remap.tsv'


def project_map_path(build: str, data_type: str, data_source: str, version: int) -> str:
    """
    Path to project mapping file which needs to be stored in the callset's bucket prior to launching pipeline
    Set is_test to false because resource files do not change during testing.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test=False)
    return filepath + 'resources/project_map.tsv'


def missing_metrics_path(build: str, data_type: str, data_source: str, version: int) -> str:
    """
    Path to file containing sample IDs that do not have any seq metrics.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test=False)
    return filepath + 'missing_metrics/samples_missing_metrics.tsv'


def ped_path(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Return the pedigree filepath
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test)
    return filepath + 'resources/callset.ped'


def sex_check_path(build: str, data_type: str, data_source: str, version: int) -> str:
    """
    Return file path for sex check. Set is_test to false because resource files do not change during testing.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test=False)
    return filepath + 'sex_check/imputed_sex.tsv'


def rdg_gnomad_pop_pca_loadings_ht_path(build: str) -> str:
    """Return the PCA loadings from joint RDG and gnomAD PCA"""
    return f'gs://seqr-datasets/v02/GRCh{build}/sample_qc/population_assignment/ancestry_pca_loadings.ht'


def pop_RF_fit_path(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Return the RandomForestClassifier path
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test)
    return filepath + 'pop_imputation/RF_population_fit.pkl'


def mt_temp_path(build: str, data_type: str,  data_source: str, version: int, test: bool) -> str:
    """
    Temporary MatrixTable path that should be deleted once pipeline is complete
    """
    filepath = get_filepath(build, data_type, data_source, version, test)
    return filepath + 'temp/mt_temp.mt'


def ht_to_tsv_path(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Returns seqr sample qc Table: can be exomes or genomes, internal or external data.
    Contains vcf sample ID, seqr sample ID, callrate, chimera, contamination, coverage, pedigree information,
    provided gender and affected status, computed gender, platform, ancestry, inbreeding coefficient.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test)
    return filepath + f'final_output/seqr_sample_qc.tsv'


def sample_qc_ht_path(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Returns seqr sample qc hail table: can be exomes or genomes, internal or external data.
    Contains vcf sample ID, seqr sample ID, callrate, chimera, contamination, coverage, pedigree information,
    provided gender and affected status, computed gender, platform, ancestry, inbreeding coefficient.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test)
    return filepath + f'final_output/seqr_sample_qc.ht'


def val_noncoding_ht_path(build):
    """
    HT of noncoding variants for validating callset data type. HT written using
    hail-elasticsearch-pipelines/download_and_create_reference_datasets/v02/hail_scripts/write_dataset_validation_ht.py
    """
    if build == '37':
        return 'gs://seqr-reference-data/GRCh37/validate_ht/common_noncoding_variants.grch37.ht'
    else:
        return 'gs://seqr-reference-data/GRCh38/validate_ht/common_noncoding_variants.grch38.ht'


def val_coding_ht_path(build):
    """
    HT of coding variants for validating callset data type. HT written using
    hail-elasticsearch-pipelines/download_and_create_reference_datasets/v02/hail_scripts/write_dataset_validation_ht.py
    """
    if build == '37':
        return 'gs://seqr-reference-data/GRCh37/validate_ht/common_coding_variants.grch37.ht'
    else:
        return 'gs://seqr-reference-data/GRCh38/validate_ht/common_coding_variants.grch38.ht'
