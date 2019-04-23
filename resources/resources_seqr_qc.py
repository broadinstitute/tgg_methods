class DataException(Exception):
    pass


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
        return f'gs://seqr-datasets/GRCh{build}/RDG_{data_type}_Broad_{data_source}/v{version}/RDG_{data_type}_Broad_{data_source}/sample_qc/test/' # TODO revisit where test should be - want to use same resource folder to avoid copying over agian
    else:
        return f'gs://seqr-datasets/GRCh{build}/RDG_{data_type}_Broad_{data_source}/v{version}/RDG_{data_type}_Broad_{data_source}/sample_qc/'


def callset_vcf_path(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Returns callset vcf path. Can be internal or external, exomes or genomes.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test)
    if is_test:
        filepath = filepath[:-5]
    if data_type == 'WGS' and data_source == 'internal':
        filepath = filepath[:-11] + f'/sharded-vcfs/RDG_{data_type}_Broad_Internal.filtered.*'
    else:
        filepath = filepath[:-11] + '.vcf.gz'
    return filepath


def mt_path(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Returns MatrixTable for seqr sample qc purposes: can be exomes or genomes, internal or external data.
    Generated from vcf.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test)
    return filepath + f'resources/callset.mt'


def liftover_mt_path(build: str, data_type: str, data_source: str, version: int, is_test: bool) -> str:
    """
    Returns seqr sample qc MatrixTable: can be exomes or genomes, internal or external data.
    Contains vcf sample ID, seqr sample ID, callrate, chimera, contamination, coverage, pedigree information,
    provided gender and affected status, computed gender, platform, ancestry, inbreeding coefficient.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test)
    return filepath + f'resources/callset_37_liftover.mt'


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
    return filepath + 'resources/remap.txt'


def project_map_path(build: str, data_type: str, data_source: str, version: int) -> str:
    """
    Path to project mapping file which needs to be stored in the callset's bucket prior to launching pipeline
    Set is_test to false because resource files do not change during testing.
    """
    filepath = get_filepath(build, data_type, data_source, version, is_test=False)
    return filepath + 'resources/project_map.txt'


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
    return filepath + 'sex_check/imputed_sex.txt'


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
