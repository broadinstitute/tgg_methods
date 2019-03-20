def get_filepath(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Returns string with file path to callset version
    :param data_type: String of WES or WGS
    :param is_external: boolean whether data is external or not
    :param version: int version of callset
    :param test: Used if running pipeline on designated test data
    :return:
    """
    origin = "External" if is_external else "Internal"
    build = "GRCh38" if data_type == "WGS" and not is_external else "GRCh37"
    if test:
        return f'gs://seqr-datasets/{build}/RDG_{data_type}_Broad_{origin}/test/v{version}/RDG_{data_type}_Broad_{origin}'
    else:
        return f'gs://seqr-datasets/{build}/RDG_{data_type}_Broad_{origin}/v{version}/RDG_{data_type}_Broad_{origin}'


def callset_vcf_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Returns callset vcf path. Can be internal or external, exomes or genomes.
    """
    filepath = get_filepath(data_type, is_external, version, test)
    if data_type == "WGS" and not is_external:
        filepath = filepath + f'/sharded-vcfs/RDG_{data_type}_Broad_Internal.filtered.*'
    return filepath + '.vcf.gz'


def mt_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Returns MatrixTable for seqr sample qc purposes: can be exomes or genomes, internal or external data.
    Generated from vcf.
    """
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '.mt'


def qc_mt_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Returns seqr sample qc MatrixTable: can be exomes or genomes, internal or external data.
    Contains vcf sample ID, seqr sample ID, callrate, chimera, contamination, coverage, pedigree information,
    provided gender and affected status, computed gender, platform, ancestry, inbreeding coefficient.
    """
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_seqr_qc.mt'


def qc_ht_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Returns seqr sample qc Table: can be exomes or genomes, internal or external data.
    Contains vcf sample ID, seqr sample ID, callrate, chimera, contamination, coverage, pedigree information,
    provided gender and affected status, computed gender, platform, ancestry, inbreeding coefficient.
    """
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_seqr_qc_samples.ht'


def metadata_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Path to metadata file associated with samples in the callset for seqr sample QC
    """
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_metadata.txt'


def remap_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Path to remap ID file which needs to be stored in the callset's bucket prior to launching pipeline
    """
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_remap.txt'


def project_map_path(data_type: str, is_external: bool, version: int, test: bool)->str:
    """
    Path to project mapping file which needs to be stored in the callset's bucket prior to launching pipeline
    """
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_project_map.txt'


def sex_check_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Return file path for sex check
    """
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '/sex.txt'


def ped_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '.ped'


class DataException(Exception):
    pass


def exome_callrate_mt_path(is_external: bool, version: int, test: bool) -> str:
    data_type = "WES"
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_callrate_mt_for_pca.mt'


def exome_platform_callrate_scores_ht_path(is_external: bool, version: int, test: bool) -> str:
    data_type = "WES"
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_platform_callrate_pca_scores.ht'


def population_assignments_ht_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '/pop_imputation/RF_pop_assignments.txt.bgz'


def population_assignments_tsv_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '/pop_imputation/RF_pop_assignments.tsv'


def population_RF_fit_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '/pop_imputation/RF_population_fit.pkl'


def qc_temp_data(data_type: str, is_external: bool, version: int, test: bool) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_platform_pca/' #TODO:Do I use temp anywhere else?

def qc_mt_temp_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + f'_temp/qc_mt_temp.mt' #TODO:Do I use temp anywhere else?

def platform_labels_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_platform_pca/platforms.txt'

def metric_MAD_data_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '.platform_pop_MAD.txt.bgz'

def test():
    print('In resources test method')