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
        return f'gs://seqr-datasets/methods_dev/test_data/GRCh37/topf'  # TODO rename test vcf
    else:
        return f'gs://seqr-datasets/{build}/RDG_{data_type}_Broad_{origin}/{version}/RDG_{data_type}_Broad_{origin}_{version}'


def callset_vcf_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Returns callset vcf path. Can be internal or external, exomes or genomes.
    """
    filepath = get_filepath(data_type, is_external, version, test)
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


def metadata_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Path to metadata file associated with samples in the callset for seqr sample QC
    """
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_metadata.txt'


def remap_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """
    Path to remap ID file
    """
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_remap.txt'


def sex_check_path(data_type: str, is_external: bool, version: int, test: bool) -> str:
    """Return file path for sex check"""
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '.sex_check.txt.bgz'


def sex_check_name():
    return "test"


def ped_path():
    return "gs://seqr-dbgap/cmg_topf_cms_wes.ped"  # TODO Use Ben's API to pull ped of all individuals in callset from seqr


class DataException(Exception):
    pass


def exome_callrate_mt_path(is_external: bool, version: int, test: bool) -> str:
    data_type = "WES"
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_callrate_mt_for_pca.mt'


def exome_callrate_scores_ht_path(is_external: bool, version: int, test: bool) -> str:
    data_type = "WES"
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_callrate_pca_scores.ht'


def qc_temp_data(data_type: str, is_external: bool, test: bool, version: int) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + '_platform_pca/'


def platform_labels_path(data_type: str, is_external: bool, test: bool, version: int) -> str:
    filepath = get_filepath(data_type, is_external, version, test)
    return filepath + 'platforms.txt'


class DataException(Exception):
    pass