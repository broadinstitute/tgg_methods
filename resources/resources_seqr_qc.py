def callset_vcf_path(data_type: str, external: bool, version: int, test: bool) -> str:
    """
    Returns callset vcf path. Can be internal or external, exomes or genomes.
    """
    #Need to work version in because if overwrite - need to reload every project everytime so return .../v' + version + "/callset_Vcf"
    #Need callsets to be named the same everytime...should it just grab whatever vcf is in the directory? or fix naming with mv in bucket?
    if test:
        return f'gs://seqr-datasets/methods_dev/test_data/GRCh37/CMG_MYOSEQ.vcf.gz'
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as either 'exomes' or genomes'")
    if data_type == 'exomes' and external:
        return f'gs://seqr-datasets/GRCh37/MacArthur_WGS_Broad_External/v{version}/MacArthur_WGS_Broad_External.vcf.gz'
    if data_type == 'exomes' and not external:
        return f'gs://seqr-datasets/GRCh37/CMG_MYOSEQ/v{version}/CMG_MYOSEQ.vcf.gz'
    if data_type == 'genomes' and external:
        return f'gs://seqr-datasets/GRCh37/MacArthur_WGS_Broad_External/v{version}/MacArthur_WGS_Broad_External.vcf.gz'
    if data_type == 'genomes' and not external:
        return f'gs://seqr-datasets/GRCh38/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets/v{version}/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets_v{version}.vcf.gz'


def mt_path(data_type: str, external: bool, version: int, test: bool) -> str:
    """
    Returns MatrixTable for seqr sample qc purposes: can be exomes or genomes, internal or external data.
    Generated from vcf.
    """
    if test:
        return f'gs://seqr-datasets/methods_dev/test_data/GRCh37/topf_test_seqr_qc.mt'
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as either 'exomes' or genomes'")
    if data_type == 'exomes' and external:
        return f'gs://seqr-datasets/GRCh37/MacArthur_WGS_Broad_External/v{version}/MacArthur_WGS_Broad_External.mt'
    if data_type == 'exomes' and not external:
        return f'gs://seqr-datasets/GRCh37/CMG_MYOSEQ/v{version}/CMG_MYOSEQ.mt'
    if data_type == 'genomes' and external:
        return f'gs://seqr-datasets/GRCh37/MacArthur_WGS_Broad_External/v{version}/MacArthur_WGS_Broad_External.mt'
    if data_type == 'genomes' and not external:
        return f'gs://seqr-datasets/GRCh38/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets/v{version}/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callset.mt'


def qc_mt_path(data_type: str, external: bool, version: int, test: bool) -> str:
    """
    Returns seqr sample qc MatrixTable: can be exomes or genomes, internal or external data.
    Contains vcf sample ID, seqr sample ID, callrate, chimera, contamination, coverage, pedigree information,
    provided gender and affected status, computed gender, platform, ancestry, inbreeding coefficient.
    """
    if test:
        return f'gs://seqr-datasets/methods_dev/test_data/GRCh37/topf_test_seqr_qc.mt'
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as either 'exomes' or genomes'")
    if data_type == 'exomes' and external:
        return f'gs://seqr-datasets/GRCh37/MacArthur_WGS_Broad_External/v{version}/MacArthur_WGS_Broad_External_seqr_qc.mt'
    if data_type == 'exomes' and not external:
        return f'gs://seqr-datasets/GRCh37/CMG_MYOSEQ/v{version}/CMG_MYOSEQ_seqr_qc.mt'
    if data_type == 'genomes' and external:
        return f'gs://seqr-datasets/GRCh37/MacArthur_WGS_Broad_External/v{version}/MacArthur_WGS_Broad_External_seqr_qc.mt'
    if data_type == 'genomes' and not external:
        return f'gs://seqr-datasets/GRCh38/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets/v{version}/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets_v{version}_seqr_qc.mt'


def metadata_path(data_type: str, external: bool, version: int, test: bool) -> str:
    """
    Path to metadata file associated with samples in the callset for seqr sample QC
    """
    if test:
        return f'gs://seqr-datasets/methods_dev/test_data/GRCh37/topf_meta.txt'
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as either 'exomes' or genomes'")
    if data_type == 'exomes' and external:
        return f'gs://seqr-datasets/GRCh37/MacArthur_WGS_Broad_External/v{version}/MacArthur_WGS_Broad_External_metadata.txt'
    if data_type == 'exomes' and not external:
        return f'gs://seqr-datasets/GRCh37/CMG_MYOSEQ/v{version}/CMG_MYOSEQ_metadata.txt'
    if data_type == 'genomes' and external:
        return f'gs://seqr-datasets/GRCh37/MacArthur_WGS_Broad_External/v{version}/MacArthur_WGS_Broad_External_metadata.txt'
    if data_type == 'genomes' and not external:
        return f'gs://seqr-datasets/GRCh38/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets/v{version}/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets_v{version}_metadata.txt'


def remap_path(data_type: str, external: bool, test: bool) -> str:
    """
    Path to remap ID file
    """
    if test:
        return f'gs://seqr-datasets/methods_dev/test_data/GRCh37/topf_remap.txt'
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as either 'exomes' or genomes'")
    if data_type == 'exomes' and external:
        return f'gs://seqr-datasets/callset_id_remap_files/MacArthur_WGS_Broad_External_remap.txt'
    if data_type == 'exomes' and not external:
        return f'gs://seqr-datasets/callset_id_remap_files/CMG_MYOSEQ_remap.txt'
    if data_type == 'genomes' and external:
        return f'gs://seqr-datasets/callset_id_remap_files/MacArthur_WGS_Broad_External_remap.txt'
    if data_type == 'genomes' and not external:
        return f'gs://seqr-datasets/callset_id_remap_files/CMG_RGP_Broad_MacArthur_RareDisease_WGS_Callsets_remap.txt'


def sex_check_name():
    return "test"


# Using full ped right now? Do we need sex check concordance right away or can it happen after subset?
def ped_path():
    return "gs://seqr-dbgap/cmg_topf_cms_wes.ped"


class DataException(Exception):
    pass