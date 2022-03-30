class VCFDataTypeError(Exception):
    def __init__(self, message):
        super().__init__(message)


def temp_sample_qc_folder(
    build: int, data_type: str, data_source: str, version: int, is_test=False
) -> str:
    """
    Returns string with file path to callset version
    :param build: The genome build: "38" or "37"
    :param data_type: String of WES or WGS
    :param data_source: Data is External or Internal
    :param version: int version of callset
    :param is_test: Used if running pipeline on designated test data
    :return:
    """
    if is_test:
        return f"gs://seqr-temp/v02/GRCh{build}/RDG_{data_type}_Broad_{data_source}/v{version}/sample_qc/test"  # TODO revisit where test should be - want to use same resource folder to avoid copying over agian
    else:
        return f"gs://seqr-temp/v02/GRCh{build}/RDG_{data_type}_Broad_{data_source}/v{version}/sample_qc"


def final_sample_qc_folder(
    build: int, data_type: str, data_source: str, version: int
) -> str:
    """
    Returns string with file path to callset version
    :param build: The genome build: "38" or "37"
    :param data_type: String of WES or WGS
    :param data_source: Data is External or Internal
    :param version: int version of callset
    :param is_test: Used if running pipeline on designated test data
    :return:
    """
    return f"gs://seqr-datasets/v02/GRCh{build}/RDG_{data_type}_Broad_{data_source}/v{version}/sample_qc"


def mt_path(
    build: int, data_type: str, data_source: str, version: int, is_test: bool
) -> str:
    """
    Returns MatrixTable for seqr sample qc purposes: can be exomes or genomes, internal or external data.
    Generated from vcf.
    """
    return f"{temp_sample_qc_folder(build, data_type, data_source, version, is_test)}/resources/callset.mt"


def seq_metrics_path(build: int, data_type: str, data_source: str, version: int) -> str:
    """
    Path to metadata file associated with samples in the callset for seqr sample QC.
    """
    return f"{final_sample_qc_folder(build, data_type, data_source, version)}/resources/callset_seq_metrics.tsv"


def remap_path(build: int, data_type: str, data_source: str, version: int) -> str:
    """
    Path to remap ID file which needs to be stored in the callset's bucket prior to launching pipeline
    """
    return f"{final_sample_qc_folder(build, data_type, data_source, version)}/resources/remap.tsv"


def missing_metrics_path(
    build: int, data_type: str, data_source: str, version: int
) -> str:
    """
    Path to file containing sample IDs that do not have any seq metrics.
    """
    return f"{temp_sample_qc_folder(build, data_type, data_source, version)}/missing_metrics/samples_missing_metrics.tsv"


def mt_temp_path(
    build: int, data_type: str, data_source: str, version: int, is_test: bool
) -> str:
    """
    Temporary MatrixTable path that should be deleted once pipeline is complete
    """
    return f"{temp_sample_qc_folder(build, data_type, data_source, version, is_test)}/temp/mt_temp.mt"


def sample_qc_tsv_path(
    build: int, data_type: str, data_source: str, version: int, is_test: bool
) -> str:
    """
    Returns seqr sample qc Table: can be exomes or genomes, internal or external data.
    Contains vcf sample ID, seqr sample ID, callrate, chimera, contamination, coverage, pedigree information,
    provided gender and affected status, computed gender, platform, ancestry, inbreeding coefficient.
    """
    return f"{final_sample_qc_folder(build, data_type, data_source, version, is_test)}/final_output/seqr_sample_qc.tsv"


def sample_qc_ht_path(
    build: int, data_type: str, data_source: str, version: int, is_test: bool
) -> str:
    """
    Returns seqr sample qc hail table: can be exomes or genomes, internal or external data.
    Contains vcf sample ID, seqr sample ID, callrate, chimera, contamination, coverage, pedigree information,
    provided gender and affected status, computed gender, platform, ancestry, inbreeding coefficient.
    """
    return f"{temp_sample_qc_folder(build, data_type, data_source, version, is_test)}/final_output/seqr_sample_qc.ht"


def rdg_gnomad_pop_pca_loadings_ht_path(build: int) -> str:
    """Return the precomputed PCA loadings from joint RDG and gnomAD PCA"""
    if build == 37:
        return "gs://seqr-datasets/sample_qc_resources/population_assignment/37_ancestry_pca_loadings.ht"
    else:
        return (
            "gs://gcp-public-data--gnomad/release/3.1/pca/gnomad.v3.1.pca_loadings.ht"
        )


def val_noncoding_ht_path(build):
    """
    HT of noncoding variants for validating callset data type. HT written using
    hail-elasticsearch-pipelines/download_and_create_reference_datasets/v02/hail_scripts/write_dataset_validation_ht.py
    """
    return f"gs://seqr-reference-data/GRCh{build}/validate_ht/common_noncoding_variants.grch{build}.ht"


def val_coding_ht_path(build):
    """
    HT of coding variants for validating callset data type. HT written using
    hail-elasticsearch-pipelines/download_and_create_reference_datasets/v02/hail_scripts/write_dataset_validation_ht.py
    """
    return f"gs://seqr-reference-data/GRCh{build}/validate_ht/common_coding_variants.grch{build}.ht"
