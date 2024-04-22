import argparse
import logging
from typing import Dict, List, Optional, Union
from copy import deepcopy

import hail as hl
from gnomad.resources.grch37.gnomad import public_release as v2_public_release
from gnomad.resources.grch38.gnomad import public_release as v3_public_release
from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.relatedness import UNRELATED, get_relationship_expr
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.vep import CSQ_ORDER

import gnomad_qc.v2.resources as v2
import gnomad_qc.v2.resources.annotations as v2_annotations
import gnomad_qc.v3.resources as v3
import gnomad_qc.v4.resources as v4

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sanna_variant_script")
logger.setLevel(logging.INFO)

########################################################################################
# Global variables that need to be updated when new gnomAD versions are added.
########################################################################################
VERSIONS = {"v2_exomes", "v2_genomes", "v3", "v3.1", "v4_exomes", "v4_genomes"}
"""Set of gnomAD versions."""

V2_VERSIONS = {"v2_exomes", "v2_genomes"}
"""Set of gnomAD v2 versions."""
V3_VERSIONS = {"v3", "v3.1"}
"""Set of gnomAD v3 versions."""
V4_VERSIONS = {"v4_exomes", "v4_genomes"}
"""Set of gnomAD v4 versions."""

VERSION_GENOME_BUILD = {v: "GRCh37" if "v2" in v else "GRCh38" for v in VERSIONS}
"""Mapping of gnomAD versions to their respective genome builds."""

GENETIC_ANCESTRY_LABEL = {
    "v2_exomes": {"gen_anc": "pop", "grpmax": "popmax"},
    "v2_genomes": {"gen_anc": "pop", "grpmax": "popmax"},
    "v3": {"gen_anc": "pop", "grpmax": "popmax"},
    "v3.1": {"gen_anc": "pop", "grpmax": "popmax"},
    "v4_exomes": {"gen_anc": "gen_anc", "grpmax": "grpmax"},
    "v4_genomes": {"gen_anc": "gen_anc", "grpmax": "grpmax"},
}
"""Mapping of gnomAD versions to their respective genetic ancestry and grpmax labels."""

VERSION_RESOURCE_MAP = {
    "v2_exomes": {"data_type": "exomes"},
    "v2_genomes": {"data_type": "genomes"},
    "v3": {
        "version": "3",
        "release_version": "3.0",
        "freq_version": "3.0",
        "vep_version": "3.1.1",
        "data_type": "genomes",
    },
    "v3.1": {
        "version": "3.1",
        "vep_version": "3.1.1",
        "data_type": "genomes",
    },
    "v4_exomes": {
        "version": v4.constants.CURRENT_VERSION,
        "freq_version": v4.constants.CURRENT_FREQ_VERSION["exomes"],
        "data_type": "exomes",
        "filter_version": v4.constants.CURRENT_VARIANT_QC_RESULT_VERSION["exomes"],
        "release_version": v4.constants.CURRENT_RELEASE,
        "vep_version": v4.constants.CURRENT_ANNOTATION_VERSION,
        "info_version": v4.constants.CURRENT_ANNOTATION_VERSION,
        "meta_version": v4.constants.CURRENT_ANNOTATION_VERSION,
    },
    "v4_genomes": {
        "version": v4.constants.CURRENT_VERSION,
        "freq_version": v4.constants.CURRENT_FREQ_VERSION["genomes"],
        "data_type": "genomes",
        "info_version": "3.1",
        "filter_version": v4.constants.CURRENT_VARIANT_QC_RESULT_VERSION["genomes"],
        "release_version": v4.constants.CURRENT_RELEASE,
        "vep_version": v4.constants.CURRENT_ANNOTATION_VERSION,
        "meta_version": v4.constants.CURRENT_ANNOTATION_VERSION,
    },
}
"""Mapping of gnomAD versions to their respective resource datatypes and versions."""

GENE_REGION_GTF_MAP = {
    "v2_exomes": "gs://gcp-public-data--gnomad/resources/grch37/gencode/gencode.v19.annotation.gtf.gz",
    "v2_genomes": "gs://gcp-public-data--gnomad/resources/grch37/gencode/gencode.v19.annotation.gtf.gz",
    "v3": "gs://hail-common/references/gencode/gencode.v29.annotation.gtf.bgz",
    "v3.1": "gs://gnomad/resources/gencode/gencode.v35.annotation.gtf.gz",
    "v4_exomes": "gs://gcp-public-data--gnomad/resources/grch38/gencode/gencode.v39.annotation.gtf.gz",
    "v4_genomes": "gs://gcp-public-data--gnomad/resources/grch38/gencode/gencode.v39.annotation.gtf.gz",
}

SAMPLE_META_MAPPING = {
    "pop": {
        "v2_exomes": "pop",
        "v2_genomes": "pop",
        "v3": "pop",
        "v3.1": "pop",
        "v4_exomes": "pop",
        "v4_genomes": "pop",
    },
    # SG Added columns
    # release info
    "release": {
        "v2_exomes": "release",
        "v2_genomes": "release",
        "v3": "release",
        "v3.1": "release",
        "v4_exomes": "release",
        "v4_genomes": "release",
    },
    # cohort
    # neuro
    "neuro": {
        "v2_exomes": "neuro",
        "v2_genomes": "neuro",
        "v3": "v2_neuro",
        "v3.1": "non_neuro",
    },
    # control
    "control": {
        "v2_exomes": "control",
        "v2_genomes": "control",
        "v3": "v2_control",
        "v3.1": "controls_and_biobanks",
    },
    # topmed
    "topmed": {
        "v2_exomes": "topmed",
        "v2_genomes": "topmed",
        "v3": "v2_topmed",
        "v3.1": "non_topmed",
    },
    # non_v2
    "non_v2": {"v3.1": "non_v2"},
    # non_cancer
    "non_cancer": {
        "v2_exomes": "non_cancer",
        "v2_genomes": "non_cancer",
        "v3": "non_cancer",
        "v3.1": "non_cancer",
    },
    # tgp
    "tgp": {"v3.1": "tgp", "v4_genomes": "tgp"},
    # hgdp
    "hgdp": {"v3.1": "hgdp", "v4_genomes": "hgdp"},
    # project_id
    "project_id": {
        "v2_exomes": "project_id",
        "v2_genomes": "project_id",
        "v3": "project_id",
        "v3.1": "project_id",
        "v4_exomes": "project",
        "v4_genomes": "project_id",
    },
    # investigator or PI
    "investigator": {
        "v2_exomes": "investigator",
        "v2_genomes": "investigator",
        "v3": "contact_pi",
        "v3.1": "contact_pi",
        "v4_exomes": "investigator",
        "v4_genomes": "contact_pi",
    },
    # name of the research project
    "project_name": {
        "v2_genomes": "project_name",
        "v3": "research_project",
        "v3.1": "research_project",
        "v4_genomes": "research_project",
    },
    # name of research project extended
    "project_description": {
        "v2_exomes": "project_description",
        "v2_genomes": "project_description",
        "v3": "title",
        "v3.1": "title",
        "v4_exomes": "cohort",
        "v4_genomes": "title",
    },
    # high quality
    "high_quality": {
        "v2_exomes": "high_quality",
        "v2_genomes": "high_quality",
        "v3": "high_quality",
        "v3.1": "high_quality",
        "v4_exomes": "high_quality",
        "v4_genomes": "high_quality",
    },
    # sex
    "sex": {
        "v2_exomes": "sex",
        "v2_genomes": "sex",
        "v3": "sex",
        "v3.1": "sex_karyotype",
        "v4_exomes": "sex_karyotype",
        "v4_genomes": "sex_karyotype",
    },
    # age
    "age": {
        "v2_exomes": "age",
        "v2_genomes": "age",
        "v3": "v2_age",
        "v3.1": "age",
        "v4_exomes": "age",
        "v4_genomes": "age",
    },
}
"""Mapping of sample meta fields to their respective fields in the gnomAD versions."""

########################################################################################
# Global variables that are relevant to specific versions -- shouldn't need to be
# updated with version additions.
########################################################################################
V3_EXCLUDE_INFO_FIELDS = {
    "AS_SB_TABLE",
    "AS_RAW_MQ",
    "AS_MQ_DP",
}
"""
Fields to exclude from the V3 info HT.

These weren't used for filtering directly -- only for computing other fields.
"""

V4_FILTERS_INFO_FIELDS = [
    "singleton",
    "transmitted_singleton",
    "omni",
    "mills",
    "monoallelic",
    "only_het",
]
"""Info fields to include from the V4 filter HT."""

V2_INCLUDE_INFO_FIELDS = [
    "variant_type",
    "allele_type",
    "n_alt_alleles",
    "was_mixed",
    "has_star",
    "qd",
    "pab_max",
    "info_MQRankSum",
    "info_SOR",
    "info_InbreedingCoeff",
    "info_ReadPosRankSum",
    "info_FS",
    "info_QD",
    "info_MQ",
    "info_DP",
    "transmitted_singleton",
    "fail_hard_filters",
    "info_POSITIVE_TRAIN_SITE",
    "info_NEGATIVE_TRAIN_SITE",
    "omni",
    "mills",
    "tp",
    "rf_train",
    "rf_label",
    "rf_probability",
    "rank",
    "was_split",
]
"""Info fields to include from the V2 info HT."""

########################################################################################
# Global variables -- shouldn't need to be updated with version additions.
########################################################################################
ENTRY_FIELDS_TO_KEEP = {
    "GT": hl.expr.types.tcall,
    "GQ": hl.expr.types.tint32,
    "PID": hl.expr.types.tstr,
    "PGT": hl.expr.types.tcall,
    "AD": hl.expr.types.tarray(hl.expr.tint32),
    "PL": hl.expr.types.tarray(hl.expr.tint32),
    "DP": hl.expr.types.tint32,
    "adj": hl.expr.types.tbool,
}
"""Entry fields to keep in the gnomAD MT."""

VEP_INCLUDE_FIELDS = [
    "most_severe_consequence",
    "motif_feature_consequences",
    "regulatory_feature_consequences",
    "transcript_consequences",
]
"""VEP fields to include from the VEP HT."""

TSV_COLUMN_ORDER = [
    "s",
    "v1_GRCh38_chrom",
    "v1_GRCh38_pos",
    "v1_GRCh38_ref",
    "v1_GRCh38_alt",
    "v1_GRCh37_chrom",
    "v1_GRCh37_pos",
    "v1_GRCh37_ref",
    "v1_GRCh37_alt",
    "v2_GRCh38_chrom",
    "v2_GRCh38_pos",
    "v2_GRCh38_ref",
    "v2_GRCh38_alt",
    "v2_GRCh37_chrom",
    "v2_GRCh37_pos",
    "v2_GRCh37_ref",
    "v2_GRCh37_alt",
    "v1_GT",
    "v1_GQ",
    "v1_PID",
    "v1_PGT",
    "v1_AD",
    "v1_PL",
    "v1_DP",
    "v1_adj",
    "v2_GT",
    "v2_GQ",
    "v2_PID",
    "v2_PGT",
    "v2_AD",
    "v2_PL",
    "v2_DP",
    "v2_adj",
]
"""
Column order for the TSV export.

Columns specified here will appear first and in given order. All other columns are 
appended after.
"""


########################################################################################
# Functions that need to be updated when new gnomAD versions are added.
########################################################################################
def get_vep_ht(gnomad_version: str, intervals: hl.tarray = None) -> hl.Table:
    """
    Get the VEP HT for the given gnomAD version.

    :param gnomad_version: gnomAD version.
    :param intervals: Intervals to filter to.
    :return: VEP HT.
    """
    data_type, version = gnomad_version_to_resource_version(
        gnomad_version, vep_version=True
    )

    if gnomad_version in V3_VERSIONS:
        vep_ht = v3.release.release_sites(public=True).versions[version].ht()
        vep_ht = vep_ht.select("vep")

        # Load VEP. There are two extra fields in v3 VEP that need to be dropped so
        # that it can be merged with v2.
        vep_ht = vep_ht.annotate(
            vep=vep_ht.vep.annotate(
                transcript_consequences=vep_ht.vep.transcript_consequences.map(
                    lambda x: x.drop("appris", "tsl")
                )
            )
        )
    elif gnomad_version in V2_VERSIONS:
        vep_ht = v2_public_release(data_type).ht().select("vep")
    elif gnomad_version in V4_VERSIONS:
        vep_ht = v4.annotations.get_vep(data_type=data_type).versions[version].ht()
    else:
        raise DataException(f"Version {gnomad_version} is not supported for VEP HT")

    vep_ht = hl.filter_intervals(vep_ht, intervals) if intervals else vep_ht

    return vep_ht


def get_freq_ht(gnomad_version: str, intervals: hl.tarray = None) -> hl.Table:
    """
    Get the frequency HT for the given gnomAD version.

    :param gnomad_version: gnomAD version.
    :param intervals: Intervals to filter to.
    :return: Frequency HT.
    """
    data_type, version = gnomad_version_to_resource_version(
        gnomad_version, freq_version=True
    )

    if gnomad_version in V3_VERSIONS:
        freq_ht = v3_public_release(data_type).versions[version].ht()
        freq_ht = freq_ht.select("freq", "popmax").select_globals("freq_meta")
        if gnomad_version == "v3.1":
            freq_ht = freq_ht.annotate(popmax=freq_ht.popmax.drop("faf95"))
    elif gnomad_version in V2_VERSIONS:
        freq_ht = (
            v2_public_release(data_type)
            .ht()
            .select("freq", "popmax")
            .select_globals("freq_meta")
        )
    elif gnomad_version in V4_VERSIONS:
        freq_ht = v4.annotations.get_freq(version=version, data_type=data_type).ht()
        if data_type == "exomes":
            freq_ht = freq_ht.transmute(grpmax=freq_ht.grpmax.gnomad)
    else:
        raise DataException(
            f"Version {gnomad_version} is not supported for frequency HT"
        )

    freq_ht = hl.filter_intervals(freq_ht, intervals) if intervals else freq_ht

    return freq_ht


def get_filter_ht(gnomad_version: str, intervals: hl.tarray = None) -> hl.Table:
    """
    Get the filter HT for the given gnomAD version.

    :param gnomad_version: gnomAD version.
    :param intervals: Intervals to filter to.
    :return: Filter HT.
    """
    data_type, version = gnomad_version_to_resource_version(
        gnomad_version, filter_version=True
    )

    if gnomad_version in V3_VERSIONS:
        try:
            filter_ht = (
                v3.annotations.get_vqsr_filters(model_id="vqsr_alleleSpecificTrans")
                .versions[version]
                .ht()
            )
        except FileNotFoundError:
            try:
                filter_ht = (
                    v3.annotations.get_vqsr_filters(
                        model_id="vqsr_alleleSpecificTrans", split=False
                    )
                    .versions[version]
                    .ht()
                )
                filter_ht = split_vqsr(filter_ht)
            except FileNotFoundError:
                raise DataException(
                    "There is no available split or unsplit info HT for use!"
                )
    elif gnomad_version in V2_VERSIONS:
        filter_ht = v2_annotations.rf(data_type).ht().select("filters")
    elif gnomad_version in V4_VERSIONS:
        filter_ht = (
            v4.variant_qc.final_filter(data_type=data_type).versions[version].ht()
        )
        filter_ht = filter_ht.transmute(**filter_ht.truth_sets)
    else:
        raise DataException(f"Version {gnomad_version} is not supported for filter HT")

    filter_ht = hl.filter_intervals(filter_ht, intervals) if intervals else filter_ht

    return filter_ht


def get_info_ht(gnomad_version: str, intervals: hl.tarray = None) -> hl.Table:
    """
    Get the info HT for the given gnomAD version.

    :param gnomad_version: gnomAD version.
    :param intervals: Intervals to filter to.
    :return: Info HT.
    """
    data_type, version = gnomad_version_to_resource_version(
        gnomad_version, info_version=True
    )

    if gnomad_version in {"v3", "v3.1", "v4_genomes"}:
        try:
            info_ht = v3.annotations.get_info().versions[version].ht()
        except FileNotFoundError:
            try:
                info_ht = v3.annotations.generate_qc_annotations.split_info()
                info_ht = info_ht.drop("old_locus", "old_alleles")
            except FileNotFoundError:
                raise DataException(
                    "There is no available split or unsplit info HT for use!"
                )
    elif gnomad_version in V2_VERSIONS:
        info_ht = v2_annotations.rf(data_type).ht()
    elif gnomad_version in {"v4_exomes"}:
        info_ht = v4.annotations.get_info().versions[version].ht()
    else:
        raise DataException(f"Version {gnomad_version} is not supported for info HT")

    info_ht = hl.filter_intervals(info_ht, intervals) if intervals else info_ht

    return info_ht


def get_gnomad_raw_data(
    version: str,
    intervals: Union[List[hl.Interval], hl.tarray] = None,
    samples_ht: Optional[hl.Table] = None,
    ref: bool = False,
    densify: bool = False,
    loci_filter_ht: Optional[hl.Table] = None,
) -> hl.MatrixTable:
    """
    Get the raw gnomAD data for the given version.

    This will also filter to the given intervals and split multi-allelics and annotate
    entries with adj where needed.

    .. warning::

        Data after v2 is represented as a Hail VariantDataset, and by default the
        returned MT will only be the filtered variant data. Since the VariantDataset
        in not densified, the MT will not have valid reference genotype information.

    :param version: gnomAD version.
    :param intervals: Intervals to filter to.
    :param samples: List of samples to filter to.
    :param ref: Whether include reference genotypes. For all versions other than v2 this
        requires a densify.
    :param densify: Whether to densify the MT even if not only getting ref GTs.
    :return: gnomAD MT.
    """
    data_type, _ = gnomad_version_to_resource_version(version)

    if isinstance(intervals, list):
        intervals = hl.literal(intervals)

    gnomad_mt = None
    vds = None
    if version in V2_VERSIONS:
        # Need to use hardcoded path because the v2 resources have been moved and not
        # updated.
        if ref:
            gnomad_mt_path = f"gs://gnomad_v2/hardcalls/hail-0.2/mt/{data_type}/gnomad.{data_type}.mt"
        else:
            gnomad_mt_path = f"gs://gnomad_v2/non_refs_only/hail-0.2/mt/{data_type}/gnomad.{data_type}.mt"

        gnomad_mt = hl.read_matrix_table(gnomad_mt_path)
        # Drop outdated GATK fields.
        gnomad_mt = gnomad_mt.drop("info")
        if intervals is not None:
            gnomad_mt = hl.filter_intervals(gnomad_mt, intervals)

        if samples_ht is not None:
            gnomad_mt = gnomad_mt.semi_join_cols(samples_ht)

        if loci_filter_ht is not None:
            gnomad_mt = gnomad_mt.filter_rows(
                hl.is_defined(loci_filter_ht[gnomad_mt.locus])
            )

    elif version in V3_VERSIONS or version == "v4_genomes":
        vds = v3.basics.get_gnomad_v3_vds(remove_hard_filtered_samples=False)
    elif version in V4_VERSIONS:
        vds = v4.basics.get_gnomad_v4_vds(
            remove_hard_filtered_samples=False, remove_dead_alleles=False
        )
    else:
        raise DataException(f"Version {version} is not supported for variants samples")

    if gnomad_mt is None and vds is not None:
        # If requested, densify the MT, otherwise use only the variant data.
        if ref or densify:
            # Filter to intervals.
            if intervals is not None:
                vds = hl.vds.filter_intervals(vds, intervals)

            if samples_ht is not None:
                vds = hl.vds.filter_samples(vds, samples_ht, keep=True)

            if loci_filter_ht is not None:
                vmt = vds.variant_data
                vds = hl.vds.VariantDataset(
                    vds.reference_data,
                    vmt.filter_rows(
                        hl.is_defined(loci_filter_ht[vmt.locus])
                    )
                )

            # Split multi-allelics.
            vds = hl.vds.split_multi(vds)

            gnomad_mt = hl.vds.to_dense_mt(vds)
        else:
            gnomad_mt = vds.variant_data

            if intervals is not None:
                gnomad_mt = hl.filter_intervals(gnomad_mt, intervals)

            if samples_ht is not None:
                gnomad_mt = gnomad_mt.semi_join_cols(samples_ht)

            if loci_filter_ht is not None:
                gnomad_mt = gnomad_mt.filter_rows(
                    hl.is_defined(loci_filter_ht[gnomad_mt.locus])
                )

            gnomad_mt = hl.experimental.sparse_split_multi(
                gnomad_mt, filter_changed_loci=True
            )
            gnomad_mt = gnomad_mt.checkpoint(
                hl.utils.new_temp_file("mt_after_split", "mt")
            )

            # This filter is required to filter refs resulting from the split.
            gnomad_mt = gnomad_mt.filter_entries(gnomad_mt.GT.is_non_ref())

        # Annotate adj.
        # TODO: DO WE NEED TO adjust ploidy?
        gnomad_mt = annotate_adj(gnomad_mt)

    if ref:
        gnomad_mt = gnomad_mt.filter_entries(gnomad_mt.GT.is_hom_ref())
        gnomad_mt = gnomad_mt.annotate_entries(
            **{
                k: hl.null(v)
                for k, v in ENTRY_FIELDS_TO_KEEP.items()
                if k not in gnomad_mt.entry
            }
        )

    return gnomad_mt


def get_gnomad_regions_mt(
    gnomad_version: str,
    samples_ht: Optional[hl.Table] = None,
    row_regions: Optional[Union[List[hl.Interval], hl.tarray]] = None,
    variant_ht: Optional[hl.Table] = None,
    ref: bool = False,
    pops: List[str] = None,
    least_consequence: str = None,
    max_af: float = None,
) -> hl.MatrixTable:
    """
    Get the gnomAD regions MT for the given version.

    :param gnomad_version: gnomAD version.
    :param samples: Optional list of samples to filter to.
    :param row_regions: Optional list of regions to filter to.
    :param variant_ht: Optional variant HT to use for filtering.
    :param ref: Whether to include reference genotypes. Default is False.
    :param pops: List of populations to filter to.
    :return: gnomAD regions MT.
    """
    if variant_ht is not None and row_regions is not None:
        raise DataException("Only one of variant_ht and row_regions can be provided.")

    if variant_ht is not None:
        # Get the intervals for the variants.
        row_regions = variant_ht.aggregate(
            hl.agg.collect(
                hl.interval(
                    start=variant_ht.locus, end=variant_ht.locus, includes_end=True
                )
            )
        )
        print("Number of input variants: ", len(row_regions))
        for region in variant_ht.regions.collect():
            row_regions.extend(region)

        print("Number of total intervals to filter to: ", len(row_regions))

    filter_ht = None
    loci_filter_ht = None
    if least_consequence is not None or max_af is not None:
        filter_ht = get_freq_ht(gnomad_version, row_regions)
        filter_ht = filter_freq_and_csq(
            filter_ht,
            get_vep_ht(gnomad_version, row_regions),
            filter_ht,
            max_af,
            least_consequence,
            variant_ht,
        )
        loci_filter_ht = filter_ht.key_by("locus")

    if gnomad_version in V3_VERSIONS:
        gnomad_mt = get_gnomad_v3_regions_mt(
            gnomad_version, samples_ht, row_regions, ref, pops, loci_filter_ht=loci_filter_ht
        )
    elif gnomad_version in V2_VERSIONS:
        gnomad_mt = get_gnomad_v2_regions_mt(
            gnomad_version, samples_ht, row_regions, ref, loci_filter_ht=loci_filter_ht
        )
    elif gnomad_version in V4_VERSIONS:
        gnomad_mt = get_gnomad_v4_regions_mt(
            gnomad_version, samples_ht, row_regions, ref, pops, loci_filter_ht=loci_filter_ht
        )
    else:
        raise DataException(
            f"Version {gnomad_version} is not supported for gnomAD regions HT"
        )

    if filter_ht is not None:
        gnomad_mt = gnomad_mt.semi_join_rows(filter_ht)

    return gnomad_mt


def get_variant_annotations_expr(
    table_key: hl.expr.StructExpression,
    gnomad_version: str,
    intervals: hl.tarray = None,
) -> Dict[str, hl.expr.Expression]:
    """
    Get the variant annotations expressions for the given gnomAD version.

    :param table_key: Table key.
    :param gnomad_version: gnomAD version.
    :param intervals: Intervals to filter to.
    :return: Dict of Variant annotations expressions.
    """
    # Load resources based on the gnomad versions.
    data_type, _ = gnomad_version_to_resource_version(gnomad_version)

    vep_ht = get_vep_ht(gnomad_version, intervals=intervals)
    freq_ht = get_freq_ht(gnomad_version, intervals=intervals)
    filter_ht = get_filter_ht(gnomad_version, intervals=intervals)
    info_ht = get_info_ht(gnomad_version, intervals=intervals)

    keyed_vep = vep_ht[table_key]
    keyed_freq = freq_ht[table_key]
    keyed_filtering = filter_ht[table_key]
    keyed_info = info_ht[table_key]

    freq_meta = freq_ht.freq_meta.collect()[0]

    annotations_expr = {"filters": keyed_filtering.filters}
    if gnomad_version in V3_VERSIONS:
        annotations_expr.update(
            {
                **{
                    f"info_{field}": keyed_filtering.info[field]
                    for field in keyed_filtering.info
                },
                **{
                    f"info_{field}": keyed_info.info[field]
                    for field in keyed_info.info
                    if field.startswith("AS") and field not in V3_EXCLUDE_INFO_FIELDS
                },
            }
        )
    if gnomad_version in V4_VERSIONS:
        filter_info_fields = V4_FILTERS_INFO_FIELDS[:]
        if data_type == "exomes":
            filter_info_fields += ["sibling_singleton"]
            # For more information on the selected compute_info_method, please see the
            # run_compute_info in:
            # gnomad_qc.v4.annotations.generate_variant_qc_annotations.py
            compute_info_method = hl.eval(
                filter_ht.filtering_model_specific_info.compute_info_method
            )
            info_struct = hl.struct(
                **keyed_info.site_info, **keyed_info[f"{compute_info_method}_info"]
            )
        else:
            # v3 info HT has no SOR or AS_SOR fields. They are computed by VQSR, so we
            # can grab them from the filters HT.
            info_struct = hl.struct(
                **keyed_info.info,
                SOR=keyed_filtering.SOR,
                AS_SOR=keyed_filtering.features.AS_SOR,
            )
        score_name = hl.eval(filter_ht.filtering_model.score_name)
        annotations_expr.update(
            {
                **{f"info_{score_name}": keyed_filtering[f"{score_name}"]},
                **{
                    f"info_{field}": keyed_filtering[field]
                    for field in filter_info_fields
                },
                **{
                    f"info_{field}": info_struct[field]
                    for field in info_struct
                    if field not in V3_EXCLUDE_INFO_FIELDS
                },
            }
        )
    elif gnomad_version in V2_VERSIONS:
        annotations_expr.update(
            {
                f"{'info_' if not f.startswith('info') else ''}{f}": keyed_info[f]
                for f in V2_INCLUDE_INFO_FIELDS
            }
        )
    else:
        raise DataException(
            f"Version {gnomad_version} is not supported for annotations"
        )

    # We only want to annotate overall and global pops frequencies.
    keys_to_keep = hl.set(["group", GENETIC_ANCESTRY_LABEL[gnomad_version]["gen_anc"]])

    # Annotate VEP and frequencies.
    annotations_expr.update(
        {
            "vep": keyed_vep.vep.select(*VEP_INCLUDE_FIELDS),
            "popmax": keyed_freq[GENETIC_ANCESTRY_LABEL[gnomad_version]["grpmax"]],
            # Combine frequency information with their metadata and filter down to
            # overall and global pops frequencies.
            "freq": hl.zip(freq_meta, keyed_freq.freq).filter(
                lambda x: x[0].keys().all(lambda k: keys_to_keep.contains(k))
            ),
        }
    )

    return annotations_expr


########################################################################################
# Functions that are specific to the gnomAD versions -- shouldn't need to be updated.
########################################################################################
def get_v2_exomes_v3_rel_ht() -> Union[hl.Table, hl.Table]:
    """
    Get the V2 exomes and V3 relatedness HTs.

    :return: V2 exomes and V3 relatedness HTs.
    """
    rel_ht = v3.sample_qc.v2_v3_relatedness.versions["3"].ht()
    rel_ht = rel_ht.key_by()
    rel_ht = rel_ht.select(
        relationship=get_relationship_expr(
            rel_ht.kin, rel_ht.ibd0, rel_ht.ibd1, rel_ht.ibd2
        ),
        s_v2_exomes=hl.if_else(
            rel_ht.i.data_type == "v2_exomes", rel_ht.i.s, rel_ht.j.s
        ),
        s_v3=hl.if_else(rel_ht.i.data_type == "v3_genomes", rel_ht.i.s, rel_ht.j.s),
    )
    rel_ht = rel_ht.filter(hl.or_else(rel_ht.relationship != UNRELATED, False))

    v3_rel_ht = rel_ht.group_by("s_v3").aggregate(
        v2_exomes_rel=hl.agg.collect(
            hl.struct(
                s=rel_ht.s_v2_exomes,
                rel=rel_ht.relationship,
            )
        )
    )

    v2_rel_ht = rel_ht.group_by("s_v2_exomes").aggregate(
        v3_rel=hl.agg.collect(
            hl.struct(
                s=rel_ht.s_v3,
                rel=rel_ht.relationship,
            )
        )
    )

    return v2_rel_ht, v3_rel_ht


def get_v4_rel_ht() -> hl.Table:
    """
    Get the V4 exomes and V4 genomes relatedness HT.

    :return: V4 exomes and V4 genomes relatedness HT.
    """
    rel_ht = v4.sample_qc.relatedness().ht()
    rel_ht = rel_ht.key_by()
    rel_ht = rel_ht.select(
        "relationship",
        s_v4_exomes=hl.if_else(rel_ht.i.data_type == "exomes", rel_ht.i.s, rel_ht.j.s),
        s_v4_genomes=hl.if_else(
            rel_ht.i.data_type == "genomes", rel_ht.i.s, rel_ht.j.s
        ),
    )
    rel_ht = rel_ht.filter(hl.or_else(rel_ht.relationship != UNRELATED, False))
    v4_exomes_rel_ht = rel_ht.group_by("s_v4_exomes").aggregate(
        v4_genomes_rel=hl.agg.collect(
            hl.struct(
                s=rel_ht.s_v4_genomes,
                rel=rel_ht.relationship,
            )
        )
    )
    v4_genomes_rel_ht = rel_ht.group_by("s_v4_genomes").aggregate(
        v4_exomes_rel=hl.agg.collect(
            hl.struct(
                s=rel_ht.s_v4_exomes,
                rel=rel_ht.relationship,
            )
        )
    )

    return v4_exomes_rel_ht, v4_genomes_rel_ht


def get_gnomad_v3_regions_mt(
    gnomad_version: str,
    samples_ht: Optional[hl.Table] = None,
    regions: Optional[Union[List[hl.Interval], hl.tarray]] = None,
    ref: bool = False,
    pops: List[str] = None,
    loci_filter_ht: Optional[hl.Table] = None,
) -> hl.MatrixTable:
    """
    Get the gnomAD v3 regions MT for the given version.

    :param gnomad_version: gnomAD version.
    :param samples: Optional list of samples to filter to.
    :param regions: Optional list of regions to filter to.
    :param ref: Whether to include reference genotypes. Default is False.
    :param pops: List of populations to filter to.
    :return: gnomAD v3 regions MT.
    """
    data_type, version = gnomad_version_to_resource_version(gnomad_version)

    meta = v3.meta.versions[version].ht()
    if gnomad_version == "v3.1":
        meta = meta.annotate(
            **meta.project_meta,
            **meta.subsets,
            **meta.sex_imputation,
            **meta.population_inference,
        )
        meta = meta.filter(~meta.sample_filters.hard_filtered)
    elif gnomad_version == "v3":
        meta = meta.filter(hl.len(meta.hard_filters) == 0)
        meta = meta.annotate(non_cancer=~meta.s.contains("TCGA"))

    if pops:
        logger.info(f"Filtering samples to {pops}")
        meta = meta.filter(hl.set(pops).contains(meta.pop))

    logger.info("Loading gnomAD MT")
    samples_keep_ht = meta.select()
    if samples_ht is not None:
        samples_keep_ht = samples_keep_ht.semi_join(samples_ht)

    gnomad_mt = get_gnomad_raw_data(
        gnomad_version,
        intervals=regions,
        samples_ht=samples_keep_ht,
        densify=ref,
        loci_filter_ht=loci_filter_ht,
    )

    logger.info("Adding sample metadata to MT")
    gnomad_mt = gnomad_mt.annotate_cols(meta=meta[gnomad_mt.col_key])

    return gnomad_mt


def get_gnomad_v4_regions_mt(
    gnomad_version: str,
    samples_ht: Optional[hl.Table] = None,
    regions: Optional[Union[List[hl.Interval], hl.tarray]] = None,
    ref: bool = False,
    pops: List[str] = None,
    loci_filter_ht: Optional[hl.Table] = None,
) -> hl.MatrixTable:
    """
    Get the gnomAD v3 regions MT for the given version.

    :param gnomad_version: gnomAD version.
    :param samples: List of samples to filter to.
    :param regions: List of regions to filter to.
    :param ref: Whether to include reference genotypes.
    :param pops: List of populations to filter to.
    :return: gnomAD v3 regions MT.
    """
    data_type, meta_version = gnomad_version_to_resource_version(
        gnomad_version, meta_version=True
    )

    meta = v4.meta(meta_version, data_type).ht()
    select_expr = {
        **meta.sex_imputation,
        **meta.population_inference,
    }
    if gnomad_version == "v4_genomes":
        select_expr = {
            **meta.project_meta,
            **select_expr,
            **meta.subsets,
        }
    elif gnomad_version == "v4_exomes":
        select_expr = {
            **meta.project_meta.drop("v2_meta", "ukb_meta"),
            **select_expr,
        }

    meta = meta.select("sample_filters", "high_quality", "release", **select_expr)
    meta = meta.filter(~meta.sample_filters.hard_filtered)

    if pops:
        logger.info(f"Filtering samples to {pops}")
        meta = meta.filter(hl.set(pops).contains(meta.pop))

    logger.info("Loading gnomAD MT")
    samples_keep_ht = meta.select()
    if samples_ht is not None:
        samples_keep_ht = samples_keep_ht.semi_join(samples_ht)

    gnomad_mt = get_gnomad_raw_data(
        gnomad_version,
        intervals=regions,
        samples_ht=samples_keep_ht,
        densify=ref,
        loci_filter_ht=loci_filter_ht,
    )

    logger.info("Adding sample metadata to MT")
    gnomad_mt = gnomad_mt.annotate_cols(meta=meta[gnomad_mt.col_key])

    return gnomad_mt


def get_gnomad_v2_regions_mt(
    gnomad_version: str,
    samples_ht: Optional[hl.Table] = None,
    regions: Optional[Union[List[hl.Interval], hl.tarray]] = None,
    ref: bool = False,
    loci_filter_ht: Optional[hl.Table] = None,
) -> hl.MatrixTable:
    """
    Get the gnomAD v2 regions MT for the given version.

    :param gnomad_version: gnomAD version.
    :param samples: List of samples to filter to.
    :param regions: List of regions to filter to.
    :param ref: Whether to include reference genotypes.
    :return: gnomAD v2 regions MT.
    """
    data_type, version = gnomad_version_to_resource_version(gnomad_version)

    logger.info("Loading and retrieving gnomAD v2 regions of interest")

    meta_ht = v2.basics.get_gnomad_meta(data_type, full_meta=True)
    meta_ht = meta_ht.annotate(non_cancer=~meta_ht.s.contains("TCGA"))

    logger.info(f"Filtering samples to samples that pass hard filters")
    meta_ht = meta_ht.filter((hl.len(meta_ht.hard_filters) == 0))

    logger.info("Loading gnomAD MT")
    samples_keep_ht = meta_ht.select()
    if samples_ht is not None:
        samples_keep_ht = samples_keep_ht.semi_join(samples_ht)

    gnomad_mt = get_gnomad_raw_data(
        gnomad_version,
        intervals=regions,
        samples_ht=samples_keep_ht,
        ref=ref,
        loci_filter_ht=loci_filter_ht,
    )

    logger.info("Adding sample metadata to MT")
    gnomad_mt = gnomad_mt.annotate_cols(meta=meta_ht[gnomad_mt.col_key])

    return gnomad_mt


########################################################################################
# Functions that shouldn't need to be updated with version additions.
########################################################################################
def split_interval_ht(ht: hl.Table, interval_field="regions") -> hl.Table:
    ht = ht.explode(ht[interval_field])
    interval_type = ht[interval_field].dtype
    intervals = ht.aggregate(
        hl.agg.collect(
            (ht[interval_field], hl.struct(locus=ht.locus, alleles=ht.alleles))
        )
    )

    rg = intervals[0][0].start.reference_genome
    contigs = rg.contigs
    sorted_intervals = sorted(
        intervals,
        key=lambda i: (
            contigs.index(i[0].start.contig),
            i[0].start.position,
            contigs.index(i[0].end.contig),
            i[0].end.position,
        ),
    )

    sorted_intervals = [[i, [s]] for i, s in sorted_intervals]
    split_intervals = deepcopy(sorted_intervals[:1])

    for interval, i in sorted_intervals:
        previous_interval = split_intervals[-1][0]
        previous_i = split_intervals[-1][1][:]
        if previous_interval.start.contig == interval.start.contig:
            if previous_interval.end.position < interval.end.position:
                if interval.start.position < previous_interval.end.position:
                    if interval.start.position != previous_interval.start.position:
                        split_intervals[-1][0] = hl.Interval(
                            previous_interval.start, interval.start
                        )
                        split_intervals.append(
                            [
                                hl.Interval(interval.start, previous_interval.end),
                                previous_i + i
                            ]
                        )
                    else:
                        split_intervals[-1][1].extend(i)

                    split_intervals.append(
                        [hl.Interval(previous_interval.end, interval.end), i[:]]
                    )
                else:
                    split_intervals.append(
                        [hl.Interval(interval.start, interval.end), i[:]])
            elif interval.start.position > previous_interval.end.position:
                if interval.start.position != previous_interval.start.position:
                    split_intervals[-1][0] = hl.Interval(
                        previous_interval.start, interval.start
                    )
                split_intervals.append(
                    [
                        hl.Interval(interval.start, interval.end),
                        previous_i + i
                    ]
                )
                split_intervals.append(
                    [hl.Interval(interval.end, previous_interval.end), previous_i]
                )
            else:
                if interval.start.position == previous_interval.start.position:
                    split_intervals[-1][1].extend(i)

        else:
            split_intervals.append([hl.Interval(interval.start, interval.end), i[:]])

    ht = hl.Table.parallelize(
        [hl.struct(region=r, variants=v) for r, v in split_intervals],
        schema=hl.tstruct(
            region=interval_type,
            variants=hl.tarray(
                hl.tstruct(locus=interval_type.point_type, alleles=hl.tarray(hl.tstr))
            ),
        ),
        key=["region"],
    )
    ht = ht.checkpoint(
        hl.utils.new_temp_file("split_intervals", "ht"), overwrite=True
    )

    return ht


def unify_sample_meta(ht: hl.Table, gnomad_version: str) -> hl.Table:
    """
    Unify the sample meta for the given gnomAD version.

    :param ht: Table to unify the sample meta for.
    :param gnomad_version: gnomAD version.
    :return: Table with unified sample meta.
    """
    _map = {
        field: SAMPLE_META_MAPPING[field][gnomad_version]
        for field in SAMPLE_META_MAPPING
        if gnomad_version in SAMPLE_META_MAPPING[field]
    }
    ht = ht.annotate(
        **{f"s_{field}": ht.meta[meta_field] for field, meta_field in _map.items()},
        meta=ht.meta.drop(*list(_map.values())),
    )

    return ht.rename({"meta": f"meta_{gnomad_version}"})


def export_to_tsv(ht: hl.Table, out_file: str, builds: List[str]) -> None:
    """
    Export the table to TSV.

    :param ht: Table to export.
    :param out_file: Output file.
    :param builds: Builds in export.
    :return: None.
    """

    def get_variant_flattened_expr(v: str) -> Dict[str, hl.expr.Expression]:
        flatten_expr = {}
        for build in builds:
            flatten_expr.update(
                {
                    f"{v}_{build}_chrom": ht[f"{v}_{build}_locus"].contig,
                    f"{v}_{build}_pos": ht[f"{v}_{build}_locus"].position,
                    f"{v}_{build}_ref": ht[f"{v}_{build}_alleles"][0],
                    f"{v}_{build}_alt": ht[f"{v}_{build}_alleles"][1],
                }
            )

        flatten_expr.update(
            {
                f"{v}_filters": hl.delimit(ht[f"{v}_filters"]),
                f"{v}_most_severe_consequence": ht[f"{v}_vep"].most_severe_consequence,
                f"{v}_motif_feature_consequences": hl.json(
                    ht[f"{v}_vep"].motif_feature_consequences
                ),
                f"{v}_regulatory_feature_consequences": hl.json(
                    ht[f"{v}_vep"].regulatory_feature_consequences
                ),
                f"{v}_transcript_consequences": hl.json(
                    ht[f"{v}_vep"].transcript_consequences
                ),
            }
        )

        flatten_expr.update(
            {
                f"{v}_{freq_field}_{freq_component}": ht[f"{v}_{freq_field}"][
                    freq_component
                ]
                for freq_component in ["AC", "AN", "AF", "homozygote_count"]
                for freq_field in ["popmax", "global_freq", "sample_pop_freq"]
            }
        )

        flatten_expr[f"{v}_filters"] = hl.delimit(ht[f"{v}_filters"])

        return flatten_expr

    # Flatten variant expressions.
    ht = ht.transmute(
        **get_variant_flattened_expr("v1"), **get_variant_flattened_expr("v2")
    )

    # Convert the samples meta to json and into a single column.
    meta_cols = [col for col in ht.row if col.startswith("meta_")]
    ht = ht.transmute(
        meta=hl.coalesce(
            *[hl.or_missing(hl.is_defined(ht[x]), hl.json(ht[x])) for x in meta_cols]
        )
    )

    # Drop some annotations that aren't useful.
    ht = ht.drop(*[f"{v}_a_index" for v in ["v1", "v2"] if f"{v}_a_index" in ht.row])

    # Order columns a little nicer.
    ht = ht.select(
        # First get all the columns for which order was specified.
        *[col for col in TSV_COLUMN_ORDER if col in ht.row],
        # Then get all the other ones.
        *[col for col in ht.row if col not in TSV_COLUMN_ORDER],
    )

    # Flatten whatever is left and export.
    ht.flatten().export(out_file)


def split_vqsr(vqsr_ht: hl.Table) -> hl.Table:
    """
    Split the VQSR HT.

    :param vqsr_ht: VQSR HT to split.
    :return: Split VQSR HT.
    """
    vqsr_ht = hl.split_multi(vqsr_ht)
    vqsr_ht = vqsr_ht.annotate(
        info=vqsr_ht.info.select(
            "NEGATIVE_TRAIN_SITE",
            "POSITIVE_TRAIN_SITE",
            SOR=vqsr_ht.info.AS_SOR[vqsr_ht.a_index - 1],
            VQSLOD=vqsr_ht.info.AS_VQSLOD[vqsr_ht.a_index - 1],
            culprit=vqsr_ht.info.AS_culprit[vqsr_ht.a_index - 1],
        )
    )

    return vqsr_ht


def gnomad_version_to_resource_version(
    gnomad_version,
    release_version=False,
    freq_version=False,
    vep_version=False,
    info_version=False,
    filter_version=False,
    meta_version=False,
):
    """
    Get the resource version for the given gnomAD version.

    :param gnomad_version: gnomAD version.
    :param release_version: Whether to get the release version.
    :param freq_version: Whether to get the frequency version.
    :param vep_version: Whether to get the VEP version.
    :param info_version: Whether to get the info version.
    :param filter_version: Whether to get the filter version.
    :param meta_version: Whether to get the meta version.
    :return: Resource version.
    """
    version_types = {
        "release_version": release_version,
        "freq_version": freq_version,
        "vep_version": vep_version,
        "info_version": info_version,
        "filter_version": filter_version,
        "meta_version": meta_version,
    }
    if sum(version_types.values()) > 1:
        raise ValueError(
            "Only one of release_version and freq_version can be set to True"
        )

    version = VERSION_RESOURCE_MAP.get(gnomad_version)
    data_type = version["data_type"]

    default_version = version.get("version", gnomad_version)
    version_type = [k for k, v in version_types.items() if v]

    if version_type:
        version = version.get(version_type[0], default_version)
    else:
        version = default_version

    return data_type, version


def get_liftover_variants_expr(
    locus: hl.expr.LocusExpression,
    alleles: hl.expr.ArrayExpression,
    destination_ref: hl.ReferenceGenome,
) -> Dict[str, hl.expr.Expression]:
    """
    Get the liftover variants expressions.

    :param locus: Locus expression.
    :param alleles: Alleles expression.
    :param destination_ref: Destination reference genome.
    :return: Liftover variants expressions.
    """
    logger.info("Running liftover for variant...")
    liftover_result = hl.liftover(locus, destination_ref, include_strand=True)
    lifted_over_locus = liftover_result.result
    lifted_over_alleles = alleles.map(
        lambda a: hl.if_else(
            liftover_result.is_negative_strand, hl.reverse_complement(a), a
        )
    )
    return {
        f"{locus.dtype.reference_genome.name}_locus": locus,
        f"{locus.dtype.reference_genome.name}_alleles": alleles,
        f"{destination_ref.name}_locus": lifted_over_locus,
        f"{destination_ref.name}_alleles": lifted_over_alleles,
    }


def add_liftover_annotations(var_ht: hl.Table, output_genome: str) -> hl.Table:
    """
    Adds lifted over variants in `lift_variants`.

    In addition, if the output genome is different from the var_ht genome,
    then regions are lifted over and `locus` and `alleles` are set to match the output
    genome.

    :param var_ht: Table with variants.
    :param output_genome: Output genome.
    :return: Table with lifted over variants.
    """
    _, destination_ref = get_liftover_genome(get_reference_genome(var_ht.locus))

    # Create lifted over variants. Make sure to assign.
    lift_variants = get_liftover_variants_expr(
        var_ht.locus, var_ht.alleles, destination_ref
    )
    lift_expr = dict(
        locus=lift_variants[f"{output_genome}_locus"],
        alleles=lift_variants[f"{output_genome}_alleles"],
        **lift_variants,
    )

    # If the destination isn't the same as the input, liftover the regions and assign
    # the lifted-over variants as primary locus, alleles.
    if destination_ref.name == output_genome:
        lift_expr["regions"] = var_ht.regions.map(
            lambda x: hl.if_else(
                (x.end.position - x.start.position) > 1,
                hl.liftover(x, destination_ref),
                hl.bind(
                    lambda l_x: hl.interval(
                        l_x, l_x, includes_start=True, includes_end=False
                    ),
                    hl.liftover(x.start, destination_ref),
                ),
            )
        )

    return var_ht.annotate(**lift_expr)


# TODO: This was copied from utils as the utils version needs some fixing.
def get_liftover_genome(source: hl.ReferenceGenome) -> List[hl.ReferenceGenome]:
    """
    Get the liftover genome for the given source.

    :param source: Source reference genome.
    :return: List of source build (with liftover chain added) and destination build
        (with sequence loaded).
    """

    logger.info(
        "Loading reference genomes, adding chain file, and loading fasta sequence for"
        " destination build"
    )
    if source.name == "GRCh38":
        target = hl.get_reference("GRCh37")
        chain = "gs://hail-common/references/grch38_to_grch37.over.chain.gz"
        if not target.has_sequence():
            target.add_sequence(
                "gs://hail-common/references/human_g1k_v37.fasta.gz",
                "gs://hail-common/references/human_g1k_v37.fasta.fai",
            )
    else:
        target = hl.get_reference("GRCh38")
        chain = "gs://hail-common/references/grch37_to_grch38.over.chain.gz"
        if not target.has_sequence():
            target.add_sequence(
                "gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz",
                "gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai",
            )
    if not source.has_liftover(target):
        source.add_liftover(chain, target)

    return [source, target]


def get_variants_samples(version: str, var_ht: hl.Table) -> hl.Table:
    """
    Get Table of samples carrying variants in a variant Table for specified gnomAD version.

    :param version: gnomAD version.
    :param var_ht: Table of variants to find carriers of.
    :return: Table of variant carriers.
    """
    # Convert variants to interval for pushdown.
    logger.info("Retrieving samples that have the variant...")
    intervals = var_ht.aggregate(
        hl.agg.collect(
            hl.interval(start=var_ht.locus, end=var_ht.locus, includes_end=True)
        )
    )
    mt = get_gnomad_raw_data(version, intervals)

    # Now filter gnomAD to keep only the correct alleles of interest after splitting.
    mt = mt.semi_join_rows(var_ht.key_by("locus", "alleles"))

    # Create a table of variants -> sample(s).
    # This Table contains one entry per variant and sample.
    # Remove all row/columns annotations -- we don't need them in this table.
    mt = mt.select_rows().select_cols().select_entries(*ENTRY_FIELDS_TO_KEEP)

    # Only keep non-reference genotypes.
    return mt.filter_entries(mt.GT.is_non_ref()).entries()


def get_samples_ht_from_gnomad_mt(
    gnomad_mt: hl.MatrixTable,
    var_ht: hl.Table
) -> hl.Table:
    """
    Get a Table of sample non ref genotypes per variant from a gnomAD MT.

    :param gnomad_mt: gnomAD MT.
    :param var_ht: Variant HT.
    :return: Table of sample non ref genotypes per variant.
    """
    gnomad_mt = gnomad_mt.select_cols()

    variant_mt = gnomad_mt.semi_join_rows(var_ht.key_by("locus", "alleles"))
    sample_mt = variant_mt.filter_entries(variant_mt.GT.is_non_ref())
    sample_ht = sample_mt.annotate_rows(
        samples=hl.agg.collect(hl.struct(**sample_mt.col, **sample_mt.entry))
    ).rows()

    return sample_ht


def get_regions_ht_for_one_variant(
    var_ht: hl.Table,
    variant_samples_ht: hl.Table,
    gnomad_mt: hl.MatrixTable,
    row: hl.struct,
    gnomad_version: str,
):
    """
    Get the regions HT for one variant.

    :param var_ht: Variant HT.
    :param variant_samples_ht: Variant samples HT.
    :param gnomad_mt: gnomAD MT.
    :param row: Row to get regions HT for.
    :param gnomad_version: gnomAD version.
    :return: Regions HT for one variant.
    """
    # Prefix all variants and genotype fields with 'v2'.
    gnomad_mt = gnomad_mt.rename({x: f"v2_{x}" for x in gnomad_mt.row})
    gnomad_mt = gnomad_mt.rename({x: f"v2_{x}" for x in gnomad_mt.entry})

    # Add the target variant information. Prefix the fields with "v1_".
    gnomad_mt = gnomad_mt.annotate_rows(
        **{
            f"v1_{k}": hl.literal(row[k], dtype=var_ht[k].dtype)
            for k in row
            if k != "regions"
        }
    )

    # Add the variants annotations. Prefix the fields with "v2_".
    intervals = gnomad_mt.aggregate_rows(
        # NOTE: Added this here in hopes to get past stuck job.
        hl.agg.collect(
            hl.interval(
                start=gnomad_mt.v2_locus, end=gnomad_mt.v2_locus, includes_end=True
            )
        )
    )

    gnomad_mt = gnomad_mt.transmute_rows(
        **{
            f"v2_{k}": v
            for k, v in get_variant_annotations_expr(
                gnomad_mt.row_key,
                gnomad_version,
                # NOTE: Added intervals here in hopes to get past stuck job.
                intervals,
            ).items()
        }
    )

    # Create a table with all variants / samples in the region.
    region_ht = gnomad_mt.entries()

    # Key table by locus, alleles and sample.
    region_ht = region_ht.key_by("v1_locus", "v1_alleles", "s")

    # Add all genotype information about the target variant for this sample. Prefix the
    # fields with "v1_".
    keyed_variant_samples = variant_samples_ht[region_ht.key]
    region_ht = region_ht.annotate(
        **{f"v1_{f}": keyed_variant_samples[f] for f in keyed_variant_samples}
    )

    return region_ht


def get_regions_ht_for_all_variants(
    var_ht: hl.Table,
    variant_samples_ht: hl.Table,
    gnomad_mt: hl.MatrixTable,
    gnomad_version: str,
    gnomad_meta_ht: hl.Table,
):
    """
    Get the regions HT for all variants.

    :param var_ht: Variant HT.
    :param variant_samples_ht: Variant samples HT.
    :param gnomad_mt: gnomAD MT.
    :param gnomad_version: gnomAD version.
    :param gnomad_meta_ht: gnomAD meta HT.
    :return: Regions HT for all variants.
    """
    gnomad_mt = gnomad_mt.rename({x: f"v2_{x}" for x in gnomad_mt.row})

    variant_ht = split_interval_ht(var_ht)
    variant_keyed = variant_ht[gnomad_mt.v2_locus]
    region_mt = gnomad_mt.annotate_rows(v1_variants=variant_keyed.variants)
    region_mt = region_mt.filter_entries(region_mt.GT.is_non_ref()).select_cols()
    region_mt = region_mt.explode_rows(region_mt.v1_variants)
    region_mt = region_mt.annotate_rows(
        **{f"v1_{k}": region_mt.v1_variants[k] for k in region_mt.v1_variants}
    ).drop("v1_variants")

    region_keyed = variant_samples_ht[region_mt.v1_locus, region_mt.v1_alleles]
    region_ht = region_mt.annotate_rows(
        **{f"v1_{k}": region_keyed[k] for k in region_keyed if k != "samples"},
        v1_entries=region_keyed.samples,
        v2_entries=hl.agg.collect(
            hl.struct(**region_mt.col, **region_mt.entry)
        ),
    ).rows()
    region_ht = region_ht.checkpoint(hl.utils.new_temp_file("int_region_ht", "ht"))
    region_ht = region_ht.annotate(
        samples=hl.set(region_ht.v1_entries.map(lambda x: x.s)).intersection(
            hl.set(region_ht.v2_entries.map(lambda x: x.s))
        )
    )
    region_ht = region_ht.filter(hl.len(region_ht.samples) > 0)
    region_ht = region_ht.annotate(
        v1_entries=region_ht.v1_entries.filter(
            lambda x: region_ht.samples.contains(x.s)
        ),
        v2_entries=region_ht.v2_entries.filter(
            lambda x: region_ht.samples.contains(x.s)
        )
    ).drop("samples")
    region_ht = region_ht.annotate(
        **{
            f"v1_{k}": v
            for k, v in get_variant_annotations_expr(
                hl.struct(v1_locus=region_ht.v1_locus, v1_alleles=region_ht.v1_alleles),
                gnomad_version,
            ).items()
        },
        **{
            f"v2_{k}": v
            for k, v in get_variant_annotations_expr(
                region_ht.key,
                gnomad_version,
            ).items()
        }
    )
    region_ht = region_ht.checkpoint(hl.utils.new_temp_file("region_ht", "ht"))
    region_ht = region_ht.explode(region_ht.v1_entries)
    region_ht = region_ht.transmute(
        v1_sample=region_ht.v1_entries,
        v2_sample=region_ht.v2_entries.filter(
            lambda x: x.s == region_ht.v1_entries.s
        )[0],
        **gnomad_meta_ht[region_ht.v1_entries.s],
    )
    # Add 'v1' and 'v2' prefixes.
    # Add the variants annotations. Prefix the fields with "v2_".
    region_ht = region_ht.transmute(
        **region_ht.v1_sample.rename(
            {x: f"v1_{x}" for x in region_ht.v1_sample}
        ),
        **region_ht.v2_sample.rename(
            {x: f"v2_{x}" for x in region_ht.v2_sample}
        ),
    )

    return region_ht


def create_regions_ht(
    var_ht: hl.Table,
    gnomad_mt: hl.MatrixTable,
    gnomad_version: str,
    variant_samples_ht: hl.Table = None,
    row: hl.struct = None,
    af_max: float = None,
    needs_liftover: bool = False,
):
    """
    Create a regions HT for the given variant HT.

    :param var_ht: Variant HT.
    :param gnomad_mt: gnomAD MT.
    :param gnomad_version: gnomAD version.
    :param variant_samples_ht: Variant samples HT.
    :param row: Row to annotate.
    :param af_max: Maximum allele frequency.
    :param needs_liftover: Whether the variants need liftover.
    :return: Regions HT.
    """
    gnomad_meta_ht = gnomad_mt.cols().checkpoint(
        hl.utils.new_temp_file("gnomad_meta_ht", "ht")
    )
    # Select only entry fields to keep.
    gnomad_mt = gnomad_mt.select_entries(*ENTRY_FIELDS_TO_KEEP)

    if needs_liftover:
        # Add the lifted-over variants for v2 and annotate the lifted-over variants.
        _, destination_ref = get_liftover_genome(get_reference_genome(gnomad_mt.locus))
        gnomad_mt = gnomad_mt.annotate_rows(
            **get_liftover_variants_expr(
                gnomad_mt.locus, gnomad_mt.alleles, destination_ref
            )
        )
    else:
        gnomad_mt = gnomad_mt.annotate_rows(
            **{
                f"{gnomad_mt.locus.dtype.reference_genome.name}_locus": gnomad_mt.locus,
                f"{gnomad_mt.locus.dtype.reference_genome.name}_alleles": gnomad_mt.alleles,
            }
        )

    if row is not None:
        region_ht = get_regions_ht_for_one_variant(
            var_ht=var_ht,
            variant_samples_ht=variant_samples_ht,
            gnomad_mt=gnomad_mt,
            row=row,
            gnomad_version=gnomad_version,
        )
    else:
        if variant_samples_ht is None:
            variant_samples_ht = get_samples_ht_from_gnomad_mt(
                gnomad_mt, var_ht
            ).checkpoint(hl.utils.new_temp_file("sample_ht", "ht"))

        region_ht = get_regions_ht_for_all_variants(
            var_ht=var_ht,
            variant_samples_ht=variant_samples_ht,
            gnomad_mt=gnomad_mt,
            gnomad_version=gnomad_version,
            gnomad_meta_ht=gnomad_meta_ht,
        )

    # Add the gnomad version.
    region_ht = region_ht.annotate(gnomad_version=gnomad_version)

    # Remove the generic locus / alleles (build-specific ones are left).
    region_ht = region_ht.key_by().drop(
        "v1_locus", "v1_alleles", "v2_locus", "v2_alleles"
    )

    # Only keep freq information for the overall pop + sample global pop.
    for v in ["v1", "v2"]:
        region_ht = region_ht.annotate(
            **{
                f"{v}_global_freq": region_ht[f"{v}_freq"][0][1],
                f"{v}_sample_pop_freq": region_ht[f"{v}_freq"].find(
                    lambda x: hl.coalesce(x[0].get("pop"), x[0].get("gen_anc"))
                    == region_ht.meta.pop
                )[1],
            }
        )
        region_ht = region_ht.drop(f"{v}_freq")

    if af_max:
        logger.info(
            "Filtering v2 variants to variants with a global AF less than or equal to"
            f" {af_max}"
        )
        region_ht = region_ht.filter(region_ht.v2_global_freq.AF <= af_max)

    # Unify samples meta.
    logger.info("Unifying sample meta...")
    region_ht = unify_sample_meta(region_ht, gnomad_version)

    return region_ht


def import_variants_regions_ht(
    input_tsv: str,
    input_genome: str,
    gnomad_version: str,
    use_genes_as_regions: bool = False,
) -> hl.Table:
    """
    This function imports the variants and regions from a TSV file.

    :param input_tsv: Input TSV file.
    :param input_genome: Input genome.
    :param gnomad_version: gnomAD version.
    :param use_genes_as_regions: Whether to use genes as regions.
    :return: Variants and regions HT.
    """
    logger.info("Importing variant table...")
    var_ht = hl.import_table(input_tsv, no_header=True, min_partitions=12)

    # Make sure that the dimensions are what is expected: chrom, pos, ref, alt, regions.
    if len(var_ht.row) != 5:
        raise Exception(
            f"Unexpected number of columns in input TSV file: {len(var_ht.row)}. There"
            " should be 5 columns: chrom, pos, ref, alt, region[s]."
        )
    # Return a Table with the standard hail types.
    var_ht = var_ht.transmute(
        locus=hl.locus(
            contig=var_ht.f0, pos=hl.int(var_ht.f1), reference_genome=input_genome
        ),
        # Alleles in hail are represented as an array of string, with the first one
        # being the ref.
        alleles=[var_ht.f2, var_ht.f3],
        regions=hl.if_else(
            var_ht.f4 == "",
            hl.missing(
                hl.tarray(
                    hl.expr.types.tinterval(
                        hl.expr.types.tlocus(reference_genome=input_genome)
                    )
                )
            ),
            # Here the split means that it will support a comma-delimited list of
            # regions. It also means that the variable regions will be an array.
            var_ht.f4.split(",").map(
                lambda x: hl.parse_locus_interval(x, reference_genome=input_genome)
            ),
        ),
    )

    if use_genes_as_regions:
        intervals = var_ht.aggregate(
            hl.agg.collect(
                hl.interval(start=var_ht.locus, end=var_ht.locus, includes_end=True)
            )
        )

        vep_ht = get_vep_ht(gnomad_version, intervals)
        vep_ht = vep_ht.select(
            gene_ids=hl.set(vep_ht.vep.transcript_consequences.map(lambda x: x.gene_id))
        )
        var_ht = var_ht.annotate(gene_id=vep_ht[var_ht.locus, var_ht.alleles].gene_ids)
        var_ht = var_ht.explode(var_ht.gene_id)
        gene_ids = var_ht.aggregate(hl.agg.collect_as_set(var_ht.gene_id))
        gencode_ht = hl.experimental.import_gtf(
            GENE_REGION_GTF_MAP[gnomad_version],
            reference_genome=input_genome,
            skip_invalid_contigs=True,
            min_partitions=500,
            force_bgz=True,
        )
        gencode_ht = gencode_ht.annotate(
            gene_id=gencode_ht.gene_id.split('\\.')[0],
            transcript_id=gencode_ht.transcript_id.split('\\.')[0]
        )
        gencode_ht = gencode_ht.filter(
            hl.any(
                lambda y: (gencode_ht.feature == 'gene') & (gencode_ht.gene_id == y.split('\\.')[0]),
                gene_ids,
            )
        )
        gencode_ht = gencode_ht.group_by(gencode_ht.gene_id).aggregate(
            intervals=hl.agg.collect(gencode_ht.interval),
            gene_symbols=hl.agg.collect_as_set(gencode_ht.gene_name),
        )
        var_ht = var_ht.annotate(
            regions=gencode_ht[var_ht.gene_id].intervals,
            gene_symbols=gencode_ht[var_ht.gene_id].gene_symbols,
        )
        var_ht = var_ht.group_by(var_ht.locus, var_ht.alleles).aggregate(
            regions=hl.agg.explode(lambda x: hl.agg.collect_as_set(x), var_ht.regions),
            gene_symbols=hl.agg.explode(lambda x: hl.agg.collect_as_set(x), var_ht.gene_symbols),
            gene_ids=hl.agg.collect_as_set(var_ht.gene_id),
        )

    # Add lifted over variants and set locus, alleles and regions to match the gnomAD
    # build.
    output_genome = VERSION_GENOME_BUILD[gnomad_version]
    if input_genome != output_genome:
        logger.warning(
            f"Variants in the input TSV will be lifter over from {input_genome} to"
            f" {output_genome} in order to use gnomAD {gnomad_version}."
        )
        var_ht = add_liftover_annotations(var_ht, output_genome)
    else:
        var_ht = var_ht.annotate(
            **{
                f"{output_genome}_locus": var_ht.locus,
                f"{output_genome}_alleles": var_ht.alleles,
            }
        )

    return var_ht.key_by("locus", "alleles")


def vep_genes_expr(
    vep_expr: hl.expr.StructExpression,
    least_consequence: str
) -> hl.expr.SetExpression:
    vep_consequences = hl.literal(
        set(CSQ_ORDER[0:CSQ_ORDER.index(least_consequence) + 1])
    )
    return (
                hl.set(
                    vep_expr.transcript_consequences
                        .filter(
                        lambda tc: (tc.biotype == 'protein_coding') &
                                   (tc.consequence_terms.any(lambda c: vep_consequences.contains(c)))
                    )
                        .map(lambda x: x.gene_id)
                )
            )


def filter_freq_and_csq(
    t: Union[hl.Table, hl.MatrixTable],
    vep_ht: hl.Table = None,
    freq_ht: hl.Table = None,
    max_freq: float = None,
    least_consequence: str = None,
    variant_ht: hl.Table = None,
) -> hl.MatrixTable:
    """
    Filters MatrixTable to include variants that:
    1. Have a global AF <= `max_freq`
    2. Have a consequence at least as severe as `least_consequence` (based on ordering from CSQ_ORDER)

    :param t: Input HT/MT
    :param max_freq: Max. AF to keep
    :param least_consequence: Least consequence to keep.
    :return: Filtered MT
    """
    vep_filter = vep_ht is not None and least_consequence is not None
    freq_filter = freq_ht is not None and max_freq is not None
    if not vep_filter and not freq_filter:
        logger.info(
            "No VEP HT and least_consequence or Freq HT and max_freq were provided, "
            "so no filtering will be done."
        )
        return t

    filter_expr = True

    if isinstance(t, hl.MatrixTable):
        t_key = t.row_key
    else:
        t_key = t.key

    if vep_filter:
        filter_expr &= hl.len(vep_genes_expr(vep_ht[t_key].vep, least_consequence)) > 0

    if freq_filter:
        af_expr = freq_ht[t_key].freq[0].AF
        filter_expr &= hl.is_missing(af_expr) | (af_expr <= max_freq)

    if variant_ht is not None:
        filter_expr |= hl.is_defined(variant_ht[t_key])

    if isinstance(t, hl.MatrixTable):
        return t.filter_rows(filter_expr)
    else:
        return t.filter(filter_expr)


def main(args):
    hl.init(
        log="sanna_gnomad_query.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day/sanna_gnomad_query",
    )
    result_tables = []
    if "v2_exomes" in args.gnomad or "v3" in args.gnomad:
        v2_rel_ht, v3_rel_ht = get_v2_exomes_v3_rel_ht()
    if "v4_exomes" in args.gnomad or "v4_genomes" in args.gnomad:
        v4_exomes_rel_ht, v4_genomes_rel_ht = get_v4_rel_ht()

    all_builds = [VERSION_GENOME_BUILD[b] for b in args.gnomad] + [args.input_genome]
    needs_liftover = len(set(all_builds)) > 1

    # Loop through all the desired gnomAD versions.
    for gnomad_version in args.gnomad:
        logger.info(f"Looking up variants in gnomaD {gnomad_version}")
        data_type, _ = gnomad_version_to_resource_version(gnomad_version)

        # First read in the variants and regions provided (lifting over if necessary).
        var_ht = import_variants_regions_ht(
            input_tsv=args.input_tsv,
            input_genome=args.input_genome,
            gnomad_version=gnomad_version,
            use_genes_as_regions=args.use_genes_as_regions,
        )
        var_ht = var_ht.repartition(min(var_ht.count(), 5000)).checkpoint(
            hl.utils.new_temp_file("input_var_ht", "ht")
        )
        var_ht.show(200)

        variant_samples_ht = None
        row = None

        # Because we store refs only in the hardcalls for v2 (which doesn't have all
        # genotype quality info). We'll get the non-ref genotypes from the non-ref
        # version of gnomAD and the ref genotypes from the hardcalls. The `gnomad_refs`
        # array allows us to loop through these as necessary.
        if args.ref and gnomad_version == "v2_genomes":
            logger.warning(
                "Reference genotypes for v2 genomes are in cold storage and"
                " will not appear in the results"
            )
            gnomad_refs = [False]
        elif args.ref and gnomad_version == "v2_exomes":
            gnomad_refs = [False, True]
        else:
            gnomad_refs = [args.ref]

        for gnomad_ref in gnomad_refs:
            # Get the variants in the region of interest.
            gnomad_mt = get_gnomad_regions_mt(
                gnomad_version,
                variant_ht=var_ht,
                ref=gnomad_ref,
                pops=args.pops,
                least_consequence=args.least_consequence,
                max_af=args.af_max,
            ).checkpoint(hl.utils.new_temp_file("gnomad_mt", "mt"))

            if gnomad_version == "v3":
                # Add ID / relationship to sample in v2 if any.
                gnomad_mt = gnomad_mt.annotate_cols(
                    v2_exomes_rel=v3_rel_ht[gnomad_mt.s].v2_exomes_rel,
                )
            # If exomes: add ID / relationship to sample in v3 if any.
            elif gnomad_version == "v2_exomes":
                gnomad_mt = gnomad_mt.annotate_cols(
                    v3_rel=v2_rel_ht[gnomad_mt.s].v3_rel
                )
            elif gnomad_version == "v4_exomes":
                gnomad_mt = gnomad_mt.annotate_cols(
                    v4_genomes_rel=v4_exomes_rel_ht[gnomad_mt.s].v4_genomes_rel
                )
            elif gnomad_version == "v4_genomes":
                gnomad_mt = gnomad_mt.annotate_cols(
                    v4_exomes_rel=v4_genomes_rel_ht[gnomad_mt.s].v4_exomes_rel
                )

            region_ht = create_regions_ht(
                var_ht=var_ht,
                gnomad_mt=gnomad_mt,
                gnomad_version=gnomad_version,
                variant_samples_ht=variant_samples_ht,
                row=row,
                af_max=args.af_max,
                needs_liftover=needs_liftover,
            )
            region_ht = region_ht.annotate(type=data_type)

            logger.info(
                    "Checkpointing variant table and appending to final result table"
                    " list"
            )
            result_tables.append(
                region_ht.checkpoint(
                    f"{args.checkpoint_dir}/part_{len(result_tables)}.ht",
                    overwrite=True,
                )
            )

    # Concatenate all the variants / regions /gnomAD builds and export.
    if len(result_tables) < 1:
        logger.warning("No results returned ?!?")
    else:
        # Get the first result table.
        final_result_ht = result_tables.pop()
        if len(result_tables) > 0:
            final_result_ht = final_result_ht.union(*result_tables, unify=True)
        if args.out_ht:
            final_result_ht = final_result_ht.checkpoint(args.out_ht, overwrite=True)
        if args.out_tsv:
            export_to_tsv(final_result_ht, args.out_tsv, list(all_builds))


# This part of python handling of script arguments.
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_tsv",
        help=(
            "Tab-separated input file containing variants and regions. If multiple"
            " regions correspond to a single variant, they should be entered as a"
            " comma-separated list."
        ),
    )
    parser.add_argument(
        "--input_genome",
        choices=["GRCh37", "GRCh38"],
        help=(
            "Build for coordinates/ref/alt alleles in the input file. IMPORTANT: when"
            ' using GRCh37, contigs should NOT start with "chr" or "chrom". When using'
            ' GRCh38, they NEED to start with "chr"'
        ),
    )
    parser.add_argument(
        "--use_genes_as_regions",
        help=(
            "If specified, the regions examined for each variant will be the gene "
            "interval(s) associated with the variant."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--gnomad",
        choices=VERSIONS,
        nargs="+",
        help="Version(s) of gnomAD to search.",
    )
    parser.add_argument(
        "--ref",
        help="If added, reference genotypes are also exported for all gnomAD variants.",
        action="store_true",
    )
    parser.add_argument(
        "--checkpoint_dir",
        help=(
            "Base directory for checkpointing MTs. Default is gs://gnomad-tmp/sanna,"
            " which is cleaned-up every 7 days."
        ),
        default="gs://gnomad-tmp/sanna",
    )

    parser.add_argument(
        "--out_ht",
        help="Output Hail Table. If specified, results are also saved as a hail table.",
    )
    parser.add_argument(
        "--out_tsv",
        help=(
            "If specified, results are output as TSV. Depending on the extension"
            " provided, the result may be uncompressed (.tsv) or compressed (.tsv.gz)"
        ),
    )
    parser.add_argument(
        "--af_max",
        default=None,
        type=float,
        help=(
            "If specified, results are filtered to variants with a global AF less than"
            " or equal to the passed value."
        ),
    )
    parser.add_argument(
        '--least_consequence',
        help=f'Least consequence for the second variant to be included in the output.',
        default=None,
        choices=CSQ_ORDER,
    )
    parser.add_argument(
        "--pops",
        help=(
            "If specified, results are filtered to passed continental pops. Default is"
            " an unfiltered dataset"
        ),
        nargs="+",
        choices=["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "mid", "sas", "oth"],
        default=[],
    )
    # Make sure we're outputting something.
    args = parser.parse_args()
    if not args.out_ht and not args.out_tsv:
        logger.error("At least one of --out_ht or --out_tsv needs to be specified.")
    else:
        main(args)
