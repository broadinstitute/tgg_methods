"""Utility functions to compare sample pairs that share a rare doubleton to related samples."""
import logging
from typing import List, Set, Tuple

import hail as hl

from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.filtering import filter_to_adj
from gnomad.utils.vcf import SPARSE_ENTRIES

from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import interval_qc_path, relatedness_ht_path
from ukbb_qc.resources.variant_qc import info_ht_path, NA12878, SYNDIP

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("doubleton_utils")
logger.setLevel(logging.INFO)


TEMP_PATH = "gs://gnomad-tmp/kc"
"""
Path to bucket for temporary data.
"""

VDS_PATH = "gs://ukbb-pharma-exome-analysis/500k_temp/500k.vds"
"""
Path to Variant Dataset (VDS) that contains high quality samples from the final UK Biobank callset.

Generated using the following commands:

vds = hl.vds.read_vds("gs://gnomad/raw/exomes/4.0/gnomad_v4.0.vds")
meta_ht = hl.read_table("gs://broad-ukbb/broad.freeze_7/sample_qc/meta.ht")
meta_ht = meta_ht.filter(meta_ht.sample_filters.high_quality)
vds = hl.vds.filter_samples(vds, meta_ht, remove_dead_alleles=True)
call_stats_ht = hl.read_table(ukb.var_annotations_ht_path('ukb_freq', *TRANCHE_DATA[CURRENT_TRANCHE]))
freq_index = get_cohort_index(call_stats_ht)
var = vds.variant_data.annotate_rows(call_stats=call_stats_ht[vds.variant_data.row_key].freq[freq_index])
vds = hl.vds.variant_dataset.VariantDataset(vds.reference_data, var)
vds.write(ukb_exomes_path, overwrite=args.overwrite)

Full VDS size: 8.69 TiB
VDS variant data size: 561.97 GiB
"""

TRANCHE_DATA = ("broad", CURRENT_FREEZE)
"""
UKB tranche data (data source and data freeze number)
"""

# Remove END from SPARSE_ENTRIES (not present in VDS variant data entries)
SPARSE_ENTRIES.remove("END")


def get_doubleton_sites(
    vds_path: str = VDS_PATH,
    temp_path: str = TEMP_PATH,
    tranche_data: Tuple[str, int] = TRANCHE_DATA,
    sparse_entries: List[str] = SPARSE_ENTRIES,
) -> hl.Table:
    """
    Filter UKB VDS to bi-allelic, autosomal sites in interval QC pass regions with an allele count of two.

    Produce Table of bi-allelic, autosomal SNPs in interval QC pass regions that have an allele count of two
    with no homozygotes (calculated on adj genotypes).

    :param str vds_path: Path to UKB 455k VDS. Default is VDS_PATH.
    :param str temp_path: Path to bucket to store Table and other temporary data. Default is TEMP_PATH.
    :param Tuple[str, int] tranche_data: UKB tranche data (data source and data freeze number). Default is TRANCHE_DATA.
    :param List[str] sparse_entries: List of fields to select from VDS. Default is SPARSE_ENTRIES.
    :return: Table of high quality sites with doubletons.
    """
    logger.info("Reading in VDS and filtering to bi-allelic SNPs...")
    mt = hl.vds.read_vds(vds_path).variant_data
    # Drop unnecessary annotations
    mt = mt.select_rows().select_entries(*sparse_entries)
    mt = mt.filter_rows(bi_allelic_expr(mt) & hl.is_snp(mt.alleles[0], mt.alleles[1]))

    logger.info("Filter to autosomes and splitting multiallelics...")
    mt = mt.filter_rows(mt.locus.in_autosome())
    mt = hl.experimental.sparse_split_multi(mt)

    logger.info("Removing AS_lowqual sites...")
    info_ht = hl.read_table(info_ht_path(*tranche_data, split=True))
    mt = mt.filter_rows(~info_ht[mt.row_key].AS_lowqual)

    logger.info("Filtering to interval QC pass regions...")
    interval_ht = hl.read_table(interval_qc_path(*tranche_data, "autosomes"))
    mt = mt.filter_rows(hl.is_defined(interval_ht[mt.locus]))

    logger.info("Filtering to adj and calculating allele count...")
    mt = filter_to_adj(mt)
    mt = mt.annotate_rows(call_stats=hl.agg.call_stats(mt.GT, mt.alleles))
    # Get AC at allele index 1 (call_stats includes a count for each allele, including reference)
    mt = mt.transmute_rows(
        ac=mt.call_stats.AC[1], n_hom=mt.call_stats.homozygote_count[1]
    )

    logger.info("Filtering to an allele count of two and returning...")
    ht = mt.rows()
    ht = ht.filter((ht.ac == 2) & (ht.n_hom == 0))
    ht = ht.checkpoint(f"{temp_path}/high_quality_sites.ht", overwrite=True)
    return ht


def get_doubleton_samples(
    vds_path: str = VDS_PATH,
    temp_path: str = TEMP_PATH,
    control_samples: Set[str, str] = {NA12878, SYNDIP},
):
    """
    Get IDs of sample that share a rare doubleton in the 455k VDS.

    Filter VDS variant data to sites present in specified input Table, collect sample IDs, annotate IDs onto rows, and write
    HT to temporary path.

    :param str vds_path: Path to UKB 455k VDS. Default is VDS_PATH.
    :param str temp_path: Path to bucket to store Table and other temporary data. Default is TEMP_PATH.
    :param Tuple[str, str] control_samples: Tuple of control sample IDs to remove. Default is (NA12878, SYNDIP).
    :return: Table keyed by sample IDs that share a rare doubleton.
    """
    logger.info("Getting IDs of samples that share a rare doubleton...")
    mt = hl.vds.read_vds(vds_path).variant_data
    mt = hl.experimental.sparse_split_multi(mt)
    mt = mt.filter_cols(~hl.literal(control_samples).contains(mt.s))
    ht = get_doubleton_sites(vds_path)
    mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
    mt = mt.annotate_rows(pair=hl.agg.filter(mt.GT.is_het(), hl.agg.collect(mt.s)))
    ht = mt.annotate_rows(s1=mt.pair[0], s2=mt.pair[1]).rows()
    ht = ht.checkpoint(f"{temp_path}/doubletons.ht", overwrite=True)

    logger.info("Filtering to unique sample pairs and returning...")
    ht = ht.key_by("s1", "s2").select().distinct()
    ht = ht.checkpoint(f"{temp_path}/doubletons_uniq.ht", overwrite=True)
    logger.info("Found %i unique sample pairs", ht.count())
    return ht


def compare_doubletons_to_related(
    tranche_data: Tuple[str, int] = TRANCHE_DATA, temp_path: str = TEMP_PATH
) -> None:
    """
    Get sample pairs that share doubletons and compare these pairs to samples in 455k relatedness Table.

    :param Tuple[str, int] tranche_data: UKB tranche data (data source and data freeze number). Default is TRANCHE_DATA.
    :param str temp_path: Path to bucket to store Table and other temporary data. Default is TEMP_PATH.
    :return: None; function prints information to stdout.
    """
    ht = get_doubleton_samples()
    rel_ht = hl.read_table(relatedness_ht_path(*tranche_data))
    rel_ht = rel_ht.key_by(i=rel_ht.i.s, j=rel_ht.j.s)

    logger.info("Annotating the doubleton sample pairs with relatedness information...")
    ht = ht.annotate(
        rel_def=(
            hl.is_defined(rel_ht[ht.s1, ht.s2]) | hl.is_defined(rel_ht[ht.s2, ht.s1])
        ),
        kin=hl.coalesce(
            rel_ht[ht.s1, ht.s2].kin,
            rel_ht[ht.s2, ht.s1].kin,
        ),
        relationship=hl.coalesce(
            rel_ht[ht.s1, ht.s2].relationship,
            rel_ht[ht.s2, ht.s1].relationship,
        ),
    )
    ht = ht.checkpoint(f"{temp_path}/doubletons_uniq_rel.ht", overwrite=True)
    ht.show()

    def _get_agg_struct(ht: hl.Table) -> hl.expr.StructExpression:
        """
        Aggregate input Table and return StructExpression describing doubleton pairs.

        Return count of pairs present in relatedness HT, kinship distribution stats, and
        dictionary counting relationship types.

        Assumes Table is annotated with:
            - `rel_def`: Boolean for whether pair was present in relatedness Table.
            - `kin`: Kinship value for sample pair.
            - `relationship`: Relationship of sample pair (if found in relatedness Table).

        :param hl.Table ht: Input Table.
        :return: StructExpression describing doubleton pairs.
        """
        return ht.aggregate(
            hl.struct(
                pair_in_relatedness_ht=hl.agg.count_where(ht.rel_def),
                kin_stats=hl.agg.stats(ht.kin),
                rel_counter=hl.agg.counter(ht.relationship),
            )
        )

    logger.info(
        "Results from HT aggregate before removing 'unrelated' relationships: %s",
        _get_agg_struct(ht),
    )

    ht = ht.filter((hl.is_missing(ht.relationship)) | (ht.relationship != "unrelated"))
    logger.info(
        "Results from HT aggregate after removing 'unrelated' relationships: %s",
        _get_agg_struct(ht),
    )
