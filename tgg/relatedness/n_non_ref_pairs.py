"""Utility functions to calculate shared sites among sample pairs"""
import argparse
import logging
from typing import Set, Tuple

import hail as hl

from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.filtering import filter_to_adj
from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import interval_qc_path
from ukbb_qc.resources.variant_qc import info_ht_path, NA12878, SYNDIP

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("n_non_ref_utils")
logger.setLevel(logging.INFO)


TEMP_PATH = "gs://gnomad-tmp/n_non_ref"
"""
Path to bucket for temporary data.
"""

VDS_PATH = "gs://ukbb-pharma-exome-analysis/500k_temp/500k.vds"
"""
Path to Variant Dataset (VDS) that contains high quality samples from the final UK Biobank callset.

Generated using the following commands:

# NOTE: If ever need to rerun, should read v4 VDS using `get_gnomad_v4_vds`
# (in gnomad_qc/v4/resoures/basics.py)
vds = hl.vds.read_vds("gs://gnomad/raw/exomes/4.0/gnomad_v4.0.vds")
meta_ht = hl.read_table("gs://broad-ukbb/broad.freeze_7/sample_qc/meta.ht")
meta_ht = meta_ht.filter(meta_ht.sample_filters.high_quality)
vds = hl.vds.filter_samples(vds, meta_ht, remove_dead_alleles=True)
call_stats_ht = hl.read_table(ukb.var_annotations_ht_path('ukb_freq', *TRANCHE_DATA[CURRENT_TRANCHE]))
freq_index = get_cohort_index(call_stats_ht)
var = vds.variant_data.annotate_rows(call_stats=call_stats_ht[vds.variant_data.row_key].freq[freq_index])
vds = hl.vds.VariantDataset(vds.reference_data, var)
vds.write(ukb_exomes_path, overwrite=args.overwrite)

Full VDS size: 8.69 TiB
VDS variant data size: 561.97 GiB
"""

TRANCHE_DATA = ("broad", CURRENT_FREEZE)
"""
UKB tranche data (data source and data freeze number)
"""


def get_n_non_ref_sites(
    vds_path: str = VDS_PATH,
    temp_path: str = TEMP_PATH,
    tranche_data: Tuple[str, int] = TRANCHE_DATA,
    autosomes_only: bool = True,
    snp_only: bool = False,
    het_only: bool = False,
    adj_only: bool = True,
    interval_qc_regions: bool = True,
    no_AS_lowqual: bool = True,
    non_ref_samples: int = 3,
    control_samples={NA12878, SYNDIP},
) -> hl.Table:
    """
    Filter VDS to autosomal sites in interval QC pass regions with an adj non reference allele count of 
    n and annotate site non ref samples.

    :param vds_path: Path to VDS. Default is VDS_PATH.
    :param temp_path: Path to bucket to store Table and other temporary data. Default is TEMP_PATH.
    :param tranche_data: UKB tranche data (data source and data freeze number). Default is TRANCHE_DATA.
    :param autosomes_only: Filter VDS to autosomes. Defaults to True.
    :param snp_only: Filter VDS to bi-allelic SNPs. Defaults to False.
    :param het_only: Filter to only het samples. Defaults to False.
    :param adj_only: Filter GT to adj. Defaults to True.
    :param interval_qc_regions: Filter to interval QC regions. Defaults to True.
    :param no_AS_lowqual: Remove AS_lowqual sites. Defaults to True.
    :param non_ref_samples: Desired number of non ref samples for each variant, e.g. for tripletons, n_samples=3. Defaults to 3.
    :param control_samples: Set of control samples to remove. Defaults to {NA12878, SYNDIP}. 
    :return: Table of high quality sites filtered to variants with specified number of non-reference samples.
    """
    logger.warning(
        "Reading in VDS. This script assumes the VDS multi-allelics have been split..."
    )
    vds = hl.vds.read_vds(vds_path)
    mt = vds.variant_data

    mt = mt.filter_cols(~hl.literal(control_samples).contains(mt.s))

    if autosomes_only:
        logger.info("Filtering to autosomes...")
        mt = mt.filter_rows(mt.locus.in_autosome())

    if snp_only:
        logger.info("Filtering to bi-allelic SNPs...")
        mt = mt.filter_rows(
            bi_allelic_expr(mt) & hl.is_snp(mt.alleles[0], mt.alleles[1])
        )

    if no_AS_lowqual:
        logger.info("Removing AS_lowqual sites...")
        info_ht = hl.read_table(info_ht_path(*tranche_data, split=True))
        mt = mt.filter_rows(~info_ht[mt.row_key].AS_lowqual)

    if interval_qc_regions:
        logger.info("Filtering to interval QC pass regions...")
        interval_ht = hl.read_table(interval_qc_path(*tranche_data, "autosomes"))
        mt = mt.filter_rows(hl.is_defined(interval_ht[mt.locus]))

    if adj_only:
        logger.info("Filtering to adj...")
        mt = filter_to_adj(mt)

    if "call_stats" in mt.row:
        logger.info("Changing old callstat annotation...")
        mt = mt.transmute_rows(original_call_stats=mt.call_stats)

    logger.info("Recalculating call stats post-quick QC...")
    mt = mt.annotate_rows(hail_call_stats=hl.agg.call_stats(mt.GT, mt.alleles))

    logger.info(
        "Get AC for alt allele and annotate site with the number of hom_var, het, hom_alt, and non ref samples..."
    )
    mt = mt.annotate_rows(
        # Get AC at allele index 1 (hail call_stats includes a count for each allele, including reference)
        ac=mt.hail_call_stats.AC[1],
        n_het=hl.agg.count_where(hl.is_defined(mt["GT"]))
        - hl.sum(mt.hail_call_stats.homozygote_count),
        n_hom_ref=mt.hail_call_stats.homozygote_count[0],
        n_hom_alt=mt.hail_call_stats.homozygote_count[1],
    )
    mt = mt.annotate_rows(
        n_non_ref=hl.agg.count_where(hl.is_defined(mt["GT"])) - mt.n_hom_ref
    )

    logger.info("Filtering sites to where n_non_ref = %i...", non_ref_samples)
    mt = mt.filter_rows(mt.n_non_ref == non_ref_samples)

    if het_only:
        mt = mt.filter_rows(mt.n_hom_alt == 0)

    logger.info("Finding non_ref samples at each site...")
    mt = mt.annotate_rows(
        samples=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.collect(mt.s))
    )

    ht = mt.rows().checkpoint(
        f"{temp_path}/n_non_ref_{non_ref_samples}_sample_high_quality_sites.ht",
        overwrite=True,
    )
    logger.info("Found %i sites where n_non_ref == %i", ht.count(), non_ref_samples)
    return ht


def get_and_count_sample_pairs(ht: hl.Table) -> hl.Table:
    """
    Return the number of shared n non_ref sites per pair.

    :param ht: Table to compute pairs on.
    :return ht: Table with sample pairs and number of variants shared per pair
    """
    logger.info("Collecting and counting sample pairs...")
    ht = ht.annotate(
        sample_pairs=hl.range(0, ht.samples.length()).flatmap(
            lambda i: hl.range(i + 1, ht.samples.length()).map(
                lambda j: hl.set([ht.samples[i], ht.samples[j]])
            )
        )
    )
    ht = ht.explode(ht.sample_pairs)
    logger.info("Aggregating shared sites per pair...")
    ht = ht.group_by(ht.sample_pairs).aggregate(n_non_ref_sites_shared=hl.agg.count())
    return ht


def get_samples_n_non_ref(
    vds_path: str = VDS_PATH,
    temp_path: str = TEMP_PATH,
    control_samples: Set[str] = {NA12878, SYNDIP},
    non_ref_samples: int = 3,
    het_only: bool = False,
) -> None:
    """
    Get number of shared non ref sites per sample pair in the 455k VDS.

    Filter VDS variant data to sites present in specified input Table, collect sample IDs, annotate IDs onto rows,
    explode on sample pairs, count pairs, and write HT to temporary path.

    :param vds_path: Path to UKB 455k VDS. Default is VDS_PATH.
    :param temp_path: Path to bucket to store Table and other temporary data. Default is TEMP_PATH.
    :param control_samples: Set of control sample IDs to remove. Default is {NA12878, SYNDIP}.
    :param non_ref_samples: Number of non_ref samples per site to filter to. Defaults to 3.
    :param het_only: Filter the non ref sample sites to those with only het sample calls. Defaults to False.
    """
    logger.info(
        "Retrieving pairwise shared %i non ref sample sites...", non_ref_samples
    )
    ht = get_n_non_ref_sites(
        vds_path=vds_path,
        temp_path=temp_path,
        non_ref_samples=non_ref_samples,
        het_only=het_only,
        control_samples=control_samples,
    )
    ht = get_and_count_sample_pairs(ht)
    ht.write(
        f"{temp_path}/pairwise_{non_ref_samples}_non_ref_shared_sites.ht",
        overwrite=True,
    )


def main(args):
    """Find number of requested non ref sites per sample pair in high quality sites."""
    try:
        hl.init(log="/n_non_ref.log", default_reference="GRCh38")
        get_samples_n_non_ref(
            vds_path=args.vds_path,
            temp_path=args.temp_path,
            non_ref_samples=args.non_ref_samples,
            het_only=args.het_only,
        )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(args.temp_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        """
        This script calculates the number of shared non ref sites per sample pair in the 455k VDS.
        """
    )
    parser.add_argument(
        "--slack-channel", help="Send message to Slack channel/user.",
    )
    parser.add_argument(
        "--slack-token",
        help="Token to authenticate slack. Must be specified if --slack-channel is also set.",
    )
    parser.add_argument(
        "--vds-path", help="Path to 455k UKB VDS.", default=VDS_PATH,
    )
    parser.add_argument(
        "--temp-path",
        help="Path to temporary bucket to store hail logs.",
        default=TEMP_PATH,
    )
    parser.add_argument(
        "--non-ref-samples",
        help="Number of samples per site with non-reference alleles to filter to. Defaults to 3.",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--het-only",
        help="Whether to only count sites where all samples are het.",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        from slack_creds import slack_token

        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
