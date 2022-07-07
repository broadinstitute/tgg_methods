"""Utility functions to calculate singletons and plot distribution."""
import argparse
import logging
from typing import List, Set, Tuple

import hail as hl

from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.filtering import filter_to_adj
from gnomad.utils.vcf import SPARSE_ENTRIES
from gnomad.utils.slack import slack_notifications

from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import interval_qc_path, relatedness_ht_path
from ukbb_qc.resources.variant_qc import info_ht_path, NA12878, SYNDIP

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("singleton_utils")
logger.setLevel(logging.INFO)


TEMP_PATH = "gs://gnomad-tmp/singleton_analysis"
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

# Remove END from SPARSE_ENTRIES (not present in VDS variant data entries)
SPARSE_ENTRIES.remove("END")


def get_singleton_sites(
    vds_path: str = VDS_PATH,
    temp_path: str = TEMP_PATH,
    tranche_data: Tuple[str, int] = TRANCHE_DATA,
) -> hl.Table:
    """
    Filter UKB VDS to bi-allelic, autosomal sites in interval QC pass regions with an adj allele count of two and no homozygotes.

    :param vds_path: Path to UKB 455k VDS. Default is VDS_PATH.
    :param temp_path: Path to bucket to store Table and other temporary data. Default is TEMP_PATH.
    :param tranche_data: UKB tranche data (data source and data freeze number). Default is TRANCHE_DATA.
    :param sparse_entries: List of fields to select from VDS. Default is SPARSE_ENTRIES.
    :return: Table of high quality sites with doubletons.
    """
    logger.info("Reading in VDS and filtering to bi-allelic SNPs...")
    mt = hl.vds.read_vds(vds_path).variant_data
    # Drop unnecessary annotations
    mt = mt.select_rows()
    mt = mt.filter_rows(bi_allelic_expr(mt) & hl.is_snp(mt.alleles[0], mt.alleles[1]))

    logger.info("Filter to autosomes...")
    mt = mt.filter_rows(mt.locus.in_autosome())

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
    mt = mt.transmute_rows(ac=mt.call_stats.AC[1])

    logger.info("Filtering to an allele count of one and returning...")
    ht = mt.rows()
    ht = ht.filter(ht.ac == 1)
    ht = ht.checkpoint(f"{temp_path}/high_quality_singleton_sites.ht", overwrite=True)
    return ht


def get_samples_n_singletons(
    vds_path: str = VDS_PATH,
    temp_path: str = TEMP_PATH,
    control_samples: Set[str] = {NA12878, SYNDIP},
):
    """
    Get number of singletons per samples in the 455k VDS.

    Filter VDS variant data to sites present in specified input Table, collect sample IDs, annotate IDs onto rows, and write
    HT to temporary path.

    :param vds_path: Path to UKB 455k VDS. Default is VDS_PATH.
    :param temp_path: Path to bucket to store Table and other temporary data. Default is TEMP_PATH.
    :param control_samples: Set of control sample IDs to remove. Default is {NA12878, SYNDIP}.
    :return: Table keyed by sample IDs and their number of singletons.
    """
    logger.info("Counting singletons per sample...")
    mt = hl.vds.read_vds(vds_path).variant_data
    mt = mt.filter_cols(~hl.literal(control_samples).contains(mt.s))
    ht = get_singleton_sites(vds_path)
    mt = mt.annotate_rows(ac=ht[mt.row_key].ac)
    mt = mt.filter_rows(hl.is_defined(mt.ac))
    mt = mt.annotate_cols(
        singletons=hl.agg.count_where(mt.GT.is_non_ref() & (mt.ac == 1)),
        hail_n_singletons=hl.agg.sum(
            hl.rbind(
                mt["GT"],
                lambda gt: hl.sum(
                    hl.range(0, gt.ploidy).map(
                        lambda i: hl.rbind(gt[i], lambda gti: (gti != 0) & (mt.ac == 1))
                    )
                ),
            )
        ),
    )
    ht = mt.select_cols("singletons", "hail_n_singletons").cols()
    ht.write(f"{temp_path}/singletons_per_sample.ht", overwrite=True)
    return ht


def main(args):
    """Find number of singletons per sample in high quality sites."""
    try:
        hl.init(log="/singletons.log", default_reference="GRCh38")
        get_samples_n_singletons()

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(args.temp_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        """
        This script calculates the number of singletons in individuals in the 455k VDS.
        """
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
    )
    parser.add_argument(
        "--slack-token",
        help="Token to authenticate slack. Must be specified if --slack-channel is also set.",
    )
    parser.add_argument(
        "--temp-path",
        help="Path to temporary bucket to store hail logs.",
    )
    args = parser.parse_args()

    if args.slack_channel:
        from slack_creds import slack_token
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
