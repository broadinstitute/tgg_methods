## THIS SCRIPT NEEDS SSDs!!!!
import argparse

import hail as hl
from gnomad.utils.slack import try_slack


def main(args):

    hl.init(log="/vcf_write.log", default_reference=f"{args.build}")

    mt = hl.import_vcf(
        f"{args.genotypes_vcf}", force_bgz=True, find_replace=("nul", ".")
    )
    filters_ht = hl.import_vcf(
        f"{args.filtered_vcf}", force_bgz=True, find_replace=("nul", ".")
    ).rows()
    mt = mt.annotate_rows(filters=filters_ht[mt.row_key].filters)
    hl.export_vcf(
        mt,
        f"gs://seqr-datasets/v02/{args.build}/RDG_{args.data_type}_Broad_{args.data_source}/v{args.version}/RDG_{args.data_type}_Broad_{args.data_source}.vcf.bgz",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genotypes-vcf",
        help="Path to genotypes vcf, i.e. 'gs://fc-a5707a81-6c5d-4ddb-aa81-5eb31a484335/RDG_Broad_WES_Internal_Mar2020/part_one_outputs/chr*/RDG_Broad_WES_Internal_chr*.*.hard_filtered_with_genotypes.vcf.gz'",
        required=True,
    )
    parser.add_argument(
        "--filtered-vcf",
        help="Path to filtered VCF, i.e. 'gs://fc-a5707a81-6c5d-4ddb-aa81-5eb31a484335/RDG_Broad_WES_Internal_Mar2020/part_two_outputs/RDG_Broad_WES_Internal.filtered.*.vcf.gz'",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--callset-version",
        help="Version of callset vcf",
        type=int,
        required=True,
    )
    parser.add_argument(
        "--data-source",
        help="Data source (Internal or External)",
        choices=["Internal", "External"],
        required=True,
    )
    parser.add_argument(
        "--data-type",
        help="Sequencing data type (WES or WGS)",
        choices=["WES", "WGS"],
        required=True,
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
