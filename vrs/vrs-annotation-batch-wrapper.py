"""
This is a batch script which adds VRS IDs to a Hail Table by creating a sharded-VCF, running a vrs-annotation script on each shard, and annotating the merged results back onto the original Hail Table.
To be run using Query-On-Batch (https://hail.is/docs/0.2/cloud/query_on_batch.html#:~:text=Hail%20Query%2Don%2DBatch%20uses,Team%20at%20our%20discussion%20forum.)
usage: python3 vrs-annotation-batch-wrapper.py --billing-project gnomad-vrs \
    --tmp-dir gnomad-vrs \
    --bucket-mount gnomad-vrs \
    --image us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten-vrs-image-v2-031323 \
    --ht-input gs://gnomad-vrs/ht-inputs/ht-1k-TESTING-ONLY.ht \
    --ht-output example-ht-out.ht \
    --user-partitions 10 \
    --header gs://gnomad-vrs/header-fix.txt
"""

import argparse
import datetime
import errno
import logging
import os
import sys

import hail as hl
import hailtop.batch as hb
from tgg.batch.batch_utils import init_job


def init_job_with_gcloud(
    batch,
    name: str = None,
    image: str = None,
    cpu: float = None,
    memory: float = None,
    disk_size: float = None,
    mount: str = None,
):
    """
    Create job and initialize glcoud authentication and gsutil commands.
    Wraps Ben Weisburd's init_job (https://github.com/broadinstitute/tgg_methods/blob/master/tgg/batch/batch_utils.py#L160) with additional gcloud steps.
    Parameters passed through to init_job:
        :param batch: Batch object.
        :param name: Job label which will show up in the Batch web UI.
        :param image: Docker image name (eg. "us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten-vrs-image-v2-031323").
        :param cpu: Number of CPUs (between 0.25 to 16).
        :param memory: Amount of RAM in Gb (eg. 3.75).
        :param disk_size: Amount of disk in Gb (eg. 50).
        :param mount: Name of GCP Bucket to mount using cloudfuse.
        :return: New job object.
    """
    job = init_job(batch, name, image, cpu, memory, disk_size)
    job.command(f"gcloud -q auth activate-service-account --key-file=/gsa-key/key.json")
    job.command(f"curl -sSL broad.io/install-gcs-connector | python3")
    job.cloudfuse(mount, "/local-seqrepo")
    return job


def main(args):

    # prefix and datetime are used to construct unique names for temporary files
    prefix = (args.ht_output).split(".")[0]

    prefix += datetime.datetime.now().strftime("%m:%d:%Y:%H:%M:%S")

    # Reads in Hail Table, partitions, and exports to sharded VCF within the folder to-be-mounted
    ht0 = hl.read_table(f"{args.ht_input}")

    # NOTE: SELECT IS CURRENTLY USED TO REDUCE COSTS AND TIME DURING TESTING
    ht0 = ht0.select()

    # NOTE: REPARTITION DOES NOT OUTPUT EVENLY? LAST SHARD IS 5-TIMES THE SIZE OF OTHERS (STORAGE) AND MANY TIMES THE NUMBER OF VARIANTS
    # REACHED OUT TO HAIL TEAM ABOUT ISSUE, IN PROGRESS
    ht1 = ht0.repartition(args.partitions_for_vcf_export)

    hl.export_vcf(
        ht1,
        f"gs://{args.mount}/vrs-temp/shard-{prefix}.vcf.bgz",
        append_to_header=args.header_path,
        parallel="header_per_shard",
    )

    # Create a list of all shards of VCF
    ret_dict = hl.utils.hadoop_ls(
        f"gs://{args.mount}/vrs-temp/shard-{prefix}.vcf.bgz/part-*.bgz"
    )

    ret_list = [ret_item["path"].split("/")[-1] for ret_item in ret_dict]

    # Creating backend and batch for coming annotation batch jobs
    backend = hb.ServiceBackend(args.billing_project, args.tmp_dir)

    batch001 = hb.Batch(name="vrs-annotation", backend=backend)

    for vcf_index in ret_list:

        # Setting up job
        new_job = init_job_with_gcloud(
            batch=batch001,
            name=f"VCF_job_{vcf_index.split('.')[0].split('-')[:2]}",
            image=args.image,
            mount=args.mount,
        )

        # Script path for annotation is not expected to change for GA4GH if functionality is installed using pip in the DockerFile
        # ...unless GA4GH team changes their file structure and a new Image is used
        vrs_script_path = (
            "/usr/local/lib/python3.9/dist-packages/ga4gh/vrs/extras/vcf_annotation.py"
        )
        # Print which file is being annotated, create directory for annotated shard, and then perform annotation
        new_job.command(f"echo at {vcf_index}")
        new_job.command("mkdir /temp-vcf-annotated/")
        new_job.command(
            f"python3 {vrs_script_path} --vcf_in /local-seqrepo/vrs-temp/shard-{prefix}.vcf.bgz/{vcf_index} --vcf_out /temp-vcf-annotated/annotated-{vcf_index.split('.')[0]}.vcf --seqrepo_root_dir /local-seqrepo/seqrepo/2018-11-26/"
        )

        # Copy shard to its appropriate place in Google Bucket
        new_job.command(
            f"gsutil cp /temp-vcf-annotated/annotated-{vcf_index.split('.')[0]}.vcf gs://{args.mount}/vrs-temp/annotated-{prefix}.vcf/"
        )

    # Execute all jobs in Batch
    batch001.run()

    # Import all annotated shards
    ht2 = hl.import_vcf(
        f"gs://gnomad-vrs/vrs-temp/annotated-{prefix}.vcf/*.vcf",
        reference_genome="GRCh38",
    ).make_table()

    # Checkpoint (write) resulting annotated table
    ht2 = ht2.checkpoint(f"gs://{args.mount}/{args.ht_output}", overwrite=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--billing-project", help="Project to bill.", type=str)
    parser.add_argument(
        "--tmp-dir", help="Name of google bucket for storing temp files.", type=str
    )
    parser.add_argument(
        "--image",
        default="us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten-vrs-image-v2-031323",
        help="Image in a GCP Artifact Registry repository.",
        type=str,
    )
    parser.add_argument(
        "--mount",
        help="GCP Bucket to mount image to",
        default="gnomad-vrs",
        type=str,
    )

    parser.add_argument(
        "--ht-input",
        help="FULL PATH of Hail Table Input",
        default="gs://gnomad-vrs/ht-inputs/ht-1k-TESTING-ONLY.ht",
    )
    parser.add_argument(
        "--partitions-for-vcf-export",
        help="Number of partitions to use when exporting the Table to a sharded VCF (each partition is exported to a separate VCF). This value determines the breakdown of jobs that will be ran (the VRS annotation script is ran in parallel on the VCFs).",
        default=100,
        type=int,
    )
    parser.add_argument(
        "--ht-output", help="Name of HT output RELATIVE to mounted bucket.", type=str
    )
    parser.add_argument(
        "--header-path",
        help="Full path of txt file containing lines to append to VCF headers for fields that maybe be missing when exporting the Table to VCF.",
        type=str,
        default="gs://gnomad-vrs/header-fix.txt",
    )

    args = parser.parse_args()

    main(args)
