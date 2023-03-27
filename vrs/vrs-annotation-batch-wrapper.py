"""
This is a batch script which adds VRS IDs to a Hail Table by creating a sharded-VCF, running a vrs-annotation script on each shard, and annotating the merged results back onto the original Hail Table.
To be run using Query-On-Batch (https://hail.is/docs/0.2/cloud/query_on_batch.html#:~:text=Hail%20Query%2Don%2DBatch%20uses,Team%20at%20our%20discussion%20forum.)
usage: python3 vrs-annotation-batch-wrapper.py \
    --billing-project gnomad-vrs \
    --working-bucket gnomad-vrs \
    --image us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten-vrs-image-v2-031323 \
    --version test_v3_1k \
    --prefix marten_prelim_test \
    --partitions-for-vcf-export 20 \
    --header-path gs://gnomad-vrs/header-fix.txt
"""

import argparse
import datetime
import errno
import logging
import os
import sys

import hail as hl
import hailtop.batch as hb
from gnomad.resources.grch38.gnomad import public_release
from gnomad_qc.v3.resources.annotations import vrs_annotations as v3_vrs_annotations
from tgg.batch.batch_utils import init_job

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_vrs_test_dataset")
logger.setLevel(logging.INFO)


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
    job.cloudfuse(mount, "/local-vrs-mount")
    return job


def main(args):

    # Working bucket definition
    working_bucket = args.working_bucket 

    # prefix to create custom named versions of outputs 
    prefix = args.prefix + args.version 

    # Input paths (fixed) - note that the Hail Tables are altered outside of the script and may partition oddly
    input_paths_dict = {
        "3.1.2": public_release("genomes").path ,
        "test_v3_1k": "gs://gnomad-vrs/ht-inputs/ht-1k-TESTING-ONLY.ht",
        "test_v3_10k": "gs://gnomad-vrs/ht-inputs/ht-10k-TESTING-ONLY.ht",
        "test_v3_100k": "gs://gnomad-vrs/ht-inputs/ht-100k-TESTING-ONLY.ht",
    }

    output_paths_dict = {
        "3.1.2": v3_vrs_annotations.path,
        "test_v3_1k": f"gs://gnomad-vrs/vrs-temp/outputs/{prefix}-Full-ht-1k-output-full-annotated.ht",
        "test_v3_10k": f"gs://gnomad-vrs/vrs-temp/outputs/{prefix}-Full-ht-10k-output-full-annotated.ht",
        "test_v3_100k": f"gs://gnomad-vrs/vrs-temp/outputs/{prefix}-Full-ht-100k-output-full-annotated.ht",
    }

    # Reads in Hail Table, partitions, and exports to sharded VCF within the folder to-be-mounted
    ht_original = hl.read_table(input_paths_dict[args.version])

    # Option to downsample for testing, if you want to test on v3.1.2 but not all of it
    if args.downsample < 1.00 and args.version == "3.1.2":
        ht_original = ht_original.sample(args.downsample)

    # Select() removed all non-key rows - VRS-Allele here is then added back to original table
    ht = (
        ht_original.select()
    )  # QUESTION: is this okay or would this be cluttering the namespace ?

    # NOTE: REPARTITION DOES NOT OUTPUT EVENLY? LAST SHARD IS 5-TIMES THE SIZE OF OTHERS (STORAGE) AND MANY TIMES THE NUMBER OF VARIANTS
    # REACHED OUT TO HAIL TEAM ABOUT ISSUE, IN PROGRESS
    # 03-21 Response from Tim Poterba: an issue since the TESTING TABLES were sampled and subsetted randomly
    # --> Release Data is not affected in the same way
    # UPDATE: Tim has diagnosed problem, not yet fixed it yet on Query-On-Batch
    ht = ht.repartition(args.partitions_for_vcf_export)

    hl.export_vcf(
        ht,
        f"gs://{working_bucket}/vrs-temp/shard-{prefix}.vcf.bgz",
        append_to_header=args.header_path,
        parallel="header_per_shard",
    )

    # Create a list of all shards of VCF
    file_dict = hl.utils.hadoop_ls(
        f"gs://{working_bucket}/vrs-temp/shard-{prefix}.vcf.bgz/part-*.bgz"
    )

    # Create a list of all file names to later annotate in parallel
    file_list = [file_item["path"].split("/")[-1] for file_item in file_dict]

    # Creating backend and batch for coming annotation batch jobs
    backend = hb.ServiceBackend(args.billing_project, working_bucket)

    batch_vrs = hb.Batch(name="vrs-annotation", backend=backend)

    for vcf_index in file_list:

        # Setting up job
        new_job = init_job_with_gcloud(
            batch=batch_vrs,
            name=f"VCF_job_{vcf_index.split('.')[0].split('-')[1]}",
            image=args.image,
            mount=working_bucket,
        )

        # Script path for annotation is not expected to change for GA4GH if functionality is installed using pip in the DockerFile
        # ...unless GA4GH team changes their file structure and a new Image is used
        vrs_script_path = (
            "/usr/local/lib/python3.9/dist-packages/ga4gh/vrs/extras/vcf_annotation.py"
        )
        # Print which file is being annotated, create directory for annotated shard, and then perform annotation
        new_job.command(f"echo now on: {vcf_index}")
        new_job.command("mkdir /temp-vcf-annotated/")
        new_job.command(
            f"python3 {vrs_script_path} --vcf_in /local-vrs-mount/vrs-temp/shard-{prefix}.vcf.bgz/{vcf_index} --vcf_out /temp-vcf-annotated/annotated-{vcf_index.split('.')[0]}.vcf --seqrepo_root_dir /local-vrs-mount/seqrepo/2018-11-26/"
        )

        # Copy shard to its appropriate place in Google Bucket
        new_job.command(
            f"gsutil cp /temp-vcf-annotated/annotated-{vcf_index.split('.')[0]}.vcf gs://{working_bucket}/vrs-temp/annotated-{prefix}.vcf/"
        )

    # Execute all jobs in Batch
    batch_vrs.run()

    # Import all annotated shards
    ht_annotated = hl.import_vcf(
        f"gs://gnomad-vrs/vrs-temp/annotated-{prefix}.vcf/*.vcf",
        reference_genome="GRCh38",
    ).make_table()
    logger.info("Annotated table constructed")  # turn into logger statement

    # When you select something, even in the 'info' field, it is no longer in that struct!
    ht_annotated = ht_annotated.select(ht_annotated.info.VRS_Allele)

    # Checkpoint (write) resulting annotated table
    ht_annotated = ht_annotated.checkpoint(
        f"gs://{working_bucket}/vrs-temp/outputs/VRS-{prefix}", overwrite=args.overwrite
    )
    logger.info("Annotated Hail Table checkpointed")

    if "test" not in args.version:
        # NOTE: this construction is very slow on a large scale, maybe consider a Join ?
        logger.info("Constructing a final table")
        ht_final = ht_original.annotate(
            info=ht_original.info.annotate(
                VRS_Allele=ht_annotated[
                    ht_original.locus, ht_original.alleles
                ].VRS_Allele
            )
        )
    else:
        logger.info("Constructing a test table")
        ht_final = ht_original.annotate(
            VRS_Allele=ht_annotated[ht_original.locus, ht_original.alleles].VRS_Allele
        )

    logger.info(f"Outputting final table at: {output_paths_dict[args.version]}")
    ht_final.write(output_paths_dict[args.version],overwrite=args.overwrite)
    # hl.write_table(ht_final,output_paths_dict[args.version])


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--billing-project", help="Project to bill.", type=str)
    parser.add_argument(
        "--image",
        default="us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten-vrs-image-v2-031323",
        help="Image in a GCP Artifact Registry repository.",
        type=str,
    )
    parser.add_argument(
        "--working-bucket",
        help="Name GCP Bucket to mount for access",
        default="gnomad-vrs",
        type=str,
    )
    parser.add_argument(
        "--version", help="Version of HT to read in", default="test_v3_1k", type=str
    )
    parser.add_argument(
        "--partitions-for-vcf-export",
        help="Number of partitions to use when exporting the Table to a sharded VCF (each partition is exported to a separate VCF). This value determines the breakdown of jobs that will be ran (the VRS annotation script is ran in parallel on the VCFs).",
        default=100,
        type=int,
    )
    parser.add_argument(
        "--prefix",
        help="Prefix to append to names of all saved temp or output files.",
        type=str,
        default="",
    )
    parser.add_argument(
        "--header-path",
        help="Full path of txt file containing lines to append to VCF headers for fields that maybe be missing when exporting the Table to VCF.",
        type=str,
        default="gs://gnomad-vrs/header-fix.txt",
    )
    parser.add_argument(
        "--downsample",
        help="If reading in the whole Release Table, option to downsample",
        default=1.00,
        type=float,
    )
    parser.add_argument(
        "--overwrite",
        help="Boolean to pass to ht.write(overwrite=_____)",
        action='store_true'
    )

    args = parser.parse_args()

    main(args)
