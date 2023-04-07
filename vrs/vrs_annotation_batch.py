"""
This is a batch script which adds VRS IDs to a Hail Table by creating a sharded-VCF, running a vrs-annotation script on each shard, and annotating the merged results back onto the original Hail Table.
To be run using Query-On-Batch (https://hail.is/docs/0.2/cloud/query_on_batch.html#:~:text=Hail%20Query%2Don%2DBatch%20uses,Team%20at%20our%20discussion%20forum.)
usage: python3 vrs_annotation_batch.py \
    --billing-project gnomad-vrs \
    --working-bucket gnomad-vrs-io-finals \
    --image us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten-vrs-image-v2-031323 \
    --version test_v3_1k \
    --prefix marten_prelim_test \
    --partitions-for-vcf-export 20 \
    --downsample 0.1 \
    --header-path gs://gnomad-vrs-io-finals/header-fix.txt
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
from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.resources.annotations import vrs_annotations as v3_vrs_annotations
from tgg.batch.batch_utils import init_job


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("annotate_vrs_ids")
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

    hl.init(
        backend="batch",
        gcs_requester_pays_configuration="gnomad-vrs",  # NOTE: still pending permissions fix for writing to requester-pays bucket
        default_reference="GRCh38",
    )

    working_bucket = args.working_bucket
    version = args.version

    # Prefix to create custom named versions of outputs
    prefix = args.prefix + version

    input_paths_dict = {
        "v3.1.2": public_release("genomes").path,
        "test_v3.1.2": public_release("genomes").path,
        "test_v3_1k": "gs://gnomad-vrs-io-finals/ht-inputs/ht-1k-TESTING-ONLY-repartition-10p.ht",
        "test_v3_10k": "gs://gnomad-vrs-io-finals/ht-inputs/ht-10k-TESTING-ONLY-repartition-50p.ht",
        "test_v3_100k": "gs://gnomad-vrs-io-finals/ht-inputs/ht-100k-TESTING-ONLY-repartition-100p.ht",
    }

    output_paths_dict = {
        "v3.1.2": v3_vrs_annotations.path,
        "test_v3.1.2": f"gs://gnomad-vrs-io-finals/ht-outputs/{prefix}-Full-ht-release-output.ht",
        "test_v3_1k": f"gs://gnomad-vrs-io-finals/ht-outputs/{prefix}-Full-ht-1k-output.ht",
        "test_v3_10k": f"gs://gnomad-vrs-io-finals/ht-outputs/{prefix}-Full-ht-10k-output.ht",
        "test_v3_100k": f"gs://gnomad-vrs-io-finals/ht-outputs/{prefix}-Full-ht-100k-output.ht",
    }

    check_resource_existence(
        output_step_resources={
            "vrs-annotation-batch-wrapper.py": [
                f"gs://{working_bucket}/vrs-temp/outputs/VRS-{prefix}",
                output_paths_dict[version],
            ]
        },
        overwrite=args.overwrite,
    )

    # Read in Hail Table, partition, and export to sharded VCF within the folder to be mounted
    ht_original = hl.read_table(input_paths_dict[version])

    # Option to downsample for testing, if you want to annotate part of a Hail Table but not all of it
    if args.downsample < 1.00:
        logger.info("Downsampling Table...")
        ht_original = ht_original.sample(args.downsample)

    # Use 'select' to remove all non-key rows - VRS-Allele is added back to original table based on just locus and allele
    ht = ht_original.select()

    logger.info("Table read in and selected")

    """
    # Delete this block comment before merging - I do want to leave it in for now (as of 03/29) to keep track of progress at least
    # NOTE: REPARTITION DOES NOT OUTPUT EVENLY? LAST SHARD IS 5-TIMES THE SIZE OF OTHERS (STORAGE) AND MANY TIMES THE NUMBER OF VARIANTS
    # REACHED OUT TO HAIL TEAM ABOUT ISSUE, IN PROGRESS
    # 03-21 Response from Tim Poterba: an issue since the TESTING TABLES were sampled and subsetted randomly
    # --> Release Data is not affected in the same way
    # UPDATE: Tim has diagnosed problem, not yet fixed it yet on Query-On-Batch
    # UPDATE 03-28: All I had to do is Repartition them in a separate Notebook running Spark and re-output them
    # It's fine AND everything matches the trush set from the VRS team!!! woo!!!
    """

    # Repartition the Hail Table if requested. In the following step, the Hail Table is exported to a sharded VCF with one shard per partition
    if args.partitions_for_vcf_export:
        if args.partitions_for_vcf_export < ht.n_partitions():
            logger.info("Repartitioning Table by Naive Coalsece")
            ht = ht.naive_coalesce(args.partitions_for_vcf_export)
        elif args.partitions_for_vcf_export > ht.n_partitions():
            logger.info(
                "Repartitioning Table by Repartition, NOTE this results in a shuffle"
            )
            ht = ht.repartition(args.partitions_for_vcf_export)

    logger.info("Exporting sharded VCF")

    hl.export_vcf(
        ht,
        f"gs://{working_bucket}/vrs-temp/shard-{version}.vcf.bgz",
        append_to_header=args.header_path,
        parallel="header_per_shard",
    )

    logger.info(
        f"Gathered list of files in gs://{working_bucket}/vrs-temp/shard-{version}.vcf.bgz using Hail's Hadoop method"
    )

    # Create a list of all shards of VCF
    file_dict = hl.utils.hadoop_ls(
        f"gs://{working_bucket}/vrs-temp/shard-{version}.vcf.bgz/part-*.bgz"
    )

    # Create a list of all file names to later annotate in parallel
    file_list = [file_item["path"].split("/")[-1] for file_item in file_dict]

    logger.info("File list created... getting ready to start Batch Jobs")

    # Create backend and batch for coming annotation batch jobs
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
        new_job.command("mkdir /temp-vcf-annotated/")
        new_job.command(
            f"python3 {vrs_script_path} --vcf_in /local-vrs-mount/vrs-temp/shard-{version}.vcf.bgz/{vcf_index} --vcf_out /temp-vcf-annotated/annotated-{vcf_index.split('.')[0]}.vcf --seqrepo_root_dir /local-vrs-mount/seqrepo/2018-11-26/"
        )

        # Copy shard to its appropriate place in Google Bucket
        new_job.command(
            f"gsutil cp /temp-vcf-annotated/annotated-{vcf_index.split('.')[0]}.vcf gs://{working_bucket}/vrs-temp/annotated-{version}.vcf/"
        )

    # Execute all jobs in Batch
    batch_vrs.run()

    logger.info("Batch Jobs executed, preparing to read in sharded VCF from prior step")

    # Import all annotated shards
    ht_annotated = hl.import_vcf(
        f"gs://{working_bucket}/vrs-temp/annotated-{version}.vcf/*.vcf",
        reference_genome="GRCh38",
    ).make_table()
    logger.info("Annotated table constructed")

    # When you select something, even in the 'info' field, it is no longer in that struct!
    ht_annotated = ht_annotated.select(ht_annotated.info.VRS_Allele)

    # Checkpoint (write) resulting annotated table
    ht_annotated = ht_annotated.checkpoint(
        f"gs://{working_bucket}/vrs-temp/outputs/annotated-checkpoint-VRS-{prefix}.ht",
        overwrite=args.overwrite,
    )
    logger.info("Annotated Hail Table checkpointed")

    # NOTE: how to generalize this when we're not working on v3.1.2?
    # Output final Hail Tables with VRS annotations
    if "3.1.2" in version:
        # NOTE: could performance be improved?
        logger.info("Adding VRS IDs to original Table")
        ht_final = ht_original.annotate(
            info=ht_original.info.annotate(
                VRS_Allele=ht_annotated[
                    ht_original.locus, ht_original.alleles
                ].VRS_Allele
            )
        )
    else:
        logger.info("Constructing a test Table")
        ht_final = ht_original.annotate(
            VRS_Allele=ht_annotated[ht_original.locus, ht_original.alleles].VRS_Allele
        )

        logger.info(f"Outputting final table at: {output_paths_dict[version]}")
        ht_final.write(output_paths_dict[version], overwrite=args.overwrite)

    # Deleting temporary files saves a great deal of space and keeps CloudFuse costs down
    logger.info("Preparing to delete temporary files and sharded VCFs generated")
    delete_temps = hb.Batch(name="delete_temps", backend=backend)
    d1 = init_job_with_gcloud(
        batch=delete_temps,
        name=f"del_files",
        image=args.image,
        mount=working_bucket,
    )
    d1.command(f"gsutil -m -q rm -r gs://{working_bucket}/vrs-temp/")
    delete_temps.run()


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
        help="Name of GCP Bucket to mount for access.",
        default="gnomad-vrs-io-finals",
        type=str,
    )
    parser.add_argument(
        "--version", help="Version of HT to read in", default="test_v3_1k", type=str
    )
    parser.add_argument(
        "--partitions-for-vcf-export",
        help="Number of partitions to use when exporting the Table to a sharded VCF (each partition is exported to a separate VCF). This value determines the breakdown of jobs that will be ran (the VRS annotation script is ran in parallel on the VCFs).",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--prefix",
        help="Prefix to append to names of all saved temp files and test version output files.",
        type=str,
        default="",
    )
    parser.add_argument(
        "--header-path",
        help="Full path of txt file containing lines to append to VCF headers for fields that maybe be missing when exporting the Table to VCF.",
        type=str,
        default="gs://gnomad-vrs-io-finals/header-fix.txt",
    )
    parser.add_argument(
        "--downsample",
        help="Proportion to which to downsample the original Hail Table input.",
        default=1.00,
        type=float,
    )
    parser.add_argument(
        "--overwrite",
        help="Boolean to pass to ht.write(overwrite=_____) determining whether or not to overwrite existing output for the final Table and checkpointed files.",
        action="store_true",
    )

    args = parser.parse_args()

    main(args)
