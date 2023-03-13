"""
This is a batch script which runs a vrs-annotation script on a VCF
Currently only outputs a NON-ZIPPED sharded VCF
TODO: possible to read from VCFs which are not only 100 partitions that need to be changed manually
TODO: add i to zip output 
usage: python3 vrs-annotation-batch-wrapper.py --billing-project gnomad-vrs \
    --tmp-dir gnomad-vrs \
    --bucket-mount gnomad-vrs \
    --image us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten-vrs-image:0221 \
    --vcf-in vcf-10k-sharded-0228.vcf.bgz \
    --vcf-out gs://gnomad-vrs/vcf-test-0302.vcf.bgz
"""

import argparse
import datetime
import errno
import logging
import os
import sys

import hail as hl
import hailtop.batch as hb
from batch_utils import init_job


def init_job_with_gcloud(
    batch,
    name: str = None,
    image: str = None,
    cpu: float = None,
    memory: float = None,
    disk_size: float = None,
    mount: str = None,
    mount_to: str = None,
):
    """
    Create job and initialize glcoud authentication and gsutil commands.
    Wraps Ben Weisburd's init_job (https://github.com/macarthur-lab/methods/blob/master/batch/batch_utils.py#L101) with additional gcloud steps.
    Parameters passed through to init_job (descriptions taken from github documentation for the function):
        :param batch: Batch object
        :param name: job label which will show up in the Batch web UI
        :param image: docker image name (eg. "weisburd/image-name@sha256:aa19845da5")
        :param cpu: number of CPUs (between 0.25 to 16)
        :param memory: amount of RAM in Gb (eg. 3.75)
        :param disk_size: amount of disk in Gb (eg. 50)
        :return: new job object
    """
    job = init_job(batch, name, image, cpu, memory, disk_size)
    job.command(f"gcloud -q auth activate-service-account --key-file=/gsa-key/key.json")
    job.command(f"curl -sSL broad.io/install-gcs-connector | python3")
    return job


def main(args):

    backend = hb.ServiceBackend(args.billing_project, args.tmp_dir)

    batch001 = hb.Batch(name="vrs", backend=backend)

    index_list = [str(yi).zfill(3) for yi in range(100)]
    # TODO: this only works for sharded VCFs with 100 partitions
    # Needs fixing and edits in future to be more universifiable

    for vcf_index in index_list:

        new_job = init_job_with_gcloud(
            batch=batch001, name=f"VCF_{vcf_index}_job", image=args.image
        )
        new_job.cloudfuse(f"{args.seqrepo_bucket}", "/local-seqrepo")
        # NOTE: the script gs://gnomad-vrs/marten-vrs-annotation.py has minor edits from 'vcf_annotation.py' from the GA4GH Team
        # namely, 'click' functionality for arguments has been changed to ArgParser
        # and some custom errors they included could not be updated
        # THESE WILL BE CHANGED IN THE FUTURE ! 
        new_job.command(
            f"python3 /local-seqrepo/marten-vrs-annotation.py --vcf-in /local-seqrepo/{args.vcf_in}/part-{vcf_index}*.bgz --vcf-out vcf-output-{vcf_index}.vcf --seqrepo-path /local-seqrepo/seqrepo/2018-11-26/"
        )
        new_job.command(f"gsutil cp vcf-output-{vcf_index}.vcf {args.vcf_out}/")

    batch001.run()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--billing-project", help="Project to bill.", type=str)
    parser.add_argument(
        "--tmp-dir",
        help="Name of google bucket for storing temp files.", type=str
    )
    parser.add_argument(
        "--image",
        default="us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten-vrs-image:0221",
        help="Image in a GCP Artifact Registry repository.",
        type=str
    )
    parser.add_argument(
        "--seqrepo-bucket", required=True, help="Name of bucket to mount containing download of seqrepo database in 'seqrepo/2018-11-26/'.", type=str
    )
    parser.add_argument(
        "--vcf-in",
        help="Path of SHARDED .VCF.BGZ to pass to annotator, RELATIVE to mounted bucket.",
        default="vcf-10k-sharded-0228.vcf.bgz", type=str
    )
    parser.add_argument("--vcf-out", help="Path and name of output VCF.", type=str)

    args = parser.parse_args()

    main(args)
