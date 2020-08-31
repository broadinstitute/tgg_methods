"""
Shared methods for Batch pipelines.

Batch docs:  https://hail.is/docs/batch/api/batch/hailtop.batch.job.Job.html#hailtop.batch.job.Job
"""

import configargparse
import contextlib
import os

import hailtop.batch as hb
from hailtop.batch.job import Job


class HG38_REF_PATHS:
    fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    fai = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"

class HG37_REF_PATHS:
    fasta = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"
    fai = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai"
    dict = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict"


def init_arg_parser(
    default_billing_project="tgg-rare-disease",
    default_temp_bucket="macarthurlab-cromwell",
    default_cpu=1,
    default_memory=3.75,
    parser=configargparse.ArgumentParser(formatter_class=configargparse.ArgumentDefaultsRawHelpFormatter),
    gsa_key_file=None,
):
    """Initializes and returns an argparse instance with common pipeline args pre-defined."""

    local_or_cluster_grp = parser.add_mutually_exclusive_group(required=True)
    local_or_cluster_grp.add_argument("--local", action="store_true", help="Batch: run locally")
    local_or_cluster_grp.add_argument("--cluster", action="store_true", help="Batch: submit to cluster")
    parser.add_argument("-r", "--raw", action="store_true", help="Batch: run directly on the machine, without using a docker image")

    parser.add_argument("--gsa-key-file", default=gsa_key_file, help="Batch: path of gcloud service account .json "
        "key file. If provided, Batch will mount this file into the docker image so gcloud commands can run as this service account.")
    parser.add_argument("--batch-billing-project", default=default_billing_project, help="Batch: this billing project will be "
        "charged when running jobs on the Batch cluster. To set up a billing project name, contact the hail team.")
    parser.add_argument("--batch-name", nargs="+", help="Batch: label for the current Batch run")
    parser.add_argument("--batch-temp-bucket", default=default_temp_bucket, help="Batch: bucket where it stores temp "
        "files. The batch service-account must have Admin permissions for this bucket. These can be added by running "
        "gsutil iam ch serviceAccount:[SERVICE_ACCOUNT_NAME]:objectAdmin gs://[BUCKET_NAME]")
    parser.add_argument("-t", "--cpu", type=float, default=default_cpu, choices=[0.25, 0.5, 1, 2, 4, 8, 16], help="Batch: number of CPUs (eg. 0.5)")
    parser.add_argument("-m", "--memory", type=float, default=default_memory, help="Batch: memory in gigabytes (eg. 3.75)")
    parser.add_argument("-f", "--force", action="store_true", help="Recompute and overwrite cached or previously computed data")
    parser.add_argument("--start-with", type=int, help="Start from this step in the pipeline")
    parser.add_argument("--dry-run", action="store_true", help="Don't run commands, just print them.")
    parser.add_argument("--verbose", action="store_true", help="Verbose log output.")
    return parser


@contextlib.contextmanager
def run_batch(args, batch_name=None):
    """Wrapper around creating, running, and then closing a Batch run.

    :param args: Parsed args from the ArgumentParser created via the init_arg_parser method
    :param batch_name: (optional) batch label which will show up in the Batch web UI

    Usage:
        with run_batch(args) as batch:
            ... batch job definitions ...
    """

    if args.local:
        backend = hb.LocalBackend() if args.raw else hb.LocalBackend(gsa_key_file=args.gsa_key_file)
    else:
        backend = hb.ServiceBackend(billing_project=args.batch_billing_project, bucket=args.batch_temp_bucket)

    try:
        batch = hb.Batch(backend=backend, name=" ".join(args.batch_name) if args.batch_name else batch_name)

        yield batch  # returned to with ... as batch:

        # run on end of with..: block
        batch.run(dry_run=args.dry_run, verbose=args.verbose)

    finally:
        if isinstance(backend, hb.ServiceBackend):
            backend.close()


def init_job(
        batch,
        name: str = None,
        image: str = None,
        cpu: float = None,
        memory: float = None,
        disk_size: float = None,
):
    """Common job init steps

    :param batch: Batch object
    :param name: job label which will show up in the Batch web UI
    :param image: docker image name (eg. "weisburd/image-name@sha256:aa19845da5")
    :param cpu: number of CPUs (between 0.25 to 16)
    :param memory: amount of RAM in Gb (eg. 3.75)
    :param disk_size: amount of disk in Gb (eg. 50)
    :return: new job object
    """

    j = batch.new_job(name=name)
    if image:
        j.image(image)

    if cpu:
        if cpu < 0.25 or cpu > 16:
            raise ValueError(f"CPU arg is {cpu}. This is outside the range of 0.25 to 16 CPUs")

        j.cpu(cpu)  # Batch default is 1

    if memory:
        if memory < 0.1 or memory > 60:
            raise ValueError(f"Memory arg is {memory}. This is outside the range of 0.1 to 60 Gb")

        j.memory(f"{memory}Gi")  # Batch default is 3.75G

    if disk_size:
        if disk_size < 1 or disk_size > 1000:
            raise ValueError(f"Disk size arg is {disk_size}. This is outside the range of 1 to 1000 Gb")

        j.storage(f'{disk_size}Gi')

    j.command("set -x")  # shell should print commands before running them - this can be useful for debugging

    return j


def switch_gcloud_auth_to_user_account(
    batch_job: Job,
    gs_path_of_gcloud_credentials: str,
    gcloud_user_account: str,
    gcloud_project: str = None,
):
    """This method adds some shell commands to your Batch job to switch gcloud auth from the Batch-provided service
    account to your user account.

    This can be used to access all your google buckets without having to first grant access to the Batch service account

    For this to work, you must first
    1) create a google bucket that only you have access to - for example: gs://weisburd-gcloud-secrets/
    2) on your local machine, make sure you're logged in to gcloud by running
		   gcloud auth login
	3) copy your local ~/.config directory (which caches your gcloud auth credentials) to the secrets bucket from step 1
		   gsutil -m cp -r ~/.config/  gs://weisburd-gcloud-secrets/
    4) grant your default Batch service-account read access to your secrets bucket so it can download these credentials
       into each docker container.
    5) make sure gcloud & gsutil are installed inside the docker images you use for your Batch jobs
    6) call this method at the beginning of your batch job:

    Example:
          switch_gcloud_auth_to_user_account(
            batch_job,
            "gs://weisburd-gcloud-secrets",
            "weisburd@broadinstitute.org",
            "seqr-project")

    :param batch_job: Batch job object
    :param gs_path_of_gcloud_credentials: google bucket path that contains your .config folder
    :param gcloud_user_account: user account to activate
    :param gcloud_project: (optional) set this as the default gcloud project
    :return:
    """

    batch_job.command(f"gcloud auth list")
    batch_job.command(f"gcloud auth activate-service-account --key-file /gsa-key/key.json")
    batch_job.command(f"gsutil -m cp -r {os.path.join(gs_path_of_gcloud_credentials, '.config')} /tmp/")
    batch_job.command(f"rm -rf ~/.config")
    batch_job.command(f"mv /tmp/.config ~/")
    batch_job.command(f"gcloud config set account {gcloud_user_account}")
    if gcloud_project:
        batch_job.command(f"gcloud config set project {gcloud_project}")
    batch_job.command(f"gcloud auth list")

