"""
Shared methods for Batch pipelines.

Batch docs:  https://hail.is/docs/batch/api/batch/hailtop.batch.job.Job.html#hailtop.batch.job.Job
"""

import collections
import configargparse
import contextlib
import itertools
import logging
import os
import subprocess

import hailtop.batch as hb
from hailtop.batch.job import Job
from typing import Union

_GCLOUD_PROJECT = None

class HG38_REF_PATHS:
    fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    fai = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
    dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    gencode_v36_gtf = "gs://macarthurlab-rnaseq/ref/gencode.v36.annotation.gtf"

class HG37_REF_PATHS:
    fasta = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"
    fai = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai"
    dict = "gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict"


def set_gcloud_project(gcloud_project):
    global _GCLOUD_PROJECT
    _GCLOUD_PROJECT = gcloud_project


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
        batch = hb.Batch(backend=backend, name=batch_name)

        batch.batch_utils_temp_bucket = args.batch_temp_bucket

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

    j.command("set -euxo pipefail")  # set bash options for easier debugging and to make command execution more robust

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
    if gcloud_project or _GCLOUD_PROJECT:
        batch_job.command(f"gcloud config set project {gcloud_project or _GCLOUD_PROJECT}")

    batch_job.command(f"gcloud auth list")  # print auth list again to show that 'gcloud config set account' succeeded.

    # attach credentials to batch_job object
    batch_job._batch_utils_gs_path_of_gcloud_credentials = gs_path_of_gcloud_credentials
    batch_job._batch_utils_gcloud_user_account = gcloud_user_account

class StorageBucketRegionException(Exception):
    pass


def check_storage_bucket_region(google_storage_paths: Union[str, list], gcloud_project: str = None, verbose: bool = True):
    """Checks whether the given google storage path(s) are stored in US-CENTRAL1 - the region where the hail Batch
    cluster is located. Localizing data from other regions will be slower and result in egress charges.

    :param google_storage_paths: a gs:// path or a list of gs:// paths to check.
    :param gcloud_project: (optional) if specified, it will be added to the gsutil command with the -u arg.
    :raises StorageRegionException: If the given path(s) is not stored in the same region as the Batch cluster.

    """
    if isinstance(google_storage_paths, str):
        google_storage_paths = [google_storage_paths]

    buckets = set([path.split("/")[2] for path in google_storage_paths])
    for bucket in buckets:
        gsutil_command = f"gsutil"
        if gcloud_project or _GCLOUD_PROJECT:
            gsutil_command += f" -u {gcloud_project or _GCLOUD_PROJECT}"

        output = subprocess.check_output(f"{gsutil_command} ls -L -b gs://{bucket}", shell=True, encoding="UTF-8")
        for line in output.split("\n"):
            if "Location constraint:" in line:
                location = line.strip().split()[-1]
                break
        else:
            raise StorageBucketRegionException(f"ERROR: Couldn't determine gs://{bucket} bucket region.")

        if location not in {"US", "US-CENTRAL1"}:
            raise StorageBucketRegionException(f"ERROR: gs://{bucket} is located in {location}. This may cause egress "
                f"charges when copying files to the Batch cluster which is in US-CENTRAL.")

        if verbose:
            print(f"Confirmed gs://{bucket} is in {location}")


# dictionary that maps a job id to the set of buckets that have been gcsfuse-mounted into this job, to avoid mounting
# the same bucket 2x
_GCSFUSE_MOUNTED_BUCKETS_PER_JOB = collections.defaultdict(set)


def localize_file(job, google_storage_path: str, gcloud_project: str = None, use_gcsfuse: bool = False) -> str:
    """Copies a file from a google bucket to the local filesystem and returns the new absolute local path.
    Requires gsutil to exist inside the docker container.

    :param job: batch Job object
    :param google_storage_path: gs:// path of file to localize
    :param gcloud_project: (optional) if specified, it will be added to the gsutil command with the -u arg.
    :param use_gcsfuse: instead of copying the file, use gcsfuse to mount the bucket containing this file.
    :returns: Local file path after localization.
    """

    path = google_storage_path.replace("gs://", "")
    dirname = os.path.dirname(path)
    bucket_name = path.split("/")[0]

    if use_gcsfuse:
        root_dir = "/gcsfuse_mounts"
        local_bucket_dir = os.path.join(root_dir, bucket_name)
        local_file_path = os.path.join(root_dir, path)
        job_hash = hash(job)
        if bucket_name not in _GCSFUSE_MOUNTED_BUCKETS_PER_JOB[job_hash]:
            job.command(f"mkdir -p {local_bucket_dir}")
            job.gcsfuse(bucket_name, local_bucket_dir, read_only=True)
            _GCSFUSE_MOUNTED_BUCKETS_PER_JOB[job_hash].add(bucket_name)
    else:
        gsutil_command = f"gsutil"
        if gcloud_project or _GCLOUD_PROJECT:
            gsutil_command += f" -u {gcloud_project or _GCLOUD_PROJECT}"

        root_dir = "/localized"
        local_dir = os.path.join(root_dir, dirname)
        local_file_path = os.path.join(root_dir, path)
        job.command(f"mkdir -p '{local_dir}'; time {gsutil_command} -m cp -r '{google_storage_path}' '{local_file_path}'")

    job.command(f"ls -lh '{local_file_path}'")  # make sure file exists

    return local_file_path


# dictionary that maps a job id to the set of paths that have been localized via temp buckets. This avoids localizing
# the same path 2x and creating a race condition for deletion of this path from the temp bucket.
_PATHS_LOCALIZED_VIA_TEMP_BUCKET_PER_JOB = collections.defaultdict(set)


def localize_via_temp_bucket(
    job,
    google_storage_path: str,
    gcloud_project: str = None,
    temp_bucket: str = None,
    gs_path_of_gcloud_credentials: str = None,
    gcloud_user_account: str = None,
):
    """This method provides a work-around to allow using gcs_fuse on a bucket where you can't grant read access to the
    Batch service account. It works by
    1. copying the file to a temp bucket that you control (and where you have granted read access to the Batch service account)
    2. localizing that file to your job using gcs_fuse
    3. deleting the temp copy created in step 1

    Requires gsutil to exist inside the docker container.

    :param job: batch Job object
    :param google_storage_path: gs:// path of file to localize
    :param gcloud_project: (optional) if specified, it will be added to the gsutil command with the -u arg.
    :param use_gcsfuse: instead of copying the file, use gcsfuse to mount the bucket containing this file.
    :param temp_bucket: (optional) temp bucket to use.
    :param gs_path_of_gcloud_credentials: (not needed if switch_gcloud_auth_to_user_account was called on job)
        this arg will be forwarded to switch_gcloud_auth_to_user_account (see that method for more details)
    :param gcloud_user_account: (not needed if switch_gcloud_auth_to_user_account was called on job)
        this arg will be forwarded to switch_gcloud_auth_to_user_account (see that method for more details)
    :returns: Local file path after localization.
    """
    gsutil_command_prefix = f"gsutil"
    if gcloud_project or _GCLOUD_PROJECT:
        gsutil_command_prefix += f" -u {gcloud_project or _GCLOUD_PROJECT}"

    batch = job._batch  # the _batch attribute is set within Batch when the Job object is created.
    batch_id = abs(hash(batch)) % 10**9
    job_id = abs(hash(job)) % 10**9

    temp_bucket = temp_bucket or getattr(batch, "batch_utils_temp_bucket", None)
    if not temp_bucket:
        raise ValueError("temp bucket not specified.")

    temp_file_path = f"gs://{temp_bucket}/batch_{batch_id}/job_{job_id}/" + google_storage_path.replace("gs://", "")

    if temp_file_path in _PATHS_LOCALIZED_VIA_TEMP_BUCKET_PER_JOB[job_id]:
        raise ValueError(f"{google_storage_path} has already been localized via temp bucket.")
    _PATHS_LOCALIZED_VIA_TEMP_BUCKET_PER_JOB[job_id].add(temp_file_path)

    job.command(f"{gsutil_command_prefix} -m cp -r {google_storage_path} {temp_file_path}")

    # temp file cleanup job
    cleanup_job_name = (f"{job.name} " if job.name else "") + f"cleanup {os.path.basename(temp_file_path)}"
    cleanup_job = init_job(batch, cleanup_job_name, job._image, cpu=0.25, memory=0.25*3.75)
    cleanup_job.always_run()

    gs_path_of_gcloud_credentials = gs_path_of_gcloud_credentials or getattr(job, "_batch_utils_gs_path_of_gcloud_credentials", None)
    gcloud_user_account = gcloud_user_account or getattr(job, "_batch_utils_gcloud_user_account", None)
    if not gs_path_of_gcloud_credentials:
        raise ValueError("gs_path_of_gcloud_credentials not specified")
    if not gcloud_user_account:
        raise ValueError("gcloud_user_account not specified")
    switch_gcloud_auth_to_user_account(cleanup_job, gs_path_of_gcloud_credentials, gcloud_user_account)

    cleanup_job.command(f"gsutil -m rm -r {temp_file_path}")
    cleanup_job.depends_on(job)

    return localize_file(job, temp_file_path, gcloud_project=gcloud_project, use_gcsfuse=True)


def generate_path_to_file_size_dict(glob):
    """Runs "gsutil ls -l {glob}" and returns a dictionary that maps each gs:// file to
    its size in bytes. This appears to be faster than running hl.hadoop_ls(..).
    """
    logging.info(f"Listing {glob}")
    try:
        gsutil_output = subprocess.check_output(
            f"gsutil -m ls -l {glob}",
            shell=True,
            stderr=subprocess.STDOUT,
            encoding="UTF-8")
    except subprocess.CalledProcessError as e:
        if "One or more URLs matched no objects." in e.output:
            return {}
        else:
            raise e

    records = [r.strip().split("  ") for r in gsutil_output.strip().split("\n") if not r.startswith("TOTAL: ")]
    return {r[2]: int(r[0]) for r in records}  # map path to size in bytes


def batch_iter(iterable, batch_size=1):
    """Takes a list, set, tuple or any other iterable and breaks it up into batches.

    Args:
        iterable: a list, tuple or any other object that can be looped over
        batch_size (int): size of batches to return

    Yields:
         batches of size 'batch_size'.
    """

    if batch_size < 1:
        raise ValueError(f"batch_size={{batch_size}}. It needs to be a positive integer.")

    it = iter(iterable)
    while True:
        batch = tuple(itertools.islice(it, batch_size))
        if not batch:
            break
        yield batch


_PREV_MEMORY_BYTES = 0


def print_memory_stats(message="", run_gc=False):
    """Prints current memory usage, as well as change in usage since this method was last called.

    Args:
        message (str): Print this message along with the memory usage.
        run_gc (bool): If True, calls gc.collect() before printing memory use.
    """
    import gc
    import psutil

    global _PREV_MEMORY_BYTES
    if message:
        message = " - " + message

    if run_gc:
        gc.collect()

    memory_bytes = psutil.Process(os.getpid()).memory_info().rss

    logging.info(
        f"memory used {message}: {memory_bytes//10**6} Mb    delta: {(memory_bytes - _PREV_MEMORY_BYTES)//10**6} Mb")

    _PREV_MEMORY_BYTES = memory_bytes

