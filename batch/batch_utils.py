"""
Shared methods for Batch pipelines.

Batch docs:  https://hail.is/docs/batch/api/batch/hailtop.batch.job.Job.html#hailtop.batch.job.Job
"""

import os
from hailtop.batch.job import Job


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
        if memory < 0.1 or memory > 64:
            raise ValueError(f"Memory arg is {memory}. This is outside the range of 0.1 to 64 Gb")

        j.memory(f"{memory}G")  # Batch default is 3.75G

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

    batch_job.command(f"gcloud auth activate-service-account --key-file /gsa-key/key.json")
    batch_job.command(f"gsutil -m cp -r {os.path.join(gs_path_of_gcloud_credentials, '.config')} /tmp/")
    batch_job.command(f"rm -rf ~/.config")
    batch_job.command(f"mv /tmp/.config ~/")
    batch_job.command(f"gcloud config set account {gcloud_user_account}")
    if gcloud_project:
        batch_job.command(f"gcloud config set project {gcloud_project}")

