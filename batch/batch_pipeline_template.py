"""
This is intended as a starting point for Batch pipelines.

Batch docs:  https://hail.is/docs/batch/api/batch/hailtop.batch.job.Job.html#hailtop.batch.job.Job
"""

import argparse
import hail as hl  # used for hadoop file utils
import hailtop.batch as hb
import logging
import os
from batch.batch_utils import switch_gcloud_auth_to_user_account

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

OUTPUT_BUCKET = "gs://some-bucket"          # TODO edit this
DOCKER_IMAGE = "weisburd/some-image:latest" # TODO edit this

def parse_args():
    p = argparse.ArgumentParser()
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--cluster", action="store_true", help="Batch: submit to cluster")
    grp.add_argument("--local", action="store_true", help="Batch: run locally")
    p.add_argument("-r", "--raw", action="store_true", help="Batch: run directly on the machine, without using a docker image")
    p.add_argument("--batch-billing-project", default="tgg-rare-disease", help="Batch: billing project. Required if submitting to cluster.")
    p.add_argument("--batch-job-name", help="Batch: (optional) job name")
    p.add_argument("--cpu", type=float, help="Batch: (optional) number of CPUs to request", default=1, choices=[0.25, 0.5, 1, 2, 4, 8, 16])
    p.add_argument("--memory", type=float, help="Batch: (optional) memory in gigabytes", default=3.75)
    p.add_argument("--force", action="store_true", help="Recompute and overwrite cached or previously computed data")
    args = p.parse_args()

    return args


def main():

    args = parse_args()

    # see https://hail.zulipchat.com/#narrow/stream/223457-Batch-support/topic/auth.20as.20user.20account for more details
    working_dir = ""
    if args.local:
        if args.raw:
            backend = hb.LocalBackend()
            working_dir = os.path.abspath(os.path.dirname(__name__))
        else:
            backend = hb.LocalBackend(gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    else:
        backend = hb.ServiceBackend(args.batch_billing_project)

    b = hb.Batch(backend=backend, name=args.batch_job_name)

    # define workflow inputs
    if args.local:
        genes_gtf = b.read_input("gencode.v26.annotation.gff3", extension=".gff3")
    else:
        genes_gtf = b.read_input("gs://macarthurlab-rnaseq/ref/gencode.v26.annotation.GRCh38.gff3", extension=".gff3")

    data = []  # TODO edit input data

    # define parallel execution for samples
    for sample_id, bam_path, bai_path in data:

        # set job inputs & outputs
        input_read_data = b.read_input_group(bam=bam_path, bai=bai_path)

        output_dir = OUTPUT_BUCKET
        output_file_path = os.path.join(output_dir, f"some_file.tar.gz")

        # check if output file already exists
        if args.local:
            bam_size = None
        else:
            import hail as hl  # used for hadoop file utils
            if hl.hadoop_is_file(output_file_path) and not args.force:
                logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
                continue

            file_stats = hl.hadoop_stat(bam_path)
            bam_size = int(round(file_stats['size_bytes']/10.**9))

        # define majiq build commands for this sample
        j = b.new_job(name=args.batch_job_name)
        if not args.raw:
            j.image("weisburd/some-image:latest")
        if bam_size:
            j.storage(f'{bam_size*2}Gi')
        j.cpu(args.cpu)  # Batch default is 1
        if args.memory:
            j.memory(f"{args.memory}G")  # Batch default is 3.75G
        else:
            j.memory(f"{3.75*args.cpu}G")  # Batch default is 3.75G
        logger.info(f'Requesting: {j._storage or "default"} storage, {j._cpu or "default"} CPU, {j._memory or "default"} memory')

        # switch to user account
        switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

        # run commands
        j.command("set -x")
        if not args.raw:
            #j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp gs://gtex-resources/GENCODE/gencode.v26.GRCh38.ERCC.genes.collapsed_only.gtf .")
            j.command(f"mv {genes_gtf} gencode.gff3")
            j.command(f"mv {input_read_data.bam} {sample_id}.bam")
            j.command(f"mv {input_read_data.bai} {sample_id}.bam.bai")
        j.command(f"pwd >> {j.logfile}")
        j.command(f"ls >> {j.logfile}")
        j.command(f"date >> {j.logfile}")

        j.command(f"majiq build gencode.gff3 -c majiq_build.cfg -j 1 -o majiq_build_{sample_id} >> {j.logfile}")

        j.command(f"tar czf majiq_build_{sample_id}.tar.gz majiq_build_{sample_id}")
        j.command(f"cp majiq_build_{sample_id}.tar.gz {j.output_tar_gz}")
        j.command(f"echo Done: {output_file_path} >> {j.logfile}")
        j.command(f"date >> {j.logfile}")

        # copy output
        b.write_output(j.output_tar_gz, output_file_path)
        print("Output file path: ", output_file_path)

        log_file_path = os.path.join(output_dir, f"{sample_id}.log")
        b.write_output(j.logfile, log_file_path)
        print("Log file path: ", log_file_path)

    b.run()

    if isinstance(backend, hb.ServiceBackend):
        backend.close()


if __name__ == "__main__":
    main()
