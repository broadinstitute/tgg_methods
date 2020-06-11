"""
This is intended as a starting point for Batch pipelines.

Batch docs:  https://hail.is/docs/batch/api/batch/hailtop.batch.job.Job.html#hailtop.batch.job.Job
"""

import argparse
import hail as hl  # used for hadoop file utils
import hailtop.batch as hb
import logging
import os
from batch.batch_utils import init_job, switch_gcloud_auth_to_user_account

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

OUTPUT_BUCKET = "gs://some-bucket"          # TODO edit this
DOCKER_IMAGE = "weisburd/some-image:latest" # TODO edit this


def main():
    p = batch_utils.init_arg_parser()
    args = p.parse_args()

    # see https://hail.zulipchat.com/#narrow/stream/223457-Batch-support/topic/auth.20as.20user.20account for more details
    working_dir = ""
    if args.local:
        if args.raw:
            backend = hb.LocalBackend()
            working_dir = os.path.abspath(os.path.dirname(__name__))
        else:
            backend = hb.LocalBackend(gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    else:
        backend = hb.ServiceBackend(billing_project=args.batch_billing_project, bucket=args.batch_temp_bucket)

    b = hb.Batch(backend=backend, name=args.batch_name)

    logger.info("Working dir: " + working_dir)

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
            if hl.hadoop_is_file(output_file_path) and not args.force:
                logger.info(f"{sample_id} output file already exists: {output_file_path}. Skipping...")
                continue

            file_stats = hl.hadoop_stat(bam_path)
            bam_size = int(round(file_stats['size_bytes']/10.**9))

        # init job commands
        logger.info(f'Requesting: {j._storage or "default"} storage, {j._cpu or "default"} CPU, {j._memory or "default"} memory')
        j = init_job(b, name=args.batch_job_name, cpu=args.cpu, memory=args.memory, disk_size=bam_size*2, image=DOCKER_IMAGE if not args.raw else None)
        switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

        # run commands
        j.command("set -x")
        if not args.raw:
            #j.command(f"gsutil -u {GCLOUD_PROJECT} -m cp gs://gtex-resources/GENCODE/gencode.v26.GRCh38.ERCC.genes.collapsed_only.gtf .")
            j.command(f"mv {genes_gtf} gencode.gff3")
            j.command(f"mv {input_read_data.bam} {sample_id}.bam")
            j.command(f"mv {input_read_data.bai} {sample_id}.bam.bai")
        j.command(f"ls")
        j.command(f"date")

        j.command(f"majiq build gencode.gff3 -c majiq_build.cfg -j 1 -o majiq_build_{sample_id}")
        j.command(f"tar czf majiq_build_{sample_id}.tar.gz majiq_build_{sample_id}")
        j.command(f"cp majiq_build_{sample_id}.tar.gz {j.output_tar_gz}")
        j.command(f"echo Done: {output_file_path}")
        j.command(f"date")

        # copy output
        b.write_output(j.output_tar_gz, output_file_path)
        print("Output file path: ", output_file_path)
    b.run()

    if isinstance(backend, hb.ServiceBackend):
        backend.close()


if __name__ == "__main__":
    main()
