"""
This is intended as a starting point for Batch pipelines.

Batch docs:  https://hail.is/docs/batch/api/batch/hailtop.batch.job.Job.html#hailtop.batch.job.Job
"""

import hail as hl   # used for hadoop file utils
import logging
import os
import pandas as pd
from batch import batch_utils

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "seqr-project"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

DOCKER_IMAGE = "weisburd/base-bam-tools@sha256:83f91e81dbf9562e425e2afc20b6efc8a70e9caee081e40e9dd9a5b1874cb79f"
OUTPUT_DIR = "gs://rgp-fastq-files-for-illumina"


def main():
    p = batch_utils.init_arg_parser(default_cpu=8, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("-R", "--reference", choices={"37", "38"}, default="38")
    p.add_argument("--dont-sort", action="store_true", help="The bam/cram needs to be sorted by read name prior to "
        "fastq conversion. If it's already sorted this way, use this option to skip the sorting step.")
    p.add_argument("tsv_path", help="Table with 1 header line and columns: sample_id, cram_path")
    p.add_argument("sample_id", nargs="*", help="(optional) 1 or more sample_id prefixes to process. If not specified, all rows in the .tsv will be processed.")
    args = p.parse_args()

    df = pd.read_table(args.tsv_path)
    if {"sample_id", "cram_path"} - set(df.columns):
        p.error(f"{args.tsv_path} must contain a 'sample_id' and 'cram_path' columns")

    if not args.force:
        hl.init(log="/dev/null", quiet=True)

    if args.sample_id:
        missing = 0
        for sample_id in args.sample_id:
            if not any(df.sample_id.str.startswith(sample_id)):
                missing += 1
                logger.error(f"Couldn't find any matches for {sample_id}")
        if missing:
            p.error(f"Couldn't find matches for {missing} of the {len(args.sample_id)} args")
                
    # process samples
    with batch_utils.run_batch(args, batch_name=f"cram => fastq") as batch:
        for _, row in df.iterrows():
            if args.sample_id and not any([row.sample_id.startswith(prefix) for prefix in set(args.sample_id)]):
                continue

            input_filename = os.path.basename(row.cram_path)
            prefix = input_filename.replace(".bam", "").replace(".cram", "")

            output_filename_r1 = f"{prefix}__R1.fastq.gz"
            output_filename_r2 = f"{prefix}__R2.fastq.gz"
            output_path_r1 = os.path.join(OUTPUT_DIR, output_filename_r1)
            output_path_r2 = os.path.join(OUTPUT_DIR, output_filename_r2)

            if not args.force and hl.hadoop_is_file(output_path_r1) and hl.hadoop_is_file(output_path_r2):
                logger.info(f"Output files exist (eg. {output_filename_r1}). Skipping {input_filename}...")
                continue

            j = batch_utils.init_job(batch, f"cram => fastq: {row.sample_id}", DOCKER_IMAGE if not args.raw else None, args.cpu, args.memory)
            batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

            # copy inputs
            REF_PATHS = batch_utils.HG38_REF_PATHS if args.reference == "38" else batch_utils.HG37_REF_PATHS
            j.command(f"""gsutil -m -u {GCLOUD_PROJECT} cp {row.cram_path} .""")
            j.command(f"""gsutil -m cp {REF_PATHS.fasta} {REF_PATHS.fai} {REF_PATHS.dict} .""")

            if not args.dont_sort:
                j.command(f"samtools sort -n -@ {max(1, int(args.cpu))} -m 3G  {input_filename} -o {prefix}.sorted.cram")
            j.command(f"samtools fastq -1 {output_filename_r1} -2 {output_filename_r2} {prefix}.sorted.cram")
            j.command(f"""gsutil -m cp {output_filename_r1} {output_filename_r2} {OUTPUT_DIR}""")

            logger.info(f"Output: {output_filename_r1} {output_filename_r2}")


if __name__ == "__main__":
    main()

