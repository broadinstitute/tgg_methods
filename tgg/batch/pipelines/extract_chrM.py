import hail as hl   # used for hadoop file utils
import logging
import os
import pandas as pd
from batch import batch_utils
from urllib import parse

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)


GCLOUD_PROJECT = "broad-mpg-gnomad"
GCLOUD_USER_ACCOUNT = "weisburd@broadinstitute.org"
GCLOUD_CREDENTIALS_LOCATION = "gs://weisburd-misc/creds"

DOCKER_IMAGE = "weisburd/gatk@sha256:db9a3ce53be4f734f0c5c35043a51bc0b05ceb6c4e2f4ed59fb4e7b2edf69ee8"
OUTPUT_DIR = "gs://fc-secure-2d4b04e6-f83f-481b-a967-b211aad522c1/v3.1_chrM"


def main():
    p = batch_utils.init_arg_parser(default_cpu=0.5, default_memory=1.75, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    p.add_argument("tsv_path", help="Table with header: sample_id, cram_path, crai_path")
    p.add_argument("sample_id", nargs="*", help="(optional) 1 or more sample_ids to process. If not specified, all rows in the .tsv will be processed.")
    args = p.parse_args()

    df = pd.read_table(args.tsv_path)
    if {"sample_id", "cram_path", "crai_path"} - set(df.columns):
        p.error(f"{args.tsv_path} must contain a 'sample_id', 'cram_path', 'crai_path' columns")

    if not args.force:
        hl.init(log="/dev/null", quiet=True)

    # process samples
    with batch_utils.run_batch(args, batch_name=f"extract chrM") as batch:
        for _, row in df.iterrows():
            if args.sample_id and row.sample_id not in set(args.sample_id):
                continue

            input_filename = os.path.basename(row.cram_path)
            prefix = input_filename.replace(".bam", "").replace(".cram", "")

            output_cram_path = os.path.join(OUTPUT_DIR, f"{prefix}.chrM.cram")
            output_crai_path = os.path.join(OUTPUT_DIR, f"{prefix}.chrM.cram.crai")

            if not args.force and hl.hadoop_is_file(output_cram_path) and hl.hadoop_is_file(output_crai_path):
                logger.info(f"Output files exist (eg. {output_cram_path}). Skipping {input_filename}...")
                continue

            j = batch_utils.init_job(batch, f"chrM: {row.sample_id}", DOCKER_IMAGE if not args.raw else None, args.cpu, args.memory)
            batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)

            # copy inputs
            REF_PATHS = batch_utils.HG38_REF_PATHS
            fasta_filename = os.path.basename(parse.urlparse(REF_PATHS.fasta).path)

            j.command(f"""set -ex
                env
                gsutil -m cp {REF_PATHS.fasta} {REF_PATHS.fai} {REF_PATHS.dict} .
                java -Xms2g -jar /gatk.jar PrintReads \
                    -R {fasta_filename} \
                    -I {row.cram_path} \
                    --read-index {row.crai_path} \
                    -L chrM \
                    --gcs-project-for-requester-pays broad-mpg-gnomad \
                    -O {prefix}.chrM.bam
                        
                samtools view -C -T {fasta_filename} {prefix}.chrM.bam > {prefix}.chrM.cram
                samtools index {prefix}.chrM.cram {prefix}.chrM.cram.crai
                
                gsutil -m cp {prefix}.chrM.cram.crai {output_crai_path}
                gsutil -m cp {prefix}.chrM.cram {output_cram_path}
            """)

            logger.info(f"Submitted {row.sample_id}: {output_cram_path}")


if __name__ == "__main__":
    main()

