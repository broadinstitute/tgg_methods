
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


def main():
    p = batch_utils.init_arg_parser(default_cpu=1, gsa_key_file=os.path.expanduser("~/.config/gcloud/misc-270914-cb9992ec9b25.json"))
    args = p.parse_args()

    # process samples
    with batch_utils.run_batch(args, "test") as batch:
        for cpu in (0.25, 0.5, 1, 2):
            args.cpu = cpu
            j = batch.new_job(f"test - {args.cpu} cpu")
            j.image(DOCKER_IMAGE)
            j.cpu(args.cpu)
            j.memory(args.cpu*3.75)
            #batch_utils.switch_gcloud_auth_to_user_account(j, GCLOUD_CREDENTIALS_LOCATION, GCLOUD_USER_ACCOUNT, GCLOUD_PROJECT)
            #j.command(f"yes > data.txt || true")
            j.command(f"ls -lh")
            j.command(f"df -kh")
            j.command(f"sleep 3600") # sleep for 0.5 hour
            j.command(f"free -h")

if __name__ == "__main__":
    main()
