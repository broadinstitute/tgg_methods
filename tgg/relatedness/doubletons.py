"""
This script tests the utility of finding doubleton pairs in relatedness inference.

The current method of relatedness inference using pc-relate will be too inefficient
on large datasets. Konrad proposed an analysis to see whether filtering to pairs of
samples that share a rare doubleton could help reduce the search space for relatedness inference.
"""
import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.utils.slack import slack_notifications

from tgg.relatedness.doubleton_utils import compare_doubletons_to_related


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("test_doubletons_relatedness")
logger.setLevel(logging.INFO)


def main(args):
    """Find doubleton pairs and compare to related pairs."""
    try:
        hl.init(log="/test_doubletons_relatedness.log", default_reference="GRCh38")
        compare_doubletons_to_related()

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(args.temp_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        """
        This script extracts doubletons from the 455k VDS and
        compares them to related sample pairs from the 455k callset.
        """
    )
    parser.add_argument(
        "--slack-channel",
        help="Send message to Slack channel/user.",
    )
    parser.add_argument(
        "--slack-token",
        help="Token to authenticate slack. Must be specified if --slack-channel is also set.",
    )
    parser.add_argument(
        "--temp-path",
        help="Path to temporary bucket to store hail logs.",
    )
    args = parser.parse_args()

    if args.slack_channel:
        if not args.slack_token:
            raise DataException("Must specify slack_token if --slack-channel is set.")
        with slack_notifications(args.slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
