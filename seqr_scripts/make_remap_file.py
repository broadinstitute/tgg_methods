import argparse
import csv
import logging
from os import replace

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Make remap file")
logger.setLevel(logging.INFO)


def read_file(file_path: str):
    '''
    Read and parse a two line text file containing:
        - vcf ids separated by a commma and a space (", ") on the first line
        - corresponding seqr ids separated by a comma and a space (", ") on the second line
    Returns:
        python zip (generator) object, with tuples of (vcf id, corresponding seqr id)

    :param file_path: string path of the text file to read.
    '''
    with open(file_path, newline="") as infile:
        reader = csv.reader((line.replace(", ", ",") for line in infile), delimiter=",")
        vcf_ids = next(reader)
        seqr_ids = next(reader)
        if len(vcf_ids)!= len(seqr_ids):
            raise ValueError("Number of vcf ids does not match number of seqr ids.")
    return zip(vcf_ids, seqr_ids)

def main(args):
    with open(args.out, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n", quoting = csv.QUOTE_NONE, quotechar="")
        writer.writerow(["s", "seqr_id"])
        writer.writerows(read_file(args.input))
    logger.info(f"File written to {args.out}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Make project remap tsv from a given text file."
    )
    parser.add_argument(
        "-i"
        "-f",
        "--input",
        required=True,
        help="Input file path.",
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Name for the remap tsv.",
    )
    args = parser.parse_args()
    main(args)