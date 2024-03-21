import argparse
import csv
import logging
from os import replace

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Make remap file")
logger.setLevel(logging.INFO)


def read_file(file_path: str):
    """
    Read and parse a two line text file containing:
        - vcf ids separated by a commma and a space (", ") on the first line
        - corresponding seqr ids separated by a comma and a space (", ") on the second line

    This format for the text file was considered given how vcf ids and seqr ids are provided in airtable; it's easiest to copy each group into a single line, instead of a column for each group. Makes formatting seqr remap files much easier.
    :param file_path: string path of the text file to read.
    :return: zip (generator) object, with tuples of (vcf id, corresponding seqr id)
    """
    with open(file_path, newline="") as infile:
        #csv doesn't parse multi-character separators
        reader = csv.reader((line.replace(", ", ",") for line in infile), delimiter=",")
        vcf_ids = next(reader)
        seqr_ids = next(reader)
        if len(vcf_ids) != len(seqr_ids):
            raise ValueError("Number of vcf ids does not match number of seqr ids.")
    return zip(vcf_ids, seqr_ids)


def main(args):
    with open(args.out, "w", newline="") as outfile:
        writer = csv.writer(
            outfile,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
            quotechar="",
        )
        # write the header
        writer.writerow(["s", "seqr_id"])
        # get rows from generator, and write each row as a line
        writer.writerows(read_file(args.input))
    logger.info(f"File written to {args.out}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Make seqr project remap tsv from a given text file."
    )
    parser.add_argument(
        "-i" "-f",
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
