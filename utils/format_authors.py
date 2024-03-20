#!/usr/bin/env python
import argparse
import logging

from collections import defaultdict
from subprocess import Popen, PIPE

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)

logger = logging.getLogger("Format gnomAD authors")
logger.setLevel(logging.INFO)

EXPECTED_HEADER = "First Name and Initials\tLast Name\tRole\tEmail\tPrimary Affiliation\tSecondary Affiliation\tTertiary Affiliation\tFunding acknowledgment\tCOI"
# NOTE: At the moment only allow for three affiliations per author


def collect_author_metadata(input_file: str) -> [list, dict]:
    """
    Add all authors to a list with a unique author id, and add author metadata (first name and initial, affiliation, funding, COI) to a dictionary.

    :param input_file: Path to infile_file, which should consist of the author spreadsheet copy and pasted into a txt file to preserve special characters
    :return: List of authors in format of last_name:author_id; dictionary containing metadata for authors with author_id as key
    """
    # Create an author list, which will contain a list of authors formatted as last_name:author_id, where author_id is a unqiue number associated with that author
    author_list = []

    # Create an author dictionary which will contain the metadata for a given author (such as first name and affliation), keyed by the author_id
    authors = defaultdict(dict)
    author_id = 0

    with open(input_file, "r") as f:
        header = next(f)
        if not header.startswith(EXPECTED_HEADER):
            logger.warning("Header does not match expected header")
        for line in f:
            line = line.rstrip()
            items = line.split("\t")

            # For each column, check if there is a value; if not, set the value to NA
            try:
                first_and_initials = items[0]
            except:
                first_and_initials = "NA"

            try:
                last_name = items[1]
            except:
                last_name = "NA"

            try:
                role = items[2]
            except:
                role = "NA"

            try:
                email = items[3]
            except:
                email = "NA"

            try:
                primary_af = items[4]
            except:
                primary_af = "NA"

            try:
                secondary_af = items[5]
            except:
                secondary_af = "NA"

            try:
                tertiary_af = items[6]
            except:
                tertiary_af = "NA"

            try:
                funding = items[7]
            except:
                funding = "NA"

            try:
                coi = items[8]
            except:
                coi = "NA"

            # Add the authors to a list, and their metadata to a dictionary
            author_list.append(f"{last_name}:{author_id}")
            authors["last"][author_id] = last_name
            authors["first_and_initials"][author_id] = first_and_initials
            authors["primary_af"][author_id] = primary_af
            authors["secondary_af"][author_id] = secondary_af
            authors["tertiary_af"][author_id] = tertiary_af
            authors["funding"][author_id] = funding
            authors["coi"][author_id] = coi

            author_id += 1

    return author_list, authors


def format_authors(author_list: list, authors: dict, output_docx: str) -> None:
    """
    Format the authors and corresponding metadata for publication.

    Final output is a word document. Authors are output by first name and initial, followed by last name, followed by affiliation numbers in superscript. The list of affiliations are output below the authors.

    :param author_list: List of authors in format of last_name:author_id
    :param authors: Dictionary of author metadata (keyed first by 'first_and_initials', 'primary_af', 'secondary_af', 'tertiary_af', 'funding', or 'coi' and then by author_id)
    :param output_docx: Path to docx file where output should be written
    :return: None
    """
    # Keep track of seen affilations so that the same affilations will be assigned to the same number across all authors
    seen_affliations = {}
    affliation_counter = 0
    formatted_author_list = []
    funding_list = []
    coi_list = []

    # Output information to a markdown file that will later be converted to a word doc to allow formatting such as superscripts
    outfile = open("formatted_authors.md", "w")
    # Sort the author list, which means authors will be sorted alphabetically by last name
    author_list = sorted(author_list)
    for i in author_list:
        last_name, author_id = i.split(":")
        author_id = int(author_id)
        first_and_initials = authors["first_and_initials"][author_id]
        primary_af = authors["primary_af"][author_id]
        secondary_af = authors["secondary_af"][author_id]
        tertiary_af = authors["tertiary_af"][author_id]
        funding = authors["funding"][author_id].rstrip().rstrip("\.")
        coi = authors["coi"][author_id]

        # Add all of the given author's affilaitons to a list
        affs = [primary_af, secondary_af, tertiary_af]
        author_aff_nos = []

        for affiliation in affs:
            affiliation = affiliation.lstrip().rstrip().rstrip("\.")

            # Only proceed if the affiliation isn't NA or blank
            if affiliation != "NA" and affiliation != "":
                # If the affiliation has not yet been seen, increment the counter and assign the affiliation to that number, otherwise grab the number from the seen_affiliations dictionary
                if affiliation not in seen_affliations:
                    affliation_counter += 1
                    affliation_no = affliation_counter
                    seen_affliations[affiliation] = affliation_no
                else:
                    affliation_no = seen_affliations[affiliation]

                author_aff_nos.append(affliation_no)

        # Add authors funding information to a funding list
        if funding != "NA" and funding != "":
            funding_list.append(f"{first_and_initials} {last_name}: {funding}")
        # Add authors conflicts of interest to a coi list
        if coi != "NA" and coi != "" and coi != "No COIs" and coi != "None":
            coi_list.append(f"{first_and_initials} {last_name}: {coi}")

        # Sort the affilations numerically for a given author and join them together by comma
        author_aff_nos = map(str, sorted(author_aff_nos))
        author_aff_nos = ",".join(author_aff_nos)

        # Output the authors first name and initial, followed by last name, followed by affiliation numbers in superscript
        formatted_author_list.append(
            f"{first_and_initials} {last_name}^{author_aff_nos}^"
        )

    outfile.write(", ".join(formatted_author_list))
    outfile.write("\n\n")

    # Sort affiliations by value (number) so that they will be written out in the correct order below the list of authors
    sorted_affiliations = {
        k: v for k, v in sorted(seen_affliations.items(), key=lambda item: item[1])
    }
    for k, v in sorted_affiliations.items():
        outfile.write(f"^{v}^{k}\n")

    outfile.write("\n\n")

    # Output funding information
    if len(funding_list) > 0:
        outfile.write("Authors received funding as follows:\n")
        outfile.write("\n".join(funding_list))
    outfile.write("\n\n")

    # Output conflict of interest information
    if len(coi_list) > 0:
        outfile.write("Conflicts of interest are as follow:\n")
        outfile.write("\n".join(coi_list))
    else:
        outfile.write("No conflicts of interest to declare")

    outfile.close()

    # Convert the markdown file to a word doc
    cmd = f"""pandoc -o {output_docx} -f markdown+hard_line_breaks  -t docx formatted_authors.md"""
    result, err = Popen([cmd], stdout=PIPE, stderr=PIPE, shell=True).communicate()
    logger.info(result)
    logger.info(err)

    # Output the total count of authors
    logger.info("A total of %d authors were output", len(author_list))


def main(args):
    # Collect the author metadata
    author_list, authors = collect_author_metadata(args.input_file)
    # Ouput the authors and the metadata in format suitable for publication
    format_authors(author_list, authors, args.output_docx)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script formats and reorders the gnomAD author list in preparation for publication"
    )
    parser.add_argument(
        "-i",
        "--input-file",
        help="Path to infile_file, which should consist of the author spreadsheet copy and pasted into a txt file",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output-docx",
        help="Path to docx file where output should be written",
        required=True,
    )

    args = parser.parse_args()

    main(args)
