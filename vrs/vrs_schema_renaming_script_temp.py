"""
 Script to update column names from GA4GH VRS-Python Branch 0.8.1 to 0.8.2 convention
 See: https://github.com/ga4gh/vrs-python/releases/tag/0.8.2
 Note that this script will be run on 05-04-23 or 05-05-23 as a temporary fix.
These changes will be permenantly added to the script which adds the VRS Annotations in a later update.
 Note that this code may not be fully productionized or up-to-date with the rest of the team's code standards. 
"""

import argparse
import hail as hl


def main(args):

    import_path = args.import_ht
    ht = hl.read_table(import_path)

    # Code which actually changes naming schema from VRS-Python Branch 0.8.1 to 0.8.2

    ht = ht.annotate(
        info=ht.info.annotate(
            vrs=ht.info.vrs.rename(
                {
                    "VRS_Allele": "VRS_Allele_IDs",
                    "VRS_Start": "VRS_Starts",
                    "VRS_End": "VRS_Ends",
                    "VRS_Alt": "VRS_States",
                }
            )
        )
    )

    ht_reformat = ht.annotate(
        info=ht.info.annotate(
            vrs=ht.info.vrs.annotate(
                VRS_Allele_IDs=ht.info.vrs.VRS_Allele_IDs.split(","),
                VRS_Starts=ht.info.vrs.VRS_Starts.split(",").map(lambda x: hl.int(x)),
                VRS_Ends=ht.info.vrs.VRS_Ends.split(",").map(lambda x: hl.int(x)),
                VRS_States=ht.info.vrs.VRS_States.split(","),
            )
        )
    )

    ht_reformat.info.vrs.show()

    # Code to add global annotation for which VRS Version (not GitHub Branch Number) is used

    ht = ht.annotate_globals(VRS_Version="v1.3")

    # Outputting updated table with path appended to include and state the update

    ht_output_path = args.export_ht
    print("Outputting to: ", ht_output_path)
    ht_reformat.write(ht_output_path, overwrite=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--import-ht",
        help="Hail Table path to reformat according to GA4GH VRS-Python Branch 0.8.2.",
        type=str,
    )
    parser.add_argument(
        "--export-ht", help="Hail Table path to export reformatted HT to.", type=str
    )

    args = parser.parse_args()

    main(args)
