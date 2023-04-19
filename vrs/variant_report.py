# Daniel Marten Variant Report
# Sumbmitted: 04-19-2023

import argparse
import datetime
import errno
import logging
import os
import sys
import json
import pickle

import hail as hl
import hailtop.batch as hb
from gnomad.resources.grch38.gnomad import public_release
from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.resources.annotations import vrs_annotations as v3_vrs_annotations
from tgg.batch.batch_utils import init_job


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("annotate_vrs_ids")
logger.setLevel(logging.INFO)

def main(args):

    hl.init(
        backend="spark"
    )

    ht_path = public_release('genomes').path
    if args.non_ref_path:
        ht_path = args.non_ref_path

    [usr_contig,usr_position,usr_ref,usr_alt] = args.user_allele.split('-')

    ht = hl.read_table(ht_path)

    # NOTE: this is changed from the original, which had hl.filter_intervals(ht,parse_locus_intervals())
    # This came more naturally and was easier, but I wonder if it performs better or worse ? 
    ht_l = ht.filter(
        ht.locus == hl.locus(usr_contig,int(usr_position),reference_genome="GRCh38")
    )

    ht_la = ht_l.filter(
        ht_l.alleles == [usr_ref,usr_alt]
    )

    if args.kc_reread: 
        print('running KCs rereading code!!!')
        # Not sure the function of this but I trust KC's code!!!
        intervals = ht_la._calculate_new_partitions(1)

        # Re-read in HT with new intervals generated above and select vrs struct only
        ht_la = hl.read_table(ht_path, _intervals=intervals)

    ## NEW: CONSTRUCT VRS STRUCT
    ht_2 = ht_la.annotate(
            info = ht_la.info.annotate(
                vrs_struct = hl.struct(
                    VRS_Allele = ht_la.info.VRS_Allele,
                    VRS_Start = ht_la.info.VRS_Start,
                    VRS_End = ht_la.info.VRS_End,
                    VRS_Alt = ht_la.info.VRS_Alt,
                )
            )
        )

    ht_2_vrs = ht_2.select(ht_2.info.vrs_struct)
    output_struct = ht_2_vrs.aggregate(hl.agg.take(ht_2_vrs.vrs_struct, 1)[0]) 

    # Naive & Simple implemention for converting that Struct into a JSON and outputting it
    # Is there supposed to be more refined/better code for this someplace from the GA4GH Team? 
    vrs_dict = {}
    for vrs_key,vrs_value in output_struct.items():
        vrs_dict[vrs_key] = vrs_value

    json_output = json.dumps(vrs_dict)

    print('Desired JSON output as: \n',json_output)

    if args.vrs_txt_out:
        with open(args.vrs_txt_out, 'w', encoding='utf-8') as f:
            json.dump(vrs_dict, f, ensure_ascii=False, indent=4)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--user-allele",
        help="Variant to produce JSON for in format: chr#-position-alt-ref",
        default="chr1-783006-A-G",
        type=str
    )
    parser.add_argument(
        "--non-ref-path",
        help="If arg passed, will read from non reference path provided",
        default=None,
        type=str
    )
    parser.add_argument(
        "--vrs-txt-out",
        help="If supplied, will output JSON as txt File",
        default=None,
        type=str
    )
    parser.add_argument(
        "--kc-reread",
        help="If passed, use KCs rereading code",
        action="store_true"
    )


    args = parser.parse_args()

    main(args)