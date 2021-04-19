import base64
import datetime
import json
import pandas as pd
import requests

INPUT_TABLE_HEADER = [
    'mim_number',                   # 0
    'phenotype_mim_number',         # 1
    'phenotypic_series_number',     # 2
    'phenotype_inheritance',        # 3

    'locus_size',                   # 4
    'locus',                        # 5
    'locus_hg19',                   # 6

    'cyto',                         # 7
    'gene_symbols',                 # 8
    'gene_id',                      # 9
    'gene_description',             # 10
    'phenotype_description',        # 11
    'date_created',                 # 12
    'date_updated',                 # 13
    'mouse_gene_id',                # 14

    'oe_lof_upper',                 # 15
    'pLI',                          # 16
    'mis_z',                        # 17

    'text',                         # 15
    'comments',                     # 16

    'xstart',
    'xend',
    'xstart_hg19',
    'xend_hg19',
    'phenotype_map_method',
    #'liftover_to_hg19_failed',
]

OUTPUT_COLUMNS = [
    'chrom',
    'start',
    'end',
    'mim_number',
    'phenotype_mim_number',
    'phenotypic_series_number',
    'phenotype_inheritance',
    'gene_symbols',
    'gene_id',
    'gene_description',
    'phenotype_description',
    'date_created',
    'date_updated',
    'mouse_gene_id',
    'oe_lof_upper',
    'pLI',
    'mis_z',
    'text',
    'comments',
]
MAX_GENE_SIZE = 5*10**6  # 5 Mbases (dystrophin is 2.3Mb)

"""
Example row:
$1                      chrom : 1
$2                      start : 7784284
$3                        end : 7845180
$4                 mim_number : 603427
$5       phenotype_mim_number : 616882
$6   phenotypic_series_number :
$7      phenotype_inheritance : Autosomal dominant
$8               gene_symbols : FASPS3, PER3
$9                    gene_id : ENSG00000049246
$10          gene_description : Period circadian regulator 3
$11     phenotype_description : ?Advanced sleep phase syndrome, familial, 3
$12              date_created :
$13              date_updated :
$14             mouse_gene_id : Per3 (MGI:1277134)
$15              oe_lof_upper : 8.4100e-01
$16                       pLI : 1.3661e-17
$17                     mis_z : -1.0855e-01
$18                      text :
$19                  comments : mutation identified in 1 FASPS3 family
"""


def get_omim_table():
    """Retrieves the latest OMIM table"""

    r = requests.get("https://broadinstitute.github.io/omim-search-p/d")
    if not r.ok:
        raise Exception(f"Failed to download latest OMIM json from omim-search-p: {r}")

    json_string = r.json()
    json_obj = json.loads(base64.b64decode(json_string[0]))
    omim_df = pd.DataFrame(columns=INPUT_TABLE_HEADER, data=json_obj["data"])

    #print(omim_df.columns)
    #print(omim_df["locus"])

    omim_df[["chrom", "interval"]] = omim_df["locus"].str.split(":", expand=True)
    omim_df[["start", "end"]] = omim_df["interval"].str.split("-", expand=True).astype("int32")

    omim_df = omim_df[omim_df["locus_size"] < MAX_GENE_SIZE]
    omim_df = omim_df[omim_df["start"] > 1]

    return omim_df[OUTPUT_COLUMNS]


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_omim_table()
    print(df)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
    output_path = f"omim_{timestamp}.tsv"
    df.to_csv(output_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(df)} records to {output_path}")
