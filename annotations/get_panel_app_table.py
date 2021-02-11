import pandas as pd
import requests

from annotations.cache_utils import cache_data_table

URL = "https://panelapp.genomicsengland.co.uk/api/v1/genes/"

@cache_data_table
def get_panel_app_table():
    """Download the HGNC table from https://www.genenames.org/download/custom/ and return it as a pandas DataFrame"""
    data = {
        'next': URL
    }
    rows = []
    while data.get('next'):
        url = data['next']

        r = requests.get(url)
        if not r.ok:
            raise Exception(f"Failed to download {url}: {r}")

        data = r.json()

        for r in data["results"]:
            ensembl_genes = r["gene_data"]["ensembl_genes"]

            # gene id may not be specified for some results
            gene_id = ""
            if ensembl_genes and 'GRch38' in ensembl_genes:
                gene_id = next(iter(ensembl_genes['GRch38'].values())).get("ensembl_id")

            rows.append({
                "hgnc": r["gene_data"]["hgnc_id"],
                "gene name": r["gene_data"]["gene_name"],
                "biotype": r["gene_data"]["biotype"],
                "gene id": gene_id,
                "confidence": r["confidence_level"],
                "penetrance": r["penetrance"],
                "mode_of_pathogenicity": r["mode_of_pathogenicity"],
                "mode_of_inheritance": r["mode_of_inheritance"],
                "publications": ", ".join(r["publications"]).replace("\t", "  "),
                "evidence": ", ".join(r["evidence"]).replace("\t", "  "),
                "phenotypes": ", ".join(r["phenotypes"]).replace("\t", "  "),
                "panel name": r["panel"]["name"],
            })

        print(f"Retrieved {url}  total rows: {len(rows)}")

    return pd.DataFrame(rows)


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_panel_app_table()
    print(df)

