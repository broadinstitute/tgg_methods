import pandas as pd
import requests
from io import StringIO

from annotations.cache_utils import cache_data_table


def _get_clingen_table(url):
    """Download one of the ClinGen .csv tables and return it as a pandas DataFrame

    Args:
        url (str): for example "https://search.clinicalgenome.org/kb/gene-validity/download"
    Return:
        pandas DataFrame
    """
    r = requests.get(url)
    if not r.ok:
        raise Exception(f"Failed to download {url}: {r}")

    table_contents = r.content.decode('UTF-8')
    lines = table_contents.split("\n")
    header_line = lines[4]
    #print(header_line)
    table_contents = "\n".join([header_line] + lines[6:])
    return pd.read_csv(StringIO(table_contents))


@cache_data_table
def get_clingen_gene_disease_validity_table():
    "Download ClinGen gene-disease validity table and return it as a pandas DataFrame"
    return _get_clingen_table("https://search.clinicalgenome.org/kb/gene-validity/download")

@cache_data_table
def get_clingen_dosage_sensitivity_table():
    return _get_clingen_table("https://search.clinicalgenome.org/kb/gene-dosage/download")


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_clingen_gene_disease_validity_table()
    print("Gene-Disease Validity Table columns:")
    print(df)

    df = get_clingen_dosage_sensitivity_table()
    print("Dosage Sensitivity Table columns:")
    print(df)
