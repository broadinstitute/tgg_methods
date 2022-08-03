import collections

import pymysql


# SQL db schema @ http://useast.ensembl.org/info/docs/api/core/core_schema.html
# To get list of available dbs, connect and run SHOW DATABASES


CURRENT_ENSEMBL_DATABASE = "homo_sapiens_core_107_38"

"""
homo_sapiens_cdna_102_38
homo_sapiens_cdna_103_38
homo_sapiens_core_102_38
homo_sapiens_core_103_38
homo_sapiens_funcgen_102_38
homo_sapiens_funcgen_103_38
homo_sapiens_otherfeatures_102_38
homo_sapiens_otherfeatures_103_38
homo_sapiens_rnaseq_102_38
homo_sapiens_rnaseq_103_38
homo_sapiens_variation_102_38
homo_sapiens_variation_103_38
"""


def _get_gene_id_to_transcript_ids(database=CURRENT_ENSEMBL_DATABASE, only_canonical_transcripts=False):
    """Internal method for retrieving a dictionary containing gene_id => transcript ids

    Args:
        database (str): The Ensembl database name (eg. "homo_sapiens_core_107_38")
        only_canonical_transcripts (bool): If True, only canonical transcripts will be returned

    Return:
        dict: mapping ENSG id string to a list of ENST id strings
    """

    gene_id_to_transcript_id = collections.defaultdict(list)
    with pymysql.connect(host="useastdb.ensembl.org", user="anonymous", database=database) as conn:
        with conn.cursor() as c:
            if only_canonical_transcripts:
                join_keys = "canonical_transcript_id = transcript_id"
            else:
                join_keys = "transcript.gene_id = gene.gene_id"

            columns = ["gene.stable_id", "transcript.stable_id"]
            columns_str = ", ".join(columns)
            c.execute(f"SELECT {columns_str} FROM gene LEFT JOIN transcript ON {join_keys}")

            for row in c:
                d = dict(zip(columns, row))
                gene_id = d['gene.stable_id']
                transcript_id = d['transcript.stable_id']
                gene_id_to_transcript_id[gene_id].append(transcript_id)

    return gene_id_to_transcript_id


def get_all_transcripts(database=CURRENT_ENSEMBL_DATABASE):
    """Returns a dictionary mapping each Ensembl gene_id => transcript ids

    Args:
        database (str): The Ensembl database name (eg. "homo_sapiens_core_107_38")

    Return:
        dict: mapping ENSG id string to a list of ENST id strings
    """

    return _get_gene_id_to_transcript_ids(database=database, only_canonical_transcripts=False)


def get_canonical_transcripts(database=CURRENT_ENSEMBL_DATABASE):
    """Returns a dictionary mapping each Ensembl gene_id => canonical transcript id

    Args:
        database (str): The Ensembl database name (eg. "homo_sapiens_core_107_38")

    Return:
        dict: mapping ENSG id string to the canonical ENST id string
    """

    gene_id_to_transcript_ids = _get_gene_id_to_transcript_ids(database=database, only_canonical_transcripts=True)

    gene_id_to_canonical_transcript_id = {}
    for gene_id, transcript_ids in gene_id_to_transcript_ids.items():
        if len(transcript_ids) > 1:
            raise Exception(f"{gene_id} has more than 1 canonical transcript: " + ", ".join(transcript_ids))
        if len(transcript_ids) == 0:
            raise Exception(f"{gene_id} has 0 canonical transcripts")

        gene_id_to_canonical_transcript_id[gene_id] = transcript_ids[0]

    return gene_id_to_canonical_transcript_id


def get_gene_created_modified_dates(database=CURRENT_ENSEMBL_DATABASE):
    gene_id_to_dates = {}
    with pymysql.connect(host="useastdb.ensembl.org", user="anonymous", database=database) as conn:
        with conn.cursor() as c:
            columns = ["stable_id", "created_date", "modified_date"]
            columns_str = ", ".join(columns)
            c.execute(f"SELECT {columns_str} FROM gene")

            for row in c:
                d = dict(zip(columns, row))
                gene_id = d['stable_id']
                gene_id_to_dates[gene_id] = d['created_date'], d['modified_date']

    return gene_id_to_dates


def get_canonical_transcript_created_modified_dates(database=CURRENT_ENSEMBL_DATABASE):
    transcript_id_to_dates = {}
    with pymysql.connect(host="useastdb.ensembl.org", user="anonymous", database=database) as conn:
        with conn.cursor() as c:
            columns = ["transcript.stable_id", "transcript.created_date", "transcript.modified_date"]
            columns_str = ", ".join(columns)
            c.execute(f"SELECT {columns_str} FROM gene LEFT JOIN transcript ON canonical_transcript_id = transcript_id")

            for row in c:
                d = dict(zip(columns, row))
                gene_id = d['transcript.stable_id']
                transcript_id_to_dates[gene_id] = d['transcript.created_date'], d['transcript.modified_date']

    return transcript_id_to_dates
