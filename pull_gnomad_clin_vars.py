import argparse
import hail as hl
from gnomad_qc.v2.resources.variant_qc import (
    release_ht_path,
    RELEASE_VERSION,
    get_gnomad_public_data,
)
from gnomad_hail.utils.slack import try_slack


CONSEQUENCES = {
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "transcript_ablation",
}

SIG = {
    "Pathogenic",
    "Pathogenic/Likely_pathogenic",
    "Likely_pathogenic",
    "Conflicting_interpretations_of_pathogenicity",
    "DM",
    "DM?",
}

E_POPS = {
    "",  # Used for total
    "nfe_seu",
    "nfe_bgr",
    "nfe_onf",
    "afr",
    "sas",
    "amr",
    "eas",
    "nfe_swe",
    "nfe_nwe",
    "eas_jpn",
    "eas_kor",
    "eas_oea",
    "nfe_est",
    "nfe",
    "fin",
    "asj",
    "oth",
}

G_POPS = {
    "",  # Used for total
    "nfe_seu",
    "afr",
    "nfe_onf",
    "amr",
    "eas",
    "nfe_nwe",
    "nfe_est",
    "nfe",
    "fin",
    "asj",
    "oth",
}


def consequence_filter(ht: hl.Table, csq_terms: set, genes: set) -> hl.Table:
    """
    Explode on transcript consequences and filter variant table down to PTV within gene of interest
    :param ht:
    :param consequences:
    :param genes:
    :return: Table containing only PTVs
    """
    csq_terms = hl.literal(csq_terms)
    genes = hl.literal(genes)
    ht = ht.explode(ht.vep.transcript_consequences)
    ht = ht.explode(ht.vep.transcript_consequences.consequence_terms)
    ht = ht.filter(
        (csq_terms.contains(ht.vep.transcript_consequences.consequence_terms))
        & (genes.contains(ht.vep.transcript_consequences.gene_symbol))
    )
    ht = ht.annotate(
        lof=ht.vep.transcript_consequences.lof,
        csq_term=csq_terms.find(
            lambda x: ht.vep.transcript_consequences.consequence_terms == x
        ),
    )
    return ht


def filter_export_to_gene_list(ht, genes):
    """Take a matrix table and return a table filtered down to a set of genes

    :param Table ht:
    :param list of str or set of str genes: Genes of interest to which to filter table
    :return: Filtered table
    :rtype: Table
    """
    gene_names = hl.literal(genes)
    ht = ht.annotate(
        gene_of_interest=gene_names.find(
            lambda x: ht.vep.transcript_consequences.gene_symbol == x
        )
    )
    ht = ht.filter(
        hl.is_defined(ht.gene_of_interest)
        & hl.is_defined(ht.vep.transcript_consequences.hgvsc)
    )  #  TODO: check on hgvsc utility
    return ht


def filter_variant_summary_to_genes(ht, genes, assembly):
    """Take a matrix table and return a table filtered down to a set of genes
    :param Table ht:
    :param list of str or set of str genes: Genes of interest to which to filter table
    :param String assembly: Genome build
    :return: Filtered table
    :rtype: Table
    """
    gene_names = hl.literal(genes)
    ht = ht.filter((ht.Chromosome == "17") & (ht.Start == "83920497"), keep=False)
    ht = ht.filter(ht.Assembly == assembly)
    ht = ht.annotate(
        gene_of_interest=gene_names.find(lambda x: ht.GeneSymbol.split(";").contains(x))
    )
    ht = ht.filter(hl.is_defined(ht.gene_of_interest))
    return ht


def filter_clinvar_ht_to_sigs(ht: hl.Table, sig: set) -> hl.Table:
    """Filter table of variants down to clinically significant using
    a set of clinical significance terms
    :param Table ht:
    :param list of str or set of str sig: Clinical significance terms
    """
    sig = hl.literal(sig)
    ht = ht.explode(ht.clinvar_clin_sig)
    ht = ht.transmute(clin_sig=sig.contains(ht.clinvar_clin_sig))
    ht = ht.filter(hl.is_defined(ht.clin_sig))
    return ht


def hgmd_path_filter(ht: hl.Table, sig: set) -> hl.Table:
    """Filter table of variants down to clinically significant using
    a set of clinical significance terms
    :param Table ht:
    :param list of str or set of str sig: Clinical significance terms
    """
    sig = hl.literal(sig)
    ht = ht.select(hgmd_sig=sig.find(lambda x: ht.info.CLASS == x), hgmd_rsid=ht.rsid)
    ht = ht.filter(hl.is_defined(ht.hgmd_sig))
    return ht


def make_ucsc_url(ht: hl.Table, window_size: hl.int) -> hl.Table:
    """Create UCSC url and annotate ht with it"""

    build = "19" if ht.locus.dtype.reference_genome == "GRCh37" else "38"
    ucsc_url = f"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg{build}&position=chr"

    return f"{ucsc_url}{ht.vep.seq_region_nam}%3A{hl.str(ht.vep.start - 1)}-{hl.str(ht.vep.end - 1)}&position=chr{ht.vep.seq_region_name}%3A{hl.str(ht.vep.start - (window_size / 2))}-{hl.str(ht.vep.end + (window_size / 2))}"


def make_clinvar_url(ht: hl.Table) -> hl.Table:
    """Create clinvar variant URL and annotate ht with it"""
    # http://www.ncbi.nlm.nih.gov/clinvar?term=rs75822236%5BVariant%20ID%5D
    return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{hl.str(ht.VariationID)}"


def make_hgmd_url(ht: hl.Table) -> hl.Table:
    """Create HGMD pro link and annotate file with it"""
    return f"https://portal.biobase-international.com/hgmd/pro/mut.php?acc={hl.str(ht.hgmd_rsid)}"


def get_datatype_pop(ht, data_type, pop, data_type_pops, stat):
    return (
        ht[f"{data_type}_{gnomad_pop_expr(pop, stat)}"]
        if pop in data_type_pops[data_type]
        else 0
    )


def gnomad_pop_expr(pop, stat) -> hl.str:
    if stat.upper() in {"AC", "AF", "AN"}:
        return f"gnomad_{stat}_adj_{pop}" if pop != "" else f"gnomad_{stat}_adj"
    else:
        ValueError(
            f"{stat} is not an accepted value. Please choose from 'AC', 'AN', or 'AF'"
        )


def gnomad_ac_dict(ht, all_pops) -> hl.dict:
    return {
        f'{gnomad_pop_expr(pop,"AC")}': (
            hl.or_else(get_datatype_pop(ht, "exomes", pop, E_POPS, "AC"), 0)
            + hl.or_else(get_datatype_pop(ht, "genomes", pop, G_POPS, "AC"), 0)
        )
        for pop in all_pops
    }


def gnomad_an_dict(ht, all_pops) -> hl.dict:
    return {
        f'{gnomad_pop_expr(pop,"AN")}': (
            hl.or_else(get_datatype_pop(ht, "genomes", pop, E_POPS, "AN"), 0)
            + hl.or_else(get_datatype_pop(ht, "genomes", pop, G_POPS, "AN"), 0)
        )
        for pop in all_pops
    }


def gnomad_af_dict(ht, all_pops) -> hl.dict:
    return {
        gnomad_pop_expr(pop, "AF"): hl.or_else(
            ht[f'{gnomad_pop_expr(pop,"AC")}'] / ht[f'{gnomad_pop_expr(pop,"AN")}'], 0
        )
        for pop in all_pops
    }


def calc_gnomad_allele_stats(ht: hl.Table, all_pops) -> hl.Table:
    ht = ht.annotate(**gnomad_ac_dict(ht, all_pops), **gnomad_an_dict(ht, all_pops),)
    ht = ht.annotate(**gnomad_af_dict(ht, all_pops))
    return ht


def filter_gnomad_to_genes_consq(data_type, genes, consequences) -> hl.Table:
    ht = get_gnomad_public_data(data_type, split=True)
    ht = hl.filter_intervals(ht, hl.experimental.get_gene_intervals(gene_symbols=genes))
    consequence_ht = consequence_filter(ht, consequences, genes)
    ht = ht.annotate(
        lof=consequence_ht[ht.key].lof, csq_term=consequence_ht[ht.key].csq_term
    )
    return ht


def prep_clinvar_data(genes, sig, output_path) -> hl.Table:

    clinvar_ht = hl.read_table(
        "gs://seqr-reference-data/GRCh37/clinvar/clinvar.GRCh37.ht"
    )
    clinvar_ht = hl.filter_intervals(
        clinvar_ht, hl.experimental.get_gene_intervals(gene_symbols=genes)
    )  # TODO: Add build arg
    clinvar_ht = clinvar_ht.select(
        clinvar_rsid=clinvar_ht.rsid,
        clinvar_info=clinvar_ht.info,
        allele_id=hl.str(clinvar_ht.info.ALLELEID),
        clinvar_clin_sig=clinvar_ht.info.CLNSIG,
    )
    clinvar_ht = filter_clinvar_ht_to_sigs(clinvar_ht, sig)
    clinvar_ht = clinvar_ht.checkpoint(output_path + "clinvar.ht", overwrite=True)
    clinvar_vs = hl.read_table(
        "gs://seqr-datasets/methods_dev/test_data/genzyme/variant_summary_09162019.ht"
    )
    clinvar_vs = clinvar_vs.select(
        "VariationID",
        "Chromosome",
        "Start",
        "Assembly",
        "NumberSubmitters",
        "ClinicalSignificance",
        "LastEvaluated",
        "PhenotypeList",
        "GeneSymbol",
        "RS# (dbSNP)",
    )
    clinvar_vs = clinvar_vs.rename(
        {"#AlleleID": "AlleleID", "RS# (dbSNP)": "clinvar_rsid"}
    )
    clinvar_vs = filter_variant_summary_to_genes(clinvar_vs, genes, "GRCh37")
    clinvar_ht = clinvar_ht.annotate(**clinvar_vs[clinvar_ht.allele_id])  # TODO: Some variants duplicated
    return clinvar_ht


def prep_hgmd_data(genes, sig) -> hl.Table:
    hgmd = hl.read_matrix_table(
        "gs://seqr-datasets/methods_dev/test_data/genzyme/hgmd.mt"
    ).rows()  # TODO: get_hgmd_path
    hgmd = hl.filter_intervals(
        hgmd, hl.experimental.get_gene_intervals(gene_symbols=genes)
    )
    hgmd = hgmd_path_filter(hgmd, sig)
    return hgmd


def annotate_ht_w_all_data(
    ht, clinvar_ht, hgmd_ht, output_path, data_type, genes
) -> hl.Table:
    ht = ht.annotate(**clinvar_ht[ht.key])
    ht = ht.annotate(**hgmd_ht[ht.key])
    ht = ht.filter(
        hl.is_defined(ht.clin_sig)
        | hl.is_defined(ht.hgmd_sig)
        | hl.is_defined(ht.csq_term)
    )
    ht = ht.checkpoint(
        f"{output_path}_gnomad_{data_type}_w_clinvar_hgmd.ht", overwrite=True
    )
    ht = ht.annotate(
        ucsc_url=make_ucsc_url(ht, 50),
        clinvar_url=make_clinvar_url(ht),
        hgmd_url=make_hgmd_url(ht),
    )
    ht = ht.annotate(Reference_Allele=ht.alleles[0], Alternate_Allele=ht.alleles[1])
    variants = ht.count()
    print(f"pre-explode variants: {variants}")  # TODO: switch to logger
    ht = ht.explode(ht.vep.transcript_consequences)
    variants = ht.count()
    print(f"post-explode variants: {variants}")
    ht = filter_export_to_gene_list(ht, genes)
    variants = ht.count()
    print(f"post filter variants: {variants}")

    ht = ht.annotate(gnomad_g_pop_max=ht.popmax[0])
    ht_export = ht.select(
        ht.vep.seq_region_name,
        ht.vep.start,
        ht.vep.end,
        ht.Reference_Allele,
        ht.Alternate_Allele,
        ht.vep.allele_string,
        ht.vep.transcript_consequences.gene_symbol,
        ht.vep.transcript_consequences.canonical,
        ht.vep.transcript_consequences.hgvsc,
        ht.vep.transcript_consequences.hgvsp,
        ht.vep.variant_class,
        ht.VariationID,
        ht.NumberSubmitters,
        ht.ClinicalSignificance,
        ht.hgmd_sig,
        ht.LastEvaluated,
        ht.PhenotypeList,
        ht.clinvar_url,
        ht.hgmd_url,
        ht.ucsc_url,
        ht.csq_term,
        ht.lof,
    )
    ht_export_flat = ht_export.flatten()
    ht_export_flat.write(
        f"{output_path}_gnomad_{data_type}_w_clinvar_hgmd_flat_export.ht",
        overwrite=True,
    )
    ht_export_flat.export(output_path + f"gnomad_{data_type}_w_clinvar_hgmd_export.tsv")
    return ht_export


def get_ac_an(data_type, pops, genes) -> hl.Table:
    """
    Retreives AC and AN for given populations within genes
    :param data_type: exome or genome
    :param pops: population in dataset
    :param genes: genes to filter HT to
    :return:
    """
    ht = hl.read_table(
        release_ht_path(
            f"{data_type}", release_tag=RELEASE_VERSION, nested=False, with_subsets=True
        )
    )
    ht = hl.filter_intervals(ht, hl.experimental.get_gene_intervals(gene_symbols=genes))
    stats = [f"gnomad_AC_adj_{pop}" if pop != "" else "gnomad_AC_adj" for pop in pops]
    stats.extend(
        [f"gnomad_AN_adj_{pop}" if pop != "" else "gnomad_AN_adj" for pop in pops]
    )
    ht = ht.select(*stats)
    rename_dict = {stat: f"{data_type}_{stat}" for stat in stats}
    ht = ht.rename(rename_dict)
    return ht


def main(args):
    genes = args.genes
    output_path = args.output_path
    print(genes)
    all_pops = E_POPS.union(G_POPS)

    clinvar_ht = prep_clinvar_data(genes, SIG, output_path)
    hgmd_ht = prep_hgmd_data(genes, SIG)

    e_ht = filter_gnomad_to_genes_consq("exomes", genes, CONSEQUENCES)
    e_ht = annotate_ht_w_all_data(
        e_ht, clinvar_ht, hgmd_ht, output_path, "exomes", genes
    )

    g_ht = filter_gnomad_to_genes_consq("genomes", genes, CONSEQUENCES)
    g_ht = annotate_ht_w_all_data(
        g_ht, clinvar_ht, hgmd_ht, output_path, "genomes", genes
    )

    ht = e_ht.union(g_ht, unify=True)
    ht = ht.key_by(
        "locus", "alleles", "hgvsc"
    ).distinct()  # TODO upstream combine any variants where multiple sigs
    ht = ht.key_by("locus", "alleles")

    ht = ht.annotate(**get_ac_an("exomes", E_POPS, genes)[ht.key])
    ht = ht.annotate(**get_ac_an("genomes", G_POPS, genes)[ht.key])
    ht = calc_gnomad_allele_stats(ht, all_pops)
    ht = ht.select_globals()
    ht = ht.annotate(
        ENST_hgvsc=hl.cond(
            hl.is_defined(ht.hgvsc), ht.hgvsc.split(":")[0], hl.null(hl.str)
        ),
        hgvsc_split=hl.cond(
            hl.is_defined(ht.hgvsc), ht.hgvsc.split(":")[1], hl.null(hl.str)
        ),
        ENST_hgvsp=hl.cond(
            hl.is_defined(ht.hgvsp), ht.hgvsp.split(":")[0], hl.null(hl.str)
        ),
        hgvsp_split=hl.cond(
            hl.is_defined(ht.hgvsp), ht.hgvsp.split(":")[1], hl.null(hl.str)
        ),
    )
    ht = ht.drop("hgvsp", "hgvsc").rename(
        {"hgvsp_split": "hgvsp", "hgvsc_split": "hgvsc"}
    )  # TODO change to transmute

    ht = ht.select(
        ht.seq_region_name,
        ht.start,
        ht.end,
        ht.Reference_Allele,
        ht.Alternate_Allele,  # TODO supply the string instead of ht.
        ht.allele_string,
        ht.gene_symbol,
        ht.canonical,
        ht.ENST_hgvsc,
        ht.hgvsc,
        ht.ENST_hgvsp,
        ht.hgvsp,
        ht.variant_class,
        ht.VariationID,
        ht.NumberSubmitters,
        ht.ClinicalSignificance,
        ht.hgmd_sig,
        ht.LastEvaluated,
        ht.PhenotypeList,
        ht.clinvar_url,
        ht.hgmd_url,
        ht.ucsc_url,
        ht.csq_term,
        ht.lof,
        ht.gnomad_AF_adj,
        ht.gnomad_AC_adj,
        ht.gnomad_AN_adj,
        ht.gnomad_AF_adj_afr,
        ht.gnomad_AC_adj_afr,
        ht.gnomad_AN_adj_afr,
        ht.gnomad_AF_adj_amr,
        ht.gnomad_AC_adj_amr,
        ht.gnomad_AN_adj_amr,
        ht.gnomad_AF_adj_asj,
        ht.gnomad_AC_adj_asj,
        ht.gnomad_AN_adj_asj,
        ht.gnomad_AF_adj_eas,
        ht.gnomad_AC_adj_eas,
        ht.gnomad_AN_adj_eas,
        ht.gnomad_AF_adj_eas_jpn,
        ht.gnomad_AC_adj_eas_jpn,
        ht.gnomad_AN_adj_eas_jpn,
        ht.gnomad_AF_adj_eas_kor,
        ht.gnomad_AC_adj_eas_kor,
        ht.gnomad_AN_adj_eas_kor,
        ht.gnomad_AF_adj_eas_oea,
        ht.gnomad_AC_adj_eas_oea,
        ht.gnomad_AN_adj_eas_oea,
        ht.gnomad_AF_adj_fin,
        ht.gnomad_AC_adj_fin,
        ht.gnomad_AN_adj_fin,
        ht.gnomad_AF_adj_nfe,
        ht.gnomad_AC_adj_nfe,
        ht.gnomad_AN_adj_nfe,
        ht.gnomad_AF_adj_nfe_bgr,
        ht.gnomad_AC_adj_nfe_bgr,
        ht.gnomad_AN_adj_nfe_bgr,
        ht.gnomad_AF_adj_nfe_est,
        ht.gnomad_AC_adj_nfe_est,
        ht.gnomad_AN_adj_nfe_est,
        ht.gnomad_AF_adj_nfe_nwe,
        ht.gnomad_AC_adj_nfe_nwe,
        ht.gnomad_AN_adj_nfe_nwe,
        ht.gnomad_AF_adj_nfe_onf,
        ht.gnomad_AC_adj_nfe_onf,
        ht.gnomad_AN_adj_nfe_onf,
        ht.gnomad_AF_adj_nfe_seu,
        ht.gnomad_AC_adj_nfe_seu,
        ht.gnomad_AN_adj_nfe_seu,
        ht.gnomad_AF_adj_nfe_swe,
        ht.gnomad_AC_adj_nfe_swe,
        ht.gnomad_AN_adj_nfe_swe,
        ht.gnomad_AF_adj_oth,
        ht.gnomad_AC_adj_oth,
        ht.gnomad_AN_adj_oth,
        ht.gnomad_AF_adj_sas,
        ht.gnomad_AC_adj_sas,
        ht.gnomad_AN_adj_sas,
    )

    ht = ht.checkpoint(
        f"{output_path}{hl.eval(hl.delimit(genes,delimiter='_'))}_final_export_w_freqs.ht",
        overwrite=True,
    )
    ht.export(
        f"{output_path}{hl.eval(hl.delimit(genes,delimiter='_'))}_final_export_w_freqs.tsv"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genes", help="Genes to run pipeline on", nargs="+", required=True
    )  # TODO Add arg to work with list and change output path, date stamp run
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite", help="Overwrite previous paths", action="store_true"
    )
    parser.add_argument(
        "--output-path",
        help="Output path for script results and intermediate files",
        required=True,
    )

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
