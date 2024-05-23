import argparse
import datetime
from copy import deepcopy
import logging

import hail as hl

from gnomad.resources.grch38.gnomad import POPS
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_prevalence_calculator")
logger.setLevel(logging.INFO)

LOF_CONSEQUENCES = {
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "transcript_ablation",
}

MISSENSE_CONSEQUENCES = {
    "missense_variant",
}

SIG = {
    "Pathogenic",
    "Pathogenic/Likely_pathogenic",
    "Pathogenic/Likely_pathogenic/Established_risk_allele",
    "Pathogenic/Likely_pathogenic/Likely_risk_allele",
    "Pathogenic/Likely_risk_allele",
    "Likely_pathogenic",
    "Likely_pathogenic/Likely_risk_allele",
    "Conflicting_interpretations_of_pathogenicity",
    "Conflicting_classifications_of_pathogenicity",
    "DM",
}

POPS = {
    "exomes": deepcopy(POPS["v4"]["exomes"]),
    "genomes": deepcopy(POPS["v4"]["genomes"]),
}
ALL_POPS = set().union(*POPS.values())
# Insert "" to account for total dataset
ALL_POPS.add("")

TEMP_PATH = "gs://gnomad-tmp-4day/mwilson/prevalence/"
RESULT_PATH = "gs://gnomad-mwilson/prevalence/results/"
CLINVAR_HT_PATH = (
    "gs://gnomad-mwilson/prevalence/clinvar/GRCh38/clinvar/clinvar.GRCh38.ht"
)
CLINVAR_VAR_SUMMARY_PATH = (
    "gs://gnomad-mwilson/prevalence/clinvar_variant_summary/variant_summary.ht"
)
HGMD_HT_PATH = "gs://gnomad-mwilson/prevalence/hgmd/HGMD_Pro_2023.1_hg38.ht"


def create_lof_consequence_ht(
    ht: hl.Table, csq_terms: set, gene: hl.str, keep_missense=True
) -> hl.Table:
    """
    Explode on transcript consequences and filter variant table down to PTV within gene of interest
    :param ht:
    :param consequences:
    :param genes:
    :return: Table containing only PTVs
    """
    csq_terms = hl.literal(csq_terms)
    ht = ht.explode(ht.vep.transcript_consequences.consequence_terms)
    ht = ht.filter(
        (csq_terms.contains(ht.vep.transcript_consequences.consequence_terms))
        & (
            (ht.vep.transcript_consequences.gene_symbol == gene)
            & (
                (
                    ht.vep.transcript_consequences.canonical == 1
                )  # TODO: Confirm this behavior is what we want -- should this be futher filtered to where mane select exist, if it doesnt, default to canonical?
                | hl.is_defined(ht.vep.transcript_consequences.mane_select)
            )
        )
    )
    ht = ht.annotate(
        lof=ht.vep.transcript_consequences.lof,
        csq_term=csq_terms.find(
            lambda x: ht.vep.transcript_consequences.consequence_terms == x
        ),
    )
    ht = ht.filter(ht.lof == "HC")
    return ht


def create_missense_revel_ht(
    ht: hl.Table, missense_csq_terms: set, gene: hl.str
) -> hl.Table:
    """
    Explode on transcript consequences and filter variant table down to PTV within gene of interest
    :param ht:
    :param csq_terms:
    :param genes:
    :return: Table containing missense variants with high revel scores
    """
    missense_csq_terms = hl.literal(missense_csq_terms)
    ht = ht.explode(ht.vep.transcript_consequences.consequence_terms)
    ht = ht.select(
        missense_high_revel=(
            missense_csq_terms.contains(
                ht.vep.transcript_consequences.consequence_terms
            )
        )
        & (
            (ht.vep.transcript_consequences.gene_symbol == gene)
            & (
                (
                    ht.vep.transcript_consequences.canonical == 1
                )  # TODO: Confirm this behavior is what we want -- should this be futher filtered to where mane select exist, if it doesnt, default to canonical?
                | hl.is_defined(ht.vep.transcript_consequences.mane_select)
            )
        )
        & (ht.in_silico_predictors.revel_max >= 0.932)
    )
    ht = ht.filter(ht.missense_high_revel)
    return ht


def filter_variant_summary_to_genes(ht: hl.Table, gene: str, assembly: str):
    """
    Take a matrix table and return a table filtered down to a set of genes

    :param Table ht:
    :param list of str or set of str genes: Genes of interest to which to filter table
    :param String assembly: Genome build
    :return: Filtered table
    :rtype: Table
    """
    gene_name = hl.literal(gene)
    ht = ht.filter(
        (ht.Chromosome == "17") & (ht.Start == "83920497"), keep=False
    )  # Weird locus in clinvar table
    ht = ht.filter(
        (ht.Assembly == assembly) & (ht.GeneSymbol.split(";").contains(gene_name))
    )
    return ht


def filter_clinvar_ht_to_sigs(ht: hl.Table, sig: set = SIG) -> hl.Table:
    """
    Filter table of variants down to clinically significant using a set of clinical significance terms.

    :param Table ht: Clinvar Table
    :param list of str or set of str sig: Clinical significance terms
    """
    sig = hl.literal(sig)
    ht = ht.explode(ht.clinvar_clin_sig)
    ht = ht.filter(sig.contains(ht.clinvar_clin_sig))
    ht = ht.filter(
        hl.if_else(
            hl.is_defined(ht.clinvar_clin_conf),
            (
                hl.any(
                    lambda x: x.startswith("Pathogenic")
                    | x.startswith("Likely_pathogenic"),
                    ht.clinvar_clin_conf,
                )
            ),
            True,
        )
    )
    return ht


def make_ucsc_url(ht: hl.Table) -> hl.Table:
    """Create UCSC url and annotate ht with it"""
    # https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&highlight=hg38.chr6%3A145686297-145686297&position=chr6%3A145686272-145686322
    ucsc_url = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&highlight=hg38."

    return hl.format(
        "%s%s%s%s%s%s%s%s%s%s%s%s",
        ucsc_url,
        ht.vep.seq_region_name,
        "%3A",
        hl.str(ht.vep.start - 1),
        "-",
        hl.str(ht.vep.end - 1),
        "&position=",
        ht.vep.seq_region_name,
        "%3A",
        hl.str(ht.vep.start - 25),
        "-",
        hl.str(ht.vep.end + 25),
    )


def make_clinvar_url(ht: hl.Table) -> hl.Table:
    """Create clinvar variant URL and annotate ht with it"""
    # http://www.ncbi.nlm.nih.gov/clinvar?term=rs75822236%5BVariant%20ID%5D
    return hl.format(
        "%s%s", "https://www.ncbi.nlm.nih.gov/clinvar/variation/", ht.VariationID
    )


def make_hgmd_url(ht: hl.Table) -> hl.Table:
    """Create HGMD pro link and annotate file with it"""
    return hl.if_else(
        hl.is_defined(ht.hgmd_rsid),
        hl.format(
            "%s%s",
            "https://portal.biobase-international.com/hgmd/pro/mut.php?acc=",
            ht.hgmd_rsid,
        ),
        "NA",
    )


def get_datatype_pop(ht, data_type, pop, data_type_pops, stat):
    """"""
    return (
        ht[f"{data_type}_{gnomad_pop_expr(pop, stat)}_adj"]
        if pop in data_type_pops
        else 0
    )


def gnomad_pop_expr(pop, stat) -> hl.str:
    if stat.upper() in {"AC", "AF", "AN"}:
        return f"{stat}_{pop}" if pop != "" else f"{stat}"
    else:
        ValueError(
            f"{stat} is not an accepted value. Please choose from 'AC', 'AN', or 'AF'"
        )


def gnomad_ac_dict(ht, all_pops=ALL_POPS) -> hl.dict:
    logger.info("Building AC Dict")
    return {
        f'{gnomad_pop_expr(pop,"AC")}': (
            hl.or_else(get_datatype_pop(ht, "exomes", pop, POPS["exomes"], "AC"), 0)
            + hl.or_else(get_datatype_pop(ht, "genomes", pop, POPS["genomes"], "AC"), 0)
        )
        for pop in all_pops
    }


def gnomad_an_dict(ht, all_pops=ALL_POPS) -> hl.dict:
    logger.info("Building AN dict")
    return {
        f'{gnomad_pop_expr(pop,"AN")}': (
            hl.or_else(get_datatype_pop(ht, "exomes", pop, POPS["exomes"], "AN"), 0)
            + hl.or_else(get_datatype_pop(ht, "genomes", pop, POPS["genomes"], "AN"), 0)
        )
        for pop in all_pops
    }


def gnomad_af_dict(ht, all_pops=ALL_POPS) -> hl.dict:
    logger.info("Building AF dict")
    return {
        gnomad_pop_expr(pop, "AF"): hl.if_else(
            ht[f'{gnomad_pop_expr(pop,"AN")}'] == 0,
            0,
            ht[f'{gnomad_pop_expr(pop,"AC")}'] / ht[f'{gnomad_pop_expr(pop,"AN")}'],
        )
        for pop in all_pops
    }


def calc_gnomad_allele_stats(ht: hl.Table, all_pops=ALL_POPS) -> hl.Table:
    """
    Calculate combined AC, AN, AF in gnomAD exomes and genomes
    :param ht:
    :param all_pops:
    :return:
    """
    ht = ht.annotate(
        **gnomad_ac_dict(ht, all_pops),
        **gnomad_an_dict(ht, all_pops),
    )
    ht = ht.annotate(**gnomad_af_dict(ht, all_pops))
    return ht


def filter_to_gene_annotate_consq(
    data_type: str,
    gene: str,
    gene_interval,
    lof_csq_terms: list,
    missense_csq_terms: list,
) -> hl.Table:
    """
    Filter to genes of interest and annotate with consequence terms.

    :param data_type: exome or genome
    :param gene: Gene to filter to
    :param gene_interval: Gene interval to filter to
    :param consequences: List of consequence terms to filter to
    :return: Filtered Table
    """
    ht = release_sites(data_type).ht()
    ht = hl.filter_intervals(ht, gene_interval)
    gene = hl.literal(gene)
    ht = ht.explode(ht.vep.transcript_consequences)
    # Filter to gene and Ensembl transcripts, v4 has both ensemble and refseq so rows are duplicated but we only need one
    ht = ht.filter(
        (ht.vep.transcript_consequences.gene_symbol == gene)
        & (ht.vep.transcript_consequences.source == "Ensembl")
    )

    # Filter filtered variants
    ht = ht.filter(
        ht.filters.contains("AC0") | ht.filters.contains("AS_VQSR"), keep=False
    )

    # Create HT with just LOF consequences within gene of interest
    consequence_ht = create_lof_consequence_ht(ht, lof_csq_terms, gene)
    ht = ht.annotate(
        lof=consequence_ht[ht.key].lof, csq_term=consequence_ht[ht.key].csq_term
    )
    # Create HT with just missense consequences with high revel score within gene of interest
    missense_ht = create_missense_revel_ht(ht, missense_csq_terms, gene)
    ht = ht.annotate(**missense_ht[ht.key])
    return ht


def prep_clinvar_data(gene, gene_interval, sig=SIG) -> hl.Table:
    """
    Filter clinvar table to genes of interest and annotate with clinvar URL.

    :param gene: Gene to filter to
    :param gene_interval: Gene interval to filter to
    :param sig: List of significance terms to filter to
    :return: Filtered clinvar Table with clinvar_url annotation
    """
    clinvar_ht = hl.read_table(CLINVAR_HT_PATH)
    # Filter to gene of interest
    clinvar_ht = hl.filter_intervals(clinvar_ht, gene_interval)
    clinvar_ht = clinvar_ht.select(
        clinvar_rsid=clinvar_ht.rsid,
        allele_id=hl.str(clinvar_ht.info.ALLELEID),
        clinvar_clin_sig=clinvar_ht.info.CLNSIG,
        clinvar_clin_conf=clinvar_ht.info.CLNSIGCONF,
    )
    # Filter to significant conseqeuence terms
    clinvar_ht = filter_clinvar_ht_to_sigs(clinvar_ht, sig)
    clinvar_ht = clinvar_ht.checkpoint(TEMP_PATH + "clinvar.ht", overwrite=True)

    # Filter clinvars variant summary to genes of interest
    clinvar_vs = hl.read_table(CLINVAR_VAR_SUMMARY_PATH).key_by(
        "#AlleleID"
    )  # Update to GRCh38
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
    clinvar_vs = filter_variant_summary_to_genes(clinvar_vs, gene, "GRCh38")
    clinvar_ht = clinvar_ht.annotate(
        **clinvar_vs[clinvar_ht.allele_id]
    )  # TODO: Some variants duplicated
    clinvar_ht = clinvar_ht.annotate(clinvar_url=make_clinvar_url(clinvar_ht))
    return clinvar_ht


def prep_hgmd_data(gene_interval, sig) -> hl.Table:
    """
    Filter HGMD table to genes of interest and annotate with HGMD URL.

    :param gene: Gene to filter to
    :param gene_interval: Gene interval to filter to
    :param sig: List of significance terms to filter to
    :return: Filtered HGMD Table with hgmd_url annotation
    """
    hgmd = hl.read_table(HGMD_HT_PATH)
    hgmd = hl.filter_intervals(hgmd, gene_interval)
    # Filter table of variants down to clinically significant using
    # a set of clinical significance terms
    sig = hl.literal(sig)
    hgmd = hgmd.select(
        hgmd_sig=sig.find(lambda x: hgmd.info.CLASS == x), hgmd_rsid=hgmd.rsid
    )
    hgmd = hgmd.filter(hl.is_defined(hgmd.hgmd_sig))

    hgmd = hgmd.annotate(hgmd_url=make_hgmd_url(hgmd))
    return hgmd


def annotate_ht_w_all_data(ht, clinvar_ht, hgmd_ht, data_type) -> hl.Table:
    ht = ht.annotate(**clinvar_ht[ht.key])
    ht = ht.annotate(**hgmd_ht[ht.key])

    ht = ht.checkpoint(
        f"{TEMP_PATH}gnomad_{data_type}_w_clinvar_hgmd.ht", overwrite=True
    )
    ht = ht.annotate(
        ucsc_url=make_ucsc_url(ht),
    )
    ht = ht.annotate(Reference_Allele=ht.alleles[0], Alternate_Allele=ht.alleles[1])

    POPS[data_type].insert(0, "")
    data_type_ac_ans = get_ac_an(ht, data_type, POPS[data_type])

    ht_export = ht.select(
        ht.filters,
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
        ht.missense_high_revel,
        mane_select=ht.vep.transcript_consequences.mane_select,
        revel_max=ht.in_silico_predictors.revel_max,
        **data_type_ac_ans[0],
        **data_type_ac_ans[1],
    )
    return ht_export


def get_ac_an(ht, data_type, pops) -> hl.Table:
    """
    Retreives AC and AN for given populations within genes
    :param ht: gnomAD HT
    :param data_type: exome or genome
    :param pops: population in dataset
    :return:
    """
    pop_indices = [f"{pop}_adj" if pop != "" else "adj" for pop in pops]
    AC_callstat_dict = {
        f"{data_type}_AC_{idx}": ht.freq[ht.freq_index_dict[idx]].AC
        for idx in pop_indices
    }
    AN_callstat_dict = {
        f"{data_type}_AN_{idx}": ht.freq[ht.freq_index_dict[idx]].AN
        for idx in pop_indices
    }
    return (AC_callstat_dict, AN_callstat_dict)


def main(args):
    hl.init(
        default_reference="GRCh38",
        log="gnomad_prevalence.log",
        tmp_dir="gs://gnomad-tmp-4day/mwilson/prevalence/tmp",
    )
    # Open genes file and iterate through each gene printing itå
    with hl.hadoop_open(args.genes, "r") as f:
        for gene in f:
            gene = gene.strip()
            gene_interval = hl.experimental.get_gene_intervals(
                gene_symbols=[gene], reference_genome="GRCh38"
            )
            logger.info("Running the gnomaAD prevalence script on %s...", gene)
            output_path = (
                RESULT_PATH + f"{datetime.datetime.now().strftime('%Y-%m-%d')}/{gene}/"
            )
            clinvar_ht = prep_clinvar_data(gene, gene_interval, SIG)
            hgmd_ht = prep_hgmd_data(gene_interval, SIG)

            e_ht = filter_to_gene_annotate_consq(
                "exomes",
                gene,
                gene_interval,
                lof_csq_terms=LOF_CONSEQUENCES,
                missense_csq_terms=MISSENSE_CONSEQUENCES,
            )
            e_ht = annotate_ht_w_all_data(e_ht, clinvar_ht, hgmd_ht, "exomes")

            g_ht = filter_to_gene_annotate_consq(
                "genomes",
                gene,
                gene_interval,
                lof_csq_terms=LOF_CONSEQUENCES,
                missense_csq_terms=MISSENSE_CONSEQUENCES,
            )
            g_ht = annotate_ht_w_all_data(g_ht, clinvar_ht, hgmd_ht, "genomes")
            ht = e_ht.union(g_ht, unify=True)
            ht = ht.key_by(
                "locus", "alleles", "hgvsc"
            ).distinct()  # TODO upstream combine any variants where multiple sigs
            ht = ht.key_by("locus", "alleles")

            ht = calc_gnomad_allele_stats(ht)
            ht = ht.select_globals()
            ht = ht.annotate(
                ENST_hgvsc=hl.if_else(
                    hl.is_defined(ht.hgvsc), ht.hgvsc.split(":")[0], hl.missing(hl.tstr)
                ),
                hgvsc_split=hl.if_else(
                    hl.is_defined(ht.hgvsc), ht.hgvsc.split(":")[1], hl.missing(hl.tstr)
                ),
                ENST_hgvsp=hl.if_else(
                    hl.is_defined(ht.hgvsp), ht.hgvsp.split(":")[0], hl.missing(hl.tstr)
                ),
                hgvsp_split=hl.if_else(
                    hl.is_defined(ht.hgvsp), ht.hgvsp.split(":")[1], hl.missing(hl.tstr)
                ),
            )
            ht = ht.drop("hgvsp", "hgvsc").rename(
                {"hgvsp_split": "hgvsp", "hgvsc_split": "hgvsc"}
            )
            ht = ht.annotate(
                gnomAD_ID=ht.locus.contig.replace("^chr", "")
                + "-"
                + hl.str(ht.locus.position)
                + "-"
                + ht.alleles[0]
                + "-"
                + ht.alleles[1],
                source=hl.if_else(
                    (
                        hl.is_defined(ht.ClinicalSignificance)
                        | hl.is_defined(ht.hgmd_sig)
                    )
                    & hl.is_defined(ht.lof),
                    "Both LoF + ACMG",
                    hl.if_else(
                        (
                            hl.is_defined(ht.ClinicalSignificance)
                            | hl.is_defined(ht.hgmd_sig)
                        )
                        & hl.is_defined(ht.missense_high_revel),
                        "Both revel + ACMG",
                        hl.if_else(
                            (
                                hl.is_defined(ht.ClinicalSignificance)
                                | hl.is_defined(ht.hgmd_sig)
                            )
                            & hl.is_missing(ht.lof)
                            & hl.is_missing(ht.missense_high_revel),
                            "ACMG",
                            hl.if_else(
                                hl.is_defined(ht.lof),
                                "LoF",
                                hl.if_else(
                                    hl.is_defined(ht.missense_high_revel),
                                    "high revel",
                                    hl.missing(hl.tstr),
                                ),
                            ),
                        ),
                    ),
                ),
                build=str(ht.locus.dtype.reference_genome),
                VCI_link="",
                ACMG_curations_verdict="",
                Second_review="",
                Notes="",
                LoF_curation_verdict="",
                Final_classification="",
                Reason_for_removal="",
            )

            ht = ht.filter(
                hl.is_defined(ht.ClinicalSignificance)
                | hl.is_defined(ht.hgmd_sig)
                | hl.is_defined(ht.csq_term)
                | hl.is_defined(
                    ht.lof
                )  # TODO: Confirm this is wanted -- just added 12/4, was not there before v4 udpates
                | hl.is_defined(ht.missense_high_revel)
            )
            ht = ht.filter(
                ht.filters.contains("AC0") | ht.filters.contains("RF"), keep=False
            )

            logger.info("Selecting export fields")
            ht = ht.select(
                ht.gnomAD_ID,
                ht.source,
                ht.build,
                ht.seq_region_name,
                ht.start,
                ht.end,
                ht.Reference_Allele,
                ht.Alternate_Allele,
                ht.allele_string,
                ht.gene_symbol,
                ht.canonical,
                ht.mane_select,
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
                ht.missense_high_revel,
                ht.revel_max,
                ht.VCI_link,
                ht.ACMG_curations_verdict,
                ht.Second_review,
                ht.Notes,
                ht.LoF_curation_verdict,
                ht.Final_classification,
                ht.Reason_for_removal,
                ht.AF,
                ht.AC,
                ht.AN,
                ht.AF_afr,
                ht.AC_afr,
                ht.AN_afr,
                ht.AF_amr,
                ht.AC_amr,
                ht.AN_amr,
                ht.AF_asj,
                ht.AC_asj,
                ht.AN_asj,
                ht.AF_eas,
                ht.AC_eas,
                ht.AN_eas,
                ht.AF_fin,
                ht.AC_fin,
                ht.AN_fin,
                ht.AF_nfe,
                ht.AC_nfe,
                ht.AN_nfe,
                ht.AF_remaining,
                ht.AC_remaining,
                ht.AN_remaining,
                ht.AF_sas,
                ht.AC_sas,
                ht.AN_sas,
            )
            if args.all_transcripts:
                ht = ht.checkpoint(
                    f"{output_path}{hl.eval(hl.delimit(gene,delimiter='_'))}_final_export_w_freqs_all.ht",
                    overwrite=True,
                )
                ht.export(
                    f"{output_path}{hl.eval(hl.delimit(gene,delimiter='_'))}_all_variants_final_export_w_freqs.tsv"
                )

            # Filter to mane select transcripts when it is defined else filter to canonical
            # TODO: Confirm this works as anticipated
            ht = ht.filter(
                hl.if_else(
                    hl.is_defined(ht.mane_select),
                    hl.is_defined(ht.mane_select),
                    ht.canonical == 1,
                )
            )
            ht.export(f"{output_path}{gene}_final_export_w_freqs.tsv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", help="Genes to run pipeline on", required=True)
    parser.add_argument(
        "--overwrite", help="Overwrite previous paths", action="store_true"
    )
    parser.add_argument(
        "--all-transcripts",
        help="Export all variants, regardless of clinvar, hgmd, and LOF status.",
        action="store_true",
    )

    args = parser.parse_args()

    main(args)
