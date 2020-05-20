#!/usr/bin/env python

import hail as hl
from gnomad.sample_qc.relatedness import (
    get_duplicated_samples,
    get_duplicated_samples_ht,
    infer_families,
    get_relationship_expr,
    explode_duplicate_samples_ht,
    PARENT_CHILD,
    SECOND_DEGREE_RELATIVES,
    DUPLICATE_OR_TWINS,
    UNRELATED,
    SIBLINGS
)
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.liftover import lift_data, get_liftover_genome

from gnomad_qc.v2.resources.sample_qc import qc_mt_path



import re
import io
import logging
import argparse
import itertools
from collections import Counter, defaultdict
from typing import Tuple, List
from os.path import dirname
import matplotlib.pyplot as plt

logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("check pedigree")
logger.setLevel(logging.INFO)



def check_subset(input_mt: hl.MatrixTable,
                 pedigree: hl.Pedigree,
                 output_dir: str,
                 output_name: str
                ) -> Tuple[hl.MatrixTable, list, list]:
    """
    Filters the MatrixTable to biallelics SNVs on the autosomes
    :param MatrixTable input_mt: MatrixTable
    :param hl.Pedigree pedigree: pedigree file
    :return: mt subsetted to the samples given in the pedigree, list of samples in the pedigree, list of samples in the vcf
    :rtype: Tuple[hl.MatrixTable, list, list]
    """
    
    # take sample names to subset from the pedigree
    samples_to_subset =  hl.literal(pedigree.Individual_ID.collect())
    mt_subset = input_mt.filter_cols(samples_to_subset.contains(input_mt['s']))

    # filter to variants that have at least one alt call after the subsetting
    mt_subset = mt_subset.filter_rows(hl.agg.any(mt_subset.GT.is_non_ref()))

    # check that the samples in the pedigree are present in the vcf subset
    OUT = hl.hadoop_open("{output_dir}/{output_name}_missing_samples_in_subset.txt".format(**locals()),'w')
    expected_samples = pedigree.Individual_ID.collect()
    vcf_samples = mt_subset.s.collect()
    
    missings = set(expected_samples) - set(vcf_samples)
    for i in missings:
        OUT.write(i + '\n')
    OUT.close()


    
    return(mt_subset, expected_samples, vcf_samples)


def remap_samples(original_mt_path: str,
                 input_mt: hl.MatrixTable,
                 pedigree: hl.Pedigree,
                 kin_ht: hl.Table,
                 sample_sexes: dict
                ) -> Tuple[hl.MatrixTable, hl.Table]:
    """
    Renames s col in the MatrixTable and i and j cols in the kinship ht based on seqr remap files
    :param str original_mt_path: path to original MatrixTable location
    :param hl.MatrixTable input_mt: MatrixTable 
    :param hl.Pedigree pedigree: pedigree file
    :param hl.Table kin_ht: kinship table
    :param dict sample_sexes: dictionary of sample sexes
    :return: mt, kinship ht, and sex dictionary with sample names remapped
    :rtype: Tuple[hl.MatrixTable, hl.Table, dict]
    """


    base_path = "/".join(dirname(original_mt_path).split("/")[:-1]) + ("/base/projects")
    project_list = list(set(pedigree.Project.collect()))

    remap_hts = []
    for i in project_list:
        remap = f'{base_path}/{i}/{i}_remap.tsv'
        if hl.hadoop_is_file(remap):
            remap_ht = hl.import_table(remap)
            remap_ht = remap_ht.key_by('s', 'seqr_id')
            remap_hts.append(remap_ht)
        

    if len(remap_hts) > 0:
        ht = remap_hts[0]
        for next_ht in remap_hts[1:]:
            ht = ht.join(next_ht, how = "outer")

        rename_dict = dict(hl.tuple([ht.s, ht.seqr_id]).collect())

        ht = ht.key_by('s')
        input_mt = input_mt.annotate_cols(seqr_id = ht[input_mt.col_key].seqr_id)
        input_mt = input_mt.key_cols_by(s = hl.if_else(hl.is_missing(input_mt.seqr_id), input_mt.s, input_mt.seqr_id))
        kin_ht = kin_ht.annotate(seqr_id_i = ht[kin_ht.i].seqr_id)
        kin_ht = kin_ht.annotate(i = hl.if_else(hl.is_missing(kin_ht.seqr_id_i), kin_ht.i, kin_ht.seqr_id_i))
        kin_ht = kin_ht.annotate(seqr_id_j = ht[kin_ht.j].seqr_id)
        kin_ht = kin_ht.annotate(j = hl.if_else(hl.is_missing(kin_ht.seqr_id_j), kin_ht.j, kin_ht.seqr_id_j))
        kin_ht = kin_ht.drop('seqr_id_j', 'seqr_id_i')

        sexes_renamed = {}
        for sample,sex in sample_sexes.items():
            if sample in rename_dict:
                new_sample_id = rename_dict[sample]
                sexes_renamed[new_sample_id] = sex
            else:
                sexes_renamed[sample] = sex

    return(input_mt, kin_ht, sexes_renamed)



def filter_to_biallelic_auto_snvs(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filters the MatrixTable to biallelics SNVs on the autosomes
    :param MatrixTable input_mt: MatrixTable
    :return: mt containing only snvs, biallelic variants, high callrate variants, and no rare variants
    :rtype: hl.MatrixTable
    """
    
    # filter to biallelic SNVs
    input_mt = input_mt.filter_rows((hl.len(input_mt.alleles) == 2) & hl.is_snp(input_mt.alleles[0], input_mt.alleles[1]))
    input_mt = hl.variant_qc(input_mt) 
    # only take variants with a high callrate
    input_mt = input_mt.filter_rows(input_mt.variant_qc.call_rate > 0.99)

    # recalculate AF since removed samples
    # could this just be pulled from variant_qc?
    input_mt = input_mt.annotate_rows(info = input_mt.info.annotate(AC = hl.agg.sum(input_mt.GT.n_alt_alleles()), 
                                           AN = 2 * hl.agg.count_where(hl.is_defined(input_mt.GT))))
    input_mt = input_mt.annotate_rows(info = input_mt.info.annotate(AF = input_mt.info.AC/input_mt.info.AN))
    # remove variants with low AF
    input_mt = input_mt.filter_rows(input_mt.info.AF > 0.001)
    input_mt = filter_to_autosomes(input_mt)
    
    return(input_mt)
    


def ld_prune(input_mt: hl.MatrixTable, build: str, gnomad_ld: bool) -> hl.MatrixTable:
    """
    LD prunes the MatrixTable
    :param MatrixTable input_mt: MatrixTable
    :param str input_mt: build for the input MatrixTable
    :return: ld-pruned MatrixTable
    :rtype: hl.MatrixTable
    """
    
    if gnomad_ld == False:
        mm_pruned = hl.ld_prune(input_mt.GT, r2=0.1)
        input_mt = input_mt.filter_rows(hl.is_defined(mm_pruned[input_mt.row_key]))
    else:
        # borrow from gnomad ld pruning
        pruned_mt = hl.read_matrix_table(qc_mt_path('joint', ld_pruned=True))

        if build == "GRCh38":
            pruned_mt_path = "gs://seqr-kml/pedigree_inference_resources/pruned_mt_lift.mt"
            pruned_mt = hl.read_matrix_table(pruned_mt_path)

        input_mt = input_mt.filter_rows(hl.is_defined(pruned_mt.index_rows(input_mt.row_key)))
        
        
    return(input_mt)



def filter_kin_ht(
    ht: hl.Table,
    out_summary: io.TextIOWrapper,
    out_review: io.TextIOWrapper,
    first_degree_pi_hat = .40,
    grandparent_pi_hat = .20,
    grandparent_ibd1 = 0.25,
    grandparent_ibd2 = 0.15
) -> hl.Table:
    
    """
    Filters the kinship table to relationships of grandparents and above
    :param ht: hl.Table
    :param first_degree_pi_hat: float = 0.40
    :param grandparent_pi_hat: float = 0.20
    :param grandparent_ibd1: float = 0.25
    :param grandparent_ibd2: float = 0.15
    :return: ht containing only relationships of grandparents and above
    :rtype: hl.Table
    """
    
    # round the values here?

    # filter to anything above the relationship of a grandparent
    ht = ht.filter((ht.pi_hat > first_degree_pi_hat) | ((ht.pi_hat > grandparent_pi_hat) & (ht.ibd1 > grandparent_ibd1) & (ht.ibd2 < grandparent_ibd2)))
    ht = ht.annotate(pair = hl.sorted([ht.i, ht.j])) #better variable name

    out_summary.write("NOTE: kinship table was filtered to:\n(kin > {first_degree_pi_hat}) or ".format(**locals()))
    out_summary.write("(kin > {grandparent_pi_hat} and ibd1 > {grandparent_ibd1} and ibd2 > {grandparent_ibd2})\n".format(**locals()))
    out_summary.write("relationships not meeting this critera were not evaluated\n\n")

    return(ht)



def separate_double_cousins_and_uncles(
    ht: hl.Table,
    ibd1_second_degree_threshold = 0.40,
    ibd0_second_degree_threshold = 0.40
) -> hl.Table:

    """
    Distinguish between uncle/grandparent and double-first cousins, setting the latter to unrelated
    Remove any relationship below uncle/grandparent from the second degree relatives category
    :param ht: hl.Table
    :param ibd1_second_degree_threshold: float = 0.40, ibd1 threshold below which will be call unrelated
    :param ibd0_second_degree_threshold: float = 0.40, ibd0 threshold below which will be call unrelated
    :return: ht with more specific purning to help ensure double cousins are filtered out
    :rtype: hl.Table
    """

    ht = ht.annotate(relationship = hl.if_else(
        ((ht.relationship == SECOND_DEGREE_RELATIVES) & 
         ((ht.ibd1 < ibd1_second_degree_threshold) | (ht.ibd0 < ibd0_second_degree_threshold))),
        UNRELATED,
        ht.relationship))
    
    return(ht)



def create_rank_table(
    ht: hl.Table,
) -> hl.Table:

    """
    Create rank table
    :param ht: hl.Table
    :return: ht with added annotation for rank
    :rtype: hl.Table
    """
    
    samples = hl.tuple([ht.i, ht.j]).collect()
    merged = set(list(itertools.chain(*samples)))
    v = 0
    ranked = []
    for i in merged:
        ranked.append({'s':i, 'rank':v})
        v += 1

    rank_table = hl.Table.parallelize(hl.literal(ranked, 'array<struct{s: str, rank: int}>'))
    rank_table = rank_table.key_by('s')
    
    return(rank_table)



def evaluate_dups(
    ht: hl.Table,
    master_kin_ht: hl.Table,
    dup_ht: hl.Table,
    out_summary: io.TextIOWrapper
) -> Tuple[hl.Table, hl.Table, list]:

    """
    Annotate master kinship table with duplicate decisions
    :param ht: hl.Table
    :param master_ht: hl.Table
    :param dup_ht: hl.Table
    :return: ht with duplicate samples removed and master kinship ht with marked duplicate pairs/annotations on if the pairwise relationship consists of a removed duplicate sample or if the pair was a duplicate, dups to remove
    :rtype: Tuple[hl.Table, hl.Table, list]
    """
    
    dup_exploded_ht = explode_duplicate_samples_ht(dup_ht)
    dups_to_remove = dup_exploded_ht.aggregate(hl.agg.filter(dup_exploded_ht.dup_filtered, hl.agg.collect(dup_exploded_ht.s)))
    hl_dups_to_remove = hl.literal(dups_to_remove) if dups_to_remove else hl.empty_array(hl.tstr)

    # create new kin_ht with dups removed
    ht = ht.filter((hl_dups_to_remove.contains(ht.i))|(hl_dups_to_remove.contains(ht.j)), keep=False)


    # dups_list = [sorted(list(i)) for i in dups2]
    dups_list = dup_exploded_ht.aggregate(hl.agg.filter(~dup_exploded_ht.dup_filtered, hl.agg.collect(dup_exploded_ht.s)))
    hl_dups_list = hl.literal(dups_list) if dups_list else hl.empty_array(hl.tarray(hl.tstr))


    # annotate master ht with duplicate annotations
    master_kin_ht = master_kin_ht.annotate(call = hl.if_else(master_kin_ht.relationship == DUPLICATE_OR_TWINS, "duplicates", "unknown"),
                                          final_evaluation = hl.if_else(master_kin_ht.relationship == DUPLICATE_OR_TWINS, "inferred_only_duplicates", "unknown")
                                          )
    # annotate master ht with duplicate removal annotations
    master_kin_ht = master_kin_ht.annotate(call = hl.if_else(((hl_dups_to_remove.contains(master_kin_ht.i)) & ((hl_dups_to_remove.contains(master_kin_ht.j))==False)) & (master_kin_ht.call != "duplicates")|
                                                  ((hl_dups_to_remove.contains(master_kin_ht.j)) & ((hl_dups_to_remove.contains(master_kin_ht.i))==False)) & (master_kin_ht.call != "duplicates")|
                                                  ((hl_dups_to_remove.contains(master_kin_ht.j)) & ((hl_dups_to_remove.contains(master_kin_ht.i)))),
                                                 "dup_removal",
                                                 master_kin_ht.call),
                                          final_evaluation = hl.if_else(((hl_dups_to_remove.contains(master_kin_ht.i)) & ((hl_dups_to_remove.contains(master_kin_ht.j))==False)) & (master_kin_ht.call != "duplicates")|
                                                  ((hl_dups_to_remove.contains(master_kin_ht.j)) & ((hl_dups_to_remove.contains(master_kin_ht.i))==False)) & (master_kin_ht.call != "duplicates")|
                                                  ((hl_dups_to_remove.contains(master_kin_ht.j)) & ((hl_dups_to_remove.contains(master_kin_ht.i)))),
                                                 "contains_removed_dup",
                                                 master_kin_ht.final_evaluation))


    out_summary.write("There were " +  str(len(dups_to_remove)) + " duplicate samples removed" + '\n')
    for i in dups_to_remove:
        out_summary.write(str(i) + '\n')
    out_summary.write('\n')
    
    return(ht, master_kin_ht, dups_to_remove)



def get_given_ped_rels(
    input_pedigree: str,
    vcf_samples: list,
    output_dir: io.TextIOWrapper,
    output_name: str
) -> Tuple[hl.Pedigree, list, dict, dict]:

    """
    Get the trios and duos given by the pedigree, also write out a new pedigree that does not contain missing samples
    :param input_pedigree: str
    :param vcf_samples: dict
    :return: trios given in the pedigree, duos given in the pedigree, and dictionary of each sample's family ID, dictionary of project IDs for each sample
    :rtype: Tuple[hl.Pedigree, list, dict, dict]
    """

    ped_trios = []
    ped_duos = []

    given_families = {}
    seqr_projects = defaultdict(str)
    
    
    out_new_ped = hl.hadoop_open("{output_dir}/{output_name}_functioning_pedigree.ped".format(**locals()),'w')
    
    
    with hl.hadoop_open(input_pedigree, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.rstrip('\n')
            items = line.split('\t')
            Project, Family_ID, Individual_ID, Paternal_ID, Maternal_ID, Sex, Affected = items[0:7]
        
            given_families[Individual_ID] = Family_ID
        
            if Individual_ID not in vcf_samples:
                Individual_ID = "."
            if Paternal_ID not in vcf_samples:
                Paternal_ID = "."
            if Maternal_ID not in vcf_samples:
                Maternal_ID = "."
            
            
            # remove line from pedigree if the proband is missing
            if Individual_ID != ".":
                seqr_projects[Individual_ID] = Project
                out_new_ped.write("{Family_ID}\t{Individual_ID}\t{Paternal_ID}\t{Maternal_ID}\t{Sex}\t{Affected}\n".format(**locals()))
         
            # for comparision purposes to ibd output, set sex and family to none (sex is already checked in a separate step)
            Sex = None
            Family_ID = None
         
            if Individual_ID in vcf_samples:
                if Paternal_ID in vcf_samples and Maternal_ID in vcf_samples:
                    ped_trios.append(hl.Trio(Individual_ID, Family_ID, Paternal_ID, Maternal_ID,Sex))
                elif Paternal_ID in vcf_samples and Maternal_ID not in vcf_samples:
                    ped_duos.append(sorted([Individual_ID,Paternal_ID]))
                elif Paternal_ID not in vcf_samples and Maternal_ID in vcf_samples:
                    ped_duos.append(sorted([Individual_ID,Maternal_ID]))
                

    given_trios = hl.Pedigree(ped_trios)
    out_new_ped.close()
    
    return(given_trios, ped_duos, given_families, seqr_projects)



def annotate_parent_decisions(
    trios: hl.Pedigree,
    duos: list,
    ht: hl.Table
) -> hl.Table:

    """
    Annotates who the parents are for trios and duos
    :param input_pedigree: str
    :param vcf_samples: dict
    :return: ht with annotated decision columns for who the parent(s) of the sample were
    :rtype: hl.Table
    """
    
    # filter out confirmed siblings
    parent_decisions = {}

    for trio in trios:
        parent_decisions[trio.s] = [trio.pat_id, trio.mat_id]  # two parents
    
    for duo in duos:
        s1 = duo[0]
        s2 = duo[1]
        parent_decisions[s1] = [str(s2)]  # duos
        parent_decisions[s2] = [str(s1)]  # duos
    
    
    hl_final_decisions = hl.literal(parent_decisions)

    
    ht = ht.annotate(decision_i = hl_final_decisions.get(ht.i),
                                    decision_j = hl_final_decisions.get(ht.j))


    # set NAs to no decision
    dec_exprj = hl.if_else(hl.is_missing(ht.decision_j), ["no_decision"],  ht.decision_j)
    ht = ht.annotate(decision_j = dec_exprj)

    dec_expri = hl.if_else(hl.is_missing(ht.decision_i), ["no_decision"],  ht.decision_i)
    ht = ht.annotate(decision_i = dec_expri)
    
    return(ht)



def evaluate_trios(
    ht: hl.Table,
    master_kin_ht: hl.Table,
    out_summary: io.TextIOWrapper,
    out_review: io.TextIOWrapper,
    trios: hl.Pedigree,
    given_trios: hl.Pedigree
) -> Tuple[hl.Table, hl.Table, List[hl.Trio], hl.Trio, List[hl.Trio]]:

    """
    Determine if trios are confirmed, inferred only, or given only, annotate master kinship table with trio decisions
    :param ht: hl.Table
    :param master_ht: hl.Table
    :return: ht with trio samples removed, master kinship ht with annotated trio decisions, given-only trios, complete inferred trios, inferred only trios
    :rtype: Tuple[hl.Table, hl.Table, List[hl.Trio], hl.Trio, List[hl.Trio]]
    """
    
    # set fam_id and is_female of inferred ped to none
    inf_trios = []
    for trio in trios.complete_trios():
        fam_id = None
        is_female = None
        inf_trios.append(hl.Trio(trio.s, fam_id, trio.pat_id, trio.mat_id, is_female))

    inf_trios = hl.Pedigree(inf_trios)

    # find the compelete trios
    complete_inferred_trios = inf_trios.complete_trios()
    complete_given_trios = given_trios.complete_trios()

    # find differences between given and inferred trios
    confirmed_trios = []
    inferred_only_trios = []
    given_only_trios = []
    inferred_only_trio_pairs = []
    confirmed_trio_pairs = []

    # find confirmed and inferred only trios
    for trio in complete_inferred_trios:
        if trio in complete_given_trios:
            confirmed_trios.append(trio)
            confirmed_trio_pairs.append([trio.s, trio.pat_id]) 
            confirmed_trio_pairs.append([trio.s, trio.mat_id])
        else:
            inferred_only_trios.append(trio)
            inferred_only_trio_pairs.append([trio.s, trio.pat_id])
            inferred_only_trio_pairs.append([trio.s, trio.mat_id])
            out_review.write("inferred_only_trio\t[s:{trio.s}, pat_id:{trio.pat_id}, mat_id:{trio.mat_id}]\t.\n".format(**locals()))


    inferred_only_trio_pairs = [sorted(x) for x in inferred_only_trio_pairs]
    confirmed_trio_pairs = [sorted(x) for x in confirmed_trio_pairs]

    # find given only trios
    for trio in complete_given_trios:
        if trio not in complete_inferred_trios:
            given_only_trios.append(trio)
            out_review.write("given_only_trio\t[s:{trio.s}, pat_id:{trio.pat_id}, mat_id:{trio.mat_id}]\t.\n".format(**locals()))

    trios_accounted = inferred_only_trios + confirmed_trios  #(given only won't be in ht so don't have to account for those)

    accounted_relations = []
    for trio in trios_accounted:
        accounted_relations.append(sorted([trio.s, trio.pat_id]))
        accounted_relations.append(sorted([trio.s, trio.mat_id]))
        
   
    out_summary.write("Number of confirmed trios: " + str(len(confirmed_trios)) + '\n')
    out_summary.write("Number of inferred only trios: " + str(len(inferred_only_trios)) + '\n')
    out_summary.write("Number of given only trios: " + str(len(given_only_trios)) + '\n\n')

    # removed the inferred trios that have been accounted for already from the next iteration of the kinship ht
    ht = ht.filter(hl.literal(accounted_relations).contains(ht.pair), keep=False)

    # annotate master ht with trio annotations
    master_kin_ht = master_kin_ht.annotate(call = hl.if_else(hl.literal(accounted_relations).contains(master_kin_ht.pair),"parent_offspring", master_kin_ht.call))
    master_kin_ht = master_kin_ht.annotate(final_evaluation = hl.if_else(hl.literal(inferred_only_trio_pairs).contains(master_kin_ht.pair), "inferred_only_trio_component", master_kin_ht.final_evaluation))
    master_kin_ht = master_kin_ht.annotate(final_evaluation = hl.if_else(hl.literal(confirmed_trio_pairs).contains(master_kin_ht.pair), "confirmed_trio_component", master_kin_ht.final_evaluation))

    return(ht, master_kin_ht, given_only_trios, complete_inferred_trios, inferred_only_trios)



def evaluate_duos(
    ht: hl.Table,
    master_kin_ht: hl.Table,
    out_summary: io.TextIOWrapper,
    out_review: io.TextIOWrapper,
    ped_duos: list
) -> Tuple[hl.Table, hl.Table, list, list]:

    """
    Determine if duos are confirmed, inferred only, or given only, annotate master kinship table with duo decisions
    :param ht: hl.Table
    :param master_ht: hl.Table
    :return: ht with duos removed and master kinship ht with annotated duos decisions, inferred duos, and given-only duos
    :rtype: Tuple[hl.Table, hl.Table, list, list]
    """
    
    parent_child_ht = ht.filter(ht.relationship == PARENT_CHILD)
    parent_child_ht_i = parent_child_ht.aggregate(hl.agg.group_by(parent_child_ht.i, hl.agg.collect(parent_child_ht.j)))
    parent_child_ht_j = parent_child_ht.aggregate(hl.agg.group_by(parent_child_ht.j, hl.agg.collect(parent_child_ht.i)))

    # for at least one of the samples in the parent-offspring pair, that sample must appear in a parent offspring pair only once.
    potential_duos = {**parent_child_ht_i , **parent_child_ht_j}  # combine the two dictionaries


    duos = []
    for sample,possible_parents in potential_duos.items():
        if len(possible_parents) == 1:
            duos.append(sorted([sample,possible_parents[0]]))
    duos = sorted(set(map(tuple, duos)))  # get only uniq duos
    duos = [list(x) for x in duos]

    inferred_duos = duos
    given_duos = ped_duos
                    
    # find differences between given and inferred trios
    confirmed_duos = []
    inferred_only_duos = []
    given_only_duos = []
                              
    for duo in inferred_duos:
        if duo in given_duos:
            confirmed_duos.append(duo)
        else:
            inferred_only_duos.append(duo)

    for duo in given_duos:
        if duo not in inferred_duos:
            given_only_duos.append(duo)
            duo_out = re.sub("'", "", str(duo)) 
            out_review.write("given_only_duo\t{duo_out}\t.\n".format(**locals()))


    duos_accounted = inferred_only_duos + confirmed_duos

    hl_duos_accounted  = hl.literal(duos_accounted) if duos_accounted else hl.empty_array(hl.tarray(hl.tstr))
    hl_inferred_only_duos  = hl.literal(inferred_only_duos) if inferred_only_duos else hl.empty_array(hl.tarray(hl.tstr))
    hl_confirmed_duos  = hl.literal(confirmed_duos) if confirmed_duos else hl.empty_array(hl.tarray(hl.tstr))

    # remove duos from temp kin ht
    ht = ht.filter(hl_duos_accounted.contains(ht.pair), keep=False)

    # annotate master kin with dup annotations
    master_kin_ht = master_kin_ht.annotate(call = hl.if_else(hl_duos_accounted.contains(master_kin_ht.pair), "parent_offspring", master_kin_ht.call))
    master_kin_ht = master_kin_ht.annotate(final_evaluation = hl.if_else(hl_inferred_only_duos.contains(master_kin_ht.pair), "inferred_only_duo", master_kin_ht.final_evaluation))
    master_kin_ht = master_kin_ht.annotate(final_evaluation = hl.if_else(hl_confirmed_duos.contains(master_kin_ht.pair), "confirmed_duo", master_kin_ht.final_evaluation))

    out_summary.write("Number of confirmed duos: " + str(len(confirmed_duos)) + '\n')
    out_summary.write("Number of inferred only duos: " + str(len(inferred_only_duos)) + '\n')
    out_summary.write("Number of given only duos: " + str(len(given_only_duos)) + '\n\n')
    
    return(ht, master_kin_ht, duos, given_only_duos)



def evaluate_siblings(
    ht: hl.Table,
    master_kin_ht: hl.Table,
    out_summary: io.TextIOWrapper
) -> Tuple[hl.Table, hl.Table, hl.expr.ArrayExpression]:

    """
    Determine if siblings are confirmed or inferred only, annotate master kinship table with sibling decisions
    :param ht: hl.Table
    :param master_ht: hl.Table
    :return: ht with sibling samples removed and master kinship ht with annotated sibling decisions, sibling pairs
    :rtype: Tuple[hl.Table, hl.Table, hl.expr.ArrayExpression]
    """
    
    sibling_ht = ht.filter(ht.relationship == SIBLINGS)
    
    # confirmed sibs are those where the final decision (parents) matched
    # need to check they are not both none/no decision because that's not really a match
    confirmed_sibs = sibling_ht.filter((sibling_ht.decision_i == sibling_ht.decision_j) & (sibling_ht.decision_i != ["no_decision"]))
    confirmed_sibs_list = confirmed_sibs.pair.collect()
    hl_confirmed_sibs_list = hl.literal(confirmed_sibs_list) if confirmed_sibs_list else hl.empty_array(hl.tarray(hl.tstr))

    unconfirmed_sibs = sibling_ht.filter(hl_confirmed_sibs_list.contains(sibling_ht.pair), keep=False)
    unconfirmed_sibs_list = unconfirmed_sibs.pair.collect()
    unconfirmed_sibs_list = [sorted(x) for x in unconfirmed_sibs_list]
    hl_unconfirmed_sibs_list = hl.literal(unconfirmed_sibs_list) if unconfirmed_sibs_list else hl.empty_array(hl.tarray(hl.tstr))


    sibs_accounted = sibling_ht.pair.collect()
    hl_sibs_accounted = hl.literal(sibs_accounted) if sibs_accounted else hl.empty_array(hl.tarray(hl.tstr))

    # now can remove all sibling relations from the temp ht since they been accounted for
    ht = ht.filter(hl_sibs_accounted.contains(ht.pair), keep=False)
    sibs_accounted_list = [sorted(x) for x in sibs_accounted]

    hl_sibs_accounted_list = hl.literal(sibs_accounted_list) if sibs_accounted_list else hl.empty_array(hl.tarray(hl.tstr))

    # annotate master ht with sibling annotations
    master_kin_ht = master_kin_ht.annotate(call = hl.if_else(hl_sibs_accounted_list.contains(master_kin_ht.pair), "siblings", master_kin_ht.call))
    master_kin_ht = master_kin_ht.annotate(final_evaluation = hl.if_else(hl_confirmed_sibs_list.contains(master_kin_ht.pair), "confirmed_siblings", master_kin_ht.final_evaluation))
    master_kin_ht = master_kin_ht.annotate(final_evaluation = hl.if_else(hl_unconfirmed_sibs_list.contains(master_kin_ht.pair), "inferred_only_siblings", master_kin_ht.final_evaluation))

    out_summary.write("Number of confirmed siblings: " + str(confirmed_sibs.count()) + '\n')
    out_summary.write("Number of inferred only siblings: " + str(unconfirmed_sibs.count()) + '\n\n')
    
    return(ht, master_kin_ht, hl_sibs_accounted_list)


def evaluate_grandparents(
    ht: hl.Table,
    hl_sibs_accounted_list: hl.expr.ArrayExpression,
    master_kin_ht: hl.Table,
    out_summary: io.TextIOWrapper
) -> Tuple[hl.Table, hl.Table]:
    
    """
    Determine if grandparents are confirmed or inferred only, annotate master kinship table with grandparent decisions
    :param ht: hl.Table
    :param master_ht: hl.Table
    :return: ht with grandparent samples removed and master kinship ht with annotated grandparent decisions
    :rtype: Tuple[hl.Table, hl.Table]
    """

    # filter kinship ht to the grandparents
    grandparent_ht = ht.filter(ht.relationship == SECOND_DEGREE_RELATIVES)


    # test if either parents have the paired sample as a sib in the unconfirmed sib list, then it would make sense to see grand/nephew relations
    # nephew/uncle: check in siblings ht for presence of the parent of the nephew and the grandparent
    # can use [0] and [-1] indices cause there will never be more than two samples here because they are only parent decisions

    expr = hl.if_else(hl_sibs_accounted_list.contains(hl.sorted([grandparent_ht.decision_j[0], grandparent_ht.i]))|
                   hl_sibs_accounted_list.contains(hl.sorted([grandparent_ht.decision_j[-1], grandparent_ht.i]))|
                   hl_sibs_accounted_list.contains(hl.sorted([grandparent_ht.decision_i[0], grandparent_ht.j]))|
                   hl_sibs_accounted_list.contains(hl.sorted([grandparent_ht.decision_i[-1], grandparent_ht.j])),
                   True,
                   False)



    grandparent_ht = grandparent_ht.annotate(uncle_evidence = expr)


    confirmed_grandparents = grandparent_ht.filter(((grandparent_ht.decision_i == grandparent_ht.decision_j) & (grandparent_ht.decision_i != ["no_decision"]))|
                                                 grandparent_ht.uncle_evidence)

    confirmed_grandparents_list = confirmed_grandparents.pair.collect()
    hl_confirmed_grandparents_list = hl.literal(confirmed_grandparents_list) if confirmed_grandparents_list else hl.empty_array(hl.tarray(hl.tstr))

    unconfirmed_grandparents = grandparent_ht.filter(hl_confirmed_grandparents_list.contains(grandparent_ht.pair), keep=False)
    unconfirmed_grandparents_list = unconfirmed_grandparents.pair.collect()
    unconfirmed_grandparents_list = [sorted(x) for x in unconfirmed_grandparents_list]

    # maybe can check confirmed with other method? check if way to pull from ped/families? will this be added in notes column?
    grandparents_accounted = grandparent_ht.pair.collect() 
    hl_grandparents_accounted = hl.literal(grandparents_accounted) if grandparents_accounted else hl.empty_array(hl.tarray(hl.tstr))
    hl_unconfirmed_grandparents_list = hl.literal(unconfirmed_grandparents_list) if unconfirmed_grandparents_list else hl.empty_array(hl.tarray(hl.tstr))

    # now can remove all grandparent relations from the ht since they been accounted for
    ht = ht.filter(hl_grandparents_accounted.contains(ht.pair),keep=False)

    master_kin_ht = master_kin_ht.annotate(call = hl.if_else(hl_grandparents_accounted.contains(master_kin_ht.pair), "grandparent_or_uncle", master_kin_ht.call))
    master_kin_ht = master_kin_ht.annotate(final_evaluation = hl.if_else(hl_confirmed_grandparents_list.contains(master_kin_ht.pair), "confirmed_grandparent_or_uncle", master_kin_ht.final_evaluation))
    master_kin_ht = master_kin_ht.annotate(final_evaluation = hl.if_else(hl_unconfirmed_grandparents_list.contains(master_kin_ht.pair), "inferred_only_grandparent_or_uncle", master_kin_ht.final_evaluation))

    n_unknown = master_kin_ht.filter(master_kin_ht.call == "unknown").count()
    out_summary.write("Number of confirmed grandparents_uncles: " + str(confirmed_grandparents.count()) + '\n')
    out_summary.write("Number of inferred only grandparents_uncles: " + str(unconfirmed_grandparents.count()) + '\n')
    out_summary.write("Remaining unknown relationships: " + str(n_unknown) + '\n\n')
    
    return(ht, master_kin_ht)

def plot_ibd(ht, output_file):
    df3 = ht.to_pandas()


    colors = {'duplicates':'yellowgreen', 
              'grandparent_or_uncle':'blue', 
              'parent_offspring':'green', 
              'siblings':'black',
              'dup_removal':'orange',
              'unknown':'magenta',
              'not_included':'cyan'}


    groups = df3.groupby('call')
    plt.rc('axes', labelsize=35)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=20)    # fontsize of the tick labels

    params = {'legend.fontsize': 17}
    plt.rcParams.update(params)

    fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1)

    for name, group in groups:
        ax1.plot(group.ibd0, group.ibd1, alpha = .5, marker='.', linestyle='', ms=18, label=name,color=colors[name])
        ax2.plot(group.ibd1, group.ibd2, alpha = .5, marker='.', linestyle='', ms=18, label=name,color=colors[name])
        ax3.plot(group.ibd0, group.ibd2, alpha = .5, marker='.', linestyle='', ms=18, label=name,color=colors[name])
        ax4.plot(group.pi_hat, group.ibd2, alpha = .5, marker='.', linestyle='', ms=18, label=name,color=colors[name])
    
    ax1.set_xlabel('ibd0')
    ax1.set_ylabel('ibd1')
    ax2.set_xlabel('ibd1')
    ax2.set_ylabel('ibd2')
    ax3.set_xlabel('ibd0')
    ax3.set_ylabel('ibd2')
    ax4.set_xlabel('pi_hat')
    ax4.set_ylabel('ibd2')
    ax4.legend(bbox_to_anchor=(.575, -.2))

    plt.show()
    fig.set_size_inches(20, 30)

    with hl.hadoop_open(output_file, 'wb') as out:
        fig.savefig(out)



def main(args):

    output_dir = args.output_dir
    output_name = args.output_name
    sex_check = args.sex_check
    input_mt = args.input_mt
    input_pedigree = args.input_pedigree

    gnomad_ld = args.gnomad_ld
    run_ibd = args.run_ibd
    first_degree_pi_hat = args.first_degree_pi_hat
    grandparent_pi_hat = args.grandparent_pi_hat
    grandparent_ibd1 = args.grandparent_ibd1

    grandparent_ibd2 = args.grandparent_ibd2
    ibd1_second_degree_threshold = args.ibd1_second_degree_threshold
    ibd0_second_degree_threshold = args.ibd0_second_degree_threshold
    second_degree_min_kin = args.second_degree_min_kin
    ibd1_100_min = args.ibd1_100_min

    ibd2_0_max = args.ibd2_0_max
    ibd0_0_max = args.ibd0_0_max
    first_degree_kin_thresholds_lower = args.first_degree_kin_thresholds_lower
    first_degree_kin_thresholds_upper = args.first_degree_kin_thresholds_upper

    first_degree_kin_thresholds = (first_degree_kin_thresholds_lower, first_degree_kin_thresholds_upper)


    # read in MT and pedigree
    mt = hl.read_matrix_table(input_mt)
    pedigree = hl.import_table(input_pedigree, impute=True)
    build = get_reference_genome(mt.locus).name  # infer build of the mt


    mt = filter_to_biallelic_auto_snvs(mt)

    mt = ld_prune(mt, build, gnomad_ld)

    out_mt = "{output_dir}/{output_name}_processed_mt.mt".format(**locals())


    # IBD
    ibd_results_ht = hl.identity_by_descent(mt, maf=mt.info.AF, min=0.10, max=1.0)
    ibd_results_ht = ibd_results_ht.annotate(ibd0 = ibd_results_ht.ibd.Z0,
                                ibd1 = ibd_results_ht.ibd.Z1,
                                ibd2 = ibd_results_ht.ibd.Z2,
                                pi_hat = ibd_results_ht.ibd.PI_HAT).drop("ibs0", "ibs1", "ibs2", "ibd")



    # write out matrix and hail tables
    if run_ibd:
        logger.info("Running IBD...")

        # ibd
        out_ht = "{output_dir}/{output_name}_ibd_kinship.tsv".format(**locals())
        ibd_results_ht.export(out_ht)

        mt.write(out_mt, overwrite=True)

    else:
        logger.warn("Skipping IBD - using previous calculations...")



    logger.info('Reading in mt...')
    mt = hl.read_matrix_table("{output_dir}/{output_name}_processed_mt.mt".format(**locals()))
    mt_path = input_mt

    logger.info('Reading in kinship ht...')
    kin_ht = hl.import_table("{output_dir}/{output_name}_ibd_kinship.tsv".format(**locals()), impute=True)

    logger.info('Reading sex check into dictionary...')  # code here is based on old sex check, need new sex check that works on dense data
    sample_sexes = {}
    with hl.hadoop_open(sex_check,'r') as infile:
        next(infile)
        for line in infile:
            line = line.rstrip()
            items = line.split('\t')
            sample = items[0]
            sex = items[6]
            if sex == "female":
                sample_sexes[sample] = True
            elif sex == "male":
                sample_sexes[sample] = False


    logger.info('Remapping sample names...')
    mt, kin_ht, sample_sexes = remap_samples(mt_path, mt, pedigree, kin_ht, sample_sexes)

    # subset mt to the samples in the pedigree
    mt_subset, expected_samples, vcf_samples = check_subset(mt, pedigree, output_dir, output_name)
    num_ped = len(expected_samples)
    num_vcf = len(vcf_samples)

    # subset ht to the samples in the pedigree
    subset = hl.literal(expected_samples)
    kin_ht = kin_ht.filter(subset.contains(kin_ht.i) | subset.contains(kin_ht.j))

    # output_name = "WK_TEST"  # temp, remove this later
    kin_ht = kin_ht.key_by("i", "j")
    original_kin_ht = kin_ht  # save the original kin_ht as it will be filtered downstream

    # setup output files
    out_summary = hl.hadoop_open("{output_dir}/{output_name}_ped_check_summary.txt".format(**locals()), 'w')
    out_review = hl.hadoop_open("{output_dir}/{output_name}_needs_review.txt".format(**locals()), 'w')
    out_review.write("status\tsamples\tkin_results[kin,idb0,ibd1,ibd2]\n")

   
    logger.info('Filtering kinship table to remove unrelateds from analysis...')
    kin_ht = filter_kin_ht(kin_ht, out_summary, out_review)
    # output basic stats
    out_summary.write("Number individuals in pedigree: " + str(num_ped) + '\n')
    out_summary.write("Number individuals in subset from the VCF: " + str(num_vcf) + '\n')
    out_summary.write("Number of relationships in the kinship table: " + str(kin_ht.count()) + '\n\n')


    logger.info('Defining pairwise relationships in the kinship ht...')
    kin_ht = kin_ht.annotate(relationship = get_relationship_expr(kin_expr = kin_ht.pi_hat, 
                        ibd0_expr = kin_ht.ibd0, 
                        ibd1_expr = kin_ht.ibd1, 
                        ibd2_expr = kin_ht.ibd2,
                        first_degree_kin_thresholds = (.4, .75),
                        second_degree_min_kin = 0.195,
                        ibd1_100_min = .70,
                        ibd2_0_max = .30,
                        ibd0_0_max = .15))

    master_kin_ht = kin_ht  # maintain master_kin_ht as ht with ped-check relationships

    logger.info('Distinguishing between double cousins and uncles/grandparents...')
    kin_ht = separate_double_cousins_and_uncles(kin_ht)

    logger.info('Finding duplicate samples...')
    dups = get_duplicated_samples(kin_ht)

    logger.info('Generating rank table...')
    rank_table = create_rank_table(kin_ht)
    # decide which dups to keep based on rank
    dup_ht = get_duplicated_samples_ht(dups,rank_table)

    logger.info('Inferring trios...')
    trios = infer_families(relationship_ht = kin_ht,
                sex = sample_sexes,
                duplicate_samples_ht = dup_ht)


    logger.info('Evaluating duplicates...')
    kin_ht, master_kin_ht, dups_to_remove = evaluate_dups(kin_ht, master_kin_ht, dup_ht, out_summary)

    logger.info('Evaluating trios...')
    
    # find relationships given in the pedigree
    given_trios, given_duos, given_families, seqr_projects = get_given_ped_rels(input_pedigree, vcf_samples, output_dir, output_name)
    # determine if inferred trios and given trios match
    kin_ht, master_kin_ht, given_only_trios, complete_inferred_trios, inferred_only_trios = evaluate_trios(kin_ht, master_kin_ht, out_summary, out_review, trios, given_trios)


    logger.info('Evaluating duos (offspring with a single parent)...')
    kin_ht, master_kin_ht, duos, given_only_duos = evaluate_duos(kin_ht, master_kin_ht, out_summary, out_review, given_duos)


    logger.info('Annotating parents in the kinship ht...')
    kin_ht = annotate_parent_decisions(trios, duos, kin_ht)


    logger.info('Evaluating siblings...')
    kin_ht, master_kin_ht, hl_sibs_accounted_list = evaluate_siblings(kin_ht, master_kin_ht, out_summary)


    logger.info('Evaluating grandparents...')
    kin_ht, master_kin_ht = evaluate_grandparents(kin_ht, hl_sibs_accounted_list, master_kin_ht, out_summary)



    # add annnotation for whether or not pairs in the same family based on the given pedigree and for review status
    hl_given_families = hl.literal(given_families)

    # check if families of sample i and sample j match
    master_kin_ht = master_kin_ht.annotate(same_given_family = hl.if_else(hl_given_families.get(master_kin_ht.i) == hl_given_families.get(master_kin_ht.j),
                                                True,
                                                False))

    # annotate whether or not the relationships need review (if it was only inferred or unknown)
    master_kin_ht = master_kin_ht.annotate(status = hl.if_else((master_kin_ht.final_evaluation == "unknown") | (master_kin_ht.final_evaluation.startswith('inferred_only')),
                                        "needs_review",
                                        "ok")
                                        )

    # count number of relationships needing review
    num_review = master_kin_ht.filter(master_kin_ht.status == "needs_review").count()
       
    
    master_kin_ht.export("{output_dir}/{output_name}_annotated_kin.txt".format(**locals()))
    out_summary.write("Inferred only relationships needing review: " + str(num_review) + '\n')

    # should two samples in same fam without parents count as given only sibs/something else?
    all_given_only = given_only_trios + given_only_duos

    out_summary.write("Given only relationships needing review: " + str(len(all_given_only)) + '\n')
    out_summary.close()


    # number inferred needing review won't match lines in out review file becuase trio-components are collapsed into trios, and given only are also output
    review = master_kin_ht.filter((master_kin_ht.status == "needs_review") & (master_kin_ht.final_evaluation != "inferred_only_trio_component"))

    df = review.to_pandas()
    df = df.sort_values(by=['final_evaluation'])
    for index, row in df.iterrows():
        out_review.write(row['final_evaluation'] + '\t' + str(re.sub("'","",str(row['pair']))) + '\t[' + ','.join(map(str,([row['pi_hat'],row['ibd0'],row['ibd1'],row['ibd2']]))) + "]\n")
    out_review.close()

    master_kin_ht.export("{output_dir}/{output_name}_annotated_kin.txt".format(**locals()))
    master_kin_ht.write("{output_dir}/{output_name}_annotated_kin.ht".format(**locals()), overwrite=True)


    # count how many of each discrepancy or confirmed
    hl_seqr_projects = hl.literal(seqr_projects)

    # add annotation for seqr projects of sample i and sample j
    master_kin_ht = master_kin_ht.annotate(seqr_id_i = hl_seqr_projects.get(master_kin_ht.i),
                                        seqr_id_j = hl_seqr_projects.get(master_kin_ht.j))

    # count samples belonging to each project in the MT
    mt_subset = mt_subset.annotate_cols(seqr_id = hl_seqr_projects.get(mt_subset.s))
    mt_subset_cols = mt_subset.cols()
    n_seqr_ids = mt_subset_cols.group_by(mt_subset_cols.seqr_id).aggregate(n=hl.agg.count())

    # count how many individuals in each project in the ped and make dictionary
    seqr_project_counter = Counter(seqr_projects.values())


    # add annotations from master ht back to the original ht and set anything not in the master ht to "not_included"
    original_kin_ht = original_kin_ht.key_by(original_kin_ht.i, original_kin_ht.j)
    master_kin_ht = master_kin_ht.key_by(master_kin_ht.i, master_kin_ht.j)
    joiner_ht = master_kin_ht.select(master_kin_ht.call)
    original_kin_ht =  original_kin_ht.annotate(**joiner_ht[original_kin_ht.i, original_kin_ht.j])
    original_kin_ht = original_kin_ht.annotate(call = hl.if_else(hl.is_defined(original_kin_ht.call), original_kin_ht.call, "not_included"))
    original_kin_ht = original_kin_ht.annotate(seqr_id_i = hl_seqr_projects.get(original_kin_ht.i),
                                        seqr_id_j = hl_seqr_projects.get(original_kin_ht.j))




    # output results per project
    OUT_PROJECT_OVERVIEW = hl.hadoop_open("{output_dir}/{output_name}_project_stats.txt".format(**locals()),'w')
    OUT_PROJECT_OVERVIEW.write("project\tinfer_only_to_review\tgiven_only_to_review\n")

    for project in set(seqr_projects.values()):  # is there a faster way to do this?
        OUT_PROJECT_SUMMARY = hl.hadoop_open("{output_dir}/{project}/{output_name}_{project}_ped_check_summary.txt".format(**locals()),'w')
        OUT_PROJECT_REVIEW = hl.hadoop_open("{output_dir}/{project}/{output_name}_{project}_needs_review.txt".format(**locals()),'w')
        OUT_PROJECT_REVIEW.write("status\tsamples\tkin_results[pi_hat,idb0,ibd1,ibd2]\n")
        
        project_ht = master_kin_ht.filter((master_kin_ht.seqr_id_i == project) | (master_kin_ht.seqr_id_j == project))
        project_ht_org = original_kin_ht.filter((original_kin_ht.seqr_id_i == project) | (original_kin_ht.seqr_id_j == project))

        n_given_only_trios = 0
        for trio in given_only_trios:
            if seqr_projects[trio.s] == project or seqr_projects[trio.pat_id] == project or seqr_projects[trio.mat_id] == project:
                n_given_only_trios += 1
                OUT_PROJECT_REVIEW.write("given_only_trio\t[s:{trio.s}, pat_id:{trio.pat_id}, mat_id:{trio.mat_id}]\t.\n".format(**locals()))
        
        n_confirmed_trios = 0
        for trio in complete_inferred_trios:
            if seqr_projects[trio.s] == project or seqr_projects[trio.pat_id] == project or seqr_projects[trio.mat_id] == project:
                n_confirmed_trios += 1

        n_inferred_only_trios = 0
        for trio in inferred_only_trios: # should change this variable name to include trio later!
            if seqr_projects[trio.s] == project or seqr_projects[trio.pat_id] == project or seqr_projects[trio.mat_id] == project:
                n_inferred_only_trios += 1
                OUT_PROJECT_REVIEW.write("inferred_only_trio\t[s:{trio.s}, pat_id:{trio.pat_id}, mat_id:{trio.mat_id}]\t.\n".format(**locals()))
        
        n_given_only_duos = 0
        for duo in given_only_duos:
            duo_out = re.sub("'","",str(duo)) 
            s1 = duo[0]
            s2 = duo[1]
            if seqr_projects[s1] == project or seqr_projects[s2] == project:
                n_given_only_duos += 1
                OUT_PROJECT_REVIEW.write("given_only_duo\t{duo_out}\t.\n".format(**locals()))
                
        # compute struct with one aggregation, hl.struct?
        n_confirmed_trio_component = project_ht.aggregate(hl.agg.count_where(project_ht.final_evaluation == "confirmed_trio_component"))
        n_inferred_only_trio_component = project_ht.aggregate(hl.agg.count_where(project_ht.final_evaluation == "inferred_only_trio_component"))
        n_confirmed_duos = project_ht.aggregate(hl.agg.count_where(project_ht.final_evaluation == "confirmed_duo"))
        n_inferred_only_duos = project_ht.aggregate(hl.agg.count_where(project_ht.final_evaluation == "inferred_only_duo"))
        n_confirmed_siblings = project_ht.aggregate(hl.agg.count_where(project_ht.final_evaluation == "confirmed_siblings"))
        n_inferred_only_siblings = project_ht.aggregate(hl.agg.count_where(project_ht.final_evaluation == "inferred_only_siblings"))
        n_confirmed_grandparent_or_uncle = project_ht.aggregate(hl.agg.count_where(project_ht.final_evaluation == "confirmed_grandparent_or_uncle"))
        n_inferred_only_grandparent_or_uncle = project_ht.aggregate(hl.agg.count_where(project_ht.final_evaluation == "inferred_only_grandparent_or_uncle"))
        n_unknown = project_ht.aggregate(hl.agg.count_where(project_ht.final_evaluation == "unknown"))

        n_relationships = project_ht.count()
        n_infer_review = project_ht.filter(project_ht.status == "needs_review").count()
        n_given_review = n_given_only_trios + n_given_only_duos
        n_samples_in_vcf = n_seqr_ids.filter(n_seqr_ids.seqr_id == project).n.collect()[0]
        n_samples_in_ped = seqr_project_counter[project]

        project_dups_to_remove = []
        for sample in dups_to_remove:
            if seqr_projects[sample] == project:
                project_dups_to_remove.append(sample)
                 
        OUT_PROJECT_SUMMARY.write("NOTE: kinship table was filtered to:\n(kin > {first_degree_pi_hat}) or ".format(**locals()))
        OUT_PROJECT_SUMMARY.write("(kin > {grandparent_pi_hat} and ibd1 > {grandparent_ibd1} and ibd2 > {grandparent_ibd2})\n".format(**locals()))
        OUT_PROJECT_SUMMARY.write("relationships not meeting this critera were not evaluated\n\n")
        OUT_PROJECT_SUMMARY.write("Number individuals in pedigree: " + str(n_samples_in_ped) + '\n')
        OUT_PROJECT_SUMMARY.write("Number individuals in subset from the VCF: " + str(n_samples_in_vcf) + '\n')
        OUT_PROJECT_SUMMARY.write("Number of relationships in the kinship table: " + str(n_relationships) + '\n\n')
           
        OUT_PROJECT_SUMMARY.write("There were " +  str(len(project_dups_to_remove)) + " duplicate samples removed" + '\n')

        for i in project_dups_to_remove:
            OUT_PROJECT_SUMMARY.write(str(i) + '\n')
        OUT_PROJECT_SUMMARY.write('\n')
            
        OUT_PROJECT_SUMMARY.write("Number of confirmed trios: " + str(n_confirmed_trios) + '\n')
        OUT_PROJECT_SUMMARY.write("Number of inferred only trios: " + str(n_inferred_only_trios) + '\n')
        OUT_PROJECT_SUMMARY.write("Number of given only trios: " + str(n_given_only_trios) + '\n\n')

        OUT_PROJECT_SUMMARY.write("Number of confirmed duos: " + str(n_confirmed_duos) + '\n')
        OUT_PROJECT_SUMMARY.write("Number of inferred only duos: " + str(n_inferred_only_duos) + '\n')
        OUT_PROJECT_SUMMARY.write("Number of given only duos: " + str(n_given_only_duos) + '\n\n')

        OUT_PROJECT_SUMMARY.write("Number of confirmed siblings: " + str(n_confirmed_siblings) + '\n')
        OUT_PROJECT_SUMMARY.write("Number of inferred only siblings: " + str(n_inferred_only_siblings) + '\n\n')

        OUT_PROJECT_SUMMARY.write("Number of confirmed grandparents_uncles: " + str(n_confirmed_grandparent_or_uncle) + '\n')
        OUT_PROJECT_SUMMARY.write("Number of inferred only grandparents_uncles: " + str(n_inferred_only_grandparent_or_uncle) + '\n')
        OUT_PROJECT_SUMMARY.write("Remaining unknown relationships: " + str(n_unknown) + '\n\n')

        OUT_PROJECT_SUMMARY.write("Inferred only relationships needing review: " + str(n_infer_review) + '\n')
        OUT_PROJECT_SUMMARY.write("Given only relationships needing review: " + str(n_given_review) + '\n')

        OUT_PROJECT_SUMMARY.close()
        
        OUT_PROJECT_OVERVIEW.write("{project}\t{n_infer_review}\t{n_given_review}\n".format(**locals()))
         
        out_review = project_ht.filter((project_ht.status == "needs_review") & (project_ht.final_evaluation != "inferred_only_trio_component"))

        df = out_review.to_pandas()
        df = df.sort_values(by=['final_evaluation'])
        for index, row in df.iterrows():
            OUT_PROJECT_REVIEW.write(row['final_evaluation'] + '\t' + str(re.sub("'", "", str(row['pair']))) + '\t[' + ','.join(map(str,([row['pi_hat'],row['ibd0'],row['ibd1'],row['ibd2']]))) + "]\n")
        OUT_PROJECT_REVIEW.close() 
        
        
        project_ht.export("{output_dir}/{project}/{output_name}_{project}_annotated_kin.txt".format(**locals()))
        plot_ibd(project_ht_org,"{output_dir}/{project}/{output_name}_{project}_ibd_org.png".format(**locals()))

    OUT_PROJECT_OVERVIEW.close() 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script checks a given pedigree to against ibd kinship output')
    parser.add_argument('-d', '--output_dir', help='path to directory to output results')
    parser.add_argument('-n', '--output_name', help='output prefix to use for results')
    parser.add_argument('-s', '--sex_check', help='path to sex check output')
    parser.add_argument('-i', '--input_mt', help='path to input MatrixTable')
    parser.add_argument('-p', '--input_pedigree', help='path to input pedigree')
    parser.add_argument('--gnomad_ld', help='use variants from gnomAD for ld prune', action='store_true')
    parser.add_argument('--run_ibd', help='run IBD', action='store_true')
    
    parser.add_argument('--first_degree_pi_hat', help='minimum pi_hat threshold to use to filter the kinship table to first degree relatives', nargs='?', const=1, type=float, default=.40)
    parser.add_argument('--grandparent_pi_hat', help='minimum pi_hat threshold to use to filter the kinship table to grandparents', nargs='?', const=1, type=float, default=.20)
    parser.add_argument('--grandparent_ibd1', help='minimum ibd1 threshold to use to filter the kinship table to grandparents', nargs='?', const=1, type=float, default=.25)
    parser.add_argument('--grandparent_ibd2', help='maximum ibd2 threshold to use to filter the kinship table to grandparents', nargs='?', const=1, type=float, default=.15)
    parser.add_argument('--ibd1_second_degree_threshold', help='ibd1 value below which second degree relatives will be called unrelated', nargs='?', const=1, type=float, default=.40)
    parser.add_argument('--ibd0_second_degree_threshold', help='ibd1 value below which second degree relatives will be called unrelated', nargs='?', const=1, type=float, default=.40)

    parser.add_argument('--second_degree_min_kin', help='min kinship threshold for 2nd degree relatives', nargs='?', const=1, type=float, default=0.195)
    parser.add_argument('--ibd1_100_min', help='min IBD1 threshold for 1.0 IBD1 sharing', nargs='?', const=1, type=float, default=0.70)
    parser.add_argument('--ibd2_0_max', help='max IBD2 threshold for 0 IBD2 sharing', nargs='?', const=1, type=float, default=0.30)
    parser.add_argument('--ibd0_0_max', help='max IBD0 threshold for 0 IBD0 sharing', nargs='?', const=1, type=float, default=0.15)

    parser.add_argument('--first_degree_kin_thresholds_lower', help='lower bound for kinship threshold for 1st degree relatives', nargs='?', const=1, type=float, default=0.40)
    parser.add_argument('--first_degree_kin_thresholds_upper', help='upper bound for kinship threshold for 1st degree relatives', nargs='?', const=1, type=float, default=0.75)



    args = parser.parse_args()
    
    main(args)


