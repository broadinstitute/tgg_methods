import argparse
import json
import logging
import pickle

import hail as hl
from gnomad.sample_qc.ancestry import assign_population_pcs, pc_project
from gnomad.sample_qc.filtering import compute_stratified_metrics_filter
from gnomad.sample_qc.pipeline import filter_rows_for_qc
from gnomad.sample_qc.platform import (
    assign_platform_from_pcs,
    compute_callrate_mt,
    run_platform_pca,
)
from gnomad.utils import slack
from gnomad.utils.filtering import filter_to_autosomes
from hail.utils.misc import new_temp_file
from resources.resources_seqr_qc import (  # missing_metrics_path,; mt_path,; remap_path,; sample_qc_ht_path,; sample_qc_tsv_path,; seq_metrics_path,
    VCFDataTypeError,
    rdg_gnomad_pop_pca_loadings_path,
    rdg_gnomad_rf_model_path,
    rdg_gnomad_v4_pop_pca_loadings_path,
    rdg_gnomad_v4_rf_model_path,
    val_coding_ht_path,
    val_noncoding_ht_path,
)

# from gnomad_qc.v4.resources.sample_qc import per_pop_min_rf_probs_json_path
from gnomad_qc.v4.sample_qc.assign_ancestry import assign_pop_with_per_pop_probs

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("seqr_rf_spike_generation")
logger.setLevel(logging.INFO)

V4_LOADINGS = "gs://gcp-public-data--gnomad/release/4.0/pca/gnomad.v4.0.pca_loadings.ht"
# V2_LOADINGS = "gs://seqr-datasets/sample_qc_resources/population_assignment/38_ancestry_pca_loadings.ht"
GLEESON_PEDIGREE_PATH = 'gs://seqr-datasets/v02/GRCh38/RDG_WES_Broad_Internal/base/projects/R0285_cmg_gleeson_exomes/R0285_cmg_gleeson_exomes_pedigree.tsv'
V2_POP_PATH =  'gs://marten-seqr-sandbox-storage/ancestry/wes_v23_gnomadv2_qc_pop_ht.tsv'

def my_rf_spike_pop(
    mt: hl.MatrixTable,
    ht_pop_path: str = V2_POP_PATH,
    pop_col: str = 'v2_qc_pop',
    pedigree_path: str = GLEESON_PEDIGREE_PATH,
    pop: str = 'mid',
) -> hl.MatrixTable: 
    
    ht_pedigree = hl.import_table(pedigree_path).key_by('Individual_ID')
    
    # Filter to GLE samples
    mt_subset = mt.filter_cols(
        hl.is_defined(ht_pedigree[mt.col_key])
    )
    
    # Filter to MID samples
    if '.tsv' in ht_pop_path: 
        ht_pop = hl.import_table(ht_pop_path)
    else: 
        ht_pop = hl.read_table(ht_pop_path)
        
    ht_pop = ht_pop.key_by('s')
    
    mt_subset = mt_subset.annotate_cols(
        qc_pop = ht_pop[mt_subset.col_key][f"{pop_col}"]
    )
        
    mt_subset = mt_subset.filter_cols(
        mt_subset.qc_pop == pop
    )
    
    return mt_subset

def do_v4_pc_project(
    mt: hl.MatrixTable,
    build: int,
    ht_pop_path=V2_POP_PATH,
    pop_col='v2_qc_pop',
    pop='mid',
    pedigree_path=GLEESON_PEDIGREE_PATH,
    num_pcs: int = 20,
    v2_loadings: bool = False
) -> hl.Table:
    
    loadings_path = V4_LOADINGS
    if v2_loadings: 
        loadings_path = rdg_gnomad_pop_pca_loadings_path(build=38)

    loadings = hl.read_table(loadings_path)
    
    logger.info(f"With num_pcs: {num_pcs}")
    
    mt = my_rf_spike_pop(mt,
                        ht_pop_path=ht_pop_path,
                        pop_col=pop_col,
                        pedigree_path=pedigree_path,
                        pop=pop)
    
    mt = mt.select_entries("GT")
    
    logger.info(f"With num_pcs: {num_pcs}")
    
    scores = pc_project(mt, loadings)
    scores = (
        scores.annotate(scores=scores.scores[:num_pcs], known_pop="Unknown")
        .key_by("s")
        .checkpoint(new_temp_file("scores_temp", extension="ht"))
    )
    
    return scores
    
def main(args):
    
    mt = hl.read_matrix_table(
       args.callset_path
    )  
    
    if args.test:
        mt = mt.sample_cols(0.01)
        mt = mt.filter_rows(mt.locus.contig=="chr22")
    
    build = args.build
    ht_pop_path= args.ht_pop_path
    pop_col=args.pop_col
    training_pop=args.training_pop
    pedigree_path=args.pedigree_path
    num_pcs = args.num_pcs
    v2_loadings = args.v2_loadings
    output_bucket = args.bucket_path 

    
    scores_output = do_v4_pc_project(
        mt = mt,
        build = build,
        ht_pop_path=ht_pop_path,
        pop_col=pop_col,
        pop=training_pop,
        pedigree_path=pedigree_path,
        num_pcs = num_pcs,
        v2_loadings = v2_loadings,
        
    )
    
    temp = new_temp_file().split('/')[-1]

    out_path = f'{output_bucket}_custom_loadings_{temp}.ht'

    logger.info(f'Printing to {out_path}')
    
    scores_output = scores_output.checkpoint(f'{output_bucket}_custom_loadings_{training_pop}_{temp}.ht',overwrite=True)
    
    scores_output_meta = scores_output.annotate(
        v4_race_ethnicity = f'v2_{training_pop}_training',
        data_type="exomes",
        hard_filters = hl.empty_set(hl.tstr),
        broad_external="broad",
        releasable=True)
    
    scores_output_meta = scores_output_meta.checkpoint(f'{output_bucket}custom_loadings_meta_{training_pop}_{temp}.ht',overwrite=True)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--callset-path",
        help="Path to callset as mt",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--pedigree-path",
        help="Path to seqr project pedigree tsv tsv",
        type=str,
        default=GLEESON_PEDIGREE_PATH
    )
    parser.add_argument(
        "--ht-pop-path",
        help="Path to hail table or tsv containing genetic ancestry assignments",
        type=str,
        default=V2_POP_PATH
    )

    parser.add_argument(
        "-b",
        "--build",
        help="Reference build, 37 or 38",
        type=int,
        choices=[37, 38],
        default=38,
    )

    parser.add_argument(
        "--test",
        help="To run a test of the pipeline using test files and directories",
        action="store_true",
    )
    parser.add_argument(
        "--v2-loadings",
        help="Choose to use gnomAD v2 loadings to genetate PCs",
        action="store_true"
    )
    parser.add_argument(
        "--pop-col",
        help="Column containing genetic ancestry group information",
        type=str,
        default="v2_qc_pop"
    )
    parser.add_argument(
        "--training_pop",
        help="Genetic ancestry group to create extra training set for",
        type=str,
        default="mid"
    )
    parser.add_argument(
        "--num-pcs",
        help="Number of PCs to project.",
        type=int,
        default=20
    )

    parser.add_argument(
        "--overwrite", help="Overwrite previous paths", action="store_true"
    )
    parser.add_argument(
        "--bucket-path",
        help="Path to bucket to output results to",
        default="gs://marten-seqr-sandbox-tmp-4day/"
    )




    args = parser.parse_args()

    main(args)
