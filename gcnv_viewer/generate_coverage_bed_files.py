from collections import defaultdict
from firecloud import api
import hail as hl
import os
import pandas as pd
import re
import tqdm

hl.init(log="/dev/null")

#%%

entities = api.get_entities("cmg-exomes-gcnv", "cmg_gcnv", "sample_set").json()

#%%

# copied from https://stackabuse.com/python-how-to-flatten-list-of-lists/
def flatten(list_of_lists):
    if len(list_of_lists) == 0:
        return list_of_lists
    if isinstance(list_of_lists[0], list):
        return flatten(list_of_lists[0]) + flatten(list_of_lists[1:])
    return list_of_lists[:1] + flatten(list_of_lists[1:])

#%%


samples_counter = 0
gcnv_cluster_to_sample_bed_paths = {}

for e in entities:
    cluster_name = e["name"]
    if not cluster_name.startswith("cluster_"):
        continue
    if "attributes" not in e:
        raise ValueError(f"Entity dictionary is missing an 'attributes' key: {e}")

    attributes = e["attributes"]
    if "denoised_copy_ratios" not in attributes:
        print(f"WARNING: {cluster_name} is missing a 'denoised_copy_ratios' key in attributes dict. Skipping this cluster...")  #: {attributes}")
        continue

    denoised_copy_ratios = e["attributes"]["denoised_copy_ratios"]
    if isinstance(denoised_copy_ratios, list):
        denoised_copy_ratios = flatten(denoised_copy_ratios)
    elif isinstance(denoised_copy_ratios, dict):
        if "items" not in denoised_copy_ratios:
            raise ValueError(f"{cluster_name} is missing an 'items' key in denoised_copy_ratios dict: {denoised_copy_ratios}")
        denoised_copy_ratios = flatten(denoised_copy_ratios["items"])
    else:
        print(f"WARNING: unexpected denoised_copy_ratios type for {cluster_name}. Skipping this cluster...")  # " + str(type(denoised_copy_ratios)))
        continue

    non_string_values = [path for path in denoised_copy_ratios if not isinstance(path, str)]
    if non_string_values:
        raise ValueError(f"{cluster_name} has non-string values among denoised_copy_ratios: {non_string_values}")

    print(f"{len(denoised_copy_ratios):5d} samples in {cluster_name}")
    samples_counter += len(denoised_copy_ratios)

    gcnv_cluster_to_sample_bed_paths[cluster_name] = denoised_copy_ratios

print(f"Parsed {samples_counter} sample paths")

#%%

# group clusters like "cluster_1_CASE", "cluster_1_COHORT", "cluster_1_last_call", etc. into a single "cluster_1"
grouped_gcnv_cluster_to_sample_bed_paths = defaultdict(list)
for cluster_name, sample_bed_paths in gcnv_cluster_to_sample_bed_paths.items():
    cluster_name_pattern = "^cluster_([0-9]+)_"
    match = re.search(cluster_name_pattern, cluster_name)
    if not match:
        print(f"WARNING: {cluster_name} doesn't match {cluster_name_pattern}. It has {len(sample_bed_paths)} samples. Will keep it as a separate group.")
        grouped_cluster_name = cluster_name
    else:
        grouped_cluster_name = f"cluster_{int(match.group(1)):2d}".replace(" ", "_")

    grouped_gcnv_cluster_to_sample_bed_paths[grouped_cluster_name].extend(sample_bed_paths)

#%%

for cluster_name, sample_bed_paths in sorted(grouped_gcnv_cluster_to_sample_bed_paths.items(), key=lambda t: t[0]):
    print(f"{len(sample_bed_paths):5d} samples in {cluster_name}")
    for sample_id in "UWA_LAI11009_LAID15-1282_1", "UWA_LAI15392_LAID19-1235_1":
        if any(sample_id in path for path in sample_bed_paths):
            print(f"{sample_id} is in {cluster_name}")



#%%

def run(command):
    print(command)
    os.system(command)

#grouped_gcnv_cluster_to_sample_bed_paths

#%%


for cluster_name, paths in sorted(grouped_gcnv_cluster_to_sample_bed_paths.items(), key=lambda t: len(t[1])):  #, reverse=True):
    print(f"Processing {cluster_name} which has {len(paths)} samples")
    #for i in range(len(paths)//250):
    #    cluster_df = None
    #    for path in tqdm.tqdm(paths[i*250:(i+1)*250], unit=" paths"):
    cluster_bed_bucket_path = f"gs://seqr-datasets-gcnv/GRCh38/RDG_WES_Broad_Internal/v3/beds/{cluster_name}.bed.gz"
    if hl.hadoop_is_file(cluster_bed_bucket_path):
        print(f"{cluster_bed_bucket_path} already exists. Skipping..")
        continue

    cluster_df = None
    paths.sort(key=lambda path: os.path.basename(path))
    for path in tqdm.tqdm(paths, unit=" paths"):
            sample_name = os.path.basename(path).replace("denoised_copy_ratios-", "")
            sample_name = re.sub(".tsv$", "", sample_name)

            column_name = sample_name
            if cluster_df is not None:
                while column_name in set(cluster_df.columns):
                    print(f"WARNING: Duplicate sample name: {column_name}  {path}")
                    column_name += "_2"

            with hl.hadoop_open(path) as f:
                # skip header
                while next(f).startswith("@"):
                    continue

                # read table
                df = pd.read_table(f, names=["chr", "start", "end", column_name]).set_index(["chr", "start", "end"])
                df = df.round({column_name: 2})

            if cluster_df is None:
                cluster_df = df
            else:
                df_length = len(df)
                if len(cluster_df) != df_length:
                    raise ValueError(f"Table {path} has {df_length} rows instead of the expected {len(cluster_df)}")

                cluster_df = cluster_df.join(df, how="outer")
                if len(cluster_df) != df_length:
                    raise ValueError(f"Outer join with {path} resulted in {len(cluster_df)} rows instead of the expected {df_length}")

            #if len(cluster_df.columns) > 5:
            #    break

    print(f"Writing to {cluster_name}.bed.gz")
    cluster_df.reset_index().to_csv(f"{cluster_name}.temp.bed.gz", header=True, index=False, sep="\t", compression="gzip")
    cluster_df = None
    run(f"gunzip -c {cluster_name}.temp.bed.gz | bgzip > {cluster_name}.bed.gz")
    run(f"rm {cluster_name}.temp.bed.gz")
    run(f"tabix -f -S 1 -s 1 -b 2 -e 3 {cluster_name}.bed.gz")
    run(f"gsutil -m cp {cluster_name}.bed.gz {cluster_bed_bucket_path}")
    run(f"gsutil -m cp {cluster_name}.bed.gz.tbi {cluster_bed_bucket_path}.tbi")
    run(f"rm {cluster_name}.bed.gz {cluster_name}.bed.gz.tbi")
