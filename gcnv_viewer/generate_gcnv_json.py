import collections
import hail as hl
import json
import os

#%%
os.chdir("/Users/weisburd/code/methods/gcnv_viewer")
print(os.getcwd())

#%%

google_storage_dir = "gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs"

assert hl.hadoop_is_dir(google_storage_dir)

#%%

batch_name_to_path_and_samples = {}

for result in hl.hadoop_ls(google_storage_dir):

    if not result['path'].endswith('.bed.gz') and not result['path'].endswith('.bed'):
        continue

    if result['size_bytes'] < 1000:
        print(f"ERROR: file size of {result['path']} is too small: {result['size_bytes']}")

    with hl.hadoop_open(result['path'], 'r') as f:
        line = f.readline()
        fields = line.rstrip("\n").split("\t")
        sample_ids = fields[3:]

    batch_name = os.path.basename(result['path']).replace(".dcr", "").replace(".bed", "").replace(".gz", "")
    batch_name_to_path_and_samples[batch_name] = (result['path'], sample_ids)

#%%
rows = []
for batch_name, (path, sample_ids) in batch_name_to_path_and_samples.items():
    rows.append({
        'name': batch_name,
        'data': [{
            'type': 'gcnv_bed',
            'url': path,
            'samples': sample_ids,
        }],
    })

#%%
# output settings

gcnv_rows = json.dumps(rows)

settings_json = """
{
    "genome": "hg38",
    "locus": "chr15:92,835,700-93,031,800",
    "bamOptions": {
	    "trackHeight": 200,
        "viewAsPairs": false,
        "showSoftClips": true,
        "alignmentShading": "strand"
    },
    "sjOptions": {
        "trackHeight": 170,
        "colorBy": "strand",
	    "colorByNumReadsThreshold": 5,
        "thicknessBasedOn": "numUniqueReads",
        "bounceHeightBasedOn": "random",
        "labelUniqueReadCount": true,
        "labelMultiMappedReadCount": false,
        "labelTotalReadCount": false,
        "labelMotif": false,
        "labelAnnotatedJunction": false,
        "labelAnnotatedJunctionValue": " [A]",
	    "showOnlyPlusStrand": false,
	    "showOnlyMinusStrand": false,
        "hideAnnotated": false,
        "hideUnannotated": false,
        "minUniquelyMappedReads": 0,
        "minTotalReads": 1,
        "maxFractionMultiMappedReads": 1,
        "minSplicedAlignmentOverhang": 0
    },
    "vcfOptions": {
	    "displayMode": "EXPANDED"
    },
    "rowsInCategories": [
        {
            "categoryName": "gCNV Batches",
            "rows": %(gcnv_rows)s
        }
    ]
}
""" % locals()

#%%
with open("gcnv_settings.json", "wt") as f:
    json.dump(json.loads(settings_json), f, indent=3)
