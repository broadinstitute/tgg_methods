import collections
import json
import os


#%%
files = """gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_10_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_11_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_12_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_13_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_13_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_14_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_15_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_15_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_15_3.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_16_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_17_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_18_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_18_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_19_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_1_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_1_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_20_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_21_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_21_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_22_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_22_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_23_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_24_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_24_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_24_3.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_24_4.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_25_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_26_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_26_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_27_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_28_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_29_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_29_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_29_3.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_29_4.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_2_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_30_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_30_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_31_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_31_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_31_3.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_3_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_3_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_3_3.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_4_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_4_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_4_3.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_4_4.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_5_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_5_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_6_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_6_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_7_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_7_2.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_7_3.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_8_1.dcr.bed.gz
gs://fc-secure-e2c5f2a5-2e76-4c01-a264-419262b2c7c8/dcr_tabs/cc_9_1.dcr.bed.gz"""

path_lookup = {f.split("/")[-1].replace(".dcr.bed.gz", ""): f for f in files.split()}


batch_to_samples = collections.defaultdict(list)
for path in files.split():
    for line in open("batch_info.txt"):
        sample_id, batch_name = line.strip().split()
        batch_to_samples[batch_name].append(sample_id)

        # check that there's a .bed file for this batch name
        path = path_lookup[batch_name]

rows = []
for batch_name, sample_ids in batch_to_samples.items():
    rows.append({
        'name': batch_name,
        'data': [{
            'type': 'gcnv_bed',
            'url': path_lookup[batch_name],
            'samples': sample_ids,
        }],
    })

    #header = f.readline()
    #fields = header.split()
    #print(fields)

#%%
# output settings

gcnv_rows = json.dumps(rows)

settings_json = """
{
    "genome": "hg38",
    "locus": "chr15:92,835,700-93,031,800",
    "bamOptions": {
	"showBams": false,
        "trackHeight": 200,
        "viewAsPairs": false,
        "showSoftClips": true,
        "alignmentShading": "strand"
    },
    "sjOptions": {
        "showCoverage": true,
        "showJunctions": true,
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
	"showVcfs": false,
        "displayMode": "EXPANDED"
    },
    "rowsInCategories": [
        {
            "categoryName": "Reference Data",
            "rows": []
        },
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
