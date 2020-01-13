
# MyoSeq reports

This series of scripts prepares PDF reports for the MyoSeq group led by Volker Straub.

## Google cloud set up
 1. Download and install the Google cloud [sdk](https://cloud.google.com/sdk/). Make sure you have a billing project set up.
 2. Create a Google bucket. 
 3. Create a VM instance. Make sure the boot disk is large enough to store all of the crams. Example command: `gcloud beta compute --project=cmg-analysis instances create coverage --zone=us-central1-b --machine-type=n1-standard-64 --subnet=default --network-tier=PREMIUM --maintenance-policy=MIGRATE --service-account=422849809944-compute@developer.gserviceaccount.com --scopes=https://www.googleapis.com/auth/devstorage.read_only,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/trace.append --image=debian-9-stretch-v20191210 --image-project=debian-cloud --boot-disk-size=500GB --boot-disk-type=pd-standard --boot-disk-device-name=coverage --reservation-affinity=any`. Be sure to shut down the VM when not in use, as Google will charge for its uptime.
 4. Install samtools, bgzip, Anaconda (Python3.7), and their necessary dependencies on your VM instance. See [samtools documentation](http://www.htslib.org/download/) and [Anaconda](https://www.anaconda.com/distribution/).
 5. Clone this repo or copy the scripts in the `coverage` folder to the VM.
 6. Ensure the Google service account associated with the VM or a personal Gmail has access to the crams. To login on the VM with a personal Gmail, use `gcloud auth login`.
 7. Install (hail)[https://hail.is/docs/0.2/getting_started.html] locally.
 8. Clone the gnomAD hail utilities [repo](https://github.com/macarthur-lab/gnomad_hail). Set `HAIL_SCRIPTS` (in `~/.bashrc`) equal to the path to the gnomAD hail repo for imports to work properly.

## Directory set up
On your local machine, set up the following directory structure:
```
reports/
    seqr/
    inference/
    per_sample/
    pdf/
    summary/ (create this only if samples have SMA or CNV results)
```

## Part 1: Preparation
 1. Create BED files for all MyoSeq gene lists if necessary. For instructions, see README in `beds/` directory. Copy BED files to the Google bucket created above. 
 2. Create lists of samples with candidate genes, without candidate genes, and with candidate CNV or SMA findings.
 3. Check which _seqr_ projects contain the samples of interest.
 4. Create an analysis group for the samples of interest in their relevant project(s).
 5. Download all variants tagged "REPORT" (saved variants page). Select only lines containing samples of interest.
 6. Download all rare variants across MyoSeq gene list (use project-wide search on the analysis group(s) created in step 4 for ClinVar pathogenic/likely pathogenic, nonsense, essential splice site, frameshift, missense, or synonymous variants in the MyoSeq gene list. Filter to AF 0.01 across all reference population databases, and do not apply any quality filters).
 7. Copy all TSVs downloaded from _seqr_ to Google bucket.
 8. Download pedigree for samples from _seqr_ (in the project(s) of interest, click "Download Table", then ".tsv" under "Individuals") and remove quotations from downloaded files.
 9. Locate all cram files for samples of interest and copy them to the Google bucket.
 10. Install [pdflatex](https://www.tug.org/applications/pdftex/).

## Part 2: Coverage
The scripts necessary for part 1 are in the `coverage` folder. This step requires access to the sample cram files. All scripts in this section are run on the Google VM created above.

 1. ssh into VM instance. Example command: `gcloud beta compute --project "cmg-analysis" ssh --zone "us-central1-b" "coverage"`
 2. **`get_coverage.sh`**: Calculates sample coverage across all MyoSeq gene list genes. Outputs a bgzipped TSV per gene, with coverage per position in the gene. The columns in this file are samples in the same order as provided (i.e., in the same order as in the input crams list).
 3. **`summarize_coverage.py`**: Summarizes coverage across MyoSeq genes for each sample and compares each sample to the other samples in the batch. Creates one TSV per sample with gene, batch mean coverage, sample mean coverage, percent of callable sites for the batch, percent of callable sites for the sample, number of uncallable sites for the batch, and number of uncallable sites for the sample.
 4. Copy the files output by step 3 first to the Google bucket and then to local storage. Don't forget to shut down the VM.

## Part 3: Ancestry and sexcheck
The scripts necessary for part 2 are in the `inference` folder. This step requires three files: one with inferred ancestry, one with reported sex, and another with imputed sex. The reported ancestry is always recorded as European. Note that `make_myoseq_report.py` assumes the output files are named `MYOSEQ_sex.tsv` and `MYOSEQ_pop.tsv` and stored in `inference/`

The two files with inferred information (`seqr_sample_qc.tsv` and `sex.txt`) are generated via the _seqr_ sample QC pipeline. 

 1. **`ancestry.py`**: Parses `seqr_sample_qc.tsv` for each sample's inferred ancestry. Outputs a file with two columns: sample ID and ancestry.
 2. **`sex_check.py`**: Parses `sex.txt` for each sample's inferred sex and pedigree file (downloaded from _seqr_) to compare inferred and reported sex. Outputs a file with four columns: sample ID, reported sex, inferred sex, f-stat, and whether the sexes match (`CONCORD`/`CONFLICT`).

## Part 4: Candidate variant wrangling
The scripts necessary for part 3 are in the `seqr` folder. This step requires TSV files with REPORT and all MyoSeq gene list variants downloaded from [_seqr_]([https://seqr.broadinstitute.org/dashboard](https://seqr.broadinstitute.org/dashboard)). Note that `make_myoseq_report.py` assumes the output files are stored in `seqr/variants/`

 1. Using `hailctl` (or preferred method), spin up a Google compute cluster to run the popmax script. Don't forget to shut down the cluster when popmax is complete. Example command: `hailctl dataproc start kc --master-machine-type n1-highmem-8 --worker-machine-type n1-highmem-8 --num-workers 10 --init gs://gnomad-public/tools/inits/master-init.sh --max-idle 20m --worker-boot-disk-size=100 --project cmg-analysis --properties=spark:spark.executor-memory=25g,spark:spark.speculation=true,spark:spark.speculation.quantile=0.9,spark:spark.speculation.multiplier=3`
 2. **`get_popmax.py`**: Joins seqr downloaded TSVs to gnomAD _hail_ tables to get gnomAD global allele frequency (AF), highest population AF in gnomAD exomes, and population with highest AF in gnomAD exomes.


## Part 5: Report generation
The scripts necessary for part 4 are in the top level `MyoSeq_reports` folder. 

1. **`make_myoseq_report.py`**: Generates .tex file for one sample at a time.
2. **`make_myoseq_report.sh`**: Calls `make_myoseq_report.py` for all samples to generate all .tex files. Example command: `bash ~/code/methods/MyoSeq_reports/pdf/make_myoseq_report.sh -l sample_lists/candidates.list -d ~/Documents/MYOSEQ/Reports/Nov_2019 -r ~/code/methods/MyoSeq_reports/resources/tex/ -o pdf/ -p ~/code/methods/MyoSeq_reports/pdf/`
3. **`make_pdfs.sh`**: Uses [`pdflatex`]([https://www.tug.org/applications/pdftex/](https://www.tug.org/applications/pdftex/)) to generate PDF files for all patients. **NOTE:** Run this twice to ensure proper formatting.

## Resources
This folder contains `.tex` and `.bed` files necessary to properly format reports.
