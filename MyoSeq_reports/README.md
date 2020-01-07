
# MyoSeq reports

This series of scripts prepares PDF reports for the MyoSeq group led by Volker Straub.

## Google cloud set up
 1. Download and install the Google cloud [sdk](https://cloud.google.com/sdk/). Make sure you have a billing project set up.
 2. Create a Google bucket. 
 3. Create a VM instance. Make sure the boot disk is large enough to store all of the crams. Example command: `gcloud beta compute --project=cmg-analysis instances create coverage --zone=us-central1-b --machine-type=n1-standard-64 --subnet=default --network-tier=PREMIUM --maintenance-policy=MIGRATE --service-account=422849809944-compute@developer.gserviceaccount.com --scopes=https://www.googleapis.com/auth/devstorage.read_only,https://www.googleapis.com/auth/logging.write,https://www.googleapis.com/auth/monitoring.write,https://www.googleapis.com/auth/servicecontrol,https://www.googleapis.com/auth/service.management.readonly,https://www.googleapis.com/auth/trace.append --image=debian-9-stretch-v20191210 --image-project=debian-cloud --boot-disk-size=500GB --boot-disk-type=pd-standard --boot-disk-device-name=coverage --reservation-affinity=any`. Be sure to shut down the VM when not in use, as Google will charge for its uptime.
 4. Install samtools, bgzip, and their necessary dependencies on your VM instance. See [samtools documentation] (http://www.htslib.org/download/).

## Part 1: Preparation
 1. Create BED files for all MyoSeq gene lists if necessary. For instructions, see README in `beds/` directory. Copy BED files to the Google bucket created above. 
 2. Create lists of samples with candidate genes, without candidate genes, and with candidate CNV or SMA findings.
 3. Download all variants tagged "REPORT" in seqr (saved variants page), and all rare variants across MyoSeq gene list (use project-wide search).
 4. Locate all cram files for samples of interest and copy them to the Google bucket.

## Part 2: Coverage
The scripts necessary for part 1 are in the `coverage` folder. This step requires access to the sample cram files. All scripts in this section are run on the Google VM created above.

 1. ssh into VM instance. Example command: `gcloud beta compute --project "cmg-analysis" ssh --zone "us-central1-b" "coverage"`
 2. **`get_coverage.sh`**: Calculates sample coverage across all MyoSeq gene list genes. Outputs a bgzipped TSV per gene.
 3. **`gene_coverage.py`**: Script to use when running `dsub` to summarize coverage calculations. (This might also work on a VM?) Creates one tsv per gene with mean coverage, median coverage, and coverage above each threshold (1, 5, 10, 15, 20, 25, 30, 50, 100) for each position.
 4. **`summarize_coverage.py`**: Summarizes sample coverage over each MyoSeq gene. Outputs bgzipped tsv with gene name, gene size, mean coverage, median coverage, % callable bases, and uncallable bases.

## Part 3: Ancestry and sexcheck
The scripts necessary for part 2 are in the `inference` folder. This step requires three files: one with inferred ancestry, one with reported sex, and another with imputed sex. The reported ancestry is always recorded as European.

The two files with inferred information (`seqr_sample_qc.tsv` and `sex.txt`) are generated via the seqr sample QC pipeline. 

 1. **`ancestry.py`**: Parses `seqr_sample_qc.tsv` for each sample's inferred ancestry. Outputs a file with two columns: sample ID and ancestry.
 2. **`sex_check.py`**: Parses `sex.txt` for each sample's inferred sex and pedigree file (downloaded from *seqr*) to compare inferred and reported sex. Outputs a file with four columns: sample ID, reported sex, inferred sex, f-stat, and whether the sexes match (`CONCORD`/`CONFLICT`).

## Part 4: Candidate variant wrangling
The scripts necessary for part 3 are in the `seqr` folder. This step requires TSV files with REPORT and candidate variants downloaded from [_seqr_]([https://seqr.broadinstitute.org/dashboard](https://seqr.broadinstitute.org/dashboard)).

1. **`bigquery.py`**: Prepares TSVs downloaded from _seqr_ for upload to BigQuery. Variants are uploaded to BigQuery to lookup [gnomAD]([http://gnomad.broadinstitute.org](http://gnomad.broadinstitute.org/)) popmax allele frequency/popmax population; _seqr_ does not store this information.
2. **`reformat.py`**: Combines JSON files with frequency downloaded from BigQuery and TSVs downloaded from seqr for report generation.

## Part 5: Report generation
The scripts necessary for part 4 are in the top level `MyoSeq_reports` folder. 

1. **`make_myoseq_report.py`**: Generates .tex file for one sample at a time.
2. **`make_myoseq_report.sh`**: Calls `make_myoseq_report.py` for all samples to generate all .tex files.
3. **`make_pdfs.sh`**: Uses [`pdflatex`]([https://www.tug.org/applications/pdftex/](https://www.tug.org/applications/pdftex/)) to generate PDF files for all patients.

## Resources
This folder contains .tex files necessary to properly format reports.
