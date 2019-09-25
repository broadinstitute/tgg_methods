
# MyoSeq reports

This series of scripts prepares PDF reports for the MyoSeq group led by Volker Straub.


## Part 1: Coverage
The scripts necessary for part 1 are in the `coverage` folder. This step requires access to the sample cram files.

 1. **`prepare_coverage.sh`**: Prepares batch of samples for coverage calculations across MyoSeq gene list (169 genes). Outputs a batch file for `dsub`.
 2. **`coverage.sh`**: Script to use when running `dsub` for coverage calculations. Calculates sample coverage across regions specified in batch file produced in step 1. Outputs one bgzipped tsv file per gene, with coverage calculated across all samples specified in step 1.
 3. **`gene_coverage.py`**: Script to use when running `dsub` to summarize coverage calculations. (This might also work on a VM?) Creates one tsv per gene with mean coverage, median coverage, and coverage above each threshold (1, 5, 10, 15, 20, 25, 30, 50, 100) for each position.
 4. **`summarize_coverage.py`**: Summarizes sample coverage over each MyoSeq gene. Outputs bgzipped tsv with gene name, gene size, mean coverage, median coverage, % callable bases, and uncallable bases.

## Part 2: Ancestry and sexcheck
The scripts necessary for part 2 are in the `inference` folder. This step requires three files: one with inferred ancestry, one with reported sex, and another with imputed sex. The reported ancestry is always recorded as European.

The two files with inferred information (`seqr_sample_qc.tsv` and `sex.txt`) are generated via the seqr sample QC pipeline. 

 1. **`ancestry.py`**: Parses `seqr_sample_qc.tsv` for each sample's inferred ancestry. Outputs a file with two columns: sample ID and ancestry.
 2. **`sex_check.py`**: Parses `sex.txt` for each sample's inferred sex and pedigree file (downloaded from *seqr*) to compare inferred and reported sex. Outputs a file with four columns: sample ID, reported sex, inferred sex, f-stat, and whether the sexes match (`CONCORD`/`CONFLICT`).

## Part 3: Candidate variant wrangling
The scripts necessary for part 3 are in the `seqr` folder. This step requires TSV files with REPORT and candidate variants downloaded from [_seqr_]([https://seqr.broadinstitute.org/dashboard](https://seqr.broadinstitute.org/dashboard)).

1. **`bigquery.py`**: Prepares TSVs downloaded from _seqr_ for upload to BigQuery. Variants are uploaded to BigQuery to lookup [gnomAD]([http://gnomad.broadinstitute.org](http://gnomad.broadinstitute.org/)) popmax allele frequency/popmax population; _seqr_ does not store this information.
2. **`reformat.py`**: Combines JSON files with frequency downloaded from BigQuery and TSVs downloaded from seqr for report generation.

## Part 4: Report generation
The scripts necessary for part 4 are in the top level `MyoSeq_reports` folder. 

1. **`make_myoseq_report.py`**: Generates .tex file for one sample at a time.
2. **`make_myoseq_report.sh`**: Calls `make_myoseq_report.py` for all samples to generate all .tex files.
3. **`make_pdfs.sh`**: Uses [`pdflatex`]([https://www.tug.org/applications/pdftex/](https://www.tug.org/applications/pdftex/)) to generate PDF files for all patients.

## Resources
This folder contains .tex files necessary to properly format reports.
