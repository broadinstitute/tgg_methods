# MyoSeq reports

Stuff


## Part 1: Coverage

 1. **`prepare_coverage.sh`**: Prepares batch of samples for coverage calculations across MyoSeq gene list (169 genes). Outputs a batch file for `dsub`.
 2. **`coverage.sh`**: Script to use when running `dsub` for coverage calculations. Calculates sample coverage across regions specified in batch file produced in step 1. Outputs one bgzipped tsv file per gene, with coverage calculated across all samples specified in step 1.
 3. **`gene_coverage.py`**: Script to use when running `dsub` to summarize coverage calculations. (This might also work on a VM?) Creates one tsv per gene with mean coverage, median coverage, and coverage above each threshold (1, 5, 10, 15, 20, 25, 30, 50, 100) for each position.
 4. **`summarize_coverage.py`**: Summarizes sample coverage over each MyoSeq gene. Outputs bgzipped tsv with gene name, gene size, mean coverage, median coverage, % callable bases, and uncallable bases.

## Part 2: Ancestry and sexcheck
Mike generates a file called something like `seqr_sample_qc.tsv` as part of his seqr sample QC pipeline, and Kristen generates a file called something like `v(callset version)_sex.txt` with imputed sex.

 1. **`ancestry.py`**: Parses `seqr_sample_qc.tsv` for each sample's inferred ancestry. Outputs a file with two columns: sample ID and ancestry.
 2. **`sex_check.py`**: Parses `sex.txt` for each sample's inferred sex and pedigree file (downloaded from *seqr*) to compare inferred and reported sex. Outputs a file with four columns: sample ID, reported sex, inferred sex, f-stat, and whether the sexes match (`CONCORD`/`CONFLICT`).
