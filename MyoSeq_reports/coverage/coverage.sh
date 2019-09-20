#!/bin/bash

# gs://broad-ukbb/resources/Homo_sapiens_assembly38.fasta.bgz

# old commands from Monkol
# ~mlek/dev/samtools/samtools depth -b $bedfile -f myoseq.bam.list -q 10 -Q 20 -r $region | bgzip > myoseq.$gene.txt.gz
# ~mlek/dev/samtools/samtools depth -b ACTA1.bed -f myoseq.bam.list -q 10 -Q 20 -r 1:229567238-229568872 | bgzip > myoseq.ACTA1.txt.gz
# ACTA1.bed
#1	229567238	229567399	ACTA1_ENSG00000143632.10_ENST00000366684.3
#1	229567457	229567659	ACTA1_ENSG00000143632.10_ENST00000366684.3
#1	229567730	229567942	ACTA1_ENSG00000143632.10_ENST00000366684.3
#1	229568006	229568188	ACTA1_ENSG00000143632.10_ENST00000366684.3
#1	229568292	229568637	ACTA1_ENSG00000143632.10_ENST00000366684.3
#1	229568723	229568872	ACTA1_ENSG00000143632.10_ENST00000366684.3

# TSV should be named gene.txt.bgz
samtools depth -f ${BAMLIST} -r ${REGION} -q 10 -Q 20 -a --reference ${FASTA} | bgzip > ${TSV}
