#!/bin/bash

# Pre-requisites / dependencies: wget, samtools, R

# download the Nanopore cfDNA dataset from zenodo (~8Gb) and extract content
wget -c https://zenodo.org/record/6642503/files/Katsman-data-files.140622.megalodonModBams.tar.gz

# check md5sum
echo "ffe9916a38fead7c6db78ea1003ed89f Katsman-data-files.140622.megalodonModBams.tar.gz" > chk.tmp
md5sum -c chk.tmp
rm chk.tmp

# unpack the tar ball
tar -xvzf Katsman-data-files.140622.megalodonModBams.tar.gz

# for each bam file, get the length of each read from the absolute value in TLEN field
for IN_FILE in $(ls megalodonModBams/*bam)
do
    OUT_FILE=$(basename ${IN_FILE} .meg242.remora1.edgefix.mod_mappings.sorted.hg38.bam)_frglen.txt
    samtools view -F 3852 ${IN_FILE} \
    | awk '{contig = $3=="chrM" ? "chrM" : "gDNA"; print contig, sqrt($9**2)}' \
    | sort -k1,1 -k2,2n \
    | uniq -c \
    | awk 'BEGIN{OFS="\t"; print "freq\tchr\tlength"} {print $1, $2, $3}' \
    > ${OUT_FILE}
done

# then, manually retrieve sample metadata
# go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185307
# note the library names and the sample names
# click SRA Run Selector
# download Metadata for all samples as a csv file
# match sample names with "disease_state" using the "Library Name" column
# concatenate all SAMPLE_frglen.txt into 1 csv file named "mtDNA_modeling_data.csv"
