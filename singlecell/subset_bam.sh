#!/usr/bin/bash

# Ref: https://bioinformatics.stackexchange.com/questions/12868/is-there-a-command-line-tool-to-split-a-sam-bam-file-by-cb-cell-barcode-tag

# $1 is the raw bam file
# $2 is the barcode list (one barcode per line) to filter
# $3 output bam file
# $4 temp dir

tmpdir=$4
temp=$(mktemp -p ${tmpdir} -q -u)
bamhead="${temp}_head"
filterSamBody="${temp}_filter_sam_body"
filterSam="${filterSamBody}.sam"

# Save the header lines
samtools view -H $1 > $bamhead

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $1 | LC_ALL=C grep -F -f $2 > ${filterSamBody}

# Combine header and body
cat ${bamhead} ${filterSamBody} > ${filterSam}

# Convert filtered.sam to BAM format
samtools view -b ${filterSam} > $3

# remove temp files
rm $bamhead $filterSamBody $filterSam

