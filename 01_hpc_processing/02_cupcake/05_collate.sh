#! /bin/bash

###############
# DESCRIPTION #
###############
# Collates fusion information and filters candidates.
# NOTE: this script must be run from the directory containing the output from step 03.
#
#############
# VARIABLES #
#############
# FUSION_COLLATE_INFO_PY=[path to `fusion_collate_info.py` executable]
# OUTPUT_PREFIX=[OUTPUT_PREFIX from step 03]
# CLASSIFICATION_FILE=[<OUTPUT_PREFIX>_classification.txt output from step 04]
# GENCODE_GTF=[path to gene annotation file gencode.<version>.annotation.gtf]
# REF_GENOME=[path to GRCh38.<version>.genome.fa]

sed -i -e '/^#/d' "$OUTPUT_PREFIX.abundance.txt"
$FUSION_COLLATE_INFO_PY $OUTPUT_PREFIX $CLASSIFICATION_FILE $GENCODE_GTF --genome $REF_GENOME