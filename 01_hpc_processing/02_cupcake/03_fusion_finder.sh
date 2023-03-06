#! /bin/bash

###############
# DESCRIPTION #
###############
# Runs Cupcake fusion finder script.
#
#############
# VARIABLES #
#############
# FUSION_FINDER_PY=[path to `fusion_finder.py` executable]
# HQ_ISOFORMS=[path to <prefix>hq.fasta.gz]
# SORTED_SAM_FILE=[path to file output from step 02]
# CLUSTER_REPORT=[path to <prefix>.cluster_report.csv from step 01]
# OUTPUT_PREFIX

cp $HQ_ISOFORMS .
gunzip $HQ_ISOFORMS
unzipped=${HQ_ISOFORMS%".gz"}
$FUSION_FINDER_PY --input $unzipped -s $SORTED_SAM_FILE --cluster_report $CLUSTER_REPORT -o $OUTPUT_PREFIX -d 100000
