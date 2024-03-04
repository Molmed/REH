#!/usr/bin/env bash

# Expected input in "../../data/sv_callsets" directory:
# - illumina.tiddit.vcf
# - ont.sniffles.vcf
# - pb.sniffles.vcf

# PATHS TO BINARIES
# Expected binary file: SURVIVOR Version: 1.0.7
# https://github.com/fritzsedlazeck/SURVIVOR
SURVIVOR="../bin/SURVIVOR"

# PATHS
OUT_DIR="out"
DATA_DIR="../../data/sv_callsets"
TMP_DIR="/tmp"

# SETTINGS
MAX_DISTANCE=1000
MIN_BP=100

### GENERATE CONSENSUS SET ###
VCF_LIST="$OUT_DIR/vcfs"
echo "$DATA_DIR/illumina.tiddit.vcf" > $VCF_LIST
echo "$DATA_DIR/ont.sniffles.vcf" >> $VCF_LIST
echo "$DATA_DIR/pb.sniffles.vcf" >> $VCF_LIST

$SURVIVOR merge $VCF_LIST $MAX_DISTANCE 3 1 1 0 $MIN_BP "$OUT_DIR/REH_consensus_3x.vcf"

# GENERATE LONG-READ CONSENSUS SET ###
VCF_LIST="$OUT_DIR/longread_vcfs"
echo "$DATA_DIR/ont.sniffles.vcf" > $VCF_LIST
echo "$DATA_DIR/pb.sniffles.vcf" >> $VCF_LIST

$SURVIVOR merge $VCF_LIST $MAX_DISTANCE 2 1 1 0 $MIN_BP "$OUT_DIR/REH_consensus_longread.vcf"
