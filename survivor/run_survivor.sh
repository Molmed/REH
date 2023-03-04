#!/usr/bin/env bash

# Expected binary file in "bin":
# SURVIVOR Version: 1.0.7
# https://github.com/fritzsedlazeck/SURVIVOR
SURVIVOR="bin/SURVIVOR"

# Paths for docker
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DOCKER_MNT="/REH"

# Path to binaries
BCFTOOLS="docker run -v $SCRIPT_DIR:$DOCKER_MNT themariya/bcftools:latest bcftools"
CONDA_BIN_DIR="/home/mariya/miniconda3/envs/reh/bin/"
BGZIP="$CONDA_BIN_DIR/bgzip"
TABIX="$CONDA_BIN_DIR/tabix"

# Expected input in "data" directory:
# - illumina.tiddit.vcf
# - ont.sniffles.vcf
# - pb.sniffles.vcf

MAX_DISTANCE=1000
MIN_BP=100
SMALL=1000
MEDIUM=10000
LARGE=1000000

OUT_DIR="out"
DATA_DIR="data"
TMP_DIR="/tmp"
BLACKLIST="ref/hg38-blacklist.v2.bed"

VCF_LIST="$OUT_DIR/vcfs"
echo "$DATA_DIR/illumina.tiddit.vcf" >> $VCF_LIST
echo "$DATA_DIR/ont.sniffles.vcf" >> $VCF_LIST
echo "$DATA_DIR/pb.sniffles.vcf" >> $VCF_LIST

### MERGE ALL SVS ###
$SURVIVOR merge $VCF_LIST $MAX_DISTANCE 1 1 1 0 $MIN_BP "$OUT_DIR/REH_merged.tmp.vcf"
$SURVIVOR filter "$OUT_DIR/REH_merged.tmp.vcf" $BLACKLIST $MIN_BP -1 0 -1 "$OUT_DIR/REH_merged_all.vcf"

### CONSENSUS SETS ###
$SURVIVOR merge $VCF_LIST $MAX_DISTANCE 2 1 1 0 $MIN_BP "$OUT_DIR/REH_consensus_2x.vcf" 
$SURVIVOR merge $VCF_LIST $MAX_DISTANCE 3 1 1 0 $MIN_BP "$OUT_DIR/REH_consensus_3x.vcf"

### SORT, COMPRESS AND INDEX CONSENSUS SETS ###
$BCFTOOLS sort "$DOCKER_MNT/$OUT_DIR/REH_consensus_2x.vcf" -T $TMP_DIR -Oz -o "$DOCKER_MNT/$OUT_DIR/REH_consensus_2x.sorted.vcf"
$BGZIP "$OUT_DIR/REH_consensus_2x.sorted.vcf"
$TABIX -p vcf "$OUT_DIR/REH_consensus_2x.sorted.vcf.gz"

$BCFTOOLS sort "$DOCKER_MNT/$OUT_DIR/REH_consensus_3x.vcf" -T $TMP_DIR -Oz -o "$DOCKER_MNT/$OUT_DIR/REH_consensus_3x.sorted.vcf"
$BGZIP "$OUT_DIR/REH_consensus_3x.sorted.vcf"
$TABIX -p vcf "$OUT_DIR/REH_consensus_3x.sorted.vcf.gz"

### NON-TRANSLOCATION STATISTICS ###
# Strip out BNDs
grep -v "SVTYPE=BND" "$DATA_DIR/illumina.tiddit.vcf" > "$OUT_DIR/illumina.nobnd.vcf"
grep -v "SVTYPE=BND" "$DATA_DIR/ont.sniffles.vcf" > "$OUT_DIR/ont.nobnd.vcf"
grep -v "SVTYPE=BND" "$DATA_DIR/pb.sniffles.vcf" > "$OUT_DIR/pb.nobnd.vcf"

echo "$OUT_DIR/illumina.nobnd.vcf" > $VCF_LIST
echo "$OUT_DIR/ont.nobnd.vcf" >> $VCF_LIST
echo "$OUT_DIR/pb.nobnd.vcf" >> $VCF_LIST

# Merge Non-translocations only
$SURVIVOR merge $VCF_LIST $MAX_DISTANCE 1 1 1 0 $MIN_BP "$OUT_DIR/REH_merged_nobnd.tmp.vcf"
$SURVIVOR filter "$OUT_DIR/REH_merged_nobnd.tmp.vcf" $BLACKLIST $MIN_BP -1 0 -1 "$OUT_DIR/REH_merged_nobnd.vcf"

# Merged non-translocations by size
$SURVIVOR filter "$OUT_DIR/REH_merged_nobnd.vcf" $BLACKLIST $MIN_BP $SMALL 0 -1 "$OUT_DIR/REH_merged_small.vcf"
$SURVIVOR filter "$OUT_DIR/REH_merged_nobnd.vcf" $BLACKLIST $SMALL $MEDIUM  0 -1 "$OUT_DIR/REH_merged_medium.vcf"
$SURVIVOR filter "$OUT_DIR/REH_merged_nobnd.vcf" $BLACKLIST $MEDIUM -1 0 -1 "$OUT_DIR/REH_merged_large.vcf"

### MATRICES FOR VENN DIAGRAMS ###
NORMALIZE=0
# All SVs
$SURVIVOR genComp "$OUT_DIR/REH_merged_all.vcf" $NORMALIZE "$OUT_DIR/REH_merged_all.mat.txt"
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' "$OUT_DIR/REH_merged_all.vcf" | sed -e 's/\(.\)/\1 /g' > "$OUT_DIR/REH_merged_all_overlap.txt"

# Small
$SURVIVOR genComp "$OUT_DIR/REH_merged_small.vcf" $NORMALIZE "$OUT_DIR/REH_merged_small.mat.txt"
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' "$OUT_DIR/REH_merged_small.vcf" | sed -e 's/\(.\)/\1 /g' > "$OUT_DIR/REH_merged_small_overlap.txt"

# Medium
$SURVIVOR genComp "$OUT_DIR/REH_merged_medium.vcf" $NORMALIZE "$OUT_DIR/REH_merged_medium.mat.txt"
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' "$OUT_DIR/REH_merged_medium.vcf" | sed -e 's/\(.\)/\1 /g' > "$OUT_DIR/REH_merged_medium_overlap.txt"

# Large
$SURVIVOR genComp "$OUT_DIR/REH_merged_large.vcf" $NORMALIZE  "$OUT_DIR/REH_merged_large.mat.txt"
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' "$OUT_DIR/REH_merged_large.vcf" | sed -e 's/\(.\)/\1 /g' > "$OUT_DIR/REH_merged_large_overlap.txt"
