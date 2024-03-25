#!/bin/bash

#SBATCH -A snic2022-22-63
#SBATCH -t 5-00:00:00
#SBATCH -p core -n 1
#SBATCH -J seurat

module load bioinfo-tools snakemake

#snakemake $1 --cluster-config seurat.json --cluster "sbatch -A {cluster.A} -t {cluster.t} -p {cluster.p} -n {cluster.n}" -s seurat.smk.py --latency-wait 10 --use-singularity --singularity-args "\-\-cleanenv" --use-envmodules
snakemake $1 --cluster-config baf_logr.json --cluster "sbatch -A {cluster.A} -t {cluster.t} -p {cluster.p} -n {cluster.n}" -s baf_logr.smk.py --latency-wait 10 --use-singularity --use-envmodules
