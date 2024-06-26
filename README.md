# REH cell line multi-omics

## Citation
This is a public repository containing scripts described in the publication:

Lysenkova Wiklander et al.
"A multiomic characterization of the leukemia cell line REH using short- and long-read sequencing" https://doi.org/10.26508/lsa.202302481

## Data
Sequencing data and BAM files for this project are available at NCBI/SRA under the BioProject accession numbers PRJNA600820 and PRJNA834955. These data have been analyzed on a high performance cluster (HPC) using the commands in `01_hpc_processing`. The resulting analysis datasets are available at https://doi.org/10.5281/zenodo.108757278.

## Instructions
The scripts are numbered in the order they should be executed.

### HPC Bash Scripts
For HPC scripts, the full paths to source files have been omitted for simplicity. You will need to adjust them to your environment.

### Other scripts
To run the other scripts in this repository, you will need to do the following.

Install:
- R 4.2.1 and an integrated environment, e.g. RStudio
- R packages: chromoMap, RColorBrewer, VennDiagram
- Python 3.8
- Necessary Python packages: `pip install -r requirements.txt`
- SURVIVOR v1.0.7 built from: https://github.com/fritzsedlazeck/SURVIVOR and placed in the directory `02_sv_callset_analysis/bin`

Download the files from Zenodo and place them in the following directories:

data/coverage:
- copycat.ont.coverage.10kb.csv
- copycat.pb.coverage.10kb.csv
- copycat.pcrfree.coverage.10kb.csv

data/sv_callsets:
- illumina.tiddit.vcf
- ont.sniffles.vcf
- pb.sniffles.vcf

data/fusion_callsets/long_read:
- cupcake.long.csv
- cupcake.std.csv
- jaffa_results.csv

data/fusion_callsets/short_read:
- GM12878.fusionreport.txt
- illumina.all.txt
- illumina.filtered.csv
- REH.arriba.fusions.tsv
- REH.fusioncatcher.fusion-genes.txt
- REH.pizzly.txt
- REH.squid.fusions.annotated.txt
- REH.starfusion.abridged.tsv
