# REH cell line digital karyotyping

## Citation
This is a public repository containing scripts described in the publication "A complete digital karyotype of the B-cell leukemia REH cell line resolved by long-read sequencing" (Manuscript)

## Instructions
The scripts are numbered in the order they should be executed. For HPC scripts, the full paths to source files have been omitted for simplicity.

## Pre-requisites
- R 4.2.1 and an integrated environment, e.g. RStudio
- R packages: chromoMap, RColorBrewer, VennDiagram
- Python 3.8 
- Install necessary Python packages: `pip install -r requirements.txt`
- SURVIVOR v1.0.7 built from: https://github.com/fritzsedlazeck/SURVIVOR and placed in the directory `02_sv_callset_analysis/bin`

## Data
Primary data for this project is available at NCBI/SRA under the BioProject accession numbers PRJNA600820 and PRJNA834955.

Analysis datasets are available at [ZENODO DOI].

To run the scripts in this directory, download the analysis files from Zenodo and place them in the following directories:

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