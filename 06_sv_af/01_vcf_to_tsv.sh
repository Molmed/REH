sed -n '/#CHROM/,$p' ../02_sv_callset_analysis/01_consensus_sets/out/REH_consensus_3x.vcf > out/trimmed.tmp.tsv
tail -n +2 out/trimmed.tmp.tsv > out/trimmed.tsv
rm out/trimmed.tmp.tsv
