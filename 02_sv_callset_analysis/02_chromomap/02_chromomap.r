# Usage: run from an integrated environment such as RStudio

library(chromoMap)

chr_file = paste("data", "chrom_lengths.tsv", sep = "/")
variants_file = paste("out", "bins.tsv", sep = "/")

chromoMap(
  chr_file,
  variants_file,
  data_based_color_map = T,
  data_type = "numeric",
  legend = T,
  data_colors = list(c("#FFFACD","#008080")),
  aggregate_func = "sum",
  export.options = T
)