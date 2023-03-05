library(VennDiagram)
library(RColorBrewer)

colors <- brewer.pal(3, "Pastel2")

a_val=0.5
alpha=rep(a_val, 3)

fontfamily="sans"
fontface="bold"

text_size_numbers = .5
text_size_labels = 0

draw_venn <- function(overlap_file, output_file) {
  input = paste("out/", overlap_file, sep = "")
  output = paste("out/", output_file, sep = "")
  t=read.table(input,header=F)
  venn_data=list(Illumina=which(t[,1]==1), ONT=which(t[,2]==1), PB=which(t[,3]==1))
  venn.diagram(venn_data,  
               fill = colors, 
               alpha = alpha, 
               filename = output, 
               output=TRUE,
               
               # Output features
               imagetype="png" ,
               height = 480 , 
               width = 480 , 
               resolution = 300,
               compression = "lzw",
               
               # Circles
               lwd = 2,
               lty = 'blank',
               
               # Numbers
               cex = text_size_numbers,
               fontface = fontface,
               fontfamily = fontfamily,
               
               # Set names
               cat.cex = text_size_labels,
               cat.fontface = fontface,
               cat.default.pos = "outer",
               cat.pos = c(-27, 27, 135),
               cat.dist = c(0.055, 0.055, 0.085),
               cat.fontfamily = fontfamily,
               rotation = 1);
}

draw_venn("REH_merged_all_overlap.txt", "venn_all.png")
draw_venn("REH_merged_small_overlap.txt", "venn_small.png")
draw_venn("REH_merged_medium_overlap.txt", "venn_medium.png")
draw_venn("REH_merged_large_overlap.txt", "venn_large.png")
