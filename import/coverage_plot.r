#install.packages("ggplot2")
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
packageVersion('ggplot2')
path1 = getwd()
print(path1)
sample = args[1]
print(sample)
bed_file <- list.files(path = path1, pattern = "\\.cov.bed$")
cov_plot <- function(BED, plot_title){
    cov_df = read.csv(BED, sep = "\t", header=F)
    names(cov_df) <- c("Ref", "Locus", "Depth")    
    library(ggplot2)
    options(repr.plot.width=20, repr.plot.height=10)
    p <- ggplot(cov_df, aes(Locus, Depth)) + 
        geom_line(color= 'red') + scale_x_continuous(n.breaks = 20, guide = guide_axis(check.overlap = TRUE)) 
    plot = p + facet_grid(cols = vars(Ref), scales = "free_x") + 
        theme_bw() +
        theme(text = element_text(size = 20), plot.title = element_text(size = 30, hjust =0.5, face = "bold")) +
        ggtitle(sample) +  theme(strip.text.x = element_text(angle = 90, size=10, face = "bold")) +theme(axis.text.x = element_text(angle =45, hjust = 1, size = 10,face = "bold" )) +theme(plot.title = element_text(size=15))
    plot_png = paste(plot_title,"_coverage.png", sep="")
    png(plot_png, width=1500, height=600)
    print(plot)
    dev.off()
}
for (file in bed_file){
  ref = strsplit(file, ".cov.bed")
  ref = gsub(sample, "", ref)
  print(ref)
  cov_plot(file, ref)
}