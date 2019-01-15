library(Cairo)
source("R/heatmap_fdr.R")
source("scripts/benchmarks-gene.R")
heatmap.fdr(all.df, "gene", c("CC2", "ScreenBEAM", "PBNPA", "PinAPL-py", "sgRSEA", "HitSelect", "MAGeCK")) -> plot.fdr
save_plot("figures/heatmap_with_name.pdf", plot.fdr, device=cairo_pdf, base_width = 8, base_height = 12)
