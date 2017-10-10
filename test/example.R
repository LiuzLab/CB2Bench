#library(CC2Sim)
library(dplyr)
library(tidyr)
#source("R/utils.R")
#source("R/expr.R")
load(system.file("extdata", "nature-biotech.Rdata", package = "CC2Sim"))
RT112 <- dataset$RT112

methods <- list(
  mageck = run.mageck,
  screenBEAM = run.ScreenBEAM,
  DESeq2 = run.DESeq2,
  edgeR = run.edgeR,
  sgRSEA = run.sgRSEA,
  PBNPA = run.PBNPA,
  CC2 = run.mbttest,
  RSA = run.RSA
)

RT112.ret <- run(RT112, methods)

UMUC3 <- dataset$UMUC3
UMUC3.ret <- run(UMUC3, methods)

RT112.1 <- plot.all(RT112.ret, "ROC", "RTT12 dataset benchmark (AUCROC)")
RT112.2 <- plot.all(RT112.ret, "PR", "RTT12 dataset benchmark (AUCPRC)")
final <- plot_grid(RT112.1, RT112.2)
final


UMUC3.1 <- plot.all(UMUC3.ret, "ROC", "UMUC3 dataset benchmark (AUCROC)")
UMUC3.2 <- plot.all(UMUC3.ret, "PR", "UMUC3 dataset benchmark (AUCPRC)")
final <- plot_grid(UMUC3.1, UMUC3.2)
final

final <- plot_grid(RT112.1, RT112.2, UMUC3.1, UMUC3.2, ncol=2)
save_plot(filename = "CC2-NatBioTech.PDF", final, base_aspect_ratio = 1.8, base_height = 8)

rank.heatmap <- function(gene.ret, filename) {
  gene.rank <- gene.ret[,-c(1,ncol(gene.ret))]
  for(i in 1:6) {
    gene.rank[,i] <- rank(-gene.ret[,i+1])
  }
  row.names(gene.rank) <- gene.ret$gene
  library(pheatmap)
  library(RColorBrewer)

  df.anno <- data.frame(label=ifelse(gene.ret$label,"essential", "non essential"))
  row.names(df.anno) <- gene.ret$gene
  order(rowSums(gene.rank))
  pheatmap(t(gene.rank[order(rowSums(gene.rank)),]), cluster_cols = F,
           color = colorRampPalette((brewer.pal(n = 11, name ="RdYlBu")))(100),
           annotation_col = df.anno, filename = filename)
}
rank.heatmap(RT112.ret$gene, "RT112-gene-heatmap.pdf")
rank.heatmap(UMUC3.ret$gene, "UMUC3-gene-heatmap.pdf")
