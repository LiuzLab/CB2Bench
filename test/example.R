library(CC2Sim)

load(system.file("extdata", "nature-biotech.Rdata", package = "CC2Sim"))
RT112 <- dataset$RT112

methods <- list(
  mageck = run.mageck,
  DESeq2 = run.DESeq2,
  edgeR = run.edgeR,
  sgRSEA = run.sgRSEA,
  PBNA = run.PBNPA,
  #screenBEAM = run.ScreenBEAM,
  CC2 = run.mbttest,
  RSA = run.RSA
)

RT112.ret <- run(RT112, methods)


RT112.1 <- plot.all(RT112.ret, "ROC", "RTT12 dataset benchmark (AUCROC)")
RT112.2 <- plot.all(RT112.ret, "PR", "RTT12 dataset benchmark (AUCPRC)")
final <- plot_grid(RT112.1, RT112.2)
final
