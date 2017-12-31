selector <- read_delim("inst/extdata/negative-expr.tsv", delim="\t")
load("inst/extdata/nature-biotech.Rdata")

methods = list(
  MAGeCK = run.mageck,
  DESeq2 = run.DESeq2,
  edgeR = run.edgeR,
  sgRSEA = run.sgRSEA,
  PBNPA = run.PBNPA,
  ScreenBEAM = run.ScreenBEAM,
  CC2 = run.CC2,
  CC2py = run.CC2py
)

ret <- list()
for(d in names(dataset)) {
  df.dat <- dataset[[d]]
  cache.dir <- file.path("/Users/hwan/Sandbox/CC2Sim/cache/nature-biotech", d )
  dir.create(cache.dir, showWarnings = T, recursive = T, mode = "0777")
  ret[[d]] <- run(df.dat, methods, selector, cache.dir)
}
