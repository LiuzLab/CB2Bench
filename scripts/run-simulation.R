methods = list(
  MAGeCK = run.mageck,
  DESeq2 = run.DESeq2,
  edgeR = run.edgeR,
  sgRSEA = run.sgRSEA,
  PBNPA = run.PBNPA,
  #ScreenBEAM = run.ScreenBEAM,
  CC2 = run.CC2,
  CC2py = run.CC2py
)

df.tmp <- expand.grid(depth = c(100, 500),
                      noise = c(0.1, 0.5, 1.0),
                      effect = c(0.1, 0.2),
                      facs = c(0.1, 0.25))

results <- tibble()

selector <- read_delim("inst/extdata/twosided-expr.tsv", delim = "\t")

for (i in 1:nrow(df.tmp)) {
  facs <- df.tmp[i,]$facs
  depth <- df.tmp[i,]$depth
  noise <- df.tmp[i,]$noise
  effect <- df.tmp[i,]$effect
  df.dat <-
    load.sim(
      depth = depth,
      facs = facs,
      noise = noise,
      effect = effect
    )
  cache.dir <-
    file.path(sprintf("cache/sim/FACS_%.2f_DEPTH_%d_NOISE_%.2f_EFFECT_%.2f/", facs, depth, noise, effect))

  if(!dir.exists(cache.dir)) {
    dir.create(cache.dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  }
  results <- bind_rows(results,
                       tibble(
                         noise = noise,
                         effect = effect,
                         ret = list(run(df.dat, methods, selector, cache.dir))))
}
