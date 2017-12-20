methods = list(
  MAGeCK = run.mageck,
  DESeq2 = run.DESeq2,
  edgeR = run.edgeR,
  sgRSEA = run.sgRSEA,
  PBNPA = run.PBNPA,
  #ScreenBEAM = run.ScreenBEAM,
  CC2 = run.CC2
)

df.tmp <- expand.grid(noise = seq(1, 6, 1),
                      effect = c(0.1, 0.2))

results <- tibble()

selector <- read_delim("inst/extdata/twosided-expr.tsv", delim = "\t")

for (i in 1:nrow(df.tmp)) {
  noise <- df.tmp[i,]$noise
  effect <- df.tmp[i,]$effect
  df.dat <-
    load.sim(
      depth = 300,
      facs = 0.25,
      noise = noise,
      effect = effect
    )
  cache.dir <-
    file.path(sprintf("/Users/hwan/Sandbox/CC2Sim/cache/sim/FACS_0.25_DEPTH_100_NOISE_%.2f_EFFECT_%.2f/", noise, effect))

  if(!dir.exists(cache.dir)) {
    dir.create(cache.dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
  }
  results <- bind_rows(results,
                       tibble(
                         noise = noise,
                         effect = effect,
                         ret = list(run(df.dat, methods, selector, cache.dir))))
}
