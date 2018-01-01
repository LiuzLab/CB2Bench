m_id <- list()
m_id$CC2py_gene.csv <- "p_value"
m_id$CC2_gene.csv <- "gpvalue"
m_id$MAGeCK_gene.csv <- "twosided.p.value"
m_id$PBNPA_gene.csv <- "p.value.twosides"
m_id$sgRSEA_gene.csv <- "p.value.twosides"
m_id$DESeq2_gene.csv <- "pvalue"

all.df <- NULL

file.name <- "cache/sim/FACS_*/*_gene.csv"

for(f in Sys.glob(file.name)) {
  m <- basename(f)
  print(m)
  dset <- basename(dirname(f))
  df <- read.csv(f)
  fdr <- df[,m_id[[m]]]
  m <- strsplit(m,"\\_")[[1]][1]
  if(m != "CC2" && m != "CC2py") {
    fdr <- p.adjust(fdr, method="fdr")
  }
  new.df <- data.frame(dataset=dset, method=m, gene=df[,1], fdr=fdr)
  colnames(new.df) <- c("dataset", "method", "gene", "fdr")

  param <- strsplit(gsub("_"," ",gsub("=", " ", dset))," ")[[1]]
  depth <- as.integer(param[4])
  noise <- as.double(param[6])
  effect <- as.double(param[8])
  facs <- as.double(param[2])

  dataset <- load.sim(depth, facs, noise, effect)

  ess <- dataset %>% mutate(essential=ifelse(class=="decreasing"|class=="increasing",1,0)) %>%
    group_by(gene) %>%
    summarise(essential=mean(essential)) %>%
    select(gene, essential)
  new.df <- left_join(new.df, ess, by=c("gene"="gene"))
  new.df$noise <- noise
  new.df$effect <- effect
  new.df$facs <- facs
  new.df$depth <- depth
  #print(head(new.df))

  if(is.null(all.df)) {
    all.df <- new.df
  }
  else {
    all.df <- rbind(all.df, new.df)
  }
}
all.df$fdr[is.na(all.df$fdr)] <- 1



ct <- c(0.2, 0.1, 0.05, 0.01, 0.005, 0.001)
df.prof <- tibble()
for (dset in unique(all.df$dataset)) {
  for (mat in unique(all.df$method)) {
    tmp <- all.df %>% filter(dataset == dset, method == mat)
    print(head(tmp))
    facs <- tmp$facs[1]
    depth <- tmp$depth[1]
    noise <- tmp$noise[1]
    effect <- tmp$effect[1]

    for (fdr in ct) {
      #fdr <- 10^fdr
      TP <- sum((tmp$fdr <= fdr) & (tmp$essential == 1))
      FP <- sum((tmp$fdr <= fdr) & (tmp$essential == 0))
      FN <- sum((tmp$fdr > fdr) & (tmp$essential == 1))

      precision <- TP / max(1, (TP + FP))
      recall <- TP / max(1, (TP + FN))
      fmeasure <- 2 * (precision * recall) / (precision + recall)
      if (is.na(fmeasure))
        fmeasure <- 0
      df.prof <-
        bind_rows(
          df.prof,
          tibble(
            dataset = dset,
            method = mat,
            FDR = fdr,
            value = precision,
            measure = "precision",
            facs = facs,
            depth = depth,
            noise = noise,
            effect = effect,
            TP = TP,
            FP = FP,
            FN = FN
          )
        )
      df.prof <-
        bind_rows(
          df.prof,
          tibble(
            dataset = dset,
            method = mat,
            FDR = fdr,
            value = recall,
            measure = "recall",
            facs = facs,
            depth = depth,
            noise = noise,
            effect = effect,
            TP = TP,
            FP = FP,
            FN = FN
          )
        )
      df.prof <-
        bind_rows(
          df.prof,
          tibble(
            dataset = dset,
            method = mat,
            FDR = fdr,
            value = fmeasure,
            measure = "F-measure",
            facs = facs,
            depth = depth,
            noise = noise,
            effect = effect,
            TP = TP,
            FP = FP,
            FN = FN
          )
        )
    }
  }
}


(pt <- df.prof %>% mutate(FDR = factor(FDR, levels = ct)) %>%
    filter(measure=="F-measure", facs==0.25, depth==100, noise < 1.0) %>%
  ggplot(aes(x = FDR, y = value)) +
  geom_point(aes(colour = method, shape =
                   method), alpha = 0.5) +
  geom_line(aes(colour = method, group = method), alpha = 0.5) +
  facet_grid(sprintf("ES = %.1f",effect)~sprintf("Noise = %f", noise)) +
  xlab("FDR") + ylab("F-measure") + scale_colour_d3() +
    theme(axis.text.x = element_text(angle = 90)))
#save_plot("figures/fig5-fdr-curve-sim-gene.tiff",pt, base_height = 8, base_aspect_ratio = 1.6)
