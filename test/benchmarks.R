library(readr)

m_id <- list()

m_id$CC2_gene.csv <- "gpvalue"
m_id$MAGeCK_gene.csv <- "neg.fdr"
m_id$PBNPA_gene.csv <- "neg.fdr"
m_id$ScreenBEAM_gene.csv <- "FDR"
m_id$sgRSEA_gene.csv <- "FDR.neg"

all.df <- NULL


load("inst/extdata/nature-biotech.Rdata")

file.name <- "cache/nature-biotech/*/*_gene.csv"
for(f in Sys.glob(file.name)) {
  m <- basename(f)
  dset <- basename(dirname(f))
  df <- read.csv(f)
  fdr <- df[,m_id[[m]]]
  if(m=="CC2_gene.csv") fdr[df$gtvalue<0] <- 1
  if(m=="ScreenBEAM_gene.csv") fdr[df$beta>0] <- 1
  m <- strsplit(m,"\\_")[[1]][1]
  new.df <- data.frame(dataset=dset, method=m, gene=df[,1], fdr=fdr)
  colnames(new.df) <- c("dataset", "method", "gene", "fdr")
  ess <- dataset[[dset]] %>% mutate(essential=ifelse(class=="decreasing",1,0)) %>%
    group_by(gene) %>%
    summarise(essential=mean(essential)) %>%
    select(gene, essential)
  new.df <- left_join(new.df, ess, by=c("gene"="gene"))

  if(is.null(all.df)) all.df <- new.df
  else {
    all.df <- rbind(all.df, new.df)
  }
}


df.gene.plot <- tibble()
for(dset in unique(all.df$dataset)) {
  for(mat in unique(all.df$method)) {
    tmp <- all.df %>% filter(dataset==dset, method==mat)
    for(fdr in seq(-10,-1,1)) {
      fdr <- 10^fdr
      TP <- sum((tmp$fdr <= fdr) & (tmp$essential == 1))
      FP <- sum((tmp$fdr <= fdr) & (tmp$essential == 0))
      FN <- sum((tmp$fdr > fdr) & (tmp$essential == 1))

      precision <- TP / max(1,(TP+FP))
      recall <- TP / max(1,(TP+FN))
      fmeasure <- 2*(precision*recall)/(precision+recall)
      if(is.na(fmeasure)) fmeasure <- 0
      #df.gene.plot <- bind_rows(df.gene.plot, tibble(dataset=d, method=mat, FDR=log10(fdr), PR=precision, RC=recall, TP=TP, FP=FP, FN=FN))
      df.gene.plot <- bind_rows(df.gene.plot, tibble(dataset=dset, method=mat, FDR=log10(fdr), value=precision, measure="precision", TP=TP, FP=FP, FN=FN))
      df.gene.plot <- bind_rows(df.gene.plot, tibble(dataset=dset, method=mat, FDR=log10(fdr), value=recall, measure="recall", TP=TP, FP=FP, FN=FN))
      df.gene.plot <- bind_rows(df.gene.plot, tibble(dataset=dset,method=mat, FDR=log10(fdr), value=fmeasure, measure="F-measure", TP=TP, FP=FP, FN=FN))
    }
  }
}

pt <- df.gene.plot %>% filter(measure=="F-measure") %>%
  ggplot( aes(x=-FDR, y=value)) + geom_point(aes(colour=method, shape=method), alpha=0.5) +
  geom_line(aes(colour=method), alpha=0.5) + facet_grid(.~dataset) + ylim(0,1) + scale_x_reverse(breaks = c(1,5,10))+ xlab("FDR") + ylab("F-measure")
save_plot(pt, filename = "nat-biotech.pdf", base_height = 4, base_width = 8)


for(dset in unique(all.df$dataset)) {
  ess <- dataset[[dset]] %>% mutate(essential=ifelse(class=="decreasing",1,0)) %>%
    group_by(gene) %>%
    summarise(essential=mean(essential)) %>%
    select(gene, essential)

  x <- all.df %>% filter(dataset==dset) %>% select(-1,-5) %>% spread(method, fdr) %>%
    remove_rownames() %>% left_join(ess, by="gene")

  write.table(x, file = sprintf("%s.tsv", dset), sep="\t", row.names = F)
}

