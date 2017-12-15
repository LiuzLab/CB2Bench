library(readr)

m_id <- list()

m_id$CC2_sgRNA.csv <- "pvalue"
m_id$MAGeCK_sgRNA.csv <- "p.low"
m_id$CC1_sgRNA.csv <- "pvalue"
m_id$DESeq2_sgRNA.csv <- "padj"
m_id$edgeR_sgRNA.csv <- "PValue"
all.df <- NULL


load("inst/extdata/nature-biotech.Rdata")

file.name <- "cache/nature-biotech/*/*_sgRNA.csv"
for(f in Sys.glob(file.name)) {
  m <- basename(f)
  dset <- basename(dirname(f))
  df <- read.csv(f)
  print(m)
  print(head(df))
  fdr <- df[,m_id[[m]]]
  if(m=="CC2_sgRNA.csv") {
    df <- df[,c(-1,-2)]
    fdr[df$tvalue<0] <- 1
  } else if(m=="MAGeCK_sgRNA.csv") {
    fdr <- p.adjust(fdr, method="fdr")
  } else if(m=="CC1_sgRNA.csv"||m=="DESeq2_sgRNA.csv"){
    fdr <- p.adjust(fdr, method="fdr")
    fdr[df$fc>0] <- 1
  } else if(m=="edgeR_sgRNA.csv") {
    fdr <- p.adjust(fdr, method="fdr")
    fdr[df$logFC>0] <- 1
  }

  m <- strsplit(m,"\\_")[[1]][1]
  new.df <- data.frame(dataset=dset, method=m, sgRNA=df[,1], fdr=fdr)
  colnames(new.df) <- c("dataset", "method", "sgRNA", "fdr")
  ess <- dataset[[dset]] %>% mutate(essential=ifelse(class=="decreasing",1,0)) %>%
    select(sgRNA, essential)
  new.df <- left_join(new.df, ess, by=c("sgRNA"="sgRNA"))

  if(is.null(all.df)) all.df <- new.df
  else {
    all.df <- rbind(all.df, new.df)
  }
}

all.df$fdr[is.na(all.df$fdr)] <- 1

df.sgRNA.plot <- tibble()
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
      #df.sgRNA.plot <- bind_rows(df.sgRNA.plot, tibble(dataset=d, method=mat, FDR=log10(fdr), PR=precision, RC=recall, TP=TP, FP=FP, FN=FN))
      df.sgRNA.plot <- bind_rows(df.sgRNA.plot, tibble(dataset=dset, method=mat, FDR=log10(fdr), value=precision, measure="precision", TP=TP, FP=FP, FN=FN))
      df.sgRNA.plot <- bind_rows(df.sgRNA.plot, tibble(dataset=dset, method=mat, FDR=log10(fdr), value=recall, measure="recall", TP=TP, FP=FP, FN=FN))
      df.sgRNA.plot <- bind_rows(df.sgRNA.plot, tibble(dataset=dset,method=mat, FDR=log10(fdr), value=fmeasure, measure="F-measure", TP=TP, FP=FP, FN=FN))
    }
  }
}

pt <- df.sgRNA.plot %>% filter(measure==measure) %>% filter(FDR>-6) %>%
  ggplot( aes(x=-FDR, y=value)) + geom_point(aes(colour=method, shape=method), alpha=0.5) +
  geom_line(aes(colour=method), alpha=0.5) + facet_grid(measure~dataset) + ylim(0,1) + scale_x_reverse(breaks = seq(1,5,1))+ xlab("FDR") + ylab("measure")
save_plot(pt, filename = "nat-biotech-sgRNA-FDR.pdf", base_height = 4, base_width = 8)

library(pheatmap)
heatmap <- list()
for(dset in unique(all.df$dataset)) {
  ess <- dataset[[dset]] %>% mutate(essential=ifelse(class=="decreasing",1,0)) %>%
    group_by(sgRNA) %>%
    summarise(essential=mean(essential)) %>%
    select(sgRNA, essential)
  x <- all.df %>% filter(dataset==dset) %>% select(-1,-5) %>% spread(method, fdr) %>%
    remove_rownames()

  x[x<1e-10] <- 1e-10
  x[,-1] <- floor(-log10(x[,-1]))
  x <- x[order(-rowSums(x[,-1])),]
  x <- x %>% left_join(ess, by="sgRNA")
  x <- x[order(-x$essential),]
  #x[x$essential==0,] <- x[x$essential==0,][order(rowSums(x[x$essential==0,c(-1,-ncol(x))])),]

  x <- x %>% remove_rownames()
  x <- x[,c(1,3,2,4:ncol(x))]
  heatmap[[dset]] <- pheatmap(t(column_to_rownames(x, "sgRNA")), scale = "none",
                              cluster_cols = F, cluster_rows = F, main = dset,
                              show_colnames=F)$gtable
}

save_plot(plot_grid(plotlist = heatmap,  ncol = 1), filename = "nat-biotech-sgRNA-heatmap.svg", base_height=10, base_width=12)
