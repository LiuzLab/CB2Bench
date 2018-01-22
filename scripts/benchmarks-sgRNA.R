m_id <- list()

m_id$CC2_sgRNA.csv <- "p_value"
m_id$MAGeCK_sgRNA.csv <- "p.low"
m_id$DESeq2_sgRNA.csv <- "padj"
m_id$edgeR_sgRNA.csv <- "PValue"
all.df <- NULL


load("inst/extdata/nature-biotech.Rdata")

file.name <- "cache/nature-biotech/*/*_sgRNA.csv"
for(f in Sys.glob(file.name)) {
  m <- basename(f)
  dset <- basename(dirname(f))
  df <- read.csv(f)
  fdr <- df[,m_id[[m]]]
  if(m=="CC2_sgRNA.csv") {
    df <- df[,-1]
    fdr[df$t_value<0] <- 1
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

(pt1 <- heatmap.fdr(all.df, "sgRNA", c("CC2", "MAGeCK", "edgeR", "DESeq2")))
save_plot("figures/fig3-heatmap-sgRNA.tiff",pt1, base_height = 8)
(pt2 <- lineplot.fdr(all.df))
save_plot("figures/fig4-fdr-curve-sgRNA.tiff",pt2, base_height = 8, base_aspect_ratio = 1.6)
