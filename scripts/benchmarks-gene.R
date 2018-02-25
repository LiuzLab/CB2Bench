m_id <- list()

m_id$MAGeCK_gene.csv <- "neg.fdr"
m_id$PBNPA_gene.csv <- "neg.fdr"
m_id$ScreenBEAM_gene.csv <- "FDR"
m_id$sgRSEA_gene.csv <- "FDR.neg"
m_id$CC2_gene.csv <- "fdr_neg"
all.df <- NULL

file.name <- "cache/nature-biotech/*/*_gene.csv"
for(f in Sys.glob(file.name)) {
  m <- basename(f)
  if(m=="DESeq2_gene.csv") next()
  dset <- basename(dirname(f))
  df <- read.csv(f)
  fdr <- df[,m_id[[m]]]
  if(m=="ScreenBEAM_gene.csv") {
    fdr[df$beta>0] <- 1
  }

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

all.df %>% write_csv("../CC2Sim-shiny-report/data/gene.csv")
all.df <- all.df %>%
  filter(dataset!="shRNA.RT112", dataset!="shRNA.UMUC3")

(pt1 <- heatmap.fdr(all.df, "gene", c("CC2", "ScreenBEAM", "PBNPA", "MAGeCK")))
save_plot("figures/fig1-heatmap-gene.tiff",pt1, base_height = 8)
(pt2 <- lineplot.fdr(all.df %>% filter(method != "DESeq2")))
save_plot("figures/fig2-fdr-curve-gene.tiff",pt2, base_height = 8, base_aspect_ratio = 1.6)

(pt3 <- lineplot.fdr.f1(all.df %>% filter(method != "DESeq2")))
save_plot("figures/fig2-fdr-curve-gene.tiff",pt2, base_height = 8, base_aspect_ratio = 1.6)


save_plot("figures/fig.tiff",plot_grid(pt1, pt3, ncol=1, labels = "auto"), base_height = 6)
