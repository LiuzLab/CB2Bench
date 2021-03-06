m_id <- list()

m_id$MAGeCK_gene.csv <- "neg.fdr"
m_id$PBNPA_gene.csv <- "neg.fdr"
m_id$ScreenBEAM_gene.csv <- "FDR"
#m_id$sgRSEA_gene.csv <- "FDR.neg"
m_id$CC2_gene.csv <- "fdr_neg"
m_id$HitSelect_gene.csv <- "fdr"
m_id[["PinAPL-py_gene.csv"]] <- "p.value.combined"
m_id$CRISPhieRmix_gene.csv <- "FDR"
all.df <- NULL

file.name <- "cache/nature-biotech/*/*_gene.csv"
for(f in Sys.glob(file.name)) {
  m <- basename(f)
  if(m=="DESeq2_gene.csv") next()
  if(m=="sgRSEA_gene.csv") next()
  dset <- basename(dirname(f))
  df <- read.csv(f)
  fdr <- df[,m_id[[m]]]
  if(m=="ScreenBEAM_gene.csv") {
     fdr[df$beta>0] <- 1
  }

  m <- strsplit(m,"\\_")[[1]][1]
  new.df <- data.frame(dataset=dset, method=m, gene=df[,1], fdr=fdr)
  colnames(new.df) <- c("dataset", "method", "gene", "fdr")
  ess <- dataset[[dset]] %>%
    group_by(gene) %>%
    summarise(essential=mean(class)) %>%
    select(gene, essential)
  new.df <- left_join(new.df, ess, by=c("gene"="gene"))

  if(is.null(all.df)) all.df <- new.df
  else {
    all.df <- rbind(all.df, new.df)
  }
}

all.df %>% write_csv("inst/extdata/gene_fdr.csv")

(pt1 <- heatmap.fdr(all.df, "gene", c("CC2", "CRISPhieRmix", "ScreenBEAM", "PBNPA", "PinAPL-py", "HitSelect", "MAGeCK" )))
save_plot("figures/fig1-heatmap-gene.tiff",pt1, base_height = 8)
(pt2 <- lineplot.fdr(all.df))
save_plot("figures/fig-fdr-curve-gene.pdf",pt2, base_height = 8, base_aspect_ratio = 1.6)

#(pt3 <- lineplot.fdr.f1(all.df %>% filter(method != "DESeq2")))
#save_plot("figures/fig2-fdr-curve-gene.tiff",pt2, base_height = 8, base_aspect_ratio = 1.6)

