library(ggsci)
m_id <- list()

m_id$MAGeCK_gene.csv <- "neg.score"
m_id$PBNPA_gene.csv <- "neg.pvalue"
m_id$ScreenBEAM_gene.csv <- "bsta"
m_id$sgRSEA_gene.csv <- "p.value.neg"
m_id$CC2_gene.csv <- "p_value_neg"
m_id$HitSelect_gene.csv <- "effect_size"
m_id[["PinAPL-py_gene.csv"]] <- "p.value.combined"
all.df <- NULL

file.name <- "cache/nature-biotech/*/*_gene.csv"
for(f in Sys.glob(file.name)) {
  m <- basename(f)
  if(m=="DESeq2_gene.csv") next()
  dset <- basename(dirname(f))
  df <- read.csv(f)
  score <- df[,m_id[[m]]]
  if(m=="MAGeCK_gene.csv") {
    score <- 1-score
  } else if(m=="ScreenBEAM_gene.csv") {
    score <- score
  } else if (m=="sgRSEA_gene.csv") {
    score <- 1-score
  } else if(m=="CC2_gene.csv") {
    score <- 1-score
  } else if(m=="PBNPA_gene.csv") {
    score <- 1-score
  } else if(m=="PinAPL-py_gene.csv") {
    score <- 1-score
  }

  m <- strsplit(m,"\\_")[[1]][1]
  new.df <- data.frame(dataset=dset, method=m, gene=df[,1], score=score)
  colnames(new.df) <- c("dataset", "method", "gene", "score")
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

all.df <- all.df %>%
  filter(dataset!="shRNA.RT112", dataset!="shRNA.UMUC3")
all.df %>% write_csv("inst/extdata/gene_score.csv")


plots <- list()
all.df <- all.df %>%
  filter(dataset!="shRNA.RT112", dataset!="shRNA.UMUC3")
for(dset in unique(all.df$dataset)) {
  p <- ggplot()
  for(mat in rev(unique(all.df$method))) {
    tmp <- all.df %>% filter(dataset==dset, method==mat)
    if(nrow(tmp)==0) {
      next
    }
    AUC <- calc_auc(ggplot()+geom_roc(data=tmp, aes(d = essential, m = score, color=method)))$AUC
    tmp$method <- paste(mat, sprintf("(AUC %.2f)", AUC))
    p <- p + geom_roc(data=tmp, aes(d = essential, m = score, color=method), labels = FALSE, size=0.5)
  }
  p <- p + ggtitle(dset) + scale_color_npg() +
    xlab("False Positive Rate") +
    ylab("True Positive Rate") +
    theme(legend.position = c(0.3,0.25)) +
    coord_equal(ratio = 1)
  plots[[dset]] <- p
}
plot_grid(plotlist=plots, nrow=1)
