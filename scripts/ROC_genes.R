m_id <- list()

m_id$MAGeCK_gene.csv <- "neg.score"
m_id$PBNPA_gene.csv <- "neg.pvalue"
m_id$ScreenBEAM_gene.csv <- "bsta"
m_id$sgRSEA_gene.csv <- "NScore.y"
m_id$CC2_gene.csv <- "p_value_neg"
m_id$HitSelect_gene.csv <- "effect_size"
m_id[["PinAPL-py_gene.csv"]] <- "Fisher.Statistic"
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
    score <- -score 
  } else if(m=="CC2_gene.csv") {
    score <- 1-score
  } else if(m=="PBNPA_gene.csv") {
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


p <- ggplot() 
for(mat in rev(unique(all.df$method))) {
  tmp <- all.df %>% filter(dataset=="CRISPR.RT112", method==mat) 
  AUC <- calc_auc(ggplot()+geom_roc(data=tmp, aes(d = essential, m = score, color=method)))$AUC
  cat(mat, AUC, "\n")                  
  p <- p + geom_roc(data=tmp, aes(d = essential, m = score, color=method), labels = FALSE, 
                    size=0.5) 
}
(p1 <- p + ggtitle("CRISPR.RT112"))

p <- ggplot() 
for(mat in rev(unique(all.df$method))) {
  tmp <- all.df %>% filter(dataset=="CRISPR.UMUC3", method==mat) 
  AUC <- calc_auc(ggplot()+geom_roc(data=tmp, aes(d = essential, m = score, color=method)))$AUC
  cat(mat, AUC, "\n")                  
  p <- p + geom_roc(data=tmp, aes(d = essential, m = score, color=method), labels = FALSE, size=0.5) 
}
p2 <- p + ggtitle("CRISPR.UMUC3")

p <- ggplot() 
for(mat in rev(unique(all.df$method))) {
  tmp <- all.df %>% filter(dataset=="CRISPRi.RT112", method==mat) 
  AUC <- calc_auc(ggplot()+geom_roc(data=tmp, aes(d = essential, m = score, color=method)))$AUC
  cat(mat, AUC, "\n")                  
  p <- p + geom_roc(data=tmp, aes(d = essential, m = score, color=method), labels = FALSE, size=0.5) 
}
p3 <- p + ggtitle("CRISPRi.RT112")

plot_grid(p1, p2, p3)

