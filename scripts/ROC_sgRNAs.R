library(ggplot2)
library(plotROC)
m_id <- list()

m_id$CC2_sgRNA.csv <- "t_value"
m_id$MAGeCK_sgRNA.csv <- "score"
m_id[["PinAPL-py_sgRNA.csv"]] <- "p.value"
all.df <- NULL


load("inst/extdata/nature-biotech.Rdata")

file.name <- "cache/nature-biotech/*/*_sgRNA.csv"
for(f in Sys.glob(file.name)) {
  m <- basename(f)
  dset <- basename(dirname(f))
  print(m)
  df <- read.csv(f)
  print(head(df))
  score <- df[,m_id[[m]]]
  if(m=="CC2_sgRNA.csv") {
    df[,c(1,2)] <- df[,c(2,1)]
    score <- -score
  } else if(m=="MAGeCK_sgRNA.csv") {
    0
  } else if(m=="PinAPL-py_sgRNA.csv") {
    score <- 1-score
  } else {
    next
  }
  m <- strsplit(m,"\\_")[[1]][1]
  new.df <- data.frame(dataset=dset, method=m, sgRNA=df[,1], score=score)
  colnames(new.df) <- c("dataset", "method", "sgRNA", "score")
  ess <- dataset[[dset]] %>% mutate(essential=ifelse(class=="decreasing",1,0)) %>%
    select(sgRNA, essential)
  new.df <- left_join(new.df, ess, by=c("sgRNA"="sgRNA"))
  
  if(is.null(all.df)) all.df <- new.df
  else {
    all.df <- rbind(all.df, new.df)
  }
}

p <- ggplot() 
for(mat in rev(unique(all.df$method))) {
  tmp <- all.df %>% filter(dataset=="CRISPR.RT112", method==mat) 
  AUC <- calc_auc(ggplot()+geom_roc(data=tmp, aes(d = essential, m = score, color=method)))$AUC
  cat(mat, AUC, "\n")                  
  p <- p + geom_roc(data=tmp, aes(d = essential, m = score, color=method), labels = FALSE, alpha=0.5) 
}
p1 <- p + ggtitle("CRISPR.RT112")

p <- ggplot() 
for(mat in rev(unique(all.df$method))) {
  tmp <- all.df %>% filter(dataset=="CRISPR.UMUC3", method==mat) 
  AUC <- calc_auc(ggplot()+geom_roc(data=tmp, aes(d = essential, m = score, color=method)))$AUC
  cat(mat, AUC, "\n")                  
  p <- p + geom_roc(data=tmp, aes(d = essential, m = score, color=method), labels = FALSE, alpha=0.5) 
}
p2 <- p + ggtitle("CRISPR.UMUC3")

p <- ggplot() 
for(mat in rev(unique(all.df$method))) {
  tmp <- all.df %>% filter(dataset=="CRISPRi.RT112", method==mat) 
  AUC <- calc_auc(ggplot()+geom_roc(data=tmp, aes(d = essential, m = score, color=method)))$AUC
  cat(mat, AUC, "\n")                  
  
  p <- p + geom_roc(data=tmp, aes(d = essential, m = score, color=method), labels = FALSE, alpha=0.5) 
}
p3 <- p + ggtitle("CRISPRi.RT112")

plot_grid(p1, p2, p3)
