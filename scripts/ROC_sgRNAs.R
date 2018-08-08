library(ggplot2)
library(plotROC)
library(cowplot)
library(tidyverse)
library(ggsci)
m_id <- list()

m_id$CC2_sgRNA.csv <- "p_value_neg"
m_id$MAGeCK_sgRNA.csv <- "p.low"
m_id[["PinAPL-py_sgRNA.csv"]] <- "p.value"
all.df <- NULL


load("inst/extdata/nature-biotech.Rdata")

file.name <- "cache/nature-biotech/*/*_sgRNA.csv"
for(f in Sys.glob(file.name)) {
  m <- basename(f)
  dset <- basename(dirname(f))
  df <- read.csv(f)
  score <- df[,m_id[[m]]]
  if(m=="CC2_sgRNA.csv") {
    df[,c(1,2)] <- df[,c(2,1)]
    score <- -score
  } else if(m=="MAGeCK_sgRNA.csv") {
    score <- 1-score
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
