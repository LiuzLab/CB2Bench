m_id <- list()

m_id$CC2_sgRNA.csv <- "fdr_neg"
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
    df[,c(1,2)] <- df[,c(2,1)]
  } else if(m=="MAGeCK_sgRNA.csv") {
    fdr <- p.adjust(fdr, method="fdr")
  } else if(m=="DESeq2_sgRNA.csv"){
    fdr[df$stat>0] <- 1
  } else if(m=="edgeR_sgRNA.csv") {
    fdr <- p.adjust(fdr, method="fdr")
    fdr[df$logFC>0] <- 1
  }

  m <- strsplit(m,"\\_")[[1]][1]
  new.df <- data.frame(dataset=dset, method=m, sgRNA=df[,1], fdr=fdr)
  colnames(new.df) <- c("dataset", "method", "sgRNA", "fdr")
  new.df <- new.df %>% left_join(dataset[[dset]] %>% select(sgRNA, gene), by="sgRNA")
  ess <- dataset[[dset]] %>% mutate(essential=ifelse(class=="decreasing",1,0)) %>%
    select(sgRNA, essential)
  new.df <- left_join(new.df, ess, by=c("sgRNA"="sgRNA"))

  if(is.null(all.df)) all.df <- new.df
  else {
    all.df <- rbind(all.df, new.df)
  }
}
all.df %>% select(dataset, method, sgRNA, gene, fdr, essential) %>%
  write_csv( "../CC2Sim-shiny-report/data/sgRNA.csv")

