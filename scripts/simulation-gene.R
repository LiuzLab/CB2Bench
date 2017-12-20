m_id <- list()
m_id$CC2_gene.csv <- "gpvalue"
m_id$MAGeCK_gene.csv <- "twosided.p.value"
m_id$PBNPA_gene.csv <- "p.value.twosides"
m_id$sgRSEA_gene.csv <- "p.value.twosides"
m_id$DESeq2_gene.csv <- "pvalue"

all.df <- NULL

file.name <- "cache/sim/*/*_gene.csv"

for(f in Sys.glob(file.name)) {
  m <- basename(f)
  print(m)
  dset <- basename(dirname(f))
  df <- read.csv(f)
  fdr <- df[,m_id[[m]]]
  m <- strsplit(m,"\\_")[[1]][1]
  if(m != "CC2") {
    fdr <- p.adjust(fdr, method="fdr")
  }
  new.df <- data.frame(dataset=dset, method=m, gene=df[,1], fdr=fdr)
  colnames(new.df) <- c("dataset", "method", "gene", "fdr")

  param <- strsplit(gsub("_"," ",gsub("=", " ", dset))," ")[[1]]
  depth <- as.integer(param[4])
  noise <- as.double(param[6])
  effect <- as.double(param[8])
  facs <- as.double(param[2])

  dataset <- load.sim(depth, facs, noise, effect)

  ess <- dataset %>% mutate(essential=ifelse(class=="decreasing"|class=="increasing",1,0)) %>%
    group_by(gene) %>%
    summarise(essential=mean(essential)) %>%
    select(gene, essential)
  new.df <- left_join(new.df, ess, by=c("gene"="gene"))
  new.df$noise <- noise
  new.df$effect <- effect

  print(head(new.df))

  if(is.null(all.df)) {
    all.df <- new.df
  }
  else {
    all.df <- rbind(all.df, new.df)
  }
}
all.df$fdr[is.na(all.df$fdr)] <- 1

lineplot.fdr(all.df)
