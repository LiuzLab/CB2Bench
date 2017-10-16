library(readr)

m_id <- list()

m_id$CC2_gene.csv <- "gpvalue"
m_id$MAGeCK_gene.csv <- "twosided.p.value"
m_id$PBNPA_gene.csv <- "p.value.twosides"
m_id$sgRSEA_gene.csv <- "p.value.twosides"

all.df <- NULL

file.name <- "cache/sim/*_F=*/*_gene.csv"

for(f in Sys.glob(file.name)) {
  m <- basename(f)
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
  depth <- as.integer(param[2])
  noise <- as.double(param[4])
  effect <- as.double(param[6])
  facs <- as.double(param[8])
  
  dataset <- load.sim(depth, facs, noise, effect) 
    
  ess <- dataset %>% mutate(essential=ifelse(class=="decreasing"|class=="increasing",1,0)) %>%
    group_by(gene) %>%
    summarise(essential=mean(essential)) %>%
    select(gene, essential)
  new.df <- left_join(new.df, ess, by=c("gene"="gene"))
  
  if(is.null(all.df)) all.df <- new.df
  else {
    all.df <- rbind(all.df, new.df)
  }
}
all.df$fdr[is.na(all.df$fdr)] <- 1

df.gene.plot <- tibble()
for(dset in unique(all.df$dataset)) {
  param <- strsplit(gsub("_"," ",gsub("=", " ", dset))," ")[[1]]
  depth <- as.integer(param[2])
  noise <- as.double(param[4])
  effect <- as.double(param[6])
  facs <- as.double(param[8])
  
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
      #df.gene.plot <- bind_rows(df.gene.plot, tibble(dataset=d, method=mat, FDR=log10(fdr), PR=precision, RC=recall, TP=TP, FP=FP, FN=FN))
      df.gene.plot <- bind_rows(df.gene.plot, tibble(depth=depth, noise=noise, facs=facs, effect=effect, method=mat, FDR=log10(fdr), value=precision, measure="precision", TP=TP, FP=FP, FN=FN))
      df.gene.plot <- bind_rows(df.gene.plot, tibble(depth=depth, noise=noise, facs=facs, effect=effect, method=mat, FDR=log10(fdr), value=recall, measure="recall", TP=TP, FP=FP, FN=FN))
      df.gene.plot <- bind_rows(df.gene.plot, tibble(depth=depth, noise=noise, facs=facs, effect=effect, method=mat, FDR=log10(fdr), value=fmeasure, measure="F-measure", TP=TP, FP=FP, FN=FN))
    }
  }
}
#df.gene.plot$FDR <- -df.gene.plot$FDR

plots <- list()
for(mm in c("recall")) {
  for(ff in c(0.1, 0.25)) { 
    for(ee in c(0.1, 0.2) ) {
      p <- df.gene.plot %>% filter(measure == mm) %>% filter(facs==ff, effect==ee) %>%
        ggplot(aes(x=-FDR, y=value)) + geom_point(aes(colour=method, shape=method), alpha=0.5) +
        geom_line(aes(colour=method), alpha=0.5) + facet_grid(depth~noise) + ylim(0,1) +
        scale_x_discrete(limits=as.character(seq(10,1,-1))) + xlab("FDR") + ylab(mm) + 
        ggtitle(sprintf("FACS = %.2f & effect size = %.1f", ff, ee)) + theme(legend.position="none")
      key <- sprintf("%s_%.2f_%.1f", mm, ff, ee)
      plots[[key]] <- p
    }
  }
}

#grobs <- ggplotGrob(plots$`F-measure_0.10_0.1`)$grobs
#legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


pt <- plot_grid(plot_grid(plotlist=plots, ncol=2) , legend, rel_widths = c(1, .1))
save_plot(pt, filename = "sim-ret-rc.pdf", base_height = 12, base_width = 14)
