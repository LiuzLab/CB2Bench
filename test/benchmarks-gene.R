library(readr)

m_id <- list()

m_id$CC2_gene.csv <- "gpvalue"
m_id$MAGeCK_gene.csv <- "neg.fdr"
m_id$PBNPA_gene.csv <- "neg.fdr"
m_id$ScreenBEAM_gene.csv <- "FDR"
m_id$sgRSEA_gene.csv <- "FDR.neg"
m_id$CC1_gene.csv <- "pvalue"
m_id$DESeq2_gene.csv <- "padj"
all.df <- NULL


load("inst/extdata/nature-biotech.Rdata")

file.name <- "cache/nature-biotech/*/*_gene.csv"
for(f in Sys.glob(file.name)) {
  m <- basename(f)
  dset <- basename(dirname(f))
  df <- read.csv(f)
  print(m)
  print(df)
  fdr <- df[,m_id[[m]]]
  if(m=="CC2_gene.csv") {
    fdr[df$gtvalue<0] <- 1
  } else if(m=="ScreenBEAM_gene.csv") {
    fdr[df$beta>0] <- 1
  } else if(m=="CC1_gene.csv"||m=="DESeq2_gene.csv"){
    fdr <- p.adjust(fdr, method="fdr")
    fdr[df$fc>0] <- 1
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

all.df$fdr[is.na(all.df$fdr)] <- 1

df.gene.plot <- tibble()
for(dset in unique(all.df$dataset)) {
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
      df.gene.plot <- bind_rows(df.gene.plot, tibble(dataset=dset, method=mat, FDR=log10(fdr), value=precision, measure="precision", TP=TP, FP=FP, FN=FN))
      df.gene.plot <- bind_rows(df.gene.plot, tibble(dataset=dset, method=mat, FDR=log10(fdr), value=recall, measure="recall", TP=TP, FP=FP, FN=FN))
      df.gene.plot <- bind_rows(df.gene.plot, tibble(dataset=dset,method=mat, FDR=log10(fdr), value=fmeasure, measure="F-measure", TP=TP, FP=FP, FN=FN))
    }
  }
}

pt <- df.gene.plot %>% filter(measure==measure) %>% filter(FDR>-6) %>%
  ggplot( aes(x=-FDR, y=value)) + geom_point(aes(colour=method, shape=method), alpha=0.5) +
  geom_line(aes(colour=method), alpha=0.5) + facet_grid(measure~dataset) + ylim(0,1) + scale_x_reverse(breaks = seq(1,5,1))+ xlab("FDR") + ylab("measure")
save_plot(pt, filename = "nat-biotech-gene-FDR.pdf", base_height = 4, base_width = 8)

library(pheatmap)
library(RColorBrewer)

heatmap <- list()

(col.pal <- RColorBrewer::brewer.pal(5, "Reds"))
col.pal[1] <- "#FFFFFF"

for(dset in unique(all.df$dataset)) {
  ess <- dataset[[dset]] %>% mutate(essential=ifelse(class=="decreasing",1,0)) %>%
    group_by(gene) %>%
    summarise(essential=mean(essential)) %>%
    select(gene, essential)
  x <- all.df %>% filter(dataset==dset, method != "DESeq2") %>% select(-1,-5) %>% spread(method, fdr) %>%
    remove_rownames()

  x[x<1e-8] <- 1e-8
  x[,-1] <- floor(-log10(x[,-1]))
  x <- x[order(-rowSums(x[,-1])),]
  x <- x %>% left_join(ess, by="gene")
  x <- x[order(-x$essential),]
  x <- x %>% remove_rownames()
  #x <- x[,c(1,3,2,4:ncol(x))]
  x$essential[x$essential==1] <- 8
  tmp <- x %>% mutate(Essential = ifelse(essential>0, "Essential", "Non-essential") ) %>% select(gene, Essential) %>% column_to_rownames("gene")

  hm <- pheatmap(column_to_rownames(x, "gene") %>% #select(-essential) %>%
                   select(CC2, MAGeCK, PBNPA, ScreenBEAM, sgRSEA), scale = "none",
                 cluster_cols = F, cluster_rows = F, main = dset, color = col.pal, #border_color = NA,
                 fontsize_row=8, display_numbers = F, number_format = "%d",
                 fontsize_number = 6, legend = F, show_rownames=F,  annotation_row = tmp, annotation_legend = F,
                 annotation_colors = list("Essential"=c("Essential" = "#000000", "Non-essential" = "#ffffff")),
                 gaps_row = sum(x$essential>0))
  heatmap[[dset]] <- hm$gtable
}


hm <- pheatmap(column_to_rownames(x, "gene") %>% select(-essential) %>%
                 select(CC2, MAGeCK, PBNPA, ScreenBEAM, sgRSEA), scale = "none",
               cluster_cols = F, cluster_rows = F, main = dset, color = col.pal, border_color = NA,
               fontsize_row=8, display_numbers = F, number_format = "%d",
               fontsize_number = 6, legend = T, show_rownames=F, annotation_row = tmp,
               annotation_colors = list("essential"=c("Essential" = "#000000", "Non-essential" = "#ffffff")),
               gaps_row = sum(x$essential>0))

#heatmap$legend <- hm$gtable[[1]][[7]]
(p <- plot_grid(plotlist = heatmap,  ncol = 4))
#save_plot(p, filename = "nat-biotech-gene-heatmap.pdf", base_height=10, base_width=12)



