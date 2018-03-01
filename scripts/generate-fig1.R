pt.f1 <- list()
all.df <- all.df %>% filter(method!="DESeq2") %>%
  filter(dataset!="shRNA.RT112", dataset!="shRNA.UMUC3")
ct <- c(0.1, 0.05, 0.01, 0.005, 0.001)
for(dset in unique(all.df$dataset)) {
  df.prof <- tibble()
  for(mat in unique(all.df$method)) {
    tmp <- all.df %>% filter(dataset==dset, method==mat)
    for(fdr in ct) {
      #fdr <- 10^fdr
      TP <- sum((tmp$fdr <= fdr) & (tmp$essential == 1))
      FP <- sum((tmp$fdr <= fdr) & (tmp$essential == 0))
      FN <- sum((tmp$fdr > fdr) & (tmp$essential == 1))

      precision <- TP / max(1,(TP+FP))
      recall <- TP / max(1,(TP+FN))
      fmeasure <- 2*(precision*recall)/(precision+recall)
      if(is.na(fmeasure)) fmeasure <- 0
      df.prof <- bind_rows(df.prof, tibble(dataset=dset,method=mat, FDR=fdr, value=fmeasure, measure="F-measure", TP=TP, FP=FP, FN=FN))
      if(dset=="shRNA.RT112" && mat=="CC2") {
        print(c(mat, dset, fdr))
        print(table(tmp$fdr <= fdr, tmp$essential))
        print(c(precision, recall, fmeasure))
      }
    }
  }

  pt.f1[[dset]] <- df.prof %>% mutate(FDR = factor(FDR, levels=ct)) %>%
      ggplot( aes(x=FDR, y=value)) + geom_point(aes(colour=method), size=2) +
      geom_line(aes(colour=method, group=method), size=0.8)+
      xlab("FDR") + ylab("F1-score") + ylim(0,1) +
    scale_color_npg() + theme(axis.text.x = element_text(angle=90)) +
    theme(strip.text.x = element_text(face="bold")) +
    theme(legend.position = "none")

}

legend <- get_legend(pt.f1$CRISPR.RT112 + theme(legend.position = "bottom", legend.justification="center"))
plot_grid(plot_grid(plotlist = pt.f1, nrow=1), legend, ncol=1, rel_heights = c(4,1))


prof.level <- "gene"
order.methods <- c("CC2", "ScreenBEAM", "PBNPA", "PinAPL-py", "sgRSEA", "HitSelect", "MAGeCK")
all.df$fdr[is.na(all.df$fdr)] <- 1
heatmap <- list()
(col.pal <- RColorBrewer::brewer.pal(9, "Reds"))
col.pal[1] <- "#FFFFFF"
for (dset in unique(all.df$dataset)) {

  if(prof.level == "gene") {
    ess <-
      dataset[[dset]] %>% mutate(essential = ifelse(class == "decreasing", 1, 0)) %>%
      group_by(gene) %>%
      summarise(essential = mean(essential)) %>%
      select(gene, essential)
  }
  else {
    ess <-
      dataset[[dset]] %>% mutate(essential = ifelse(class == "decreasing", 1, 0)) %>%
      group_by(sgRNA) %>%
      summarise(essential = mean(essential)) %>%
      select(sgRNA, essential)
  }
  x <-
    all.df %>% filter(dataset == dset) %>% select(-1,-5) %>% spread(method, fdr) %>%
    remove_rownames()

  x[x < 1e-8] <- 1e-8
  x[,-1] <- floor(-log10(x[,-1]))
  x <- x[order(-rowSums(x[,-1])),]
  x <- x %>% left_join(ess, by = prof.level)
  x <- x[order(-x$essential),]
  x <- x %>% remove_rownames()
  x <- x[, c(1, 3, 2, 4:ncol(x))]
  x$essential[x$essential == 1] <- 8
  tmp <-
    x %>% mutate(Essentiality = ifelse(essential > 0, "Essential", "Non-essential")) %>%
    select(prof.level, Essentiality) %>%
    column_to_rownames(prof.level)

  if(is.null(order.methods)) {
    order.methods <- 1:(ncol(x)-1)
  }

  hm <-
    pheatmap(
      column_to_rownames(x, prof.level) %>%
        select(order.methods),
      scale = "none",
      cluster_cols = F,
      cluster_rows = F,
      main = dset,
      color = col.pal,
      legend = F,
      show_rownames = F,
      annotation_row = tmp,
      annotation_legend = F,
      annotation_colors = list(
        "Essentiality" = c("Essential" = "#000000", "Non-essential" = "#ffffff")
      ),
      gaps_row = sum(x$essential > 0)
    )
  heatmap[[dset]] <- hm$gtable
}

legend_heatmap <- generate_heatmap_legend()
pt.merged <- list()
for(d in unique(all.df$dataset)) {
  x <- plot_grid(heatmap[[d]]) + theme(plot.margin = unit(c(12,24,6,40), "pt"))
  if(d == "CRISPR.RT112") {
    pt.merged[[d]] <- plot_grid(x, pt.f1[[d]], ncol=1, rel_heights = c(6,4),
                                labels = "AUTO",
                                label_size = 18)
  } else {
    pt.merged[[d]] <- plot_grid(x, pt.f1[[d]], ncol=1, rel_heights = c(6,4))
  }
}

plot_grid(plotlist = pt.merged, nrow=1)

legend_f1 <- get_legend(pt.f1$CRISPR.RT112 + theme(legend.position = "left"))


fig1 <- plot_grid(plot_grid(plotlist = pt.merged, nrow=1),
                  plot_grid(legend_heatmap, legend_f1, ncol=1),
                  rel_widths = c(8.5,1.5))

save_plot(filename = "figures/fig1-heatmap-f1.tiff", fig1, base_width = 10, base_height = 8)

