heatmap.fdr <- function(all.df, prof.level, order.methods = NULL) {

  all.df$fdr[is.na(all.df$fdr)] <- 1
  heatmap <- list()
  (col.pal <- RColorBrewer::brewer.pal(5, "Reds"))
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
      x %>% mutate(Essential = ifelse(essential > 0, "Essential", "Non-essential")) %>%
      select(prof.level, Essential) %>%
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
          "Essential" = c("Essential" = "#000000", "Non-essential" = "#ffffff")
        ),
        gaps_row = sum(x$essential > 0)
      )
    heatmap[[dset]] <- hm$gtable
  }

  (p <- plot_grid(plotlist = heatmap,  ncol = length(heatmap)))
}
