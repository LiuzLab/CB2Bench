cor.heatmap <- function(df, dset) {
  df %>%
    filter(dataset==dset) %>%
    group_by(method) %>%
    mutate(rank = rank(fdr)) %>%
    ungroup %>%
    as.data.frame %>%
    select(-fdr) %>%
    spread(method, rank) %>%
    select(-dataset, -essential, -gene) %>%
    cor(method = "spearman") %>%
    pheatmap(display_numbers = T,
             main = dset,
             legend = F,
             breaks = seq(0.7, 1.0, length.out = 100)) %>%
    .$gtable
}

plot_grid(
  all.df %>% cor.heatmap("CRISPR.RT112"),
  all.df %>% cor.heatmap("CRISPR.UMUC3"),
  all.df %>% cor.heatmap("CRISPRi.RT112"),
  nrow = 1,
  labels = "AUTO")
