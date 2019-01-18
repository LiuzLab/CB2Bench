library(tidyverse)
library(ggsci)
library(pheatmap)
library(cowplot)

draw <- function(screen) {
  sprintf("cache/nature-biotech/%s/CC2_gene.csv", screen) %>% read_csv() %>%
    mutate(method = "CC2", rank = rank(p_value_neg)) %>% select(method, gene, rank) -> df_rank

  sprintf("cache/nature-biotech/%s/ScreenBEAM_gene.csv", screen) %>% read_csv() %>%
    mutate(method = "ScreenBEAM", rank = rank(zstat)) %>% select(method, gene, rank) %>%
    bind_rows(df_rank) -> df_rank


  sprintf("cache/nature-biotech/%s/PinAPL-Py_gene.csv", screen) %>% read_csv() %>%
    mutate(method = "PinAPL-Py", rank = rank(-`Fisher Statistic`)) %>% select(method, gene, rank) %>%
    bind_rows(df_rank) -> df_rank

  sprintf("cache/nature-biotech/%s/HitSelect_gene.csv", screen) %>% read_csv() %>%
    mutate(method = "HitSelect", rank = rank(fdr)) %>% select(method, gene, rank) %>%
    bind_rows(df_rank) -> df_rank

  sprintf("cache/nature-biotech/%s/MAGeCK_gene.csv", screen) %>% read_csv() %>%
    mutate(method = "MAGeCK", rank = neg.rank) %>% select(method, gene=id, rank) %>%
    bind_rows(df_rank) -> df_rank


  sprintf("cache/nature-biotech/%s/PBNPA_gene.csv", screen) %>% read_csv() %>%
    mutate(method = "PBNPA", rank = rank(neg.pvalue)) %>% select(method, gene=Gene, rank) %>%
    bind_rows(df_rank) -> df_rank


  sprintf("cache/nature-biotech/%s/CRISPhieRmix_gene.csv", screen) %>% read_csv() %>%
    mutate(method = "CRISPhieRmix", rank = rank(FDR)) %>% select(method, gene, rank) %>%
    bind_rows(df_rank) -> df_rank

  df_rank %>% spread(method, rank) -> df_matrix


  df_matrix %>% select(-gene) %>% cor %>% pheatmap(display_numbers = T, silent=F, fontsize_number = 12) %>% .$gtable -> hm

  load("~/Projects/InProgress/CC2-manuscript-materials/CC2Bench/inst/extdata/nature-biotech.Rdata")

  dataset[[screen]] %>%
    select(gene, class) %>%
    group_by(gene) %>%
    summarise(essential = mean(class)) %>%
    ungroup() %>%
    mutate(essential = c("non-essential", "essential")[essential+1]) -> df_info

  df_rank %>% left_join(df_info, by = "gene") %>%
    ggplot(aes(x=method, y=rank)) +
    geom_boxplot(aes(fill=essential), width = 0.25, alpha=0.7) +
    geom_jitter(aes(fill=essential, color=essential), width = 0.25, alpha=0.5) + scale_y_reverse() +
    scale_fill_npg() + scale_color_npg() + ggtitle(screen) +
    theme(axis.title.x=element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) -> gg_bo

  if(screen == "CRISPR.RT112")
    plot_grid(gg_bo, hm, ncol=1, labels ="AUTO")
  else
    plot_grid(gg_bo, hm, ncol=1)
}


plot_grid(
  "CRISPR.RT112" %>% draw,
  "CRISPR.UMUC3" %>% draw,
  "CRISPRi.RT112" %>% draw,
  nrow = 1
) %>% save_plot(filename = "figures/rank-box.png", ., base_width = 16, base_height = 9)
