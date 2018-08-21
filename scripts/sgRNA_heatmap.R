library(tidyverse)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(eulerr)

CRISPR.RT112 <- list()
CRISPR.UMUC3 <- list()
CRISPRi.RT112 <- list()
read_csv("inst/extdata/CRISPR.RT112.csv") %>% rename(gene=Gene) -> CRISPR.RT112$readcount
read_csv("inst/extdata/CRISPR.UMUC3.csv") -> CRISPR.UMUC3$readcount
read_csv("inst/extdata/CRISPRi.RT112.csv") -> CRISPRi.RT112$readcount

rename_cols <- function(obj) {
  obj$readcount %>%
    rename("T0 1" = A1,
           "T0 2" = A2,
           "T0 3" = A3,
           "T1 1" = B1,
           "T1 2" = B2,
           "T1 3" = B3) -> obj$readcount
  obj
}

CRISPR.RT112 %>% rename_cols() -> CRISPR.RT112
CRISPR.UMUC3 %>% rename_cols() -> CRISPR.UMUC3
CRISPRi.RT112 %>% rename_cols() -> CRISPRi.RT112

read_and_merge <- function(obj, screen) {
  sprintf("cache/nature-biotech/%s/CC2_sgRNA.csv", screen) %>% read_csv() -> obj$cc2_sg
  sprintf("cache/nature-biotech/%s/MAGeCK_sgRNA.csv", screen) %>%
    read_csv() -> obj$mageck_sg
  sprintf("cache/nature-biotech/%s/PinAPL-py_sgRNA.csv", screen) %>% read_csv() -> obj$pp_sg
  obj$cc2_sg %>% left_join(obj$mageck_sg,
                           by = c("sgRNA"="sgrna")) %>%
    left_join(obj$pp_sg, by = "sgRNA") -> obj$merged_sg
  obj
}


CRISPR.RT112 %>% read_and_merge("CRISPR.RT112") -> CRISPR.RT112
CRISPR.UMUC3 %>% read_and_merge("CRISPR.UMUC3") -> CRISPR.UMUC3
CRISPRi.RT112 %>% read_and_merge("CRISPRi.RT112") -> CRISPRi.RT112

normalize <- function(obj) {
  for(i in 3:ncol(obj$readcount)) {
    s <- sum(obj$readcount[,i])
    obj$readcount[,i] <- obj$readcount[,i] / s * 10^6
  }
  obj
}

CRISPR.RT112 %>% normalize() -> CRISPR.RT112
CRISPR.UMUC3 %>% normalize() -> CRISPR.UMUC3
CRISPRi.RT112 %>% normalize() -> CRISPR.RT112i

plot_heatmap <- function(obj, df_sg, main_title) {
  if(nrow(df_sg) == 0) {
    return(NULL)
  }
  dplyr::select(df_sg, sgRNA) %>%
    left_join(obj$readcount, by = "sgRNA") %>%
    as.data.frame %>%
    column_to_rownames("sgRNA") %>%
    dplyr::select(-gene) -> df_rc
  df_rc %>%
    `+`(1) %>% log2 %>%
    pheatmap(
      scale = "row",
      #cluster_cols = F,
      #cluster_rows = F,
      show_rownames = F,
      border_color = NA,
      legend = F,
      #annotation_row = df_an,
      clustering_method = "average",
      #clustering_distance_cols = "correlation",
      treeheight_row = 0,
      treeheight_col = 0,
      main = main_title,
      silent = T) #%>% .$gtable
}

generate_figure <- function(obj, cutoff = 0.05) {
  obj$merged_sg %>%
    dplyr::filter(fdr_twosided < cutoff, FDR > cutoff, `p-value (adj.)` > cutoff)  %>%
    plot_heatmap(obj, ., "CC2") -> hm.cc
  obj$merged_sg %>%
    dplyr::filter(fdr_twosided > cutoff, FDR < cutoff, `p-value (adj.)` > cutoff)  %>%
    plot_heatmap(obj, ., "MAGeCK") -> hm.mg
  obj$merged_sg %>%
    dplyr::filter(fdr_twosided > cutoff, FDR > cutoff, `p-value (adj.)` < cutoff)  %>%
    plot_heatmap(obj, ., "PinAPL-Py") -> hm.pp

  obj$merged_sg %>%
    select(CC2 = fdr_twosided,
           MAGeCK = FDR,
           `PinAPL-Py` = `p-value (adj.)`) %>%
    mutate(CC2 = CC2 < cutoff,
           MAGeCK = MAGeCK < cutoff,
           `PinAPL-Py` = `PinAPL-Py` < cutoff) %>% write_csv("tmp.csv")
  obj$merged_sg %>%
    select(CC2 = fdr_twosided,
           MAGeCK = FDR,
           `PinAPL-Py` = `p-value (adj.)`) %>%
    mutate(CC2 = CC2 < cutoff,
           MAGeCK = MAGeCK < cutoff,
           `PinAPL-Py` = `PinAPL-Py` < cutoff) %>%
    euler(shape = "ellipse") %>%
    plot(quantities=T,
         fill = c("#db3236", "#3cba54", "#f4c20d"),
         edges = c("#ff0000", "#00ff00", "#ffff00"), alpha=0.2) -> vd

  plot_grid(vd,
            plot_grid(hm.cc$gtable, hm.mg$gtable,hm.pp$gtable, nrow=1),
            rel_widths = c(1.5,1), scale=c(0.8, 1), nrow = 1)
}

CRISPRi.RT112 %>% generate_figure()

plot_grid(
  CRISPR.RT112 %>% generate_figure(),
  CRISPR.UMUC3 %>% generate_figure(),
  CRISPRi.RT112 %>% generate_figure(),
  ncol=1, labels = "AUTO") %>%
  save_plot("figures/Figure-Heatmap-Evers.png",., base_width = 8, base_height = 12)

dev.off()
par(mfrow=c(2,1))
CRISPR.RT112$merged_sg %>%
  filter(fdr_twosided < 0.01, FDR > 0.01) %>%
  select(log2FC) %>% unlist %>% hist(main = "CC2")

CRISPR.RT112$merged_sg %>%
  filter(fdr_twosided > 0.01, FDR < 0.01) %>%
  select(log2FC) %>% unlist %>% hist(main = "MAGeCK")

CRISPR.UMUC3$merged_sg %>%
  filter(fdr_twosided < 0.01, FDR > 0.01) %>%
  select(log2FC) %>% unlist %>% hist(main = "CC2")

CRISPR.UMUC3$merged_sg %>%
  filter(fdr_twosided > 0.01, FDR < 0.01) %>%
  select(log2FC) %>% unlist %>% hist(main = "MAGeCK")

CRISPRi.RT112$merged_sg %>%
  filter(fdr_twosided < 0.01, FDR > 0.01) %>%
  select(log2FC) %>% unlist %>% hist(main = "CC2")

CRISPRi.RT112$merged_sg %>%
  filter(fdr_twosided > 0.01, FDR < 0.01) %>%
  select(log2FC) %>% unlist %>% hist(main = "MAGeCK")
