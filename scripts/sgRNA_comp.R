library(tidyverse)
library(pheatmap)
library(eulerr)

read_count <- read_csv("inst/extdata/CRISPR.RT112.csv")

cc <- read_csv("cache/nature-biotech/CRISPR.RT112/CC2_sgRNA.csv")
mg <- read_csv("cache/nature-biotech/CRISPR.RT112/MAGeCK_sgRNA.csv")

merge <- left_join(cc, mg, by = c("sgRNA" = "sgrna"))
merge$fdr_low <- merge$p.low
merge %>% filter(p_value_neg < 0.1 | fdr_low < 0.1) %>% filter(p_value_neg < 0.1, fdr_low > 0.1) %>% .$sgRNA -> cc_only
merge %>% filter(p_value_neg < 0.1 | fdr_low < 0.1) %>% filter(p_value_neg > 0.1, fdr_low < 0.1) %>% .$sgRNA -> mg_only

for(i in 3:8) {
  s <- sum(read_count[,i])
  read_count[,i] <- read_count[,i] / s * 10^6
}

read_count %>% filter(sgRNA %in% cc_only) %>% select(-Gene) %>% column_to_rownames("sgRNA") %>%
  pheatmap(cluster_cols = F, scale="row", main = "MAGeCK Only", cellwidth = 25) -> hm.cc

read_count %>% filter(sgRNA %in% mg_only) %>% select(-Gene) %>% column_to_rownames("sgRNA") %>%
  pheatmap(cluster_cols = F, scale="row", main = "MAGeCK Only", cellwidth = 25) -> hm.mg


merge %>% filter(p_value_neg < 0.1 | fdr_low < 0.1) %>% select(CC2 = p_value_neg, MAGeCK = fdr_low) %>%
  mutate(CC2 = CC2 < 0.1,
         MAGeCK = MAGeCK < 0.1) %>% euler %>% plot(quantities=T, lty=1) -> plot.venn

plot_grid(plot.venn, hm.mg$gtable)

