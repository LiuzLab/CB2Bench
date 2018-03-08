gene_stat <- list()

for(dset in unique(all.df$dataset)) {
  gene_stat[[dset]] <-
    all.df %>% filter(dataset==dset) %>% select(-dataset) %>% spread(method, fdr) %>% arrange(-essential, CC2)
}

library(openxlsx)
write.xlsx(gene_stat, "inst/extdata/gene-level-stats.xlsx")
