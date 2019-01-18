library(CRISPhieRmix)
library(tidyverse)
library(DESeq2)
library(glue)
load("inst/extdata/nature-biotech.Rdata")

run_CRISPhieRmix <- function(dat_name) {
  dat <- dataset[[dat_name]]
  dat %>% dplyr::select(-gene, -X, -class)%>% column_to_rownames("sgRNA") -> df_count
  df_coldata <- data.frame(name = colnames(df_count), group = c(rep("ctl", 3), rep("trt", 3)))

  dat$gene <- factor(dat$gene, levels = unique(dat$gene))
  es_genes <- dat %>% filter(class=="decreasing") %>% pull(gene) %>% unique

  dds <- DESeqDataSetFromMatrix(countData = df_count, colData = df_coldata, design = ~group)
  dds <- DESeq(dds)
  df_ret <- results(dds) %>% as.data.frame
  all((df_ret %>% rownames()) == dat$sgRNA)

  ret <- CRISPhieRmix(x = df_ret$log2FoldChange,
                      geneIds = dat$gene,
                      VERBOSE = TRUE,
                      PLOT = TRUE)

  data.frame(gene=ret$genes, FDR = ret$FDR) %>%
    write_csv(glue("cache/nature-biotech/{dat_name}/CRISPhieRmix_gene.csv"))
}

sapply( c("CRISPR.RT112", "CRISPR.UMUC3", "CRISPRi.RT112"), function(i) run_CRISPhieRmix(i) )
