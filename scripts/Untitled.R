library(CRISPhieRmix)
library(readxl)
library(tidyverse)
library(DESeq2)
df_crisprn_a375_count <- read_xlsx("~/Projects/InProgress/CB2_Revision/data/Brunello.xlsx", sheet = 2, skip = 1)
df_crisprn_a375_info <- read_xlsx("~/Projects/InProgress/CB2_Revision/data/Brunello.xlsx", sheet = 3)


left_join(df_crisprn_a375_count, df_crisprn_a375_info %>% select(1,2), by = "sgRNA Sequence") %>%
  select(1, 6, 2:5) %>%
  dplyr::rename(sgRNA = "sgRNA Sequence") %>%
  dplyr::rename(Gene = "Annotated Gene Symbol") -> df_raw


es_genes <- scan("~/Projects/InProgress/CB2_Revision/data/essential-genes.txt", what = character())

left_join(df_crisprn_a375_count, df_crisprn_a375_info %>% select(1,2), by = "sgRNA Sequence") %>%
  dplyr::rename(sgRNA = "sgRNA Sequence") %>%
  dplyr::rename(Gene = "Annotated Gene Symbol") %>%
  mutate(Gene = str_replace_all(Gene, "_", "-")) %>%
  unite("sgRNA", c("Gene", "sgRNA"), sep = "_") %>%
  as.data.frame %>%
  column_to_rownames("sgRNA") -> df_crisprn_a375_count


df_coldata <- data.frame(name = colnames(df_crisprn_a375_count), group = c(rep("ctl", 1), rep("trt", 3)))


dds <- DESeqDataSetFromMatrix(countData = df_crisprn_a375_count, colData = df_coldata, design = ~group)
dds <- DESeq(dds)
df_ret <- results(dds) %>% as.data.frame
df_ret$log2FoldChange[is.na(df_ret$log2FoldChange)] <- 0
neg_ctrl <- df_ret %>% row.names() %>% startsWith("NO-")
ret <- CRISPhieRmix(
  x = df_ret$log2FoldChange[!neg_ctrl],
  geneIds = df_raw$Gene[!neg_ctrl] %>% as.factor,
  negCtrl = df_ret$log2FoldChange[neg_ctrl],
  mu = 0,
  VERBOSE = TRUE, PLOT = TRUE)

intersect(
  es_genes,
  ret$genes[ret$FDR<0.01]
)

pROC::roc(ret$genes %in% es_genes, ret$FDR, auc=TRUE)

sum(ret$FDR<=0.1)
