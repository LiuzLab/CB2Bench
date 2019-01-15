library(CRISPhieRmix)
library(tidyverse)
library(DESeq2)
load("inst/extdata/nature-biotech.Rdata")

RT112 <- dataset$CRISPR.RT112 %>% as.data.frame
RT112 %>% select(-gene, -X, -class) %>% as.data.frame() %>% column_to_rownames("sgRNA") -> df_count
df_coldata <- data.frame(name = colnames(df_count), group = c(rep("ctl", 3), rep("trt", 3)))

dds <- DESeqDataSetFromMatrix(countData = df_count, colData = df_coldata, design = ~group)
dds <- DESeq(dds)
df_ret <- results(dds) %>% as.data.frame
all((df_ret %>% rownames()) == RT112$sgRNA)

ret <- CRISPhieRmix(x = df_ret$log2FoldChange, geneIds = RT112$gene %>% as.factor, VERBOSE = TRUE, PLOT = TRUE, max_iter = 1000)

essential <- RT112 %>% filter(class == "decreasing") %>% pull(gene) %>% unique
all_genes <- RT112 %>% pull(gene) %>% unique

library(precrec)
library(ggplot2)

mmdat <- mmdata(ret$FDR, ret$genes %in% essential)

evalmod(mmdat) %>% autoplot()

sum(ret$FDR < 0.1)

gamma.CRISPhieRmixROC = pROC::roc(ret$genes %in% essential,
                                  ret$FDR, auc = TRUE)
gamma.CRISPhieRmixROC
plot(gamma.CRISPhieRmixROC)


data.frame(fdr=ret$FDR, ess = ret$genes %in% essential) %>% ggplot(aes(x=ess,y=-log10(1e-10+fdr))) + geom_boxplot(aes(color=ess))  + geom_jitter()
