library(DESeq2)
library(tidyverse)
library(stringr)

tidy_count <- function(df_count) {
  df_count %>%
    gather(sample_id, read_count, -sgRNA, -gene) %>%
    mutate(group = str_sub(sample_id, 1,2)) %>%
    select(gene, sgRNA, group, sample_id, read_count)
}

size_factor_norm_count <- function(csv_path) {
  df_count <- read_csv(csv_path)
  df_deseq2 <- df_count %>%
    select(-gene) %>%
    as.data.frame %>%
    column_to_rownames("sgRNA")
  colnames(df_deseq2) <-
    c("T0_1", "T0_2", "T0_3", "T1_1", "T1_2", "T1_3")
  df_condition <- data.frame(
    condition = c(rep("T0", 3), rep("T1", 3)),
    sample_id = colnames(df_deseq2)
  )

  dds <-
    DESeqDataSetFromMatrix(df_deseq2,
                           df_condition,
                           ~ condition)
  dds <- DESeq(dds)
  df_norm <- counts(dds, normalized=TRUE) %>%
    as.data.frame %>%
    rownames_to_column("sgRNA")

  select(df_count, sgRNA, gene) %>%
    left_join(df_norm, by="sgRNA") %>% tidy_count
}

raw_count <- function(csv_path) {
  df_count <- read_csv(csv_path)
  colnames(df_count)[3:8] <-
    c("T0_1", "T0_2", "T0_3", "T1_1", "T1_2", "T1_3")
  df_count %>% tidy_count
}

cpm_norm_count  <- function(csv_path) {
  df_count <- read_csv(csv_path)
  colnames(df_count)[3:8] <-
    c("T0_1", "T0_2", "T0_3", "T1_1", "T1_2", "T1_3")
  t_count <- df_count %>% tidy_count
  t_count %>%
    left_join(
      t_count %>% group_by(sample_id) %>%
        summarise(total_count=sum(read_count)),
      by="sample_id") %>%
    mutate(read_count = read_count * 1. / total_count * 1e6) %>%
    select(-total_count)
}

count_functions <-
  list(
    "sizeFactor" = size_factor_norm_count,
    "countPerMil" = cpm_norm_count,
    "raw" = raw_count
  )

df_count <- tibble()
for(csv_path in Sys.glob("inst/extdata/*.csv")) {
  dataset <- basename(csv_path) %>% str_remove(".csv")
  for(cnt_func in names(count_functions)) {
    df_count <- bind_rows(
      df_count,
      cbind(data.frame(dataset=dataset, count_type=cnt_func),
            count_functions[[cnt_func]](csv_path))
    )
  }
}

df_count %>% write_csv("../CC2Sim-shiny-report/data/read_count.csv")
