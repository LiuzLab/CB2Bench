df_count <- read_csv("inst/extdata/CRISPR.RT112_essentiality.csv") %>%
  select(-X) %>% gather(sample_id, read_count, -gene, -class, -sgRNA)

ggplot(df_count, aes(read_count)) +
  geom_histogram(aes(fill=class)) +
  facet_wrap(~sample_id, ncol=1)

ggplot(df_count, aes(read_count)) +
  geom_histogram() +
  facet_wrap(~sample_id, ncol=1)

ggplot(df_count, aes(read_count)) +
  geom_histogram(aes(fill=sample_id)) +
  facet_wrap(~class)

library(GGally)

read_csv("inst/extdata/CRISPR.RT112_essentiality.csv") %>%
  filter(class == "inactive") %>%
  select(A1,A2,A3,B1,B2,B3) %>% log10() %>% ggpairs()


df_data <- read_csv("inst/extdata/CRISPR.RT112_essentiality.csv")
df_sg <- read_csv("cache/nature-biotech/CRISPR.RT112/CC2_sgRNA.csv") %>% select(-gene)
df_data %>% left_join(df_sg, by="sgRNA") %>% filter(class=="decreasing") %>%
  ggplot(aes(p_value_neg)) + geom_histogram()

sdf_count %>% group_by(sample_id) %>% summarise(library_size=sum(read_count))
