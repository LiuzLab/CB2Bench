library(precrec)
# The ggplot2 package is required
library(ggplot2)

dset <- "CRISPR.RT112"

scores <- join_scores(
  all.df %>% filter(dataset==dset, method == "CC2") %>% mutate(fdr = 1-fdr) %>% pull(fdr),
  all.df %>% filter(dataset==dset, method == "MAGeCK") %>% mutate(fdr = 1-fdr) %>% pull(fdr),
  all.df %>% filter(dataset==dset, method == "ScreenBEAM")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  all.df %>% filter(dataset==dset, method == "PBNPA")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  all.df %>% filter(dataset==dset, method == "sgRSEA")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  all.df %>% filter(dataset==dset, method == "HitSelect")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  all.df %>% filter(dataset==dset, method == "PinAPL-py")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  chklen = FALSE
)

labels <- join_labels(
  all.df %>% filter(dataset==dset, method == "CC2") %>% pull(essential),
  all.df %>% filter(dataset==dset, method == "MAGeCK") %>% pull(essential),
  all.df %>% filter(dataset==dset, method == "ScreenBEAM") %>% pull(essential),
  all.df %>% filter(dataset==dset, method == "PBNPA")%>% pull(essential),
  all.df %>% filter(dataset==dset, method == "sgRSEA")%>% pull(essential),
  all.df %>% filter(dataset==dset, method == "HitSelect")%>% pull(essential),
  all.df %>% filter(dataset==dset, method == "PinAPL-py")%>% pull(essential),

  chklen = FALSE

)

mmdat <- mmdata(scores, labels, modnames = c("CC2", "MAGeCK", "ScreenBEAM", "PBNPA", "sgRSEA", "HitSelect", "PinAPL-py"))

evalmod(mmdat) %>% autoplot(curvetype="PRC", alpha=0.1) + scale_color_npg()
evalmod(mmdat) %>% autoplot(curvetype="ROC", alpha=0.1) + scale_color_npg()

#+ facet_grid(modname~.)
