library(precrec)
# The ggplot2 package is required
library(ggplot2)
library(gghighlight)

draw_curve <- function(eval_ret, ctype = "PRC") {
  df_auprc <- auc(eval_ret) %>% filter(curvetypes == ctype) %>%
    dplyr::rename(modname = modnames) %>%
    mutate(method = sprintf("%s (%.3f)", modname, aucs))

  df_pr <- eval_ret %>% as.data.frame %>% filter(type==ctype) %>% left_join(df_auprc, by="modname")
  df_pr %>%
    ggplot(aes(x=x,y=y)) +
    geom_line(aes(color=method), size = 1) +
    geom_line(data = df_pr %>% filter(type==ctype, startsWith(method, "CB2")) ,
              aes(x=x,y=y,color=method), size=0.1) +
    theme_minimal() + scale_color_npg()
}

dset <- "CRISPR.UMUC3"

scores <- join_scores(
  all.df %>% filter(dataset==dset, method == "CC2") %>% mutate(fdr = 1-fdr) %>% pull(fdr),
  all.df %>% filter(dataset==dset, method == "MAGeCK") %>% mutate(fdr = 1-fdr) %>% pull(fdr),
  all.df %>% filter(dataset==dset, method == "CRISPhieRmix") %>% mutate(fdr = 1-fdr) %>% pull(fdr),
  all.df %>% filter(dataset==dset, method == "ScreenBEAM")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  all.df %>% filter(dataset==dset, method == "PBNPA")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  #all.df %>% filter(dataset==dset, method == "sgRSEA")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  all.df %>% filter(dataset==dset, method == "HitSelect")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  all.df %>% filter(dataset==dset, method == "PinAPL-py")%>% mutate(fdr = 1-fdr) %>%  pull(fdr),
  chklen = FALSE
)

labels <- join_labels(
  all.df %>% filter(dataset==dset, method == "CC2") %>% pull(essential),
  all.df %>% filter(dataset==dset, method == "MAGeCK") %>% pull(essential),
  all.df %>% filter(dataset==dset, method == "CRISPhieRmix") %>% pull(essential),
  all.df %>% filter(dataset==dset, method == "ScreenBEAM") %>% pull(essential),
  all.df %>% filter(dataset==dset, method == "PBNPA")%>% pull(essential),
  #all.df %>% filter(dataset==dset, method == "sgRSEA")%>% pull(essential),
  all.df %>% filter(dataset==dset, method == "HitSelect")%>% pull(essential),
  all.df %>% filter(dataset==dset, method == "PinAPL-py")%>% pull(essential),

  chklen = FALSE

)

mmdat <- mmdata(scores, labels, modnames = c("CC2", "MAGeCK", "CRISPhieRmix", "ScreenBEAM", "PBNPA",
                                             #"sgRSEA",
                                             "HitSelect", "PinAPL-py"))

eval_ret <- evalmod(mmdat)


p1 <- draw_curve(eval_ret, "PRC") + xlab("Recall") + ylab("Precision") + gghighlight(use_direct_label=F) + facet_wrap(~method, nrow=1) + ggtitle("PR Curves")
p2 <- draw_curve(eval_ret, "ROC") + xlab("1- Specificity") + ylab("Sensitivity") + gghighlight(use_direct_label=F) + facet_wrap(~method, nrow=1) + ggtitle("ROC Curves")

title <- ggdraw() + draw_label("Evers et al.'s CRISPR library (UMUC3)", fontface='bold')

plot_grid(title, plot_grid(p1, p2, ncol=1, labels = "auto"), rel_heights = c(0.1,1), ncol=1)



