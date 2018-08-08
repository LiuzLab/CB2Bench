library(cowplot)
library(ggsci)

plot.dot.bar <- function(gene.name,
                         dataset.name,
                         ncol) {
  df_fdr <- read_csv("inst/extdata/gene_fdr.csv")
  df_count <- read_csv("inst/extdata/read_count.csv")
  plot.bar <-
    df_fdr %>%
    filter(gene==gene.name,
           dataset==dataset.name) %>%
    mutate(FDR=-log10(fdr)) %>%
    ggplot(aes(x=method, y=FDR)) +
    geom_bar(stat = "identity", aes(fill=method)) +
    ylab("-log10(FDR)") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = -log10(0.05), color = "black", alpha=0.5) +
    scale_fill_npg()


  plot.dots <-
    df_count %>%
    filter(gene==gene.name,
           dataset==dataset.name,
           count_type=="countPerMil") %>%
    mutate(CPM=read_count) %>%
    ggplot(aes(x=group, y=CPM)) +
    geom_dotplot(aes(fill=group, color=group),
                 binaxis='y',
                 stackdir='center',
                 stackratio=1.5,
                 dotsize=1.2) +
    facet_wrap(~sgRNA, scales = "free_y", ncol=ncol) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_color_npg() +
    scale_fill_npg() +
    ggtitle(str_c(dataset.name, "-", gene.name))
  plot_grid(plot.dots, plot.bar, ncol = 1)
}

gene.name <- "RPL5"
plots <- list()
for(dataset in c("CRISPR.RT112", "CRISPR.UMUC3", "CRISPRi.RT112")) {
  plots[[dataset]] <- plot.dot.bar(gene.name, dataset, 3)
}

plot_merged <- plot_grid(plotlist = plots, labels = "AUTO", nrow=1)
save_plot("figures/RPL5.png", plot_merged, base_width = 14, base_height = 12)
#save_plot("figures/RPL5.tiff", plot_merged, base_width = 14, base_height = 8)
