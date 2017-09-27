library(CC2Sim)
library(tidyverse)
library(cowplot)
library(plotROC)
results <- expand.grid(noise = seq(4, 6, 0.5),
                       effect = c(0.1, 0.2)) %>%
  rowwise() %>%
  mutate(ret = list(run.simulation(depth=50, facs = 0.10, noise, effect)))
# I need to check that without list is fine for above code

results.all <- data.frame()
for (i in 1:nrow(results)) {
  row <- results[i, ]
  noise <- row$noise
  effect <- row$effect
  tidy <-
    cbind(data.frame(noise = noise, effect = effect, row$ret[[1]]$tidy))
  results.all <- rbind(results.all, tidy)
}

results.all.gene <- data.frame()
for (i in 1:nrow(results)) {
  row <- results[i, ]
  noise <- row$noise
  effect <- row$effect
  tidy <-
    cbind(data.frame(noise = noise, effect = effect, row$ret[[1]]$tidy.gene))
  results.all.gene <- rbind(results.all.gene, tidy)
}


p <- (ggplot(results.all, aes(m = -pvalue, d = label, color = methods)) +
        geom_roc(labels = FALSE, pointsize = 0, linealpha=0.7, pointalpha=0.5, size=0.5) +
        facet_grid(sprintf("beta = %.2f", effect)~sprintf("noise = %.2f", noise)) +
        ylab("TPR") + xlab("FPR") + theme(legend.position="bottom") + ylim(0,1) +
        scale_color_brewer(palette = "Set1"))
p

p2 <- (ggplot(results.all.gene, aes(m = -pvalue, d = label, color = methods)) +
         geom_roc(labels = FALSE, pointsize = 0, linealpha=0.7, pointalpha=0.5, size=0.5) +
         facet_grid(sprintf("beta = %.2f", effect)~sprintf("noise = %.2f", noise)) +
         ylab("TPR") + xlab("FPR") + theme(legend.position="bottom") + ylim(0,1) +
         scale_color_brewer(palette = "Set1"))
p2

save_plot("inst/figures/simulation_seq_depth_50_facs_10_sgRNA.PDF", p, base_height=6, base_aspect_ratio = 2.0)
save_plot("inst/figures/simulation_seq_depth_50_facs_10_gene.PDF", p2, base_height=6, base_aspect_ratio = 2.0)

save_plot("benchmarking/simulation_ROC.pdf", p, base_height=6, base_aspect_ratio = 2.0)
save(file="benchmarking/results.rData", results)

# An example code for investigation the correlation between method A and B
# results.sgRNA %>%
#   left_join(select(sim.dat,sgRNA, class), by="sgRNA") %>%
#   mutate(label = ifelse(class %in% c("increasing", "decreasing"), 1, 0)) %>%
#   #filter( class != "increasing" & class != "decreasing") %>%
#   ggplot(aes(x=-log10(CC2), y=-log10(mageck), colour=as.factor(label))) + geom_point(alpha=0.5)
