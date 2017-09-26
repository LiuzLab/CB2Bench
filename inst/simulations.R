library(CC2Sim)
library(tidyverse)
library(cowplot)
library(plotROC)
results <- expand.grid(noise = seq(4, 6, 0.5),
                       effect = c(0.05, 0.1, 0.2)) %>%
  rowwise() %>%
  mutate(ret = list(run.simulation(depth=10, facs = 0.1, noise, effect)))
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

results.all.tmp <- filter(results.all, noise >= 3 & effect < 0.2)
p <- (ggplot(results.all.tmp, aes(m = -pvalue, d = label, color = methods)) +
        geom_roc(labels = FALSE, pointsize = 0, linealpha=0.7, pointalpha=0.5, size=0.5) +
        facet_grid(sprintf("beta = %.2f", effect)~sprintf("noise = %.2f", noise)) +
        ylab("TPR") + xlab("FPR") + theme(legend.position="bottom") + ylim(0,1) +
        scale_color_brewer(palette = "Set1"))
save_plot("benchmarking/simulation_ROC.pdf", p, base_height=6, base_aspect_ratio = 2.0)
save(file="benchmarking/results.rData", results)

# An example code for investigation the correlation between method A and B
# results.sgRNA %>%
#   left_join(select(sim.dat,sgRNA, class), by="sgRNA") %>%
#   mutate(label = ifelse(class %in% c("increasing", "decreasing"), 1, 0)) %>%
#   #filter( class != "increasing" & class != "decreasing") %>%
#   ggplot(aes(x=-log10(CC2), y=-log10(mageck), colour=as.factor(label))) + geom_point(alpha=0.5)
