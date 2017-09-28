library(CC2Sim)
library(tidyverse)
library(cowplot)
library(plotROC)
df.tmp <- expand.grid(noise = seq(1, 6, 1),
                       effect = c(0.1, 0.2))

results <- tibble()
for(i in 1:nrow(df.tmp)) {
  noise <- df.tmp[i,]$noise
  effect <- df.tmp[i,]$effect
  results <- bind_rows(results,
                       tibble(noise=noise,
                       effect=effect,
                       ret=list(run.simulation(depth=100, facs = 0.25, noise, effect))))
}

curve <- tibble()
pr.sgRNA <- tibble()
for (i in 1:nrow(results)) {
  row <- results[i, ]
  noise <- row$noise
  effect <- row$effect
  print(c(noise, effect))
  tidy <- row$ret[[1]]$tidy
  for(m in unique(tidy$methods)) {
    x <- tidy %>% filter(methods == m)
    x[is.na(x$pvalue),"pvalue"] <- 1
    fg <- 1-x$pvalue[x$label == 1]
    bg <- 1-x$pvalue[x$label == 0]
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

    y <- tibble(
      noise=noise,
      effect=effect,
      method=m,
      x=pr$curve[,1], y=pr$curve[,2], density=pr$curve[,3]
    )
    curve <- bind_rows(curve, y)
    pr.sgRNA <- bind_rows(pr.sgRNA, tibble(noise=noise,
                                         effect=effect,
                                         method=m,
                                         AUPRC=pr$auc.integral))
  }
}
p1 <- ggplot(pr.sgRNA, aes(x=method, y=AUPRC)) +
  geom_bar(stat = "identity", aes(fill=method)) +
  facet_grid(effect~noise) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p1

p <- ggplot(curve, aes(x=x, y=y)) +
  geom_line(aes(colour=method)) +
  facet_grid(effect~noise) +
  ylab("Precision") + xlab("Recall")
p

curve.gene <- tibble()
pr.gene <- tibble()
for (i in 1:nrow(results)) {
  row <- results[i, ]
  noise <- row$noise
  effect <- row$effect
  print(c(noise, effect))
  tidy <- row$ret[[1]]$tidy.gene
  for(m in unique(tidy$methods)) {
    x <- tidy %>% filter(methods == m)
    x[is.na(x$pvalue),"pvalue"] <- 1
    fg <- 1-x$pvalue[x$label == 1]
    bg <- 1-x$pvalue[x$label == 0]
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    y <- tibble(
      noise=noise,
      effect=effect,
      method=m,
      x=pr$curve[,1], y=pr$curve[,2], density=pr$curve[,3]
    )
    curve.gene <- bind_rows(curve.gene, y)
    pr.gene <- bind_rows(pr.gene, tibble(noise=noise,
                                          effect=effect,
                                          method=m,
                                          AUPRC=pr$auc.integral))
  }
}

p22 <- ggplot(pr.gene, aes(x=method, y=AUPRC)) +
  geom_bar(stat = "identity", aes(fill=method)) +
  facet_grid(effect~noise) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2 <- ggplot(curve.gene, aes(x=x, y=y)) +
  geom_line(aes(colour=method)) +
  facet_grid(effect~noise) + ylab("Precision") + xlab("Recall")
p2

plots <- plot_grid(p1, p22, p, p2, ncol=2, labels=c("A","B","C", "D"))

save_plot("inst/simulation_depth_100.pdf", plots,
          base_height=8, base_aspect_ratio = 3
          )

# An example code for investigation the correlation between method A and B
# results.sgRNA %>%
#   left_join(select(sim.dat,sgRNA, class), by="sgRNA") %>%
#   mutate(label = ifelse(class %in% c("increasing", "decreasing"), 1, 0)) %>%
#   #filter( class != "increasing" & class != "decreasing") %>%
#   ggplot(aes(x=-log10(CC2), y=-log10(mageck), colour=as.factor(label))) + geom_point(alpha=0.5)
