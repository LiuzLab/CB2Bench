plot.performance <- function(tidy, type="PR") {

  fun.curve <- roc.curve
  if(type=="PR") fun.curve <- pr.curve

  curve <- tibble()
  auc <- tibble()
  for(m in unique(tidy$methods)) {
    x <- tidy %>% filter(methods == m)
    x[is.na(x$score),"score"] <- -1e8
    fg <- x$score[x$label == 1]
    bg <- x$score[x$label == 0]
    pr <- fun.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    y <- tibble(
      method=m,
      x=pr$curve[,1], y=pr$curve[,2], density=pr$curve[,3]
    )
    curve <- bind_rows(curve, y)
    auc <- bind_rows(auc, tibble(method=m,
                                 AUC=ifelse(is.null(pr$auc.integral), pr$auc, pr$auc.integral)))
  }
  #curve <- arrange(curve, desc(y))

  ret <- list()
  ret$bar.plot <- ggplot(auc, aes(x=method, y=AUC)) +
    geom_bar(stat = "identity", aes(fill=method)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + ylim(0,1)

  ret$curve.plot <- ggplot(curve, aes(x=x, y=y)) +
    #geom_line(aes(colour=method)) +
    #geom_smooth(aes(colour=method),method="loess") + ylim(0,1)
    geom_point(aes(colour=method), alpha=0.5, size=0.5) + ylim(0, 1)
    #geom_path(aes(colour=method), alpha=0.5, size=1) + ylim(0,1)

  ret
}

plot.all <- function(ret, type="PR", title = NULL) {
  lab.x <- "False Positive Rate"
  lab.y <- "True Postivie Rate"
  lab.measure <- "AUCROC"
  if(type=="PR") {
    lab.x <- "Recall"
    lab.y <- "Precision"
    lab.measure <- "AUCPRC"
  }

  gene.plots <- plot.performance(ret$tidy.gene, type)
  sgRNA.plots <- plot.performance(ret$tidy.sgRNA, type)
  p <- plot_grid(sgRNA.plots$bar.plot + labs(title="sgRNA") +
                   labs(y=lab.measure),
            gene.plots$bar.plot + labs(title="gene") +
              theme(axis.title.y= element_blank()),
            sgRNA.plots$curve.plot + labs(x=lab.x, y=lab.y) ,
            gene.plots$curve.plot + labs(x=lab.x, y=lab.y) +
              theme(axis.title.y= element_blank()))
  if(!is.null(title)) {
    title <- ggdraw() + draw_label(title, fontface='bold')
    p <- plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  }
  p
}


# RT112.1 <- plot.all(RT112.ret, "ROC", "RTT12 dataset benchmark (AUCROC)")
# RT112.2 <- plot.all(RT112.ret, "PR", "RTT12 dataset benchmark (AUCPRC)")
#
# UMUC3.1 <- plot.all(UMUC3.ret, "ROC", "UMUC3 dataset benchmark (AUCROC)")
# UMUC3.2 <- plot.all(UMUC3.ret, "PR", "UMUC3 dataset benchmark (AUCPRC)")
#
