plot.AUPRC <- function(tidy) {
  curve <- tibble()
  prv <- tibble()
  for(m in unique(tidy$methods)) {
    x <- tidy %>% filter(methods == m)
    x[is.na(x$pvalue),"pvalue"] <- 1
    fg <- 1-x$pvalue[x$label == 1]
    bg <- 1-x$pvalue[x$label == 0]
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    y <- tibble(
      method=m,
      x=pr$curve[,1], y=pr$curve[,2], density=pr$curve[,3]
    )
    print(y)
    curve <- bind_rows(curve, y)
    prv <- bind_rows(prv, tibble(method=m,
                                 AUPRC=pr$auc.integral))
  }

  #curve <- curve %>% arrange(desc(y))

  ret <- list()
  ret$bar.plot <- ggplot(prv, aes(x=method, y=AUPRC)) +
    geom_bar(stat = "identity", aes(fill=method)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  ret$curve.plot <- ggplot(curve, aes(x=x, y=y)) +
    geom_step(aes(colour=method)) +
    ylab("Precision") + xlab("Recall") + ylim(0,1)

  ret
}


run <- function(dat, methods = list(
    mageck = run.mageck,
    DESeq2 = run.DESeq2,
    edgeR = run.edgeR,
    sgRSEA = run.sgRSEA,
    PBNA = run.PBNPA,
    #screenBEAM = run.ScreenBEAM,
    CC2 = run.mbttest
  )) {

  sim.dat <- dat
  results.sgRNA <- NULL
  results.gene <- NULL
  for (i in names(methods)) {
    cat("Running", i, "...", "\n")
    df.ret <- methods[[i]](sim.dat)
    #df.ret$sgRNA$pvalue <- p.adjust(df.ret$sgRNA$pvalue, method="fdr")
    if(!is.null(df.ret$sgRNA)) {
      if (is.null(results.sgRNA)) {
        results.sgRNA <- df.ret$sgRNA
      } else {
        results.sgRNA <- dplyr::left_join(results.sgRNA, df.ret$sgRNA, by = "sgRNA")
      }
      nc <- ncol(results.sgRNA)
      colnames(results.sgRNA)[nc] <- i
    }

    if (!is.null(df.ret$gene)) {
      if (is.null(results.gene)) {
        #df.ret$gene$pvalue <- p.adjust(df.ret$gene$pvalue, method="fdr")
        results.gene <- df.ret$gene
      }
      else {
        results.gene <- dplyr::left_join(results.gene, df.ret$gene, by = "gene")
      }
      nc <- ncol(results.gene)
      colnames(results.gene)[nc] <- i
    }
  }

  ret <- list()
  ret$df <- (
    df.sgRNA.summary <- results.sgRNA %>%
      dplyr::left_join(select(sim.dat, sgRNA, class), by = "sgRNA") %>%
      mutate(label = ifelse(
        class %in% c("increasing", "decreasing"), 1, 0
      )) %>%
      select(-class)
  )

  ret$tidy.sgRNA <- (tidy.sgRNA.summary <- df.sgRNA.summary %>%
                       gather(methods, pvalue,-sgRNA,-label))

  ret$plot.roc.sgRNA <-
    ggplot(tidy.sgRNA.summary, aes(
      m = -pvalue,
      d = label,
      color = methods
    )) +
    geom_roc(labels = FALSE)

  # Not think about gene level, yet
  tmp <- sim.dat %>% group_by(gene) %>%
    filter(row_number()==1) %>% select(gene, class)
  df.gene.summary <- results.gene %>%
    left_join(select(tmp, gene, class), by="gene") %>%
    mutate(label = ifelse(class %in% c("increasing", "decreasing"), 1, 0)) %>%
    select(-class)

  tidy.gene.summary <- df.gene.summary %>%
    gather(methods, pvalue, -gene, -label)
  ret$tidy.gene <- tidy.gene.summary

  ret$plot.roc.gene <- ggplot(tidy.gene.summary, aes(m=-pvalue, d=label, color= methods)) +
    geom_roc(labels=FALSE)

  AUPRC.sgRNA <- plot.AUPRC(ret$tidy.sgRNA)
  AUPRC.gene <- plot.AUPRC(ret$tidy.gene)
  ret$plot.prc.sgRNA <- AUPRC.sgRNA$curve.plot
  ret$plot.prc.gene <- AUPRC.gene$curve.plot
  ret
}
