run.simulation <- function(depth = 10,
                           facs = 0.1,
                           noise = 2.0,
                           effect = 0.05) {
  cat(
    "current parameter:",
    sprintf("depth = %d, facs = %f, noise = %f, effect = %f\n", depth, facs, noise, effect)
  )
  methods <- list(
    mageck = run.mageck,
    DESeq2 = run.DESeq2,
    edgeR = run.edgeR,
    sgRSEA = run.sgRSEA,
    #screenBEAM = run.ScreenBEAM,
    CC2 = run.mbttest
  )
  sim.dat <- load.sim(depth, facs, noise, effect)
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

  ret$tidy <- (tidy.sgRNA.summary <- df.sgRNA.summary %>%
                 gather(methods, pvalue,-sgRNA,-label))

  ret$plot <-
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

  ret$plot.gene <- ggplot(tidy.gene.summary, aes(m=-pvalue, d=label, color= methods)) +
    geom_roc(labels=FALSE)

  ret$tidy.gene <- tidy.gene.summary
  ret
}

