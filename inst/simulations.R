library(CC2Sim)
library(tidyverse)
library(cowplot)
library(PRROC)

load.sim <- function(depth, facs, noise, effect) {
  raw.url <- "https://raw.githubusercontent.com/hyunhwaj/Crispulator.jl/master/simulation/matrix/scenario_%d_%.2f_%.2f_%.2f.csv"
  url <- sprintf(raw.url, depth, facs, noise, effect)
  read.csv(url)
}


read_delim("inst/extdata/twosided-expr.tsv", delim="\t")

plot.AUPRC <- function(tidy) {
  curve <- tibble()
  prv <- tibble()
  for(m in unique(tidy$methods)) {
    x <- tidy %>% filter(methods == m)
    x[is.na(x$score),"score"] <- -1e+8
    fg <- x$score[x$label == 1]
    bg <- x$score[x$label == 0]
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
        results.gene <- df.ret$gene
      }
      else {
        results.gene <- dplyr::left_join(results.gene, df.ret$gene, by = "gene")
      }
      nc <- ncol(results.gene)
      colnames(results.gene)[nc] <- i
    }
    print(head(results.gene))
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
  ret$sgRNA <- df.sgRNA.summary

  ret$tidy.sgRNA <- (tidy.sgRNA.summary <- df.sgRNA.summary %>%
                       gather(methods, score,-sgRNA,-label))

  ret$plot.roc.sgRNA <-
    ggplot(tidy.sgRNA.summary, aes(
      m = score,
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
    gather(methods, score, -gene, -label)
  ret$tidy.gene <- tidy.gene.summary
  ret$gene <- df.gene.summary

  ret$plot.roc.gene <- ggplot(tidy.gene.summary, aes(m=score, d=label, color= methods)) +
    geom_roc(labels=FALSE)

  AUPRC.sgRNA <- plot.AUPRC(ret$tidy.sgRNA)
  AUPRC.gene <- plot.AUPRC(ret$tidy.gene)
  ret$plot.prc.sgRNA <- AUPRC.sgRNA$curve.plot
  ret$plot.prc.gene <- AUPRC.gene$curve.plot
  ret
}



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
    print(m)
    x <- tidy %>% filter(methods == m)
    print(head(x))
    x[is.na(x$score),"score"] <- 1
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
