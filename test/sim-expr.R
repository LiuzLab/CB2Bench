library(CC2Sim)
run.mageck <- function(dat) {
  tmp.fname <- tempfile(pattern = "file", tmpdir = tempdir())

  dat.mageck <- dat[, c(3,2, 5:ncol(dat))]
  write.table(file = tmp.fname, dat.mageck, sep="\t", quote = F, row.names = F)
  mageck.cmd <- "mageck"

  nx <- ncol(dat)-4
  control.samples <- colnames(dat)[5:(5+nx/2-1)]
  case.samples <- colnames(dat)[(5+nx/2):(5+nx-1)]

  treatment.id <- paste0( case.samples, collapse="," )
  control.id <- paste0( control.samples, collapse="," )

  out.dir <- paste0(tempdir(),"/", "mageck")
  cmd  <- paste(mageck.cmd, "test", "-k", tmp.fname, "-t", treatment.id, "-c", control.id,
                "-n", out.dir)


  system(cmd)
  df.gene <- read.table(paste0(out.dir, ".gene_summary.txt"),
                        sep="\t", row.names = NULL, head=T) %>%
    mutate(twosided.p.value=pmin(1,pmin(neg.p.value, pos.p.value)*2))



  df.sgRNA <- read.table(paste0(out.dir, ".sgrna_summary.txt"),
                         sep="\t", row.names = NULL, head=T)
  list("sgRNA"=df.sgRNA, "gene"=df.gene)
}

run.mbttest <- function(dat) {
  nx <- (ncol(dat)-4)
  df.sgRNA <- mbetattest(X=dat, nci=4, na=nx/2, nb=nx/2, alpha=0.05, level="sgRNA")
  df.gene <- mbetattest(X=dat, nci=4, na=nx/2, nb=nx/2, alpha=0.05, level="gene")
  list("gene"=df.gene, "sgRNA"=df.sgRNA)
}

run.DESeq2<-function(dat){
  df.deseq2<-dat[,-(1:4)]

  nx <- ncol(df.deseq2)
  col.data <- data.frame(row.names = colnames(df.deseq2),
                         condition=c(rep("T0", nx/2),rep("T1", nx/2)))
  col.data$condition <- factor(col.data$condition, levels = c("T0", "T1"))

  row.names(df.deseq2)<-dat$sgRNA

  dds <- DESeqDataSetFromMatrix(countData = df.deseq2, colData = col.data, design = ~ condition)
  dds <- DESeq(dds)
  res <- rownames_to_column(as.data.frame(results(dds)), "sgRNA")
  list("sgRNA"=res)
}

run.edgeR <- function(dat) {
  df.edgeR <-dat[,-(1:4)]
  nx <- ncol(df.edgeR)
  group <- c(rep("T0", nx/2),rep("T1", nx/2))
  res1<- DGEList(counts=df.edgeR, group=group)
  res2<- estimateCommonDisp(res1)
  res3<- estimateTagwiseDisp(res2)
  res4<- exactTest(res3)
  list("sgRNA"=cbind(data.frame(sgRNA=dat$sgRNA), data.frame(res4$table)))
}

run.ScreenBEAM <- function(dat) {
  tmp.name <- tempfile()
  save(dat, file=tmp.name)
  tmp.outname <- tempfile()
  cmd <- paste("Rscript", system.file("extdata", "ScreenBEAM.R", package="CC2Sim"), tmp.name, tmp.outname)
  system(cmd)
  load(tmp.outname)
  colnames(df.ret) <- c("gene", "nSGRNA", "nSGRNA.good", "beta", "zstat", "pval", "FDR", "bsta")
  list("gene"=df.ret)
}

run.sgRSEA <- function(dat) {
  nx <- ncol(dat)-4
  control.samples <- colnames(dat)[5:(5+nx/2-1)]
  case.samples <- colnames(dat)[(5+nx/2):(5+nx-1)]

  dat <- UQnormalize(dat, trt=case.samples, ctrl=control.samples)
  results <- sgRSEA(dat=dat, multiplier=30)

  pos <- results$gene.pos
  neg <- results$gene.neg

  df.gene <- dplyr::left_join(rownames_to_column(as.data.frame(pos), "gene"),
                              rownames_to_column(as.data.frame(neg), "gene"), by="gene") %>%
    dplyr::mutate(p.value.twosides=pmin(1,pmin(p.value.pos , p.value.neg )*2))


  ret <- list("gene"=df.gene)
  ret
}

run.PBNPA <- function(dat) {
  dat$gene <- as.character(dat$gene)
  dat$sgRNA <- as.character(dat$sgRNA)
  nx <- (ncol(dat)-4)/2

  datlist <- list()
  for(i in 1:nx) {
    datlist[[i]] <- data.frame(sgRNA = dat$sgRNA,
                               Gene = dat$gene,
                               initial.count = dat[,i+4],
                               final.count = dat[,i+4+nx])
  }
  result <- PBNPA(datlist)$final.result %>%
    dplyr::mutate(p.value.twosides=pmin(1,pmin(pos.pvalue , neg.pvalue )*2))

  ret <- list("gene"=result)
}

run.RSA <- function(dat) {
  nx <- ncol(dat)-4
  ctrl.median <- dat[,5:(5+nx/2-1)] %>%
    apply(1, median, na.rm = TRUE)

  test.median <- dat[,(5+nx/2):(5+nx-1)] %>%
    apply(1, median, na.rm = TRUE)

  FC <- test.median / ctrl.median
  df.RSA <- data.frame(Gene_ID=dat$gene,
                       Well_ID=dat$sgRNA,
                       Score=FC)

  ret.RSA <- RSA(df.RSA, LB=0, UB=1e8)
  ret.gene <- ret.RSA %>% dplyr::group_by(Gene_ID) %>%
    dplyr::summarise(score = mean(LogP))
  ret <- list("gene"=ret.gene)
}

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
    x[is.na(x$score),"score"] <- -1e8
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

run <- function(dat, methods, selector, cache.dir = NULL) {
  sim.dat <- dat
  results.sgRNA <- NULL
  results.gene <- NULL
  selector <- as.data.frame(selector)
  for (i in names(methods)) {
    cat("Running", i, "...", "\n")

    cache.sgRNA <- NULL
    cache.gene <- NULL

    if(!is.null(cache.dir)) {
      cache.sgRNA <- file.path(cache.dir, paste0(i, "_", "sgRNA.csv"))
      cache.gene <- file.path(cache.dir, paste0(i, "_", "gene.csv"))
    }

    df.ret <- list()

    if((is.null(cache.sgRNA)&&is.null(cache.gene)) ||
       (!file.exists(cache.sgRNA) && !file.exists(cache.gene)) ) {
      df.ret <- methods[[i]](dat)
      if(!is.null(df.ret$sgRNA)) {
        write_csv(df.ret$sgRNA, cache.sgRNA)
      }
      if(!is.null(df.ret$gene)) {
        write_csv(df.ret$gene, cache.gene)
      }
    } else {
      if(file.exists(cache.sgRNA)) {
        df.ret$sgRNA <- read_csv(cache.sgRNA)
      }
      if(file.exists(cache.gene)) {
        df.ret$gene <- read_csv(cache.gene)
      }
    }

    gene_id <- selector[selector$name==i,]$gene_id
    sgrna_id <- selector[selector$name==i,]$sgrna_id
    tar.sgrna.column <- selector[selector$name==i, "sgRNA.column"]
    tar.sgrna.func <- selector[selector$name==i, "sgRNA.func"]
    tar.gene.column <- selector[selector$name==i, "gene.column"]
    tar.gene.func <- selector[selector$name==i, "gene.func"]

    if(!is.na(sgrna_id)) {
      tmp.sgRNA <- data.frame("sgRNA"=df.ret$sgRNA[,sgrna_id],
                              "score"=sapply(df.ret$sgRNA[,tar.sgrna.column],
                                           get(tar.sgrna.func)))
      colnames(tmp.sgRNA) <- c("sgRNA", "score")
      if (is.null(results.sgRNA)) {
        results.sgRNA <- tmp.sgRNA
      } else {
        results.sgRNA <- dplyr::left_join(results.sgRNA, tmp.sgRNA, by = "sgRNA")
      }
      nc <- ncol(results.sgRNA)
      colnames(results.sgRNA)[nc] <- i
    }
    if (!is.na(gene_id)) {
      tmp.gene <- data.frame("gene"=df.ret$gene[,gene_id],
                             "score"=sapply(df.ret$gene[,tar.gene.column],
                                          get(tar.gene.func)))
      colnames(tmp.gene) <- c("gene", "score")

      if (is.null(results.gene)) {
        results.gene <- tmp.gene
      }
      else {
        results.gene <- dplyr::left_join(results.gene, tmp.gene, by = "gene")
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

dat <- load.sample()

methods = list(
  MAGeCK = run.mageck,
  DESeq2 = run.DESeq2,
  edgeR = run.edgeR,
  sgRSEA = run.sgRSEA,
  PBNPA = run.PBNPA,
  ScreenBEAM = run.ScreenBEAM,
  CC2 = run.mbttest
)

library(tidyverse)
library(cowplot)
library(edgeR)
library(DESeq2)
library(sgRSEA)
library(PBNPA)
library(CRISPRCloud2)
library(plotROC)
library(PRROC)
selector <- read_delim("inst/extdata/negative-expr.tsv", delim="\t")
load("inst/extdata/nature-biotech.Rdata")

cache.dir <- "/Users/hwan/Sandbox/CC2Sim/cache/RT112"
dir.create(cache.dir, showWarnings = T, recursive = T, mode = "0777")
ret <- run(dataset$RT112, methods, selector, cache.dir)
rtt.x <- plot.all(ret, "ROC", "RTT12 dataset benchmark (AUCROC)")
rtt.y <- plot.all(ret, "PR", "RTT12 dataset benchmark (AUCPRC)")

cache.dir <- "/Users/hwan/Sandbox/CC2Sim/cache/UMUC3"
dir.create(cache.dir, showWarnings = T, recursive = T, mode = "0777")
ret2 <- run(dataset$UMUC3, methods, selector, cache.dir)
umuc.x  <- plot.all(ret2, "ROC", "UMUC3 dataset benchmark (AUCROC)")
umuc.y  <- plot.all(ret2, "PR", "UMUC3 dataset benchmark (AUCPRC)")
both <- plot_grid(rtt.x, umuc.x, rtt.y, umuc.y)

rank.heatmap <- function(gene.ret, filename) {
  gene.ret <- ret$gene
  gene.rank <- gene.ret[,-c(1,ncol(gene.ret))]
  for(i in 1:(ncol(gene.ret)-2)) {
    gene.rank[,i] <- rank(gene.ret[,i+1])
  }
  row.names(gene.rank) <- gene.ret$gene
  library(pheatmap)
  library(RColorBrewer)

  df.anno <- data.frame(label=ifelse(gene.ret$label,"essential", "non essential"))
  row.names(df.anno) <- gene.ret$gene
  order(rowSums(gene.rank))
  pheatmap(t(gene.rank[order(gene.rank$CC2),]), cluster_cols = F,
           color = colorRampPalette((brewer.pal(n = 11, name ="RdYlBu")))(100),
           annotation_col = df.anno)
}
rank.heatmap(ret$gene, "RT112-gene-heatmap.pdf")
rank.heatmap(UMUC3.ret$gene, "UMUC3-gene-heatmap.pdf")
