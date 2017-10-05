Report real data analysis
================

Introduction
------------

This is the R markdown report of the analysis for the CRISPR dataset of [Evers et al., Nat. Bio. Tech. 2016](http://www.nature.com/nbt/journal/v34/n6/full/nbt.3536.html)

Download dataset & Data Preprocessing
-------------------------------------

``` r
library(CRISPRCloud2)
library(CC2Sim)
library(readxl)
library(tidyverse)
library(knitr)
library(plotROC)
library(cowplot)
library(PRROC)
system("wget http://www.nature.com/nbt/journal/v34/n6/extref/nbt.3536-S3.xlsx")
nat <- readxl::read_excel("nbt.3536-S3.xlsx")
```

See I selected proper positions in the sheet correctly.

``` r
RT112 <- select(nat, 2, 3, 1, 4, 5:7, 8:10)  
RT112
```

    ## # A tibble: 961 x 10
    ##                Sequence   Gene LibraryID Essential `t=0 Rep. 1\r\nRT112`
    ##                   <chr>  <chr>     <chr>     <dbl>                 <dbl>
    ##  1 TTGACGGTATCCACCCATTC  RPS19  RPS19-10         1                 10180
    ##  2 CAGTACACGGCCCGCTTGTA  NUP93   NUP93-4         1                  9073
    ##  3 AGAAGGGTCGTACATACTTC PSMD11  PSMD11-2         1                  9408
    ##  4 CTTACTGGTTCCTTGGGTCC  RPS3A   RPS3A-3         1                  8922
    ##  5 ATGGGTGACCGGCTGTACAT  PSMB3   PSMB3-3         1                  7434
    ##  6 GAGGGCAAGCGGATTCCATT POLR2A  POLR2A-5         1                  7672
    ##  7 GCTCGCACCCACCTTGGTGT   RPS8    RPS8-8         1                  6999
    ##  8 CATGGCGCGCCGCTCGTACG  RPL36   RPL36-3         1                  5988
    ##  9 CCAGCGCGCTACTTACAGTG  RPS13   RPS13-4         1                  5763
    ## 10 AGATGGCGGACATTCAGGTG  RPS11   RPS11-2         1                  9444
    ## # ... with 951 more rows, and 5 more variables: `t=0 Rep.
    ## #   2\r\nRT112` <dbl>, `t=0 Rep. 3\r\nRT112` <dbl>, `t=1 Rep.
    ## #   1\r\nRT112` <dbl>, `t=1 Rep. 2\r\nRT112` <dbl>, `t=1 Rep.
    ## #   3\r\nRT112` <dbl>

For UMUC3, too.

``` r
UMUC3 <- select(nat, 2, 3, 1, 4, 17:22)
UMUC3
```

    ## # A tibble: 961 x 10
    ##                Sequence   Gene LibraryID Essential `t=0 Rep. 1\r\nUMUC3`
    ##                   <chr>  <chr>     <chr>     <dbl>                 <dbl>
    ##  1 TTGACGGTATCCACCCATTC  RPS19  RPS19-10         1                 10693
    ##  2 CAGTACACGGCCCGCTTGTA  NUP93   NUP93-4         1                  9073
    ##  3 AGAAGGGTCGTACATACTTC PSMD11  PSMD11-2         1                  9695
    ##  4 CTTACTGGTTCCTTGGGTCC  RPS3A   RPS3A-3         1                  9600
    ##  5 ATGGGTGACCGGCTGTACAT  PSMB3   PSMB3-3         1                  7575
    ##  6 GAGGGCAAGCGGATTCCATT POLR2A  POLR2A-5         1                  8467
    ##  7 GCTCGCACCCACCTTGGTGT   RPS8    RPS8-8         1                  7857
    ##  8 CATGGCGCGCCGCTCGTACG  RPL36   RPL36-3         1                  6367
    ##  9 CCAGCGCGCTACTTACAGTG  RPS13   RPS13-4         1                  6284
    ## 10 AGATGGCGGACATTCAGGTG  RPS11   RPS11-2         1                 10211
    ## # ... with 951 more rows, and 5 more variables: `t=0 Rep.
    ## #   2\r\nUMUC3` <dbl>, `t=0 Rep. 3\r\nUMUC3` <dbl>, `t=1 Rep.
    ## #   1\r\nUMUC3` <dbl>, `t=1 Rep. 2\r\nUMUC3` <dbl>, `t=1 Rep.
    ## #   3\r\nUMUC3` <dbl>

Need to set proper names for each column, and change the labels.

``` r
dataset <- list("RT112" = RT112, "UMUC3" = UMUC3)
for(d in names(dataset)) {
  colnames(dataset[[d]]) <- c("X", "gene", "sgRNA", "class", "B1", "B2", "B3", "A1", "A2", "A3")
  dataset[[d]][dataset[[d]]==0] <- "inactive"
  dataset[[d]][dataset[[d]]==1] <- "increasing"

}
dataset
```

    ## $RT112
    ## # A tibble: 961 x 10
    ##                       X   gene    sgRNA      class    B1    B2    B3    A1
    ##                   <chr>  <chr>    <chr>      <chr> <dbl> <dbl> <dbl> <dbl>
    ##  1 TTGACGGTATCCACCCATTC  RPS19 RPS19-10 increasing 10180  9768  9406  1005
    ##  2 CAGTACACGGCCCGCTTGTA  NUP93  NUP93-4 increasing  9073  8598  8363   688
    ##  3 AGAAGGGTCGTACATACTTC PSMD11 PSMD11-2 increasing  9408  9573  9384  1014
    ##  4 CTTACTGGTTCCTTGGGTCC  RPS3A  RPS3A-3 increasing  8922  8965  8779  1266
    ##  5 ATGGGTGACCGGCTGTACAT  PSMB3  PSMB3-3 increasing  7434  6958  6871   450
    ##  6 GAGGGCAAGCGGATTCCATT POLR2A POLR2A-5 increasing  7672  7537  7259   647
    ##  7 GCTCGCACCCACCTTGGTGT   RPS8   RPS8-8 increasing  6999  6966  6775   673
    ##  8 CATGGCGCGCCGCTCGTACG  RPL36  RPL36-3 increasing  5988  5893  5498   266
    ##  9 CCAGCGCGCTACTTACAGTG  RPS13  RPS13-4 increasing  5763  5546  5131   187
    ## 10 AGATGGCGGACATTCAGGTG  RPS11  RPS11-2 increasing  9444  9621  8835  1900
    ## # ... with 951 more rows, and 2 more variables: A2 <dbl>, A3 <dbl>
    ## 
    ## $UMUC3
    ## # A tibble: 961 x 10
    ##                       X   gene    sgRNA      class    B1    B2    B3    A1
    ##                   <chr>  <chr>    <chr>      <chr> <dbl> <dbl> <dbl> <dbl>
    ##  1 TTGACGGTATCCACCCATTC  RPS19 RPS19-10 increasing 10693 11662  7010  1517
    ##  2 CAGTACACGGCCCGCTTGTA  NUP93  NUP93-4 increasing  9073  9115  6162   514
    ##  3 AGAAGGGTCGTACATACTTC PSMD11 PSMD11-2 increasing  9695 10711  6476   923
    ##  4 CTTACTGGTTCCTTGGGTCC  RPS3A  RPS3A-3 increasing  9600  9622  5971  3901
    ##  5 ATGGGTGACCGGCTGTACAT  PSMB3  PSMB3-3 increasing  7575  8571  5262   680
    ##  6 GAGGGCAAGCGGATTCCATT POLR2A POLR2A-5 increasing  8467  8619  5644  1285
    ##  7 GCTCGCACCCACCTTGGTGT   RPS8   RPS8-8 increasing  7857  7964  5314   706
    ##  8 CATGGCGCGCCGCTCGTACG  RPL36  RPL36-3 increasing  6367  6452  4452   271
    ##  9 CCAGCGCGCTACTTACAGTG  RPS13  RPS13-4 increasing  6284  6109  3943   297
    ## 10 AGATGGCGGACATTCAGGTG  RPS11  RPS11-2 increasing 10211  9993  6170  3438
    ## # ... with 951 more rows, and 2 more variables: A2 <dbl>, A3 <dbl>

Running test
------------

edit functions from `CC2Sim` package for this simulations

``` r
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
  
  curve <- curve %>% arrange(desc(y))
  
  ret <- list()
  ret$bar.plot <- ggplot(prv, aes(x=method, y=AUPRC)) +
    geom_bar(stat = "identity", aes(fill=method)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  ret$curve.plot <- ggplot(curve, aes(x=x, y=y)) +
    geom_line(aes(colour=method)) +
    ylab("Precision") + xlab("Recall") + ylim(0,1)
  
  ret
}
run <- function(dat) {
  methods <- list(
    mageck = run.mageck,
    DESeq2 = run.DESeq2,
    edgeR = run.edgeR,
    sgRSEA = run.sgRSEA,
    #screenBEAM = run.ScreenBEAM,
    CC2 = run.mbttest
  )
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
```

Run RT112 dataset

``` r
RT112.ret <- run(dataset$RT112)
```

    ## Running mageck ... 
    ## Running DESeq2 ...

    ## Warning: Setting row names on a tibble is deprecated.

    ## Running edgeR ... 
    ## Running sgRSEA ... 
    ## Running CC2 ...

    ## Warning: Column `sgRNA` joining factor and character vector, coercing into
    ## character vector

    ## Warning: Column `gene` joining factor and character vector, coercing into
    ## character vector

    ## # A tibble: 538 x 4
    ##    method         x         y density
    ##     <chr>     <dbl>     <dbl>   <dbl>
    ##  1 mageck 1.0000000 0.5057232 0.00140
    ##  2 mageck 1.0000000 0.5057232 0.00140
    ##  3 mageck 0.9979424 0.5052083 0.03206
    ##  4 mageck 0.9979424 0.5057351 0.03885
    ##  5 mageck 0.9979424 0.5062630 0.04370
    ##  6 mageck 0.9979424 0.5067921 0.09556
    ##  7 mageck 0.9979424 0.5073222 0.10023
    ##  8 mageck 0.9979424 0.5078534 0.10087
    ##  9 mageck 0.9958848 0.5073375 0.12518
    ## 10 mageck 0.9958848 0.5078699 0.12753
    ## # ... with 528 more rows
    ## # A tibble: 521 x 4
    ##    method         x         y    density
    ##     <chr>     <dbl>     <dbl>      <dbl>
    ##  1 DESeq2 1.0000000 0.5057232 0.01551059
    ##  2 DESeq2 1.0000000 0.5057232 0.01551059
    ##  3 DESeq2 0.9979424 0.5052083 0.02826186
    ##  4 DESeq2 0.9979424 0.5057351 0.04230840
    ##  5 DESeq2 0.9979424 0.5062630 0.07668917
    ##  6 DESeq2 0.9979424 0.5067921 0.07901892
    ##  7 DESeq2 0.9979424 0.5073222 0.09265575
    ##  8 DESeq2 0.9979424 0.5078534 0.10481658
    ##  9 DESeq2 0.9979424 0.5083857 0.12327689
    ## 10 DESeq2 0.9958848 0.5078699 0.12655657
    ## # ... with 511 more rows
    ## # A tibble: 373 x 4
    ##    method         x         y   density
    ##     <chr>     <dbl>     <dbl>     <dbl>
    ##  1  edgeR 1.0000000 0.5057232 0.1074991
    ##  2  edgeR 1.0000000 0.5057232 0.1074991
    ##  3  edgeR 0.9979424 0.5052083 0.1521552
    ##  4  edgeR 0.9958848 0.5046924 0.1765458
    ##  5  edgeR 0.9938272 0.5041754 0.1879511
    ##  6  edgeR 0.9917695 0.5036573 0.2289082
    ##  7  edgeR 0.9897119 0.5031381 0.2338115
    ##  8  edgeR 0.9876543 0.5026178 0.3067432
    ##  9  edgeR 0.9855967 0.5020964 0.4026057
    ## 10  edgeR 0.9835391 0.5015740 0.4243014
    ## # ... with 363 more rows
    ## # A tibble: 495 x 4
    ##    method     x         y     density
    ##     <chr> <dbl>     <dbl>       <dbl>
    ##  1    CC2     1 0.5057232 0.000369413
    ##  2    CC2     1 0.5057232 0.000369413
    ##  3    CC2     1 0.5062500 0.002493961
    ##  4    CC2     1 0.5067779 0.004607585
    ##  5    CC2     1 0.5073069 0.006553482
    ##  6    CC2     1 0.5078370 0.008879870
    ##  7    CC2     1 0.5083682 0.017415824
    ##  8    CC2     1 0.5089005 0.018313099
    ##  9    CC2     1 0.5094340 0.024979438
    ## 10    CC2     1 0.5099685 0.025312515
    ## # ... with 485 more rows
    ## # A tibble: 88 x 4
    ##    method         x         y density
    ##     <chr>     <dbl>     <dbl>   <dbl>
    ##  1 mageck 1.0000000 0.4946237 0.00000
    ##  2 mageck 1.0000000 0.4946237 0.00000
    ##  3 mageck 0.9782609 0.5113636 0.15784
    ##  4 mageck 0.9565217 0.5057471 0.16420
    ##  5 mageck 0.9565217 0.5116279 0.28960
    ##  6 mageck 0.9347826 0.5058824 0.32120
    ##  7 mageck 0.9347826 0.5119048 0.35490
    ##  8 mageck 0.9130435 0.5060241 0.39860
    ##  9 mageck 0.8913043 0.5000000 0.49326
    ## 10 mageck 0.8913043 0.5061728 0.57134
    ## # ... with 78 more rows
    ## # A tibble: 51 x 4
    ##    method         x         y   density
    ##     <chr>     <dbl>     <dbl>     <dbl>
    ##  1 sgRSEA 1.0000000 0.4946237 0.3715514
    ##  2 sgRSEA 1.0000000 0.4946237 0.3715514
    ##  3 sgRSEA 0.9782609 0.4891304 0.4016482
    ##  4 sgRSEA 0.9565217 0.4835165 0.5012540
    ##  5 sgRSEA 0.9347826 0.4777778 0.5449660
    ##  6 sgRSEA 0.9130435 0.4719101 0.5757793
    ##  7 sgRSEA 0.8913043 0.4659091 0.5786456
    ##  8 sgRSEA 0.8913043 0.4712644 0.6495880
    ##  9 sgRSEA 0.8695652 0.4651163 0.6911501
    ## 10 sgRSEA 0.8695652 0.4705882 0.7155142
    ## # ... with 41 more rows
    ## # A tibble: 51 x 4
    ##    method     x         y     density
    ##     <chr> <dbl>     <dbl>       <dbl>
    ##  1    CC2     1 0.4946237 0.006611907
    ##  2    CC2     1 0.4946237 0.006611907
    ##  3    CC2     1 0.5000000 0.014266939
    ##  4    CC2     1 0.5054945 0.023658520
    ##  5    CC2     1 0.5111111 0.079843203
    ##  6    CC2     1 0.5168539 0.107330912
    ##  7    CC2     1 0.5227273 0.127124635
    ##  8    CC2     1 0.5287356 0.127337937
    ##  9    CC2     1 0.5348837 0.129916705
    ## 10    CC2     1 0.5411765 0.160793205
    ## # ... with 41 more rows

``` r
plot_grid(RT112.ret$plot.roc.sgRNA, RT112.ret$plot.roc.gene, RT112.ret$plot.prc.sgRNA, RT112.ret$plot.prc.gene)
```

![](real_data_analysis_report_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-7-1.png)

Run UMUC3 dataset

``` r
UMUC3.ret <- run(dataset$UMUC3)
```

    ## Running mageck ... 
    ## Running DESeq2 ...

    ## Warning: Setting row names on a tibble is deprecated.

    ## Running edgeR ... 
    ## Running sgRSEA ... 
    ## Running CC2 ...

    ## Warning: Column `sgRNA` joining factor and character vector, coercing into
    ## character vector

    ## Warning: Column `gene` joining factor and character vector, coercing into
    ## character vector

    ## # A tibble: 507 x 4
    ##    method         x         y density
    ##     <chr>     <dbl>     <dbl>   <dbl>
    ##  1 mageck 1.0000000 0.5057232 0.00001
    ##  2 mageck 1.0000000 0.5057232 0.00001
    ##  3 mageck 0.9979424 0.5052083 0.00638
    ##  4 mageck 0.9958848 0.5046924 0.01025
    ##  5 mageck 0.9958848 0.5052192 0.01291
    ##  6 mageck 0.9958848 0.5057471 0.02335
    ##  7 mageck 0.9958848 0.5062762 0.03620
    ##  8 mageck 0.9938272 0.5057592 0.06231
    ##  9 mageck 0.9938272 0.5062893 0.07818
    ## 10 mageck 0.9917695 0.5057712 0.09461
    ## # ... with 497 more rows
    ## # A tibble: 636 x 4
    ##    method         x         y     density
    ##     <chr>     <dbl>     <dbl>       <dbl>
    ##  1 DESeq2 1.0000000 0.5057232 0.007332787
    ##  2 DESeq2 1.0000000 0.5057232 0.007332787
    ##  3 DESeq2 1.0000000 0.5062500 0.024709861
    ##  4 DESeq2 0.9979424 0.5057351 0.032101098
    ##  5 DESeq2 0.9958848 0.5052192 0.034922649
    ##  6 DESeq2 0.9938272 0.5047022 0.043008426
    ##  7 DESeq2 0.9938272 0.5052301 0.052273584
    ##  8 DESeq2 0.9938272 0.5057592 0.053396558
    ##  9 DESeq2 0.9938272 0.5062893 0.077409702
    ## 10 DESeq2 0.9938272 0.5068206 0.081724017
    ## # ... with 626 more rows
    ## # A tibble: 639 x 4
    ##    method         x         y    density
    ##     <chr>     <dbl>     <dbl>      <dbl>
    ##  1  edgeR 1.0000000 0.5057232 0.01661292
    ##  2  edgeR 1.0000000 0.5057232 0.01661292
    ##  3  edgeR 1.0000000 0.5062500 0.04380816
    ##  4  edgeR 1.0000000 0.5067779 0.04706000
    ##  5  edgeR 0.9979424 0.5062630 0.06283521
    ##  6  edgeR 0.9958848 0.5057471 0.06304816
    ##  7  edgeR 0.9938272 0.5052301 0.06343362
    ##  8  edgeR 0.9938272 0.5057592 0.06623456
    ##  9  edgeR 0.9938272 0.5062893 0.07999450
    ## 10  edgeR 0.9938272 0.5068206 0.10227404
    ## # ... with 629 more rows
    ## # A tibble: 860 x 4
    ##    method         x         y      density
    ##     <chr>     <dbl>     <dbl>        <dbl>
    ##  1    CC2 1.0000000 0.5057232 0.0004109738
    ##  2    CC2 1.0000000 0.5057232 0.0004109738
    ##  3    CC2 0.9979424 0.5052083 0.0005750083
    ##  4    CC2 0.9979424 0.5057351 0.0007038608
    ##  5    CC2 0.9958848 0.5052192 0.0099259465
    ##  6    CC2 0.9938272 0.5047022 0.0134857245
    ##  7    CC2 0.9938272 0.5052301 0.0303216326
    ##  8    CC2 0.9917695 0.5047120 0.0408500330
    ##  9    CC2 0.9897119 0.5041929 0.0410724366
    ## 10    CC2 0.9876543 0.5036726 0.0435300053
    ## # ... with 850 more rows
    ## # A tibble: 84 x 4
    ##    method         x         y density
    ##     <chr>     <dbl>     <dbl>   <dbl>
    ##  1 mageck 1.0000000 0.4946237 0.00000
    ##  2 mageck 1.0000000 0.4946237 0.00000
    ##  3 mageck 1.0000000 0.4946237 0.00000
    ##  4 mageck 0.9782609 0.5037313 0.00000
    ##  5 mageck 0.9565217 0.5136187 0.00000
    ##  6 mageck 0.9347826 0.5243902 0.00000
    ##  7 mageck 0.9347826 0.5243902 0.19422
    ##  8 mageck 0.9347826 0.5308642 0.34802
    ##  9 mageck 0.9347826 0.5375000 0.40960
    ## 10 mageck 0.9347826 0.5443038 0.41106
    ## # ... with 74 more rows
    ## # A tibble: 65 x 4
    ##    method         x         y   density
    ##     <chr>     <dbl>     <dbl>     <dbl>
    ##  1 sgRSEA 1.0000000 0.4946237 0.1078466
    ##  2 sgRSEA 1.0000000 0.4946237 0.1078466
    ##  3 sgRSEA 0.9782609 0.4891304 0.1744894
    ##  4 sgRSEA 0.9565217 0.4835165 0.3880330
    ##  5 sgRSEA 0.9347826 0.4777778 0.4353278
    ##  6 sgRSEA 0.9130435 0.4719101 0.4869223
    ##  7 sgRSEA 0.8913043 0.4659091 0.5463991
    ##  8 sgRSEA 0.8695652 0.4597701 0.5471157
    ##  9 sgRSEA 0.8695652 0.4651163 0.5521319
    ## 10 sgRSEA 0.8695652 0.4705882 0.5643139
    ## # ... with 55 more rows
    ## # A tibble: 88 x 4
    ##    method         x         y    density
    ##     <chr>     <dbl>     <dbl>      <dbl>
    ##  1    CC2 1.0000000 0.4946237 0.06986857
    ##  2    CC2 1.0000000 0.4946237 0.06986857
    ##  3    CC2 0.9782609 0.4891304 0.20890706
    ##  4    CC2 0.9565217 0.4835165 0.29009046
    ##  5    CC2 0.9565217 0.4888889 0.29778659
    ##  6    CC2 0.9347826 0.4831461 0.31910239
    ##  7    CC2 0.9130435 0.4772727 0.34381437
    ##  8    CC2 0.9130435 0.4827586 0.37393165
    ##  9    CC2 0.9130435 0.4883721 0.39800792
    ## 10    CC2 0.8913043 0.4823529 0.56428488
    ## # ... with 78 more rows

``` r
plot_grid(UMUC3.ret$plot.roc.sgRNA, UMUC3.ret$plot.roc.gene, UMUC3.ret$plot.prc.sgRNA, UMUC3.ret$plot.prc.gene)
```

![](real_data_analysis_report_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-9-1.png)
