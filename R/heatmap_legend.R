generate_heatmap_legend <- function() {
  # df <- data.frame(
  #   x = 1:5,
  #   y = 1:5,
  #   FDR = c("≥8", "≥6", "≥4", "≥2", "≥0"),
  #   w = 1:5,
  #   stringsAsFactors = F
  # )

  df <- data.frame(
    x = 1:5,
    y = 1:5,
    FDR = c("\u2265 8", "\u2265 6", "\u2265 4", "\u2265 2", "\u2265 0"),
    #FDR = sprintf("\u2265%d", seq(8,0)),
    w = 1:5,
    stringsAsFactors = F
  )

  (col.pal <- RColorBrewer::brewer.pal(9, "Reds"))
  col.pal[1] <- "#FFFFFF"

  test <- ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = FDR), colour = "grey50") +
    scale_fill_manual(
      name = "-log10(FDR)",
      #values = rev(col.pal)[c(1,3,5,7,9)],
      values = rev(col.pal)[c(1,3,5,7,9)],
      limits = df$FDR
    )

  legend_heatmap_1 <- get_legend(test)


  df <- data.frame(
    x = 1:2,
    y = 1:2,
    essential = c("Essential", "Non-essential"),
    w = 1:2,
    stringsAsFactors = F
  )

  test <- ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = essential), colour = "grey50") +
    scale_fill_manual(
      name = "Essentiality",
      values = c("#000000", "#ffffff"),
      limits = df$essential
    )

  legend_heatmap_2 <- get_legend(test)

  plot_grid(
    NULL,
    legend_heatmap_1,
    legend_heatmap_2,
    NULL,
    ncol = 1,
    rel_heights = c(0.5, 0.8, 0.8, 0.5)
  )
}


generate_heatmap_legend_no_essential <- function() {
  df <- data.frame(
    x = 1:5,
    y = 1:1,
    FDR = c("8", "6", "4", "2", "0"),
    #FDR = sprintf("\u2265%d", seq(8,0)),
    w = 1:5,
    stringsAsFactors = F
  )
  print(df)
  (col.pal <- RColorBrewer::brewer.pal(9, "Reds"))
  col.pal[1] <- "#FFFFFF"

  test <- ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = FDR), colour = "grey50") +
    scale_fill_manual(
      name = "-log10(FDR)",
      #values = rev(col.pal)[c(1,3,5,7,9)],
      values = rev(col.pal)[c(1,3,5,7,9)],
      limits = df$FDR
    ) + geom_text(aes(label=FDR)) +
    theme(legend.position = "none") + ylab("-log10(FDR)") +
    theme(axis.title.y = element_text(angle=0)) +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank()) + xlab("")
  test
}

generate_heatmap_legend_no_essential()
