suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(patchwork)
  library(ggrepel)
  library(ggpubr)
  library(Nebulosa)
  library(grid)
})
if (requireNamespace("extrafont", quietly = TRUE)) extrafont::loadfonts(quiet = TRUE)
# F6G
plot_F6G <- function(tpm, ann_col, ann_colors, output_dir = ".") {
  p = pheatmap::pheatmap(tpm, show_rownames = T, cluster_cols = F, border = FALSE, scale = "row",
    angle_col = 45, annotation_col = ann_col, annotation_colors = ann_colors, gaps_col = 3,
    color = colorRampPalette(c("#516da0", "#fcf8f8", "#a94c4a"))(100), cellwidth = 25,
    cellheight = 17, fontsize = 13, silent = TRUE)
  ggsave(file.path(output_dir, "F6G.pdf"), plot = p$gtable, width = 6, height = 5)
}
# F6I
plot_F6I <- function(umap, stype_color, output_dir = ".") {
  dens_layers <- umap %>%
    split(.$subtype) %>%
    imap(~ {
      bins_g <- round(scales::rescale(log10(nrow(.x)), to = c(2, 5)))
      h_umap <- c(MASS::bandwidth.nrd(.x$UMAP_1), MASS::bandwidth.nrd(.x$UMAP_2)) * 1.18
      stat_density_2d(data = .x, aes(x = UMAP_1, y = UMAP_2, fill = subtype, color = subtype),
        geom = "polygon", alpha = 0.2, bins = bins_g, size = 0.7, h = h_umap)
    })
  p <- ggplot(umap, aes(UMAP_1, UMAP_2)) +
    ggrastr::geom_point_rast(color = "#d4d4d4", size = 0.3,
      alpha = 0.1) + dens_layers + scale_color_manual(values = stype_color, na.value = NA) +
    scale_fill_manual(values = stype_color, na.value = NA) +
    theme_void(base_size = 14) +
    theme(legend.title = element_blank(), plot.margin = margin(r = 10))
  ggsave(file.path(output_dir, "F6I.pdf"), plot = p, width = 4.2, height = 2.7)
}
# F6J
plot_F6J <- function(umap, output_dir = ".") {
  p1 <- ggplot(umap, aes(UMAP_1, UMAP_2)) +
    ggrastr::geom_point_rast(color = "#d4d4d4", size = 0.3,
      alpha = 0.1) +
    stat_density2d(aes(fill = treat, color = treat), size = 0.7, geom = "polygon",
      alpha = 0.3, bins = 4) +
    scale_color_manual(values = c(CTRL = "#548adb", CKO = "#cf5b53"),
      na.value = NA) +
    scale_fill_manual(values = c(CTRL = "#548adb", CKO = "#cf5b53"), na.value = NA) +
    scale_x_continuous(limits = c(-4, 21)) +
    theme_void(base_size = 14) +
    theme(legend.title = element_blank(),
      plot.margin = margin(r = 10))
  ggsave(file.path(output_dir, "F6J.pdf"), plot = p1, width = 3.9, height = 2.7)
}
# F6K
plot_F6K <- function(sce.T, output_dir = ".") {
  p = plot_density(sce.T, reduction = "umap", features = c("Pdcd1", "Gzmk", "Foxp3", "Cxcr3"),
    size = 0.1, pal = "inferno")
  p = lapply(p, function(x) x + theme_void() +
    theme(plot.title = element_text(hjust = 0.57,
      size = 6.5, margin = margin(b = 0.5)), plot.margin = margin(b = 6, r = 10), text = element_text(size = 6,
      family = "Arial"), legend.key.size = unit(0.2, "cm"), legend.margin = margin(l = 0.5),
    legend.text = element_text(margin = margin(l = 2))))
  p = patchwork::wrap_plots(p)
  ggsave(file.path(output_dir, "F6K.pdf"), plot = p, width = 2.9, height = 2)
}
# F6L
plot_F6L <- function(de_markers, anti_tumor_Tmarkers, pro_tumor_Tmarkers, output_dir = ".") {
  library(ggrepel)
  de_markers <- de_markers %>%
    mutate(group = case_when(gene %in% anti_tumor_Tmarkers ~ "Anti-tumor", gene %in% pro_tumor_Tmarkers ~
      "Pro-tumor", TRUE ~ "Other"))
  fs = 6
  p = ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj + 9.99988867182683e-321))) +
    ggrastr::geom_point_rast(color = "grey", alpha = 0.3, size = 0.3) +
    geom_point(data = de_markers %>%
      filter(group != "Other" & -log10(p_val_adj + 9.99988867182683e-321) > 2), aes(color = group),
    size = 0.4) +
    scale_color_manual(values = c(`Anti-tumor` = "tomato", `Pro-tumor` = "steelblue",
      Other = "grey")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey",
      size = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey",
      size = 0.2) +
    geom_text_repel(seed = 1, fontface = "italic", size = fs / .pt, family = "Arial",
      data = de_markers %>%
        filter(group != "Other"), aes(label = gene), segment.size = 0.15, segment.alpha = 0.4,
      max.overlaps = 20, min.segment.length = 0.3, force_pull = 10, force = 20, point.padding = 0.1) +
    theme_bw() +
    theme(legend.title = element_blank(), text = element_text(family = "Arial",
      size = fs), panel.grid = element_blank(), legend.position = "none", axis.ticks = element_line(size = 0.2),
    axis.text = element_text(color = "black")) +
    labs(x = expression(paste("Average ",
      log[2], "FC")), y = expression(paste(-log[10], " adjusted P-value")))
  ggsave(file.path(output_dir, "F6L.pdf"), plot = p, width = 2.4, height = 2)
}
# S6G
plot_S6G <- function(sce.T, subtype_markers, output_dir = ".") {
  p = DotPlot(sce.T, subtype_markers, group.by = "subtype", dot.scale = 2.3) +
    RotatedAxis() +
    theme(text = element_text(size = 6, family = "Arial"), axis.text = element_text(size = 6,
      family = "Arial"), axis.ticks = element_line(size = 0.2), legend.key.size = unit(0.25,
      "cm"), panel.border = element_rect(color = "black", size = 0.25), panel.grid.minor = element_line(color = "gray90",
      size = 0.25), panel.spacing = unit(1, "mm"), strip.text = element_text(margin = margin(b = 1,
      unit = "mm"), size = 6), strip.placement = "outlet", axis.line = element_blank()) +
    labs(x = "", y = "") +
    scale_color_gradient2(low = "#1f273a", mid = "#fdf8f8", high = "#892424",
      midpoint = 0)
  ggsave(file.path(output_dir, "S6G.pdf"), plot = p, width = 6.4, height = 2.1)
}
# S6H
plot_S6H <- function(sce.T, output_dir = ".") {
  cols11 <- c("#1b9e77", "#ee781e", "#8984c7", "#c56f9b", "#66a61e", "#eebf3c", "#ca8b17",
    "#a89191", "#4c93c2", "#b2df8a", "#33a02c")
  p = DimPlot(sce.T, group.by = "subtype", cols = cols11, pt.size = 6, raster = T, raster.dpi = c(1024,
    1024)) +
    labs(title = "") +
    theme_void() +
    theme(text = element_text(size = 13, family = "Arial"),
      legend.key.size = unit(0.5, "cm"))
  ggsave(file.path(output_dir, "S6H.pdf"), plot = p, width = 4.2, height = 3)
}
# S6I
plot_S6I <- function(da_results_se, output_dir = ".") {
  p = ggplot(da_results_se, aes(y = reorder(subtype, log2FC), x = log2FC, color = color)) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin = log2FC - se_log2FC, xmax = log2FC +
      se_log2FC), width = 0.2) +
    geom_vline(xintercept = 0, linetype = "solid") +
    geom_vline(xintercept = c(-1,
      1), linetype = "dashed") +
    theme_classic(base_size = 13) +
    xlab("log2FC of mean proportion (CKO vs CTRL)") +
    ylab(NULL) +
    theme(axis.text = element_text(color = "black"), panel.border = element_rect(fill = NA,
      color = "black", linewidth = 1), axis.line = element_blank(), axis.text.y = element_text(size = 11),
    ) +
    scale_color_identity()
  ggsave(file.path(output_dir, "S6I.pdf"), plot = p, width = 6, height = 4)
}
# S6J
plot_S6J <- function(sce.T, output_dir = ".") {
  p = plot_density(sce.T, reduction = "umap", features = c("Cd8a", "Gzma", "Ccl5", "Lag3",
    "Havcr2", "Cd4", "Tigit", "Isg15", "Tcf7", "Trdc"), size = 0.1, pal = "inferno")
  p = lapply(p, function(x) x + theme_void() +
    theme(plot.title = element_text(hjust = 0.56,
      size = 6.5, margin = margin(b = 0.5)), plot.margin = margin(b = 6, r = 10), text = element_text(size = 6,
      family = "Arial"), legend.key.size = unit(0.2, "cm"), legend.margin = margin(l = 0.5),
    legend.text = element_text(margin = margin(l = 2))))
  p = patchwork::wrap_plots(p, ncol = 5)
  ggsave(file.path(output_dir, "S6J.pdf"), plot = p, width = 6.5, height = 2)
}
# S6K
plot_S6K <- function(sce.T, output_dir = ".") {
  plots_violins = VlnPlot(sce.T, group.by = "treat", features = c("Pdcd1", "Tox", "Il2ra",
    "Ikzf2", "Gzmb", "Ccl5"), pt.size = 0, ncol = 3, combine = FALSE, cols = c("#8bbae4",
    "#e29a92"))
  for (i in 1:length(plots_violins)) {
    plots_violins[[i]] <- plots_violins[[i]] + theme(panel.border = element_rect(fill = NA,
      color = "black", linewidth = 1), plot.title = element_text(face = "plain", size = 15),
    axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text = element_text(color = "black")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      ggpubr::stat_compare_means(label.x.npc = 0.4,
        vjust = 1, label = "p.signif") +
      labs(x = NULL) +
      NoLegend()
  }
  p = patchwork::wrap_plots(plots_violins)
  ggsave(file.path(output_dir, "S6K.pdf"), plot = p, width = 6, height = 4)
}
# S6L
plot_S6L <- function(gse_res2, output_dir = ".") {
  p <- gse_res2 %>%
    mutate(logFDR = -log10(p.adjust), logq = -log10(qvalues)) %>%
    ggplot(aes(x = NES, y = reorder(Description, NES), size = logq, color = logFDR)) +
    geom_point() +
    scale_color_gradient(low = "#d8bfab", high = "#d06a1b") +
    scale_size(range = c(3,
      5.5)) +
    scale_x_continuous(limits = c(2.45, 3.15)) +
    labs(x = "NES", y = "", color = "-log10(FDR)",
      size = "-log10(qvalue)") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 11), plot.margin = margin(t = 10))
  ggsave(file.path(output_dir, "S6L.pdf"), plot = p, width = 5.5, height = 3.5)
}
# S6M
plot_S6M <- function(df1, df2, df3, df4, output_dir = ".") {
  plot_go = function(df, gocolor, gotitle) {
    df$LogP = -log10(df$p.adjust)
    df$labelx = rep(0.02 * max(df$LogP), nrow(df))
    df$labely = seq(nrow(df), 1) - 0.2
    df$Description = paste0(stringr::str_wrap(df$Description, width = 35), "\n")
    ggplot(data = df, aes(LogP, reorder(Description, LogP))) +
      geom_bar(stat = "identity",
        alpha = 0.3, fill = gocolor, width = 0.8) +
      geom_text(aes(x = labelx, y = labely,
        label = Description), size = 6 / .pt, family = "Arial", hjust = 0) +
      theme_classic() +
      theme(text = element_text(size = 6, family = "Arial"), plot.title = element_text(size = 6.5,
        family = "Arial", margin = margin(b = 0.5)), axis.text.y = element_blank(),
      axis.line.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line.x = element_line(colour = "black", size = 0.25), axis.text.x = element_text(colour = "black"),
      axis.ticks.x = element_line(colour = "black", size = 0.2), axis.ticks.length.x = unit(0.6,
        "mm")) +
      xlab("-log10(adj.p)") + ggtitle(gotitle) +
      scale_x_continuous(expand = c(0,
        0))
  }
  plot_go(df1, "#63aed6", "Cd8_Tex") -> p1
  plot_go(df2, "#f66f5d", "Cd8_Tem") -> p2
  plot_go(df3, "#496297", "Cd4_Treg") -> p3
  plot_go(df4, "#a55361", "Cd4_Tm") -> p4
  p = p1 + p2 + p3 + p4 + plot_layout(ncol = 4)
  ggsave(file.path(output_dir, "S6M.pdf"), plot = p, width = 7, height = 2)
}
