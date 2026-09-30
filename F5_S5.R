suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(patchwork)
  library(ggrepel)
  library(ggpubr)
  library(ggalluvial)
  library(Nebulosa)
  library(grid)
})
if (requireNamespace("extrafont", quietly = TRUE)) extrafont::loadfonts(quiet = TRUE)

# F5G
plot_F5G <- function(resM_CTSL, genelists, output_dir = ".") {
  Volcano <- function(deg, title, l1 = "up", l2 = "down", label = NULL) {
    deg <- deg[, c("log2FoldChange", "padj", "sig", "gene")] %>%
      na.omit()
    FC_t <- 1
    p_adj_t <- 0.05
    k1 = (deg$padj < p_adj_t) & (deg$log2FoldChange < -FC_t)
    k2 = (deg$padj < p_adj_t) & (deg$log2FoldChange > FC_t)
    if (!"sig" %in% colnames(deg)) {
      sig = ifelse(k1, l2, ifelse(k2, l1, "non-sig"))
      deg$sig <- factor(sig, levels = c(l2, "non-sig", l1))
    }
    if (is.null(label)) {
      label <- deg %>%
        filter(k1) %>%
        arrange(log2FoldChange) %>%
        head(15) %>%
        rbind(deg %>%
          filter(k2) %>%
          arrange(log2FoldChange) %>%
          tail(15)) %>%
        rbind(deg %>%
          filter(k1) %>%
          arrange(padj) %>%
          head(10)) %>%
        rbind(deg %>%
          filter(k2) %>%
          arrange(padj) %>%
          head(10)) %>%
        unique.data.frame()
    } else {
      label = deg[deg$gene %in% label, ]
    }
    label$sig = ifelse(label$log2FoldChange < 0, paste0(l2, "_label"), paste0(l1, "_label"))
    p <- ggplot(data = deg, aes(x = log2FoldChange, y = -log10(padj + 1e-300))) +
      ggrastr::geom_point_rast(alpha = 0.5,
        size = 1, aes(color = sig)) +
      ylab("-log10(Padj)") +
      scale_color_manual(values = c(up = "#d88984",
        up_label = "#ce2b2b", `non-sig` = "#cfcdcd", down = "#8abddc", down_label = "#2863a3")) +
      geom_point(size = 1.02, data = label, aes(color = sig)) +
      ggrepel::geom_text_repel(seed = 1,
        aes(label = label$gene), data = label, point.padding = 0.2, force = 10, max.overlaps = 20,
        color = "black", size = 4, label.size = 0.001, segment.size = 0.3, label.padding = 0.1,
        box.padding = 0.3) +
      theme_classic() +
      labs(x = expression(paste(Log[2], " Fold Change")),
        y = expression(paste(-Log[10], "Padj"))) +
      theme(legend.position = "none") +
      theme(panel.border = element_rect(fill = NA,
        color = "black", linewidth = 1.3), axis.line = element_blank())
  }
  p = Volcano(resM_CTSL, "M_CTSL_KO", label = genelists)

  ggsave(file.path(output_dir, "F5G.pdf"), plot = p, width = 4.5, height = 4)
}

# F5H
plot_F5H <- function(data1, data2, output_dir = ".") {
  p1 <- ggplot(data1, aes(x = x, y = runningScore, color = Description)) + ggrastr::rasterize(geom_line(size = 1,
    alpha = 0.95), dpi = 300) +
    geom_hline(yintercept = 0, lty = 3, col = "black", lwd = 0.6) +
    scale_color_manual(values = c("#a45f5d", "#71b28a", "#f5b386", "#8abddc")) +
    theme(legend.position = "right") +
    theme_classic() +
    labs(x = "", y = "Enrichment score") +
    theme(panel.border = element_rect(fill = NA,
      color = "black", linewidth = 1.2), axis.line = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank(), plot.margin = margin(t = 0, r = 0, b = -20, l = 0),
    axis.text = element_text(color = "black")) +
    guides(color = guide_legend(title = ""))
  p2 <- ggplot(data2, aes(x = x, ymin = -ymin, ymax = -ymax, color = Description)) + ggrastr::rasterise(geom_linerange(alpha = 1,
    size = 0.15), dpi = 300) +
    scale_color_manual(values = c("#a45f5d", "#71b28a", "#f5b386",
      "#8abddc")) +
    theme_classic() +
    labs(x = "") +
    theme(panel.border = element_rect(fill = NA,
      color = "black", linewidth = 1.2), axis.line = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), plot.margin = margin(t = -20, r = 0, b = 0, l = 0),
    axis.text = element_text(color = "black"), legend.position = "none")
  library(patchwork)
  combined_plot <- p1 / p2 + plot_layout(heights = c(7, 3)) +
    plot_layout(guides = "collect")

  ggsave(file.path(output_dir, "F5H.pdf"), plot = combined_plot, width = 5.2, height = 3.3)
}

# F5I
plot_F5I <- function(Cellratio, cols, output_dir = ".") {
  p <- Cellratio %>%
    ggplot(aes(x = Group, y = Ratio, fill = Celltype, stratum = Celltype, alluvium = Celltype)) +
    geom_col(width = 0.5, color = NA) +
    geom_flow(width = 0.5, alpha = 0.3, knot.pos = 0.35) +
    scale_fill_manual(values = cols) +
    theme_void(base_size = 14) +
    theme(axis.text.x = element_text(size = 14,
      vjust = 5), legend.position = "right")
  ggsave(file.path(output_dir, "F5I.pdf"), plot = p, width = 4, height = 3.5)
}

# F5J
plot_F5J <- function(sce.macro, output_dir = ".") {

  p = DimPlot(sce.macro, group.by = "subtype", cols = paletteer::paletteer_d("RColorBrewer::Set3")[-1],
    raster = T, pt.size = 3) +
    labs(title = "")

  ggsave(file.path(output_dir, "F5J.pdf"), plot = p, width = 6.2, height = 4.3)
}

# F5K
plot_F5K <- function(umap, output_dir = ".") {
  p1 <- ggplot(umap, aes(UMAP_1, UMAP_2)) +
    ggrastr::geom_point_rast(color = "#d4d4d4", size = 0.3,
      alpha = 0.1) +
    stat_density2d(data = filter(umap, treat == "CTRL"), aes(fill = treat,
      color = treat), size = 0.7, geom = "polygon", alpha = 0.3, bins = 4) +
    scale_color_manual(values = c(CTRL = "#113264",
      CKO = "#cf5b53"), na.value = NA) +
    scale_fill_manual(values = c(CTRL = "#113264", CKO = "#cf5b53"),
      na.value = NA) +
    theme_void(base_size = 14) +
    theme(legend.title = element_blank(),
      plot.margin = margin(r = 10))
  p <- p1 + stat_density2d(data = filter(umap, treat == "CKO"), aes(fill = treat, color = treat),
    size = 0.7, geom = "polygon", alpha = 0.3, bins = 4)
  ggsave(file.path(output_dir, "F5K.pdf"), plot = p, width = 4.8, height = 3.3)
}

# F5L
plot_F5L <- function(df1, df2, output_dir = ".") {
  plot_go = function(df, gocolor, gotitle) {
    df$LogP = -log10(df$p.adjust)
    df$labelx = rep(0.02 * max(df$LogP), nrow(df))
    df$labely = seq(nrow(df), 1) - 0.2
    df$Description = paste0(stringr::str_wrap(df$Description, width = 35), "\n")
    ggplot(data = df, aes(LogP, reorder(Description, LogP))) +
      geom_bar(stat = "identity",
        alpha = 0.3, fill = gocolor, width = 0.8) +
      geom_text(aes(x = labelx, y = labely,
        label = Description), size = 6 / .pt, family = "Arial", hjust = 0, vjust = 0.55) +
      theme_classic() +
      theme(text = element_text(size = 6, family = "Arial"), plot.title = element_text(size = 6.5,
        family = "Arial", margin = margin(b = 0.5)), axis.text.y = element_blank(), axis.line.y = element_blank(),
      axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line(colour = "black",
        size = 0.25), axis.text.x = element_text(colour = "black"), axis.ticks.x = element_line(colour = "black",
        size = 0.2), axis.ticks.length.x = unit(0.6, "mm")) +
      xlab("-log10(adj.p)") +
      ggtitle(gotitle) +
      scale_x_continuous(expand = c(0, 0))
  }

  plot_go(df1, "#6ac0bf", "Macro_Mrc1") -> p1
  plot_go(df2, "#e0c988", "Macro_Cxcl9") -> p2
  p = p1 + p2 + plot_layout(ncol = 2)

  ggsave(file.path(output_dir, "F5L.pdf"), plot = p, width = 3.5, height = 2.1)
}

# F5M
plot_F5M <- function(sce.macro, output_dir = ".") {

  plots_violins = VlnPlot(sce.macro, group.by = "treat", features = c("Cxcl9", "Gbp2", "Cx3cr1",
    "Tgfbr1"), pt.size = 0, ncol = 2, combine = FALSE, cols = c("#8bbae4", "#e29a92"))
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

  ggsave(file.path(output_dir, "F5M.pdf"), plot = p, width = 5, height = 5)
}

# F5N
plot_F5N <- function(joined_results, output_dir = ".") {

  interesting_factors <- c("STAT6", "STAT3", "STAT1")
  joined_results <- joined_results %>%
    mutate(interesting_factor = ifelse(factor %in% interesting_factors, factor, ifelse(log10.UP.mouse >
      6 & log10.DOWN.mouse > 12.5, "highlight", "other")))
  palette <- c(STAT6 = "#64dfbe", STAT3 = "#9359d6", STAT1 = "#e99165", highlight = "#e7cd86",
    other = "gray80")
  p = ggplot() +
    geom_point(data = joined_results %>%
      filter(interesting_factor %in% c("highlight", "other")), aes(x = log10.DOWN.mouse,
      y = log10.UP.mouse, color = interesting_factor), alpha = 0.6, size = 1.5) +
    geom_point(data = joined_results %>%
      filter(!interesting_factor %in% c("highlight", "other")), aes(x = log10.DOWN.mouse,
      y = log10.UP.mouse, color = interesting_factor), size = 2.5) +
    ggrepel::geom_text_repel(seed = 1,
      data = joined_results %>%
        filter((log10.UP.mouse > 6 & log10.DOWN.mouse > 12.5) | factor == "STAT6"), aes(x = log10.DOWN.mouse,
        y = log10.UP.mouse, label = factor)) +
    geom_vline(xintercept = 12.5, linetype = "dashed",
      color = "gray50") +
    geom_hline(yintercept = 6, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = palette) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none",
      axis.text = element_text(color = "black"), panel.grid = element_blank()) +
    labs(x = "Pred. factors for down-reg. genes (-log10p)",
      y = "Pred. factors for up-reg. genes (-log10p)")

  ggsave(file.path(output_dir, "F5N.pdf"), plot = p, width = 4.1, height = 3.7)
}

# S5J
plot_S5J <- function(sce.macro, subtype_markers, output_dir = ".") {

  p = DotPlot(sce.macro, subtype_markers, group.by = "subtype") +
    RotatedAxis() +
    theme(panel.border = element_rect(color = "black"),
      panel.grid.minor = element_line(color = "gray90", size = 0.25), panel.spacing = unit(1,
        "mm"), strip.text = element_text(margin = margin(b = 3, unit = "mm")), strip.placement = "outlet",
      axis.line = element_blank(), ) +
    labs(x = "", y = "") +
    scale_color_gradient2(low = "#1f273a",
      mid = "#fafafa", high = "#790202", midpoint = -0.2)

  ggsave(file.path(output_dir, "S5J.pdf"), plot = p, width = 15.3, height = 4.5)
}

# S5K
plot_S5K <- function(umap, stype_color, output_dir = ".") {
  dens_layers <- umap %>%
    split(.$subtype) %>%
    purrr::imap(~ {
      bins_g <- 5
      stat_density_2d(data = .x, aes(x = UMAP_1, y = UMAP_2, fill = subtype, color = subtype),
        geom = "polygon", alpha = 0.25, bins = bins_g, size = 0.6)
    })
  p <- ggplot(umap, aes(UMAP_1, UMAP_2)) +
    ggrastr::geom_point_rast(color = "#d4d4d4", size = 0.6,
      alpha = 0.1) + dens_layers + scale_color_manual(values = stype_color, na.value = NA) +
    scale_fill_manual(values = stype_color, na.value = NA) +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(), legend.title = element_blank(), axis.text = element_text(size = 9,
      color = "black"), axis.title = element_text(size = 10)) +
    labs(x = "UMAP 1", y = "UMAP 2")

  ggsave(file.path(output_dir, "S5K.pdf"), plot = p, width = 5.8, height = 3.2)
}

# S5L
plot_S5L <- function(Cellratio, cols, output_dir = ".") {
  Cellratio %>%
    ggplot(aes(x = Group, y = Ratio, fill = Celltype, stratum = Celltype, alluvium = Celltype)) +
    geom_col(width = 0.5, color = NA) +
    geom_flow(width = 0.5, alpha = 0.3, knot.pos = 0.35) +
    scale_fill_manual(values = cols) +
    theme_void(base_size = 13) +
    theme(axis.text.x = element_text(size = 12,
      vjust = 5), legend.position = "right", legend.title = element_blank(), plot.margin = margin(r = 10)) ->
  p

  ggsave(file.path(output_dir, "S5L.pdf"), plot = p, width = 4, height = 3.5)
}

# S5M
plot_S5M <- function(sce.macro, output_dir = ".") {

  p = plot_density(sce.macro, reduction = "umap", features = c("Mrc1", "Apoe", "Trem2", "Hspa1a",
    "Hsph1", "Bag3", "Cxcl9", "Gbp4", "Ccl5"), pal = "inferno", raster = TRUE)
  p = lapply(p, function(x) x + theme_void() +
    theme(plot.title = element_text(hjust = 0.4),
      plot.margin = margin(l = 15, b = 20)))
  p = patchwork::wrap_plots(p)

  ggsave(file.path(output_dir, "S5M.pdf"), plot = p, width = 9, height = 7)
}

# S5N
plot_S5N <- function(de_markers, anti_tumor_TAM_markers, pro_tumor_TAM_markers, output_dir = ".") {

  library(ggrepel)
  de_markers <- de_markers %>%
    mutate(group = case_when(gene %in% anti_tumor_TAM_markers ~ "Anti-tumor", gene %in%
      pro_tumor_TAM_markers ~ "Pro-tumor", TRUE ~ "Other"))
  fs = 6
  p = ggplot(de_markers, aes(x = avg_log2FC, y = -log10(p_val_adj + 9.99988867182683e-321))) +
    ggrastr::geom_point_rast(color = "grey", alpha = 0.25, size = 0.3) +
    geom_point(data = de_markers %>%
      filter(group != "Other" & -log10(p_val_adj + 9.99988867182683e-321) > 2), aes(color = group),
    size = 0.4) +
    scale_color_manual(values = c(`Anti-tumor` = "tomato", `Pro-tumor` = "steelblue",
      Other = "grey")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey",
      size = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey",
      size = 0.2) +
    geom_text_repel(fontface = "italic", data = de_markers %>%
      filter(group != "Other"), aes(label = gene), segment.size = 0.15, segment.alpha = 0.4,
    size = fs / .pt, family = "Arial", max.overlaps = 20, min.segment.length = 0.3, force_pull = 10,
    force = 20, point.padding = 0.1) +
    theme_bw() +
    theme(legend.title = element_blank(),
      text = element_text(family = "Arial", size = fs), panel.grid = element_blank(), legend.position = "none",
      axis.ticks = element_line(size = 0.2), axis.text = element_text(color = "black")) +
    labs(x = "Average log2FC", y = "-log10 adjusted P-value")

  ggsave(file.path(output_dir, "S5N.pdf"), plot = p, width = 2.8, height = 2.1)
}

# S5O
plot_S5O <- function(gse_res2, output_dir = ".") {
  p <- gse_res2 %>%
    mutate(logFDR = -log10(p.adjust), logq = -log10(qvalues)) %>%
    ggplot(aes(x = NES, y = reorder(Description, NES), size = logq, color = logFDR)) +
    geom_point() +
    scale_color_gradient(low = "#e2c3cb", high = "#a0304c") +
    scale_size(range = c(3,
      5)) +
    scale_x_continuous(limits = c(2.57, 2.83)) +
    labs(x = "NES", y = "", color = "-log10(FDR)",
      size = "-log10(qvalue)") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 12, vjust = 0.7), plot.margin = margin(t = 10))

  ggsave(file.path(output_dir, "S5O.pdf"), plot = p, width = 6.5, height = 3.5)
}

# S5Q
plot_S5Q <- function(df1, output_dir = ".") {
  plot_babble_asinh = function(df, targetname) {
    df$Pval = 10^(-abs(df$Pval))
    df$Regulation <- ifelse(df$Coef > 0, "Upregulated", "Downregulated")
    df <- df %>%
      arrange(Coef)
    ggplot(df, aes(x = reorder(Gene, Coef), y = asinh(Coef / 0.01), size = -log10(Pval),
      color = Regulation)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c(Upregulated = "firebrick",
        Downregulated = "steelblue")) +
      scale_size_continuous(limits = c(0, 90), range = c(1.5,
        2.8), breaks = c(30, 60, 90)) +
      theme_linedraw(base_size = 7, base_family = "Arial") +
      theme(panel.grid = element_blank()) +
      labs(y = paste0("Scaled perturbation effects of ",
        targetname), x = "") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1,
        size = 6.5), legend.position = "right", legend.key.size = unit(0.2, "cm")) +
      geom_hline(yintercept = 0,
        linetype = "dashed", color = "darkgray")
  }
  p = plot_babble_asinh(df1, "STAT3")
  ggsave(file.path(output_dir, "S5Q.pdf"), plot = p, width = 2.8, height = 1.8)
}

# S5S
plot_S5S <- function(df_focus, label_genes, rho_rev, slope_rev, lim_rev, pct_without, pct_reversed,
    output_dir = ".") {
  cell_cols <- c(`Stronger CKO effect without STAT3i` = "#B4556C", `Stronger CKO effect with STAT3i` = "#4C6A92",
    `Direction reversed` = "#9D7FA8")
  cell_labels <- c(`Direction reversed` = "Reversed", `Stronger CKO effect with STAT3i` = "Stat3i stronger",
    `Stronger CKO effect without STAT3i` = "no Stat3i stronger")
  group_cols <- c(CTRL = "#2A86CC", CKO = "#E96858", `CTRL + STAT3i` = "#19A18D", `CKO + STAT3i` = "#A688D8")
  theme_cell <- theme_classic(base_size = 6, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.65),
      axis.ticks = element_line(linewidth = 0.6), axis.text = element_text(color = "black",
        size = 6), axis.title = element_text(color = "black", size = 6), legend.title = element_blank(),
      legend.text = element_text(size = 6), legend.key.width = grid::unit(13, "pt"), plot.title = element_text(face = "bold",
        size = 6), plot.margin = margin(4, 5, 4, 4))
  p_fc <- ggplot(df_focus, aes(FC_inh, FC_CTSL)) +
    geom_hline(yintercept = 0, color = "#D9D9D9",
      linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "#D9D9D9", linewidth = 0.4) +
    geom_abline(slope = 1, linetype = "22", color = "#A6A6A6", linewidth = 0.55) +
    geom_abline(slope = slope_rev,
      color = "#262626", linewidth = 0.9) +
    geom_point(aes(fill = state), shape = 21, size = 1.35,
      stroke = 0, alpha = 0.72) +
    geom_text_repel(data = df_focus %>%
      filter(gene %in% label_genes), aes(label = gene), size = 2.1, fontface = "bold.italic",
    color = "#111111", box.padding = 0.42, point.padding = 0.32, segment.size = 0.4, segment.color = "#666666",
    max.overlaps = Inf, seed = 1) +
    scale_fill_manual(values = cell_cols, labels = cell_labels) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + coord_equal(xlim = c(-lim_rev,
      lim_rev), ylim = c(-lim_rev, lim_rev), clip = "off") +
    labs(x = expression(log[2] *
      FC ~ "(Ctsl CKO vs CTRL, Il10-Stat3i)"), y = expression(log[2] * FC ~ "(Ctsl CKO vs CTRL, Il10)"),
    subtitle = sprintf("Ctsl-CKO genes: rho = %.2f; slope = %.2f\n%.1f%% stronger without Stat3i; %.1f%% reversed",
      rho_rev, slope_rev, pct_without, pct_reversed)) +
    theme_cell + theme(legend.position = "bottom",
      legend.justification = "left", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0,
        0, 0, 0), legend.key.height = grid::unit(6, "pt"), legend.spacing.y = grid::unit(0,
        "pt"), plot.subtitle = element_text(size = 6, hjust = 0, margin = margin(b = 2)))
  ggsave(file.path(output_dir, "S5S.pdf"), plot = p_fc, width = 4.4, height = 4.6)
}

# S5T
plot_S5T <- function(df_focus, output_dir = ".") {
  theme_cell <- theme_classic(base_size = 6, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.65),
      axis.ticks = element_line(linewidth = 0.6), axis.text = element_text(color = "black",
        size = 6), axis.title = element_text(color = "black", size = 6), legend.title = element_blank(),
      legend.text = element_text(size = 6), legend.key.width = grid::unit(13, "pt"), plot.title = element_text(face = "bold",
        size = 6), plot.margin = margin(4, 5, 4, 4))
  p_retention <- ggplot(df_focus %>%
    filter(retention < 3), aes(retention)) +
    geom_histogram(aes(fill = after_stat(x) >=
      1), bins = 40, color = "white", linewidth = 0.35, alpha = 0.9) +
    geom_vline(xintercept = 1,
      linetype = "22", linewidth = 0.8, color = "#262626") +
    geom_vline(xintercept = median(df_focus$retention,
      na.rm = TRUE), linewidth = 0.9, color = "#444444") +
    scale_fill_manual(values = c(`FALSE` = "#B4556C",
      `TRUE` = "#4C6A92"), guide = "none") + annotate("text", x = 1.03, y = Inf, label = "Equal CKO effect",
      hjust = 0, vjust = 1.05, size = 2.1, fontface = "bold") + annotate("text", x = 0.04,
      y = Inf, label = "no Stat3i\nstronger", hjust = 0, vjust = 3, size = 2.1, color = "#B4556C",
      fontface = "bold", lineheight = 0.9) + annotate("text", x = 2.35, y = Inf, label = "Stat3i\nstronger",
      hjust = 1, vjust = 3, size = 2.1, color = "#4C6A92", fontface = "bold", lineheight = 0.9) +
    labs(x = "|log2FC| (Stat3i) / |log2FC| (no Stat3i)", y = "Number of genes") +
    theme_cell +
    theme(legend.position = "none", axis.title = element_text(size = 6), axis.text = element_text(size = 6))
  ggsave(file.path(output_dir, "S5T.pdf"), plot = p_retention, width = 3, height = 2.6)
}

# S5U
plot_S5U <- function(df_focus, padj_label_genes, output_dir = ".") {
  cell_cols <- c(`Stronger CKO effect without STAT3i` = "#B4556C", `Stronger CKO effect with STAT3i` = "#4C6A92",
    `Direction reversed` = "#9D7FA8")
  cell_labels <- c(`Direction reversed` = "Reversed", `Stronger CKO effect with STAT3i` = "Stat3i stronger",
    `Stronger CKO effect without STAT3i` = "no Stat3i stronger")
  group_cols <- c(CTRL = "#2A86CC", CKO = "#E96858", `CTRL + STAT3i` = "#19A18D", `CKO + STAT3i` = "#A688D8")
  theme_cell <- theme_classic(base_size = 6, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.65),
      axis.ticks = element_line(linewidth = 0.6), axis.text = element_text(color = "black",
        size = 6), axis.title = element_text(color = "black", size = 6), legend.title = element_blank(),
      legend.text = element_text(size = 6), legend.key.width = grid::unit(13, "pt"), plot.title = element_text(face = "bold",
        size = 6), plot.margin = margin(4, 5, 4, 4))
  sig_line <- -log10(0.05)
  p_padj <- ggplot(df_focus, aes(padj_inh_cap, padj_CTSL_cap)) +
    geom_hline(yintercept = sig_line,
      color = "#D9D9D9", linewidth = 0.4) +
    geom_vline(xintercept = sig_line, color = "#D9D9D9",
      linewidth = 0.4) +
    geom_abline(slope = 1, linetype = "22", color = "#9E9E9E", linewidth = 0.55) +
    geom_point(aes(fill = state), shape = 21, size = 0.9, stroke = 0, alpha = 0.52) +
    geom_text_repel(data = df_focus %>%
      filter(gene %in% padj_label_genes), aes(label = gene), size = 2.1, fontface = "bold.italic",
    color = "#111111", box.padding = 0.35, point.padding = 0.25, segment.size = 0.35, segment.color = "#666666",
    max.overlaps = Inf, seed = 1) +
    scale_fill_manual(values = cell_cols, labels = cell_labels) +
    coord_equal(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
    labs(x = expression(-log[10] *
      "(padj, CtslCKO-Il10-Stat3i)"), y = expression(-log[10] * "(padj, CtslCKO-Il10)"),
    caption = sprintf("%.1f%% more significant\nwithout Stat3i", mean(df_focus$padj_CTSL_log >
      df_focus$padj_inh_log, na.rm = TRUE) * 100)) +
    theme_cell + theme(legend.position = "none",
      plot.caption = element_text(size = 6, hjust = 0, margin = margin(t = -3, b = 0)), axis.title = element_text(size = 6),
      axis.text = element_text(size = 6))
  ggsave(file.path(output_dir, "S5U.pdf"), plot = p_padj, width = 3.5, height = 4)
}

# S5V
plot_S5V <- function(ora_show, pct_without, output_dir = ".") {
  theme_cell <- theme_classic(base_size = 6, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.65),
      axis.ticks = element_line(linewidth = 0.6), axis.text = element_text(color = "black",
        size = 6), axis.title = element_text(color = "black", size = 6), legend.title = element_blank(),
      legend.text = element_text(size = 6), legend.key.width = grid::unit(13, "pt"), plot.title = element_text(face = "bold",
        size = 6), plot.margin = margin(4, 5, 4, 4))
  p_go <- ggplot(ora_show, aes(signed_score, Description, fill = Regulation)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, color = "#555555", linewidth = 0.5) +
    scale_x_continuous(limits = c(-5,
      2.5), breaks = c(-5, -2.5, 0, 2.5), labels = ~ abs(.x), oob = scales::squish) +
    scale_y_discrete(labels = ~ stringr::str_wrap(.x,
      32)) +
    scale_fill_manual(values = c(`Down in CKO` = "#3E8FA3", `Up in CKO` = "#E06965"),
      labels = c(`Down in CKO` = "Down", `Up in CKO` = "Up"), drop = FALSE) +
    labs(x = "-log10(adjusted P)",
      y = NULL, title = sprintf("Stronger without\nStat3i (%.1f%%)", pct_without)) +
    theme_cell +
    theme(legend.position = "top", legend.justification = "center", legend.key.height = grid::unit(5,
      "pt"), legend.key.width = grid::unit(8, "pt"), legend.spacing.x = grid::unit(2,
      "pt"), legend.margin = margin(3, 0, 1, 0), axis.text.y = element_text(size = 6),
    axis.text.x = element_text(size = 6), axis.title.x = element_text(size = 6), plot.title = element_text(hjust = 0.5,
      size = 6, margin = margin(b = 4)))
  ggsave(file.path(output_dir, "S5V.pdf"), plot = p_go, width = 4, height = 4.5)
}

# S5W
plot_S5W <- function(gene_long, output_dir = ".") {
  cell_cols <- c(`Stronger CKO effect without STAT3i` = "#B4556C", `Stronger CKO effect with STAT3i` = "#4C6A92",
    `Direction reversed` = "#9D7FA8")
  cell_labels <- c(`Direction reversed` = "Reversed", `Stronger CKO effect with STAT3i` = "Stat3i stronger",
    `Stronger CKO effect without STAT3i` = "no Stat3i stronger")
  group_cols <- c(CTRL = "#2A86CC", CKO = "#E96858", `CTRL + STAT3i` = "#19A18D", `CKO + STAT3i` = "#A688D8")
  theme_cell <- theme_classic(base_size = 6, base_family = "sans") +
    theme(axis.line = element_line(linewidth = 0.65),
      axis.ticks = element_line(linewidth = 0.6), axis.text = element_text(color = "black",
        size = 6), axis.title = element_text(color = "black", size = 6), legend.title = element_blank(),
      legend.text = element_text(size = 6), legend.key.width = grid::unit(13, "pt"), plot.title = element_text(face = "bold",
        size = 6), plot.margin = margin(4, 5, 4, 4))
  make_gene_plot <- function(g, show_x = FALSE, show_y = FALSE) {
    d <- gene_long %>%
      filter(Gene == g)
    ggplot(d, aes(Group, Expression, fill = Group)) +
      geom_boxplot(width = 0.67, alpha = 0.25,
        outlier.shape = NA, color = "#303030", linewidth = 0.4) +
      geom_point(position = position_jitter(width = 0.045,
        height = 0, seed = 11), shape = 21, size = 1.2, stroke = 0.1) +
      stat_compare_means(comparisons = list(c("CTRL + STAT3i",
        "CKO + STAT3i"), c("CTRL", "CTRL + STAT3i"), c("CTRL", "CKO")), method = "t.test",
      label = "p.signif", size = 2.1, fontface = "bold", bracket.size = 0.65, tip.length = 0.015) +
      scale_fill_manual(values = group_cols, guide = "none") +
      scale_x_discrete(labels = c(CTRL = "CTRL",
        CKO = "CKO", `CTRL + STAT3i` = "CTRL + Stat3i", `CKO + STAT3i` = "CKO + Stat3i")) +
      scale_y_continuous(expand = expansion(mult = c(0.08, 0.3))) +
      labs(x = NULL, y = if (show_y)
        expression(log[2](TPM + 1)) else NULL, title = g) +
      theme_cell + theme(plot.title = element_text(face = "bold.italic",
        hjust = 0.5, size = 6), axis.text.x = if (show_x)
        element_text(angle = 42, hjust = 1, size = 5.3) else element_blank(), axis.ticks.x = if (show_x)
        element_line(linewidth = 0.6) else element_blank(), axis.title.y = element_text(size = 6), axis.text.y = element_text(size = 6),
      plot.margin = margin(1, 2, 1, 2))
  }
  p_arg2 <- make_gene_plot("Arg2", show_x = FALSE, show_y = TRUE)
  p_ccl8 <- make_gene_plot("Ccl8", show_x = FALSE, show_y = TRUE)
  p_mmp14 <- make_gene_plot("Mmp14", show_x = TRUE, show_y = TRUE)
  p_genes <- p_arg2 / p_ccl8 / p_mmp14 + plot_layout(heights = c(1, 1, 1))
  ggsave(file.path(output_dir, "S5W.pdf"), plot = p_genes, width = 2.3, height = 6.5)
}
