suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ggrepel)
  library(corrplot)
  library(pheatmap)
})
if (requireNamespace("extrafont", quietly = TRUE)) {
  extrafont::loadfonts(quiet = TRUE)
}

# F2B
plot_F2B <- function(TME0, module_list, file = "F2B.pdf") {
  pp0 <- list()
  for (i in 1:length(module_list)) {
    module_name <- names(module_list)[i]
    df <- data.frame(TME0@meta.data, TME0@reductions$umap@cell.embeddings)
    df$score <- df[[module_name]]
    limits_max <- c(0.6, 0.6, 0.8, 1, 0.9, 1.4, 1.2, 0.8, 1)
    interp_data <- MBA::mba.surf(
      cbind(df$UMAP_1, df$UMAP_2, df$score),
      no.X = 30, no.Y = 30, extend = TRUE
    )$xyz.est

    contour_df <- expand.grid(x = interp_data$x, y = interp_data$y)
    contour_df$z <- as.vector(interp_data$z)

    pp0[[module_name]] <- ggplot(df, aes(UMAP_1, UMAP_2)) +
      ggrastr::geom_point_rast(aes(color = score), size = 0.4) +
      viridis::scale_color_viridis(
        option = "G", limits = c(min(df$score), limits_max[i])
      ) + geom_contour(
        data = contour_df, aes(x = x, y = y, z = z),
        color = "white", bins = 3, size = 0.7, alpha = 0.6
      ) + theme_void() + theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 12), legend.text = element_text(size = 9)
      ) + guides(color = guide_colorbar(barwidth = 0.8, barheight = 4)) +
      labs(title = gsub("_", " ", module_name))
  }
  p_final <- wrap_plots(pp0)
  ggsave(file, p_final, width = 9, height = 6)
}

# F2C
plot_F2C <- function(cor.guide, tb, file = "F2C.pdf") {
  pdf(file, width = 5, height = 5)
  corrplot(
    cor.guide[tb$genes, tb$genes], method = "square", tl.pos = "n", type = "lower",
    cl.cex = 0.8, col = colorRampPalette(c("#1E3163", "#99eff7", "#f8edba",
      "#f35a13"))(10),
    mar = c(4, 4, 4, 4), cl.pos = "b", cl.length = 5, cl.ratio = 0.1
  ) %>%
    corrRect(c(1, 21, 35, 47, 64, 80), lwd = 2)
  gene_anno <- rev(tb$.)
  anno_colors <- c(
    `1` = "#A2B8A2", `2` = "#e0b6db", `3` = "#f5df91", `4` = "#bbaa9a", `5` = "#f5d3c4"
  )
  anno_vector <- anno_colors[gene_anno]
  for (i in 1:length(anno_colors)) {
    class_genes <- which(gene_anno == names(anno_colors)[i])
    rect(
      -1, min(class_genes) - 0.5, 0, max(class_genes) +
        0.5, col = anno_colors[names(anno_colors)[i]],
      border = NA
    )
    text(
      -2.5, mean(class_genes), labels = names(anno_colors)[i], col = "black", cex = 1
    )
  }
  dev.off()
}

# F2D
plot_F2D <- function(SPP1_df, file = "f2_TME2_M6_SPP1_coef_vln_asinh_color2.pdf") {
  anno_colors_dark <- rev(c("#559cc0", "#e0848e", "#897abb", "#d8b36e", "#6da87e"))
  p = ggplot(SPP1_df, aes(x = group, y = asinh(coef / 0.01))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#868686", size = 0.5) +
    geom_violin(
      aes(group = group, fill = factor(group)),
      width = 0.6, alpha = 0.1, color = "#9c9c9c", size = 0.4
    ) + geom_jitter(
      width = 0.1, aes(color = factor(group)),
      shape = 21, fill = "white", size = 0.6, stroke = 0.6, alpha = 1
    ) + scale_fill_manual(values = anno_colors_dark) +
    scale_color_manual(values = anno_colors_dark) +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) +
    labs(x = "", y = "Scaled effects on MP6_SPP1") + theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6),
      plot.margin = margin(t = 5, r = 20, b = 5, l = -2),
      legend.position = "none", axis.line = element_blank(),
    ) + coord_flip()
  pdf(file = file, width = 1.4, height = 2.4)
  print(p)
  dev.off()
}

# F2F
plot_F2F <- function(prtb_ntc_df, plot_eg_df, x_axis = 2.7, file = "F2F.pdf") {
  p1 = prtb_ntc_df %>%
    ggplot(.) +
    geom_density(aes(x = effects), fill = "grey", linewidth = 0.2, bw = 0.07) +
    xlim(c(-x_axis, x_axis)) + theme_linedraw(base_size = 6.5, base_family = "Arial") +
    theme(plot.margin = margin(0, 0, 0, 0)) + theme(
      panel.grid = element_blank(), panel.border = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(), axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  p1.1 <- prtb_ntc_df %>%
    ggplot(., aes(x = effects)) +
    geom_vline(aes(xintercept = effects), alpha = 0.5, color = "grey", size = 0.3) +
    xlim(c(-x_axis, x_axis)) + theme_linedraw(base_size = 6.5, base_family = "Arial") +
    theme(plot.margin = margin(0, 0, 0, 0)) + theme(
      axis.title.x = element_blank(), axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(), panel.border = element_rect(size = 0.4,
        color = "black"),
      panel.spacing = unit(0, "pt")
    )
  p2 = plot_eg_df %>%
    ggplot(., aes(x = effects)) + geom_vline(
      aes(xintercept = effects, color = line_color), size = 0.75, alpha = 0.7
    ) + xlim(c(-x_axis, x_axis)) + theme_bw(base_size = 6.5, base_family = "Arial") +
    theme(panel.background = element_rect(color = "black")) +
    facet_wrap(. ~ guide_group, ncol = 1, strip.position = "left") +
    theme(plot.margin = margin(r = 10), panel.spacing = unit(0, "pt")) + theme(
      axis.text = element_text(color = "black"),
      strip.background = element_blank(), strip.placement = "outside",
      strip.text.y.left = element_text(size = 6, family = "Arial", angle = 0,
        color = "black"),
      panel.grid = element_blank(), panel.border = element_rect(size = 0.4,
        color = "black"),
      panel.spacing = unit(0.1, "lines"), legend.position = "none"
    ) + xlab("") + scale_y_continuous(expand = c(0, 0)) +
    scale_color_manual(values = c("#ec936d", "#dfdfdf", "#1b536e"))
  p <- p1 / p1.1 / p2 + plot_layout(heights = c(1.9, 0.5, 3.4))
  ggsave(file, p, width = 2.3, height = 1.8)
}

# F2G
plot_F2G <- function(t1, xlab = "Scaled perturbation effects on MP6_SPP1",
    file = "F2G.pdf") {
  p <- ggplot() + geom_point(
    data = t1, mapping = aes(x = effects, y = guide_group, color = max),
    size = 1.2, alpha = 0.7
  ) + geom_line(
    data = t1, mapping = aes(x = effects, y = guide_group, color = max),
    size = 0.9, alpha = 0.5
  ) + scale_color_gradientn(
    colours = c("#15638a", "#4f94b4", "#dfdfdf", "#e69877", "#c76a43"),
    values = scales::rescale(c(min(t1$max), -1.2, 0, 1.2, max(t1$max)))
  ) + theme_linedraw(base_size = 6.5, base_family = "Arial") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) + theme(
      axis.text.y = element_text(size = 6),
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#eeeeee", size = 0.02),
      legend.position = "none", plot.margin = margin(0, 0, 0, 0),
      panel.border = element_rect(size = 0.4, color = "black"),
      panel.spacing = unit(0, "pt")
    ) + xlim(-2.7, 2.7) + labs(x = xlab, y = "")
  ggsave(file, p, width = 2.3, height = 2.25)
}

# F2H
plot_F2H <- function(volcano_df, coef_cut_off = 1, fdr_cut_off = 5,
    file = "F2H.pdf") {
  volcano_df$guide_cutoff_class <- ifelse(
    volcano_df$mlg10fdr > fdr_cut_off & abs(volcano_df$coef) > coef_cut_off,
    "effect", "no effect"
  )
  p1 = ggplot(volcano_df) + ggrastr::geom_point_rast(
    aes(x = coef, y = mlg10fdr, color = guide_cutoff_class), size = 0.7, alpha = 0.6
  ) + geom_vline(
    xintercept = coef_cut_off, linetype = "dashed", color = "#424141", alpha = 0.6,
    size = 0.2
  ) + geom_vline(
    xintercept = -coef_cut_off, linetype = "dashed", color = "#424141", alpha = 0.6,
    size = 0.2
  ) + geom_hline(
    yintercept = fdr_cut_off, linetype = "dashed", color = "#424141", alpha = 0.6,
    size = 0.2
  ) + ggrepel::geom_text_repel(
    data = filter(volcano_df, guide_cutoff_class == "effect"),
    mapping = aes(x = coef, y = mlg10fdr, label = guide_group),
    size = 6 / .pt, family = "Arial", max.overlaps = Inf, segment.size = 0.1,
    min.segment.length = 1
  ) + scale_color_manual(values = c("#bb610d", "grey")) +
    theme_linedraw(base_size = 6.2, base_family = "Arial") +
    theme(panel.grid = element_blank(), legend.position = "none") + xlim(c(-2.5, 1.8)) +
    xlab("Scaled perturbation effects on MP6_SPP1") + ylab("-log10(FDR)")
  ggsave(file, p1, width = 1.8, height = 1.8)
}

# F2L
plot_F2L <- function(plot_df, file = "f2_TME2_narrowdown_asinh.pdf") {
  p = ggplot(plot_df) + geom_point(
    aes(x = asinh(pg1 / 0.01), y = asinh(pg2 / 0.01), color = class2),
    size = 1, alpha = 0.6
  ) + ggrepel::geom_text_repel(
    data = filter(plot_df, class2 != "others"), mapping = aes(
      x = asinh(pg1 / 0.01), y = asinh(pg2 / 0.01), color = class2, label = gene
    ), size = 6 / .pt, family = "Arial", segment.size = 0.2, box.padding = 0.3,
    point.padding = 0.3, force = 2, max.overlaps = Inf
  ) + scale_colour_manual(values = c("#a5336a", "grey")) +
    geom_abline(intercept = -0.005, slope = -1, linetype = "dashed", alpha = 0.3,
      size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3, size = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3, size = 0.3) +
    theme_linedraw(base_size = 6.5, base_family = "Arial") +
    xlab("Guide1 effects on MP6_SPP1") + ylab("Guide2 effects on MP6_SPP1") +
    theme(panel.grid = element_blank(), legend.position = "none")
  ggsave(file, p, width = 1.8, height = 1.75)
}

# F2I
plot_F2I <- function(TREM2_GOdf, file = "F2I.pdf") {
  TREM2_GOdf %>%
    mutate(
      x = -log10(p.adjust) * type, Description = stringr::str_wrap(
        sub("^([a-z])", "\\U\\1", Description, perl = TRUE), width = 28
      ), margin = ifelse(type == 1, -0.2, 0.2), hjust = ifelse(type == 1, 1, 0)
    ) %>%
    arrange(x) %>%
    mutate(Description = factor(Description, levels = Description)) %>%
    ggplot(aes(x, Description, fill = as.character(type))) + geom_col(alpha = 0.6) +
    geom_text(
      aes(x = margin, label = Description, hjust = hjust),
      lineheight = 0.7, show.legend = FALSE, colour = "black", size = 6.5 / .pt
    ) + geom_hline(yintercept = 0) + labs(x = "-log10(p.adjust)", y = NULL) +
    theme_classic(base_size = 7) + scale_fill_manual(
      name = "", values = c("#306682", "#d1493a"), labels = c("down", "up")
    ) + scale_x_continuous(labels = abs, limits = c(-8.2, 7.2)) + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.line.y = element_blank(), plot.margin = margin(t = 0, r = 5, b = 10,
        l = -10),
      panel.border = element_blank(), axis.line.x = element_line(color = "black"),
      legend.key.size = unit(0.3, "cm"),
      legend.position = "top", legend.text = element_text(size = 6.5)
    ) -> p1
  ggsave(file, p1, width = 2.9, height = 2.7)
}

# F2J
plot_F2J <- function(coef_TREM2, file = "f2_TREM2_perturb_top_dotplot_asinh.pdf") {
  ggplot(
    coef_TREM2, aes(x = reorder(gene, -coef), y = asinh(coef / 0.01), color = class)
  ) + coord_cartesian(xlim = c(-220000, nrow(coef_TREM2) + 22000)) +
    ggrastr::geom_point_rast(size = 0.1, alpha = 0.5) +
    geom_point(data = coef_TREM2 %>% filter(class != "NS"), size = 0.7, alpha = 0.8) +
    scale_color_manual(
      values = c(NS = "#c5c5c5d0", neg_perturbed = "#2b6b88", pos_perturbed = "#cd5a4d")
    ) + geom_text_repel(
      data = coef_TREM2 %>%
        filter(class != "NS"),
      aes(label = gene),
      size = 5.5 / .pt, family = "Arial", fontface = "bold", max.overlaps = 20,
      force = 8, box.padding = 0.1, point.padding = 0.01, min.segment.length = 0.2,
      segment.size = 0.2
    ) + labs(x = "Genes", y = "Scaled perturbation effects of TREM2") +
    geom_hline(yintercept = 0, color = "#353434d0", linetype = "dashed", size = 0.3) +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.position = "none"
    ) -> p1
  ggsave(file, p1, width = 1.3, height = 2.3)
}

# F2K
plot_F2K <- function(df, file = "f2_TREM2_mice_compare_KO_WT_SPP1_module.pdf") {
  p = ggplot(df, aes(x = condition, y = module_score)) +
    geom_violin(aes(fill = condition), alpha = 0.7) +
    geom_boxplot(width = 0.25, alpha = 0.7) +
    ggpubr::stat_compare_means(label.x = 1.3, label.y = 1.15) +
    ggpubr::theme_pubr() + ylab("M6_SPP1 score") + labs(x = "") + theme(
      panel.border = element_rect(fill = NA, color = "black"),
      axis.line = element_blank(), legend.key.size = unit(0.4, "cm"),
      legend.position = "top", legend.justification = "right",
      axis.text.x = element_blank()
    ) + scale_y_continuous(expand = expansion(mult = c(0, 0.08))) + scale_fill_manual(
      values = c("#92aa9d", "#f1dca0"), name = ""
    )
  ggsave(p, file = file, width = 3.3, height = 5.2)
}

# S2D
plot_S2D <- function(c2, file = "S2D.pdf") {
  p <- DimPlot(c2, label = T, group.by = "anno")
  p1 = p <- DimPlot(
    c2, group.by = "anno", cols = paletteer::paletteer_d("miscpalettes::pastel",
      direction = -1)
  ) + theme_void() + ggtitle("")
  ggsave(file, p, width = 5, height = 3.6)
}

# S2E
plot_S2E <- function(avg_expression, file = "s2_C2_clustered_heatmap_new.pdf") {
  annotation_col <- data.frame(Cluster = factor(colnames(avg_expression$RNA)))
  rownames(annotation_col) <- colnames(avg_expression$RNA)
  anno_colors <- list(
    Cluster = setNames(
      paletteer::paletteer_d("miscpalettes::pastel",
        direction = -1)[1:length(unique(colnames(avg_expression$RNA)))],
      unique(colnames(avg_expression$RNA))
    )
  )
  p2 = pheatmap(
    avg_expression$RNA, cluster_rows = FALSE, cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE, scale = "row", color = colorRampPalette(c("#4575b4",
      "#f3f2f2", "#d89c99"))(50),
    border_color = NA, annotation_col = annotation_col, annotation_colors = anno_colors,
    annotation_legend = FALSE, fontsize = 13, angle_col = 315, silent = TRUE
  )
  pdf(file, width = 3.6, height = 6.3)
  print(p2)
  dev.off()
}

# S2F
plot_S2F <- function(sig_mat, file = "s2_noTMEvsC2_module_score_new.pdf") {
  p3 <- sig_mat %>%
    pivot_longer(!label) %>%
    ggplot(.) +
    geom_boxplot(aes(x = name, y = value, color = label), outlier.size = 0.5) +
    ggpubr::stat_compare_means(
      aes(x = name, y = value, group = label), method = "t.test", label = "p.signif"
    ) + ggthemes::theme_few() + scale_color_manual(
      values = c("#334A45", "#F3BFA9"), name = ""
    ) + ylab("Module score") + xlab("") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) + theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      legend.position = "top"
    )
  ggsave(file, p3, width = 4.8, height = 4.2)
}

# S2G
plot_S2G <- function(TME, df_freq, file = "s2_TME2_perturbed_umap.pdf") {
  p = DimPlot(TME, group.by = "if_Perturbed", cols = c("#f57f6a", "#bae0e2")) +
    theme(legend.position = c(0.8, 6)) + theme_void() + ggtitle("") + inset_element(
      ggplot(df_freq) +
        geom_col(aes(x = "", y = Freq, fill = Var1), width = 1, alpha = 0.8) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = c("#f57f6a", "#bae0e2")) +
        theme_void() + theme(legend.position = "none"),
      align_to = "plot", left = 0.95, bottom = -0.7, right = 1.4, top = 0.8
    ) + geom_text(
      data = subset(df_freq, Var1 == "Perturbed"), aes(
        x = 1, y = cumsum(Freq)[1] -
          Freq[1] / 2, label = scales::percent(percentage, accuracy = 0.1)
      ), size = 3.5
    ) + theme(plot.margin = margin(0, 1, 1.5, 0.2, "cm"))
  ggsave(file, p, width = 4.9, height = 3.4)
}

# S2H
plot_S2H <- function(hist_df, column = "sg", file = "S2H.pdf") {
  hist_df_mod <- hist_df
  hist_df_mod[hist_df_mod[, column] > 100, column] <- 100
  breaks = seq(-0.5, 102, 3)
  p1 <- ggplot(data = hist_df_mod) + geom_histogram(
    breaks = breaks, alpha = 0.5, linewidth = 0.3, aes(x = hist_df_mod[, column],
      y = stat(count) / sum(count)),
    fill = "#becac0", color = "#2c2c2c"
  ) + ylab("Relative frequency") +
    xlab("Number of detected sgRNA categories per cell") +
    ggthemes::theme_few() + scale_x_continuous(
      breaks = c(seq(0, 80, 20), 100), labels = c(seq(0, 80, 20), ">100")
    )
  xpars <- sapply(c("mean", "median"), do.call, list(x = hist_df[, column])) %>%
    round(., 2) %>%
    paste(names(.), ., sep = " = ") %>%
    paste(., collapse = "\n")
  p1 <- p1 + geom_vline(xintercept = median(hist_df[, column]), linetype = "dashed") +
    annotate(geom = "text", x = 30, y = 0.1, label = xpars, hjust = 0)
  ggsave(file, p1, width = 4.5, height = 2.1)
}

# S2I
plot_S2I <- function(hist_df, column = "gene", file = "S2I.pdf") {
  breaks = seq(-0.5, min(2000, max(hist_df[, column]) + 1), 2)
  p1 <- ggplot(data = hist_df) + geom_histogram(
    breaks = breaks, alpha = 0.6, linewidth = 0.3, aes(x = hist_df[, column],
      y = stat(count) / sum(count)),
    fill = "#e4dfce", color = "#2c2c2c"
  ) + ylab("Relative frequency") +
    xlab("Number of perturbed gene categories per cell") + ggthemes::theme_few()
  xpars <- sapply(c("mean", "median"), do.call, list(x = hist_df[, column])) %>%
    round(., 2) %>%
    paste(names(.), ., sep = " = ") %>%
    paste(., collapse = "\n")
  p1 <- p1 + geom_vline(xintercept = median(hist_df[, column]), linetype = "dashed") +
    annotate(
      geom = "text", x = 3 * median(hist_df[, column]),
      y = 0.08, label = xpars, hjust = 0
    )
  ggsave(file, p1, width = 4.5, height = 2.1)
}

# S2J
plot_S2J <- function(target_effects, file = "S2J.pdf") {
  p <- target_effects %>%
    ggplot(data = .) + geom_point(
      aes(x = reorder(prtb_genes, number), y = coef, color = color),
      size = 1, alpha = 0.8
    ) + scale_color_manual(values = c("#835d70", "darkgrey")) + geom_text_repel(
      data = subset(target_effects, color == "A"),
      aes(x = reorder(prtb_genes, number), y = coef, label = prtb_genes),
      size = 4, nudge_x = 2, max.overlaps = 15, force = 2
    ) + geom_hline(yintercept = 0, color = "grey") + ggthemes::theme_few() + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.position = "none"
    ) + xlab("") + ylab("Pertrurb effects on targets")
  ggsave(file, p, width = 5, height = 2.1)
}

# S2K
plot_S2K <- function(target_control_long, file = "S2K.pdf") {
  p = target_control_long %>%
    ggplot(aes(x = value)) + geom_histogram(
      aes(fill = variable),
      alpha = 0.6, position = "identity", color = "black", bins = 50, linewidth = 0.3
    ) + geom_density(aes(color = variable), size = 0.7) + scale_color_manual(
      values = c("#835d70", "#757574"), name = ""
    ) + ggthemes::theme_few() + geom_vline(xintercept = 0, color = "black") +
    scale_fill_discrete(
      type = c("#835d70", "#fcf4e6"), name = ""
    ) + ylab("Number of targets") + xlab("Pertrurb effects on targets")
  ggsave(file, p, width = 5, height = 2.1)
}

# S2L
plot_S2L <- function(flt_df, file = "S2L.pdf") {
  p1 = ggplot(flt_df, aes(x = magnitude, y = variance)) +
    geom_point(color = "#8f3f8f", size = 2, alpha = 0.7) +
    theme_linedraw() + theme(panel.grid = element_blank()) +
    labs(x = "Magnitude (Mean absolute scaled effects)", y = "Variance") +
    ggrepel::geom_text_repel(aes(label = module), point.padding = 1, size = 3)
  ggsave(file, p1, width = 3.5, height = 4)
}

# S2M
plot_S2M <- function(flt_df, file = "S2M.pdf") {
  p2 = ggplot(flt_df, aes(x = homogeneity, y = heterogeneity)) + geom_point(
    aes(size = magnitude), color = "#419672", alpha = 0.7
  ) + theme_linedraw() + theme(panel.grid = element_blank()) +
    labs(x = "Intra-CM homogeneity", y = "Inter-CM heterogeneity") +
    ggrepel::geom_text_repel(aes(label = module), point.padding = 5, size = 3)
  ggsave(file, p2, width = 3.5, height = 4)
}

# S2N
plot_S2N <- function(coef_means,
    file = "s2_TME2_modulexgroup_meancoef_heatmap_asinh.pdf") {
  anno_colors_dark <- c("#7a9b7a", "#3f71b3", "#d1b226", "#8a8383", "#df9555")
  names(anno_colors_dark) = paste0("CM", c(1:5))
  color_used <- colorRamp2(
    c(-0.8, 0, 0.8), c("#554c7a", "#F0F0F0", "#a87839")
  )
  annotation <- HeatmapAnnotation(
    Module_Group = paste0("CM", c(1:5)), col = list(Module_Group = anno_colors_dark),
    show_annotation_name = FALSE, show_legend = FALSE
  )
  p <- Heatmap(
    coef_means, col = color_used,
    heatmap_legend_param = list(title = "Scaled\neffects"),
    top_annotation = annotation, cluster_rows = F, cluster_columns = F,
    show_row_names = T,
    show_column_names = T, column_names_centered = TRUE, column_names_rot = 0
  )
  pdf(file, width = 3.7, height = 4)
  print(p)
  dev.off()
}

# S2O
plot_S2O <- function(HIF1A_GOdf, file = "s2_HIF1A_GO.pdf") {
  HIF1A_GOdf %>%
    mutate(
      x = -log10(p.adjust) *
        type, Description = stringr::str_wrap(Description, width = 28),
      margin = ifelse(type == 1, -0.2, 0.2), hjust = ifelse(type == 1, 1, 0)
    ) %>%
    arrange(x) %>%
    mutate(Description = factor(Description, levels = Description)) %>%
    ggplot(aes(x, Description, fill = as.character(type))) + geom_col(alpha = 0.6) +
    geom_text(
      aes(x = margin, label = Description, hjust = hjust),
      lineheight = 0.7, show.legend = FALSE, colour = "black", size = 6.5 / .pt
    ) + geom_hline(yintercept = 0) + labs(x = "-log10(p.adjust)", y = NULL) +
    theme_classic(base_size = 7) + scale_fill_manual(
      name = "", values = c("#306682", "#d1493a"), labels = c("down", "up")
    ) + scale_x_continuous(labels = abs, limits = c(-8.2, 7.2)) + theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.line.y = element_blank(), plot.margin = margin(t = 0, r = 5, b = 10,
        l = -10),
      panel.border = element_blank(), axis.line.x = element_line(color = "black"),
      legend.key.size = unit(0.3, "cm"),
      legend.position = "top", legend.text = element_text(size = 6.5)
    ) -> p1
  ggsave(file, p1, width = 2.7, height = 2.8)
}

# S2P
plot_S2P <- function(coef_HIF1A, file = "s2_HIF1A_perturb_top_dotplot_asinh.pdf") {
  ggplot(
    coef_HIF1A, aes(x = reorder(gene, -coef), y = asinh(coef / 0.01), color = class)
  ) + coord_cartesian(xlim = c(-23000, nrow(coef_HIF1A) + 22000)) +
    ggrastr::geom_point_rast(size = 0.1, alpha = 0.5) +
    geom_point(data = coef_HIF1A %>% filter(class != "NS"), size = 0.7, alpha = 0.8) +
    scale_color_manual(
      values = c(NS = "#c5c5c5d0", neg_perturbed = "#2b6b88", pos_perturbed = "#cd5a4d")
    ) + geom_text_repel(
      data = coef_HIF1A %>%
        filter(class != "NS"),
      aes(label = gene),
      size = 5.5 / .pt, family = "Arial", fontface = "bold", max.overlaps = 20,
      force = 10, box.padding = 0.1, point.padding = 0.01, min.segment.length = 0.4,
      segment.size = 0.2
    ) + labs(x = "Genes", y = "Scaled perturbation effects of HIF1A") +
    geom_hline(yintercept = 0, color = "#353434d0", linetype = "dashed", size = 0.3) +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) + theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      legend.position = "none"
    ) -> p1
  ggsave(file, p1, width = 1.5, height = 2.3)
}
