library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggvenn)
library(ggradar)
library(pheatmap)
if (requireNamespace("extrafont", quietly = TRUE)) extrafont::loadfonts(quiet = TRUE)

# F3F
plot_F3F <- function(pca_data, output_dir = ".") {

  color_palette <- c("#E9C46A", "#2A9D8F", "#457B9D", "#E76F51")
  fs = 6
  pca_plot = ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 0.9, alpha = 0.8) +
    stat_ellipse(
      aes(group = Group), linetype = 2, type = "norm", level = 0.9, size = 0.5
    ) +
    labs(title = "", x = "PC1", y = "PC2") +
    scale_color_manual(values = color_palette) +
    theme_classic() + theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.line = element_blank(), axis.text = element_text(color = "black"),
      legend.position = "top", legend.title = element_blank(),
      plot.margin = margin(r = 3)
    ) +
    theme(
      axis.text = element_text(size = fs), plot.title = element_text(size = fs),
      axis.title = element_text(size = fs), legend.text = element_text(size = fs),
      legend.key.size = unit(0.3, "cm"), legend.margin = margin(t = -10, b = -15),
      axis.ticks = element_line(size = 0.2)
    )
  boxplot_y = ggplot(pca_data, aes(x = Group, y = PC2, fill = Group)) +
    geom_violin(alpha = 0.6, width = 0.8, scale = "width", color = NA) +
    labs(y = NULL, x = NULL) +
    theme_classic() + scale_fill_manual(values = color_palette) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.line = element_blank(), axis.text = element_text(color = "black")
    ) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), legend.position = "none"
    ) +
    theme(
      axis.text = element_text(size = fs), plot.title = element_text(size = fs),
      axis.title = element_text(size = fs), legend.text = element_text(size = fs),
      legend.title = element_text(size = fs), legend.key.size = unit(0.2, "cm")
    )
  p = pca_plot + boxplot_y + plot_layout(widths = c(4, 1.1))

  ggsave(file.path(output_dir, "F3F.pdf"), p, width = 2.6, height = 1.9)
}

# F3G
plot_F3G <- function(tpm_mean_all, column_groups, output_dir = ".") {

  ann_col = data.frame(group = c(column_groups))
  rownames(ann_col) = colnames(tpm_mean_all)
  ann_col$group = factor(ann_col$group, levels = unique(ann_col$group))
  group_col = c("#cfc3a6", "#7fb7b1", "#7692a3", "#e0a596")
  names(group_col) = unique(ann_col$group)
  fs = 6
  hc_row = hclust(dist(tpm_mean_all), method = "complete")
  ha = HeatmapAnnotation(
    group = ann_col$group, col = list(group = group_col),
    show_annotation_name = F, simple_anno_size = unit(1.5, "mm"),
    annotation_legend_param = list(
      title = NULL, labels_gp = gpar(fontsize = fs), grid_height = unit(2.7, "mm"),
      grid_width = unit(2.7, "mm")
    )
  )
  col_fun = colorRampPalette(c("#A3B7DC", "#F5F8FA", "#FFA899"))(100)
  ht = Heatmap(
    tpm_mean_all, name = "Expr. Zscore", col = col_fun, cluster_columns = FALSE,
    cluster_rows = T, show_column_names = F, show_row_names = T,
    row_names_gp = gpar(fontsize = fs),
    column_names_gp = gpar(fontsize = fs), row_dend_gp = gpar(lwd = 0.5),
    column_split = ann_col$group, column_title = NULL, column_gap = unit(0.7, "mm"),
    top_annotation = ha, heatmap_legend_param = list(
      labels_gp = gpar(fontsize = fs), at = c(-4, 0, 4), labels = c("-4", "0", "4"),
      title_gp = gpar(fontsize = fs), legend_height = unit(7, "mm"),
      grid_width = unit(2, "mm")
    )
  )

  pdf(file.path(output_dir, "F3G.pdf"), width = 5.6, height = 1.3)
  draw(ht, merge_legends = TRUE)
  dev.off()
}

# F3H
plot_F3H <- function(radar_data_scaled, output_dir = ".") {
  p = ggradar(
    radar_data_scaled, grid.min = 0.1, grid.mid = 0.7, grid.max = 1,
    values.radar = c("", "", ""),
    grid.line.width = 0.2, gridline.min.linetype = "solid",
    gridline.mid.linetype = "solid",
    gridline.max.linetype = "solid", gridline.min.colour = "grey",
    gridline.mid.colour = "grey",
    gridline.max.colour = "black", axis.label.size = 4.3, group.line.width = 0.5,
    group.point.size = 3, background.circle.colour = "white",
    group.colours = c(tam2402 = "#b17e3cff", tam0201 = "#f5ae22ff", N = "#425cabff",
    T = "#d81616ff"),
    fill.alpha = 0.1, fill = TRUE
  ) +
    theme(legend.position = "top")

  ggsave(file.path(output_dir, "F3H.pdf"), p, width = 5.5, height = 4)
}

# F3I
plot_F3I <- function(coef_2M, output_dir = ".") {
  p1 = ggplot(coef_2M, aes(y = asinh(MP1_ISG / 0.01), x = asinh(MP6_SPP1 / 0.01))) +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "dark grey",
    linewidth = 0.3) +
    geom_vline(xintercept = -0.1, linetype = "dashed", color = "dark grey",
    linewidth = 0.3) +
    geom_point(aes(color = color), size = 0.5, alpha = 0.7) +
    scale_color_manual(values = c(darkred = "#7c2727", lightpink = "orange",
    black = "#424242")) +
    ggrepel::geom_text_repel(
      seed = 1, aes(label = gene, color = color),
      fontface = "bold", family = "Arial", size = 6 / .pt, segment.size = 0.2,
      min.segment.length = 0.5, point.padding = 0.05, max.overlaps = 20
    ) +
    xlim(c(-1.55, 0.15)) +
    ylim(c(-0.1, 0.75)) +
    ylab("Scaled effects on MP1_ISG") +
    xlab("Scaled effects on MP6_SPP1 ") +
    ggtitle("") +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "F3I.pdf"), p1, width = 1.8, height = 1.9)
}

# F3J
plot_F3J <- function(gene_lists, output_dir = ".") {
  p <- ggvenn(
    gene_lists, show_elements = TRUE, label_sep = "\n", fill_color = c("#b8cb7d",
    "#e8bd80", "#acc0e8", "#e5c0be"),
    fill_alpha = 0.5, stroke_size = 0, set_name_size = 6, text_size = 6
  )

  ggsave(file.path(output_dir, "F3J.pdf"), p, width = 7.5, height = 7)
}

# F3K
plot_F3K <- function(df, output_dir = ".") {
  colors = colorRampPalette(c("#2e538a", "#2a4670", "#f8f8fc", "#f5b1b0",
  "#a15f5e"))(50)
  breaks = seq(-0.5, 0.52, length.out = 51)

  p = pheatmap::pheatmap(
    df, silent = TRUE, cluster_rows = F, cluster_cols = F, scale = "none",
    show_rownames = T,
    show_colnames = T, fontsize = 10, border_color = "#504f4f", color = colors,
    breaks = breaks, cellwidth = 35, cellheight = 35, fontsize_row = 12,
    fontsize_col = 12,
    angle_col = 45, main = "CTSL Scaled effects"
  )

  ggsave(file.path(output_dir, "F3K.pdf"), p, width = 5, height = 3.5)
}

# S3F
plot_S3F <- function(cor_mat, output_dir = ".") {

  library(ComplexHeatmap)
  library(circlize)
  symbol_mat <- matrix(
    "", nrow = nrow(cor_mat), ncol = ncol(cor_mat), dimnames = dimnames(cor_mat)
  )
  for (i in seq_len(nrow(cor_mat))) {
    symbol_mat[i, which.max(cor_mat[i, ])] <- "*"
  }
  p = Heatmap(
    cor_mat, name = "pearson correlation", column_title = "CRC", row_title = "TAM",
    col = colorRamp2(
      c(-0.7, 0, 0.4, 0.7), c("#5b71a3", "#f3f7ff", "#f5bda7", "#e26736")
    ),
    cluster_rows = FALSE, cluster_columns = FALSE, cell_fun = function(j, i, x, y,
    width, height, fill) {
      if (symbol_mat[i, j] == "*") {
        grid.text("*", x, y, gp = gpar(fontsize = 15, fontface = "bold"))
      }
    }, column_names_rot = 45
  )

  pdf(file.path(output_dir, "S3F.pdf"), width = 6.5, height = 5)
  draw(p)
  dev.off()
}

# S3G
plot_S3G <- function(histograms, output_dir = ".") {
  draw_histogram <- function(target_control_long) {
    p = target_control_long %>%
      ggplot(aes(x = value)) +
      geom_histogram(
        aes(fill = variable),
        alpha = 0.6, position = "identity", color = "black", bins = 50, linewidth = 0.3
      ) +
      geom_density(aes(color = variable), size = 0.7) +
      scale_color_manual(
        values = c("#835d70", "#757574"), name = ""
      ) +
      ggthemes::theme_few() + geom_vline(xintercept = 0, color = "black") +
      scale_fill_discrete(
        type = c("#835d70", "#fcf4e6"), name = ""
      ) +
      ylab("Number of targets") +
      xlab("Pertrurb effects on targets")
    print(p)
  }
  p = draw_histogram(histograms$MDM) / draw_histogram(histograms$TAM)
  ggsave(file.path(output_dir, "S3G.pdf"), p, width = 4.3, height = 3.6)
}

# S3H
plot_S3H <- function(scaled_effects, output_dir = ".") {

  color_used <- colorRamp2(
    c(-0.9, 0, 0.9), c("#4a6285", "#fffbfb", "#743634")
  )
  p = Heatmap(
    scaled_effects, name = "Scaled\nperturbation\neffects", cluster_rows = F,
    cluster_columns = T, col = color_used, show_row_names = T, show_column_names = T,
  )

  pdf(file.path(output_dir, "S3H.pdf"), width = 6, height = 3)
  draw(p)
  dev.off()
}

# S3I
plot_S3I <- function(inputs, output_dir = ".") {

  {
    t1 <- inputs$M6_SPP1$t1
    prtb_ntc_df <- inputs$M6_SPP1$prtb_ntc_df
    plot_eg_df <- inputs$M6_SPP1$plot_eg_df
    x_axis <- 1.4
    xlab <- "Scaled perturbation effects on MP6_SPP1"
    p <- ggplot() + geom_point(
      data = t1, mapping = aes(x = effects, y = guide_group, color = max),
      size = 1.2, alpha = 0.7
    ) +
      geom_line(
        data = t1, mapping = aes(x = effects, y = guide_group, color = max),
        size = 0.9, alpha = 0.5
      ) +
      scale_color_gradientn(
        colours = c("#15638a", "#4f94b4", "#dfdfdf", "#e69877", "#c76a43"),
        values = scales::rescale(
          c(min(-max(t1$max), min(t1$max)), -0.5, 0, 0.5, max(t1$max))
        )
      ) +
      theme_linedraw(base_size = 6.5, base_family = "Arial") +
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
      theme(
        axis.text.y = element_text(size = 6),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#eeeeee", size = 0.02),
        legend.position = "none", plot.margin = margin(0, 0, 0, 0),
        panel.border = element_rect(size = 0.4, color = "black"),
        panel.spacing = unit(0, "pt")
      ) +
      xlim(-1.4, 1.4) +
      labs(x = xlab, y = "")
    p1 = prtb_ntc_df %>%
      ggplot(.) +
      geom_density(aes(x = effects), fill = "grey", linewidth = 0.2, bw = 0.07) +
      xlim(c(-x_axis, x_axis)) +
      theme_linedraw(base_size = 6.5, base_family = "Arial") +
      theme(plot.margin = margin(0, 0, 0, 0)) +
      theme(
        panel.grid = element_blank(), panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      )
    p1.1 <- prtb_ntc_df %>%
      ggplot(., aes(x = effects)) +
      geom_vline(aes(xintercept = effects), alpha = 0.5, color = "grey", size = 0.3) +
      xlim(c(-x_axis, x_axis)) +
      theme_linedraw(base_size = 6.5, base_family = "Arial") +
      theme(plot.margin = margin(0, 0, 0, 0)) +
      theme(
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.border = element_rect(size = 0.4,
        color = "black"),
        panel.spacing = unit(0, "pt")
      )
    p2 = plot_eg_df %>%
      ggplot(., aes(x = effects)) +
      geom_vline(
        aes(xintercept = effects, color = line_color), size = 0.75, alpha = 0.7
      ) +
      xlim(c(-x_axis, x_axis)) +
      theme_bw(base_size = 6.5, base_family = "Arial") +
      theme(panel.background = element_rect(color = "black")) +
      facet_wrap(. ~ guide_group, ncol = 1, strip.position = "left") +
      theme(plot.margin = margin(r = 10), panel.spacing = unit(0, "pt")) +
      theme(
        axis.text = element_text(color = "black"),
        strip.background = element_blank(), strip.placement = "outside",
        strip.text.y.left = element_text(size = 6, family = "Arial", angle = 0,
        color = "black"),
        panel.grid = element_blank(), panel.border = element_rect(size = 0.4,
        color = "black"),
        panel.spacing = unit(0.1, "lines"), legend.position = "none"
      ) +
      xlab("") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_color_manual(values = c("#e47e52", "#dfdfdf", "#1b536e"))
    pc = p1 / p1.1 / p2 / p + plot_layout(heights = c(1.9, 0.5, 3.4, 9)) +
      plot_layout(guides = "collect")

    ggsave(file.path(output_dir, "S3I_M6_SPP1.pdf"), pc, width = 2.3, height = 3.2)
  }

  {
    t1 <- inputs$M1_ISG$t1
    prtb_ntc_df <- inputs$M1_ISG$prtb_ntc_df
    plot_eg_df <- inputs$M1_ISG$plot_eg_df
    x_axis <- 1.4
    xlab <- "Scaled perturbation effects on MP1_ISG"
    p <- ggplot() + geom_point(
      data = t1, mapping = aes(x = effects, y = guide_group, color = max),
      size = 1.2, alpha = 0.7
    ) +
      geom_line(
        data = t1, mapping = aes(x = effects, y = guide_group, color = max),
        size = 0.9, alpha = 0.5
      ) +
      scale_color_gradientn(
        colours = c("#dfdfdf", "#e69877", "#c76a43"),
        values = scales::rescale(c(0, 0.2, max(t1$max)))
      ) +
      theme_linedraw(base_size = 6.5, base_family = "Arial") +
      geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
      theme(
        axis.text.y = element_text(size = 6),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "#eeeeee", size = 0.02),
        legend.position = "none", plot.margin = margin(0, 0, 0, 0),
        panel.border = element_rect(size = 0.4, color = "black"),
        panel.spacing = unit(0, "pt")
      ) +
      xlim(-1.4, 1.4) +
      labs(x = xlab, y = "")
    p1 = prtb_ntc_df %>%
      ggplot(.) +
      geom_density(aes(x = effects), fill = "grey", linewidth = 0.2, bw = 0.07) +
      xlim(c(-x_axis, x_axis)) +
      theme_linedraw(base_size = 6.5, base_family = "Arial") +
      theme(plot.margin = margin(0, 0, 0, 0)) +
      theme(
        panel.grid = element_blank(), panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      )
    p1.1 <- prtb_ntc_df %>%
      ggplot(., aes(x = effects)) +
      geom_vline(aes(xintercept = effects), alpha = 0.5, color = "grey", size = 0.3) +
      xlim(c(-x_axis, x_axis)) +
      theme_linedraw(base_size = 6.5, base_family = "Arial") +
      theme(plot.margin = margin(0, 0, 0, 0)) +
      theme(
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(), panel.border = element_rect(size = 0.4,
        color = "black"),
        panel.spacing = unit(0, "pt")
      )
    p2 = plot_eg_df %>%
      ggplot(., aes(x = effects)) +
      geom_vline(
        aes(xintercept = effects, color = line_color), size = 0.75, alpha = 0.7
      ) +
      xlim(c(-x_axis, x_axis)) +
      theme_bw(base_size = 6.5, base_family = "Arial") +
      theme(panel.background = element_rect(color = "black")) +
      facet_wrap(. ~ guide_group, ncol = 1, strip.position = "left") +
      theme(plot.margin = margin(r = 10), panel.spacing = unit(0, "pt")) +
      theme(
        axis.text = element_text(color = "black"),
        strip.background = element_blank(), strip.placement = "outside",
        strip.text.y.left = element_text(size = 6, family = "Arial", angle = 0,
        color = "black"),
        panel.grid = element_blank(), panel.border = element_rect(size = 0.4,
        color = "black"),
        panel.spacing = unit(0.1, "lines"), legend.position = "none"
      ) +
      xlab("") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_color_manual(values = c("#ec936d", "#dfdfdf", "#1b536e"))
    pc = p1 / p1.1 / p2 / p + plot_layout(heights = c(1.9, 0.5, 3.4, 9)) +
      plot_layout(guides = "collect")

    ggsave(file.path(output_dir, "S3I_M1_ISG.pdf"), pc, width = 2.3, height = 3.2)
  }

}

# S3J
plot_S3J <- function(volcano_data, output_dir = ".") {
  draw_volcano <- function(volcano_df, out_gene, coef_cut_off, min) {
    fdr_cut_off = 5
    max = 1.1
    p1 = ggplot(volcano_df) +
      ggrastr::geom_point_rast(
        aes(x = coef, y = mlg10fdr, color = guide_cutoff_class), size = 0.7, alpha = 0.6
      ) +
      geom_vline(
        xintercept = coef_cut_off, linetype = "dashed", color = "#424141",
        alpha = 0.6, size = 0.2
      ) +
      geom_vline(
        xintercept = -coef_cut_off, linetype = "dashed", color = "#424141",
        alpha = 0.6, size = 0.2
      ) +
      geom_hline(
        yintercept = fdr_cut_off, linetype = "dashed", color = "#424141",
        alpha = 0.6, size = 0.2
      ) +
      ggrepel::geom_text_repel(
        seed = 1, data = filter(volcano_df, guide_cutoff_class == "effect"),
        mapping = aes(x = coef, y = mlg10fdr, label = guide_group),
        size = 6 / .pt, family = "Arial", segment.size = 0.1, min.segment.length = 1
      ) +
      scale_color_manual(values = c("#bb610d", "grey")) +
      theme_linedraw(base_size = 6.2, base_family = "Arial") +
      theme(panel.grid = element_blank(), legend.position = "none") +
      xlim(c(min, max)) +
      xlab(paste0("Scaled perturbation effects on MP", gsub("^M", "", out_gene))) +
      ylab("-log10(FDR)")
    print(p1)
  }
  p1 = draw_volcano(volcano_data$M6_SPP1, "M6_SPP1", 0.8, -2)
  p2 = draw_volcano(volcano_data$M1_ISG, "M1_ISG", 0.3, -0.5)
  ggsave(file.path(output_dir, "S3J.pdf"), p1 + p2, width = 3.5, height = 1.8)
}

# S3K
plot_S3K <- function(coef_2M, output_dir = ".") {
  p2 = ggplot(coef_2M, aes(y = asinh(ISG15 / 0.01), x = asinh(SPP1 / 0.01))) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "dark grey",
    linewidth = 0.3) +
    geom_vline(xintercept = -1, linetype = "dashed", color = "dark grey",
    linewidth = 0.3) +
    geom_point(aes(color = color), size = 0.5, alpha = 0.7) +
    scale_color_manual(values = c(darkred = "#7c2727", lightpink = "orange",
    black = "#424242")) +
    ggrepel::geom_text_repel(
      seed = 1, aes(label = gene, color = color),
      fontface = "bold", family = "Arial", size = 6 / .pt, segment.size = 0.2,
      min.segment.length = 0.5, point.padding = 0.05, max.overlaps = 20
    ) +
    ylim(c(-0.8, 2)) +
    xlim(c(-3.4, 0.5)) +
    ylab("Scaled effects on ISG15") +
    xlab("Scaled effects on SPP1") +
    ggtitle("") +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "S3K.pdf"), p2, width = 1.8, height = 2)
}

# S3L
plot_S3L <- function(comparisons, output_dir = ".") {
  comp_df = comparisons$M6_SPP1
  out_gene = "M6_SPP1"

  p1 <- ggplot(comp_df) +
    geom_point(aes(x = x, y = y, color = class), size = 2.5, alpha = 0.8) +
    ggrepel::geom_text_repel(
      seed = 1, data = filter(comp_df, class != "non-sig"),
      mapping = aes(x = x, y = y, label = gene), size = 3, max.overlaps = 50
    ) +
    geom_abline(intercept = -0.4, linetype = "dashed", alpha = 0.3) +
    geom_abline(intercept = 0.4, linetype = "dashed", alpha = 0.3) +
    geom_segment(
      aes(y = -0.4, yend = -0.4, x = -Inf, xend = -0.4), color = "darkgrey", alpha = 0.3
    ) +
    geom_segment(
      aes(x = -0.4, xend = -0.4, y = -Inf, yend = -0.4), color = "darkgrey", alpha = 0.3
    ) +
    scale_color_manual(
      values = c(Consistent = "#c55650", TAM = "#fcac3d", MDM = "#89a536",
      `non-sig` = "grey"),
      name = ""
    ) +
    theme_linedraw() + theme(panel.grid = element_blank()) +
    xlab(paste0("TAM perturb effects on MP", gsub("^M", "", out_gene))) +
    ylab(paste0("MDM perturb effects on MP", gsub("^M", "", out_gene)))

  comp_df = comparisons$M1_ISG
  out_gene = "M1_ISG"
  cut_off = 0.3

  p2 <- ggplot(comp_df) +
    geom_point(aes(x = x, y = y, color = class), size = 2.5, alpha = 0.8) +
    ggrepel::geom_text_repel(
      seed = 1, data = filter(comp_df, class != "non-sig"),
      mapping = aes(x = x, y = y, label = gene), size = 3, max.overlaps = 50
    ) +
    geom_segment(
      aes(y = cut_off, yend = cut_off, x = Inf, xend = cut_off),
      color = "darkgrey", alpha = 0.3
    ) +
    geom_segment(
      aes(x = cut_off, xend = cut_off, y = Inf, yend = cut_off),
      color = "darkgrey", alpha = 0.3
    ) +
    scale_color_manual(
      values = c(Consistent = "#c55650", TAM = "#fcac3d", MDM = "#89a536",
      `non-sig` = "grey"),
      name = ""
    ) +
    theme_linedraw() + theme(panel.grid = element_blank()) +
    xlab(paste0("TAM perturb effects on MP", gsub("^M", "", out_gene))) +
    ylab(paste0("MDM perturb effects on MP", gsub("^M", "", out_gene)))

  ggsave(
    file.path(output_dir, "S3L.pdf"), p1 + p2 + plot_layout(guides = "collect"),
    width = 9.2, height = 3.7
  )
}

# S3M
plot_S3M <- function(coef_CTSL, genes_l, output_dir = ".") {

  fs = 6
  p = ggplot(coef_CTSL, aes(x = MDM, y = TAM)) +
    ggrastr::geom_point_rast(aes(color = class), alpha = 0.4, size = 0.5, shape = 1) +
    geom_point(
      data = coef_CTSL[coef_CTSL$gene %in% genes_l, ], aes(color = class),
      size = 1, shape = 20
    ) +
    geom_hline(yintercept = 0, color = "darkgrey", size = 0.2) +
    geom_vline(xintercept = 0, color = "darkgrey", size = 0.2) +
    annotate(
      "segment", x = 0.06, y = 0, xend = 0, yend = 0.06, linetype = "dashed", size = 0.3
    ) +
    annotate(
      "segment", x = -0.06, y = 0, xend = 0, yend = -0.06, linetype = "dashed",
      size = 0.3
    ) +
    scale_color_manual(values = c("#3f8a7d", "grey", "#d96d38")) +
    ggrepel::geom_text_repel(
      seed = 1, data = coef_CTSL[coef_CTSL$gene %in% genes_l, ], size = fs / .pt,
      aes(label = coef_CTSL[coef_CTSL$gene %in% genes_l, ]$gene),
      segment.size = 0.2, segment.alpha = 0.4, min.segment.length = 1,
      max.overlaps = 100,
      force = 20, force_pull = 0.5
    ) +
    scale_x_continuous(limits = c(-0.15, 0.16)) +
    scale_y_continuous(limits = c(-0.165, 0.15)) +
    xlab("Perturbation effects of CTSL in MDM") +
    ylab("Perturbation effects of CTSL in TAM") +
    theme_bw() + theme(
      panel.grid = element_blank(), axis.ticks = element_line(size = 0.2),
      axis.text = element_text(size = fs, color = "black"),
      axis.title = element_text(size = fs), legend.position = "none"
    )

  ggsave(file.path(output_dir, "S3M.pdf"), p, width = 2.3, height = 1.9)
}
