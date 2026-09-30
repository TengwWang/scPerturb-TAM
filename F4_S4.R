library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggpubr)
library(Seurat)
if (requireNamespace("extrafont", quietly = TRUE)) extrafont::loadfonts(quiet = TRUE)

# F4A
plot_F4A <- function(mat, coef_module, guide_modules, output_dir = ".") {
  gene_names <- rownames(mat)
  module_colors <- setNames(
    c(
      "#889BCB", "#7FB0A3", "#b992a4", "#97BC68", "#E6B85C", "#E89A90", "#C49A6C",
      "#7C7C7C", "#A6D1C5"
    )[1:8], rownames(coef_module)[1:8]
  )
  top_anno <- HeatmapAnnotation(
    MP = rownames(coef_module)[1:8], col = list(MP = module_colors),
    show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(1.7,
    "mm")
  )
  group_levels <- unique(as.character(guide_modules$GuideGroup))
  group_colors <- setNames(
    c("#7f8ca0", "#295f61", "#8a8998", "#374940", "#bebfbf",
    "#8c6f67")[1:length(group_levels)],
    group_levels
  )
  left_anno <- rowAnnotation(
    Cluster = as.character(guide_modules$GuideGroup),
    col = list(Cluster = group_colors),
    annotation_legend_param = list(title = "Guide Group"),
    show_annotation_name = TRUE, show_legend = FALSE, simple_anno_size = unit(2, "mm")
  )
  gene_highlight <- c("EREG", "TIMP2", "CTSL", "CCL2", "CXCL8", "FCGRT", "TNFSF12",
  "MAP1LC3B")
  highlight_idx <- which(gene_names %in% gene_highlight)
  right_anno <- rowAnnotation(
    gene = anno_mark(
      at = highlight_idx, labels = gene_names[highlight_idx], side = "right",
      labels_gp = gpar(fontsize = 9), link_width = unit(0.5, "cm"),
      extend = unit(0.5, "cm")
    )
  )

  color_used <- colorRamp2(
    seq(-2, 2, length = 50),
    colorRampPalette(c("#204d6c", "#93a5cd", "white", "#ca3755", "#a52d45"))(50)
  )
  p <- Heatmap(
    asinh(mat / 0.01),
    name = "Scaled\nperturbation\neffects", col = color_used, cluster_rows = FALSE,
    cluster_columns = FALSE, show_row_names = FALSE, show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 10),
    top_annotation = top_anno, left_annotation = left_anno,
    right_annotation = right_anno,
    row_split = guide_modules[, 2], column_split = 1:8, gap = unit(1.4, "mm"),
    column_gap = unit(1.4, "mm"),
    column_title = NULL,
    heatmap_legend_param = list(title = "Scaled\nperturbation\neffects"),
    height = unit(nrow(mat) * 0.15, "mm")
  )

  pdf(file.path(output_dir, "F4A.pdf"), width = 4.5, height = 4.5)
  draw(p)
  dev.off()
}

# F4B
plot_F4B <- function(SPP1_df, output_dir = ".") {
  anno_colors_dark <- c("#7f8ca0", "#295f61", "#8a8998", "#374940", "#bebfbf",
  "#8c6f67")
  p = ggplot(
    SPP1_df, aes(x = factor(group, levels = rev(unique(group))), y = asinh(coef / 0.01))
  ) +
    geom_hline(yintercept = -0.4, linetype = "dashed", color = "#868686", size = 0.5) +
    geom_boxplot(
      aes(group = group, fill = factor(group)),
      color = "#515151", width = 0.3, outlier.size = -1, alpha = 0.3, size = 0.35
    ) +
    geom_jitter(width = 0.07, aes(color = factor(group)), size = 0.2, alpha = 0.6) +
    scale_fill_manual(values = anno_colors_dark) +
    scale_color_manual(values = anno_colors_dark) +
    theme_classic() + labs(x = "Cluster", y = "Scaled effects on MP6_SPP1") +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6),
      axis.text = element_text(color = "black"),
      legend.position = "none", axis.line = element_blank(),
    ) +
    coord_flip()

  pdf(file.path(output_dir, "F4B.pdf"), width = 2.1, height = 3.9)
  withr::with_seed(1, print(p))
  dev.off()
}

# F4C
plot_F4C <- function(cor.gene, gene_modules, output_dir = ".") {
  color_used <- colorRamp2(
    seq(-0.45, 0.45, length = 11), rev(RColorBrewer::brewer.pal(11, "RdBu"))
  )
  library(scico)
  library(ggsci)
  colors_d3 <- pal_igv("default")(23)[-c(4, 9)]
  anno_colors <- setNames(colors_d3, unique(gene_modules$GeneGroup))
  top_anno <- HeatmapAnnotation(
    clusters = as.character(gene_modules$GeneGroup), col = list(clusters = anno_colors),
    simple_anno_size = unit(2.5, "mm"), show_annotation_name = FALSE, show_legend = F
  )
  gene_highlight <- c(
    "CXCL10", "IFI44", "IFIT1", "IFIT2", "ISG15", "ISG20", "CD74", "HLA-DPA1",
    "HLA-DPB1", "HLA-DQA1", "HLA-DRA", "MALAT1", "NEAT1", "ZEB2", "FOXN3", "AHRR",
    "CCL2", "CCL7", "CCL8", "S100A9", "GPNMB", "CTSB", "CTSD"
  )
  gene_names = gene_modules$GeneName
  highlight_idx <- which(gene_names %in% gene_highlight)
  right_anno <- rowAnnotation(
    gene = anno_mark(
      at = highlight_idx, labels = gene_names[highlight_idx], side = "right",
      labels_gp = gpar(fontsize = 9.5), link_width = unit(0.5, "cm"),
      extend = unit(2, "cm")
    )
  )
  p1 <- Heatmap(
    cor.gene[gene_names, gene_names], name = "Gene Correlation", col = color_used,
    cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F,
    column_names_gp = gpar(fontsize = 5),
    use_raster = F, row_names_gp = gpar(fontsize = 5),
    top_annotation = top_anno, right_annotation = right_anno, gap = unit(0.2, "mm"),
    column_gap = unit(0.2, "mm"),
    row_split = gene_modules$GeneGroup, column_split = gene_modules$GeneGroup,
    row_title = NULL, column_title = NULL, heatmap_width = unit(12, "cm"),
    heatmap_height = unit(10, "cm"),
    layer_fun = function(j, i, x, y, width, height, fill, slice_r, slice_c) {
      if (slice_r == slice_c) {
        grid.rect(gp = gpar(lwd = 1.7, col = "black", fill = NA))
      }
    }
  )
  pdf(file.path(output_dir, "F4C.pdf"), width = 6.5, height = 4.5)
  draw(p1)
  dev.off()
}

# F4D
plot_F4D <- function(rank_4gp, output_dir = ".") {

  slct_row <- rank_4gp[rank_4gp$target %in% c("CTSL", "FCGRT", "BRD2", "TNFSF12",
  "CCL4", "CD59", "COX20", "TREM2"),
  ]
  p = ggplot(data = rank_4gp) +
    geom_point(
      aes(x = reorder(target, rank_sum), y = score),
      color = "#bfbfbf", size = 0.2, alpha = 0.6
    ) +
    geom_point(
      data = slct_row, aes(x = target, y = score), color = "#c6715a", size = 1
    ) +
    ggrepel::geom_text_repel(
      seed = 1, data = slct_row, point.padding = 0.5, min.segment.length = 1,
      force = 10, aes(x = target, y = score, label = target),
      color = "#c6715a", size = 3, fontface = "bold", max.overlaps = 10
    ) +
    ggthemes::theme_few() + scale_x_discrete(expand = expansion(add = 25)) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10), plot.title = element_text(size = 10),
      legend.position = "none",
    ) +
    xlab("Targets") +
    ylab("Pro-tumor target score")

  pdf(file.path(output_dir, "F4D.pdf"), width = 2.5, height = 4)
  withr::with_seed(1, print(p))
  dev.off()
}

# F4E
plot_F4E <- function(go_rank1_asinh, go_rank2_asinh, output_dir = ".") {
  draw_rank_list = function(go_rank, CTSL_col) {
    p_rank_go_list = list()
    for (n in names(go_rank)) {
      top_rank = go_rank[[n]]
      top_rank$color[top_rank$target == "CTSL"] = "B"
      ctsl_row <- top_rank[top_rank$target == "CTSL", ]
      other_rows <- top_rank[top_rank$target != "CTSL", ]
      p_rank_go_list[[n]] <- ggplot(data = top_rank) +
        ggrastr::rasterise(
          geom_segment(
            aes(
              x = reorder(target, rank), xend = reorder(target, rank),
              y = 0, yend = effects
            ), colour = "lightgrey", linewidth = 0.08, alpha = 0.3
          ), dpi = 300
        ) +
        ggrastr::geom_point_rast(
          aes(x = reorder(target, rank), y = effects),
          color = "#cccccc", size = 0.2, alpha = 0.6
        ) +
        geom_segment(
          data = ctsl_row, aes(x = target, xend = target, y = 0, yend = effects),
          colour = CTSL_col, linewidth = 0.6, alpha = 0.8
        ) +
        geom_point(
          data = ctsl_row, aes(x = target, y = effects), color = CTSL_col, size = 1
        ) +
        geom_hline(yintercept = 0, color = "grey") +
        ggrepel::geom_text_repel(
          seed = 1, data = ctsl_row, point.padding = 20, force = 3, nudge_y = 0.001,
          aes(x = target, y = effects, label = paste0("CTSL (rank: ", rank, ")")),
          color = CTSL_col, size = 3.6, fontface = "bold"
        ) +
        theme_bw() + scale_x_discrete(expand = expansion(add = 25)) +
        theme(
          axis.text = element_text(color = "black", size = 8),
          panel.grid = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          plot.title = element_text(size = 10, hjust = 0.5, vjust = 0.2),
          legend.position = "none"
        ) +
        xlab("Targets") +
        ylab(paste0("pertrurb effects")) +
        ggtitle(stringr::str_wrap(paste0("GO: ", n), width = 24))
    }
    return(p_rank_go_list)
  }
  p_rank_go_list1 = draw_rank_list(go_rank1_asinh, "#ac826a")
  p_rank_go_list2 = draw_rank_list(go_rank2_asinh, "#5e887a")
  p_rank_go_list = c(p_rank_go_list1, p_rank_go_list2)
  p_rank_go_list[[1]] = p_rank_go_list[[1]] + ylab("Mean scaled effects (Donor 1)")
  p_rank_go_list[[4]] = p_rank_go_list[[4]] + ylab("Mean scaled effects (Donor 2)")
  for (n in c(2, 3, 5, 6)) p_rank_go_list[[n]] = p_rank_go_list[[n]] +
    ylab("") +
    theme(plot.margin = margin(l = 0))
  for (n in c(1:3)) p_rank_go_list[[n]] = p_rank_go_list[[n]] +
    theme(axis.title.x = element_blank()) +
    theme(plot.margin = margin(b = 0))
  for (n in c(4:6)) p_rank_go_list[[n]] = p_rank_go_list[[n]] +
    ggtitle("") +
    theme(plot.margin = margin(t = 0))

  p = wrap_plots(p_rank_go_list, ncol = 3)

  pdf(file.path(output_dir, "F4E.pdf"), width = 6, height = 4.7)
  withr::with_seed(1, print(p))
  dev.off()
}

# F4F
plot_F4F <- function(coef_2M, genes, output_dir = ".") {

  p = ggplot(coef_2M, aes(y = M1_ISG, x = M6_SPP1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "dark grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "dark grey") +
    ggrastr::geom_point_rast(color = "#b3b3b3", alpha = 1, size = 0.4) +
    geom_point(
      data = coef_2M[coef_2M$gene %in% genes, ], aes(y = M1_ISG, x = M6_SPP1),
      color = "#C02A1B", size = 1
    ) +
    ggrepel::geom_text_repel(
      seed = 1, data = coef_2M[coef_2M$gene %in% genes, ], aes(label = gene),
      color = "#C02A1B", max.overlaps = 50, force = 60, force_pull = 5,
      point.padding = 0.001,
      min.segment.length = 1.4, size = 3.2, fontface = "bold", segment.size = 0.5,
      segment.alpha = 0.4
    ) +
    ylab("Mean scaled effects on MP1_ISG") +
    xlab("Mean scaled effects on MP6_SPP1") +
    theme_bw(base_size = 11.5) +
    theme(
      axis.text = element_text(color = "black"),
      panel.grid = element_blank(), legend.position = "none"
    )

  pdf(file.path(output_dir, "F4F.pdf"), width = 4.5, height = 4.2)
  withr::with_seed(1, print(p))
  dev.off()
}

# F4G
plot_F4G <- function(
    tmp, all_genes, angiogenesis, immunosuppressive, inflammation, antigen_presentation,
    T_activation, colors, breaks, output_dir = "."
    ) {
  gene_categories <- c(
    rep("Angiogenesis", length(angiogenesis)),
    rep("Immunosuppressive", length(immunosuppressive)),
    rep("Inflammation", length(inflammation)),
    rep("Antigen presentation", length(antigen_presentation)),
    rep("T cell activation", length(T_activation))
  )
  annotation_row <- data.frame(Category = gene_categories)
  rownames(annotation_row) <- all_genes
  annotation_colors <- list(
    Category = c(
      Angiogenesis = "#b8a9c9", Inflammation = "#a8d6d6", Immunosuppressive = "#d6c1ab",
      `Antigen presentation` = "#a1b5d8", `T cell activation` = "#e4aa90"
    )
  )
  p = pheatmap::pheatmap(
    silent = TRUE, asinh(tmp / 0.01),
    cluster_rows = F, cluster_cols = F, scale = "none", annotation_row = annotation_row,
    annotation_colors = annotation_colors, show_rownames = T, show_colnames = T,
    fontsize = 9, color = colors, breaks = breaks, cellwidth = 22, cellheight = 9,
    border_color = NA, angle_col = 45
  )

  pdf(file.path(output_dir, "F4G.pdf"), width = 5, height = 6)
  withr::with_seed(1, print(p))
  dev.off()
}

# F4H
plot_F4H <- function(resM_CTSL, genelists, output_dir = ".") {
  Volcano <- function(deg, title, l1 = "up", l2 = "down", label = NULL) {
    deg <- deg[, c("log2FoldChange", "padj", "sig", "gene")] %>%
      na.omit()
    FC_t <- 1
    p_adj_t <- 0.05
    k1 <- (deg$padj < p_adj_t) & (deg$log2FoldChange < -FC_t)
    k2 <- (deg$padj < p_adj_t) & (deg$log2FoldChange > FC_t)
    if (!"sig" %in% colnames(deg)) {
      sig <- ifelse(k1, l2, ifelse(k2, l1, "non-sig"))
      deg$sig <- factor(sig, levels = c(l2, "non-sig", l1))
    }
    if (is.null(label)) {
      label <- deg %>%
        filter(k1) %>%
        arrange(log2FoldChange) %>%
        head(15) %>%
        rbind(deg %>% filter(k2) %>% arrange(log2FoldChange) %>% tail(15)) %>%
        rbind(deg %>% filter(k1) %>% arrange(padj) %>% head(10)) %>%
        rbind(deg %>% filter(k2) %>% arrange(padj) %>% head(10)) %>%
        unique.data.frame()
    } else label <- deg[deg$gene %in% label, ]
    label$sig <- ifelse(
      label$log2FoldChange < 0, paste0(l2, "_label"), paste0(l1, "_label")
    )
    p <- ggplot(deg, aes(log2FoldChange, -log10(padj + 1e-300))) +
      ggrastr::geom_point_rast(aes(color = sig), alpha = 0.5, size = 0.3) +
      geom_point(data = label, aes(color = sig), size = 0.4) +
      ggrepel::geom_text_repel(
        seed = 1, data = label, aes(label = gene),
        size = 6.5 / .pt, color = "black", min.segment.length = 1, point.padding = 0.05,
        force = 20, force_pull = 2, max.overlaps = 200, segment.size = 0.2,
        segment.alpha = 0.5
      ) +
      scale_color_manual(
        values = c(
          up = "#d88984", up_label = "#ce2b2b", `non-sig` = "#cfcdcd", down = "#8abddc",
          down_label = "#2863a3"
        )
      ) +
      labs(
        x = expression(paste(Log[2], " Fold Change")),
        y = expression(paste(-Log[10], " Padj"))
      ) +
      theme_classic(base_size = 7) +
      theme(
        legend.position = "none", axis.title = element_text(size = 7),
        axis.text = element_text(size = 6.5, color = "black"),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 0.4),
        axis.line = element_blank()
      )
    return(p)
  }

  p <- Volcano(resM_CTSL, "M_CTSL_KO", label = genelists)

  pdf(file.path(output_dir, "F4H.pdf"), width = 3, height = 3)
  withr::with_seed(1, print(p))
  dev.off()
}

# F4I
plot_F4I <- function(tpm_gsva_type, tpm, output_dir = ".") {

  ann_col = data.frame(group = rep(c("AAVS1KO", "CTSLKO"), each = 3))
  rownames(ann_col) = colnames(tpm)
  colnames(tpm_gsva_type) = paste0(
    rep(c("AAVS1KO", "CTSLKO"), each = 3), rep(c("-1", "-2", "-3"), times = 2)
  )
  group_col = c("#9d9d9d", "#b47a7a")
  names(group_col) = c("AAVS1KO", "CTSLKO")
  fs = 6.5
  hc_row = hclust(dist(tpm_gsva_type), method = "complete")
  ha = HeatmapAnnotation(
    group = ann_col$group, col = list(group = group_col),
    show_annotation_name = F, simple_anno_size = unit(2, "mm"),
    annotation_legend_param = list(
      title = NULL, direction = "horizontal", grid_height = unit(3, "mm"),
      grid_width = unit(3, "mm"), labels_gp = gpar(fontsize = fs)
    )
  )
  col_fun = colorRampPalette(c("#A3B7DC", "#F5F8FA", "#FFA899"))(100)
  ht = Heatmap(
    tpm_gsva_type, name = "GSVA score", col = col_fun, cluster_columns = FALSE,
    cluster_rows = as.dendrogram(hc_row),
    show_column_names = T, show_row_names = TRUE, column_names_rot = 45,
    row_names_gp = gpar(fontsize = fs),
    column_names_gp = gpar(fontsize = fs), row_dend_gp = gpar(lwd = 0.5),
    width = unit(3, "cm"), height = unit(4, "cm"),
    column_split = ann_col$group, column_title = NULL, column_gap = unit(0.7, "mm"),
    top_annotation = ha, heatmap_legend_param = list(
      direction = "horizontal", labels_gp = gpar(fontsize = fs), at = c(-1, 0, 1),
      labels = c("-1", "0", "1"), grid_height = unit(2.5, "mm"),
      title_gp = gpar(fontsize = fs), legend_width = unit(1, "cm")
    )
  )

  pdf(file.path(output_dir, "F4I.pdf"), width = 2.5, height = 3.2)
  draw(
    ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
    merge_legends = TRUE
  )
  dev.off()
}

# S4B
plot_S4B <- function(hist_df1, hist_df2, output_dir = ".") {
  # Source: s4_1374lib_QC, cell 15 (zero-based)

  plot_dist_for_guide_per_cell1 <- function(hist_df, column,
  xlab = "Number of sgRNA categories per cell") {
    hist_df_mod <- hist_df
    m = 20
    hist_df_mod[hist_df_mod[, column] > m, column] <- m
    breaks = seq(-0.5, m + 1, 1)
    p1 <- ggplot(data = hist_df_mod) +
      geom_histogram(
        breaks = breaks, alpha = 0.5, linewidth = 0.3, aes(x = hist_df_mod[, column],
        y = stat(count) / sum(count)),
        fill = "#becac0", color = "#2c2c2c"
      ) +
      ylab("Relative frequency") +
      xlab(xlab) +
      theme_linedraw(base_size = 7, base_family = "Arial") +
      theme(panel.grid = element_blank()) +
      scale_x_continuous(
        breaks = c(seq(0, m - 5, 5), m), labels = c(seq(0, m - 5, 5), paste0(">", m))
      )
    xpars <- sapply(c("mean", "median"), do.call, list(x = hist_df[, column])) %>%
      round(., 2) %>%
      paste(names(.), ., sep = " = ") %>%
      paste(., collapse = "\n")
    p1 <- p1 + geom_vline(
      xintercept = median(hist_df[, column]), linetype = "dashed", size = 0.3
    ) +
      annotate(geom = "text", x = 10, y = 0.1, label = xpars, hjust = 0, size = 7 / .pt)

  }
  p1 = plot_dist_for_guide_per_cell1(hist_df1, "sgnoflt")

  # Source: s4_1374lib_QC, cell 16 (zero-based)

  p2 = plot_dist_for_guide_per_cell1(hist_df2, "sgnoflt")

  pdf(file.path(output_dir, "S4B.pdf"), width = 2, height = 2.8)
  withr::with_seed(1, print(p1 / p2))
  dev.off()
}

# S4C
plot_S4C <- function(ef_abs_df_long1, ef_abs_df_long2, output_dir = ".") {
  p1 <- ggplot(ef_abs_df_long1, aes(x = variable, y = asinh(value / 0.01))) +
    geom_boxplot(color = "grey40", width = 0.6, alpha = 0.9, outlier.size = -1,
    size = 0.35) +
    ggrastr::geom_jitter_rast(
      width = 0.1, aes(color = variable), size = 0.1, alpha = 0.4
    ) +
    scale_color_brewer(palette = "Set2") +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1,
      size = 6)
    ) +
    scale_x_discrete(
      labels = c(target_ef_abs = "Target", ntc_ef_abs = "NTC", aavs1_ef_abs = "AAVS1")
    ) +
    labs(y = "Absolute value of scaled effects", x = "") +
    ggpubr::stat_compare_means(
      comparisons = list(
        c("target_ef_abs", "ntc_ef_abs"), c("target_ef_abs", "aavs1_ef_abs")
      ), method = "wilcox.test", size = 6 / .pt, bracket.size = 0.2
    )
  p2 <- ggplot(ef_abs_df_long2, aes(x = variable, y = asinh(value / 0.01))) +
    geom_boxplot(color = "grey40", width = 0.6, alpha = 0.9, outlier.size = -1,
    size = 0.35) +
    ggrastr::geom_jitter_rast(
      width = 0.1, aes(color = variable), size = 0.1, alpha = 0.4
    ) +
    scale_color_brewer(palette = "Set2") +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1,
      size = 6)
    ) +
    scale_x_discrete(
      labels = c(target_ef_abs = "Target", ntc_ef_abs = "NTC", aavs1_ef_abs = "AAVS1")
    ) +
    labs(y = "Absolute value of scaled effects", x = "") +
    ggpubr::stat_compare_means(
      comparisons = list(
        c("target_ef_abs", "ntc_ef_abs"), c("target_ef_abs", "aavs1_ef_abs")
      ), method = "wilcox.test", size = 6 / .pt, bracket.size = 0.2
    )
  pdf(file.path(output_dir, "S4C.pdf"), width = 2.6, height = 2.5)
  withr::with_seed(1, print(p1 + p2))
  dev.off()
}

# S4D
plot_S4D <- function(percent_data, output_dir = ".") {
  draw_percent <- function(data) {
    list2env(data, envir = environment())
    library(scales)
    p2 = ggplot(all_ef_prct_df, aes(x = percent)) +
      geom_density(fill = "#5f85cc", alpha = 0.4, linewidth = 0.3) +
      geom_vline(
        aes(xintercept = ntc_prct), color = "#F28C28", linetype = "solid", size = 0.6
      ) +
      geom_vline(
        aes(xintercept = aavs1_prct), color = "#C3A6FF", linetype = "solid", size = 0.6
      ) +
      annotate(
        "text", x = 0.06, y = 1, label = "NTC", color = "#F28C28", hjust = 1,
        size = 6 / .pt, fontface = "bold"
      ) +
      annotate(
        "text", x = 0.13, y = 1, label = "AAVS1", color = "#C3A6FF", hjust = 0,
        size = 6 / .pt, fontface = "bold"
      ) +
      annotate(
        "text", x = 0.5, y = 7, label = "Targets", color = "#5f85cc", size = 6 / .pt,
        fontface = "bold", hjust = 1
      ) +
      labs(title = "", x = "Percentage of affected genes (|coef|>0.005)",
      y = "Density") +
      scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.7)) +
      theme_linedraw(base_size = 7, base_family = "Arial") +
      theme(panel.grid = element_blank()) +
      theme(legend.position = "none")

  }
  p1 = draw_percent(percent_data$TAM1)
  p2 = draw_percent(percent_data$TAM2)
  pdf(file.path(output_dir, "S4D.pdf"), width = 2.6, height = 2.2)
  withr::with_seed(1, print(p1 / p2))
  dev.off()
}

# S4E
plot_S4E <- function(df2, df3, output_dir = ".") {
  plot_babble_asinh = function(df, targetname) {
    df$Pval = 10^(-abs(df$Pval))
    df$Regulation <- ifelse(df$Coef > 0, "Upregulated", "Downregulated")
    df <- df %>%
      arrange(Coef)
    ggplot(
      df, aes(
        x = reorder(Gene, Coef), y = asinh(Coef / 0.01), size = -log10(Pval),
        color = Regulation
      )
    ) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c(Upregulated = "firebrick",
      Downregulated = "steelblue")) +
      scale_size_continuous(
        limits = c(0, 90), range = c(1.5, 2.8), breaks = c(30, 60, 90)
      ) +
      theme_linedraw(base_size = 7, base_family = "Arial") +
      theme(panel.grid = element_blank()) +
      labs(y = paste0("Scaled perturbation effects of ", targetname), x = "") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5),
        legend.position = "right", legend.key.size = unit(0.2, "cm")
      ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray")
  }
  p = plot_babble_asinh(df2, "STAT6") +
    plot_babble_asinh(df3, "ID3") +
    patchwork::plot_layout(widths = c(3, 2))
  pdf(file.path(output_dir, "S4E.pdf"), width = 4.5, height = 2)
  withr::with_seed(1, print(p))
  dev.off()
}

# S4F
plot_S4F <- function(IntGraph, output_dir = ".") {
  # Source: s4_1374lib_QC, cell 48 (zero-based)
  library(ggpubr)
  library(ggplot2)
  p1 = ggplot(IntGraph, aes(x = CorValue, y = combined_score)) +
    ggrastr::geom_point_rast(color = "#fc8772", size = 0.4, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#d85a44",
    size = 0.6) +
    stat_cor(
      aes(label = paste(..r.label.., ..p.label.., sep = "~~~")),
      label.x = -0.2, label.y = 1200, color = "black", size = 7 / .pt
    ) +
    labs(x = "Correlation on perturbation profile", y = "STRING Score") +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank())

  # Source: s4_1374lib_QC, cell 49 (zero-based)
  IntGraph$CorGroup <- cut(
    IntGraph$CorValue, breaks = c(-Inf, 0.2, 0.5, Inf),
    labels = c("< 0.2", "0.2 - 0.5", "> 0.5")
  )
  p2 = ggplot(IntGraph, aes(x = CorGroup, y = combined_score, fill = CorGroup)) +
    geom_violin(trim = FALSE, alpha = 0.8, linewidth = 0.4, width = 1.1) +
    labs(x = "Correlation on perturbation profile", y = "STRING Score") +
    theme_minimal() + scale_fill_brewer(palette = "Set3") +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) +
    theme(legend.position = "none")

  pdf(file.path(output_dir, "S4F.pdf"), width = 3.8, height = 1.8)
  withr::with_seed(1, print(p1 + p2))
  dev.off()
}

# S4G
plot_S4G <- function(coef_consis, output_dir = ".") {
  library(ggpubr)

  library(patchwork)
  p1 <- ggplot(coef_consis, aes(x = G8_ISG, y = M1_ISG)) +
    geom_point(alpha = 0.6, color = "#7C9A92", size = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "#A35D5D", linetype = "22",
    linewidth = 0.5) +
    stat_cor(
      aes(label = after_stat(r.label)),
      label.x.npc = "left", label.y.npc = "top", size = 7 / .pt
    ) +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) +
    labs(x = "GP8_ISG", y = "MP1_ISG")
  p2 <- ggplot(coef_consis, aes(x = G19_inflam, y = M6_SPP1)) +
    geom_point(alpha = 0.6, color = "#7C9A92", size = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "#A35D5D", linetype = "22",
    linewidth = 0.5) +
    stat_cor(
      aes(label = after_stat(r.label)),
      label.x.npc = "left", label.y.npc = "top", size = 7 / .pt
    ) +
    theme_linedraw(base_size = 7, base_family = "Arial") +
    theme(panel.grid = element_blank()) +
    labs(x = "GP19_protumor_chemo", y = "MP6_SPP1")
  p1 + p2
  pdf(file.path(output_dir, "S4G.pdf"), width = 3.5, height = 1.8)
  withr::with_seed(1, print(p1 + p2))
  dev.off()
}

# S4H
plot_S4H <- function(df_long, output_dir = ".") {
  p <- ggplot(df_long, aes(x = id, y = celltype)) +
    geom_point(aes(size = `Percent expressed`, color = `Average Expression`)) +
    scale_size(range = c(1, 7)) +
    scale_color_gradientn(colours = c("#F9D4D0", "#F4A6A0", "#EA6A63", "#D73027",
    "#67000D")) +
    theme_classic() + theme(
      axis.title = element_blank(), axis.text = element_text(color = "black",
      size = 12),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      axis.line = element_blank(), legend.key.size = unit(16, "pt"),
      legend.title = element_text(size = 12), legend.text = element_text(size = 12)
    )

  pdf(file.path(output_dir, "S4H.pdf"), width = 2, height = 3.5)
  withr::with_seed(1, print(p))
  dev.off()
}

# S4I
plot_S4I <- function(panMye, output_dir = ".") {

  genelist = c("CTSL", "HIF1A", "TREM2", "PIK3CA", "PIK3CG")
  panMye$MajorCluster[panMye$MajorCluster == "Myeloid"] = "Other"
  p = DotPlot(panMye, genelist, group.by = "MajorCluster") +
    labs(x = NULL, y = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.85, hjust = 0.75),
      legend.key.size = unit(16, "pt"),
      axis.text = element_text(color = "black", size = 12),
      legend.title = element_text(size = 12), legend.text = element_text(size = 12)
    ) +
    scale_color_gradientn(colours = c("#f7f6f6", "#b5c7d4", "#46819F"))

  pdf(file.path(output_dir, "S4I.pdf"), width = 3.5, height = 3.5)
  withr::with_seed(1, print(p))
  dev.off()
}

# S4J
plot_S4J <- function(box.data.melt_ctsl, downsampled_data, fill_colors,
output_dir = ".") {
  box.data.melt_ctsl %>%
    ggplot(aes(ClusterName, Expr)) +
    stat_boxplot(
      geom = "errorbar", position = position_dodge(width = 0.2), width = 0.3
    ) +
    geom_boxplot(
      position = position_dodge(width = 0.2), width = 0.6, outlier.shape = NA
    ) +
    geom_point(
      data = downsampled_data, aes(fill = ClusterName, color = ClusterName),
      size = 0.2, alpha = 0.4, position = position_jitterdodge(jitter.width = 0.3,
      dodge.width = 0.2)
    ) +
    scale_fill_manual(values = fill_colors) +
    scale_color_manual(values = fill_colors) +
    labs(x = NULL, y = "CTSL Expression") +
    scale_x_discrete(limits = rev(levels(factor(box.data.melt_ctsl$ClusterName)))) +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
      axis.line = element_line(color = "black", size = 0.4),
      panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
      panel.background = element_blank(), axis.text.y = element_text(color = "black",
      size = 12),
      axis.text.x = element_text(color = "black", size = 12, angle = 45, vjust = 1,
      hjust = 1),
      axis.line.x.top = element_line(color = "black"),
      axis.text.x.top = element_blank(), axis.ticks.y.right = element_blank(),
      axis.text.y.right = element_blank(), axis.ticks.x.top = element_blank(),
      panel.spacing.x = unit(0, "cm"),
      panel.border = element_blank(), legend.position = "none",
      panel.spacing = unit(0, "lines")
    ) +
    guides(x.sec = "axis", y.sec = "axis") ->
  p

  pdf(file.path(output_dir, "S4J.pdf"), width = 6, height = 4)
  withr::with_seed(1, print(p))
  dev.off()
}

# S4M
plot_S4M <- function(tpm_slt, ann_col, ann_colors, output_dir = ".") {
  p = pheatmap::pheatmap(
    tpm_slt, silent = TRUE, show_rownames = TRUE, cluster_cols = FALSE,
    cluster_rows = FALSE,
    border = FALSE, scale = "row", annotation_col = ann_col,
    annotation_colors = ann_colors,
    gaps_col = 3, color = colorRampPalette(c("#A3B7DC", "#F5F8FA", "#FFA899"))(100),
    cellwidth = 20, cellheight = 20
  )
  pdf(file.path(output_dir, "S4M.pdf"), width = 5, height = 8)
  withr::with_seed(1, print(p))
  dev.off()
}

# S4N
plot_S4N <- function(resM0_CTSL, genelists, output_dir = ".") {
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
        rbind(deg %>% filter(k2) %>% arrange(log2FoldChange) %>% tail(15)) %>%
        rbind(deg %>% filter(k1) %>% arrange(padj) %>% head(10)) %>%
        rbind(deg %>% filter(k2) %>% arrange(padj) %>% head(10)) %>%
        unique.data.frame()
    } else {
      label = deg[deg$gene %in% label, ]
    }
    label$sig = ifelse(
      label$log2FoldChange < 0, paste0(l2, "_label"), paste0(l1, "_label")
    )
    p <- ggplot(data = deg, aes(x = log2FoldChange, y = -log10(padj + 1e-300))) +
      geom_point(alpha = 0.5, size = 1, aes(color = sig)) +
      ylab("-log10(Padj)") +
      scale_color_manual(
        values = c(
          up = "#d88984", up_label = "#ce2b2b", `non-sig` = "#cfcdcd", down = "#8abddc",
          down_label = "#2863a3"
        )
      ) +
      geom_point(size = 1.02, data = label, aes(color = sig)) +
      ggrepel::geom_text_repel(
        seed = 1, aes(label = label$gene),
        data = label, point.padding = 0.2, force = 10, max.overlaps = 20,
        color = "black", size = 4, label.size = 0.001, segment.size = 0.3,
        label.padding = 0.1, box.padding = 0.3
      ) +
      theme_classic() + labs(
        x = expression(paste(Log[2], " Fold Change")),
        y = expression(paste(-Log[10], "Padj"))
      ) +
      theme(legend.position = "none") +
      theme(
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1.3),
        axis.line = element_blank()
      )
  }
  p = Volcano(resM0_CTSL, "M0_CTSL_KO", label = genelists)

  pdf(file.path(output_dir, "S4N.pdf"), width = 5.1, height = 5.1)
  withr::with_seed(1, print(p))
  dev.off()
}

# S4O
plot_S4O <- function(gsego1, output_dir = ".") {
  p <- gsego1 %>%
    as_tibble() %>%
    arrange(desc(NES)) %>%
    filter(p.adjust < 0.05) %>%
    mutate(Description = sub("^([a-z])", "\\U\\1", Description, perl = TRUE)) %>%
    ggplot(aes(stats::reorder(Description, NES), NES)) +
    geom_col(aes(fill = NES), width = 0.6, alpha = 0.8) +
    theme_classic() + ggplot2:::coord_flip() + scale_fill_gradient2(low = "#8abddc",
    mid = "white", high = "#d88984") +
    labs(x = "", y = "Normalized enrichment score") +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1.2),
      axis.line = element_blank(), axis.text = element_text(color = "black")
    ) +
    theme(
      legend.key.size = unit(0.3, "cm"), legend.text = element_text(size = 8),
      legend.title = element_text(size = 8.5)
    )

  pdf(file.path(output_dir, "S4O.pdf"), width = 5.2, height = 3)
  withr::with_seed(1, print(p))
  dev.off()
}
