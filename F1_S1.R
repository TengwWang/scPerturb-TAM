
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggpubr)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(survminer)
})

# F1B
plot_F1B <- function(plot_mat, module2_vec, file = "f1_NMF_cluster_heatmap.pdf") {
  color_used <- colorRamp2(
    seq(0, 15, length = 7), rev(RColorBrewer::brewer.pal(11, "RdBu")[1:7])
  )
  module2_vec <- factor(module2_vec)
  module2_levels <- unique(module2_vec)
  module2_colors <- setNames(
    c(
      "#889BCB", "#7FB0A3", "#b992a4", "#97BC68", "#E6B85C", "#E89A90",
      "#C49A6C", "#7C7C7C", "#A6D1C5"
    ), module2_levels
  )
  top_anno <- HeatmapAnnotation(
    MP = module2_vec, col = list(MP = module2_colors),
    show_annotation_name = F, show_legend = FALSE, simple_anno_size = unit(2, "mm")
  )
  ht <- Heatmap(
    plot_mat, col = color_used, name = "Similarity", top_annotation = top_anno,
    show_row_names = FALSE, show_column_names = FALSE, cluster_rows = FALSE,
    cluster_columns = FALSE, row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    row_title = "NMF consistent programs", column_title = "NMF consistent programs",
    row_title_side = "left", column_title_side = "bottom"
  )
  pdf(file, width = 5, height = 4.5)
  draw(ht, padding = unit(c(4, 4, 5, 4), "mm"))
  for (i in seq_along(unique(module2_vec))) {
    inds <- which(module2_vec == unique(module2_vec)[i])
    start <- min(inds)
    end <- max(inds)
    decorate_heatmap_body("Similarity", {
        grid.rect(
          x = unit((start - 1) / ncol(plot_mat), "npc"),
          y = unit(1 - (start - 1) / nrow(plot_mat), "npc"),
          width = unit((end - start + 1) / ncol(plot_mat), "npc"),
          height = unit((end - start + 1) / nrow(plot_mat), "npc"),
          just = c("left", "top"), gp = gpar(col = "black", fill = NA, lwd = 1.5)
        )
    })
  }
  dev.off()
}

# F1C
plot_F1C <- function(
  dotplot_data, group_bounds, colorbar_df,
  file = "f1_NMF_module_GO_dotplot_slct.pdf"
) {
  p <- ggplot(dotplot_data, aes(y = Module, x = Description)) +
    geom_rect(
      data = group_bounds, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = NA, color = "black", linewidth = 0.5
    ) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradientn(
      name = "Adjusted\np-value", colors = c("#B2182B", "#d8431a", "#ea9036",
        "#ecc377", "#D1E5F0", "#67A9CF"),
      values = scales::rescale(c(0, 0.01, 0.03, 0.04, 0.05, 0.1)), limits = c(0, 0.1),
      oob = scales::squish
    ) +
    scale_size_continuous(
      name = "Gene Count", range = c(0.5, 3.8), limits = c(2, max(dotplot_data$Count))
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 13.5),
      plot.title = element_text(hjust = 0.5), plot.margin = margin(1, 110, 10, 30),
      legend.position = c(1.09, 0.1)
    ) +
    labs(x = "", y = "")
  group_colors = c(
    "#889BCB", "#7FB0A3", "#b992a4", "#97BC68", "#E6B85C", "#E89A90", "#C49A6C",
    "#7C7C7C", "#A6D1C5"
  )
  p_bar <- ggplot(colorbar_df, aes(x = x, y = 0.5, fill = Group)) +
    geom_tile(width = 1, height = 0.3) +
    scale_fill_manual(values = group_colors) +
    theme_void() + theme(legend.position = "none", plot.margin = margin(2, 0, 0, 0))
  final_plot <- p_bar / p + plot_layout(heights = c(0.2, 6))
  ggsave(file, final_plot, width = 10, height = 5.6)
}

# F1D
plot_F1D <- function(df_2d, file = "f1_compare2d_NMF_SPP1_TN_sc_meta.pdf") {
  my_colors <- c(
    "#5A739E", "#5A937B", "#866A7A", "#6A8E4C", "#d6ac56", "#B76F5C", "#b68959",
    "#707070", "#7DA69A"
  )
  module_names <- unique(df_2d$module)
  color_map <- setNames(my_colors[1:length(module_names)], module_names)
  p <- ggplot(df_2d[1:8, ], aes(y = mean_logFC, x = SMD)) +
    geom_errorbar(
      aes(xmin = ci.lb, xmax = ci.ub), color = "#8b8b8b", width = 0.03
    ) +
    geom_point(size = 2, aes(color = module)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ggrepel::geom_text_repel(
      aes(label = module, color = module),
      size = 4.2, max.overlaps = 10, box.padding = 0.5, segment.color = NA,
      fontface = "bold"
    ) +
    scale_color_manual(values = color_map) +
    coord_cartesian(
      ylim = c(min(df_2d$mean_logFC[1:8] - 0.01), max(df_2d$mean_logFC[1:8]) + 0.06),
      xlim = c(min(df_2d$SMD[1:8]) - 0.02, max(df_2d$SMD[1:8]) + 0.05)
    ) +
    ylab("Mean log2(T/N) module score (sc)") +
    xlab("Standardized mean difference (meta: TCGA+GTEx)") +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none", axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      axis.line = element_blank()
    )
  ggsave(file, p, width = 4.5, height = 3.8)
}

# F1E
plot_F1E <- function(crc, file = "f1_MP6SPP1_CRC_featureplot.pdf") {
  p1 = DimPlot(
    crc, group.by = "ClusterName", cols = paletteer::paletteer_d(
      "miscpalettes::pastel")[-c(8, 1)]
  ) +
    NoAxes() + labs(title = "CRC", X = NULL, Y = NULL) +
    theme(plot.title = element_text(hjust = 0.65))
  p2 = FeaturePlot(crc, "M6_SPP1") +
    viridis::scale_color_viridis() + NoAxes() + labs(title = "M6_SPP1 score",
      X = NULL, Y = NULL) +
    theme(plot.title = element_text(hjust = 0.75))
  p = p1 / p2
  ggsave(file, p, height = 6.5, width = 5)
}


# F1G
plot_F1G <- function(patient_scores, file = "f1_M6SPP1_ICB_boxplot.pdf") {
  resp_lv <- c("NR_Pre", "NR_Post", "R_Pre", "R_Post")
  colors <- c(
    NR_Pre = "#5395a8", NR_Post = "#2A5567", R_Pre = "#bd7f66", R_Post = "#8C3D2E"
  )
  plist = list()
  for (n in names(patient_scores)) {
    pseudo_bulk_data <- patient_scores[[n]]
    pseudo_bulk_data$resp <- factor(pseudo_bulk_data$resp, levels = resp_lv)
    y_max <- max(pseudo_bulk_data$M6_SPP1, na.rm = TRUE) *
      1.065
    p <- ggplot(pseudo_bulk_data, aes(x = resp, y = M6_SPP1, color = resp)) +
      scale_colour_manual(values = colors) +
      geom_boxplot(outlier.size = -1) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.8) +
      ggpubr::stat_compare_means(
        comparisons = list(c("NR_Post", "R_Post")), method = "t.test", size = 3
      ) +
      theme_classic() + scale_y_continuous(limits = c(NA, y_max)) +
      ggtitle(n) +
      labs(y = "AUCell score of M6_SPP1") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 11),
        axis.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        plot.margin = margin(r = 0, l = 20, t = 3)
      )
    if (n != "TNBC")
      p = p + theme(legend.position = "none")
    plist[[n]] = p
  }
  p_final <- wrap_plots(plist, guides = "collect", ncol = 2)
  ggsave(file, plot = p_final, width = 4.6, height = 4.4)
}

# S1A
plot_S1A <- function(panMye, file = "s1_NMF_panMarco_major_dotplot.pdf") {
  p = DotPlot(
    panMye, names(panMye@meta.data)[17:25],
    group.by = "MajorCluster", col.max = 2, col.min = -0.75
  ) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7),
      axis.line = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    guides(
      color = guide_colorbar(title = "Zscore of mean\nModule score"),
      size = guide_legend(title = "Fraction of cells\nin group (%)")
    ) +
    scale_y_discrete(
      limits = rev(c("Mast", "Mono", "Macro", "cDC1", "cDC2", "cDC3", "pDC", "Prolif."))
    ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1),
      legend.key.size = unit(12, "pt"), legend.title = element_text(size = 10),
      legend.text = element_text(size = 10)
    ) +
    scale_color_gradientn(colours = c("#fbfdfb", "#dcece0", "#bdd1c3",
      "#435c4c"))
  ggsave(file, p, height = 5.2, width = 4.5)
}

# S1B
plot_S1B <- function(GO_result_list, file = "s1_NMF_GO_barplot.pdf") {
  p_GO_list = list()
  for (i in names(GO_result_list)) {
    go_enrichment_pathway <- GO_result_list[[i]] %>%
      arrange(-p.adjust)
    go_enrichment_pathway$Description = paste0(
      stringr::str_wrap(go_enrichment_pathway$Description, width = 33), "\n"
    )
    go_enrichment_pathway$Description <- factor(
      go_enrichment_pathway$Description, levels = go_enrichment_pathway$Description
    )
    p_GO_list[[i]] = ggplot(go_enrichment_pathway, aes(x = Description, y =
      Count, fill = p.adjust)) +
      geom_bar(stat = "identity", width = 0.7) +
      coord_flip() + scale_fill_gradientn(
        colours = c("#a499aa", "#cbbecb", "#cbbecb", "#e3d9df"),
        name = "Adjusted\np-value"
      ) +
      labs(x = "", y = "Count", title = i, ) +
      ggthemes::theme_few() + theme(
        axis.text.y = element_text(size = 13.5, colour = "black", vjust = 0.7),
        axis.text.x = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 7, r = -200, b = 1, l = -200),
        legend.margin = margin(r = 0, l = 0), legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.4, "cm")
      )
  }
  ps = wrap_plots(p_GO_list, ncol = 3)
  ggsave(file, ps, height = 14.5, width = 17.5)
}

# S1C
plot_S1C <- function(seu_module, file = "s1_NMF_module_score_cyq_umap.pdf") {
  ps = FeaturePlot(
    seu_module, features = rownames(seu_module),
    reduction = "umap", cols = c("#f8eaea", "#ebb9b7", "#b1211c"),
    pt.size = 0.3, ncol = 3
  )
  for (i in 1:length(ps)) {
    ps[[i]] <- ps[[i]] + theme_void() + labs(color = "Module\nscore") +
      theme(
        legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"),
        legend.title = element_text(size = 11), legend.text = element_text(size = 9.5),
        plot.title = element_text(hjust = 0.5)
      )
  }
  ps = ps + patchwork::plot_layout(guides = "collect")
  ggsave(file, ps, height = 5.3, width = 6)
}

# S1D
plot_S1D <- function(panMacro, file = "s1_NMF_panMacro_sub_dotplot.pdf") {
  p = DotPlot(
    panMacro, names(panMacro@meta.data)[17:25],
    group.by = "ClusterName", col.max = 2, col.min = -0.75
  ) +
    theme(
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.7),
      axis.line = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    guides(
      color = guide_colorbar(title = "Zscore of mean\nModule score"),
      size = guide_legend(title = "Fraction of cells\nin group (%)")
    ) +
    scale_y_discrete(
      limits = rev(
        c(
          "Macro_ISG15", "Macro_LYVE1", "Macro_IL1B", "Macro_INHBA",
          "Macro_FN1", "Macro_GPNMB", "Macro_PPARG", "Macro_C1QC", "Macro_SPP1",
          "Macro_NLRP3", "Myeloid_MKI67", "Macro_CX3CR1", "Macro_VCAN"
        )
      )
    ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1),
      legend.key.size = unit(12, "pt"), legend.title = element_text(size = 10),
      legend.text = element_text(size = 10)
    ) +
    scale_color_gradientn(colours = c("#fffefe", "#f5e9ec", "#fdeaef",
      "#e7a8d1", "#6E016B"))
  ggsave(file, p, height = 5.1, width = 5.2)
}

# S1E
plot_S1E <- function(
  exp, comp, ylabname = "MP6_SPP1 module score",
  file = "s1_TcgaGtex_MP6_SPP1.pdf"
) {
  p1 <- ggboxplot(
    exp, x = "Tissue", y = "Gene", fill = NULL, outlier.size = 0.5, xlab = "",
    ylab = ylabname, color = "Group", palette = c("#134a77", "#85391b"),
    ggtheme = theme_classic()
  ) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(), legend.title = element_blank(),
      axis.line = element_blank(), panel.border = element_rect(color =
        "black", fill = NA, size = 0.5)
    )
  p2 <- p1 + stat_pvalue_manual(
    comp, x = "Tissue", y.position = 2.2, label = "p.signif", position =
      position_dodge(0.8)
  )
  ggsave(file, p2, width = 9.5, height = 3.2)
}

# S1F
plot_S1F <- function(cell_scores, file = "s1_M6SPP1_ICB_violinplot.pdf") {
  plist = list()
  for (n in names(cell_scores)) {
    df <- na.omit(cell_scores[[n]])
    df$resp <- factor(df$resp, levels = c("NR_Pre", "NR_Post", "R_Pre", "R_Post"))
    existing_groups <- unique(df$resp)
    all_comparisons <- list(c("NR_Post", "R_Post"), c("NR_Pre", "R_Pre"))
    valid_comparisons <- Filter(
      function(comp) all(comp %in% existing_groups), all_comparisons
    )
    if (length(valid_comparisons) == 0) {
      message(paste("Skipping", n, ": not enough groups to compare"))
      next
    }
    colors <- c(
      NR_Pre = "#7cafbd", NR_Post = "#3e6778", R_Pre = "#cda595", R_Post = "#b16b5d"
    )
    y_max <- max(df$M6_SPP1, na.rm = TRUE) *
      1.25
    if (length(valid_comparisons) == 1)
      y_max <- max(df$M6_SPP1, na.rm = TRUE) *
        1.12
    p <- ggplot(df, aes(x = resp, y = M6_SPP1, fill = resp)) +
      geom_violin(trim = FALSE, scale = "width", color = NA, alpha = 0.7) +
      scale_fill_manual(values = colors) +
      theme_classic() + labs(title = n, x = NULL, y = "AUCell score of M6_SPP1") +
      scale_y_continuous(limits = c(NA, y_max)) +
      stat_compare_means(comparisons = valid_comparisons, method = "t.test",
        size = 2.8) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 11),
        axis.text = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
        axis.line = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        plot.margin = margin(r = 0, l = 20, t = 3)
      )
    if (n != "TNBC")
      p = p + theme(legend.position = "none")
    plist[[n]] = p
  }
  p_final <- wrap_plots(plist, ncol = 2) +
    plot_layout(guides = "collect")
  ggsave(file, plot = p_final, width = 5, height = 4.3)
}

# S1G
plot_S1G <- function(
  count_1374, count_80, morandi_colors,
  file = "s1_lib_category_1374_80.pdf"
) {
  p1 = ggplot(count_1374, aes(x = 2, y = n, fill = Category)) +
    geom_bar(stat = "identity", width = 0.5, color = "white", size = 0.5) +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(values = morandi_colors) +
    theme_void() + theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right", plot.margin = margin(r = 0)
    ) +
    xlim(1, 2.5)
  p2 = ggplot(count_80, aes(x = 2, y = n, fill = Category)) +
    geom_bar(stat = "identity", width = 0.5, color = "white", size = 0.5) +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(values = morandi_colors) +
    theme_void() + theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right", plot.margin = margin(l = 0)
    ) +
    xlim(1, 2.5)
  p = p1 + p2 + plot_layout(ncol = 2, guides = "collect")
  ggsave(file, plot = p, width = 7.5, height = 6, dpi = 300)
}
