#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
input_dir <- Sys.getenv("PAPER_INPUT", file.path(script_dir, "../audit/F7_S7/inputs"))
output_dir <- Sys.getenv("PAPER_OUTPUT", file.path(script_dir, "output/F7_S7"))
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(stringr)
  library(patchwork)
  library(metafor)
  library(forestplot)
  library(grid)
})
# F7A
plot_F7A <- function() {
  inputs <- readRDS(file.path(input_dir, "F7A.rds"))
  colors <- c("#498a85", "#cc7069")
  pseudo_bulk_data <- transform(inputs, CTSL = CTSL_scaled)
  ratio_df <- pseudo_bulk_data %>%
    group_by(cancer, label) %>%
    summarise(med = median(CTSL), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = label, values_from = med) %>%
    mutate(ratio = Tumor - Normal) %>%
    arrange(desc(ratio))
  pseudo_bulk_data$cancer <- factor(pseudo_bulk_data$cancer, levels = ratio_df$cancer)
  p <- ggplot(pseudo_bulk_data, aes(x = cancer, y = CTSL, color = label, split.by = label)) +
    scale_colour_manual(values = colors) +
    geom_boxplot(outlier.size = -1, width = 0.5,
      position = position_dodge(width = 0.6), size = 0.5) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.6,
      jitter.width = 0.1), size = 0.6, alpha = 0.9) +
    theme_classic() +
    scale_y_continuous(limits = c(NA,
      2.1)) +
    labs(y = "Expression of CTSL in macrophages", x = "") +
    theme(axis.text = element_text(color = "black"),
      legend.position = "top", legend.title = element_blank(), panel.border = element_rect(color = "black",
        fill = NA, size = 0.5), axis.line = element_blank(), axis.text.x = element_text(angle = 45,
        hjust = 1, size = 9))
  comp <- compare_means(CTSL ~ label, group.by = "cancer", data = pseudo_bulk_data, symnum.args = list(cutpoints = c(0,
    1e-04, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "")), p.adjust.method = "holm")
  p2 <- p + stat_pvalue_manual(comp, x = "cancer", y.position = 1.9, label = "p.signif",
    position = position_dodge(0.6))
  ggsave(file.path(output_dir, "F7A.pdf"), p2, width = 4.8, height = 3.8)
}
# F7D
plot_F7D <- function() {
  inputs <- readRDS(file.path(input_dir, "F7D.rds"))
  cox_df <- inputs
  colnames(cox_df)[1] <- "cancer"
  cox_df1 = cox_df %>%
    filter(!cancer %in% c("DLBC", "UVM", "SARC"))
  cox_df1$logHR <- log(cox_df1$HR)
  cox_df1$SE <- (log(cox_df1$upper95) - log(cox_df1$lower95)) / (2 * 1.96)
  meta <- rma(yi = logHR, sei = SE, data = cox_df1, method = "REML")
  pan <- data.frame(cancer = "Pan-cancer", HR = exp(meta$b), lower95 = exp(meta$ci.lb), upper95 = exp(meta$ci.ub),
    p_value = meta$pval)
  df <- rbind(cox_df1[, c("cancer", "HR", "lower95", "upper95", "p_value")], pan)
  df$star <- cut(df$p_value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***",
    "**", "*", ""))
  df$col <- ifelse(df$HR > 1, "#ca7377", "#a5bbe6")
  tabletext <- cbind(c("Cancer type", df$cancer), c("Sig.", as.character(df$star)))
  p = forestplot(labeltext = tabletext, mean = c(NA, df$HR), lower = c(NA, df$lower95), upper = c(NA,
    df$upper95), is.summary = c(TRUE, rep(FALSE, nrow(df) - 1), TRUE), xlog = TRUE, xlab = "Hazard Ratio (95% CI)",
  txt_gp = fpTxtGp(label = list(gpar(cex = 1), gpar(cex = 1.3)), ticks = gpar(cex = 0.8),
    xlab = gpar(cex = 1, fontface = "bold")), shapes_gp = fpShapesGp(summary = "#bb4f43",
    box = lapply(c("black", df$col), function(x) gpar(fill = x, col = x)), lines = gpar(col = "black"),
    zero = gpar(col = "#7a7a7a", lty = 2)))
  pdf(file.path(output_dir, "F7D.pdf"), width = 6.5, height = 5.5)
  print(p)
  dev.off()
}
# F7E
plot_F7E <- function() {
  inputs <- readRDS(file.path(input_dir, "F7E.rds"))
  resp_lv <- c("NR_Pre", "NR_Post", "R_Pre", "R_Post")
  colors <- c(NR_Pre = "#5395a8", NR_Post = "#2A5567", R_Pre = "#bd7f66", R_Post = "#8C3D2E")
  fs = 6
  plist = list()
  macro_slct = inputs
  names(macro_slct) = c("CRC", "NSCLC")
  for (n in names(macro_slct)) {
    obj <- macro_slct[[n]]
    pseudo_bulk_data <- obj %>%
      mutate(resp = str_replace(resp, "-", "_"), resp = factor(resp, levels = resp_lv),
        CTSL = CTSL_scaled)
    y_max <- max(pseudo_bulk_data$CTSL, na.rm = TRUE) + 0.5
    p <- ggplot(pseudo_bulk_data, aes(x = resp, y = CTSL, color = resp)) +
      scale_colour_manual(values = colors) +
      geom_boxplot(outlier.size = -1, size = 0.4, width = 0.5) +
      geom_jitter(width = 0.1,
        size = 0.5, alpha = 0.8) +
      ggpubr::stat_compare_means(comparisons = list(c("NR_Post",
        "R_Post")), method = "t.test", size = fs / .pt, label.y = y_max - 0.4, bracket.size = 0.2) +
      theme_classic() +
      scale_y_continuous(limits = c(NA, y_max)) + ggtitle(n) +
      labs(y = "Mean expression of CTSL in TAM") +
      theme(plot.title = element_text(hjust = 0.5, size = fs), axis.text = element_text(color = "black",
        size = fs), axis.title = element_text(size = fs), axis.ticks = element_line(linewidth = 0.2),
      panel.border = element_rect(color = "black", fill = NA, size = 0.4), legend.text = element_text(size = fs),
      legend.title = element_blank(), legend.key.size = unit(3, "mm"), legend.box.spacing = unit(1,
        "mm"), axis.line = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
      axis.title.x = element_blank(), plot.margin = margin(r = 20, t = 3))
    if (n != "CRC")
      p = p + theme(legend.position = "none")
    plist[[n]] = p
  }
  p2 = wrap_plots(plist, guides = "collect", ncol = 2)
  ggsave(file.path(output_dir, "F7E.pdf"), p2, width = 3.6, height = 1.5)
}
# S7A
plot_S7A <- function() {
  inputs <- readRDS(file.path(input_dir, "S7A.rds"))
  p = VlnPlot(subset(inputs, cancer %in% c("NPC", "LYM"), invert = T), features = "CTSL",
    split.by = "label", cols = c("#498a85", "#cc7069"), group.by = "cancer", pt.size = 0) +
    xlab("") +
    stat_compare_means(aes(group = split), method = "wilcox.test", label = "p.signif")
  ggsave(file.path(output_dir, "S7A.pdf"), p, width = 10, height = 4)
}
panels <- c("F7A", "F7D", "F7E", "S7A")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) {
  cat("Available panels:", paste(panels, collapse = " "), "\n")
} else {
  if (identical(args, "all"))
    args <- panels
  stopifnot(all(args %in% panels))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  for (panel in args) {
    set.seed(1)
    get(paste0("plot_", panel))()
  }
}
