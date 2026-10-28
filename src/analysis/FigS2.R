#---- Libraries ------------
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(GO.db)
source("src/functions/plotting_fun.R")
source("src/functions/ggplot2_theme.R")

#---- Settings -----
PARAMS <- list(
  dc_cutoff = 10,
  gene_count_cutoff = 3,
  padj_cutoff = 0.01,
  padj_eva = 0.1,
  log2OR_cutoff = 0,
  log2OR_cutoff_eva = 1,
  pseudocount = 1
)

#---- Parse data ------------

# AID screen data
Counts <- fread("db/count_tables/AID_screen_counts.csv")

#---- Supplementary Figure 2D: Correlation between samples ------------
corr_mat <- Counts %>%
  mutate_if(is.numeric, ~log2(.x+1)) %>%
  dplyr::select(-c(gene_id, gene_name)) %>%
  cor(use = "everything") 

GFPlow_corr_mat <- corr_mat[grepl("GFPlow", rownames(corr_mat)), grepl("GFPlow", colnames(corr_mat))]
colnames(GFPlow_corr_mat) <- gsub("_GFPlow", "", colnames(GFPlow_corr_mat))

Input_corr_mat <- corr_mat[grepl("Input", rownames(corr_mat)), grepl("Input", colnames(corr_mat))]
colnames(Input_corr_mat) <- gsub("_Input", "", colnames(Input_corr_mat))


ht_opt(
  heatmap_row_names_gp    = gpar(fontsize = 6, fontfamily = "Helvetica"),
  heatmap_column_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
  heatmap_row_title_gp    = gpar(fontsize = 7, fontfamily = "Helvetica"),
  heatmap_column_title_gp = gpar(fontsize = 7, fontfamily = "Helvetica")
)

GFPlow_cor_hm <- mv_hm(GFPlow_corr_mat, col_labels = c("reporter", "frame", "replicate"),
                       col_title = "GFPlow samples",
                       colors = list(
                         reporter = wong_colors[1:2],
                         frame = wong_colors[3:5],
                         replicate = wong_colors[6:7]
                       ),
                       value_fontsize = 5, names_fontsize = 7)

Input_cor_hm <- mv_hm(Input_corr_mat, col_labels = c("reporter", "frame", "replicate"), col_title = "Input samples",colors = list(
  reporter = wong_colors[1:2],
  frame = wong_colors[3:5],
  replicate = wong_colors[6:7]
),
value_fontsize = 5, names_fontsize = 7)

mm2in <- function(x) x / 25.4

pdf("results/FigS2B.pdf", width = mm2in(120), height = mm2in(60))

draw(GFPlow_cor_hm + Input_cor_hm)

dev.off()