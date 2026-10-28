#=================================================#
# Libraries
#=================================================#
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(circlize)
  library(ComplexHeatmap)
})

devtools::load_all("/groups/stark/vloubiere/vlite/")
source("src/functions/Genomic_intervals.R")
#=================================================#
# Parse input
#=================================================#
Pprc1_PROseq <- fread("db/FC_tables/Pprc1_PROseq_DESeq2_FC_Pprc1Bound_TCT.csv")

# Bin genes by base expression levels
Pprc1_PROseq[, expression_bin := cut(log2(baseMean +1), breaks = quantile(log2(baseMean +1), probs = seq(0, 1, length.out = 21)),
                include.lowest = TRUE, labels = FALSE)]
mm10.blacklist <- fread("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz")
setnames(mm10.blacklist, c("seqnames", "start", "end", "category"))
mm10.blacklist[, start := start+1] # change to 1-based

Pprc1_PROseq[, Blacklisted := .I %in% Pprc1_PROseq[mm10.blacklist, on = .(seqnames, tss >= start, tss <= end), which = TRUE, nomatch = NULL]]

tracks <- list(
  "KAT2A"      = "db/bw/ChIP/CutNRun_exp2_Kat2bKO-Kat2a.AID_V5_noTreatment_merged.bw",
  "GFP"        = "db/bw/ChIP/CutNRun_exp2_Kat2bKO-Kat2a.AID-rescue.GFP_GFP_IAADox_merged.bw",
  "WT"   = "db/bw/ChIP/CutNRun_exp2_Kat2bKO-Kat2a.AID-rescue.GFP.Kat2a.wt_GFP_IAADox_merged.bw",
  "dN"   = "db/bw/ChIP/CutNRun_exp2_Kat2bKO-Kat2a.AID-rescue.GFP.Kat2a.dN_GFP_IAADox_merged.bw"
)


# Subset direct targets and expression matched controls
subset <- Pprc1_PROseq[(DirectTarget)]
controls <- Pprc1_PROseq[diff != "Down-regulated" & !Pprc1_bound & !Blacklisted, {
  
  ctrl_distrib <- table(subset[, expression_bin])
  
  .res <- lapply(names(ctrl_distrib), function(x) {
    pool <- .SD[expression_bin == as.numeric(x)]
    n <- min(ctrl_distrib[x], nrow(pool))
    set.seed(123)
    pool[sample(.N, n)]
  })
  
  rbindlist(.res)
  
}]

# combine
comb_dt <- rbind(subset, controls)

## Read in binding
extend = 5000
comb_dt <- vlite::resizeBed(comb_dt, center = "start", upstream = extend-1, downstream = extend, ignore.strand = F)

comb_dt_binned <- binBed(comb_dt, bins.width = 50)
comb_dt_binned[, names(tracks) := lapply(tracks, vlite::bwCoverage, bed = .SD)]

#--- draw heatmap ---
mat_list <- lapply(names(tracks), function(vc) {
  wide <- dcast(comb_dt_binned, gene_id ~ bin.idx, value.var = vc)
  m <- as.matrix(wide[, -"gene_id", with = FALSE])
  rownames(m) <- wide$gene_id
  m
})
names(mat_list) <- names(tracks)

# transform
mat_list <- lapply(mat_list, function(x){log2(x+1)})

# Sort
sort_ids <- names(sort(rowMeans(mat_list[[1]]), decreasing = T))
mat_list <- lapply(mat_list, function(x) x[sort_ids, , drop = FALSE])
# Get Groups
group_vector <- comb_dt[match(sort_ids, gene_id), fifelse(DirectTarget, "PPRC1 target", "control")]

# Build each heatmap; only show row split/labels on the first one
v5_range <- quantile(mat_list[["KAT2A"]], c(0.01,0.5, 0.99), na.rm = TRUE)
GFP_range <- quantile(mat_list[["WT"]], c(0.01,0.5, 0.99), na.rm = TRUE)
col_list <- list(
  "KAT2A"= colorRamp2(v5_range, c("white", "#fff3ef", "#d85f43")),
  "GFP"   = colorRamp2(GFP_range, c("white", "#f0fafa", "#4f9b97")),
  "WT"    = colorRamp2(GFP_range, c("white","#f0fafa", "#4f9b97")),
  "dN"    = colorRamp2(GFP_range, c("white","#f0fafa", "#4f9b97"))
)

# hm options
ht_opt(
  heatmap_column_title_gp = gpar(fontsize = 7),
  legend_title_gp = gpar(fontsize = 7),
  legend_labels_gp = gpar(fontsize = 5),
  legend_gap         = unit(1, "mm"),   # space between stacked legends
  COLUMN_ANNO_PADDING = unit(0, "mm")
  )
brks = c(0, 0.4, 0.8)
hm_list <- NULL
for(i in names(mat_list)){
  
  hm <- Heatmap(
    mat_list[[i]],
    col = col_list[[i]],
    name = i,
    row_split = group_vector,
    show_heatmap_legend = i %in% c("KAT2A", "GFP"),
    heatmap_legend_param = list(
      at     = brks,
      labels = format(round(brks, 1), nsmall = 1),
      legend_height = unit(3, "mm"),   # vertical legend
      grid_width    = unit(1, "mm")    # width of the color bar
    ),
    use_raster = TRUE,
    bottom_annotation = HeatmapAnnotation(
      tick = anno_mark(
        at = 101,                         # column index
        labels = "",               
        which = "column",
        side = "bottom",
        labels_gp = gpar(fontsize = 6),
        link_height = unit(1.5, "mm")     # tick length
      ),
      show_annotation_name = FALSE
    ),
    raster_quality = 10,
    raster_by_magick = TRUE,
    raster_resize_mat = TRUE,
    row_labels = FALSE,
    row_title_gp = gpar(fontsize = 7),
    row_gap = unit(1, "mm"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    column_title = i,
    border = TRUE
  )
  hm_list <- if (is.null(hm_list)) hm else hm_list + hm
}
pdf("results/Fig5E.pdf", width = 3.1, height = 1.9)
draw(hm_list, ht_gap = unit(1, "mm"), merge_legend = TRUE)
dev.off()
