#=================================================#
# Libraries
#=================================================#
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(circlize)
})
devtools::load_all("/groups/stark/vloubiere/vlite/")
source("src/functions/Genomic_intervals.R")
source("src/functions/plotting_fun.R")

#=================================================#
# Parse input
#=================================================#
meta_ChIP <- readRDS("Rdata/ChIP_metadata_peakCalling.rds")
meta_ChIP[, c("assay", "experiment", "cell_line", "antibody", "treatment") := tstrsplit(condition, "_")]

# Get all Kat2a peaks
bed_paths <- unique(meta_ChIP[grepl("Kat2a", condition), peaks])
Kat2a.beds <- lapply(bed_paths, vlite::importBed)
Kat2a.peaks <- rbindlist(Kat2a.beds, idcol = "condition")

# Center the peaks on the summit
extend_summit = 1000
Kat2a.peaks.centered <- copy(Kat2a.peaks)
Kat2a.peaks.centered[, summit := start + peak]
Kat2a.peaks.centered[, start := summit - extend_summit+1]
Kat2a.peaks.centered[, end := summit +extend_summit]

# Collapse peaks and take the coordinates for the one with the highest score
Kat2a.peaks.centered[, contig := vlite::collapseBed(Kat2a.peaks.centered, return.idx.only = T, ignore.strand = T)]
setorder(Kat2a.peaks.centered, -score)
Kat2a.peaks.collapse <- Kat2a.peaks.centered[, .SD[1], contig]

# Remove peaks overlapping blacklisted regions
mm10.blacklist <- fread("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz")
setnames(mm10.blacklist, c("seqnames", "start", "end", "category"))
mm10.blacklist[, start := start+1] # change to 1-based
Kat2a.peaks.collapse <- Kat2a.peaks.collapse[intersectRanges(Kat2a.peaks.collapse,mm10.blacklist, mode = "nonoverlap",return_index = T, ignore.strand = T )]


##--- Flow data ----
flow_data <- readRDS("/groups/stark/sabath/Bioinfo/ReporterInduction/gfp_filtered.rds")
flow_data_filtered <- flow_data[cellline %in% c("FNR5", "FNR6") & BL1.A > quantile(flow_data$BL1.A, .001),]
flow_data_filtered[, cellline := fcase(
  cellline == "FNR5", "Rpl5",
  cellline == "FNR6", "Rps24"
)]
setnames(flow_data_filtered, c("name.y", "BL1.A"), c("fusion", "GFP"))

 flow_stats <- flow_data_filtered %>%
  dplyr::as_tibble() %>%
  dplyr::filter(fusion %in% c("VPR", "DsRed")) %>%
  dplyr::mutate(fusion = factor(fusion, levels = c("VPR", "DsRed"))) %>%  # estimate = VPR - DsRed
  dplyr::group_by(cellline) %>%
  dplyr::group_modify(~ broom::tidy(wilcox.test(GFP ~ fusion,
                                                data = .x
                                                alternative = "greater"))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p.adj = p.adjust(p.value, method = "BH"))

flow_data_filtered %>%
  ggplot(aes(cellline, GFP)) +
  geom_violin(aes(fill = fusion),colour = NA,
              alpha = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(aes(fill = fusion), outlier.shape = NA,
               linewidth = 0.235,
               median.linewidth = 0.235,
               width = 0.2,
               position = position_dodge(0.9)) +
  scale_fill_discrete(palette = c("brown", "dodgerblue")) +
  labs(x = NULL, y = "GFP intensity (log2)", fill = "TetR fusion") +
  theme(legend.position = "bottom")

##--- Pprc1 targets ----
Pprc1_PROseq <- fread("db/FC_tables/Pprc1_PROseq_DESeq2_FC_Pprc1Bound_TCT.csv")
Pprc1.filt <- fread("db/bed/Pprc1_confident_peaks_1kb_centered.txt")

Pprc1_DirTargetGenes <- Pprc1_PROseq[(DirectTarget), unique(gene_id)]

##--- Hcfc1 PROseq ----
Hcfc1_PROseq <- fread("/groups/stark/nemcko/Hcfc1_paper/Hcfc1/db/FC_tables/PROseq/transcript/Hcfc1.depl.N_transcript_spikeIn_norm_DESeq2_FC.txt")

#=================================================#
# Fig. S8D: Effect of Hcfc1 depletion on Pprc1 targets ----
#=================================================#
FigS8D <- Hcfc1_PROseq[gene_id %in% Pprc1_DirTargetGenes & condition == "IAA_3h"] %>%
  ggplot(aes("", y = log2FoldChange)) +
  geom_boxplot(fill = "#c2ae8d",
               linewidth = 0.235,
               median.linewidth = 0.235,
               outlier.size = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.235,
             linetype = "dashed") +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = NULL, y = "log2 Fold change\n(-/+ Hcfc1 PRO-seq)")
ggsave(filename = "results/FigS8D.pdf",
       plot = FigS8D,
       width = 25, height = 30, units = "mm" )

#=================================================#
# Fig. S8H: Kat2a PCC ----
#=================================================#
bws <- meta_ChIP[grepl("Kat2a|Parental", condition) & !grepl("merged", bw),]

bws[, cell_line := fcase(
  cell_line == "Parental", "Control",
   cell_line == "Kat2bKO-Kat2a.AID", "KAT2A",
   cell_line == "Kat2bKO-Kat2a.AID-rescue.GFP", "GFP",
   cell_line == "Kat2bKO-Kat2a.AID-rescue.GFP.Kat2a.dN", "KAT2A-dN",
  cell_line == "Kat2bKO-Kat2a.AID-rescue.GFP.Kat2a.wt", "KAT2A"
    
  )]
bws[, name := paste(cell_line, antibody, treatment, replicate, sep = "_")]

Kat2a.peaks.collapse[, (bws$name) := lapply(bws$bw, function(x) log10(vlite::bwCoverage(Kat2a.peaks.collapse, x) + 0.1))]

clust_pcc <- as.matrix(Kat2a.peaks.collapse[, .SD,
                                              .SDcols = bws$name]) |>
  Hmisc::varclus(similarity = "pearson", trans = "none")

pdf("results/FigS8H.pdf", width = 47/25.4, height = 60/25.4, pointsize = 6)

par(
  mar = c(3, 5, 2, 2),
  tcl = -0.2,
  mgp = c(3, 0.5, 0),
  las = 1,
  cex = 0.8, lwd = 0.5
)

plot(
  clust_pcc,
  lwd = 0.5,
  hang = -1,
  xlab = NA,
  ylab = "PCC",
  las = 1,        # rotate leaf labels vertically so they fit narrow spacing
  cex = 1       # shrink leaf label font relative to pointsize = 6
)

labs <- factor(gsub("_rep[12]", "", colnames(clust_pcc$sim)))
Cc <- rainbow(nlevels(labs))[labs]

points(
  seq_along(clust_pcc$hclust$order),
  rep(0, length(clust_pcc$hclust$order)),
  col = Cc[clust_pcc$hclust$order],
  pch = 19,
  cex = 0.8
)

dev.off()

#=================================================#
# Fig. S8I: Kat2a genome browser ----
#=================================================#
# Rps1 locus
FigS8I_region <- "chr18:60,703,180-60,866,287"
tracks <- c(
  "Kat2a"      = "/groups/stark/sabath/CutandRun/260618_Kat2a_Rescue/db/bw/CUTNRUN/V5_cGW5cl1_exp2_-_merged.bw",
  "Kat2a +IAA" = "/groups/stark/sabath/CutandRun/260618_Kat2a_Rescue/db/bw/CUTNRUN/V5_cGW5cl1_exp2_IAA_merged.bw",
  "Input"   = "/groups/stark/sabath/CutandRun/260618_Kat2a_Rescue/db/bw/CUTNRUN/V5_FNR26_exp2_-_merged.bw",
  "Kat2a wt"   = "/groups/stark/sabath/CutandRun/260618_Kat2a_Rescue/db/bw/CUTNRUN/GFP_cKeS241_exp2_IAADox_merged.bw",
  "Kat2a dN"   = "/groups/stark/sabath/CutandRun/260618_Kat2a_Rescue/db/bw/CUTNRUN/GFP_cKeS240_exp2_IAADox_merged.bw",
  "GFP"        = "/groups/stark/sabath/CutandRun/260618_Kat2a_Rescue/db/bw/CUTNRUN/GFP_cKeS239_exp2_IAADox_merged.bw"
)
track.cols <- c("#d85f43","#d85f43","#babdbc","#4f9b97","#4f9b97","#babdbc")

# Annotations
genes <- rtracklayer::import('/groups/stark/vojtek/Resources/01_annotations/mm10_GRCm38_v102/Mus_musculus.GRCm38.102.gtf.gz', format = "gtf") %>%
  as.data.table()
genes[, transcript_version := as.numeric(transcript_version)]
genes[, transcript_support_level := as.numeric(transcript_support_level)]
setorder(genes, -transcript_version, na.last = T)

# Get one transcript per gene
Select_tx <- genes[ tag == "basic" & type == "transcript" & !is.na(transcript_support_level), transcript_id[1], gene_id][,V1]
genes <- genes[transcript_id %in% Select_tx] 
genes[, seqnames := paste0("chr", seqnames)]
genes[, tss := fifelse(strand == "+", start, end)]

# Get Tx id for RPS14
RPS14_tx <- genes[transcript_name == "Rps14-201" & type == "transcript", transcript_id]

# BW tracks
FigS8I_signalTrack <- plot_signal_tracks(bw = tracks, region = FigS8I_region, free_y = F, fill_color = track.cols, fill_alpha = 0.8) + 
  theme_publish() +
  annotate("rect",
           xmin = genes[gene_name == "Rps14" &  type == "transcript", tss +-1000],
           xmax = genes[gene_name == "Rps14" &  type == "transcript", tss + 1000], ymin = -Inf, ymax = Inf, fill = "cornflowerblue", alpha = 0.3) +
  annotate("rect",
           xmin = genes[gene_name == "Ndst1" &  type == "transcript", tss - 1000],
           xmax = genes[gene_name == "Ndst1" &  type == "transcript", tss + 1000], ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.3) +
annotate("rect",
         xmin = genes[gene_name == "Cd74" &  type == "transcript", tss - 1000],
         xmax = genes[gene_name == "Cd74" &  type == "transcript", tss + 1000], ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.3) +
  annotate("rect",
           xmin =  genes[gene_name == "Tcof1" &  type == "transcript", tss - 1000],
           xmax = genes[gene_name == "Tcof1" &  type == "transcript", tss + 1000], ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.3)

FigS8I_annoTrack <- plot_annotation_track(genes, FigS8I_region, arrow_bin = 20000, min_intron_length = 1000, arrow_width = 500, arrow_height = 0.1 ) +
  theme(plot.margin = margin(t = -5, r = 0, b = 0, l = 0))
FigS8I <- cowplot::plot_grid(FigS8I_signalTrack,
                             FigS8I_annoTrack, ncol = 1, align = "v",axis = "lr", rel_heights = c(1, 0.2))
ggsave(filename = "results/FigS8I.pdf", plot = FigS8I, width = 43, height = 50, units = "mm")
