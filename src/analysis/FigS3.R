#=====================================================================#
#--- Libraries ------
#=====================================================================#
library(ComplexHeatmap)
devtools::load_all("/groups/stark/vloubiere/vlite/")
source("src/functions/plotting_fun.R")
source("src/functions/ggplot2_theme.R")
source("src/functions/find_elbow.R")
source("src/functions/Genomic_intervals.R")

#=====================================================================#
#--- Parse input ------
#=====================================================================#
##-- EPD ----
EPD <- fread("db/annotations/GRCm38.v102.EPDnew_promoters.txt.gz")

##-- ChIP ----
Pprc1_peaks <- fread("db/bed/Pprc1_all_peaks_1kb_centered.txt.gz")
Pprc1_peaks[, rep1_CPM := vlite::bwCoverage(.SD, bw = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_noTreatment_exp5_rep1.bw")]
Pprc1_peaks[, rep2_CPM := vlite::bwCoverage(.SD, bw = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_noTreatment_exp5_rep2.bw")]

##-- PROseq ----
PROseq_down_geneIDs <- fread("db/FC_tables/PROseq/Pprc1.depl_transcript_spikeIn_norm_DESeq2_FC.txt")[diff == "Down-regulated",tstrsplit(gene_id, "__",keep = 1)[[1]] ]
PROseq_dds <- readRDS("db/dds/PROseq/Pprc1.depl_transcript_spikeIn_norm_DESeq2.dds")

##-- MS data ----
MS.dat <-  readxl::read_excel("/groups/stark/sabath/MS/241202_PRM_PPRC1/data1.xlsx")
MS.dat <- as.data.table(MS.dat)
# Keep only needed columns
MS.dat <- MS.dat[, c("Genes", "WT")]

# Convert types
MS.dat[, Genes := as.character(Genes)]
MS.dat[, WT := as.numeric(WT)]
MS.dat <- MS.dat[!is.na(WT) & !is.na(Genes), ]
MS.dat <- MS.dat[WT >0] # Keep only positive values for log10 transformation
MS.dat[, WT_log10 := log10(WT)] # Log10 transform

# Define highlighting groups
MS.dat[, group := NA]
exclude_ribosomal <- c("Rps6kc1","Rps6ka3","Rps6kb2","Rps6kb1","Rps6ka2","Rps6ka5","Rps6ka6","Rps6ka4","Rps6kc1", "Rps6ka1", "Rps19bp1")

MS.dat[, group := fcase(
  grepl("^(Rpl|Rps)", Genes) & !(Genes %in% exclude_ribosomal), "Ribosomal",
  Genes == "Pprc1", "PPRC1"
)]

#=====================================================================#
#--- Figure S3A: Pprc1 depletion MS ------
#=====================================================================#
ggplot(MS.dat,
       aes(x = "", y = WT_log10)) +
  geom_boxplot(width = 0.3) +
  geom_point(data = . %>% dplyr::filter(Genes == "Pprc1"), color = "#144e25", size = 2) +
  geom_text(data = . %>% dplyr::filter(Genes == "Pprc1"), aes(label = Genes),  color = "#144e25", hjust = -0.5, size = 7/.pt) +
  labs(x= NULL, y = "Relative protein abundance (log10)") +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_discrete(expand = expansion(add = 0.3), continuous.limits = c(1,1.3)) +
  scale_y_continuous(expand = c(0,0)) +
  annotate(geom = "text", x = 1.3, y = 0, size = 7/.pt,
           label = paste0("n=", format(nrow(MS.dat), big.mark = ",")))
ggsave("results/FigS3A.pdf", width = 36, height = 45, units = "mm")


#=====================================================================#
#--- Figure S3C: PROseq PCC dendo ------
#=====================================================================#
# Change names
PROseq_dds@colData@rownames <- gsub("Pprc1\\.AID_IAA_", "", PROseq_dds@colData@rownames)
pdf(file = "results/FigS3C.pdf", width = 3, height = 4)
Dendroplot_vsaNorm(PROseq_dds)
dev.off()

#=====================================================================#
#--- Figure S3D: Pprc1 ChIP-seq PCC ------
#=====================================================================#
pcc <- cor(log10(Pprc1_peaks$rep1_CPM), log10(Pprc1_peaks$rep2_CPM), method = "pearson")

FigS3D <- Pprc1_peaks %>%
  ggplot(aes(x = log10(rep1_CPM), y = log10(rep2_CPM))) +
  ggpointdensity::geom_pointdensity(size = 0.2) +
  viridis::scale_color_viridis(name = "Density") +
  geom_smooth(method = "lm", se = FALSE, color = "tomato", alpha = 0.8, linewidth = 0.25) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5,
           label = paste0("Pearson's r = ", round(pcc, 3)), size = 6 / .pt ) +
  labs(x = expression(log[10]("CPM rep1")),
       y = expression(log[10]("CPM rep2"))) +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0))

#=====================================================================#
#--- Figure S3E: Confident peaks cutoff ------
#=====================================================================#
ChIP_cutoff <- find_elbow(Pprc1_peaks$Pprc1_merged_CPM)
setorder(Pprc1_peaks, -Pprc1_merged_CPM )
ChIP_cutoff.idx <- which(Pprc1_peaks$Pprc1_merged_CPM <ChIP_cutoff)[1]

FigS3E <- Pprc1_peaks[, .(idx = .I, Pprc1_merged_CPM)] %>% 
  ggplot(aes(x = idx, y = Pprc1_merged_CPM)) + 
  geom_line() +
  geom_vline(xintercept = ChIP_cutoff.idx, linetype = "dashed", linewidth = 0.25) +
  labs(y = "PPRC1 binding (CPM)",
       x = "rank")

cowplot::plot_grid(FigS3D,FigS3E, nrow = 1)
ggsave("results/FigS3D,E.pdf", units = "mm", width = 90, height = 45)
