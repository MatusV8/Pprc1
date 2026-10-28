#==============================================================#
# Libraries
#==============================================================#
source("src/functions/plotting_fun.R")
source("src/functions/ggplot2_theme.R")
devtools::load_all("/groups/stark/vloubiere/vlite/")
library(data.table)
library(tidyverse)

#==============================================================#
# Input ----
#==============================================================#
## Bigwigs -----
bw <- list(
  "Pprc1" = "db/bw/ChIP/ChIPseq_exp6_Pprc1.AID_V5_noTreatment_merged.bw",
  "Taf1" = "/groups/stark/sabath/ChIP/251028_ChIP_Hcfc1Depletion2/db/bw/ChIP/Taf1_FNR94_noTreatment_exp2_merged.bw"
)

bw_treatment <- list(
  "noTreatment" = "db/bw/ChIP/ChIPseq_exp6_Pprc1.AID_V5_noTreatment_merged.bw",
  "Flavopiridol" = "db/bw/ChIP/ChIPseq_exp6_Pprc1.AID_V5_Flavopiridol_merged.bw",
  "RNAse" = "db/bw/ChIP/ChIPseq_exp6_Pprc1.AID_V5_RNase_merged.bw",
  "Triptolide" = "db/bw/ChIP/ChIPseq_exp6_Pprc1.AID_V5_Triptolide_rep1.bw"
)

## annotations ----
mm10_gtf <- rtracklayer::import('/groups/stark/vojtek/Resources/01_annotations/mm10_GRCm38_v102/Mus_musculus.GRCm38.102.gtf.gz', format = "gtf") %>% as.data.table()
mm10_gtf[, seqnames := paste0("chr", seqnames)]

Pprc1_PROseq <- fread("db/FC_tables/Pprc1_PROseq_DESeq2_FC_Pprc1Bound_TCT.csv")

#==============================================================#
# Fig 3A: Genome browser track of Rpl3 (Pprc1, Rpb2, Taf1) ----
#==============================================================#
Fig3A_region <-  "chr15:80077452-80084029"

# BW tracks
Fig3A_signalTrack <- plot_signal_tracks(bw = bw, region = Fig3A_region, free_y = F, fill_color = c("#2bbcb6", "gray30"), fill_alpha = 0.8) + 
  theme_publish() +
  scale_x_reverse(expand = c(0,0), position = "top") 

# Annotations
Fig3_annoTrack <- plot_annotation_track(mm10_gtf[transcript_id == "ENSMUST00000081650"], Fig3A_region) +
  scale_x_reverse(expand = c(0,0), limits = c(80084029, 80077452)) +
  theme(plot.margin = margin(t = -5, r = 0, b = 0, l = 0))
cowplot::plot_grid(Fig3A_signalTrack,
                   Fig3_annoTrack, ncol = 1, align = "v",axis = "lr", rel_heights = c(1, 0.2))

#==============================================================#
# Fig 3B: Metaplot around TSS of Pprc1 target genes for (Pprc1, Rpb2, Taf1) ----
#==============================================================#
TSS_up = 500
TSS_down = 2500

Pprc1_targets_TSS <- Pprc1_PROseq[(DirectTarget), .(seqnames, strand, start = tss)]

# Scaled max
bwPlotAverage(Pprc1_targets_TSS, tracks = bw, center = "center", upstream = TSS_up, downstream = TSS_down, ignore.strand = F,
              colors = c(c("#2bbcb6", "gray30")), normalisation  = "max") +
  theme_publish() +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25) +
  labs(x = "Distance from TSS (bp)",
       y = "Normalised binding",
       color = NULL) +
  guides(fill = "none") +
  theme(legend.position = "inside", legend.position.inside = c(0.8, 0.8),
        legend.justification = c(0.5, 0.5),
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0, "pt"))
ggsave("results/Fig3B_maxsNorm.pdf", units = "mm", width = 40, height = 35)

#=====================================================================#
#----- Fig 3C: RNAse treatment ----
#=====================================================================#
Pprc1_PROseq[, start := tss]
Fig3C <- bwPlotAverage(Pprc1_PROseq[(DirectTarget)],
                       tracks = bw_treatment,
                       colors = c("#2bbcb6", "cornflowerblue", "brown", "brown2"),
                       center = "center", region.start = "TSS", region.end = "TES", upstream = 750, downstream = 750) +
  theme_publish() +
  scale_y_continuous() +
  scale_x_continuous(breaks = c(-750, 0, 750)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25) +
  labs(x = "Distance from TSS (bp)",
       y = "Pprc1 binding (CPM)",
       color = "treatment") +
  guides(fill = "none") +
  theme(legend.position = "inside", legend.position.inside = c(0.3, 0.7),
        legend.justification = c(0.5, 0.5),
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0, "pt"))
ggsave(filename = "results/Fig3C.pdf", plot = Fig3C, units = "mm", width = 40, height = 35)
