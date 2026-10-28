#=====================================================================#
#---- Libraries ------------
#=====================================================================#
library(data.table)
library(tidyverse)
library(DESeq2)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggseqlogo)
devtools::load_all("/groups/stark/vloubiere/vlite/")
source("src/functions/ggplot2_theme.R")
source("src/functions/GO_functions.R")
source("src/functions/plotting_fun.R")

#=====================================================================#
#---- Parse data -----
#=====================================================================#
## --- Annotations ----
prom <- fread("db/annotations/GRCm38.v102.EPDnew_promoters.txt.gz")

## --- PROseq ----
Pprc1_PROseq <- fread("db/FC_tables/db/Pprc1_PROseq_DESeq2_FC_Pprc1Bound_TCT.csv")

## -- ChIP-seq ----
Pprc1_peaks <- fread("db/bed/Pprc1_confident_peaks_1kb_centered.txt")

Fig <- list()

#=====================================================================#
#---- Figure 2A: Growth curve  ------------
#=====================================================================#
CompetitonAssay <- fread("../db/Pprc1_competitionAssay.txt")

Fig$Fig2A <- CompetitonAssay %>%
  ggplot(aes(x = day, y = mCherry/100, color = auxin)) +
  geom_point(shape = 16, alpha = 0.5, size = 2, show.legend = F) +
  stat_summary(geom = "line", fun = "mean", show.legend = F) +
  labs(y = "Fraction of living cells",
       x = "Time (days)") +
  scale_color_manual(values = c("gray", "tomato"))

#=====================================================================#
#--- Figure 2B: Pprc1 MA-plot ------
#=====================================================================#
Fig$Fig2B <- Pprc1_PROseq %>%
  mutate(diff = paste0(diff, " (", dplyr::n(), ")"), .by = diff) %>%
  mutate(alpha_val = ifelse(grepl("Unaffected", diff), 0.2, 0.5),
         size_val = ifelse(grepl("Unaffected", diff), 0.5, 0.8)) %>%
  ggplot(aes(x = log10(baseMean), y = log2FoldChange, color = diff)) +
  ggrastr::rasterise(
    geom_point(shape = 16, aes(alpha = alpha_val, size = size_val)),
    dpi = 300
  ) +
  scale_color_manual(values = c("cornflowerblue", "gray", "tomato")) +
  scale_alpha_identity() +
  scale_size_identity() +
  guides(size = "none", alpha = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 1),
        legend.justification = c(1, 1)) +
  labs(x = "log10 Mean expression",
       y = "log2 Fold change",
       color = NULL)

#=====================================================================#
#--- Figure 2C: GO analysis ------
#=====================================================================#
GOs <- getGOs()
GOs <- GOs[, .SD[1], by = .(ensembl_id, go_id, TERM)] # deduplicate
Pprc1_down_GO <- go_enrichment(
  genes = Pprc1_PROseq$gene_id,
  hit = Pprc1_PROseq[, diff == "Down-regulated"],
  GO_table = GOs,
  min_term_size = 10,
  max_term_size = 200)

Fig$Fig2C <- Pprc1_down_GO %>% 
  dplyr::filter(padj < 0.05, term_size < 200, Ontology == "CC") %>%
  arrange(padj) %>%
  dplyr::slice_head(n = 5, by = Ontology) %>%
  mutate(description = paste0(TERM, " (", GO_hit, "/", term_size, ")")) %>%
  ggplot(aes(x = log2OR, y = reorder(description, log2OR))) +
  geom_bar(stat = "identity", fill = "cornflowerblue", alpha = 0.8, width = 0.8) +
  geom_text(
    aes(label = sprintf("%.2e", padj)),
    hjust = 1.2,
    size = 5, size.unit = "pt",
    color = "black"
  ) +
  scale_x_continuous(expand = c(0,0)) +
  theme(legend.position = "none", axis.line.y = element_blank()) +
  labs(y = NULL, x = "log2 Odds ratio")

#=====================================================================#
#--- Figure 1D: Ribosome protein coding genes ------
#=====================================================================#
ribosome_subunits <- c(
  "Fau",
  "Rpl10a",
  "Rpl13",
  "Rpl13a",
  "Rpl14",
  "Rpl18a",
  "Rpl24",
  "Rpl27",
  "Rpl27a",
  "Rpl28",
  "Rpl36",
  "Rps11",
  "Rps23",
  "Rps25",
  "Rpl15",
  "Rpl18",
  "Rpl19",
  "Rpl29",
  "Rpl34",
  "Rpl36a",
  "Rpl4",
  "Rpl6",
  "Rpl7",
  "Rpl7a",
  "Rpl8",
  "Rps8",
  "Rack1",
  "Rpl11",
  "Rpl23",
  "Rpl37a",
  "Rpl38",
  "Rpl39",
  "Rps10",
  "Rps12",
  "Rps14",
  "Rps15",
  "Rps19",
  "Rps24",
  "Rps27l",
  "Rps28",
  "Rpsa",
  "Rpl10",
  "Rpl12",
  "Rpl35a",
  "Rpl9",
  "Rplp0",
  "Rplp2",
  "Rps13",
  "Rps15a",
  "Rps16",
  "Rps18",
  "Rps2",
  "Rps3",
  "Rps4x",
  "Rps5",
  "Rps9",
  "Rpl10l",
  "Rpl22l1",
  "Rpl3",
  "Rpl39l",
  "Rpl3l",
  "Rplp1",
  "Rpl17",
  "Rpl21",
  "Rpl22",
  "Rpl23a",
  "Rpl26",
  "Rpl30",
  "Rpl31",
  "Rpl32",
  "Rpl35",
  "Rpl37",
  "Rpl5",
  "Rps17",
  "Rps20",
  "Rps21",
  "Rps26",
  "Rps27",
  "Rps27a",
  "Rps29",
  "Rps3a1",
  "Rps6",
  "Rps7",
  "Uba52"
)

# Scatter plot
Fig$Fig2D <- Pprc1_PROseq[(gene_name %in% ribosome_subunits)] %>%
  ggplot(aes(x = "Ribosomal\nprotein genes", y = log2FoldChange, color = diff)) +
  geom_jitter(width = 0.4, shape = 16, alpha = 0.7, size = 0.5) +
  theme_publish() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c( "cornflowerblue", "gray")) +
  annotate("text", x = -Inf, y = 0.6,
           label = paste0("n = ", length(ribosome_subunits)),
           hjust = -0.1, vjust = 1,
           size = 6 / .pt, color = "black") +
  labs(x = NULL, y = "log2 Fold change") +
  theme(legend.position = "none", axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_x_discrete(position = "bottom", expand = c(0.5, 0.5))

cowplot::plot_grid(
  Fig$Fig2A,
  Fig$Fig2B,
  Fig$Fig2C,
  Fig$Fig2D,
  nrow = 1,
  rel_widths = c(0.8, 1, 1.6, 0.5)
)
ggsave("results/Fig2_topPanel.pdf", units = "mm", width = 180, height = 35)

#=====================================================================#
#--- Figure 2E: genome-browser ----
#=====================================================================#
# Rps1 locus
Fig2E_region <-  "chr7:45095286-45160468"
Fig2_tracks <- c(
  "Pprc1" = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_noTreatment_exp6_merged.bw",
  "noAB" = "../db/bw/ChIP/ChIPseq_Pprc1.AID_noAb_noTreatment_exp6_merged.bw")

# BW tracks
Fig2E_signalTrack <- plot_signal_tracks(bw = Fig2_tracks, region = Fig2E_region, free_y = F, fill_color = "#2bbcb6", fill_alpha = 0.8) + 
  theme_publish()  +
  annotate("rect", xmin = 45102567, xmax = 45103657, ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.3) +
  annotate("rect", xmin = 45123733, xmax = 45124823, ymin = -Inf, ymax = Inf, fill = "cornflowerblue", alpha = 0.3) +
  annotate("rect", xmin = 45128211, xmax = 45129301, ymin = -Inf, ymax = Inf, fill = "cornflowerblue", alpha = 0.3) +
  annotate("rect", xmin = 45135870, xmax = 45136960, ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.3) +
  annotate("rect", xmin = 45153823, xmax = 45154913, ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.3)

# Annotations
genes <- rtracklayer::import('/groups/stark/vojtek/Resources/01_annotations/mm10_GRCm38_v102/Mus_musculus.GRCm38.102.gtf.gz', format = "gtf") %>%
  as.data.table()
genes[, seqnames := paste0("chr", seqnames)]
select_tx <- genes[transcript_support_level == 1 & transcript_biotype == "protein_coding" & gene_name %in% c("Rps11", "Rpl13a", "Flt3l", "Aldh16a1", "Fcgrt")][order(-gene_version)][!duplicated(gene_name), transcript_id]

Fig2E_annoTrack <- plot_annotation_track(genes[transcript_id %in% select_tx], Fig2E_region, arrow_bin = 5000, min_intron_length = 500) +
  theme(plot.margin = margin(t = -5, r = 0, b = 0, l = 0))
Fig$Fig2E <- cowplot::plot_grid(Fig2E_signalTrack,
                                Fig2E_annoTrack, ncol = 1, align = "v",axis = "lr", rel_heights = c(1, 0.2))

#=====================================================================#
#--- Figure 2F: Pprc1-ChIP promoter frequency ----
#=====================================================================#
# Set promoter length 
prom[, prom_start := tss-99]
prom[, prom_end := tss+100]
Pprc1_peaks[, promoter := .I %in% Pprc1_peaks[prom, on = .(seqnames, start <= prom_end, end >= prom_start), which = TRUE,nomatch = NULL]]

Fig$Fig2F <- Pprc1_peaks[, .N, promoter] %>% 
  mutate(prop = N*100/sum(N)) %>%
  ggplot(aes(x = "", y = prop, fill = promoter)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, width = 4, show.legend = F) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c("gray", "cornflowerblue")) +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) +
  annotate("text", y = Pprc1_peaks[, mean(promoter)*100]/2, x = 1,
           label = Pprc1_peaks[, sprintf("%.1f%%", 100*mean(promoter))],
           color = "black", size = 6 / ggplot2::.pt) +
  labs(x = NULL,
       y = "% Peaks overlapping a promoter")

#=====================================================================#
#--- Figure 2G: Prop. DEGs bound  ------
#=====================================================================#
# Bar
Fig$Fig2G <- Pprc1_PROseq[diff == "Down-regulated", .(freq = .N/nrow(Pprc1_PROseq[diff == "Down-regulated"])), Pprc1_bound] %>%
  ggplot(aes(x = "", y = freq*100, fill = Pprc1_bound)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.8, width = 4, show.legend = F) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(axis.line.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("gray", "cornflowerblue")) +
  annotate("text", y = Pprc1_PROseq[diff == "Down-regulated", 100*mean(Pprc1_bound)/2], x = 1,
           label = sprintf("%.1f%%", 100*Pprc1_PROseq[diff == "Down-regulated", mean(Pprc1_bound)]) ,
           color = "black", size = 6 / ggplot2::.pt)+
  labs(x = NULL,
       fill = "Bound",
       y = "% PRO-seq down and bound")

#=====================================================================#
#--- Figure 2H: Proportion of Pprc1 direct targets with a TCT motif ----
#=====================================================================#

Fig$Fig2H <- Pprc1_PROseq[(DirectTarget), .N, hTCT] %>%
  mutate(prop = N * 100 / sum(N)) %>%
  {
    ggplot(., aes(x = "", y = prop, fill = hTCT)) +
      geom_bar(stat = "identity", position = "stack", alpha = 0.8, show.legend = FALSE) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_fill_manual(values = c("gray", "cornflowerblue")) +
      theme(axis.line.x = element_blank(),
            axis.ticks.x = element_blank()) +
      annotate("text",
               y = .$prop[.$hTCT == TRUE] / 2,
               x = 1,
               label = sprintf("%.1f%%", .$prop[.$hTCT == TRUE]),
               color = "black", size = 6 / ggplot2::.pt) +
      labs(x = NULL,
           y = "% PPRC1 targets with a TCT-motif")
  }

#=====================================================================#
#----- Figure 2I: +1 nucleotide ----
#=====================================================================#
Pprc1_PROseq[, start := tss]
Pprc1_PROseq[, plus1nt := as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, seqnames, start = start, end = start, strand = strand))] 

Fig$Fig2I <- bind_rows(
  mutate(Pprc1_PROseq, group = "all"),
  mutate(Pprc1_PROseq[(DirectTarget)], group = "Pprc1\ntarget")
) %>%
  summarise(N = dplyr::n(), .by = c(group,plus1nt)) %>%
  mutate(freq = N*100/sum(N), .by =c(group) ) %>%
  ggplot(aes(x = group,y = freq, fill = plus1nt)) +
  geom_bar(stat = "identity", position = "stack", alpha = 1) +
  scale_fill_manual(values = c("#0a9648", "#265d9a", "#f8b42c", "#d7293b")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x = NULL, y = "% +1 nucleotide", fill = NULL) +
  theme(axis.line.x = element_blank())

#=====================================================================#
#----- Figure 2J: PPRC1 TSS logo ----
#=====================================================================#
# Get -1 to +7 seq seq
Pprc1_PROseq[strand == "+", TSS_seq := as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, seqnames, start = tss-1, end = tss+6, strand = strand))]
Pprc1_PROseq[strand == "-", TSS_seq := as.character(getSeq(BSgenome.Mmusculus.UCSC.mm10, seqnames, start = tss-6, end = tss+1, strand = strand))]

Fig$Fig2J <- ggseqlogo(Pprc1_PROseq[(DirectTarget), toupper(TSS_seq)], method = "bits") +
  theme_publish() +
  theme(text = element_text(size = 6), axis.ticks.length.x = unit(0, "pt"),
        plot.margin = margin(0, 0, 0, 0)) +
  scale_x_continuous(expand = c(0,0),
                     breaks = 1:8,
                     labels = c("-1", paste0("+", 1:7))) +
  scale_y_continuous(expand = c(0,0))

#=====================================================================#
#----- Figure 2 bottom row ----
#=====================================================================#
cowplot::plot_grid(
  Fig$Fig2E,
  cowplot::plot_grid(
    Fig$Fig2F +coord_flip() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line()),
    Fig$Fig2G + coord_flip() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line()),
    Fig$Fig2H + coord_flip() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.line.x = element_line()), ncol = 1
  ),
  Fig$Fig2I,
  cowplot::plot_grid(NULL, Fig$Fig2J, NULL, rel_heights = c(0.2, 0.6, 0.2), nrow = 3), 
  nrow = 1, rel_widths = c(1.4,1, 1, 0.8), scale = 0.9, rel_heights = c(1,1,1, 0.5))

ggsave("results/Fig2_bottomPanel_2.pdf", units = "mm", width = 180, height = 42)
