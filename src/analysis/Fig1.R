#---- Libraries ------------
library(data.table)
library(tidyverse)
library(GO.db)
source("src/functions/GO_functions.R")
source("src/functions/databases.R")
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
int <- fread("Rdata/AIDscreen_Integrations.txt.gz")
int[, c("reporter","selection_type", "frame", "fraction", "rep") := tstrsplit(sampleID, "_" )]
Counts <- fread("db/count_tables/AID_screen_counts.csv")
enr <- fread("db/FC_tables/AIDscreen_CMH_enrichment.csv")
enr_EvA <- fread("db/FC_tables/AIDscreen_EF1a_vs_Atp5k_WaldTest.csv")
enr_EvA[, diff := fcase(
  padj_wald < PARAMS$padj_eva & `log2OR_EF1A-Atp5k` >= PARAMS$log2OR_cutoff_eva, "EF1a_specific",
  padj_wald < PARAMS$padj_eva & `log2OR_EF1A-Atp5k` <= -PARAMS$log2OR_cutoff_eva, "Atp5k_specific")
]

Atp5k_hits <- enr[reporter == "Atp5k" & hit, gene_id]
EF1a_hits <- enr[reporter == "EF1a" & hit, gene_id]

Atp5k_specific_hits <- enr_EvA[diff == "Atp5k_specific", gene_id]
EF1a_specific_hits <- enr_EvA[diff == "EF1a_specific", gene_id]

Shared_hits <- Atp5k_hits[(Atp5k_hits %in% EF1a_hits) & !(Atp5k_hits %in% c(Atp5k_specific_hits,EF1a_specific_hits ))]

# Gene annotations
mm10_anno <- get_ensdb("mm10")
genes <- genes(mm10_anno) %>% as.data.table()
prot_coding_genes <-genes[gene_biotype == "protein_coding", unique(gene_id)]
Intron_containing_genes <- fread("db/annotations/GRCm38.102_GeneTrapCompatibleCodingIntrons.txt.gz")$gene_id
Intronless_genes <- prot_coding_genes[!(prot_coding_genes %in% unique(Intron_containing_genes))]

## RNAseq data
RNAseq <- fread("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE99971&format=file&file=GSE99971%5FRNAseq%5Ftpm%2Etxt%2Egz")
RNAseq <- RNAseq[Geneid %in% genes$gene_id] # Remove genes which are not found in the used annotation
RNAseq[, mean_TPM := rowMeans(.SD), .SDcols = grepl("mESC", names(RNAseq))]

## GOs
GOs <- getGOs()
GOs <- GOs[, .(go_id, gene_id = ensembl_id, Ontology, description = TERM)]
GOs <- unique(GOs, by = c("gene_id", "go_id"))

# Pol II transcription
Tx_terms <- c("GO:0006366",
              "GO:0045893",
              "GO:0016591", # RNAPII holoenzyme
              "GO:0071339", # MLL1
              "GO:0006357"  #regulation of transcription by RNA polymerase II
)

Tx_neg_terms <- c("GO:0000122",
                  "GO:0045892",
                  "GO:0045943", # Pol I
                  "GO:0045945" # Pol III
)

# translation
Translation_terms <- c(
  "GO:0005840", # ribosome
  "GO:0045727", # positive regulation of translation
  "GO:0006414", # translational elongation
  "GO:0006413", # translational initiation
  "GO:0006415", # translational termination
  "GO:0030533", # triplet codon-amino acid adaptor activity
  "GO:0006418" # tRNA aminoacylation for protein translation
)
Translation_neg_terms <- "GO:0017148"

# mRNA processing
mRNA_process_terms <- c("GO:0006397", "GO:0048024")
mRNA_process_neg_terms <-c("GO:0050686", "GO:0048025")

# mRNA export
mRNA_export_terms <- c("GO:0006406", "GO:0046831")
mRNA_export_neg_terms <- c("GO:0046832")

GO_category_ids <- list(
  "transcription" = unique(
      c(
      setdiff(
        GOs[go_id %in% Tx_terms, gene_id],
        GOs[go_id %in% Tx_neg_terms, gene_id]
      ),
      "ENSMUSG00000033323") # Ctdp1
    ),
  
  "mRNA processing" = unique(
    setdiff(
    GOs[go_id %in% c(mRNA_process_terms, mRNA_export_terms), gene_id],
    GOs[go_id %in% c(mRNA_process_neg_terms, mRNA_export_neg_terms), gene_id]
    
  )),
  
  "mRNA cap\nrelated" = unique(c(
    GOs[go_id == "GO:0006370", gene_id], # 7-methylguanosine mRNA capping
    GOs[go_id == "GO:0000339", gene_id], # RNA cap binding 
    GOs[go_id == "GO:0098808", gene_id], # mRNA Cap binding
   "ENSMUSG00000029167", # Ppargc1a (PMID: 29654059)
   "ENSMUSG00000033871", # Ppargc1b paralog
    "ENSMUSG00000055491" # Pprc1 paralog
  )),
  
  "translation" = unique(setdiff(
    GOs[go_id %in% Translation_terms, gene_id],
    GOs[go_id %in% Translation_neg_terms, gene_id]
  )
  )
)

GO_anno <- as.data.table(stack(GO_category_ids))[, .(GO_category = ind, gene_id = values)]

# Set only 1 annotation per gene using hierarchical assignment.
GO_enr <- GO_anno[gene_id %in% enr$gene_id]
GO_enr[, GO_category := factor(GO_category, levels = c("mRNA cap\nrelated", "transcription", "mRNA processing", "translation"))]
GO_enr[gene_id == "ENSMUSG00000016554", GO_category := "translation"] # Eif3d as translation
GO_enr <- GO_enr[order(GO_category)][!duplicated(gene_id)]

GO_enr <- rbindlist(list(
  GO_enr,
  data.table(
    GO_category = "other",
    gene_id = setdiff(enr$gene_id, GO_anno$gene_id)
  )
))

#---- Figure 1C: Number of trapped genes -----
Prop_trapped.p <- int[dc >= PARAMS$dc_cutoff & !grepl("(AN3-12|test)", sampleID) & InFrame,.N, by = gene_id][
  RNAseq[Geneid %in% prot_coding_genes & mean_TPM >= 1,.(gene_id = Geneid)], on = "gene_id"][,
                                                                                             .(Category = factor(fcase(
                                                                                               N >2 , ">2",
                                                                                               N >0, "1-2",
                                                                                               gene_id %in% Intronless_genes, "Intronless",
                                                                                               default = "Untrapped"), levels = c("Intronless", "Untrapped", "1-2", ">2")
                                                                                             ))][, .N, by = Category] %>%
  mv_pie(., category = "Category", value = "N", fill_colors =c("gray","gray65", "lightskyblue", "lightskyblue2"), label_position = "inside") +
  theme(legend.position = "bottom")
ggsave("results/Fig1C_Screen_trappedGenes_pie.pdf", Prop_trapped.p, width = 35, height = 30, units = "mm")

#---- Figure 1E: Volcano plots for EF1a and Atp5k screens -----
Fig1E <- enr %>%
  mutate(
    padj = pmax(padj, 10^-220),
    category = case_when(
      gene_name %in% c("Gtf2f2", "Brd4") ~ "Transcription",
      gene_name %in% c("Eif2b3", "Eif3c") ~ "Translation",
      gene_name %in% c("Snrnp70", "Prpf4b") ~ "mRNA processing",
      hit ~ "Enriched",
      .default = "NS")
  ) %>%
  mutate(reporter = paste0(reporter, " (", sum(hit), " hits)"), .by = reporter) %>%
  ggplot(aes(x = log2OR, y = -log10(padj), color = category)) +
  ggrastr::geom_point_rast(alpha = 0.5, shape = 16, size = 0.5, raster.dpi = 300) +
  geom_hline(yintercept = -log10(PARAMS$padj_cutoff), linetype = "dashed", linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25) +
  scale_color_manual(
    values = c(
      "NS"        = "gray",
      "Enriched"         = "tomato",
      "Transcription"  = "steelblue",
      "Translation"    = "pink4",
      "mRNA processing"       = "orange"
    ),
    breaks = c("NS","Enriched","Transcription", "Translation", "mRNA processing")
  ) +
  facet_wrap(~ reporter, axes = "all_y") +
  labs(x = "log2 Odds ratio (GFPlow/Input)",
       y = "-log10 FDR") +
  #theme(legend.position = "none") +
  ggrepel::geom_text_repel(
    data = . %>%
      dplyr::filter(gene_name %in% c(
        "Brd4", "Gtf2f2",
        "Eif2b3", "Eif3c",
        "Snrnp70", "Prpf4b")),
    aes(label = gene_name, color = category),
    size = 2,
    show.legend = FALSE,
    min.segment.length = 0.1,       # skip very short segments (less line clutter)
    max.overlaps = Inf,
    force = 3,                       # stronger label repulsion from each other
    force_pull = 0.5,                # weaker pull toward the point (labels drift further away)
    point.padding = 0.5,             # space between label and its own point
    box.padding = 0.4,               # space between label boxes
    segment.size = 0.25,             # thinner connecting lines
    segment.alpha = 0.5,             # semi-transparent lines (less visual noise)
    bg.color = "white",              # white halo behind text for readability
    bg.r = 0.1,                      # halo radius
    nudge_x = 0.5,                   # push labels toward the right margin
    direction = "both",               # allow repulsion in both x and y
    xlim = c(-3, 0)
  ) +
  theme(legend.position = "top")
ggsave("results/Fig1E.pdf",plot =  Fig1E, width =95, height = 50, units = "mm")


#---- Figure 1F: log OR EF1a vs Atp5k screens -----
Fig1F <- enr_EvA %>%
  dplyr::filter(gene_id %in% c(Atp5k_hits, EF1a_hits))  %>%
  mutate(diff = ifelse(is.na(diff), "NS", diff)) %>%
  ggplot(aes(log2OR_Atp5k, log2OR_EF1a, color = diff)) +
  geom_point(alpha = 0.7, shape = 16, size = 0.7) +
  scale_x_continuous(limits = c(-1.2, 6)) +
  scale_y_continuous(limits = c(-1.2, 6)) +
  scale_color_manual(values = c(
    "NS" = "gray",
    "EF1a_specific" = "brown",
    "Atp5k_specific" = "steelblue"
  )) +
  geom_abline(intercept = 0, linewidth = 0.25) +
  geom_abline(intercept = c(-1,1), linetype = "dashed", linewidth = 0.25) +
  ggrepel::geom_text_repel(
    data = . %>%
      dplyr::filter(gene_name %in% c(
        "Kansl1", "Kansl2","Kansl3",
        "Pprc1", "Ncbp1", "Ncbp2")),
    aes(label = gene_name),
    size = 2,
    show.legend = FALSE,
    min.segment.length = 0,       # skip very short segments (less line clutter)
    max.overlaps = Inf,
    force = 3,                       # stronger label repulsion from each other
    force_pull = 0.1,                # weaker pull toward the point (labels drift further away)
    point.padding = 0.3,             # space between label and its own point
    box.padding = 0.4,               # space between label boxes
    segment.size = 0.5,             # thinner connecting lines
    segment.alpha = 0.5,             # semi-transparent lines (less visual noise)
    bg.color = "white",              # white halo behind text for readability
    bg.r = 0.1,                      # halo radius
    seed = 123,
    max.iter = 20000,                 # more iterations → better convergence
    max.time = 3,                     # give it more wall-clock time to resolve
    direction = "both"               # allow repulsion in both x and y
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c(0, 1)) +
  labs(x = "Atp5k log2 OR",
       y = "EF1a log2 OR",
       color = NULL)

#---- Figure 1G: GO enrichment -----
eval_genes <- unique(enr[(max_GFPlow_count > 2 | max_Input_count > 2), gene_id])
enr_SpecificHits_selectGO <- map_dfr(levels(GO_anno$GO_category), ~{
  .ft <- fisher.test(eval_genes %in% GO_enr[GO_category == .x, gene_id],
                     eval_genes %in% c(Atp5k_specific_hits, EF1a_specific_hits),
                     alternative = "two.sided")
  
  data.table(
    hit_category = "specific",
    category = .x,
    estimate = .ft$estimate,
    p.value  = .ft$p.value)
})

enr_UnionHits_selectGO <- map_dfr(levels(GO_anno$GO_category), ~{
  .ft <- fisher.test(eval_genes %in% GO_enr[GO_category == .x, gene_id],
                     eval_genes %in% c(Atp5k_hits, EF1a_hits),
                     alternative = "two.sided")
  
  data.table(
    hit_category = "union",
    category = .x,
    estimate = .ft$estimate,
    p.value  = .ft$p.value)
})

GOenr_UnionVSspecific <- rbind(enr_UnionHits_selectGO, enr_SpecificHits_selectGO)
GOenr_UnionVSspecific[, padj := p.adjust(p.value, method = "fdr"), hit_category]

# Plot
Fig1G <- GOenr_UnionVSspecific %>%
  dplyr::filter(category != "other") %>%
  mutate(category = factor(category,
                           levels = c("mRNA cap\nrelated",
                                      "transcription",
                                      "mRNA processing",
                                      "translation")),
    hit_category = factor(hit_category, levels = c("union", "specific"))) %>%
  ggplot(aes(x = hit_category,
             y = category,
             fill = log2(estimate+0.01),
             size = -log10(padj))) +
  geom_point(aes(color = padj < 0.01),
             shape = 21,
             stroke = 0.5) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  scale_fill_gradient2(low = "blue", high = "red",mid = "white", midpoint = 0) +
  scale_size_continuous(
    breaks = c(2, 5, 10),
    range  = c(0.1, 6)
  ) +
  guides(
    fill  = guide_colorbar(order = 1),
    size  = guide_legend(order = 2),
    color = guide_legend(order = 3)
  ) +
  scale_y_discrete(expand = expansion(add = 0.2)) +
  scale_x_discrete(expand = expansion(add  = 0.3)) +
  coord_cartesian(clip = "off") +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        legend.spacing.y = unit(4, "pt"),      # Removes space between legend boxes
        legend.margin = margin(t = 0, b = 0) # Removes internal padding around each legend) 
  ) +  
  labs(x = "hit category", title = "",
       y = "gene category", fill = "log2 OR", color = "padj < 0.01")

#--- Fig1H: mRNA cap related ---- 
Fig1H <- enr[gene_id %in% GO_anno[GO_category == "mRNA cap\nrelated",gene_id], .SD[sum(hit) > 0], gene_id] %>% 
  mutate(padj = pmax(padj, 10^-5)) %>%
  ggplot(aes(x = reporter, y = gene_name, fill = log2OR, size = -log10(padj))) + 
  geom_point(aes(color = padj < 0.01), shape = 21, stroke = 0.5) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  scale_fill_gradient2(low = "blue",mid =  "white", high = "red", midpoint = 0) +
  scale_size_continuous(
    breaks = c(0, 2, 5),
    range  = c(0.1, 6)
  ) +
  coord_cartesian(clip = "off") +
  guides(
    fill  = guide_colorbar(order = 1),
    size  = guide_legend(order = 2),
    color = guide_legend(order = 3)
  )  +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.spacing.y = unit(4, "pt"),      # Removes space between legend boxes
        legend.margin = margin(t = 0, b = 0) # Removes internal padding around each legend) ) 
  ) +
  scale_y_discrete(expand = expansion(add = 0.2)) +
  scale_x_discrete(expand = expansion(add  = 0.3)) +
  labs(title = "mRNA cap related hits",
       x = "Reporter", y = "Gene", fill = "log2 OR", color = "padj < 0.01")

cowplot::plot_grid(Fig1F + theme(plot.margin = margin(r = 20)),
                   Fig1G + theme(plot.margin = margin(r = 20)),
                   Fig1H,
                   rel_widths = c(0.8, 0.7,0.6), nrow = 1, align = "b") 
ggsave("results/Fig1F-H.pdf", width = 180, height = 45, units = "mm")
