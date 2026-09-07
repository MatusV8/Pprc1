#================================#
# Setup ----
#================================#
# Libraries
library(data.table)
library(tidyverse)
source("/groups/stark/vojtek/matools/mv_OddsRatio.R")
source("src/functions/Fisher_grouped.R")
source("src/functions/Rowwise_fisher_test.R")
source("/groups/stark/vojtek/matools/mv_hm.R")
source("/groups/stark/vojtek/matools/mv_pie.R")
library(ComplexHeatmap)

PARAMS <- list(
  dc_cutoff = 10,
  gene_count_cutoff = 3,
  padj_cutoff = 0.05,
  l2OR_cutoff = 0
)


#================================#
# Parse input ----
#================================#
Int <- fread("Rdata/AIDscreen_Integrations.csv")
Int[, c("reporter","selection_type", "frame", "fraction", "rep") := tstrsplit(sampleID, "_" )]

Counts <- fread("db/count_tables/AID_screen_counts.csv")

#================================#
# Correlation between samples ----
#================================#

# Reshape data for heatmap format
corr_mat <- Counts %>%
  mutate_if(is.numeric, ~log2(.x+1)) %>%
  dplyr::select(-c(gene_id, gene_name)) %>%
  cor(use = "everything") 

GFPlow_corr_mat <- corr_mat[grepl("GFPlow", rownames(corr_mat)), grepl("GFPlow", colnames(corr_mat))]
colnames(GFPlow_corr_mat) <- gsub("_GFPlow", "", colnames(GFPlow_corr_mat))

Input_corr_mat <- corr_mat[grepl("Input", rownames(corr_mat)), grepl("Input", colnames(corr_mat))]
colnames(Input_corr_mat) <- gsub("_Input", "", colnames(Input_corr_mat))

GFPlow_cor_hm <- mv_hm(GFPlow_corr_mat, col_labels = c("Promoter", "frame", "replicate"), col_title = "GFPlow samples",colors = list(
  Promoter = c("#66C2A5", "#FC8D62"),
  frame = c("#E41A1C", "#377EB8", "#4DAF4A"),
  replicate = c("#A65628", "#F781BF")
))
for (x in names(GFPlow_cor_hm@top_annotation@anno_list)) {
  GFPlow_cor_hm@top_annotation@anno_list[[x]]@name_param$show <- FALSE
}

Input_cor_hm <- mv_hm(Input_corr_mat, col_labels = c("Promoter", "frame", "replicate"), col_title = "Input samples",colors = list(
  Promoter = c("#66C2A5", "#FC8D62"),
  frame = c("#E41A1C", "#377EB8", "#4DAF4A"),
  replicate = c("#A65628", "#F781BF")
))

pdf("results/Screen_replicates_PCC.pdf", width = 12, height = 5)

grid::grid.newpage()

# Left heatmap
pushViewport(viewport(x = 0.25, y = 0.5, width = 0.5, height = 1))
draw(Input_cor_hm, newpage = FALSE)
popViewport()

# Right heatmap
pushViewport(viewport(x = 0.75, y = 0.5, width = 0.5, height = 1))
draw(GFPlow_cor_hm, newpage = FALSE)
popViewport()

dev.off()
#================================#
# Number of trapped genes  ----
#================================#

# Intronless genes
source("/groups/stark/vojtek/matools/get_EsDb.R")
mm10_anno <- get_ensdb("mm10")

prot_coding_genes <- as.data.table(genes(mm10_anno, 
                                         filter = ~ gene_biotype == "protein_coding"))[, unique(gene_id)]
Intron_containing_genes <- rtracklayer::import("/groups/stark/vojtek/Resources/01_annotations/mm10_GRCm38_v102/Mus_musculus.GRCm38.102_CodingIntrons.gtf.gz")$gene_id
Intronless_genes <- prot_coding_genes[!(prot_coding_genes %in% unique(Intron_containing_genes))]

# Get IDs of trapped genes
Trapped_genes <- Int[dc >= PARAMS$dc_cutoff & !grepl("AN3-12", sampleID) & CompatibleIntron, unique(gene_id)] 

## Get a set of expressed genes from RNAseq data
RNAseq <- fread("/groups/stark/nemcko/work/ORFtrap/general/GSE99971_RNAseq_tpm.txt.gz")
RNAseq[, mean_TPM := rowMeans(.SD), .SDcols = grepl("mESC", names(RNAseq))]
RNAseq[, quantile_group := {
  q <- quantile(mean_TPM[mean_TPM >= 1], probs = c(0.25, 0.5, 0.75, 1))
  fcase(
    mean_TPM < 1,              0L,
    mean_TPM <= q[1],          1L,
    mean_TPM <= q[2],          2L,
    mean_TPM <= q[3],          3L,
    mean_TPM <= q[4],          4L
  )
}]

# Calculate fraction of trapped genes per quantile
trapped_ExpQuant <- RNAseq[, {
  Expressed_genes <- unique(Geneid)
  taggable_genes <- prot_coding_genes[(prot_coding_genes %in% Expressed_genes)]
  
  data.table(
   trapped = taggable_genes %in% Trapped_genes)
}, by = quantile_group]

# Trapped genes per quantile of expressed
trapped_ExpQuant %>%
  ggplot(aes(x = quantile_group, fill = trapped)) +
  theme_minimal() +
  scale_fill_manual(values = c("gray", "tomato")) +
  labs(y = "Fraction of trapped genes",
       x = "Expression quantile") +
  scale_y_continuous(expand = c(0,0 )) +
  geom_bar(position = "fill", stat = "count")
ggsave("results/Screen_trappedGenes_byQuantile.pdf", width = 3, height = 3)


# Pie chart for TPM > 1
Prop_trapped.p <- RNAseq[Geneid %in% prot_coding_genes & quantile_group != 0,
                         .(Category = fcase(
                           Geneid %in% Trapped_genes,   "Trapped",
                           Geneid %in% Intronless_genes, "Intronless",
                           default = "Untrapped"
                         ))][, .N, by = Category] %>%
  mv_pie(., category = "Category", value = "N", fill_colors =c("gray", "#4DBBD5", "gray65"), label_position = "inside")


# Pie chart for top 50% expressed
Prop_trapped_q05.p <- RNAseq[Geneid %in% prot_coding_genes & quantile_group %in% c(3,4),
                         .(Category = fcase(
                           Geneid %in% Trapped_genes,   "Trapped",
                           Geneid %in% Intronless_genes, "Intronless",
                           default = "Untrapped"
                         ))][, .N, by = Category] %>%
  mv_pie(., category = "Category", value = "N", fill_colors =c("gray", "#4DBBD5", "gray65"), label_position = "outside")

#================================#
# Prop correct  ----
#================================#
Int %>%
  dplyr::filter(dc >= PARAMS$dc_cutoff, !grepl("(AN3-12|test)", sampleID)) %>%
  mutate(Category = case_when(
    !is.na(phase) ~ paste0("Frame ", phase+1),
    is.na(phase) & IntronOrient ~ "Wrong frame",
    Intron & !IntronOrient ~ "Wrong orientation",
    Intragenic & !Intron ~ "Not intron",
    !Intragenic ~ "Not gene")) %>%
  dplyr::select(sampleID, Category, frame) %>%
  dplyr::summarise(N = dplyr::n(), frame = unique(frame), .by = c(sampleID, Category)) %>%
  mutate(freq = N / sum(N), .by = sampleID) %>%
  mutate(frame = paste0("ORF ", frame)) %>%
  ggplot(aes(y = reorder(sampleID, -freq), x = freq, fill = as.character(Category))) +
  geom_bar(position="stack", stat="identity") +
  scale_x_continuous(expand = c(0,0), breaks = seq(from = 0, to = 1, by = 0.1)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(fill = "Annotation",
       x = "Proportion of unique integrations",
       y = element_blank()) +
  theme_bw() +
  geom_text(aes(label = round(freq,2 )), position = position_stack(vjust = 0.5)) +
  facet_wrap(~ frame, ncol = 1, scales = "free_y")
#ggsave("results/TEST_2reps_IntegrationPosition.pdf", width =9, height = 6)

# Test samples
## Fraction plot
Int %>%
  dplyr::filter(dc >= PARAMS$dc_cutoff, grepl("AN3-12", sampleID)) %>%
  mutate(Category = case_when(
    !is.na(phase) ~ paste0("Frame ", phase+1),
    is.na(phase) & IntronOrient ~ "Wrong frame",
    Intron & !IntronOrient ~ "Wrong orientation",
    Intragenic & !Intron ~ "Not intron",
    !Intragenic ~ "Not gene"),
    sampleID = gsub(".*_(Intein|Classic).*_(rep[1-9])", "\\1_\\2", sampleID)) %>%
  dplyr::select(sampleID, Category, frame) %>%
  dplyr::summarise(N = dplyr::n(), frame = unique(frame), .by = c(sampleID, Category)) %>%
  mutate(freq = N / sum(N), .by = sampleID) %>%
  mutate(frame = paste0("ORF ", frame)) %>%
  ggplot(aes(y = sampleID, x = freq, fill = as.character(Category))) +
  geom_bar(position="stack", stat="identity") +
  scale_x_continuous(expand = c(0,0), breaks = seq(from = 0, to = 1, by = 0.1)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(fill = "Annotation",
       x = "Proportion of unique integrations",
       y = element_blank()) +
  theme_bw() +
  geom_text(aes(label = round(freq,2 )), position = position_stack(vjust = 0.5)) 

## simple plot
Int[dc >= PARAMS$dc_cutoff & grepl("AN3-12", sampleID),
    .(
      selection = gsub(".*(Intein|Classic).*", "\\1", .BY),
      freq = sum(CompatibleIntron, na.rm = T)/.N
      ), by = sampleID] %>%
  ggplot(aes(x = selection, y = freq)) +
  theme_classic() +
  geom_jitter(width = 0.1, alpha = 0.7, size = 3) +
  labs(y = "Corret integration / all",
       x = NULL)

