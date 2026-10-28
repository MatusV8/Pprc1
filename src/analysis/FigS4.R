#--------------------------------------------------#
# Libraries
#--------------------------------------------------#
devtools::load_all("/groups/stark/vloubiere/vlite/")
library(BSgenome.Mmusculus.UCSC.mm10)
library(tidyverse)
library(Hmisc)
library(ggpointdensity)
library(ggrastr)
library(viridis)
source("src/functions/ggplot2_theme.R")
source("src/functions/Genomic_intervals.R")

#--------------------------------------------------#
# Parameters ----
#--------------------------------------------------#
# IP-MS thresholds
fc_thresh  <- 4
fdr_thresh <- 0.005

#--------------------------------------------------#
# Parse Input ----
#--------------------------------------------------#
# Blacklisted regions
mm10.blacklist <- fread("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz",
                        col.names = c("seqnames", "start", "end", "name"))

#Pprc1 peaks
exp6_beds <- list.files(path = "../db/bed/ChIP_peaks/", pattern = ".*V5.*exp6.*narrowPeak", full.names = T)

exp6_peaks <- lapply(exp6_beds, function(x) {
  
  rtracklayer::import(x, format = "narrowPeak") |>
    as.data.table()
  
  })
exp6_peaks <- rbindlist(exp6_peaks)
exp6_CollapsedPeaks <- vlite::collapseBed(exp6_peaks) # Collapse
exp6_CollapsedPeaks <- intersectRanges(exp6_CollapsedPeaks, mm10.blacklist, mode = "nonoverlap") # remove blacklisted regions

# Taf1 peaks
Taf1_files <- list.files("/groups/stark/sabath/ChIP/251028_ChIP_Hcfc1Depletion2/db/peaks/ChIP/",
                         pattern = "Taf1.*noTreatment.*exp2.*narrowPeak", full.names = T)


Taf1_peaks <- lapply(Taf1_files, function(x) {
  
  rtracklayer::import(x, format = "narrowPeak") |>
    as.data.table()
  
})

Taf1_peaks <- rbindlist(Taf1_peaks)
Taf1_CollapsedPeaks <- vlite::collapseBed(Taf1_peaks) # collapse
Taf1_CollapsedPeaks <- intersectRanges(Taf1_CollapsedPeaks, mm10.blacklist,
                                       mode = "nonoverlap") # remove blacklisted regions

# Bigwigs
bws <- list(
  "Taf1_rep1" = "/groups/stark/sabath/ChIP/251028_ChIP_Hcfc1Depletion2/db/bw/ChIP/Taf1_FNR94_noTreatment_exp2_rep1.bw",
  "Taf1_rep2" = "/groups/stark/sabath/ChIP/251028_ChIP_Hcfc1Depletion2/db/bw/ChIP/Taf1_FNR94_noTreatment_exp2_rep2.bw",
  "noTreatment_rep1" = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_noTreatment_exp6_rep1.bw",
  "noTreatment_rep2" = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_noTreatment_exp6_rep2.bw",
  "Flavopiridol_rep1" = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_Flavopiridol_exp6_rep1.bw",
  "Flavopiridol_rep2" = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_Flavopiridol_exp6_rep2.bw",
  "RNAse_rep1" = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_RNase_exp6_rep1.bw",
  "RNAse_rep2" = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_RNase_exp6_rep2.bw",
  "Triptolide_rep1" = "../db/bw/ChIP/ChIPseq_Pprc1.AID_V5_Triptolide_exp6_rep1.bw"
)

# Check if files exist
all(sapply(bws, file.exists))

# Quantify signal over peaks
treatment_names <- grep("Taf1", names(bws), invert = T, value = T)
exp6_CollapsedPeaks[, (treatment_names) := lapply(bws[treatment_names], function(x) vlite::bwCoverage(exp6_CollapsedPeaks, x))]

Taf1_CollapsedPeaks[, "Taf1_rep1" := vlite::bwCoverage(Taf1_CollapsedPeaks, bws$Taf1_rep1)]
Taf1_CollapsedPeaks[, "Taf1_rep2" := vlite::bwCoverage(Taf1_CollapsedPeaks, bws$Taf1_rep2)]

## ------- IP-MS ----------------------#
# Import proteomics results
dat <- fread("/groups/stark/sabath/MS/260117_IPMS_PPRC1dAD_triplicates/db/volcano_20241011_E3_NEO7_Sabath_Stark_IMP_ID1828_onbead_dig_8x_direct_normsum_limma_Proteins.txt")
dat <- dat[, .(
  gene_name = `Gene Symbol`, 
  mouse_ensembl = `Ensembl Gene ID`,
  mouse_uniprot = Accession, 
  
  PPRC1_wt_rep1 = `Imputed Normalised apQuant Area 20241011_E3_NEO7_Sabath_Stark_IMP_ID1742_onbead_dig_5per_repeat_control_FLAG_GFP_PPRC_1_rep1`,
  PPRC1_wt_rep2 = `Imputed Normalised apQuant Area 20241011_E3_NEO7_Sabath_Stark_IMP_ID1742_onbead_dig_5per_repeat_control_FLAG_GFP_PPRC_1_rep2`,
  PPRC1_wt_rep3 = `Imputed Normalised apQuant Area 20241011_E3_NEO7_Sabath_Stark_IMP_ID1828_onbead_dig_10per_control_FLAG_GFP_PPRC_1_rep3`,
  
  control_rep1 = `Imputed Normalised apQuant Area 20241011_E3_NEO7_Sabath_Stark_IMP_ID1742_onbead_dig_5per_repeat_control_FLAG_GFP_rep1`,
  control_rep2 = `Imputed Normalised apQuant Area 20241011_E3_NEO7_Sabath_Stark_IMP_ID1742_onbead_dig_5per_repeat_control_FLAG_GFP_rep2`,
  control_rep3 = `Imputed Normalised apQuant Area 20241011_E3_NEO7_Sabath_Stark_IMP_ID1828_onbead_dig_10per_control_FLAG_GFP_rep3`,
  
  log2FC_wt = `log2Fold change sample  control`, 
  FDR_wt = `Adjusted p-Value normalised sample  control`
)]

dat <- dat[!grep("cont_", mouse_uniprot)]

#--------------------------------------------------#
# Figure S4A: PCC Taf1 ChIP ----
#--------------------------------------------------#
Taf1_CollapsedPeaks %>%
  ggplot(aes(log10(Taf1_rep1+0.1), log10(Taf1_rep2+0.1))) +
  ggpointdensity::geom_pointdensity(size = 0.2)+
  viridis::scale_color_viridis(name = "Density") +
  geom_smooth(method = "lm", se = FALSE, alpha = 1, linewidth = 0.5) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1,
           label = paste0("Pearson's r = ", round(cor(log10(Taf1_CollapsedPeaks$Taf1_rep1 + 0.1),
                                                      log10(Taf1_CollapsedPeaks$Taf1_rep2 +0.1)), 2)), size = 6 / .pt ) +
  labs(x = "rep1 log10 CPM",
       y = "rep2 log10 CPM",
       title = "Taf1 ChIP-seq")
ggsave("results/FigS4A.pdf", width = 50, height = 40, units = "mm")  

#--------------------------------------------------#
# Figure S4B: Pprc1 treatment PCC ----
#--------------------------------------------------#
treatment_pcc <- (as.matrix(exp6_CollapsedPeaks[, .SD,
                                           .SDcols = grep("Pprc1|Taf1", names(bws), invert = TRUE, value = TRUE)]) +0.1) |>
  log10() |>
  Hmisc::varclus(similarity = "pearson", trans = "none")

pdf("results/FigS4B.pdf", width = 35/25.4, height = 40/25.4, pointsize = 6)

par(
  mar = c(3, 5, 2, 2),
  tcl = -0.2,
  mgp = c(3, 0.5, 0),
  las = 1,
  cex = 0.8, lwd = 0.5
)

plot(
  treatment_pcc,
  lwd = 0.5,
  hang = -1,
  xlab = NA,
  ylab = "PCC",
  las = 1,        # rotate leaf labels vertically so they fit narrow spacing
  cex = 1       # shrink leaf label font relative to pointsize = 6
)

labs <- factor(gsub("_rep[12]", "", colnames(treatment_pcc$sim)))
Cc <- rainbow(nlevels(labs))[labs]

points(
  seq_along(treatment_pcc$hclust$order),
  rep(0, length(treatment_pcc$hclust$order)),
  col = Cc[treatment_pcc$hclust$order],
  pch = 19,
  cex = 0.8
)

dev.off()

#--------------------------------------------------#
# Figure S4C: Pprc1 IP-MS PCC ----
#--------------------------------------------------#
ip <- dat[, .(PPRC1_wt_rep1, PPRC1_wt_rep2, PPRC1_wt_rep3,
              control_rep1, control_rep2, control_rep3)]

pdf("results/FigS4C.pdf", width = 35/25.4, height = 40/25.4, pointsize = 6)
par(
  mar = c(3, 5, 2, 2),
  tcl = -0.2,
  mgp = c(3, 0.5, 0),
  las = 1,
  cex = 0.8,lwd = 0.5)

ip[, {
  # Import normalized counts
  counts <- as.matrix(.SD)
  samples <- factor(gsub("^(.*)_rep.*", "\\1", colnames(counts)))
  Cc <- viridis::viridis(length(levels(samples)))
  Cc <- Cc[samples]
  
  # Prepare matrix and log-transform
  mat <- as.matrix(counts)
  mat <- log10(mat + 1)
  
  # Compute correlation clustering
  cor <- varclus(mat, similarity = "pearson", trans = "none")
  
  # Plot dendrogram
  plot(
    cor,
    lwd = 0.5,
    hang = -1,
    xlab = NA,
    las = 1,
    ylab = "PCC"
  )
  
  
  # Add colored points for each sample
  points(
    1:ncol(mat),
    rep(0, ncol(mat)),
    col = Cc[cor$hclust$order],
    pch = 19,
    cex = 0.8
  )
},]

dev.off()


#--------------------------------------------------#
# Figure S4D: Pprc1 IP-MS volcano ----
#--------------------------------------------------#
gene_groups <- list(
  "PPRC1"   = c("Q6NZN1;PPRC1"),   
  "ATAC"    = c("Q99LQ1", "Q8CHV6", "Q8CID0", "Q9DA08", "Q9JHD2", "Q3TUF7", "Q8R0L9", "Q6KAQ7" ),
  "PAF1C" = c("Q62018", "Q8K2T8","Q5XJE5", "Q8JZM7"),
  "CBC" = c("Q3UYV9", "Q9CQ49"),
  "NELF" = c("Q8BG30",  "Q8C4Y3", "P19426"),
  "HCFC1" = c("Q61191")
)

group_colors <- c(
  "PPRC1"   = "#1B9E77",
  "ATAC"    = "#D95F02",
  "PAF1C"   = "#7570B3",
  "CBC"     = "#1F78B4",
  "NELF"    = "#E7298A", 
  "HCFC1"   = "#00C853",
  "Other" = "grey80"
)

Uniprot_lookup <- setNames(rep(names(gene_groups), lengths(gene_groups)), unlist(gene_groups))
dat[, group := Uniprot_lookup[mouse_uniprot]]
dat[is.na(group), group := "Other"]
dat[, group := factor(group,
                          levels = rev(names(group_colors)))]

dat[, hit := abs(log2FC_wt) >= fc_thresh & FDR_wt <= fdr_thresh]

dat %>%
  ggplot(aes(log2FC_wt, -log10(FDR_wt), color = group)) +
  geom_point(data = . %>% dplyr::filter(group == "Other"), shape = 16, alpha = 0.5, size =0.4) +
  geom_point(data = . %>% dplyr::filter(group != "Other"), shape = 16, alpha = 1, size =1) +
  ggrepel::geom_text_repel(
    data = . %>% dplyr::filter(group != "Other", FDR_wt < 10^-2.5),
    aes(label = gene_name),
    max.overlaps = 30,
    size = 2,
    bg.color = "white",   # outline/halo color around the text
    bg.r = 0.15,           # thickness of the outline
    segment.size = 0.3,
    box.padding = 0.1,
    force_pull = 1, max.iter = 10000,max.time = 2, force = 1,
    show.legend = FALSE
  ) +
  #scale_y_continuous(limits = c(0, 5)) +
  theme(legend.position = "top",
        legend.key.spacing.x = unit(3, "pt"),
        legend.text = element_text(margin = margin(l = -1, unit = "pt"))) +
  labs( x = expression(log[2]("Fold change")),
        y = expression(-log[10]("FDR")),
        color = NULL) +
  scale_color_manual(values = group_colors)
ggsave("results/FigS4D.pdf", width = 65, height = 48, units = "mm")
