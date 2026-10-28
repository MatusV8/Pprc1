#-----------------------------------#
# Libraries
#-----------------------------------#
library(data.table)

#-----------------------------------#
# Input
#-----------------------------------#
## --- PROseq ----
Pprc1_PROseq <- fread("db/FC_tables/PROseq/Pprc1.depl_transcript_spikeIn_norm_DESeq2_FC.txt")
Pprc1_PROseq[, c("gene_id", "gene_name", "coord") := tstrsplit(gene_id, "__")]
Pprc1_PROseq[, c("seqnames", "tss", "strand") := tstrsplit(coord, "_")]
Pprc1_PROseq[, tss := as.numeric(tss)]

## --- ChIP ----
Pprc1_peaks <- fread("db/bed/Pprc1_confident_peaks_1kb_centered.txt")

# Add Peaks to PROSeq
Pprc1_PROseq[, Pprc1_bound := .I %in% Pprc1_PROseq[Pprc1_peaks,
                                                  on = .(seqnames, tss >= start, tss <= end),
                                                  which = TRUE, nomatch = NULL]]
Pprc1_PROseq[, DirectTarget := diff == "Down-regulated" & Pprc1_bound]

# Add TCT annotations to PROseq table
prom <- fread("db/annotations/GRCm38.v102.EPDnew_promoters.txt.gz")
tct_tss <- unique(prom[!is.na(hTCT) & hTCT != "no",
                       .(seqnames, strand, tss)])

Pprc1_PROseq[, hTCT := FALSE]
Pprc1_PROseq[tct_tss, hTCT := TRUE, on = .(seqnames, strand, tss)]

#-----------------------------------#
# Export table
#-----------------------------------#
fwrite(Pprc1_PROseq[, .(
  gene_id,
  gene_name,
  seqnames,
  tss,
  strand,
  baseMean,
  log2FoldChange,
  lfcSE,
  stat,
  pvalue,
  padj,
  diff,
  Pprc1_bound,
  DirectTarget,
  hTCT
)], "db/FC_tables/Pprc1_PROseq_DESeq2_FC_Pprc1Bound_TCT.csv")