#---------------------------------------------------------#
# Dependencies
#---------------------------------------------------------#
library(data.table)
devtools::load_all("/groups/stark/vloubiere/vlite/")
source("src/functions/find_elbow.R")
source("src/functions/Genomic_intervals.R")

#---------------------------------------------------------#
# Params
#---------------------------------------------------------#
extend_summit = 500

#---------------------------------------------------------#
# Parse inputs
#---------------------------------------------------------#
meta <- readRDS("Rdata/ChIP_metadata_peakCalling.rds")
meta[, c("method", "experiment", "cell_line", "antibody", "treatment") := tstrsplit(condition, "_")]
meta[, sampleID := paste(condition, replicate, sep = "_")]
meta <- unique(meta, by = "sampleID")

Pprc1.peaks <- meta[antibody == "V5" & cell_line == "Pprc1.AID" & treatment == "noTreatment" & experiment == "exp5", vlite::importBed(peaks), sampleID]

# Center the peaks on the summit
Pprc1.peaks.centered <- copy(Pprc1.peaks)
Pprc1.peaks.centered[, summit := start + peak]
Pprc1.peaks.centered[, start := summit - extend_summit+1]
Pprc1.peaks.centered[, end := summit +extend_summit]

# Collapse peaks and take the coordinates for the one with the highest score
Pprc1.peaks.centered[, contig := vlite::collapseBed(Pprc1.peaks.centered, return.idx.only = T, ignore.strand = T)]
setorder(Pprc1.peaks.centered, -score)
Ppprc1.peaks.collapse <- Pprc1.peaks.centered[, .SD[1], contig]

# Remove peaks overlapping blacklisted regions
mm10.blacklist <- fread("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz")
setnames(mm10.blacklist, c("seqnames", "start", "end", "category"))
mm10.blacklist[, start := start+1] # change to 1-based
Ppprc1.peaks.collapse <- intersectRanges(Ppprc1.peaks.collapse,mm10.blacklist, mode = "nonoverlap", ignore.strand = T )

#---------------------------------------------------------#
# Define confident peaks ----
#---------------------------------------------------------#
setorder(Ppprc1.peaks.collapse, seqnames, start, end, strand)
Ppprc1.peaks.collapse[, name := paste0("Peak_", .I)]
Ppprc1.peaks.collapse[,
                    merged := vlite::bwCoverage(.SD, "db/bw/ChIP/ChIPseq_exp5_Pprc1.AID_V5_noTreatment_merged.bw")
]

ChIP_cutoff <- find_elbow(Ppprc1.peaks.collapse$merged)
Pprc1.peaks.confident <- Ppprc1.peaks.collapse[merged >= ChIP_cutoff]

#---------------------------------------------------------#
# Save files ----
#---------------------------------------------------------#
# All peaks
rtracklayer::export.bed(Ppprc1.peaks.collapse, "db/bed/Pprc1_all_peaks_1kb_centered.bed")
fwrite(x = Ppprc1.peaks.collapse[, .(seqnames, start, end, strand, name, score, Pprc1_merged_CPM = merged)],
       file = "db/bed/Pprc1_all_peaks_1kb_centered.txt.gz", compress = "gzip")

# High confidence peaks
rtracklayer::export.bed(Pprc1.peaks.confident, "db/bed/Pprc1_confident_peaks_1kb_centered.bed")
fwrite(x = Pprc1.peaks.confident[, .(seqnames, start, end, strand, name, score, Pprc1_merged_CPM = merged)],
       file = "db/bed/Pprc1_confident_peaks_1kb_centered.txt")
