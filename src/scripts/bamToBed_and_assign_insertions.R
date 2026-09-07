#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is correct number of args: if not, return an error
if (length(args)!=9) {
  stop("Please specify:\n
       [required] 1/ Collapsed bam file\n
       [required] 2/ Path to a gtf file containing genome annotation\n
       [required] 3/ Path to a txt file containing compatible introns used for assignment\n
       [required] 4/ Phase of the sample (1, 2 or 3)\n
       [required] 5/ Output bed file (.bed)\n
       [required] 6/ Output Slippage filtered bed file (.bed)\n
       [required] 7/ Output annotated bed file (.txt)\n
       [required] 8/ Output Compatible only bed file (.bed)\n
       [required] 9/ Output filtering and annotation stats (.txt)\n")

}

#========================================#
# Packages
#========================================#
require(tidyverse)
require(Rsamtools)
require(rtracklayer)
require(GenomicRanges)
require(data.table)
source("src/functions/CollapseBam_functions.R")
source("src/functions/Genomic_intervals.R")
source("src/functions/annotate_insertions.R")

#========================================#
# Input and output paths ####
#========================================#

# Input
bam <- args[1]
anno_path <- args[2]
introns_path <- args[3]
Phase <- as.numeric(args[4])-1 # Assuming that Phase takes values 1 to 3 in the metadata.txt

bam <- "/groups/stark/vojtek/Bioinfo_project/PPRC1_paper/db/bam_unique/AIDscreen/EF1a_Intein_1_Input_rep1_q30_unique.bam"
anno_path <- "/groups/stark/vojtek/Resources/01_annotations/mm10_GRCm38_v102/Mus_musculus.GRCm38.102.gtf.gz"
introns_path <- "/groups/stark/vojtek/Bioinfo_project/PPRC1_paper/Pprc1_git/db/annotations/GRCm38.102_GeneTrapCompatibleCodingIntrons.txt.gz"
Phase <- 0

# Output
o.bed <- args[5]
o.filtBed <- args[6]
o.AnnoBed <- args[7]
o.CompatibleBed <- args[8]
o.Stats <- args[9]

#========================================#
# Parse files ####
#========================================#
# Annotations
Annotations_mm10 <- rtracklayer::import(anno_path, format  = "gtf") %>% as.data.table()
Annotations_mm10[ , transcript_support_level := as.numeric(gsub(" .*", "", transcript_support_level))]
Annotations_mm10 <- Annotations_mm10[
  seqnames %in% c(1:19, "X", "Y") & # only standard chromosomes
    tag == "basic" & # keep only complete annotations
    transcript_support_level %in% c(1,2) # keep only transcripts with the highest support level
]

Compatible_introns <- fread(introns_path)
Compatible_introns <- Compatible_introns[seqnames %in% c(1:19, "X", "Y")]

# Filter reads
Raw <- readBam(bam, tags = "dc")
Collapsed <- collapseReadStarts(Raw)
Filtered <- SlippageFilter(Collapsed, hamming_dist = 3, mingap =1,fragment_size = 30, seq_length = 15)

# Get positions of insertions (1-based).
# Note: strand needs to be reversed if iPCR cassette is in a reverse orientation
ins <- copy(Filtered[, .(seqnames, pos, strand, width, dc)])
ins[, start := fifelse(strand == "+", pos-1, pos+width)]
ins[, strand := fifelse(strand == "+", "-", "+")]
ins <- ins[, .(seqnames, start, end = start, strand, width, dc )]

# Annotate
Annotated <- annotateInsertions(bed.dt = ins,
                                gtf = Annotations_mm10,
                                compatible_introns = Compatible_introns,
                                Phase = Phase)
#========================================#
# Export files ####
#========================================#
# Annotation file
write_tsv(Annotated, file = o.AnnoBed, col_names = T)

# Export stats
tibble(
  "N_reads" = sum(Raw$dc),
  "N_integrations" = nrow(Raw),
  "N_collapsed" = nrow(Collapsed),
  "N_collapsed_filtered" = nrow(Filtered),
  "N_collapsed_filtered_annotated" = nrow(Annotated)) %>%
  write_tsv(., file = o.Stats, col_names = T)
