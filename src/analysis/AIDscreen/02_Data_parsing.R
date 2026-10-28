setwd("/groups/stark/vojtek/Bioinfo_project/PPRC1_paper/Pprc1_git/")
## Libraries ----
library(tidyverse)
library(data.table)
#===================================#
# Variables ----
#===================================#

META_PATH = "Rdata/AIDscreen_processed_metadata.rds"

PARAMS <- list(
  dc_cutoff = 10
)

#===================================#
# Compile annotated insertions ----
#===================================#
# Read metadata
META <- readRDS(META_PATH)

Annotated_insertions <- META[, {
  fread(paste0("db/annotated_integrations/Pprc1_AIDscreen/", .BY, "_filt_anno.txt.gz"))
}, by = sampleID ]

# Export data
data.table::fwrite(Annotated_insertions, file = "Rdata/AIDscreen_Integrations.txt.gz", sep = "\t", compress = "gzip")

#===================================#
# DC filter and count  ----
#===================================#

# DC filtering and keep only compatible introns
cor_int <- Annotated_insertions[(InFrame) & dc >= PARAMS$dc_cutoff]

# Split by experiment
intTest_int <- cor_int[grepl("test", sampleID)]
screen_int <- cor_int[!grepl("test", sampleID)]
screen_int[, sampleID := gsub("Intein_", "", sampleID)]

# create a count table for each sample
count.tbl <- copy(screen_int)
count.tbl <- unique(count.tbl, by = c("sampleID", "seqnames", "start", "strand", "gene_id", "gene_name")) # Make sure that no integrations are duplicated
count.tbl <- dcast(count.tbl, gene_id + gene_name ~ sampleID, fun.aggregate = length, fill = 0)

# create a count table for merged replicates, merge inputs together
count.tbl.merg <- copy(screen_int)
count.tbl.merg[, sampleID := gsub("_rep[1-9]", "", sampleID)]
count.tbl.merg[, sampleID := gsub(".*_([1-3])_Input", "Input_\\1", sampleID)]
count.tbl.merg <- unique(count.tbl.merg, by = c("sampleID", "seqnames", "start", "strand", "gene_id")) # Make sure that no integrations are duplicated
count.tbl.merg <- dcast(count.tbl.merg, gene_id +gene_name ~ sampleID, fun.aggregate = length, fill = 0)

# Save count table
if(!dir.exists("db/count_tables")){
  dir.create("db/count_tables")
}

fwrite(count.tbl, "db/count_tables/AID_screen_counts.csv")
fwrite(count.tbl.merg, "db/count_tables/AID_screen_counts_merged.csv")
