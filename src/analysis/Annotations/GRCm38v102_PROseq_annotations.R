# Libraries ----
library(data.table)
devtools::load_all("/groups/stark/vloubiere/vlite")
source("src/functions/databases.R")

#===================================================#
# Settings ----
#===================================================#
Promoter_up = 250
Promoter_down = 250
gene_extension = 250
TES_distance = 500

#===================================================#
# Parse data ----
#===================================================#
EPD <- fread("db/annotations/GRCm38.v102.EPDnew_promoters.txt.gz")
EPD[,start := tss]
EPD[,end := tss]
EPD[strand == "+", CAGE_max_TSS := vlite::bwCoverage(.SD, "db/public/ctssMaxCounts_FANTOM5_CAGE_TSS_mm10.fwd.bw")]
EPD[strand == "-", CAGE_max_TSS := vlite::bwCoverage(.SD, "db/public/ctssMaxCounts_FANTOM5_CAGE_TSS_mm10.rev.bw")]
EPD[is.na(CAGE_max_TSS), CAGE_max_TSS:= 0]

# Get ensembl DB
ensdb <- get_ensdb("mm10")
ensdb_genes <- as.data.table(genes(ensdb))
ensdb_genes[,seqnames := paste0("chr", seqnames)]
ensdb_genes[,entrezid := NULL]
ensdb_genes[, seqnames := as.character(seqnames)]
ensdb_genes[, gene_version := as.numeric(gsub(".*\\.(.*)$", "\\1", gene_id_version))]

# Add TES information
EPD[strand == "+", `:=`(start = tss, end = ensdb_genes[.SD, on = "gene_id", end])]
EPD[strand == "-", `:=`(start = ensdb_genes[.SD, on = "gene_id", start], end = tss)]
EPD[, width := end-start+1]

#==========================================================#
# Protein coding and lincRNA only ----
#==========================================================#
EPD_ProtLincRNA_reduced <- EPD[gene_biotype %in% c("protein_coding", "lincRNA")][
  order(-CAGE_max_TSS),
  .SD[1, .(seqnames, start, end, strand, cluster.id = paste0(gene_id, "__", gene_name, "__", paste(seqnames,tss,strand, sep = "_")))],
  gene_name]

# Save data
purrr::walk(c("transcript", "geneBody", "promoter"), ~{
  
  .dt <- copy(EPD_ProtLincRNA_reduced)
  
  # Resize based on the parameters around TSS
  if(.x == "transcript"){
    .dt[, `:=` (start = start-gene_extension, end = end +gene_extension)] # Extend gene
    .dt[, width:= end-start+1]
  } else if(.x == "geneBody"){
    .dt[strand=="+", `:=` (start = start+Promoter_down, end = end -TES_distance)]
    .dt[strand=="-", `:=` (start = start+TES_distance, end = end - Promoter_down)]
    .dt[, width:= end-start+1]
  } else if(.x == "promoter"){
    .dt[strand=="+", `:=` (start = start-Promoter_up, end = start+ Promoter_down)]
    .dt[strand=="-", `:=` (start = end-Promoter_down, end = end+Promoter_up)]
    .dt[, width:= end-start+1]
  } else {
      warning("Category not found. Maybe misspelled?")
  }
 
  .f_name = paste0("db/annotations/mm10.GRCm38.v102.EPDnew.ProtCodingLincRNAonly.CAGEcorrected.", Promoter_up, "up", Promoter_down, "down_", .x, ".rds")
  saveRDS(object = .dt,
          file = .f_name)
  
})
