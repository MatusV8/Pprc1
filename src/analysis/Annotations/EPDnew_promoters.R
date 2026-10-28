#---------------------------------#
# Dependencies ------
#---------------------------------#
devtools::load_all("/groups/stark/vloubiere/vlite/")
source("src/functions/databases.R")
source("src/functions/parse_EPDnew.R")
library(data.table)
library(tidyverse)

#---------------------------------#
# Retrieve EPD ------
#---------------------------------#
dat.file <- "db/annotations/Mm_EPDnew.txt.gz"
if(!file.exists(dat.file)){
  # Download file
  temp_file = "db/annotations/EPDnew.dat"
  download.file(url = "https://epd.expasy.org/ftp/epdnew/M_musculus/current/Mm_EPDnew.dat",
                destfile = temp_file)
  
  # Parse file
  EPD <- read_epd_dat(temp_file)
  
  # Save the file
  fwrite(EPD, file = dat.file, compress = "gzip")
  
  # Remove the raw file
  file.remove(temp_file)
  
}
  

ElemenT_file <- "db/annotations/mm10_EPDnew_minPlus100__start101smooth10.txt"
if(!file.exists(ElemenT_file)){
  ElemeNT_ulr <- "https://www.juven-gershonlab.org/wp-content/uploads/2023/11/CoreFiles_vCore.zip"
  download.file(ElemeNT_ulr, destfile = "db/annotations/ElemeNT_vCore.zip")
  unzip("db/annotations/ElemeNT_vCore.zip", exdir = "db/annotations")
  file.copy("db/annotations/CoreFiles/mm10_EPDnew_minPlus100__start101smooth10.txt", ElemenT_file)
  
  # remove unnecessary files
  unlink("db/annotations/CoreFiles/", recursive = TRUE)
  file.remove("db/annotations/ElemeNT_vCore.zip")
}


#---------------------------------#
#---- Parse data ----
#---------------------------------#
motifs <- fread(ElemenT_file)
motifs[name == "S-adenosylmethionin", name := "S-adenosylmethioninedecarboxylase246_1"] # Fix so the names match
setnames(motifs, "name", "entry_name")

EPD <- fread(dat.file)
EPD_clean <- EPD[!is.na(chr),.(
  entry_name,
  promoter_id = tstrsplit(entry_name, "_", keep = 2)[[1]],
  dr_ensemblt,
  dr_ensembl,
  dr_gene_symbol,
  chr,
  strand,
  epd_tss = ifelse(strand == "+", tss_position+1,tss_position), # in EPD the + positions are shifted
  idx = tss_position
)]

# Separate the rows where multiple ensembl gene IDs are assigned
EPD_clean <- EPD_clean[, .(dr_ensembl = trimws(unlist(strsplit(dr_ensembl, " \\| "))),
                           dr_gene_symbol = trimws(unlist(strsplit(dr_gene_symbol, " \\| ")))), by = setdiff(names(EPD_clean), c("dr_ensembl", "dr_gene_symbol"))]
# remove duplicates
EPD_clean <- unique(EPD_clean)

# Get Ensembl data
ensdb <- get_ensdb("mm10")
ensdb_genes <- as.data.table(genes(ensdb))
ensdb_genes[, chr := paste0("chr", seqnames)]
ensdb_genes[, `:=`(
  entrezid = NULL,
  seqnames = NULL) ]

biotype_levels <- c(
  # Protein-coding (highest significance)
  "protein_coding",
  
  # Immunoglobulin genes
  "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "IG_LV_gene",
  
  # T-cell receptor genes
  "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
  
  # Long non-coding RNAs
  "lincRNA", "macro_lncRNA", "bidirectional_promoter_lncRNA", "lncRNA",
  
  # Other functional ncRNAs
  "antisense", "sense_intronic", "sense_overlapping", "3prime_overlapping_ncRNA",
  "processed_transcript", "TEC",
  
  # Small/structural RNAs
  "miRNA", "snRNA", "snoRNA", "scaRNA", "misc_RNA", "sRNA", "scRNA",
  "rRNA", "Mt_tRNA", "ribozyme",
  
  # Pseudogenes (lowest biological significance / non-functional)
  "processed_pseudogene", "unprocessed_pseudogene",
  "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene",
  "transcribed_unitary_pseudogene", "unitary_pseudogene",
  "polymorphic_pseudogene", "pseudogene",
  "IG_C_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "IG_V_pseudogene",
  "TR_J_pseudogene", "TR_V_pseudogene",
  "translated_unprocessed_pseudogene"
)

ensdb_genes[, gene_biotype := factor(gene_biotype, levels = biotype_levels, ordered = TRUE)]
setorder(ensdb_genes, gene_biotype, na.last = TRUE)

ensdb_tx <- as.data.table(transcripts(ensdb))
ensdb_tx[, chr := paste0("chr", seqnames)]
ensdb_tx[, seqnames := NULL]
ensdb_tx[, ensdb_TSS := ifelse(strand == "+", start, end)]
ensdb_tx[, ensdb_TES := ifelse(strand == "+", end, start)]
setnames(ensdb_tx, c("start", "end"), c("tx_start", "tx_end"))
ensdb_tx <- ensdb_genes[,.(gene_id, gene_name)][ensdb_tx, on = "gene_id"] # Add gene_name

transcript_biotype_levels <- c(
  # Protein-coding (highest significance)
  "protein_coding",
  
  # Immunoglobulin genes
  "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "IG_LV_gene",
  
  # T-cell receptor genes
  "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
  
  # Long non-coding RNAs
  "lincRNA", "macro_lncRNA", "bidirectional_promoter_lncRNA", "lncRNA",
  
  # Other functional ncRNAs / transcript categories
  "antisense", "sense_intronic", "sense_overlapping", "3prime_overlapping_ncRNA",
  "processed_transcript", "retained_intron", "TEC",
  
  # Transcripts flagged for aberrant/degraded processing
  "nonsense_mediated_decay", "non_stop_decay",
  
  # Small/structural RNAs
  "miRNA", "snRNA", "snoRNA", "scaRNA", "misc_RNA", "sRNA", "scRNA",
  "rRNA", "Mt_rRNA", "Mt_tRNA", "ribozyme",
  
  # Pseudogenes (lowest biological significance / non-functional)
  "processed_pseudogene", "unprocessed_pseudogene",
  "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene",
  "transcribed_unitary_pseudogene", "unitary_pseudogene",
  "translated_processed_pseudogene", "translated_unprocessed_pseudogene",
  "polymorphic_pseudogene", "pseudogene",
  "IG_C_pseudogene", "IG_D_pseudogene", "IG_pseudogene", "IG_V_pseudogene",
  "TR_J_pseudogene", "TR_V_pseudogene"
)
ensdb_tx[, tx_biotype := factor(tx_biotype, levels = transcript_biotype_levels, ordered = TRUE)]
setorder(ensdb_tx, tx_biotype, na.last = TRUE)

# Set gene IDs to NA if not part of the Ensembl gene database
EPD_clean[!(dr_ensembl %in% ensdb_genes$gene_id), dr_ensembl := NA]

# Set tx IDs to NA if not part of the Ensembl tx database
EPD_clean[!(dr_ensemblt %in% ensdb_tx$tx_id), dr_ensemblt := NA]

# Fill in gene names if ids present (if strand and chromosome matches)
# Note: gene IDs are unique ensdb_genes[, .SD[.N>1], gene_id]

# Fill gene_id from gene_symbol (chr/strand match), take the gene_id with the highest priority if multiple matches
EPD_clean[ensdb_genes[, .SD[1], .(gene_name)], dr_ensembl := fcoalesce(dr_ensembl, i.gene_id),
          on = .(chr, strand, dr_gene_symbol = gene_name)]

# Fill gene_id from transcript_id (chr/strand match), tx_ids are unique 
EPD_clean[ensdb_tx, dr_ensembl := fcoalesce(dr_ensembl, i.gene_id),
          on = .(chr, strand, dr_ensemblt = tx_id)]

# Fill gene_symbol from gene_id (chr/strand match), gene_ids are unique
EPD_clean[ensdb_genes, dr_gene_symbol := fcoalesce(dr_gene_symbol, i.gene_name),
          on = .(chr, strand, dr_ensembl = gene_id)]

# Add ensembl gene information 
EPD_clean <- ensdb_genes[,.(gene_chr = chr, gene_strand = strand, gene_start = start, gene_end = end, gene_id, gene_name, gene_biotype)][EPD_clean, on = .(gene_id = dr_ensembl)]


#---------------------------------------------#
#---- Evaluate and correct annotations ----
#---------------------------------------------#

# I consider an annotation correct if the gene IDs, chromosome and strand match and the TSS is upstream of the gene end
EPD_correct <- unique(EPD_clean[!is.na(gene_name) & chr == gene_chr & strand == gene_strand & ifelse(strand == "+", epd_tss < gene_end, epd_tss > gene_start)])
out_cols <- c("entry_name", "promoter_id", "gene_id","dr_ensembl", "gene_name",
              "dr_gene_symbol", "dr_ensemblt", "chr", "strand", "epd_tss")

# Select Promoters which need to be fixed: 126
EPD_toFix <- unique(EPD_clean[!(entry_name %in% EPD_correct$entry_name)])[,.(entry_name, promoter_id, strand, chr= chr, epd_tss,dr_ensembl = gene_id, dr_gene_symbol, dr_ensemblt, idx = epd_tss)]
length(unique(EPD_toFix$entry_name))

# Match based on transcript IDs: 3
EPD_FixedBYtx <- ensdb_tx[EPD_toFix[!is.na(dr_ensemblt)], on = .(tx_id = dr_ensemblt, strand, chr)][, dr_ensemblt := tx_id][!is.na(gene_id) & ifelse(strand == "+", epd_tss < tx_end, epd_tss > tx_start), ..out_cols]
fixed_entries <- EPD_FixedBYtx$entry_name
length(fixed_entries)

# Match based on gene name (gene_ids are wrong?): 7
EPD_FixedBYname <- ensdb_genes[EPD_toFix[!is.na(dr_gene_symbol)], on = .(gene_name = dr_gene_symbol, strand, chr)][, dr_gene_symbol := gene_name][!is.na(gene_id) & ifelse(strand == "+", epd_tss < end, epd_tss > start), ..out_cols]
fixed_entries <- c(fixed_entries, EPD_FixedBYname$entry_name)
unique(EPD_FixedBYname$entry_name) %>% length()

# If still doesn't match assign the closest annotated TSS within 500 bp: 30
setkey(EPD_toFix, chr, strand, idx)
ensdb_tx[, idx := ensdb_TSS]
setkey(ensdb_tx, chr, strand, idx)
EPD_FixedBYprox <- ensdb_tx[EPD_toFix[!(entry_name %in% fixed_entries)], roll = "nearest"][ifelse(strand == "+", epd_tss < tx_end, epd_tss > tx_start)][order(abs(epd_tss- ensdb_TSS)), .SD[1], entry_name][abs(epd_tss- ensdb_TSS) <= 500, ..out_cols]
fixed_entries <- c(fixed_entries, EPD_FixedBYprox$entry_name)
length(EPD_FixedBYprox$entry_name)

# Remove the TSS which gene annotations are missing in the ensdb: 0 promoters
EPD_missing_geneIDs <- EPD_toFix[!(entry_name %in% fixed_entries) & !is.na(dr_ensembl) & !(dr_ensembl %in% ensdb_genes$gene_id), entry_name]
fixed_entries <- c(fixed_entries, EPD_missing_geneIDs)

# Check how many entries were not resolved: result 93
sum(!(EPD_toFix$entry_name %in% fixed_entries))
EPD_toFix[!(entry_name %in% fixed_entries)]

# Assign to a transcript with the closest downstream TES: 40
EPD_FixedBYclosest <- EPD_toFix %>%
  dplyr::filter(!(entry_name %in% fixed_entries)) %>%
  dplyr::left_join(
    ensdb_tx,
    by = c("chr", "strand"),
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(if_else(strand == "+", epd_tss < ensdb_TES, epd_tss > ensdb_TES)) %>%
  dplyr::mutate(dist = abs(epd_tss - ensdb_TSS)) %>% # distance between the start of the transcript and EPD TSS
  dplyr::slice_min(dist, by = entry_name, n = 1, with_ties = FALSE) %>%
  dplyr::filter(dist <= 10000) %>%
  dplyr::select(all_of(out_cols))

fixed_entries <- c(fixed_entries, EPD_FixedBYclosest$entry_name)
# join together
EPD_fixed <- rbindlist(list(
  EPD_correct[ , `:=`(dr_ensembl = gene_id, chr = chr)][, ..out_cols],
  EPD_FixedBYname,
  EPD_FixedBYprox,
  EPD_FixedBYtx,
  EPD_FixedBYclosest
), fill = T)

# How many were not assigned? 49
setdiff(EPD_toFix$entry_name, EPD_fixed$entry_name) %>% length()
EPD_toFix[entry_name %in% setdiff(EPD_toFix$entry_name, EPD_fixed$entry_name)]

#--- Fill promoters for missing gene names from Ensdb ----
Ensdb_specific_prom <- ensdb_tx[!(gene_name %in% EPD_fixed$gene_name) & !is.na(tx_support_level) & tx_biotype %in% c("protein_coding", "lincRNA") & chr %in% paste0("chr", c(1:19, "X", "Y"))][
    order(tx_start)][, .(entry_name = paste(gene_name, 1:.N, sep = "_"),
                       promoter_id = 1:.N,
                       gene_id,
                       dr_ensembl =NA,
                       dr_ensemblt = tx_id,
                       gene_name,
                       dr_gene_symbol = NA,
                       chr, strand, epd_tss = ensdb_TSS), by = gene_name]

# Bind promoter and select the TSS with the highest signal
# Bin TSSs and compute signal per nt
Ensdb_specific_prom_bin <- Ensdb_specific_prom[, {
  start <- seq(epd_tss-150, epd_tss+150)
  if(strand=="-")
    start <- rev(start)
  .(start= start)
}, .(seqnames =chr, gene_id, strand, entry_name)]
Ensdb_specific_prom_bin[, end:= start]

# Compute CAGE signal
Ensdb_specific_prom_bin[strand=="+", CAGE_max:= vlite::bwCoverage(.SD,"db/public/ctssMaxCounts_FANTOM5_CAGE_TSS_mm10.fwd.bw")]
Ensdb_specific_prom_bin[strand=="-", CAGE_max:= vlite::bwCoverage(.SD, "db/public/ctssMaxCounts_FANTOM5_CAGE_TSS_mm10.rev.bw")]

# take the best match
Ensdb_specific_prom <- Ensdb_specific_prom_bin[order(-CAGE_max)][CAGE_max > 0, .SD[1], entry_name][, .(entry_name, corrected_tss = start)][Ensdb_specific_prom, on = "entry_name"]
Ensdb_specific_prom[is.na(corrected_tss), corrected_tss := epd_tss] # if all are 0 then return the original

# Combine EPD and Ensdb promoters
prom <- rbind(EPD_fixed,Ensdb_specific_prom[, .(entry_name,promoter_id, gene_id, gene_name, chr, strand, epd_tss = corrected_tss)], fill = T)

# Remove duplicates
prom <- unique(prom, by = c("gene_id", "gene_name", "chr", "strand", "epd_tss"))

# Check that all TSSs are upstream of TESs
ensdb_genes[prom, on = "gene_id"][ifelse(strand == "+", epd_tss >= end, epd_tss <= start)]

#--- Add CAGE signal ----
up <- 150
down <- 150

prom[, start := ifelse(strand == "+", epd_tss - up, epd_tss - down )]
prom[, end := ifelse(strand == "+", epd_tss + down, epd_tss + up)]
setnames(prom, "chr", "seqnames")
prom[strand == "+", CAGE_max := vlite::bwCoverage(.SD, "db/public/ctssMaxCounts_FANTOM5_CAGE_TSS_mm10.fwd.bw")]
prom[strand == "-", CAGE_max := vlite::bwCoverage(.SD, "db/public/ctssMaxCounts_FANTOM5_CAGE_TSS_mm10.rev.bw")]

# Add biotype
prom[ensdb_genes, gene_biotype := i.gene_biotype, on = .(gene_id, gene_name)]


# Add motifs
prom <- motifs[prom, on = "entry_name"]
setnames(prom, "epd_tss", "tss" )
# Export 
out_file = "db/annotations/GRCm38.v102.EPDnew_promoters.txt.gz"
select_cols <- c(
  "entry_name","gene_id", "gene_name", "gene_biotype",
  "seqnames", "strand", "tss", "CAGE_max",
  names(motifs)[-1]
)
fwrite(x =  prom[, ..select_cols],
       file = out_file, sep = "\t", compress = "gzip")
