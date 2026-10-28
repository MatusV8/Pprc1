
# Load packages
devtools::load_all("/groups/stark/vloubiere/vlite/")

#==============================================================================
# IMPORT METADATA
#==============================================================================

# Load CUT&RUN sample metadata (skip first 4 rows)
meta <- fread("metadata_ChIP.txt")
meta[, sampleID := paste(method, experiment, cell_line, antibody, treatment, replicate, sep = "_")]
meta[, condition := paste(method, experiment, cell_line, antibody, treatment, sep = "_")]
meta[, input := as.character(input)]
meta[method == "ChIPseq", input := gsub("V5", "noAb", sampleID)]
meta[method == "CutNRun", input := fifelse(grepl("rep1", sampleID), "CutNRun_exp2_Parental_V5_noTreatment_rep1", "CutNRun_exp2_Parental_V5_noTreatment_rep2")]

#==============================================================================
# DEFINE I7 COLUMN BASED ON SEQUENCER
#==============================================================================

meta[, i7.column := fcase(sequencer %in% c("NextSeq", "Aviti"), 14)]

#==============================================================================
# DEMULTIPLEXING (ONLY IF FASTQ FILES ARE MISSING)
#==============================================================================

meta <- meta[, {
  # Create demultiplexing command(s)
  cmd <- cmd_demultiplexVBCfile(
    vbcFile = bam_path,
    layout = layout,
    i7 = barcode,
    i7.column = i7.column,
    i5 = i5, 
    fq.output.folder = "db/fq/ChIP/"
  )
  # Submit demultiplexing jobs
  vl_submit(
    cmd,
    logs = "db/logs/ChIP/demultiplexing/",
    execute =T,
    overwrite = F,
    time = "08:00:00"
  )
  
  # Collect output file names
  output.files <- cmd[, .(path = paste0(path, collapse = ",")), file.type]
  output.files <- data.table::transpose(output.files, make.names = 1)
  if (!"fq2" %in% names(output.files))
    output.files[, fq2 := as.character(NA)]
  
  # Return metadata and output file names
  c(.SD, output.files)
}, names(meta)]

#==============================================================================
# SAVE UPDATED METADATA
#==============================================================================

saveRDS(meta, "Rdata/ChIP_metadata_demultiplexed.rds")
