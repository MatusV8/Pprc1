# Set working directory
setwd("/groups/stark/vojtek/Bioinfo_project/PPRC1_paper/Pprc1_git/")

# Load packages
devtools::load_all("/groups/stark/vloubiere/vlite/")

#==============================================================================
# LOAD DEMULTIPLEXED METADATA
#==============================================================================

meta <- readRDS("Rdata/ChIP_metadata_demultiplexed.rds")

#==============================================================================
# PROCESS FASTQ FILES FOR CUT&RUN
#==============================================================================

meta.processed <- meta[, {
  # Create processing command for each sample
  cmd <- vlite::cutnrunProcessing(
    fq1= fq1,
    fq2 = fq2,
    output.prefix = sampleID,
    genome = genome, # The genome version to use for alignment (here mm10)
    max.ins= 1000, # The maximum insert size for Bowtie2. Default= 1000.
    fq.output.folder = "db/fq/ChIP", # Output folder for trimmed fq files
    bam.output.folder = "db/bam/ChIP", # Output folder for aligned bam files
    alignment.stats.output.folder = "db/stats/ChIP", # Output folder for alignment statistics
    cores = 8 # Number of cores to use for the job
  )
  
  # Submit processing job
  vl_submit(
    cmd,
    logs = "db/logs/ChIP/processing/",
    overwrite = FALSE,
    execute = T
  )
  
  # Collect output file names
  output <- cmd[, .(path = paste0(path, collapse = ",")), file.type]
  output <- data.table::transpose(output, make.names = 1)
  
  # Add empty fq2.trim for single-end files
  if (!"fq2.trim" %in% names(output)) {
    output[, fq2.trim := NA_character_]
    setcolorder(output, c("fq1.trim", "fq2.trim"))
  }
  
  # Return metadata and output file names
  cbind(.SD, output)
}, .(sampleID, genome)]

#==============================================================================
# SAVE UPDATED METADATA
#==============================================================================

saveRDS(meta.processed, "Rdata/ChIP_metadata_processed.rds")
