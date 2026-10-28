
# Load packages
devtools::load_all("/groups/stark/vloubiere/vlite/")

#==============================================================================
# LOAD PROCESSED METADATA
#==============================================================================

meta <- readRDS("Rdata/ChIP_metadata_processed.rds")

# Map input BAM files to each sample
meta$input.bam <- meta$bam[match(meta$input, meta$sampleID)]

#==============================================================================
# CALL PEAKS ON INDIVIDUAL REPLICATES
#==============================================================================

meta.peak <- meta[, {
  # Create peak calling command for each replicate
  cmd <- cmd_peakCalling(
    bam = unique(bam),
    bam.input = unique(input.bam),
    layout = layout,
    output.prefix = sampleID,
    keep.dup = 1,
    compute.model = FALSE,
    extsize = 300,
    shift = 0,
    broad = FALSE,
    genome.macs2 = "mm",
    genome = genome,
    peaks.output.folder = "db/bed/ChIP_peaks/",
    bw.output.folder = "db/bw/ChIP/",
    Rpath = "Rscript"
  )
  
  # Submit peak calling job
  vl_submit(
    time = '01:30:00',
    cmd,
    logs = "db/logs/ChIP/peakCalling/",
    overwrite = T,
    execute = T
  )
  
  # Collect output file names
  output <- cmd[, .(path = paste0(path, collapse = ",")), file.type]
  output <- data.table::transpose(output, make.names = 1)
  
  # Return metadata and output file names
  c(.SD, output)
}, .(sampleID, genome, layout)]

#==============================================================================
# CALL PEAKS ON MERGED REPLICATES (PER CONDITION)
#==============================================================================

meta.peak.merged <- meta[, {
  # Create peak calling command for merged replicates
  cmd <- cmd_peakCalling(
    bam = unique(bam),
    bam.input = unique(input.bam),
    layout = layout,
    output.prefix = paste0(condition, "_merged"),
    keep.dup = 1,
    compute.model = FALSE,
    extsize = 300,
    shift = 0,
    broad = FALSE,
    genome.macs2 = "mm",
    genome = genome,
    peaks.output.folder = "db/bed/ChIP_peaks/",
    bw.output.folder = "db/bw/ChIP/",
    Rpath = "Rscript"
  )
  
  # Submit peak calling job
  vl_submit(
    time = '01:30:00',
    cmd,
    logs = "db/logs/ChIP/peakCalling/",
    overwrite = T,
    execute = T
  )
  
  # Collect output file names and append ".merged" to column names
  output <- cmd[, .(path = paste0(path, collapse = ",")), file.type]
  output <- data.table::transpose(output, make.names = 1)
  #setnames(output, function(x) paste0(x, ".merged"))
  
  # Return metadata and output file names
  c(.SD, output)
}, .(condition, genome, layout)]

#==============================================================================
# SAVE UPDATED METADATA WITH PEAK CALLING OUTPUTS
#==============================================================================
value_cols <- c("condition", "replicate", "mapq.stats", "align.stats", "bam", "peaks", "bedgraph", "bw")


meta.peak.merged.sub <- meta.peak.merged[, ..value_cols]
meta.peak.merged.sub[, replicate := "merged"]

meta.merge <- rbind(unique(meta.peak[, ..value_cols]), meta.peak.merged.sub)

saveRDS(unique(meta.merge), "Rdata/ChIP_metadata_peakCalling.rds")

