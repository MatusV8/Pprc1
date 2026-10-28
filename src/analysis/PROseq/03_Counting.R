devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import the demultiplexed fastq metadata generated earlier
meta.processed <- readRDS("Rdata/PROseq_processed_files.rds")

features <- list(
  transcript = "db/annotations/mm10.GRCm38.v102.EPDnew.ProtCodingLincRNAonly.CAGEcorrected.250up250down_transcript.rds",
  geneBody   = "db/annotations/mm10.GRCm38.v102.EPDnew.ProtCodingLincRNAonly.CAGEcorrected.250up250down_geneBody.rds",
  promoter   = "db/annotations/mm10.GRCm38.v102.EPDnew.ProtCodingLincRNAonly.CAGEcorrected.250up250down_promoter.rds"
)

# At this point, each sample should correspond to a unique entry in the metadata
stopifnot(uniqueN(meta.processed$sampleID)==nrow(meta.processed))

results_list <- lapply(names(features), function(feature) {
  
  # Validate annotation file exists
  if (!file.exists(features[[feature]])) {
    stop("Annotation file not found for feature: ", feature)
  }
  
  # Add feature type to metadata
  local_meta <- copy(meta.processed)[, feature_type := feature]
  
  # Generate and submit count commands for each sample
  local_meta[, {
    
    cmd <- vlite::cmd_countPROseqReads(
      umi.count.file = umi.counts.ref, # UMI-collapsed reads
      output.prefix = sampleID,
      annotation.file = features[[feature]], # Full transcript annot
      feature = feature, # The name of the features
      count.tables.output.folder = "db/count_tables/PROseq/", # Folder where count tables should be saved
      Rpath = "Rscript" # The R executable used to run the R script
    )
    
    # Submit to the clusters
    vl_submit(
      cmd,
      cores= 2,
      mem= 16,
      time= '00:15:00',
      overwrite = T, # In the case where output files exist, should they be overwritten?
      execute = T, # Excute the command
      logs = "db/logs/PROseq/counting", # Where to save the logs
    )
    out <- transpose(cmd[, .(file.type, path)], make.names = T)
    out[, feature_type := feature]
  }, by = names(meta.processed)]
  
})

# Combine results for all features
meta <- rbindlist(results_list)

#==============================================================================
# SAVE UPDATED METADATA WITH COUNT TABLE PATHS
#==============================================================================

saveRDS(meta, "Rdata/PROseq_count_files.rds")
