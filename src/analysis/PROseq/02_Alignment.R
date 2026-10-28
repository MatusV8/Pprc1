# Libraries
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Check the help function
#?vlite::proseqProcessing

# Import the demultiplexed fastq metadata generated earlier
meta.demultiplex <- readRDS("Rdata/PROseq_demultiplexed_fastqs.rds")

# Generate the commands and submit to the server
meta.processed <- meta.demultiplex[, {
  # Create the command
  cmd <- vlite::proseqProcessing(
    fq1= fq1, # Only the first read is used
    ref.genome = genome, # Genome to which the reads should first be aligned
    output.prefix = sampleID, # Output file name prefix
    spike.genome = spikein_genome, # The read that did not align to the reference genome will then be aligned to the spike in genome
    fq.output.folder = "db/fq/PROseq/", # Output folder where trimmed fq files will be saved
    bam.output.folder = "db/bam/PROseq/", # Output folder where aligned bam files will be saved
    alignment.stats.output.folder = "db/alignment_stats/PROseq/", # Output folder where alignment statistics will be saved  
    counts.output.folder = "db/count_tables/PROseq/", # Output folder where alignment statistics will be saved  
    counts.stats.output.folder = "db/stats/PROseq/", # Output folder where alignment statistics will be saved
    Rpath = "Rscript",
    bw.output.folder = "db/bw/PROseq/" # Output folder for bigwig files
  )

  # Submit to the cluster
  vl_submit(
    cmd,
    time = "02:00:00",
    overwrite = F, # In the case where output file exist, should the analyses still be ran again?
    execute = T, # Execute the command
    logs = "db/logs/processing/PROseq/", # Where to save the logs
  )
  
  # Return commands and output files
  transpose(cmd[, .(path= paste0(path, collapse= ",")), file.type], make.names = T)
}, by= .(user, batch, experiment, cell_line, treatment, replicate, sampleID,
         DESeq2_condition, DESeq2_control,
         genome, spikein_genome)] # For each sampleID

# Save the metadata containing processed files
saveRDS(meta.processed, "Rdata/PROseq_processed_files.rds")
