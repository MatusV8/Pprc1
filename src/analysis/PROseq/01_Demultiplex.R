#=============================================================#
# Requirements ----
#=============================================================#
setwd("/groups/stark/vojtek/Bioinfo_project/PPRC1_paper/Pprc1_git/")
require(data.table)
devtools::load_all("/groups/stark/vloubiere/vlite/")

#=============================================================#
# Parameters ----
#=============================================================#

# Load metadata
meta <- fread("metadata_PROseq.txt")

# Loop over the sampleID values (see the 'by=' argument on the last line)
meta.demultiplex <- meta[, {
  cmd <- cmd_demultiplexVBCfile(
    vbcFile = bam_path,
    layout = layout,
    i7 = barcode, # i7 barcode
    i7.column = 14, # This specifies in which column of the bam file the i7 sequences are stored. By default it is set to 14 but should sometimes be set to 12, depending on the sequencer. If you are not sure in which column are your indexes, check the head of your bam file using samtools view.
    output.prefix= NULL, # The prefix that will be appended to the output file's name
    fq.output.folder= "db/fq/PROseq", # The output folder where the output fastq files will be saved
    proseq.eBC = eBC, # eBC sequence (only used for PRO-Seq)
    proseq.umi.length = 10, # The length of the PRO-Seq eBC sequence (will not be used as proseq.eBC is set to NULL)
    cores = 8 # The number of processors to use
  )

  # Submit to the cluster
  vl_submit(
    cmd,
    execute = TRUE, # Should the command be execute? Default= TRUE.
    overwrite = FALSE, # If set to FALSE (default), the command will only be submitted if output files don't exist. If set to TRUE, they are overwritten
    mem = 16, # The memory for the job (in Gb)
    logs = "db/logs/PROseq/demultiplex/", # Chose the folder where log files will be saved
    job.name= sampleID
  )
  # Return the output files paths
  transpose(cmd[, .(file.type, path)], make.names = T)
}, by= names(meta)] # For each line in the file

# Create a folder to save RDS files
if(!dir.exists("Rdata")){
  dir.create("Rdata", recursive = T, showWarnings = F)
}


# Save the file containing demultiplexed fastq paths
meta.demultiplex.file <- "Rdata/PROseq_demultiplexed_fastqs.rds"
saveRDS(meta.demultiplex,
        meta.demultiplex.file)
