# Libs
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import the counts metadata
meta.counts <- readRDS("Rdata/PROseq_count_files.rds")

# Edit the DESeq2_condition and control columns 
meta.counts[,DESeq2_condition := treatment]
meta.counts[, DESeq2_control := gsub("[1-9]{1,2}", "0", treatment)]
meta.counts[, DESeq2_name := gsub("_rep[1-9]", "", sampleID)]

# Differential analysis
meta.FC <- lapply(c("spikeIn", "libSize"), function(x){
  meta.counts[, {
  # Generate the command
  cmd <- vlite::cmd_DESeq2_PROseq(
    umi.count.tables = count.table, # Tables containing the UMI counts per gene feature (one per line) and per sample (one per column)
    output.prefix = experiment, # The name of the output files
    sample.names = sampleID, # The names of the sample. 
    conditions = DESeq2_condition,  # The condition to which each sample belongs. Here it will be: control, IAA_3h, control, IAA_3h. This is how DESeq2 knows that each condition has two replicates.
    ctl.conditions = DESeq2_control, # The control condition (of course, only the comparisons where the DESeq2_condition and the DESeq2_control are different will be computed).
    ref.genome.stat.files = umi.stats.ref, # The file containing the statistics for the reference genome.
    spikein.stat.files = umi.stats.spike, # The file containing the statistics for the spike-in genome. Required for spike-in normalization! See next column.
    normalization = x, # We want to use spike-in normalization. Other possible values include 'default' and 'libsize'. See ?vlite::cmd_DESeq2_PROseq
    feature = feature_type, # The name of the feature (which will be added ot the output file names)
    dds.output.folder = "db/dds/PROseq/", # The folder where .dds DESeq2 object should be saved.
    FC.tables.output.folder = "db/FC_tables/PROseq/", # The folder where FC tables should be saved.
    MAplots.output.folder = "db/FC_tables/PROseq/", # The folder where MA plots (.pdf) will be saved.
    Rpath = "Rscript" # The R path executable that will be used to run the Rscript. DESeq2 should be installed there!
  )
  
  # Submit to the cluster
  vl_submit(
    cmd,
    cores= 2,
    mem= 8,
    time= '00:40:00',
    overwrite = T, # In the case where output files exist, should they be overwritten?
    execute = T, # Excute the command
    logs = "db/logs/PROseq_DESeq2/", # Where to save the logs
  )
  
  # Return output files paths
  out <- data.table::transpose(cmd[, .(file.type, path)], make.names = T)
  out[, normalisation := x]
  }, by = .(experiment, feature_type)] # If you want to analyse several experiments separately, use the "by" argument to split the samples accordingly and create separated objects.
})

saveRDS(rbindlist(meta.FC), file = "Rdata/PROseq_DESeq2_mm10.GRCm38.v102.EPDnew.CAGEcorrected.rds")