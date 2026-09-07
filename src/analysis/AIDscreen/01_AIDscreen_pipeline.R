#-------------------------------------------#
# Dependencies
#-------------------------------------------#
require(data.table)
require(parallel)
devtools::load_all("/groups/stark/vloubiere/vlite")
library(here)

#-------------------------------------------#
# Working directory
#-------------------------------------------#
# Set working directory
wd <- "/groups/stark/vojtek/Bioinfo_project/PPRC1_paper/Pprc1_git/"
setwd(wd)

experiment_prefix <- "Pprc1_AIDscreen"
#-------------------------------------------#
# select directory for intermediate files
#-------------------------------------------#
# tmp_dir = "db"
tmp_dir = paste("/scratch-cbe/users/matus.vojtek", gsub(".*/(.*?)$", "\\1", wd), experiment_prefix, sep = "/")

# Create tmp dir if it doesn't exists
if(!dir.exists(tmp_dir)){
  dir.create(tmp_dir, recursive = T)
}

#-------------------------------------------#
# Create directories for tmp files
#-------------------------------------------#
dirs_ls <- list(
  "db/fastq",
  "db/bam",
  "db/bam_unique",
  "db/bed",
  "db/annotated_integrations",
  "logs"
)

# Create directories if not existing
lapply(dirs_ls, function(x){

  # Path to the directory
  .tmp_path = paste(tmp_dir, x, experiment_prefix, sep = "/")

  # If doesn't exists create
  if(!dir.exists(.tmp_path)){
    dir.create(.tmp_path, recursive = T)
  }

})

#-------------------------------------------#
# Check directories in the wd
#-------------------------------------------#
dirs_keep_ls <- list(
  "Rdata",
  "db/bam",
  "db/bam_unique",
  "db/bed",
  "db/annotated_integrations"
)

# Create directories if not existing
lapply(dirs_keep_ls, function(x){

  # Path to the directory
  .tmp_path = paste(wd, x, experiment_prefix, sep = "/")

  # If doesn't exists create
  if(!dir.exists(.tmp_path)){
    dir.create(.tmp_path, recursive = T)
  }

})

#-------------------------------------------#
# Load metadata and create unique IDs
#-------------------------------------------#
# Read in metadata
meta <- fread("metadata_AIDscreen.txt", header = T, na.strings = c("",NA))


# Replace "_" with "-"
meta <- meta[, c(lapply(.SD, function(x) gsub("_", "-", x)), .(bam_path = bam_path)), .SDcols = setdiff(names(meta), "bam_path")]

# convert to numeric
meta$frame <- as.numeric(meta$frame)

# Create IDs
InfoCols <- c("expName", "reporter", "sample", "frame", "auxin", "fraction", "enzyme", "replicate")

InfoToUse <- sapply(meta[, ..InfoCols],
                    function(x) length(unique(x)) > 1)

meta$sampleID <- do.call(paste, c(meta[, .SD, .SDcols = InfoCols[InfoToUse]], sep="_"))

# Check if sampleIDs are unique
if(any(duplicated(meta$sampleID))){
  warning("some sampleIDs are duplicated\n\nCreating unique IDs...")
  meta$sampleID <- make.unique(meta$sampleID)
}


#-------------------------------------------#
# Output files
#-------------------------------------------#
meta[, fq1:= paste0(tmp_dir,"/db/fastq/",experiment_prefix, "/", sampleID, "_1.fq.gz")]
meta[(layout=="PAIRED"), fq2:= paste0(tmp_dir, "/db/fastq/",experiment_prefix, "/", sampleID, "_2.fq.gz")]
meta[, fq1_trim:= gsub(".fq.gz$", "_trimmed.fq.gz", fq1)]
meta[, bam:= paste0(tmp_dir, "/db/bam/",experiment_prefix, "/",  sampleID, ".bam"),]
meta[, bam_stats:= gsub(".bam$", "_stats.txt", bam)]
meta[, bam_unique:= paste0(tmp_dir, "/db/bam_unique/",experiment_prefix, "/", sampleID, "_q30_unique.bam"),]
meta[, bed_file:= paste0(tmp_dir, "/db/bed/",experiment_prefix, "/", sampleID, ".bed"),]
meta[, bed_file_SlippageFilt:= paste0(tmp_dir, "/db/bed/",experiment_prefix, "/", sampleID, "_filt.bed"),]
meta[, bed_anno:= paste0(tmp_dir, "/db/annotated_integrations/",experiment_prefix, "/", sampleID, "_filt_anno.bed"),]
meta[, bed_compatible:= paste0(tmp_dir, "/db/bed/",experiment_prefix, "/", sampleID, "_compatible.bed"),]
meta[, bed_anno_stats:= paste0(tmp_dir, "/db/annotated_integrations/",experiment_prefix, "/", sampleID, "_filt_anno_stats.txt"),]

#-------------------------------------------#
# Job parameters and modules
#-------------------------------------------#
cores <- 12
mem <- 32
conda_env <- "/groups/stark/conda_envs/.conda/envs/vlite_pipelines/"

#-------------------------------------------#
# Demultiplex VBC bam file
#-------------------------------------------#
meta[, demultiplex_cmd:= {
  if(!file.exists(fq1)) # fq2 is not used anyway
    if(T)
    {
      BC <- paste0("'^BC:Z:", gsub("\\|", "\\|^BC:Z:", barcodes), "'")
      fq_prefix <- normalizePath(gsub("_1.fq.gz$", "", fq1), mustWork = F)
      .f <- fcase(layout=="SINGLE" & sequencer=="NextSeq", "src/functions/demultiplex_se_12.pl", # se reads, BC is in column 12
                  layout=="PAIRED" & sequencer=="NextSeq550", "src/functions/demultiplex_pe_12.pl", # pe reads, BC is in column 12 (typically what we use)

                  layout=="PAIRED" & sequencer=="MiSeq", "src/functions/demultiplex_pe_12.pl", # pe reads, BC is in column 12 (typically what we use)

                  layout=="PAIRED" & sequencer=="NextSeq2000", "src/functions/demultiplex_pe_12.pl") # pe reads, BC is in column 12
      cmd <- paste("samtools view -@", cores-1, bam_path,
                   "| perl",
                   normalizePath(.f),
                   BC,
                   fq_prefix,
                   "; gzip",
                   paste0(fq_prefix, "_1.fq"))
      if(!is.na(fq2))
        cmd <- paste0(cmd, "; gzip ", fq_prefix, "_2.fq")
      cmd
    }
}, .(layout, barcodes, fq1, fq2)]


## Trim reads
meta[, trim_cmd:= {
  if(!file.exists(fq1_trim))
  {
    cmd <- paste0("trim_galore --gzip -o ", dirname(fq1_trim), "/ ", fq1)
    cmd
  }
}, fq1_trim]

#-------------------------------------------#
# Alignment
#-------------------------------------------#
meta[, aln_cmd:= {
  if(!file.exists(bam))
  {
    #### BOWTIE 2
    cmd <- paste("bowtie2 -p", cores,
                 "-U", paste0(fq1_trim, collapse= ","),
                 "-x", normalizePath("/groups/stark/vojtek/Resources/03_index/Bowtie2Index/GRCm38/GRCm38", mustWork = F),
                 "| samtools sort -@", cores-1, "-o", bam)
  }
}, bam]

#-------------------------------------------#
# Get bam stats
#-------------------------------------------#
meta[, stats_cmd:= {
  if(!file.exists(bam_stats))
    paste("samtools stats -@", cores-1, bam, "| grep ^SN | cut -f 2- >", bam_stats)
}, .(bam, bam_stats)]

#-------------------------------------------#
# Collapsed bam files
#-------------------------------------------#
meta[, collapse_cmd:= {
  if(!file.exists(bam_unique))
    paste("samtools sort -n -@", cores-1,  bam,
          "| samtools fixmate -m - - | samtools sort -@", cores-1,
          "| samtools markdup --duplicate-count -r - - | samtools view -q 30 -b -o",  bam_unique)
}, .(bam, bam_unique)]

#-------------------------------------------#
# Compute insertions
#-------------------------------------------#
meta[, insertions_cmd:= {
  if(any(!file.exists(c(bed_file, bed_anno, bed_file_SlippageFilt))))
  {
    paste("Rscript",
          normalizePath("src/scripts/bamToBed_and_assign_insertions.R"),
          bam_unique,
          normalizePath("/groups/stark/vojtek/Resources/01_annotations/mm10_GRCm38_v102/Mus_musculus.GRCm38.102.gtf.gz"),
          normalizePath("/groups/stark/vojtek/Bioinfo_project/PPRC1_paper/Pprc1_git/db/annotations/GRCm38.102_GeneTrapCompatibleCodingIntrons.txt.gz"),
          frame,
          bed_file,
          bed_file_SlippageFilt,
          bed_anno,
          bed_compatible,
          bed_anno_stats
    )
  }
}, .(bam_unique, bed_file, bed_anno)]


#-------------------------------------------#
# Copy from tmp_dir
#-------------------------------------------#
cols_to_copy <- c("bam_unique", "bam_stats", "bed_file", "bed_file_SlippageFilt",
                  "bed_anno", "bed_compatible", "bed_anno_stats")

# Create new columns with cp_ prefix
meta[, paste0("cp_", cols_to_copy, "_cmd") :=  lapply(.SD, function(x) paste("cp",
                                                                             x,
                                                                             gsub(tmp_dir, wd, x),
                                                                             sep = " ")),
     .SDcols = cols_to_copy]

#-------------------------------------------#
# CHECK AND SUBMIT
#-------------------------------------------#
meta[, cmd:= paste(
  c(
    paste0("source activate ", conda_env), # Activate conda environment
    na.omit(unique(unlist(.SD)))), # execute commands
  collapse= "; "),
  by = sampleID, .SDcols= patterns("_cmd$")]

meta[, {
  vlite::bsub(cmd,
              cores= cores,
              mem = mem,
              name = "AIDscreen",
              time = '06:00:00',
              wdir = wd,
              logs= paste0(tmp_dir, "/logs/"),
              execute = T)
}, cmd]

# -------------------------------------------#
# Save
# -------------------------------------------#
saveRDS(meta[, grep("cmd$", invert = T, names(meta)), with= F],
        paste0(wd, "/Rdata/AIDscreen_processed_metadata.rds"))

