library(data.table)
library(tidyverse)


#------ Degron-based screens -----
meta_screen <- readRDS("Rdata/AIDscreen_processed_metadata.rds")
meta_screen[, sampleID := gsub("_([123])_", "_frame\\1_", sampleID)]

meta_screen %>%
  transmute(
    sampleID,
    fraction,
    frame = paste0("frame", frame),
    title = sampleID,
    lib_strategy = "Degron-based screen",
    organism = "Mus musculus",
    cell_line = case_when(
      reporter == "EF1a" ~ "cMV13c2",
      reporter == "Atp5k" ~ "cMV18c5",
      reporter == "test" ~ "AN3-12"),
    cell_type = "Mouse embryonic stem cells (haploid)",
    molecule = "genomic DNA",
    layout = "paired-end",
    sequencer,
    additional_file = "Degron-based-screen_counts.csv",
    fq1 = basename(fq1),
    fq2 = basename(fq2)) %>%
  group_by(sampleID) %>%
  summarise(across(everything(), ~ list(unique(na.omit(.x)))), .groups = "drop") %>%
  tidyr::unnest_wider(-sampleID, names_sep = "_") %>%
  write_tsv(., file = "Screen_geo.txt")

# Count table
AIDscreen_counts <- fread("db/count_tables/AID_screen_counts.csv")
colnames(AIDscreen_counts) <- gsub("(EF1a|Atp5k|test)_([1-3])_(GFPlow|Input|GFPhigh)_(rep[12])", "\\1_Intein_frame\\2_\\3" , colnames(AIDscreen_counts))
fwrite(AIDscreen_counts, "../geo/DegronScreen/AIDscreen_counts.csv.gz", compress = "gzip")

# FC table
paste0(
  "cp ",
  "db/FC_tables/AIDscreen_CMH_enrichment.csv",
  " ../geo/DegronScreen/"
)

# Integrations
int <- fread("Rdata/AIDscreen_Integrations.txt.gz")
int[, sampleID := gsub("(EF1a|Atp5k|test)_(Intein|Classic)_([1-3])_(GFPlow|Input|GFPhigh)_(rep[12])", "\\1_\\2_frame\\3_\\4" , sampleID)]
fwrite(int, "../geo/DegronScreen/AIDscreen_Integrations.txt.gz", sep = "\t", compress = "gzip")

#------ ChIP-seq -----
meta_ChIP_<- readRDS("Rdata/ChIP_metadata_demultiplexed.rds")

ChIP_processed <- readRDS("Rdata/ChIP_metadata_processed.rds")
ChIP_bw <- readRDS("Rdata/ChIP_metadata_peakCalling.rds")

# Copy fq files
ChIP_fqs <- ChIP_processed[method == "ChIPseq", paste(c(fq1, fq2), collapse = " ")]

paste("cp ",
      ChIP_fqs,
      " ../geo/ChIPseq/")


# Copy bw files
ChIP_bws <- ChIP_bw[ grepl("ChIPseq", condition) & !grepl("noAb", condition), unique(bw)]
paste0(
  "cp ",
  paste(ChIP_bws, collapse = " "),
  " ../geo/ChIPseq/"
)

# Copy bed files
ChIP_narrowPeak <- ChIP_bw[ grepl("ChIPseq", condition) & !grepl("noAb", condition), unique(peaks)]
paste0(
  "cp ",
  paste(ChIP_narrowPeak, collapse = " "),
  " ../geo/ChIPseq/"
)

#------ Cut & Run -----
# Copy fq files
CutNRun_fqs <- ChIP_processed[method == "CutNRun", paste(c(fq1, fq2), collapse = " ")]

paste("cp ",
      CutNRun_fqs,
      " ../geo/CutNRun/")


# Copy bw files
CutNRun_bws <- ChIP_bw[ grepl("CutNRun", condition) & !grepl("Parental", condition), unique(bw)]
paste0(
  "cp ",
  paste(CutNRun_bws, collapse = " "),
  " ../geo/CutNRun/"
)

# Copy bed files
CutNRun_narrowPeak <- ChIP_bw[ grepl("CutNRun", condition) & !grepl("Parental", condition), unique(peaks)]
paste0(
  "cp ",
  paste(CutNRun_narrowPeak, collapse = " "),
  " ../geo/CutNRun/"
)


#------ PROseq -----
PROseq_meta <- readRDS("Rdata/PROseq_demultiplexed_fastqs.rds")
PROseq_processed <- readRDS("Rdata/PROseq_processed_files.rds")

# FQs
PROseq_fqs <- PROseq_meta[, unique(c(fq1, fq2))]
paste("cp ",
      paste(PROseq_fqs, collapse = " "),
      "../geo/PROseq/")

#BWs
PROseq_bws<- PROseq_processed[, unique(c(ps.bw, ns.bw))]
paste("cp ",
      paste(PROseq_bws, collapse = " "),
      "../geo/PROseq/")

#DESeq FC table
file.copy("db/FC_tables/PROseq/Pprc1.depl_transcript_spikeIn_norm_DESeq2_FC.txt", to = "../geo/PROseq/")

