#----------------------------------------------#
# Dependencies
#----------------------------------------------#
library(data.table)
library(rtracklayer)

#----------------------------------------------#
# Parse inputs
#----------------------------------------------#
# Import blaclisted regions UCSC
blacklist <- fread("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz")
setnames(blacklist, c("seqnames", "start", "end", "name"))
blacklist[, start:= start+1]
blacklist <- rbind(
  data.table::copy(blacklist)[, strand:= "+"],
  data.table::copy(blacklist)[, strand:= "-"]
)

# Import tRNAs gtf file
# tRNAs from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.tRNAs.gtf.gz
tRNAs <- rtracklayer::import("/groups/stark/vloubiere/genomes/Mus_musculus/GENCODE/gencode.vM25.tRNAs.gtf.gz")
tRNAs <- as.data.table(tRNAs)

# Import genes gtf file
others <- rtracklayer::import("/groups/stark/vojtek/Resources/01_annotations/mm10_GRCm38_v102/Mus_musculus.GRCm38.102.gtf.gz")
others <- as.data.table(others)
other_mask_types <- c(
  "rRNA","Mt_rRNA",
  "Mt_tRNA",
  "snRNA","snoRNA","scaRNA",
  "miRNA",
  "scRNA","sRNA",
  "misc_RNA",
  "ribozyme"
)
others <- others[gene_biotype %in% other_mask_types & type == "gene"]

# Combine
mask <- rbind(others[, .(seqnames, start, end, strand, cluster.id= gene_biotype)],
              tRNAs[, .(seqnames, start, end, strand, cluster.id= gene_type)],
              blacklist[, .(seqnames, start, end, strand, cluster.id= name)])
mask[, start:= as.integer(start)]
mask[, cluster.id:= paste0(cluster.id, "__", seqnames, ":", start, "-", end, ":", strand)]
mask[, width:= as.integer(end-start+1)]

# Order
setcolorder(mask, "width", after= "end")
setorderv(mask, c("seqnames", "start", "end"))

# SAVE
#exportBed(mask, "db/annotations/mm10_blacklisted_regions.bed")
saveRDS(mask, "db/annotations/mm10_PROseq_blacklisted_regions.rds")