# Retrieve genomic coordinates ----
file <- "db/public/ctssMaxCounts_FANTOM5_CAGE_TSS_mm10.fwd.bw"
if(!file.exists(file))
  download.file(url = "https://hgdownload.soe.ucsc.edu/gbdb/mm10/fantom5/ctssMaxCounts.fwd.bw",
                destfile = file)

file <- "db/public/ctssMaxCounts_FANTOM5_CAGE_TSS_mm10.rev.bw"
if(!file.exists(file))
  download.file(url = "https://hgdownload.soe.ucsc.edu/gbdb/mm10/fantom5/ctssMaxCounts.rev.bw",
                destfile = file)
