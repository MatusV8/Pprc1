binBed <- function(bed, bins.width = 1) {
  .dt <- copy(bed)
  .dt[, line.idx := .I]
  
  res <- .dt[, {
    starts <- seq(start, end, by = bins.width)
    ends <- pmin(starts + bins.width - 1, end)
    .(start = starts, end = ends, seqnames = seqnames, strand = strand)
  }, by = setdiff(names(.dt), c("seqnames", "start", "end", "strand"))]
  
  res[, bin.idx := rowid(line.idx)]
  res[strand == "-", bin.idx := rev(bin.idx), line.idx]
  res[, .SD, .SDcols = c(names(bed), "bin.idx")]
}
