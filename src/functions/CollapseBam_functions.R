#' Read BAM files as data.tables
#'
#' @description
#' Reads one  BAM files using \pkg{Rsamtools}, extracts the standard
#' alignment fields (reference name, position, strand, width, sequence)
#' together with any requested optional BAM tags.
#'
#' @param bam_path Path to one or more BAM file as a character.
#' @param tags Character vector. Names of optional BAM tags to extract
#'   (e.g. \code{"NM"}, \code{"MD"}). Defaults to \code{character(0)}
#'   (no extra tags are extracted).
#'
#' @return A \code{data.table} with columns \code{seqnames}, \code{start},
#'   \code{strand}, \code{width}, \code{sequence}, plus one column for each
#'   tag requested in \code{tags}.
#'
#' @details
#' Internally this calls \code{\link[Rsamtools]{scanBam}} with a
#' \code{ScanBamParam} requesting \code{rname}, \code{strand}, \code{pos},
#' \code{qwidth}, \code{cigar}, \code{seq}, and the fields listed in \code{tags}.
#' To consider possible indels,
#'  \code{\link[GenomicAlignments]{cigarWidthAlongReferenceSpace}}
#'  is used to correct the \code{qwidth}.
#' @examples
#' \dontrun{
#' readBam("sample.bam")
#' lapply(c("sample1.bam", "sample2.bam"), readBam, tags = c("NM", "MD"))
#' }
#'
#' @importFrom Rsamtools scanBam ScanBamParam
#' @importFrom GenomicAlignments cigarWidthAlongReferenceSpace
#' @import data.table
#' @export
readBam <- function(bam_path, tags = character(0)) {

  # --- Validations -----------------------------------------------------
  if (!is.character(bam_path) || length(bam_path) != 1) {
    stop("`bam_path` must be a string.")
  }
  if (!is.character(tags)) {
    stop("`tags` must be a character vector.")
  }

  if (!file.exists(bam_path)) {
    stop("The following file do not exist:\n  ",
         bam_path)
  }

  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    stop("Package 'Rsamtools' is required but not installed.")
  }
  
  if (!requireNamespace("GenomicAlignments", quietly = TRUE)) {
    stop("Package 'GenomicAlignments' is required but not installed.")
  }

  # --- Read, format and merge BAM files ---------------------------------
  raw <- Rsamtools::scanBam(
      bam_path,
      param = Rsamtools::ScanBamParam(
        what = c("rname", "strand", "pos", "qwidth","cigar", "seq"),
        tag  = tags
      )
    )[[1]]

    raw$seq <- as.character(raw$seq)

    core <- data.table::as.data.table(raw[c("rname", "pos", "strand", "qwidth", "seq")])
    data.table::setnames(core, "rname", "seqnames")
    # Alignment span on the genome (accounts for insertions/deletions in the CIGAR)
    core[, width := GenomicAlignments::cigarWidthAlongReferenceSpace(raw$cigar)]
    if (length(tags) > 0) {
      core[, (tags) := raw$tag[tags]]
    }

    core[]
}

collapseReadStarts <- function(bam.dt){
  require(data.table)

  # assumes seqnames pos strand width seq and dc columns
  # Collapse reads starting at the same position
  # If the reads have the same start position (pos+width-1 for - strand),
  # on the same chromosome and strand
  # Keep only one with the highest number of duplicates and add the duplicates of others

  .dt <- copy(bam.dt)
  .dt[, tmp_start := fifelse(strand == "+", pos, pos + width - 1L)] # 1-based closed annotations

  # Sum of dc per group (vectorized, no .SD construction)
  .dt[, total_dc := sum(dc), by = .(seqnames, strand, tmp_start)]

  # Row index of the max-dc read per group (vectorized via .I)
  idx <- .dt[, .I[which.max(dc)], by = .(seqnames, strand, tmp_start)]$V1

  # Subset to just those rows, once
  .res <- .dt[idx]
  .res[, dc := total_dc]
  .res[, c("tmp_start", "total_dc") := NULL]

  .res[]
}

#' Group sequences by Hamming-distance similarity
#'
#' Assigns sequences to groups using a greedy, representative-based clustering
#' algorithm. The first ungrouped sequence becomes the representative of a new
#' group, and all remaining sequences within `max_hamming_dist` of that
#' representative are assigned to the same group.
#'
#' @param sequences A character vector containing sequences of equal length.
#'   Missing values are not permitted.
#' @param max_hamming_dist A non-negative integer specifying the maximum Hamming
#'   distance allowed between a sequence and the representative of its group.
#' @param nthread A positive integer specifying the number of threads passed to
#'   [stringdist::stringdist()]. Defaults to `2`.
#'
#' @details
#' This is a greedy clustering procedure. Group assignments can therefore depend on
#' the order of `sequences`. Recommend to sort reads by number of reads so the most
#' sequenced read is the representative for that group.
#'
#' Hamming distance is defined only for sequences of equal length. Consequently,
#' this function requires all supplied sequences to have the same number of
#' characters.
#'
#' @return An integer vector with one group label for each element of
#'   `sequences`. Labels are consecutive integers starting at `1`.
#'
#' @examples
#' sequences <- c("AAAA", "AAAT", "TTTT", "TTTA")
#' .GroupBySeqSimilarity(
#'   sequences,
#'   max_hamming_dist = 1,
#'   nthread = 1
#' )
#'
#' @keywords internal
.GroupBySeqSimilarity <- function(
    sequences,
    max_hamming_dist,
    nthread = 2L
) {
  if (!requireNamespace("stringdist", quietly = TRUE)) {
    stop(
      "Package 'stringdist' is required. Please install it first.",
      call. = FALSE
    )
  }

  if (!is.character(sequences)) {
    stop("`sequences` must be a character vector.", call. = FALSE)
  }

  if (anyNA(sequences)) {
    stop("`sequences` must not contain missing values.", call. = FALSE)
  }

  if (
    length(max_hamming_dist) != 1L ||
    is.na(max_hamming_dist) ||
    !is.numeric(max_hamming_dist) ||
    max_hamming_dist < 0 ||
    max_hamming_dist != as.integer(max_hamming_dist)
  ) {
    stop(
      "`max_hamming_dist` must be a single non-negative integer.",
      call. = FALSE
    )
  }

  if (
    length(nthread) != 1L ||
    is.na(nthread) ||
    !is.numeric(nthread) ||
    nthread < 1 ||
    nthread != as.integer(nthread)
  ) {
    stop("`nthread` must be a single positive integer.", call. = FALSE)
  }

  if (length(sequences) == 0L) {
    return(integer())
  }

  sequence_lengths <- nchar(sequences, type = "chars")

  if (length(unique(sequence_lengths)) != 1L) {
    stop(
      "All elements of `sequences` must have the same length.",
      call. = FALSE
    )
  }

  max_hamming_dist <- as.integer(max_hamming_dist)
  nthread <- as.integer(nthread)

  group_labels <- integer(length(sequences))
  ungrouped <- seq_along(sequences)
  current_group <- 0L

  while (length(ungrouped) > 0L) {
    current_group <- current_group + 1L
    representative_index <- ungrouped[[1L]]

    if (length(ungrouped) == 1L) {
      is_similar <- TRUE
    } else {
      distances <- stringdist::stringdist(
        sequences[representative_index],
        sequences[ungrouped[-1L]],
        method = "hamming",
        nthread = nthread
      )

      is_similar <- c(TRUE, distances <= max_hamming_dist)
    }

    matched_indices <- ungrouped[is_similar]
    group_labels[matched_indices] <- current_group
    ungrouped <- ungrouped[!is_similar]
  }

  group_labels
}

SlippageFilter <- function(bam.dt,
                           seq_length = 15,
                           hamming_dist = 3,
                           mingap = 0,
                           fragment_size = 30){

  # Dependencies
  devtools::load_all("/groups/stark/vloubiere/vlite/")

  N_allow_mismatch = hamming_dist
  N_bp_for_haming = seq_length

  # Remove reads which are shorter than the sequence which will be compared
  .c <- copy(bam.dt)
  .c <- .c[qwidth >= seq_length]

  # Resize to fragment size
  if(!is.na(fragment_size)){
    if(is.numeric(fragment_size) & fragment_size > 0){
      .c[strand == "+", `:=`(start = pos, end = pos+fragment_size-1) ]
      .c[strand == "-", `:=`(start = (pos+width-1)-fragment_size+1, end = pos+width-1) ]
    }
  }

  # Find contigs
  .c[, contig :=  vlite::collapseBed(.c, ignore.strand = F,  return.idx.only = T, min.gapwidth = mingap)]

  # contigs with more than 1 read
  .contigs_to_assess <- .c[, .N, by = contig][N > 1, contig]

  # Change the orientation of the - reads because Rsamtools:scanBam returns reverse complement for
  # the - strand
  .fix <- .c[contig %in% .contigs_to_assess]
  .fix[, tag := seq]
  .fix[strand == "-", tag := stringi::stri_reverse(seq)]
  .fix[, tag := stringi::stri_sub(tag, 1, N_bp_for_haming)]
  # Greedy clustering of sequences, take the most sequence read per group as representative
  .fix[order(-dc), Similarity_group := .GroupBySeqSimilarity(tag, N_allow_mismatch), contig]

  # Collapse
  .fix[, dc_sum := sum(dc), .(contig, Similarity_group)]
  idx <- .fix[, .I[which.max(dc)], by = .(contig, Similarity_group)]$V1
  .fix <- .fix[idx]
  .fix[, dc := dc_sum]
  .fix[, c("Similarity_group", "dc_sum", "tag") := NULL]

  # combine
  .out <- rbind(.fix, .c[!(contig %in% .contigs_to_assess)])

  .out[]
}
