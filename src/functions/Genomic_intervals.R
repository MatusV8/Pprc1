#' Intersect Two Sets of Genomic Ranges
#'
#' @description
#' Finds overlaps between two sets of genomic ranges (\code{a} and \code{b})
#' and returns either the overlapping/non-overlapping subset of \code{a}, or
#' every overlapping (a, b) pair, depending on \code{mode}.
#'
#' @param a A \code{GRanges} object, or a \code{data.frame}/\code{data.table}/
#'   \code{tibble} coercible to \code{GRanges} via
#'   \code{\link[GenomicRanges]{makeGRangesFromDataFrame}} (needs
#'   \code{seqnames}, \code{start}, \code{end}, and optionally \code{strand}).
#' @param b Same as \code{a}; the "subject" ranges compared against \code{a}.
#' @param mode One of \code{"overlap"}, \code{"nonoverlap"}, or \code{"all"}:
#'   \describe{
#'     \item{\code{"overlap"}}{Rows of \code{a} overlapping >=1 range in \code{b}.}
#'     \item{\code{"nonoverlap"}}{Rows of \code{a} overlapping no range in \code{b}.}
#'     \item{\code{"all"}}{Every overlapping (a, b) pair (one row per hit,
#'       inner-join style). Shared column names in b get \code{"b."}
#'       prefixes.}
#'   }
#' @param return_index Logical; if \code{TRUE} return integer indices
#'   (\code{"overlap"}/\code{"nonoverlap"}) or an \code{a_index}/\code{b_index}
#'   \code{data.table} (\code{"all"}) instead of subsetted data.
#' @param ignore.strand Logical; see \code{\link[GenomicRanges]{findOverlaps}}
#' @param ... Passed to \code{\link[GenomicRanges]{findOverlaps}} /
#'   \code{\link[GenomicRanges]{countOverlaps}} (e.g. \code{maxgap},
#'   \code{minoverlap}, \code{type}, \code{ignore.strand}).
#'
#' @return A \code{data.table}, or integer vector / index \code{data.table}
#'   when \code{return_index = TRUE}.
#'
#' @details
#' Overlap semantics are fully controlled via \code{...}. Note
#' \code{ignore.strand} defaults to \code{FALSE}; pass
#' \code{ignore.strand = TRUE} explicitly if strand should be ignored.
#'
#' @export
intersectRanges <- function(a, b,
                            mode = c("overlap", "nonoverlap", "all"),
                            return_index = FALSE, ignore.strand = FALSE, ...) {

  mode <- match.arg(mode)

  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required but not installed.")
  }

  .asGRanges <- function(x, argname) {
    if (methods::is(x, "GRanges")) return(x)
    if (is.data.frame(x)) {
      return(GenomicRanges::makeGRangesFromDataFrame(
        as.data.frame(x), keep.extra.columns = TRUE
      ))
    }
    stop(sprintf("`%s` must be a GRanges or a data.frame/data.table/tibble.", argname))
  }

  # Fast GRanges -> data.table conversion (avoids slow as.data.frame.GRanges).
  # Always returns a fresh, safely-mutable data.table.
  .toDT <- function(x) {
    if (methods::is(x, "GRanges")) {
      dt <- data.table::as.data.table(S4Vectors::mcols(x))
      dt[, `:=`(
        seqnames = as.character(GenomicRanges::seqnames(x)),
        start    = GenomicRanges::start(x),
        end      = GenomicRanges::end(x),
        strand   = as.character(GenomicRanges::strand(x))
      )]
      data.table::setcolorder(dt, c("seqnames", "start", "end", "strand"))
      dt
    } else if (data.table::is.data.table(x)) {
      data.table::copy(x)          # avoid mutating caller's data.table
    } else {
      data.table::as.data.table(x)
    }
  }

  a_gr <- .asGRanges(a, "a")
  b_gr <- .asGRanges(b, "b")

  # Match seqlevel styles
  if (!any(GenomeInfoDb::seqlevelsStyle(a_gr) %in% GenomeInfoDb::seqlevelsStyle(b_gr))) GenomeInfoDb::seqlevelsStyle(b_gr) <- GenomeInfoDb::seqlevelsStyle(a_gr)

  if (mode %in% c("overlap", "nonoverlap")) {
    has_hit <- GenomicRanges::countOverlaps(a_gr, b_gr, ignore.strand = ignore.strand, ...) > 0
    idx <- if (mode == "overlap") which(has_hit) else which(!has_hit)

    if (isTRUE(return_index)) return(idx)
    return(.toDT(a)[idx])
  }

  # mode == "all" -> every a/b pair that overlaps
  hits   <- GenomicRanges::findOverlaps(a_gr, b_gr, ignore.strand = ignore.strand, ...)
  a_hits <- S4Vectors::queryHits(hits)
  b_hits <- S4Vectors::subjectHits(hits)

  if (isTRUE(return_index)) {
    return(data.table::data.table(a_index = a_hits, b_index = b_hits))
  }

  a_dt <- .toDT(a)
  b_dt <- .toDT(b)

  shared_names <- intersect(names(a_dt), names(b_dt))
  if (length(shared_names) > 0) {
    data.table::setnames(b_dt, old = shared_names, new = paste0("b.", shared_names))
  }

  cbind(a_dt[a_hits], b_dt[b_hits])
}
