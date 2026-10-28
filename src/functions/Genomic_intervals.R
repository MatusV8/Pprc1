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
  a_style <- GenomeInfoDb::seqlevelsStyle(a_gr)[1]
  b_style <- GenomeInfoDb::seqlevelsStyle(b_gr)
  
  if (!any(a_style %in% b_style)) {
    tryCatch(
      GenomeInfoDb::seqlevelsStyle(b_gr) <- a_style,
      error = function(e) {
        warning("Could not harmonize seqlevels style between `a` and `b`: ", conditionMessage(e),
                ". Proceeding without renaming; make sure seqnames match manually.")
      }
    )
  }

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

#' Split genomic intervals into fixed-width bins
#'
#' `binBed()` takes a set of genomic intervals (a BED-like table) and splits
#' every interval into consecutive, fixed-width bins. It accepts a
#' [data.table::data.table()], a [base::data.frame()], or a
#' `GenomicRanges::GRanges` object and always returns a
#' [data.table::data.table()].
#'
#' Each input interval of width `w = end - start + 1` is divided into
#' `ceiling(w / bins.width)` bins (at least one). Bin `k` (1-based) spans
#' `start + (k - 1) * bins.width` to
#' `min(start + k * bins.width - 1, end)`, so the final bin of an interval is
#' truncated at the original `end` and may be narrower than `bins.width`.
#'
#' A `bin.idx` column is added giving the 1-based position of each bin within
#' its source interval. For intervals on the negative strand (`strand == "-"`)
#' the bin index is reversed, so `bin.idx == 1` always corresponds to the 5'
#' end of the feature. If the input has no `strand` column, all bins are
#' numbered in ascending coordinate order.
#'
#' @param bed A `data.table`, `data.frame`, or `GRanges` object describing
#'   genomic intervals. It must contain (at least) integer `start` and `end`
#'   columns. For `GRanges` input the `seqnames` column is renamed to `chrom`
#'   and the `width` column is dropped. All other columns are carried through
#'   to the output unchanged (repeated for every bin of the interval).
#' @param bins.width A single positive number giving the desired bin width in
#'   base pairs.
#'
#' @return A `data.table` with the same columns as the (normalised) input plus
#'   a `bin.idx` column, containing one row per bin. The input object is never
#'   modified (copy semantics).
#'
#' @section Overwrite behaviour:
#' If the input already contains a `bin.idx` column it is overwritten with the
#' computed bin index. Any other columns are preserved.
#'
#' @examples
#' library(data.table)
#'
#' bed <- data.table(
#'   chrom  = c("chr1", "chr2"),
#'   start  = c(1L, 100L),
#'   end    = c(250L, 149L),
#'   strand = c("+", "-")
#' )
#' binBed(bed, bins.width = 100L)
#'
#' # data.frame input works too
#' binBed(as.data.frame(bed), bins.width = 100L)
#'
#' @importFrom data.table as.data.table copy is.data.table set setcolorder
#' @export
binBed <- function(bed, bins.width) {
  ## ---- validate bins.width -------------------------------------------------
  if (!is.numeric(bins.width) || length(bins.width) != 1L ||
      is.na(bins.width) || bins.width < 1) {
    stop("`bins.width` must be a single positive number.", call. = FALSE)
  }
  
  ## ---- normalise input to a fresh data.table -------------------------------
  if (inherits(bed, "GRanges")) {
    .dt <- data.table::as.data.table(as.data.frame(bed))
    if ("seqnames" %in% names(.dt)) {
      data.table::setnames(.dt, "seqnames", "chrom")
    }
    if ("width" %in% names(.dt)) {
      data.table::set(.dt, j = "width", value = NULL)
    }
    ## as.data.frame() turns strand into a factor; use plain characters so the
    ## downstream comparison strand == "-" behaves predictably.
    if ("strand" %in% names(.dt)) {
      data.table::set(.dt, j = "strand", value = as.character(.dt$strand))
    }
  } else if (data.table::is.data.table(bed)) {
    .dt <- data.table::copy(bed)
  } else if (is.data.frame(bed)) {
    .dt <- data.table::as.data.table(bed)
  } else {
    stop("`bed` must be a data.table, data.frame, or GRanges object.",
         call. = FALSE)
  }
  
  if (!all(c("start", "end") %in% names(.dt))) {
    stop("`bed` must contain `start` and `end` columns.", call. = FALSE)
  }
  
  ## Column names *after* normalisation (Bug 2 fix): never rely on names(bed).
  orig.names <- names(.dt)
  has.strand <- "strand" %in% orig.names
  
  n.src <- nrow(.dt)
  if (n.src == 0L) {
    ## Nothing to bin: return an empty table with the expected columns.
    out.cols <- c(setdiff(orig.names, "bin.idx"), "bin.idx")
    if (!"bin.idx" %in% names(.dt)) {
      data.table::set(.dt, j = "bin.idx", value = integer(0))
    }
    data.table::setcolorder(.dt, out.cols)
    return(.dt[])
  }
  
  ## ---- compute per-interval bin counts and expand --------------------------
  n.bins <- pmax(1L, as.integer(ceiling((.dt$end - .dt$start + 1) / bins.width)))
  idx     <- rep.int(seq_len(n.src), n.bins)   # source row for each output row
  bin.idx <- sequence(n.bins)                  # 1..n1, 1..n2, ... (ascending)
  
  res <- .dt[idx]
  
  ## Add bin.idx as a genuine column *first* so the subsequent start/end
  ## expressions unambiguously reference columns (Bug 4 fix).
  data.table::set(res, j = "bin.idx", value = bin.idx)
  data.table::set(res, j = "start",
                  value = res$start + (res$bin.idx - 1L) * bins.width)
  data.table::set(res, j = "end",
                  value = pmin(res$start + bins.width - 1L, res$end))
  
  ## ---- reverse bin index on the negative strand (vectorised) ---------------
  if (has.strand) {
    neg <- which(res$strand == "-")
    if (length(neg)) {
      n.bins.exp <- rep.int(n.bins, n.bins)  # per output row: bins in its interval
      data.table::set(
        res, i = neg, j = "bin.idx",
        value = n.bins.exp[neg] - res$bin.idx[neg] + 1L
      )
    }
  }
  
  ## ---- assemble output -----------------------------------------------------
  out.cols <- c(setdiff(orig.names, "bin.idx"), "bin.idx")
  data.table::setcolorder(res, out.cols)
  ## Drop any helper columns that are not part of the requested output.
  extra <- setdiff(names(res), out.cols)
  if (length(extra)) {
    for (col in extra) data.table::set(res, j = col, value = NULL)
  }
  res[]
}

