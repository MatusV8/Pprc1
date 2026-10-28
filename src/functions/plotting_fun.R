mv_pie <- function(tib, category, value, fill_colors = c("gray","#4DBBD5", "#E64B35" , "#00A087", "#3C5488", "#F39B7F", "#8491B4"),
                   alpha = 1, bar_color = NULL,
                   label = "b", label_position = "inside", label_size = 7, legend.key.size.pt = 7, label_units = "pt") {
  
  require(dplyr)
  require(ggplot2)
  
  # Validate input data structure
  if (!inherits(tib, c("data.frame", "tbl_df", "tbl", "data.table"))) {
    stop("'tib' must be a data.frame, tibble, or data.table")
  }
  
  # Validate category parameter
  if (!is.character(category) || length(category) != 1) {
    stop("'category' must be a character vector of length 1")
  }
  if (!category %in% names(tib)) {
    stop(sprintf("Column '%s' not found in data", category))
  }
  
  # Validate value parameter
  if (!is.character(value) || length(value) != 1) {
    stop("'value' must be a character vector of length 1")
  }
  if (!value %in% names(tib)) {
    stop(sprintf("Column '%s' not found in data", value))
  }
  if (!is.numeric(tib[[value]])) {
    stop(sprintf("Column '%s' must be numeric", value))
  }
  
  # Validate label parameter
  valid_labels <- c("v", "p", "b")
  if (!label %in% valid_labels) {
    stop(sprintf("'label' must be one of: %s", paste(valid_labels, collapse = ", ")))
  }
  
  # Validate and set label position
  if (is.numeric(label_position)) {
    lab_pos <- label_position
  } else if (is.character(label_position) && length(label_position) == 1) {
    lab_pos <- switch(label_position,
                      "inside" = 1,
                      "outside" = 1.7,
                      stop(sprintf("'label_position' must be 'inside', 'outside', or a numeric value")))
  } else {
    stop("'label_position' must be 'inside', 'outside', or a numeric value")
  }
  
  # Calculate percentages
  tib_plot <- tib %>%
    mutate(
      .value = .data[[value]],
      .category = .data[[category]],
      .percentage = .value / sum(.value) * 100
    )
  
  # Create label based on label parameter
  tib_plot <- tib_plot %>%
    mutate(.label = case_when(
      label == "v" ~ as.character(.value),
      label == "p" ~ sprintf("%.1f%%", .percentage),
      label == "b" ~ sprintf("%d\n(%.1f%%)", .value, .percentage)
    ))
  
  # Create plot
  p <- tib_plot %>%
    ggplot(aes(x = "", y = .value, fill = .category)) +
    geom_bar(stat = "identity", width = 1, color = bar_color, alpha = alpha) +
    theme_void() +
    geom_text(aes(label = .label, x = lab_pos), 
              position = position_stack(vjust = 0.5),
              size = label_size, size.unit = label_units) +
    coord_polar("y", start = 0) +
    labs(fill = category) +
    theme(plot.margin = margin(-10, 10, -10, -10),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -10),  # negative value pulls legend closer
          legend.box.spacing = unit(0, "pt"),
          legend.key.size = unit(legend.key.size.pt, "pt")) +
    scale_y_continuous(expand = expansion(0, 0)) +
    scale_x_discrete(expand = expansion(0, 0))
  
  # Apply custom colors if provided
  if (!is.null(fill_colors)) {
    p <- p + scale_fill_manual(values = fill_colors)
  }
  
  return(p)
}

#' Create a heatmap with optional column and/or row metadata annotations
#'
#' @param mat         A numeric matrix. Row and column names are used to derive
#'                    annotations when \code{show_row_anno} or \code{show_col_anno}
#'                    are \code{TRUE}.
#' @param col_labels  Character vector of annotation label names for columns.
#'                    Must match the number of fields produced by splitting column
#'                    names on \code{sep}. If \code{NULL} (default), labels are
#'                    set to \code{"1"}, \code{"2"}, ... automatically.
#' @param row_labels  Character vector of annotation label names for rows.
#'                    Must match the number of fields produced by splitting row
#'                    names on \code{sep}. If \code{NULL} (default), labels are
#'                    set to \code{"1"}, \code{"2"}, ... automatically.
#' @param colors      Optional. Either:
#'                    \itemize{
#'                      \item A named list of character vectors, one per label.
#'                      \item A single character vector recycled across all labels.
#'                      \item \code{NULL} (default) to use a built-in colour palette.
#'                    }
#' @param sep         Character used to split row/column names into annotation
#'                    fields. Default: \code{"_"}.
#' @param hm_name     Legend title for the heatmap colour scale. Default: \code{"PCC"}.
#' @param col_breaks  Numeric vector of length 3 defining the colour scale breakpoints.
#'                    Default: \code{c(0.4, 0.7, 1)}.
#' @param col_colors  Character vector of length 3 defining the colours at each
#'                    breakpoint. Default: \code{c("blue", "white", "red")}.
#' @param row_title   Row title string. Default: \code{"Samples"}.
#' @param col_title   Column title string. Default: \code{"Input samples"}.
#' @param show_col_anno Logical. Whether to build and display a top annotation
#'                    by splitting column names on \code{sep}. Default: \code{TRUE}.
#' @param show_row_anno Logical. Whether to build and display a left annotation
#'                    by splitting row names on \code{sep}. Default: \code{FALSE}.
#' @param show_values Logical. Whether to print cell values inside the heatmap.
#'                    Default: \code{TRUE}.
#' @param value_digits Integer. Number of decimal places for cell value labels.
#'                    Default: \code{1}.
#' @param value_fontsize Numeric. Font size for cell value labels. Default: \code{8}.
#' @param names_fontsize Numeric. Font size for row/column name labels. Default: \code{10}.
#' @param ...         Additional arguments passed to \code{ComplexHeatmap::Heatmap()}.
#'
#' @return A \code{Heatmap} object from the \pkg{ComplexHeatmap} package.
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar grid.text
#' @importFrom data.table as.data.table tstrsplit
mv_hm <- function(
    mat,
    col_labels       = NULL,
    row_labels       = NULL,
    colors           = NULL,
    sep              = "_",
    hm_name          = "PCC",
    col_breaks       = c(0.4, 0.7, 1),
    col_colors       = c("blue", "white", "red"),
    row_title        = NULL,
    col_title        = NULL,
    show_col_anno    = TRUE,
    show_row_anno    = FALSE,
    show_values      = TRUE,
    cluster_rows        = TRUE,
    cluster_columns     = TRUE,
    show_row_names      = FALSE,
    show_column_names   = FALSE,
    value_digits     = 1,
    value_fontsize   = 8,
    names_fontsize   = 10,
    ...
) {
  
  # ---- Input checks ----------------------------------------------------------
  
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("`mat` must be a numeric matrix.")
  }
  if (!is.null(col_labels) && !is.character(col_labels)) {
    stop("`col_labels` must be a character vector or NULL.")
  }
  if (!is.null(row_labels) && !is.character(row_labels)) {
    stop("`row_labels` must be a character vector or NULL.")
  }
  if (!is.logical(show_col_anno) || length(show_col_anno) != 1) {
    stop("`show_col_anno` must be a single logical value.")
  }
  if (!is.logical(show_row_anno) || length(show_row_anno) != 1) {
    stop("`show_row_anno` must be a single logical value.")
  }
  if (show_col_anno && is.null(colnames(mat))) {
    stop("`mat` must have column names when `show_col_anno = TRUE`.")
  }
  if (show_row_anno && is.null(rownames(mat))) {
    stop("`mat` must have row names when `show_row_anno = TRUE`.")
  }
  if (length(col_breaks) != 3 || length(col_colors) != 3) {
    stop("`col_breaks` and `col_colors` must each be of length 3.")
  }
  if (!is.numeric(col_breaks) || is.unsorted(col_breaks)) {
    stop("`col_breaks` must be a sorted numeric vector.")
  }
  
  # ---- Default colour palette ------------------------------------------------
  
  default_palette <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62"
  )
  
  # ---- Helper: resolve NULL labels from first element of dim names -----------
  
  .default_labels <- function(dim_names) {
    n <- length(data.table::tstrsplit(dim_names[[1]], sep, fixed = TRUE))
    as.character(seq_len(n))
  }
  
  # ---- Resolve labels --------------------------------------------------------
  
  if (show_col_anno && is.null(col_labels)) {
    col_labels <- .default_labels(colnames(mat))
  }
  if (show_row_anno && is.null(row_labels)) {
    row_labels <- .default_labels(rownames(mat))
  }
  
  # ---- Helper: split dim names into an annotation data.table -----------------
  
  .build_anno_dt <- function(dim_names, dim_label, labels) {
    dt <- data.table::as.data.table(
      data.table::tstrsplit(dim_names, sep, fixed = TRUE)
    )
    if (ncol(dt) != length(labels)) {
      stop(
        "Splitting ", dim_label, " names by '", sep, "' yields ", ncol(dt),
        " fields, but `", dim_label, "_labels` has length ", length(labels), "."
      )
    }
    names(dt) <- labels
    dt
  }
  
  # ---- Helper: resolve colors for a given annotation data.table --------------
  
  .resolve_colors <- function(anno, colors, labels) {
    if (is.null(colors)) {
      col_l <- lapply(labels, function(x) {
        lvls <- unique(anno[[x]])
        if (length(lvls) > length(default_palette)) {
          stop("More unique values in '", x, "' than available default colours. ",
               "Please supply a `colors` argument.")
        }
        setNames(default_palette[seq_along(lvls)], lvls)
      })
      names(col_l) <- labels
      
    } else if (is.list(colors)) {
      if (!all(labels %in% names(colors))) {
        stop("When `colors` is a list it must be named with all entries in the relevant labels.")
      }
      col_l <- colors[labels]
      for (x in labels) {
        lvls <- unique(anno[[x]])
        # Auto-assign level names if the vector is unnamed
        if (is.null(names(col_l[[x]]))) {
          if (length(col_l[[x]]) < length(lvls)) {
            stop("Color vector for label '", x, "' has ", length(col_l[[x]]),
                 " colours but ", length(lvls), " levels.")
          }
          col_l[[x]] <- setNames(col_l[[x]][seq_along(lvls)], lvls)
        } else {
          missing_lvls <- setdiff(lvls, names(col_l[[x]]))
          if (length(missing_lvls) > 0) {
            stop("Color list for label '", x, "' is missing entries for: ",
                 paste(missing_lvls, collapse = ", "), ".")
          }
        }
      }
      
    } else if (is.character(colors)) {
      col_l <- lapply(labels, function(x) {
        lvls <- unique(anno[[x]])
        if (length(lvls) > length(colors)) {
          stop("Not enough colors supplied for label '", x, "': need ",
               length(lvls), ", got ", length(colors), ".")
        }
        setNames(colors[seq_along(lvls)], lvls)
      })
      names(col_l) <- labels
      
    } else {
      stop("`colors` must be NULL, a named list, or a character vector.")
    }
    col_l
  }
  
  # ---- Helper: build a HeatmapAnnotation from a data.table ------------------
  
  .build_heatmap_anno <- function(anno, col_l, labels, which = c("column", "row")) {
    which     <- match.arg(which)
    anno_args <- lapply(labels, function(x) anno[[x]])
    names(anno_args) <- labels
    anno_args$col                  <- col_l
    anno_args$which                <- which
    anno_args$annotation_name_side <- if (which == "column") "left" else "bottom"
    anno_args$annotation_name_gp       <- grid::gpar(fontsize = names_fontsize)
    anno_args$annotation_legend_param  <- list(
      labels_gp = grid::gpar(fontsize = names_fontsize),
      title_gp  = grid::gpar(fontsize = names_fontsize+1)
    )
    anno_args$simple_anno_size <- unit(names_fontsize * 0.4, "mm")
    do.call(ComplexHeatmap::HeatmapAnnotation, anno_args)
  }
  
  # ---- Build column annotation -----------------------------------------------
  
  top_anno <- NULL
  if (show_col_anno) {
    col_anno_dt <- .build_anno_dt(colnames(mat), "column", col_labels)
    
    if(any(is.na(col_anno_dt))){
      stop("NAs are present in some of the categories.")
    }
    
    col_col_l   <- .resolve_colors(col_anno_dt, colors, col_labels)
    top_anno    <- .build_heatmap_anno(col_anno_dt, col_col_l, col_labels, which = "column")
  }
  
  # ---- Build row annotation --------------------------------------------------
  
  left_anno <- NULL
  if (show_row_anno) {
    row_anno_dt <- .build_anno_dt(rownames(mat), "row", row_labels)
    row_col_l   <- .resolve_colors(row_anno_dt, colors, row_labels)
    left_anno   <- .build_heatmap_anno(row_anno_dt, row_col_l, row_labels, which = "row")
  }
  
  # ---- Heatmap colour scale --------------------------------------------------
  
  col_fun <- circlize::colorRamp2(col_breaks, col_colors)
  
  # ---- Cell label function ---------------------------------------------------
  
  fmt <- paste0("%.", value_digits, "f")
  
  cell_fn <- if (show_values) {
    function(j, i, x, y, width, height, fill) {
      grid::grid.text(
        sprintf(fmt, mat[i, j]), x, y,
        gp = grid::gpar(fontsize = value_fontsize, col = "black")
      )
    }
  } else {
    NULL
  }
  
  # ---- Draw heatmap ----------------------------------------------------------
  
  hm <- ComplexHeatmap::Heatmap(
    mat,
    name                = hm_name,
    col                 = col_fun,
    cluster_rows        = cluster_rows,
    cluster_columns     = cluster_columns,
    show_row_names      = show_row_names,
    show_column_names   = show_column_names,
    row_names_gp        = grid::gpar(fontsize = names_fontsize),
    column_names_gp     = grid::gpar(fontsize = names_fontsize),
    row_title           = row_title,
    column_title        = col_title,
    top_annotation      = top_anno,
    left_annotation     = left_anno,
    cell_fun            = cell_fn,
    heatmap_legend_param = list(
      labels_gp = grid::gpar(fontsize = names_fontsize),
      title_gp  = grid::gpar(fontsize = names_fontsize+1)
    ),
    ...
  )
  
  return(hm)
}

#====================================================================#
# Helper functions ----
#====================================================================#
#' Parse a genomic region specification
#'
#' Accepts a region in one of several convenient forms and normalizes it to a
#' plain list with numeric coordinates.
#'
#' @param region One of:
#'   * a \code{GRanges} object of length 1,
#'   * a single string \code{"chr:start-end"} (commas in the numbers are
#'     tolerated, e.g. \code{"chr1:1,000-2,000"}),
#'   * a length-3 vector \code{c(chr, start, end)}.
#'
#' @return A list with elements \code{seqnames} (character scalar),
#'   \code{start} (numeric scalar) and \code{end} (numeric scalar).
#'
#' @keywords internal
#' @noRd
.parse_region <- function(region) {
  # ---- GRanges input --------------------------------------------------------
  if (methods::is(region, "GRanges")) {
    if (length(region) != 1L) {
      stop("`region` GRanges must have exactly one range.", call. = FALSE)
    }
    return(list(
      seqnames = as.character(GenomicRanges::seqnames(region)),
      start    = as.numeric(GenomicRanges::start(region)),
      end      = as.numeric(GenomicRanges::end(region))
    ))
  }
  
  # ---- Single "chr:start-end" string ---------------------------------------
  if (is.character(region) && length(region) == 1L) {
    # Strip thousands separators for forgiving parsing.
    clean <- gsub(",", "", region, fixed = TRUE)
    m <- regmatches(clean, regexec("^([^:]+):([0-9]+)-([0-9]+)$", clean))[[1]]
    if (length(m) != 4L) {
      stop(
        "Could not parse region string. Expected 'chr:start-end', got: ",
        region, call. = FALSE
      )
    }
    return(list(
      seqnames = m[2],
      start    = as.numeric(m[3]),
      end      = as.numeric(m[4])
    ))
  }
  
  # ---- Length-3 vector c(chr, start, end) ----------------------------------
  if (length(region) == 3L) {
    return(list(
      seqnames = as.character(region[[1]]),
      start    = as.numeric(region[[2]]),
      end      = as.numeric(region[[3]])
    ))
  }
  
  stop(
    "`region` must be a GRanges, a 'chr:start-end' string, or c(chr,start,end).",
    call. = FALSE
  )
}

#' Greedy interval packing into rows (levels)
#'
#' Assigns each interval the lowest integer "row" such that it does not overlap
#' any interval already placed in that row (allowing a small \code{gap}). This
#' is used to lay out transcripts on as few rows as possible.
#'
#' @param starts Integer/numeric vector of interval start coordinates.
#' @param ends   Integer/numeric vector of interval end coordinates.
#' @param gap    Minimum horizontal gap required between two intervals sharing a
#'   row. Defaults to 0.
#'
#' @return An integer vector (same length as \code{starts}) of 1-based row levels.
#'
#' @keywords internal
#' @noRd
.pack_intervals <- function(starts, ends, gap = 0) {
  n <- length(starts)
  if (n == 0L) return(integer(0))
  
  # Process intervals left-to-right by start coordinate for a stable, tight pack.
  ord <- order(starts, ends)
  levels_out <- integer(n)
  # `row_ends[k]` = right-most occupied coordinate on row k.
  row_ends <- numeric(0)
  
  for (i in ord) {
    placed <- FALSE
    if (length(row_ends) > 0L) {
      # First row whose last interval ends before this one starts (minus gap).
      free <- which(row_ends + gap < starts[i])
      if (length(free) > 0L) {
        k <- free[1L]
        row_ends[k] <- ends[i]
        levels_out[i] <- k
        placed <- TRUE
      }
    }
    if (!placed) {
      # Open a new row.
      row_ends <- c(row_ends, ends[i])
      levels_out[i] <- length(row_ends)
    }
  }
  
  levels_out
}

#====================================================================#
# Main functions ----
#====================================================================#

plot_signal_tracks <- function(bw,
                               region,
                               nbins = 1000L,
                               fill_color = "black",
                               fill_alpha = 0.4,
                               free_y     = TRUE,
                               show_x_axis = TRUE,
                               strip_position = "right",
                               x_axis_position = "top",
                               y_title = "CPM",
                               y_expand = 0.1,
                               show_baseline = T,
                               lwd = 0.25,
                               show_legend = FALSE
)
{
  # Parse genomic ranges
  r <- .parse_region(region)
  
  # Create x-axis label
  x_title = paste0(r$seqnames,":", r$start,"-", r$end)
  
  # Create a gr object
  .gr <- GenomicRanges::GRanges(
    seqnames = r$seqnames,
    ranges   = IRanges::IRanges(start = r$start, end = r$end)
  )
  
  # Load BW
  .dt <- purrr::map_dfr(bw, ~{
    .bwf <- rtracklayer::BigWigFile(.x)
    summary(.bwf, .gr, size = nbins, type = "mean")[[1]] %>% as.data.table()
  }, .id = "sample")
  
  setDT(.dt)
  
  # Keep the order of the samples
  .dt[, sample := factor(sample, names(bw))]
  
  # Get plotting position
  .dt[, position := floor((end+start)/2)]
  
  # Check colors
  if (length(fill_color) == 1) {
    .fill_color <- rep(fill_color, length(bw))
  } else if (length(fill_color) == length(bw)) {
    .fill_color <- fill_color
  } else {
    stop("Number of colors doesn't match number of samples.", call. = FALSE)
  }
  
  # plot
  .p <- .dt %>%
    ggplot(aes(x = position, y = score, fill = sample)) +
    theme_classic() +
    theme(strip.background = element_blank()) +
    geom_area(alpha = fill_alpha, show.legend = show_legend) +
    scale_fill_manual(values = .fill_color) +
    scale_y_continuous(expand = expansion(mult = c(0, y_expand)),
                       breaks = function(limits) c(0, limits[2]),
                       labels = scales::label_number(drop0trailing = TRUE)) +
    scale_x_continuous(expand = c(0,0),
                       labels = scales::label_comma(),
                       position = x_axis_position) +
    facet_wrap(~ sample,ncol = 1,  scales = ifelse(free_y, "free_y", "fixed"), strip.position = strip_position) +
    labs(x = x_title,
         y = y_title)
  
  # Add horizontal line 
  if(show_baseline){
    .p <- .p +
      geom_hline(yintercept = 0, linewidth = lwd)
  }
  
  .p
}

plot_annotation_track <- function(gtf_dt,
                                  region,
                                  #color_by       = NULL,
                                  colors         = "black",
                                  track_name     = NULL,
                                  show_labels    = TRUE,
                                  label_size     = 2,
                                  label_padding = 0.5,
                                  padding_nrow  = 2,
                                  exon_height    = 0.5,
                                  body_lwd       = 0.25,
                                  show_x_axis    = FALSE,
                                  arrow_height    = "auto",
                                  arrow_width = 100,
                                  arrow_bin      = 1000,
                                  min_intron_length = 200,
                                  base_size      = 11) {
  
  
  if(arrow_height == "auto"){
    arrow_height <- exon_height*0.6
  } else if(!is.numeric(arrow_height)){
    stop("arrow_height must be 'auto' or numeric")
  }
  
  # Check chromosome formating 
  reg <- .parse_region(region)
  if(!(reg$seqnames %in% gtf_dt$seqnames)){
    stop("The chromosome doesn't match the reference.")
  }
  
  .dt <- copy(gtf_dt)
  
  # 1. Keep only features overlapping the region -----------------------------
  .dt[, seqnames := as.character(seqnames)]
  .dt <- .dt[seqnames == reg$seqnames & end > reg$start & start < reg$end]
  
  # 2. Clip coordinates to the region ----------------------------------------
  .dt[, `:=`(start = pmax(start, reg$start), end = pmin(end, reg$end))]
  
  # 3. Add row levels  ----------------------------------------
  tx_rows <- .dt[type == "transcript"]
  tx_rows[, row_level := .pack_intervals(start, end, gap = 0)]
  .dt[tx_rows, row_level := i.row_level, on = "transcript_id"]
  .dt[, row_level := min(row_level), transcript_id ] # join back
  
  # 4. Calculate hights for exons and UTRs ----------------------------------------
  .dt[, h := fcase(
    type == "CDS", exon_height,
    type == "exon", exon_height * 0.6, default = NA )
  ]
  .dt[, `:=`(ymin = row_level - h / 2, ymax = row_level + h / 2)]
  
  # 5. Get introns ----------------------------------------
  setorder(.dt, start, end) # resort annotation
  introns <- .dt[type == "exon", .(start = end[-(.N)]+1, end = start[-1]-1), by = .(seqnames,strand, transcript_id, gene_name, gene_id, row_level)]
  
  # 6. Calculate arrow polygons for introns ----------------------------------------
  arrow_dt <- introns[(end - start) > min_intron_length, {
    
    if(.N >0){
      pos <- seq(start, end, by = arrow_bin)
      pos <- pos[-c(1L, length(pos))]   # trim endpoints (was 'post <- post[...]')
      if (length(pos) == 0L) {
        pos <- (start + end) / 2      # midpoint if no interior positions
      }
      .tipx <- if (strand[1] == "+") pos + arrow_width else pos - arrow_width
      .tail <- if (strand[1] == "+") pos - arrow_width else pos + arrow_width
      # 3 vertices per triangle: top-tail, tip, bottom-tail
      data.table::data.table(
        grp = rep(seq_along(pos), each = 3L),
        x   = as.vector(rbind(.tail, .tipx, .tail)),
        y   = as.vector(rbind(row_level + arrow_height/2, row_level, row_level - arrow_height/2))
      )
    } else(
      data.table::data.table(
        grp = NA,
        x   = NA,
        y   = NA
      )
    )
    
    
  }, by = names(introns)]
  arrow_dt[, poly_id := .GRP, by = .(transcript_id, grp, start, end)]
  
  # 7. Labels
  lab <- tx_rows[,.SD[1], transcript_id]
  lab[, xmid := (start + end) / 2]
  
  .dt %>% 
    ggplot() +
    theme_void() +
    scale_x_continuous(limits = c(reg$start, reg$end), expand = c(0,0)) +
    # Add arrows
    ggplot2::geom_path(data = arrow_dt,
                       ggplot2::aes(x = x, y = y, group = poly_id), linewidth = body_lwd*0.6) +
    # Add introns
    ggplot2::geom_segment(data = . %>% dplyr::filter(type == "transcript"),
                          ggplot2::aes(x = start, xend = end, y = row_level, yend = row_level), linewidth = body_lwd) +
    # Add CDS and UTRs
    ggplot2::geom_rect(data = . %>% dplyr::filter(type != "transcript"),
                       ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax), fill = colors) +
    
    # Add gene names
    #shadowtext::geom_shadowtext(data = lab,
    #                ggplot2::aes(x = xmid, y = row_level, label = gene_name),
    #                size = label_size, 
    #                colour = "grey20",    # The text color
    #                bg.colour = "white",  # The glow/halo color
    #                bg.r = 0.1) +          # The radius of the glow
    
    
    geom_text(data = lab,
              ggplot2::aes(x = xmid, y = row_level -label_padding, label = gene_name),
              size = label_size, color = "black", vjust = 1) +
    scale_y_continuous(limits = c(1-exon_height-label_padding-padding_nrow, max(.dt$row_level)+exon_height/2), expand = c(0,0))
  
}

bwPlotAverage <- function(bed,
                          by= NULL,
                          tracks,
                          center= "center",
                          center.name= ifelse(center!="region", "Center", "Start"),
                          region.start= "Start",
                          region.end= "End",
                          upstream= 5000L,
                          downstream= 5000L,
                          nbins= if(center=="region") c(50L, 150L, 50L) else 251L,
                          ignore.strand= FALSE,
                          colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4"),
                          plot= TRUE,
                          normalisation = "none",
                          xlab= "Genomic distance",
                          ylab= "Mean signal",
                          legend_title = "sample",
                          cleanup.cache= FALSE)
{
  
  # Checks ----
  if(!(normalisation %in% c("none", "max", "sum"))){
    stop("Invalid normalisation. Pick one of 'none' 'max' 'sum'")
  }
  
  names <- if(is.null(names(tracks))){
    gsub("(.*)[.].*$", "\\1", basename(tracks))
  } else {
    names(tracks)
  }
  
  if(any(names == "")) {
    stop("Track names cannot be empty. Please provide valid names.")
  }
  
  if(!is.null(by))
    names <- paste0(rep(names, each= length(unique(by))), " @ ", sort(unique(by)))
  
  # Compute signal if necessary (cache file) ----
  signal <- bwBinnedCoverage(bed = bed,
                             tracks = tracks,
                             center = center,
                             upstream = upstream,
                             downstream = downstream,
                             nbins = nbins,
                             ignore.strand= ignore.strand,
                             cleanup.cache = cleanup.cache)
  signal <- lapply(signal, as.data.table)
  
  # Scale to 1, max of means across positions
  if(normalisation == "max"){
    signal <- map(signal, function(df) {
      col_means <- apply(df, 2, mean, na.rm = TRUE)
      df / max(col_means, na.rm = TRUE)
    })
  }
  
  if(normalisation == "sum"){
    signal <- map(signal, function(df) {
      col_means <- apply(df, 2, mean, na.rm = TRUE)
      df / sum(col_means, na.rm = TRUE)
    })
  }
  
  # Split groups using by ----
  if(!is.null(by)) {
    signal <- lapply(signal, function(x) split(x, by))
    signal <- unlist(signal, recursive = FALSE, use.names = FALSE)
  }
  
  # Compute mean signal, standard error, bin.x.pos and x.pos ----
  Nregions <- sapply(signal, nrow)
  signalMean <- lapply(signal, apply, 2, mean, na.rm= TRUE)
  signalSe <- lapply(signal, function(x) {
    apply(x, 2, function(y) {
      sd(y, na.rm = TRUE)/sqrt(length(na.omit(y)))
    })
  })
  bin.x.pos <- lapply(signal, function(x) as.numeric(colnames(x)))
  
  # Format data table ----
  dat <- data.table(Nregions,
                    signalMean,
                    signalSe,
                    bin.x.pos)
  
  # Add names and colors ----
  dat[, name := factor(names, names)]
  
  # Melt data ----
  dat <- dat[, lapply(.SD, unlist), name]
  
  .p <- dat %>% 
    ggplot(aes(x = bin.x.pos)) + 
    theme_bw() +
    geom_ribbon(aes(ymin = signalMean - signalSe,
                    ymax = signalMean + signalSe, fill = name), alpha = 0.3) +
    geom_line(aes(y = signalMean, color = name)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = xlab, y = ylab, color = legend_title, fill = legend_title) 
  
  if(center == "region"){
    .p <- .p +
      scale_x_continuous(breaks = c(-upstream, 0, dat[nbins[1]+nbins[2], bin.x.pos], max(dat$bin.x.pos)),
                         labels = c(-upstream, region.start, region.end, downstream))
  }
  .p
  
}

Dendroplot_vsaNorm <- function(dds){
  vst_mat <- DESeq2::vst(dds, blind = TRUE) # Stabilize variance
  mat <- SummarizedExperiment::assay(vst_mat)
  cor <- Hmisc::varclus(mat, similarity = "pearson", trans = "none")
  
  plot(
    cor,
    hang = -1,
    xlab = NA,
    las = 1,
    ylab = "PCC"
  )
  
  # Add colored points for each sample
  labs <- factor(gsub("(.*)_[1-9]$", "\\1", colnames(mat)))
  Cc <- rainbow(nlevels(labs))[labs]
  points(
    1:ncol(mat),
    rep(0, ncol(mat)),
    col = Cc[cor$hclust$order],
    pch = 19
  )
}