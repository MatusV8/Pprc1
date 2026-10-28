#---- Function to extract EPD .dat data ----
read_epd_dat <- function(file) {

  lines <- readLines(file)
  codes <- substr(lines, 1L, 2L)
  n     <- length(lines)

  # --- Entry boundaries ---
  end_idx <- which(trimws(lines) == "//")
  n_entries <- length(end_idx)
  start_idx <- c(1L, head(end_idx + 1L, -1L))

  # --- Pre-index all fields at once ---
  field_idx <- function(code) which(codes == code)
  field_val <- function(idx) trimws(substr(lines[idx], 4L, nchar(lines[idx])))

  idx  <- lapply(c("ID","GN","DE","OS","ME","DT","SE","FP","DR"), field_idx)
  vals <- lapply(idx, field_val)
  names(idx) <- names(vals) <- c("ID","GN","DE","OS","ME","DT","SE","FP","DR")

  # --- Pre-parse DR lines ---
  dr_raw   <- vals$DR
  dr_db <- trimws(sub(";.*", "", dr_raw))
  dr_value <- sub("\\.$", "", trimws(sub("^[^;]+;\\s*", "", dr_raw)))

  known_dbs <- c("EnsemblT","GENOME","UCSC","Ensembl","RefSeq","Gene Symbol","EPD")

  # --- Helper: get line indices for entry i using binary search (O(log n)) ---
  entry_range <- function(field_line_nums, s, e) {
    lo <- findInterval(s - 1L, field_line_nums) + 1L
    hi <- findInterval(e,       field_line_nums)
    if (lo > hi) integer(0L) else lo:hi
  }

  first_val <- function(field, s, e) {
    pos <- entry_range(idx[[field]], s, e)
    if (length(pos) == 0L) NA_character_ else vals[[field]][pos[1L]]
  }

  all_vals <- function(field, s, e) {
    pos <- entry_range(idx[[field]], s, e)
    if (length(pos) == 0L) NA_character_ else vals[[field]][pos]
  }

  # --- Main loop ---
  result <- vector("list", n_entries)

  for (i in seq_len(n_entries)) {
    s <- start_idx[i]
    e <- end_idx[i]

    # ID
    id_raw   <- first_val("ID", s, e)
    id_split   <- strsplit(trimws(id_raw), "\\s{2,}")[[1]]   # 2+ spaces separate name from classification
    entry_name <- id_split[1L]
    if (length(id_split) >= 2L) {
      classif      <- trimws(strsplit(gsub("\\.$", "", id_split[2L]), ";")[[1]])
      data_class   <- if (length(classif) >= 1L) trimws(classif[1L]) else NA_character_
      tss_type     <- if (length(classif) >= 2L) trimws(classif[2L]) else NA_character_
      tax_division <- if (length(classif) >= 3L) trimws(classif[3L]) else NA_character_
    } else {
      data_class <- tss_type <- tax_division <- NA_character_
    }

    # GN
    gn_raw    <- first_val("GN", s, e)
    gene_name <- if (!is.na(gn_raw)) gsub(";$", "", sub("^Name=", "", gn_raw)) else NA_character_

    # Simple fields
    de_all      <- all_vals("DE", s, e)
    description <- if (!all(is.na(de_all))) paste(de_all, collapse = " ") else NA_character_
    organism    <- first_val("OS", s, e)
    method      <- first_val("ME", s, e)
    date        <- first_val("DT", s, e)
    sequence    <- first_val("SE", s, e)
    fp          <- first_val("FP", s, e)

    # DR fields — O(log n) range, then simple vector ops on small subset
    dr_pos    <- entry_range(idx$DR, s, e)
    sub_db    <- dr_db[dr_pos]
    sub_val   <- dr_value[dr_pos]

    get_db <- function(db) {
      v <- sub_val[sub_db == db]
      if (length(v) == 0L) NA_character_ else paste(unique(v), collapse = " | ")
    }

    dr_ensemblt    <- get_db("EnsemblT")
    dr_genome      <- get_db("GENOME")
    dr_ucsc        <- get_db("UCSC")
    dr_ensembl     <- get_db("Ensembl")
    dr_refseq      <- get_db("RefSeq")
    dr_gene_symbol <- get_db("Gene Symbol")
    dr_epd         <- get_db("EPD")
    dr_other       <- { v <- sub_val[!sub_db %in% known_dbs]
    if (length(v) == 0L) NA_character_ else paste(v, collapse = " | ") }

    # Coordinates from UCSC DR
    chr <- strand <- NA_character_
    tss_position  <- NA_integer_
    if (!is.na(dr_ucsc)) {
      p <- strsplit(dr_ucsc, "\\s+")[[1]]
      if (length(p) >= 4L) {
        chr          <- p[2L]
        strand       <- p[3L]
        tss_position <- suppressWarnings(as.integer(p[4L]))
      }
    }

    result[[i]] <- list(
      entry_name     = entry_name,
      data_class     = data_class,
      tss_type       = tss_type,
      tax_division   = tax_division,
      gene_name      = gene_name,
      description    = description,
      organism       = organism,
      method         = method,
      date           = date,
      chr            = chr,
      strand         = strand,
      tss_position   = tss_position,
      sequence       = sequence,
      fp             = fp,
      dr_ensemblt    = dr_ensemblt,
      dr_genome      = dr_genome,
      dr_ucsc        = dr_ucsc,
      dr_ensembl     = dr_ensembl,
      dr_refseq      = dr_refseq,
      dr_gene_symbol = dr_gene_symbol,
      dr_epd         = dr_epd,
      dr_other       = dr_other
    )
  }

  rbindlist(result, fill = TRUE)
}
