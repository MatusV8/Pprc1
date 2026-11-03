
# Get GO annotations
getGOs <- function() {

  go2all <- as.data.table(AnnotationDbi::toTable(org.Mm.egGO2ALLEGS))
  ens    <- unique(as.data.table(AnnotationDbi::toTable(org.Mm.egENSEMBL)))
  
  GO <- merge(ens, go2all, by = "gene_id")
  
  go_desc <- as.data.table(AnnotationDbi::select(
    GO.db::GO.db, keys = unique(GO$go_id), keytype = "GOID",
    columns = c("TERM", "DEFINITION")))
  GO <- merge(GO, go_desc, by.x = "go_id", by.y = "GOID", all.x = TRUE)
  
  evidence_levels <- c(
    # Experimental
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
    # High-throughput experimental
    "HTP", "HDA", "HMP", "HGI", "HEP",
    # Phylogenetic / computational, curator-reviewed
    "IBA", "IBD", "IKR", "IRD",
    "ISS", "ISO", "ISA", "ISM", "IGC", "RCA",
    # Author / curator statements
     "IC", "NAS","TAS",
    # Electronic and no data
    "IEA", "ND"
  )
  
  GO[, Evidence := factor(Evidence, levels = evidence_levels, ordered = TRUE)]
  setorder(GO, Evidence)
  GO[]
}

go_enrichment <- function(genes, hit, GO_table, min_term_size = 10, max_term_size = 200, padj_method = "BH" ) {
  
  # Restrict the universe to the intersect
  keep  <- genes %in% unique(GO_table$ensembl_id)
  genes <- genes[keep]
  hit   <- hit[keep]
  
  annotation_genes = GO_table$ensembl_id
  annotation_go = GO_table$go_id
  
  genes <- as.character(genes)
  if (!is.logical(hit) || anyNA(hit)) stop("hit must be logical without NA.")
  if (anyNA(genes) || anyDuplicated(genes)) stop("genes must be unique and non-missing.")
  
  N <- length(genes)
  H <- sum(hit)
  U <- N - H
  
  idx <- match(as.character(annotation_genes), genes)
  go  <- as.character(annotation_go)
  ok  <- !is.na(idx) & !is.na(go) & nzchar(go)
  pairs <- unique(data.frame(idx = idx[ok], GO = go[ok]))
  
  f <- factor(pairs$GO)
  t <- tabulate(f)
  a <- rowsum(as.integer(hit[pairs$idx]), f)[, 1]
  
  b <- H - a
  c <- t - a
  d <- U - c
  
  p <- phyper(a - 1, m = t, n = N - t, k = H, lower.tail = FALSE)
  
  res <- data.table(
    go_id = levels(f),
    term_size = t,
    GO_hit = a,
    log2OR    = log2((a * d) / (b * c)),                       # Inf if b or c = 0, NaN if 0/0
    log2OR_hc = log2(((a + .5) * (d + .5)) / ((b + .5) * (c + .5))),  # Haldane-Anscombe
    p.value = p,
    row.names = NULL
  )
  res[term_size >= min_term_size & term_size <= max_term_size, padj := p.adjust(p.value, method = padj_method)]
  go_info <- unique(GO_table[, !"ensembl_id"], by = "go_id")
  res <- go_info[res, on = "go_id"]
  setorder(res, padj, na.last = T)
  res[]
}


get_descendant_genes <- function(GO_ID, GOs_df) {
  requireNamespace("GO.db")
  GOs_df <- as.data.table(GOs_df)
  
  # Get all descendant GO terms (children, grandchildren, etc.)
  descendants <- c(
    tryCatch(GO.db::GOBPOFFSPRING[[GO_ID]],  error = function(e) NULL), # Use GOMFOFFSPRING or GOCCOFFSPRING for MF/CC
    tryCatch(GO.db::GOMFOFFSPRING[[GO_ID]], error = function(e) NULL),
             tryCatch(GO.db::GOCCOFFSPRING[[GO_ID]], error = function(e) NULL)
  )
  
  # Combine parent term with all descendants
  all_GO_terms <- c(GO_ID, descendants)
  
  # Remove NAs
  all_GO_terms <- all_GO_terms[!is.na(all_GO_terms)]
  
  # Get all unique gene_ids for these GO terms
  gene_ids <- GOs_df[GOs_df$ID %in% all_GO_terms, gene_id] %>%
    unique()
  
  return(gene_ids)
}