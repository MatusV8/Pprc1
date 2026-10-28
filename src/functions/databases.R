#' Retrieve EnsDb object from AnnotationHub
#'
#' @param species Character string specifying species (e.g., "Mus musculus", "Homo sapiens")
#' @param assembly Character string specifying genome assembly (e.g., "mm10", "mm39", "hg38", "hg19")
#' @param version Integer specifying Ensembl version (e.g., 102, 79)
#' @return An EnsDb object
#' @examples
#' # Get latest mouse EnsDb
#' ensdb <- get_ensdb(species = "Mus musculus")
#'
#' # Get latest EnsDb for mm10
#' ensdb <- get_ensdb(assembly = "mm10")
#'
#' # Get specific version for mouse
#' ensdb <- get_ensdb(species = "Mus musculus", version = 102)
get_ensdb <- function(assembly = NULL, species = NULL, version = NULL) {

  # Load AnnotationHub
  if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
    stop("Package 'AnnotationHub' is required. Install with: BiocManager::install('AnnotationHub')")
  }

  # Load EnsDb metadata
  ah <- AnnotationHub::AnnotationHub()
  ensdb <- AnnotationHub::query(ah, "EnsDb")

  # Check if assembly or species are valid
  if(!is.character(assembly) && !is.character(species)){
    stop("assembly or species need to be specified.")
  }


  # Assembly to species/version mapping
  assembly_map <- list(
    # Mouse assemblies
    "mm39" = list(species = "Mus musculus", assembly = "GRCm39", version = NULL),
    "mm10" = list(species = "Mus musculus", assembly = "GRCm38", version = 102),
    "mm9" = list(species = "Mus musculus", assembly = "MGSCv37", version = 67),

    # Human assemblies
    "hg38" = list(species = "Homo sapiens", assembly = "GRCh38", version = NULL),
    "hg19" = list(species = "Homo sapiens", assembly = "GRCh37", version = 75),

    # Rat assemblies
    "rn7" = list(species = "Rattus norvegicus", assembly = "mRatBN7.2", max_version = NULL),
    "rn6" = list(species = "Rattus norvegicus", assembly = "Rnor_6.0", max_version = 103)
  )

  # Handle assembly first
  assembly_info <- list(species = NULL, assembly = NULL, max_version = NULL)
  if (!is.null(assembly)) {
    if (assembly %in% names(assembly_map)) {
      assembly_info <- assembly_map[[assembly]]
    } else if(!(assembly %in% AnnotationHub::mcols(ensdb)$genome)){
      stop(sprintf("Assembly '%s' not recognized. Supported: %s or ensembl assembly IDs e.g. 'GRCm38'",
                   assembly, paste(names(assembly_map), collapse = ", ")))
    } else{
      assembly_info$assembly <- assembly
      assembly_info$species <- AnnotationHub::mcols(ensdb)$species[AnnotationHub::mcols(ensdb)$genome == assembly][1]
    }

    if (!is.null(species) && !identical(species, assembly_info$species)) {
      warning(sprintf("species = '%s' ignored; using species implied by assembly '%s' ('%s')",
                      species, assembly, assembly_info$species))
    }
  } else {
    # Else we use species only information
    assembly_info$species <- species
  }

  # Validate that species is a valid objext
  if(!(assembly_info$species %in% AnnotationHub::mcols(ensdb)$species )){
    stop("Invalid species name")
  }

  # Build query terms
  query_terms <- c(
    if(is.null(assembly_info$assembly)) assembly_info$species else assembly_info$assembly,
       assembly_info$version)

      # Get query ID
      query_result <- AnnotationHub::query(ensdb, query_terms)

    if(length(names(query_result)) == 1){
      ensdb_id <- names(query_result)
    } else if(length(names(query_result)) > 1){

      .versions<- as.integer(gsub(".*\\b (\\d{2,3}) \\b.*", "\\1", query_result$title))
      names(.versions) <- names(query_result)
      ensdb_id <- names(.versions)[which.max(.versions)]
    } else if (length(query_result) == 0) {
      stop(sprintf("No EnsDb found for: %s", paste(query_terms, collapse = ", ")))
    }
    message(sprintf("Retrieving EnsDb: %s (%s)",
                    ensdb_id,
                    query_result[ensdb_id]$title))

    # Retrieve and return the EnsDb object
    ensdb <- ah[[ensdb_id]]

    return(ensdb)
}

#============================================================================#
# GO related functions
#============================================================================#


#' Get GO Annotations for a Given Organism
#'
#' Retrieves Gene Ontology (GO) term annotations for genes in a specified
#' organism annotation database (\code{OrgDb}), joins them with human-readable
#' GO term descriptions from \code{GO.db}, and returns a tidy data frame.
#'
#' @param orgdb An \code{OrgDb} object (e.g. \code{org.Mm.eg.db::org.Mm.eg.db},
#'   \code{org.Hs.eg.db::org.Hs.eg.db}, or an \code{OrgDb} object retrieved via
#'   \code{AnnotationHub}). Must support the specified \code{keytype} and
#'   provide the \code{"GO"} and \code{"ONTOLOGY"} columns.
#' @param keytype Character string specifying the key type to use for querying
#'   \code{orgdb} (e.g. \code{"ENSEMBL"}, \code{"ENTREZID"}, \code{"SYMBOL"}).
#'   Defaults to \code{"ENSEMBL"}. Must be one of the values returned by
#'   \code{AnnotationDbi::keytypes(orgdb)}.
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{gene_id}{Gene identifier, named according to \code{keytype}.}
#'     \item{ID}{GO term identifier (GOID).}
#'     \item{type}{GO ontology category (\code{"BP"}, \code{"CC"}, or \code{"MF"}).}
#'     \item{description}{Human-readable description of the GO term.}
#'   }
#'
#' @details
#' The function validates that \code{orgdb} supports the requested
#' \code{keytype} and contains the required \code{"GO"} and \code{"ONTOLOGY"}
#' columns before querying. Rows with missing values are dropped, and duplicate
#' rows are removed. GO term descriptions are retrieved separately from
#' \code{GO.db} and merged in via a left join on GO identifiers.
#'
#' @examples
#' \dontrun{
#' library(org.Mm.eg.db)
#' GO_mouse <- get_GOs(org.Mm.eg.db, keytype = "ENSEMBL")
#'
#' library(org.Hs.eg.db)
#' GO_human <- get_GOs(org.Hs.eg.db, keytype = "ENSEMBL")
#' }
#'
#' @seealso \code{\link[AnnotationDbi]{select}}, \code{\link[AnnotationDbi]{keytypes}},
#'   \code{\link[AnnotationDbi]{columns}}
#'
#' @importFrom AnnotationDbi select keys keytypes columns
#' @importFrom dplyr select distinct rename left_join
#' @importFrom tidyr drop_na
#'
#' @export
get_GOs <- function(orgdb, keytype = "ENSEMBL") {
  require(AnnotationDbi)
  require(GO.db)
  require(dplyr)
  require(tidyr)

  # Validate that orgdb supports the requested keytype and columns
  valid_keytypes <- AnnotationDbi::keytypes(orgdb)
  if (!keytype %in% valid_keytypes) {
    stop(sprintf(
      "keytype '%s' not supported by this orgdb. Available keytypes: %s",
      keytype, paste(valid_keytypes, collapse = ", ")
    ))
  }

  valid_columns <- AnnotationDbi::columns(orgdb)
  required_cols <- c("GO", "ONTOLOGY")
  missing_cols <- setdiff(required_cols, valid_columns)
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "orgdb is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  GO <- AnnotationDbi::select(
    orgdb,
    keys = AnnotationDbi::keys(orgdb, keytype = keytype),
    columns = c(keytype, "GO", "ONTOLOGY", "SYMBOL"),
    keytype = keytype
  )

  # Remove EVIDENCE, drop NAs, dedupe
  GO <- GO %>%
    dplyr::select(-EVIDENCE) %>%
    tidyr::drop_na() %>%
    dplyr::distinct() %>%
    dplyr::rename(gene_id = !!keytype,
                  gene_name = SYMBOL,
                  ID = GO,
                  type = ONTOLOGY)

  # Retrieve GO term descriptions
  extra <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = GO$ID,
    keytype = "GOID",
    columns = "TERM"
  ) %>%
    tidyr::drop_na() %>%
    dplyr::distinct() %>%
    dplyr::rename(description = TERM)

  # Merge GO terms with descriptions
  GO <- dplyr::left_join(GO, extra, by = c("ID" = "GOID")) %>%
    as_tibble()
  return(GO)
}

#' Compute GO Term Background Statistics for a Gene Set
#'
#' Filters a GO annotation table to a given set of gene IDs and computes,
#' for each GO term (\code{type} x \code{description} combination), the
#' number of genes in the input set annotated with that term, alongside
#' the total number of unique genes queried.
#'
#' @param gene_ids Character vector of gene identifiers to filter on. Should
#'   match the values in the \code{gene_id} column of \code{GO_table}.
#' @param GO_table A data frame or tibble containing GO annotations, such as
#'   the output of \code{\link{getGOs}}. Must contain the columns
#'   \code{gene_id}, \code{type}, and \code{description}.
#'
#' @return A tibble filtered to rows where \code{gene_id \%in\% gene_ids},
#'   with two additional columns:
#'   \describe{
#'     \item{anno_count}{Number of unique genes (from \code{gene_ids}) annotated
#'       with each \code{type}/\code{description} combination.}
#'     \item{genes_count}{Total number of unique genes supplied in \code{gene_ids}.}
#'   }
#'
#' @details
#' \code{anno_count} and \code{genes_count} can be used together to compute
#' the proportion of input genes annotated with a given GO term (e.g. for
#' background/foreground comparisons in enrichment analysis).
#'
#' @examples
#' \dontrun{
#' GO_table <- getGOs(org.Mm.eg.db, keytype = "ENSEMBL")
#' bg <- get_GO_bg(gene_ids = c("ENSMUSG00000001", "ENSMUSG00000002"), GO_table)
#' }
#'
#' @importFrom dplyr filter mutate
#' @importFrom tibble as_tibble
#'
#' @export
get_GO_bg <- function(gene_ids, GO_table) {

  required_cols <- c("gene_id", "type", "description")
  missing_cols <- setdiff(required_cols, colnames(GO_table))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "GO_table is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  gene_ids <- unique(gene_ids)

  GO_table %>%
    tibble::as_tibble() %>%
    dplyr::filter(gene_id %in% gene_ids) %>%
    dplyr::mutate(
      anno_count = length(unique(gene_id)),
      genes_count = length(gene_ids),
      .by = c(type, description)
    )
}
