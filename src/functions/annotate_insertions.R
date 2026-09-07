annotateInsertions <- function(bed.dt,
                               gtf,
                               compatible_introns,
                               Phase = NA){

  #-----------------------------------------------------------#
  # The function returns annotated data.table with these columns:

  # 1. Index
  # 2. Chromosome name
  # 3. Start
  # 4. End
  # 5. Strand
  # 6. dc: Number of PCR duplicates
  # 7. element type
  # 8. gene type
  # 9. Sense
  # 10. Intron InFrame: logical if based on the phase of the construct integrations is compatible
  # 11. Phase
  # 12. transcript_id
  # 13. gene_id
  # 14. gene_name


  #-----------------------------------------------------------#


  # There might be these issues:
  # 1. duplicated rows (multiple introns, same gene, same phase)
  #     Resolved by keeping unique rows

  # 2. different phase (multiple introns, same gene, different phase)
  #     Resolved by filtering rows with only correct phase (single phase samples) or ignoring phase and
  #      removing duplicated rows

  # 3. different gene (overlapping introns, different gene, same phase)
  #     Resolved by keeping all

  # 4. different gene different phase (overlapping introns, different gene, different phase)
  #     Resolved by filtering rows with only correct phase (single phase samples) or
  #     by sorting genes by transcript support level and keeping all


  #--------------------------------------#
  # Read in file
  #--------------------------------------#
  BED <- copy(bed.dt)
  Compatible_introns <- copy(compatible_introns)
  anno <- copy(gtf)
  biotype_levels <- c(
    "protein_coding",
    "IG_C_gene",
    "TR_C_gene",
    "polymorphic_pseudogene",
    "lincRNA",
    "macro_lncRNA",
    "bidirectional_promoter_lncRNA",
    "3prime_overlapping_ncRNA",
    "antisense",
    "sense_intronic",
    "sense_overlapping",
    "processed_transcript",
    "transcribed_unprocessed_pseudogene",
    "transcribed_processed_pseudogene",
    "transcribed_unitary_pseudogene",
    "unprocessed_pseudogene",
    "processed_pseudogene"
  )

  anno[, gene_biotype :=  factor(gene_biotype, levels = biotype_levels, ordered = TRUE)]
  anno <- anno[type %in% c("transcript", "exon", "five_prime_utr", "trhee_prime_utr")]
  anno[, type := factor(type, levels = c("exon", "five_prime_utr", "three_prime_utr", "transcript"), ordered = TRUE)]

  #--------------------------------------#
  # Output tibble with annotated integrations
  #--------------------------------------#
  # Add annotation columns
  BED[, idx := paste(seqnames, start, strand, sep = ":")]

  overlaps <- intersectRanges(BED, Compatible_introns, mode = "all", ignore.strand = T)
  overlaps <- unique(overlaps, by = c("idx","phase", "gene_id")) # get rid of redundant annotations (intron for the same gene, same phase, same integration)
  overlaps[, Sense := strand == b.strand]
  overlaps[, InFrame := if(is.na(Phase)) {Sense} else {(Phase == phase) & Sense}] # If phase is NA, than every phase is InFrame

  # Select all In-frame annotations (if not clear which gene to assign, then assign to both)
  setorder(overlaps,idx, -InFrame, -Sense, transcript_support_level)
  InframeIns <- overlaps[(InFrame), .SD[1], .(idx, gene_id)]

  # For out of phase insertion keep only one gene name if multiple match (prioritise in-frame and sense annotation and then pick the one with the best transcript support level)
  InframeIns <- rbind(InframeIns, overlaps[!(idx %in% InframeIns$idx), .SD[1], .(idx)])

  # The rest
  others <- intersectRanges(BED[!(idx %in% InframeIns$idx)],
                            anno,
                            mode = "all", return_index = F, ignore.strand = T)
  others[,  Sense := strand == b.strand]
  # Keep only one annotation per insertion
  setorder(others, type, -Sense, gene_biotype, transcript_support_level)

  others <- others[, .SD[1], .(idx)]
  others[, type := ifelse(type == "transcript", "non_coding_intron", type)]

  # Combine
  bed.anno <- rbindlist(
    list(
      InframeIns[, .(idx,
                     seqnames,
                     start,
                     end,
                     strand,
                     dc,
                     element = "intron",
                     gene_type = "protein_coding",
                     Sense,
                     InFrame,
                     phase,
                     transcript_id,
                     gene_id,
                     gene_name)],

      others[, .(idx,
                     seqnames,
                     start,
                     end,
                     strand,
                     dc,
                     element = type,
                     gene_type = gene_biotype,
                     Sense,
                     InFrame = F,
                     phase,
                     transcript_id,
                     gene_id,
                     gene_name)],

      BED[!(idx %in% c(others$idx, InframeIns$idx)),.(idx,
                                                      seqnames,
                                                      start,
                                                      end,
                                                      strand,
                                                      dc,
                                                      element = "intergenic",
                                                      gene_type = NA,
                                                      Sense = F,
                                                      InFrame = F,
                                                      phase = NA,
                                                      transcript_id = NA,
                                                      gene_id = NA,
                                                      gene_name = NA) ]
    )
  )



  return(bed.anno)
}
