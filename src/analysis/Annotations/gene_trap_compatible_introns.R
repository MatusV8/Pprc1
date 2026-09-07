#=========================================#
# Libraries       ====
#=========================================#
require(rtracklayer)
require(GenomicRanges)
require(data.table)

#=========================================#
# Download GRCm38 v102 GTF file ====
#=========================================#
#GTF_file <- "Mus_musculus.GRCm38.102.gtf.gz" # name your file
#mm10v102_url <- "https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"
#
#if(!file.exist(GTF_file)){
#
#  download.file(mm10v102_url, GTF_file)
#
#}

#=========================================#
# Get tagging compatible introns ====
#=========================================#
GeneTrap_Intron_path = "db/annotations/GRCm38.102_GeneTrapCompatibleCodingIntrons.txt.gz"
if(!file.exists(GeneTrap_Intron_path)) {

  # Read in annotations for the assembly
  Annotation_GR <- rtracklayer::import(anno_path) %>% as.data.table()
  Annotation_GR[ ,transcript_support_level := as.numeric(gsub(" .*", "", transcript_support_level))] # Fix transcript support level

  # get coding exons
  coding_exons <- Annotation_GR[tag == "basic" & type == "CDS" & seqnames %in% c(1:19, "X", "Y")] |> unique( by = c("seqnames", "strand", "start", "end", "protein_id", "gene_name"))
  coding_exons[, exon_number := as.numeric(coding_exons$exon_number)]

  # get introns
  setorderv(coding_exons, "exon_number", 1)
  # Compute coordinates between exons for each transcript isoform = introns
  ps <- coding_exons[strand=="+",
                     .(phase= phase[-1],
                       start= end[-(.N)]+1,
                       end= start[-1]-1,
                     intron_number = exon_number[-(.N)]),
                     .(seqnames, strand, transcript_id, gene_id, gene_name, protein_id, transcript_support_level)]

  setorderv(coding_exons, "exon_number", -1)
  ns <- coding_exons[strand=="-", .(phase= phase[-(.N)],
                                    start= end[-(.N)]+1,
                                    end= start[-1]-1,
                     intron_number = exon_number[-1]),
                     .(seqnames, strand, transcript_id, gene_id, gene_name, protein_id, transcript_support_level)]
  introns <- rbind(ps, ns)
  introns[, width:= end-start+1]

  # Collapse introns per gene
  setorder(introns, transcript_support_level)
  introns[, contig :=  vlite::collapseBed(.SD, ignore.strand = F,  return.idx.only = T), by = phase] # Find overlapping introns
  introns_collapsed <- introns[, .(
    start = min(start),
    end = max(end),
    transcript_id = stringi::stri_paste(unique(transcript_id), collapse = ","),
    transcript_support_level = stringi::stri_paste(unique(transcript_support_level), collapse = ","),
    gene_name = stringi::stri_paste(unique(gene_name), collapse = ","),
    protein_id = stringi::stri_paste(unique(protein_id), collapse = ",")
  ), keyby = .(seqnames, strand, phase, gene_id, contig)]
  introns_collapsed[, contig := NULL]
  setorder(introns_collapsed, seqnames, start,end,strand)

  # Save as a text file, tab separated
  fwrite(introns_collapsed, GeneTrap_Intron_path, sep = "\t", compress = "gzip")
}

