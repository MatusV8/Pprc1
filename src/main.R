if(F){

  # 01. Annotations
  file.edit("src/analysis/Annotations/gene_trap_compatible_introns.R")
  file.edit("src/analysis/Annotations/Retrieve_FANTOM5_CAGE.R")
  file.edit("src/analysis/Annotations/EPDnew_promoters.R")
  file.edit("src/analysis/Annotations/GRCm38v102_PROseq_annotations.R")
  
  # 02. Degron screen analysis
  file.edit("src/analysis/AIDscreen/01_AIDscreen_pipeline.R")
  file.edit("src/analysis/AIDscreen/02_Data_parsing.R")
  file.edit("src/analysis/AIDscreen/03_Enrichment_CMH.R")

  # 03. ChIP-seq, Cut&Run analysis
  file.edit("src/analysis/ChIP/01_ChIP_demultiplexing.R")
  file.edit("src/analysis/ChIP/02_ChIP_processing.R")
  file.edit("src/analysis/ChIP/03_ChIP_peakCalling.R")
  file.edit("src/analysis/ChIP/04_ChIP_Pprc1confidentPeaks.R")
  
  # 04. PRO-seq analysis
  file.edit("src/analysis/PROseq/01_Demultiplex.R")
  file.edit("src/analysis/PROseq/02_Alignment.R")
  file.edit("src/analysis/PROseq/03_Counting.R")
  file.edit("src/analysis/PROseq/04_DESeq_enrichment.R")
  file.edit("src/analysis/PROseq/05_PROseq_CNR_combined.R")
  
  # 05. Figures
  file.edit("src/analysis/Fig1.R")
  file.edit("src/analysis/Fig2.R")
  file.edit("src/analysis/Fig3.R")
  file.edit("src/analysis/Fig5.R")
  
  file.edit("src/analysis/FigS1.R")
  file.edit("src/analysis/FigS2.R")
  file.edit("src/analysis/FigS3.R")
  file.edit("src/analysis/FigS4.R")
  file.edit("src/analysis/FigS8.R")

  # 06. GEO submission
  file.edit("src/GEO_submission.R")

}
