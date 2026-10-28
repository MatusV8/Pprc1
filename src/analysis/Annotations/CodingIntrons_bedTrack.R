GeneTrap_Intron_path <- fread("db/annotations/GRCm38.102_GeneTrapCompatibleCodingIntrons.txt.gz")

GeneTrap_Intron_path[, .(
  seqnames,
  start,
  end,
  name = as.character(phase),
  score = 1L,
  strand,
  thickStart = start,
  thickEnd = end,
  itemRgb = data.table::fcase(
    phase == 0, "#3E8DED",  # 62,141,237
    phase == 1, "#298C15",  # 41,140,21
    phase == 2, "#ED9911",  # 237,153,17
    default = "#808080"
  )
)] |>
  rtracklayer::export(
    "db/annotations/GRCm38.102_GeneTrapCompatibleCodingIntrons.bed",
    format = "BED"
  )
