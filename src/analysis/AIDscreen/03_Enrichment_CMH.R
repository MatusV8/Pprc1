#================================#
# Setup ----
#================================#
# Libraries
library(data.table)
library(tidyverse)
source("src/functions/CMH.R")

PARAMS <- list(
  dc_cutoff = 10,
  gene_count_cutoff = 1,
  padj_cutoff = 0.01,
  l2OR_cutoff = 0
)

#================================#
# Parse input ----
#================================#
Counts <- fread("db/count_tables/AID_screen_counts.csv")

#=====================================#
# Enrichment over Input using CMH ----
#=====================================#
cmh_enr <- data.table(sample = colnames(Counts)[grepl("GFPlow", colnames(Counts))])
cmh_enr[, input := gsub("GFPlow", "Input", sample)]
cmh_enr <- cmh_enr[, {
  .dt <- Counts[, c("gene_id","gene_name", .BY$sample, .BY$input), with = FALSE]
  names(.dt) <- c("gene_id","gene_name", "GFPlow_count", "Input_count")
  
  .dt[, GFPlow_total := sum(GFPlow_count)]
  .dt[, Input_total := sum(Input_count)]
}, by = .(sample, input)]

cmh_enr[, reporter := gsub("(Atp5k|EF1a).*", "\\1", sample)]

# remove empty measurements
cmh_enr <- cmh_enr[GFPlow_count >= PARAMS$gene_count_cutoff | Input_count >= PARAMS$gene_count_cutoff]

cmh_enr1 <-  cmh_enr[, {
  K <- .N
  
  if (K == 0) {
    list(
      n_strata         = 0L,
      max_GFPlow_count = NA_real_,
      max_Input_count  = NA_real_,
      logOR            = NA_real_,
      SE_logOR         = NA_real_,
      p_value          = NA_real_,
      method           = NA_character_
    )
    
  } else {
    # Array construction and OR/SE estimation are identical for k >= 1
    .arr <- build_array(GFPlow_count, GFPlow_total, Input_count, Input_total)
    est  <- mh_rbg(.arr)
    
    # P-value source depends only on whether there is one stratum or many
    if (K == 1) {
      method <- "Fisher_exact"
      p <- tryCatch(
        fisher.test(.arr[, , 1], alternative = "greater")$p.value,
        error = function(e) {
          print(paste("Fisher's exact test failed at:", .BY$reporter, .BY$gene_name))
          NA_real_
        })
      if (is.na(p)) method <- NA_character_
      
    } else {
      method <- "CMH_exact"
      p <- tryCatch(
        mantelhaen.test(.arr, alternative = "greater", exact = TRUE)$p.value,
        error = function(e) NULL)
      
      if (is.null(p)) {
        method <- "CMH_approx"
        p <- tryCatch(
          mantelhaen.test(.arr, alternative = "greater", exact = FALSE)$p.value,
          error = function(e2) {
            print(paste("mantelhaen test failed at:", .BY$reporter, .BY$gene_name))
            NA_real_
          })
        if (is.na(p)) method <- NA_character_
      }
    }
    
    list(
      n_strata         = K,
      max_GFPlow_count = max(.arr[1, 1, ]),
      max_Input_count  = max(.arr[1, 2, ]),
      logOR            = est$logOR,
      SE_logOR         = est$SE,
      p_value          = as.numeric(p),
      method           = method
    )
  }
  
}, .(reporter, gene_id, gene_name)]

cmh_enr1[,padj := p.adjust(p_value, method = "fdr"), reporter]
cmh_enr1[, log2OR := logOR/log(2)]
cmh_enr1[, log2SE := SE_logOR/log(2)]
cmh_enr1[, hit := padj < PARAMS$padj_cutoff & log2OR > PARAMS$l2OR_cutoff]
fwrite(cmh_enr1, "db/FC_tables/AIDscreen_CMH_enrichment.csv")

#=======================#
# CMH, compare Promoters, Wald Test ----
#=======================#
wald_reporters <- dcast(cmh_enr1[is.finite(log2OR) & (max_GFPlow_count >= PARAMS$gene_count_cutoff | max_Input_count >=PARAMS$gene_count_cutoff)],
                        gene_id + gene_name ~ reporter,
                        value.var = c("log2OR", "log2SE", "padj")) %>%
  drop_na() %>%
  mutate(`log2OR_EF1A-Atp5k` = log2OR_EF1a -log2OR_Atp5k,
         SE_delta = sqrt(log2SE_EF1a^2 + log2SE_Atp5k^2),
         z = `log2OR_EF1A-Atp5k` / SE_delta,
         p_wald = 2 * pnorm(-abs(z)),
         .by = c(gene_id, gene_name)
  )
setDT(wald_reporters)

# FDR correction across genes
wald_reporters[gene_name %in% cmh_enr1[(hit), gene_name] & !is.na(p_wald), padj_wald := p.adjust(p_wald, method = "BH")]
fwrite(wald_reporters, "db/FC_tables/AIDscreen_EF1a_vs_Atp5k_WaldTest.csv")

# Gm42715 = Srcap
# Gm49804 = Nudt3
