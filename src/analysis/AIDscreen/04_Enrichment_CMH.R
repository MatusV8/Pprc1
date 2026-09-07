#================================#
# Setup ----
#================================#
# Libraries
library(data.table)
library(tidyverse)
source("/groups/stark/vojtek/matools/mv_OddsRatio.R")
source("src/functions/Fisher_grouped.R")
source("src/functions/Rowwise_fisher_test.R")
source("/groups/stark/vojtek/matools/mv_hm.R")
library(ComplexHeatmap)

PARAMS <- list(
  dc_cutoff = 10,
  gene_count_cutoff = 1,
  padj_cutoff = 0.05,
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

cmh_enr1 <- cmh_enr[, {
  K <- .N
  
  if (K == 0) {
    list(
      n_strata = 0L,
      max_GFPlow_count =NA_real_,
      max_Input_count = NA_real_,
      OR      = NA_real_,
      CI_low  = NA_real_,
      CI_high = NA_real_,
      p_value = NA_real_,
      SE_logOR = NA_real_,
      method = NA_character_
    )
  } else if (K == 1) {
    # Single stratum — use Fisher's exact test instead
    mat <-  matrix(c(GFPlow_count,
                     GFPlow_total - GFPlow_count,
                     Input_count, Input_total  - Input_count), nrow = 2)
    
    tryCatch({
      test     <- fisher.test(mat, alternative = "greater")
      test_two <- fisher.test(mat, alternative = "two.sided")
      
      # If estimate is Inf then use a pseudocount to estimate
      #if(is.infinite(test$estimate)){
      #  test$estimate     <- fisher.test(mat + 1, alternative = "greater")$estimate
      #  test_two$estimate          <- fisher.test(mat + 1, alternative = "two.sided")$estimate
      #}
      
      # SE of logOR back-calculated from two-sided 95% CI
      se_logOR <- (log(test_two$conf.int[2]) - log(test_two$conf.int[1])) / (2 * 1.96)
      
      list(
        n_strata = 1L,
        max_GFPlow_count = as.numeric(mat[1,1]),
        max_Input_count = as.numeric(mat[1,2]),
        OR      = as.numeric(test_two$estimate),
        CI_low  = as.numeric(test_two$conf.int[1]),
        CI_high = as.numeric(test_two$conf.int[2]),
        p_value = as.numeric(test$p.value),
        SE_logOR      = as.numeric(se_logOR),
        method  = "Fisher's_exact"
      )
    }, error = function(e) {
      print(paste("Fisher's exact test failed at:",.BY$group, .BY$gene_name, sep = " "))
      list(
        n_strata = 1L,
        max_GFPlow_count = mat[1,1],
        max_Input_count = mat[1,2],
        OR      = NA_real_,
        CI_low  = NA_real_,
        CI_high = NA_real_,
        p_value = NA_real_,
        SE_logOR = NA_real_,
        method = NA_character_
      )
    })
  } else {
    
    # Construct array
    .arr <- array(NA, dim = c(2, 2, K),
                  dimnames = list(
                    c("in_gene", "not_in_gene"),
                    c("GFPlow", "Input")
                  ))
    for (k in 1:K) {
      a = as.numeric(GFPlow_count[k])
      b = as.numeric(Input_count[k])
      c = as.numeric(GFPlow_total[k] - GFPlow_count[k])
      d = as.numeric(Input_total[k]  - Input_count[k])
      
      m <- matrix(c(a, c,
                    b, d), nrow = 2)
      
      .arr[, , k] <- m
    }
    
    
    # Perform test
    tryCatch({
      test     <- mantelhaen.test(.arr, alternative = "greater",   exact = TRUE)
      test_two <- mantelhaen.test(.arr, alternative = "two.sided", exact = TRUE)
      
      # SE of logOR back-calculated from two-sided 95% CI
      se_logOR <- (log(test_two$conf.int[2]) - log(test_two$conf.int[1])) / (2 * 1.96)
      
      list(
        n_strata = K,
        max_GFPlow_count = max(.arr[1,1,]),
        max_Input_count = max(.arr[1,2,]),
        OR      = as.numeric(unname(test_two$estimate)),
        CI_low  = as.numeric(test_two$conf.int[1]),
        CI_high = as.numeric(test_two$conf.int[2]),
        p_value = as.numeric(test$p.value),
        SE_logOR      = as.numeric(se_logOR),
        method  = "CMH_exact"
      )
    }, error = function(e) {
      # exact = TRUE failed — fall back to asymptotic
      tryCatch({
        test     <- mantelhaen.test(.arr, alternative = "greater",   exact = FALSE)
        test_two <- mantelhaen.test(.arr, alternative = "two.sided", exact = FALSE)
        
        # SE of logOR back-calculated from two-sided 95% CI
        se_logOR <- (log(test_two$conf.int[2]) - log(test_two$conf.int[1])) / (2 * 1.96)
        
        list(
          n_strata = K,
          max_GFPlow_count = max(.arr[1,1,]),
          max_Input_count = max(.arr[1,2,]),
          OR      = as.numeric(unname(test_two$estimate)),
          CI_low  = as.numeric(test_two$conf.int[1]),
          CI_high = as.numeric(test_two$conf.int[2]),
          p_value = as.numeric(test$p.value),
          SE_logOR      = as.numeric(se_logOR),
          method  = "CMH_aprox"
        )
      }, error = function(e2) {
        print(paste("mantelhaen test failed at:",.BY$reporter, .BY$gene_name, sep = " "))
        list(
          n_strata = K,
          max_GFPlow_count = max(.arr[1,1,]),
          max_Input_count = max(.arr[1,2,]),
          OR      = NA_real_,
          CI_low  = NA_real_,
          CI_high = NA_real_,
          p_value = NA_real_,
          SE_logOR = NA_real_,
          method  = "CMH_aprox"
        )
      })
    })
  }
  
}, .(reporter, gene_id, gene_name)]

cmh_enr1[,padj := p.adjust(p_value, method = "fdr"), reporter]
cmh_enr1[, log2OR := log2(OR)]
cmh_enr1[, log2SE := SE_logOR/log(2)]
cmh_enr1[, hit := padj < 0.01 & log2OR >0]
fwrite(cmh_enr1, "db/FC_tables/AIDscreen_CMH_enrichment.csv")

## Use pseudocount of 1 to all
cmh_enr_pc1 <- cmh_enr[, {
  K <- .N
  
  if (K == 0) {
    list(OR_mh = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_, SE_logOR = NA_real_)
  } else if (K == 1) {
    # Single stratum — use Fisher's exact test instead
    mat <-  matrix(c(GFPlow_count,
                     GFPlow_total - GFPlow_count,
                     Input_count, Input_total  - Input_count), nrow = 2) +1
    
    tryCatch({
      test     <- fisher.test(mat, alternative = "greater")
      test_two <- fisher.test(mat, alternative = "two.sided")
      
      # If estimate is Inf then use a pseudocount to estimate
      if(is.infinite(test$estimate)){
        test$estimate     <- fisher.test(mat + 1, alternative = "greater")$estimate
        test_two$estimate          <- fisher.test(mat + 1, alternative = "two.sided")$estimate
      }
      
      # SE of logOR back-calculated from two-sided 95% CI
      se_logOR <- (log(test_two$conf.int[2]) - log(test_two$conf.int[1])) / (2 * 1.96)
      
      list(
        n_strata = 1L,
        max_GFPlow_count = mat[1,1],
        max_Input_count = mat[1,2],
        OR      = as.numeric(test_two$estimate),
        CI_low  = as.numeric(test_two$conf.int[1]),
        CI_high = as.numeric(test_two$conf.int[2]),
        p_value = as.numeric(test$p.value),
        SE_logOR      = as.numeric(se_logOR),
        method  = "Fisher's_exact"
      )
    }, error = function(e) {
      print(paste("Fisher's exact test failed at:",.BY$group, .BY$gene_name, sep = " "))
      list(
        n_strata = 1L,
        max_GFPlow_count = mat[1,1],
        max_Input_count = mat[1,2],
        OR      = NA_real_,
        CI_low  = NA_real_,
        CI_high = NA_real_,
        p_value = NA_real_,
        SE_logOR = NA_real_,
        method = NA_real_
      )
    })
  } else {
    
    # Construct array
    .arr <- array(NA, dim = c(2, 2, K),
                  dimnames = list(
                    c("in_gene", "not_in_gene"),
                    c("GFPlow", "Input")
                  ))
    for (k in 1:K) {
      a = as.numeric(GFPlow_count[k])
      b = as.numeric(Input_count[k])
      c = as.numeric(GFPlow_total[k] - GFPlow_count[k])
      d = as.numeric(Input_total[k]  - Input_count[k])
      
      m <- matrix(c(a, c,
                    b, d), nrow = 2) +1
      
      .arr[, , k] <- m
    }
    
    
    # Perform test
    tryCatch({
      test     <- mantelhaen.test(.arr, alternative = "greater",   exact = TRUE)
      test_two <- mantelhaen.test(.arr, alternative = "two.sided", exact = TRUE)
      
      # SE of logOR back-calculated from two-sided 95% CI
      se_logOR <- (log(test_two$conf.int[2]) - log(test_two$conf.int[1])) / (2 * 1.96)
      
      list(
        n_strata = K,
        max_GFPlow_count = max(.arr[1,1,]),
        max_Input_count = max(.arr[1,2,]),
        OR      = as.numeric(unname(test_two$estimate)),
        CI_low  = as.numeric(test_two$conf.int[1]),
        CI_high = as.numeric(test_two$conf.int[2]),
        p_value = as.numeric(test$p.value),
        SE_logOR      = as.numeric(se_logOR),
        method  = "CMH_exact"
      )
    }, error = function(e) {
      # exact = TRUE failed — fall back to asymptotic
      tryCatch({
        test     <- mantelhaen.test(.arr, alternative = "greater",   exact = FALSE)
        test_two <- mantelhaen.test(.arr, alternative = "two.sided", exact = FALSE)
        
        # SE of logOR back-calculated from two-sided 95% CI
        se_logOR <- (log(test_two$conf.int[2]) - log(test_two$conf.int[1])) / (2 * 1.96)
        
        list(
          n_strata = K,
          max_GFPlow_count = max(.arr[1,1,]),
          max_Input_count = max(.arr[1,2,]),
          OR      = as.numeric(unname(test_two$estimate)),
          CI_low  = as.numeric(test_two$conf.int[1]),
          CI_high = as.numeric(test_two$conf.int[2]),
          p_value = as.numeric(test$p.value),
          SE_logOR      = as.numeric(se_logOR),
          method  = "CMH_aprox"
        )
      }, error = function(e2) {
        print(paste("mantelhaen test failed at:",.BY$reporter, .BY$gene_name, sep = " "))
        list(
          n_strata = K,
          max_GFPlow_count = max(.arr[1,1,]),
          max_Input_count = max(.arr[1,2,]),
          OR      = NA_real_,
          CI_low  = NA_real_,
          CI_high = NA_real_,
          p_value = NA_real_,
          SE_logOR = NA_real_,
          method  = "CMH_aprox"
        )
      })
    })
  }
  
}, .(reporter, gene_id, gene_name)]
cmh_enr_pc1[,padj := p.adjust(p_value, method = "fdr"), reporter]
cmh_enr_pc1[, log2OR := log2(OR)]
cmh_enr_pc1[, log2SE := SE_logOR/log(2)]
cmh_enr_pc1[, hit := padj < 0.05 & log2OR >0]

#=======================#
# CMH, compare Promoters, Wald Test ----
#=======================#
wald_reporters <- dcast(cmh_enr1[is.finite(log2OR) & max_GFPlow_count >= PARAMS$gene_count_cutoff | max_Input_count >=PARAMS$gene_count_cutoff], gene_id + gene_name ~ reporter, value.var = c("log2OR", "log2SE")) %>%
  drop_na() %>%
  mutate(EF1a_Atp5k_log2OR_delta = log2OR_EF1a -log2OR_Atp5k,
         SE_delta = sqrt(log2SE_EF1a^2 + log2SE_Atp5k^2),
         z = EF1a_Atp5k_log2OR_delta / SE_delta,
         p_wald = 2 * pnorm(-abs(z)),
         .by = c(gene_id, gene_name)
  )
setDT(wald_reporters)

# FDR correction across genes
wald_reporters[gene_name %in% cmh_enr1[padj < 0.01, gene_name] & !is.na(p_wald), p_wald_fdr := p.adjust(p_wald, method = "BH")]

wald_reporters[!is.na(p_wald_fdr)] %>%
  ggplot(aes(x = log2OR_Atp5k, y = log2OR_EF1a, label = gene_name)) +
  theme_classic() +
  geom_point(alpha = 0.5, shape = 16) +
  geom_abline(intercept = c(-1, 1), linetype = "dashed") +
  ggrepel::geom_text_repel(data = . %>% dplyr::slice_max(order_by = abs(EF1a_Atp5k_log2OR_delta), n = 10))
fwrite(wald_reporters, "db/FC_tables/AIDscreen_EF1a_vs_Atp5k_WaldTest.csv")
       
#----#
# GO analysis
#---#
## GOs
source("/groups/stark/vojtek/matools/GO_functions.R")
library(GO.db)
GOs <- getGOs()

GO_bg <- GO_bg(unique(cmh_enr1$gene_id), GOs)

setDT(GOs)

# Continue here
# get_genes_with_descendants()

# Pol II transcription
Tx_terms <- unique(unlist(mget(c("GO:0006366",
                                 "GO:0045893",
                                 "GO:0006357" #regulation of transcription by RNA polymerase II
), GOBPOFFSPRING)))
Tx_terms <- c("GO:0006366", "GO:0045893",Tx_terms)
Tx_neg_terms <- unique(unlist(mget(c("GO:0000122",
                                     "GO:0045892",
                                     "GO:0045943", # Pol I
                                     "GO:0045945", # Pol III
                                     "GO:1904281" # Pol V
), GOBPOFFSPRING))) # Exclude Pol I and Pol III
Tx_neg_terms <- c("GO:0000122",
                  "GO:0045892",
                  "GO:0045943", # Pol I
                  "GO:0045945", # Pol III
                  "GO:1904281", Tx_neg_terms)
Tx_terms <- setdiff(Tx_terms, Tx_neg_terms)

# translation
Translation_terms <- unique(unlist(mget(c("GO:0006412", "GO:0006417"), GOBPOFFSPRING)))
Translation_terms <- c("GO:0006412", "GO:0006417", Translation_terms)
Translation_neg_terms <- unique(unlist(mget("GO:0017148", GOBPOFFSPRING)))
Translation_neg_terms <- c("GO:0017148")
Translation_terms <- setdiff(Translation_terms, Translation_neg_terms)
Translation_terms <- setdiff(Translation_terms, "GO:1903912")

# mRNA processing
mRNA_process_terms <- unique(unlist(mget(c("GO:0006397", "GO:0048024"), GOBPOFFSPRING)))
mRNA_process_terms <- c("GO:0006397", "GO:0048024", mRNA_process_terms)
mRNA_process_neg_terms <- unique(unlist(mget(c("GO:0050686", "GO:0048025"), GOBPOFFSPRING)))
mRNA_process_neg_terms <- c("GO:0050686", "GO:0048025", mRNA_process_neg_terms)
mRNA_process_terms <- setdiff(mRNA_process_terms, mRNA_process_neg_terms)

# mRNA export
mRNA_export_terms <- unique(unlist(mget(c("GO:0006406", "GO:0046831"), GOBPOFFSPRING)))
mRNA_export_terms <- c("GO:0006406", "GO:0046831", mRNA_export_terms)
mRNA_export_neg_terms <- unique(unlist(mget(c("GO:0046832"), GOBPOFFSPRING)))
mRNA_export_terms <- setdiff(mRNA_export_terms, mRNA_export_neg_terms)

GO_category_ids <- list(
  "transcription" = unique(c(
    GOs[ID %in% Tx_terms, gene_id], # Positive regulation of PolII
    GOs[ID %in% unique(unlist(mget("GO:0016591", GOCCOFFSPRING))) , gene_id], # RNPII and GTFs
    GOs[ID %in% "GO:0071339", gene_id], # MLL1
    "ENSMUSG00000033323" # Ctdp1
  )),
  
  "mRNA processing" = unique(c(
    GOs[ID %in% mRNA_export_terms, gene_id],
    GOs[ID %in% mRNA_process_terms, gene_id]
    #GOs[ID %in% as.character(GOBPOFFSPRING[["GO:0006397"]]) & grepl("positive", description), unique(ID)], # positive regulation of mRNA processing
    #GOs[ID == "GO:0005681", unique(gene_id)], # Spliceosome components
    #GOs[ID == "GO:0006406", unique(gene_id)] # mRNA export from nucleus
  )),
  
  "mRNA cap" = unique(c(
    GOs[ID == "GO:0006370", gene_id], # 7-methylguanosine mRNA capping
    GOs[ID == "GO:0000339", gene_id], # RNA cap binding 
    GOs[ID == "GO:0098808", gene_id], # mRNA Cap binding
    "ENSMUSG00000055491" # Pprc1
  )),
  
  "translation" = unique(
    GOs[ID %in% Translation_terms, gene_id]
    #GOs[ID %in% as.character(GOBPOFFSPRING[["GO:0006412"]]) & grepl("positive", description), unique(gene_id)], # Positive regulation of translation
    #GOs[ID %in% as.character(GOBPOFFSPRING[["GO:0006418"]]), unique(gene_id)], # tRNA aminoacylation
    #GOs[ID == "GO:0003735", unique(gene_id)], # Ribosome
    #GOs[grepl("eukaryotic translation initiation factor", description), unique(gene_id)], # translation initiation factors
    #GOs[ID == "GO:0003746", unique(gene_id)], # translation elongation
    #GOs[ID == "GO:0008079", unique(gene_id)] # translation termination
    #"cytoplasmic translation" = "GO:0002181",
    #"positive regulation of translation" = "GO:0045727",
    #"translation factor activity, RNA binding" = "GO:0008135",
    #"translational elongation" = "GO:0006414",
    #"translational initiation" = "GO:0006413",
    #"translational termination" = "GO:0006415",
    #"triplet codon-amino acid adaptor activity" = "GO:0030533",
    #"tRNA aminoacylation for protein translation" = "GO:0006418",
    #"aminoacyl-tRNA synthetase multienzyme complex" = "GO:0017101"
  )
)

GO_anno <- as.data.table(stack(GO_category_ids))[, .(GO_category = ind, gene_id = values)]
#GO_anno <- data.table(
#  GO_category = names(GO_category_ids),
#  gene_id = map(GO_category_ids, ~{
#    
#    map(.x, function(GO_ID){
#      
#      # Get gene id's including for the offspring
#      get_genes_with_descendants(GO_ID, GOs)
#    }) %>%
#      unlist() %>%
#      unique()
#  })
#)[, .(gene_id = unlist(gene_id)), by = GO_category]

# Set the hierarchy of annotations and remove duplicates
GO_anno[, GO_category := factor(GO_category, levels = rev(c("translation","mRNA processing", "transcription","mRNA cap" )))]
GO_anno <- GO_anno[order(GO_category), .SD[1], by  = gene_id]
#GO_anno[gene_id == "ENSMUSG00000055491", GO_category := "mRNA cap"] # Pprc1 as mRNA cap
GO_anno[gene_id == "ENSMUSG00000016554", GO_category := "translation"] # Eif3d as translation

Ensbl_id_to_names <- function(ensbl_ids) {
  # Connect to the Ensembl Mart for mouse genes (use "hsapiens_gene_ensembl" for human)
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  
  # Query Ensembl for gene symbols based on Ensembl IDs
  gene_map <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = ensbl_ids,
    mart = mart
  )
  
  # Return matched gene names; unmapped entries will become NA
  mapped <- gene_map$external_gene_name[match(ensbl_ids, gene_map$ensembl_gene_id)]
  return(mapped)
}
#GO_anno[, gene_name := Ensbl_id_to_names(gene_id)]



# Test enrichment of custom GOs in the anno table
Test_custom_GOs <- GO_anno[, {
  
  GO_enrichement_test(gene_id, GO_bg, group_by_cols = c("type", "description"))
  
},GO_category]

# Specific hits


# Annotate FC tables
GO_enr <- GO_anno[cmh_enr1, on = "gene_id"]
GO_enr[is.na(GO_category), GO_category := "other"]
GO_enr <-distinct(GO_enr, gene_id, .keep_all = T)
GO_enr_EvA <- GO_anno[wald_reporters, on = "gene_id"]
GO_enr_EvA[is.na(GO_category), GO_category := "other"]

Atp5k_specific_hits <- wald_reporters[p_wald_fdr < 0.05 & EF1a_Atp5k_log2OR_delta <= -1, gene_id]
EF1a_specific_hits <- wald_reporters[p_wald_fdr < 0.05 & EF1a_Atp5k_log2OR_delta >= 1, gene_id]

Shared_hits <- unique(cmh_enr1[,.SD[all(hit)], gene_id][!(gene_id %in% c(Atp5k_specific_hits, EF1a_specific_hits)), gene_id])
Union_hits <- cmh_enr1[(hit), unique(gene_id)]

## Enrichment in promoter specific hits
enr_SpecificHits_selectGO <- map_dfr(unique(GO_enr$GO_category), ~{
  .ft <- fisher.test(GO_enr$GO_category == .x, GO_enr$gene_id %in% c(Atp5k_specific_hits, EF1a_specific_hits), alternative = "two.sided")
  
  data.table(
    hit_category = "specific",
    category = .x,
    estimate = .ft$estimate,
    p.value  = .ft$p.value)
})

enr_SharedHits_selectGO <- map_dfr(unique(GO_enr$GO_category), ~{
  .ft <- fisher.test(GO_enr$GO_category == .x, GO_enr$gene_id %in% Shared_hits, alternative = "two.sided")
  
  data.table(
    hit_category = "shared",
    category = .x,
    estimate = .ft$estimate,
    p.value  = .ft$p.value)
})

enr_Union_selectGO <- map_dfr(unique(GO_enr$GO_category), ~{
  .ft <- fisher.test(GO_enr$GO_category == .x, GO_enr$gene_id %in% Union_hits, alternative = "two.sided")
  
  data.table(
    hit_category = "shared",
    category = .x,
    estimate = .ft$estimate,
    p.value  = .ft$p.value)
})

GOenr_sharedVSspecific <- rbind(enr_SpecificHits_selectGO, enr_SharedHits_selectGO)[category != "other", padj := p.adjust(p.value, method = "fdr"), hit_category]

GOenr_unionVSspecific <- rbind(enr_SpecificHits_selectGO, enr_Union_selectGO)[category != "other", padj := p.adjust(p.value, method = "fdr"), hit_category]
# Plot
GOenr_sharedVSspecific %>%
  dplyr::filter(category != "other") %>%
  mutate(padj = ifelse(padj < 10^-5, 10^-5, padj)) %>%
  ggplot(aes(x = hit_category, y= category, color = padj < 0.001, size = log2(estimate))) +
  geom_point(shape = 16) +
  scale_color_manual(values = c("TRUE" = "tomato", "FALSE" = "gray")) +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "hit category", y = "gene category", color = "padj < 0.001")

GOenr_unionVSspecific %>%
  dplyr::filter(category != "other") %>%
  mutate(padj = ifelse(padj < 10^-5, 10^-5, padj)) %>%
  ggplot(aes(x = hit_category, y= category, color = padj < 0.001, size = log2(estimate))) +
  geom_point(shape = 16) +
  scale_color_manual(values = c("TRUE" = "tomato", "FALSE" = "gray")) +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "hit category", y = "gene category", color = "padj < 0.001")


# Plot mRNA cap related
GO_anno[cmh_enr1, on = "gene_id"][GO_category == "mRNA cap", .SD[sum(hit) > 0], gene_id] %>% 
  mutate(padj = ifelse(padj < 10^-3, 10^-3, padj)) %>%
  #mutate(cluster = tidy_kmeans(value = log2OR, category = reporter, id = gene_id, centers = 4, scale_cols = F,scale_rows = F)) %>%
  ggplot(aes(x = reporter, y = gene_name, color = log2OR, size = -log10(padj))) + 
  geom_point() + 
  scale_color_gradient2(low = "blue",mid =  "white", high = "red", midpoint = 0) +
  scale_size_continuous(
    breaks = c(0, 1, 2, 3),   # show only these values in legend
    range  = c(1, 3)           # min and max dot size in plot (in pts)
  ) +
  #scale_x_discrete(expand = c(0,0)) +
  #scale_y_discrete(expand = c(0,0)) +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Reporter", y = "Gene")
