#---- Libraries ------------
library(data.table)
library(tidyverse)
source("src/functions/plotting_fun.R")
source("src/functions/ggplot2_theme.R")

#---- Settings -----
PARAMS <- list(
  dc_cutoff = 10,
  gene_count_cutoff = 3,
  padj_cutoff = 0.01,
  padj_eva = 0.1,
  log2OR_cutoff = 0,
  log2OR_cutoff_eva = 1,
  pseudocount = 1
)

#---- Parse data ------------

# AID screen data
int <- fread("Rdata/AIDscreen_Integrations.txt.gz")
int[, c("reporter","selection_type", "frame", "fraction", "rep") := tstrsplit(sampleID, "_" )]
setorder(int, -InFrame)
# Xon data
Xon.dat <- fread("/groups/stark/vojtek/FACS/Xon_analysis/Combined_data.csv")
Xon.dat[, rep := seq_len(.N), by = .(Sample, Drug, `Time[h]`)]

#---- Supplementary Figure 1C: Classical vs intein approach ------------

# Background as reverse strand intronic
Classic_counts <- int[grepl("Classic", sampleID) & dc >= PARAMS$dc_cutoff] %>%
  dplyr::mutate(InFrame = (!Sense & element == "intron")) %>%
  dplyr::arrange(-InFrame) %>%
  distinct(sampleID, idx, .keep_all = T) %>%
  dplyr::count(InFrame, sampleID, phase) %>% 
  mutate(prop = n/sum(n),
         selection = "Classic",
         frame = as.character(phase + 1),
         .by = sampleID) %>%
  dplyr::filter(InFrame) %>%
  dplyr::select(-phase)
 
Intein_counts <- int %>%
  dplyr::filter(dc >= PARAMS$dc_cutoff, grepl("test_Intein", sampleID)) %>% 
  distinct(sampleID, idx, .keep_all = T) %>%
  dplyr::count(InFrame, sampleID) %>% 
  mutate(prop = n/sum(n),
         selection = "SIFT",
         frame = gsub(".*(Intein|Classic)_([1-3]).*", "\\2", sampleID),
         .by = sampleID) %>%
  dplyr::filter(InFrame)
  
FigS1B <- rbind(Intein_counts, Classic_counts) %>%
  { 
    df <- .
    
    # Summary: mean prop per frame x selection
    summary_df <- df %>%
      summarise(Mean = mean(prop),
                Min = min(prop) -0.05,
                Max = max(prop)+0.05, .by = c(frame, selection))
    
    # Arrow + fold change data: one row per frame
    arrow_df <- summary_df %>%
      tidyr::pivot_wider(names_from = selection, values_from = c(Mean, Min, Max)) %>%
      mutate(fc = Mean_SIFT/ Mean_Classic,
             y = Max_Classic,
             y_end = Min_SIFT,
             fc_label = paste0(round(fc, 1), "x"))
    
    ggplot(df, aes(x = frame, y = prop, color = selection)) +
      geom_point(alpha = 0.7, size = 2, shape = 16,
                 position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0.01, dodge.width = 0)) +
      
      # Double-sided arrow between the two group means
      geom_segment(data = arrow_df,
                   aes(x = as.numeric(frame) ,
                       xend = as.numeric(frame),
                       y = y, yend = y_end),
                   color = "gray67",
                   arrow = arrow(ends = "both", length = unit(6, "pt"), type = "closed"),
                   inherit.aes = FALSE) +
      # Fold change label
      geom_text(data = arrow_df,
                aes(x = as.numeric(frame) + 0.2,
                    y = (Mean_Classic + Mean_SIFT) / 2,
                    label = fc_label),
                color = "black", size = 7, size.unit = "pt",
                inherit.aes = FALSE)
  } +
  theme(legend.position = "top") +
  scale_color_manual(values = c("gray", "brown")) +
  labs(y = "Correct integrations / all",
       x = "Frame")
ggsave("results/FigS1B.pdf",FigS1B, width = 50, height =40, units = "mm")

#----- Supl. Fig. 1F: Xon vs dGFP -----
FigS1F <- Xon.dat %>%
  dplyr::filter(cell_line %in% c("13", "12"), Drug %in% c("LMI070", "Triptolide")) %>%
  dplyr::summarise(GFPpos = mean(GFPpos), .by = c(label, cell_line,clone, Drug,`Time[h]`)) %>%
  mutate(label = gsub("EF1a-", "", label),
         GFPpos = case_when(
           !grepl("Xon", label) ~ 100-GFPpos,
           .default = GFPpos
         )) %>%
  ggplot() +
  aes(x = `Time[h]`, y = GFPpos, color = label, group = label, fill = label) +
  geom_point(size =1, alpha = 0.6, shape = 16,show.legend = F) +  
  stat_summary(
    geom = "smooth",
    fun.data = mean_cl_normal,
    fun.args = list(conf.int = 0.95),
    alpha = .2,
    size = 0.5,
    se = F
  ) +
  labs(y = "% conversion",
       x ="Time [h]",
       fill = NULL,
       color = NULL) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.25) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray60", linewidth = 0.25) +
  scale_color_manual(values = c("brown", "black")) +
  scale_fill_manual(values =  c("brown", "black")) +
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.83, 0.23))


#----- Supl. Fig. 1G: Dynamic range -----
FigS1G_dat <- Xon.dat %>%
  dplyr::filter(cell_line %in% c("13", "12"), Drug %in% c("LMI070", "Triptolide")) %>%
  dplyr::summarise(GFP_intensity = mean(GFP_intensity),
                                .by = c(label, Sample, cell_line,clone, Drug,`Time[h]`)) %>%
  summarise(DR = log10(max(GFP_intensity)/min(GFP_intensity)),
         .by = label) %>%
  mutate(label = gsub("EF1a-", "", label))

DR_FC <- 10^(FigS1G_dat %>% dplyr::filter(label == "Xon-EGFP") %>% pull(DR) %>% mean() -
               FigS1G_dat %>% dplyr::filter(label == "mNG(d2)") %>% pull(DR) %>% mean())

FigS1G <-   FigS1G_dat %>%
  ggplot(aes(y = label, x = DR)) +
  geom_jitter(alpha = 0.6, size = 2, height = 0.1, width = 0, shape = 16) +
  labs(y = NULL,
       x = "Dynamic range (log10)") +
  annotate("text",
           x = 2,
           y = 1.5,  # position between the two labels, adjust as needed
           label = paste0(round(DR_FC, 1), "x"),
           size = 2, fontface = "bold", color = "brown") +
  annotate("segment",
           x = 1.5, xend = 2.388,
           y = 2, yend = 2,  # adjust to match your label positions
           arrow = arrow(ends = "both", length = unit(0.2, "cm"), type = "closed"), color = "brown")

cowplot::plot_grid(FigS1F, FigS1G, nrow = 2, align = "hv", rel_heights = c(1,0.7))
ggsave("results/FigS1F,G.pdf", width = 50, height = 55, units = "mm")
