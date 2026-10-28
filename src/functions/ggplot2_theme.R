require(ggplot2)
theme_set(
  theme_classic() +
    theme(
      axis.text        = element_text(size = 6, color = "black"),
      axis.title       = element_text(size = 7, color = "black"),
      legend.text      = element_text(size = 6, color = "black"),
      legend.title     = element_text(size = 7, color = "black"),
      strip.text       = element_text(size = 7, color = "black"),
      plot.title       = element_text(size = 7, color = "black"),
      axis.line        = element_line(linewidth = 0.25),
      axis.ticks       = element_line(linewidth = 0.25),
      panel.grid       = element_blank(),
      text             = element_text(family = "Helvetica"),
      
      legend.margin      = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0, "mm"),
      legend.box.margin  = margin(4, 4, 4, 4),   # margin around the entire legend box
      legend.key.size      = unit(2, "mm"),   # size of the legend key/symbol
      legend.key.spacing.y = unit(0, "mm"),   # vertical space between items
      legend.key.spacing.x = unit(0, "mm"),    # horizontal space between items
      
      plot.background  = element_blank(),
      panel.background = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key       = element_blank(),
      strip.background = element_blank()
    )
)


theme_publish <- function(){
    theme_classic() +
      theme(
        axis.text        = element_text(size = 6, color = "black"),
        axis.title       = element_text(size = 7, color = "black"),
        legend.text      = element_text(size = 6, color = "black"),
        legend.title     = element_text(size = 7, color = "black"),
        strip.text       = element_text(size = 7, color = "black"),
        plot.title       = element_text(size = 7, color = "black"),
        axis.line        = element_line(linewidth = 0.25),
        axis.ticks       = element_line(linewidth = 0.25),
        panel.grid       = element_blank(),
        text             = element_text(family = "Helvetica"),
        
        plot.background  = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key       = element_blank(),
        strip.background = element_blank()
      )
}

# Wong palette
wong_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

# Set as default ggplot2 color/fill scales
options(
  ggplot2.discrete.colour = wong_colors,
  ggplot2.discrete.fill   = wong_colors
)