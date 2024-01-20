# ------------------------------------------------
# Utils
# Chris Hoover, Jan 2022
# ------------------------------------------------

library(ggplot2)

# Plotting Utils --------------------------------
scale_color_hpiq <- function(){
  scale_color_manual(values = c("royalblue", "lightblue","lightgreen","darkgreen"),
                     breaks = c(1:4))
}

scale_fill_hpiq <- function(){
  scale_fill_manual(values = c("royalblue", "lightblue","lightgreen","darkgreen"),
                    breaks = c(1:4))
}


theme_here <- function(){
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
}
