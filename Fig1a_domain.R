# Fig. 1a: Stacked bar plot of domain-level 16S rRNA composition.
# Input : Fig1a_domain.xlsx
# Output: Fig1a_domain.pdf
# Run with the current folder as the working directory.

library(tidyverse)
library(readxl)

rRNA16s <- read_xlsx("Fig1a_domain.xlsx")
rRNA16s$group <- rep(c("ApoEpi", "ApoGas", "SymEpi", "SymGas"), each = 4)
colnames(rRNA16s)[1] <- "sample"

data_long <- rRNA16s %>%
  pivot_longer(cols = -c(sample, group), names_to = "taxa", values_to = "percentage")

color <- c("Archaea"   = "#3c2b67",
           "Bacteria"  = "#5d4c89",
           "Eukaryota" = "#9693b6",
           "unknown"   = "#efebf0")

ggplot(data_long, aes(x = sample, y = percentage, fill = taxa)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  facet_grid(~group, scales = "free_x") +
  scale_fill_manual(values = color) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text       = element_text(color = "black")
  ) +
  scale_y_continuous(expand = c(0, 0))

ggsave("Fig1a_domain.pdf", width = 8, height = 4, units = "in")
