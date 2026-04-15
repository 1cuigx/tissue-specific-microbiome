# Fig. 1b: Stacked bar plot of phylum-level 16S rRNA composition.
# Input : Fig1b_phyla.xlsx
# Output: Fig1b_phyla.pdf
# Run with this folder as the working directory.

library(tidyverse)
library(readxl)

phyla16s <- read_xlsx("Fig1b_phyla.xlsx")
phyla16s <- phyla16s[1:16, ]
phyla16s$group <- rep(c("ApoEpi", "ApoGas", "SymEpi", "SymGas"), each = 4)

data_long <- phyla16s %>%
  pivot_longer(cols = -c(sample, group), names_to = "phylum", values_to = "percentage")

color <- c(
  "Actinobacteria"  = "#5e9f8e",
  "Bacteroidetes"   = "#366c3f",
  "Planctomycetes"  = "#85568c",
  "Verrucomicrobia" = "#633e27",
  "Firmicutes"      = "#d2c67e",
  "Proteobacteria"  = "#b2656f",
  "Other_Phyla"     = "#90908d"
)

ggplot(data_long, aes(x = sample, y = percentage, fill = phylum)) +
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

ggsave("Fig1b_phyla.pdf", width = 8, height = 4, units = "in")
