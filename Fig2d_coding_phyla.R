# Fig. 2d: Stacked bar plot of phylum composition from coding mRNA data,
# with two-way ANOVA, Tukey HSD post-hoc and non-parametric analyses.
#
# Input : Fig2d_coding_phyla.xlsx
# Output: Fig2d_coding_phyla.pdf
#         Supplementary_Table_assumption_checks_mRNA_Phyla.csv
#         Supplementary_Table_ANOVA_FDR_mRNA_Phyla.csv
#         Supplementary_Table_TukeyHSD_mRNA_Phyla.csv
#         Supplementary_Table_NonParametric_mRNA_Phyla.csv
# Run with this folder (DataCode) as the working directory.

library(tidyverse)
library(readxl)
library(car)
library(emmeans)
library(multcomp)
library(rstatix)

phyla <- read_xlsx("Fig2d_coding_phyla.xlsx")
colnames(phyla)[7] <- "TOther_Phyla"

data_long <- phyla %>%
  pivot_longer(cols = -c(sample, group), names_to = "phylum", values_to = "percentage") %>%
  mutate(
    symbiosis = ifelse(grepl("Sym|sym", group), "Symbiotic", "Aposymbiotic"),
    tissue    = ifelse(grepl("Epi|epi", group), "Epidermis", "Gastrodermis")
  )

# ---- 1. ANOVA assumption checks ----
assumption_checks <- data_long %>%
  group_by(phylum) %>%
  summarize(
    shapiro_p = {
      mod <- aov(percentage ~ tissue * symbiosis, data = cur_data())
      shapiro.test(residuals(mod))$p.value
    },
    levene_p = {
      car::leveneTest(percentage ~ interaction(tissue, symbiosis),
                      data = cur_data())$`Pr(>F)`[1]
    },
    .groups = "drop"
  )
write.csv(assumption_checks, "Supplementary_Table_assumption_checks_mRNA_Phyla.csv",
          row.names = FALSE)

# ---- 2. Two-way ANOVA + FDR correction ----
full_anova_results <- data_long %>%
  group_by(phylum) %>%
  reframe({
    mod <- aov(percentage ~ tissue * symbiosis, data = cur_data())
    at  <- summary(mod)[[1]]
    tibble(
      Source       = rownames(at),
      Df           = at$Df,
      Sum_Sq       = round(at$`Sum Sq`, 4),
      Mean_Sq      = round(at$`Mean Sq`, 4),
      F_value      = round(at$`F value`, 3),
      p_value      = at$`Pr(>F)`,
      Variance_pct = round(at$`Sum Sq` / sum(at$`Sum Sq`) * 100, 1)
    )
  })

pval_summary <- full_anova_results %>%
  filter(!is.na(p_value)) %>%
  group_by(Source) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write.csv(pval_summary, "Supplementary_Table_ANOVA_FDR_mRNA_Phyla.csv",
          row.names = FALSE)

# ---- 3. Tukey HSD post-hoc ----
posthoc_results <- data_long %>%
  group_by(phylum) %>%
  reframe({
    mod <- aov(percentage ~ group, data = cur_data())
    overall_p <- summary(mod)[[1]]$`Pr(>F)`[1]
    if (!is.na(overall_p) && overall_p < 0.05) {
      tukey <- TukeyHSD(mod, "group")$group
      tibble(Comparison = rownames(tukey),
             Difference = round(tukey[, "diff"], 4),
             Lower_CI = round(tukey[, "lwr"], 4),
             Upper_CI = round(tukey[, "upr"], 4),
             p_adjusted = round(tukey[, "p adj"], 6),
             Significant = tukey[, "p adj"] < 0.05)
    } else {
      tibble(Comparison = "ANOVA not significant",
             Difference = NA_real_, Lower_CI = NA_real_,
             Upper_CI = NA_real_, p_adjusted = NA_real_, Significant = NA)
    }
  })
write.csv(posthoc_results, "Supplementary_Table_TukeyHSD_mRNA_Phyla.csv",
          row.names = FALSE)

# ---- 4. Non-parametric alternative ----
nonparam_results <- data_long %>%
  group_by(phylum) %>%
  reframe({
    kt <- kruskal.test(percentage ~ group, data = cur_data())
    if (kt$p.value < 0.05) {
      dunn <- rstatix::dunn_test(cur_data(), percentage ~ group,
                                 p.adjust.method = "BH")
      tibble(Test = "Kruskal-Wallis + Dunn",
             KW_chi2 = round(kt$statistic, 3),
             KW_df = kt$parameter, KW_p = kt$p.value,
             Comparison = dunn$.y.,
             Group1 = dunn$group1, Group2 = dunn$group2,
             Dunn_p_adj = dunn$p.adj, Dunn_signif = dunn$p.adj.signif)
    } else {
      tibble(Test = "Kruskal-Wallis",
             KW_chi2 = round(kt$statistic, 3),
             KW_df = kt$parameter, KW_p = kt$p.value,
             Comparison = NA_character_,
             Group1 = NA_character_, Group2 = NA_character_,
             Dunn_p_adj = NA_real_, Dunn_signif = NA_character_)
    }
  })
write.csv(nonparam_results, "Supplementary_Table_NonParametric_mRNA_Phyla.csv",
          row.names = FALSE)

# ---- 5. Stacked bar plot ----
color <- c(
  "Actinobacteria" = "#5e9f8e",
  "Bacteroidetes"  = "#366c3f",
  "Cyanobacteria"  = "#909247",
  "Firmicutes"     = "#d2c67e",
  "Proteobacteria" = "#b2656f",
  "TOther_Phyla"   = "#cbcccb"
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

ggsave("Fig2d_coding_phyla.pdf", width = 8, height = 4, units = "in")
