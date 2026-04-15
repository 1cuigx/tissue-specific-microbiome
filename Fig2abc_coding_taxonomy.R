# Fig. 2a-c: mRNA-derived coding taxonomy boxplots with two-way ANOVA,
# Tukey HSD post-hoc tests and non-parametric alternatives.
#
# Input : Fig2abc_coding_taxonomy.xlsx
# Output: Fig2abc_coding_taxonomy.pdf
#         Supplementary_Table_assumption_checks_Fig2.csv
#         Supplementary_Table_ANOVA_Fig2.csv
#         Supplementary_Table_ANOVA_FDR_Fig2.csv
#         Supplementary_Table_TukeyHSD_Fig2.csv
#         Supplementary_Table_NonParametric_Fig2.csv
# Run with this folder (DataCode) as the working directory.

library(tidyverse)
library(readxl)
library(ggpubr)
library(car)
library(emmeans)
library(multcomp)
library(rstatix)

df <- read_xlsx("Fig2abc_coding_taxonomy.xlsx")

tissue_col <- c("ApoEpi" = "#A0CAE0",
                "ApoGas" = "#3082BC",
                "SymEpi" = "#FDAF69",
                "SymGas" = "#E55409")

# ---- Derive tissue/symbiosis factors from the four-level Group column ----
df <- df %>%
  mutate(group = Group) %>%
  separate(Group, into = c("symbiosis", "tissue"), sep = 3)

sample_sizes <- df %>% count(Taxon, group, name = "n")

# ---- 1. ANOVA assumption checks ----
assumption_checks <- df %>%
  group_by(Taxon) %>%
  summarize(
    shapiro_p = {
      mod <- aov(Percentage ~ tissue * symbiosis, data = cur_data())
      shapiro.test(residuals(mod))$p.value
    },
    levene_p = {
      car::leveneTest(Percentage ~ interaction(tissue, symbiosis),
                      data = cur_data())$`Pr(>F)`[1]
    },
    .groups = "drop"
  )
write.csv(assumption_checks, "Supplementary_Table_assumption_checks_Fig2.csv",
          row.names = FALSE)

# ---- 2. Full two-way ANOVA tables ----
full_anova_results <- df %>%
  group_by(Taxon) %>%
  reframe({
    mod <- aov(Percentage ~ tissue * symbiosis, data = cur_data())
    at  <- summary(mod)[[1]]
    total_ss <- sum(at$`Sum Sq`)
    tibble(
      Source       = rownames(at),
      Df           = at$Df,
      Sum_Sq       = round(at$`Sum Sq`, 4),
      Mean_Sq      = round(at$`Mean Sq`, 4),
      F_value      = round(at$`F value`, 3),
      p_value      = at$`Pr(>F)`,
      Variance_pct = round(at$`Sum Sq` / total_ss * 100, 1)
    )
  })
write.csv(full_anova_results, "Supplementary_Table_ANOVA_Fig2.csv", row.names = FALSE)

# ---- 3. FDR correction across taxa ----
fdr_corrected <- full_anova_results %>%
  filter(!is.na(p_value)) %>%
  group_by(Source) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write.csv(fdr_corrected, "Supplementary_Table_ANOVA_FDR_Fig2.csv", row.names = FALSE)

# ---- 4. Tukey HSD post-hoc + compact letter display ----
posthoc_results <- df %>%
  group_by(Taxon) %>%
  reframe({
    mod <- aov(Percentage ~ group, data = cur_data())
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
write.csv(posthoc_results, "Supplementary_Table_TukeyHSD_Fig2.csv", row.names = FALSE)

cld_results <- df %>%
  group_by(Taxon) %>%
  reframe({
    mod <- aov(Percentage ~ group, data = cur_data())
    overall_p <- summary(mod)[[1]]$`Pr(>F)`[1]
    if (!is.na(overall_p) && overall_p < 0.05) {
      emm <- emmeans(mod, "group")
      cld_out <- cld(emm, Letters = letters, adjust = "tukey")
      tibble(group = as.character(cld_out$group),
             emmean = round(cld_out$emmean, 3),
             letter = trimws(cld_out$.group))
    } else {
      tibble(group = c("ApoEpi", "ApoGas", "SymEpi", "SymGas"),
             emmean = NA_real_, letter = "ns")
    }
  })

# ---- 5. Non-parametric alternative ----
nonparam_results <- df %>%
  group_by(Taxon) %>%
  reframe({
    kt <- kruskal.test(Percentage ~ group, data = cur_data())
    if (kt$p.value < 0.05) {
      dunn <- rstatix::dunn_test(cur_data(), Percentage ~ group,
                                 p.adjust.method = "BH")
      tibble(Test = "Kruskal-Wallis + Dunn",
             KW_chi2 = round(kt$statistic, 3),
             KW_df = kt$parameter, KW_p = kt$p.value,
             Group1 = dunn$group1, Group2 = dunn$group2,
             Dunn_p_adj = dunn$p.adj, Dunn_signif = dunn$p.adj.signif)
    } else {
      tibble(Test = "Kruskal-Wallis",
             KW_chi2 = round(kt$statistic, 3),
             KW_df = kt$parameter, KW_p = kt$p.value,
             Group1 = NA_character_, Group2 = NA_character_,
             Dunn_p_adj = NA_real_, Dunn_signif = NA_character_)
    }
  })
write.csv(nonparam_results, "Supplementary_Table_NonParametric_Fig2.csv",
          row.names = FALSE)

# ---- 6. Figure ----
anova_labels <- df %>%
  group_by(Taxon) %>%
  summarize(
    label = {
      mod <- aov(Percentage ~ tissue * symbiosis, data = cur_data())
      at  <- summary(mod)[[1]]
      total_ss <- sum(at$`Sum Sq`)
      paste0(
        "Tissue: ", round(at$`Sum Sq`[1]/total_ss*100, 1),
        "%, F(", at$Df[1], ",", at$Df[4], ") = ", round(at$`F value`[1], 2),
        ", p = ", signif(at$`Pr(>F)`[1], 3),
        "\nSymbiosis: ", round(at$`Sum Sq`[2]/total_ss*100, 1),
        "%, F(", at$Df[2], ",", at$Df[4], ") = ", round(at$`F value`[2], 2),
        ", p = ", signif(at$`Pr(>F)`[2], 3),
        "\nTissue \u00D7 Symbiosis: ", round(at$`Sum Sq`[3]/total_ss*100, 1),
        "%, F(", at$Df[3], ",", at$Df[4], ") = ", round(at$`F value`[3], 2),
        ", p = ", signif(at$`Pr(>F)`[3], 3)
      )
    },
    max_y = max(Percentage, na.rm = TRUE),
    .groups = "drop"
  )

plot_letters <- cld_results %>%
  left_join(df %>%
              group_by(Taxon, group) %>%
              summarize(y_pos = max(Percentage, na.rm = TRUE), .groups = "drop"),
            by = c("Taxon", "group"))

p <- ggplot(df, aes(x = group, y = Percentage, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  facet_wrap(~Taxon, scales = "free_y") +
  geom_text(data = anova_labels,
            aes(x = 1, y = max_y * 0.9, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.5, lineheight = 0.9) +
  geom_text(data = plot_letters,
            aes(x = group, y = y_pos * 1.08, label = letter),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  labs(y = "Percentage", x = NULL,
       caption = paste0("n = ", unique(sample_sizes$n)[1],
                        " biological replicates per group. ",
                        "Two-way ANOVA (tissue \u00D7 symbiosis); ",
                        "letters indicate Tukey HSD post-hoc groups (p < 0.05).")) +
  scale_fill_manual(values = tissue_col) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text       = element_text(color = "black"),
    legend.position  = "bottom",
    plot.caption     = element_text(hjust = 0, size = 8)
  )

ggsave("Fig2abc_coding_taxonomy.pdf", plot = p,
       width = 12, height = 8, units = "in")
