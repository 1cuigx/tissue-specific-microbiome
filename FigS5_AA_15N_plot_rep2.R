# Fig. S5: 15N enrichment of amino acids (replicate 2), with two-way ANOVA
# (Sym x Germ), Tukey HSD post-hoc, non-parametric alternatives, and
# cumulative pairwise t-tests.
#
# Input : FigS5_AA_rep2.xlsx
# Output: FigS5_AA_delta15N_boxplot_rep2.pdf
#         Supplementary_Table_assumption_checks_FigS5.csv
#         Supplementary_Table_ANOVA_FigS5.csv
#         Supplementary_Table_ANOVA_FDR_FigS5.csv
#         Supplementary_Table_TukeyHSD_FigS5.csv
#         Supplementary_Table_NonParametric_FigS5.csv
#         Supplementary_Table_Cumulative_Pairwise_FigS5.csv
# Run with this folder as the working directory.

library(tidyverse)
library(readxl)
library(ggpubr)
library(car)
library(emmeans)
library(multcomp)
library(rstatix)

df <- read_xlsx("FigS5_AA_rep2.xlsx")

df <- df %>%
  mutate(
    Sym  = ifelse(str_starts(Animal, "H"), "Sym", "Apo"),
    Germ = case_when(
      str_detect(Animal, "R")  ~ "Germ-present",
      str_detect(Animal, "GF") ~ "Germ-depleted",
      TRUE ~ NA_character_
    ),
    ratio = N15 / N14
  )

normalized_df <- df %>%
  group_by(Animal, AAs) %>%
  mutate(control_avg = mean(ratio[Treatment == "Control"], na.rm = TRUE),
         delta_15N   = ((ratio - control_avg) / control_avg) * 1000) %>%
  ungroup()


# ---- 1. Sample sizes ----
sample_sizes_figS5 <- normalized_df %>% count(AAs, Animal, name = "n")

# ---- 2. Assumption checks ----
assumption_checks_figS5 <- normalized_df %>%
  group_by(AAs) %>%
  summarize(
    shapiro_p = {
      mod <- aov(delta_15N ~ Sym * Germ, data = cur_data())
      shapiro.test(residuals(mod))$p.value
    },
    levene_p = {
      car::leveneTest(delta_15N ~ interaction(Sym, Germ),
                      data = cur_data())$`Pr(>F)`[1]
    },
    .groups = "drop"
  )
write.csv(assumption_checks_figS5,
          "Supplementary_Table_assumption_checks_FigS5.csv", row.names = FALSE)

# ---- 3. Full two-way ANOVA ----
full_anova_figS5 <- normalized_df %>%
  group_by(AAs) %>%
  reframe({
    mod <- aov(delta_15N ~ Sym * Germ, data = cur_data())
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
write.csv(full_anova_figS5, "Supplementary_Table_ANOVA_FigS5.csv",
          row.names = FALSE)

# ---- 4. FDR correction ----
fdr_corrected_figS5 <- full_anova_figS5 %>%
  filter(!is.na(p_value)) %>%
  group_by(Source) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "BH")) %>%
  ungroup()
write.csv(fdr_corrected_figS5, "Supplementary_Table_ANOVA_FDR_FigS5.csv",
          row.names = FALSE)

# ---- 5. Tukey HSD ----
posthoc_figS5 <- normalized_df %>%
  group_by(AAs) %>%
  reframe({
    mod <- aov(delta_15N ~ Animal, data = cur_data())
    overall_p <- summary(mod)[[1]]$`Pr(>F)`[1]
    if (!is.na(overall_p) && overall_p < 0.05) {
      tukey <- TukeyHSD(mod, "Animal")$Animal
      tibble(Comparison = rownames(tukey),
             Difference = round(tukey[, "diff"], 4),
             Lower_CI   = round(tukey[, "lwr"], 4),
             Upper_CI   = round(tukey[, "upr"], 4),
             p_adjusted = round(tukey[, "p adj"], 6),
             Significant = tukey[, "p adj"] < 0.05)
    } else {
      tibble(Comparison = "ANOVA not significant",
             Difference = NA_real_, Lower_CI = NA_real_,
             Upper_CI = NA_real_, p_adjusted = NA_real_, Significant = NA)
    }
  })
write.csv(posthoc_figS5, "Supplementary_Table_TukeyHSD_FigS5.csv",
          row.names = FALSE)

cld_figS5 <- normalized_df %>%
  group_by(AAs) %>%
  reframe({
    mod <- aov(delta_15N ~ Animal, data = cur_data())
    overall_p <- summary(mod)[[1]]$`Pr(>F)`[1]
    if (!is.na(overall_p) && overall_p < 0.05) {
      emm <- emmeans(mod, "Animal")
      cld_out <- cld(emm, Letters = letters, adjust = "tukey")
      tibble(Animal = as.character(cld_out$Animal),
             emmean = round(cld_out$emmean, 3),
             letter = trimws(cld_out$.group))
    } else {
      tibble(Animal = c("AGF", "AR", "HGF", "HR"),
             emmean = NA_real_, letter = "ns")
    }
  })

# ---- 6. Non-parametric alternative ----
nonparam_figS5 <- normalized_df %>%
  group_by(AAs) %>%
  reframe({
    kt <- kruskal.test(delta_15N ~ Animal, data = cur_data())
    if (kt$p.value < 0.05) {
      dunn <- rstatix::dunn_test(cur_data(), delta_15N ~ Animal,
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
write.csv(nonparam_figS5, "Supplementary_Table_NonParametric_FigS5.csv",
          row.names = FALSE)

# ---- 7. Cumulative Δδ15N across amino acids ----
cumulative_df <- normalized_df %>%
  group_by(Animal, AAs) %>%
  mutate(Rep = row_number()) %>%
  ungroup() %>%
  group_by(Animal, Rep) %>%
  summarize(cumulative_delta15N = sum(delta_15N, na.rm = TRUE),
            Sym  = first(Sym),
            Germ = first(Germ),
            .groups = "drop")

cumulative_anova <- aov(cumulative_delta15N ~ Sym * Germ, data = cumulative_df)
summary(cumulative_anova)

comparisons <- list(c("AGF", "AR"),
                    c("AGF", "HGF"),
                    c("AR", "HR"),
                    c("HGF", "HR"))

cumulative_pairwise <- compare_means(
  cumulative_delta15N ~ Animal,
  data = cumulative_df,
  method = "t.test",
  p.adjust.method = "BH",
  comparisons = comparisons
)
write.csv(cumulative_pairwise,
          "Supplementary_Table_Cumulative_Pairwise_FigS5.csv", row.names = FALSE)

# ---- 8. Figure ----
anova_labels_figS5 <- normalized_df %>%
  group_by(AAs) %>%
  summarize(
    label = {
      mod <- aov(delta_15N ~ Sym * Germ, data = cur_data())
      at  <- summary(mod)[[1]]
      total_ss <- sum(at$`Sum Sq`)
      paste0(
        "Symbiosis: ",     round(at$`Sum Sq`[1]/total_ss*100, 1),
        "%, F(", at$Df[1], ",", at$Df[4], ") = ", round(at$`F value`[1], 2),
        ", p = ", signif(at$`Pr(>F)`[1], 3),
        "\nGerm: ",        round(at$`Sum Sq`[2]/total_ss*100, 1),
        "%, F(", at$Df[2], ",", at$Df[4], ") = ", round(at$`F value`[2], 2),
        ", p = ", signif(at$`Pr(>F)`[2], 3),
        "\nGerm \u00D7 Symbiosis: ", round(at$`Sum Sq`[3]/total_ss*100, 1),
        "%, F(", at$Df[3], ",", at$Df[4], ") = ", round(at$`F value`[3], 2),
        ", p = ", signif(at$`Pr(>F)`[3], 3)
      )
    },
    max_y = max(delta_15N, na.rm = TRUE),
    .groups = "drop"
  )

plot_letters_figS5 <- cld_figS5 %>%
  left_join(normalized_df %>%
              group_by(AAs, Animal) %>%
              summarize(y_pos = max(delta_15N, na.rm = TRUE), .groups = "drop"),
            by = c("AAs", "Animal"))

col <- c("AGF" = "#A0CAE0", "AR" = "#3082BC",
         "HGF" = "#f9b070", "HR" = "#e35f39")

p_revised <- ggplot(normalized_df, aes(x = Animal, y = delta_15N, fill = Animal)) +
  geom_boxplot(outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~ AAs, scales = "free_y") +
  geom_text(data = anova_labels_figS5,
            aes(x = 1, y = max_y * 0.85, label = label),
            inherit.aes = FALSE, hjust = 0, size = 2.5, lineheight = 0.9) +
  geom_text(data = plot_letters_figS5,
            aes(x = Animal, y = y_pos * 1.08, label = letter),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = col) +
  labs(y = expression(Delta * delta^15 * N ~ "(\u2030, relative to control)"),
       x = NULL,
       caption = "n = 3 biological replicates per group. Two-way ANOVA (Sym \u00D7 Germ); letters indicate Tukey HSD post-hoc groups (p < 0.05).") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text       = element_text(color = "black"),
        plot.caption     = element_text(hjust = 0, size = 8))

ggsave("FigS5_AA_delta15N_boxplot_rep2.pdf", plot = p_revised,
       width = 16, height = 14, units = "in")
