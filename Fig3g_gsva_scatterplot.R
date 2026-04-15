# Fig. 3g: GSVA pairwise limma contrasts (Sym vs Apo, Gas vs Epi) rendered
# as a scatterplot of t-statistics.
#
# Inputs : kallisto/                        (pseudoalignment outputs)
#          aipBac.kegg.gs.all.RData         (KEGG gene-set object `gs`)
#          map_des_category.xlsx            (pathway descriptors + category)
# Output : Fig3g_gsva_scatterplot.pdf
#          Fig3g_gsva_Sym_vs_Apo_limma.csv
#          Fig3g_gsva_Gas_vs_Epi_limma.csv
#          Fig3g_gsva_Sym_vs_Apo_and_Gas_vs_Epi_t.csv
# Run with this folder as the working directory.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(sleuth)
library(GSVA)
library(tidyverse)
library(limma)
library(ComplexHeatmap)
library(ggrepel)

base_dir <- "."

sample_id <- dir(file.path(base_dir, "kallisto"))
kal_dirs  <- sapply(sample_id, function(id) file.path(base_dir, "kallisto", id))

s2c <- data.frame(
  sample    = sample_id,
  condition = rep(c("ApoEpi", "ApoGas", "SymEpi", "SymGas"), each = 4),
  state     = rep(c("Apo", "Sym"), each = 8),
  tissue    = rep(rep(c("Epi", "Gas"), each = 4), times = 2)
)
s2c_kal <- dplyr::mutate(s2c, path = kal_dirs)

so_ansn <- sleuth::sleuth_prep(s2c_kal, ~ condition)
so_ansn <- sleuth::sleuth_fit(so_ansn)
so_ansn <- sleuth::sleuth_wt(so_ansn, "conditionSymGas") #no actual DE, only TPM

sva_norm_kt <- kallisto_table(so_ansn, use_filtered = FALSE,
                              normalized = TRUE, include_covariates = FALSE)

expr <- sva_norm_kt %>% select(target_id, sample, tpm) %>% spread(sample, tpm)
aip.expr <- expr[, 2:17]
row.names(aip.expr) <- expr$target_id

load("aipBac.kegg.gs.all.RData")   # provides `gs`

# ---- GSVA ----
gsvapar <- gsvaParam(as.matrix(aip.expr), gs, maxDiff = FALSE)
gsva_es <- gsva(gsvapar)
gsva_es <- as.data.frame(gsva_es)
gsva_es$ko <- rownames(gsva_es)

ko_des <- readxl::read_xlsx("map_des_category.xlsx",
                            col_names = c("ko", "name", "show", "category"))

gsva_snan <- left_join(gsva_es, ko_des, by = "ko")
gsva_t <- gsva_snan[, 1:16]
rownames(gsva_t) <- gsva_snan$name

# ---- limma: Sym vs Apo ----
group_list <- data.frame(sample = colnames(gsva_t),
                         group = rep(c("Apo", "Sym"), each = 8))
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_t)
contrast.matrix <- makeContrasts(Sym - Apo, levels = design)

fit  <- lmFit(gsva_t, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
write.csv(x, "Fig3g_gsva_Sym_vs_Apo_limma.csv", row.names = TRUE)
df <- data.frame(ID = rownames(x), sva = x$t)

# ---- limma: Gas vs Epi ----
group_list <- data.frame(sample = colnames(gsva_t),
                         group = rep(c("Epi", "Gas", "Epi", "Gas"), each = 4))
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_t)
contrast.matrix <- makeContrasts(Gas - Epi, levels = design)

fit  <- lmFit(gsva_t, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
write.csv(x, "Fig3g_gsva_Gas_vs_Epi_limma.csv", row.names = TRUE)
df1 <- data.frame(ID = rownames(x), nvc = x$t)

df_t <- left_join(df, df1, by = "ID")
write.csv(df_t, "Fig3g_gsva_Sym_vs_Apo_and_Gas_vs_Epi_t.csv", row.names = FALSE)

df_t1 <- left_join(df_t, ko_des, by = c("ID" = "name")) %>%
  mutate(category = ifelse(is.na(category), "F", category),
         category = factor(category, levels = c("F", "B", "V", "C", "N")))

cat_color <- c("F" = "#d9d9d9",
               "V" = "#41ae76", "N" = "#f768a1",
               "C" = "#008fcd", "B" = "#984ea3")

ggplot() +
  geom_point(data = df_t1 %>% filter(category == "F"),
             aes(x = sva, y = nvc),
             size = 3, color = "#d9d9d9", alpha = 0.5) +
  geom_point(data = df_t1 %>% filter(category != "F"),
             aes(x = sva, y = nvc, color = category),
             size = 3, alpha = 1) +
  theme_classic() +
  scale_color_manual(values = cat_color) +
  labs(x = "Sym/Apo", y = "Gas/Epi")

ggsave("Fig3g_gsva_scatterplot.pdf", width = 8, height = 6)
