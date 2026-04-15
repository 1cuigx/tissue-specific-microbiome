# Fig. 3f: GSVA pathway activity heatmap with two-way factorial limma analysis.
#
# Inputs : kallisto/                   (pseudoalignment outputs)
#          aipBac.kegg.gs.all.RData    (KEGG gene-set object `gs`)
#          map_des.tsv                 (pathway descriptors with `show` flag)
# Output : Fig3f_gsva_heatmap.pdf
#          Fig3f_gsva_scores.csv
#          Supplementary_Table_GSVA_Limma_full.csv
# Run with this folder as the working directory.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(sleuth)
library(GSVA)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(limma)

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

so_all <- sleuth::sleuth_prep(s2c_kal, ~ condition)
so_all <- sleuth::sleuth_fit(so_all)
so_all <- sleuth::sleuth_wt(so_all, "conditionSymGas") # no actual DE, only TPM

sva_norm_kt <- kallisto_table(so_all, use_filtered = FALSE,
                              normalized = TRUE, include_covariates = FALSE)

load("aipBac.kegg.gs.all.RData")  # provides `gs`

expr <- sva_norm_kt %>% select(target_id, sample, tpm) %>% spread(sample, tpm)
aip.expr <- expr[, 2:17]
row.names(aip.expr) <- expr$target_id

# ---- GSVA pathway-level activity ----
gsvapar <- gsvaParam(as.matrix(aip.expr), gs, maxDiff = FALSE)
gsva_es <- gsva(gsvapar)
gsva_es <- as.data.frame(gsva_es)
gsva_es$ko <- rownames(gsva_es)

ko_des <- read_tsv("map_des.tsv", col_names = c("ko", "name", "show"))

gsva_snan <- left_join(gsva_es, ko_des, by = "ko")
gsva_t <- gsva_snan[, 1:16]
rownames(gsva_t) <- gsva_snan$name

tShow <- ko_des %>% filter(show == TRUE)
tShow_label <- tShow$name

# ---- Heatmap ----
Condition <- rep(c("ApoEpi", "ApoGas", "SymEpi", "SymGas"), each = 4)

ha <- HeatmapAnnotation(type = Condition,
                        col = list(type = c(
                          "ApoEpi" = "#A0CAE0",
                          "ApoGas" = "#3082BC",
                          "SymEpi" = "#FDAF69",
                          "SymGas" = "#E55409"
                        )),
                        annotation_name_side = "left")

ch <- Heatmap(gsva_t,
              show_row_names        = FALSE,
              top_annotation        = ha,
              clustering_method_rows = "complete",
              cluster_columns       = FALSE,
              name                  = "GSVA score",
              border                = FALSE,
              show_column_names     = FALSE) +
  rowAnnotation(foo = anno_mark(
    at     = which(rownames(gsva_t) %in% tShow_label),
    labels = rownames(gsva_t)[rownames(gsva_t) %in% tShow_label],
    labels_gp = gpar(fontsize = 8),
    padding   = unit(2.5, "mm")))

pdf("Fig3f_gsva_heatmap.pdf", width = 7, height = 10)
draw(ch, merge_legends = TRUE, heatmap_legend_side = "right")
dev.off()

write.csv(gsva_snan, "Fig3f_gsva_scores.csv", row.names = FALSE)

# ---- Two-way factorial limma (~ state * tissue) ----
s2c_limma <- s2c
s2c_limma$state  <- factor(s2c_limma$state,  levels = c("Apo", "Sym"))
s2c_limma$tissue <- factor(s2c_limma$tissue, levels = c("Epi", "Gas"))

design <- model.matrix(~ state * tissue, data = s2c_limma)
fit <- lmFit(gsva_t, design)
fit <- eBayes(fit)

# Joint F-test across main effects + interaction
limma_joint <- topTable(fit, coef = 2:4, number = Inf, adjust.method = "BH") %>%
  tibble::rownames_to_column("Pathway") %>%
  select(Pathway, F_joint = F, FDR_joint = adj.P.Val)

res_symbiosis <- topTable(fit, coef = "stateSym",
                          number = Inf, adjust.method = "BH") %>%
  tibble::rownames_to_column("Pathway") %>%
  select(Pathway, logFC_Symbiosis = logFC, t_Symbiosis = t,
         P_Symbiosis = P.Value, FDR_Symbiosis = adj.P.Val)

res_tissue <- topTable(fit, coef = "tissueGas",
                       number = Inf, adjust.method = "BH") %>%
  tibble::rownames_to_column("Pathway") %>%
  select(Pathway, logFC_Tissue = logFC, t_Tissue = t,
         P_Tissue = P.Value, FDR_Tissue = adj.P.Val)

res_interaction <- topTable(fit, coef = "stateSym:tissueGas",
                            number = Inf, adjust.method = "BH") %>%
  tibble::rownames_to_column("Pathway") %>%
  select(Pathway, logFC_Interaction = logFC, t_Interaction = t,
         P_Interaction = P.Value, FDR_Interaction = adj.P.Val)

limma_full <- limma_joint %>%
  left_join(res_symbiosis,   by = "Pathway") %>%
  left_join(res_tissue,      by = "Pathway") %>%
  left_join(res_interaction, by = "Pathway")

write.csv(limma_full, "Supplementary_Table_GSVA_Limma_full.csv", row.names = FALSE)

cat("Joint F-test  FDR<0.05:", sum(limma_joint$FDR_joint           < 0.05, na.rm = TRUE), "\n")
cat("Symbiosis main FDR<0.05:", sum(res_symbiosis$FDR_Symbiosis    < 0.05, na.rm = TRUE), "\n")
cat("Tissue main    FDR<0.05:", sum(res_tissue$FDR_Tissue          < 0.05, na.rm = TRUE), "\n")
cat("Interaction    FDR<0.05:", sum(res_interaction$FDR_Interaction < 0.05, na.rm = TRUE), "\n")
