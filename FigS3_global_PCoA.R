# Fig. S3: Global PCoA and two-way PERMANOVA (State x Layer) on all 16 samples,
# plus PERMDISP dispersion check.
#
# Inputs : kallisto/                   (pseudoalignment outputs)
# Output : FigS3_global_PCoA.pdf
# Run with this folder as the working directory.

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(sleuth)
library(UpSetR)
library(pheatmap)
library(vegan)

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

sva_fil_norm_kt  <- kallisto_table(so_all, use_filtered = TRUE,
                                   normalized = TRUE, include_covariates = FALSE)
sva_fil_norm_tpm <- sva_fil_norm_kt %>%
  select(gene_id = target_id, sample, tpm) %>%
  spread(gene_id, tpm)

expr_mat <- sva_fil_norm_tpm
rownames(expr_mat) <- expr_mat$sample
expr_mat$sample    <- NULL
expr_mat <- as.matrix(expr_mat)

# log2(TPM + 1) stabilizes the heavy right tail of expression data
expr_log <- log2(expr_mat + 1)

meta <- data.frame(
  sample = rownames(expr_log),
  State  = factor(rep(c("Apo", "Sym"), each = 8),          levels = c("Apo", "Sym")),
  Layer  = factor(rep(rep(c("Epi", "Gas"), each = 4), 2),  levels = c("Epi", "Gas"))
)
stopifnot(all(meta$sample == rownames(expr_log)))

# ---- Distance matrices ----
dist_bc <- vegdist(expr_mat, method = "bray")       # compositional
dist_eu <- dist(expr_log,   method = "euclidean")   # log-Euclidean

# ---- PERMANOVA: two-way with interaction ----
set.seed(42)
permanova_bc <- adonis2(dist_bc ~ State * Layer, data = meta,
                        permutations = 9999, by = "terms")
set.seed(42)
permanova_eu <- adonis2(dist_eu ~ State * Layer, data = meta,
                        permutations = 9999, by = "terms")
print(permanova_bc)
print(permanova_eu)

# ---- PERMDISP: within-group dispersions ----
bd_state <- betadisper(dist_eu, meta$State); print(anova(bd_state))
bd_layer <- betadisper(dist_eu, meta$Layer); print(anova(bd_layer))

# ---- PCoA visualization ----
pcoa <- cmdscale(dist_eu, k = 2, eig = TRUE)
ve   <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)
pcoa_df <- data.frame(PCo1  = pcoa$points[, 1],
                      PCo2  = pcoa$points[, 2],
                      State = meta$State,
                      Layer = meta$Layer)

p_global <- ggplot(pcoa_df, aes(PCo1, PCo2, colour = State, shape = Layer)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#67a9cf", "#ef8a62")) +
  xlab(paste0("PCo1 (", ve[1], "%)")) +
  ylab(paste0("PCo2 (", ve[2], "%)")) +
  theme_bw() +
  theme(panel.grid   = element_blank(),
        panel.border = element_blank(),
        axis.line    = element_line(colour = "black"))

ggsave("FigS3_global_PCoA.pdf", p_global, width = 6, height = 4, dpi = 300)
