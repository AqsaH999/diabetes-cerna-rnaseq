# ==========================================================
# 10_cerna_correlation.R — filter ceRNA axes by expression correlation
#   Test 1: across all 57 donors (confounded by ND/T2D group difference)
#   Test 2: within T2D donors only (removes the group effect - stricter)
# ==========================================================
setwd("~/Desktop/diabetes-cerna-rnaseq")   # remove during Day 12 cleanup
library(data.table)

triads <- fread("results/tables/cerna_triads.csv")
vst    <- readRDS("results/tables/vst_matrix.rds")
meta   <- read.csv("data/processed/metadata.csv", row.names = 1)
meta   <- meta[colnames(vst), , drop = FALSE]

# ---- map gene symbols -> Ensembl IDs (vst is keyed by Ensembl) ----
ann <- fread("data/reference/gene_annotation.tsv", header = FALSE,
             col.names = c("gene_id", "gene_name", "gene_biotype"))
ann[, gene_id := sub("\\..*$", "", gene_id)]
ann <- ann[gene_name != "" & gene_id %in% rownames(vst)]
ann[, key := toupper(gene_name)]
sym2ens <- setNames(ann$gene_id, ann$key)[!duplicated(ann$key)]

pairs <- unique(triads[, .(lncRNA, target_symbol)])
pairs[, lnc_id  := sym2ens[toupper(lncRNA)]]
pairs[, mrna_id := sym2ens[toupper(target_symbol)]]
cat("lncRNA-mRNA pairs:", nrow(pairs), "| mappable:",
    sum(!is.na(pairs$lnc_id) & !is.na(pairs$mrna_id)), "\n")
pairs <- pairs[!is.na(lnc_id) & !is.na(mrna_id)]

# ---- helper: Spearman correlation over a given sample set ----
corr_over <- function(mat) {
  rho <- numeric(nrow(pairs)); pv <- numeric(nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    ct <- suppressWarnings(cor.test(mat[pairs$lnc_id[i], ],
                                    mat[pairs$mrna_id[i], ], method = "spearman"))
    rho[i] <- unname(ct$estimate); pv[i] <- ct$p.value
  }
  list(rho = rho, padj = p.adjust(pv, method = "BH"))
}

# ---- Test 1: all 57 donors ----
a <- corr_over(vst)
pairs[, `:=`(rho_all = a$rho, padj_all = a$padj)]

# ---- Test 2: within T2D donors only ----
t2d <- rownames(meta)[meta$condition == "T2D"]
cat("T2D donors used for within-group test:", length(t2d), "\n")
b <- corr_over(vst[, t2d])
pairs[, `:=`(rho_t2d = b$rho, padj_t2d = b$padj)]

fwrite(pairs, "results/tables/cerna_pair_correlation.csv")

# ---- Report both ----
bt_all <- binom.test(sum(pairs$rho_all > 0), nrow(pairs), 0.5)
bt_t2d <- binom.test(sum(pairs$rho_t2d > 0), nrow(pairs), 0.5)

cat("\n=== TEST 1: all 57 donors (ND + T2D) ===\n")
cat("Positive rho:", sum(pairs$rho_all > 0), "of", nrow(pairs),
    sprintf("(%.0f%%, binomial p = %.2g)\n",
            100*mean(pairs$rho_all > 0), bt_all$p.value))
cat("Positive AND padj<0.05:", sum(pairs$rho_all > 0 & pairs$padj_all < 0.05),
    "| median rho:", round(median(pairs$rho_all), 3), "\n")

cat("\n=== TEST 2: within T2D donors only (group effect removed) ===\n")
cat("Positive rho:", sum(pairs$rho_t2d > 0), "of", nrow(pairs),
    sprintf("(%.0f%%, binomial p = %.2g)\n",
            100*mean(pairs$rho_t2d > 0), bt_t2d$p.value))
cat("Positive AND padj<0.05:", sum(pairs$rho_t2d > 0 & pairs$padj_t2d < 0.05),
    "| median rho:", round(median(pairs$rho_t2d), 3), "\n")

cat("\nAgreement between the two tests (Spearman of rho_all vs rho_t2d):",
    round(cor(pairs$rho_all, pairs$rho_t2d, method = "spearman"), 3), "\n")

# ---- Build supported sets under each criterion ----
kept_all <- pairs[rho_all > 0 & padj_all < 0.05]
kept_t2d <- pairs[rho_t2d > 0 & padj_t2d < 0.05]
kept_both <- pairs[rho_all > 0 & padj_all < 0.05 & rho_t2d > 0]

fwrite(kept_all,  "results/tables/cerna_pairs_supported.csv")
fwrite(kept_t2d,  "results/tables/cerna_pairs_supported_t2d.csv")

mk_tri <- function(k) merge(triads, k[, .(lncRNA, target_symbol)],
                            by = c("lncRNA", "target_symbol"))
tri_all  <- mk_tri(kept_all)
tri_t2d  <- mk_tri(kept_t2d)
tri_both <- mk_tri(kept_both)

fwrite(tri_all,  "results/tables/cerna_triads_supported.csv")
fwrite(tri_t2d,  "results/tables/cerna_triads_supported_t2d.csv")
fwrite(tri_both, "results/tables/cerna_triads_supported_both.csv")

cat("\n=== AXES RETAINED (of", nrow(triads), ") ===\n")
cat("All-donor criterion :", nrow(tri_all),
    sprintf("| %d lncRNA, %d miRNA, %d mRNA\n",
            uniqueN(tri_all$lncRNA), uniqueN(tri_all$mature_mirna_id),
            uniqueN(tri_all$target_symbol)))
cat("Within-T2D criterion:", nrow(tri_t2d),
    sprintf("| %d lncRNA, %d miRNA, %d mRNA\n",
            uniqueN(tri_t2d$lncRNA), uniqueN(tri_t2d$mature_mirna_id),
            uniqueN(tri_t2d$target_symbol)))
cat("Both criteria       :", nrow(tri_both),
    sprintf("| %d lncRNA, %d miRNA, %d mRNA\n",
            uniqueN(tri_both$lncRNA), uniqueN(tri_both$mature_mirna_id),
            uniqueN(tri_both$target_symbol)))

cat("\n--- Top 15 pairs by within-T2D rho ---\n")
print(head(kept_t2d[order(-rho_t2d), .(lncRNA, target_symbol,
                                       rho_all = round(rho_all, 3),
                                       rho_t2d = round(rho_t2d, 3),
                                       padj_t2d = signif(padj_t2d, 3))], 15))