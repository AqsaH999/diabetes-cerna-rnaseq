# ==========================================================
# 06_enrichment.R — GO + KEGG enrichment on DE mRNAs
# ==========================================================
setwd("~/Desktop/diabetes-cerna-rnaseq")
library(clusterProfiler); library(org.Hs.eg.db); library(data.table); library(ggplot2)

options(timeout = 600)   # KEGG's server is slow; the 60s default often fails

de  <- fread("results/tables/DE_mRNA.csv")
res <- fread("results/tables/deseq2_T2D_vs_ND.csv"); setnames(res, 1, "gene_id")

# --- Map Ensembl IDs -> Entrez IDs (what the databases use) ---
eg   <- bitr(de$gene_id,  fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
uni  <- bitr(res$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
up   <- bitr(de[log2FoldChange > 0]$gene_id, "ENSEMBL", "ENTREZID", org.Hs.eg.db)
down <- bitr(de[log2FoldChange < 0]$gene_id, "ENSEMBL", "ENTREZID", org.Hs.eg.db)

# Relaxed down-set: padj<0.05 only (no fold-change cutoff).
# The strict |log2FC|>1 filter leaves too few down genes for statistical power.
down_relaxed <- bitr(res[!is.na(padj) & padj < 0.05 & log2FoldChange < 0]$gene_id,
                     "ENSEMBL", "ENTREZID", org.Hs.eg.db)

cat("DE mRNAs mapped:", nrow(eg), "| up:", nrow(up), "| down:", nrow(down),
    "| down (relaxed):", nrow(down_relaxed),
    "| background universe:", nrow(uni), "\n")

# --- GO Biological Process ---
ego         <- enrichGO(eg$ENTREZID,   universe = uni$ENTREZID, OrgDb = org.Hs.eg.db,
                        ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ego_up      <- enrichGO(up$ENTREZID,   universe = uni$ENTREZID, OrgDb = org.Hs.eg.db,
                        ont = "BP", qvalueCutoff = 0.05, readable = TRUE)
ego_down    <- enrichGO(down$ENTREZID, universe = uni$ENTREZID, OrgDb = org.Hs.eg.db,
                        ont = "BP", qvalueCutoff = 0.05, readable = TRUE)
ego_down_rx <- enrichGO(down_relaxed$ENTREZID, universe = uni$ENTREZID, OrgDb = org.Hs.eg.db,
                        ont = "BP", qvalueCutoff = 0.05, readable = TRUE)

fwrite(as.data.frame(ego),         "results/tables/GO_BP.csv")
fwrite(as.data.frame(ego_up),      "results/tables/GO_BP_up.csv")
fwrite(as.data.frame(ego_down),    "results/tables/GO_BP_down.csv")
fwrite(as.data.frame(ego_down_rx), "results/tables/GO_BP_down_relaxed.csv")

# --- KEGG (needs internet; wrapped so a server timeout can't halt the script) ---
kk <- tryCatch({
  k <- enrichKEGG(gene = eg$ENTREZID, universe = uni$ENTREZID,
                  organism = "hsa", qvalueCutoff = 0.05)
  setReadable(k, org.Hs.eg.db, keyType = "ENTREZID")
}, error = function(e) { message("KEGG unavailable (skipped): ", e$message); NULL })

if (!is.null(kk)) fwrite(as.data.frame(kk), "results/tables/KEGG.csv")

# --- Figures (only drawn if there are results) ---
if (nrow(as.data.frame(ego)) > 0)
  ggsave("results/figures/GO_dotplot.pdf",      dotplot(ego,         showCategory = 20), width = 9, height = 8)
if (nrow(as.data.frame(ego_up)) > 0)
  ggsave("results/figures/GO_dotplot_up.pdf",   dotplot(ego_up,      showCategory = 20), width = 9, height = 8)
if (nrow(as.data.frame(ego_down_rx)) > 0)
  ggsave("results/figures/GO_dotplot_down.pdf", dotplot(ego_down_rx, showCategory = 20), width = 9, height = 8)
if (!is.null(kk) && nrow(as.data.frame(kk)) > 0)
  ggsave("results/figures/KEGG_dotplot.pdf",    dotplot(kk,          showCategory = 20), width = 9, height = 8)

# --- Console summary ---
cat("\n--- Top GO terms (UP in T2D) ---\n")
print(head(as.data.frame(ego_up)[, c("Description","GeneRatio","p.adjust","Count")], 15))

cat("\n--- Top GO terms (DOWN in T2D, relaxed padj<0.05) ---\n")
print(head(as.data.frame(ego_down_rx)[, c("Description","GeneRatio","p.adjust","Count")], 15))

if (!is.null(kk)) {
  cat("\n--- Top KEGG pathways ---\n")
  print(head(as.data.frame(kk)[, c("Description","GeneRatio","p.adjust","Count")], 15))
}

cat("\nGO terms: all =", nrow(as.data.frame(ego)),
    "| up =", nrow(as.data.frame(ego_up)),
    "| down (strict) =", nrow(as.data.frame(ego_down)),
    "| down (relaxed) =", nrow(as.data.frame(ego_down_rx)), "\n")

# ==========================================================
# Targeted check: are beta-cell / ECM terms present but ranked lower?
# (Enrichment favours large gene sets, so small but real programs can
#  sit far down the list rather than being absent.)
# ==========================================================
d <- as.data.frame(ego_down_rx)
u <- as.data.frame(ego_up)

cat("\n--- DOWN terms matching insulin / secretion / glucose / hormone / ion transport ---\n")
hit_d <- grep("insulin|secretion|glucose|hormone|ion transport", d$Description, ignore.case = TRUE)
if (length(hit_d) > 0) {
  out_d <- data.frame(Rank = hit_d, d[hit_d, c("Description","p.adjust","Count")])
  print(out_d, row.names = FALSE)
} else cat("(none found)\n")

cat("\n--- UP terms matching extracellular matrix / collagen / complement ---\n")
hit_u <- grep("extracellular matrix|collagen|complement", u$Description, ignore.case = TRUE)
if (length(hit_u) > 0) {
  out_u <- data.frame(Rank = hit_u, u[hit_u, c("Description","p.adjust","Count")])
  print(out_u, row.names = FALSE)
} else cat("(none found)\n")