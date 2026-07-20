# ==========================================================
# 04b_figures.R — polished volcano + PCA (run AFTER 04_deseq2.R)
# ==========================================================
setwd("~/Desktop/diabetes-cerna-rnaseq")
library(EnhancedVolcano); library(ggplot2); library(data.table)

# ---- DE results + GTF-based gene names ----
res <- read.csv("results/tables/deseq2_T2D_vs_ND.csv", row.names = 1)
res$ensembl <- sub("\\..*$", "", rownames(res))

ann <- fread("data/reference/gene_annotation.tsv", header = FALSE,
             col.names = c("gene_id","gene_name","gene_biotype"))
ann$gene_id <- sub("\\..*$", "", ann$gene_id)
name_map <- setNames(ann$gene_name, ann$gene_id)

res$symbol <- name_map[res$ensembl]
res$label  <- ifelse(is.na(res$symbol) | res$symbol == "", res$ensembl, res$symbol)

sig    <- which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
sig    <- sig[order(res$padj[sig])]
topLab <- res$label[head(sig, 25)]
up <- sum(res$padj < 0.05 & res$log2FoldChange >  1, na.rm = TRUE)
dn <- sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)

# ---- clean volcano ----
pdf("results/figures/volcano_clean.pdf", width = 8, height = 7)
print(EnhancedVolcano(res,
                      lab = res$label, selectLab = topLab,
                      x = "log2FoldChange", y = "padj",
                      pCutoff = 0.05, FCcutoff = 1.0,
                      pointSize = 1.6, labSize = 3.6, colAlpha = 0.55,
                      drawConnectors = TRUE, widthConnectors = 0.4,
                      maxoverlapsConnectors = 25, arrowheads = FALSE,
                      legendLabels = c("Not sig", "Fold-change only", "p-value only", "Significant"),
                      title = "T2D vs ND - differential expression",
                      subtitle = paste0(up, " up in T2D  \u00b7  ", dn, " down   (padj<0.05, |log2FC|>1)"),
                      caption = "GSE164416 (islets) \u00b7 DESeq2",
                      legendPosition = "top",
                      gridlines.major = FALSE, gridlines.minor = FALSE))
dev.off()

# ---- clean PCA (from saved vst matrix) ----
vst  <- readRDS("results/tables/vst_matrix.rds")
meta <- read.csv("data/processed/metadata.csv", row.names = 1)
meta <- meta[colnames(vst), , drop = FALSE]
v    <- apply(vst, 1, var)
sel  <- order(v, decreasing = TRUE)[1:500]
pc   <- prcomp(t(vst[sel, ]))
pv   <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
df   <- data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], group = meta$condition)

pdf("results/figures/pca_clean.pdf", width = 7, height = 5)
print(ggplot(df, aes(PC1, PC2, color = group)) +
        geom_point(size = 3, alpha = 0.85) +
        stat_ellipse(level = 0.8, linewidth = 0.6) +
        scale_color_manual(values = c(ND = "#E64B35", T2D = "#4DBBD5")) +
        labs(x = paste0("PC1: ", pv[1], "% variance"),
             y = paste0("PC2: ", pv[2], "% variance"),
             title = "PCA - ND vs T2D islets", color = "Group") +
        theme_bw(base_size = 13) +
        theme(panel.grid.minor = element_blank(),
              plot.margin = margin(10, 20, 10, 10)))
dev.off()

cat("Saved volcano_clean.pdf and pca_clean.pdf\n")