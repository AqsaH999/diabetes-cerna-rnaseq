# ==========================================================
# 05_biotype_split.R — label DE genes by biotype, split mRNA vs lncRNA, heatmap
# ==========================================================
setwd("~/Desktop/diabetes-cerna-rnaseq")
library(data.table); library(pheatmap)

# --- 1. DESeq2 results ---
res <- read.csv("results/tables/deseq2_T2D_vs_ND.csv", row.names = 1)
res$gene_id <- rownames(res)

# --- 2. Annotation straight from the local Ensembl GTF (gene rows only) ---
gtf  <- fread(cmd = "awk -F'\t' '$3==\"gene\"' data/reference/Homo_sapiens.GRCh38.110.gtf",
              sep = "\t", header = FALSE)
a    <- gtf$V9
ann  <- data.table(
  gene_id      = sub('.*gene_id "([^"]+)".*', "\\1", a),
  gene_name    = ifelse(grepl('gene_name "', a),
                        sub('.*gene_name "([^"]+)".*', "\\1", a), NA_character_),
  gene_biotype = sub('.*gene_biotype "([^"]+)".*', "\\1", a)
)
ann <- unique(ann, by = "gene_id")

# --- 3. Merge annotation onto results ---
res <- merge(as.data.table(res), ann, by = "gene_id", all.x = TRUE)

# --- 4. Significant genes, split by biotype ---
sig       <- res[!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1]
de_mrna   <- sig[gene_biotype == "protein_coding"][order(padj)]
de_lncrna <- sig[gene_biotype == "lncRNA"][order(padj)]

fwrite(de_mrna,   "results/tables/DE_mRNA.csv")
fwrite(de_lncrna, "results/tables/DE_lncRNA.csv")

cat("\nTotal significant:", nrow(sig),
    "| DE mRNAs:", nrow(de_mrna),
    "| DE lncRNAs:", nrow(de_lncrna), "\n\nBiotype breakdown:\n")
print(sort(table(sig$gene_biotype), decreasing = TRUE))

# --- 5. Heatmap of the top 50 DE genes ---
vst  <- readRDS("results/tables/vst_matrix.rds")
meta <- read.csv("data/processed/metadata.csv", row.names = 1)
meta <- meta[colnames(vst), , drop = FALSE]

top <- head(sig[order(padj)]$gene_id, 50)
mat <- vst[rownames(vst) %in% top, , drop = FALSE]

lab <- ann$gene_name[match(rownames(mat), ann$gene_id)]      # use symbols where available
lab[is.na(lab) | lab == ""] <- rownames(mat)[is.na(lab) | lab == ""]
rownames(mat) <- lab

mat <- t(scale(t(mat)))                                       # z-score each gene
ann_col <- data.frame(Group = meta$condition); rownames(ann_col) <- rownames(meta)

pdf("results/figures/heatmap_top50.pdf", width = 9, height = 10)
pheatmap(mat, annotation_col = ann_col, show_colnames = FALSE,
         main = "Top 50 DE genes — T2D vs ND", fontsize_row = 7)
dev.off()

cat("\nSaved: DE_mRNA.csv, DE_lncRNA.csv, heatmap_top50.pdf\n")