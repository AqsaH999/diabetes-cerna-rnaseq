# ==========================================================
# 17_readme_pngs.R - regenerate heatmap + GO dotplots as high-res PNG,
#                    white background (readable on GitHub dark theme)
# ==========================================================
library(data.table); library(ggplot2)
dir.create("results/figures/png", showWarnings = FALSE)

## ---------------- Heatmap (native PNG, white background) ----------------
library(pheatmap)
res <- read.csv("results/tables/deseq2_T2D_vs_ND.csv", row.names = 1)
res$ens <- sub("\\..*$", "", rownames(res))
ann <- fread("data/reference/gene_annotation.tsv", header = FALSE,
             col.names = c("gene_id","gene_name","gene_biotype"))
ann$gene_id <- sub("\\..*$", "", ann$gene_id)
res$symbol  <- setNames(ann$gene_name, ann$gene_id)[res$ens]

sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
sig <- sig[order(sig$padj), ]
top <- head(rownames(sig), 50)

vst  <- readRDS("results/tables/vst_matrix.rds")
meta <- read.csv("data/processed/metadata.csv", row.names = 1)
meta <- meta[colnames(vst), , drop = FALSE]

mat <- vst[rownames(vst) %in% top, , drop = FALSE]
lab <- res$symbol[match(rownames(mat), res$ens)]
lab[is.na(lab) | lab == ""] <- rownames(mat)[is.na(lab) | lab == ""]
rownames(mat) <- lab
mat <- t(scale(t(mat)))
ann_col <- data.frame(Group = meta$condition); rownames(ann_col) <- rownames(meta)

png("results/figures/png/heatmap.png", width = 2000, height = 2400, res = 220, bg = "white")
pheatmap(mat, annotation_col = ann_col, show_colnames = FALSE,
         main = "Top 50 DE genes (T2D vs ND)", fontsize_row = 7,
         annotation_colors = list(Group = c(ND = "#E64B35", T2D = "#4DBBD5")))
dev.off()

## ---------------- GO dotplots from saved tables (white background) -------
go_dot <- function(csv, title, out) {
  d <- fread(csv)
  d <- head(d[order(p.adjust)], 15)
  d[, gr := sapply(strsplit(GeneRatio, "/"),
                   function(x) as.numeric(x[1]) / as.numeric(x[2]))]
  d[, Description := factor(Description, levels = rev(Description))]
  p <- ggplot(d, aes(gr, Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#E64B35", high = "#4477AA", trans = "log10") +
    scale_size(range = c(2, 7)) +
    labs(x = "GeneRatio", y = NULL, title = title, color = "p.adjust", size = "Count") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12))
  ggsave(out, p, width = 8, height = 6, dpi = 220, bg = "white")
}

up_csv <- if (file.exists("results/tables/GO_BP_up_simplified.csv"))
            "results/tables/GO_BP_up_simplified.csv" else "results/tables/GO_BP_up.csv"
dn_csv <- if (file.exists("results/tables/GO_BP_down_simplified.csv"))
            "results/tables/GO_BP_down_simplified.csv" else "results/tables/GO_BP_down_relaxed.csv"

go_dot(up_csv, "GO Biological Process - up in T2D",   "results/figures/png/GO_up.png")
go_dot(dn_csv, "GO Biological Process - down in T2D", "results/figures/png/GO_down.png")

cat("Regenerated heatmap.png, GO_up.png, GO_down.png (high-res, white background)\n")
