# ==========================================================
# 16_export_pngs.R — re-render headline figures as PNG for the README
# ==========================================================
library(data.table); library(igraph)
dir.create("results/figures/png", showWarnings = FALSE)

# ---- 1. Validation scatter (Track A vs Track B) ----
# (the Python script already made a PNG; copy it if present)
if (file.exists("results/figures/validation_scatter.png"))
  file.copy("results/figures/validation_scatter.png",
            "results/figures/png/validation_scatter.png", overwrite = TRUE)

# ---- 2. Volcano ----
if (file.exists("results/tables/deseq2_T2D_vs_ND.csv")) {
  library(EnhancedVolcano)
  res <- read.csv("results/tables/deseq2_T2D_vs_ND.csv", row.names = 1)
  res$ensembl <- sub("\\..*$", "", rownames(res))
  ann <- fread("data/reference/gene_annotation.tsv", header = FALSE,
               col.names = c("gene_id","gene_name","gene_biotype"))
  ann$gene_id <- sub("\\..*$", "", ann$gene_id)
  res$label <- setNames(ann$gene_name, ann$gene_id)[res$ensembl]
  res$label[is.na(res$label) | res$label == ""] <- res$ensembl[is.na(res$label) | res$label == ""]
  sig <- which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
  topLab <- res$label[sig[order(res$padj[sig])][1:25]]
  png("results/figures/png/volcano.png", width = 2000, height = 1750, res = 220)
  print(EnhancedVolcano(res, lab = res$label, selectLab = topLab,
                        x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 1,
                        pointSize = 1.6, labSize = 3.6, drawConnectors = TRUE,
                        legendLabels = c("Not sig","Fold-change only","p-value only","Significant"),
                        title = "T2D vs ND", subtitle = NULL, legendPosition = "top",
                        gridlines.major = FALSE, gridlines.minor = FALSE))
  dev.off()
}

# ---- 3. Expression-supported ceRNA network (layered) ----
edges <- fread("results/tables/cerna_edges_supported.csv")
nodes <- fread("results/tables/cerna_hubs_supported.csv")
g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$class  <- nodes$class[match(V(g)$name, nodes$node)]
V(g)$degree <- nodes$degree[match(V(g)$name, nodes$node)]
pal  <- c(lncRNA = "#E64B35", miRNA = "#F0B429", mRNA = "#4DBBD5")
xpos <- c(lncRNA = 0, miRNA = 1, mRNA = 2)

lm_e <- edges[type == "lncRNA-miRNA"]
lm_e[, mir_deg := nodes$degree[match(to, nodes$node)]]
setorder(lm_e, from, -mir_deg)
sel_mir <- unique(lm_e[, head(.SD, 4), by = from]$to)
core <- edges[(type=="lncRNA-miRNA" & to %in% sel_mir) | (type=="miRNA-mRNA" & from %in% sel_mir)]
gc2 <- graph_from_data_frame(core, directed = FALSE)
V(gc2)$class  <- nodes$class[match(V(gc2)$name, nodes$node)]
V(gc2)$degree <- degree(gc2)
V(gc2)$color  <- pal[V(gc2)$class]
V(gc2)$size   <- pmin(pmax(3 + 1.6*sqrt(V(gc2)$degree), 7), 12)
V(gc2)$label  <- V(gc2)$name
lay <- matrix(0, vcount(gc2), 2)
for (cc in names(xpos)) {
  idx <- which(V(gc2)$class==cc); idx <- idx[order(V(gc2)$degree[idx], decreasing=TRUE)]
  lay[idx,1] <- xpos[[cc]]
  lay[idx,2] <- if (length(idx)==1) 0 else seq(-1,1,length.out=length(idx))
}
lab_deg <- ifelse(V(gc2)$class=="lncRNA", pi, 0)
lab_dst <- ifelse(V(gc2)$class=="miRNA", 0, 1.1)
png("results/figures/png/cerna_network.png", width = 2200, height = 2600, res = 200)
par(mar = c(5,1,3,1))
plot(gc2, layout=lay, rescale=TRUE, vertex.frame.color="grey45",
     vertex.label.cex=0.5, vertex.label.color="black",
     vertex.label.dist=lab_dst, vertex.label.degree=lab_deg,
     edge.color=adjustcolor("grey65", alpha.f=0.4), edge.width=0.5, edge.curved=0.1,
     main="Expression-supported ceRNA network")
legend("bottom", legend=names(pal), pt.bg=pal, pch=21, pt.cex=1.3, cex=0.8,
       bty="n", horiz=TRUE, inset=c(0,-0.04), xpd=TRUE)
dev.off()

# ---- 4. NEAT1 - miR-29 axis ----
mir29 <- grep("^hsa-miR-29", nodes$node, value = TRUE)
tgt   <- unique(edges[from %in% mir29 & type=="miRNA-mRNA"]$to)
gs <- induced_subgraph(g, which(V(g)$name %in% c("NEAT1", mir29, tgt)))
gs <- delete_vertices(gs, V(gs)[V(gs)$class=="lncRNA" & V(gs)$name!="NEAT1"])
V(gs)$color <- pal[V(gs)$class]
V(gs)$size  <- ifelse(V(gs)$class=="lncRNA",22, ifelse(V(gs)$class=="miRNA",15,10))
V(gs)$label <- V(gs)$name
set.seed(7)
png("results/figures/png/NEAT1_miR29_axis.png", width = 1800, height = 1700, res = 200)
par(mar = c(4,1,3,1))
plot(gs, layout=layout_with_fr(gs, niter=3000), vertex.frame.color="grey30",
     vertex.label.cex=0.65, vertex.label.color="black", vertex.label.dist=0,
     edge.color="grey70", edge.width=1, main="NEAT1 - miR-29 - ECM axis")
legend("bottom", legend=names(pal), pt.bg=pal, pch=21, pt.cex=1.3, cex=0.8,
       bty="n", horiz=TRUE, inset=c(0,-0.03), xpd=TRUE)
dev.off()

cat("PNGs written to results/figures/png/:\n")
print(list.files("results/figures/png"))