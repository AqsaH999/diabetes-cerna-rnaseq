# ==========================================================
# 11_network_supported.R — network + figures from correlation-supported axes
# ==========================================================
setwd("~/Desktop/diabetes-cerna-rnaseq")   # remove during Day 12 cleanup
library(data.table); library(igraph)

tri <- fread("results/tables/cerna_triads_supported.csv")

# ---- rebuild edges from the supported axes ----
edges <- unique(rbind(
  tri[, .(from = lncRNA,          to = mature_mirna_id, type = "lncRNA-miRNA")],
  tri[, .(from = mature_mirna_id, to = target_symbol,   type = "miRNA-mRNA")]
))
fwrite(edges, "results/tables/cerna_edges_supported.csv")

g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$class  <- ifelse(V(g)$name %in% tri$lncRNA, "lncRNA",
                      ifelse(grepl("^hsa-", V(g)$name), "miRNA", "mRNA"))
V(g)$degree <- degree(g)

nodes <- data.table(node = V(g)$name, class = V(g)$class, degree = V(g)$degree)
nodes[, betweenness := betweenness(g, v = node, directed = FALSE)]
setorder(nodes, -degree)
fwrite(nodes, "results/tables/cerna_hubs_supported.csv")

pal  <- c(lncRNA = "#E64B35", miRNA = "#F0B429", mRNA = "#4DBBD5")
xpos <- c(lncRNA = 0, miRNA = 1, mRNA = 2)
V(g)$color <- pal[V(g)$class]

add_legend <- function() {
  legend("bottom", legend = names(pal), pt.bg = pal, pch = 21,
         pt.cex = 1.3, cex = 0.8, bty = "n", horiz = TRUE,
         inset = c(0, -0.05), xpd = TRUE)
}

# ---------- Figure A: layered, all supported lncRNAs ----------
lm_e <- edges[type == "lncRNA-miRNA"]
lm_e[, mir_deg := nodes$degree[match(to, nodes$node)]]
setorder(lm_e, from, -mir_deg)
sel_mir <- unique(lm_e[, head(.SD, 4), by = from]$to)

core <- edges[(type == "lncRNA-miRNA" & to   %in% sel_mir) |
                (type == "miRNA-mRNA"   & from %in% sel_mir)]

gc2 <- graph_from_data_frame(core, directed = FALSE)
V(gc2)$class  <- nodes$class[match(V(gc2)$name, nodes$node)]
V(gc2)$degree <- degree(gc2)
V(gc2)$color  <- pal[V(gc2)$class]
V(gc2)$size   <- pmin(pmax(3 + 1.6 * sqrt(V(gc2)$degree), 7), 12)
V(gc2)$label  <- V(gc2)$name

lay <- matrix(0, nrow = vcount(gc2), ncol = 2)
for (cc in names(xpos)) {
  idx <- which(V(gc2)$class == cc)
  idx <- idx[order(V(gc2)$degree[idx], decreasing = TRUE)]
  lay[idx, 1] <- xpos[[cc]]
  lay[idx, 2] <- if (length(idx) == 1) 0 else seq(-1, 1, length.out = length(idx))
}
lab_deg <- ifelse(V(gc2)$class == "lncRNA", pi, 0)
lab_dst <- ifelse(V(gc2)$class == "miRNA", 0, 1.1)

pdf("results/figures/cerna_network_supported.pdf", width = 12, height = 14)
par(mar = c(6, 1, 3, 1))
plot(gc2, layout = lay, rescale = TRUE,
     vertex.frame.color = "grey45",
     vertex.label.cex = 0.5, vertex.label.color = "black",
     vertex.label.dist = lab_dst, vertex.label.degree = lab_deg,
     edge.color = adjustcolor("grey65", alpha.f = 0.4), edge.width = 0.5,
     edge.curved = 0.1,
     main = "Expression-supported ceRNA network:  lncRNA  ->  miRNA  ->  mRNA")
add_legend()
mtext(sprintf("%d nodes, %d edges | axes retained after lncRNA-mRNA correlation filtering (rho>0, padj<0.05)",
              vcount(gc2), ecount(gc2)), side = 1, line = 3, cex = 0.8)
dev.off()

# ---------- Figure B: NEAT1 - miR-29 axis (supported edges only) ----------
mir29     <- grep("^hsa-miR-29", nodes$node, value = TRUE)
mir29_tgt <- unique(edges[from %in% mir29 & type == "miRNA-mRNA"]$to)
keep      <- unique(c("NEAT1", mir29, mir29_tgt))
gs <- induced_subgraph(g, vids = which(V(g)$name %in% keep))
gs <- delete_vertices(gs, V(gs)[V(gs)$class == "lncRNA" & V(gs)$name != "NEAT1"])

V(gs)$color <- pal[V(gs)$class]
V(gs)$size  <- ifelse(V(gs)$class == "lncRNA", 22,
                      ifelse(V(gs)$class == "miRNA", 15, 10))
V(gs)$label <- V(gs)$name

set.seed(7)
pdf("results/figures/cerna_NEAT1_miR29_supported.pdf", width = 9, height = 8.5)
par(mar = c(5, 1, 3, 1))
plot(gs, layout = layout_with_fr(gs, niter = 3000),
     vertex.frame.color = "grey30",
     vertex.label.cex = 0.65, vertex.label.color = "black",
     vertex.label.dist = 0,
     edge.color = "grey70", edge.width = 1,
     main = "NEAT1 - miR-29 family - ECM targets (expression-supported)")
add_legend()
dev.off()

# ---------- Report ----------
cat("=== EXPRESSION-SUPPORTED ceRNA NETWORK ===\n")
cat("Nodes:", vcount(g), "| Edges:", ecount(g), "| Axes:", nrow(tri), "\n")
print(table(nodes$class))
cat("\n--- lncRNAs ---\n");            print(nodes[class == "lncRNA"])
cat("\n--- Top 10 miRNA hubs ---\n");  print(head(nodes[class == "miRNA"], 10))
cat("\n--- Top 10 mRNA hubs ---\n");   print(head(nodes[class == "mRNA"], 10))
cat("\nSaved: cerna_network_supported.pdf, cerna_NEAT1_miR29_supported.pdf,",
    "cerna_edges_supported.csv, cerna_hubs_supported.csv\n")