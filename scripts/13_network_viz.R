# ==========================================================
# 09_network_viz.R — ceRNA network figures + hub analysis
# (run AFTER 08_network.R)
# ==========================================================
library(data.table); library(igraph)

edges <- fread("results/tables/cerna_edges.csv")
nodes <- fread("results/tables/cerna_node_degree.csv")

g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$class  <- nodes$class[match(V(g)$name, nodes$node)]
V(g)$degree <- nodes$degree[match(V(g)$name, nodes$node)]

pal  <- c(lncRNA = "#E64B35", miRNA = "#F0B429", mRNA = "#4DBBD5")
xpos <- c(lncRNA = 0, miRNA = 1, mRNA = 2)
V(g)$color <- pal[V(g)$class]
V(g)$size  <- pmin(2 + 1.8 * sqrt(V(g)$degree), 15)

top_mir  <- head(nodes[class == "miRNA"][order(-degree)]$node, 8)
top_mrna <- head(nodes[class == "mRNA"][order(-degree)]$node, 8)
keep_lab <- c(nodes[class == "lncRNA"]$node, top_mir, top_mrna)
V(g)$label <- ifelse(V(g)$name %in% keep_lab, V(g)$name, NA)

# legend drawn OUTSIDE the plot, in the bottom margin
add_legend <- function() {
  legend("bottom", legend = names(pal), pt.bg = pal, pch = 21,
         pt.cex = 1.3, cex = 0.8, bty = "n", horiz = TRUE,
         inset = c(0, -0.05), xpd = TRUE)
}

# ---------- Figure 1a: full network, force-directed ----------
set.seed(42)
lay <- layout_with_fr(g, niter = 3000)

pdf("results/figures/cerna_network.pdf", width = 14, height = 14)
par(mar = c(6, 1, 3, 1))
plot(g, layout = lay, rescale = TRUE,
     vertex.frame.color = "grey45",
     vertex.label.cex = 0.65, vertex.label.color = "black",
     vertex.label.font = 2,
     vertex.label.dist = 1.4, vertex.label.degree = -pi/2,
     edge.color = adjustcolor("grey70", alpha.f = 0.45), edge.width = 0.5,
     main = "T2D vs ND ceRNA network  (lncRNA - miRNA - mRNA)")
add_legend()
mtext(sprintf("%d nodes, %d edges | labels: all lncRNAs + top 8 miRNA/mRNA hubs",
              vcount(g), ecount(g)), side = 1, line = 3, cex = 0.85)
dev.off()

# ---------- Figure 1b: layered view, ALL 10 lncRNAs represented ----------
lm_e <- edges[type == "lncRNA-miRNA"]
lm_e[, mir_deg := nodes$degree[match(to, nodes$node)]]
setorder(lm_e, from, -mir_deg)
sel_mir <- unique(lm_e[, head(.SD, 4), by = from]$to)     # top 4 miRNAs per lncRNA

core_edges <- edges[(type == "lncRNA-miRNA" & to   %in% sel_mir) |
                      (type == "miRNA-mRNA"   & from %in% sel_mir)]

gc2 <- graph_from_data_frame(core_edges, directed = FALSE)
V(gc2)$class  <- nodes$class[match(V(gc2)$name, nodes$node)]
V(gc2)$degree <- degree(gc2)
V(gc2)$color  <- pal[V(gc2)$class]
V(gc2)$size   <- pmin(pmax(3 + 1.6 * sqrt(V(gc2)$degree), 7), 12)
V(gc2)$label  <- V(gc2)$name

lay3 <- matrix(0, nrow = vcount(gc2), ncol = 2)
for (cc in names(xpos)) {
  idx <- which(V(gc2)$class == cc)
  idx <- idx[order(V(gc2)$degree[idx], decreasing = TRUE)]
  lay3[idx, 1] <- xpos[[cc]]
  lay3[idx, 2] <- if (length(idx) == 1) 0 else seq(-1, 1, length.out = length(idx))
}

# miRNA labels centred ON the node; lncRNA labels left; mRNA labels right
lab_deg <- ifelse(V(gc2)$class == "lncRNA", pi, 0)
lab_dst <- ifelse(V(gc2)$class == "miRNA", 0, 1.4)

pdf("results/figures/cerna_network_layered.pdf", width = 12, height = 14)
par(mar = c(6, 1, 3, 1))
plot(gc2, layout = lay3, rescale = TRUE,
     vertex.frame.color = "grey45",
     vertex.label.cex = 0.5, vertex.label.color = "black",
     vertex.label.dist = lab_dst, vertex.label.degree = lab_deg,
     edge.color = adjustcolor("grey65", alpha.f = 0.4), edge.width = 0.5,
     edge.curved = 0.1,
     main = "ceRNA network:  all 10 DE lncRNAs  ->  miRNA  ->  mRNA")
add_legend()
mtext(sprintf("%d nodes, %d edges | top 4 miRNA partners shown per lncRNA",
              vcount(gc2), ecount(gc2)), side = 1, line = 3, cex = 0.85)
dev.off()
cat("lncRNAs shown in layered figure:", sum(V(gc2)$class == "lncRNA"), "of 10\n")

# ---------- Figure 2: focused NEAT1 - miR-29 axis ----------
mir29     <- grep("^hsa-miR-29", nodes$node, value = TRUE)
mir29_tgt <- unique(edges[from %in% mir29 & type == "miRNA-mRNA"]$to)
keep      <- unique(c("NEAT1", mir29, mir29_tgt))
gs        <- induced_subgraph(g, vids = which(V(g)$name %in% keep))
gs        <- delete_vertices(gs, V(gs)[V(gs)$class == "lncRNA" & V(gs)$name != "NEAT1"])

V(gs)$size  <- ifelse(V(gs)$class == "lncRNA", 22,
                      ifelse(V(gs)$class == "miRNA", 15, 10))
V(gs)$label <- V(gs)$name

set.seed(7)
pdf("results/figures/cerna_NEAT1_miR29_axis.pdf", width = 9, height = 8.5)
par(mar = c(5, 1, 3, 1))
plot(gs, layout = layout_with_fr(gs, niter = 3000),
     vertex.frame.color = "grey30",
     vertex.label.cex = 0.65, vertex.label.color = "black",
     vertex.label.dist = 0,
     edge.color = "grey70", edge.width = 1,
     main = "NEAT1 - miR-29 family - ECM target axis")
add_legend()
dev.off()

# ---------- Hub analysis table ----------
nodes[, betweenness := betweenness(g, v = node, directed = FALSE)]
setorder(nodes, -degree)
fwrite(nodes, "results/tables/cerna_hubs.csv")

cat("\nSaved: cerna_network.pdf, cerna_network_layered.pdf,",
    "cerna_NEAT1_miR29_axis.pdf, cerna_hubs.csv\n\n")
cat("--- Top 10 hubs (degree) ---\n");   print(head(nodes, 10))
cat("\n--- Top 10 (betweenness) ---\n"); print(head(nodes[order(-betweenness)], 10))
cat("\nNEAT1-miR-29 axis:", vcount(gs), "nodes |", ecount(gs), "edges\n")