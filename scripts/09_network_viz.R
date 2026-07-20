# ==========================================================
# 09_network_viz.R — ceRNA network figures + hub analysis
# ==========================================================
setwd("~/Desktop/diabetes-cerna-rnaseq")
library(data.table); library(igraph)

edges <- fread("results/tables/cerna_edges.csv")
nodes <- fread("results/tables/cerna_node_degree.csv")

g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$class  <- nodes$class[match(V(g)$name, nodes$node)]
V(g)$degree <- nodes$degree[match(V(g)$name, nodes$node)]

pal <- c(lncRNA = "#E64B35", miRNA = "#F0B429", mRNA = "#4DBBD5")
V(g)$color <- pal[V(g)$class]
V(g)$size  <- 2.5 + 2.2 * sqrt(V(g)$degree)
V(g)$label <- ifelse(V(g)$degree >= 5 | V(g)$class == "lncRNA", V(g)$name, NA)

# ---------- Figure 1: full network ----------
set.seed(42)
lay <- layout_with_fr(g)

pdf("results/figures/cerna_network.pdf", width = 12, height = 12)
par(mar = c(1, 1, 3, 1))
plot(g, layout = lay,
     vertex.frame.color = "grey40",
     vertex.label.cex = 0.5, vertex.label.color = "black",
     vertex.label.dist = 0.35, vertex.label.family = "sans",
     edge.color = "grey85", edge.width = 0.6,
     main = "T2D vs ND ceRNA network  (lncRNA - miRNA - mRNA)")
legend("bottomleft", legend = names(pal), pt.bg = pal, pch = 21,
       pt.cex = 2, bty = "n", title = "Node class")
mtext(sprintf("%d nodes, %d edges | labels shown for hubs (degree >= 5) and all lncRNAs",
              vcount(g), ecount(g)), side = 1, cex = 0.8)
dev.off()

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
pdf("results/figures/cerna_NEAT1_miR29_axis.pdf", width = 9, height = 8)
par(mar = c(1, 1, 3, 1))
plot(gs, layout = layout_with_fr(gs),
     vertex.frame.color = "grey30",
     vertex.label.cex = 0.65, vertex.label.color = "black",
     vertex.label.family = "sans",
     edge.color = "grey70", edge.width = 1,
     main = "NEAT1 - miR-29 family - ECM target axis")
legend("bottomleft", legend = names(pal), pt.bg = pal, pch = 21,
       pt.cex = 2, bty = "n", title = "Node class")
dev.off()

# ---------- Hub analysis table ----------
nodes[, betweenness := betweenness(g, v = node, directed = FALSE)]
setorder(nodes, -degree)
fwrite(nodes, "results/tables/cerna_hubs.csv")

cat("Saved: cerna_network.pdf, cerna_NEAT1_miR29_axis.pdf, cerna_hubs.csv\n\n")
cat("--- Top 10 hubs overall (degree) ---\n"); print(head(nodes, 10))
cat("\n--- Top 10 by betweenness (bridging nodes) ---\n")
print(head(nodes[order(-betweenness)], 10))
cat("\nNEAT1-miR-29 axis:", vcount(gs), "nodes |", ecount(gs), "edges\n")
cat("miR-29 family targets:\n"); print(sort(mir29_tgt))