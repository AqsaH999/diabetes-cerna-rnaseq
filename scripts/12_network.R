# ==========================================================
# 08_network.R — assemble the lncRNA-miRNA-mRNA ceRNA network
# ==========================================================
library(data.table); library(igraph)

lm <- fread("results/tables/lncRNA_miRNA_edges_final.csv")   # lncRNA -- miRNA
mm <- fread("results/tables/miRNA_mRNA_edges_final.csv")     # miRNA  -- mRNA

# fix symbol-casing duplicates (miRTarBase mixes e.g. "FAS" and "Fas")
mm[, target_symbol := toupper(target_symbol)]
mm <- unique(mm)

# --- 1. Unified edge list with layer labels ---
edges <- rbind(
  lm[, .(from = lncRNA,          to = mature_mirna_id, type = "lncRNA-miRNA")],
  mm[, .(from = mature_mirna_id, to = target_symbol,   type = "miRNA-mRNA")]
)
edges <- unique(edges)
fwrite(edges, "results/tables/cerna_edges.csv")

# --- 2. Build the graph and label node classes ---
g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$class <- ifelse(V(g)$name %in% lm$lncRNA, "lncRNA",
                     ifelse(grepl("^hsa-", V(g)$name), "miRNA", "mRNA"))
V(g)$degree <- degree(g)

nodes <- data.table(node = V(g)$name, class = V(g)$class, degree = V(g)$degree)[order(-degree)]
fwrite(nodes, "results/tables/cerna_node_degree.csv")

# --- 3. Complete ceRNA axes (lncRNA -> miRNA -> mRNA triads) ---
triads <- merge(lm[, .(lncRNA, mature_mirna_id)],
                mm[, .(mature_mirna_id, target_symbol)],
                by = "mature_mirna_id", allow.cartesian = TRUE)
setcolorder(triads, c("lncRNA", "mature_mirna_id", "target_symbol"))
fwrite(triads, "results/tables/cerna_triads.csv")

# --- 4. Report ---
comp <- components(g)
cat("=== ceRNA NETWORK ===\n")
cat("Nodes:", vcount(g), "| Edges:", ecount(g), "\n")
print(table(nodes$class))
cat("\nConnected components:", comp$no, "| largest component:", max(comp$csize), "nodes\n")
cat("Complete lncRNA-miRNA-mRNA axes:", nrow(triads), "\n")

cat("\n--- lncRNA hubs (all) ---\n"); print(nodes[class == "lncRNA"])
cat("\n--- Top 15 miRNA hubs ---\n"); print(head(nodes[class == "miRNA"], 15))
cat("\n--- Top 15 mRNA hubs ---\n");  print(head(nodes[class == "mRNA"], 15))