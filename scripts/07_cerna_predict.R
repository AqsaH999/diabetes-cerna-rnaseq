# ==========================================================
# 07_cerna_predict.R — Part 1: miRNA -> mRNA edges (multiMiR)
# ==========================================================
setwd("~/Desktop/diabetes-cerna-rnaseq")
library(multiMiR); library(data.table)

de_mrna <- fread("results/tables/DE_mRNA.csv")
syms <- unique(de_mrna[!is.na(gene_name) & gene_name != ""]$gene_name)
cat("Querying multiMiR for", length(syms), "DE mRNAs (validated interactions)...\n")

chunks   <- split(syms, ceiling(seq_along(syms) / 50))
res_list <- list()

for (i in seq_along(chunks)) {
  cat(sprintf("  chunk %d/%d (%d genes)... ", i, length(chunks), length(chunks[[i]])))
  r <- tryCatch({
    q <- get_multimir(org = "hsa", target = chunks[[i]],
                      table = "validated", summary = FALSE)
    as.data.table(q@data)
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r) && nrow(r) > 0) {
    res_list[[i]] <- r; cat(nrow(r), "interactions\n")
  } else cat("none\n")
}

mm <- rbindlist(res_list, fill = TRUE)
mir_mrna <- unique(mm[, .(mature_mirna_id, target_symbol)])
mir_mrna <- mir_mrna[!is.na(mature_mirna_id) & mature_mirna_id != ""]

fwrite(mir_mrna, "results/tables/miRNA_mRNA_edges.csv")
cat("\nSaved miRNA_mRNA_edges.csv:", nrow(mir_mrna), "unique edges |",
    uniqueN(mir_mrna$mature_mirna_id), "distinct miRNAs |",
    uniqueN(mir_mrna$target_symbol), "distinct mRNAs\n")



# ==========================================================
# 07_cerna_predict.R — Part 2: lncRNA -> miRNA edges (ENCORI)
# ==========================================================
library(data.table)

de_lnc   <- fread("results/tables/DE_lncRNA.csv")
lnc_syms <- unique(de_lnc[!is.na(gene_name) & gene_name != ""]$gene_name)
cat("Querying ENCORI for", length(lnc_syms), "named DE lncRNAs...\n")

encori_lnc <- function(gene, clip = 1) {
  url <- paste0("https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38",
                "&geneType=lncRNA&miRNA=all&clipExpNum=", clip,
                "&degraExpNum=0&pancancerNum=0&programNum=0&program=None",
                "&target=", utils::URLencode(gene, reserved = TRUE), "&cellType=all")
  txt <- tryCatch(system(sprintf("curl -s --max-time 90 '%s'", url), intern = TRUE),
                  error = function(e) NULL)
  if (is.null(txt) || length(txt) == 0) return(NULL)
  txt <- txt[!grepl("^#", txt)]                       # drop citation lines
  if (length(txt) < 2) return(NULL)                   # header only = no hits
  dt <- tryCatch(fread(text = paste(txt, collapse = "\n"), sep = "\t"),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  need <- c("geneName", "miRNAname", "clipExpNum")
  if (!all(need %in% names(dt))) return(NULL)         # skip malformed responses
  unique(dt[, .(lncRNA = geneName, mature_mirna_id = miRNAname, clipExpNum)])
}

lnc_list <- list()
for (i in seq_along(lnc_syms)) {
  g <- lnc_syms[i]
  cat(sprintf("  [%2d/%d] %-15s ... ", i, length(lnc_syms), g))
  r <- tryCatch(encori_lnc(g), error = function(e) NULL)   # belt-and-braces
  if (!is.null(r) && nrow(r) > 0) { lnc_list[[g]] <- r; cat(nrow(r), "sites\n") }
  else cat("none\n")
  Sys.sleep(0.4)                                      # be polite to the server
}

lnc_mir <- rbindlist(lnc_list, fill = TRUE)

if (nrow(lnc_mir) == 0) {
  cat("\nNo lncRNA-miRNA interactions returned. Check connectivity.\n")
} else {
  fwrite(lnc_mir, "results/tables/lncRNA_miRNA_sites_raw.csv")
  
  # collapse multiple binding sites into one edge per lncRNA-miRNA pair
  lnc_edges <- lnc_mir[, .(clipExpNum = max(clipExpNum)), by = .(lncRNA, mature_mirna_id)]
  fwrite(lnc_edges, "results/tables/lncRNA_miRNA_edges.csv")
  
  cat("\nSaved lncRNA_miRNA_edges.csv:", nrow(lnc_edges), "edges |",
      uniqueN(lnc_edges$lncRNA), "lncRNAs |",
      uniqueN(lnc_edges$mature_mirna_id), "miRNAs\n")
  
  cat("\nEdges contributed per lncRNA:\n")
  print(lnc_edges[, .N, by = lncRNA][order(-N)])
  
  cat("\nEvidence spread (clipExpNum):\n")
  print(summary(lnc_edges$clipExpNum))
}