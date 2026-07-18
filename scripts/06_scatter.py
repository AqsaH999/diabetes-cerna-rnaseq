#!/usr/bin/env python3
"""Validation scatter: our counts vs published counts for one donor.
Saves results/figures/validation_scatter.png"""
import re, math
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REP_SRR, REP_DP = "SRR13380421", "DP003"      # representative donor
OUR = "results/tables/subset_counts.txt"
PUB = "data/processed/GSE164416_counts.txt"
strip = lambda g: g.split(".")[0]

# our counts for the representative sample
our = {}
with open(OUR) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if p[0] == "Geneid":
            col = next(i for i in range(6, len(p)) if REP_SRR in p[i]); continue
        our[strip(p[0])] = int(p[col])

# published counts for the matched donor
pub = {}
with open(PUB) as f:
    hdr = f.readline().rstrip("\n").split("\t"); j = hdr.index(REP_DP)
    for line in f:
        p = line.rstrip("\n").split("\t"); pub[strip(p[0])] = int(p[j])

genes = [g for g in our if g in pub]
x = [math.log2(pub[g]+1) for g in genes]   # published
y = [math.log2(our[g]+1) for g in genes]   # ours

# Pearson r
n=len(x); mx=sum(x)/n; my=sum(y)/n
r=sum((a-mx)*(b-my) for a,b in zip(x,y))/(
   math.sqrt(sum((a-mx)**2 for a in x))*math.sqrt(sum((b-my)**2 for b in y)))

plt.figure(figsize=(6,6))
plt.hexbin(x, y, gridsize=60, bins="log", cmap="viridis")
lim=max(max(x),max(y))
plt.plot([0,lim],[0,lim], "r--", lw=1, label="y = x")
plt.xlabel("Published counts  log2(count+1)")
plt.ylabel("Our pipeline  log2(count+1)")
plt.title(f"Track A vs Track B — {REP_DP}  (Pearson r = {r:.3f})")
plt.legend(loc="upper left"); plt.colorbar(label="genes per bin (log)")
plt.tight_layout()
import os; os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/validation_scatter.png", dpi=150)
print(f"Saved results/figures/validation_scatter.png  (r={r:.3f}, {len(genes)} genes)")