#!/usr/bin/env python3
"""Validate Track-A counts (ours) against the published matrix (Track B).
For each donor, correlate per-gene log2 counts. Expect r >= ~0.95."""
import re, math
# SRR -> donor (DP), verified from GEO sample titles
srr2dp = {
 "SRR13380421":"DP003","SRR13380425":"DP010","SRR13380426":"DP011",
 "SRR13380427":"DP012","SRR13380430":"DP019","SRR13380434":"DP030",
 "SRR13380441":"DP049","SRR13380446":"DP055","SRR13380448":"DP058",
 "SRR13380458":"DP070","SRR13380466":"DP079","SRR13380471":"DP086",
}
OUR = "results/tables/subset_counts.txt"
PUB = "data/processed/GSE164416_counts.txt"
strip = lambda g: g.split(".")[0]

# --- read our count matrix ---
our = {}; our_cols = []
with open(OUR) as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if p[0] == "Geneid":                       # header row
            for i in range(6, len(p)):             # sample cols start at 7th (index 6)
                m = re.search(r"(SRR\d+)", p[i])
                if m: our_cols.append((i, m.group(1)))
            continue
        our[strip(p[0])] = {srr: int(p[i]) for i, srr in our_cols}
# --- read published matrix ---
pub = {}
with open(PUB) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    dp_idx = {hdr[i]: i for i in range(1, len(hdr))}
    need = set(srr2dp.values())
    for line in f:
        p = line.rstrip("\n").split("\t")
        pub[strip(p[0])] = {dp: int(p[dp_idx[dp]]) for dp in need if dp in dp_idx}

def pearson(x, y):
    n=len(x); mx=sum(x)/n; my=sum(y)/n
    cov=sum((a-mx)*(b-my) for a,b in zip(x,y))
    sx=math.sqrt(sum((a-mx)**2 for a in x)); sy=math.sqrt(sum((b-my)**2 for b in y))
    return cov/(sx*sy) if sx and sy else float("nan")
genes = [g for g in our if g in pub]
print(f"Shared genes compared: {len(genes)}\n")
print(f"{'Sample':15}{'Donor':7}{'Pearson r':>11}")
rs=[]
for srr, dp in sorted(srr2dp.items()):
    if dp not in dp_idx:
        print(f"{srr:15}{dp:7}   (not in published)"); continue
    xs=[]; ys=[]
    for g in genes:
        xs.append(math.log2(our[g][srr]+1)); ys.append(math.log2(pub[g][dp]+1))
    r=pearson(xs,ys); rs.append(r)
    print(f"{srr:15}{dp:7}{r:11.4f}")
print(f"\nMean Pearson r across {len(rs)} samples: {sum(rs)/len(rs):.4f}")
