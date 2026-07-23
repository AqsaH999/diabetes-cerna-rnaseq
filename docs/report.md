# ceRNA Network Reconstruction in Type 2 Diabetic Pancreatic Islets

### Reproducible RNA-seq pipeline (raw FASTQ → validated counts → lncRNA–miRNA–mRNA network) on GSE164416, uncovering a NEAT1–miR-29–collagen axis

**Dataset:** GSE164416 (Wigger et al., *Nature Metabolism* 2021) - human pancreatic islets,
laser-capture microdissected, non-diabetic (ND) vs type-2-diabetic (T2D) donors.

---

## 1. Question

Which genes behave differently in diabetic pancreatic islets, and can a regulatory
lncRNA -> miRNA -> mRNA (competing-endogenous-RNA, ceRNA) network be reconstructed that
plausibly underlies those changes? The original study did not perform a ceRNA analysis, so
this project **reproduces** its primary processing and **extends** it with the ceRNA layer.

## 2. Approach

A two-track design was used. **Track A** downloaded 6 ND + 6 T2D raw FASTQ samples and ran the
full pipeline (FastQC -> fastp -> HISAT2 -> featureCounts) to demonstrate command-line
competence. **Track B** ran the biology on the full published count matrix (18 ND vs 39 T2D)
for statistical power. The two were tied together by a validation step. Differential expression
(DESeq2) was followed by biotype splitting, functional enrichment (clusterProfiler), miRNA
prediction (multiMiR + ENCORI), network assembly (igraph), and an expression-correlation
filter.

## 3. Key findings

**Pipeline validation.** Track-A counts correlated with the published counts at
**mean Pearson r = 0.975** across 12 matched donors, confirming the from-scratch pipeline
reproduces the published processing despite different tools (HISAT2/featureCounts vs the
authors' STAR/htseq-count).

**Differential expression.** 726 genes were significant (padj < 0.05, |log2FC| > 1): 515
protein-coding mRNAs and 125 lncRNAs.

**Functional enrichment (two-sided story).**
- *Up in T2D:* immune activation, leukocyte chemotaxis, MHC-II antigen presentation, and
  extracellular-matrix (ECM) / collagen organisation - active immune recruitment plus fibrotic
  remodelling.
- *Down in T2D:* ribosome biogenesis, translation, Golgi vesicle transport, and the TCA cycle -
  the biosynthetic and secretory machinery beta cells require for glucose sensing and insulin
  production.

**ceRNA network.** A database-derived network (multiMiR gold-standard + ENCORI CLIP evidence)
of 208 nodes / 336 axes was filtered against the expression data: **92% of lncRNA-mRNA pairs
were positively correlated** (p = 5e-28), retaining 238 axes. A within-T2D robustness check
(removing the ND/T2D group effect) kept 80% positive (p = 1e-13; agreement rho = 0.93),
indicating within-disease co-variation rather than group separation alone. Hub rankings were
unchanged by filtering.

**Central axis: NEAT1 -> miR-29 -> ECM.** The dominant, filter-robust axis is the up-regulated
lncRNA **NEAT1** sponging the **miR-29 family**, whose validated targets are collagens and ECM
genes (COL1A1, COL3A1, COL4A2, SERPINH1, ELN, MMP2, PDGFRB).

## 4. Literature sanity-checks

Three headline genes were checked against independent published work:

- **SLC2A2 / GLUT2 (down in my data):** SLC2A2 mRNA is reduced by 70-90% in human T2D islets,
  part of documented beta-cell dedifferentiation (loss of GLUT2, PDX1, MafA). My result matches
  this well-established finding [Diabetologia 2014; SAGE 2023].
- **NEAT1 (up, central sponge):** NEAT1 is an established driver of diabetic inflammation
  (NLRP3/IL-1beta) and of ECM accumulation and fibrosis in diabetic complications, consistent
  with its position here [Springer 2026 review].
- **miR-29 (predicted functional reduction):** In profibrotic environments miR-29 is suppressed
  (via TGF-beta), which de-represses collagen synthesis and drives ECM accumulation - exactly
  the direction my NEAT1-sponge model predicts [Am J Physiol Cell Physiol 2022].

**Honest nuance on miR-29.** Circulating/bulk miR-29 is often reported *increased* in T2D and
metabolic disease, and miR-29a can itself impair glucose-stimulated insulin secretion. My model
proposes a *local functional reduction* of miR-29 through NEAT1 sponging, consistent with the
pro-fibrotic "miR-29-low -> collagen-high" mechanism - but because mature miRNAs were not
measured here, this remains an inference to be tested directly.

## 5. Limitations

- Mature miRNAs are not captured by standard RNA-seq; the miRNA layer is database-inferred.
- Bulk laser-captured tissue means the immune signal partly reflects cell-composition
  differences between donors rather than beta-cell-intrinsic change.
- Positive lncRNA-mRNA correlation supports, but does not prove, sponging (shared upstream
  regulation gives the same pattern).
- 59 of 125 DE lncRNAs lack gene symbols and could not be queried for miRNA interactions.

## 6. Conclusion and next steps

Diabetic islets show a coherent shift - immune infiltration and fibrotic ECM remodelling up,
beta-cell biosynthetic machinery down - and the reconstructed ceRNA network places an
experimentally-supported **NEAT1 / miR-29 / collagen** axis at the centre of the fibrotic arm,
in agreement with three independent lines of evidence (differential expression, enrichment,
expression correlation) and with prior literature. The natural next steps are **matched
small-RNA-seq** to measure the miRNA layer directly and **qRT-PCR / luciferase reporter
validation** of the NEAT1-miR-29 axis.

## References

- Wigger et al. *Nature Metabolism* 2021. PMID 34183850.
- GLUT2, glucose sensing and glucose homeostasis. *Diabetologia* 2014.
- Beta-cell dysfunction in type 2 diabetes. *SAGE / J Diabetes Res* 2023.
- lncRNAs in type 2 diabetes (review). *Discover Endocrinology and Metabolism* 2026.
- The microRNA-29 family: role in metabolism and metabolic disease. *Am J Physiol Cell Physiol*
  2022. PMID 35704699.
