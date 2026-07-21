# Methods

## Dataset
- GEO: GSE164416 (Wigger et al., Nature Metabolism 2021, PMID 34183850)
- SRA: PRJNA690574. Human pancreatic islets, laser-capture microdissected.
- Track A (pipeline demonstration): 12 samples (6 ND, 6 T2D), single-end SMART-Seq, 76 bp.
- Track B (statistics): published HTSeq count matrix, 18 ND vs 39 T2D.

## Reference
- Genome: GRCh38 (prebuilt HISAT2 index)
- Annotation: Ensembl release 110 (chosen over GENCODE v44 to match index chromosome naming)

## Tool versions
- FastQC 0.12.1
- fastp 1.3.6
- HISAT2 2.2.2
- samtools 1.x
- featureCounts (Subread) 2.1.1
- MultiQC 1.35
- sra-tools (fasterq-dump) 3.4.1
- DESeq2, clusterProfiler 4.20.0, multiMiR 2.4.0 (R/Bioconductor)

## Key parameters
- Trimming (fastp): single-end, --qualified_quality_phred 20 --length_required 25
- Alignment (HISAT2): single-end -U, -p 6
- Counting (featureCounts): -s 0 unstranded, single-end (no -p)
- DE (DESeq2): ~ condition, ND reference; padj < 0.05 and |log2FC| > 1
- Enrichment: GO BP, background universe = tested genes only; simplify(cutoff = 0.7)
- miRNA-mRNA: multiMiR, miRTarBase Functional MTI (gold standard) only
- lncRNA-miRNA: ENCORI API (hg38); miRcode is defunct
- ceRNA filtering: Spearman lncRNA-mRNA correlation across 57 donors, rho > 0 and padj < 0.05
- Robustness: correlations recomputed within T2D donors alone (80% positive, p = 1e-13; agreement rho = 0.93 with all-donor), confirming within-disease co-variation rather than group separation

## Key results
- Alignment 91.1-97.5%; featureCounts assignment ~48% (typical for total-RNA SMART-Seq)
- Pipeline validation vs published counts: mean Pearson r = 0.9746 (12 donors)
- 726 DE genes: 515 mRNAs, 125 lncRNAs
- Database ceRNA network: 208 nodes, 401 edges, 336 axes
- Expression-supported: 180 nodes, 335 edges, 238 axes; 92% of lncRNA-mRNA pairs positively correlated (p = 5e-28)
- Hub rankings unchanged after filtering (miR-29 family; IGF1/MMP2/VIM/CCN2)

## Limitations
- Mature miRNAs are not captured by standard RNA-seq; the miRNA layer is database-inferred, not measured.
- Bulk LCM tissue: immune signal partly reflects cell-composition differences between donors.
- 59 of 125 DE lncRNAs lack gene symbols and could not be queried against ENCORI.
- Positive lncRNA-mRNA correlation is consistent with ceRNA activity but also with shared upstream regulation.
