# Biology Background: What This Project Is Really Doing

A plain-English companion to the 15-day plan. No fancy words — where a real term is needed, it's explained right away. Reread this anytime you feel lost about *why* a step exists.

---

## The big picture

Your cells contain DNA — the full instruction manual for the body. But DNA just sits in the nucleus; it doesn't do anything by itself. When a cell wants to *use* an instruction, it makes a working copy of that gene. That copy is called **RNA**. So the flow is:

**DNA → RNA → protein** (proteins are the machines that do the work in a cell).

The amount of RNA a cell makes for each gene tells you how "switched on" that gene is right now. If a diabetic person's cells make a lot more RNA for gene X than a healthy person's cells, that's a clue that gene X is involved in diabetes. **Measuring all the RNA in a sample is called RNA-seq.** That is the heart of this whole project.

The basic question the project asks:

> *"Which genes behave differently in diabetic islet cells compared to healthy ones — and can I map the hidden control system making them behave that way?"*

---

## The three types of RNA you care about

Not all RNA becomes protein. Three kinds matter here:

- **mRNA (messenger RNA)** — the "normal" kind. The working copy of a gene that gets turned into a protein. When people say "this gene is expressed," they usually mean its mRNA. These are the main **workers**.
- **lncRNA (long non-coding RNA)** — long RNA that does *not* become a protein. Many of them are regulators: they help control which other genes are on or off. Think **managers**.
- **miRNA (microRNA)** — very tiny RNA pieces whose job is to *shut genes down*. A miRNA finds a matching mRNA and silences it, so no protein gets made. Think **off-switches / brakes**.

---

## The ceRNA idea (the sponge story)

miRNAs (the off-switches) work by grabbing matching RNA. Normally they grab mRNAs and silence them. But lncRNAs can have the *same* matching spots that mRNAs do. So a lncRNA can act like a **sponge** — soaking up miRNAs and keeping them busy, so those miRNAs can't reach the mRNAs they were meant to silence.

Analogy: a predator (miRNA) eats a certain prey (mRNA). Drop in decoys (lncRNA) that look like prey. The predator wastes time chasing decoys, and the real prey survives — so more of that mRNA survives and makes more protein.

This chain — **lncRNA soaks up miRNA, which frees up mRNA** — is a **ceRNA network** (ceRNA = "competing endogenous RNA," just "RNAs competing for the same miRNAs"). Three layers:

**lncRNA → (soaks up) → miRNA → (would have silenced) → mRNA**

In diabetes these control chains can go wrong. Mapping them can reveal *why* certain genes misbehave. That map — the ceRNA network — is your final deliverable, and exactly what the job asks for.

---

## Why this dataset (diabetic islets)

The samples come from **pancreatic islets** — clusters of cells in the pancreas that make **insulin**, the hormone that controls blood sugar. In type-2 diabetes these cells stop working properly. Comparing diabetic vs healthy islets is one of the most direct ways to see diabetes at the RNA level.

---

## What each step gives you, and why

**Raw data (FASTQ files).** The sequencer doesn't say "gene X = 500 copies." It hands you millions of tiny RNA fragments, each written as letters (A, T, G, C), in files called **FASTQ**. On their own they're meaningless — like a shredded book on the floor. The next steps reassemble the book and count the pages.

**Quality control + trimming (FastQC / fastp).** Some fragments are low-quality or have leftover lab chemicals ("adapters") stuck on. This checks quality and snips off the bad bits. *Yield:* clean reads + a before/after report. (Cleaning smudges off the shreds.)

**Alignment (HISAT2).** Figure out *where in the human genome* each clean fragment came from — matching each shred back to its page. Output is a **BAM file** (every read + its location). *Yield:* you now know which gene each fragment belongs to.

**Counting (featureCounts).** Tally how many fragments landed on each gene. *Yield:* a **count table** — a grid of "gene × sample" with numbers. Higher number = that gene was more active in that sample.

**Differential expression (DESeq2).** Compare counts in diabetic vs healthy and ask, statistically, *"which genes are reliably different?"* A gene must be different *consistently*, beyond random noise. *Yield:* a table of **differentially expressed genes (DEGs)**, each with a fold-change (how much it moved) and a p-value (how confident). The **volcano plot** is just a picture of this table.

**Splitting mRNA vs lncRNA (biotype step).** Sort the changed genes: which are protein-coding (mRNA) and which are lncRNA. *Yield:* two clean lists — **DE mRNAs** and **DE lncRNAs** — the outer layers of the network.

**Enrichment (clusterProfiler).** Ask *"what do these changed genes have in common?"* (e.g. "many are involved in insulin secretion"). *Yield:* the biological *story* — which pathways are disturbed.

**miRNA prediction (multiMiR / starBase).** Standard RNA-seq can't see tiny miRNAs directly, so use databases that already know "miRNA A targets mRNA B" and "lncRNA C soaks up miRNA A." *Yield:* two connection lists (the network's wires).

**Building + drawing the network (igraph / Cytoscape).** Join the two wire-lists wherever they share the same miRNA, forming full lncRNA → miRNA → mRNA chains. *Yield:* the **ceRNA network**. The most-connected dots ("hubs") are your candidate master-regulators of diabetes.

---

## One honest note for interviews

miRNAs here come from *prediction databases*, not your own sequencing. That's standard, but be ready to say:

> *"Standard RNA-seq doesn't capture mature miRNAs; they need small-RNA-seq. So I built the miRNA layer from curated interaction databases, and the natural next step would be matched small-RNA-seq to confirm them experimentally."*

---

## FAQ (your follow-up questions)

**Q1. In the counting step, "gene X has 500 reads in sample 1, 1200 in sample 2" — what is a "sample"? Is it cells?**

A "sample" is **one biological specimen** — here, the islet material taken from **one donor (one person)**. It is *not* a single cell.

This dataset is **bulk RNA-seq**: the RNA from *all* the cells in that person's islet material is mixed together in one tube and measured as a group. So the number you get is the **average activity across thousands/millions of cells** from that one donor. "Sample 1" and "sample 2" are two different people (say, one healthy and one diabetic). "Gene X = 500 in sample 1, 1200 in sample 2" means gene X was about twice as active in donor 2 as in donor 1.

(There's another kind called **single-cell RNA-seq** that measures each individual cell separately — but that's not what you're doing here. In your project, sample = donor.)

**Q2. In DESeq2, what does "differentially expressed" mean?**

Break the phrase in two:

- **"Expressed"** = how active a gene is / how much RNA it's making. High counts = highly expressed. Low counts = barely expressed.
- **"Differentially"** = the difference *between your two groups* (diabetic vs healthy).

So a **differentially expressed gene** is one whose activity level is **reliably different in diabetics than in healthy people** — and "reliably" is the key word. DESeq2 doesn't just check if the average looks higher; it checks whether the difference is *consistent* across the donors and *bigger than random noise*. A gene can be differentially expressed in either direction: **up** (more active in diabetics) or **down** (less active in diabetics). Both count.

**Q3. In enrichment, what do "300 changed genes" mean — is it the differentially expressed genes, the high ones or the low ones?**

"300 changed genes" is exactly your **differentially expressed gene list from DESeq2** — the genes that came out as reliably different. ("Changed" just means "changed compared to healthy.")

Importantly, it includes **both directions**: genes that went **up** (higher reads in diabetics) *and* genes that went **down** (lower reads in diabetics). It's not only the high ones or only the low ones — it's every gene that moved significantly, up or down.

Enrichment then takes that whole list and looks for common themes ("a lot of these genes are involved in insulin secretion"). Tip: it's often insightful to run enrichment **separately** on the up genes and the down genes, because "what's being switched on in diabetes" and "what's being switched off" can tell different stories.
