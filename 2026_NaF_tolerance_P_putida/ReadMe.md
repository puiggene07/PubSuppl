## RNAseq analysis on the adaptaton of *P. putida* towards sodium fluoride in R

- This serves as the Supplementary Information for the RNAseq data performed in the article "Adaptive evolution of *Pseudomonas putida* in the presence of fluoride exposes novel functions of a benzoate transporter" (JB00479-25R2) by Lea Ets, Heili Ilves, Lisette Juhe, Tanel Ilmjärv, Òscar Puiggené, Pablo Ivan Nikel, and Maia Kivisaar.
- Included here is RNA‑Seq differential expression analysis using DESeq2 which includes QC, normalization, PCA, correlation heatmaps, DE analysis, volcano plots, and DEG extraction.

## Notes on the RNAseq Analysis

- Raw data was publicly stored in NCBI with the GEO accession [GSE309641](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309641)

### Required input files: 

1. FPKM HTSeq counts provided here under raw_data/ folder.
2. Sample table as *sample_table.txt*
3. KT2440_gene references.csv to transform gene IDs to PP_#### numbers,

### Important methodological notes: 

- DESeq2’s median‑ratio method means normalized counts do not sum to the same value across samples.
- Genes with total normalized counts < 10 are removed.
- Multiple testing correction: Benjamini–Hochberg FDR.
- PCA uses rlog transformation.
- Dependencies are listed in the attached code.

### Pipeline overview

1. Load raw counts
2. Normalize using DESeq2 size factors
3. QC: total counts, histograms, PCA, correlation heatmap
4. Differential expression using DESeq2
5. Automated contrast generation
6. DEG filtering: |log2FC| > 2 and padj < 0.1
7. Annotation with gene reference table
8. Heatmaps and volcano plots
9. Export of DEG tables

### Output files

- PCA plot
- Correlation heatmap
- Volcano plots
- DEG tables for comparisons:
    A. Within‑genotype ±F
    B.Between genotypes with F
    C. Between genotypes without F
