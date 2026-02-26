# Data analysis of RNA-Seq Data with  DESeq

#Library loading
library(readxl)
library(tidyverse)
library(tidyselect)
library(tidyr)
library(broom)
library(ggplot2)
library(dplyr)
library(base)
library(DESeq2)
library(tibble)
library(pheatmap)
library(ashr) 
library(RColorBrewer) 
library(calibrate)
library(patchwork)
library(gplots)
library(EnhancedVolcano)

#Load sample information
sampleTable <- read.table("sample_table.txt",header=TRUE)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design=~condition)

vector	<- c("D3_1", "D3_2", "D3F_1", "D3F_2", "DC_1", "DC_2", "DCF_1", "DCF_2", "DD_1", "DD_2", "DDF_1", "DDF_2", "WT_1", "WT_2", "WTF_1", "WTF_2")

rawCounts <- counts(dds) %>%
  as_tibble() %>%
  pivot_longer(cols = all_of(vector),
               names_to = "Treatment",
               values_to = "Counts")
TotalRaw <- rawCounts %>%
  group_by(Treatment) %>%
  summarise(Total_counts = sum(Counts))
TotalRaw

posRawCounts <- rawCounts %>% 
  mutate(positiveTRC = Counts+1)

#PLOTS TOTAL RAW Counts
gg1 <- ggplot(TotalRaw, aes(y=Treatment, x=Total_counts, fill = Treatment))+ geom_col()+ ggtitle("Total RAW Count per Replicate")+xlab("Total Raw Counts")
gg2 <- ggplot(posRawCounts, aes(x=log(positiveTRC)))+ geom_histogram(bins=25)+ ggtitle("Raw Histogramm")+xlab("log(Raw Counts+1)")

(gg1 + gg2) + patchwork::plot_layout(guides = 'collect') & theme_bw() 

##Normalisation by size factor
ddsNorm <- estimateSizeFactors(dds)
ddsNorm

normCounts <- counts(ddsNorm, normalized=TRUE)
normCounts_tibble <- counts(ddsNorm, normalized=TRUE) %>%
  as_tibble() %>%
  pivot_longer(cols = all_of(vector),
               names_to = "Treatment",
               values_to = "Normalised_Counts")
TotalNorm <- normCounts_tibble %>%
  group_by(Treatment) %>%
  summarise(Total_Norm_counts = sum(Normalised_Counts))
TotalNorm

posNormCounts <- normCounts_tibble %>% 
  mutate(positiveTNC = Normalised_Counts+1)

#PLOTS TOTAL Counts AFTER NORMALIZATION
gg1N <- ggplot(TotalNorm, aes(y=Treatment, x=Total_Norm_counts, fill = Treatment))+ geom_col()+ ggtitle("Total Norm Counts per Replicate (After Normalization)")+xlab("Total Norm. Counts")
gg2N <- ggplot(posNormCounts, aes(x=log(positiveTNC)))+ geom_histogram(bins=25)+ ggtitle("Normalized Histogramm")+xlab("log(Normalized Counts+1)")

(gg1N + gg2N) + patchwork::plot_layout(guides = 'collect') & theme_bw() 

### IMPORTANT** Why are the total read counts still not identical even after normalisation? -> In the normalization process the `estimateSizeFactors` results in different total normalized counts. This is due to the process of the median-ratio method.

## Extracting differentially expressed genes

### Quality Control (Pair-wise correlation and Principal component analysis)

## Extracting differentially expressed genes

keep <- rowSums(normCounts) >= 10 #this checks whether the sum of read counts across samples is equal or higher than 10.
ddsExpressed <- ddsNorm[keep,]
ddsExpressed #Only 4 genes were not considered expressed

#Quality control: Pair-wise correlation
redblue <-rev(brewer.pal(11,"RdBu")) #red-to-blue color palette with 11 different shades

ddsExpressedCount <- counts(ddsExpressed)
ddsHeatmap <- cor((ddsExpressedCount), method = c("spearman"))

heatmap.2(ddsHeatmap, symm = TRUE, margin = c(5,5), trace="none", cexRow = 0.9, cexCol=0.8, dendrogram="column", lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(1.5, 4, 2 ) , col="heat.colors")

#Principal component analysis

rld <- rlog(ddsExpressed, blind=TRUE) #create a normalized data set compatible with the PCA function

redblue <-rev(brewer.pal(11,"RdBu"))
pcaData <- plotPCA(rld, intgroup=c("condition", "replicate"), returnData=TRUE) #create PCA object
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggPCA <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank())
#theme(legend.position = "bottom") + scale_color_manual(values=c("#4393C3", "#D6604D", "#67001F"))
ggPCA

ggsave("PCA_F_transporter_collab.png", plot=ggPCA)

#After we have made sure that the RNA-seq experiment worked in principle and that the data behaves as we expected, we can now do the actual RNA-seq analysis, ***i.e.*** ask how the transcriptional profile changes in dependence of the treatment, basically, which genes get upregulated, which ones overlap between treatments, etc. 

ddsExpressed$condition <- relevel(ddsExpressed$condition, ref = "control") ##which one is the reference/control sample i.e. glucose

#Next, we need to do the actual differential expression analysis. In DESeq2, the whole statistical analysis is wrapped in the _DESeq()_ function

ddsAnalysis <- DESeq(ddsExpressed)
ddsAnalysis

### Extraction of the results of the differential analysis into a table

#Same genotype - comparison versus F

# Define all contrasts as a list
contrast_list <- list(
  D3_vs_D3F   = c("condition","D3","D3F"),
  DC_vs_DCF   = c("condition","DC","DCF"),
  DD_vs_DDF   = c("condition","DD","DDF"),
  WT_vs_WTF   = c("condition","WT","WTF"),
  D3F_vs_WTF  = c("condition","D3F","WTF"), #Comparisons having F among genotypes
  DCF_vs_WTF  = c("condition","DCF","WTF"),
  DDF_vs_WTF  = c("condition","DDF","WTF"),
  DCF_vs_D3F  = c("condition","DCF","D3F"),
  DDF_vs_D3F  = c("condition","DDF","D3F"),
  DCF_vs_DDF  = c("condition","DCF","DDF"),
  D3_vs_WT    = c("condition","D3","WT"), #Comparisons minus F among genotypes
  DC_vs_WT    = c("condition","DC","WT"),
  DD_vs_WT    = c("condition","DD","WT"),
  DC_vs_D3    = c("condition","DC","D3"),
  DD_vs_D3    = c("condition","DD","D3"),
  DC_vs_DD    = c("condition","DC","DD")
)

# Compute all results
results_list <- lapply(contrast_list, function(ct) {
  res <- results(ddsAnalysis, contrast = ct)
  res <- as.data.frame(res) %>%
    rownames_to_column("gene_id")
  return(res)
})


# Visualizing fold changes

### Using the plotMA function from the DESeq2 package, we generate a Manhatten plot for each treatment. 

plotMA(D3_vs_D3F, main = "D3_vs_D3F")
plotMA(DC_vs_DCF, main = "DC_vs_DCF")
plotMA(DD_vs_DDF, main = "DD_vs_DDF")
plotMA(WT_vs_WTF, main = "WT_vs_WTF")

plotMA(D3F_vs_WTF, main = "D3F_vs_WTF")
plotMA(DCF_vs_WTF, main = "DCF_vs_WTF")
plotMA(DDF_vs_WTF, main = "DDF_vs_WTF")
plotMA(DCF_vs_D3F, main = "DCF_vs_D3F")
plotMA(DDF_vs_D3F, main = "DDF_vs_D3F")
plotMA(DCF_vs_DDF, main = "DCF_vs_DDF")

plotMA(D3_vs_WT, main = "D3_vs_WT")
plotMA(DC_vs_WT, main = "DC_vs_WT")
plotMA(DD_vs_WT, main = "DD_vs_WT")
plotMA(DC_vs_D3, main = "DC_vs_D3")
plotMA(DD_vs_D3, main = "DD_vs_D3")
plotMA(DC_vs_DD, main = "DC_vs_DD")

#Above, we ran the two treatment comparisons against the control sample separately, ending up with two tables that each contain the log2FoldChange and the adj. p-value. For the downstream analysis, it will be much more convenient to have these combined in a single table.

# Standardize column names for each contrast
results_mod <- lapply(names(results_list), function(name) {
  df <- results_list[[name]] %>%
    dplyr::select(gene_id, log2FoldChange, padj) %>%
    dplyr::rename(
      !!paste0("log2FC_", name) := log2FoldChange,
      !!paste0(name, "_padj")   := padj
    )
  return(df)
})
names(results_mod) <- names(results_list)

# Instead of dealing with all the genes in the dataset, we consider only those genes that are at least 2-fold up- or downregulated by at least one of the treatments and for which the adjusted p-value is < 0.1. ***DEG = Differentially Expressed Genes, LFC = log2 Fold Change***

# Helper function to merge multiple tables by gene_id
merge_all <- function(tbl_list) {
  Reduce(function(x, y) merge(x, y, by="gene_id", sort=FALSE), tbl_list)
}

combined_tbl      <- merge_all(results_mod[c("D3_vs_D3F","DC_vs_DCF","DD_vs_DDF","WT_vs_WTF")])
combined_F_tbl    <- merge_all(results_mod[c("D3F_vs_WTF","DCF_vs_WTF","DDF_vs_WTF","DCF_vs_D3F","DDF_vs_D3F","DCF_vs_DDF")])
combined_noF_tbl  <- merge_all(results_mod[c("D3_vs_WT","DC_vs_WT","DD_vs_WT","DC_vs_D3","DD_vs_D3","DC_vs_DD")])

# Generic DEG filter function
filter_deg <- function(tbl, fc_cols, padj_cols, fc_cut=2, padj_cut=0.1) {
  tbl %>%
    filter(if_any(all_of(fc_cols), ~ .x > fc_cut | .x < -fc_cut)) %>%
    filter(if_any(all_of(padj_cols), ~ .x < padj_cut))
}

## Filter DEG tables

DEG_tbl     <- filter_deg(combined_tbl,     grep("log2FC", names(combined_tbl), value=TRUE),     grep("padj", names(combined_tbl), value=TRUE))
DEG_F_tbl   <- filter_deg(combined_F_tbl,   grep("log2FC", names(combined_F_tbl), value=TRUE),   grep("padj", names(combined_F_tbl), value=TRUE))
DEG_noF_tbl <- filter_deg(combined_noF_tbl, grep("log2FC", names(combined_noF_tbl), value=TRUE), grep("padj", names(combined_noF_tbl), value=TRUE))

# Add gene names

KT2440_gene_references <- read.csv("KT2440_gene_references.csv", sep=";")

DEG_tbl_names     <- inner_join(DEG_tbl,     KT2440_gene_references, by="gene_id")
DEG_F_tbl_names   <- inner_join(DEG_F_tbl,   KT2440_gene_references, by="gene_id")
DEG_noF_tbl_names <- inner_join(DEG_noF_tbl, KT2440_gene_references, by="gene_id")

#Save outputs
write.csv(DEG_tbl_names,     "DEG_tbl_names_2.csv",     row.names = FALSE)
write.csv(DEG_F_tbl_names,   "DEG_F_tbl_names_2.csv",   row.names = FALSE)
write.csv(DEG_noF_tbl_names, "DEG_noF_tbl_names_2.csv", row.names = FALSE)

### Heatmap of DEG via fold change

#Plot heatmap
LFC_DEG <- DEG_tbl %>%
  dplyr::select(c(log2FC_D3_vs_D3F, log2FC_DC_vs_DCF, log2FC_DD_vs_DDF, log2FC_WT_vs_WTF))

LFC_F_DEG <- DEG_F_tbl %>%
  dplyr::select(c(log2FC_D3F_vs_WTF, log2FC_DCF_vs_WTF, log2FC_DDF_vs_WTF, log2FC_DCF_vs_D3F, log2FC_DDF_vs_D3F, log2FC_DCF_vs_DDF))

redblue <-rev(brewer.pal(11,"RdBu")) #red-to-blue color palette with 11 different shades
pheatmap(LFC_DEG,
         breaks = c(-5:7),
         cluster_rows=T,
         cluster_cols=F,
         col=redblue,
         angle_col=0,
         cellwidth=70)

pheatmap(LFC_F_DEG,
         breaks = c(-5:7),
         cluster_rows=T,
         cluster_cols=F,
         col=redblue,
         angle_col=0,
         cellwidth=70)

# VOLCANO PLOTS

EnhancedVolcano(D3_vs_D3F, #Change with desired comparison
                lab = rownames(D3_vs_D3F), #Change with desired comparison
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                pCutoff = 10e-16,
                FCcutoff = 2.0,
                pointSize = 3.0,
                labSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.50)
