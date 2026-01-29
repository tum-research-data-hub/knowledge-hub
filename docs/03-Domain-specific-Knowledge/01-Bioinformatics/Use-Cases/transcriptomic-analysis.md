---
title: Transcriptomic Analysis
slug: /domain-knowledge/bioinformatics/use-cases/transcriptomic-analysis
---

# **Use Case 4: Transcriptomic Analysis**

This type of analysis starts with RNA sequencing data in the form of sequencing reads, typically stored in FASTQ files, which contain both nucleotide sequences and their per-base quality scores. In transcriptomics, these reads represent RNA molecules that were present in a cell or tissue at the time of sampling and are commonly generated through RNA-seq experiments.

## **Workflow Overview**

```
Read mapping / quantification → Expression analysis → Interpretation
FASTQ (+ reference genome/transcriptome) → Expression matrix → Differential expression results
```

## **Key Concepts**

### **Read Mapping and Quantification**

Before analysis, RNA reads are typically mapped (aligned) to a reference genome or transcriptome, similar to reference-based DNA analysis. Alternatively, **alignment-free** or **pseudo-alignment** methods may be used to directly assign reads to known transcripts. 

The goal of this step is not to identify sequence variation, but to **quantify gene or transcript expression levels**.

### **Expression Quantification**

The mapping or quantification process produces a file that links sequencing reads to genes or transcripts, along with quality and confidence measures. From this information, **expression counts** or **abundance estimates** are generated, describing how frequently each gene or transcript is observed in the data.

### **Downstream Analysis**

Once expression levels have been quantified, downstream analyses can be performed:

- **Differential expression analysis** - Identifies genes whose expression differs significantly between experimental conditions
- **Clustering** - Groups samples or genes with similar expression patterns
- **Gene set enrichment** - Tests whether predefined gene sets show coordinated expression changes
- **Pathway analysis** - Identifies biological pathways affected by experimental conditions

## **Applications**

Transcriptomic analysis provides insight into the functional activity of a genome under specific conditions. It is widely used to study:

- Cellular responses to environmental changes
- Disease mechanisms
- Developmental processes
- Effects of genetic or pharmacological perturbations

Common applications:
- Understanding cellular responses to stimuli (e.g., drugs, toxins)
- Cell type identification
- Improving cellular characteristics

## **Popular Tools**

### Alignment & Quantification
- **STAR** - Fast RNA-seq aligner
- **HISAT2** - Hierarchical indexing aligner
- **Salmon** - Alignment-free quantification
- **kallisto** - Pseudo-alignment quantification
- **featureCounts** - Read counting

### Differential Expression
- **DESeq2** (R) - Statistical analysis
- **edgeR** (R) - Differential expression
- **limma** (R) - Linear models

### Visualization & Interpretation
- **ggplot2** (R) - Publication-quality plots
- **pheatmap** (R) - Heatmaps
- **clusterProfiler** (R) - Enrichment analysis

## **Code Example: RNA-seq Analysis Workflow**

<details>
<summary>**Click to expand Bash workflow with STAR**</summary>

```bash
#!/bin/bash
# RNA-seq transcriptomic analysis workflow

# Step 1: Quality control
fastqc sample_*.fastq.gz -o qc_before/
multiqc qc_before/ -o qc_before/

# Step 2: Trim adapters and low-quality bases
for sample in sample1 sample2 sample3 control1 control2 control3; do
    trim_galore --paired \
        --quality 20 \
        --length 36 \
        --fastqc \
        -o trimmed/ \
        ${sample}_R1.fastq.gz \
        ${sample}_R2.fastq.gz
done

# Step 3: Build STAR index (only needed once)
STAR --runMode genomeGenerate \
    --genomeDir star_index/ \
    --genomeFastaFiles reference_genome.fasta \
    --sjdbGTFfile annotation.gtf \
    --sjdbOverhang 100 \
    --runThreadN 16

# Step 4: Align reads with STAR
for sample in sample1 sample2 sample3 control1 control2 control3; do
    STAR --genomeDir star_index/ \
        --readFilesIn trimmed/${sample}_R1_val_1.fq.gz trimmed/${sample}_R2_val_2.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix aligned/${sample}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --runThreadN 8
done

# Step 5: Count reads per gene
featureCounts -p -T 8 \
    -a annotation.gtf \
    -o counts.txt \
    aligned/*_Aligned.sortedByCoord.out.bam

# Step 6: MultiQC report for all samples
multiqc aligned/ -o multiqc_report/

echo "Alignment and quantification complete!"
echo "Count matrix: counts.txt"
```

</details>

## **Alternative: Alignment-Free Quantification**

<details>
<summary>**Click to expand Salmon quantification workflow**</summary>

```bash
#!/bin/bash
# Salmon alignment-free quantification

# Step 1: Build Salmon index (only needed once)
salmon index \
    -t transcriptome.fasta \
    -i salmon_index \
    -k 31

# Step 2: Quantify with Salmon
for sample in sample1 sample2 sample3 control1 control2 control3; do
    salmon quant \
        -i salmon_index \
        -l A \
        -1 trimmed/${sample}_R1.fastq.gz \
        -2 trimmed/${sample}_R2.fastq.gz \
        -p 8 \
        --validateMappings \
        -o salmon_quant/${sample}
done

echo "Salmon quantification complete!"
```

</details>

## **R Example: Differential Expression Analysis with DESeq2**

<details>
<summary>**Click to expand R script with DESeq2**</summary>

```r
#!/usr/bin/env Rscript
# Differential expression analysis with DESeq2

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)

# ============================================================
# Step 1: Load count data
# ============================================================

# Read featureCounts output
counts <- read.table("counts.txt", header = TRUE, row.names = 1, skip = 1)

# Extract only count columns (remove annotation columns)
count_matrix <- counts[, grep("Aligned", colnames(counts))]

# Clean column names
colnames(count_matrix) <- gsub("_Aligned.sortedByCoord.out.bam", "", 
                                gsub("aligned.", "", colnames(count_matrix)))

# ============================================================
# Step 2: Create sample metadata
# ============================================================

sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(c(rep("treatment", 3), rep("control", 3))),
  row.names = colnames(count_matrix)
)

print("Sample information:")
print(sample_info)

# ============================================================
# Step 3: Create DESeq2 object
# ============================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# Pre-filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat(sprintf("\nRetained %d genes after filtering\n", sum(keep)))

# ============================================================
# Step 4: Run DESeq2 analysis
# ============================================================

dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "treatment", "control"))

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Summary
summary(res)

# ============================================================
# Step 5: Save results
# ============================================================

# Convert to data frame and save
res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj))

write_csv(res_df, "differential_expression_results.csv")

# Save significantly differentially expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

write_csv(sig_genes, "significant_genes.csv")

cat(sprintf("\nFound %d significantly DE genes (padj < 0.05, |log2FC| > 1)\n", 
            nrow(sig_genes)))

# ============================================================
# Step 6: Visualizations
# ============================================================

# MA plot
pdf("MA_plot.pdf", width = 8, height = 6)
plotMA(res, ylim = c(-5, 5), main = "MA Plot")
dev.off()

# Volcano plot
volcano_data <- res_df %>%
  mutate(
    significant = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("volcano_plot.png", width = 10, height = 8, dpi = 300)

# Heatmap of top 50 genes
vsd <- vst(dds, blind = FALSE)
top_genes <- head(rownames(res_ordered), 50)

pheatmap(
  assay(vsd)[top_genes, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = sample_info["condition"],
  main = "Top 50 Differentially Expressed Genes",
  filename = "heatmap_top50.png",
  width = 8,
  height = 10
)

# PCA plot
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(title = "PCA Plot") +
  theme_minimal()

ggsave("pca_plot.png", width = 8, height = 6, dpi = 300)

cat("\nAnalysis complete! Results and plots saved.\n")
```

</details>

## **Python Example: Gene Expression Clustering**

<details>
<summary>**Click to expand Python script**</summary>

```python
#!/usr/bin/env python3
"""
Clustering and visualization of gene expression data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import dendrogram, linkage

def load_expression_data(file_path):
    """
    Load normalized expression data
    
    Args:
        file_path: Path to expression matrix (genes x samples)
        
    Returns:
        DataFrame with expression values
    """
    return pd.read_csv(file_path, index_col=0)

def filter_variable_genes(expr_df, top_n=1000):
    """
    Select most variable genes
    
    Args:
        expr_df: Expression DataFrame
        top_n: Number of top variable genes to keep
        
    Returns:
        Filtered DataFrame
    """
    # Calculate variance for each gene
    gene_var = expr_df.var(axis=1)
    
    # Select top N most variable genes
    top_genes = gene_var.nlargest(top_n).index
    
    return expr_df.loc[top_genes]

def cluster_genes(expr_df, n_clusters=5):
    """
    Perform k-means clustering on genes
    
    Args:
        expr_df: Expression DataFrame
        n_clusters: Number of clusters
        
    Returns:
        Series with cluster assignments
    """
    # Standardize data
    scaler = StandardScaler()
    expr_scaled = scaler.fit_transform(expr_df)
    
    # Perform k-means
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(expr_scaled)
    
    return pd.Series(clusters, index=expr_df.index, name='cluster')

def plot_expression_heatmap(expr_df, output_file='expression_heatmap.png'):
    """
    Create clustered heatmap
    
    Args:
        expr_df: Expression DataFrame
        output_file: Output file path
    """
    # Standardize for visualization
    scaler = StandardScaler()
    expr_scaled = scaler.fit_transform(expr_df.T).T
    
    # Create heatmap
    plt.figure(figsize=(12, 10))
    sns.clustermap(
        expr_scaled,
        cmap='RdBu_r',
        center=0,
        yticklabels=False,
        figsize=(12, 10)
    )
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to {output_file}")

# Usage example
if __name__ == '__main__':
    # Load data
    expr_data = load_expression_data('normalized_counts.csv')
    
    # Filter for variable genes
    variable_genes = filter_variable_genes(expr_data, top_n=500)
    
    # Cluster genes
    gene_clusters = cluster_genes(variable_genes, n_clusters=6)
    
    # Print cluster sizes
    print("\nCluster sizes:")
    print(gene_clusters.value_counts().sort_index())
    
    # Create heatmap
    plot_expression_heatmap(variable_genes)
```

</details>

## **Expected Outputs**

A complete RNA-seq analysis produces:

- **Count matrix** - Gene × sample table of read counts
- **Quality reports** - MultiQC summary of alignment and QC metrics
- **Differential expression results** - List of significant genes with fold changes
- **Visualizations** - MA plots, volcano plots, heatmaps, PCA plots
- **Enrichment results** - GO terms and pathways associated with DE genes

### **Interpreting Results**

| Metric | Good | Acceptable | Concerning |
|--------|------|------------|------------|
| **Alignment rate** | `>85%` | 70-85% | `<70%` |
| **Assigned reads** | `>60%` | 40-60% | `<40%` |
| **DE genes (typical)** | 100-5000 | 50-100 or 5000+ | `<10` or `>10,000` |
| **Sample correlation** | `>0.9` (replicates) | 0.8-0.9 | `<0.8` |

### **Statistical Thresholds**

- **Adjusted p-value** - Typically padj < 0.05 (sometimes 0.01 for stringency)
- **Fold change** - Often |log2FC| > 1 (2-fold change), adjust based on biology
- **Base mean** - Filter low-expressed genes (e.g., baseMean > 10)

## **Computational Requirements**

### **Resource Estimates**

| Step | Samples | CPU Cores | RAM | Storage | Time |
|------|---------|-----------|-----|---------|------|
| Alignment (STAR) | 6 samples | 8 per sample | 32 GB | ~100 GB | 2-4 hours |
| Quantification (Salmon) | 6 samples | 8 per sample | 8 GB | ~50 GB | 30-60 min |
| DESeq2 analysis | Any | 1-4 | 8-16 GB | Minimal | 10-30 min |

:::tip Performance Tip
**Salmon** is 10-100x faster than alignment-based quantification and produces comparable results for most applications.
:::

## **Common Issues & Troubleshooting**

:::warning RNA-seq Problems

**Low alignment rate**
- Wrong reference genome or annotation version
- Heavy rRNA contamination (use rRNA depletion kits)
- Poor quality libraries
- Adapter contamination

**Few or no differentially expressed genes**
- Insufficient biological effect (treatments not different enough)
- Too much biological variability (need more replicates)
- Batch effects overwhelming signal
- Incorrect experimental design in DESeq2

**Too many DE genes (>50% of genome)**
- Samples from different tissues/cell types
- Major batch effects not accounted for
- Wrong comparison group
- Contamination or sample mix-up

**PCA shows outlier samples**
- Technical failure (exclude if confirmed)
- Biological outlier (may be real)
- Batch effects (correct with ComBat-seq or include in model)
- Sample swap (check metadata)

**Low assignment rate**
- Annotation doesn't match genome version
- Many reads in intergenic regions (genomic DNA contamination)
- Annotation incomplete (common in non-model organisms)

**Memory errors in DESeq2**
- Too many genes (filter low-count genes first)
- Very large sample size (subsample for testing)
- Use more RAM or run on cluster

:::

## **Sample Data & Practice**

Public RNA-seq datasets:
- **GEO (Gene Expression Omnibus)** - [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/)
- **ENA RNA-seq studies** - [EBI ENA](https://www.ebi.ac.uk/ena/browser/)
- **ENCODE** - [ENCODE Portal](https://www.encodeproject.org/)
- **GTEx** - Human tissue expression [GTEx Portal](https://gtexportal.org/)

## **Next Steps & Related Analyses**

After differential expression analysis:

### **Functional Interpretation**
- **Gene Ontology enrichment** - Find enriched biological processes (clusterProfiler, DAVID)
- **KEGG pathway analysis** - Identify affected metabolic/signaling pathways
- **Gene Set Enrichment Analysis (GSEA)** - Ranked gene list analysis
- **Network analysis** - Identify hub genes and modules (WGCNA)

### **Visualization**
- **Interactive exploration** - Use tools like iDEP, Degust, or custom Shiny apps
- **Pathway visualization** - KEGG visualization, Reactome
- **Alluvial plots** - Show gene overlap across comparisons

### **Advanced Analyses**
- **Isoform-level analysis** - Differential transcript usage (DRIMSeq, DEXSeq)
- **Allele-specific expression** - Detect imprinting or cis-regulatory effects
- **Alternative splicing** - MISO, rMATS for splice junction analysis
- **Co-expression networks** - WGCNA to identify gene modules

### **Validation**
- **qRT-PCR** - Confirm top DE genes
- **Protein validation** - Western blots, immunofluorescence
- **Functional assays** - CRISPR knockouts, overexpression studies

**Related Use Cases:**
- [Reference-Based Analysis](/domain-knowledge/bioinformatics/use-cases/reference-based-analysis) - Variant calling from RNA-seq
- [De Novo Assembly](/domain-knowledge/bioinformatics/use-cases/de-novo-analysis) - For organisms without references

## **Key Considerations**

- **Biological replicates are essential** - Minimum 3 per condition, 4-6 recommended for more power
- **Library preparation matters** - Know whether your data is stranded or unstranded
- **Normalization is critical** - DESeq2 and edgeR handle this internally (don't use RPKM/FPKM for DE)
- **Batch effects** - Account for technical variation in experimental design (date, lane, technician)
- **Multiple testing correction** - Always use adjusted p-values (padj/FDR), never raw p-values
- **Validation** - Confirm key findings with qRT-PCR or other orthogonal methods
- **Experimental design** - Plan analysis strategy before sequencing (power analysis, blocking design)
