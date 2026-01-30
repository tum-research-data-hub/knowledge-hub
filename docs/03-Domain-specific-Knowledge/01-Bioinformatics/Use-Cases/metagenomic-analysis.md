---
title: Metagenomic Analysis
slug: /domain-knowledge/bioinformatics/use-cases/metagenomic-analysis
---

# **Use Case 3: Metagenomic Analysis**

This type of analysis starts with DNA sequencing data in the form of sequencing reads, typically stored in FASTQ files, which contain both nucleotide sequences and their per-base quality scores. In metagenomics, these reads originate from a mixture of multiple organisms present in a single sample, such as environmental, clinical, or host-associated communities.

## **Workflow Overview**

```
Optional assembly → Binning / Profiling → Functional analysis
FASTQ → Contigs (FASTA) → MAGs / Taxonomic & functional profiles
```

## **Key Concepts**

### **Why Metagenomics?**

Because no single reference genome can adequately represent all organisms in a mixed sample, metagenomic analyses often rely on de novo approaches to reconstruct genomic fragments from the sequencing reads. Depending on the research question, reads may be:
- Assembled into contigs
- Analyzed directly without full assembly

### **Binning**

Following assembly, contigs can be grouped into bins, each representing a putative genome or taxonomic unit. This process, known as **binning**, uses:
- Sequence composition (e.g., GC content, tetranucleotide frequencies)
- Coverage patterns
- Similarity to known genomes

The resulting bins are commonly referred to as **metagenome-assembled genomes (MAGs)**.

### **Profiling**

In addition to genome reconstruction, metagenomic workflows frequently include:
- **Taxonomic profiling** - Estimates which organisms are present in the sample
- **Functional profiling** - Identifies metabolic pathways or gene functions encoded by the community

These analyses allow researchers to characterize both the composition and functional potential of complex microbial ecosystems.

## **Applications**

Metagenomic analysis is widely used to study microbial communities in environments where cultivation of individual organisms is difficult or impossible.

Typical applications include:
- Microbiome research (gut, skin, environmental)
- Environmental monitoring
- Pathogen discovery
- Investigation of community-level responses to environmental or experimental changes

## **Popular Tools**

### Assembly & Binning
- **MEGAHIT** - Fast metagenomic assembler
- **MetaSPAdes** - SPAdes for metagenomes
- **MetaBAT2** - Binning tool
- **MaxBin2** - Automated binning
- **CheckM** - Assess MAG quality

### Taxonomic Profiling
- **Kraken2** - Fast taxonomic classification
- **MetaPhlAn** - Marker-based profiling
- **Bracken** - Species abundance estimation

### Functional Profiling
- **HUMAnN** - Functional profiling
- **DIAMOND** - Fast protein alignment
- **eggNOG-mapper** - Functional annotation

## **Code Example: Metagenomic Analysis Workflow**

<details>
<summary>**Click to expand Bash workflow**</summary>

```bash
#!/bin/bash
# Metagenomic analysis workflow

# Step 1: Quality control
fastqc metagenome_R1.fastq.gz metagenome_R2.fastq.gz -o qc_results/

# Step 2: Remove host contamination (if applicable)
# Map reads to host genome and extract unmapped reads
bowtie2 -x host_genome_index \
    -1 metagenome_R1.fastq.gz \
    -2 metagenome_R2.fastq.gz \
    --un-conc-gz clean_reads.fastq.gz \
    -S host_mapped.sam

# Step 3: Taxonomic profiling with Kraken2
kraken2 --db minikraken2_v2 \
    --paired clean_reads_1.fastq.gz clean_reads_2.fastq.gz \
    --threads 16 \
    --report kraken_report.txt \
    --output kraken_output.txt

# Generate visualization with Bracken
bracken -d minikraken2_v2 \
    -i kraken_report.txt \
    -o bracken_species.txt \
    -r 150 \
    -l S

# Step 4: Metagenomic assembly
megahit -1 clean_reads_1.fastq.gz \
    -2 clean_reads_2.fastq.gz \
    -o assembly_output \
    --min-contig-len 1000 \
    -t 16

# Step 5: Binning - map reads back to assembly
bowtie2-build assembly_output/final.contigs.fa assembly_index
bowtie2 -x assembly_index \
    -1 clean_reads_1.fastq.gz \
    -2 clean_reads_2.fastq.gz \
    -S mapped_reads.sam \
    --threads 16

# Convert to BAM and sort
samtools view -bS mapped_reads.sam | samtools sort -o mapped_reads.sorted.bam
samtools index mapped_reads.sorted.bam

# Bin contigs using MetaBAT2
metabat2 -i assembly_output/final.contigs.fa \
    -a mapped_reads.sorted.bam \
    -o bins/bin \
    -t 16 \
    -m 1500

# Step 6: Assess MAG quality with CheckM
checkm lineage_wf bins/ checkm_output/ \
    -x fa \
    -t 16 \
    --tab_table \
    -f checkm_results.txt

# Step 7: Functional profiling with HUMAnN
humann --input clean_reads_1.fastq.gz \
    --output humann_output \
    --threads 16

echo "Metagenomic analysis complete!"
```

</details>

## **Python Example: Taxonomic Profile Analysis**

<details>
<summary>**Click to expand Python script**</summary>

```python
#!/usr/bin/env python3
"""
Analyze and visualize Kraken taxonomic profile
"""

import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

def parse_kraken_report(report_file):
    """
    Parse Kraken2 report file
    
    Args:
        report_file: Path to Kraken report
        
    Returns:
        DataFrame with taxonomic information
    """
    columns = ['percentage', 'reads_clade', 'reads_taxon', 'rank', 'taxid', 'name']
    df = pd.read_csv(report_file, sep='\t', names=columns)
    
    # Clean up name column (remove leading spaces)
    df['name'] = df['name'].str.strip()
    
    return df

def summarize_taxonomy(df, rank='S', top_n=10):
    """
    Summarize taxonomy at specific rank
    
    Args:
        df: Kraken report DataFrame
        rank: Taxonomic rank (S=species, G=genus, F=family, etc.)
        top_n: Number of top taxa to return
        
    Returns:
        DataFrame with top taxa
    """
    # Filter by rank
    rank_df = df[df['rank'] == rank].copy()
    
    # Sort by percentage
    rank_df = rank_df.sort_values('percentage', ascending=False)
    
    # Get top N
    top_taxa = rank_df.head(top_n)
    
    return top_taxa[['name', 'percentage', 'reads_clade']]

def plot_taxonomic_composition(df, rank='S', top_n=10, output_file='taxonomy_plot.png'):
    """
    Create bar plot of taxonomic composition
    
    Args:
        df: Kraken report DataFrame
        rank: Taxonomic rank
        top_n: Number of top taxa to plot
        output_file: Output file path
    """
    top_taxa = summarize_taxonomy(df, rank, top_n)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.barh(top_taxa['name'], top_taxa['percentage'])
    ax.set_xlabel('Percentage of Reads (%)')
    ax.set_ylabel('Organism')
    ax.set_title(f'Top {top_n} Most Abundant Taxa (Rank: {rank})')
    ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved to {output_file}")

def calculate_diversity_metrics(df, rank='S'):
    """
    Calculate alpha diversity metrics
    
    Args:
        df: Kraken report DataFrame
        rank: Taxonomic rank
    """
    rank_df = df[df['rank'] == rank]
    
    # Shannon diversity
    proportions = rank_df['reads_clade'] / rank_df['reads_clade'].sum()
    proportions = proportions[proportions > 0]  # Remove zeros
    shannon = -sum(proportions * np.log(proportions))
    
    # Species richness
    richness = len(rank_df[rank_df['reads_clade'] > 0])
    
    print(f"\nDiversity Metrics (Rank: {rank}):")
    print(f"Species Richness: {richness}")
    print(f"Shannon Diversity: {shannon:.3f}")

# Usage example
if __name__ == '__main__':
    import numpy as np
    
    # Parse Kraken report
    df = parse_kraken_report('kraken_report.txt')
    
    # Print summary at species level
    print("Top 10 Most Abundant Species:")
    print(summarize_taxonomy(df, rank='S', top_n=10))
    
    # Calculate diversity
    calculate_diversity_metrics(df, rank='S')
    
    # Create visualization
    plot_taxonomic_composition(df, rank='S', top_n=15)
```

</details>

## **R Example: Functional Pathway Analysis**

<details>
<summary>**Click to expand R script**</summary>

```r
#!/usr/bin/env Rscript
# Analyze HUMAnN functional profiles

library(tidyverse)

# Read HUMAnN pathway abundance file
pathways <- read_tsv("humann_output/pathways_abundance.tsv")

# Clean column names
colnames(pathways) <- gsub("_Abundance", "", colnames(pathways))

# Calculate summary statistics
pathway_summary <- pathways %>%
  summarise(
    total_pathways = n(),
    mean_abundance = mean(`sample_name`, na.rm = TRUE),
    sd_abundance = sd(`sample_name`, na.rm = TRUE)
  )

print(pathway_summary)

# Plot top pathways
top_pathways <- pathways %>%
  arrange(desc(`sample_name`)) %>%
  head(20)

ggplot(top_pathways, aes(x = reorder(`# Pathway`, `sample_name`), y = `sample_name`)) +
  geom_col(fill = "#0065BD") +
  coord_flip() +
  labs(
    title = "Top 20 Metabolic Pathways",
    x = "Pathway",
    y = "Relative Abundance"
  ) +
  theme_minimal()

ggsave("top_pathways.png", width = 12, height = 8, dpi = 300)
```

</details>

## **Expected Outputs**

A comprehensive metagenomic analysis produces:

- **Taxonomic profile** - Species abundance table and visualization
- **MAGs (Metagenome-Assembled Genomes)** - Individual genome bins
- **Quality reports** - CheckM completeness and contamination scores
- **Functional profile** - Pathway abundance and gene family counts
- **Diversity metrics** - Alpha and beta diversity measures

### **Interpreting Results**

| Metric | Good MAG | Medium MAG | Poor MAG |
|--------|----------|------------|----------|
| **Completeness** | `>90%` | 50-90% | `<50%` |
| **Contamination** | `<5%` | 5-10% | `>10%` |
| **Quality score** | `>80` | 50-80 | `<50` |

**Quality score** = Completeness - 5 × Contamination

### **Typical Taxonomic Diversity**

- **Gut microbiome** - 500-1000 species, dominated by Firmicutes and Bacteroidetes
- **Soil** - 1000s of species, highly diverse
- **Marine** - Variable, often dominated by few abundant species

## **Computational Requirements**

### **Resource Estimates**

| Analysis Type | Reads | CPU Cores | RAM | Storage | Time |
|---------------|-------|-----------|-----|---------|------|
| Taxonomic profiling only | 10M pairs | 8-16 | 16-32 GB | ~10 GB | 1-4 hours |
| Assembly + binning | 50M pairs | 32-64 | 128-256 GB | ~500 GB | 1-2 days |
| Deep metagenomics | 200M+ pairs | 64+ | 256-512 GB | ~2 TB | 3-7 days |

:::tip Resource Recommendation
Metagenomic assembly is extremely resource-intensive. Use **LRZ AI Systems** or high-memory cluster nodes. See [Infrastructure Guide](/02-General-Knowledge/infrastructure).
:::

## **Common Issues & Troubleshooting**

:::warning Metagenomic Challenges

**Low taxonomic assignment rates**
- Database may be outdated - update Kraken2/MetaPhlAn databases
- Novel organisms not in reference databases
- Poor sequencing quality - check with FastQC
- Host contamination - filter host reads first

**Poor MAG quality**
- Insufficient sequencing depth (aim for 10-50 Gb total)
- Low abundance organisms (`<1%`) difficult to bin
- Strain heterogeneity within species
- Try multiple binning tools and use DAS Tool to merge results

**Memory errors during assembly**
- Use MEGAHIT instead of MetaSPAdes (more memory efficient)
- Reduce k-mer size
- Subsample reads for pilot assembly
- Use cluster with high-memory nodes

**Too many/too few bins**
- Adjust minimum contig length (1000-2000 bp)
- Modify MetaBAT2 sensitivity parameters
- Check coverage depth across samples
- Consider differential coverage binning

**Functional profiling failures**
- HUMAnN requires substantial resources
- Try subsampling to 10M reads first
- Ensure databases are properly installed
- Use pre-computed profiles when available

:::

## **Sample Data & Practice**

Public metagenomic datasets:
- **Human Microbiome Project** - [HMP Data Portal](https://hmpdacc.org/)
- **Earth Microbiome Project** - [EMP Portal](http://www.earthmicrobiome.org/)
- **NCBI SRA** - Filter by "metagenome" datatype
- **MG-RAST** - Metagenomic analysis server with example datasets

## **Next Steps & Related Analyses**

After metagenomic analysis:
- **MAG refinement** - Manual curation and polishing
- **Metabolic reconstruction** - KEGG pathway analysis
- **Comparative metagenomics** - Compare communities across conditions
- **Integration with metadata** - Correlate with environmental/clinical variables
- **Time series analysis** - Track community changes over time
- **MAG phylogeny** - Place genomes in evolutionary context

**Related Use Cases:**
- [De Novo Assembly](/domain-knowledge/bioinformatics/use-cases/de-novo-analysis) - Assembling individual MAGs
- [Transcriptomics](/domain-knowledge/bioinformatics/use-cases/transcriptomic-analysis) - Metatranscriptomics for active genes

## **Key Considerations**

- **Sequencing depth is critical** - Deeper sequencing captures more community diversity (aim for 10-50 Gb)
- **Assembly vs. read-based** - Choose based on whether you need MAGs or just profiles
- **Database selection** - Taxonomic assignment quality depends on reference database completeness
- **MAG quality assessment** - Always check completeness and contamination with CheckM
- **Multiple profiling methods** - Different tools may give different results; use complementary approaches
- **Replication is essential** - Biological variability high in microbiomes; use multiple samples per condition
- **Controls matter** - Include negative controls to assess contamination
