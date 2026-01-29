---
title: De Novo Sequence Analysis
slug: /domain-knowledge/bioinformatics/use-cases/de-novo-analysis
---

# **Use Case 2: De Novo Sequence Analysis**

This type of analysis also starts with DNA sequencing data in the form of sequencing reads, typically stored in FASTQ files, which contain both nucleotide sequences and their per-base quality scores. In contrast to reference-based approaches, no suitable reference genome is assumed to be available, or the goal is to reconstruct the genome independently of an existing reference.

## **Workflow Overview**

```
De novo assembly → Genome annotation (→ optional comparative analysis)
FASTQ → Contigs/Scaffolds (FASTA) → Annotated genome
```

## **Key Concepts**

### **De Novo Assembly**

Instead of mapping reads to a reference, the sequencing reads are assembled **de novo** (Latin for "from the beginning"), meaning they are computationally combined based on sequence overlap and consistency to reconstruct longer contiguous sequences (contigs). These contigs may be further connected into larger structures called scaffolds, representing an approximation of the original genome sequence.

### **Assembly Output**

The result of the assembly process is typically stored in FASTA format, containing the reconstructed contig or scaffold sequences. Quality metrics associated with the assembly include:
- Contig length distributions
- Coverage statistics
- Measures of completeness
- N50 values (median contig length)

### **Genome Annotation**

Once a genome has been assembled, an additional analysis step known as genome annotation is commonly performed. During annotation, genomic features such as genes, coding sequences, and regulatory elements are predicted and assigned functional information based on:
- Sequence patterns
- Similarity to known genes
- External databases

## **Applications**

De novo analysis enables the identification of:
- Previously unknown genes
- Genomic structures
- Large-scale sequence differences not detectable using reference-based methods

It is especially important for:
- Studying non-model organisms
- Newly discovered species
- Highly diverged strains
- Genomes with substantial structural rearrangements
- Metagenomic samples (see Use Case 3)

Common applications:
- Sequencing new species
- Pathogen discovery

## **Popular Tools**

- **SPAdes** - Versatile assembler for bacterial genomes
- **Velvet** - Classic de Bruijn graph assembler
- **MEGAHIT** - Memory-efficient metagenome assembler
- **Flye** - Long-read assembler
- **Prokka** - Rapid prokaryotic genome annotation
- **BUSCO** - Assembly quality assessment

## **Code Example: De Novo Assembly Workflow**

<details>
<summary>**Click to expand Bash workflow**</summary>

```bash
#!/bin/bash
# De novo genome assembly workflow

# Step 1: Quality control and trimming
# Check raw read quality
fastqc raw_reads_R1.fastq.gz raw_reads_R2.fastq.gz -o qc_before/

# Trim low-quality bases and adapters
trimmomatic PE -threads 8 \
    raw_reads_R1.fastq.gz raw_reads_R2.fastq.gz \
    trimmed_R1.fastq.gz unpaired_R1.fastq.gz \
    trimmed_R2.fastq.gz unpaired_R2.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Check quality after trimming
fastqc trimmed_R1.fastq.gz trimmed_R2.fastq.gz -o qc_after/

# Step 2: De novo assembly using SPAdes
spades.py \
    -1 trimmed_R1.fastq.gz \
    -2 trimmed_R2.fastq.gz \
    -o assembly_output \
    --careful \
    -t 16 \
    -m 64

# Step 3: Assess assembly quality using QUAST
quast.py assembly_output/scaffolds.fasta \
    -o assembly_quality \
    --threads 8

# Step 4: Check completeness with BUSCO
busco \
    -i assembly_output/scaffolds.fasta \
    -o busco_output \
    -m genome \
    -l bacteria_odb10 \
    --cpu 8

# Step 5: Annotate genome using Prokka (for bacterial genomes)
prokka \
    --outdir annotation_output \
    --prefix my_genome \
    --kingdom Bacteria \
    --cpus 8 \
    assembly_output/scaffolds.fasta

echo "De novo assembly and annotation complete!"
echo "Assembly: assembly_output/scaffolds.fasta"
echo "Annotation: annotation_output/my_genome.gff"
```

</details>

## **Python Example: Assembly Quality Assessment**

<details>
<summary>**Click to expand Python script**</summary>

```python
#!/usr/bin/env python3
"""
Analyze assembly statistics from FASTA file
"""

def parse_fasta(fasta_file):
    """
    Parse FASTA file and return sequences
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        Dictionary of sequence_id: sequence
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Start new sequence
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def calculate_n50(lengths):
    """
    Calculate N50 statistic
    
    Args:
        lengths: List of contig lengths
        
    Returns:
        N50 value
    """
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target = total_length / 2
    
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= target:
            return length
    return 0

def assembly_statistics(fasta_file):
    """
    Calculate comprehensive assembly statistics
    
    Args:
        fasta_file: Path to assembly FASTA file
    """
    sequences = parse_fasta(fasta_file)
    lengths = [len(seq) for seq in sequences.values()]
    
    # Calculate statistics
    num_contigs = len(lengths)
    total_length = sum(lengths)
    mean_length = total_length / num_contigs if num_contigs > 0 else 0
    longest_contig = max(lengths) if lengths else 0
    shortest_contig = min(lengths) if lengths else 0
    n50 = calculate_n50(lengths)
    
    # Calculate GC content
    all_seq = ''.join(sequences.values()).upper()
    gc_count = all_seq.count('G') + all_seq.count('C')
    gc_content = (gc_count / len(all_seq) * 100) if all_seq else 0
    
    # Print report
    print("=" * 50)
    print("Assembly Statistics")
    print("=" * 50)
    print(f"Number of contigs: {num_contigs:,}")
    print(f"Total assembly length: {total_length:,} bp")
    print(f"Mean contig length: {mean_length:,.0f} bp")
    print(f"Longest contig: {longest_contig:,} bp")
    print(f"Shortest contig: {shortest_contig:,} bp")
    print(f"N50: {n50:,} bp")
    print(f"GC content: {gc_content:.2f}%")
    print("=" * 50)
    
    # Contig length distribution
    print("\nContig length distribution:")
    bins = [500, 1000, 5000, 10000, 50000, 100000]
    for i, threshold in enumerate(bins):
        count = sum(1 for l in lengths if l >= threshold)
        print(f"  >= {threshold:,} bp: {count} contigs")

# Usage example
if __name__ == '__main__':
    assembly_statistics('assembly_output/scaffolds.fasta')
```

</details>

## **Expected Outputs**

A successful de novo assembly produces:

- **Assembly FASTA file** - Contains contigs/scaffolds with reconstructed sequences
- **Quality metrics** - QUAST report with N50, total length, largest contig
- **Completeness scores** - BUSCO results showing % of conserved genes found
- **Annotation files** - GFF3 and protein FASTA files from Prokka

### **Interpreting Assembly Quality**

| Metric | Good | Acceptable | Poor |
|--------|------|------------|------|
| **N50** (bacterial genome) | `>50 kb` | 10-50 kb | `<10 kb` |
| **# Contigs** (bacterial) | `<100` | 100-500 | `>500` |
| **BUSCO completeness** | `>95%` | 80-95% | `<80%` |
| **Total length** (bacterial) | 4-6 Mb | 3-7 Mb | Outside range |

## **Computational Requirements**

### **Resource Estimates**

| Genome Type | Reads | CPU Cores | RAM | Time |
|-------------|-------|-----------|-----|------|
| Bacterial | 100x coverage | 16 | 32-64 GB | 2-6 hours |
| Fungal | 100x coverage | 16-32 | 64-128 GB | 6-24 hours |
| Small eukaryote | 100x coverage | 32+ | 128-256 GB | 1-3 days |
| Mammalian | 50x coverage | 64+ | 500+ GB | 3-7 days |

:::tip Resource Recommendation
For genomes >100 Mb, use the **LRZ Linux Cluster** high-memory nodes. See [Infrastructure Guide](/general-knowledge/infrastructure).
:::

## **Common Issues & Troubleshooting**

:::warning Assembly Problems

**Fragmented assembly (high contig count)**
- Increase sequencing coverage (aim for 100x+)
- Use longer reads or hybrid assembly (short + long reads)
- Try different k-mer values or assemblers

**Low BUSCO completeness**
- Wrong lineage database selected
- Insufficient sequencing depth
- Contamination or mixed samples
- Organism genuinely lacks some conserved genes

**Assembly too large/small**
- Check for contamination with Kraken2
- Verify expected genome size for organism
- Look for remaining adapter sequences

**Out of memory errors**
- Reduce k-mer size
- Use memory-efficient assemblers (MEGAHIT)
- Request more RAM or use cluster resources
- Pre-filter low-quality/error-prone reads

**Assembly takes too long**
- Reduce thread count (sometimes helps)
- Pre-normalize read depth
- Use faster assemblers for draft assemblies

:::

## **Sample Data & Practice**

Practice datasets:
- **E. coli K-12** - Well-characterized bacterial genome
- **SRA database** - Search for your organism of interest
- **SPAdes test data** - Included with SPAdes installation

## **Next Steps & Related Analyses**

After assembly and annotation:
- **Comparative genomics** - Compare to related species
- **Pangenome analysis** - Identify core and accessory genes
- **Phylogenetic analysis** - Construct evolutionary trees
- **Functional analysis** - Predict metabolic capabilities
- **Reference-based analysis** - Use your assembly as a new reference

**Related Use Cases:**
- [Reference-Based Analysis](/domain-knowledge/bioinformatics/use-cases/reference-based-analysis) - Compare strains to your assembly
- [Metagenomic Analysis](/domain-knowledge/bioinformatics/use-cases/metagenomic-analysis) - For mixed samples

## **Key Considerations**

- **Sequencing depth matters** - Higher coverage (100x+) generally leads to better assemblies
- **Read length impacts quality** - Longer reads (PacBio, Nanopore) produce dramatically better assemblies
- **Computational resources** - De novo assembly is memory-intensive; plan accordingly
- **Quality assessment is critical** - Always use BUSCO and QUAST to evaluate assembly completeness
- **Multiple k-mer sizes** - Tools like SPAdes use multiple k-mer sizes for optimal results
- **Hybrid assembly** - Combining short and long reads often yields best results
