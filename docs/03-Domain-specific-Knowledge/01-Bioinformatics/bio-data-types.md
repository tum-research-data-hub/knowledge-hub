---
 title: Data Types 
 sidebar\_label: Data Types & Formats
---

# **Common Bioinformatics Data Types and Formats**

Understanding the various file formats used in bioinformatics is essential for effective data management, analysis, and reproducibility. This guide covers the standard formats for sequences, alignments, variants, and genomic features.

## **🧬 Sequence Formats**

FASTA and FASTQ files are the foundational formats for next-generation sequencing (NGS). They represent the starting point for most bioinformatics pipelines, carrying the raw or processed code of life.

### **FASTA**

A **FASTA** file can contain multiple nucleotide or amino acid sequences. Each entry starts with a header line beginning with the \> symbol, followed by a unique identifier and metadata (e.g., gene name, organism). The subsequent lines contain the actual sequence data.

* **Note:** FASTA files typically do not contain information about the quality of the sequence.

### **FASTQ**

A **FASTQ** file is an extension of the FASTA format that also stores corresponding quality scores. Each read entry is recorded in four lines:

1. Sequence identifier  
2. Nucleotide sequence  
3. A \+ delimiter  
4. Quality scores (encoded as ASCII characters) indicating the confidence of each nucleotide call.

## **📍 Alignment Formats**

After sequencing, a major task is "mapping" or aligning these reads to a reference genome to determine their origin.

### **SAM (Sequence Alignment Map)**

The **SAM** format is a generic, human-readable text format for storing read alignments. It includes metadata such as genomic position, mapping quality, and a **CIGAR** string describing matches, insertions, and deletions.

### **BAM (Binary Alignment Map)**

The **BAM** format is a compressed, binary version of SAM. It is significantly smaller (30-50%) and is the standard for high-performance computing, though it requires specific tools like samtools or Picard to read.

### **STO (Stockholm)**

The **Stockholm** format is used for multiple sequence alignments (MSAs) and is the standard for databases like Pfam and Rfam. It includes detailed mark-up for sequence features and consensus structures.

## **🔍 Variants**

### **VCF (Variant Call Format)**

**VCF** is the standard for storing genome sequence variations, such as **SNPs** (Single Nucleotide Polymorphisms) and **indels**. It consists of a large header describing the data, followed by records for every position where a variation was found.

### **BCF (Binary Call Format)**

**BCF** is the binary counterpart to VCF. It allows for much faster processing and is highly recommended for large-scale population studies.

## **📊 Counts and Expression Matrices**

For transcriptomics (RNA-seq), data is often represented as an expression matrix—a table where rows represent genes and columns represent samples or individual cells.

* **Text Files (TSV/CSV):** Simple tabular formats used for smaller datasets.  
* **HDF5 (Hierarchical Data Format):** A container format designed for massive amounts of data. It is commonly used for single-cell RNA-seq because it can store count matrices and metadata in a single, fast-access file.

## **🗺️ Gene Feature Formats**

These formats define the "map" of the genome, indicating where genes, exons, and introns are located.

* **GFF3 (General Feature Format):** A robust 9-column format for storing genomic annotations.  
* **GTF (Gene Transfer Format):** A variation of GFF specifically focused on gene structures.  
* **BED (Browser Extensible Data):** A flexible, tab-separated format used primarily for defining genomic intervals and visualization on genome browsers.

## **🔗 External Specifications**

* [SAM/BAM/VCF/BCF Specifications](https://samtools.github.io/hts-specs/)  
* [GFF3 Documentation](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)  
* [UCSC BED Format Guide](http://genome.ucsc.edu/FAQ/FAQformat.html#format1)