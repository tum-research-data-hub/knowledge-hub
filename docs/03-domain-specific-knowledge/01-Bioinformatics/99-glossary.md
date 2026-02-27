---
title: Glossary
slug: /domain-knowledge/bioinformatics/glossary
---

# **Bioinformatics Glossary**

This glossary provides definitions for common terms used in bioinformatics and computational biology. Understanding these terms is essential for navigating the field and communicating effectively with collaborators.

---

## **Genetics and Biology Concepts**

### **Gene**
A region of DNA that encodes a functional product, such as a protein or functional RNA.

### **Allele**
One of multiple alternative versions of a gene at a specific genomic location.

### **Genotype**
The genetic makeup of an organism, defined by the variants it carries.

### **Phenotype**
The observable traits or characteristics of an organism resulting from its genotype and environment.

### **Mutation**
A change in DNA sequence that may affect gene function.

### **Genome Editing**
The deliberate modification of DNA sequences using targeted molecular tools (e.g., CRISPR-Cas9).

---

## **Computational Context**

### **Bioinformatics**
The application of computational methods to analyze and interpret biological data.

### **Pipeline / Workflow**
An ordered sequence of computational steps that transform raw data into biological results.

### **Method**
A computational or experimental approach used to analyze biological data.

---

## **Core Sequencing and Data Terms**

### **Sequencing**
The laboratory process of determining the order of nucleotides (A, C, G, T or U) in DNA or RNA molecules.

### **Sequencing Reads (Reads)**
Short nucleotide sequences produced by a sequencing machine that represent fragments of the original DNA or RNA.

### **FASTQ File**
A standard file format that stores sequencing reads together with a quality score for each nucleotide.

### **Quality Score (Base Quality)**
A numerical estimate of how confident the sequencing machine is that a given nucleotide was read correctly. Typically encoded as [Phred quality scores](https://en.wikipedia.org/wiki/Phred_quality_score).

### **Coverage / Depth**
The average number of reads that align to (or "cover") each base in a reference sequence. Higher coverage generally leads to more confident variant calls or assemblies.

---

## **Reference-Based Analysis Terms**

### **Reference Genome**
A curated representative genome sequence used as a coordinate system for aligning and analyzing sequencing data.

### **Mapping (Alignment)**
The process of assigning sequencing reads to their most likely positions on a reference genome.

### **Mapping Quality**
A measure of how confident the algorithm is that a sequencing read was aligned to the correct genomic location.

### **SAM / BAM File**
File formats that store aligned sequencing reads and their associated quality information. SAM (Sequence Alignment/Map) is text-based, while BAM is a compressed binary version.

### **Variant**
A difference in DNA sequence between a sample and a reference genome.

### **Single Nucleotide Variant (SNV)**
A variant affecting a single nucleotide position in the genome.

### **Single Nucleotide Polymorphism (SNP)**
An SNV that is common within a population (typically >1% frequency).

### **Insertion / Deletion (Indel)**
A variant where one or more nucleotides are inserted into or deleted from the genome.

### **Variant Calling**
The computational process of identifying genetic variants from aligned sequencing reads.

### **Variant Call Format (VCF)**
A standardized file format that stores detected genetic variants along with supporting evidence and quality metrics.

---

## **De Novo Analysis and Assembly Terms**

### **De Novo Analysis**
A type of analysis that reconstructs genome sequences directly from sequencing reads without using a reference genome.

### **Genome Assembly**
The process of combining overlapping sequencing reads into longer continuous sequences.

### **Contig**
A continuous DNA sequence assembled from overlapping sequencing reads (contiguous sequence).

### **Scaffold**
A set of ordered and oriented contigs that approximate larger genomic regions, often containing gaps.

### **N50**
A statistic indicating that 50% of the total assembly length is contained in contigs of this length or longer. Higher N50 values indicate better assemblies.

### **FASTA File**
A simple text-based file format for storing nucleotide or protein sequences without quality information.

### **Genome Annotation**
The process of identifying genes and other functional elements within a genome sequence.

---

## **Metagenomics Terms**

### **Metagenomics**
The study of genetic material obtained directly from mixed communities of organisms in a single sample.

### **Metagenome**
The combined genetic content of all organisms present in a metagenomic sample.

### **Binning**
The process of grouping assembled sequences that are likely to originate from the same organism.

### **Metagenome-Assembled Genome (MAG)**
A draft genome reconstructed from metagenomic data through assembly and binning.

### **Taxonomic Profiling**
The identification and quantification of organisms present in a biological sample.

### **Functional Profiling**
The identification of genes and metabolic pathways encoded by a biological community.

### **Microbiome**
The collection of microorganisms (bacteria, fungi, viruses) living in a specific environment, such as the human gut.

---

## **Transcriptomics Terms**

### **Transcriptomics**
The study of all RNA molecules expressed in a cell, tissue, or organism at a given time.

### **RNA-seq**
A sequencing method used to measure RNA abundance and gene expression levels.

### **Transcriptome**
The complete set of RNA transcripts produced by a genome under specific conditions.

### **Gene Expression**
The process by which genetic information is used to produce RNA and proteins. In the context of RNA-seq, it refers to the abundance of RNA transcripts.

### **Expression Quantification**
The measurement of how frequently each gene or transcript is observed in sequencing data.

### **Differential Expression Analysis**
A statistical analysis that identifies genes whose expression levels differ significantly between conditions or sample groups.

### **Expression Matrix**
A table containing gene or transcript expression values across multiple samples (rows = genes, columns = samples).

### **TPM (Transcripts Per Million)**
A normalized measure of gene expression that accounts for both sequencing depth and gene length.

### **FPKM / RPKM**
Fragments/Reads Per Kilobase of transcript per Million mapped reads. Older normalization methods for RNA-seq data.

---

## **File Formats**

### **FASTA**
Text format for nucleotide or protein sequences. Each sequence starts with a header line (beginning with `>`) followed by the sequence.

### **FASTQ**
Extension of FASTA that includes quality scores for each base. Standard output format from sequencing machines.

### **SAM/BAM**
Sequence Alignment/Map format for storing aligned reads. BAM is the compressed binary version.

### **VCF/BCF**
Variant Call Format for storing genetic variants. BCF is the compressed binary version.

### **GFF3/GTF**
Gene annotation formats describing genomic features (genes, exons, etc.).

### **BED**
Browser Extensible Data format for genomic intervals and annotations.

### **HDF5**
Hierarchical Data Format for storing large, complex datasets, commonly used in single-cell sequencing.

---

## **Statistical and Analysis Terms**

### **P-value**
The probability of observing results as extreme as those obtained, assuming the null hypothesis is true.

### **Adjusted P-value (FDR, padj)**
P-value corrected for multiple testing to control the false discovery rate.

### **Log2 Fold Change**
A measure of how much a gene's expression has changed between conditions, on a log2 scale. A value of 1 means doubled expression, -1 means halved.

### **Principal Component Analysis (PCA)**
A dimensionality reduction technique used to visualize similarities and differences between samples.

### **Clustering**
Grouping samples or genes with similar characteristics or expression patterns.

---

## **Ontologies and Databases**

### **Ontology**
A controlled vocabulary of terms and their relationships, used to standardize biological knowledge.

### **Gene Ontology (GO)**
A standardized system for describing gene function, cellular component, and biological process.

### **BLAST**
Basic Local Alignment Search Tool - algorithm for comparing sequences to find regions of similarity.

### **NCBI**
National Center for Biotechnology Information - major database of biological sequences and literature.

### **UniProt**
Universal Protein Resource - comprehensive protein sequence and functional information database.

### **Ensembl**
Genome browser and annotation database for vertebrate and other eukaryotic genomes.

---

## **Quality Control Terms**

### **PHRED Score**
Quality score representing the probability that a base call is incorrect. Score of 20 = 1% error rate, 30 = 0.1% error rate.

### **Adapter Contamination**
Sequencing artifacts where adapter sequences from library preparation are present in reads.

### **PCR Duplicates**
Identical reads arising from PCR amplification during library preparation rather than independent sequencing events.

### **GC Content**
The percentage of guanine (G) and cytosine (C) bases in a sequence, which can affect sequencing quality.

---

## **Additional Resources**

For more detailed information on bioinformatics terminology:
- [NCBI Glossary](https://www.ncbi.nlm.nih.gov/glossary/)
- [EMBL-EBI Training](https://www.ebi.ac.uk/training/)
- [Galaxy Training Materials](https://training.galaxyproject.org/)

---

**Need clarification on a term?** Contact us at [rdhub@mdsi.tum.de](mailto:rdhub@mdsi.tum.de) or explore our [Use Cases](/domain-knowledge/bioinformatics/use-cases/reference-based-analysis) for practical examples.
