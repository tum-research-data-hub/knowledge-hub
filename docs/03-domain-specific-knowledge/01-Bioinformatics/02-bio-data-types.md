---
 title: Data Types 
 slug: /domain-knowledge/bioinformatics/bio-data-types
---

# **Common Bioinformatics Data Types and Formats**

Understanding the various file formats used in bioinformatics is essential for effective data management, analysis, and reproducibility. This guide covers the standard formats for sequences, alignments, variants, and genomic features.

## **🧬 Sequence Formats**

FASTA and FASTQ files are the foundational formats for next-generation sequencing (NGS). They represent the starting point for most bioinformatics pipelines, carrying the raw or processed code of life.

### **FASTA**

A [**FASTA**](https://en.wikipedia.org/wiki/FASTA_format) file can contain multiple nucleotide or amino acid sequences. Each entry starts with a header line beginning with the \> symbol, followed by a unique identifier and metadata (e.g., gene name, organism). The subsequent lines contain the actual sequence data.

* **Note:** FASTA files typically do not contain information about the quality of the sequence.

### **FASTQ**

A **[FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)** file is an extension of the FASTA format that also stores corresponding quality scores. Each read entry is recorded in four lines:

1. Sequence identifier  
2. Nucleotide sequence  
3. A \+ delimiter  
4. Quality scores (encoded as ASCII characters) indicating the confidence of each nucleotide call.

## **📍 Alignment Formats**

After sequencing, a major task is "mapping" or aligning these reads to a reference genome to determine their origin.

### **SAM (Sequence Alignment Map)**

The [**Sequence Alignment Map (SAM)**](https://en.wikipedia.org/wiki/SAM_(file_format)) format is generic text format for storing the alignments of sequencing reads against reference nucleotide sequences and is typically the output of read alignment algorithms, such as Bowtie2. Each file contains one header section and alignment section with the read sequence and quality scores such as the reference chromosome, genomic position, mapping quality, and a CIGAR string describing matches/insertions/deletions. The SAM format supports short and long reads (up to 128 Mbp) produced by different sequencing platforms.

### **BAM (Binary Alignment Map)**

The [**Binary Alignment Map (BAM)**](https://en.wikipedia.org/wiki/BAM_(file_format)) format is a compressed, binary format of SAM and thus takes up much less storage space: a BAM file is typically around a quarter to a fifth of the size of the uncompressed SAM file, depending on the data. Some special tools are needed in order to make sense of BAM, such as Samtools, Picard, and Integrative Genomics Viewer (IGV).

### **CRAM**

The **CRAM** format is a further-compressed alternative to BAM. Instead of storing every read sequence in full, it stores only the differences to the reference genome (reference-based compression), which typically makes CRAM files 30-60% smaller than the corresponding BAM file. CRAM is the preferred submission format of the European Nucleotide Archive (ENA) for aligned reads and is well suited for archiving.

* **Note:** Decoding a CRAM file requires the exact reference sequence it was written against. Record the reference (name, version, checksum) alongside the file and keep it available; otherwise the reads cannot be restored. Samtools and htslib read and write CRAM (`samtools view -C -T reference.fasta`).

The detailed format specifications of the SAM/BAM/CRAM formats can be found on the official web site: https://samtools.github.io/hts-specs/. 

### **STO (Stockholm)**


The **STO (Stockholm)** format is a text file format for storing multiple sequence alignments (MSAs) and additional alignment features and metadata. They typically have .sto or .stk file extensions. The Stockholm format is used by the Pfam and Rfam databases and tools like HMMER. A STO filE consists of a header line with a format and version identifier; mark-up lines starting with "#=GF","#=GC","#=GS" or "#=GR"; alignment lines with the sequence name and aligned sequence; a "//" line indicating the end of the alignment. Alignments are shown with inserts as lower case characters, matches as upper case characters, and gaps as ' . ' or ' - '. 

A more detailed description of Stockholm files can be found here: https://rfam.xfam.org/help#tabview=tab13.


## **🔍 Variants**

### **VCF (Variant Call Format)**

The **VCF (Variant Call Format)** is a tab-delimited text file format originally developed under the 1000 Genomes Project and is the standard output of variant callers. It is used to store information sequence variations in the genome, such as Single Nucleotide Polymorphisms (SNPs) and insertions/deletions (indels). It contains meta-information lines (prefixed with ##), a header line (prefixed with #), and data lines each containing 8 fields for each variant in the genome, typically followed by additional columns showing format and sample information. 

### **BCF (Binary Call Format)**

The **Binary Call Format (BCF)** is the binary version of the VCF format, allowing for faster processing, reading, and storage for large projects. Tools like BCFtools can work with both VCF and BCF files and are often used to convert between them. 

More detailed information on the VCF and BCF files can be found here: https://samtools.github.io/hts-specs/VCFv4.2.pdf.

## **📊 Counts and Expression Matrices**

For transcriptomics (RNA-seq), data is often represented as an expression matrix—a table where rows represent genes and columns represent samples or individual cells.

* **Text Files (TSV/CSV):** Simple tabular formats used for smaller datasets.  
* The **Hierarchal Data File (HDF5)** is a file format designed to store large amounts of omics data and metadata. HDF5 uses a hierarchical directory structure that act as containers to store two types of objects: 1) groups, containing zero or more groups; and 2) datasets, multi-dimensional arrays of read counts together with supporting metadata. This format allows for storing both the count matrices and all metadata in a single file rather than having separate features, barcodes and matrix files. This makes it suitable for single-cell sequencing data.

## **🗺️ Gene Feature Formats**

These formats define the "map" of the genome, indicating where genes, exons, and introns are located.

* The **General Feature Format (GFF)** format is a standard, general-purpose format for storing genomic annotations. A GFF file generally contains the name of the sequence, the type of feature described (i.e., a whole gene, an exon or an intron), the coordinates of the feature within the overall sequence and other additional attributes. It contains genes and other features of DNA, RNA and protein sequences, each line representing a single feature and containing nine fields: seqid, source, type, start, end, score, strand, phase, and attributes. The most recent version of GFF is GFF3. 
	The official documentation for the GFF format can be found here: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

* **GTF (Gene Transfer Format):** A variation of GFF specifically focused on gene structures.  

* The **Browser Extensible Data (BED)** is a tab-separated file used to store basic genomic features. BED files are typically generated by feature detection algorithms or through manual curation of alignments. Each entry consists of one line per feature, with the first three required columns specifying the genomic region interval (contig/chromsome name, begin, and end position); the remaining nine columns contain additional track information. 

    The full BED specification can be found here: http://genome.ucsc.edu/FAQ/FAQformat.html#format1.

For how these formats are produced and quality-checked, see the [Gene & Genome Annotation](/domain-knowledge/bioinformatics/gene-annotation) guide.

## **⚖️ Mass Spectrometry Formats**

Proteomics and metabolomics instruments write vendor-specific raw files (e.g. Thermo `.raw`, Bruker `.d`) that can only be read with the vendor's software libraries. Open, XML-based community formats maintained by the [HUPO Proteomics Standards Initiative (PSI)](https://www.psidev.info/) make the data readable independent of the vendor.

### **mzML**

**mzML** is the open standard format for mass spectrometry raw data: the spectra (m/z and intensity arrays), chromatograms, and instrument and acquisition metadata described with a controlled vocabulary. Vendor files are usually converted with `msconvert` from [ProteoWizard](https://proteowizard.sourceforge.io/). Repositories such as PRIDE and MetaboLights accept mzML alongside (not instead of) the original raw files, so keep both.

Related PSI formats cover the later analysis steps: **mzIdentML** (peptide and protein identifications) and **mzTab** (tabular summary of identification and quantification results).

The mzML specification can be found here: https://www.psidev.info/mzML.

## **🧊 Macromolecular Structure Formats**

### **mmCIF (PDBx/mmCIF)**

**mmCIF** (macromolecular Crystallographic Information File) is the standard format for three-dimensional structures of proteins and nucleic acids: atomic coordinates together with experimental and model metadata, organised as named categories and key-value items. It is the master archive format of the worldwide Protein Data Bank (wwPDB) and the mandatory deposition format for crystallographic structures. It is also used for predicted models (e.g. AlphaFold DB, via the ModelCIF extension).

* **Note:** The older, fixed-column **PDB format** (`.pdb`) is still widely used by tools but is frozen and cannot represent large structures (more than 99,999 atoms or 62 chains). Use mmCIF for anything you archive or deposit.

The PDBx/mmCIF dictionary and documentation can be found here: https://mmcif.wwpdb.org/.
