# TUM Data Knowledge Hub
By the TUM Research Data Hub


# 01 About TUM Data Knowledge Hub



# 02 General Knowledge


## Data Types
Depending on your research questions, you may encounter various categories of data.
Data type describes the kind of data and defines how it is stored, what operations can be applied to it, and what analysis or models are valid.
Understanding data types is essential because it affects every step of your entire research process, and data misclassification may cause errored results and flawed conclusions.
The most common data types despite the research disciplines include:
Numeric data (including integers and float numbers)
Numeric data refer to quantities that can be measured or counted.
Numeric data is the foundation of quantitative analysis, storing numbers in the correct format is crucial for valid modelling and calculations.

Categorical data
Categorial data represent values that fall into groups or categories without any specific numerical meaning.
Example: experimental conditions group
‘High’
‘Medium’
‘Low’

Ordinal data
Ordinal data are similar to categorical data but possess an inherent order or ranking.
Example: ordered survey responses
‘Definitely agree’ --> 5
‘Mostly agree’ --> 4
‘Neither agree nor disagree’ --> 3
‘Mostly disagree’ --> 2
‘Definitely disagree’ --> 1

Text or String data
Text or string data typically includes free-format text, names, labels, and descriptions. Text is also widely used for metadata and documentation.

Binary or Boolean data
Binary or Boolean data refer to yes/no or true/false values.
Booleans are also used in filtering and subsetting datasets, building logical conditions in programming, creating binary classification models, and documenting quality-control checks.
Example: ‘Did the robot detect an obstacle?’ --> False

Date and Time data
Date and time data record temporal information and they are crucial for examining trends, sequences, intervals, or rhythms. Proper handling of temporal data is essential especially when raw data are international date formats and across time zones.
Example: Timestamps
2025-01-03 10:15:00,
2025-01-03 10:15:01, ...

Depending on your research subject and area, you might also see special data types, including:
Matrix and Vector
Matrix and vector data are structure collection of numeric data. Vector is one-dimensional array whereas matrices are two-dimensional.

Geospatial data
Geospatial data refers to the information that is linked to a specific location on or near the surface of the Earth. It includes location information (e.g. geographic coordinates), attribute information (characteristics or details about the location) and temporal information (time or lifespan of the data).

Image data

Genomic data
Link down to the domain-specific data type

This is by means a complete list.




## Repositories
## Tools/Programming languages
R (Gagneur book: https://gagneurlab.github.io/dataviz/)
Python
pip
MatLAB
SPSS
Excel
Minitab?
Power BI
Microsoft Power BI is popular visualisation tool to build interactive reports and dashboards.  It enables connection with enterprise data warehouses with live data refresh. Power BI has been widely used for organisations to discover business intelligence and empower strategic decision making.
Other similar tools: Tableau,
Columns:
Name
Usage
Type
Description
Links (internal)
Links (documentation)


Home Software - Software - BayernCollab
SQL



## Infrastructure
This page will give you an overview of the different types of computing and storage infrastructure available to researchers at the Technical University of Munich.
At the end of this page, you will find a decision guide to help you decide on a solution, also going into the different specifications of computing infrastructure (CPU, RAM, GPU) and the type of limitations you can face depending on your requirements.
Domain-specific guides can be found in their specific sub-chapters and may be linked on this site.
LRZ
Linux Cluster
Compute Cloud
AI Systems Cluster
Quantum Computing
Other GPU options
Hosting own hardware
Großgeräteanträge (Major Research Instrumentation)
Terrabyte
EOSC
EU Node
NFDI
Jupyter4NFDI
NFDI4Ing
NHR (National High Performance Computing Alliance)
Application by Doctoral Researcher
AWS
MS Azure
Google Cloud Computing
EuroHPC https://www.eurohpc-ju.europa.eu/index_en
“your own computer”
WSL
Seeing specifications
Limitations
Columns:
Name
Link
Costs
Type of Service (Cluster, VM, Jupyter, etc.)
Access requirement / Application Procedure
Support
Intended Use?
Storage/Transfer Options
Additional Information:
Decision help / How to decide?
Also limiting factors (RAM, CPUs, storage)
Recommendation for LRZ services, incl. Feedback options and Consultation offer by the Servicedesk

## Metadata
Metadata refers to the information to describe about a data point or a dataset. Metadata are used to describe the structure, content and features. It is crucial to help systems and users to organise, manage and use data more effectively.
For example, the metadata of an image file from your phone might include file format, file size, resolution, resolution, date created, etc.


# 03 Domain-specific Knowledge
## Bioinformatics/Biology
TODO: make into a table
TODO: Include a workflow/diagram

### Data Types
Sequencing counts
HDF5 format
Text files (tsv, csv)
Excel files
Sequence formats
FASTA
FASTQ
Alignment formats
SAM
BAM
Stockholm
Variants
VCF
BCF
Gene feature formats
GFF3
BED
GTF files

SEQUENCE formats
FASTA and FASTQ files are the standard output formats used by next-generation sequencing platforms, such as Illumina sequencers. <ADD MORE INTRO/MOTIVATION>.
A FASTA file is text-based format for representing nucleotide sequences or amino acid sequences, in which base pairs or amino acids are represented using single-letter codes. A FASTA file can contain multiple nucleotide or amino acid sequences. Each entry starts with a header line starting with the “>” symbol, followed by other information about the read, typically starting with a unique identifier, name of the gene, the type of sequence, and the organism or sample where the read was sequenced.  The second line of each entry contains the nucleotide or amino acid sequence. The entire sequence can be on a single long line but is often split into multiple lines of 60-100 characters each. FASTA entries typically do not contain information about the quality of the sequence.

A FASTQ file is a modified version of a traditional FASTA format, additionally storing the corresponding quality scores of nucleotide sequences. FASTQ files are commonly used to store raw sequencing reads, short fragments produced by the sequencing machine. Each read entry is recorded in four lines: a sequence identifier, the nucleotide sequence, a `+` delimiter line, and quality scores that indicate the confidence of each nucleotide call, each encoded with a single ASCII character. The per-base quality scores are important because during sequencing, some nucleotides are more likely to be incorrect than others, which could affect downstream analyses, such as variant calling or read mapping. In such cases, it may be necessary to trim or filter out low-quality sequencing reads.

ALIGNMENT FORMATS
A major task in bioinformatics is aligning reads from a sequencing machine to a reference genome. <CONTINUE MOTIVATION/INTRO>
The Sequence Alignment Map (SAM) format is generic text format for storing the alignments of sequencing reads against reference nucleotide sequences and is typically the output of read alignment algorithms, such as Bowtie2. Each file contains one header section (each header line starting with @) and alignment section, with each alignment having 11 mandatory fields containing the read sequence, reference chromosome, genomic position, mapping quality, and a CIGAR string describing matches/insertions/deletions. The SAM format supports both short and long reads (up to 128 Mbp) produced by various sequencing platforms.
The Binary Alignment Map (BAM) format is a compressed, binary format of SAM and thus take up less storage space than SAM files (usually 30-50% smaller). Some special tools are needed in order to make sense of BAM, such as Samtools, Picard, and IGV.
The detailed format specification of the SAM/BAM formats can be found on the official website: https://samtools.github.io/hts-specs/.

The STO (Stockholm) format is a text file format for storing multiple sequence alignments (MSAs) and additional alignment features and metadata (typically .sto or .stk file extensions). The STO format is used by the Pfam and Rfam databases and tools like HMMER. A STO filE consists of a header line with a format and version identifier; mark-up lines starting with "#=GF","#=GC","#=GS" or "#=GR"; alignment lines with the sequence name and aligned sequence; a "//" line indicating the end of the alignment. Alignments are shown with inserts as lower case characters, matches as upper case characters, and gaps as ' . ' or ' - '.  A more detailed description of Stockholm files can be found here: https://rfam.xfam.org/help#tabview=tab13.

VARIANTS
The VCF (Variant Call Format) is a tab-delimited text file format originally developed under the 1000 Genomes Project and is the standard output of variant callers. It is used to store information sequence variations in the genome, such as Single Nucleotide Polymorphisms (SNPs) and insertions/deletions (indels). It contains meta-information lines (prefixed with ##), a header line (prefixed with #), and data lines each containing 8 fields for each variant in the genome, typically followed by additional columns showing format and sample information.
The Binary Call Format (BCF) is the binary version of the VCF format, allowing for faster processing, reading, and storage for large projects. Tools like BCFtools can work with both VCF and BCF files and are often used to convert between them. A more detailed specification for the VCF and BCF files can be found here: https://samtools.github.io/hts-specs/VCFv4.2.pdf.

COUNTS
For example, aAgene expression matrix is a table where the rows represent the genes, the columns represent the samples, and the values in the cells represent the counts or expression levels of the genes.

TEXT files
The text format consists of tab or comma separated columns with genes on the columns and cells on the rows.  Gene expression matrices and their associated metadata are commonly stored in tabular formats, such as using TSV, CSV or TXT formats.
The Hierarchal Data File (HDF5) is a file format designed to store large amounts of omics data and metadata. HDF5 uses a hierarchical directory structure that act as containers to store two types of objects: 1) groups, containing zero or more groups; and 2) datasets, multi-dimensional arrays of read counts together with supporting metadata. This format allows for storing both the count matrices and all metadata in a single file rather than having separate features, barcodes and matrix files.  This makes it suitable for single-cell sequencing data, which are based on family of scRNA-seq storage formats

GENE FEATURE FORMATS

The General Feature Format (GFF) format is a standard, general-purpose format for storing genomic annotations. A GFF file generally contains the name of the sequence, the type of feature described (i.e., a whole gene, an exon or an intron), the coordinates of the feature within the overall sequence and other additional attributes. It  genes and other features of DNA, RNA and protein sequences, each line representing a single feature and containing nine fields: seqid, source, type, start, end, score, strand, phase, and attributes. The most recent version of GFF is GFF3. The official documentation can be found here: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
The Gene Transfer Format (GTF) format borrows from the GFF file format. For a detailed descripttion of the fields, please see http://gmod.org/wiki/GFF3.
The Browser Extensible Data (BED) is a tab-separated file used to store basic genomic features. BED files are typically generated by feature detection algorithms or through manual curation of alignments. Each entry consists of one line per feature, with the first three required columns specifying the genomic region interval (contig/chromsome name, begin, and end position); the remaining nine columns contain additional track information. BED files can be manipuluated using standard Unix tools (e.g. sed, awk, and sort) or dedicated packages, such as the bedtools suite with additional functionality.
The full BED specification can be found here:
http://genome.ucsc.edu/FAQ/FAQformat.html#format1.

### Repositories

Data centers and repositories are common tools for making research data accessible. In repositories, the data are enriched with metadata containing relevant information for findability and subsequent use. This makes the data findable in the catalog of the respective repository, but also in metadata search portals and in data search engines. Repositories are therefore both search and storage locations for accessible data. Here are some useful links for finding trustworthy research data repositories:
The Registry of Research Data Repositories (re3data): https://www.re3data.org/
Global registry of research data repositories across all academic disciplines. Although not a search engine for finding data, this resource is a starting point for finding repositories that may contain relevant searchable data.
https://fairsharing.org/: - community-maintained portal that allows searching for data and metadata standards, as well as repositories and policies
DataCite’s Repository finder: https://repositoryfinder.datacite.org/
DataCite Commons indexes datasets registered with a DOI with DataCite, CrossRef or one of six other scholarly DOI agencies. DataCite links the repository to information in re3data. The DataCite Commons support page provides guidance on how to search effectively.

Biology-related
The Expasy (Expert Protein Analysis System) portal from the SIB Swiss Institute of Bioinformatics (SIB), providing free access to over 160 databases, software tools, and datasets for life sciences research, especially in proteomics, genomics, and structural biology. It's a central hub connecting resources like SWISS-PROT (protein database) and PROSITE (protein patterns) with analysis tools for tasks like sequence alignment, structure prediction, and phylogeny, supporting researchers with in-depth protein knowledge and biological data analysis.
NIAID https://data.niaid.nih.gov/
Can be used to find datasets on Infectious and Immune-mediated Diseases (IID) across many repositories. Provides Resource Catalogs (collections of scientific information or research outputs) and Dataset Repositories (collections of data of a particular experimental type)
DatabaseCommons is a global catalog of biological databases. It provides easy access and retrieval to a full collection of worldwide biological databases, assesses the database impact by factoring both citation and age, and delivers a series of useful statistics and trends to investigate their status and impact on biomedical research

To help get you started, we provide a list of some of the major databases used in biology and bioinformatics. Note that this is a non-exhaustive list and there are many other domain-specific databases.

Omics and Sequencing repositories
NCBI Sequence Read Archive (SRA)
DNA Data Bank of Japan (DDBJ )
EMBL European Nucleotide Archive (ENA)
EBI Metagenomics
Gene Expression Omnibus (GEO)
ArrayExpress



### Tools
Programming languages
Software managers
Conda
Pip
Bioconductor
BiomaRt
### Infrastructure
Sequera (Domain-specific)
De.NBI  (maybe only in domain-specific)
Things to include:
General Use-cases and their required computing specifications
Things to look out for (bulletpoints)

### Metadata & Ontologies

Researchers need specific data in order to answer specific questions. Standardized metadata are important for enabling search based on matches to specific terms across datasets. <IMPROVE INTRO/MOTIVATION>

Ontologies are a controlled set of terms and structured relationships between these terms. Ontologies represent the shared background knowledge for a community and are useful for describing data. One of the most basic forms of an ontology is a hierarchical classification, which uses an “is a” relationship between terms. Sample metadata can also be structured to use ontologies. The structured relationships provided by ontologies make our metadata more machine readable.

Gene Ontology
The Gene Ontology is a structured, standardized representation of biological knowledge designed to describe functional characteristics of gene products across the tree of life. The GO terms are linked together by relations into a labeled directed acyclic graph.
GO is organized into three main aspects used to describe gene function:
Molecular function (MF): the activities performed by a gene product at the molecular level.
Cellular component (CC): the locations, relative to cellular structures, where MFs are performed.
Biological process (BP): a “biological program” comprising molecular activities acting in concert to achieve a particular outcome; this program can be at the cellular level or at the organism level of multicellular organisms.
Link: https://geneontology.org/

Human Phenotype Ontology (HPO)

The Human Phenotype Ontology (HPO) is a standardized vocabulary for describing human abnormalities and symptoms (phenotypes) associated with a disease (e.g., "muscle weakness," "atrial septal defect"). Developed from literature and resources like Orphanet and OMIM, HPO helps clinicians and researchers by providing consistent terms for data, enabling deep phenotpying, phenotype-driven diagnosis and integrating genomics for rare disease understanding, serving as a key part of global health initiatives like the GA4GH.

Link: https://hpo.jax.org/

Human Disease Ontology (DO)
The Human Disease Ontology (DO) is a well-established ontology of human diseases useful for disease classification. The DO uses disease-specific terms and identifiers from biomedical vocabularies, such as ICD, SNOMED CT, OMIM, MeSH and the NCBI Thesaurus. DO provides a consistent, reusable, and sustainable description of human disease terms to enable data integration and analysis across different biomedical resources and genomic databases.
Link: https://disease-ontology.org/about/

Plant Ontology (PO)
The Plant Ontology is a structured vocabulary developed by the Plant Ontology Consortium that describes plant anatomy, morphology growth and development. The PO is intended for annotating gene expression and phenotype data, plant structures and stages of plant development, using the data model adopted by the Gene Ontology.
Link: http://www.plantontology.org/.

Ontology of Host-Microbiome Interactions (OHMI)

The Ontology of Host-Microbiome Interactions aims to ontologically represent and standardize various entities and relations related to microbiomes, microbiome host organisms (e.g., human and mouse), and the interactions between the hosts and microbiomes across different conditions. OHMI incorporates established ontologies to create logically structured representations of (1) microbiomes, microbial taxonomy, host species, host anatomical entities, and HMIs under different conditions and (2) associated study protocols and types of data analysis and experimental results. OHMI comprises over 1000 terms, including terms imported from more than 10 existing ontologies together with some 500 OHMI-specific terms.
Link: https://github.com/ohmi-ontology/ohmi

Chemical Entities of Biological Interest (ChEBI)
The Chemical Entities of Biological Interest (ChEBI) is a freely accessible and manually curated database and ontology of molecular entities, primarily small chemical compounds. The chemical entities in ChEBI are either naturally occurring molecules or synthetic compounds that interfere with biological processes . The database covers the chemical nomenclatures, structures, and ontology of over 195,000 entries. The nomenclature, symbolism and terminology endorsed by the International Union of Pure and Applied Chemistry (IUPAC) and the Nomenclature Committee of the International Union of Biochemistry and Molecular Biology (NC-IUBMB).
Link: https://www.ebi.ac.uk/chebi/

Other Resources for Ontologies
More ontologies can be found here:
BioPortal: https://bioportal.bioontology.org/
https://www.ebi.ac.uk/ols4/ontologies

# 04 Use-cases
## Computing Infrastructure
Who needs computing?
Individual Researcher
Institutional
Private entity
Research Group
Students
Chair
Research Project / Collaboration
Type of Grant
National/International
Hardware requirements?
CPU
Cores / Threads
Clock Speed
Cache Size
GPU
Core Count
Memory bandwith
VRAM
Memory
Capacity
Clock Speed
Latency
What type of need?
Discipline?
Link to individual Usecases and “Starting Points”
Application
AI
Data Analysis
Visualization
Data Protection (Medical / Human data)
Storage size / Capacity
Possibilites for Sharing


## Bioinformatics/Biology
### Central terms in bioinformatics:
The use-cases on the following pages include several terms specific to biological and bioinformatics, which we want to explain shortly here for reference.
#### Genetics and biology concepts
Gene
 A region of DNA that encodes a functional product, such as a protein or functional RNA.
Allele
 One of multiple alternative versions of a gene at a specific genomic location.
Genotype
 The genetic makeup of an organism, defined by the variants it carries.
Phenotype
 The observable traits or characteristics of an organism resulting from its genotype and environment.
Mutation
 A change in DNA sequence that may affect gene function.
Genome editing
 The deliberate modification of DNA sequences using targeted molecular tools.
#### Computational context
Bioinformatics
 The application of computational methods to analyze and interpret biological data.
Pipeline / Workflow
 An ordered sequence of computational steps that transform raw data into biological results.
Method
 A computational or experimental approach used to analyze biological data.
#### Core sequencing and data terms
Sequencing
 The laboratory process of determining the order of nucleotides (A, C, G, T or U) in DNA or RNA molecules.
Sequencing reads (Reads)
 Short nucleotide sequences produced by a sequencing machine that represent fragments of the original DNA or RNA.
FASTQ file
 A standard file format that stores sequencing reads together with a quality score for each nucleotide.
Quality score (Base quality)
 A numerical estimate of how confident the sequencing machine is that a given nucleotide was read correctly.
#### Reference-based analysis terms
Reference genome
 A curated representative genome sequence used as a coordinate system for aligning and analyzing sequencing data.
Mapping (Alignment)
 The process of assigning sequencing reads to their most likely positions on a reference genome.
Mapping quality
 A measure of how confident the algorithm is that a sequencing read was aligned to the correct genomic location.
SAM / BAM file
 File formats that store aligned sequencing reads and their associated quality information, with BAM being a compressed binary version of SAM.
Variant
 A difference in DNA sequence between a sample and a reference genome.
Single Nucleotide Variant (SNV)
 A variant affecting a single nucleotide position in the genome.
Single Nucleotide Polymorphism (SNP)
 An SNV that is common within a population.
Insertion / Deletion (Indel)
 A variant where one or more nucleotides are inserted into or deleted from the genome.
Variant calling
 The computational process of identifying genetic variants from aligned sequencing reads.
Variant Call Format (VCF)
 A standardized file format that stores detected genetic variants along with supporting evidence and quality metrics.
#### De novo analysis and assembly terms
De novo analysis
 A type of analysis that reconstructs genome sequences directly from sequencing reads without using a reference genome.
Genome assembly
 The process of combining overlapping sequencing reads into longer continuous sequences.
Contig
 A continuous DNA sequence assembled from overlapping sequencing reads.
Scaffold
 A set of ordered and oriented contigs that approximate larger genomic regions, often containing gaps.
FASTA file
 A simple file format for storing nucleotide or protein sequences without quality information.
Genome annotation
 The process of identifying genes and other functional elements within a genome sequence.
#### Metagenomics terms
Metagenomics
 The study of genetic material obtained directly from mixed communities of organisms in a single sample.
Metagenome
 The combined genetic content of all organisms present in a metagenomic sample.
Binning
 The process of grouping assembled sequences that are likely to originate from the same organism.
Metagenome-assembled genome (MAG)
 A draft genome reconstructed from metagenomic data through assembly and binning.
Taxonomic profiling
 The identification and quantification of organisms present in a biological sample.
Functional profiling
 The identification of genes and metabolic pathways encoded by a biological community.
#### Transcriptomics terms
Transcriptomics
 The study of all RNA molecules expressed in a cell, tissue, or organism at a given time.
RNA-seq
 A sequencing method used to measure RNA abundance and gene expression levels.
Transcriptome
 The complete set of RNA transcripts produced by a genome under specific conditions.
Gene expression
 The process by which genetic information is used to produce RNA and proteins.
Expression quantification
 The measurement of how frequently each gene or transcript is observed in sequencing data.
Differential expression analysis
 A statistical analysis that identifies genes whose expression levels differ between conditions or sample groups.
Expression matrix
 A table containing gene or transcript expression values across multiple samples.

## Bioinformatics use-cases
### Use-case 1: Reference-based sequence analysis
This type of analysis starts with DNA sequencing data in the form of “sequencing reads”, typically stored in FASTQ-files (link), which contain both nucleotide sequences and theit per-base quality scores. These reads are then mapped (aligned) to a reference genome, a curated, representative genome sequence that serves as a coordinate system for the organism under study (rather link to a “central terms” section?).
The mapping process produces an alignment file, commonly in SAM (or its binary equivalent, BAM) format, which records the individual reads and their positions on the reference genome. In addition to the aligned sequences, these files store two kinds of quality measurements: the original base quality scores from the sequencing and a new mapping quality score, reflecting the confidence of a read being placed at the correct genomic location.
Using this information, differences between the sequenced sample and the reference genome can be identified. These differences may include naturally occuring variation as well as deliberately introduced genetic changes, such as targeted mutations or deletions genereted by genome editing experiments.
Such changes in the genetic code include single-nucleotide substitutions(Single Nucleotide Variants, SNVs; or when common in a population “Single Nucleotide Polymorphisms”, SNPs), as well as insertions or deletions (indels). Identifying these variants is a dedicated analysis step known as “variant calling”, which evaluates the aligned reads in the BAM/SAM file.
The result of variant calling is typically stored in a “Variant Call Format” (VCF) file. A VCF summarizes the detected variations, their genomic positions, the observed alleles, various quality and confidence metrics as well as how many sequencing reads support each variant.
Reference-based variant analysis is widely used to identify specific mutations or to determine the genotype of an organism associated with a particular phenotype. By locating affected genes or genomic regions, researchers can investigate how genetic changes alter biological function and lead to observable traits, diseases, or experimental outcomes.
Typical workflow:
Read mapping (alignment) -> Variant calling
FASTQ (+ reference FASTA) -> BAM/SAM -> VCF
Common applications:
Microbial strain identification
Evolutionary genetics
Genetic variant discovery
Popular Tools:
BLAST
BWA, BWA-MEM
STAR
Minimap2
samtools

### Use-case 2: De-novo sequence analysis
This type of analysis also starts with DNA sequencing data in the form of sequencing reads, typically stored in FASTQ-files (link), which contain both nucleotide sequences and their per-base quality scores. In contrast to reference-based approaches, no suitable reference genome is assumed to be available, or the goal is to reconstruct the genome independently of an existing reference.
Instead of mapping reads to a reference, the sequencing reads are assembled de novo ( latin for “from the beginning”), meaning they are computationally combined based on sequence overlap and consistency to reconstruct longer contiguous sequences (contigs). These contigs may be further connected into larger structures called scaffolds, representing an approximation of the original genome sequence.
The result of the assembly process is typically stored in FASTA format, containing the reconstructed contig or scaffold sequences. Quality metrics associated with the assembly, such as contig length distributions, coverage statistics, or measures of completeness, are usually generated alongside the assembled sequences rather than embedded directly in the FASTA file.
Once a genome has been assembled, an additional analysis step known as genome annotation is commonly performed. During annotation, genomic features such as genes, coding sequences, and regulatory elements are predicted and assigned functional information based on sequence patterns, similarity to known genes, or external databases.
De novo analysis enables the identification of previously unknown genes, genomic structures, or large-scale sequence differences, which may not be detectable using reference-based methods. It is therefore especially important for studying non-model organisms, newly discovered species, highly diverged strains, or genomes that have undergone substantial structural rearrangements.
De novo assembly workflows are also frequently used in metagenomics, where sequencing data originates from a mixture of multiple organisms and no single reference genome can adequately represent the sample.
Typical workflow:
De novo assembly -> Genome annotation (-> optional comparative analysis)
FASTQ -> Contigs/Scaffolds (FASTA) -> Annotated genome
Common applications:
Sequencing new species
Pathogen discovery

### Use-case 3: Metagenomic analysis
This type of analysis starts with DNA sequencing data in the form of sequencing reads, typically stored in FASTQ-files (link), which contain both nucleotide sequences and their per-base quality scores. In metagenomics, these reads originate from a mixture of multiple organisms present in a single sample, such as environmental, clinical, or host-associated communities.
Because no single reference genome can adequately represent all organisms in such a sample, metagenomic analyses often rely on de novo approaches (see de novo analysis chapter) to reconstruct genomic fragments from the sequencing reads. Depending on the research question, reads may be assembled into contigs or analyzed directly without full assembly.
Following assembly, contigs can be grouped into bins, each representing a putative genome or taxonomic unit. This process, known as binning, uses sequence composition, coverage patterns, and similarity to known genomes. The resulting bins are commonly referred to as metagenome-assembled genomes (MAGs).
In addition to genome reconstruction, metagenomic workflows frequently include taxonomic profiling, which estimates which organisms are present in the sample, and functional profiling, which identifies metabolic pathways or gene functions encoded by the community. These analyses allow researchers to characterize both the composition and functional potential of complex microbial ecosystems.
Metagenomic analysis is widely used to study microbial communities in environments where cultivation of individual organisms is difficult or impossible. Typical applications include microbiome research, environmental monitoring, pathogen discovery, and the investigation of community-level responses to environmental or experimental changes.
Typical workflow:
Optional assembly -> Binning / Profiling -> Functional analysis
FASTQ -> Contigs (FASTA) -> MAGs / Taxonomic & functional profiles

### Use-case 4: Transcriptomic analysis
This type of analysis starts with RNA sequencing data in the form of sequencing reads, typically stored in FASTQ-files (link), which contain both nucleotide sequences and their per-base quality scores. In transcriptomics, these reads represent RNA molecules that were present in a cell or tissue at the time of sampling and are commonly generated through RNA-seq experiments.
Before analysis, RNA reads are typically mapped (aligned) to a reference genome or transcriptome, similar to reference-based DNA analysis. Alternatively, alignment-free or pseudo-alignment methods may be used to directly assign reads to known transcripts. The goal of this step is not to identify sequence variation, but to quantify gene or transcript expression levels.
The mapping or quantification process produces a file that links sequencing reads to genes or transcripts, along with quality and confidence measures. From this information, expression counts or abundance estimates are generated, describing how frequently each gene or transcript is observed in the data.
Once expression levels have been quantified, downstream analyses can be performed, most notably differential expression analysis, which identifies genes whose expression differs significantly between experimental conditions, developmental stages, or sample groups. Additional analyses often include clustering, gene set enrichment, and pathway analysis.
Transcriptomic analysis provides insight into the functional activity of a genome under specific conditions. It is widely used to study cellular responses to environmental changes, disease mechanisms, developmental processes, and the effects of genetic or pharmacological perturbations.

Typical workflow:
Read mapping / quantification -> Expression analysis -> Interpretation
FASTQ (+ reference genome/transcriptome) -> Expression matrix -> Differential expression results

Common applications:
Understanding cellular responses to stimuli (e.g. drugs, toxins)
Cell type identification
Improving cellular characteristics
























From TUM_Services_Knoweldge_Hub.xlsx












Ethics
Ethical Data Intiative — TUM Think Tank



| Name | Usage | Type | Description | Link (TUM internal) | Link (Documentation) |
| --- | --- | --- | --- | --- | --- |
|  |  |  |  |  |  |
|  |  |  |  |  |  |
|  |  |  |  |  |  |
|  |  |  |  |  |  |
|  |  |  |  |  |  |
|  |  |  |  |  |  |
|  |  |  |  |  |  |


| Database | Description | Link |
| --- | --- | --- |
| Gene Expression Omnibus (GEO) | GEO is a major public functional genomics data repository that archives high-throughput microarray and next-generation sequence functional genomic datasets. | https://www.ncbi.nlm.nih.gov/geo/ |
| NCBI Sequence Read Archive (SRA) | The SRA is NIH's archive of high-throughput sequencing data and is part of the International Nucleotide Sequence Database Collaboration (INSDC) that includes the NCBI Sequence Read Archive (SRA), the European Bioinformatics Institute (EBI), and the DNA Database of Japan (DDBJ). Data submitted to any of the three organizations are shared among them. | https://www.ncbi.nlm.nih.gov/sra/ |
| DNA Data Bank of Japan (DDBJ ) | Bioinformation and DDBJ Center (BI-DDBJ) collects nucleotide sequence data as a member of INSDC(International Nucleotide Sequence Database Collaboration). | https://www.ddbj.nig.ac.jp/ |
| EMBL European Nucleotide Archive (ENA) | The European Nucleotide Archive (ENA) provides a comprehensive record of the world’s nucleotide sequencing information, covering raw sequencing data, sequence assembly information and functional annotation. | https://www.ebi.ac.uk/ena/browser/ |
| ArrayExpress | The functional genomics data collection (ArrayExpress), stores data from high-throughput functional genomics experiments, and provides data for reuse to the research community. In line with community guidelines, a study typically contains metadata such as detailed sample annotations, protocols, processed data and raw data. Raw sequence reads from high-throughput sequencing studies are brokered to the European Nucleotide Archive (ENA), and links are provided to download the sequence reads from ENA. | https://www.ebi.ac.uk/biostudies/arrayexpress |
| UniProt | UniProtKB is on of the most widely used protein databases. It consists of the expertly curated component of UniProtKB and non-peer reviewed. It contains hundreds of thousands of protein descriptions, including function, domain structure, subcellular location, post-translational modifications and functionally characterized variants. | http://www.ebi.ac.uk/uniprot/ |
| The European Genome-phenome Archive (EGA) | The European Genome-phenome Archive (EGA) is a service for permanent archiving and sharing of personally identifiable genetic, phenotypic, and clinical data generated for the purposes of biomedical research projects or in the context of research-focused healthcare systems. | https://ega-archive.org/ |
| GTEx (Genotype-Tissue Expression) | Maps the relationship between genotype and tissue-specific gene expression. |  |
| Genome Sequence Archive (GSA) |  | https://ngdc.cncb.ac.cn/gsa/ |
| The Cancer Genome Atlas (TCGA) | The TCGA is a collaboration between the National Cancer Institute (NCI) and National Human Genome Research Institute (NHGRI). It has generated comprehensive, multi-dimensional maps of the key genomic changes in 33 types of cancer. The TCGA dataset, 2.5 petabytes of data describing tumor tissue and matched normal tissues from more than 11,000 patients. | http://cancergenome.nih.gov/ |
| PDB |  | https://www.rcsb.org/ |


| RefSeq | The Reference Sequence (RefSeq) collection provides a comprehensive, integrated, non-redundant, well-annotated set of sequences, including genomic DNA, transcripts, and proteins. RefSeq sequences form a foundation for medical, functional, and diversity studies. They provide a stable reference for genome annotation, gene identification and characterization, mutation and polymorphism analysis (especially RefSeqGene records), expression studies, and comparative analyses. | https://www.ncbi.nlm.nih.gov/refseq/ |
| --- | --- | --- |
| Ensembl | Maintained by EMBL European Bioinformatics Institute | https://www.ensembl.org/info/about/index.html |
| GenBank |  |  |
| UCSC Genome Browser | It is an interactive website offering access to genome sequence data from a variety of vertebrate and invertebrate species and major model organisms, integrated with a large collection of aligned annotations. The Browser is a graphical viewer optimized to support fast interactive performance and is an open-source, web-based tool suite built on top of a MySQL database for rapid visualization, examination, and querying of the data at many levels. |  |
|  |  |  |


| Category | Name | Link | Responsible | Contact | About | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| Statistical Consulting |  | https://www.math.cit.tum.de/en/math/department/statistical-consulting/ |  |  |  |  |
| RDM |  | https://web.tum.de/en/researchdata/hub/ |  |  |  |  |
| Ethics Committee |  | https://www.tum.de/en/about-tum/organization/ethics-committee |  |  |  |  |


| Category | Name | Link | Responsible | Contact | About | Notes | Add? |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RDM | Metadata Crawler | https://gitlab.lrz.de/nfdi4ing/crawler | NFDI4Ing |  | The crawler is capable to read out ontologies (.owl-files) and create a dictionary of relevant properties to be filled by the user. In a second step this dictionary is read-out to create a metadata-file which accompanies the main dataset. The dictionary therein may be filled with instruction where to find the relevant inputs rather than hard-coding them directly. |  | Yes |
| Data Analysis | PredictProtein | https://predictprotein.org/ |  |  | PredictProtein at www.predictprotein.org is a service that combines diverse traditional, homology-based protein sequence analysis tools. You can directly submit a protein sequence or an alignment there and receive predictions for various features of the submitted sequence. | Check with group |  |
| High-Performance Computing (HPC) |  | https://marge.aer.ed.tum.de/ | NFDI4Ing |  | MARGE (Multi-Access Research Gateway for HPC Experts): Supports researchers in managing HPC-related data challenges |  |  |