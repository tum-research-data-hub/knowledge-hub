---
 title: Repositories 
 slug: /domain-knowledge/bioinformatics/bio-repo
---

# **Bioinformatics Repositories**

Bioinformatics resources fall into two groups that are easy to confuse:

- **Deposition archives** accept *your* data, assign stable accession numbers, and preserve the data long-term. Journals and funders expect raw data to be deposited here.
- **Knowledge bases and reference resources** integrate, curate, and re-annotate data from many sources. You use them for analysis, but you do not normally submit your own experiment to them.

## **Where Do I Submit?**

| I have ... | Submit to | Notes |
|---|---|---|
| Raw sequencing reads (any organism, non-human or fully consented open data) | [ENA](https://www.ebi.ac.uk/ena/browser/) (or SRA / DDBJ) | The three INSDC archives mirror each other, so submit only once. ENA is the European entry point. |
| Genome or transcriptome assemblies, annotated sequences | ENA (or GenBank / DDBJ) | Register a project and samples first, then submit the assembly and its annotation. |
| Functional genomics experiments (RNA-seq, ChIP-seq, methylation, arrays) with processed data | [GEO](https://www.ncbi.nlm.nih.gov/geo/) or [ArrayExpress in BioStudies](https://www.ebi.ac.uk/biostudies/arrayexpress) | Both store processed data (e.g. count matrices) together with the experimental design; raw reads are brokered to SRA or ENA. |
| Human data that can identify individuals (genomes, exomes, clinical phenotypes) | [EGA](https://ega-archive.org/); in Germany the [GHGA](https://www.ghga.de/) national node | Controlled access: a Data Access Committee decides on requests. Never submit such data to open archives. |
| Genetic variants (VCF), non-human or aggregated | [European Variation Archive (EVA)](https://www.ebi.ac.uk/eva/) | |
| Mass spectrometry proteomics | [PRIDE](https://www.ebi.ac.uk/pride/) (ProteomeXchange) | |
| Metabolomics | [MetaboLights](https://www.ebi.ac.uk/metabolights/) | |
| Microscopy and other biological images | [BioImage Archive](https://www.ebi.ac.uk/bioimage-archive/) | |
| Macromolecular structures | [wwPDB](https://www.wwpdb.org/) (via RCSB PDB, PDBe or PDBj) | |
| Anything without a dedicated archive, or supplementary files linking several datasets | [BioStudies](https://www.ebi.ac.uk/biostudies/), [Zenodo](https://zenodo.org/) or [mediaTUM](https://mediatum.ub.tum.de/) | See the general [Repository Collections](/general-knowledge/repo) page. |

:::tip Projects and samples come first
Almost every submission starts with two records that tie everything together:

- A **project** record (**BioProject**, accession `PRJEB...`, `PRJNA...` or `PRJDB...`) describes the overall study. It is the single identifier you cite in the data availability statement of a paper, and it links reads, assemblies, and other data types.
- A **sample** record (**BioSample**, accession `SAMEA...`, `SAMN...` or `SAMD...`) describes each biological sample: organism, tissue, treatment, collection details.

Register the project and samples before uploading data files, and reuse the same sample accessions when you submit to different archives (e.g. reads to ENA and proteomics to PRIDE). This is what makes a multi-omics study traceable.
:::

Practical advice:

- **Submit early.** Archives let you keep data private until a release date you choose, so you can obtain accession numbers for your manuscript in advance.
- **Use the archive's metadata checklist.** The sample checklists implement community standards; see [Metadata and Ontologies](/domain-knowledge/bioinformatics/bio-metadata-ontologies#3-metadata-standards-and-checklists).
- **Keep the submission files.** Store the sample sheet and manifest files with your project documentation so that corrections remain possible later.

## **Deposition Archives**

| Database | Description | Link |
|---|---|---|
| EMBL European Nucleotide Archive (ENA) | The European Nucleotide Archive (ENA) provides a comprehensive record of the world’s nucleotide sequencing information, covering raw sequencing data, sequence assembly information and functional annotation. | https://www.ebi.ac.uk/ena/browser/ |
| NCBI Sequence Read Archive (SRA) | The SRA is NIH's archive of high-throughput sequencing data and is part of the International Nucleotide Sequence Database Collaboration (INSDC) that includes the NCBI Sequence Read Archive (SRA), the European Bioinformatics Institute (EBI), and the DNA Database of Japan (DDBJ). Data submitted to any of the three organizations are shared among them. SRA is a primary repository for functional genomics data such as RNA-Seq, ChIP-Seq, and Methyl-Seq. | https://www.ncbi.nlm.nih.gov/sra/ |
| DNA Data Bank of Japan (DDBJ) | DDBJ collects nucleotide sequence data along with their study and sample information as a member of INSDC (International Nucleotide Sequence Database Collaboration). | https://www.ddbj.nig.ac.jp/ |
| GenBank | GenBank is an open access database of all publicly available nucleotide sequences and their protein translations and annotations | https://www.ncbi.nlm.nih.gov/genbank/  |
| Genome Sequence Archive (GSA) |  GSA is a repository for archiving raw sequence reads in the National Genomics Data Center, China National Center for Bioinformation (CNCB-NGDC). | https://ngdc.cncb.ac.cn/gsa/ |
| BioProject | BioProject records describe a research project and group all data generated within it under one accession. Projects registered at ENA, NCBI or DDBJ are exchanged between the INSDC partners, so a BioProject accession resolves at all three. | https://www.ncbi.nlm.nih.gov/bioproject/ |
| BioSamples | BioSamples stores and supplies descriptions and metadata about biological samples used in research. Sample accessions are shared between the INSDC partners and are referenced by many other EMBL-EBI archives (ENA, PRIDE, MetaboLights, EGA), which makes them the link between different data types derived from the same sample. | https://www.ebi.ac.uk/biosamples/ |
| Gene Expression Omnibus (GEO) | GEO is a major public functional genomics data repository that archives high-throughput microarray and next-generation sequence functional genomic datasets. GEO stores the raw data, processed files, and metadata | https://www.ncbi.nlm.nih.gov/geo/ |
| Biostudies | Biostudies is a database dedicated to storing and organizing multimodal data from biological studies and relevant links to other databases. Biostudies also hosts functional genomics data from ArrayExpress, including RNA and protein expression and methylation profiling data. A central concept is a `study,` which contains metadata, including annotations, protocols, processed data and raw data. | https://www.ebi.ac.uk/biostudies/ |
| The European Genome-phenome Archive (EGA) | The European Genome-phenome Archive (EGA) is a service for permanent archiving and sharing of personally identifiable genetic, phenotypic, and clinical data generated for biomedical research projects or in research-focused healthcare systems. | https://ega-archive.org/ |
| PRoteomics IDEntifications Database (PRIDE) | PRIDE is the largest public repository for mass spectrometry-based proteomics data and a founding member of the ProteomeXchange consortium. It stores raw files, peak lists, identification and quantification results together with the experimental metadata, and assigns a ProteomeXchange accession (`PXD...`) to each dataset. | https://www.ebi.ac.uk/pride/ |
| MetaboLights | MetaboLights is an open repository for metabolomics experiments and derived information. It covers metabolite structures and their reference spectra as well as raw and processed data from mass spectrometry and NMR studies, described with the ISA metadata framework. Studies receive an accession of the form `MTBLS...`. | https://www.ebi.ac.uk/metabolights/ |
| BioImage Archive | The BioImage Archive stores and distributes biological images from all imaging modalities and scales, in particular the image data underlying publications. It also provides the storage for added-value image resources such as the Image Data Resource (IDR) and EMPIAR. | https://www.ebi.ac.uk/bioimage-archive/ |
| RCSB Protein Data Bank (PDB) | PDB is a global repository of the experimentally determined 3D structures of large biological molecules, including proteins, DNA and RNA. It is managed by the Worldwide Protein Data Bank organization (wwPDB; wwpdb.org). | https://www.rcsb.org/ |

## **Knowledge Bases and Reference Resources**

| Database | Description | Link |
|---|---|---|
| Reference Sequence (RefSeq) | The RefSeq database collects a comprehensive, integrated, non-redundant, set of sequences, including genomic DNA, transcripts, and proteins. RefSeq records are subjected to an automated processing pipeline combined with manual curation component to provide a robust set of reference sequences of genomic, transcript and protein data across a wide range of taxa. | https://www.ncbi.nlm.nih.gov/refseq/ |
| Ensembl |  Ensembl is a widely used open resource that integrates and analyses publicly available genomics data from both eukaryotes and prokaryotes. It provides genome assemblies and rich annotations—such as genes, transcripts, proteins, regulatory elements, genetic variants, and comparative genomics (e.g., homology and gene trees)—that can be explored through an interactive genome browser and downloaded in standard formats.  Ensembl is maintained by EMBL European Bioinformatics Institute. | https://www.ensembl.org/info/about/index.html |
| UCSC Genome Browser | The UCSC Genome Browser is an interactive website providing access to genome sequence data from a wide range of vertebrate and invertebrate species and major model organisms, integrated with a large collection of aligned annotations. |  https://genome.ucsc.edu/ |
| UniProt | UniProtKB is one of the most widely used protein databases. It consists of an expertly curated component and a non-peer-reviewed component. It contains hundreds of thousands of protein descriptions, including function, domain structure, subcellular location, post-translational modifications, and functionally characterized variants. | https://www.uniprot.org/ |
| The Encyclopedia of DNA Elements (ENCODE) | ENCODE is a public resource containing comprehensive data stemming from efforts to map all functional DNA elements present in the human genome. ENCODE integrates and annotates data on genes, transcripts, and transcriptional regulatory regions, as well as their chromatin states and DNA methylation patterns. | https://www.encodeproject.org/ |
| Genotype-Tissue Expression (GTEx) | GTEx is a comprehensive resource and tissue biobank that aims to map the genetic effects on the transcriptome across diverse human tissues and to link these regulatory mechanisms to trait and disease. GTEx provides molecular assays, including whole-genome sequencing, whole-exome sequencing, and RNA-seq.| https://gtexportal.org/home/ |
| NCI Genomic Data Commons (GDC) | The GDC is the US National Cancer Institute's data portal for harmonised cancer genomics data. It provides uniformly re-processed genomic, transcriptomic, epigenomic and clinical data from large programmes, including The Cancer Genome Atlas (TCGA), which characterised tumour and matched normal tissue across 33 cancer types. Summary-level data are open; individual-level sequence data require an approved access request. | https://portal.gdc.cancer.gov/ |

**Looking for a resource that is not listed?** Search [FAIRsharing](https://fairsharing.org/) or [re3data](https://www.re3data.org/), or browse the [ELIXIR Deposition Databases](https://elixir-europe.org/platforms/data/elixir-deposition-databases) and the RDMkit pages for your subfield, e.g. [plant sciences](https://rdmkit.elixir-europe.org/plant_sciences), [proteomics](https://rdmkit.elixir-europe.org/proteomics), [single-cell sequencing](https://rdmkit.elixir-europe.org/single_cell_sequencing) or [bioimaging](https://rdmkit.elixir-europe.org/bioimaging_data).
