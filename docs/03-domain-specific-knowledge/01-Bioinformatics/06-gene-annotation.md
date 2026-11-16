---
title: Gene & Genome Annotation
slug: /domain-knowledge/bioinformatics/gene-annotation
---

# **Gene & Genome Annotation**

An assembled genome is a long string of letters. **Annotation** turns it into something biologists can work with: where the genes are and what they do. Almost every downstream analysis, from RNA-seq quantification to comparative genomics, inherits the quality of the annotation it is based on, and annotations of the same organism can differ considerably between pipelines, versions and databases.

This guide helps you choose an approach, check the result, and document it so that others can reproduce and compare it. It follows on from the [De Novo Sequence Analysis](/domain-knowledge/bioinformatics/use-cases/de-novo-analysis) use case.

:::tip Reuse first
Check whether a maintained annotation of your genome already exists (Ensembl, RefSeq, a community database for your organism) before producing a new one, and state exactly which release you use. If you do annotate, prefer an established pipeline over a custom chain of tools.
:::

## **1\. Structural vs. Functional Annotation**

Annotation consists of two layers that are produced by different tools and should be assessed separately.

| | Structural annotation | Functional annotation |
|---|---|---|
| **Question** | *Where* are the features? | *What* do the gene products do? |
| **Input** | Genome assembly (repeat-masked), plus evidence such as RNA-seq reads or proteins of related species | The predicted protein (or transcript) sequences |
| **Output** | Coordinates of genes, exons, coding sequences, UTRs, non-coding RNAs, repeats (GFF3) | Product names, protein domains, GO terms, EC numbers, KEGG orthologs, pathway assignments |
| **Main approaches** | *Ab initio* prediction from sequence signals; evidence-based prediction from RNA-seq and protein alignments; projection (liftover) from an annotated relative; deep learning models | Transfer by sequence similarity or orthology; domain and family signatures (HMMs); structure-based similarity; machine learning |

**Structure caps function.** Functional annotators work on the gene models they are given. A fragmented, fused, or transposon-derived gene model receives a wrong or empty function regardless of the tool. If functional annotations look poor, check the gene models first.

A typical eukaryotic annotation therefore runs in this order:

```
Assembly → repeat masking → structural annotation (+ ncRNA) → quality checks
         → functional annotation → quality checks → standard formats → submission
```

## **2\. Tools by Kingdom**

Gene structure differs fundamentally between the domains of life (introns, gene density, genetic code), so there is no single tool for everything. Using a pipeline built for the wrong kingdom produces confident nonsense. The tables list widely used, openly documented options; they are a starting point, not a ranking.

### **Plants**

Plant genomes are large, repeat-rich and often polyploid. Repeat annotation and masking (e.g. EDTA, or RepeatModeler2 with RepeatMasker) before gene prediction is essential.

| Tool | Layer | Notes |
|---|---|---|
| [**BRAKER3**](https://github.com/Gaius-Augustus/BRAKER) | Structural | Combines RNA-seq and protein evidence with GeneMark and AUGUSTUS; a common reference pipeline. GeneMark has licence restrictions. |
| **GALBA** | Structural | Protein evidence only; useful when no RNA-seq is available but a well-annotated relative exists |
| [**Helixer**](https://github.com/weberlab-hhu/Helixer) | Structural | Deep learning, no evidence data needed, runs on a GPU; tends to over-predict and is best combined with other evidence |
| **Liftoff / miniprot** | Structural | Projection of an existing annotation onto a new assembly of the same or a close species |
| [**Mercator4 / MapMan**](https://www.plabipd.de/mercator_main.html) | Functional | Assigns proteins to the curated, plant-specific MapMan bin hierarchy |
| [**eggNOG-mapper**](https://github.com/eggnogdb/eggnog-mapper) | Functional | Orthology-based transfer of GO terms, KEGG orthologs, EC numbers and more; broad but generous with GO terms |
| [**InterProScan**](https://www.ebi.ac.uk/interpro/about/interproscan/) | Functional | Protein domains and families from the InterPro member databases; conservative and the most widely used baseline |
| **EnTAP** | Functional | Wrapper combining similarity search, eggNOG and InterProScan with contaminant filtering; designed for non-model organisms |
| **AHRD** | Functional | Generates human-readable gene descriptions with a per-gene quality code |

Organelle genomes need their own tools (e.g. GeSeq for plastids and mitochondria); nuclear pipelines mis-annotate organellar contigs.

### **Animals and Human**

| Tool | Notes |
|---|---|
| **GENCODE / Ensembl** | Reference annotations for human and mouse (manual curation combined with automated pipelines). Use them; do not re-annotate these genomes. |
| **NCBI EGAPx** | Public version of the NCBI eukaryotic annotation pipeline that produces RefSeq annotations |
| **BRAKER3, GALBA, Helixer** | As for plants |
| **TOGA, Liftoff** | Projection from a well-annotated reference species; performs very well when a close relative exists |
| **StringTie, IsoQuant** | Transcript assembly from short- or long-read RNA-seq, as evidence or as an annotation in its own right |

Functional annotation uses the same general-purpose tools as in plants (InterProScan, eggNOG-mapper), complemented by orthology to well-studied model organisms.

### **Bacteria and Archaea**

Prokaryotic genes are dense and lack introns, so annotation is fast and largely automated, with structure and function delivered in one run.

| Tool | Notes |
|---|---|
| [**Bakta**](https://github.com/oschwengers/bakta) | Current community standard for isolates and MAGs; versioned database, stable cross-references to RefSeq and UniProt, standards-compliant output |
| [**PGAP**](https://github.com/ncbi/pgap) | The NCBI pipeline used for RefSeq; can be run locally; a good choice when the genome will be submitted to GenBank |
| **Prokka** | Long the default and still common in the literature, but no longer actively developed; prefer Bakta for new projects |
| **Pyrodigal / Prodigal** | Gene calling only (used inside the pipelines above) |
| Specialist tools | AMRFinderPlus (resistance genes), antiSMASH (biosynthetic gene clusters), dbCAN (carbohydrate-active enzymes), DefenseFinder (defence systems) |

Remember to set the correct genetic code (translation table 11 for most bacteria and archaea, 4 for e.g. *Mycoplasma*).

### **Viruses and Phages**

| Tool | Notes |
|---|---|
| [**Pharokka**](https://github.com/gbouras13/pharokka) | Phage annotation based on the PHROG protein families, plus tRNA, CRISPR, virulence and resistance genes |
| **Phold** | Structure-based phage annotation (protein language model + Foldseek); assigns functions to many proteins that sequence similarity misses |
| **geNomad**, **CheckV** | Identification of viral sequences in assemblies; completeness and contamination estimates |
| **VADR** | Reference-model-based annotation of eukaryotic viruses, used for GenBank submission (e.g. SARS-CoV-2, influenza) |

A large fraction of viral proteins has no known function; a high share of "hypothetical protein" is expected and not by itself a sign of poor annotation.

### **Fungi**

| Tool | Notes |
|---|---|
| [**funannotate**](https://github.com/nextgenusfs/funannotate) | Widely used pipeline tailored to fungi, covering gene prediction, functional annotation and submission-ready output |
| **BRAKER3** | Evidence-based prediction from RNA-seq and protein data (see plants and animals) |
| Specialist tools | antiSMASH (fungal mode), dbCAN, EffectorP and SignalP for secreted proteins and effectors |


### **Ready-Made Workflows**

Several community workflows bundle these steps, for example [nf-core/funcscan](https://nf-co.re/funcscan) (prokaryotic contigs: annotation plus resistance and biosynthetic gene screening), [nf-core/proteinannotator](https://nf-co.re/proteinannotator) (functional annotation of protein sets) and [nf-core/genomeqc](https://nf-co.re/genomeqc) (comparing genomes and annotations). The annotation tutorials of the [Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/genome-annotation/) let you try complete workflows without installation.

## **3\. Quality Checks**

No single number describes annotation quality. Report several complementary metrics, and always state the tool and database versions with them.

| Check | What it tells you | Tools |
|---|---|---|
| **Format validity** | Whether the file is well-formed and internally consistent (parent-child relations, coordinates, phases) | AGAT, GenomeTools `gff3validator`, the validators of ENA and NCBI |
| **Basic statistics** | Number of genes, mean gene and exon length, exons per gene, share of single-exon genes; compare with related species | AGAT (`agat_sp_statistics.pl`) |
| **Completeness** | Share of conserved single-copy genes found complete, duplicated, fragmented or missing | BUSCO (run in protein mode on the predicted proteome), compleasm |
| **Consistency and contamination** | Whether the proteome fits the expected lineage; fragments, contaminants and implausible genes | OMArk |
| **Coding potential** | Likelihood that each predicted coding sequence is a real protein-coding gene | PSAURON |
| **Evidence support** | Agreement of each gene model with RNA-seq or protein alignments | Annotation Edit Distance (AED), intron support statistics |
| **Agreement with a reference** | Sensitivity and precision at gene, transcript and exon level, if a trusted annotation exists | gffcompare |
| **Transposon-derived genes** | Gene models that are really transposable elements (frequent in plants) | TEsorter, comparison with the repeat annotation |
| **Genome-level checks** | Completeness and contamination of the underlying assembly | CheckM2 (prokaryotes), CheckV (viruses), BUSCO in genome mode |

How to read the results:

* **BUSCO covers only a small, conserved fraction of the gene set.** A high BUSCO score shows that well-known genes were found; it says little about the remaining, lineage-specific genes. Choose the most specific lineage dataset available and report its name and version.
* **Compare genome-mode and protein-mode completeness.** If the assembly contains a BUSCO gene but the annotation does not, the gene prediction missed it.
* **Functional annotation has no accepted quality score.** Useful proxies are the share of proteins with any annotation, the share named "hypothetical protein", the agreement between independent tools, and the specificity of the assigned terms.
* **More is not better.** Many GO terms per protein often indicate over-annotation by permissive transfer, not better knowledge. Unsupported specific terms are harder to detect, and more damaging, than missing ones.

## **4\. Databases and Version Pinning**

Functional annotation is a lookup against reference databases. The result therefore depends on the database release as much as on the tool.

| Database | Role | Things to know |
|---|---|---|
| **InterPro** (Pfam, PANTHER, CDD, NCBIfam and others) | Protein domains and families | Regular releases; some member analyses (SignalP, TMHMM) need separate licences |
| **UniProtKB** (Swiss-Prot, TrEMBL) | Product names by best hit | Swiss-Prot is curated but small for non-model organisms; TrEMBL is large but propagates automatic annotations |
| **Gene Ontology (GO)** | Standard vocabulary for function, process and location | Released monthly; terms are merged and made obsolete over time |
| **eggNOG** | Orthologous groups for annotation transfer | Check which database version your eggNOG-mapper installation uses |
| **KEGG Orthology** | Pathway mapping | KOfam profiles are free; bulk access to KEGG itself is licensed |
| **Specialist databases** | CAZy (carbohydrate-active enzymes), MEROPS (peptidases), CARD (antibiotic resistance), PHROG (phage proteins), MapMan bins (plants) | Use them where relevant to your organism group |

**Pin and record everything.** Two runs of the same tool with database releases a year apart can give visibly different annotations. For every annotation, keep a record of:

* the assembly: accession or file checksum, and version
* each tool: name, exact version, and the parameters that deviate from the defaults
* each database: name, release number or download date, and checksum where available
* the container image or environment file used (see [Containers](/general-knowledge/containers))
* the evidence data (accessions of RNA-seq runs and protein sets)

Workflow managers write most of this automatically ([Workflow Management Systems](/general-knowledge/wms)). Download databases once to a shared, read-only location (they often take tens of gigabytes, and compute nodes usually have no internet access) and never update them in the middle of a project. When comparing several genomes, annotate all of them with the same tool and database versions. Mind the licences of tools and databases if you want to share your workflow.

## **5\. Output Formats**

| Content | Format | Notes |
|---|---|---|
| Gene structures | **GFF3** | The standard exchange format; validate before sharing. GTF is still required by some RNA-seq tools and can be converted with AGAT or gffread. |
| Sequences | **FASTA** | Proteins, transcripts and coding sequences extracted from the genome with the GFF3 (gffread, AGAT), so that they are guaranteed to match |
| GO annotations | **GAF 2.2** | The Gene Ontology Consortium's format, including evidence codes and the source of each assignment |
| Other functional annotation | Tab-separated table | One row per gene and assignment, with tool, database entry and score; keep the tools' native output as well |
| Submission to archives | **EMBL flat file** (ENA) or **GenBank / feature table** (NCBI) | Created from GFF3 with converters such as EMBLmyGFF3 or NCBI table2asn; the archives' validators enforce naming rules |
| Provenance | Workflow reports, a JSON or YAML manifest, or an RO-Crate | Versions and parameters as listed in Section 4 |

Good practice for identifiers and files:

* Give genes **stable, systematic identifiers** (locus tags with a registered prefix) and never reuse an identifier for a different gene. When an annotation is updated, publish a mapping from old to new identifiers.
* **Version the annotation** separately from the assembly (e.g. assembly v2, annotation v2.1) and always cite both.
* Decide on a rule for the **representative transcript** per gene (e.g. longest coding sequence or a tagged canonical isoform) and apply it consistently; many statistics and downstream tools depend on it.
* Deposit the annotation together with the assembly; see [Where Do I Submit?](/domain-knowledge/bioinformatics/bio-repo#where-do-i-submit). The file formats themselves are described under [Data Types](/domain-knowledge/bioinformatics/bio-data-types).

## **6\. Common Pitfalls**

:::warning Things that go wrong

**Mixing annotations from different sources**
- Comparing genomes annotated with different pipelines or versions turns every methodological difference into an apparent biological one, for example in gene presence/absence analyses and pangenomes. Re-annotate all genomes consistently or restrict comparisons to genes that pass the same criteria.

**Mismatched assembly and annotation versions**
- Coordinates are only valid for the exact assembly they were made for. Check sequence names and lengths before combining files from different sources (e.g. Ensembl and UCSC name chromosomes differently).

**Skipping repeat masking**
- Unmasked transposable elements are predicted as genes and inflate gene counts, especially in plants.

**Wrong kingdom or genetic code**
- Contaminant contigs (bacterial sequences in a eukaryotic assembly, organelles in a nuclear assembly) are annotated with the wrong model. Screen and separate them first.

**Isoform handling**
- Completeness and duplication scores change depending on whether all transcripts or one per gene are evaluated. State which you used.

**Trusting names by best hit**
- A product name copied from the best database hit may itself have been transferred automatically several times. Prefer domain- and orthology-based evidence, and keep the evidence for every assignment.

**Lost provenance**
- Without the tool and database versions, an annotation cannot be reproduced or meaningfully compared. Record them when the annotation is made; they are almost impossible to reconstruct later.

**Polyploid and highly heterozygous genomes**
- Collapsed or duplicated gene copies distort duplication statistics and gene family analyses. Interpret BUSCO duplication with the ploidy in mind.

**Circular benchmarks**
- A reference annotation used as "truth" was often produced by similar tools and databases. Agreement with it is not proof of correctness.

:::

## **Further Information**

* [Earth BioGenome Project: annotation recommendations](https://www.earthbiogenome.org/report-on-annotation-recommended-tools) list tools for each annotation step
* [Galaxy Training Network: Genome Annotation](https://training.galaxyproject.org/training-material/topics/genome-annotation/) offers hands-on tutorials for prokaryotes and eukaryotes
* [GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) and [GAF 2.2 specification](https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/)
* Related pages: [De Novo Sequence Analysis](/domain-knowledge/bioinformatics/use-cases/de-novo-analysis), [Metadata and Ontologies](/domain-knowledge/bioinformatics/bio-metadata-ontologies), [Repositories](/domain-knowledge/bioinformatics/bio-repo), [Glossary](/domain-knowledge/bioinformatics/glossary)
