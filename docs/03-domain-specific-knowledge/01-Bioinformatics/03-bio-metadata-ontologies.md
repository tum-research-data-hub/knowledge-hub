---
title: Metadata and Ontologies
slug: /domain-knowledge/bioinformatics/bio-metadata-ontologies
---

# **Metadata and Ontologies in Bioinformatics**

In data-intensive research, raw data is only as valuable as the context provided with it. Standardised metadata and ontologies enable researchers to search, compare, and integrate datasets across different laboratories and global repositories.

## **1\. Why Ontologies Matter💡**

To answer specific biological questions, researchers must find datasets that match precise criteria. However, natural language is often ambiguous (e.g., "cancer" vs "neoplasm"). **Ontologies** solve this by providing a controlled set of terms and structured relationships between them.

An ontology represents the shared background knowledge of a community. One of the most basic forms is a **hierarchical classification**, which uses an "is a" relationship (e.g., a *Hepatocyte* **is a** *Eukaryotic Cell*). By using ontologies, metadata becomes machine-readable, allowing for automated data integration and sophisticated computational analysis.

## **2\. Core Biological Ontologies**

### **Gene Ontology (GO)**

The [Gene Ontology](https://geneontology.org/) is the most widely used resource for describing gene function across the tree of life. It organises biological knowledge into three main aspects:

* **Molecular Function (MF):** Activities performed by a gene product (e.g., *catalytic activity*).  
* **Cellular Component (CC):** The locations relative to cellular structures where these functions occur (e.g., *mitochondrion*).  
* **Biological Process (BP):** A "biological programme" comprising multiple molecular activities acting in concert (e.g., *DNA repair*).

### **Human Phenotype Ontology (HPO)**

The [HPO](https://hpo.jax.org/) provides a standardised vocabulary for describing human abnormalities and clinical symptoms (phenotypes). It is a key tool for deep phenotyping and phenotype-driven diagnosis, particularly in rare disease research and global health initiatives like the GA4GH.

### **Human Disease Ontology (DO)**

The [Human Disease Ontology](https://disease-ontology.org/) facilitates disease classification by integrating disease-specific terms and identifiers from biomedical vocabularies such as ICD, SNOMED CT, and MeSH. It enables consistent data analysis across genomic databases.

### **Plant Ontology (PO)**

Developed by the [Plant Ontology Consortium](http://www.plantontology.org/), the PO describes plant anatomy, morphology, and growth stages. It uses the same data model as the Gene Ontology to annotate gene expression and phenotype data in plant sciences.

### **Ontology of Host-Microbiome Interactions (OHMI)**

The [OHMI](https://github.com/ohmi-ontology/ohmi) standardises entities related to microbiomes, host organisms (human, mouse, etc.), and their interactions. It comprises over 1,000 terms, including specific terms for study protocols and experimental results.

### **Chemical Entities of Biological Interest (ChEBI)**

[ChEBI](https://www.ebi.ac.uk/chebi/) is a database and ontology focusing on small chemical compounds. It adheres to nomenclature and terminology endorsed by IUPAC and NC-IUBMB.

## **3\. Metadata Standards and Checklists**

Ontologies provide the *terms*; metadata standards define *which information must be reported* about an experiment and how it is structured. General-purpose schemas such as DataCite are covered on the [Metadata](/general-knowledge/metadata) page and by the [TUM Research Data Hub](https://web.tum.de/en/researchdata/support-information/store-manage-data/). The standards below are specific to the life sciences.

### **Minimum Information Standards**

"Minimum information" checklists list the details that must accompany a dataset so that others can interpret and reuse it. Public archives implement them as submission checklists, so choosing the right checklist at the start of a project saves time at submission.
Examples of Metadata Standards used in popular archives include:

| Standard | Applies to | Where you meet it |
|---|---|---|
| [**MIxS**](https://genomicsstandardsconsortium.github.io/mixs/) (Minimum Information about any (x) Sequence), incl. MIGS, MIMS, MIMARKS, MIMAG, MISAG | Genomes, metagenomes, marker-gene surveys, MAGs and single-cell genomes, with environment-specific packages (soil, water, host-associated, ...) | ENA and NCBI BioSample checklists |
| [**MINSEQE**](https://www.fged.org/projects/minseqe/) / **MIAME** | High-throughput sequencing and microarray functional genomics experiments | GEO and ArrayExpress submissions |
| [**MIAPPE**](https://www.miappe.org/) | Plant phenotyping experiments, incl. the description of the biological material | BioSamples plant checklist, FAIDARE, e!DAL |
| [**MIAPE**](https://www.psidev.info/miape) and the SDRF-Proteomics sample table | Mass spectrometry proteomics | PRIDE / ProteomeXchange |
| [**MSI** reporting standards](https://metabolomicssociety.org/) | Metabolomics | MetaboLights |
| [**REMBI**](https://www.ebi.ac.uk/bioimage-archive/rembi-help-overview/) | Biological imaging (light and electron microscopy) | BioImage Archive |
| **minSCe** | Single-cell sequencing experiments | ArrayExpress, Human Cell Atlas |

### **Structuring a Whole Study: ISA, ARC and RO-Crate**

* [**ISA (Investigation - Study - Assay)**](https://isa-tools.org/) is a general framework for describing multi-omics experiments: the *investigation* is the overall project, a *study* describes the biological material and its treatment, and an *assay* describes a measurement performed on it. MetaboLights and BioStudies use ISA-based metadata.
* The [**Annotated Research Context (ARC)**](https://arc-rdm.org/), developed by the NFDI consortium DataPLANT, packages data, ISA metadata, workflows and results of a project into one version-controlled folder structure.
* [**RO-Crate**](https://www.researchobject.org/ro-crate/) is a lightweight, domain-independent way of packaging research data with machine-readable metadata; it is used by WorkflowHub and can wrap an ARC.

### **Sample Metadata in Practice**

1. Decide early which archive your data will go to (see [Where Do I Submit?](/domain-knowledge/bioinformatics/bio-repo#where-do-i-submit)) and download its sample checklist.
2. Record the checklist fields while the samples are collected, for example in [TUM eLabFTW](https://www.ub.tum.de/en/eln), instead of reconstructing them at submission time.
3. Fill fields with ontology terms (Section 2) wherever the checklist allows it, and always use the NCBI Taxonomy identifier for the organism.

**Further information:** The ELIXIR [RDMkit](https://rdmkit.elixir-europe.org/metadata_management) explains documentation and metadata for the life sciences in general, with dedicated pages for e.g. [plant sciences](https://rdmkit.elixir-europe.org/plant_sciences), [proteomics](https://rdmkit.elixir-europe.org/proteomics), [single-cell sequencing](https://rdmkit.elixir-europe.org/single_cell_sequencing) and [bioimaging](https://rdmkit.elixir-europe.org/bioimaging_data). [FAIRsharing](https://fairsharing.org/) lists all registered standards by discipline.

## **4\. Finding Other Ontologies🔍**

If your specific niche is not covered by the examples above, you can explore comprehensive registries:

* [**BioPortal**](https://bioportal.bioontology.org/)**:** A web portal for accessing and sharing biomedical ontologies.  
* [**OLS4 (Ontology Lookup Service)**](https://www.ebi.ac.uk/ols4/)**:** A repository provided by EMBL-EBI that allows you to search for terms across hundreds of ontologies simultaneously.



:::info Tip for TUM Researchers  
When using the TUM DataTagger or uploading data to mediaTUM, check if your field has a mandated ontology. Using these standard terms early in your project will significantly increase the citability of your final data.  
:::