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

## **3\. Finding Other Ontologies🔍**

If your specific niche is not covered by the examples above, you can explore comprehensive registries:

* [**BioPortal**](https://bioportal.bioontology.org/)**:** A web portal for accessing and sharing biomedical ontologies.  
* [**OLS4 (Ontology Lookup Service)**](https://www.ebi.ac.uk/ols4/)**:** A repository provided by EMBL-EBI that allows you to search for terms across hundreds of ontologies simultaneously.



:::info Tip for TUM Researchers  
When using the TUM DataTagger or uploading data to mediaTUM, check if your field has a mandated ontology. Using these standard terms early in your project will significantly increase the citability of your final data.  
:::