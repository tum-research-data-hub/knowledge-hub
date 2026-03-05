---
title: Chemistry & Materials Metadata & Standards
slug: /domain-knowledge/chemistry-materials/metadata-standards
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Chemistry & Materials Metadata & Standards**

Reproducibility in chemistry and materials science requires thorough documentation of experimental conditions, instrument settings, and computational parameters. Several standards and best practices help ensure that data is FAIR (Findable, Accessible, Interoperable, Reusable).

---

## **Chemical Identifiers and Registries**

### **CAS Registry Numbers**
Unique numerical identifiers assigned by the Chemical Abstracts Service to chemical substances. The most widely recognized identifier in chemistry.
- **Example:** 7732-18-5 (water)
- **Note:** CAS numbers are proprietary; use InChI for open science

### **InChI / InChIKey**
IUPAC's open standard for unique molecular identification. InChIKey is a fixed-length hash for database lookups.
- **Use:** Cross-database linking, deduplication
- **Spec:** [InChI Trust](https://www.inchi-trust.org/)

### **PubChem CID / ChEBI ID**
Database-specific identifiers linking to curated compound entries with properties, synonyms, and references.

---

## **Metadata Standards**

### **ISA Framework (Investigation-Study-Assay)**
A general-purpose metadata framework for describing multi-assay experiments, widely used in life sciences and increasingly in chemistry.
- **Use:** Multi-step experimental workflows
- **Spec:** [ISA Tools](https://isa-tools.org/)

### **NeXus / HDF5 Metadata**
A standardized data format and metadata schema for neutron, X-ray, and muon experiments built on HDF5.
- **Use:** Synchrotron and neutron source experiments
- **Spec:** [NeXus Format](https://www.nexusformat.org/)

### **NOMAD Metadata Schema**
A comprehensive metadata schema for computational materials science, covering DFT, MD, and other simulation methods.
- **Use:** Archiving and sharing simulation data
- **Spec:** [NOMAD](https://nomad-lab.eu/)

### **Allotrope Data Format (ADF)**
A data standard for analytical chemistry, based on HDF5, developed by the Allotrope Foundation.
- **Use:** Pharmaceutical and analytical chemistry data
- **Spec:** [Allotrope Foundation](https://www.allotrope.org/)

---

## **Ontologies and Vocabularies**

### **ChEBI (Chemical Entities of Biological Interest)**
A freely available dictionary of molecular entities focused on chemical compounds.
- **Spec:** [ChEBI](https://www.ebi.ac.uk/chebi/)

### **ChEMBL**
A database of bioactive drug-like small molecules with activity data.
- **Spec:** [ChEMBL](https://www.ebi.ac.uk/chembl/)

### **Materials Design Ontology (MDO)**
An ontology for materials science covering materials, properties, computational methods, and provenance.
- **Spec:** [MDO](https://w3id.org/mdo/)

---

## **Best Practices**

1. **Use standard chemical identifiers** — InChI for open data, CAS numbers for literature reference
2. **Document instrument and method parameters** — Detailed enough for reproduction
3. **Record sample preparation** — Synthesis conditions, purification steps, sample handling
4. **Include software versions** — For computational work, record code version, basis sets, functionals, convergence criteria
5. **Use electronic lab notebooks** — eLabFTW or Chemotion for structured data capture at the point of creation
