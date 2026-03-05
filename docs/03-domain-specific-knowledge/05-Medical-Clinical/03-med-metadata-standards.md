---
title: Medical & Clinical Metadata & Standards
slug: /domain-knowledge/medical-clinical/metadata-standards
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Medical & Clinical Metadata & Standards**

Medical and clinical data is governed by strict regulatory requirements that dictate how data must be described, stored, shared, and protected. Metadata standards in this domain serve both scientific reproducibility and legal compliance.

---

## **Core Interoperability Standards**

### **HL7 FHIR (Fast Healthcare Interoperability Resources)**
The modern standard for health data exchange. FHIR defines a set of "resources" (Patient, Observation, Condition, Medication, etc.) with standardized fields, enabling interoperable data exchange via REST APIs.
- **Use:** EHR integration, clinical data exchange, research data access
- **Spec:** [HL7 FHIR](https://www.hl7.org/fhir/)

### **OMOP CDM (Common Data Model)**
A standardized relational model for harmonizing heterogeneous clinical data from different sources into a common format, maintained by OHDSI.
- **Use:** Multi-site observational studies, pharmacovigilance
- **Spec:** [OMOP CDM](https://ohdsi.github.io/CommonDataModel/)

### **CDISC Standards**
A suite of standards for clinical trial data:
- **CDASH** — Standardized data collection forms
- **SDTM** — Tabulation format for regulatory submission
- **ADaM** — Analysis-ready dataset format
- **Spec:** [CDISC](https://www.cdisc.org/)

---

## **Medical Terminology and Ontologies**

### **SNOMED CT**
The most comprehensive clinical terminology, covering diseases, findings, procedures, and substances. Used in electronic health records worldwide.
- **Spec:** [SNOMED International](https://www.snomed.org/)

### **ICD-10 / ICD-11**
International Classification of Diseases — the WHO standard for coding diagnoses, used for billing and epidemiology.
- **Spec:** [WHO ICD](https://icd.who.int/)

### **LOINC (Logical Observation Identifiers Names and Codes)**
A universal coding system for laboratory tests and clinical observations.
- **Use:** Standardizing lab result reporting
- **Spec:** [LOINC](https://loinc.org/)

### **MeSH (Medical Subject Headings)**
A controlled vocabulary for indexing biomedical literature, maintained by the NLM.
- **Use:** Literature search, metadata tagging
- **Spec:** [MeSH](https://www.nlm.nih.gov/mesh/)

### **HPO (Human Phenotype Ontology)**
A standardized vocabulary of phenotypic abnormalities, used in rare disease research and clinical genomics.
- **Spec:** [HPO](https://hpo.jax.org/)

---

## **Privacy and Consent Standards**

### **GDPR (General Data Protection Regulation)**
The EU regulation governing personal data processing. Requires lawful basis, purpose limitation, data minimization, and rights of data subjects.
- **Relevance:** All research involving EU patient data

### **Informed Consent Ontology (ICO)**
An ontology for representing informed consent information in research.
- **Spec:** [ICO](https://github.com/ICO-ontology/ICO)

### **GA4GH Data Use Ontology (DUO)**
A standard for machine-readable data use conditions, enabling automated access control decisions.
- **Use:** Genomic and clinical data sharing with consent restrictions
- **Spec:** [GA4GH DUO](https://github.com/EBISPOT/DUO)

---

## **DICOM Metadata**

DICOM files contain extensive metadata organized in a hierarchical structure:
- **Patient level** — Name, ID, birth date, sex
- **Study level** — Date, description, referring physician
- **Series level** — Modality, body part, protocol
- **Image level** — Pixel data, acquisition parameters

**Important:** DICOM metadata must be anonymized (de-identified) before research use. Tools like `deid` (Python) or CTP (Clinical Trial Processor) automate this process.

---

## **Best Practices**

1. **Anonymize early** — De-identify patient data as soon as possible in the research pipeline
2. **Use standard terminologies** — SNOMED CT, ICD, LOINC for consistent coding
3. **Document consent scope** — Record what the patient consented to and encode it machine-readably (DUO)
4. **Map to common data models** — Use OMOP CDM for observational research to enable cross-site analysis
5. **Version your data dictionaries** — Track changes to variable definitions over time
6. **Maintain an audit trail** — Log all data access and modifications for regulatory compliance
