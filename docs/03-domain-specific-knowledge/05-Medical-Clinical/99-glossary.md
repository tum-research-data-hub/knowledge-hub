---
title: Glossary
slug: /domain-knowledge/medical-clinical/glossary
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Medical & Clinical Data Glossary**

This glossary provides definitions for common terms used in medical data science, clinical research, and health informatics.

---

## **Clinical Research**

### **Clinical Trial**
A prospective research study that assigns participants to interventions (treatments, procedures, devices) to evaluate their effects on health outcomes.

### **Randomized Controlled Trial (RCT)**
A clinical trial in which participants are randomly assigned to treatment or control groups, considered the gold standard for evaluating interventions.

### **Cohort Study**
An observational study that follows a group of individuals over time to identify factors associated with outcomes.

### **Case-Control Study**
An observational study that compares individuals with a condition (cases) to similar individuals without it (controls) to identify risk factors.

### **Blinding (Single / Double)**
A method to prevent bias by concealing group assignments. Single-blind: participants don't know their group. Double-blind: neither participants nor investigators know.

### **Informed Consent**
The process by which a participant voluntarily confirms their willingness to participate in a study, after being informed of all relevant aspects.

### **GCP (Good Clinical Practice)**
An international ethical and scientific quality standard for designing, conducting, and reporting clinical trials (ICH E6).

### **IRB / Ethics Committee**
An institutional review board or ethics committee that reviews and approves research involving human subjects.

---

## **Data Standards and Formats**

### **DICOM**
Digital Imaging and Communications in Medicine — the universal standard for medical images, containing both pixel data and extensive metadata.

### **NIfTI**
Neuroimaging Informatics Technology Initiative — a simplified imaging format used in brain imaging research after anonymization from DICOM.

### **HL7 FHIR**
Fast Healthcare Interoperability Resources — a modern standard for exchanging electronic health records using RESTful APIs and JSON/XML resources.

### **OMOP CDM**
Observational Medical Outcomes Partnership Common Data Model — a standardized relational model for harmonizing clinical data from heterogeneous sources.

### **CDISC**
Clinical Data Interchange Standards Consortium — a family of standards for clinical trial data (CDASH, SDTM, ADaM).

### **SDTM**
Study Data Tabulation Model — a CDISC standard for organizing clinical trial data into standardized domains (DM, AE, LB, VS, etc.).

### **ADaM**
Analysis Data Model — a CDISC standard for analysis-ready datasets derived from SDTM data.

### **BIDS**
Brain Imaging Data Structure — a standard for organizing neuroimaging data in a consistent directory structure.

---

## **Medical Terminology Systems**

### **ICD (International Classification of Diseases)**
The WHO standard for coding diagnoses and health conditions. Current versions: ICD-10 (widely used) and ICD-11 (latest).

### **SNOMED CT**
Systematized Nomenclature of Medicine — Clinical Terms. The most comprehensive clinical terminology covering diseases, procedures, and substances.

### **LOINC**
Logical Observation Identifiers Names and Codes — a universal coding system for laboratory tests and clinical observations.

### **MedDRA**
Medical Dictionary for Regulatory Activities — a standardized terminology for coding adverse events and medical history in clinical trials.

### **MeSH**
Medical Subject Headings — a controlled vocabulary maintained by the NLM for indexing biomedical literature.

### **ATC**
Anatomical Therapeutic Chemical classification — a WHO system for classifying drugs by their therapeutic use and chemical properties.

### **HPO**
Human Phenotype Ontology — a standardized vocabulary for describing phenotypic abnormalities, used in rare disease research.

---

## **Privacy and Ethics**

### **GDPR**
General Data Protection Regulation — the EU regulation governing the processing of personal data, including health data.

### **Anonymization**
The irreversible process of removing all identifying information so that an individual cannot be re-identified from the data.

### **Pseudonymization**
Replacing identifying information with artificial identifiers (pseudonyms), such that re-identification is possible only with a separate key.

### **De-identification**
The general process of removing or obscuring personal identifiers from data. Encompasses both anonymization and pseudonymization.

### **Data Use Agreement (DUA)**
A legal agreement specifying the conditions under which a dataset may be accessed and used for research.

### **Federated Learning**
A machine learning approach that trains models across multiple institutions without sharing raw patient data — the algorithm travels to the data, not the other way around.

---

## **Biostatistics**

### **Primary Endpoint**
The main outcome measure of a clinical trial, pre-specified in the protocol.

### **Power Analysis**
A statistical calculation to determine the minimum sample size needed to detect a clinically meaningful effect with sufficient probability.

### **Kaplan-Meier Estimator**
A non-parametric method for estimating survival probabilities from time-to-event data, accounting for censored observations.

### **Hazard Ratio**
A measure of the relative risk of an event (e.g., death, recurrence) in the treatment group compared to the control group over time.

### **Intention-to-Treat (ITT)**
An analysis that includes all randomized participants regardless of whether they completed the study protocol, preserving the benefits of randomization.

### **Per-Protocol Analysis**
An analysis restricted to participants who completed the study according to the protocol. Less conservative than ITT.

### **Multiple Testing Correction**
Statistical adjustment for performing many simultaneous tests to control the false discovery rate. Common methods: Bonferroni, Benjamini-Hochberg.

### **Censoring**
In survival analysis, an observation is censored when the event of interest has not occurred by the end of the study period or the participant is lost to follow-up.

---

## **Medical Imaging**

### **CT (Computed Tomography)**
An imaging technique that uses X-rays to create cross-sectional images of the body. Measured in Hounsfield Units (HU).

### **MRI (Magnetic Resonance Imaging)**
An imaging technique using strong magnetic fields and radio waves to produce detailed images of soft tissues.

### **Hounsfield Unit (HU)**
The unit of measurement for CT image intensity. Water = 0 HU, air = -1000 HU, bone = +400 to +1000 HU.

### **Segmentation**
The process of delineating anatomical structures or regions of interest in a medical image.

### **PACS**
Picture Archiving and Communication System — the hospital system for storing, retrieving, and distributing medical images.

### **Radiomics**
The extraction of quantitative features (shape, texture, intensity) from medical images for use in predictive modeling.

---

**Need clarification on a term?** Contact us at [rdhub@mdsi.tum.de](mailto:rdhub@mdsi.tum.de) or explore our [Use Cases](/domain-knowledge/medical-clinical/use-cases) for practical examples.
