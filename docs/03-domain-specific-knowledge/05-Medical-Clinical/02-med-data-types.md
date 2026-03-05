---
title: Medical & Clinical Data Types & Formats
slug: /domain-knowledge/medical-clinical/data-types
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Medical & Clinical Data Types & File Formats**

Medical and clinical research relies on a broad range of data types — from imaging data and structured health records to clinical trial databases and wearable sensor streams. These formats are shaped by strict interoperability requirements and privacy regulations.

---

## **Medical Imaging Formats**

### **DICOM (.dcm)**
Digital Imaging and Communications in Medicine — the universal standard for medical images. A DICOM file contains both the image data and extensive metadata (patient info, acquisition parameters, equipment details).
- **Use:** Radiology (CT, MRI, X-ray, PET), ultrasound, nuclear medicine
- **Tools:** 3D Slicer, OsiriX/Horos, OHIF Viewer, pydicom (Python)
- **Spec:** [DICOM Standard](https://www.dicomstandard.org/)
- **Note:** DICOM files contain patient-identifiable information and must be anonymized before sharing

### **NIfTI (.nii, .nii.gz)**
Neuroimaging Informatics Technology Initiative — a simplified imaging format widely used in neuroimaging research after anonymization from DICOM.
- **Use:** Brain MRI analysis, fMRI, DTI, research image sharing
- **Tools:** FSL, FreeSurfer, ANTs, nibabel (Python), 3D Slicer
- **Spec:** [NIfTI-1 Format](https://nifti.nimh.nih.gov/)

### **NRRD (.nrrd)**
Nearly Raw Raster Data — a flexible format for multi-dimensional raster data, commonly used in medical image analysis research.
- **Use:** Segmentation results, image analysis pipelines
- **Tools:** 3D Slicer, ITK, pynrrd (Python)
- **Spec:** [NRRD Format](http://teem.sourceforge.net/nrrd/format.html)

### **Whole Slide Images (WSI) — SVS, NDPI, MRXS**
Vendor-specific formats for digitized histopathology slides. These are typically multi-resolution pyramid images reaching several gigabytes per slide.
- **Use:** Digital pathology, computational pathology, AI-based diagnosis
- **Tools:** OpenSlide, QuPath, ASAP
- **Note:** OpenSlide provides a unified API for reading various WSI formats

---

## **Clinical and Health Record Formats**

### **HL7 FHIR (Fast Healthcare Interoperability Resources)**
A modern RESTful standard for exchanging electronic health records. Data is organized as discrete resources (Patient, Observation, Condition, etc.) in JSON or XML.
- **Use:** EHR data exchange, clinical data integration, mobile health apps
- **Tools:** HAPI FHIR (Java), fhirpy (Python), SMART on FHIR
- **Spec:** [HL7 FHIR](https://www.hl7.org/fhir/)

### **HL7 v2 / CDA**
Older health information exchange standards still widely used in hospital systems. HL7v2 uses pipe-delimited messages; CDA uses XML clinical documents.
- **Use:** Lab results, admission/discharge messages, clinical summaries
- **Note:** Being gradually replaced by FHIR in new implementations

### **OMOP CDM (Observational Medical Outcomes Partnership Common Data Model)**
A standardized relational data model that maps heterogeneous clinical data into a common structure, enabling cross-institutional observational research.
- **Use:** Observational studies, pharmacovigilance, real-world evidence
- **Tools:** OHDSI tools (ATLAS, Achilles, HADES)
- **Spec:** [OMOP CDM](https://ohdsi.github.io/CommonDataModel/)

### **CDISC Standards (CDASH, SDTM, ADaM)**
A family of standards for clinical trial data:
- **CDASH** — Clinical Data Acquisition Standards Harmonization (data collection)
- **SDTM** — Study Data Tabulation Model (regulatory submission)
- **ADaM** — Analysis Data Model (statistical analysis)
- **Use:** Clinical trial data management and regulatory submissions (FDA, EMA)
- **Spec:** [CDISC](https://www.cdisc.org/)

### **REDCap Export Formats (.csv, .r, .sps)**
REDCap (Research Electronic Data Capture) exports data in CSV along with statistical software syntax files for labeling variables.
- **Use:** Survey data, clinical study databases, registries
- **Tools:** REDCap, R (REDCapR), Python (PyCap)

---

## **Biosignal and Wearable Data**

### **EDF / EDF+ (.edf)**
European Data Format — a standard for multi-channel biosignal recordings.
- **Use:** EEG, EMG, ECG, polysomnography, sleep studies
- **Tools:** MNE-Python, EDFbrowser, EEGLAB
- **Spec:** [EDF Spec](https://www.edfplus.info/)

### **WFDB (.hea, .dat)**
WaveForm DataBase format — used by PhysioNet for physiological signal data.
- **Use:** ECG, heart rate variability, physiological signal databases
- **Tools:** WFDB (Python/MATLAB), PhysioNet
- **Spec:** [PhysioNet WFDB](https://physionet.org/content/wfdb/)

### **GDT / Accelerometer CSV**
Various proprietary and open formats for wearable device data (accelerometers, gyroscopes, heart rate monitors).
- **Use:** Activity recognition, remote patient monitoring
- **Note:** No universal standard; most devices export CSV or JSON

---

## **Genomic and Molecular Data in Clinical Context**

Clinical genomics uses many formats shared with bioinformatics (VCF, FASTQ, BAM). Additionally:

### **Phenopackets (.json)**
A GA4GH standard for sharing disease and phenotype information associated with a patient or sample.
- **Use:** Rare disease diagnosis, clinical genomics, phenotype-genotype association
- **Spec:** [GA4GH Phenopackets](https://phenopacket-schema.readthedocs.io/)

---

## **Format Selection Guide**

| Use Case | Recommended Format |
|----------|-------------------|
| Medical imaging (clinical) | DICOM |
| Neuroimaging (research) | NIfTI |
| Digital pathology | SVS / OpenSlide-compatible |
| Health record exchange | HL7 FHIR |
| Observational research | OMOP CDM |
| Clinical trial data | CDISC (SDTM/ADaM) |
| Study databases | REDCap |
| EEG/ECG signals | EDF+ |
| Clinical genomics phenotypes | Phenopackets |
