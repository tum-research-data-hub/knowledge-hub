---
title: Clinical Trial Data Management
slug: /domain-knowledge/medical-clinical/use-cases/clinical-trial-data
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Use Case 2: Clinical Trial Data Management**

Clinical trials generate complex, longitudinal datasets that must comply with regulatory requirements (ICH-GCP, GDPR), follow standardized formats (CDISC), and undergo rigorous quality control. This use case covers the data management lifecycle from study design through analysis-ready datasets.

## **Workflow Overview**

```
Study design → Data collection → Cleaning & coding → Analysis datasets → Reporting
Protocol/CRF → REDCap/EDC entries → SDTM domains → ADaM datasets → Tables, Figures, Listings
```

## **Key Concepts**

### **Electronic Data Capture (EDC)**

Modern clinical trials collect data electronically using systems like:
- **REDCap** — Free, widely used in academic research
- **OpenClinica** — Open-source, GCP-compliant
- **Medidata Rave** — Commercial, industry standard
- **Castor EDC** — Cloud-based, European-hosted

### **CDISC Standards**

The Clinical Data Interchange Standards Consortium defines how clinical trial data should be structured:

| Standard | Purpose | Stage |
|----------|---------|-------|
| **CDASH** | Standardized data collection forms | Data entry |
| **SDTM** | Study Data Tabulation Model — standardized domains (DM, AE, LB, VS, etc.) | Tabulation |
| **ADaM** | Analysis Data Model — analysis-ready datasets | Statistical analysis |
| **Define-XML** | Machine-readable data dictionary describing variables | Documentation |

### **Data Quality**

Clinical data undergoes rigorous quality control:
- **Edit checks** — Automated validation rules in the EDC system
- **Query management** — Data managers raise queries for inconsistent or missing data
- **Source data verification (SDV)** — Comparing EDC entries against source documents
- **Medical coding** — Standardizing adverse events (MedDRA) and medications (WHO Drug Dictionary)

### **Regulatory Requirements**

- **ICH E6 (R2) / GCP** — Good Clinical Practice guidelines for data integrity
- **21 CFR Part 11** (FDA) — Requirements for electronic records and signatures
- **GDPR** — EU data protection for patient data
- **Data retention** — Trial data must typically be retained for 15-25 years

## **Popular Tools**

### **Data Capture**
- **REDCap** — Secure web-based EDC, strong academic community
- **OpenClinica** — Open-source, CDISC-aware
- **LimeSurvey** — Open-source for surveys and patient-reported outcomes

### **Data Management and Analysis**
- **R** (tidyverse, haven, admiral) — CDISC dataset creation and analysis
- **SAS** — Traditional pharmaceutical industry standard
- **Python** (pandas, lifelines) — Data cleaning and survival analysis
- **Pinnacle 21** — CDISC validation tool

### **Randomization and Monitoring**
- **randomizeR** (R) — Randomization procedures
- **DSMB tools** — Interim analysis and monitoring

## **Code Example: REDCap Data Processing with R**

<details>
<summary>**Click to expand R clinical data processing script**</summary>

```r
#!/usr/bin/env Rscript
# Clinical trial data processing: REDCap export to analysis-ready dataset

library(tidyverse)
library(REDCapR)
library(lubridate)
library(survival)
library(survminer)

set.seed(42)

# ============================================================
# Step 1: Import data from REDCap
# ============================================================

# Option A: API access (preferred for reproducibility)
# redcap_data <- redcap_read(
#   redcap_uri = "https://redcap.institution.de/api/",
#   token = Sys.getenv("REDCAP_TOKEN")
# )$data

# Option B: CSV export (manual download)
raw_data <- read_csv("redcap_export.csv", show_col_types = FALSE)

cat("Imported", nrow(raw_data), "records with", ncol(raw_data), "variables\n")

# ============================================================
# Step 2: Data cleaning and validation
# ============================================================

clean_data <- raw_data %>%
  # Remove test/dummy records
  filter(!str_detect(record_id, "TEST|DUMMY")) %>%
  # Parse dates
  mutate(
    enrollment_date = ymd(enrollment_date),
    birth_date = ymd(birth_date),
    age = as.numeric(difftime(enrollment_date, birth_date, units = "days")) / 365.25
  ) %>%
  # Validate ranges
  mutate(
    age_valid = between(age, 18, 100),
    weight_valid = between(weight_kg, 30, 300),
    height_valid = between(height_cm, 100, 250)
  )

# Report validation issues
validation_issues <- clean_data %>%
  filter(!age_valid | !weight_valid | !height_valid) %>%
  dplyr::select(record_id, age, age_valid, weight_kg, weight_valid,
                height_cm, height_valid)

if (nrow(validation_issues) > 0) {
  cat("WARNING:", nrow(validation_issues), "records with validation issues\n")
  write_csv(validation_issues, "data_queries.csv")
}

# ============================================================
# Step 3: Create SDTM-like domains
# ============================================================

# Demographics (DM) domain
dm <- clean_data %>%
  transmute(
    STUDYID = "TRIAL-001",
    USUBJID = paste0("TRIAL-001-", str_pad(record_id, 4, pad = "0")),
    SUBJID = record_id,
    AGE = round(age, 1),
    AGEU = "YEARS",
    SEX = case_when(
      sex == 1 ~ "M",
      sex == 2 ~ "F",
      TRUE ~ "U"
    ),
    RACE = race_label,
    ARM = treatment_arm,
    ARMCD = case_when(
      treatment_arm == "Treatment" ~ "TRT",
      treatment_arm == "Placebo" ~ "PBO",
      TRUE ~ "UNK"
    ),
    RFSTDTC = format(enrollment_date, "%Y-%m-%d")
  )

cat("DM domain:", nrow(dm), "subjects\n")
cat("  Treatment:", sum(dm$ARMCD == "TRT"), "\n")
cat("  Placebo:", sum(dm$ARMCD == "PBO"), "\n")

# ============================================================
# Step 4: Basic analysis — survival analysis example
# ============================================================

# Create analysis dataset
analysis <- dm %>%
  left_join(
    clean_data %>%
      transmute(
        USUBJID = paste0("TRIAL-001-", str_pad(record_id, 4, pad = "0")),
        event_time_days = as.numeric(event_time),
        event_status = as.integer(event_occurred)
      ),
    by = "USUBJID"
  ) %>%
  filter(!is.na(event_time_days))

# Kaplan-Meier survival analysis
surv_obj <- Surv(analysis$event_time_days, analysis$event_status)
km_fit <- survfit(surv_obj ~ ARM, data = analysis)

# Plot
pdf("kaplan_meier.pdf", width = 8, height = 6)
ggsurvplot(km_fit, data = analysis,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           xlab = "Time (days)",
           ylab = "Survival probability",
           title = "TRIAL-001: Kaplan-Meier Survival Curves",
           palette = c("#E41A1C", "#377EB8"))
dev.off()

# Log-rank test
logrank <- survdiff(surv_obj ~ ARM, data = analysis)
cat("Log-rank test p-value:", format.pval(1 - pchisq(logrank$chisq, 1)), "\n")

# ============================================================
# Step 5: Export datasets
# ============================================================

write_csv(dm, "output/sdtm_dm.csv")
write_csv(analysis, "output/adam_analysis.csv")

cat("Data processing complete!\n")
sessionInfo()
```

</details>

## **Expected Outputs**

- **SDTM datasets** — Standardized tabulation domains (DM, AE, LB, VS, etc.)
- **ADaM datasets** — Analysis-ready datasets with derived variables
- **Data quality reports** — Validation summaries, query logs
- **Define-XML** — Machine-readable data dictionary
- **Tables, Figures, Listings (TFLs)** — Statistical analysis outputs
- **Audit trail** — Complete log of data modifications

## **Computational Requirements**

| Task | CPU | RAM | Storage | Time |
|------|-----|-----|---------|------|
| REDCap data export and cleaning | 1-2 | 2-4 GB | `<1 GB` | Minutes |
| SDTM/ADaM dataset creation | 2-4 | 4-8 GB | `<1 GB` | Minutes-hours |
| Statistical analysis (standard) | 2-4 | 4-8 GB | `<1 GB` | Minutes |
| Complex modeling (Bayesian, ML) | 8-16 | 16-32 GB | 1-5 GB | Hours |

## **Common Issues & Troubleshooting**

:::warning Common Problems

**Missing data**
- Clinical trials always have missing data — document the extent and mechanism (MCAR, MAR, MNAR)
- Pre-specify imputation strategies in the Statistical Analysis Plan (SAP)
- Never impute primary endpoint data without strong justification

**Date inconsistencies**
- Different sites may use different date formats — standardize early
- Partial dates (e.g., only month/year) are common — define imputation rules in advance
- Use ISO 8601 format (YYYY-MM-DD) consistently

**Coding inconsistencies**
- Free-text entries for adverse events need medical coding (MedDRA)
- Medication names vary — map to WHO Drug Dictionary or ATC codes
- Establish a coding convention document before data collection begins

:::

## **Key Considerations**

- **Plan data management before enrollment** — Create a Data Management Plan (DMP) alongside the protocol
- **Use validated systems** — EDC systems must comply with 21 CFR Part 11 / GCP requirements
- **Separate roles** — Data managers, biostatisticians, and monitors should have distinct responsibilities
- **Lock data before analysis** — Freeze the database after cleaning and before unblinding
- **Archive everything** — Raw data, analysis code, outputs, and correspondence must be retained per regulatory requirements
- **Pre-register your analysis** — Publish the Statistical Analysis Plan before database lock to avoid bias
