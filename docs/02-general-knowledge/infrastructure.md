---
title: Computing Infrastructure
slug: /general-knowledge/infrastructure
---

# **Computing Infrastructure**

This page provides an overview of the different types of computing and storage infrastructure available to researchers at the Technical University of Munich. Whether you need high-performance computing resources, cloud-based solutions, or specialized AI systems, TUM and its partners offer a wide range of options to support your research.

At the end of this page, you will find a decision guide to help you choose the right solution based on your computing specifications (CPU, RAM, GPU) and research requirements.

## **Available Infrastructure**

| Name | Type of Service | Access Requirement | Intended Use | Storage/Transfer Options | Link |
|------|----------------|-------------------|--------------|-------------------------|------|
| **LRZ Linux Cluster** | HPC Cluster | TUM affiliation, application required | High-performance computing, parallel processing | LRZ DSS | [LRZ Linux Cluster](https://doku.lrz.de/linux-cluster-10746275.html) |
| **LRZ Compute Cloud** | Virtual Machines | TUM affiliation, project application | Flexible computing environments, web services | LRZ DSS, object storage| [LRZ Compute Cloud](https://doku.lrz.de/compute-cloud-10746035.html) |
| **LRZ AI Systems** | GPU Cluster | Application required, peer review | AI/ML training, deep learning | Specialized storage, LRZ DSS | [LRZ AI Systems](https://doku.lrz.de/ai-systems-11481865.html) |
| **LRZ Quantum Computing** | Quantum Systems Cluster | Special application | Quantum algorithm research | LRZ DSS, Globus | [LRZ Quantum](https://www.lrz.de/wir/newsletter/2021-06/#quantum) |
| **Terrabyte** | HPC Cluster | TUM account, DFL employee | Earth observational data science | LRZ DSS, Globus | [Terrabyte](https://www.lrz.de/services/compute/linux-cluster/terrabyte/) |
| **EOSC** | Virtual Machines, Jupyter Notebooks | European research affiliation | Flexible computing environments | Federated storage | [EOSC Portal](https://eosc-portal.eu/) |
| **Jupyter4NFDI** | Jupyter Notebooks | Affiliation with German University | Interactive data analysis | Integrated storage | [Jupyter4NFDI](https://jupyter.nfdi.de/) |
| **NHR** | National HPC | Doctoral researchers can apply | Large-scale simulations | High-throughput storage | [NHR Alliance](https://www.nhr-verein.de/) |
| **EuroHPC** | European HPC | Competitive application, calls open regularly | Extreme-scale computing | Parallel file systems | [EuroHPC JU](https://www.eurohpc-ju.europa.eu/index_en) |
| **AWS** | Commercial Cloud | Credit card/budget | Scalable cloud computing | S3, EBS storage | [AWS](https://aws.amazon.com/) |
| **Microsoft Azure** | Commercial Cloud | Credit card/budget | Cloud services, integration with M365 | Blob storage, disk storage | [Azure](https://azure.microsoft.com/) |
| **Google Cloud** | Commercial Cloud | Credit card/budget | ML/AI services, BigQuery | Cloud storage, persistent disks | [Google Cloud](https://cloud.google.com/) |
| **Local Workstation** | Personal Computer | Direct access | Development, small-scale analysis | Local disks, network drives | N/A |
| **WSL (Windows)** | Linux on Windows | Windows 10/11 | Linux tools on Windows | Windows file system integration | [WSL Docs](https://learn.microsoft.com/en-us/windows/wsl/) |

## **TUM-Specific Resources**

### **Leibniz Supercomputing Centre (LRZ)**

The LRZ is TUM's primary partner for high-performance computing and provides several tiers of computing and storage solutions:

* **Linux Cluster:** General-purpose HPC for parallel computing tasks
* **Compute Cloud:** Virtualized infrastructure for flexible deployments
* **AI Systems:** State-of-the-art GPU clusters for artificial intelligence research
* **Storage Solutions:** Tiered storage from "hot data" to long-term archiving

**Support:** The LRZ offers comprehensive documentation, helpdesk support, and regular training courses. 
**Overview:** The LRZ wiki contains information on all services provided: [LRZ Wiki, Overview of Computing options](https://doku.lrz.de/high-performance-computing-10613431.html)
**Consultation:** Contact the [LRZ Servicedesk](https://servicedesk.lrz.de/) for personalized advice on which service fits your needs.

### **NFDI (National Research Data Infrastructure)**

Germany's National Research Data Infrastructure provides domain-specific services:

* **Jupyter4NFDI:** Interactive computing environments
* **NFDI4Ing:** Infrastructure for engineering sciences
* **Domain-specific tools:** Specialized services for different research fields

### **Major Research Instrumentation (Großgeräteanträge)**

For research groups requiring dedicated hardware, TUM supports applications for major research instrumentation through DFG funding programs.

## **Decision Guide**

### **How to Choose the Right Infrastructure?**

Consider these key factors:

#### **1. Computational Requirements**

* **CPU-intensive tasks** (simulations, modeling): LRZ Linux Cluster, NHR
* **GPU-intensive tasks** (AI/ML, deep learning): LRZ AI Systems
* **Memory-intensive tasks** (large datasets in RAM): High-memory nodes at LRZ
* **Interactive analysis**: LRZ Compute Cloud, Jupyter4NFDI

#### **2. Data Volume and Storage**

* **< 1 TB:** Local workstation, LRZ Compute Cloud
* **1-10 TB:** LRZ Linux Cluster storage
* **> 10 TB:** Terrabyte, specialized storage solutions
* **Long-term archiving:** Data Science Archive (LRZ), mediaTUM

#### **3. Data Sensitivity**

* **Public data:** Any service
* **Sensitive/Personal data:** LRZ services with data protection agreement
* **Medical/Human data:** Specialized secure environments, consult [researchdata@tum.de](mailto:researchdata@tum.de)

#### **4. Budget Considerations**

* **Free (for TUM researchers):** LRZ services, NFDI services
* **Application-based:** NHR, EuroHPC (free but competitive)
* **Commercial:** AWS, Azure, Google Cloud (pay-as-you-go)

#### **5. Expertise Level**

* **Beginners:** Jupyter4NFDI, LRZ Compute Cloud (with GUI)
* **Intermediate:** LRZ Linux Cluster (batch systems)
* **Advanced:** Custom setups, cloud infrastructure

## **Getting Started**

1. **Assess your needs:** Determine your computational, storage, and data protection requirements
2. **Consult documentation:** Review the linked resources for each service
3. **Request access:** Follow the application procedures for your chosen service
4. **Get training:** Attend LRZ courses or contact the TUM Research Data Hub for guidance

## **Support and Feedback**

For questions about which infrastructure is right for your project:
* **General RDM support:** [researchdata@tum.de](mailto:researchdata@tum.de)
* **LRZ Services:** [LRZ Servicedesk](https://servicedesk.lrz.de/)
* **Domain-specific guidance:** Check our [Domain-Specific Knowledge](/domain-knowledge/bioinformatics/bio-data-types) sections

---

**Note:** Infrastructure availability and access procedures may change. Always check the official websites for the most current information.

