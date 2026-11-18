---
title: Metadata
slug: /general-knowledge/metadata
sidebar_position: 3
---

# **Metadata**

Metadata is often defined as "data about data." It refers to the structured information used to describe a data point or a dataset, covering its structure, content, and unique features. High-quality metadata is the key to ensuring your research remains findable and usable by both systems and other researchers (including your future self).

## **1\. Why Metadata Matters**

Metadata is the backbone of the **FAIR (Findable, Accessible, Interoperable, and Reusable)** principles. Without it, a dataset is just a collection of bits with no context. It helps you:

* **Organise:** Sort through thousands of files using searchable tags.  
* **Manage:** Track versions, authors, and modifications over time.  
* **Reproduce:** Provide the exact parameters and settings used to generate results.

## **2\. Common Metadata Examples**

Metadata exists at different levels of complexity. A simple example is an image file from your phone, which automatically includes:

* **File Metadata:** Format (.jpg), size (4MB), resolution (4032x3024).  
* **Contextual Metadata:** Date created, GPS coordinates, and camera settings (ISO, aperture).

In research, we often categorise metadata into three types:

1. **Descriptive:** For finding and identifying (e.g., Title, Author, Keywords, Abstract).  
2. **Structural:** How the data is put together (e.g., How different files in a folder relate to one another).  
3. **Administrative:** For technical management (e.g., File formats, License types, Access permissions).

## **3\. Metadata Standards and Schemas**

To make metadata "Interoperable," researchers use **Schemas**—standardised sets of fields that everyone in a field agrees to use.

* **Dublin Core:** A basic, universal set of 15 elements (Title, Creator, Date, etc.) used by repositories like **mediaTUM**.  
* **Domain-Specific Standards:** Fields like Biology or Physics often have complex requirements (e.g., the "Minimum Information" standards in genomics).  
  * 🧬 **See Detail:** [Metadata & Ontologies (Bioinformatics)](/domain-knowledge/bioinformatics/bio-metadata-ontologies)


## **4\. Data Management Plans (DMPs)**

A **Data Management Plan (DMP)** describes how you will handle your research data during and after your project. Many funding agencies now require a DMP as part of your proposal — it demonstrates that you have thought through data collection, storage, documentation, and sharing.

### **Why are DMPs important?**

Beyond funder requirements, a DMP helps you:
- Plan storage and backup needs early
- Clarify responsibilities within your team
- Identify metadata standards for your discipline
- Consider legal and ethical aspects before data collection begins

### **The DFG Checklist: What Funders Want to Know**

The [DFG provides a checklist](https://www.dfg.de/resource/blob/174736/forschungsdaten-checkliste-en.pdf) that outlines the key questions any good DMP should answer:

1. **Data description** — What types and formats of data will be generated or collected?
2. **Documentation & quality** — How will you document and ensure the quality of your data?
3. **Storage & backup** — Where will data be stored during the project and how will it be backed up?
4. **Legal & ethical aspects** — What legal requirements (e.g., data protection, licensing) apply?
5. **Sharing & archiving** — How and where will data be made accessible after the project?
6. **Responsibilities** — Who in your team is responsible for data management?

These six questions essentially define what a DMP covers. Answering them thoroughly demonstrates good research data management practice.

**Further information:** [TUM Research Data Hub — Data Management Plans](https://web.tum.de/en/researchdata/support-information/submit-application-planning/)

---

## **5\. Best Practices & Management Tools**

Capturing metadata manually can be time-consuming. TUM provides and recommends several tools to streamline this process:

* **Start Early:** Don't wait until the end of your project. Document your data as you create it.  
* **Use README files:** For simple projects, a well-structured README.txt file in your root directory acts as a manual for your metadata.  
* [**TUM eLabFTW**](https://www.ub.tum.de/en/eln)**:** The official Electronic Lab Notebook (ELN) for TUM researchers. It allows you to document experiments in a structured way, automatically capturing timestamps and linking datasets to specific experimental parameters.  
* [**TUM DataTagger**](https://www.ub.tum.de/en/datatagger)**:** An open-source tool designed to help you annotate and version your data collaboratively. It simplifies the process of adding descriptive metadata to datasets before they are shared or archived.

**Note:** Consistent metadata is what turns "raw data" into "scientific evidence." If you are unsure which standard to use for your metadata or what tools are best suited for your research needs, you can send us an email at [researchdata@tum.de](mailto:researchdata@tum.de).