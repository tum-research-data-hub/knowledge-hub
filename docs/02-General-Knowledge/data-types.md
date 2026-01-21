---
title: Data Types
slug: /general-knowledge/data-types
---

# **Data Types**

Depending on your research questions, you may encounter various categories of data. **Data type** describes the kind of data and defines how it is stored, what operations can be applied to it, and what analysis or models are valid.

Understanding data types is essential because it affects every step of your entire research process, and data misclassification may cause errored results and flawed conclusions.

## **1\. Universal Data Types**

The most common data types across almost all research disciplines include:

### **Numeric Data (Integers and Floats)**

Numeric data refer to quantities that can be measured or counted.

* **Integers:** Whole numbers (e.g., number of participants).  
* **Floats:** Numbers with decimal points (e.g., temperature, pH levels).  

Numeric data is the foundation of quantitative analysis; storing numbers in the correct format is crucial for valid modelling and calculations.


### **Categorical Data**

Represent values that fall into groups or categories without any specific numerical meaning.

* **Example:** Experimental condition groups (e.g., ‘High’, ‘Medium’, ‘Low’).


### **Ordinal Data**

Similar to categorical data but possess an inherent order or ranking.

* **Example:** Ordered survey responses:  
  * ‘Definitely agree’ --> 5  
  * ‘Mostly agree’ --> 4  
  * ‘Neither agree nor disagree’ --> 3  
  * ‘Mostly disagree’ --> 2  
  * ‘Definitely disagree’ --> 1


### **Text or String Data**

Typically includes free-format text, names, labels, and descriptions. 

Text is also widely used for metadata and documentation.

### **Binary or Boolean Data**

Refer to yes/no or true/false values. 

Booleans are essential for filtering datasets, building logical conditions in programming, and documenting quality-control checks.

* **Example:** ‘Did the robot detect an obstacle?’ --> False


### **Date and Time Data**

Record temporal information crucial for examining trends, sequences, or rhythms. Proper handling is essential when dealing with international date formats and different time zones.

* **Example:** Timestamps (2025-01-03 10:15:00)

##

## **2\. Specialised Data Types**

Depending on your research subject and area, you might also see special data types, including:

### **Matrix and Vector** 

Structured collections of numeric data. A **Vector** is a one-dimensional array, whereas **Matrices** are two-dimensional.  

### **Geospatial Data**

Information linked to a specific location on Earth. It includes location information (coordinates), attribute information (details about the location), and temporal information.  

### **Image Data** 
Visual representations (e.g., microscopy, satellite imagery) requiring specific archival formats like TIFF or DICOM.  

### **Genomic Data** 

Highly specialised sequences and variants used in bioinformatics.  
  
   🧬 **See Detail:** [Domain-specific Data Types (Bioinformatics)](/domain-knowledge/bioinformatics/bio-data-types)

**Note:** This is not an exhaustive list. Always consult your Data Steward or data expert in your area to ensure your data is classified correctly for long-term storage and analysis.