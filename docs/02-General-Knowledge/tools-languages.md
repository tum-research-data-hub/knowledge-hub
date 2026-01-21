---
title: Tools and Languages
---

# **Tools and Programming Languages**

Selecting the right tool depends on your research goals, data volume, and the need for reproducibility. At TUM, we support a range of open-source and proprietary solutions for data science and statistical analysis.

## **1\. Programming Languages**

Programming-based approaches are the 'gold standard' for reproducibility, as every step of your analysis is documented in code.

### **R**

R is a language specifically designed for statistical computing and graphics.

* **Best for:** Complex statistical modelling and publication-quality visualisations.  
* **TUM Context:** Widely used in Bioinformatics and Social Sciences.  
* **Data Visualisation with ggplot2:**  
  * The ggplot2 package is the industry standard for R visualisations.  
  * It is based on the **'Grammar of Graphics,'** which allows you to build plots layer by layer (data \+ aesthetics \+ geometric objects).  
* **Core Resource:** We highly recommend [**Gagneur Lab: Data Analysis and Visualisation in R**](https://gagneurlab.github.io/dataviz/).

### **Python**

Python is a general-purpose language that has become the dominant tool for Machine Learning, Artificial Intelligence, and Data Engineering.

* **Package Management with pip:**  
  * pip stands for **'Pip Installs Packages'** (a recursive acronym). It is the standard package manager for Python.  
  * **Best Practice:** Always use **Virtual Environments** (venv or conda) to keep your research project dependencies isolated.  
* **Essential Libraries:** Pandas (data), Scikit-learn (ML), PyTorch (AI), and Matplotlib (viz).

### **C and C++**

These are compiled languages used when performance is the highest priority.

* **Best for:** Developing high-performance software, simulations, and compute-intensive algorithms that need to run on HPC clusters.  
* **TUM Context:** Core to Engineering, Physics, and the development of new bioinformatics tools (like sequence aligners).  
* **Note:** While harder to learn than Python, they offer significantly faster execution speeds for massive datasets.

### **JavaScript**

JavaScript is the language of the web, increasingly used in research for sharing results.

* **Best for:** Creating interactive web-based visualisations and data dashboards.  
* **Key Libraries:** D3.js (complex data-driven documents) and React (building user interfaces).  
* **Usage:** Essential if you are building a web tool for other researchers to explore your findings.

### **SQL**

Structured Query Language (SQL) is the standard for managing and querying data held in relational database management systems.

* **Best for:** Interacting with large enterprise data warehouses and structured datasets.

## **2\. Statistical and Scientific Software**

* **MATLAB:** High-level language for numerical computation and engineering.  
* **SPSS:** Powerful platform used in the social sciences for survey analysis.  
* **Minitab:** Often used for quality improvement and teaching statistical concepts.

## **3\. Data Visualisation and Business Intelligence**

* **Microsoft Power BI:** Tool for building interactive reports with live data refresh.  
* **Tableau:** Intuitive 'drag-and-drop' interface for exploring large datasets.

## **4\. Spreadsheets**

* **Excel:** Common for quick entry but generally not recommended for complex, reproducible data pipelines.

## **5\. Summary Table**

| Name | Usage | Type | Description | Link (TUM Internal) | Link (Doc) |
| :---- | :---- | :---- | :---- | :---- | :---- |
| **R** | Statistics / Viz | Open Source | Statistical computing and graphics. | [TUM R-Kurs](https://www.google.com/search?q=https://www.groups.tum.de/stat/lehre/) | [R Project](https://www.r-project.org/) |
| **Python** | AI / Data Sci | Open Source | General purpose / Data Engineering. | [LRZ Python Docs](https://www.google.com/search?q=https://doku.lrz.de/python-10745934.html) | [Python.org](https://www.python.org/) |
| **C / C++** | HPC / Performance | Open Source | High-speed compiled languages. | [LRZ HPC Intro](https://www.google.com/search?q=https://doku.lrz.de/hpc-10745941.html) | [isocpp.org](https://isocpp.org/) |
| **JavaScript** | Web / Interactivity | Open Source | Interactive web visualisations. | [TUM Online](https://campus.tum.de/) | [MDN Docs](https://developer.mozilla.org/en-US/docs/Web/JavaScript) |
| **MATLAB** | Computation | Proprietary | Numerical computing and engineering. | [TUM Software](https://www.it.tum.de/it/software/) | [MathWorks](https://www.mathworks.com/) |
| **SPSS** | Social Science | Proprietary | Advanced statistical analysis. | [BayernCollab](https://www.google.com/search?q=https://bayerncollab.uni-wuerzburg.de/) | [IBM SPSS](https://www.ibm.com/spss) |
| **Power BI** | BI / Dashboards | Proprietary | Interactive reports and data refresh. | [TUM IT Services](https://www.google.com/search?q=https://www.it.tum.de/it/m365/) | [MS Learn](https://learn.microsoft.com/power-bi/) |

**Note:** For a full list of software available under campus licences, visit the [TUM Software Portal](https://www.it.tum.de/it/software/).