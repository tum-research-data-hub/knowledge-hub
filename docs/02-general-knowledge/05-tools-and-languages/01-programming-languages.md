---
title: Programming Languages
slug: /general-knowledge/tools/programming-languages
---

# **Programming Languages**

Programming-based approaches are the 'gold standard' for reproducibility, as every step of your analysis is documented in code.

## **R**

R is a language specifically designed for statistical computing and graphics.

* **Best for:** Complex statistical modelling and publication-quality visualisations.  
* **TUM Context:** Widely used in Statistics, Bioinformatics, Life Sciences and Social Sciences.  
* **Data Visualisation with ggplot2:**  
  * The ggplot2 package is the industry standard for R visualisations.  
  * It is based on the **'Grammar of Graphics'**, which allows you to build plots layer by layer (data + aesthetics + geometric objects).  
* **Core Resource:** We highly recommend [**Gagneur Lab: Data Analysis and Visualisation in R**](https://gagneurlab.github.io/dataviz/).

## **Python**

Python is a general-purpose language that has become the dominant tool for Machine Learning, Artificial Intelligence, and Data Engineering.

* **Package Management with pip and uv:**  
  * pip stands for **'Pip Installs Packages'** (a recursive acronym). It is the standard package manager for Python.
  * [uv](https://docs.astral.sh/uv/) has established itself as the superior and vastly faster alternative to pip as a package installer, venv as a virtual environment manager and other tools such as pyenv or poetry as dependency managers.
* **Package Management with conda:**
  * [Conda](https://docs.conda.io/) is a package and environment manager that is not limited to Python: it also installs R, command-line tools and the compiled libraries they depend on. This makes it the usual choice when a project mixes languages or needs non-Python software.
  * Packages come from *channels*. The most important ones for research are the community channel [conda-forge](https://conda-forge.org/) and [Bioconda](https://bioconda.github.io/), which provides several thousand bioinformatics tools.
* **Dependency management:**
  * most code is dependent on other libraries / Python packages
  * as libraries change, code that depends on these libraries might become incorrect over time. It is therefore essential to not only know on *which* libraries your code depends on, but also *what version* of these dependencies it assumes.
  * the Python Foundation has developed an open standard to denote dependencies in a file called [`pyproject.toml`](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/). It contains (alongside other metadata) a list of dependencies and their (usually minimal) version requirements
  * `pyproject.toml` does not need to be handwritten. Tools such as uv manage it for you.
  * A **conda** environment is described in an `environment.yml` file that lists the channels and packages (ideally with version numbers). `conda env create -f environment.yml` recreates the environment on another machine; keep the file under version control with your code.
  * a **lock file** is a file that contains the exact version numbers and checksums of all currently installed dependencies including all transitive dependencies. Such a file is necessary in order to exactly reproduce a Python environment and should definitely be included in scientific analysis code. The lock file uv produces is called `uv.lock`. For conda, `conda env export` or [conda-lock](https://conda.github.io/conda-lock/) serve the same purpose, because an `environment.yml` usually lists only the packages you asked for.
* **Best Practice:** Always use **Virtual Environments** (uv, venv or conda) to keep your research project dependencies isolated.
* **Essential Libraries:** Pandas (data), Polars (data, viz), Scikit-learn (ML), PyTorch (AI), and Matplotlib (viz).

## **C and C++**

These are compiled languages used when performance is the highest priority.

* **Best for:** Developing high-performance software, simulations, and compute-intensive algorithms that need to run on HPC clusters.  
* **TUM Context:** Core to Engineering, Physics, and the development of new bioinformatics tools (like sequence aligners).  
* **Note:** While harder to learn than Python, they offer significantly faster execution speeds for massive datasets. However, slower languages can still benefit from the speed of compiled C/C++ code via so-called *bindings*. Bindings allow you to make calls to compiled code from other languages such as Python or R. This is how Python libraries like numpy achieve their speed.

## **JavaScript**

JavaScript is the language of the web, increasingly used in research for sharing results.

* **Best for:** Creating interactive web-based visualisations and data dashboards.  
* **Key Libraries:** D3.js (complex data-driven documents) and React (building user interfaces).  
* **Usage:** Essential if you are building a web tool for other researchers to explore your findings.

## **SQL**

Structured Query Language (SQL) is the standard for managing and querying data held in relational database management systems.

* **Best for:** Interacting with large enterprise data warehouses and structured datasets.
