---
title: Simulation & Modeling Metadata & Standards
slug: /domain-knowledge/simulation-modeling/metadata-standards
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Simulation & Modeling Metadata & Standards**

Reproducing a simulation requires far more than just the output data — the exact software version, input parameters, mesh resolution, boundary conditions, and hardware used must all be documented. Metadata standards in this domain aim to capture this information systematically.

---

## **Metadata Standards and Frameworks**

### **MODA (Model Data)**
A semi-formal description framework for documenting simulation workflows, developed within the European Materials Modelling Council (EMMC).
- **Use:** Documenting simulation workflows in materials modeling
- **Spec:** [EMMC MODA](https://emmc.info/moda-workflow/)

### **W3C PROV (Provenance)**
A W3C standard for representing provenance information — who did what, when, and how. Can be applied to any computational workflow.
- **Use:** Recording simulation lineage, parameter history, data transformations
- **Spec:** [W3C PROV](https://www.w3.org/TR/prov-overview/)

### **CodeMeta**
A metadata schema for research software, enabling proper citation and discovery of simulation codes.
- **Use:** Documenting simulation software for citation and reproducibility
- **Spec:** [CodeMeta](https://codemeta.github.io/)

### **CGNS Metadata**
The CFD General Notation System includes rich metadata support for mesh topology, boundary conditions, flow solutions, and convergence histories.
- **Use:** CFD data exchange and archival
- **Spec:** [CGNS Standard](https://cgns.github.io/)

---

## **Essential Simulation Metadata**

Every simulation dataset should document:

### **Software Environment**
- Solver name and version
- Compiler and version (e.g., GCC 12.2, Intel oneAPI 2023)
- MPI implementation and version
- Operating system and HPC environment

### **Input Parameters**
- Governing equations and models (e.g., RANS, LES, DNS for turbulence)
- Material properties and constitutive models
- Boundary and initial conditions
- Time stepping scheme and step size

### **Mesh Information**
- Element types and count
- Mesh quality metrics (skewness, aspect ratio, orthogonality)
- Refinement strategy (uniform, adaptive, local)

### **Computational Resources**
- Number of cores/nodes used
- Wall-clock time and CPU hours
- Memory usage

### **Convergence and Validation**
- Residual history
- Grid convergence study results
- Comparison with experimental or analytical reference data

---

## **Best Practices**

1. **Version-control input files** — Use Git for solver input files, configuration, and post-processing scripts
2. **Automate metadata capture** — Script the extraction of solver version, run parameters, and compute resources
3. **Store mesh and results together** — Keep mesh, input, and output files in a coherent directory structure
4. **Document convergence** — Include residual plots and grid independence studies with published results
5. **Use DVC or similar** — For large binary files that don't fit in Git, use Data Version Control
6. **Create README files** — Each simulation case directory should contain a README describing the setup
