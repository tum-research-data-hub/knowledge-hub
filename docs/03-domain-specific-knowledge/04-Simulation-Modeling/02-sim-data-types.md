---
title: Simulation & Modeling Data Types & Formats
slug: /domain-knowledge/simulation-modeling/data-types
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Simulation & Modeling Data Types & File Formats**

Simulation and modeling workflows produce a variety of data types — from mesh definitions and solver input files to large time-dependent result datasets and post-processing outputs. This page covers the most important formats.

---

## **Mesh and Geometry Formats**

### **STL (Stereolithography, .stl)**
A simple surface mesh format using triangulated surfaces. Does not support volume elements or metadata.
- **Use:** 3D printing, CAD-to-simulation geometry transfer, surface visualization
- **Tools:** Gmsh, ParaView, MeshLab, FreeCAD
- **Note:** Available in ASCII and binary variants; binary is much more compact

### **STEP (.step, .stp)**
Standard for the Exchange of Product Model Data — an ISO standard for CAD geometry exchange.
- **Use:** CAD geometry exchange between software, simulation preprocessing
- **Tools:** FreeCAD, SALOME, CATIA, SolidWorks
- **Spec:** [ISO 10303](https://www.iso.org/standard/72237.html)

### **MSH (.msh)**
Gmsh's native mesh format supporting 1D, 2D, and 3D elements with physical groups and tags.
- **Use:** Mesh generation for FEA/CFD solvers
- **Tools:** Gmsh, FEniCS, GetDP
- **Spec:** [Gmsh File Formats](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format)

### **CGNS (.cgns)**
CFD General Notation System — a standard for CFD data including grids, flow solutions, and boundary conditions.
- **Use:** CFD mesh and solution exchange between solvers
- **Tools:** ANSYS, OpenFOAM (with converters), ParaView
- **Spec:** [CGNS Project](https://cgns.github.io/)

### **EXODUS II (.exo)**
A finite element data format built on NetCDF/HDF5, supporting unstructured meshes with time-dependent variables.
- **Use:** FEA results, multi-physics simulations
- **Tools:** ParaView, Cubit, SEACAS
- **Spec:** [Sandia SEACAS](https://sandialabs.github.io/seacas-docs/)

---

## **Solver Input and Output Formats**

### **OpenFOAM Case Directory**
A directory-based structure containing mesh (`polyMesh/`), initial/boundary conditions (`0/`), solver settings (`system/`), and time-step results. Uses ASCII dictionaries.
- **Use:** CFD simulations with OpenFOAM
- **Tools:** OpenFOAM, ParaView, pyFoam
- **Spec:** [OpenFOAM User Guide](https://www.openfoam.com/documentation/user-guide)

### **Abaqus Input File (.inp)**
A keyword-driven text file defining geometry, materials, boundary conditions, and analysis steps for the Abaqus FEA solver.
- **Use:** Structural, thermal, and multiphysics FEA
- **Tools:** Abaqus/CAE, CalculiX (partially compatible)

### **ANSYS Result Files (.rst, .rth)**
Binary result files from ANSYS Mechanical containing nodal/element solutions, reactions, and derived quantities.
- **Use:** Post-processing ANSYS simulations
- **Tools:** ANSYS Mechanical APDL, pyansys (Python)

---

## **Result and Visualization Formats**

### **VTK / VTU (.vtk, .vtu, .vtp, .vts)**
Visualization Toolkit formats for structured and unstructured mesh data with associated field variables. VTU (unstructured) is the most commonly used variant.
- **Use:** Simulation result visualization, data exchange between tools
- **Tools:** ParaView, VTK (Python), Mayavi, PyVista
- **Spec:** [VTK File Formats](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html)

### **HDF5 (.h5, .hdf5)**
Hierarchical Data Format version 5 — a flexible binary format for large, structured datasets with built-in compression and parallel I/O support.
- **Use:** Large simulation outputs, restart files, parameter studies
- **Tools:** h5py (Python), HDFView, ParaView
- **Spec:** [HDF Group](https://www.hdfgroup.org/solutions/hdf5/)

### **EnSight Gold (.case + .geo + .scl/.vec)**
A multi-file format for transient simulation results, widely supported by commercial and open-source tools.
- **Use:** CFD and FEA transient results
- **Tools:** ParaView, EnSight, ANSYS

### **CSV / TSV (.csv, .tsv)**
Plain-text tabular data for extracted quantities like time histories, convergence data, or probe measurements.
- **Use:** Monitoring data, post-processed scalar quantities, comparison with experiments
- **Tools:** Any spreadsheet or programming language

---

## **Workflow and Provenance Formats**

### **YAML / JSON (.yaml, .json)**
Human-readable configuration formats used for simulation parameters, workflow definitions, and metadata.
- **Use:** Solver configuration, parameter studies, workflow pipelines
- **Tools:** Any text editor, Python (PyYAML, json)

### **Provenance Metadata (W3C PROV)**
A set of standards for recording the provenance of data and processes, important for reproducibility.
- **Use:** Tracking simulation lineage and parameter history
- **Spec:** [W3C PROV](https://www.w3.org/TR/prov-overview/)

---

## **Format Selection Guide**

| Use Case | Recommended Format |
|----------|-------------------|
| CAD geometry exchange | STEP |
| Mesh generation | MSH (Gmsh) |
| CFD meshes and solutions | CGNS or OpenFOAM native |
| FEA results | VTU or EXODUS II |
| Large time-series data | HDF5 |
| Visualization | VTK/VTU |
| Configuration files | YAML |
| Extracted scalar data | CSV |
