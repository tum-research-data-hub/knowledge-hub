---
title: Chemistry & Materials Data Types & Formats
slug: /domain-knowledge/chemistry-materials/data-types
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Chemistry & Materials Data Types & File Formats**

Chemistry and materials science use a wide range of data formats — from molecular structure representations and spectroscopic data to crystallographic files and simulation outputs. This page provides an overview of the most important formats you will encounter.

---

## **Molecular Structure Representations**

### **SMILES (Simplified Molecular-Input Line-Entry System)**
A line notation for encoding molecular structures as ASCII strings. Compact and widely supported by cheminformatics tools.
- **Example:** `c1ccccc1` (benzene), `CC(=O)O` (acetic acid)
- **Use:** Database queries, molecular property prediction, ML model input
- **Tools:** RDKit, Open Babel, CDK
- **Spec:** [OpenSMILES](http://opensmiles.org/)

### **InChI (International Chemical Identifier)**
A non-proprietary, textual identifier for chemical substances developed by IUPAC. Designed for unique identification rather than human readability.
- **Example:** `InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H`
- **Use:** Cross-database linking, deduplication, unique substance identification
- **Spec:** [IUPAC InChI](https://www.inchi-trust.org/)

### **MOL / SDF (.mol, .sdf)**
MDL Molfile and Structure-Data File — widely used formats for 2D and 3D molecular structures. SDF extends MOL to store multiple molecules with associated data fields.
- **Use:** Chemical databases, structure exchange, property storage
- **Tools:** RDKit, Open Babel, Marvin, Avogadro
- **Spec:** [CTfile Formats (Dassault)](https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf)

### **PDB (.pdb)**
Protein Data Bank format — stores 3D coordinates of atoms in macromolecules. Also used for small molecules in some workflows.
- **Use:** Protein structures, molecular docking, visualization
- **Tools:** PyMOL, VMD, Chimera, Biopython
- **Spec:** [wwPDB Format Guide](https://www.wwpdb.org/documentation/file-format)

### **XYZ (.xyz)**
A minimal plain-text format listing atom types and Cartesian coordinates. Simple and widely supported, but carries no bonding or metadata information.
- **Use:** Quick structure exchange, QM calculation input/output
- **Tools:** Avogadro, ASE (Python), Jmol

### **CIF (Crystallographic Information File, .cif)**
The standard format for crystallographic data, including unit cell parameters, symmetry, and atomic coordinates.
- **Use:** Crystal structures from X-ray or neutron diffraction
- **Tools:** Mercury (CCDC), VESTA, Olex2, pymatgen
- **Spec:** [IUCr CIF](https://www.iucr.org/resources/cif)

---

## **Spectroscopic Data**

### **JCAMP-DX (.jdx, .dx)**
Joint Committee on Atomic and Molecular Physical Data — Data Exchange format. A vendor-neutral text format for spectroscopic data.
- **Use:** IR, Raman, UV-Vis, NMR, mass spectra
- **Tools:** Various spectroscopy software, jcamp (Python)
- **Spec:** [IUPAC JCAMP-DX](https://iupac.org/what-we-do/digital-standards/jcamp-dx/)

### **mzML (.mzML)**
An open XML-based format for mass spectrometry data, replacing proprietary vendor formats.
- **Use:** LC-MS, GC-MS, metabolomics, proteomics
- **Tools:** OpenMS, MZmine, ProteoWizard
- **Spec:** [HUPO-PSI mzML](https://www.psidev.info/mzml)

### **NMReDATA (.nmredata.sdf)**
An extension of the SDF format to include NMR spectral assignments linked to molecular structures.
- **Use:** NMR data sharing with structure-spectrum assignments
- **Spec:** [NMReDATA Initiative](https://nmredata.org/)

### **Bruker / Varian / JEOL (vendor-specific)**
Proprietary binary formats from instrument manufacturers. Typically require vendor software or conversion tools.
- **Use:** Raw instrument data
- **Tools:** MestReNova, TopSpin, nmrglue (Python)

---

## **Simulation and Computation Data**

### **Gaussian Output (.log, .chk, .fchk)**
Output files from Gaussian quantum chemistry calculations, containing energies, geometries, molecular orbitals, and vibrational frequencies.
- **Use:** DFT and ab initio calculations
- **Tools:** GaussView, cclib (Python), Avogadro

### **VASP Files (POSCAR, CONTCAR, OUTCAR, vasprun.xml)**
A set of files used by the Vienna Ab initio Simulation Package for solid-state DFT calculations.
- **Use:** Periodic DFT calculations, materials properties
- **Tools:** pymatgen (Python), ASE, VESTA

### **Trajectory Files (.xtc, .trr, .dcd, .nc)**
Binary formats storing atomic positions (and optionally velocities/forces) at each time step of a molecular dynamics simulation.
- **Use:** Molecular dynamics analysis
- **Tools:** MDAnalysis (Python), GROMACS tools, VMD

### **NOMAD Archive Format**
A standardized format used by the NOMAD repository for computational materials science data.
- **Use:** Archiving and sharing simulation data
- **Spec:** [NOMAD](https://nomad-lab.eu/)

---

## **Reaction and Process Data**

### **RXN / RD Files (.rxn, .rdf)**
Extensions of the MOL/SDF format family for storing chemical reactions, including reactants, products, and agents.
- **Use:** Reaction databases, retrosynthesis
- **Tools:** RDKit, MarvinSketch

### **AnIML (Analytical Information Markup Language, .animl)**
An XML-based format for analytical chemistry data, supporting multiple technique types.
- **Use:** Multi-technique analytical data
- **Spec:** [ASTM E1947](https://www.animl.org/)

---

## **Format Selection Guide**

| Use Case | Recommended Format |
|----------|-------------------|
| Molecular structures (2D/3D) | MOL/SDF or SMILES |
| Unique substance identification | InChI |
| Crystal structures | CIF |
| Protein/macromolecule structures | PDB / mmCIF |
| Spectroscopic data (vendor-neutral) | JCAMP-DX |
| Mass spectrometry | mzML |
| QM calculation results | Gaussian/ORCA output, VASP files |
| MD trajectories | XTC, TRR, or DCD |
| Archiving computational data | NOMAD Archive |
