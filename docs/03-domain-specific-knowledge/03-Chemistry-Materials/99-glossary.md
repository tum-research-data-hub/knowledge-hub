---
title: Glossary
slug: /domain-knowledge/chemistry-materials/glossary
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Chemistry & Materials Science Glossary**

This glossary provides definitions for common terms used in computational chemistry, materials science, and chemical data management.

---

## **Molecular Representation**

### **SMILES**
Simplified Molecular-Input Line-Entry System — a compact, ASCII-based notation for encoding molecular structures as text strings. Example: `CC(=O)O` represents acetic acid.

### **InChI**
International Chemical Identifier — an IUPAC standard for uniquely identifying chemical substances using a layered text string. Unlike SMILES, InChI is canonical (one molecule = one InChI).

### **InChIKey**
A fixed-length (27-character) hash of an InChI string, designed for database searching and indexing.

### **CAS Number**
A unique numerical identifier assigned by the Chemical Abstracts Service to every described chemical substance. Example: 7732-18-5 (water).

### **Molecular Formula**
A notation indicating the types and numbers of atoms in a molecule (e.g., C₆H₁₂O₆ for glucose).

### **Force Field**
A mathematical model describing the potential energy of a molecular system as a function of atomic positions. Used in molecular dynamics and energy minimization.

---

## **Computational Chemistry**

### **DFT (Density Functional Theory)**
A quantum mechanical method that calculates electronic structure by modeling electron density rather than individual electron wave functions. Widely used for molecules and materials.

### **Ab Initio**
Computational methods based on fundamental quantum mechanics without empirical parameters. More accurate but more expensive than semi-empirical methods.

### **Molecular Dynamics (MD)**
A simulation method that computes the trajectory of atoms over time by numerically solving Newton's equations of motion.

### **Basis Set**
A set of mathematical functions used to represent electron wave functions in quantum chemistry calculations. Larger basis sets are more accurate but more expensive. Examples: 6-31G*, cc-pVTZ.

### **Functional (DFT)**
The mathematical approximation used in DFT to describe electron exchange and correlation. Examples: B3LYP, PBE, M06-2X.

### **Geometry Optimization**
The process of finding the atomic arrangement that minimizes the total energy of a molecule or material.

### **Transition State**
A first-order saddle point on the potential energy surface, representing the highest-energy configuration along a reaction pathway.

### **Potential Energy Surface (PES)**
A mathematical surface describing the energy of a system as a function of its atomic coordinates.

---

## **Materials Science**

### **Crystal Structure**
The periodic arrangement of atoms in a crystalline solid, described by a unit cell and its symmetry.

### **Unit Cell**
The smallest repeating unit of a crystal that, when translated in all three dimensions, reproduces the entire crystal structure.

### **Lattice Parameters**
The lengths (a, b, c) and angles (α, β, γ) that define the shape and size of a unit cell.

### **Band Gap**
The energy difference between the top of the valence band and the bottom of the conduction band in a material. Determines whether a material is a conductor, semiconductor, or insulator.

### **Phase Diagram**
A graphical representation of the thermodynamic conditions (temperature, pressure, composition) under which different phases of a material exist.

### **Defect**
An irregularity in the crystal lattice, such as a vacancy (missing atom), interstitial (extra atom), or substitutional impurity.

---

## **Spectroscopy**

### **NMR (Nuclear Magnetic Resonance)**
A spectroscopic technique that exploits the magnetic properties of atomic nuclei to determine molecular structure, dynamics, and environment.

### **Chemical Shift**
The resonance frequency of a nucleus relative to a standard reference in NMR, expressed in parts per million (ppm). Reflects the electronic environment of the atom.

### **Mass Spectrometry (MS)**
An analytical technique that measures the mass-to-charge ratio (m/z) of ions, used for molecular weight determination and structural elucidation.

### **IR Spectroscopy (Infrared)**
A technique measuring the absorption of infrared light by molecular vibrations, used to identify functional groups.

### **Raman Spectroscopy**
A technique based on inelastic scattering of light by molecular vibrations, complementary to IR spectroscopy.

### **X-ray Diffraction (XRD)**
A technique that determines crystal structure by measuring the diffraction pattern produced when X-rays interact with a crystalline sample.

---

## **Data and File Formats**

### **CIF (Crystallographic Information File)**
The standard file format for crystallographic data, including unit cell parameters, symmetry, and atomic coordinates.

### **MOL / SDF**
MDL Molfile and Structure-Data File — widely used formats for 2D and 3D molecular structures.

### **JCAMP-DX**
A vendor-neutral text format for spectroscopic data exchange.

### **Trajectory File**
A file containing atomic positions (and optionally velocities) at successive time steps of an MD simulation. Common formats: XTC, TRR, DCD.

### **mzML**
An open XML-based format for mass spectrometry data, replacing proprietary vendor formats.

---

## **Cheminformatics**

### **Virtual Screening**
Computational evaluation of large compound libraries against a target (protein, material property) to identify promising candidates.

### **QSAR/QSPR**
Quantitative Structure-Activity/Property Relationship — statistical models relating molecular structure to biological activity or physical properties.

### **Molecular Fingerprint**
A binary or count vector encoding structural features of a molecule, used for similarity searching and machine learning.

### **Tanimoto Coefficient**
A similarity metric between two molecular fingerprints, ranging from 0 (no overlap) to 1 (identical).

---

**Need clarification on a term?** Contact us at [rdhub@mdsi.tum.de](mailto:rdhub@mdsi.tum.de) or explore our [Use Cases](/domain-knowledge/chemistry-materials/use-cases) for practical examples.
