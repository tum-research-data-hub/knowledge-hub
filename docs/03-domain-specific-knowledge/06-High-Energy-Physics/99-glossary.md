---
title: Glossary
slug: /domain-knowledge/high-energy-physics/glossary
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **High Energy Physics Glossary**

This glossary provides definitions for terms used elsewhere in this chapter, covering the event data model, the ROOT ecosystem, and derived analysis formats.

---

## **Event Data & Processing**

### **Event**
A single particle collision (or a fixed time window of one) and everything the detector recorded for it, including tracks, energy deposits, timing information, and reconstructed physics objects.

### **Reconstruction**
The processing stage that turns raw detector signals into physics objects (electrons, muons, jets, and similar) that an analysis can work with directly.

### **Derivation**
A further processing step that reduces reconstructed event data down to a smaller, analysis-focused format (e.g. ATLAS's DAOD derivations), trading generality for size and speed.

---

## **The ROOT Ecosystem**

### **ROOT**
A self-describing, hierarchical binary file format and the software framework built around it. Developed at CERN and used to store HEP data at essentially every processing stage.

### **TTree**
The columnar data structure historically used to store event data inside a ROOT file, where each variable is stored in its own branch and each entry corresponds to one event.

### **Branch**
A single named column within a TTree or RNTuple, holding the values of one variable across all stored events.

### **RNTuple**
ROOT's next-generation columnar event-data format, designed to succeed TTree with better throughput and compression for the High-Luminosity LHC era.

### **Ntuple**
A generic term for a flat, row-per-event data structure.

---

## **Derived Analysis Formats**

### **AOD (Analysis Object Data)**
Full reconstructed-event format, storing physics objects in experiment-specific C++ classes.

### **MiniAOD**
A reduced version of CMS's AOD, still based on CMS-specific classes but with less redundant content.

### **NanoAOD**
A flat, TTree-based CMS format using only standard C++ types, roughly 20 times smaller than MiniAOD and 200 times smaller than AOD.

### **DAOD_PHYSLITE**
An ATLAS derivation format from which the flat ROOT ntuples distributed through ATLAS Open Data are produced.

### **DecayTreeTuple**
The LHCb DaVinci tool used to write selected candidate and event information into a ROOT ntuple, the standard output of an LHCb analysis.

---

## **Software & Analysis Tools**

### **RDataFrame**
ROOT's modern, declarative interface for analyzing TTree and RNTuple data.

### **uproot**
A pure-Python library for reading and writing ROOT files without requiring a ROOT installation.

### **PyROOT**
The Python bindings for the ROOT framework, giving Python access to ROOT's full C++ API.

### **TBrowser**
ROOT's graphical browser for interactively inspecting the contents of a ROOT file, including its trees and branches.

### **CMSSW**
The CMS software framework used to process raw data into AOD/MiniAOD-level formats.

### **DaVinci**
The LHCb analysis software application used to select candidates and produce ROOT ntuples via tools like DecayTreeTuple.

---

## **Facilities & Data Sharing**

### **HL-LHC (High-Luminosity LHC)**
The upgraded phase of the LHC, expected to begin with Run 4, that will substantially increase collision rates and data volumes.

### **HEPData**
The open-access repository for the digitized data points, tables, and plots underlying HEP publications.

### **INSPIRE-HEP**
The literature database for high energy physics that HEPData records are anchored to, linking datasets to their originating publication.
