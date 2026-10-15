---
title: High Energy Physics Data Types & Formats
slug: /domain-knowledge/high-energy-physics/data-types
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **High Energy Physics Data Types & File Formats**

High energy physics (HEP) experiments study the fundamental constituents of matter and their interactions by accelerating particles to close to the speed of light and colliding them. At the LHC, for example, protons or heavy ions are smashed together tens of millions of times per second inside detectors like ATLAS, CMS, LHCb, and ALICE. Each collision can produce a shower of new particles, and it is the properties of these particles (momentum, energy, trajectory, decay products) that experiments reconstruct from detector responses.

Most HEP data is **event-based**: detectors record discrete "events," each corresponding to a single collision and the resulting particles passing through the detector. An event typically holds the tracks, energy deposits, timing information, and reconstructed physics objects (electrons, muons, jets, etc.) observed for that collision, and analyses proceed by selecting, transforming, and aggregating over large collections of independent events. Other experiments like Belle II, for instance, perform candidate-based analysis, where the primary unit is a reconstructed decay candidate (which may be one of several found within a single event) rather than the event itself, and store their data accordingly. This page focuses on event-based formats, centered on the ROOT ecosystem that dominates HEP software.

---

## **The ROOT Ecosystem**

### **ROOT (.root)**
A self-describing, hierarchical binary file format developed at CERN, and the de facto standard for storing HEP data at every processing stage — raw, reconstructed, and analysis-ready. A ROOT file behaves like a directory tree that can hold histograms, trees, and arbitrary C++ objects, which the ROOT framework can serialize and deserialize automatically.
- **Use:** Raw, reconstructed, and analysis-level data across ATLAS, CMS, LHCb, ALICE, and most other HEP and nuclear physics experiments
- **Tools:** ROOT (C++/Python), PyROOT, RDataFrame, uproot (pure-Python reader, no ROOT installation required)
- **Spec:** [ROOT File Format](https://root.cern/manual/root_files/)

### **TTree — the classic event container**
The columnar data structure historically used to store event data inside a ROOT file. Each variable is stored in its own "branch," so an analysis can read only the columns it needs rather than loading whole events, and each entry in the tree corresponds to one event. TTree underlies more than two exabytes of LHC Run 1–3 data as well as the CERN Open Data releases.
- **Use:** Primary event-data container in AOD/MiniAOD/NanoAOD-style formats and virtually all pre-Run 4 analysis ntuples
- **Tools:** ROOT `TTree`/`TChain` API, RDataFrame, uproot
- **Note:** The format is over 25 years old; its I/O throughput, limited multicore scalability, and inability to fully exploit modern storage (fast NVMe drives, object stores) become a bottleneck at HL-LHC data rates

### **RNTuple — the HL-LHC successor**
ROOT's next-generation columnar event-data format, a redesign of both the on-disk layout and the C++ API that is not compatible with TTree (RNTuple itself guarantees forward and backward compatibility across ROOT versions). It targets the throughput and storage requirements of the High-Luminosity LHC (Run 4 and beyond), with read throughput several times higher than TTree, better write performance and multicore scaling, and file sizes typically 10–50% smaller. RNTuple data is still stored inside a standard `.root` file.
- **Use:** Intended to replace TTree as the primary storage format for raw and derived event data starting with Run 4; the on-disk format was finalized in ROOT 6.34 (2024) and the core API became production-ready in ROOT 6.36 (2025), with initial support already in ATLAS, CMS, and LHCb software
- **Tools:** ROOT (`RNTupleReader`/`RNTupleWriter`), RDataFrame, uproot (reading supported in uproot 5)
- **Spec:** [RNTuple Binary Format Specification](https://root.cern.ch/doc/v636/md_tree_2ntuple_2doc_2BinaryFormatSpecification.html), [First release of the RNTuple on-disk format](https://root.cern.ch/blog/rntuple-binary-format/), [RNTuple: where are we now?](https://root.cern/blog/rntuple-update/)

---

## **Derived and Analysis-Ready Ntuples**

Raw and fully reconstructed event data are too large and too specific to analyze directly, so experiments derive progressively slimmer, flatter ROOT-based formats for physics analysis and public release.

### **CMS AOD / MiniAOD / NanoAOD**
A tiered reduction chain used by CMS. AOD contains full reconstructed physics objects in CMS-specific C++ classes; MiniAOD trims this down while remaining class-based; NanoAOD flattens everything into a plain TTree of standard C++ types (floats, ints, and arrays), making it roughly 20 times smaller than MiniAOD and about 200 times smaller than AOD.
- **Use:** CMS physics analyses and the bulk of CMS releases on the CERN Open Data Portal
- **Tools:** ROOT, uproot, RDataFrame for NanoAOD; CMSSW for AOD/MiniAOD

### **ATLAS DAOD_PHYSLITE ntuples**
ATLAS Open Data distributes flat ROOT ntuples (a single TTree with several dozen branches, e.g. lepton `pt`, `eta`, `phi`, `charge`, and energy per event) derived from the experiment's DAOD_PHYSLITE derivation format.
- **Use:** ATLAS Open Data educational and outreach analyses
- **Tools:** ROOT, uproot

### **LHCb DecayTreeTuple ntuples**
LHCb's DaVinci analysis application writes selected candidate and event information into a ROOT ntuple using the `DecayTreeTuple` tool, the standard end product of an LHCb analysis and ready for inspection with `TBrowser` or further processing.
- **Use:** LHCb physics analyses, including the workflow taught in the LHCb Starterkit
- **Tools:** ROOT (`TTree`, `TBrowser`), DaVinci

---

## **Open Data Portals & Learning Resources**

- **[CERN Open Data Portal](https://opendata.cern.ch/)** — the central repository for released LHC datasets (raw, AOD, MiniAOD, NanoAOD) and associated software, spanning ATLAS, CMS, LHCb, ALICE, and OPERA
- **[ATLAS Open Data](https://opendata.atlas.cern/)** — ATLAS-specific datasets and documentation, including a dedicated guide to its [ROOT ntuple format](https://opendata.atlas.cern/docs/documentation/data_format/ntuple)
- **[LHCb Starterkit (Learning Kit)](https://lhcb.github.io/starterkit-lessons/)** — hands-on lessons that walk through the LHCb software stack from raw data to a finished ROOT ntuple, including how to inspect it with `TBrowser`
