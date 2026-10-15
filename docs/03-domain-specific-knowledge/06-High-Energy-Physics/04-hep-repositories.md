---
title: High Energy Physics Repositories
slug: /domain-knowledge/high-energy-physics/repositories
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **High Energy Physics Repositories**

HEP research data is shared through two distinct kinds of repositories:
- publication-level results: the data points, tables, and cross-sections behind a paper's plots
- event-level open data: the underlying raw and derived collision datasets released by the collaborations themselves (see also the [format details](/domain-knowledge/high-energy-physics/data-types#open-data-portals--learning-resources) for these releases).

Both are tied together by INSPIRE-HEP, the community's literature and metadata hub, which links papers to their authors, experiments, and associated data records.

---

## **Publication-Level Data Repository**

| Repository | Description | Link |
|------------|-------------|------|
| **HEPData** | The primary open-access repository for digitized data points, tables, and plots underlying HEP publications, with records going back to the 1970s. Hosted on CERN infrastructure and maintained by Durham University, it is anchored to INSPIRE-HEP and arXiv identifiers, so datasets are directly linked to their originating paper. It does not host raw detector data, simulation samples, or analysis code. | [hepdata.net](https://www.hepdata.net/) |

---

## **Event-Level Open Data Repository**

| Repository | Description | Link |
|------------|-------------|------|
| **CERN Open Data Portal** | The central repository for released LHC datasets — raw, AOD, MiniAOD, and NanoAOD-level event data — and associated software, spanning ATLAS, CMS, LHCb, ALICE, and OPERA. | [opendata.cern.ch](https://opendata.cern.ch/) |

---

## **Literature & Metadata Hub**

| Repository | Description | Link |
|------------|-------------|------|
| **INSPIRE-HEP** | The community information platform for HEP scholarly content, run jointly by CERN, DESY, Fermilab, IHEP, IN2P3, and SLAC, and the successor of SPIRES, which served the field for almost 50 years. It consists of interlinked databases covering literature, authors, institutions, experiments, conferences, journals, jobs, and data, with citation tracking and BibTeX export. Researchers log in with ORCID to claim their papers, and a public REST API provides programmatic access. INSPIRE record IDs are the anchors HEPData uses to link datasets to their papers. It indexes publications and metadata rather than hosting research data itself. | [inspirehep.net](https://inspirehep.net/) |
