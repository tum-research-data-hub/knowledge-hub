---
title: Publishing Research Software
slug: /general-knowledge/publishing-research-software
sidebar_position: 2
---

:::caution Draft
This page is a first draft. If you notice errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Publishing Research Software**

Research software is any code written in the course of research: analysis scripts, notebooks, workflows, simulation codes, libraries and tools. "Publishing" it means more than making a repository visible. Published software has a **fixed version**, a **licence**, a **persistent identifier**, and enough **documentation** for someone else to run it.

## **1\. Why Publish Research Software?**

* **Results cannot be verified without it.** A methods section describes what the code was meant to do; only the code shows what was actually done. Data without the code that processed it is rarely reproducible.
* **It is part of good research practice.** The [DFG Code of Conduct](https://wissenschaftliche-integritaet.de/en/) counts software among the research results that should be made accessible and archived (guidelines 13 and 17), and funders and journals increasingly ask for a code availability statement.
* **It avoids duplicated work.** Published code can be reused and improved instead of being rewritten by the next group, or by the next doctoral researcher in your own group.
* **It preserves knowledge.** Code that lives only on a personal laptop or account is lost when its author leaves.

The [FAIR Principles for Research Software (FAIR4RS)](https://doi.org/10.15497/RDA00068) translate the FAIR principles to software and are the common reference for what "well published" means.

## **2\. How Authors Benefit**

* **Citable credit.** A software release with a DOI can be cited like a paper and listed in your publication record and ORCID profile. Research assessment is moving towards recognising software and data as outputs in their own right.
* **Visibility and collaboration.** Users report bugs, contribute improvements, and become collaborators and co-authors.
* **Better code.** Code that is written to be read by others tends to be cleaner and better documented, and errors are found earlier.
* **Less support effort.** A documented, installable release answers many of the e-mails you would otherwise receive.
* **Your own future use.** The most frequent reuser of your code is you, two years later.
* **An additional publication.** Substantial software can be described in a software paper (see below).

## **3\. How Research Software Is Usually Published**

What you publish depends on what kind of software it is.

| Situation | Usual way of publishing |
|---|---|
| **Analysis code for one paper** (scripts, notebooks) | One repository per paper. Tag the exact version used for the submitted results, archive that release with a DOI, and cite the DOI in the paper's code availability statement. |
| **Reusable tool or library** | Its own repository with versioned releases, each archived with a DOI. Distributed through the package index of its language (PyPI, CRAN, Bioconductor, conda-forge, Bioconda). Often accompanied by a software paper. |
| **Workflow** (Snakemake, Nextflow, Galaxy) | Repository plus registration in [WorkflowHub](https://workflowhub.eu/), which assigns DOIs and records the workflow type; see [Workflow Management Systems](/general-knowledge/wms). |
| **Large simulation or community code** | Project website and repository with a governance model, releases and a reference publication. |

### **With or Without Data?**

Publish code and data **separately, and link them**:

* **Code** goes into a code repository and, for each release, into an archive that issues a DOI.
* **Data** goes into a data repository, preferably a discipline-specific one; see [Repository Collections](/general-knowledge/repo).
* **Link them in both directions**: the dataset record names the code DOI ("is compiled by"), the code record names the dataset DOI ("requires" or "is supplement to"), and the paper cites both.

Keeping them apart has practical reasons: Git is not made for large files, code and data usually need different licences, and data may have access restrictions that the code does not. Small example or test datasets are the exception; they belong with the code, so that others can try it immediately.

Do not publish code only as a supplementary file of the article. Supplements are hard to find, cannot be versioned, and carry no licence.

### **The Usual Steps**

1. **Develop under version control** from the start ([Git & Version Control](/general-knowledge/git)).
2. **Add the essential files**: a `README` (what it does, how to install and run it, an example), a `LICENSE`, and a [`CITATION.cff`](https://citation-file-format.github.io/) file that states how the software should be cited. GitHub and Zenodo read this file automatically.
3. **Record the dependencies** with exact versions; see [Environments and Dependency Management](/general-knowledge/containers#first-step-environments-and-dependency-management).
4. **Check what must not be published**: passwords and access tokens (also in the Git history), personal data, internal paths, and third-party code or data whose licence does not allow redistribution.
5. **Create a release** with a version number (e.g. `v1.0.0`).
6. **Archive the release and obtain a DOI.** The most common route is the [Zenodo integration for GitHub](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content): once it is switched on for a repository, every release is archived on Zenodo automatically. Zenodo issues one DOI per version and one "concept DOI" for the software as a whole. Cite the version DOI for the exact code behind a result.
7. **Cite it** in the paper, including the version, and add the software to your ORCID record.

Optionally, describe substantial software in a **software paper**, for example in the [Journal of Open Source Software (JOSS)](https://joss.theoj.org/), SoftwareX, or a journal of your discipline, and register it in a discipline catalogue such as [bio.tools](https://bio.tools/).

## **4\. Licences**

Without a licence, nobody is legally allowed to use, modify or redistribute your code, even if it is publicly visible. Choose a licence when you create the repository.

| Licence | Type | In short |
|---|---|---|
| **MIT**, **BSD 3-Clause** | Permissive | Anyone may use, modify and redistribute, also in closed-source products, if the copyright notice is kept. The most common choice for research software. |
| **Apache 2.0** | Permissive | As above, with an explicit patent grant and clearer rules for contributions. A good default for larger projects. |
| **GPL v3** | Copyleft | Derived works that are distributed must be published under the same licence. Ensures that improvements stay open. |
| **LGPL v3** | Weak copyleft | For libraries: changes to the library stay open, but software that merely uses it can have any licence. |
| **AGPL v3** | Strong copyleft | Like GPL, and also applies when the software is offered as a web service. |

Points to keep in mind:

* **Do not use Creative Commons licences for code.** They are meant for text, images and data. Use them for the documentation or the example data in the repository if you wish.
* **Your licence must be compatible with your dependencies.** Including GPL-licensed code generally means your combined work must be GPL as well.
* **Check that you are allowed to decide.** The economic rights to software written as part of employment generally lie with the employer, and project contracts with industry partners or funders can contain their own rules. If the software has commercial potential or was developed in a cooperation, talk to TUM ForTe before publishing.
* **Make the licence machine-readable**: a `LICENSE` file with the unmodified licence text in the top folder of the repository. The [REUSE](https://reuse.software/) specification shows how to label every file.

[choosealicense.com](https://choosealicense.com/) helps with the choice. The TUM handout [Information on Licenses for Published Scientific Data and Software Programs (PDF)](https://mediatum.ub.tum.de/doc/1711141/1711141.pdf) explains the options in detail.

## **5\. Where to Publish**

Development platforms and archives do different jobs, and you need both.

| Platform | Role | Notes |
|---|---|---|
| [**GitHub**](https://github.com) | Development and visibility | Where most open-source research software lives and is found. Not an archive: repositories can be renamed, made private or deleted at any time, and there is no DOI. |
| [**LRZ GitLab**](https://gitlab.lrz.de) | Development | Hosted in Munich, with TUM login. Well suited for internal and unpublished work. Projects can be made public, but external contributors cannot simply sign up. |
| [**Zenodo**](https://zenodo.org/) | Archive with DOI | Long-term archive operated by CERN; versioned DOIs; integrates with GitHub. The standard choice for archiving releases. |
| [**Software Heritage**](https://www.softwareheritage.org/) | Source code archive | Archives public repositories including their full history and provides identifiers (SWHIDs) down to a single line of code. You can trigger archiving of your repository with "Save code now". |
| **Package indexes** (PyPI, CRAN, Bioconductor, conda-forge, Bioconda) | Distribution | Make installation easy for users. They complement an archive and do not replace it. |
| [**WorkflowHub**](https://workflowhub.eu/) | Registry for workflows | DOIs and rich metadata for workflows. |

### **Is a Private GitHub Repository Enough?**

For **development**: yes, a private repository on GitHub or, preferably, on LRZ GitLab is a normal place for work in progress.

As a **publication**: no. A private repository is not accessible, not citable and not archived. "Code available from the authors on request" does not meet the expectations of most journals and funders, and in practice such requests often go unanswered once people have moved on. The same applies to a public repository without a licence and without an archived release.

If you want to keep code closed until the paper is accepted, there are better options than waiting: archive the release on Zenodo with restricted access or an embargo, which reserves the DOI for the manuscript and lets you share an access link with reviewers, then open the record on publication. Remember that making a repository public later exposes its entire history; clean it up beforehand.

### **Should a Research Group Have Its Own GitHub Organisation?**

Yes, in most cases. A GitHub **organisation** (or a **group** on LRZ GitLab) is a shared space that belongs to the research group, not to an individual.

* **Continuity.** Repositories in personal accounts leave with the person. In an organisation they stay with the group, and repositories can be transferred into it without breaking links.
* **Shared responsibility.** Several people can administer it. Appoint at least two owners, one of them permanent staff, e.g. the group leader.
* **Visibility.** All software of the group can be found in one place and linked from the group's website.
* **Access management.** Teams and roles control who can read or change which repository, including external collaborators.

Organisations are free of charge for public repositories. Everybody joins with their personal account: do not share one login among several people, which is against the platform's terms and makes it impossible to see who changed what. Agree on a few conventions (naming, a README and licence for every repository, what happens when someone leaves), and keep code that must stay internal on LRZ GitLab.

## **6\. Minimum Checklist**

- [ ] The code is under version control and the published version is tagged as a release
- [ ] `README` with purpose, installation, usage and a small example
- [ ] `LICENSE` file with a standard open-source licence, and the right to publish is clarified
- [ ] `CITATION.cff` with authors (and ORCID iDs), title and version
- [ ] Dependencies are listed with versions (lock file, environment file or container)
- [ ] No passwords, tokens or personal data in the files or the history
- [ ] The release is archived with a DOI (e.g. Zenodo)
- [ ] Code, data and paper reference each other by their identifiers
- [ ] A named contact or issue tracker for questions

A few [tests](/general-knowledge/testing-ci) and an example dataset make the difference between code that is available and code that others can actually use.

## **Further Information**

* [TUM Research Data Hub: Publish & Share Data](https://web.tum.de/en/researchdata/support-information/publish-share-data/)
* Chue Hong N. P. et al. (2022): [FAIR Principles for Research Software (FAIR4RS Principles)](https://doi.org/10.15497/RDA00068)
* [The Turing Way: Guide for Reproducible Research](https://book.the-turing-way.org/reproducible-research/reproducible-research/)
* [Citation File Format (CFF)](https://citation-file-format.github.io/) and [CodeMeta](https://codemeta.github.io/) for software metadata
* [de-RSE](https://de-rse.org/en/), the society for research software engineering in Germany
