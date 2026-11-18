---
title: Git & Version Control
slug: /general-knowledge/git
sidebar_position: 2
---

# **Git & Version Control**

## **What is Version Control?**

A **Version Control System (VCS)** tracks changes to files over time, allowing you to recall specific versions, compare changes, and collaborate without overwriting each other's work. Think of it as an "unlimited undo" combined with a complete audit trail of who changed what and when.

## **What is Git?**

**Git** is the most widely used distributed VCS, created by Linus Torvalds in 2005 for Linux kernel development. Unlike older centralized systems (CVS, SVN), every Git user has a complete copy of the repository history, enabling offline work and robust backup.

Git excels at:
- **Tracking changes** to code, scripts, manuscripts, and configuration files
- **Branching and merging** for parallel development
- **Collaboration** across teams and institutions
- **Reproducibility** by linking analysis code to specific data versions

---

## **GitHub vs. GitLab**

Both are web platforms built around Git, but they serve different needs.

**GitHub** is the largest public code hosting platform, home to most open-source projects. It's ideal for public collaboration, sharing code with the community, and building a visible portfolio. However, free accounts have limited private repository features, and data is stored on US servers.

**GitLab** can be self-hosted, which is why many institutions run their own instances. **LRZ GitLab** ([gitlab.lrz.de](https://gitlab.lrz.de)) is the recommended choice for TUM researchers working on non-public projects — your data stays within the European research infrastructure, and you get fine-grained access control for team collaboration.

**In practice:** Use GitHub for open-source and public projects; use LRZ GitLab for internal research, unpublished work, or anything requiring institutional data governance.

**Key links:**
- **GitHub:** [github.com](https://github.com)
- **LRZ GitLab:** [gitlab.lrz.de](https://gitlab.lrz.de) — Use your LRZ/TUM credentials; recommended for non-public research projects

---

## **Additional Git functionalities: Git Pages & Wikis**

### **Git Pages** (GitHub Pages / GitLab Pages)
Host static websites directly from your repository — perfect for:
- Project documentation
- Lab websites
- Research data landing pages
- This Knowledge Hub itself!

### **Git Wikis**
Every repository can have an associated wiki for:
- Extended documentation
- Meeting notes
- Onboarding guides for new team members

Both are version-controlled and free with your repository.

---

## **Should You Use Git?**

**Git makes sense for:**
- Writing code, scripts, or analysis pipelines
- Collaborative writing in text-based formats (LaTeX, Markdown, R Markdown)
- Configuration files, workflow definitions, and documentation
- Any project where you need to track changes over time or collaborate with others

**Git requires extra consideration for:**
- **Large binary files** (images, videos, compressed archives) — Git stores every version, so repositories grow quickly. Consider Git LFS or external storage for large assets.
- **Sensitive or unpublished data** — Never put research data on public GitHub. Use LRZ GitLab or institutional storage with appropriate access controls.

**Git may not be the right tool for:**
- **Purely experimental lab work** with no computational component — a well-organized folder structure and regular backups may be sufficient.
- **Non-technical collaborators** who work exclusively in Word/Excel and are unfamiliar with version control — the learning curve may outweigh the benefits for short-term collaborations. (Though for long-term projects, it's often worth the investment.)
- **Very large datasets** — Git is designed for code, not data. Use dedicated research data repositories (Zenodo, institutional storage) and reference data via DOIs instead.

:::tip Rule of Thumb
If your work involves text-based files that change over time and you want to track who changed what and when — Git is likely worth learning.
:::

---

## **Getting Started**

This is not a tutorial — excellent resources exist:
- [Pro Git Book](https://git-scm.com/book/en/v2) (free, comprehensive)
- [GitHub Skills](https://skills.github.com/) (interactive)
- [LRZ GitLab Documentation](https://doku.lrz.de/display/PUBLIC/GitLab)

For research data management, the key insight is: **Git tracks your analysis code; your data repository tracks your data.** Link them via DOIs and version tags for full reproducibility. How to turn a repository into a citable publication is described under [Publishing Research Software](/general-knowledge/publishing-research-software).
