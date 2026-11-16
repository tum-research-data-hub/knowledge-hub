---
title: Containers
slug: /general-knowledge/containers
---

# **Containers**

## **What Are Containers?**

A **container** packages software and all its dependencies into a single, portable unit that runs consistently across different computing environments. Unlike virtual machines, containers share the host operating system's kernel, making them lightweight and fast to start.

**Why containers matter for research:**
- **Reproducibility** — Your analysis runs the same way on your laptop, a colleague's machine, or an HPC cluster
- **Dependency isolation** — No more "it works on my machine" problems
- **Portability** — Move workflows between local systems, cloud, and HPC without modification
- **Version control** — Pin exact software versions for long-term reproducibility

---

## **Container Systems You Should Know**

### **Docker**
The most popular container platform; the de facto standard for development and cloud environments.
- Excellent documentation and ecosystem
- Requires root privileges (problematic on shared HPC systems)
- Best for: local development, cloud deployment, CI/CD pipelines

### **Enroot**
[Enroot](https://doku.lrz.de/4-enroot-2750921686.html) is NVIDIA's container runtime designed for HPC environments with GPU workloads.
- Lightweight, unprivileged container solution
- Native GPU support for deep learning and CUDA applications
- Can import Docker and Singularity images
- **Available on LRZ systems for GPU-accelerated workloads**

### **Singularity / Apptainer**
[Apptainer](https://doku.lrz.de/apptainer-singularity-at-lrz-400558542.html) (formerly Singularity) is designed specifically for HPC and scientific computing.
- Runs without root privileges
- Directly executes on HPC job schedulers (SLURM)
- Can convert Docker images to Apptainer format
- **This is the primary container runtime on LRZ HPC systems**

### **Kubernetes** (Orchestration)
Not a container runtime itself, but a system for managing containers at scale. Relevant if you're deploying services, less so for typical research workflows.

---

## **Containers at LRZ**

**Apptainer is the supported container runtime on LRZ HPC clusters** (SuperMUC-NG, Linux Cluster).

| Task | Command |
|------|---------|
| Pull a Docker image | `apptainer pull docker://python:3.11` |
| Run a container | `apptainer exec python_3.11.sif python script.py` |
| Build from definition | `apptainer build myimage.sif mydef.def` |

**LRZ Documentation:**
- [Apptainer on LRZ Systems](https://doku.lrz.de/apptainer-singularity-at-lrz-400558542.html)
- [Enroot on LRZ Systems](https://doku.lrz.de/enroot-10746586.html)
- [Container Best Practices](https://doku.lrz.de/container-best-practices-116098071.html)

---

## **Quick Decision Guide**

| Environment | Recommended Container System |
|-------------|------------------------------|
| Local development | Docker |
| LRZ HPC (SuperMUC-NG, Linux Cluster) | **Apptainer** |
| LRZ HPC with GPUs | **Enroot** or Apptainer |
| Cloud (AWS, GCP, Azure) | Docker |
| CI/CD pipelines | Docker |

:::tip For Reproducible Research
Build your container image with a **definition file** (Dockerfile, Apptainer .def) and store it in version control alongside your analysis code. This ensures anyone can rebuild your exact environment years later.
:::
