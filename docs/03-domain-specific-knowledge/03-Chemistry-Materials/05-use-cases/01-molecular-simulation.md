---
title: Molecular Dynamics Simulation
slug: /domain-knowledge/chemistry-materials/use-cases/molecular-simulation
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Use Case 1: Molecular Dynamics Simulation**

Molecular dynamics (MD) simulations model the physical movements of atoms and molecules over time by numerically solving Newton's equations of motion. These simulations generate trajectory data — snapshots of atomic positions at each time step — that can reveal protein folding, material properties, or chemical reaction mechanisms.

## **Workflow Overview**

```
System setup → Energy minimization → Equilibration → Production run → Analysis
Structure files (PDB/GRO) + Force field → Trajectory (XTC/TRR) → Properties & statistics
```

## **Key Concepts**

### **Force Fields**

A force field defines how atoms interact through mathematical functions describing:
- **Bonded interactions** — Bond stretching, angle bending, dihedral torsions
- **Non-bonded interactions** — Van der Waals (Lennard-Jones) and electrostatic (Coulomb) forces

Common force fields:
- **AMBER** — Widely used for biomolecules
- **CHARMM** — Proteins, lipids, nucleic acids
- **OPLS** — Organic molecules and liquids
- **ReaxFF** — Reactive force field for bond breaking/forming

### **Trajectory Data**

An MD simulation produces a trajectory file containing atomic coordinates (and optionally velocities and forces) at regular intervals. A typical trajectory:
- Time step: 1-2 femtoseconds
- Save interval: every 1-10 picoseconds
- Total duration: nanoseconds to microseconds
- File size: gigabytes to terabytes

### **Analysis Types**

- **RMSD / RMSF** — Structural deviation and flexibility
- **Radial distribution function (RDF)** — Packing and structure in liquids/materials
- **Hydrogen bond analysis** — Interaction networks
- **Free energy calculations** — Binding affinities, solvation energies
- **Diffusion coefficients** — Mean square displacement analysis

## **Popular Tools**

- **GROMACS** — High-performance MD engine, excellent for biomolecules
- **LAMMPS** — Flexible MD code for materials science
- **NAMD** — Parallel MD for large biomolecular systems
- **OpenMM** — GPU-accelerated MD with Python API
- **MDAnalysis** (Python) — Trajectory analysis library
- **VMD** — Molecular visualization and analysis

## **Code Example: GROMACS Workflow**

<details>
<summary>**Click to expand Bash workflow**</summary>

```bash
#!/bin/bash
# GROMACS molecular dynamics simulation workflow

# Step 1: Prepare topology from PDB structure
gmx pdb2gmx -f protein.pdb -o processed.gro -water tip3p -ff amber99sb

# Step 2: Define simulation box and add solvent
gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

# Step 3: Add ions to neutralize the system
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
echo "SOL" | gmx genion -s ions.tpr -o neutral.gro -p topol.top -pname NA -nname CL -neutral

# Step 4: Energy minimization
gmx grompp -f minim.mdp -c neutral.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# Step 5: NVT equilibration (100 ps)
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# Step 6: NPT equilibration (100 ps)
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# Step 7: Production MD (100 ns)
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md

# Step 8: Basic analysis
echo "Backbone Backbone" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
echo "Backbone" | gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg

echo "Simulation and basic analysis complete!"
```

</details>

## **Code Example: Trajectory Analysis with Python**

<details>
<summary>**Click to expand Python analysis script**</summary>

```python
#!/usr/bin/env python3
"""Analyze MD trajectory using MDAnalysis."""

import MDAnalysis as mda
from MDAnalysis.analysis import rms, diffusionmap
import matplotlib.pyplot as plt
import numpy as np

# Load trajectory
u = mda.Universe("md.tpr", "md.xtc")

# RMSD analysis
rmsd_analysis = rms.RMSD(u, select="backbone")
rmsd_analysis.run()

# Plot RMSD over time
time_ns = rmsd_analysis.results.rmsd[:, 1] / 1000  # ps to ns
rmsd_values = rmsd_analysis.results.rmsd[:, 2]      # RMSD in Angstrom

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(time_ns, rmsd_values, linewidth=0.5)
ax.set_xlabel("Time (ns)")
ax.set_ylabel("RMSD (Å)")
ax.set_title("Backbone RMSD")
plt.tight_layout()
plt.savefig("rmsd.pdf", dpi=300)
plt.savefig("rmsd.png", dpi=300)

# RMSF per residue
rmsf_analysis = rms.RMSF(u.select_atoms("name CA"))
rmsf_analysis.run()

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(rmsf_analysis.results.rmsf)
ax.set_xlabel("Residue Index")
ax.set_ylabel("RMSF (Å)")
ax.set_title("Per-Residue Flexibility (Cα RMSF)")
plt.tight_layout()
plt.savefig("rmsf.pdf", dpi=300)
plt.savefig("rmsf.png", dpi=300)

print("Analysis complete!")
```

</details>

## **Expected Outputs**

- **Trajectory files** — XTC/TRR (coordinates over time, typically the largest output)
- **Energy files** — EDR (thermodynamic properties over time)
- **Checkpoint files** — CPT (for restarting simulations)
- **Analysis plots** — RMSD, RMSF, RDF, energy profiles (PDF + PNG)
- **Extracted properties** — CSV/XVG files with time series data

## **Computational Requirements**

| System Size | Atoms | CPU Cores | GPU | RAM | Time (100 ns) |
|-------------|-------|-----------|-----|-----|----------------|
| Small protein in water | ~50,000 | 8-16 | 1 | 4-8 GB | 1-3 days |
| Large protein complex | ~500,000 | 32-64 | 1-2 | 16-32 GB | 1-2 weeks |
| Membrane system | ~200,000 | 16-32 | 1 | 8-16 GB | 3-7 days |

## **Common Issues & Troubleshooting**

:::warning Common Problems

**Simulation crashes (LINCS/SHAKE warnings)**
- Usually indicates bad initial geometry — ensure proper energy minimization
- Reduce time step (try 1 fs instead of 2 fs)
- Check for steric clashes in the starting structure

**System not equilibrated**
- Monitor temperature, pressure, and energy during equilibration
- Extend equilibration if properties haven't stabilized

**Trajectory too large**
- Save frames less frequently (every 10 ps instead of every 1 ps)
- Use compressed trajectory format (XTC over TRR)
- Store only atoms of interest using output groups

:::

## **Key Considerations**

- **Document your force field choice** and any custom parameters
- **Always equilibrate** before production runs — check that temperature, pressure, and density are stable
- **Store input files alongside results** — MDP files, topology, and starting structure are essential for reproducibility
- **Use version control** for analysis scripts and parameter files
