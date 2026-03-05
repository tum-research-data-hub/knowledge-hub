---
title: Computational Fluid Dynamics (CFD) Simulation
slug: /domain-knowledge/simulation-modeling/use-cases/cfd-simulation
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Use Case 1: Computational Fluid Dynamics (CFD) Simulation**

CFD simulations solve the Navier-Stokes equations numerically to predict fluid flow, heat transfer, and related phenomena. These simulations generate large, time-dependent datasets that require careful data management for reproducibility, post-processing, and archival.

## **Workflow Overview**

```
Geometry → Meshing → Solver setup → Simulation → Post-processing
CAD (STEP) → Mesh (MSH/polyMesh) → Case setup → Time-step results → Visualization & analysis
```

## **Key Concepts**

### **Mesh Generation**

The computational domain must be discretized into cells (finite volumes) or elements:
- **Structured mesh** — Regular, hexahedral cells; efficient but limited to simple geometries
- **Unstructured mesh** — Tetrahedral/polyhedral cells; flexible for complex geometries
- **Hybrid mesh** — Combination of structured boundary layers and unstructured core
- **Adaptive mesh refinement (AMR)** — Dynamically refines the mesh during simulation in regions of interest

### **Turbulence Modeling**

Most engineering flows are turbulent. Modeling approaches (in order of increasing cost):
- **RANS** (Reynolds-Averaged Navier-Stokes) — Time-averaged; cheapest; most common in industry
- **LES** (Large Eddy Simulation) — Resolves large scales, models small scales; 10-100x more expensive than RANS
- **DNS** (Direct Numerical Simulation) — Resolves all scales; extremely expensive; research only

### **Data Volume Challenges**

A typical transient CFD simulation produces:
- Mesh: 1-100 million cells
- Fields per time step: velocity (3D), pressure, temperature, turbulence quantities
- Time steps: thousands to millions
- **Total size: tens of GB to several TB per simulation**

## **Popular Tools**

### **Solvers**
- **OpenFOAM** — Open-source, widely used in academia
- **ANSYS Fluent / CFX** — Commercial, industry standard
- **SU2** — Open-source, optimization-focused
- **COMSOL Multiphysics** — Multiphysics, user-friendly GUI

### **Meshing**
- **Gmsh** — Open-source mesh generator
- **snappyHexMesh** (OpenFOAM) — Automatic mesh generation from STL surfaces
- **ANSYS Meshing / ICEM CFD** — Commercial mesh tools

### **Post-Processing**
- **ParaView** — Open-source visualization for large datasets
- **Tecplot** — Commercial visualization
- **Python** (PyVista, matplotlib) — Scripted post-processing

## **Code Example: OpenFOAM Workflow**

<details>
<summary>**Click to expand OpenFOAM case setup and run**</summary>

```bash
#!/bin/bash
# OpenFOAM CFD simulation workflow: flow over a backward-facing step

# Step 1: Create case directory structure
mkdir -p backwardStep/{0,constant,system}

# Step 2: Generate mesh
cd backwardStep
blockMesh          # Generate structured mesh from blockMeshDict

# Step 3: Check mesh quality
checkMesh          # Reports mesh statistics, non-orthogonality, skewness

# Step 4: Decompose for parallel run
decomposePar       # Split domain across processors

# Step 5: Run solver in parallel
mpirun -np 8 simpleFoam -parallel > log.simpleFoam 2>&1

# Step 6: Reconstruct parallel results
reconstructPar

# Step 7: Post-process
# Calculate wall shear stress
simpleFoam -postProcess -func wallShearStress -latestTime

# Sample data along a line for comparison with experiments
postProcess -func "sampleDict" -latestTime

# Step 8: Convert to VTK for ParaView visualization
foamToVTK

echo "Simulation complete! Open backwardStep/VTK/ in ParaView"
```

</details>

## **Code Example: Post-Processing with Python**

<details>
<summary>**Click to expand Python post-processing script**</summary>

```python
#!/usr/bin/env python3
"""Post-process CFD results: extract velocity profiles and plot convergence."""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Step 1: Read OpenFOAM log file for residual convergence
def parse_residuals(log_file: Path) -> dict:
    """Parse OpenFOAM log file for residual history."""
    residuals = {"Ux": [], "Uy": [], "p": [], "iteration": []}
    iteration = 0

    with open(log_file, "r") as f:
        for line in f:
            if "Solving for Ux" in line:
                res = float(line.split("Final residual = ")[1].split(",")[0])
                residuals["Ux"].append(res)
                iteration += 1
                residuals["iteration"].append(iteration)
            elif "Solving for Uy" in line:
                res = float(line.split("Final residual = ")[1].split(",")[0])
                residuals["Uy"].append(res)
            elif "Solving for p" in line:
                res = float(line.split("Final residual = ")[1].split(",")[0])
                residuals["p"].append(res)

    return residuals

residuals = parse_residuals(Path("log.simpleFoam"))

# Step 2: Plot convergence
fig, ax = plt.subplots(figsize=(10, 5))
for field in ["Ux", "Uy", "p"]:
    ax.semilogy(residuals["iteration"][:len(residuals[field])],
                residuals[field], label=field)
ax.set_xlabel("Iteration")
ax.set_ylabel("Residual")
ax.set_title("Solver Convergence")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("convergence.pdf", dpi=300)
plt.savefig("convergence.png", dpi=300)

# Step 3: Read sampled velocity profile
data = np.loadtxt("postProcessing/sampleDict/0/line_U.csv",
                   delimiter=",", skiprows=1)
y = data[:, 0]
Ux = data[:, 1]

fig, ax = plt.subplots(figsize=(6, 8))
ax.plot(Ux, y, "b-", linewidth=2)
ax.set_xlabel("Velocity Ux (m/s)")
ax.set_ylabel("y (m)")
ax.set_title("Velocity Profile")
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("velocity_profile.pdf", dpi=300)
plt.savefig("velocity_profile.png", dpi=300)

print("Post-processing complete!")
```

</details>

## **Expected Outputs**

- **Field data** — Velocity, pressure, temperature at each cell for each saved time step
- **Residual history** — Convergence log showing how solver residuals decrease
- **Derived quantities** — Wall shear stress, drag/lift coefficients, heat flux
- **Visualizations** — Contour plots, streamlines, velocity profiles (PDF + PNG)
- **Validation plots** — Comparison of simulation results with experimental data

## **Computational Requirements**

| Simulation Type | Mesh Size | CPU Cores | RAM | Storage | Time |
|----------------|-----------|-----------|-----|---------|------|
| 2D steady RANS | 50K cells | 4 | 2-4 GB | `<1 GB` | Minutes |
| 3D steady RANS | 5M cells | 32-64 | 16-64 GB | 5-20 GB | Hours |
| 3D transient LES | 20M cells | 128-512 | 64-256 GB | 100 GB - 1 TB | Days-weeks |
| DNS (research) | 1B+ cells | 1000+ | 1+ TB | 10+ TB | Weeks-months |

## **Common Issues & Troubleshooting**

:::warning Common Problems

**Solution divergence**
- Reduce under-relaxation factors or CFL number
- Improve mesh quality (check non-orthogonality, skewness)
- Start with first-order schemes, then switch to second-order after partial convergence

**Mesh quality issues**
- Non-orthogonality above 70° causes instability — refine problematic regions
- Ensure adequate boundary layer resolution (y+ appropriate for turbulence model)
- Use `checkMesh` (OpenFOAM) to identify problems before running

**Storage overflow**
- Save only every Nth time step for transient simulations
- Write only fields of interest (not all default fields)
- Use binary output format instead of ASCII
- Compress old time steps and archive to long-term storage

:::

## **Key Considerations**

- **Version-control your case setup** — Input files (dictionaries, boundary conditions) should be in Git
- **Document mesh convergence** — Run at least 3 mesh resolutions to demonstrate grid independence
- **Save restart files** — For long simulations, save checkpoint files at regular intervals
- **Validate against experiments** — Simulation results without validation have limited scientific value
- **Plan storage early** — Estimate total output size before starting and ensure sufficient storage
