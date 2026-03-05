---
title: Finite Element Analysis (FEA)
slug: /domain-knowledge/simulation-modeling/use-cases/finite-element-analysis
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Use Case 2: Finite Element Analysis (FEA)**

Finite element analysis discretizes a physical domain into small elements and solves governing equations (elasticity, heat transfer, etc.) to predict structural behavior. FEA is the workhorse of mechanical and civil engineering, used to simulate stress, deformation, vibration, fatigue, and failure in components and structures.

## **Workflow Overview**

```
Geometry → Meshing → Material & BCs → Solve → Post-process
CAD (STEP/IGES) → FE mesh (INP/MSH) → Loads & constraints → Displacement/stress fields → Reports & plots
```

## **Key Concepts**

### **Element Types**

- **1D elements** — Beams, trusses, rods
- **2D elements** — Shells (thin structures), plane stress/strain
- **3D elements** — Tetrahedra, hexahedra, prisms for solid bodies
- **Higher-order elements** — Quadratic elements with mid-side nodes for better accuracy

### **Analysis Types**

- **Static linear** — Small deformation, linear material; most common
- **Static nonlinear** — Large deformations, contact, plasticity
- **Modal analysis** — Natural frequencies and mode shapes
- **Transient dynamic** — Time-dependent loads (impact, vibration)
- **Thermal analysis** — Heat conduction, convection, radiation
- **Fatigue** — Predict component lifetime under cyclic loading

### **Boundary Conditions**

- **Constraints** — Fixed, pinned, symmetry, or prescribed displacement
- **Loads** — Forces, pressures, temperatures, body forces (gravity)
- **Contact** — Surface-to-surface interaction between parts

## **Popular Tools**

### **Solvers**
- **Abaqus** — Industry-standard for nonlinear structural analysis
- **ANSYS Mechanical** — Comprehensive structural and multiphysics
- **CalculiX** — Open-source, Abaqus-compatible input format
- **FEniCS** (Python) — Open-source, flexible PDE solver
- **Code_Aster** — Open-source structural mechanics solver (EDF)

### **Pre-Processing**
- **SALOME** — Open-source CAD + meshing platform
- **Gmsh** — Open-source mesh generator
- **HyperMesh** (Altair) — Commercial pre-processor

### **Post-Processing**
- **ParaView** — Open-source visualization
- **Abaqus/Viewer** — Abaqus-native post-processing
- **PyVista** (Python) — Scripted 3D visualization

## **Code Example: FEniCS Linear Elasticity**

<details>
<summary>**Click to expand Python FEniCS workflow**</summary>

```python
#!/usr/bin/env python3
"""Linear elasticity: cantilever beam under tip load using FEniCS."""

from fenics import *
import matplotlib.pyplot as plt

# Step 1: Create mesh (beam: 10 x 1 x 1)
mesh = BoxMesh(Point(0, 0, 0), Point(10, 1, 1), 40, 4, 4)

# Step 2: Define function space (vector-valued for displacement)
V = VectorFunctionSpace(mesh, "Lagrange", 1)

# Step 3: Material properties (steel)
E = 210e9      # Young's modulus [Pa]
nu = 0.3       # Poisson's ratio
mu = E / (2 * (1 + nu))
lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))

# Step 4: Boundary conditions — fix left end
def left_boundary(x, on_boundary):
    return near(x[0], 0.0) and on_boundary

bc = DirichletBC(V, Constant((0, 0, 0)), left_boundary)

# Step 5: Define variational problem
def epsilon(u):
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)

def sigma(u):
    return lmbda * nabla_div(u) * Identity(3) + 2 * mu * epsilon(u)

u = TrialFunction(V)
v = TestFunction(V)

# Body force (gravity) and surface traction (tip load)
f = Constant((0, 0, -9810))  # gravity [N/m^3] for steel
T = Constant((0, 0, -1e6))   # tip traction [Pa]

# Mark right end for tip load
class RightEnd(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 10.0) and on_boundary

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
RightEnd().mark(boundaries, 1)
ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

a = inner(sigma(u), epsilon(v)) * dx
L = dot(f, v) * dx + dot(T, v) * ds(1)

# Step 6: Solve
u_sol = Function(V, name="Displacement")
solve(a == L, u_sol, bc)

# Step 7: Post-process
# Calculate von Mises stress
s = sigma(u_sol) - (1.0 / 3) * tr(sigma(u_sol)) * Identity(3)
von_Mises = project(sqrt(3.0 / 2 * inner(s, s)), FunctionSpace(mesh, "DG", 0))

# Maximum displacement
max_disp = u_sol.vector().max()
print(f"Maximum displacement: {max_disp:.6e} m")

# Save results
File("displacement.pvd") << u_sol
File("von_mises.pvd") << von_Mises

print("FEA complete! Open .pvd files in ParaView for visualization.")
```

</details>

## **Code Example: CalculiX Workflow**

<details>
<summary>**Click to expand CalculiX Bash workflow**</summary>

```bash
#!/bin/bash
# CalculiX FEA workflow

# Step 1: Generate mesh with Gmsh
gmsh -3 beam.geo -o beam.msh -format msh2

# Step 2: Convert mesh to Abaqus format (CalculiX uses .inp)
# (Gmsh can directly export .inp, or use a converter)
gmsh -3 beam.geo -o beam.inp -format inp

# Step 3: Edit/prepare the CalculiX input file
# beam.inp should contain:
# - *NODE and *ELEMENT definitions (from mesh)
# - *MATERIAL (e.g., steel: E=210000, nu=0.3)
# - *BOUNDARY (fixed support)
# - *CLOAD or *DLOAD (applied loads)
# - *STEP, *STATIC
# - *NODE FILE, U (request displacement output)
# - *EL FILE, S (request stress output)
# - *END STEP

# Step 4: Run CalculiX solver
ccx beam

# Step 5: Convert results to VTK for ParaView
ccx2paraview beam.frd vtk

# Step 6: Extract max displacement
grep "DISP" beam.dat | sort -k5 -n | tail -1

echo "FEA complete! Open beam.vtk in ParaView"
```

</details>

## **Expected Outputs**

- **Displacement fields** — Deformed shape with displacement magnitudes
- **Stress fields** — Von Mises stress, principal stresses, shear stresses
- **Reaction forces** — Forces at constrained nodes
- **Natural frequencies** (modal) — Eigenvalues and mode shapes
- **Result files** — VTU/VTK for visualization, DAT for extracted values

## **Computational Requirements**

| Analysis Type | Mesh Size | CPU Cores | RAM | Storage | Time |
|---------------|-----------|-----------|-----|---------|------|
| Linear static (single part) | 100K elements | 4-8 | 4-8 GB | `<1 GB` | Minutes |
| Nonlinear static (contact) | 1M elements | 16-32 | 16-64 GB | 1-10 GB | Hours |
| Modal analysis (large assembly) | 5M elements | 32-64 | 64-128 GB | 5-20 GB | Hours |
| Transient dynamic (crash) | 10M+ elements | 64-256 | 128-512 GB | 50-500 GB | Days |

## **Common Issues & Troubleshooting**

:::warning Common Problems

**Convergence failure in nonlinear analysis**
- Reduce load increment size (more substeps)
- Improve contact definitions (adjust penalty stiffness, stabilization)
- Check for rigid body modes — ensure sufficient constraints

**Stress singularities**
- Sharp corners and point loads cause artificially infinite stress
- Use fillets in geometry or interpret results away from singularity
- Stress at singularities does NOT converge with mesh refinement — this is expected

**Mesh-dependent results**
- Always perform a mesh convergence study (refine mesh until results stabilize)
- Use at least 3 mesh densities (coarse, medium, fine)
- Report mesh convergence alongside results

:::

## **Key Considerations**

- **Mesh convergence is mandatory** — Results from a single mesh resolution are insufficient
- **Verify before validate** — First check against analytical solutions, then against experiments
- **Document units consistently** — FEA codes often don't enforce unit systems; mixing units is a common error
- **Store input files with results** — The .inp or .py input file is the key to reproducibility
- **Report element types and formulation** — Reduced integration vs. full integration affects accuracy
