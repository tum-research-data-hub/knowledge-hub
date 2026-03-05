---
title: Glossary
slug: /domain-knowledge/simulation-modeling/glossary
---

:::caution Work in Progress
This page is currently under development. Content may be incomplete or contain inaccuracies. If you notice any errors or have suggestions, please [contact us](mailto:rdhub@mdsi.tum.de).
:::

# **Simulation & Modeling Glossary**

This glossary provides definitions for common terms used in computational simulation, numerical modeling, and related data management.

---

## **Numerical Methods**

### **Finite Element Method (FEM / FEA)**
A numerical technique that divides a domain into small elements and solves governing equations (e.g., elasticity, heat transfer) approximately on each element, then assembles a global solution.

### **Finite Volume Method (FVM)**
A numerical method that divides the domain into control volumes and ensures conservation of quantities (mass, momentum, energy) across each volume. The basis of most CFD solvers.

### **Finite Difference Method (FDM)**
A numerical method that approximates derivatives using differences between values at discrete grid points. Simple but limited to structured grids.

### **Navier-Stokes Equations**
The fundamental partial differential equations describing the motion of viscous fluids. Solving these equations (in simplified or full form) is the core of CFD.

### **PDE (Partial Differential Equation)**
A differential equation involving partial derivatives of an unknown function with respect to multiple independent variables. Most physics-based simulations solve PDEs.

---

## **Mesh and Discretization**

### **Mesh (Grid)**
A discretization of the computational domain into small cells (finite volumes) or elements (finite elements) on which equations are solved.

### **Structured Mesh**
A mesh with regular connectivity (each interior cell has the same number of neighbors), typically using hexahedral or quadrilateral cells. Efficient but limited to simple geometries.

### **Unstructured Mesh**
A mesh with irregular connectivity, typically using tetrahedral or polyhedral cells. Flexible for complex geometries but less computationally efficient.

### **Mesh Convergence (Grid Independence)**
The process of demonstrating that simulation results do not change significantly when the mesh is refined. Essential for credible results.

### **Adaptive Mesh Refinement (AMR)**
A technique that dynamically refines (or coarsens) the mesh during simulation in regions where higher (or lower) resolution is needed.

### **Boundary Layer Mesh**
A region of refined, thin cells near walls to resolve the steep gradients in velocity and temperature that occur in boundary layers.

### **y+ (y-plus)**
A dimensionless wall distance used to determine whether the near-wall mesh resolution is appropriate for the chosen turbulence model.

---

## **CFD-Specific Terms**

### **Turbulence**
Chaotic, irregular fluid motion characterized by eddies, vortices, and rapid fluctuations. Most engineering flows are turbulent.

### **RANS (Reynolds-Averaged Navier-Stokes)**
A turbulence modeling approach that solves time-averaged equations with turbulence models (e.g., k-ε, k-ω SST) to close the system.

### **LES (Large Eddy Simulation)**
A turbulence approach that directly resolves large-scale eddies and models only the smallest scales. More accurate than RANS but 10-100x more expensive.

### **DNS (Direct Numerical Simulation)**
Solving the Navier-Stokes equations without any turbulence modeling, resolving all scales of motion. Extremely expensive; limited to low Reynolds numbers and research.

### **CFL Number (Courant-Friedrichs-Lewy)**
A dimensionless number relating time step size, cell size, and flow velocity. Must be below a threshold (typically 1) for stability in explicit time-stepping schemes.

### **Residual**
A measure of how well the discrete equations are satisfied at a given iteration. Decreasing residuals indicate convergence toward a solution.

### **Under-Relaxation**
A technique that limits the change in a variable between iterations to improve solver stability, at the cost of slower convergence.

---

## **FEA-Specific Terms**

### **Degrees of Freedom (DOF)**
The number of independent values that can vary in the system. In structural FEA, each node typically has 3 translational and 3 rotational DOF.

### **Von Mises Stress**
A scalar stress measure derived from the stress tensor that predicts yielding in ductile materials under complex loading.

### **Eigenvalue / Natural Frequency**
In modal analysis, the eigenvalues of the structural system correspond to the squares of the natural frequencies at which the structure vibrates.

### **Mode Shape**
The characteristic deformation pattern of a structure at a natural frequency.

### **Contact**
The interaction between surfaces of different bodies that come into contact during deformation. Adds significant nonlinearity to the analysis.

### **Stress Singularity**
A point where stress theoretically approaches infinity (e.g., sharp re-entrant corner). An artifact of the mathematical model, not physical reality.

---

## **General Simulation Terms**

### **Steady-State**
A solution that does not change with time. The system has reached equilibrium.

### **Transient (Unsteady)**
A time-dependent simulation where the solution evolves over time.

### **Boundary Condition**
A constraint applied at the boundary of the computational domain (e.g., fixed wall, inlet velocity, prescribed temperature).

### **Initial Condition**
The state of the system at the start of a transient simulation.

### **Convergence**
The process by which an iterative solver approaches the solution. A simulation has converged when residuals are sufficiently small and monitored quantities have stabilized.

### **Validation**
Comparing simulation results against experimental measurements to assess accuracy.

### **Verification**
Confirming that the numerical implementation correctly solves the mathematical model (code correctness), typically against analytical or benchmark solutions.

### **Checkpoint / Restart File**
A saved state of the simulation at a specific time, allowing the simulation to be resumed if interrupted.

### **Post-Processing**
The extraction, visualization, and analysis of results after a simulation has completed.

### **HPC (High-Performance Computing)**
The use of parallel computing resources (clusters, supercomputers) to run large-scale simulations.

### **MPI (Message Passing Interface)**
A standard for parallel computing that allows processes running on different nodes to communicate. Used by most simulation solvers for parallel execution.

---

**Need clarification on a term?** Contact us at [rdhub@mdsi.tum.de](mailto:rdhub@mdsi.tum.de) or explore our [Use Cases](/domain-knowledge/simulation-modeling/use-cases) for practical examples.
