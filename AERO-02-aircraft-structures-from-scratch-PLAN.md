# PLAN.md — AERO-02: Aircraft Structures From Scratch

## Overview
A from-scratch implementation of aerospace structural analysis, progressing from basic stress/strain and beam theory through a self-built finite element solver to composite laminate analysis and fatigue life prediction. Second of five projects covering the classical pillars of aerospace engineering.

## Objective
Build a structural analysis toolkit capable of predicting stress, deflection, buckling loads, and fatigue life for representative aerospace structural elements (beams, panels, laminates), validated against closed-form solutions.

## Why This Project
Structures is arguably the most "engineering" of the classical pillars — it's where aerospace engineering meets solid mechanics and materials science directly. A self-built FEA solver (even a simple 1D/2D one) demonstrates you understand what commercial tools (NASTRAN, ANSYS) are actually solving under the hood, which is rare even among practicing engineers who only operate the software.

## Knowledge Covered
1. **Stress/Strain Fundamentals**: implement stress transformation (Mohr's circle), principal stresses, von Mises equivalent stress.
2. **Beam Theory**: Euler-Bernoulli beam bending and torsion — deflection, shear/moment diagrams for statically determinate and indeterminate beams.
3. **Finite Element Method (from scratch)**: implement a 1D/2D FEA solver for truss and beam elements — stiffness matrix assembly, boundary conditions, solving for displacements and reaction forces. Validate against closed-form beam deflection formulas.
4. **Buckling Analysis**: Euler column buckling, validate critical buckling loads against classical formulas for various boundary conditions.
5. **Composite Laminate Theory**: implement Classical Lamination Theory (CLT) — compute the ABD stiffness matrix for a laminate stack-up, predict ply-by-ply stress/strain under load, apply a failure criterion (e.g., Tsai-Wu).
6. **Fatigue Analysis**: implement S-N curve-based life prediction and Miner's rule for cumulative damage under variable amplitude loading — directly relevant to aircraft structural life management.

## Prerequisites / Background
- Mechanics of materials: stress, strain, Hooke's law
- Linear algebra (stiffness matrix assembly and solving)
- Basic differential equations (beam bending ODEs)

## Data / Dataset
- **Closed-form benchmarks**: beam, torsion, and buckling solutions from standard structures textbooks — the primary validation data for the FEA solver.
- **Metallic material properties — MIL-HDBK-5J** (2003, Distribution Statement A — approved for public release, free on DTIC, EverySpec, and the Internet Archive): design allowables (A- and B-basis) and fatigue S-N data for aerospace aluminum (e.g., 2024-T3, 7075-T6), titanium, and steel alloys. It was the last edition before the handbook moved to the FAA-maintained, paid MMPDS; MIL-HDBK-5J is equivalent to MMPDS-01. Values may since have been revised in later MMPDS editions, which is acceptable for a portfolio project but not for certification work. As U.S. government work, its tables can be reproduced in the repository with citation.
- **Composite ply properties — MIL-HDBK-17-2F** (2002, Volume 2: polymer matrix composite material properties; free on DTIC — confirm the Distribution A statement on the cover page) and **NCAMP** (National Center for Advanced Materials Performance, Wichita State University): statistically based allowables for aerospace prepregs, free after registration on the NCAMP portal. Cite NCAMP values rather than bulk-republishing them, since redistribution terms are not stated.
- **Fatigue data**: S-N curves from MIL-HDBK-5J and from NACA/NASA technical notes on the NASA Technical Reports Server (NTRS, free); laminate worked examples from published composite mechanics textbooks.
- **MMPDS is deliberately not used**: it is a paid subscription.

## Tech Stack
- Python, NumPy, SciPy (sparse linear solvers for FEA)
- Matplotlib for stress contour plots, mode shapes, S-N curves
- Pytest for validation against closed-form solutions

## Repository Structure
```
aircraft-structures-from-scratch/
  src/
    stress_strain/
    beam_theory/
    fea_solver/
    buckling/
    composite_laminates/
    fatigue/
  validation/
  notebooks/
  docs/
  tests/
```

## Execution Plan
- [ ] **Phase 1 — Stress/strain fundamentals**: implement stress transformation and von Mises criterion, validate with textbook examples.
- [ ] **Phase 2 — Beam theory**: implement Euler-Bernoulli bending/torsion for standard load cases, validate against closed-form deflection formulas.
- [ ] **Phase 3 — FEA solver**: build a 1D truss/beam finite element solver from scratch (stiffness matrix assembly, boundary conditions, solve), validate against Phase 2 closed-form results.
- [ ] **Phase 4 — Buckling**: implement Euler buckling analysis for various end conditions, validate against classical critical load formulas.
- [ ] **Phase 5 — Composite laminates**: implement Classical Lamination Theory, compute ABD matrix, apply Tsai-Wu failure criterion to a representative laminate stack-up.
- [ ] **Phase 6 — Fatigue**: implement S-N curve fitting and Miner's rule cumulative damage, apply to a simulated variable-amplitude load spectrum (connects to the ML Structural Health Monitoring project).
- [ ] **Phase 7 — Documentation and publishing**: theory notes per module, visualizations of stress fields, mode shapes, and fatigue life curves.

## Validation Strategy
- Beam deflection/stress results checked against closed-form solutions (cantilever, simply supported, various load cases).
- FEA solver validated by convergence study (mesh refinement should converge to closed-form result).
- Buckling loads validated against classical Euler formula for pinned-pinned, fixed-free, fixed-fixed conditions.
- Laminate analysis validated against published composite mechanics textbook examples.

## Deliverables
- Public repository with a working, documented FEA solver and composite/fatigue analysis modules
- Validation report comparing all modules against closed-form/published results
- Visualizations: stress contours, buckling mode shapes, ply-by-ply laminate stress, S-N/fatigue life curves

## Portfolio Differentiators
- A hand-built FEA solver is a strong, concrete signal of structural mechanics fluency — most portfolios show FEA *usage*, not FEA *implementation*.
- The fatigue module connects directly to real airline reliability work (S-N/Miner's rule underlies real aircraft structural life limits), tying this project back to your day-to-day role at GOL.

## Possible Extensions
- Extend the FEA solver to 2D plate/shell elements
- Add a probabilistic fatigue model (accounting for load scatter) instead of deterministic S-N/Miner's rule
