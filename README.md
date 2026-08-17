# 1D FEM Solver for the Viscous Burgers Equation

MATLAB implementation of a finite element solver for the 1D viscous 
Burgers equation, built as part of my bachelor thesis in numerical 
mathematics at Charles University (MFF UK).

## What it does

- Solves the viscous Burgers equation using the finite element method, 
  implemented from scratch (no black-box solvers)
- Implements and compares four time-stepping schemes: implicit Euler, 
  explicit Euler, Crank-Nicolson, and IMEX
- Verifies correctness using manufactured solutions and the Cole-Hopf 
  transformation (which maps the nonlinear Burgers equation to the 
  linear heat equation)
- Includes convergence rate analysis across all four schemes

## Why

The Burgers equation is a standard testbed for numerical methods 
because it combines nonlinear advection with diffusion — the same 
mathematical structure that shows up in more complex fluid models. 
Comparing explicit, implicit, and IMEX schemes highlights the 
stability/accuracy tradeoffs that matter in real PDE solvers.

## Status

Work in progress — part of an ongoing bachelor thesis (2025–present), 
supervised by doc. Václav Kučera.
