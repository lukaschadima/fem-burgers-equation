# 1D FEM Solver for the Viscous Burgers' Equation

A finite element solver for the 1D viscous Burgers' equation, implemented in MATLAB
as part of my bachelor's thesis at Charles University (Faculty of Mathematics and
Physics, General Mathematics / Numerical Analysis).

The project implements a spatial discretization with Lagrange finite elements of
arbitrary order, combined with a range of time-stepping schemes — from simple
Euler methods up to high-order fully implicit and IMEX Runge-Kutta schemes — and
includes convergence and stability studies comparing them.

## The problem

The viscous Burgers' equation is a nonlinear PDE that combines advection and
diffusion:

```
u_t + u * u_x = nu * u_xx + f(x, t),   x in (0, L), t in (0, T]
```

with Dirichlet boundary conditions `u(0,t) = u(L,t) = 0`. It's a standard testbed
for numerical methods because it's simple to state but exhibits real nonlinear
and (for small `nu`) stiff behavior, similar in flavor to the Navier-Stokes
equations.

Two test problems are included, both with known exact solutions used to measure
convergence:

- **Example 1 (Cole-Hopf):** an exact solution derived via the Cole-Hopf
  transformation, which linearizes Burgers' equation into the heat equation.
- **Example 2 (Manufactured solution):** `u(x,t) = sin(pi*x)*sin(omega*t)`,
  with a forcing term `f` derived to make this an exact solution — useful for
  isolating and testing time-integration accuracy independently of spatial error.

## What's implemented

**Spatial discretization (`src/core/`)**
- Lagrange finite elements of order `p = 1, 2, 3, 4` on a uniform mesh
- Gauss-Legendre quadrature (closed-form for low order, Newton's method on
  Legendre polynomials for higher order)
- Assembly of the mass matrix, stiffness matrix, nonlinear advection term, and
  its Jacobian (needed for Newton iterations in implicit schemes)

**Time integration (`src/solvers/`)**
- Explicit Euler, Implicit Euler, Crank-Nicolson (baseline methods)
- Explicit Runge-Kutta: ERK1–ERK4
- Fully implicit SDIRK2 / SDIRK3 (Newton iteration on the full nonlinear system
  at every stage)
- IMEX (implicit-explicit) ARS schemes: ARS(1,1,1), ARS(1,2,2), ARS(2,3,2),
  ARS(3,4,3) — diffusion treated implicitly (linear solve only, no Newton),
  advection treated explicitly

The IMEX schemes are the practically interesting part: they avoid the harsh
`dt ~ h^2` stability restriction that explicit methods have on the diffusion
term, without paying for a full nonlinear solve at every stage the way fully
implicit DIRK methods do.

**Analysis (`analysis/`)**
- Spatial convergence study across polynomial order `p` and mesh size `n`
- Time-stepping stability and convergence-order comparison across all methods
- Accuracy-vs-CPU-time efficiency comparison (which method gives the best
  accuracy per unit of compute)

## Repository structure

```
fem-burgers-equation/
├── setup_paths.m              <- run this first (see below)
├── src/
│   ├── core/                  <- FEM assembly, basis functions, quadrature
│   ├── solvers/                <- time-stepping schemes, Butcher tables
│   └── problems/               <- test problem definitions (IC, exact sol., source)
├── analysis/                   <- convergence and stability studies
└── results/                    <- example output plots (PDF)
```

## Running it

Requires MATLAB (developed on R2023a+, no special toolboxes needed).

```matlab
setup_paths                  % adds src/ and analysis/ to the MATLAB path

% Example: run the full stability + convergence comparison across all
% time-stepping schemes (produces PDF plots in the current folder)
test_stability_new

% Example: pure spatial convergence study (p = 1,2,3 vs mesh size)
test_spatial_convergence

% Example: solve directly and inspect the solution
params = burgers_params(2, 0.05, 10);   % example 2, nu=0.05, omega=10
[x, U_history, t] = solve_burgers_equation(2, 32, 0.001, params, 'imex_ars232');
```

## Notes

- `analysis/test_stability_legacy.m` is an earlier, simpler version of the
  stability comparison (3 methods); `test_stability_new.m` is the current,
  more complete version (10 methods, cleaner plots). Kept both since the
  legacy version is smaller and easier to read as a starting point.
- This is bachelor's thesis work in progress — the codebase will keep
  evolving as the thesis is finalized.

## License

MIT — see [LICENSE](LICENSE).
