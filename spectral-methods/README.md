# Spectral collocation — an alternative to the FEM solver

Before settling on the finite element solver in `src/`, I explored solving the
same viscous Burgers' equation with a **pseudospectral (Chebyshev collocation)**
discretization instead. This folder holds that exploration. It's not part of
the thesis's final solver, but it's the same PDE, verified against the same two
exact solutions, so it's a useful side-by-side of a different numerical
approach to the same problem.

## The two test problems (same as `src/problems/`)

- **Cole-Hopf** — `u_t + u·u_x = ν·u_xx`, no source term, exact solution via the
  Cole-Hopf transformation.
- **Manufactured solution** — `u(x,t) = sin(πx)·sin(ωt)`, with a forcing term
  `f(x,t)` derived so this is exact; isolates time-integration error from
  spatial error.

## Progression

Five scripts, in the order I actually wrote them, each one answering a
question the previous one raised:

1. **`cole_hopf_verification_pde15s.m`** and
   **`manufactured_solution_pde15s.m`** — call Chebfun's built-in `pde15s`
   solver directly. These don't implement anything themselves; they exist to
   check that both exact solutions are actually correct and that the problem
   is well-posed, before spending time hand-rolling a solver for it.
2. **`cole_hopf_imex_euler.m`** and
   **`manufactured_solution_imex_euler.m`** — my own Chebyshev collocation:
   differentiation matrices from `diffmat`, first-order IMEX Euler in time
   (diffusion implicit, advection explicit), Dirichlet boundary rows
   overwritten by hand. This is the first version I wrote myself rather than
   delegating to a library solver.
3. **`manufactured_solution_imex_ars343.m`** — the same spatial discretization,
   but with a proper 4-stage, 3rd-order IMEX Runge-Kutta scheme
   (ARS(3,4,3), Ascher–Ruuth–Spiteri 1997) instead of first-order Euler, coded
   as a general-purpose stage loop rather than four separate hard-coded
   updates.

## Why FEM in the end

Chebyshev collocation converges faster per degree of freedom on smooth
solutions, but it's a global method — every basis function is nonzero across
the whole domain, so the differentiation matrices are dense and local
mesh refinement isn't possible. FEM's local, element-wise basis is what the
thesis needed: refining near the initial boundary layer, and (later) treating
h- and p-refinement on equal footing. The Chebyshev experiments above were
what convinced me of that trade-off rather than assuming it.

## Running it

Requires MATLAB plus the [Chebfun](https://www.chebfun.org/) toolbox (all
five scripts use `chebfun`, `chebpts`, and/or `diffmat`). Each script is
standalone — run it directly, it plots the solution against the exact
reference and reports the max-norm error.

`manufactured_solution_imex_ars343.m` takes noticeably longer to run than the
others: it uses `dt = 1e-4` over `T = 10.67`, i.e. ~107,000 time steps, each
solving several linear systems.

## Note on this folder specifically

Unlike the rest of this repository, these five scripts haven't been re-run
during the 2026 GitHub cleanup — I don't have MATLAB/Chebfun in the
environment I used to review the rest of this repo, so this folder is
reviewed for mathematical correctness (the source term for the manufactured
solution checks out by hand, the ARS(3,4,3) Butcher tableau matches the
published coefficients) but not for a clean run. Worth actually re-running
before relying on it for anything beyond reading the code.
