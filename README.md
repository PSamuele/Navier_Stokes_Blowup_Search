# AI-Assisted CFD on the Cloud

**A verification case study: a vortex ring pushed into a narrowing axisymmetric domain**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue)](LICENSE)
![Stack](https://img.shields.io/badge/Stack-FEniCSx%20%7C%20gmsh%20%7C%20PETSc%20%7C%20MPI-orange)
![Compute](https://img.shields.io/badge/Compute-AWS%20EC2-yellow)
![Tests](https://img.shields.io/badge/Tests-50%20regression-green)
![Environment](https://img.shields.io/badge/Environment-conda-lightgrey)

## What this project is about

I started this project with a question about method, not about physics. How far
can one engineer get on a hard fluid dynamics problem when most of the code is
written with an AI assistant (Claude) and the computing power is rented by the
hour on Amazon Web Services (AWS)? And what does it take to actually trust the
answer?

As a test case I picked a vortex ring pushed into a narrowing, rotationally
symmetric domain. The setup is inspired by research on whether the Navier–Stokes
equations, the equations that describe how a liquid or a gas moves, can produce
infinite values in a finite time. That research question is explained
[below](#the-physical-question).

What I found, in short:

- **The tools made code and computing power cheap.** The final three-grid study
  used 15.7 hours of solver time on one AWS machine, about 11 US dollars at the
  hourly price I used.
- **The hard part moved to checking the results.** My first two attempts
  reported a dramatic "blow-up" that was not real. Both ran to the end without
  a single error message. I later found five separate defects behind it.
- **This repository is the third attempt.** It adds checks that compare what
  the program produces against things that are known independently: an exact
  formula, an exact volume, a physical law that must hold. With these checks
  the simulation is numerically consistent on the two finer grids.
- **On the physics, the result is limited.** The vortex never reached the
  narrow tip of the domain, so the effect the setup was built to test was
  never actually tested. I explain why in
  [The vortex never reached the pole](#the-vortex-never-reached-the-pole).

---

## Contents

| | |
| :-- | :-- |
| [Key terms](#key-terms) | plain definitions of every technical word used here |
| [How the work was done](#how-the-work-was-done) | who did what, the role of the AI, the AWS setup |
| [The physical question](#the-physical-question) | where the idea comes from, and what a simulation can and cannot say |
| [The setup](#the-setup) | domain, initial flow, parameters |
| [Three attempts](#three-attempts) | what happened in each run |
| [Five defects that raised no error](#five-defects-that-raised-no-error) | what went wrong in the first two runs |
| [Five checks built into the code](#five-checks-built-into-the-code) | how the third run guards against the same problems |
| [The convergence study](#the-convergence-study) | how much the results change with the mesh, and why |
| [What the results mean](#what-the-results-mean) | what can and cannot be concluded |
| [Limitations](#limitations) and [next steps](#next-steps) | |
| [Reproducing it](#reproducing-it) | commands |
| [References](#references) | |

---

## Key terms

I use these words throughout. Each one is explained here once, in plain terms.

**Physics**

| term | meaning |
| :-- | :-- |
| Navier–Stokes equations | The equations for the motion of a fluid that has viscosity. They say how velocity and pressure change in time. |
| Euler equations | The same equations with viscosity set to zero (an ideal fluid with no internal friction). |
| Viscosity (`ν`) | The internal friction of a fluid. It smooths out sharp differences in velocity. Here `ν = 0.001`. |
| Incompressible | The fluid cannot be squeezed: its density stays constant. Mathematically, the **divergence** of the velocity is zero (as much fluid flows into any small volume as flows out of it). |
| Axisymmetric | Symmetric under rotation around an axis. The 3D flow is fully described by what happens in one half plane, with coordinates `r` (distance from the axis) and `z` (height). This turns a 3D problem into a 2D one. |
| Swirl, `u_θ` | The part of the velocity that goes around the axis, like water spinning in a drain. |
| Vortex ring | A ring of spinning fluid, like a smoke ring. |
| Vorticity `ω` | A measure of how fast the fluid is locally spinning. Very large vorticity in a very small region is what a singularity would look like. |
| Circulation `Γ = r·u_θ` | Swirl multiplied by the distance from the axis. Without viscosity it stays constant as a fluid particle moves. So if a particle moves closer to the axis (smaller `r`), its swirl `u_θ = Γ/r` must grow. This is the amplification effect the whole setup relies on. |
| Kinetic energy `E` | The total energy of motion of the fluid, `½∫|u|²`. |
| Enstrophy | The integral of `|ω|²` over the domain. It grows when small, intense vortices form. |
| Finite-time singularity (blow-up) | The case where a quantity such as vorticity becomes infinite at a finite time, starting from a smooth, well-behaved flow. |
| Reynolds number `Re` | `Re = U·L/ν` (typical speed × typical size / viscosity). It compares inertia with viscosity. Large `Re` means viscosity matters only in thin regions. Here `U ≈ 23`, `L = 1`, so `Re ≈ 2.3 × 10⁴`. |
| No-slip | The fluid touching a solid wall has zero velocity. |
| Boundary layer | The thin layer next to a wall where the velocity drops to zero because of no-slip. Its thickness is roughly `δ ≈ L/√Re`, here about `0.0066`. |

**Numerics**

| term | meaning |
| :-- | :-- |
| Mesh (grid), cell, `h` | The domain is split into small triangles (cells). `h` is the size of a cell. Smaller `h` means more detail and more cost. |
| Finite element method (FEM) | The numerical method used here. Inside each triangle the solution is approximated by a simple polynomial. |
| P2 / P1 | The polynomial degree: velocity is quadratic (P2) inside each triangle, pressure is linear (P1). This is a standard stable pairing. |
| DOF (degree of freedom) | One unknown number the solver computes. The fine grid has 4.5 million velocity DOFs. |
| Time step `dt` | The simulation advances in small jumps of time `dt`. |
| CFL number | `u·dt/h`: how many cells the fluid crosses in one time step. If it gets too large the computation becomes unstable. I keep it at 0.5. |
| IPCS | Incremental Pressure Correction Scheme. A standard way to solve incompressible flow. Each time step does three simpler solves: (1) predict the velocity using the old pressure, (2) solve for a pressure correction, (3) correct the velocity so it is divergence free. The version here is first order in time (the error from the time stepping shrinks in proportion to `dt`). |
| Linear solver, preconditioner | Each of the three steps above is a large system of linear equations. PETSc solves them iteratively. A preconditioner is a helper that makes the iterations converge faster. |
| Round-off | The tiny errors a computer makes because it stores numbers with limited precision (about 16 digits). |
| FEniCSx, gmsh, PETSc, MPI | FEniCSx: the finite element library. gmsh: the mesh generator. PETSc: the linear solver library. MPI: the standard that lets one simulation run on many processor cores at once (here 8). |

**Checking the results**

| term | meaning |
| :-- | :-- |
| Verification | Checking that the program solves the equations correctly (as opposed to checking that the equations describe reality). |
| Convergence study | Running the same case on meshes of different sizes. If the results stop changing as the mesh gets finer, they can be trusted. |
| Observed order `p` | How fast the error shrinks when `h` is reduced. With `p = 2`, halving `h` divides the error by 4. With `p = 1`, by 2. |
| Richardson extrapolation | Using results from three meshes to estimate the value on an infinitely fine mesh. |
| GCI (Grid Convergence Index) | An error bar on the fine-mesh value, computed from the differences between meshes and the observed order, with a safety factor of 1.25 (Roache's method). |
| Monotone fraction | The fraction of time samples in which the coarse, medium and fine values move steadily in one direction. If they zigzag, the three meshes are not yet in the range where the error behaves predictably, and the GCI and the extrapolation mean nothing. |
| BKM criterion | Beale–Kato–Majda (1984): a smooth solution of the Euler equations can only break down at time `T*` if the integral of the maximum vorticity, `∫₀^T* max|ω| dt`, becomes infinite. Similar criteria exist for Navier–Stokes. |

---

## How the work was done

**Who did what.** I chose the problem and the geometry, defined the physics and
the parameters, sized the cloud machines, ran the three campaigns and read the
output. Claude wrote most of the code. In the first two attempts it also wrote
the code that contained the five defects listed
[below](#five-defects-that-raised-no-error). Later, following my direction, it
did much of the analysis that found them, and it wrote the 50 regression tests
(automatic tests that fail if a fixed problem comes back).

**How I worked with the AI.** I did not take the first answer to any
substantive question. I pushed back on it until it either held or broke, and
several times it broke. An AI assistant gives a confident, well written answer
to almost any question, including the ones where it is wrong, and the wrong
answers look exactly like the right ones. What the tool does not provide is
judgement about the physics: deciding that a number is wrong before any analysis
says so, choosing which checks are worth building, and deciding when the work is
done.

**What actually caught the errors.** None of the five defects raised an error,
and none would have been found by reading the code. Each one was found by
comparing the program's output against something known independently:

- the initial flow, which is given by an exact formula;
- the volume of the domain, which can be computed exactly;
- the element size gmsh actually delivered, compared with the size I asked for;
- a physical law that must hold (kinetic energy cannot grow in this setup).

This is the main lesson of the project for me: writing code became almost free,
but checking it did not. Without checks, the same tools just produce wrong
results faster.

**The AWS setup.** The production runs used one `c6i.4xlarge` machine on AWS EC2
(Amazon's rented virtual machines): 16 virtual CPUs, which are 8 physical cores.
I ran 8 MPI processes, one per physical core. Before committing to the full
study, a `--benchmark` mode runs a few steps on the real machine and estimates
wall time, memory and cost. Setup scripts and cost control are in
[`deploy/`](deploy/).

| mesh level | solver wall time |
| :-- | --: |
| coarse | 5.1 min |
| medium | 1.6 h |
| fine | 14.0 h |
| **total** | **15.7 h**, about $11 at $0.68/hour |

The cost covers solver time only, at the on-demand hourly price stored in
[`configs/aws_production.json`](configs/aws_production.json). It does not include
setup time or storage.

*A different use of AI in this field.* Wang, Lai, Gómez-Serrano and Buckmaster
(2023) used neural networks as the numerical method itself, to find blow-up
profiles. That is not what happened here. Here an AI assistant wrote ordinary
finite element code.

---

## The physical question

**The open problem.** Whether smooth solutions of the 3D incompressible
Navier–Stokes equations can develop a finite-time singularity is one of the
Clay Millennium Prize Problems. It is unsolved.

**What a simulation can say.** Nothing that settles it. The official problem is
about the whole space, a periodic box or smooth domains, not a specially built
geometry like this one. And a computer simulation, with finite precision and
finite resolution, is evidence at best, never a proof. Whenever a result here
looks suggestive, the correct reading is "the data are consistent with", never
"this shows".

**Where the idea comes from.** Thomas Hou, Jiajie Chen and collaborators have
spent about ten years turning "can it blow up?" into something that can be
computed.

| work | equations | domain | viscosity | outcome |
| :-- | :-- | :-- | :-- | :-- |
| Luo & Hou 2014 | Euler | cylinder with a solid wall | none | blow-up on a ring on the wall |
| Hou 2022 | Euler | interior (no wall) | none | possible singularity at the origin |
| Hou & Huang 2023 | Euler / Navier–Stokes | — | viscosity that goes to zero near the singular point | possible self-similar singularity |
| Hou 2023 | Navier–Stokes | — | constant | possibly singular behaviour, vorticity grows by 10⁷ |
| Chen & Hou 2025 | Euler | smooth data and boundary | none | computer-assisted proof of blow-up |

Two points from this table shaped how I read my own result.

1. **Constant viscosity does not rule it out by itself.** Hou (2023) reports
   possibly singular behaviour with ordinary, constant viscosity. So "viscosity
   always wins" is not a valid explanation for a null result.
2. **The initial flow and the resolution are what matter.** In that work the
   singular behaviour appears only when the starting flow is carefully built to
   match a *self-similar* profile (a shape that keeps the same form while
   shrinking in size), and only with *dynamic rescaling*, a technique that keeps
   zooming the mesh into the shrinking region. That reaches resolutions no fixed
   mesh can reach.

My run uses a generic starting flow, a fixed mesh and `Re ≈ 2.3 × 10⁴`. Nothing
in that regime is expected to become singular. So even a clean negative result
would only confirm what the literature already suggests.

---

## The setup

**The domain.** The fluid fills the solid of revolution of

```
f(z) = R₀ · cos(πz / 2H) · exp(−k·z²),   R₀ = 1,  H = 2,  k = 0.5
```

around the `z` axis. `f(z)` is the radius of the domain at height `z`. It is 1
at the middle (`z = 0`) and goes to zero at the two poles (`z = ±2`). The outer
boundary is a solid wall with no-slip, and nothing pushes the fluid from
outside: no inflow, no outflow, no external force.

![The domain](media/domain/domain-3d.png)

**The initial flow** has three parts, all given by exact formulas in
[`src/ic.py`](src/ic.py):

- A **background flow**, described by a *streamfunction* (a single function
  from which the velocity in the `r–z` plane is computed, which guarantees the
  flow is incompressible): `ψ_jet = J·r²·(f² − r²)²`, with `J = 10`.
- A **vortex ring** centred at `r = 0.5`, `z = 0`, with size `σ = 0.2`.
- **Swirl** on the ring: circulation `Γ₀ = S·r²·exp(−((r−0.5)² + z²)/σ²)`,
  with `S = 20`.

**The idea.** If the ring is carried toward a pole, where the domain narrows, its
fluid is forced closer to the axis. Circulation `Γ` is nearly conserved, so the
swirl `u_θ = Γ/r` grows as `r` shrinks. In the equations for axisymmetric flow,
the term that creates new vorticity from swirl is proportional to
`(1/r³)·∂(Γ²)/∂z`, so the same difference in swirl has a much stronger effect
close to the axis. The question was whether this amplification could beat
viscosity.

**Parameters.** `ν = 0.001`, final time `T = 0.55`, CFL number 0.5.

### The domain is a cone, not a cusp

In 3D renders the domain looks like it ends in a sharp spike, and I originally
called it a cusp. It is not.

- A **cusp** is a tip where the two sides meet with the same tangent, so the
  opening angle shrinks to zero at the tip. Near the tip the radius behaves like
  `r ≈ c·(H − z)^α` with `α > 1`.
- A **cone** has a fixed opening angle all the way down. That is `α = 1`.

Near `z = H` the profile goes to zero linearly (`α = 1`). Its slope at the pole
is `f'(H) = −(πR₀/2H)·exp(−kH²) = −0.106`, which is a cone with a half-angle of
**6.07°**. The `exp(−kz²)` factor does not change this: at `z = H` it is just a
constant, not zero, so it scales the shape but does not change how it narrows.

![The domain narrows linearly: a 6.07° cone, not a cusp](media/domain/cone-not-cusp.png)

A 6° cone and a cusp look the same at normal zoom, so this has to be checked with
numbers. The test `tests/test_mesh.py::test_domain_is_a_cone_not_a_cusp_at_the_poles`
checks it.

**A cone is actually the better choice**, for three reasons:

1. **The mathematics still works.** The standard existence theory for
   Navier–Stokes assumes the boundary is *Lipschitz*: roughly, locally it looks
   like the graph of a function with bounded slope. A cone of fixed angle is
   Lipschitz. A cusp is not, so the standard theory does not apply to it.
2. **A result can be interpreted.** A cusp is already singular as a shape. Any
   singular behaviour found at its tip could simply come from the shape itself,
   not from the fluid.
3. **It can be meshed in a controlled way.** A convergence study needs a set of
   meshes that are all refined by the same factor. Near a cusp tip the cells
   would have to get more and more stretched without limit, so that is not
   possible.

To be honest about this: this was not my original reasoning. The shape was
inherited from the earlier runs, where I believed it was a cusp, and I kept it
so all three runs stay comparable. The reasons above came afterwards.

The production mesh at the pole, with the exact cone drawn on top:

![The Run 3 mesh at the pole](media/domain/pole-zoom.png)

---

## Three attempts

<table>
<tr>
<td width="33%"><img src="archive/run_01/results_R1/media_R1/vortex_blowup_R1.gif" alt="Run 1"></td>
<td width="33%"><img src="archive/run_02/results_R2/media_R2/vortex_blowup_R2.gif" alt="Run 2"></td>
<td width="33%"><img src="results/convergence_aws/fine/vortex_blowup_fine.gif" alt="Run 3, fine grid"></td>
</tr>
<tr>
<td><b>Run 1.</b> On a local machine, one core, to <code>T = 0.55</code>. Final
max |u| = 11.85, which is plausible. Final max |ω| = 3.12 × 10⁶, which is not:
on that mesh the largest vorticity the cells can represent is about 10³.</td>
<td><b>Run 2.</b> On AWS, 16 processes, time step adapted to the flow. Meant to
be the high-resolution version of Run 1. It reached |u| = 9.4 × 10²¹ without
stopping. This "blow-up" is not real. The mesh was in fact the same as Run 1.</td>
<td><b>Run 3, fine mesh.</b> 752,803 cells, 4.5 million velocity DOFs, 14 hours
on 8 cores. Reaches <code>T = 0.55</code>. Kinetic energy goes down in every
one of the 550 intervals between its 551 recorded samples.</td>
</tr>
</table>

Run 2 produces the more convincing animation, and it is the wrong one. Run 3 is
this repository: a rewritten solver, a rewritten mesh generator, diagnostics that
handle the symmetry axis correctly, 50 regression tests and a three-mesh
convergence study with error bars. Runs 1 and 2 are kept unchanged in
[`archive/`](archive/), each with a note on what is wrong with it.

---

## Five defects that raised no error

None of these produced an error message. Each produced numbers that looked
reasonable. Full details and measurements are in [`docs/findings.md`](docs/findings.md).

| # | what went wrong | measured effect |
| --: | :-- | :-- |
| 1 | **The mesh was never refined.** The element size was set by a gmsh formula that used the coordinate `z`. But the 2D geometry was drawn in the x–y plane, so gmsh's `z` was zero everywhere and the formula gave the same size everywhere. | asked for `h = 1e-4` at the poles, got `0.015`: **150 times too coarse**. Runs 1 and 2 used the same mesh file. |
| 2 | **The vorticity was computed wrongly on the axis.** Part of the vorticity is `u_θ/r`. On the axis both are zero, and the true ratio is finite. The code computed `u_θ/(r + 10⁻¹⁴)` at points lying exactly on the axis, where `u_θ` was not exactly zero but a tiny round-off value. Dividing that by `10⁻¹⁴` multiplies it by `10¹⁴`. | **1373 times** too large, already on the exact initial flow, before any time step |
| 3 | **The "start of instability" time was noise.** A script flagged the first sample where the growth rate of the velocity passed a threshold. | it fired on a **0.665 %** change in velocity over one 25 µs sample, while vorticity was *decreasing* |
| 4 | **Part of the equations was not updated.** The matrix of step 1 depends on the current velocity, but it was rebuilt only when the time step changed. | the matrix was out of date for **62.8 %** of the run |
| 5 | **The velocity correction ignored the boundary conditions.** Step 3 of the scheme was solved without imposing no-slip on the wall or the conditions on the axis. | this is what produced the non-zero `u_θ` on the axis in defect 2 |

What they have in common: reading the code does not reveal any of them. Measuring
does. Reading the mesh script does not show that `z` is zero. Measuring the size
of the cells that came out shows it immediately.

![The mesh that was asked for versus the mesh that was delivered](media/domain/mesh-comparison.png)

---

## Five checks built into the code

Each check below is useful in general, not only as a fix for one past mistake.
More detail in [`docs/methodology.md`](docs/methodology.md).

**1. Energy check.** In a closed container with no-slip walls and no external
force, the exact solution loses kinetic energy through viscosity and can never
gain any: `dE/dt = −2ν·∫|D(u)|² ≤ 0`, where `D(u)` is the rate of deformation of
the fluid. So if the computed energy goes up, the computation is wrong, by
definition. The check needs no tuning and no knowledge of the flow. The run stops
at the first sample that is more than 1 % above the lowest energy seen so far,
and the time of that lowest energy is recorded as the last reliable time. A CFL
check cannot do this job: when the mesh is too coarse for the flow, the CFL
number can stay perfectly on target while the results go wrong.

Important: this is a check on the numerics, not on the physics. It tells me the
computation is still consistent. It says nothing about whether the true solution
would blow up (see [What the results mean](#what-the-results-mean)).

**2. Vorticity computed away from the axis.** Instead of adding a small number
to `r`, the code evaluates vorticity at the *quadrature points* inside each
triangle (the points where the method evaluates integrals). These points lie
strictly inside the triangle, so they can never be on the axis, and `r > 0`
there. No small number is needed. `Diagnostics.audit()` checks this on the actual
mesh and refuses to run if it fails. On the exact initial flow it gives **352.59**
against the exact **351.64**, an error of 0.27 %.

**3. Mesh check.** After building the mesh, the generator measures the cell size
it actually got and stops with an error if the size at the poles is far from what
was asked. Run on the old recipe, it fails immediately.

**4. Checked linear solves.** When PETSc's preconditioner fails, PETSc does not
stop. It returns a negative "converged reason" code and can leave infinite values
in the solution. This happened during the rewrite: the matrix of step 3 has
diagonal entries as small as **6.7 × 10⁻¹⁵** in cells touching the axis, the
original preconditioner broke down there, and the solve returned code **−11**
(preconditioner failed) with an infinite velocity field. The function `check_ksp`
now stops the run on any negative code.

**5. Local CFL.** The usual time step estimate divides the smallest cell in the
mesh by the highest velocity anywhere, even if the two are far apart. The code
instead takes, for each cell, that cell's size divided by the velocity in that
cell, and uses the smallest result. On these meshes this gives a time step
**6.1 times larger** for the same real CFL number: 15 steps instead of 92, 10.0 s
instead of 52.7 s on the same test.

**An extra check I did afterwards on the recorded data.** With viscosity and zero
swirl on the walls, circulation `Γ` cannot create a new maximum: its largest
value can only stay the same or go down. On the fine mesh the maximum of `|Γ|`
never increases in any of the 550 intervals (5.749 → 5.168). On the coarse mesh
it wobbles by at most 0.1 % early on, then jumps from 5.28 at `t = 0.2631` to
10.17 at `t = 0.267`. The jump starts exactly where the energy check marks the
last reliable time (`t = 0.2631`), so the two independent checks agree. This one
is not yet built into the code as an automatic stop.

---

## The convergence study

**The idea.** I ran the same case on three meshes. Each is twice as fine as the
previous one in every direction (all lengths divided by 2). If a quantity
changes less and less from coarse to medium to fine, the error can be estimated
and the result can be trusted. The actual ratios measured on the delivered meshes
are **1.985** and **1.969**, close to the target 2. The analysis uses these
measured ratios, not the target.

Common settings: `ν = 0.001`, `T = 0.55`, IPCS, 8 MPI processes on AWS. Full
details in [`docs/convergence.md`](docs/convergence.md).

| level | cells | velocity DOFs | `h` at poles | outcome | reliable up to `t =` | wall time |
| :-- | --: | --: | --: | :-- | --: | --: |
| coarse | 47,134 | 287,205 | 2.137e-3 | stopped by energy check | **0.2631** | 5.1 min |
| medium | 188,462 | 1,139,571 | 1.077e-3 | reached `T` | 0.5500 | 1.6 h |
| fine | 752,803 | 4,534,410 | 5.467e-4 | reached `T` | 0.5500 | 14.0 h |

<table>
<tr>
<td width="50%"><img src="results/convergence_aws/coarse/vortex_blowup_coarse.gif" alt="coarse grid"></td>
<td width="50%"><img src="results/convergence_aws/fine/vortex_blowup_fine.gif" alt="fine grid"></td>
</tr>
<tr>
<td><b>Coarse mesh.</b> The energy check stops it at <code>t = 0.268</code>.
The last reliable state is at <code>t = 0.2631</code>.</td>
<td><b>Fine mesh.</b> Runs the full time with energy going down the whole way.
The energy check never fires.</td>
</tr>
</table>

![Kinetic energy decay and the energy check](media/convergence/energy-decay.png)

The right panel shows what the energy check looks at. The coarse mesh crosses
zero at `t ≈ 0.263`. Medium and fine never do.

**Error bars.** A three-mesh error estimate needs all three meshes to be valid,
so it exists only up to the coarse mesh's reliable time, `t ≤ 0.2631`:

| quantity | fine value | observed `p` | GCI | monotone | verdict |
| :-- | --: | --: | --: | --: | :-- |
| kinetic energy | 152.34 | 1.39 | **0.16 %** | **100 %** | converged |
| max \|Γ\| | 5.4036 | 1.00 | 0.86 % | 92 % | converged |
| enstrophy | 91,232 | 1.82 | 0.40 % | **52 %** | **not converged** |
| max \|ω\| | 2,156.5 | 1.15 | 6.93 % | 75 % | **not converged** |
| BKM integral | 636.18 | 0.98 | 8.85 % | 92 % | low order, use with care |

How to read it:

- A monotone fraction of 52 % means the three meshes move in a consistent
  direction only about half the time. For enstrophy and maximum vorticity, the
  meshes are not fine enough yet, and the GCI and extrapolated values for them
  should not be trusted.
- The two quantities that converge (energy and circulation) are *integrals* over
  the whole flow. The two that do not are *peak values* at the sharpest point of
  the flow. The next section explains why.
- For the quantities that converge, the observed orders are between 1.0 and
  1.4. This is what I would expect from the time stepping: the scheme is first order in time, and the time step is tied to the
  cell size by the CFL condition, so the time error shrinks only in proportion
  to `h`. That hides the higher accuracy of the P2 elements in space. Measuring
  the spatial accuracy on its own would need a smaller time step or a second
  order time scheme.

### The mesh was refined in the wrong place

The mesh is 15 times finer at the poles than at the middle, because I expected
the interesting physics to happen at the poles. It did not.

![Where the vorticity maximum actually lives](media/convergence/vorticity-location.png)

Across the 551 recorded samples on the fine mesh, the point of maximum vorticity
is:

| region | share of samples |
| :-- | --: |
| on the wall (`r/f(z) > 0.95`) | **90.2 %** |
| near the axis (`r < 0.05`) | 1.1 % |
| near the poles (`\|z\| > 1.8`) | **0.0 %** |

Its median height is `z = 0.617`. The maximum vorticity comes from the boundary
layer on the wall around mid-height, and that is where the mesh is coarsest:

![Resolution against the flow it has to carry](media/convergence/resolution-vs-flow.png)

The boundary layer is about `δ ≈ 0.0066` thick. Compared with the cell size at
the equator, it is only **0.2, 0.4 and 0.9 cells** thick on the three meshes:
less than one cell, even on the fine mesh. P2 elements help a little because they can bend
inside a cell, but a boundary layer normally needs around ten points across it.
This is why maximum vorticity and enstrophy do not converge, while energy and
circulation, which depend mostly on the bulk of the flow, do.

On the fine mesh, 12,162 cells went to the pole regions, where the flow does
almost nothing. The monotone column in the table above (52 %) already showed a
problem before I looked at where the vorticity was.

### The vortex never reached the pole

This is the most important physical finding, and it limits what the study can
say.

The setup assumed the ring would be carried into the narrowing tip, where
compression would amplify its swirl. The recorded data show this did not happen:

- **The velocity never grew.** On the fine mesh `max|u|` stays between 18.6 and
  30.0 for the whole run, starting from 25.2. If the ring, with `Γ ≈ 5`, had been
  pushed to `r = 0.05`, its swirl alone would be `u_θ = Γ/r ≈ 100`.
- **The circulation only decreased**, from 5.749 to 5.168, as viscosity slowly
  wore it down.
- **The vorticity maximum never went beyond `|z| = 0.87`**, while the poles are at
  `|z| = 2`.

The likely cause is the background flow. Its streamfunction
`ψ_jet = J·r²·(f² − r²)²` is zero both on the axis and on the wall. That makes it
a closed loop: fluid goes up along the axis and comes back down along the wall.
On the axis its upward speed is `2J·f(z)⁴`, which drops very quickly with height:

| height `z` | 0 | 0.5 | 1.0 | 1.5 |
| :-- | --: | --: | --: | --: |
| upward speed on the axis | 20 | 8.8 | 0.68 | 0.005 |

So the background flow has practically stopped well before the pole, and nothing
carries the ring into the tip. I worked this out from the formula and the
recorded samples; I have not tracked the ring through the saved 3D fields.

In short: the null result is about a vortex ring that stayed in the middle of
the domain. It is not a test of a vortex being compressed in a cone.

### Medium versus fine over the full run

Medium and fine both reached `T = 0.55`, so they can be compared over the whole
run, without error bars.

![Medium versus fine](media/convergence/medium-vs-fine-drift.png)

- **Integral quantities agree.** Kinetic energy differs by 0.62 % at
  `t = 0.2631` and the difference grows slowly to 4.66 % at `T = 0.55`.
- **Peak values do not.** Maximum vorticity differs by up to about 200 %.
- **The BKM integral** ends at 1850.3 (medium) and 1782.3 (fine), 3.8 % apart,
  but the difference reaches 12.3 % during the run. It is built from the
  maximum vorticity, which is not converged, so this agreement is weak.
- **Divergence.** The pointwise error in incompressibility falls from 0.42 to
  0.12 to 0.032 across the three meshes, close to second order. The weak
  version, which is what IPCS actually enforces, is `1.6 × 10⁻⁶` on the fine
  mesh.

---

## What the results mean

**What the numerics show.** On the medium and fine meshes the simulation stays
consistent with the physics for the whole run: energy goes down in every
interval, circulation never makes a new maximum on the fine mesh, and the two
meshes agree on kinetic energy and circulation to within 5 %. On the coarse mesh the checks caught the
point where the mesh stopped being fine enough (`t = 0.2631`) and stopped the
run. So the checking apparatus works, and it would have caught the kind of
failure that fooled Runs 1 and 2.

**What they do not show.** They do not show that "there is no singularity",
for two reasons:

1. **Energy going down is not evidence against blow-up.** Leray (1934) and
   Hopf (1951) proved that Navier–Stokes always has solutions in a weaker,
   averaged sense, and that their energy never grows. The open question is whether they stay smooth.
   A solution can have decreasing energy and still blow up. So the energy
   decrease proves that the *computation* is sound, not that the *solution*
   stays smooth.
2. **The BKM integral staying finite is weak evidence here.** Any finite
   computation over a finite time gives a finite number. What would matter is
   the trend, for example the maximum vorticity growing like `1/(T* − t)`. And
   here the maximum vorticity sits in the unresolved wall boundary layer, so the
   BKM integral mostly measures the boundary layer, not any possible singular
   point inside the flow.

**The honest statement:**

> A checked solver, run on a generic starting flow at `Re ≈ 2.3 × 10⁴` on a
> fixed mesh, shows no singular behaviour up to `T = 0.55`. But the vortex ring
> never reached the narrowing part of the domain, so the compression mechanism
> the setup was designed around was not tested. This is also the regime where
> the literature expects no singularity. What the study does show is that the
> checks detect a loss of numerical validity when it happens: they did, on the
> coarse mesh, at `t = 0.2631`.

**A secondary finding: a resolution threshold.** From the two small validation
meshes I fitted a rule `t_breakdown ~ h^−0.35` for when a mesh stops being fine
enough. It predicted the coarse production mesh within 8 %, and it was completely
wrong about the medium mesh, which never broke down. So the change is not
gradual: the threshold lies between `h_pole = 2.1e-3` and `1.1e-3`. Three meshes
are thin evidence for a threshold, and I would not claim more than that.

---

## Limitations

- **The vortex never reached the pole.** The compression effect was not tested.
- **The domain is a 6.07° cone, not a cusp.**
- **The wall boundary layer is not resolved** (less than one cell across, even
  on the fine mesh). That is where the maximum vorticity is, so maximum
  vorticity and enstrophy are not converged.
- **The time stepping is first order**, which limits the accuracy the
  convergence study can show.
- **Axisymmetry** rules out any instability that breaks the rotational symmetry.
- **One set of parameters:** one viscosity, one starting flow, one geometry.
- **Error bars exist only up to `t = 0.2631`**, because the coarse mesh stopped
  there.
- **The code's order of accuracy has not been verified independently.** There is
  no test with the Method of Manufactured Solutions (a test where you choose an
  exact solution, add the forcing needed to make it solve the equations, and
  check that the error shrinks at the expected rate).
- **A simulation is not a proof**, and this one does not bear on the Clay
  problem.

## Next steps

In order of value:

1. **Change the starting flow so the vortex actually reaches the pole.** This can
   be checked cheaply on a coarse mesh by recording where the swirl is largest
   over time, before spending any cloud time.
2. **Refine the mesh where the flow is:** the wall boundary layer, not the poles.
3. **Use a second-order time scheme**, so the convergence study can measure the
   accuracy in space.
4. **Method of Manufactured Solutions**, to verify the code's order of accuracy.
   It runs locally in minutes.
5. **A real cusp** (`α > 1`), if the geometric question is revisited.
6. **A range of viscosities.**

---

## Reproducing it

Create the environment:

```bash
conda env create -f environment.yml && conda activate fenicsx-env
```

Run the 50 regression tests (about 40 seconds):

```bash
python -m pytest tests -q
```

Some of the tests that check the defects directly:
`test_polar_refinement_is_actually_delivered` fails by 150× on the old mesh
recipe; `test_vorticity_matches_the_analytic_initial_condition` checks 352.59
against 351.64; `test_the_r2_formula_reproduces_the_bug` shows the old vorticity
formula really does blow up on round-off; `test_energy_guard_stops_an_unphysical_run`
checks the energy check.

Run a small three-mesh study locally:

```bash
python scripts/run_convergence.py --config configs/local_validation.json
```

Recompute the convergence analysis and the figures from the recorded AWS data
(these two need only numpy and matplotlib):

```bash
python scripts/analyze_convergence.py --study results/convergence_aws
```

```bash
python scripts/make_figures.py
```

Running on AWS, with the benchmark step and cost control, is described in
[`deploy/`](deploy/).

Meshes are not stored in the repository. They are regenerated: `src/mesh.py`
always produces the same mesh and checks it. The large HDF5 velocity files are
not stored either, because of their size. The recorded diagnostics that every
number in this README comes from are in [`results/`](results/).

## Further documentation

- [`docs/findings.md`](docs/findings.md): every defect in Runs 1 and 2, with measurements.
- [`docs/methodology.md`](docs/methodology.md): what the rewrite does and why.
- [`docs/convergence.md`](docs/convergence.md): the full convergence analysis.
- [`docs/how-this-was-built.md`](docs/how-this-was-built.md): the AI-assisted workflow in detail.

---

## References

**The Hou–Chen research programme**

1. G. Luo, T. Y. Hou, *Potentially singular solutions of the 3D axisymmetric Euler equations*, PNAS **111**(36), 12968–12973 (2014). [doi:10.1073/pnas.1405238111](https://doi.org/10.1073/pnas.1405238111) · [arXiv:1310.0497](https://arxiv.org/abs/1310.0497)
2. G. Luo, T. Y. Hou, *Toward the finite-time blowup of the 3D axisymmetric Euler equations: a numerical investigation*, Multiscale Model. Simul. **12**(4), 1722–1776 (2014).
3. T. Y. Hou, *Potentially Singular Behavior of the 3D Navier–Stokes Equations*, Found. Comput. Math. **23**, 2251–2299 (2023). [doi:10.1007/s10208-022-09578-4](https://doi.org/10.1007/s10208-022-09578-4) · [arXiv:2107.06509](https://arxiv.org/abs/2107.06509)
4. T. Y. Hou, *Potential singularity of the 3D Euler equations in the interior domain*, Found. Comput. Math. (2022). [arXiv:2107.05870](https://arxiv.org/abs/2107.05870)
5. T. Y. Hou, D. Huang, *Potential Singularity Formation of Incompressible Axisymmetric Euler Equations with Degenerate Viscosity Coefficients*, Multiscale Model. Simul. **21**(1), 218–268 (2023). [doi:10.1137/22M1470906](https://doi.org/10.1137/22M1470906) · [arXiv:2102.06663](https://arxiv.org/abs/2102.06663)
6. J. Chen, T. Y. Hou, *Finite Time Blowup of 2D Boussinesq and 3D Euler Equations with C^{1,α} Velocity and Boundary*, Comm. Math. Phys. **383**(3), 1559–1667 (2021). [arXiv:1910.00173](https://arxiv.org/abs/1910.00173). See also the published [correction](https://doi.org/10.1007/s00220-022-04548-x).
7. J. Chen, T. Y. Hou, *Stable nearly self-similar blowup of the 2D Boussinesq and 3D Euler equations with smooth data*, I: Analysis [arXiv:2210.07191](https://arxiv.org/abs/2210.07191); II: Rigorous Numerics, Multiscale Model. Simul. **23**(1), 25–130 (2025). [arXiv:2305.05660](https://arxiv.org/abs/2305.05660)
8. J. Chen, T. Y. Hou, *Singularity formation in 3D Euler equations with smooth initial data and boundary*, PNAS **122**(27), e2500940122 (2025). [doi:10.1073/pnas.2500940122](https://doi.org/10.1073/pnas.2500940122)

**Theory and criteria**

9. J. T. Beale, T. Kato, A. Majda, *Remarks on the breakdown of smooth solutions for the 3-D Euler equations*, Comm. Math. Phys. **94**(1), 61–66 (1984).
10. J. Leray, *Sur le mouvement d'un liquide visqueux emplissant l'espace*, Acta Math. **63**, 193–248 (1934).
11. E. Hopf, *Über die Anfangswertaufgabe für die hydrodynamischen Grundgleichungen*, Math. Nachr. **4**, 213–231 (1951).
12. T. M. Elgindi, *Finite-time singularity formation for C^{1,α} solutions to the incompressible Euler equations on ℝ³*, Annals of Math. **194**(3), 647–727 (2021). [doi:10.4007/annals.2021.194.3.2](https://doi.org/10.4007/annals.2021.194.3.2)
13. C. L. Fefferman, *Existence and smoothness of the Navier–Stokes equation*, Clay Mathematics Institute Millennium Problem description.

**Background reading**

14. D. Barkley, *A fluid mechanic's analysis of the teacup singularity*, Proc. R. Soc. A **476**(2240), 20200348 (2020). [doi:10.1098/rspa.2020.0348](https://doi.org/10.1098/rspa.2020.0348). The most accessible explanation of the Luo–Hou result.
15. Y. Wang, C.-Y. Lai, J. Gómez-Serrano, T. Buckmaster, *Asymptotic Self-Similar Blow-Up Profile for Three-Dimensional Axisymmetric Euler Equations Using Neural Networks*, Phys. Rev. Lett. **130**, 244002 (2023). [doi:10.1103/PhysRevLett.130.244002](https://doi.org/10.1103/PhysRevLett.130.244002)

**Verification and validation**

16. P. J. Roache, *Verification and Validation in Computational Science and Engineering*, Hermosa (1998).
17. ASME V&V 20-2009, *Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer*. The procedure for unequal refinement ratios used here.
18. W. L. Oberkampf, C. J. Roy, *Verification and Validation in Scientific Computing*, Cambridge (2010).

---

## Layout

```text
├── src/         solver, mesh generator, initial flow, diagnostics
├── scripts/     convergence study driver, error analysis, plots, figures
├── configs/     settings for the local validation study and the AWS study
├── tests/       50 regression tests, one per defect found
├── deploy/      AWS setup, launch and cost control
├── results/     recorded diagnostics and run details for every mesh level
├── docs/        defects, methodology, convergence, AI workflow
├── media/       figures and animations used above
└── archive/     Runs 1 and 2, unchanged, kept as evidence
```

Licensed under the terms in [LICENSE](LICENSE).
