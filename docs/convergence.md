# Convergence — what converged, what did not, and why

The production study in detail. Two of the five quantities I tracked converge,
two do not, and one converges only at a low order. This document says which, by
how much, and why. It also explains two limits on what the study can conclude:
the vortex ring never reached the pole, and the time stepping is first order.

A short reminder of the terms used here (all are defined in plain words in the
[main README](../README.md#key-terms)): `h` is the cell size, `p` the observed
order (how fast the error shrinks as `h` shrinks), GCI the error bar on the
fine-grid value, and the monotone fraction the share of samples in which the
three grids move steadily in one direction.

Raw data and per-level provenance are in [`../results/`](../results/). Regenerate
everything here with:

```bash
python scripts/analyze_convergence.py --study results/convergence_aws
```

---

## The grids

Three levels with **all lengths scaled by the same factor**, so a single `h`
identifies each — without that, Richardson extrapolation is not defined.

| level | cells | velocity DOFs | h_pole asked | h_pole achieved | h_equator | outcome | final t | trustworthy to |
| :-- | --: | --: | --: | --: | --: | :-- | --: | --: |
| coarse | 47,134 | 287,205 | 2.000e-3 | 2.137e-3 | 3.206e-2 | `energy_growth` | 0.26804 | **0.26310** |
| medium | 188,462 | 1,139,571 | 1.000e-3 | 1.077e-3 | 1.566e-2 | **`completed`** | 0.55000 | 0.55000 |
| fine | 752,803 | 4,534,410 | 5.000e-4 | 5.467e-4 | 7.777e-3 | **`completed`** | 0.55000 | 0.55000 |

Refinement ratios measured **from the delivered meshes**, not assumed:
`coarse/medium = 1.9848`, `medium/fine = 1.9692`, against a nominal 2.0. The
observed order is solved for these two ratios using the ASME V&V 20 procedure for
a non-uniform ratio rather than forcing a constant one, so the deviation does not
bias the result.

Wall time: 5.1 min, 1.6 h and 14.0 h respectively on 8 MPI ranks of an AWS
`c6i.4xlarge` — 15.7 h of solver time in total.

## The three-grid window

Each level is truncated at its `t_at_energy_min`, not at its last written row: the
samples between the energy minimum and the guard's stop are already contaminated.
The window in which all three levels are valid is therefore set by the coarsest:

> **`t ∈ [0.00024, 0.26310]`, 200 sample points.**

Formal error bars exist only there. Medium and fine remain comparable beyond it,
but as a two-grid pair without GCI.

## Results

`p` is the observed order; `GCI` the fine-grid Grid Convergence Index (Roache,
safety factor 1.25); `monotone` the fraction of sample times at which the three
grids form a monotone sequence — the precondition for the extrapolation to mean
anything at all.

| quantity | fine value | Richardson h→0 | `p` | GCI | monotone | verdict |
| :-- | --: | --: | --: | --: | --: | :-- |
| kinetic energy | 152.337 | 152.491 | 1.39 | **0.16 %** | **100 %** | converged |
| max \|Γ\| = r·u_θ | 5.40361 | 5.43324 | 1.00 | 0.86 % | 92 % | converged |
| enstrophy | 91,231.6 | 91,372.1 | 1.82 | 0.40 % | **52 %** | **not asymptotic** |
| max \|ω\| | 2,156.47 | 2,155.35 | 1.15 | 6.93 % | 75 % | **not asymptotic** |
| BKM integral | 636.183 | 637.514 | 0.98 | 8.85 % | 92 % | order low |

![Kinetic energy convergence](../media/convergence/gci-kinetic-energy.png)

**A monotone fraction of 52 % is a coin flip.** For enstrophy the three grids do
not form an ordered sequence more than half the time, which means they are not in
the asymptotic range and the extrapolated value is an artefact, not a better
answer than the fine grid. The same applies, less severely, to peak vorticity at
75 %.

The two quantities that do converge, kinetic energy (100 % monotone, 0.16 % GCI)
and circulation (92 %, 0.86 %), are integrals over the whole flow, so they depend
mostly on the bulk. The two that do not are peak values at the sharpest feature
in the domain. The section after next explains why.

### Why the observed orders are close to 1

The time stepping in step 1 of IPCS is implicit Euler, which is first order in
time: the error it adds shrinks in proportion to `dt`. The time step is set by
the CFL condition, so `dt` shrinks in proportion to `h` (measured ratio 1.98
between grids differing by 2). The time error therefore shrinks like `h`, and it
hides the higher spatial accuracy of the P2 elements. That is consistent with
`p = 1.00` for circulation, `0.98` for the BKM integral and `1.39` for kinetic
energy. Measuring the spatial order on its own would need a time step well below
the CFL limit, or a second-order time scheme (Crank–Nicolson or BDF2).

## Why: the refinement went to the wrong place

The mesh is graded **15:1 toward the poles**, because that is where the geometric
feature of interest is. The flow disagreed.

![Where the vorticity maximum actually lives](../media/convergence/vorticity-location.png)

Across all 551 recorded samples on the fine grid, the location of `max|ω|`:

| region | share of samples |
| :-- | --: |
| on the wall (`r/f(z) > 0.95`) | **90.2 %** |
| near the axis (`r < 0.05`) | 1.1 % |
| near the poles (`\|z\| > 1.8`) | **0.0 %** |

Median height of the maximum: `z = 0.617`. Range: `−0.558` to `+0.862`. **It is
never at the poles, in any sample.** The extreme vorticity is a wall
boundary-layer feature at mid-latitude.

And the grid is at its coarsest precisely there:

![Resolution against the flow](../media/convergence/resolution-vs-flow.png)

| level | h at poles | h at equator wall | ratio | cells across `δ ≈ 0.0066` |
| :-- | --: | --: | --: | --: |
| coarse | 2.137e-3 | 3.206e-2 | 15.0 | **0.2** |
| medium | 1.077e-3 | 1.566e-2 | 14.6 | **0.4** |
| fine | 5.467e-4 | 7.777e-3 | 14.2 | **0.9** |

With `U ≈ 23` and `ν = 1e-3`, `Re ≈ 2.3 × 10⁴` and the laminar wall layer is
`δ ~ L/√Re ≈ 0.0066`. Even the fine grid puts **under one element** across it. P2
elements give some sub-cell resolution, so this is not catastrophic, but by any
standard criterion — where one wants of order ten points across a boundary layer —
the layer is unresolved.

**12,162 cells went to the pole region on the fine grid.** Nothing happens there,
and the next section shows why.

## The vortex never reached the pole

This is the most important limit on the physical conclusion. The setup assumed
the ring would be carried into the narrowing tip, where compression would amplify
its swirl (`u_θ = Γ/r` grows as `r` shrinks, because circulation `Γ` is nearly
conserved). The recorded data show this did not happen:

| signal on the fine grid | recorded | what compression to `r = 0.05` would give |
| :-- | :-- | :-- |
| `max|u|` | stays between 18.6 and 30.0 (starts at 25.2) | about 100 from the swirl alone (`Γ ≈ 5`, `r = 0.05`) |
| `max|Γ|` | only decreases, 5.749 → 5.168 | unchanged or slowly decreasing |
| height of `max|ω|` | never above `|z| = 0.86` | near `|z| = 2` |

The likely cause is the background flow. Its streamfunction
`ψ_jet = J·r²·(f² − r²)²` is zero on the axis and on the wall, so it is a closed
loop: up along the axis, back down along the wall. On the axis the upward speed
is `u_z = 2J·f(z)⁴`:

| `z` | 0 | 0.5 | 1.0 | 1.5 |
| :-- | --: | --: | --: | --: |
| `u_z` on the axis | 20 | 8.8 | 0.68 | 0.005 |

The background flow has practically stopped well before the pole, so nothing
carries the ring into the tip. This comes from the formula and the recorded
samples; I have not tracked the ring through the saved 3D fields.

So the null result is about a ring that stayed in the middle of the domain. It is
not a test of a vortex compressed in a cone.

This is the mechanism behind the table above: `max|ω|` and enstrophy are boundary
layer quantities measured where the mesh is coarse, so they do not converge;
kinetic energy and circulation are bulk quantities, so they do.

It is worth being explicit that **the GCI table diagnosed this before the location
analysis did.** `monotone = 52 %` is the signal. A reader who takes the
extrapolated enstrophy at face value and skips the monotone column gets a
confident wrong number.

## Two grids over the full interval

Medium and fine both reached `T = 0.55`, so they can be compared directly beyond
the three-grid window — without error bars, but across the whole run.

![Medium versus fine](../media/convergence/medium-vs-fine-drift.png)

| quantity | median difference | at `t = 0.2631` | at `T = 0.55` | max |
| :-- | --: | --: | --: | --: |
| kinetic energy | 0.69 % | 0.62 % | **4.66 %** | 4.66 % |
| max \|Γ\| | 0.68 % | — | 2.90 % | 2.90 % |
| BKM integral | 6.60 % | — | **3.81 %** | 12.3 % |
| enstrophy | 3.23 % | — | 9.08 % | 17.2 % |
| max \|ω\| | 8.70 % | 4.60 % | 24.1 % | **198 %** |

Two readings, both worth stating:

- **The integral quantities agree.** Kinetic energy within 0.62 % at the
  three-grid horizon, BKM within 3.81 % at `T`. The null result rests on these.
- **The pointwise peaks do not.** `max|ω|` differs by up to 198 % (at `t = 0.29`,
  medium 9,622 against fine 3,223). Any claim built on peak vorticity at late
  times would be unsupported by this data.

The kinetic-energy difference grows monotonically with time — 0.05 % → 0.24 % →
0.63 % → 1.78 % → 4.66 % — which is ordinary error accumulation, but it means
"two grids agree over the whole interval" is true to about 5 %, not tightly.

## Energy: the check the study relies on

![Energy decay and the guard](../media/convergence/energy-decay.png)

In a closed unforced domain `dE/dt ≤ 0` exactly, so any rise is numerical by
construction. Counting rising samples:

| level | E(0) | E_min | rising samples | max ΔE |
| :-- | --: | --: | --: | --: |
| coarse | 175.597 | 144.562 | **4 / 267** | +0.4974 |
| medium | 175.680 | 106.883 | **0 / 550** | −0.0297 |
| fine | 175.704 | 112.109 | **0 / 550** | −0.0236 |

Medium and fine are **strictly monotone over the entire run**. The guard never
fired on either. On the coarse grid it fired at `t = 0.2680`, having recorded the
energy minimum at `t = 0.2631`, which becomes that level's reliable horizon and
therefore the three-grid window above.

**What this shows and what it does not.** Decreasing energy shows the
*computation* is consistent with the physics. It is not evidence that the
*solution* stays smooth: Leray (1934) and Hopf (1951) showed that weak solutions
of Navier–Stokes always exist and never gain energy, and whether they can still
become singular is exactly the open question. A blowing-up solution could have
decreasing energy too.

**A second, independent check: circulation.** With viscosity and zero swirl on
the boundary, `Γ` obeys a maximum principle: its largest value can only stay the
same or go down. Counting intervals in which `max|Γ|` rises:

| level | rises | largest rise | comment |
| :-- | --: | --: | :-- |
| coarse | 92 / 267 | +61 % | wobbles below 0.1 % early, then jumps 5.28 → 10.17 from `t = 0.2631` to `0.267` |
| medium | 10 / 550 | +0.03 % | sampling noise |
| fine | **0 / 550** | — | 5.749 → 5.168 |

The coarse jump starts exactly at the energy guard's reliable horizon, so the two
checks agree. This one is not yet coded as an automatic stop.

**The BKM integral** ends at **1782.3** (fine) and **1850.3** (medium), 3.8 %
apart, with differences up to 12.3 % along the way. A singularity at `T*`
requires `∫‖ω‖_∞ dt` to diverge, and it does not. But this is weak evidence here,
for two reasons. Any finite computation over a finite time gives a finite
number; what would matter is the trend, for example `‖ω‖_∞` growing like
`1/(T* − t)`. And `‖ω‖_∞` sits in the unresolved wall boundary layer and is not
converged (up to 198 % apart between medium and fine), so the integral mostly
measures the boundary layer, not a possible singular point inside the flow.

## The divergence constraint

| level | strong `div_u_rel` at `T` | weak residual at `T` |
| :-- | --: | --: |
| coarse | 0.4236 | 2.71e-4 |
| medium | 0.1179 | 3.17e-5 |
| fine | **0.0317** | **1.60e-6** |

The strong residual converges at close to second order (ratios 3.6 and 3.7 for a
mesh ratio of ~2). The weak residual — which is what a projection method actually
enforces — is at machine-adjacent levels. A 3 % pointwise divergence on the fine
grid is normal for IPCS and is reported rather than omitted.

## The resolution threshold

From the two validation grids alone, a fitted scaling `t_breakdown ~ h^−0.35`
predicted breakdown at `t = 0.243` (coarse) and `t = 0.309` (medium). It got the
coarse level right to within 8 % and was **completely wrong about the medium
level**, which never broke down at all.

**The degradation is not a smooth power law.** There is a *resolution threshold*:
below it the grid loses the flow at some point; above it the problem disappears
entirely, at least out to `T = 0.55`. For this study the threshold sits between
`h_pole = 2.1e-3` and `1.1e-3`, now confirmed by the fine level as well.

This is a genuine finding, and it comes with a caveat that should be stated
whenever it is quoted: **three points is thin evidence for a threshold.** What can
be said with confidence is the negative half — extrapolating a breakdown horizon
from two coarse points and then trusting it does not work here.

Practically, this is why `T` is left at 0.55 and the guard decides: it needs no
such extrapolation and bounds the cost automatically.

## What would change these numbers

In descending order of expected value:

1. **Change the starting flow so the vortex actually reaches the pole.** Without
   this, the compression mechanism stays untested and the other items only
   improve error bars. It can be checked cheaply on a coarse mesh by recording
   where the swirl is largest over time, before spending cloud time.
2. **Refine the wall rather than the poles.** This addresses the cause of the
   non-convergence of `max|ω|` and enstrophy. A two-level study costs about 1.7 h
   of the same instance.
3. **Second-order time stepping**, so the observed order reflects the spatial
   accuracy and not the time error.
4. **Method of Manufactured Solutions.** Establishes the formal order of the
   implementation, which is currently unknown, giving the observed `p` values a
   reference to be judged against. Runs locally in minutes.
5. A true cusp profile with `α > 1`, if the geometric question is to be revisited.
6. A parameter sweep in `ν`.

Items 2 to 6 would change the error bars. Only item 1 would change what the study
can say about the physics.
