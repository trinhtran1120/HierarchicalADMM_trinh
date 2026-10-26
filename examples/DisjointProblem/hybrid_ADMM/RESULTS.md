# Hybrid ADMM comparison results

Saved results snapshot: October 25, 2026. The default comparison uses its saved CSV exports. The all-16-configurations comparison was rerun on October 25, 2026; its table and linked files reflect that run.

## Findings

All three methods converged on the same fixed disjoint problem. The centralized reference objective is **22436.245603277854**.

The default hybrid uses global fADMM with M1=hADMM, M2=fADMM, and M3=hADMM (`F:H/F/H`). It takes **55 outer iterations**, while standalone hADMM takes **86 whole-tree sweeps**. These iteration units represent different work. Including its inner solves, the hybrid uses **10.49 times** as many node updates and **6.94 times** as many transmitted scalars as standalone hADMM. Thus, the 36.0% smaller outer iteration count does not indicate lower computational or communication cost for this setup.

Both distributed solutions agree closely with the centralized reference: maximum coordinate errors are **5.953745e-07** for hADMM and **2.607303e-07** for the default hybrid. Their feasible reconstructed objective gaps are below **1e-14**. In the recorded sequential run, the centralized solver and standalone hADMM are faster than the default hybrid. Runtime values are individual measurements, rather than averages over repeated trials.

## Problem and solver settings

- Fixed Candidate B tree: **16 physical nodes, 15 physical edges, five levels**.
- Three modules with **4, 7, and 4 nodes**; the global root is outside the modules.
- Scalar leaves and nine root coordinates; parent decisions concatenate their physical children in increasing index order.
- Node IDs are assigned in lexicographic `(depth, index)` order.
- Each node uses `f_v(x) = 0.5 * (sum((x_i - ID_v)^2))^2 + 2 * norm(x, 1)`.
- Distributed primal variables start at one and scaled duals at zero. Central decision starts are one.
- Outer and inner penalties: `rho = inner_rho = 1000`.
- Outer stopping tolerance: `1e-4`; inner tolerance: `1e-6`.
- Outer and inner iteration budgets: `1000` each.
- Centralized Ipopt tolerance: `1e-10`.
- Module decisions and internal duals are warm-started between hybrid outer iterations.

## Accuracy and feasibility

| Metric | Centralized | hADMM | Hybrid F:H/F/H |
| --- | --- | --- | --- |
| Converged | true | true | true |
| Objective at current node copies | 22436.245603277854 | 22436.245591916100 | 22436.245599862170 |
| Relative objective gap | 0.000000e+00 | 5.064017e-10 | 1.522395e-10 |
| Feasible reconstructed objective | 22436.245603277854 | 22436.245603277941 | 22436.245603277879 |
| Feasible relative objective gap | 0.000000e+00 | 3.891538e-15 | 1.135032e-15 |
| Maximum solution error | 0.000000e+00 | 5.953745e-07 | 2.607303e-07 |
| Relative solution error (stacked 2-norm) | 0.000000e+00 | 4.432412e-08 | 1.808331e-08 |
| Physical constraint residual | 0.000000e+00 | 8.081505e-09 | 2.128280e-08 |
| ADMM primal residual | N/A | 8.081505e-09 | 2.717241e-08 |
| ADMM dual residual | N/A | 8.722596e-05 | 8.388735e-05 |

The objective gap is `abs(J - Jcentral) / max(1, abs(Jcentral))`. The current physical node copies can be slightly inconsistent, so their raw objective can lie below the centralized feasible objective. The feasible objective reconstructs every node decision from the root vector before evaluating the original cost.

Maximum solution error is the largest absolute coordinate error over all physical nodes. Relative solution error uses the stacked physical decisions: `norm(x - xcentral, 2) / max(1, norm(xcentral, 2))`. Physical residual is the largest infinity norm of the original physical-edge constraint residuals. The two distributed methods report their own ADMM primal/dual stopping residuals; their splittings differ.

## Work, communication, and runtime

**Selection rule:** “Best” here means the converged hybrid with the **lowest total transmitted scalars** among all 16 cases. `F:F/H/H` (global fADMM; M1=fADMM, M2=hADMM, M3=hADMM) achieves **37,194 scalars**, tied with `H:F/H/H`. The global fADMM choice is shown as the representative of that tie. It also ties for the lowest total node updates and maximum per-node communication among the tested hybrids.

All three columns below use [results/all/summary.csv](results/all/summary.csv) from the same comparison batch, including its baseline timings. The accuracy and root-solution sections retain the default `F:H/F/H` comparison.

| Metric | Centralized | hADMM | Best hybrid F:F/H/H |
| --- | --- | --- | --- |
| Iterations | 8 | 86 | 55 |
| Iteration type | Ipopt barrier iteration | Whole-tree sweep | Outer ADMM iteration |
| Total inner sweeps across modules | N/A | N/A | 2,711 |
| Total node proximal updates | N/A | 1,376 | 13,920 |
| Maximum updates at one node | N/A | 86 | 1,007 |
| Total transmitted scalars | 0 | 5,934 | 37,194 |
| Maximum scalars transmitted by one node | 0 | 1,118 | 5,035 |
| Wall time after warm-up (seconds) | 0.002800 | 0.013852 | 0.441572 |

The fastest recorded hybrid is `F:H/H/H` at **0.147545 s**, using **39,726 transmitted scalars**. Thus, the selected hybrid minimizes communication rather than runtime; this table reports all metrics from one selected configuration instead of combining the best value from different cases.

Node updates count node proximal solves, including all hybrid inner sweeps; PANOC's internal numerical steps are excluded. Centralized Ipopt iterations, whole-tree hADMM sweeps, and hybrid outer iterations must be interpreted using their listed units.

Communications count scalar payloads on every physical hop, including relays and stopping reports. Root-to-M2 queries and replies follow `(1,1) <-> (2,2) <-> (3,4)`. Centralized optimization communication is zero under the assumption that problem data already resides at the solver. Static setup, model upload, solution distribution, message headers, and offline diagnostics are excluded.

Wall time includes model/solver setup and diagnostics, after a full warm-up of each method. Plotting, writing output, warm-up, and the independent reformulation check are excluded. Execution is a sequential simulation, so these timings do not estimate parallel network wall time.

## Root solution

The root vector orders leaf decisions by recursive physical child order. Values below are rounded to 12 decimal places; the CSV retains full precision and every physical node's coordinates.

| Coordinate | Centralized | hADMM | Hybrid F:H/F/H |
| --- | --- | --- | --- |
| 1 | 1.173236469987 | 1.173235874613 | 1.173236405800 |
| 2 | 1.380400439356 | 1.380400025877 | 1.380400404647 |
| 3 | 1.664571717249 | 1.664571550278 | 1.664571712438 |
| 4 | 5.663760834805 | 5.663760825469 | 5.663760865933 |
| 5 | 6.199250247931 | 6.199250258046 | 6.199250267494 |
| 6 | 3.518320902647 | 3.518320928573 | 3.518320669239 |
| 7 | 6.830511241584 | 6.830511250231 | 6.830511254486 |
| 8 | 4.547359439807 | 4.547359447530 | 4.547359446942 |
| 9 | 5.047166063029 | 5.047166071854 | 5.047166070081 |

Full coordinate results: [solutions.csv](results/solutions.csv).

## Figures

- [Convergence](results/convergence.pdf): feasible objective gap versus iterations and communication, solution error, and physical feasibility.
- [Resource comparison](results/resources.pdf): total communication, maximum node communication, node updates, and runtime.

Curves retain the observed histories without monotone smoothing. Values below `1e-16` are clipped only for logarithmic display. The central reference has no ADMM trajectory.

## All 16 hybrid assignments

`F` means fADMM and `H` means hADMM. The configuration order is `global:M1/M2/M3`. All listed configurations converged. Global fADMM and hADMM implement equivalent updates on this two-level computational star; their matching numerical/communication results do not represent two distinct global algorithms.

| Global:M1/M2/M3 | Outer iterations | Relative objective gap | TotalVals | MaxNodeVals | Node updates | Seconds |
| --- | --- | --- | --- | --- | --- | --- |
| `F:F/F/F` | 55 | 1.406835e-10 | 41,098 | 7,945 | 14,959 | 1.571190 |
| `H:F/F/F` | 55 | 1.406835e-10 | 41,098 | 7,945 | 14,959 | 1.575720 |
| `F:H/F/F` | 55 | 1.406835e-10 | 43,630 | 7,945 | 14,959 | 1.275546 |
| `H:H/F/F` | 55 | 1.406835e-10 | 43,630 | 7,945 | 14,959 | 1.291792 |
| `F:F/H/F` | 55 | 1.333573e-10 | 39,638 | 5,952 | 14,448 | 0.816784 |
| `H:F/H/F` | 55 | 1.333573e-10 | 39,638 | 5,952 | 14,448 | 0.808079 |
| `F:H/H/F` | 55 | 1.333572e-10 | 42,170 | 5,952 | 14,448 | 0.538675 |
| `H:H/H/F` | 55 | 1.333572e-10 | 42,170 | 5,952 | 14,448 | 0.507076 |
| `F:F/F/H` | 55 | 1.522395e-10 | 38,654 | 7,945 | 14,431 | 1.207996 |
| `H:F/F/H` | 55 | 1.522395e-10 | 38,654 | 7,945 | 14,431 | 1.212897 |
| `F:H/F/H` | 55 | 1.522395e-10 | 41,186 | 7,945 | 14,431 | 0.920849 |
| `H:H/F/H` | 55 | 1.522395e-10 | 41,186 | 7,945 | 14,431 | 0.917281 |
| `F:F/H/H` | 55 | 1.449233e-10 | 37,194 | 5,035 | 13,920 | 0.441572 |
| `H:F/H/H` | 55 | 1.449233e-10 | 37,194 | 5,035 | 13,920 | 0.443501 |
| `F:H/H/H` | 55 | 1.449230e-10 | 39,726 | 5,035 | 13,920 | 0.147545 |
| `H:H/H/H` | 55 | 1.449230e-10 | 39,726 | 5,035 | 13,920 | 0.158224 |

These are from the separate `--all` run, so their runtimes can differ from the default run above. Full results: [all/summary.csv](results/all/summary.csv), [all/history.csv](results/all/history.csv), [all/solutions.csv](results/all/solutions.csv), [all convergence](results/all/convergence.pdf), and [all resources](results/all/resources.pdf).

## Sources, validation, and reproduction

- [Default summary](results/summary.csv), [history](results/history.csv), and [solutions](results/solutions.csv).
- [Implementation and metric definitions](README.md).
- The previously completed validation suite passed **343 assertions**, including the centralized/reformulated comparison, module-oracle agreement, all 16 configurations, physical relay accounting, and nondefault standalone penalty propagation.

From the repository root:

```sh
# Centralized, standalone hADMM, and the default mixed hybrid
julia --project=. examples/DisjointProblem/hybrid_ADMM/main.jl

# Centralized and hADMM baselines with all 16 hybrid assignments
julia --project=. examples/DisjointProblem/hybrid_ADMM/main.jl --all

# Validation suite
julia --project=. examples/DisjointProblem/hybrid_ADMM/runtests.jl
```

This Markdown report is a snapshot. Rerunning these commands refreshes the generated CSV/PDF files; it does not automatically refresh this report.
