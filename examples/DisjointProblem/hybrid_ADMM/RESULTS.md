# Hybrid ADMM comparison results

Results snapshot: October 25, 2026. This run compares **centralized, standalone hADMM, standalone fADMM, the provided nested hybrid, and all 16 Candidate B assignments**. Both hybrid methods use inner tolerance **`1e-4`**. The default CSVs contain the four baselines and Candidate B `F:H/F/H`; `results/all/` contains the four baselines and all 16 Candidate B cases. Both exports come from the same warmed-up benchmark run.

## Findings

All **20 cases** in the full benchmark converged. The centralized reference objective is **22436.245603277854**. Standalone hADMM took **86 whole-tree sweeps**, standalone fADMM took **87 flattened iterations**, the nested hybrid took **50 outer iterations**, and every Candidate B assignment took **55 outer iterations**. These iteration units represent different work.

The nested hybrid used **7,878 node proximal updates** and **30,994 transmitted scalars**. The communication-best Candidate B `F:F/H/H` used **8,793 updates** and **24,580 scalars**. Thus, the nested method used fewer node updates but more communication than the best Candidate B case. Standalone hADMM and fADMM still used fewer node updates and scalar transmissions than either hybrid.

The nested hybrid's relative objective gap is **7.524748e-10**, maximum solution error is **6.727496e-07**, and physical constraint residual is **4.126058e-08**. Its recorded wall time after warm-up was **1.480088 s**.

The new `1e-4` inner tolerance is a looser accuracy requirement than `1e-6`. It reduces inner work while retaining convergence in all tested configurations. The before/after results below compare identical configurations. Reported runtimes are individual sequential measurements after warm-up, rather than averages over repeated trials.

## Problem and solver settings

- Fixed Candidate B tree: **16 physical nodes, 15 physical edges, five levels**.
- Candidate B modules M1, M2, and M3 contain **4, 7, and 4 nodes**; the global root is outside the modules.
- The provided nested hybrid uses two complete root-child subtrees containing **4 and 11 nodes**.
- Scalar leaves and nine root coordinates; parent decisions concatenate their physical children in increasing index order.
- Legacy numeric IDs follow lexicographic `(depth, index)` order.
- Node cost: `f_v(x) = 0.5 * (sum((x_i - ID_v)^2))^2 + 2 * norm(x, 1)`.
- Distributed primal variables start at one and scaled duals at zero; central decision starts are one.
- Outer penalty `rho=1000`; hybrid inner penalty `inner_rho=1000`.
- Outer stopping tolerance `1e-4`; hybrid inner stopping tolerance **`1e-4`**.
- Outer and inner iteration budgets: `1000` each.
- Centralized Ipopt solver tolerance: `1e-10`.
- Both hybrid methods warm-start module decisions and internal duals between outer iterations.

## Accuracy and feasibility

This table retains the default hybrid `F:H/F/H` (global fADMM; M1=hADMM, M2=fADMM, M3=hADMM).

| Metric | Centralized | hADMM | fADMM | Nested hybrid F:H/H | Default Candidate B F:H/F/H |
| --- | --- | --- | --- | --- | --- |
| Converged | true | true | true | true | true |
| Objective at current node copies | 22436.245603277854 | 22436.245591916100 | 22436.245601354065 | 22436.245586395144 | 22436.245593972712 |
| Relative objective gap | 0.000000e+00 | 5.064017e-10 | 8.574468e-11 | 7.524748e-10 | 4.147370e-10 |
| Feasible reconstructed objective | 22436.245603277854 | 22436.245603277941 | 22436.245603277945 | 22436.245603277977 | 22436.245603277908 |
| Feasible relative objective gap | 0.000000e+00 | 3.891538e-15 | 4.053685e-15 | 5.513011e-15 | 2.432211e-15 |
| Maximum solution error | 0.000000e+00 | 5.953745e-07 | 6.077959e-07 | 6.727496e-07 | 3.961805e-07 |
| Relative solution error (stacked 2-norm) | 0.000000e+00 | 4.432412e-08 | 4.510647e-08 | 5.492744e-08 | 3.112212e-08 |
| Physical constraint residual | 0.000000e+00 | 8.081505e-09 | 4.852807e-09 | 4.126058e-08 | 2.700873e-08 |
| ADMM primal residual | N/A | 8.081505e-09 | 5.318737e-09 | 4.126058e-08 | 3.201131e-08 |
| ADMM dual residual | N/A | 8.722596e-05 | 8.859194e-05 | 9.390993e-05 | 8.710389e-05 |

Objective gap is `abs(J-Jcentral)/max(1,abs(Jcentral))`, using actual node costs. Current physical node copies can remain slightly inconsistent and their raw objective can lie below the central feasible objective. The feasible objective reconstructs every physical decision from the root vector before evaluating costs.

Maximum solution error is the largest absolute coordinate error over all physical nodes. Relative solution error is the stacked 2-norm error divided by `max(1,norm(xcentral,2))`. Physical residual is the largest infinity norm of the original edge constraints. ADMM primal/dual residuals follow each method's computational splitting and are not identical definitions across methods.

## Work, communication, and runtime

**Selection rule:** “Best” means the converged hybrid with the **lowest total transmitted scalars** among all 16 Candidate B cases. `F:F/H/H` achieves **24,580 scalars**, tied with `H:F/H/H`; global fADMM is shown as the representative. It also ties for the lowest node-update count and maximum per-node communication among the tested Candidate B assignments. All columns below use the same comparison batch.

| Metric | Centralized | hADMM | fADMM | Nested hybrid F:H/H | Best Candidate B F:F/H/H |
| --- | --- | --- | --- | --- | --- |
| Iterations | 8 | 86 | 87 | 50 | 55 |
| Iteration type | Ipopt barrier iteration | Whole-tree sweep | Flattened ADMM iteration | Nested outer iteration | Outer ADMM iteration |
| Total inner sweeps across modules | N/A | N/A | N/A | 1,047 | 1,694 |
| Total node proximal updates | N/A | 1,376 | 1,392 | 7,878 | 8,793 |
| Maximum updates at one node | N/A | 86 | 87 | 527 | 654 |
| Total transmitted scalars | 0 | 5,934 | 9,918 | 30,994 | 24,580 |
| Maximum scalars transmitted by one node | 0 | 1,118 | 3,132 | 4,701 | 3,270 |
| Wall time after warm-up (seconds) | 0.002813 | 0.014485 | 0.837991 | 1.480088 | 0.304428 |

The fastest recorded Candidate B assignment is `F:H/H/H` at **0.103694 s**, with **26,062 scalars**. The selected Candidate B column minimizes communication rather than runtime, and every entry in that column comes from one configuration.

Node updates count node proximal solves, including hybrid inner sweeps; PANOC's internal numerical steps are excluded. Ipopt barrier iterations, whole-tree sweeps, flattened ADMM iterations, and hybrid outer iterations represent different work.

Communication counts scalar payloads on every physical hop, including relays and stopping reports. Root-to-M2 hybrid traffic uses `(1,1) <-> (2,2) <-> (3,4)`. Standalone fADMM coordinates every non-root node from the global root and counts its physical relay routes. Its coordinator can calculate stopping residuals from received decisions, so no extra stopping payload is required. hADMM counts one subtree stopping scalar per child-to-parent reply.

Centralized communication is zero during optimization with data already resident. Static setup, model upload, solution distribution, message headers, and offline diagnostics are excluded for all methods. Runtime includes setup and diagnostics, after a full warm-up of each method; plotting, output, warm-up, and the independent reformulation check are excluded. Execution is sequential and does not estimate parallel network wall time.

## Nested hybrid baseline and partition

The benchmark directly includes [`hybrid_nested_admm/hybrid_ADMM.jl`](../../hybrid_nested_admm/hybrid_ADMM.jl) and calls its `hybridADMM` implementation. The random-tree experiment in that folder's `main.jl` is not used: this baseline runs on the same fixed 16-node tree, numerical node costs, initialization, penalties, tolerances, budgets, and warm-up policy as the other benchmark methods.

`Nested hybrid F:H/H` means global fADMM with local hADMM in two complete child subtrees:

| Nested module root | Physical nodes | Relation to Candidate B |
| --- | --- | --- |
| `(2,1)` | 4 | Same nodes as M1 |
| `(2,2)` | 11 | M2 and M3 combined |

The physical problem is identical; the partition and coordination schedule differ. The source sends the outer target again on every module inner sweep, and these transmissions are retained and counted. The adapter additionally counts one scalar subtree stopping report per internal return edge and three report scalars (inner residual, iterations, convergence flag) per module-to-global response. Forwarding does not count as an optimization update.

The nested method completed **1,047 inner sweeps across its two modules**, with **527** maximum node updates and **4,701** maximum transmitted scalars at one node. Its decisions, histories, resource metrics, and curves appear in the same CSV and PDF exports as the other methods.

## Effect of changing the hybrid inner tolerance

The previous `1e-6` summaries were preserved before rerunning. The default and communication-best Candidate B configurations compare as follows (the previous archive predates the nested baseline):

| Hybrid | Metric | Previous inner_tol=1e-6 | New inner_tol=1e-4 | Reduction |
| --- | --- | --- | --- | --- |
| `F:H/F/H` | Total transmitted scalars | 41,186 | 26,742 | 35.07% |
| `F:H/F/H` | Total node updates | 14,431 | 9,031 | 37.42% |
| `F:H/F/H` | Inner sweeps | 2,784 | 1,728 | 37.93% |
| `F:F/H/H` | Total transmitted scalars | 37,194 | 24,580 | 33.91% |
| `F:F/H/H` | Total node updates | 13,920 | 8,793 | 36.83% |
| `F:F/H/H` | Inner sweeps | 2,711 | 1,694 | 37.51% |

Both listed configurations take 55 outer iterations at either tolerance. The communication-best hybrid's maximum solution error changes from **2.607544e-07** to **3.966154e-07**; its raw relative objective gap changes from **1.449233e-10** to **1.181538e-09**. Both remain within the benchmark's accuracy checks. Its current feasible objective gap is **2.594358e-15**.

Previous summaries: [default](results/inner_tol_1e-6/summary.csv) and [all cases](results/inner_tol_1e-6/all/summary.csv). There is no previous standalone fADMM measurement in that archive.

## Root solution

Coordinates use recursive physical child order and are rounded to 12 decimal places. The hybrid here is the default `F:H/F/H`; [all/solutions.csv](results/all/solutions.csv) also contains the communication-best hybrid and every other assignment.

| Coordinate | Centralized | hADMM | fADMM | Nested hybrid F:H/H | Default Candidate B F:H/F/H |
| --- | --- | --- | --- | --- | --- |
| 1 | 1.173236469987 | 1.173235874613 | 1.173235866439 | 1.173235842917 | 1.173236228569 |
| 2 | 1.380400439356 | 1.380400025877 | 1.380400025281 | 1.380399912897 | 1.380400286991 |
| 3 | 1.664571717249 | 1.664571550278 | 1.664571550762 | 1.664571420781 | 1.664571676945 |
| 4 | 5.663760834805 | 5.663760825469 | 5.663760829075 | 5.663760856764 | 5.663760890385 |
| 5 | 6.199250247931 | 6.199250258046 | 6.199250259041 | 6.199250266300 | 6.199250281812 |
| 6 | 3.518320902647 | 3.518320928573 | 3.518320918154 | 3.518320921111 | 3.518320560190 |
| 7 | 6.830511241584 | 6.830511250231 | 6.830511251369 | 6.830511257023 | 6.830511263200 |
| 8 | 4.547359439807 | 4.547359447530 | 4.547359447509 | 4.547359454784 | 4.547359452552 |
| 9 | 5.047166063029 | 5.047166071854 | 5.047166071775 | 5.047166075170 | 5.047166074737 |

Full default-coordinate results: [solutions.csv](results/solutions.csv).

## Figures

- [Default convergence](results/convergence.pdf): feasible gap versus iterations/communication, solution error, and physical feasibility, including standalone fADMM and the nested hybrid.
- [Default resource comparison](results/resources.pdf): total communication, peak node communication, node updates, and runtime.
- [All-case convergence](results/all/convergence.pdf) and [resources](results/all/resources.pdf).

Curves retain observed histories without monotone smoothing. Values below `1e-16` are clipped only for logarithmic display. The centralized reference has no ADMM trajectory.

## All 16 hybrid assignments

`F` means fADMM; `H` means hADMM. Configuration order is `global:M1/M2/M3`. Every listed case converged with inner tolerance `1e-4`. Global fADMM and hADMM share equivalent updates on this two-level computational star, explaining matching numerical/communication results for each module assignment.

| Global:M1/M2/M3 | Outer iterations | Relative objective gap | TotalVals | MaxNodeVals | Node updates | Seconds |
| --- | --- | --- | --- | --- | --- | --- |
| `F:F/F/F` | 55 | 3.668864e-10 | 26,418 | 5,201 | 9,235 | 1.131690 |
| `H:F/F/F` | 55 | 3.668864e-10 | 26,418 | 5,201 | 9,235 | 1.119157 |
| `F:H/F/F` | 55 | 3.668862e-10 | 27,900 | 5,201 | 9,235 | 0.924157 |
| `H:H/F/F` | 55 | 3.668862e-10 | 27,900 | 5,201 | 9,235 | 0.942526 |
| `F:F/H/F` | 55 | 1.142550e-09 | 25,738 | 3,582 | 8,997 | 0.576679 |
| `H:F/H/F` | 55 | 1.142550e-09 | 25,738 | 3,582 | 8,997 | 0.577013 |
| `F:H/H/F` | 55 | 1.142550e-09 | 27,220 | 3,582 | 8,997 | 0.372759 |
| `H:H/H/F` | 55 | 1.142550e-09 | 27,220 | 3,582 | 8,997 | 0.381826 |
| `F:F/F/H` | 55 | 4.147372e-10 | 25,260 | 5,201 | 9,031 | 0.856723 |
| `H:F/F/H` | 55 | 4.147372e-10 | 25,260 | 5,201 | 9,031 | 0.872116 |
| `F:H/F/H` | 55 | 4.147370e-10 | 26,742 | 5,201 | 9,031 | 0.648375 |
| `H:H/F/H` | 55 | 4.147370e-10 | 26,742 | 5,201 | 9,031 | 0.660428 |
| `F:F/H/H` | 55 | 1.181538e-09 | 24,580 | 3,270 | 8,793 | 0.304428 |
| `H:F/H/H` | 55 | 1.181538e-09 | 24,580 | 3,270 | 8,793 | 0.308855 |
| `F:H/H/H` | 55 | 1.181538e-09 | 26,062 | 3,270 | 8,793 | 0.103694 |
| `H:H/H/H` | 55 | 1.181538e-09 | 26,062 | 3,270 | 8,793 | 0.106965 |

## Sources, validation, and reproduction

- Default exports: [summary](results/summary.csv), [history](results/history.csv), and [solutions](results/solutions.csv).
- All-case exports: [summary](results/all/summary.csv), [history](results/all/history.csv), and [solutions](results/all/solutions.csv).
- [Implementation and metric definitions](README.md).
- The updated suite passed **378 assertions**, including all 16 hybrid assignments, both standalone baselines, nested-hybrid convergence and exact payload accounting, physical-hop communication counts, inner-oracle feasible objective agreement, and solution/feasibility checks at the new tolerance.

From the repository root:

```sh
# Four baselines plus the default mixed hybrid
julia --project=. examples/DisjointProblem/hybrid_ADMM/main.jl

# Four baselines plus all 16 hybrid assignments
julia --project=. examples/DisjointProblem/hybrid_ADMM/main.jl --all

# Validation suite
julia --project=. examples/DisjointProblem/hybrid_ADMM/runtests.jl
```

This Markdown report is a snapshot. The commands refresh CSV/PDF exports without automatically refreshing this report.
