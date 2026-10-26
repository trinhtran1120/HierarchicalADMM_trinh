# Candidate B hybrid ADMM

Run from the repository root with the existing Julia environment:

```sh
julia --project=. examples/DisjointProblem/hybrid_ADMM/main.jl
julia --project=. examples/DisjointProblem/hybrid_ADMM/main.jl --all
julia --project=. examples/DisjointProblem/hybrid_ADMM/runtests.jl
```

The first command compares the centralized solver, standalone physical-tree
hADMM, and one mixed hybrid configuration on the same fixed problem. `--all`
compares the central and hADMM baselines with all 16 global/module assignments.
A budget-limited run reports `converged=false`.
The tests check the formulation, centralized models, exact-block reference,
individual module solvers, all assignments, and physical relay accounting.

```julia
include("examples/DisjointProblem/hybrid_ADMM/HybridADMM.jl")
using .HybridADMM

problem = build_problem()
reference = centralized_original(problem)
config = SolverConfig(
    global_algorithm = :fADMM,
    module_algorithms = (:hADMM, :fADMM, :hADMM), # M1, M2, M3
    rho = 1000.0,
    inner_rho = 1000.0,
    tol = 1e-4,
    inner_tol = 1e-6,
    max_iter = 1000,
    inner_max_iter = 1000,
    record_messages = false,
)
result = solve_hybrid(problem; config, reference=reference.objective)
last(result.history)

report = compare_solvers(problem; configs=[config])
print_comparison(report)
save_comparison(report, "comparison_results")
```

## Comparison metrics and saved results

The default script writes `results/summary.csv`, `results/history.csv`,
`results/solutions.csv`, `results/convergence.pdf`, and `results/resources.pdf`
beside `main.jl`. Generated results are ignored by Git. `solutions.csv` contains
every coordinate at every physical node for every method, its central reference,
and its absolute error. Histories include the distributed methods; the central
reference is a single solution without an ADMM trajectory.

| Metric | Definition |
| --- | --- |
| `objective_gap` | `abs(J-Jcentral)/max(1,abs(Jcentral))`, using actual node costs |
| `feasible_gap` | Same gap after reconstructing all decisions from the current root vector |
| `solution_error` | Maximum coordinate error over all physical node decisions |
| `relative_solution_error` | Stacked solution error in 2-norm divided by `max(1,norm(reference))` |
| `physical` | Maximum infinity norm of all original physical-edge constraint residuals |
| `primal`, `dual` | Each method's ADMM stopping residuals; central N/A |
| `iterations` | Ipopt barrier iterations, hADMM full-tree sweeps, or hybrid outer iterations, as indicated by `iteration_kind` |
| `inner_sweeps` | Sum of hybrid inner sweeps over all modules and outer iterations |
| `TotalNodeUpdates`, `MaxNodeUpdates` | Total/maximum per-node proximal updates, including hybrid inner iterations; central N/A |
| `TotalVals`, `MaxNodeVals` | Total/maximum per-node scalar transmissions, including physical relays and stopping reports |
| `seconds` | Sequential wall time after a separate full warm-up of each solver |

The central objective is evaluated from its returned physical solution, avoiding
small differences from Ipopt's L1 epigraph slack tolerances. The epigraph model
objective is also available as `reference.solver_objective`. Feasible objective
gaps are used in the convergence plots so transient infeasible node copies do
not suggest a better feasible objective. Curves are plotted without monotone
smoothing; values below `1e-16` are clipped only for logarithmic display.

All distributed runs start with primal variables at one and scaled duals at
zero, use the same outer `rho`, `tol`, and iteration budget, and solve the same
costs. Central decision starts are also one. Ipopt uses `1e-10` solver tolerance
for the reference. Standalone hADMM calls the existing `hierarchicalADMM!` with
`lambda=1/rho`; that routine now propagates the requested penalty to every node.
It counts one subtree stopping scalar on each child-to-parent response, matching
the local hybrid hADMM accounting.

Ipopt iterations, whole-tree sweeps, and hybrid outer iterations represent
different work; use node updates and communications alongside them. Node updates
count proximal solves rather than PANOC's internal numerical steps. Central
communication is zero during optimization with problem data already resident;
model upload and solution distribution are excluded for all methods.

Runtime includes model/solver setup and diagnostics, excludes plotting, output,
and warm-up, and measures this sequential simulation rather than parallel
network execution. Call `compare_solvers(...; warmup=false)` to skip warm-up;
those times include compilation. The centralized reformulation check is outside
timed runs.

Each solve starts from ones with zero scaled duals. Module decisions and internal
duals are warm-started between outer iterations. History entries include each
module's residual, iterations, and convergence flag. `reference` is optional and
is used only to calculate the relative objective gap. Stopping requires outer
primal/dual and internal feasibility residuals within `tol`, and all inner solves
within `inner_tol`. Inner budget exhaustion is reported; it never counts as an
accurate block minimization. No convergence guarantee for arbitrary inexact
iterations is assumed.

## Problem and existing code reuse

`HybridADMM.jl` isolates the existing example files in a Julia module, leaving
their standalone behavior unchanged. It uses:

- `TreeGeneration.jl` / `src/linknodes.jl`: node construction, ordered child
  vectors, `assign!`, quartic costs, PANOC solvers, resets, scalar transmission
  counters, and legacy result helpers.
- `CentralSolution.jl`: an independent leaf-only centralized reference.
- `hADMM.jl`: the analytic `CostFunc` proximal interface and `Proximal_Iteration`.
  The local hierarchical driver follows its preorder primal and return-time
  edge-dual update schedule.
- `fADMM.jl`: physical path/block-position helpers and `update_leaf_flat` for
  flattened module non-root updates.

The random tree generator is replaced by exactly the fixed tree in `guide.md`.
Physical nodes use `(depth,index)` keys; legacy numeric IDs are assigned in that
lexicographic order. Leaves have dimension one, as in the existing examples.
Internal dimensions are the sum of physical child dimensions. Every node uses
the existing cost

\[
f_v(x)=\tfrac12\left(\sum_i(x_i-\mathrm{ID}_v)^2\right)^2+2\|x\|_1.
\]

The partition, physical edge operators `A/B/c`, module `D/d/S`, physical routes,
and solver choices are separate. Module stacks use lexicographic node order;
child blocks use increasing physical child index. M3 retains its full `(2,2)`
root vector, including the unconstrained internal component corresponding to
M2. Its internal `(2,2)--(3,5)` equation selects only the `(3,5)` block.

`centralized_reformulation` builds (9) explicitly, with all physical node
variables, internal constraints, auxiliary copies, and all three boundary
equations. `feasible_decisions` and `formulation_residuals` expose the projection
equivalence checks independently of any iterative solver.

## Global fADMM and hADMM derivation

For this disjoint problem, the global feasible set can be eliminated exactly.
Write `y=x_(1,1)` and let `T_m` select the physical module-root block from `y`.
Then every `Z` in `G` has `z_m=T_m*y`. In particular, `T_2` selects the M2 block
inside `T_3`, so `A_(2,24)*z_3=z_2` holds identically. The resulting computational
star has root block `y`, module blocks `X_m`, and edges

\[
S_mX_m-T_my=0.
\]

The root owns its cost; module `m` owns `F_m + indicator(C_m)`. Group all module
blocks into the second ADMM block. The root-first updates are exactly

\[
\begin{aligned}
y^+ &= \arg\min_y f_{1,1}(y)
       +\tfrac\rho2\sum_m\|S_mX_m-T_my+u_m\|^2,\\
X_m^+ &= \arg\min_{X_m\in C_m}F_m(X_m)
       +\tfrac\rho2\|S_mX_m-T_my^++u_m\|^2,\\
u_m^+ &= u_m+S_mX_m^+-T_my^+.
\end{aligned}
\]

There is only one module depth and no edges between computational module blocks.
Consequently the global flattened and hierarchical root-first schedules have
the same updates and communication routes. Both global config labels explicitly
share `global_update`, and results report `global_updates_equivalent=true`.
They do not claim distinct global algorithms or different convergence behavior.
All three original boundary equations remain enforced through `T_m`.

`solve_oracle` executes these updates with high-accuracy centralized Ipopt solves
of each module subproblem. It is a diagnostic reference and reports no
distributed communication costs. The normal solver always uses the selected
distributed module solvers; unsupported labels throw an error.

The reported outer primal residual is `max_m norm(S_m*X_m-z_m, Inf)`. For this
Z-first splitting, the dual residual is
`rho*max_m norm(S_m*(X_m-X_m_old), Inf)`. Internal and original physical-edge
residuals are also recorded. Objective values use the current physical node
copies, which can be slightly infeasible before convergence.

## Inner updates

For local fADMM, flatten only the module's internal nodes. With module root `r`
and root-relative block selectors `R_v`, use constraints `R_v*x_r=x_v` for all
other module nodes. Transitive physical disjoint constraints make this equivalent
to `D_m*X_m=0`. Update the root with its own cost, the outer penalty, and all
flattened internal penalties. Then update every other node independently and
update the flattened duals at the coordinator. Requests/replies follow physical
module paths, counting relays.

For local hADMM, store one scaled dual for each internal physical edge. At each
node, minimize its cost plus its incident internal quadratic penalties; the
module root additionally receives `rho/2*norm(x_r-target)^2`. Traverse from root
to children, then update each edge dual on return. Edge selectors are retained,
so M3's external M2 block is included in the outer penalty but excluded from
internal penalties. These are the same signed queries as the existing hADMM:
parent-to-child `A*x_parent+u`, child-to-parent `x_child`, with parent update
target `x_child-u`.

Both solvers use `target=z_m-u_m`. Their stopping residual is the maximum
internal consensus residual and `inner_rho` times node changes during a sweep.
This includes root changes. It is an iterative stopping measure, not a claim of
exact solution or a theorem for mixed/inexact updates. Tests compare both local
methods against the independent penalized module oracle.

## Communication accounting

`TotalVals` sums transmitted scalars over senders and physical hops.
`MaxNodeVals` is the maximum sender total. Root-to-M2 requests and replies use
`(1,1) <-> (2,2) <-> (3,4)`; `(2,2)` forwarding increments communication counters
without incrementing its optimization counter. `record_messages=true` retains
one `Transmission` per physical hop for inspection; it is off by default to
avoid storing a large trace.

Payloads counted on every hop are:

| Exchange | Payload |
| --- | --- |
| Global query | Module-root target vector |
| Global response | Module-root decision, inner residual, iteration count, convergence flag |
| Local fADMM query / response | Target vector / node decision vector |
| Local hADMM query / response | Target vector / node decision vector and subtree stopping residual |

The coordinator can compute stopping quantities from these payloads and cached
previous decisions. Static topology, solver configuration, initial values,
message headers, and synchronization signals are excluded. No dynamic solver
residual payload is silently omitted. Centralized validation, objective/error
diagnostics, and history collection run in the simulation harness and are
excluded from distributed optimization traffic, matching the existing examples.
