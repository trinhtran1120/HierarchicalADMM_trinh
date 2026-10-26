# Codex implementation guide — Candidate B hybrid ADMM (disjoint problem)

## Objective and scope

Implement the **same optimization problem** on the fixed 16-node tree below, but permit **independent algorithm selection** (`fADMM` or `hADMM`) at the global coordinator and within each module. The partition and auxiliary-variable formulation must **not** depend on the selected algorithm. Start with the disjoint problem from `hADMM.md`; do not generalize to arbitrary affine coupling unless separately requested.

Success means: (1) the reformulated feasible set projects exactly onto the original disjoint feasible set, (2) objective values match an independently solved centralized reference, (3) every physical communication event uses a real tree edge, including relays, and (4) all supported global/module algorithm combinations can be selected by configuration rather than edits to mathematical or partitioning code.

**Reference:** `hADMM.md`, especially “Hierarchical optimization over tree networks,” “Disjoint problems over tree networks,” “Flattened ADMM,” and “Hierarchical ADMM algorithm.” The *hybrid splitting and auxiliary copies below are proposed design choices*, not an existing theorem from that document.

## 1. Immutable physical tree and Candidate B partition

Each node is indexed `(d,i)`, with `d=1,...,5` and root `(1,1)`.

```text
(1,1)                       [global root]
├── (2,1)                   [M1]
│   ├── (3,1)               [M1]
│   ├── (3,2)               [M1]
│   └── (3,3)               [M1]
└── (2,2)                   [M3]
    ├── (3,4)               [M2]
    │   ├── (4,1)           [M2]
    │   │   ├── (5,1)       [M2]
    │   │   └── (5,2)       [M2]
    │   ├── (4,2)           [M2]
    │   └── (4,3)           [M2]
    │       └── (5,3)       [M2]
    └── (3,5)               [M3]
        ├── (4,4)           [M3]
        └── (4,5)           [M3]
```

Exactly 16 nodes, 15 original edges. The global root is **not a member** of any module.

- `V1 = {(2,1),(3,1),(3,2),(3,3)}` (4 nodes), coordinator `(2,1)`.
- `V2 = {(3,4),(4,1),(4,2),(4,3),(5,1),(5,2),(5,3)}` (7 nodes), coordinator `(3,4)`.
- `V3 = {(2,2),(3,5),(4,4),(4,5)}` (4 nodes), coordinator `(2,2)`.

There are exactly **12 internal edges** (3, 6, 3 respectively) and **3 boundary edges**:

1. `(1,1)--(2,1)` (global--M1)
2. `(1,1)--(2,2)` (global--M3)
3. `(2,2)--(3,4)` (M3--M2)

A module need not contain every descendant of its physical root: M3 omits the branch rooted at `(3,4)`.

## 2. Original disjoint optimization problem

For each node `v=(d,i)`, define local decision vector `x_{d,i}` of dimension `n_{d,i}` and proper convex local cost `f_{d,i}(x_{d,i})`.

The **disjoint** condition is that every parent variable is the ordered concatenation of its physical children:

\[
 x_{d,i}=[x_{d+1,j}]_{j\in C_{d,i}},\qquad d=1,\ldots,4.
\]

Consequently, `n_{d,i} = sum_{j in C_{d,i}} n_{d+1,j}`. Set the child order by increasing `j`; never infer order from Python sets. For leaves, dimensions and convex costs are problem data. The original problem is

\[
\begin{aligned}
\min_{\{x_{d,i}\}}\quad&J=\sum_{d=1}^5\sum_{i\in P_d\cup L_d} f_{d,i}(x_{d,i})\\
\text{s.t.}\quad&x_{d,i}=[x_{d+1,j}]_{j\in C_{d,i}},\quad \forall\text{ physical parent nodes}.
\end{aligned}
\tag{1}
\]

In the general notation of `hADMM.md`, each physical edge `(d,i)--(d+1,j)` carries

\[
 A_{d,ij}x_{d,i}+B_{d,ij}x_{d+1,j}=c_{d,ij},
\tag{2}
\]

where, **for disjoint problems only**, `A_{d,ij}` selects the child `j` block from the parent, `B_{d,ij}=-I_{n_{d+1,j}}`, `c_{d,ij}=0`. Thus (2) asserts `A_{d,ij}x_{d,i}=x_{d+1,j}`. Store these matrices or equivalent linear operators for each actual physical edge, with dimensions tested.

## 3. Module-level stacking and internal constraints

Let `X_m=[x_{d,i}]_{(d,i) in V_m}` in a **deterministic documented order**. Let

\[
 F_m(X_m)=\sum_{(d,i)\in V_m}f_{d,i}(x_{d,i}),\qquad
 \mathcal C_m=\{X_m:D_mX_m=d_m\}.
\tag{3}
\]

Build `D_m X_m=d_m` by stacking (2) for **internal edges only**. Boundary edges must not be inserted in `D_m`.

Internal edges:

- M1: `(2,1)--(3,1)`, `(2,1)--(3,2)`, `(2,1)--(3,3)`.
- M2: `(3,4)--(4,1)`, `(3,4)--(4,2)`, `(3,4)--(4,3)`, `(4,1)--(5,1)`, `(4,1)--(5,2)`, `(4,3)--(5,3)`.
- M3: `(2,2)--(3,5)`, `(3,5)--(4,4)`, `(3,5)--(4,5)`.

**Critical point for the disjoint case:** M3 does *not* internally require `x_{2,2}=x_{3,5}`; it requires only the *selected* `(3,5)` block of `x_{2,2}` to equal `x_{3,5}`. Its `(3,4)` block remains a boundary component. This is why the original `A_{2,24}` must be retained.

## 4. Global auxiliary copies / virtual coordinator connections

Define module-root selection operators

\[
 S_1X_1=x_{2,1},\qquad S_2X_2=x_{3,4},\qquad S_3X_3=x_{2,2}.
\tag{4}
\]

Introduce **global auxiliary copies** `z_1`, `z_2`, `z_3` with the matching dimensions and consistency constraints

\[
 \boxed{S_mX_m=z_m,\quad m=1,2,3.}
\tag{5}
\]

In particular, `z_2=x_{3,4}` is the **virtual root-to-M2 consistency constraint**. It is a *computational link*, not a new physical edge and not a substitute for the original M3--M2 physical coupling.

The three boundary constraints now involve only `x_{1,1}` and the global copies:

\[
\begin{aligned}
 A_{1,11}x_{1,1}+B_{1,11}z_1&=c_{1,11},\\
 A_{1,12}x_{1,1}+B_{1,12}z_3&=c_{1,12},\\
 A_{2,24}z_3+B_{2,24}z_2&=c_{2,24}.
\end{aligned}
\tag{6}
\]

For disjoint constraints these reduce to the selected-block equalities:

\[
\begin{aligned}
 A_{1,11}x_{1,1}&=z_1,\\
 A_{1,12}x_{1,1}&=z_3,\\
 A_{2,24}z_3&=z_2.
\end{aligned}
\tag{7}
\]

Set `Z=(x_{1,1},z_1,z_2,z_3)` and define the global feasible set

\[
\mathcal G=\{Z:\text{all three equations in (6) hold}\}.
\tag{8}
\]

The root does not directly physically connect to `(3,4)`. To send a query to M2 or receive its reply, route it along `(1,1) <-> (2,2) <-> (3,4)`. During that communication, `(2,2)` is a **relay**, not the solver of M2's subproblem.

## 5. Equivalent hybrid optimization formulation (algorithm-independent)

\[
\boxed{
\begin{aligned}
\min_{Z,X_1,X_2,X_3}\quad&f_{1,1}(x_{1,1})+\sum_{m=1}^3 F_m(X_m)\\
\text{s.t.}\quad& Z\in\mathcal G,\\
& X_m\in\mathcal C_m,\quad m=1,2,3,\\
& S_mX_m=z_m,\quad m=1,2,3.
\end{aligned}}
\tag{9}
\]

Equivalently, with indicator functions,

\[
 H(Z)=f_{1,1}(x_{1,1})+\delta_{\mathcal G}(Z),\qquad
 \widetilde F_m(X_m)=F_m(X_m)+\delta_{\mathcal C_m}(X_m),
\]

\[
 \min_{Z,X_1,X_2,X_3} H(Z)+\sum_{m=1}^3\widetilde F_m(X_m)
 \quad\text{s.t. } S_mX_m-z_m=0\ (m=1,2,3).
\tag{10}
\]

**Equivalence proof (implementation invariant):** Any original feasible collection `x` gives a feasible point in (9) by stacking each module and setting `z_m=S_mX_m`. Conversely, (5) substitutes the original boundary variables into (6), restoring precisely the three omitted original constraints; internal constraints are unchanged. Thus feasible physical node decisions and objective values coincide with the original problem.

## 6. Algorithm selection — independent of the formulation

Expose a configuration conceptually equivalent to:

```yaml
global_algorithm: fADMM  # or hADMM
module_algorithms:
  M1: hADMM              # or fADMM
  M2: hADMM              # or fADMM
  M3: hADMM              # or fADMM
```

The values above are **an example, not defaults that may be hardcoded into partition/model code**. Support both algorithm choices at all four positions. Keep `physical_tree`, `partition`, `edge_constraints`, `auxiliary_copies`, and `solver_assignment` separate.

Each solver needs a well-defined problem/communication interface. The outer solver coordinates module blocks; each inner solver receives its current boundary query (which includes the auxiliary consistency penalty or equivalent ADMM information), returns a module-boundary response, and reports convergence/accuracy. A local solver must include its assigned external boundary term in the module-root update rather than solving `F_m` alone.

### Concrete global fADMM-compatible splitting

One supported outer-solver interpretation of (10) is a two-block, scaled-form ADMM with group blocks `Z` and `(X1,X2,X3)`:

\[
 Z^{k+1}=\arg\min_{Z\in\mathcal G}
 f_{1,1}(x_{1,1})+\frac\rho2\sum_{m=1}^3
 \|S_mX_m^k-z_m+u_m^k\|_2^2,
\tag{11}
\]

\[
 X_m^{k+1}=\arg\min_{X_m\in\mathcal C_m}
 F_m(X_m)+\frac\rho2\|S_mX_m-z_m^{k+1}+u_m^k\|_2^2,
\tag{12}
\]

\[
 u_m^{k+1}=u_m^k+S_mX_m^{k+1}-z_m^{k+1}.
\tag{13}
\]

All three `X_m` minimizations separate. For a local `hADMM`, solve/approximate (12) via tree-aware internal node updates on the module's actual internal edges; for a local `fADMM`, flatten *that module's internal computational communication* through relays in its own physical tree. Do **not** silently replace (12) with an unrelated proximal operator.

### Global hADMM selection requires an explicit update design

Equations (9)–(10) remain unchanged if `global_algorithm=hADMM`, but **the hADMM of `hADMM.md` cannot simply be pasted onto the original physical level indices**: the global virtual computational topology is a star of module blocks, and the global cost `H` includes *cross-module boundary equations* from (6). A global hADMM implementation must explicitly define valid global computational blocks/edge constraints and corresponding updates for (9), and verify it actually enforces the M2--M3 boundary equation. Treat this as a solver-specific derivation and test requirement; do not claim the two solver choices are automatically interchangeable mathematically or convergence-equivalent.

When the global virtual topology is only a two-level star, distinguish fADMM and hADMM *implementations* only where their actual update/communication schedules differ; don't invent a difference solely to satisfy a config label. If a particular combination is not yet implemented, fail explicitly as unsupported rather than silently falling back to another algorithm.

### Inner accuracy and convergence caveat

Equation (12) states an **exact** block minimization. Running only a fixed number of hADMM/fADMM inner iterations gives an **inexact** outer update. Expose an inner tolerance/iteration budget and report its residual. Do not assume convergence of an arbitrary mixed or inexact hybrid from the convergence theorem for standalone hADMM in `hADMM.md`.

## 7. Communication accounting

Maintain **two distinct graphs**:

- **Computational graph:** the global root coordinates abstract module blocks; `z_m=S_mX_m` creates virtual root-to-module links.
- **Physical graph:** only the 15 edges in Section 1 may transmit; all virtual communications are routed along unique physical paths.

For M2 the global-root route is `(1,1)--(2,2)--(3,4)` (two physical hops each way). Count the relay transmission at `(2,2)` in `TotalVals` and `MaxNodeVals`. For local hADMM use module-internal parent--child exchanges; for local fADMM count any module-internal relays. Count actual transmitted **scalar values** (query and response lengths can differ), not just message counts. Keep optimization updates separate from pure forwarding actions. Exclude no data silently; document any omitted metadata or stopping checks.

## 8. Minimal implementation sequence and checks

1. **Build tree and deterministic indices** → assert 16 unique nodes, 15 edges, each non-root has exactly one parent, and child-block dimensions concatenate correctly.
2. **Build modules and edge partition** → assert disjoint membership, global root excluded, module sizes `(4,7,4)`, internal edge counts `(3,6,3)`, and exactly the three boundary edges above.
3. **Build linear operators** → check every `A/B/c`, `D_m/d_m`, and `S_m` has correct dimensions and picks the expected node slice; verify `A_{2,24}z_3=z_2` in the disjoint case.
4. **Check equivalence before solving** → construct feasible leaf variables, recursively concatenate parent variables, build `X_m,z_m`, verify all (9) residuals are zero; perturb a boundary copy and verify consistency/boundary residual becomes nonzero.
5. **Solve a small convex test problem centrally** → use fixed convex node objectives and leaf dimensions. Compare centralized original (1) against centralized reformulation (9): feasible physical variables, objective, and residuals must agree within stated numerical tolerance.
6. **Implement algorithms behind interfaces** → run an exact or high-accuracy oracle for (11)–(13) first, then implement local fADMM/hADMM individually with inner residual tracking. Add global hADMM only with its explicitly justified block formulation.
7. **Exercise configuration** → test all 2^4 algorithm selections once both solver choices actually exist at every level; report objective gap, primal/dual residuals, iterations, `TotalVals`, and `MaxNodeVals`. Do not count a test as passed merely because a config is accepted.
8. **Verify routes** → no virtual one-hop `(1,1)--(3,4)` messages; every root-to-M2 request and reply is forwarded at `(2,2)` and counted there.

## Out of scope for this first implementation

Do not change the physical tree, optimize partitions, add learned response models, introduce non-disjoint affine constraints, or assert a novel convergence guarantee. Prioritize exact reformulation, configurable solver architecture, and explicit numerical checks.
