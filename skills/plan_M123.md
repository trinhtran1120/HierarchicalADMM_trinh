# plan.md

# Current Plan: Run the \(M_1,M_2,M_3\) Hybrid ADMM Example

## 1. Immediate goal

Replace the old six-node toy fixture with the exact physical tree and module layout used in the current EPCN hybrid-ADMM discussion.

The implementation target is:

\[
\boxed{
\text{global fADMM over peer modules }
M_1,M_2,M_3
}
\]

with

```text
M1 = local fADMM
M2 = local hADMM
M3 = local fADMM
```

and with communication routed over the physical tree.

The first objective is to make the repository reproduce the correct:

- physical tree;
- module partition;
- logical module graph;
- physical routes;
- per-module local update events;
- communication accounting.

A full mixed numerical correctness claim remains blocked only by the module-level fADMM-to-\(M_2\) interface identified in `formulation.md` Section 10.

---

# 2. Exact physical tree

Use:

```text
                               R=(1,1)
                              /       \
                         (2,1)         (2,2)
                        /  |  \        /    \
                    (3,1)(3,2)(3,3) (3,4)  (3,5)
                                      / | \    / \
                                  (4,1)(4,2)(4,3)(4,4)(4,5)
                                   / \          \
                               (5,1)(5,2)       (5,3)
```

Physical edges:

```text
(1,1) -> (2,1)
(1,1) -> (2,2)

(2,1) -> (3,1)
(2,1) -> (3,2)
(2,1) -> (3,3)

(2,2) -> (3,4)
(2,2) -> (3,5)

(3,4) -> (4,1)
(3,4) -> (4,2)
(3,4) -> (4,3)

(4,1) -> (5,1)
(4,1) -> (5,2)
(4,3) -> (5,3)

(3,5) -> (4,4)
(3,5) -> (4,5)
```

Add an exact topology test.

---

# 3. Exact module partition

Use:

| Module | Physical nodes | Module root | Algorithm |
|---|---|---|---|
| `M1` | `(2,1),(3,1),(3,2),(3,3)` | `(2,1)` | `fadmm` |
| `M2` | `(3,4),(4,1),(4,2),(4,3),(5,1),(5,2),(5,3)` | `(3,4)` | `hadmm` |
| `M3` | `(2,2),(3,5),(4,4),(4,5)` | `(2,2)` | `fadmm` |

The global root `(1,1)` belongs to no local module.

Add a module-membership test that fails if any physical node belongs to the wrong module.

---

# 4. Algorithmic graph

Represent separately from the physical tree:

```text
R -> M1
R -> M2
R -> M3
```

This is the global fADMM star.

Do **not** represent:

```text
R -> M3 -> M2
```

as an algorithmic graph.

That chain exists only as a communication consequence for routes to `M2`.

---

# 5. Physical routes to module roots

Expected shortest physical paths:

```text
R -> M1:
(1,1) -> (2,1)

R -> M3:
(1,1) -> (2,2)

R -> M2:
(1,1) -> (2,2) -> (3,4)
```

Reverse response paths are the reverse of these routes.

Add exact route tests.

For `R -> M2`, `(2,2)` is a relay on the second physical send.

It must not receive an optimization event solely because it forwards the global query.

---

# 6. Repository inspection

Before modifying implementation, report:

```text
hADMM implementation:
fADMM implementation:
old hybrid implementation:
physical-tree class:
module class:
routing/path function:
communication tracker:
existing tests:
```

For fADMM also report:

```text
root/coordinator state:
query:
local prox/update:
response:
dual update:
payload dimensions:
one-step ordering:
```

Reuse existing solvers where possible.

---

# 7. Data structures

Use or minimally adapt:

```python
HybridModule(
    id,
    physical_nodes,
    module_root,
    algorithm,     # "hadmm" | "fadmm"
)
```

The global root is represented separately.

The solver must store both:

```text
physical_tree
algorithmic_module_graph
```

Do not derive one from the other.

---

# 8. Communication events

For every physical transmission record:

```text
iteration
logical_origin
logical_destination
physical_sender
physical_receiver
message_type
payload_dimension
is_relay
```

Required communication counters:

```text
TotalMessages
TotalScalarValues
NodeSentMessages(v)
MaxNodeMessages
NodeSentScalars(v)
MaxNodeScalarValues
RelayMessages
```

A logical message of dimension `p` on the path

```text
(1,1) -> (2,2) -> (3,4)
```

contributes:

```text
2 physical messages
2*p scalar values
1 relay transmission
```

for the query.

The response contributes the same again.

Therefore an `R <-> M2` scalar round trip contributes:

```text
4 physical messages
4 scalar values
2 relay transmissions
```

assuming payload dimension 1 in both directions.

---

# 9. Local M1 test

`M1` uses fADMM.

Topology:

```text
        (2,1)
       /  |  \
   (3,1)(3,2)(3,3)
```

This is a local star physically and algorithmically.

Use the existing fADMM solver/update logic.

Verify one local M1 step produces the expected sequence:

```text
LOCAL_ROOT_UPDATE(M1)
QUERY(3,1)
QUERY(3,2)
QUERY(3,3)
PROX/UPDATE targets
RESPONSES
LOCAL_DUAL_UPDATE(M1)
```

All local M1 query/response paths are one physical hop.

---

# 10. Local M2 test

`M2` uses hADMM.

Internal topology:

```text
             (3,4)
           /   |   \
       (4,1) (4,2) (4,3)
        / \           \
    (5,1)(5,2)       (5,3)
```

Only these internal edges belong to the local hADMM:

```text
(3,4)->(4,1)
(3,4)->(4,2)
(3,4)->(4,3)
(4,1)->(5,1)
(4,1)->(5,2)
(4,3)->(5,3)
```

The cross-module physical edge

```text
(2,2)->(3,4)
```

must not be treated as an ordinary local hADMM parent edge.

Verify one local hADMM sweep over M2 uses:

```text
PRIMAL(node) once per local node
QUERY(edge) once per internal M2 edge
RESPONSE(edge) once per internal M2 edge
DUAL(edge) once per internal M2 edge
```

subject to the module-root incoming query interface in `formulation.md`.

---

# 11. Local M3 test

`M3` uses fADMM.

Physical topology:

```text
       (2,2)
          \
          (3,5)
          /   \
      (4,4)  (4,5)
```

Algorithmically local fADMM flattens its internal targets under module root `(2,2)`.

Therefore logical queries to `(4,4)` and `(4,5)` are routed through `(3,5)`.

Required path tests:

```text
(2,2) -> (4,4):
(2,2)->(3,5)->(4,4)

(2,2) -> (4,5):
(2,2)->(3,5)->(4,5)
```

Verify relay events at `(3,5)`.

---

# 12. Global scheduler

One global `step()` must expose the following conceptual events:

```text
GLOBAL_ROOT_UPDATE

GLOBAL_QUERY(M1)
LOCAL_STEP(M1)
GLOBAL_RESPONSE(M1)

GLOBAL_QUERY(M3)
LOCAL_STEP(M3)
GLOBAL_RESPONSE(M3)

GLOBAL_QUERY(M2)
ROUTE via (2,2)
LOCAL_STEP(M2)
GLOBAL_RESPONSE(M2)
ROUTE via (2,2)

GLOBAL_DUAL_UPDATE
```

The exact order of independent peer modules may follow the existing fADMM implementation.

The scheduler must never run an entire local ADMM solve to convergence before returning.

---

# 13. Deterministic scalar fixture

Use one scalar variable per physical node.

For pure/local hADMM validation use

\[
x_{\text{parent}}-x_{\text{child}}=0
\]

on applicable internal physical edges.

Use local objectives

\[
f_v(x)=\frac12 a_v(x-b_v)^2.
\]

Coefficients:

| Node | \(a\) | \(b\) |
|---|---:|---:|
| (1,1) | 1.00 | 1.00 |
| (2,1) | 1.10 | -0.50 |
| (2,2) | 0.90 | 2.00 |
| (3,1) | 1.20 | 0.00 |
| (3,2) | 1.30 | 1.25 |
| (3,3) | 0.80 | -1.00 |
| (3,4) | 1.40 | 0.75 |
| (3,5) | 1.00 | -0.25 |
| (4,1) | 1.15 | 1.50 |
| (4,2) | 0.95 | -1.25 |
| (4,3) | 1.25 | 0.50 |
| (4,4) | 1.05 | 2.25 |
| (4,5) | 0.85 | -0.75 |
| (5,1) | 1.35 | 1.00 |
| (5,2) | 0.75 | -1.50 |
| (5,3) | 1.45 | 0.25 |

Initialization:

```text
rho = 1.0
x0 = 0
dual0 = 0
seed = 0
```

Do not force the general-affine scalar fixture into an fADMM implementation that uses a different consensus-copy formulation.

---

# 14. Pure-mode regression on the full 16-node physical tree

Replace the old six-node pure-mode fixture with this 16-node tree.

## Pure hADMM

Use one module containing all non-root and root nodes and run the existing hADMM on the full physical tree.

Compare the hybrid-framework all-hADMM configuration with standalone hADMM after fixed `K`.

Default target:

```text
K = 10
max iterate error <= 1e-10
```

if closed-form/deterministic local solves permit it.

## Pure fADMM

Use the existing compatible fADMM formulation on the same physical topology when possible.

Verify that logical messages to deep nodes are routed over all physical hops.

Do not force numerical equivalence between hADMM and fADMM problems if their variable-copy formulations differ.

---

# 15. Exact routing checks

At minimum test these routes:

```text
(1,1) -> (3,4):
(1,1),(2,2),(3,4)

(1,1) -> (5,1):
(1,1),(2,2),(3,4),(4,1),(5,1)

(1,1) -> (5,3):
(1,1),(2,2),(3,4),(4,3),(5,3)

(2,2) -> (4,4):
(2,2),(3,5),(4,4)

(2,2) -> (4,5):
(2,2),(3,5),(4,5)
```

Also test reverse paths.

---

# 16. Mixed-mode blocker that remains

Do not invent this step:

```text
global fADMM query q2
    -> exact incoming mathematical quantity for M2 root (3,4)
    -> one hADMM local sweep
    -> exact response xi2 returned to global fADMM
```

The routing and scheduling can be implemented now.

The exact mathematical transformation above must be defined before mixed numerical results are labeled correct.

This is narrower than the old blocker: physical/algorithmic parenthood and coordinator structure are now fixed.

---

# 17. Required deliverables

After implementation report:

```text
1. Files changed
2. Exact 16-node tree representation
3. Exact M1/M2/M3 membership
4. Algorithmic graph representation
5. R-to-M2 relay trace through (2,2)
6. Local M1 fADMM event trace
7. Local M2 hADMM event trace
8. Local M3 fADMM event trace
9. Pure hADMM regression result
10. Pure fADMM regression result
11. Communication metrics for one example iteration
12. Explicit status of the remaining global-fADMM/M2 interface
```

---

# 18. Do not implement yet

Do not add:

- automatic partition discovery;
- automatic algorithm selection;
- STEP-GP;
- wireless traffic/collision model;
- asynchronous updates;
- convergence-rate estimation;
- learned proximal operators.
