# AGENTS.md

## Purpose

This file defines how Codex should work on the EPCN Hybrid ADMM example with the fixed physical tree and the three algorithmic modules \(M_1,M_2,M_3\).

- `formulation.md` is the mathematical source of truth.
- `plan.md` is the executable implementation plan.
- Do not invent missing mathematics.
- Keep the physical communication tree separate from the algorithmic module graph.

The current example intentionally uses:

\[
\text{global algorithm}=\text{fADMM},
\]

with three peer algorithmic modules

\[
R\rightarrow\{M_1,M_2,M_3\}.
\]

The physical tree is different from this logical star.

---

## 1. Fixed example topology

Use exactly this physical tree:

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

The final edge placement is:

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

Check the existing repository representation against this tree before changing code.

---

## 2. Fixed module decomposition

The global root \(R=(1,1)\) belongs to no local module.

Use:

```text
M1 = {(2,1),(3,1),(3,2),(3,3)}
M2 = {(3,4),(4,1),(4,2),(4,3),(5,1),(5,2),(5,3)}
M3 = {(2,2),(3,5),(4,4),(4,5)}
```

Module roots:

```text
root(M1) = (2,1)
root(M2) = (3,4)
root(M3) = (2,2)
```

For the worked example use:

```text
M1.algorithm = "fadmm"
M2.algorithm = "hadmm"
M3.algorithm = "fadmm"
```

This assignment is an example configuration, not a universal selection rule.

---

## 3. Two different graphs

### Algorithmic graph

The global root treats all three modules as peers:

```text
R -> M1
R -> M2
R -> M3
```

No module is an inner ADMM of another module.

### Physical communication tree

Messages must follow the fixed physical tree in Section 1.

In particular, every message between \(R\) and \(M_2\) passes through node `(2,2)`, which belongs to `M3`.

This physical forwarding does **not** make `M3` the algorithmic parent of `M2`.

---

## 4. Important interface interpretation

Use the design choice agreed for this example:

> `(2,2)` is the physical parent of `(3,4)` only for routing/communication. It is not the algorithmic hADMM parent of `M2`.

Therefore:

- do **not** automatically apply the ordinary hADMM parent query
  \[
  q_{2,2,34}^{k+1}
  =
  A_{2,2,34}x_{2,2}^{k+1}
  +y_{2,2,34}^k/\rho
  -c_{2,2,34}
  \]
  merely because `(2,2)->(3,4)` is a physical edge;
- do **not** make the update of `M2` depend on \(x_{2,2}^{k+1}\) unless `formulation.md` explicitly defines such a mathematical dependency;
- node `(2,2)` may relay a global fADMM query to `M2` without creating a data dependency on `M3`.

This separation is central to the new hybrid concept.

---

## 5. Before coding

Inspect and report:

```text
hADMM implementation:
fADMM implementation:
old hybrid implementation:
physical-tree representation:
module representation:
routing/path utility:
communication tracker:
test files:
```

For the fADMM implementation also report:

```text
coordinator state:
query:
target update:
response:
dual/coordinator update:
payload dimension:
update order:
```

Then compare the code with `formulation.md`.

Stop before changing the mixed mathematics if the repository does not provide a compatible module-level fADMM interface.

---

## 6. Mathematical discipline

Canonical notation:

```text
node                    (d,i)
child set               C_{d,i}
node variable           x_{d,i}
physical edge           (d,i)->(d+1,j)
edge constraint         A_{d,ij}x_{d,i}+B_{d,ij}x_{d+1,j}=c_{d,ij}
```

Do not use the physical parent-child relation to infer an algorithmic parent-child relation across module boundaries.

Do not claim convergence, fixed-point equivalence, or KKT correctness for the mixed example unless those properties are actually derived or numerically checked.

---

## 7. No hidden nested ADMM

The implementation must not behave as

```text
global fADMM iteration
    -> solve M1 ADMM to convergence
    -> solve M2 ADMM to convergence
    -> solve M3 ADMM to convergence
```

One global `step()` is one global iteration.

Each scheduled module may execute only the single local update/sweep prescribed for that global iteration.

A local numerical solver may iterate internally to solve a single proximal/argmin subproblem. That is not an additional ADMM iteration.

---

## 8. Physical routing and relay semantics

For every logical message, compute the physical path.

For each physical send record:

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

`is_relay=True` iff that physical sender is neither the logical origin nor the final logical destination.

Example global query to `M2`:

```text
logical: R -> M2

physical:
(1,1) -> (2,2)
(2,2) -> (3,4)
```

The second send is a relay transmission by `(2,2)`.

The reverse response is:

```text
(3,4) -> (2,2)
(2,2) -> (1,1)
```

and the second send is a relay transmission by `(2,2)`.

Pure forwarding does not trigger an optimization update at `(2,2)`.

---

## 9. Communication accounting

Track:

```text
TotalMessages
TotalScalarValues
NodeSentMessages(v)
MaxNodeMessages
NodeSentScalars(v)
MaxNodeScalarValues
RelayMessages
```

A logical message of payload dimension `p` traversing `h` physical hops contributes:

```text
h physical messages
h*p scalar values
```

to network totals.

---

## 10. Required tests

Tests must separately verify:

1. exact physical tree;
2. exact module membership;
3. algorithmic graph `R->{M1,M2,M3}`;
4. route from `R` to every module root;
5. relay semantics for `M2`;
6. pure hADMM regression;
7. pure fADMM regression;
8. local `M1` fADMM behavior;
9. local `M2` hADMM internal behavior;
10. local `M3` fADMM behavior;
11. one-global-step event trace;
12. no optimization event is triggered by pure relay forwarding.

The mixed global numerical update may be tested only after the module-level fADMM query/response interface in `formulation.md` is fully defined.

---

## 11. Reporting format

Before coding:

```text
1. Files to modify
2. Existing components to reuse
3. Tree/module representation
4. Mathematical rules implemented
5. Remaining unresolved rule
6. Verification plan
```

After coding:

```text
1. Files changed
2. Main logic added
3. Tests added
4. Tests actually run
5. Event/message trace
6. Communication counts
7. Remaining mathematical ambiguity, if any
```

Do not report mixed-algorithm correctness while a required interface equation remains unresolved.

---

## 12. Non-goals for this example

Do not add unless explicitly requested:

- automatic partition discovery;
- automatic hADMM/fADMM selection;
- STEP-GP;
- adaptive communication;
- wireless loss/collision models;
- asynchronous ADMM;
- convergence-rate estimation.
