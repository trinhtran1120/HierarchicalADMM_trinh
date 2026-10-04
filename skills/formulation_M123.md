# formulation.md

# Hybrid ADMM Mathematical Specification for the \(M_1,M_2,M_3\) Example

## 1. Purpose and status

This document specializes the EPCN hybrid-ADMM discussion to the fixed tree with modules \(M_1,M_2,M_3\).

It separates:

1. source-established hADMM/fADMM mathematics;
2. design choices already agreed for this example;
3. the agreed boundary-variable interface for M2 and the limits of a one-sweep module response.

Section 10 defines the M2 interface. This mathematical definition does not by itself establish implementation correctness or mixed convergence.

The central design choice is:

\[
\boxed{
\text{algorithmic structure} \neq \text{physical communication structure}.
}
\]

The global algorithm uses fADMM and treats \(M_1,M_2,M_3\) as peer algorithmic modules.

---

# 2. Physical tree

The root is

\[
R=(1,1).
\]

The physical tree is

```text
                               (1,1)
                              /     \
                          (2,1)     (2,2)
                         / | \       /   \
                     (3,1)(3,2)(3,3)(3,4)(3,5)
                                      /|\    /\
                                  (4,1)(4,2)(4,3)(4,4)(4,5)
                                   / \          \
                               (5,1)(5,2)       (5,3)
```

The physical edge set is

\[
\begin{aligned}
E_{\mathrm{phys}}=\{&
(1,1)\!\to\!(2,1),\;
(1,1)\!\to\!(2,2),\\
&(2,1)\!\to\!(3,1),\;
(2,1)\!\to\!(3,2),\;
(2,1)\!\to\!(3,3),\\
&(2,2)\!\to\!(3,4),\;
(2,2)\!\to\!(3,5),\\
&(3,4)\!\to\!(4,1),\;
(3,4)\!\to\!(4,2),\;
(3,4)\!\to\!(4,3),\\
&(4,1)\!\to\!(5,1),\;
(4,1)\!\to\!(5,2),\;
(4,3)\!\to\!(5,3),\\
&(3,5)\!\to\!(4,4),\;
(3,5)\!\to\!(4,5)
\}.
\end{aligned}
\]

For a pure-tree reference problem, every physical edge may carry the standard affine coupling

\[
A_{d,ij}x_{d,i}
+
B_{d,ij}x_{d+1,j}
=
c_{d,ij}.
\]

In the mixed algorithm, the cross-module edge `(2,2)->(3,4)` is routing-only and has no hADMM coupling or hADMM dual; Section 10 supplies the global boundary coupling instead.

For the scalar numerical example, a consensus specialization may be used:

\[
A_{d,ij}=1,\qquad
B_{d,ij}=-1,\qquad
c_{d,ij}=0.
\]

---

# 3. Module partition

The global root does not belong to any local module.

Define

\[
M_1=
\{(2,1),(3,1),(3,2),(3,3)\},
\]

\[
M_2=
\{(3,4),(4,1),(4,2),(4,3),(5,1),(5,2),(5,3)\},
\]

and

\[
M_3=
\{(2,2),(3,5),(4,4),(4,5)\}.
\]

Module roots are

\[
r_1=(2,1),\qquad
r_2=(3,4),\qquad
r_3=(2,2).
\]

The worked assignment is

\[
a_1=\mathrm{fADMM},\qquad
a_2=\mathrm{hADMM},\qquad
a_3=\mathrm{fADMM}.
\]

This assignment is fixed only for this example.

---

# 4. Algorithmic graph

The global root uses fADMM over the three modules:

\[
\boxed{
G_{\mathrm{alg}}:
R\to M_1,\qquad
R\to M_2,\qquad
R\to M_3.
}
\]

Thus \(M_1,M_2,M_3\) are peer algorithmic modules.

In particular,

\[
M_2\not\subset_{\mathrm{alg}} M_3.
\]

No local module is an inner ADMM solved to convergence inside another module.

---

# 5. Communication graph and routing

Communication follows the physical tree.

The module-level routes induced by the physical tree are

\[
R\leftrightarrow M_1
\]

through

\[
(1,1)\leftrightarrow(2,1),
\]

\[
R\leftrightarrow M_3
\]

through

\[
(1,1)\leftrightarrow(2,2),
\]

and

\[
R\leftrightarrow M_2
\]

through the two-hop physical route

\[
(1,1)
\leftrightarrow
(2,2)
\leftrightarrow
(3,4).
\]

Therefore every global query/response involving \(M_2\) is physically relayed by node \((2,2)\in M_3\).

However,

\[
\boxed{
\text{relay} \neq \text{data dependency}.
}
\]

Forwarding a message through \((2,2)\) does not by itself require the message to depend on \(x_{2,2}\).

---

# 6. Agreed cross-module interpretation

The physical edge

\[
(2,2)\to(3,4)
\]

crosses the \(M_3/M_2\) module boundary.

For this example, use the following interpretation:

\[
\boxed{
(2,2)\text{ is the physical parent of }(3,4)
\text{ only for communication/routing.}
}
\]

It is **not** the algorithmic hADMM parent of \(M_2\).

Therefore the ordinary hADMM cross-edge query

\[
q_{2,2,34}^{k+1}
=
A_{2,2,34}x_{2,2}^{k+1}
+
\frac1\rho y_{2,2,34}^{k}
-
c_{2,2,34}
\]

must not be inserted automatically into the mixed algorithm.

The global root is the logical fADMM coordinator of \(M_2\), even though its messages are physically routed through \((2,2)\).

This is the key distinction between the new algorithmic structure and the physical tree.

---

# 7. Local hADMM mathematics inside \(M_2\)

Inside \(M_2\), hADMM applies only to physical edges whose two endpoints both belong to \(M_2\):

\[
E_{M_2}=
\{
(3,4)\to(4,1),
(3,4)\to(4,2),
(3,4)\to(4,3),
(4,1)\to(5,1),
(4,1)\to(5,2),
(4,3)\to(5,3)
\}.
\]

For an internal edge

\[
(d,i)\to(d+1,j)\in E_{M_2},
\]

the hADMM parent query is

\[
q_{d,ij}^{k+1}
=
A_{d,ij}x_{d,i}^{k+1}
+
\frac{1}{\rho_h} y_{d,ij}^{k}
-
c_{d,ij}.
\]

A local leaf update is

\[
x_{d,i}^{k+1}
=
\arg\min_{x_{d,i}}
f_{d,i}(x_{d,i})
+
\frac{\rho_h}{2}
\left\|
B_{d-1,pi}x_{d,i}
+
q_{d-1,pi}^{k+1}
\right\|_2^2.
\]

A local non-leaf update is

\[
\begin{aligned}
x_{d,i}^{k+1}
=
\arg\min_{x_{d,i}}
\;&f_{d,i}(x_{d,i})\\
&+
\frac{\rho_h}{2}
\sum_{j\in C_{d,i}\cap M_2}
\left\|
A_{d,ij}x_{d,i}+q_{d,ij}^{k}
\right\|_2^2\\
&+
\frac{\rho_h}{2}
\left\|
B_{d-1,pi}x_{d,i}
+
q_{d-1,pi}^{k+1}
\right\|_2^2,
\end{aligned}
\]

for internal nodes whose parent is also inside \(M_2\).

The module root \((3,4)\) is special: its incoming quantity comes from the global fADMM module interface, not automatically from the physical parent \((2,2)\).

Section 10 defines that incoming query and root update. Node \((3,4)\) updates once, followed by one sweep through the remaining internal nodes, per global iteration.

---

# 8. Local fADMM mathematics inside \(M_1\) and \(M_3\)

\(M_1\) and \(M_3\) use fADMM locally.

For each module, its module root acts as the local star coordinator for its internal nodes, unless the existing standalone fADMM implementation uses an equivalent vectorized representation.

For \(M_1\),

\[
r_1=(2,1),
\]

with local targets

\[
(3,1),(3,2),(3,3).
\]

For \(M_3\),

\[
r_3=(2,2),
\]

with internal targets

\[
(3,5),(4,4),(4,5),
\]

using the existing fADMM query/response implementation.

Physical routing inside \(M_3\) is still required. For example, a logical local-fADMM query from \((2,2)\) to \((4,4)\) follows

\[
(2,2)\to(3,5)\to(4,4).
\]

The standalone fADMM solver is authoritative for the exact local query/prox/dual equations.

---

# 9. Global fADMM over modules

The global optimization topology is a star:

\[
R
\leftrightarrow
M_1,\qquad
R
\leftrightarrow
M_2,\qquad
R
\leftrightarrow
M_3.
\]

Conceptually, each module exposes a response operator

\[
\mathcal T_m,
\qquad m\in\{1,2,3\},
\]

and the global root sends a module query

\[
q_m^{k+1}
\]

and receives a module response

\[
\xi_m^{k+1}
=
\mathcal T_m(q_m^{k+1};\text{local state at }k).
\]

One global iteration must invoke each module exactly once.

The local module response is produced by:

\[
\mathcal T_1:\text{one local fADMM update/sweep},
\]

\[
\mathcal T_2:\text{one local hADMM update/sweep},
\]

\[
\mathcal T_3:\text{one local fADMM update/sweep}.
\]

No module is solved to convergence before returning.

---

# 10. Defined global-fADMM/M2 boundary interface

## 10.1 Exposed variable and global query

The global coordinator interacts with M2 through its module-root boundary variable:

\[
\boxed{\xi_2 := x_{3,4}.}
\]

Let \(z_2\) be the global coordinator's copy of this boundary variable. Adopt the scaled global dual \(u_2\) for the constraint \(\xi_2-z_2=0\), with global penalty \(\rho_g>0\). The query is

\[
\boxed{q_2^{k+1}=z_2^{k+1}-u_2^k.}
\]

The vectors \(z_2,q_2,u_2,\xi_2\) all have the dimension of \(x_{3,4}\); for the scalar fixture, this dimension is one.

The query contains no physical-parent term \(x_{2,2}^{k+1}\). It travels physically along

\[
(1,1)\rightarrow(2,2)\rightarrow(3,4),
\]

but `(2,2)` only forwards the packet. Forwarding neither transforms the query using M3 state nor triggers an M3 optimization update.

## 10.2 One root update and one internal hADMM sweep

Upon receiving \(q_2^{k+1}\), the module root performs

\[
\boxed{
\begin{aligned}
x_{3,4}^{k+1}=\arg\min_{x_{3,4}}\;&
f_{3,4}(x_{3,4})
+\frac{\rho_g}{2}\left\|x_{3,4}-q_2^{k+1}\right\|_2^2\\
&+\frac{\rho_h}{2}
\sum_{j\in C_{3,4}\cap M_2}
\left\|A_{3,4,j}x_{3,4}+\bar q_{3,4,j}^{k}\right\|_2^2.
\end{aligned}
}
\]

Here \(\rho_h>0\) is the internal hADMM penalty. The child contribution written as \(q_{3,4,j}^{k}\) in the proposed interface is denoted \(\bar q_{3,4,j}^{k}\) here to distinguish it from the forward parent query in Section 7:

\[
\bar q_{3,4,j}^{k}
=B_{3,4,j}x_{4,j}^{k}+\frac{y_{3,4,j}^{k}}{\rho_h}-c_{3,4,j}.
\]

The old child-query terms \(q_{d,ij}^{k}\) inside the non-leaf argmin of Section 7 have this same backward-message meaning. They are not the forward query \(A_{d,ij}x_{d,i}^{k+1}+y_{d,ij}^{k}/\rho_h-c_{d,ij}\).

After this single root update, the remaining nodes execute the ordinary one-pass hADMM updates of Section 7 over the six edges in \(E_{M_2}\). Each internal node updates once; child responses return on the backward traversal and each internal edge dual updates once:

\[
y_{d,ij}^{k+1}=y_{d,ij}^{k}+\rho_h
\left(A_{d,ij}x_{d,i}^{k+1}+B_{d,ij}x_{d+1,j}^{k+1}-c_{d,ij}\right).
\]

There is no second root update after the descendants finish and no local ADMM convergence loop. A numerical solver may iterate to solve a single local argmin.

M2 then returns the boundary value

\[
\boxed{\xi_2^{k+1}=x_{3,4}^{k+1}.}
\]

Thus the interface is

\[
\boxed{
q_2^{k+1}\longrightarrow x_{3,4}^{k+1}
\longrightarrow\text{one internal hADMM sweep}
\longrightarrow\xi_2^{k+1}=x_{3,4}^{k+1}.
}
\]

## 10.3 Dual ownership and global update

The dual states are distinct:

- \(u_2\): scaled global fADMM dual between R and M2;
- \(y_{d,ij}\): unscaled hADMM duals only for \((d,i)\to(d+1,j)\in E_{M_2}\).

The physical edge `(2,2)->(3,4)` has **no hADMM interface dual** in the mixed algorithm. It is a communication edge, not a module-level optimization edge.

After the response reaches R along `(3,4)->(2,2)->(1,1)`, the global dual updates as

\[
\boxed{u_2^{k+1}=u_2^k+\xi_2^{k+1}-z_2^{k+1}.}
\]

The complete global-interface chain is

\[
\boxed{
z_2^{k+1}\rightarrow q_2^{k+1}=z_2^{k+1}-u_2^k
\rightarrow M_2\rightarrow\xi_2^{k+1}\rightarrow u_2^{k+1}.
}
\]

**Existing-code sign convention.** The standalone fADMM implementation stores a scaled dual with the opposite sign: its query is `z + dual` and its update is `dual += z - response`. Therefore reuse requires `dual = -u2`. Under this substitution the equations agree; do not combine the query from one convention with the dual update from the other. The corresponding coordinator penalty is \(\frac{\rho_g}{2}\|z_2-\xi_2^k-u_2^k\|_2^2\).

## 10.4 Stateful response, not an exact module prox

Let \(s_2^k\) denote the internal hADMM state, including the local primals and internal edge duals. Define

\[
\boxed{
\xi_2^{k+1}=\mathcal T_{M_2}^{\mathrm h}
\left(q_2^{k+1};s_2^k\right).
}
\]

The sweep also advances the internal state to \(s_2^{k+1}\). Because only one hADMM sweep is performed, this response is not necessarily the exact module proximal operator \(\operatorname{prox}_{\lambda F_{M_2}}(q_2)\). The global method is fADMM with a module response generated by one hADMM iteration. Exact-prox ADMM convergence conclusions must not be assumed for this stateful response without a separate derivation or appropriate verification.

The former M2 query/root/response ambiguity is resolved by this section. Section 10.5 specifies the approved scalar global coordinator and M1/M3 completion.

---

## 10.5 Approved scalar completion for all modules

For the scalar worked example, every module exposes \(\xi_m=x_{r_m}\).
Use one global scalar \(z\), with \(z_1=z_2=z_3=z\):

\[
z^{k+1}=\arg\min_z f_R(z)+\frac{\rho_g}{2}
\sum_{m=1}^3(z-\xi_m^k-u_m^k)^2
=\frac{a_Rb_R+\rho_g\sum_m(\xi_m^k+u_m^k)}{a_R+3\rho_g}.
\]

For every module, \(q_m^{k+1}=z^{k+1}-u_m^k\) and
\(u_m^{k+1}=u_m^k+\xi_m^{k+1}-z^{k+1}\).
M1 and M3 use scalar consensus stars over their listed internal targets, with
scaled local duals \(v_{mt}\) oriented root-minus-target. Their root update is

\[
x_{r_m}^{k+1}=\arg\min_x f_{r_m}(x)
+\frac{\rho_g}{2}(x-q_m^{k+1})^2
+\frac{\rho_f}{2}\sum_t(x-x_t^k+v_{mt}^k)^2.
\]

Each target then updates once with query \(x_{r_m}^{k+1}+v_{mt}^k\),
minimizing its local objective plus the scalar proximal penalty; afterwards
\(v_{mt}^{k+1}=v_{mt}^k+x_{r_m}^{k+1}-x_t^{k+1}\).
Return the root value without updating it a second time. Use
\(\rho_g=\rho_h=\rho_f=1\) for the deterministic fixture.

The scalar reference enforces global root/module-root consensus, M1/M3 local
star consensus, and M2 internal-edge consensus. These imply a common value
\(x^*=\sum_v a_vb_v/\sum_v a_v\), with
\(J^*=\sum_v a_v(x^*-b_v)^2/2\). Report this experiment separately from the
unchanged vector-copy disjoint comparisons; do not compare their objective gaps.

---

# 11. Scalar quadratic numerical fixture

For deterministic implementation tests, assign one scalar variable per physical node:

\[
x_{d,i}\in\mathbb R.
\]

Use consensus-like physical-edge constraints

\[
x_{d,i}-x_{d+1,j}=0
\]

when a pure-tree reference problem is required.

Use strongly convex quadratic local objectives

\[
f_{d,i}(x)
=
\frac12 a_{d,i}(x-b_{d,i})^2.
\]

Use the fixed coefficient table:

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

Use

```text
rho_g = 1.0
rho_h = 1.0
rho = 1.0  # pure-tree reference penalty
all primal initial states = 0
all applicable dual initial states = 0
random seed = 0
```

These coefficients are a reproducible engineering fixture, not source-derived data.

---

# 12. Event semantics for one global iteration

One mixed global `step()` conceptually consists of:

```text
GLOBAL_ROOT_UPDATE

GLOBAL_QUERY(M1)
ROUTE_TO_M1
LOCAL_STEP(M1, fadmm)
GLOBAL_RESPONSE(M1)

GLOBAL_QUERY(M3)
ROUTE_TO_M3
LOCAL_STEP(M3, fadmm)
GLOBAL_RESPONSE(M3)

GLOBAL_QUERY(M2)
ROUTE R -> (2,2) -> (3,4)
LOCAL_STEP(M2, hadmm)
GLOBAL_RESPONSE(M2)
ROUTE (3,4) -> (2,2) -> R

GLOBAL_DUAL_UPDATE
```

The relative order of independent peer-module queries may differ if the global fADMM equations allow it.

The route to \(M_2\) must not trigger a local `M3` optimization update merely because `(2,2)` forwards the packet.

---

# 13. Residuals and correctness claims

For internal hADMM edges, use

\[
r_{d,ij}^{k+1}
=
A_{d,ij}x_{d,i}^{k+1}
+
B_{d,ij}x_{d+1,j}^{k+1}
-
c_{d,ij}.
\]

For the global M2 interface, the primal residual is \(r_2^{k+1}=\xi_2^{k+1}-z_2^{k+1}\), consistent with Section 10. For other local/global fADMM components, use the residual definition of the existing compatible implementation and document its sign convention.

Do not compare unlike residual definitions numerically.

For the scalar mixed implementation, use the standalone solvers' termination
rule: `max(max_primal_residual, max_dual_residual) < tol`, using the same `tol`
passed from `main.jl`, with `max_iter=1000` as the cap. The primal residual
covers global module boundaries and local optimization edges. The dual
residual is the maximum of \(\rho_g|\xi_m^{k+1}-\xi_m^k|\) over modules and
\(\rho_{\mathrm{local}}|x_t^{k+1}-x_t^k|\) over local targets/children, matching
the standalone per-edge/per-target residual convention. All penalties in
this scalar fixture equal one. Primal change remains a diagnostic only.
Report the actual iteration count and stopping reason. Applying this residual
criterion does not establish a general convergence theorem for the stateful
mixed method.

The M2 interface is now defined, but this alone is not a convergence proof. Do not claim mixed convergence, exact-prox equivalence, or KKT correctness without a separate derivation or appropriate numerical verification of the implemented method.

---

# 14. Partition interpretation for this example

The right branch is initially viewed as one candidate module \(M_{23}\).

It is split into \(M_2\) and \(M_3\) because \(M_2\) forms a large internally coupled subtree attached to the rest of the right branch through a small interface at the cut

\[
(2,2)\to(3,4).
\]

This follows the optimization-separator principle:

\[
\Phi(S)
=
\frac{
\operatorname{rank}(A_{\partial S})
}{
\sum_{i\in S}\dim(x_i)
}.
\]

A low local value of \(\Phi(S)\) indicates that a large internal subproblem communicates with the rest of the network through a relatively low-dimensional interface.

For this worked example, the partition is fixed in advance; automatic partition discovery is not part of the current implementation.
