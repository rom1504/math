# Independent audit: query-restricted deterministic synchronization

**Scope.** This report audits `phase2_exposed_linkage.md`, reconstructing
EL.1--EL.3 and the strictness example. It makes no surface-file edits.

## Verdict

**Accept.** The three constants are correct, the robust carrier is logically
noncircular, and the rare-matching example proves that carrier linkage is
strictly weaker than global linkage.

The result merits **concise promotion as one corollary of the existing finite
synchronization theorem**, not three new headline theorems. EL.1 is the
restricted isotonic-envelope argument, EL.2 is the already-audited local
no-crossing certificate, and EL.3 is their genuinely useful query-relative
consequence. Promotion should retain one caveat: the carrier is defined
without the unknown synchronization map, but it need not be cheap to compute
or describe from a large base-score table.

## 1. Reconstruction of EL.1

For one skeleton increment write

```math
v_j=(R_s(e_{j+1})-R_s(e_j))_s,
\qquad
\Delta_j=q(e_{j+1})-q(e_j).
```

The definition of the cancellation defect is exactly

```math
\mathfrak c(e_j,e_{j+1})
=\|v_j\|_\lambda-|\Delta_j|,
```

so a skeleton edge gives

```math
\|v_j\|_\lambda\le |\Delta_j|+\zeta.
```

For carrier endpoints oriented so that `q(e)<=q(f)`, the path has nonnegative
net increment. Hence

```math
\begin{aligned}
\mathfrak c(e,f)
&\le \sum_j\|v_j\|_\lambda-
       \left|\sum_j\Delta_j\right|\\
&\le \sum_j|\Delta_j|-\sum_j\Delta_j+\ell\zeta\\
&=2\sum_j(-\Delta_j)_++\ell\zeta\\
&\le 2\tau+D\zeta.
\end{aligned}
```

Thus the restricted carrier has uniform cancellation defect

```math
\delta=2\tau+D\zeta.
```

The isotonic envelope

```math
L_s(p)=\inf_{f\in A}
\left\{R_s(f)+{(p-q(f))_+\over\lambda_s}\right\}
```

uses only carrier comparisons. The usual two-point inequality gives, at
every attained `q(e)`,

```math
R_s(e)-{\delta\over2\lambda_s}
\le L_s(q(e))\le R_s(e).
```

It is nondecreasing and `1/lambda_s`-Lipschitz because every function in the
infimum has those same two properties. Substituting `delta` proves exactly

```math
\max_{e\in A}|R_s(e)-L_s(q(e))|
\le {\tau+D\zeta/2\over\lambda_s}.
```

There is no missing factor from a path of length below `D`: the proof uses
`ell zeta<=D zeta`.

## 2. Reconstruction of EL.2

On one skeleton edge in the pair line graph, the earlier three-state
no-crossing proof uses only the three orientations of the ultrametric
inequality for that one state triple, for each species and each distinct
pair-sum. It gives

```math
\mathfrak c(e,f)\le6\eta.
```

No triangle outside the skeleton is invoked. Taking `zeta=6 eta` in EL.1
gives

```math
{\tau+D\zeta/2\over\lambda_s}
={\tau+3D\eta\over\lambda_s},
```

so EL.13 has the correct constant. For maximal clarity, “obeys the
ultrametric inequalities on the triple” should continue to mean all required
permutations on that triple, not just one displayed orientation.

## 3. Reconstruction and noncircularity of EL.3

Let `e_0` maximize `G`. For `e` outside

```math
A_B(G)=\{e:G(e)\ge\max G-B\},
```

one has `G(e)<G(e_0)-B`. The oscillation bound implies

```math
\Psi(R(e))\le\Psi(R(e_0))+B,
```

and therefore

```math
G(e)+\Psi(R(e))<G(e_0)+\Psi(R(e_0)).
```

Every true maximizer lies in the carrier, uniformly for the declared query
class. On the carrier, EL.1 and coordinatewise Lipschitzness give pointwise
error at most

```math
(\tau+D\zeta/2)
\sum_s{\kappa_s\over\lambda_s}.
```

Taking maxima over the same carrier preserves this bound, proving EL.19.
Restricting the surrogate maximum to the carrier is essential and is done
correctly.

The construction is noncircular in the mathematical sense relevant here:
`A_B(G)` depends only on the fixed base score and a declared uniform
oscillation budget. It does not use the unknown maps `L_s`, an unknown
coupled optimizer, or the particular `Psi` later selected. It is not
automatically an information compression result, because listing or finding
the top band of an arbitrary `G` may itself require the whole base landscape.
The draft already limits its claim to models with a small or succinct robust
carrier.

## 4. Strictness example

In the rare-matching construction, choose a nonmatching edge `h` adjacent to
the exceptional matching edge `e_1`. Their profiles are

```math
R(h)=(0,0),\quad q(h)=0,
\qquad
R(e_1)=(\rho,0),\quad q(e_1)=\rho/2.
```

With equal species weights, the one skeleton edge has

```math
\mathfrak c(h,e_1)
={\rho\over2}-{\rho\over2}=0.
```

It is a length-one nondecreasing path, so `(A,H;1,0)` linkage and exact
carrier synchronization are valid. For

```math
G=0\text{ on }A,\qquad G=-B-\gamma\text{ off }A,
```

the robust band is exactly `A`, including when `B=0`.

Global exact monotone linkage fails. All matching edges have the same total
overlap `rho/2`, and they are mutually disjoint in the underlying complete
graph. A zero-backtracking path between two equal-`q` endpoints would have to
remain entirely in that exact fibre, whose pair-line subgraph is edgeless.
This proves strictness for every path-length bound, while every nonnegative
species mixture remains exactly ultrametric as in the audited matching
construction.

## 5. Promotion recommendation

Promote only a compact statement of EL.3 immediately after the global finite
synchronization theorem:

> Global linkage and global ultrametricity may be restricted to a robust
> top-score carrier and a local linking skeleton; the same zero-temperature
> error bound holds with `tau+D zeta/2`, and the matching example makes this
> restriction strict.

Keep the full EL.1--EL.2 proof decomposition and carrier strictness details in
the draft. This adds a genuine query-relative reduction without making the
surface look like three independent synchronization breakthroughs. It should
not be advertised as a natural-model synchronization mechanism until a model
supplies a succinct robust carrier and verified skeleton constants.
