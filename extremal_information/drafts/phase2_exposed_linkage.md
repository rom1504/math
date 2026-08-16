# Query-restricted deterministic synchronization

**Status:** independently audited theorem draft.  This strengthens the scope
of Theorem 9.2 without changing its constants: global linkage is replaced by
linkage between a declared query carrier, through a chosen local path
skeleton.  Its concise consequence is promoted as Corollary 9.4.

## 1. Result in one sentence

For a fixed zero-temperature query class, species synchronization is needed
only on states that can still win after the entire allowed coupling
oscillation is paid.  If those states are linked by paths with small local
cancellation and backtracking, the same uniform maximum estimate as Theorem
9.2 follows; neither linkage nor ultrametricity is needed elsewhere.

This is genuinely weaker when the base score has a robust exposed band.  It
does not help for a flat base score, where every rare fibre can be exposed.

## 2. Carrier and path-skeleton definitions

Retain the finite setup of Theorem 9.2.  Thus `E=binom(Omega,2)`, the species
profiles are `R_s:E->[0,1]`, the weights `lambda_s>0` sum to one, and

```math
q(e)=\sum_s\lambda_sR_s(e).                        \tag{EL.1}
```

Write

```math
\mathfrak c(e,f)=
\sum_s\lambda_s|R_s(e)-R_s(f)|-|q(e)-q(f)|.       \tag{EL.2}
```

Let `A` be a nonempty subset of pair labels.  It will be the carrier of all
states relevant to the declared queries.  Let `H` be any subgraph of the
pair line graph: an edge of `H` joins two pair labels that share one state.

Say that `H` has **local cancellation defect at most `zeta`** if

```math
\mathfrak c(e,f)\le\zeta
\qquad\text{for every }ef\in H.                   \tag{EL.3}
```

Say that `q` is **`(A,H;D,tau)`-linked** if, whenever `e,f in A` and
`q(e)<=q(f)`, there is an `H`-path

```math
e=e_0,e_1,\ldots,e_ell=f,
\qquad \ell\le D,                                  \tag{EL.4}
```

whose total downward variation obeys

```math
\sum_{j=0}^{\ell-1}(q(e_j)-q(e_{j+1}))_+\le\tau.  \tag{EL.5}
```

Only the endpoints must lie in `A`; intermediate labels may be auxiliary
states outside it.  Unlike global `(D,tau)`-linkage, no path obligation is
placed on a pair of endpoints that the declared query class cannot expose.

## 3. Restricted synchronization theorem

### Theorem EL.1 (carrier synchronization from a local skeleton)

Assume (EL.3)--(EL.5).  For every species `s`, there is a nondecreasing
`1/lambda_s`-Lipschitz function `L_s` on
`[min_(e in A)q(e),max_(e in A)q(e)]` such that

```math
\boxed{
\max_{e\in A}|R_s(e)-L_s(q(e))|
\le {\tau+D\zeta/2\over\lambda_s}.}               \tag{EL.6}
```

#### Proof

Take `e,f in A` with `q(e)<=q(f)` and an allowed path.  Put

```math
v_j=(R_s(e_{j+1})-R_s(e_j))_s,
\qquad
\Delta_j=q(e_{j+1})-q(e_j).                       \tag{EL.7}
```

With `||v||_lambda=sum_s lambda_s|v_s|`, definition (EL.2) gives on every
skeleton edge

```math
\|v_j\|_\lambda\le|\Delta_j|+\zeta.               \tag{EL.8}
```

The triangle inequality and nonnegative net change give

```math
\begin{aligned}
\mathfrak c(e,f)
&\le\sum_j\|v_j\|_\lambda-
       \left|\sum_j\Delta_j\right|\\
&\le\sum_j|\Delta_j|-\sum_j\Delta_j+D\zeta\\
&=2\sum_j(-\Delta_j)_++D\zeta\\
&\le2\tau+D\zeta.                                \tag{EL.9}
\end{aligned}
```

Thus the uniform cancellation defect over `A x A` is at most
`delta=2tau+D zeta`.  The isotonic-envelope proof of Proposition 9.1 uses
only comparisons between the points being approximated.  Applying it to the
restricted set `A`, namely defining

```math
L_s(p)=\inf_{f\in A}
\left\{R_s(f)+{(p-q(f))_+\over\lambda_s}\right\}, \tag{EL.10}
```

gives error `delta/(2lambda_s)`, which is (EL.6). `square`

### Corollary EL.2 (ultrametric certificate only on used triangles)

Suppose that for every skeleton edge

```math
\{x,y\}\;H\;\{x,z\},                              \tag{EL.11}
```

each individual species kernel and each sum of two distinct species obeys
the `eta`-ultrametric inequalities on the single state triple `{x,y,z}`.
Then the local no-crossing argument used in Theorem 9.2 applies on that edge
and gives

```math
\mathfrak c(\{x,y\},\{x,z\})\le6\eta.             \tag{EL.12}
```

Consequently Theorem EL.1 gives

```math
\max_{e\in A}|R_s(e)-L_s(q(e))|
\le {\tau+3D\eta\over\lambda_s}.                 \tag{EL.13}
```

Global ultrametricity is therefore unnecessary: only the triangles traversed
by the carrier-linking skeleton are used.

## 4. A noncircular zero-temperature carrier

Let `G:E->R` be a fixed base score.  Let `Q` be a declared class of coupling
potentials `Psi:[0,1]^S->R`.  Assume every `Psi in Q` has

```math
\operatorname{osc}(\Psi)
:=\sup_{r,r'}|\Psi(r)-\Psi(r')|\le B              \tag{EL.14}
```

and the common coordinatewise Lipschitz bound

```math
|\Psi(r)-\Psi(r')|
\le\sum_s\kappa_s|r_s-r_s'|.                     \tag{EL.15}
```

For example, (EL.15) on `[0,1]^S` implies (EL.14) with
`B=sum_s kappa_s`, though a smaller declared oscillation bound may be used.
Define the robust top band

```math
A_B(G)=\{e\in E:G(e)\ge\max_fG(f)-B\}.            \tag{EL.16}
```

This set depends only on the declared base query and its allowed coupling
budget, not on an unknown synchronization map or on which state actually
maximizes the coupled landscape.

### Corollary EL.3 (exposed-carrier zero-temperature control)

Assume `q` is `(A_B(G),H;D,tau)`-linked and (EL.3) holds.  For every
`Psi in Q`, put

```math
V_\Psi=\max_{e\in E}
\{G(e)+\Psi((R_s(e))_s)\},                         \tag{EL.17}
```

```math
\widetilde V_\Psi=\max_{e\in A_B(G)}
\{G(e)+\Psi((L_s(q(e)))_s)\}.                     \tag{EL.18}
```

Then

```math
\boxed{
|V_\Psi-\widetilde V_\Psi|
\le(\tau+D\zeta/2)
\sum_s{\kappa_s\over\lambda_s}.}                 \tag{EL.19}
```

Under the local ultrametric hypothesis of Corollary EL.2, replace
`D zeta/2` by `3D eta`.

#### Proof

Choose `e_0` maximizing `G`.  If `e` lies outside `A_B(G)`, then

```math
G(e)+\Psi(R(e))
<G(e_0)-B+\Psi(R(e))
\le G(e_0)+\Psi(R(e_0)).                           \tag{EL.20}
```

Thus every true maximum in (EL.17) lies in `A_B(G)`.  On this carrier,
Theorem EL.1 and (EL.15) bound the pointwise score error by the right side of
(EL.19).  The elementary inequality
`|max F-max Ftilde|<=||F-Ftilde||_infinity` on the same finite carrier proves
the claim. `square`

The restriction to `A_B(G)` in (EL.18) is essential.  Unlinked states that
cannot win the true query need not be accurately calibrated, so allowing
them back into the surrogate maximum could create a spurious winner.

## 5. Strictness: global linkage may fail arbitrarily far away

The weakening is strict even under exact ultrametricity of every nonnegative
species mixture.  Use the rare-matching construction from Proposition 5.1
with `m>=2`, weights `lambda_1=lambda_2=1/2`, and matching edges
`e_1,...,e_m`.  Choose a nonmatching edge `h` sharing one endpoint with
`e_1`.  Then

```math
(R_1(h),R_2(h),q(h))=(0,0,0),
```

```math
(R_1(e_1),R_2(e_1),q(e_1))=(\rho,0,\rho/2).       \tag{EL.21}
```

Set `A={h,e_1}` and take the skeleton consisting of their one adjacent
edge.  It has `zeta=0`, and `q` is `(A,H;1,0)`-linked.  Hence both species
synchronize exactly on `A`.

For any `B>=0` and `gamma>0`, define

```math
G(e)=\begin{cases}
0,&e\in A,\\
-B-\gamma,&e\notin A.
\end{cases}                                       \tag{EL.22}
```

Then `A_B(G)=A`, so Corollary EL.3 answers every coupling query of
oscillation at most `B` exactly.  Yet global exact monotone linkage fails:
the fibre `q=rho/2` consists of the mutually disjoint matching edges, so no
constant-`q` line-graph path joins `e_1` to `e_2`.

As `m` grows, the carrier and its synchronization state remain of size two
while the uncontrolled pair landscape grows quadratically.  This proves
that exposed-carrier linkage is not global linkage in disguise.

## 6. Scope and director judgment

The theorem is useful only when the query apparatus supplies a small robust
carrier.  If `G` is flat and the coupling class can expose any species
direction, `A_B(G)=E`; the rare-matching obstruction returns and the result
reduces to global synchronization.  Similarly, defining `A` retrospectively
as the unknown true argmax would be circular; (EL.16) avoids that by using
only `G` and a declared oscillation budget.

What has been gained is exact and finite:

1. path linkage is required only between potentially winning endpoints;
2. no-crossing/ultrametric inequalities are required only on triangles in a
   chosen linking skeleton; and
3. all unexposed states may have arbitrary synchronization error without
   affecting the declared zero-temperature responses.

This is a genuine query-relative strengthening of Theorem 9.2, but not a
general deterministic Parisi mechanism.  The strongest next question is
whether a natural model supplies a small robust carrier without an external
base gap--for example, from an exposed normal cone or a deterministic margin
theorem.  Without such a mechanism, the result remains a useful finite
localization theorem rather than a new global order parameter.
