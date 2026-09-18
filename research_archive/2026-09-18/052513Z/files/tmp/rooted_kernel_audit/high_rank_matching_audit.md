# High-rank partial-matching moving projections: exact result and ceiling

Date: 2026-08-15. Scratch research report; no project ledger or steering file
was modified.

## Verdict

There is an explicit rank-`n-1`, difference-of-matchings moving projection
whose complete (unsplit) Gram remainder is positive and whose transition
eigenvalue has the required `Theta(n^{-1/2})` scale.  Its cut-code root twirl
has an exact hafnian-square formula and a uniform parity floor.  Thus this is
a rigorous inequality outside both direct-matching Fourier support and scalar
partial transversals.

It does not produce a useful numerical bound.  Its root statistic is exactly
the scalar matching hafnian statistic multiplied by `(n-1)^2`, and every
tested nontrivial instance is vacuous.  More decisively, **every** construction
based on the add/delete partial-matching incidence graph, including arbitrary
common `S_n` irreducibles and multiplicity blocks, has transition coefficient

```math
\lambda\le \frac{4}{3\sqrt3}\frac{1+o(1)}{\sqrt n}
=0.769800\ldots\,n^{-1/2}.
```

It therefore cannot reach the conference coefficient `n^{-1/2}`, even before
paying any root-twirl loss.  Noncommuting perfect-matching association-scheme
blocks remain an exact open target for a lower-bound improvement below that
ceiling, but no closed uniform root theorem was found.

## 1. Matching layers and coordinate isometries

Let `E=binom(n,2)` and let `M_l` be the set of `l`-edge matchings of `K_n`.
Write

```math
d_l^+=\binom{n-2l}{2},\qquad
|M_l|=\binom n{2l}(2l-1)!!.
```

The unnormalized up map is

```math
(U_lf)(T)=\sum_{e\in T}f(T\setminus\{e\}),
\qquad T\in M_{l+1},
```

and its adjoint is the down map.  The coordinate maps

```math
C_{l,l+1}e_T={1\over\sqrt{l+1}}
 \sum_{e\in T}e_e\otimes e_{T\setminus e},
```

```math
C_{l+1,l}e_R={1\over\sqrt{d_l^+}}
 \sum_{e\cap V(R)=\varnothing}e_e\otimes e_{R\cup e}
```

are isometries.  At a common output layer their ranges are orthogonal: one
has `e` in the output matching and the other has `e` disjoint from it.  If
`z in {+-1}^E`, `D_z e_R=z_R e_R`, and
`ell_z=E^{-1/2}sum_e z_e e_e`, direct contraction gives

```math
C_{l,l+1}^*(\ell_z\otimes D_zf)
 ={1\over\sqrt{E(l+1)}}D_zU_lf,
```

and the analogous down identity has denominator `sqrt(E d_l^+)`.

These are the exact hypotheses used by the moving-projection proof; no
channel is bounded separately.

## 2. An explicit common standard module

Let

```math
W=\{h\in\mathbb R^n:\sum_i h_i=0\},\qquad \dim W=n-1.
```

For `R in M_l`, put `S=V(R)`, `s=2l`, and

```math
b_S=1_S-{s\over n}1,
```

so that `<b_S,h>=sum_{i in S}h_i`.  Define

```math
\gamma_l={|M_l|\,2l(n-2l)\over n(n-1)},
\qquad
(\phi_lh)(R)={\langle b_{V(R)},h\rangle\over\sqrt{\gamma_l}}.
```

Uniformity of the covered `2l`-set gives
`phi_l^*phi_l=I_W`.  Moreover every vertex of an `(l+1)`-matching occurs in
exactly `l` of its `l`-edge deletions.  Hence

```math
U_l\phi_l=s_l\phi_{l+1},
\qquad
s_l^2={l(n-2l-1)(n-2l-2)\over2}.
```

The standard representation occurs with multiplicity one in `R[M_l]`, so
adjointness also gives `D_{l+1}phi_{l+1}=s_l phi_l`; there is no omitted
multiplicity channel.

The associated symmetric path on layers `1,...,L` has edge

```math
c_l^{std}={s_l^2\over
 E\sqrt{(l+1)\binom{n-2l}{2}}}.
```

Let `lambda` and `v` be its Perron eigenvalue/vector, put
`w_l=sqrt(|M_l|)v_l`, `alpha_l=w_l/sum_j w_j`, and combine the transported
copies `D_z phi_l W` with amplitudes `sqrt(alpha_l)`.  Repeating equations
(24)--(27) of the OpenAI moving-projection construction gives rank-`n-1`
projections `P_z` and Hilbert--Schmidt vectors `Theta_z` satisfying exactly

```math
K(z,w)=\operatorname{tr}(P_zP_w)\ge0,
```

```math
\langle\Theta_z,\Theta_w\rangle_{HS}
=\bigl(\tau(zw)-\lambda\bigr)K(z,w).
```

Thus both `K` and the complete remainder `(tau-lambda)K` are Gram kernels.
This is the promised rigorous same-switch, difference-of-matchings
inequality at `Theta(n^{-1/2})` scale.

## 3. Exact cut twirl

For a signing `a`, and an even vertex set `S`, define its unsigned hafnian

```math
H_S(a)=\sum_{R\in PM(S)}\prod_{e\in R}a_e.
```

For a cut-code word `c_e=sigma x_i x_j`, one has

```math
c_R=\sigma^l\prod_{i\in V(R)}x_i.
```

Averaging over `sigma,x` kills cross layers and pairs of matchings with
different covered sets.  Since

```math
\|b_Sb_S^*\|_{HS}^2=\|b_S\|^4,
\qquad
\gamma_l={|M_l|\|b_S\|^2\over n-1},
```

the root mass is exactly

```math
{T_a\over|C_n^+|}
=(n-1)^2\sum_{l=1}^L{\alpha_l^2\over|M_l|^2}
\sum_{|S|=2l}H_S(a)^2.                 \tag{3.1}
```

Every `H_S(a)` is a sum of `(2l-1)!!` signs.  This number is odd, hence
`H_S(a)` is an odd nonzero integer.  Therefore

```math
{T_a\over|C_n^+|}\ge
(n-1)^2\sum_{l=1}^L{\alpha_l^2\binom n{2l}\over|M_l|^2}. \tag{3.2}
```

This is a uniform, algebraically closed root theorem.  It also exposes the
failure: (3.1) is precisely `(n-1)^2` times the scalar matching hafnian
twirl, so the extra projection rank did not create a new root channel.

## 4. Universal add/delete architecture ceiling

The incidence graph between `M_l` and `M_{l+1}` is biregular, with degrees
`d_l^+` and `l+1`.  Consequently

```math
\|U_l\|=\sqrt{(l+1)d_l^+},
```

with equality on the constant vectors.  Every transition block obtained by
restricting `U_l` to a common `S_n` irrep, including an arbitrary
multiplicity-space matrix, therefore has singular norm at most this value.
After the two coordinate normalizations, every symmetric representation-graph
edge block has norm at most

```math
c_l^{max}={\sqrt{(l+1)d_l^+}\over E}.
```

For a block tridiagonal operator `J`,

```math
|\langle v,Jv\rangle|
\le2\max_l\|J_{l,l+1}\|\sum_l\|v_l\|^2,
```

so `lambda_max(J)<=2 max_l c_l^max`.  Taking `l/n -> alpha` gives

```math
\sqrt n\,c_l^{max}\longrightarrow
\sqrt{2\alpha}(1-2\alpha).
```

The maximum occurs at `alpha=1/6` and equals `2/(3sqrt(3))`.  Hence

```math
\boxed{\lambda\le
 {4\over3\sqrt3}{1+o(1)\over\sqrt n}.} \tag{4.1}
```

This proof does not assume multiplicity one, scalar channels, radiality, or
separate positivity of mixed terms.  It closes the entire partial-matching
add/delete representation graph at a fixed `23.0%` shortfall from the
conference coefficient.

## 5. What remains open in noncommuting multiplicity blocks

For a general `S_n`-equivariant projection `P` on partial-matching layers,
the same cut twirl is

```math
{T_a\over|C_n^+|}
=\sum_l\sum_{|S|=2l}
h_S(a)^T\bigl(P_l[S]\circ P_l[S]\bigr)h_S(a),       \tag{5.1}
```

where `h_S(a)=(a_R)_{R in PM(S)}` and `circ` is Schur product.  Each Schur
block is positive semidefinite and lies in the perfect-matching association
scheme.  Formula (5.1) is genuinely nonconstant and retains all
alternating-even-cycle terms jointly.

No polynomially closed uniform lower bound for its nonconstant components
was obtained.  There is a structural reason that the naive spectral floor
cannot work: a fixed `S_n` irrep has dimension `exp(O(n))`, whereas a linear
matching fiber has

```math
(2l-1)!!=\exp(\Theta(n\log n))
```

states.  Thus `P_l[S] circ P_l[S]` has zero minimum eigenvalue for every
single-irrep construction at the relevant scale.  The only component forced
uniformly by the present argument is the constant component, through odd
hafnian parity, and that is exactly the scalar root collapse.

An escape would require a new theorem forcing edge-induced matching sign
vectors to have quantitative mass in selected nonconstant primitive modules.
Even if found, (4.1) limits this architecture to a possible lower-bound
improvement; it cannot prove the conference constant or convergence to it.

## 6. Reproducible finite audit

Checker:

`/home/math/quadra/tmp/rooted_kernel_audit/audit_matching_projection.py`

Output:

`/home/math/quadra/tmp/rooted_kernel_audit/audit_matching_projection.json`

The checker verifies the transition matrices, computes the code budget by
exact cut enumeration, computes all subset hafnians by integer dynamic
programming, and tests saved exact minimizer representatives.  Representative
standard-module results are:

| `n` | layers | `lambda` | normalized certificate |
|---:|---:|---:|---:|
| 5 | 1--2 | 0.040825 | -0.020412 |
| 6 | 1--2 | 0.057735 | -0.028868 |
| 7 | 1--2 | 0.063888 | -0.010648 |
| 8 | 1--2 | 0.065205 | -0.010868 |
| 8 | 1--3 | 0.082479 | -0.226337 |

These are exact finite evaluations up to floating-point diagonalization of
the displayed explicit path.  They are evidence of root loss, not the proof
of the architecture ceiling; (4.1) is fully analytic.
