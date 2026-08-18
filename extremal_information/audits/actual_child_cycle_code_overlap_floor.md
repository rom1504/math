# Cycle-code rigidity of every actual-child bridge likelihood

Status: **rigorous actual-channel identity with a sharp scope limitation**.
Every bridge pressure induced by the actual children has mandatory Eulerian
Fourier coefficients.  This rules out the pure full-parity landscape used as
the generic high-row-order falsifier.  It also gives an exact
likelihood-weighted cavity floor.  The weight is a collision likelihood,
however, so the result does not lower-bound the raw negative-tilt overlap and
does not reset `L_raw-negative-overlap` by itself.

## 1. Normalized actual-child likelihood

Fix any two children `A,D`, either relative orientation, and put

```math
t={\beta\over\sqrt N},\qquad \rho=\tanh t,\qquad d=mn.
```

Writing `I(x,y)=H_A(x)+epsilon H_D(y)`, their bridge partition is

```math
Z(B)=E_{x,y}\cosh(tI(x,y))\cosh(t x^TBy).
```

Let `Z_0=E_(x,y) cosh(tI(x,y))`, bias `(x,y)` proportionally to
`cosh(tI(x,y))`, take an independent fair `tau`, and set

```math
Q_{ij}=\tau x_i y_j.
```

Since `E_U exp(t sum_e B_eQ_e)=(cosh t)^d`, the normalized output
likelihood is

```math
p(B)={Z(B)\over Z_0(\cosh t)^d}
    =E_Q\prod_e(1+\rho B_eQ_e),
\qquad E_Up=1.                                      \tag{CY.1}
```

This is an exact representation of the actual bridge pressure up to an
additive constant.  No child optimality or surrogate construction is used.

## 2. Mandatory Eulerian Fourier coefficients

For an edge set `S` of `K_(m,n)`, (CY.1) gives

```math
\widehat p(S)=\rho^{|S|}E_Q\prod_{e\in S}Q_e.       \tag{CY.2}
```

If `S` is Eulerian, every left and right vertex has even degree.  Hence
`|S|` is even and, pointwise in the planted state,

```math
\prod_{(i,j)\in S}Q_{ij}
=\tau^{|S|}\prod_i x_i^{\deg_S(i)}
                 \prod_j y_j^{\deg_S(j)}=1.
```

Therefore every actual child pair satisfies

```math
\boxed{\widehat p(S)=\rho^{|S|}
       \quad\hbox{for every Eulerian }S.}           \tag{CY.3}
```

In particular every four-cycle has coefficient `rho^4`, independently of
the children and their orientation.  If `min(m,n)>=2` and `mn>4`, a pure
full-parity likelihood

```math
p_{\rm par}(B)=1+a\prod_eB_e
```

has zero coefficient on a proper four-cycle and therefore cannot be an
actual-child likelihood.  Thus the pure-parity falsifier in Theorem 37.51 is
genuinely excluded by the planted rank-one algebra, rather than merely by an
informal assertion that it is nonphysical.

## 3. Exact collision-weighted cavity identity

For bridge edge `e`, let

```math
p_{-e}(B_{-e})=E_{B_e}p(B)
```

and let `r_e(B_-e)` be the exact deleted-edge cavity response.  Expanding the
last channel factor in (CY.1) gives

```math
p(B)=p_{-e}(B_{-e})\{1+\rho B_er_e(B_{-e})\}.       \tag{CY.4}
```

For the half-flip derivative `D_e`,

```math
D_ep=\rho B_ep_{-e}r_e.
```

Parseval therefore proves the exact sum rule

```math
\boxed{
\rho^2\sum_eE_U[p_{-e}^2r_e^2]
=\sum_S|S|\widehat p(S)^2
\ge\sum_{S\ {m Eulerian}}|S|\rho^{2|S|}.}        \tag{CY.5}
```

Keeping only four-cycles yields the explicit consequence

```math
\boxed{
{1\over mn}\sum_eE_U[p_{-e}^2r_e^2]
\ge(m-1)(n-1)\rho^6.}                             \tag{CY.6}
```

At a comparable split and `t=beta/sqrt(N)`, (CY.6) is of order `N^-1`.
The complete right side of (CY.5) is the derivative of the cycle-space
weight enumerator of `K_(m,n)` at activity `rho^2`.  Determining whether its
complete-bipartite Ising asymptotics can survive removal of the collision
weight is a separate, presently open step; no phase-transition assertion is
used here.

## 4. Why this does not settle the negative path

The target statistic is

```math
{1\over\lambda mn}\int_{-\lambda}^0
 E_{\widehat\Pi_s}\sum_er_e^2\,ds,
\qquad d\widehat\Pi_s\propto p^s dU.
```

Equation (CY.5) instead weights each deleted-edge response by `p_-e^2`.
That collision factor may concentrate on an exponentially small set and
cannot be removed using the known child-minimality bounds.  Mandatory
four-cycles rule out a *pure* parity response, but they do not exclude a
landscape containing both the cycle-code floor and a high-row-order
component.  Converting (CY.5), or the full cycle-space transition, into a
bound under `p^s` for `s<=0` is exactly a rare-event/renormalization problem.

The reproducible finite tests of the unweighted target are in
[`../experiments/actual_child_negative_overlap_exact.py`](../experiments/actual_child_negative_overlap_exact.py)
and
[`../experiments/actual_child_negative_overlap_sample.py`](../experiments/actual_child_negative_overlap_sample.py).
