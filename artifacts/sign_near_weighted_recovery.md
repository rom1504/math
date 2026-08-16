# Globally sign-near weighted recovery rounds to exact sign recovery

Date: 2026-08-16.

Status: **verified rounding theorem**, independently audited twice.  Global
sign-nearness permits exact sign rounding with vanishing normalized
`L^infinity -> L^1`, action-profile, and Boolean-objective error.  This removes
the final integrality substep, but it also makes sign-near weighted recovery
equivalent to exact recovery as an existential recovery obligation.  The note
does not construct a recovery sequence.

## 1. Scalar rounding at the exact objective scale

For a symmetric hollow \(W=(w_{ij})\in[-1,1]^{n\times n}\), define

```math
Q(W):=\max_{x\in\{\pm1\}^n}
\left|\sum_{i<j}w_{ij}x_ix_j\right|,
\qquad
V(W):=\sum_{i<j}(1-w_{ij}^2).
```

> **Theorem (global sign-near objective rounding).** There is a symmetric
> hollow signing \(A\) such that
>
> ```math
> \boxed{
> Q(A)\le Q(W)+C\bigl(\sqrt{nV(W)}+n\bigr)}                \tag{WR.S1}
> ```
>
> for a universal constant \(C\).  Hence \(V(W)=o(n^2)\) gives
> \(Q(A)\le Q(W)+o(n^{3/2})\).

**Proof.** Independently round each upper-triangular entry with

```math
\mathbb P(A_{ij}=1)=\frac{1+w_{ij}}2,
\qquad A_{ji}=A_{ij},
\qquad A_{ii}=0.
```

For a fixed Boolean spin \(x\), the error

```math
H_{A-W}(x)=\sum_{i<j}(A_{ij}-w_{ij})x_ix_j
```

is a sum of independent centered variables, each bounded by two, with total
variance \(V(W)\).  Bernstein's inequality and a union bound give

```math
\mathbb P\!\left(\max_x|H_{A-W}(x)|\ge t\right)
\le 2^{n+1}\exp\!\left[-\frac{t^2}{2(V(W)+2t/3)}\right].   \tag{WR.S2}
```

Taking \(t=C(\sqrt{nV(W)}+n)\) with a sufficiently large absolute constant
makes the right side smaller than one.  Some supported outcome then satisfies
(WR.S1).  \(\square\)

## 1A. Bilinear rounding preserves every fixed action profile

For a symmetric matrix `E`, define its Boolean bilinear norm

```math
B(E):=\max_{x,y\in\{\pm1\}^n}|x^{\mathsf T}Ey|.
```

For the same biased rounding `E=A-W`, a fixed pair `(x,y)` gives a sum over
independent unordered edges with coefficient
`x_i y_j+x_j y_i in {-2,0,2}`.  Its variance is at most `4V(W)` and each
summand has absolute value at most four.  Bernstein and a union bound over
the `4^n` Boolean pairs give

```math
\mathbb P\{B(A-W)\ge t\}
\le2\,4^n\exp\left[-\frac{t^2}{2(4V(W)+4t/3)}\right].    \tag{WR.B1}
```

Consequently some supported exact signing satisfies

```math
\boxed{B(A-W)\le C\bigl(\sqrt{nV(W)}+n\bigr).}           \tag{WR.B2}
```

This has the exact operator interpretation

```math
\sup_{|f|\le1}\|(A-W)f\|_{\ell^1}=B(A-W),
\qquad
\|T_A-T_W\|_{L^\infty\to L^1}
=\frac{B(A-W)}{n^{3/2}}.                                 \tag{WR.B3}
```

The first identity follows by dualizing the `ell^1` norm and then using
multi-affinity to move both cube variables to Boolean vertices.  If
`V(W)=o(n^2)`, write the last quantity in (WR.B3) as `r_n=o(1)`.  For any
fixed `k` and any bounded inputs `f_1,...,f_k`, couple the two joint profile
laws by using the same uniform vertex and the same inputs.  Markov gives

```math
d_H\bigl(\mathcal S_k(T_A),\mathcal S_k(T_W)\bigr)
\le\min\{1,\sqrt{k r_n}\}.                               \tag{WR.B4}
```

Thus the standard action metric obeys

```math
d_M(T_A,T_W)
\le\left(\sum_{k\ge1}2^{-k}\sqrt k\right)\sqrt{r_n}
=o(1).                                                    \tag{WR.B5}
```

No row-variance deletion and no operator-norm bound on the rounded matrix are
needed for (WR.B4)--(WR.B5).  Directly,

```math
|\Phi(T_A)-\Phi(T_W)|\le r_n,                            \tag{WR.B6}
```

so continuity is not invoked after rounding.  The two-spin proof has only a
finite Bernstein-constant loss compared with (WR.S2), never a leading
`n^(3/2)` loss under `V=o(n^2)`.  An independent referee checked every
coefficient, variance, normalization, and action-metric tail; see
`ar_bilinear_rounding_referee.md`.

## 2. Spectral rounding theorem

Let \(W=(w_{ij})\) be a symmetric hollow order-\(n\) matrix with entries in
\([-1,1]\), and define its maximum row fractional variance by

```math
v(W):=\max_i\sum_{j\ne i}(1-w_{ij}^2).
```

> **Theorem (sign-near spectral rounding).** There is a symmetric hollow
> signing \(A\) such that
>
> ```math
> \boxed{
> \|A-W\|_{op}
> \le C\bigl(\sqrt{v(W)}+\sqrt{\log n}\bigr)}              \tag{WR.1}
> ```
>
> for a universal constant \(C\).  In particular, if \(v(W)=o(n)\), then
> \(\|A-W\|_{op}=o(\sqrt n)\).

### Proof

Independently for \(i<j\), choose

```math
\mathbb P(A_{ij}=1)=\frac{1+w_{ij}}2,
\qquad A_{ji}=A_{ij},
\qquad A_{ii}=0.
```

Then \(X=A-W\) is centered, its upper-triangular entries are independent,
\(|X_{ij}|\le2\), and

```math
\mathbb E X_{ij}^2=1-w_{ij}^2.
```

Let \(X'\) be an independent copy.  Conditional Jensen gives the standard
symmetrization inequality

```math
\mathbb E\|X\|_{op}\le\mathbb E\|X-X'\|_{op}.             \tag{WR.2}
```

The entries of \(Y=X-X'\) are symmetric, independent above the diagonal,
bounded by two, and have maximum row variance at most \(2v(W)\).
Corollary 3.6 of Bandeira--van Handel, *Sharp nonasymptotic bounds on the
norm of random matrices with independent entries*, Ann. Probab. 44 (2016),
[arXiv:1408.6185](https://arxiv.org/abs/1408.6185), applied with a fixed
\(\alpha\ge3\), gives

```math
\mathbb E\|Y\|_{op}
\le C\bigl(\sqrt{v(W)}+\sqrt{\log n}\bigr).               \tag{WR.3}
```

Explicitly, set
\(b_{ij}=(\mathbb E Y_{ij}^2)^{1/2}\) and
\(\xi_{ij}=Y_{ij}/b_{ij}\), using an auxiliary Rademacher when \(b_{ij}=0\).
Then \(\xi_{ij}\) is symmetric with unit variance, and Corollary 3.6 with
\(p=2\lceil\alpha\log n\rceil\) gives the source's explicit bound

```math
e^{2/\alpha}\left(2\sqrt{2v(W)}+28\alpha\sqrt{\log n}\right).
```

Here \(|Y_{ij}|\le2\), so its high-moment entry parameter is at most two.
The case \(n=1\) is trivial.  Equations
(WR.2)--(WR.3) imply that at least one outcome satisfies (WR.1).  The outcome
is automatically symmetric, hollow, and exactly sign-valued off the
diagonal.  \(\square\)

Recent sharp matrix-concentration results improve the leading constant and
lower-order terms, but are unnecessary for the little-\(o\) conclusion.

## 3. Consequence for the Boolean objective and action profile

On the same uniform \(n\)-point probability space, put
\(T_W=W/\sqrt n\) and \(T_A=A/\sqrt n\).  For every \(|f|\le1\),

```math
|\langle f,(T_A-T_W)f\rangle|
\le \|T_A-T_W\|_{2\to2}
=\frac{\|A-W\|_{op}}{\sqrt n}.                            \tag{WR.4}
```

Consequently, when \(v(W)=o(n)\),

```math
|\Phi(T_A)-\Phi(T_W)|=o(1).                               \tag{WR.5}
```

If \(\|T_W\|_{2\to2}=O(1)\), the rounded operators retain the same bound.
Moreover, under the identity coupling,

```math
\|(T_A-T_W)f\|_2=o(1)
```

uniformly over \(|f|\le1\).  More quantitatively, if this norm is at most
\(\epsilon\), Markov's inequality at threshold \(\epsilon^{2/3}\) and the
identity coupling give Levy--Prokhorov distance at most
\(\epsilon^{2/3}\).  Thus every one-profile of \(T_A\) is
Levy--Prokhorov \(o(1)\)-close to the corresponding profile of \(T_W\).
Thus sign-near weighted directed recovery implies exact-sign directed
recovery, not merely objective recovery.

## 4. Removing exceptional rows

The maximum-row hypothesis in Section 2 is not needed at the input.
Suppose \(V(W_n)=o(n^2)\).  Put

```math
\epsilon_n:=\frac{2V(W_n)}{n^2}=o(1),
\qquad
\tau_n:=\sqrt{\epsilon_n}.
```

Delete every row \(i\) for which

```math
\sum_{j\ne i}(1-w_{ij}^2)>\tau_n n.
```

Since the sum of all row variances is \(2V(W_n)=\epsilon_n n^2\), at
most \(\tau_n n=o(n)\) rows are deleted.  The retained principal matrix
\(U_m\), \(m=n-o(n)\), satisfies

```math
v(U_m)\le\tau_n n=o(m).                                  \tag{WR.6}
```

Deletion is harmless for the weighted Boolean objective.  For a spin on the
retained set, extend it by independent signs on the deleted vertices; all
removed terms have mean zero.  Hence

```math
Q(U_m)\le Q(W_n).                                         \tag{WR.7}
```

It is also harmless for directed action recovery under a common normalized
operator bound.  Extend a retained test function by zero.  With probability
\(m/n=1-o(1)\), a uniform point lies in the retained set and the two outputs
differ only by the normalization factor \(\sqrt{n/m}=1+o(1)\); the exceptional
points have mass \(o(1)\).  The operator bound makes the matched output
difference \(o(1)\) in \(L^2\).

More explicitly, if \(g\) is the zero extension of \(f\) and \(R\) is the
retained set, then

```math
(T_Wg)|_R=\sqrt{m/n}\,T_Uf,
\qquad
\|T_U\|_{2\to2}
\le\sqrt{n/m}\,\|T_W\|_{2\to2}.                          \tag{WR.8}
```

Couple a uniform retained point with the same point in the original space;
the coupling fails only on the deleted mass.  The displayed identity and
operator bound control the normalization error on the matched event.
Therefore

```math
\partial_1(T_{U_m},T_{W_n})=o(1).                         \tag{WR.8a}
```

Applying Section 2 to \(U_m\) produces exact signings at orders
\(m=n-o(n)\).  If a prescribed target order falls between the original and
retained sizes, principal deletion or the random \(o(n)\)-vertex completion
from `minimal_all_order_action_recovery.md` transfers the objective with
\(o(n^{3/2})\) loss.

## 5. Weighted recovery implications

For convergence alone, exact-sign recovery can be replaced by the following
strictly relaxed realization condition.

> **Objective sign-near weighted recovery.** For each member of a null
> sequence of purification tolerances, one selected liminf cluster \(T\) has
> symmetric hollow weighted matrices \(W_m\in[-1,1]^{m\times m}\) on an
> upward ratio-dense set of orders such that
>
> ```math
> \Phi(T_{W_m})\le\Phi(T)+o(1),
> \qquad
> V(W_m)=o(m^2).                                           \tag{WR.9}
> ```

For a weighted hollow matrix, separate affinity in every coordinate implies

```math
\Phi(T_{W_m})=\frac{2Q(W_m)}{m^{3/2}}.                    \tag{WR.9a}
```

Indeed, the extrema of the quadratic form over \([-1,1]^m\) occur at Boolean
vertices.  Thus the first condition in (WR.9) is exactly the required
weighted Boolean bound.

Section 1 rounds these matrices at the same orders with \(o(m^{3/2})\)
objective loss.  Principal deletion then proves convergence exactly as in
`minimal_all_order_action_recovery.md`.

If preservation with an operator-norm bound on the final exact signing is
desired, one may use the following stronger form and the spectral route:

> **Profile sign-near weighted recovery.** At the same orders, require
>
> ```math
> \|T_{W_m}\|_{2\to2}=O(1),
> \qquad
> \partial_1(T_{W_m},T)\to0,
> \qquad
> V(W_m)=o(m^2).                                           \tag{WR.10}
> ```

Delete exceptional rows as in Section 4 and apply spectral rounding.  This
produces exact signings at orders \(m-o(m)\) with normalized operator error
\(o(1)\).  Directed-profile distances compose:

```math
\partial_1(T_A,T)
\le \partial_1(T_A,T_U)
   +\partial_1(T_U,T_W)
   +\partial_1(T_W,T)=o(1).                               \tag{WR.11}
```

The perturbed order set is still upward ratio-dense.  Indeed, write
\(k(m)=m-o(m)\) for the retained order and

```math
\delta_N:=\sup_{\substack{m\in\mathcal N\\m\ge N}}
\left(1-\frac{k(m)}m\right)=o(1).
```

For a target \(N\), take the first original recovery order
\(m\ge N/(1-2\delta_N)\).  Upward ratio-density gives \(m/N=1+o(1)\), while
\(k(m)\ge m(1-\delta_N)\ge N\) and \(k(m)/N=1+o(1)\).  Thus the rounded
sequence itself satisfies \(\mathrm{AR}_{\min}^{\to}\); no padding is needed
for this profile claim.  In particular, (WR.10) proves convergence.  If only
action-profile and objective recovery are needed, Section 1A rounds directly
at the original orders and makes the deletion step unnecessary.

## 6. Strictness and boundary

This is a syntactic relaxation of the **exact-sign realization constraint**:
the entries of \(W_m\) may all be fractional, provided their total fractional
variance is \(o(m^2)\).  However, it is **not a strict existential recovery
reduction**.  Exact recovery implies weighted recovery by taking `W=A` and
`V=0`; conversely (WR.B2)--(WR.B6) turn any sign-near weighted recovery into
an exact recovery at the same orders with vanishing action and objective
error.  The theorem is a useful rounding module, but constructing its weighted
input has not removed an asymptotic recovery obligation.

The information boundary is equally sharp.  For every fixed `epsilon>0`,
`V(W)=o(n^2)` implies that all but `o(n^2)` entries satisfy
`|w_ij|>1-epsilon`; the signs of those entries reveal an almost complete sign
skeleton.  Thus a sign-near weighted state is not inherently a compressed
description of the target signing.

The variance threshold precisely excludes the archived naive blow-ups.  A
block bias of order \(k^{-1/2}\) leaves \(1-w_{ij}^2=\Theta(1)\) on a positive
fraction of entries, so \(V(W)=\Theta(n^2)\) and \(v(W)=\Theta(n)\); both
rounding estimates then allow a leading residual: \(\Theta(\sqrt n)\) in
operator norm and \(\Theta(n^{3/2})\) in the Boolean objective.  The zero
weighted matrix is the sharpest simple
falsifier: every sign rounding has

```math
\|A\|_{op}\ge\frac{\|A\|_F}{\sqrt n}=\sqrt{n-1}.
```

There is also an entropy boundary for microcanonical attempts.  Let \(\mu\)
be any law on the \(N=\binom n2\) edge signs, let \(U\) be the uniform product
law, and put \(w_e=\mathbb E_\mu A_e\).  Subadditivity of entropy and Pinsker's
inequality give, in natural logarithms,

```math
\begin{aligned}
D(\mu\|U)
&=N\log2-H(\mu)\\
&\ge\sum_e\left[\log2-h\left(\frac{1+w_e}{2}\right)\right]\\
&=\sum_eD\!\left(\operatorname{Bern}\left(\frac{1+w_e}{2}\right)
                \middle\|\operatorname{Bern}\left(\frac12\right)\right)\\
&\ge\frac12\sum_e w_e^2
=\frac12\bigl(N-V(W)\bigr).                              \tag{WR.12}
\end{aligned}
```

Thus a law whose barycenter is globally sign-near necessarily pays
\(\Theta(n^2)\) relative entropy against iid signs.  A low-entropy-cost
conditioning argument therefore cannot produce the weighted target in
(WR.9).  This does not obstruct a deterministic construction or a deliberately
quadratic-cost microcanonical ensemble.

In fact the cost is asymptotically maximal.  Put `delta_e=1-w_e^2` and let
`h` be binary entropy in natural logarithms.  The minority marginal
probability is `(1-|w_e|)/2<=delta_e/2`.  Entropy subadditivity, monotonicity
and concavity of `h`, and Jensen give

```math
\begin{aligned}
H(\mu)
&\le\sum_e h\left(\frac{1-|w_e|}{2}\right)
 \le\sum_e h\left(\frac{\delta_e}{2}\right)\\
&\le N h\left(\frac{V(W)}{2N}\right),\\
D(\mu\|U)
&\ge N\left[\log2-h\left(\frac{V(W)}{2N}\right)\right]. \tag{WR.13}
\end{aligned}
```

Therefore `V(W)=o(N)` implies

```math
D(\mu\|U)\ge N\log2-o(N),
\qquad H(\mu)=o(N).                                      \tag{WR.14}
```

At the `N=Theta(n^2)` edge scale, a sign-near barycenter identifies an
almost complete global sign phase rather than a low-cost canonical tilt.

Thus no rounding theorem can replace global sign-nearness by an unrestricted
fractional input assumption.  Sign-near rounding is now solved, but weighted
recovery itself should not be advertised as a new architecture: it is an
equivalent presentation of exact recovery unless an independent constructor
uses the fractional representation in a provably simpler way.
