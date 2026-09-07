# Wave 52: a K-functional completion tail and its exact failure profile

## Status and conclusion

- **Verified (primary source, hypotheses checked):** Montgomery-Smith's
  one-sided lower tail for a finite Rademacher sum gives a loss `12^{-H}` at
  the coefficient-profile scale `K_{1,2}(beta,sqrt(H))`; see the precise
  statement below.
- **Verified (derived here):** a one-sided Hanson--Wright cutoff for the
  outside quadratic part can be combined with that lower tail by a union
  bound.  It does not assume that the linear and quadratic parts are
  independent.  This gives a genuine far-negative completion theorem.
- **Verified but recurrence-circular:** antipodal pairing and the parent cap
  eliminate the quadratic part exactly, but the resulting sufficient
  condition forces the local oriented energy `e` to be at least the desired
  completion depth `r`.  Saved mass in this branch already proves the local
  principal recurrence.
- **Open target:** prove that a saved set of exact-minimizer local states has
  both a genuinely far negative margin and the K-functional profile required
  by the noncircular theorem.  No such abundance theorem is proved here.

The useful advance is therefore a sharp, sourced, coefficient-level
sufficient condition and an equally sharp falsification profile.  It does
not by itself prove an asymptotic restriction estimate.

## 1. Exact completion setup

Use the notation of ledger (10.1243)--(10.1244).  Thus `T=S^c`, `k=|T|`,
`q=Q(A)`, and for an oriented local state `a=(sigma,y)`,

```math
e=\sigma y^{\mathsf T}A[S]y,
\qquad
Z(w)=L(w)+Q_T(w),
```

where, for uniform `w in {+-1}^T`,

```math
L(w)=\beta^{\mathsf T}w,
\quad \beta=2\sigma A[T,S]y,
\qquad
Q_T(w)=w^{\mathsf T}Bw,
\quad B=\sigma A[T].
```

The diagonal of `B` is zero, so `E Q_T=0`, and

```math
\lVert B\rVert _F=\sqrt{k(k-1)},
\qquad
\lVert B\rVert _{\rm op}\le \lVert A\rVert _{\rm op}
\le\sqrt{2q}.
```

For a negative completion margin put

```math
r=r(S,a):=-\frac{g(S,a)}{p_2}>0,
\qquad
g=(1-p_2)e-h.
```

Then the conditional term in the annealed incidence is exactly
`P_w{Z(w)<=-r}`.

## 2. The sourced linear lower tail

For an integer `H>=1`, define

```math
\mathcal K_H(\beta)
:=K_{1,2}(\beta,\sqrt H)
=\inf_{\beta=u+v}
   \big\{\lVert u\rVert _1+\sqrt H\lVert v\rVert _2\big\}.
```

Montgomery-Smith, *The Distribution of Rademacher Sums*, Proc. Amer. Math.
Soc. **109** (1990), 517--522,
[author PDF](https://stephenmontgomerysmith.github.io/preprints/tail.pdf),
proves for every real `l_2` coefficient sequence and every `t>0` a matching
K-functional tail theorem.  In the proof, when `t^2=H` is an integer, the
constants are explicit (paper pp. 5--6):

```math
\boxed{
\Pr\!\left\{L>\frac12\mathcal K_H(\beta)\right\}
\ge 12^{-H}.}
\tag{R52.1}
```

Finite `beta` is in `l_2`; the signs here are independent Rademachers, exactly
as required by the paper.  Symmetry gives the same bound for
`L<-mathcal K_H/2`.  Thus no central-limit regularity assumption is being
inserted.

If `b_1^*>=...>=b_k^*` is the decreasing rearrangement of `|beta_i|`, put

```math
A_H(\beta)=\sum_{i\le H}b_i^*,
\qquad
D_H(\beta)=\sqrt H
 \left(\sum_{i>H}(b_i^*)^2\right)^{1/2},
\qquad
\Phi_H=A_H+D_H,
```

with the evident truncation if `H>=k`.  The Holmstedt formula quoted and used
in the same paper gives a universal `C_K>=1` such that

```math
C_K^{-1}\Phi_H(\beta)
\le \mathcal K_H(\beta)
\le C_K\Phi_H(\beta).
\tag{R52.2}
```

Thus `A_H` is the concentrated/head mechanism and `D_H` is the diffuse-tail
mechanism.  They are not two new inequalities: they are the two exact scales
seen by the single K-functional theorem.

## 3. Noncircular coupling to the quadratic completion

Theorem 1.1 of Rudelson--Vershynin, *Hanson--Wright inequality and
sub-gaussian concentration*, Electron. Commun. Probab. **18** (2013), no. 82,
[arXiv:1306.2872](https://arxiv.org/abs/1306.2872), applies because the
coordinates of `w` are independent, centered, uniformly subgaussian.  In our
normalization there is a universal `c_HW>0` such that

```math
\Pr\{Q_T>b\}
\le 2\exp\left[-c_{HW}\min\left{
 \frac{b^2}{\lVert B\rVert_F^2},
 \frac b{\lVert B\rVert_{\rm op}}
\right\}\right].
\tag{R52.3}
```

Set

```math
u_H=\frac{H\log 12+\log4}{c_{HW}},
\qquad
b_H=\lVert B\rVert_F\sqrt{u_H}
       +\lVert B\rVert_{\rm op}u_H.
\tag{R52.4}
```

Then (R52.3) gives

```math
\Pr\{Q_T>b_H\}\le\frac12\,12^{-H}.
```

On the intersection

```math
\left\{L<-\frac12\mathcal K_H(\beta)\right\}
\cap\{Q_T\le b_H\},
```

one has `Z<-r` whenever

```math
\boxed{
r+b_H\le\frac12\mathcal K_H(\beta).}
\tag{R52.5}
```

The elementary union bound `P(E cap F)>=P(E)-P(F^c)`, rather than any
independence assertion, now proves the **verified coupled tail**

```math
\boxed{
\Pr_w\{Z(w)\le-r\}\ge\frac12\,12^{-H}.}
\tag{R52.6}
```

This is the desired rigorous lower tail for the cross-linear Rademacher sum
with its actual, dependent, outside quadratic term.

## 4. Target-scale specialization

Fix `0<c<1/4` and put

```math
H=\left\lceil n^{3/4-c}\right\rceil,
\qquad T_n=n^{3/2-c}.
```

For an exact minimizer, `q=Theta(n^{3/2})`, and (R52.4) is uniformly

```math
b_H=O\!\left(n\sqrt H+\sqrt q\,H\right)
=O\!\left(n^{11/8-c/2}+n^{3/2-c}\right)
=O(T_n).
\tag{R52.7}
```

The strict inequality `c<1/4` makes the first term lower order.  Therefore,
if a local state has

```math
r(S,a)/T_n\longrightarrow\infty,
\qquad
\mathcal K_H(\beta)\ge2\{r(S,a)+b_H\},
\tag{R52.8}
```

then it contributes at least `(1/2)12^{-H}` to the completion CDF at a
margin genuinely outside the local `O(T_n)` band.

Via (R52.2), either of the following is a concrete sufficient alternative
(with a sufficiently large universal constant):

```math
\begin{array}{ll}
\text{concentrated:}&
A_H(\beta)\ge 2C_K\{r+b_H\},\\[1mm]
\text{diffuse:}&
D_H(\beta)\ge 2C_K\{r+b_H\}.
\end{array}
\tag{R52.9}
```

The diffuse condition at a far depth entails

```math
\sum_{i>H}(\beta_i^*)^2
\gg \frac{T_n^2}{H}=n^{9/4-c}.
\tag{R52.10}
```

More generally, from the definition
`mathcal K_H(beta)<=sqrt(H)||beta||_2`.  Hence every state certified by
(R52.8), concentrated or diffuse, necessarily obeys

```math
\boxed{
\lVert A[T,S]y\rVert_2^2
=\frac14\lVert\beta\rVert_2^2
\gg n^{9/4-c}.}
\tag{R52.11}
```

Thus the noncircular completion route lives precisely on a high cross-energy
population, above the project-row box scale.  High cross energy alone is not
sufficient: a coefficient profile can concentrate its `l_2` mass on too few
coordinates while having both `A_H` and `D_H` below the required depth.

## 5. Exact minimizer-specific abundance statement

Let `nu_m` denote the product law of a uniform order-`m` selector and a
uniform oriented local projective state, exactly as in (10.1244).  For fixed
`t=O(T_n)`, define `G_n` to be the states satisfying all of

```math
g<0,
\qquad r=-g/p_2,
\qquad r/T_n\to\infty,
\qquad r+b_H\le\mathcal K_H(\beta)/2.
\tag{R52.12}
```

The following is an **open minimizer-specific sufficient package**:

```math
\boxed{\nu_m(G_n)\ge\exp\{-C_0H\}.}
\tag{R52.13}
```

Indeed, averaging (R52.6) gives

```math
Z_t\ge\nu_m(G_n)\frac12\,12^{-H}
     \ge\exp\{-O(H)\}.
\tag{R52.14}
```

Ledger (10.1228)--(10.1230) then gives the bare arbitrary-cut lemma and
convergence.  Unlike the retired moment envelopes, (R52.12) explicitly
places the contributing states outside the local recurrence band.

No argument here proves (R52.13).  It is the exact new burden: simultaneous
saved abundance of far negative local margins and sufficiently rich cross
coefficient profiles.

## 6. The exact cap-pairing bound is circular

There is a second, attractive way to remove `Q_T`.  In each antipodal pair,

```math
\max\{Z(w),Z(-w)\}=Q_T(w)+|L(w)|\le q-e,
```

by the parent cap.  Consequently

```math
\min\{Z(w),Z(-w)\}=Q_T(w)-|L(w)|
\le q-e-2|L(w)|,
```

and hence the **verified exact pairing inequality** is

```math
\Pr\{Z\le-r\}
\ge\frac12\Pr\left\{|L|\ge\frac{q-e+r}{2}\right\}.
\tag{R52.15}
```

Combining (R52.15) with (R52.1) would give `12^{-H}` whenever

```math
\mathcal K_H(\beta)\ge q-e+r.
\tag{R52.16}
```

However this cannot be a new far-margin mechanism.  The two full energies
`e+Z(w)` and `e+Z(-w)` both lie in `[-q,q]`, so

```math
|L(w)|=\frac12|Z(w)-Z(-w)|\le q
\quad\hbox{for every }w.
```

Choosing `w_i=sign(beta_i)` gives

```math
\mathcal K_H(\beta)\le\lVert\beta\rVert_1
=\max_w L(w)\le q.
\tag{R52.17}
```

Thus (R52.16) forces `e>=r`.  If `r/T_n->infinity`, any saved local-state
mass used in the annealed average already has `e/T_n->infinity`; inverse
Hanson--Wright gives the desired local principal recurrence directly.  This
**falsifies the cap-pairing/K-functional combination as an independent
implementation**, although (R52.15) itself remains an exact useful identity.

The Hanson--Wright coupling in Section 3 is different: it pays only
`b_H=O(T_n)`, not the generally order-`q` slack `q-e`, and therefore does
not force `e>=r`.

## 7. Precise failure profile

For this Wave 52 implementation, failure has three distinguishable forms.

1. **No far-margin population.**  The set `g<0` with
   `-g/p_2 >> T_n` has super-saved-small `nu_m` mass.  Then this route cannot
   bypass the local band regardless of coefficient regularity.
2. **Profile failure.**  On every saved far-margin population,

   ```math
   A_H(\beta)+D_H(\beta)<2C_K\{r+b_H\}.
   ```

   Equivalently up to the universal Holmstedt constants, the linear
   K-functional cannot pay both the desired depth and the quadratic cutoff.
   A particularly clear falsifier is `Phi_H=O(T_n)` on all but
   `exp{-omega(H)}` local-state mass.
3. **Quadratic-noise domination.**  `mathcal K_H` may exceed `T_n` but fail
   to exceed `2(r+b_H)` at the actual local threshold.  Without additional
   signed dependence information, subtracting marginal probabilities cannot
   recover this case.  Generic upper concentration at a smaller cutoff is
   insufficient at probability `exp{-Theta(H)}`.

The cap-pairing alternative is already classified separately: its
far-margin success is recurrence-circular by (R52.17).  A future signed-
spectrum theorem could still improve the noncircular quadratic coupling,
but it must use structure beyond antipodal caps and marginal upper tails.

**Final research judgment:** the Hanson--Wright/union theorem
(R52.5)--(R52.6) is the live, noncircular K-functional route, conditional on
the open saved-abundance statement (R52.13).  The cap-pairing theorem
(R52.15) is verified but retired as an independent convergence mechanism,
because its K-functional sufficient condition forces `e>=r`.

## 8. Reproducibility

`tmp/far_tail_kfunctional_r52_check.py` exhaustively checks on stored small
exact minimizers (and one explicitly labeled sampled order-ten minimizer):

- the `L+Q_T` decomposition and antipodal identities;
- `||beta||_1<=q` from the parent cap;
- (R52.15) for every relevant finite threshold; and
- the dependence-safe union inequality underlying (R52.6).

Its saved output is `tmp/far_tail_kfunctional_r52_check.out`.  These finite
checks audit algebra only and are not asymptotic evidence for (R52.13).
