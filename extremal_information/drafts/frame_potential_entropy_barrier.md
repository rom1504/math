# A quadratic entropy barrier for low Rademacher frame potential

**Status.** Task-local theorem draft.  The probability and relative-entropy
theorems below are unconditional.  They use only row exposure, spectral
truncation, and the standard Hanson--Wright inequality.  The final section
separates what the theorem proves about the conference fourth-cumulant
coordinate from the still-open pressure-to-quartic implication.

## 1. Statement

Let `U_r` be the uniform law on sign matrices
`B in {+-1}^{r times r}` and put

```math
\mathcal F(B)=\|BB^T\|_F^2.
\tag{FE.1}
```

The exact mean and deterministic floor are

```math
\mathbb E_{U_r}\mathcal F(B)=2r^3-r^2,
\qquad
\mathcal F(B)\ge r^3.
\tag{FE.2}
```

The second inequality follows from `Tr(BB^T)=r^2` and Cauchy--Schwarz on
the `r` eigenvalues of `BB^T`.

### Theorem FE.1 (fixed frame deficit has speed `r^2`)

For every `delta in (0,1)` there are constants
`c_delta>0` and `r_delta<infinity` such that

```math
\boxed{
 U_r\{\mathcal F(B)\le(2-\delta)r^3\}
 \le \exp\{-c_\delta r^2\}}
\qquad(r\ge r_\delta).
\tag{FE.3}
```

Consequently, if `q_r` is any law on sign matrices and

```math
\mathbb E_{q_r}\mathcal F(B)
\le(2-\delta)r^3+o(r^3),
\tag{FE.4}
```

then there is `c'_delta>0` such that

```math
\boxed{D(q_r\|U_r)\ge c'_\delta r^2}
\tag{FE.5}
```

for all large `r`.

Thus no speed-`O(r)` change of bridge law can lower the frame potential by
a fixed fraction of its leading typical value.

## 2. A uniform one-row lower-tail lemma

Write the rows of `B` as `R_1,...,R_r in {+-1}^r`.  For `k<r`, set

```math
P_k=\left\|[R_1;\ldots;R_k][R_1;\ldots;R_k]^T\right\|_F^2,
\qquad
Z_{k+1}=\sum_{j\le k}\langle R_{k+1},R_j\rangle^2.
\tag{FE.6}
```

Then

```math
\mathcal F(B)=r^3+2\sum_{k=1}^{r-1}Z_{k+1}.
\tag{FE.7}
```

### Lemma FE.2 (late-row charge under a bounded prefix potential)

Fix `a,eta in (0,1]`.  There is `c_(a,eta)>0` such that, for every
`k in [ar,r-1]`, every deterministic prefix satisfying `P_k<=2r^3`, and
an independent uniform sign row `X`,

```math
\Pr\left\{
 \sum_{j\le k}\langle X,R_j\rangle^2
 \le(1-\eta)kr
 \right\}
\le2e^{-c_{a,\eta}r}.
\tag{FE.8}
```

**Proof.**  Put

```math
G=\sum_{j\le k}R_jR_j^T.
\tag{FE.9}
```

This is positive semidefinite and satisfies

```math
\operatorname {Tr}G=kr,
\qquad
\operatorname {Tr}G^2=P_k\le2r^3.
\tag{FE.10}
```

Let

```math
L={8\over a\eta},
\qquad
H=G\,1_{[0,Lr]}(G)
\tag{FE.11}
```

be the spectral truncation of `G`.  The discarded trace is at most

```math
\operatorname {Tr}(G-H)
\le{\operatorname {Tr}G^2\over Lr}
\le{2r^2\over L}
\le{\eta\over4}kr.
\tag{FE.12}
```

Therefore `Tr H >=(1-eta/4)kr`, while

```math
\|H\|_{op}\le Lr,
\qquad
\|H\|_F^2\le\|H\|_{op}\operatorname {Tr}H
\le Lkr^2.
\tag{FE.13}
```

Since `G-H` is positive semidefinite, the event in (FE.8) implies

```math
X^THX-\operatorname {Tr}H\le-{3\eta\over4}kr.
\tag{FE.14}
```

The Hanson--Wright inequality for a Rademacher vector gives a universal
`c_HW>0` for which the probability in (FE.14) is at most

```math
2\exp\left\{-c_{\rm HW}\min\left(
 {9\eta^2k^2r^2/16\over\|H\|_F^2},
 {3\eta kr/4\over\|H\|_{op}}
 \right)\right\}
\le
2\exp\left\{-{9c_{\rm HW}\eta^2a\over16L}r\right\}.
\tag{FE.15}
```

This proves the lemma. `square`

The spectral truncation is essential.  A raw prefix Gram matrix can have a
large top eigenvalue, but (FE.10) shows that this top spectral sector carries
too little trace to account for a fixed lower deviation.  After it is
removed, the ordinary Frobenius/operator Hanson--Wright bound has an
`exp(-Omega(r))` exponent.

## 3. Proof of the probability theorem

Fix `delta in (0,1)` and choose

```math
a=\sqrt{\delta/8},
\qquad
\eta={\delta\over8},
\qquad
\rho={\delta\over16}.
\tag{FE.16}
```

Call an index `k in {ceil(ar),...,r-1}` low if

```math
Z_{k+1}\le(1-\eta)kr.
\tag{FE.17}
```

On the event in (FE.3), every prefix satisfies

```math
P_k\le\mathcal F(B)\le2r^3.
\tag{FE.18}
```

Moreover, that event has at least `rho r` low late indices for all large
`r`.  Indeed, if fewer than `rho r` indices were low, (FE.7) would give

```math
\begin{aligned}
\sum_{k=1}^{r-1}Z_{k+1}
&\ge(1-\eta)\left(
 r\sum_{k=\lceil ar\rceil}^{r-1}k-\rho r^3
 \right)\\
&\ge{r^3\over2}(1-\eta)
 (1-a^2-2\rho-o(1)).
\end{aligned}
\tag{FE.19}
```

But

```math
(1-\eta)(1-a^2-2\rho)
=(1-\delta/8)(1-\delta/4)>1-\delta,
\tag{FE.20}
```

contradicting

```math
\sum_{k=1}^{r-1}Z_{k+1}
\le{1-\delta\over2}r^3
\tag{FE.21}
```

from (FE.3) and (FE.7).

For a specified set `I` of `t=ceil(rho r)` late indices, expose the rows in
order.  At every index in `I`, either the prefix violates (FE.18), in which
case the target event has already failed, or Lemma FE.2 charges conditional
probability at most

```math
p_r=2e^{-c_{a,\eta}r}.
\tag{FE.22}
```

Hence the probability that all indices in `I` are low while their prefixes
obey (FE.18) is at most `p_r^t`.  A union bound over the at most `2^r`
choices of `I` gives

```math
U_r\{\mathcal F(B)\le(2-\delta)r^3\}
\le2^r p_r^{\lceil\rho r\rceil}
\le e^{-c_\delta r^2}
\tag{FE.23}
```

for all large `r`.  This proves (FE.3). `square`

## 4. Relative entropy

It is enough first to assume the right side of (FE.4) is exactly
`(2-delta)r^3`.  Let

```math
E_r=\{\mathcal F(B)\le(2-\delta/2)r^3\},
\qquad p_r=q_r(E_r).
\tag{FE.24}
```

Using the deterministic floor in (FE.2) on `E_r` and the threshold in
(FE.24) on its complement,

```math
(2-\delta)r^3
\ge\mathbb E_{q_r}\mathcal F(B)
\ge\left[2-\frac\delta2-p_r\left(1-\frac\delta2\right)\right]r^3.
\tag{FE.25}
```

Thus

```math
p_r\ge{\delta\over2-\delta}.
\tag{FE.26}
```

Theorem FE.1 with `delta/2` gives
`u_r=U_r(E_r)<=exp(-c_(delta/2)r^2)`.  Data processing under the indicator
of `E_r` yields

```math
\begin{aligned}
D(q_r\|U_r)
&\ge d_{\rm bin}(p_r\|u_r)\\
&\ge p_r\log(1/u_r)-h(p_r)
\ge c'_\delta r^2.
\end{aligned}
\tag{FE.27}
```

The `o(r^3)` term in (FE.4) is absorbed by replacing `delta` with, say,
`3delta/4` for all large `r`.  This proves (FE.5). `square`

## 5. The full conference quartic coordinate also has a quadratic lower tail

The frame theorem combines with the exact conference intertwiner projection
to control every fixed lower deviation of the complete quartic coordinate,
not only its near-minimum edge.

### Theorem FE.3 (fixed quartic-coordinate deficit has speed `r^2`)

Let `A` be any symmetric conference signing of order `r`, let
`epsilon in {+-1}`, and define

```math
J_\epsilon(B)=\|BB^T\|_F^2
 +\|AB+\epsilon BA\|_F^2.
\tag{FE.28}
```

For every `delta in (0,1)` there are `c_delta>0` and `r_delta<infinity`,
uniform in `A` and `epsilon`, such that

```math
\boxed{
 U_r\{J_\epsilon(B)\le(4-\delta)r^3\}
 \le e^{-c_\delta r^2}.}
\tag{FE.29}
```

Consequently, any laws `q_r` satisfying

```math
\mathbb E_{q_r}J_\epsilon(B)
\le(4-\delta)r^3+o(r^3)
\tag{FE.30}
```

obey `D(q_r||U_r)>=c'_delta r^2` for all large `r`.

**Proof.**  If (FE.29)'s event occurs, then at least one of

```math
\|BB^T\|_F^2\le(2-\delta/2)r^3
\tag{FE.31}
```

and

```math
\|AB+\epsilon BA\|_F^2\le(2-\delta/2)r^3
\tag{FE.32}
```

occurs.  The first has probability `exp(-Omega_delta(r^2))` by FE.1.

For the second, conference symmetry gives an orthogonal projection
`P_epsilon` on `R^(r^2)` of rank `d=r^2/2` such that, for `b=vec(B)`,

```math
\|AB+\epsilon BA\|_F^2
=4(r-1)b^TP_\epsilon b.
\tag{FE.33}
```

In particular its exact mean is `4(r-1)d=2r^3-2r^2`.  For all large `r`,
(FE.32) is a lower deviation of `b^TP_epsilon b` from `Tr P_epsilon=d`
by at least `delta r^2/16`.  Hanson--Wright, together with
`||P_epsilon||_F^2=d` and `||P_epsilon||_op=1`, bounds its probability by
`2exp(-c delta^2 r^2)`.  The union bound proves (FE.29).

For the entropy statement, use `J_epsilon(B)>=r^3`.  Under (FE.30), the
event

```math
J_\epsilon(B)\le(4-\delta/2)r^3
\tag{FE.34}
```

has `q_r`-probability bounded below by a positive constant depending only
on `delta` (in the zero-error case, at least `delta/(6-delta)`).  Apply
binary data processing and (FE.29) with `delta/2`, exactly as in
(FE.24)--(FE.27). `square`

This theorem turns the typical values

```math
\mathbb E\|BB^T\|_F^2=2r^3-r^2,
\quad
\mathbb E\|AB+\epsilon BA\|_F^2=2r^3-2r^2,
\quad
\mathbb EJ_\epsilon(B)=4r^3-3r^2
\tag{FE.35}
```

into an exact entropy statement: every fixed leading improvement of the
gauge-invariant fourth-cumulant coordinate costs quadratic bridge entropy.

## 6. Relation to the conference fourth-cumulant target

For a symmetric conference child `A`, orientation `epsilon`, and completed
parent

```math
S_{\epsilon,B}=\begin{pmatrix}A&B\\B^T&\epsilon A\end{pmatrix},
\tag{FE.36}
```

the task-local quartic identity in
`conference_quartic_basin_reduction.md` is

```math
J_\epsilon(B)=\|BB^T\|_F^2
 +\|AB+\epsilon BA\|_F^2,
\tag{FE.37}
```

and

```math
\kappa_4\left({1\over2}z^TS_{\epsilon,B}z\right)
=6J_\epsilon(B)-30r^3+32r^2-10r.
\tag{FE.38}
```

Therefore any rigorous pressure-to-quartic implication of the form

```math
J_\epsilon(B)\le(2-\delta)r^3
\tag{FE.39}
```

does force the frame-potential deficit in FE.1, because the second summand
in (FE.37) is nonnegative.  In fact, the existing small-`beta`,
FMW-power-regular reduction proves the stronger

```math
J_\epsilon(B)\le(1+O(\beta^2))r^3.
\tag{FE.40}
```

For sufficiently small fixed `beta`, FE.1 therefore gives an independent
speed-`r^2` entropy charge on that regular class, even if one discards the
intertwiner channel.

The conclusion is deliberately **not** asserted for an arbitrary
target-reaching bridge.  At one fixed positive temperature, higher
cumulants can compensate for the fourth coefficient, and the current
project has not shown that an unconditioned low-pressure bridge is
FMW-power-regular.  Thus

```math
\boxed{
\text{quartic target plus the proved regularity hypothesis}
\Longrightarrow\text{quadratic frame entropy cost},
}
\tag{FE.41}
```

but the bare pressure target alone is not yet known to force (FE.31).
The new theorem closes the frame-potential entropy question; it does not
close the pressure-to-frame bridge.

## 7. Literature comparison

The only probabilistic input in the proof is Rudelson--Vershynin,
[*Hanson--Wright inequality and sub-Gaussian concentration*](https://arxiv.org/abs/1306.2872),
Theorem 1.1.

As an independent asymptotic check, Groux,
[*Asymptotic Freeness for Rectangular Random Matrices and Large Deviations
for Sample Covariance Matrices With Sub-Gaussian Tails*](https://arxiv.org/abs/1505.05733),
Theorem 1.7 in the arXiv numbering, gives the empirical sample-covariance
LDP at speed `r^(1+alpha/2)` for every fixed `alpha<2`.  Bounded Rademacher
entries lie in the paper's `S_alpha(infinity)` class, and the closed moment
sublevel `int x^2 dmu<=2-delta` excludes the Marchenko--Pastur law.  That
route gives superexponential decay at every speed `r^(2-epsilon)`.  The
elementary row argument above is stronger for this statistic: it reaches
the exact quadratic speed needed for (FE.5).
