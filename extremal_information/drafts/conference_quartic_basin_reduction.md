# A quartic basin reduction for conference bridges

**Status.** Task-local theorem draft.  The exact identities and the
Rademacher small-ball theorem are unconditional.  The pressure-to-quartic
implication is proved on an explicit strict-high-temperature power-regular
class, using the already-audited Fan--Misiakiewicz--Wang--Wen free-energy
theorem.  It does **not** prove the unconditioned bridge lower tail: a low
pressure bridge may still evade the reduction by being power-irregular.

The main structural conclusion is that the fourth pressure coefficient is
not an arbitrary trace statistic.  It is the nonnegative defect from being a
conference completion.  A fixed near-zero defect has probability
`exp(-Theta(r^2))` under a random bridge.

## 1. Setup

Let `A=A_r` be a symmetric conference signing,

```math
A^2=(r-1)I,\qquad A_{ii}=0,\qquad A_{ij}\in\{-1,1\}\quad(i\ne j),
\tag{QR.1}
```

and, for `epsilon in {+-1}` and a sign matrix `B`, put

```math
S_{\epsilon,B}=
\begin{pmatrix}A&B\\B^T&\epsilon A\end{pmatrix},
\qquad
t={\beta\over\sqrt{2r}}.
\tag{QR.2}
```

Write

```math
L_{\epsilon,B}(\beta)
=\log\left[2^{-2r}\sum_z
 \cosh\left({t\over2}z^TS_{\epsilon,B}z\right)\right]
\tag{QR.3}
```

and use the same-temperature conference target

```math
T_r(\beta)=2\log\overline Z_r(A,\beta/\sqrt r).
\tag{QR.4}
```

The gauge-invariant quartic bridge defect is

```math
\boxed{
J_\epsilon(B)=
 \|BB^T\|_F^2+\|AB+\epsilon BA\|_F^2.}
\tag{QR.5}
```

More generally, with two independently gauged children `A,C`, the second
term is `||AB+epsilon BC||_F^2`.  Under

```math
(A,C,B)\mapsto(DAD,ECE,DBE),
\tag{QR.6}
```

for diagonal signs `D,E`, both terms in (QR.5) are unchanged.  Thus this is
genuine bridge information after quotienting child switching labels.

## 2. The exact fourth cumulant is a completion defect

### Theorem QR.1 (exact trace and Boolean-cumulant identities)

For every sign bridge `B` and either orientation,

```math
\boxed{
\operatorname {Tr}S_{\epsilon,B}^4
=6r^3-8r^2+2r+2J_\epsilon(B).}
\tag{QR.7}
```

If `z` is uniform on `{+-1}^{2r}` and

```math
H={1\over2}z^TS_{\epsilon,B}z,
\tag{QR.8}
```

then

```math
\boxed{
\kappa_4(H)
=6J_\epsilon(B)-30r^3+32r^2-10r.}
\tag{QR.9}
```

Consequently the exact fourth-order pressure term is

```math
{t^4\over24}\kappa_4(H)
=\beta^4\left{
 {J_\epsilon(B)\over16r^2}
 -{5r\over16}+{1\over3}-{5\over48r}
 \right}.
\tag{QR.10}
```

The target comparison makes the role of `J-r^3` even more transparent.
The conference child expansion gives

```math
\boxed{
L_{\epsilon,B}(\beta)-T_r(\beta)
= {\beta^2\over4}
 +\beta^4\left\{
 {J_\epsilon(B)-r^3\over16r^2}
 -{1\over3}+{5\over16r}
 \right\}
 +O_{r,B}(\beta^6).}
\tag{QR.10a}
```

Thus the only bridge-dependent fourth-order contribution is the
nonnegative completion defect `(J-r^3)/(16r^2)`.  The remainder in
(QR.10a) is deliberately labelled nonuniform; Section 4 supplies the
uniform replacement on a stated asymptotic class.

Moreover,

```math
\boxed{J_\epsilon(B)\ge r^3.}
\tag{QR.11}
```

Equality holds exactly when

```math
BB^T=B^TB=rI,\qquad AB+\epsilon BA=0,
\tag{QR.12}
```

or, equivalently,

```math
S_{\epsilon,B}^2=(2r-1)I.
\tag{QR.13}
```

Thus `J-r^3` is precisely the quartic defect from an order-`2r`
conference completion.

**Proof.**  Block multiplication gives

```math
S^2=
\begin{pmatrix}
 (r-1)I+BB^T&AB+\epsilon BA\\
 (AB+\epsilon BA)^T&(r-1)I+B^TB
\end{pmatrix}.
\tag{QR.14}
```

Taking its squared Frobenius norm proves (QR.7).  For any complete hollow
signing `S` of order `N`, the even-multigraph expansion gives

```math
\kappa_4\left(\sum_{i<j}s_{ij}z_iz_j\right)
=3\operatorname {Tr}S^4-2N(N-1)(3N-4).
\tag{QR.15}
```

Substitution of `N=2r` and (QR.7) proves (QR.9)--(QR.10).
Finally, `Tr(BB^T)=r^2`, so Cauchy--Schwarz on its `r` eigenvalues gives

```math
\|BB^T\|_F^2\ge {\operatorname {Tr}(BB^T)^2\over r}=r^3.
\tag{QR.16}
```

Equality in (QR.11) forces equality in (QR.16) and zero intertwiner,
which is (QR.12).  Since a square `B` with `BB^T=rI` also has `B^TB=rI`,
(QR.14) proves the equivalence with (QR.13). `square`

The known universal double lies at the correct asymptotic edge.  For
`epsilon=-1` and `B=A+I`,

```math
J_{-}(A+I)=r^3+4r(r-1)=r^3+O(r^2).
\tag{QR.17}
```

So the near-minimum event below is nonempty even when an exact conference
completion does not exist at order `2r`.

## 3. Near-minimal quartic defect already has speed `r^2`

The next theorem needs no pressure theorem.

There is no hidden tradeoff between the two summands at this edge.  Since
the Gram term is already at least `r^3`,

```math
J_\epsilon(B)\le(1+\delta)r^3
\quad\Longrightarrow\quad
r^3\le\|BB^T\|_F^2\le(1+\delta)r^3,
\qquad
\|AB+\epsilon BA\|_F^2\le\delta r^3.
\tag{QR.17a}
```

Thus a near-minimal sum forces both approximate row orthogonality and an
approximate intertwining relation; neither channel can pay for the other.
For comparison, a uniform bridge has the exact means

```math
\mathbb E\|BB^T\|_F^2=2r^3-r^2,
\qquad
\mathbb E\|AB+\epsilon BA\|_F^2=2r^3-2r^2,
\tag{QR.17b}
```

so `E J=4r^3-3r^2`.  The target edge is a constant lower deviation in
both gauge-invariant channels, not a fluctuation on their natural standard
deviation scale.

### Theorem QR.2 (quartic completion small ball)

There is a universal `c>0` such that, for every conference `A`, either
orientation, every `0<delta<=1/2`, and a uniform Rademacher bridge `B`,

```math
\boxed{
\Pr\{J_\epsilon(B)\le(1+\delta)r^3\}
\le2e^{-cr^2}.}
\tag{QR.18}
```

**Proof.**  Conference symmetry and `Tr A=0` show that the two eigenvalues
`+-sqrt(r-1)` each have multiplicity `r/2`.  In that eigenbasis the map

```math
B\mapsto AB+\epsilon BA
\tag{QR.19}
```

vanishes on half of matrix space and has singular value `2sqrt(r-1)` on
the other half.  Hence there is an orthogonal projection `P_epsilon` on
`R^(r^2)`, of rank

```math
d={r^2\over2},
\tag{QR.20}
```

such that, for `b=vec(B)`,

```math
\|AB+\epsilon BA\|_F^2
=4(r-1)b^TP_\epsilon b.
\tag{QR.21}
```

By (QR.16), the event in (QR.18) implies

```math
b^TP_\epsilon b
\le {\delta r^3\over4(r-1)}\le {d\over2}.
\tag{QR.22}
```

Here `E b^TP_epsilon b=Tr P_epsilon=d`,
`||P_epsilon||_F^2=d`, and `||P_epsilon||_op=1`.
The Hanson--Wright inequality therefore gives

```math
\Pr\{b^TP_\epsilon b\le d/2\}
\le2\exp\{-c\min(d/4,d/2)\}
\le2e^{-c'r^2}.
\tag{QR.23}
```

This proves the claim. `square`

The only imported input in this proof is the standard Hanson--Wright
inequality in its Frobenius/operator-norm form; see Rudelson--Vershynin,
[*Hanson--Wright inequality and sub-Gaussian concentration*](https://arxiv.org/abs/1306.2872),
Theorem 1.1.

## 4. A uniform spectral remainder at small beta

For a compactly supported probability law `nu`, let `R_nu` be its
`R`-transform and put

```math
\mathfrak f(\nu)={1\over2}\int_0^1R_\nu(u)\,du.
\tag{QR.24}
```

This is the strict-high-temperature limiting normalized pressure in the
Fan--Misiakiewicz--Wang--Wen theorem already audited in
`conference_reverse_kl_fixed_temperature_obstruction.md`.

### Lemma QR.3 (uniform even fourth-order remainder)

For every `K>=1` there are explicit finite constants `C_K,beta_K>0` such
that every centered probability law `nu` supported on `[-K,K]` and every
`0<beta<beta_K` satisfy

```math
\boxed{
{\mathfrak f(\beta\nu)+\mathfrak f(-\beta\nu)\over2}
={\beta^2m_2(\nu)\over4}
 +{\beta^4\{m_4(\nu)-2m_2(\nu)^2\}\over8}
 +R_6,
\qquad |R_6|\le C_K\beta^6.}
\tag{QR.25}
```

One may take `beta_K=1/(32K)` and, with the crude cumulant count used
below, `C_K=(16K)^6/6`.

**Proof.**  Moment--free-cumulant inversion on noncrossing partitions,
`|m_j(nu)|<=K^j`, `|NC(j)|<=4^j`, and the Catalan bound on the
noncrossing Möbius function give the deliberately crude uniform estimate

```math
|r_j(\nu)|\le(16K)^j.
\tag{QR.26}
```

For `16K beta<=1/2`, the `R`-series is absolutely and uniformly
convergent on `[0,1]`, so

```math
\mathfrak f(\beta\nu)
=\sum_{j\ge1}{\beta^jr_j(\nu)\over2j}.
\tag{QR.27}
```

Averaging `beta` and `-beta` deletes odd terms.  Since `r_2=m_2` and
`r_4=m_4-2m_2^2`, the first two terms are those in (QR.25), while the
geometric tail is at most `(16K)^6 beta^6/6`. `square`

This is a genuinely uniform remainder theorem for the spectral limiting
functional.  It is not, by itself, a uniform finite-matrix expansion of
Boolean pressure; that transfer is exactly where power regularity enters.

## 5. Low regular pressure forces a near-conference quartic statistic

Call a sequence `Y_r` **FMW-regular with bound `K`** when

```math
\|Y_r\|_{op}\le K
\tag{QR.28}
```

and, for every fixed positive integer `p` and every `eta>0`,

```math
\max_i\left|(Y_r^p)_{ii}-{1\over2r}\operatorname {Tr}Y_r^p\right|
 +\max_{i\ne j}|(Y_r^p)_{ij}|
 <(2r)^{-1/2+\eta}
\tag{QR.29}
```

eventually.  This is precisely the maximum-entry power hypothesis of the
imported theorem, separated from its strict operator condition.

### Theorem QR.4 (small-beta pressure-to-completion reduction)

Fix `K>=1`.  There is `beta_0(K)>0` such that the following holds.  Let
`A_r` run through a Paley conference sequence and let `epsilon_r,B_r` be
any bridge sequence for which

```math
Y_r={S_{\epsilon_r,B_r}\over\sqrt{2r}}
\tag{QR.30}
```

is FMW-regular with bound `K`.  If `0<beta<beta_0(K)` and

```math
L_{\epsilon_r,B_r}(\beta)
\le T_r(\beta)+o(r),
\tag{QR.31}
```

then

```math
\boxed{
\limsup_{r\to\infty}{J_{\epsilon_r}(B_r)\over r^3}
\le1+32(C_K+C_1)\beta^2.}
\tag{QR.32}
```

Here `C_K,C_1` are the explicit constants from Lemma QR.3.  In particular,
`beta_0(K)` may be chosen so that the right side is at most `5/4`.

**Proof.**  Pass to a subsequence on which the empirical spectral laws
`nu_r` of `Y_r` converge to `nu` and `J/r^3` converges to `j` at its
limsup.  The laws are centered and supported on `[-K,K]`.  Directly from
(QR.7),

```math
\begin{aligned}
m_2(\nu_r)&=1-{1\over2r},\\
m_4(\nu_r)&={3\over4}+{J_\epsilon(B)\over4r^3}
             -{1\over r}+{1\over4r^2}.
\end{aligned}
\tag{QR.33}
```

Write

```math
F_r(\pm\beta)={1\over2r}\log Z_{2r}(\pm\beta Y_r),
\tag{QR.34a}
```

where `Z` is normalized by `2^(-2r)`.  Arithmetic--geometric mean gives
the exact inequality

```math
{L_{\epsilon,B}(\beta)\over2r}
\ge {F_r(\beta)+F_r(-\beta)\over2}.
\tag{QR.34}
```

Choose `beta_0` also so that `beta K<1/2`.  The audited FMW theorem and
(QR.29) identify the subsequential limits on the right with
`mathfrak f(+-beta nu)`.  Lemma QR.3 and (QR.33) therefore give

```math
\liminf {L_{\epsilon,B}(\beta)\over2r}
\ge {\beta^2\over4}
 +{\beta^4(j-5)\over32}-C_K\beta^6.
\tag{QR.35}
```

The same theorem for Paley children, or the exact Bernoulli spectral law,
and Lemma QR.3 give

```math
{T_r(\beta)\over2r}
=\psi(\beta)+o(1)
\le {\beta^2\over4}-{\beta^4\over8}+C_1\beta^6+o(1).
\tag{QR.36}
```

Combining (QR.31), (QR.35), and (QR.36) yields

```math
0\ge {\beta^4(j-1)\over32}-(C_K+C_1)\beta^6,
\tag{QR.37}
```

which is (QR.32). `square`

There is a finite, checkable version of the regular class.  Let
`q_r -> infinity`, `eta_r -> 0`, and require (QR.28) together with
(QR.29), using exponent `eta_r`, for all `p<=q_r`; for example one can
take `q_r=floor(log log r)` and `eta_r=1/sqrt(log log r)` after harmless
small-order conventions.  Denote this class by `D_r(K)`.

### Corollary QR.5 (speed-`r^2` on the power-regular class)

Choose `beta<beta_0(K)` so that the bound in (QR.32) is at most `5/4`, and
let `xi_r -> 0`.  For all sufficiently large `r`,

```math
\left\{
Y_r\in\mathcal D_r(K),\quad
L_{\epsilon,B}(\beta)\le T_r(\beta)+\xi_r r
\right\}
\subseteq
\left\{J_\epsilon(B)\le {3\over2}r^3\right\}.
\tag{QR.38}
```

Consequently, under the uniform bridge law,

```math
\boxed{
\Pr\left\{
Y_r\in\mathcal D_r(K),\quad
L_{\epsilon,B}(\beta)\le T_r(\beta)+\xi_r r
\right\}
\le2e^{-cr^2}.}
\tag{QR.39}
```

**Proof.**  Failure of the eventual containment would select a sequence
in `D_r(K)` contradicting Theorem QR.4.  Apply Theorem QR.2 with
`delta=1/2`. `square`

## 6. What the quartic statistic does and does not prove

The implication is now rigorous on a nontrivial class:

```math
\text{same-temperature pressure hit + power regularity}
\Longrightarrow
\text{near-conference quartic completion}
\Longrightarrow
\text{speed-}r^2\text{ rarity}.
\tag{QR.40}
```

This is stronger than a formal Taylor coefficient and exposes a concrete
gauge-invariant basin coordinate.  It also identifies the exact remaining
gap: prove that a target-reaching bridge cannot be power-irregular, or prove
the pressure-to-quartic lower bound directly under a weaker hypothesis such
as strict operator regularity.  Neither follows from this note.

The fourth statistic cannot control arbitrary fixed temperatures without
some small-`beta` or higher-cumulant hypothesis.  For example, the two
symmetric laws

```math
\nu_1={1\over2}\delta_0
      +{1\over4}(\delta_{\sqrt2}+\delta_{-\sqrt2}),
\tag{QR.41}
```

and

```math
\nu_2={2\over5}(\delta_{1/\sqrt2}+\delta_{-1/\sqrt2})
      +{1\over10}(\delta_{\sqrt3}+\delta_{-\sqrt3})
\tag{QR.42}
```

both have `m_2=1,m_4=2`, but have `m_6=4` and `m_6=11/2`, respectively.
Their even spectral free energies therefore first separate at order
`beta^6` (the sixth free cumulants differ by `3/2`).  This is a spectral
model warning, not a sign-bridge counterexample.  It explains why (QR.25)
must keep a quantitative remainder and why `J` alone is not an all-
temperature state.

## 7. Reproducible exact audit

The task-local program
`extremal_information/experiments/conference_quartic_identity_audit.py`
checks (QR.7), (QR.9), and (QR.21) by exact Boolean enumeration.  It
exhausts the `16` order-two bridges and uses `64` seeded order-six bridges:

```bash
.venv/bin/python \
  extremal_information/experiments/conference_quartic_identity_audit.py \
  --samples-6 64 \
  --output \
  extremal_information/experiments/results/conference_quartic_identity_audit.json
```

The frozen SHA-256 hashes are

```text
script  116f4d6b175160ecfe8dffdf5180c1bbc5353e7d6c79148f4706b5f4ad2f8b23
result  6f5cac6d45b9ccac76026da66f707be1a97773177eca2429e15d70f7bf359723
```

The computation is an identity regression only; no asymptotic claim is
inferred from its sample.
