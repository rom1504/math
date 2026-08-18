# Operator localization does not identify the conference pressure basin

**Status.** Task-local theorem/counterexample note.  This note does not edit
the canonical theory files.  Its purpose is to stratify one concrete part of
the operator-irregular conference sector and identify exactly what that
stratification can and cannot prove about the remaining lower-pressure
large-deviation problem.

## 1. Setup

Let `A=A_r` be a symmetric conference signing, fix an orientation
`epsilon`, and write

```math
S_{\epsilon,B}=
\begin{pmatrix}A&B\\B^T&\epsilon A\end{pmatrix},
\qquad
t={\beta\over\sqrt{2r}},
\tag{OL.1}
```

and

```math
f_{\epsilon,r}(B)=
\log\mathbb E_{x,y}
\cosh\left({t\over2}(x,y)^TS_{\epsilon,B}(x,y)\right).
\tag{OL.2}
```

The expectation is normalized.  If one bridge sign is flipped, the
Hamiltonian in (OL.2) changes by `2` for every spin pair.  Hence

```math
|f(B)-f(B')|\le 2t
\tag{OL.3}
```

for a one-bit flip.  Under the uniform bridge law `U_r`, the archived
conference theorem gives

```math
\mathbb E_{U_r}f_{\epsilon,r}(B)=h_\beta r+o(r)
\tag{OL.4}
```

in the strict high-temperature range used by the current conference
campaign.

## 2. A localized-information transport lemma

The following statement is independent of the conference identities.  It
is the useful coupling between a localized entropy stratum and Boolean
pressure.

### Lemma OL.1 (conditional-information diffuseness)

Order the `d=r^2` bridge coordinates by a permutation `pi`.  Under an
arbitrary law `q`, define the chain-rule entropy increments

```math
d_j^\pi=
\mathbb E_q D\left(
 q(B_{\pi(j)}\mid B_{\pi(1)},\ldots,B_{\pi(j-1)})
 \middle\| {1\over2}(\delta_{-1}+\delta_{+1})
 \right).
\tag{OL.5}
```

They obey `sum_j d_j^pi=D(q||U_r)`.  If this entropy is nonzero, put

```math
s_*(q)=\inf_\pi
 {\left(\sum_{j=1}^{r^2}\sqrt{d_j^\pi}\right)^2
  \over D(q\|U_r)};
\tag{OL.5a}
```

set `s_*(U_r)=0`.  This effective conditional-information support lies
between `1` and `r^2` whenever it is nonzero.  Then

```math
\boxed{
 |\mathbb E_qf-\mathbb E_{U_r}f|
 \le \beta\sqrt{{D(q\|U_r)s_*(q)\over r}}.}
\tag{OL.6}
```

Consequently, if `D(q||U_r)<=Cr` and

```math
\mathbb E_q f\le(h_\beta-\eta)r,
\tag{OL.7}
```

then

```math
\boxed{
 s_*(q)\ge
 \left({\eta^2\over\beta^2C}-o(1)\right)r^2.}
\tag{OL.8}
```

Thus a speed-`r` tilted law which gains a fixed linear amount of pressure
must distribute its conditional information over a positive fraction of
all bridge coordinates.

**Proof.**  For a fixed ordering, let `q^(j)` have the `q`-joint marginal
on the first `j` ordered coordinates and independent uniform signs on the
remaining coordinates.  Thus `q^(0)=U_r` and `q^(r^2)=q`.  Conditional on
a prefix of length `j-1`, average `f` over all coordinates after `j` and
call the resulting two-point function `g_j`.  Its oscillation in coordinate
`j` is at most `2t`.  Pinsker's inequality therefore gives

```math
|\mathbb E_{q^{(j)}}f-\mathbb E_{q^{(j-1)}}f|
\le 2t\,\mathbb E_q
 \sqrt{{D(q(B_{\pi(j)}\mid B_{\pi(<j)})\|U_1)\over2}}
\le\sqrt2t\sqrt{d_j^\pi}.
\tag{OL.9}
```

Sum the hybrid increments, use `sqrt(2)t=beta/sqrt(r)`, and then minimize
over `pi`:

```math
|\mathbb E_qf-\mathbb E_{U_r}f|
\le {\beta\over\sqrt r}\inf_\pi\sum_j\sqrt{d_j^\pi}
=\beta\sqrt{{D(q\|U_r)s_*(q)\over r}},
\tag{OL.10}
```

which is (OL.6).  Equations (OL.4), (OL.6), and (OL.7) imply (OL.8).
`square`

The result applies directly to the canonical negative-disorder tilt, rather
than only to hand-designed laws.  For fixed `lambda>0`, let

```math
{dq_{\lambda,r}\over dU_r}(B)
={e^{-\lambda f(B)}\over\mathbb E_{U_r}e^{-\lambda f}}.
\tag{OL.10a}
```

Because `f>=0` and Jensen's inequality gives
`-log E_U e^(-lambda f)<=lambda E_U f`, one has

```math
D(q_{\lambda,r}\|U_r)
=-\lambda\mathbb E_{q_{\lambda,r}}f
 -\log\mathbb E_{U_r}e^{-\lambda f}
\le\lambda h_\beta r+o(r).
\tag{OL.10b}
```

Consequently, if this actual tilted law satisfies
`E_(q_(lambda,r))f<=(h_beta-eta)r`, then

```math
\boxed{
s_*(q_{\lambda,r})
\ge\left({\eta^2\over\beta^2\lambda h_\beta}-o(1)\right)r^2.}
\tag{OL.10c}
```

Thus any genuine fixed-tilt phase is necessarily a diffuse conditional-
information phase.  Indeed, with
`phi(lambda)=-log E_U exp(-lambda f)`, concavity and `phi(0)=0` give

```math
\mathbb E_{q_{\lambda,r}}f
=\phi'(\lambda)
\le {\phi(\lambda)\over\lambda}
=\mathcal R_{\lambda,r}.
\tag{OL.10c'}
```

Hence a fixed linear drop in the negative-moment pressure itself already
triggers (OL.10c).

If, more specially, `q=q_C\otimes U_{C^c}` for a fixed set `C` of `m`
coordinates, order `C` first.  All later entropy increments vanish and
Cauchy--Schwarz gives the useful support form

```math
\boxed{
 |\mathbb E_qf-\mathbb E_{U_r}f|
 \le \beta\sqrt{{mD(q\|U_r)\over r}}.}
\tag{OL.10d}
```

Thus a factorized-support tilt with entropy at most `Cr` and a fixed linear
pressure gain must have
`m>=(eta^2/(beta^2C)-o(1))r^2`.  If `C` is contained in `k` complete rows,
this requires `k=Omega(r)`.

### Corollary OL.2 (fixed-support localized spectral strata are annealed in mean)

Let `I` be a fixed set of `k=o(r)` rows and let `E_I` be any event determined
only by the `kr` bridge signs in those rows.  If

```math
U_r(E_I)\ge e^{-Cr},
\tag{OL.11}
```

then

```math
\boxed{
 \mathbb E[f\mid E_I]=h_\beta r+o(r).}
\tag{OL.12}
```

More quantitatively, its difference from the unconditioned mean is at most
`beta sqrt(Ckr)`.  The same statement holds with rows replaced by columns.

**Proof.**  Conditioning on `E_I` leaves the complementary signs independent
and uniform, and its relative entropy is `log(1/U_r(E_I))`.  Apply OL.1 with
`m=kr`. `square`

This is an information-support statement, not merely a spectral estimate.
It says that a speed-`r` event supported on a sublinear set of rows may cause
operator irregularity, but that irregularity cannot itself account for a
fixed linear pressure gain.

## 3. An exact speed-`r` localized operator stratum

The preceding regime is nonempty and has precisely the entropy scale left
open by the conference phase-boundary calculation.

Fix `k>=2` and a specified row set `I` of size `k`.  Let

```math
E_I^{\rm twin}
=\{B:R_i=R_j\text{ for every }i,j\in I\}.
\tag{OL.13}
```

Then

```math
U_r(E_I^{\rm twin})=2^{-(k-1)r}.
\tag{OL.14}
```

On this event the `k by r` row submatrix has singular value exactly
`sqrt(kr)`, and hence

```math
\|B\|_{op}\ge\sqrt{kr},
\qquad
\|tS_{\epsilon,B}\|_{op}
\ge\beta\sqrt{k/2}.
\tag{OL.15}
```

Consequently, for any fixed regularity cutoff `kappa`, choosing
`k>2kappa^2/beta^2` puts this entire speed-`r` stratum outside the
operator-regular set.  Nevertheless Corollary OL.2 gives

```math
\mathbb E[f\mid E_I^{\rm twin}]
=h_\beta r+O_{\beta,k}(\sqrt r)+o(r).
\tag{OL.16}
```

Thus the most elementary localized singular spike is operator-irregular at
speed `r` but is not a favorable pressure phase in conditional mean.

There is a sharp limitation.  Under the conditioned law in (OL.13), `f` is
still a bounded-difference function.  Treat the common row as one block of
`r` independent signs and all other rows as ordinary signs.  The sum of
squared one-coordinate oscillations is `O_(beta,k)(r)`, so the elementary
lower-tail bound at a deviation `eta r` is only

```math
\Pr\{f\le\mathbb E[f\mid E_I^{\rm twin}]-\eta r
       \mid E_I^{\rm twin}\}
\le e^{-c_{\beta,k,\eta}r}.
\tag{OL.17}
```

This has exactly, rather than better than, the disorder speed under
investigation.  Localization plus bounded differences cannot exclude a
speed-`r` favorable subset inside the stratum.

## 4. Pointwise pressure cannot be controlled by singular localization

There is also a deterministic obstruction to any theorem saying that a
localized high singular mode forces typical or high pressure.

### Proposition OL.3 (localized spikes can be planted into a low-pressure bridge)

Use orientation `epsilon=-1` and the universal-double bridge

```math
B_r^0=A_r+I,
\qquad
{f_{-,r}(B_r^0)\over r}\longrightarrow\tau_\beta<h_\beta.
\tag{OL.18}
```

Fix an integer `k` and any sign row `w in {+-1}^r`.  Overwrite `k`
specified rows of `B_r^0` by `w`, producing `B_r^w`.  Then

```math
\|B_r^w\|_{op}\ge\sqrt{kr},
\qquad
d_H(B_r^w,B_r^0)\le kr,
\tag{OL.19}
```

and

```math
\boxed{
 f_{-,r}(B_r^w)=\tau_\beta r+O_{\beta,k}(\sqrt r)+o(r).}
\tag{OL.20}
```

In particular, choosing `k>2kappa^2/beta^2` gives operator-irregular,
localized-spike bridges which retain the full linear conference pressure
advantage.  There are `2^r` distinct choices of `w`.

**Proof.**  The repeated-row singular certificate proves the first part of
(OL.19).  At most `kr` bits are changed, and (OL.3) gives

```math
|f(B_r^w)-f(B_r^0)|\le 2tkr=O_{\beta,k}(\sqrt r).
\tag{OL.21}
```

Combine this with (OL.18). `square`

The family in OL.3 has only `exp(Theta(r))` members inside a cube of
`2^(r^2)` bridges, so it is quadratically rare.  It is a pointwise
falsifier, not the missing speed-`r` basin.

## 5. Singular-vector stratification and the surviving target

For completeness, standard rectangular subgaussian norm concentration gives
universal constants `C_0,c>0` such that, for fixed `I`,

```math
\Pr\{\|B_{I,*}\|_{op}\ge C_0(\sqrt r+\sqrt k)+s\}
\le2e^{-cs^2}.
\tag{OL.22}
```

At `s=Theta(sqrt r)` this is a speed-`r` event.  The twin-row construction
gives a matching speed-`r` lower example for a constant-strength localized
spike.  Hence spectral localization genuinely identifies a large part of
the exceptional operator tail; it is not merely an empty classification.

The results above give the following exact trichotomy.

1. A fixed-support localized operator event of speed `r` has typical
   pressure **in conditional mean** (OL.2).
2. Localized operator irregularity is compatible **pointwise** with the
   full low conference pressure (OL.3).
3. Elementary concentration inside such a stratum stops at speed `r`
   (OL.17).

Therefore singular-vector localization alone cannot close the conference
lower-pressure LDP.  What it does prove is a stricter information-support
condition on any proposed entropy-`O(r)` tilted bridge law:

```math
\boxed{
\text{a linear pressure gain at entropy }O(r)
\text{ requires }s_*(q)=\Theta(r^2).}
\tag{OL.23}
```

In particular, in the factorized-support case it must reweight
`Theta(r^2)` bridge coordinates.

The natural unresolved stratum is consequently not “localized versus
delocalized singular vectors” by itself.  It is a **diffuse weak tilt**:
`O(r)` total information distributed over a positive fraction of the
`r^2` bridge bits, possibly producing a delocalized rank-one or higher-rank
spike.  For example, the product law with means
`E B_ij=alpha u_i v_j/sqrt(r)` has relative entropy
`(alpha^2/2+o(1))r` and a rank-one mean of singular value
`alpha sqrt(r)` while touching all coordinates.  Whether such a diffuse
law, or a dependent analogue, can lower the conference pressure by
`Theta(r)` is not decided here.

That diffuse entropy/pressure problem is strictly narrower than the former
“deeply operator-irregular” target.  Another net stratification by singular
localization cannot decide it; a joint variational inequality for diffuse
bridge tilts is required.
