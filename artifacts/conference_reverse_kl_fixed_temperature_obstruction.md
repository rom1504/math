# Fixed-temperature conference obstruction to joint reverse-KL compensation

Status: **proved scalable obstruction for the conference implementation**.  The only
random-matrix input beyond published theorems is Lemma 3.1 below, a fixed-degree walk/chaos
lemma whose proof is included because the maximum-entry power condition needed by the July
2026 theorem is stronger than ordinary empirical spectral convergence.  The July 2026
deterministic high-temperature free-energy theorem of Fan--Misiakiewicz--Wang--Wen turns the
previous formal `beta^4/beta^8` tangent mismatch into a theorem at every sufficiently small fixed
positive temperature.  Along an infinite Paley conference sequence, the *full* joint reverse KL
is `o(r)`, while the uncompensated same-temperature interface is `gamma(beta) r+o(r)` with
`gamma(beta)>0`.

This falsifies the uniform-output reverse-KL mechanism on conference children.  It does **not**
falsify the minimizer-optimized criterion, because conference signings are not known to minimize
the finite-temperature pressure.  It also does not rule out selecting an exponentially rare bridge.

## 1. Setup and statement

For an order-`r` symmetric conference signing `A_r`,

```math
A_r^2=(r-1)I,
\qquad (A_r)_{ii}=0,
\qquad (A_r)_{ij}\in\{-1,1\}.
\tag{1.1}
```

There is an infinite such sequence: take the Paley conference matrices of order `r=q+1` along
primes `q=1 mod 4`.  Put

```math
s={\beta\over\sqrt r},
\qquad t={\beta\over\sqrt{2r}},
\qquad
S_{\epsilon,B}=
\begin{pmatrix}A_r&B\\B^{\mathsf T}&\epsilon A_r\end{pmatrix},
\tag{1.2}
```

where `epsilon` is uniform in `{-1,1}` and the entries of the `r` by `r` matrix `B` are iid
Rademacher signs.  As in the reverse-KL interface note, let

```math
\overline Z_k(C,u)=2^{-k}\sum_x\cosh\!\left({u\over2}x^{\mathsf T}Cx\right),
\tag{1.3}
```

and define

```math
\begin{aligned}
\mathcal M_r(\beta)
&=2\log\overline Z_r(A_r,s)
  -\mathbb E_{\epsilon,B}\log\overline Z_{2r}(S_{\epsilon,B},t),\\
D_r(\beta)
&=D_{\rm KL}(U\Vert\Pi_{A_r,A_r,t}),\\
G_r(\beta)
&=r^2\log\cosh t
 -2\{\log\overline Z_r(A_r,s)-\log\overline Z_r(A_r,t)\}
 -D_r(\beta).
\end{aligned}
\tag{1.4}
```

The exact joint output identity gives `G_r=-mathcal M_r`.

Define

```math
\psi(c)={1\over4}\left[
 \sqrt{1+4c^2}-1
 -\log\!\left({1+\sqrt{1+4c^2}\over2}\right)
 \right]
\tag{1.5}
```

and

```math
\gamma(\beta)
={\beta^2\over4}-2\psi(\beta)+2\psi(\beta/\sqrt2).
\tag{1.6}
```

**Theorem 1.1 (fixed-temperature linear obstruction).**  For every fixed

```math
0<\beta<{\sqrt2\over6},
\tag{1.7}
```

along the Paley conference sequence,

```math
\boxed{
{D_r(\beta)\over r}\longrightarrow0,
\qquad
{G_r(\beta)\over r}\longrightarrow\gamma(\beta)>0,
\qquad
{\mathcal M_r(\beta)\over r}\longrightarrow-\gamma(\beta)<0.}
\tag{1.8}
```

Moreover, for each orientation separately, there is `c_beta>0` such that for all sufficiently
large `r`,

```math
\Pr_B\!\left\{
 \log\overline Z_{2r}(S_{\epsilon,B},t)
 \le 2\log\overline Z_r(A_r,s)
 \right\}
\le \exp(-c_\beta r).
\tag{1.9}
```

Thus neither the mean output nor any fixed lower quantile supplies the needed cancellation.

## 2. Exact audit of the imported free-energy theorem

Corollaries 2.8 and 2.10(c) of Fan, Misiakiewicz, Wang, and Wen,
[“Dynamical mean-field limit and replica-symmetric free energy for the
orthogonally-invariant SK model”](https://arxiv.org/abs/2607.10102), state the following special
case.  Let `X_k` be a deterministic symmetric matrix sequence with `||X_k||_op<1/2`, a compact
limiting spectral law `Lambda`, and, for every fixed positive integer `j` and every `eta>0`,

```math
\max_i\left|(X_k^j)_{ii}-{1\over k}\operatorname{Tr}X_k^j\right|
 +\max_{i\ne l}|(X_k^j)_{il}|
 <k^{-1/2+\eta}
\tag{2.1}
```

eventually.  At zero external field,

```math
{1\over k}\log\sum_x
 \exp\!\left({1\over2}x^{\mathsf T}X_kx\right)
\longrightarrow
\log2+{1\over2}\int_0^1R_\Lambda(u)\,du.
\tag{2.2}
```

There is no hidden overlap parameter here.  Global spin-flip symmetry makes every one-spin mean
zero at every finite order, so the paper's two-replica parameter is exactly `q_*=0`.  Equation
(2.2) is then its replica-symmetric formula with `q_*=0`.

The normalization matches ours exactly:

```math
X=uC
\quad\Longrightarrow\quad
{1\over2}x^{\mathsf T}Xx
=u\sum_{i<j}c_{ij}x_ix_j.
\tag{2.3}
```

The paper treats the ordinary exponential partition function, whereas (1.3) uses `cosh`.  But

```math
\overline Z_k(C,u)
=2^{-k-1}\{Z_{uC}+Z_{-uC}\}.
\tag{2.4}
```

All limiting laws below are symmetric.  Therefore (2.2) gives the same limit for `uC` and
`-uC`, and (2.4) yields

```math
{1\over k}\log\overline Z_k(C,u)
\longrightarrow {1\over2}\int_0^1R_\Lambda(v)\,dv.
\tag{2.5}
```

The condition `beta<sqrt(2)/6` is a safe, explicit high-temperature range, not an optimized one.
For the child matrices the limiting operator norms are respectively `beta` and `beta/sqrt(2)`.
For the parent,

```math
\left\|{\beta\over\sqrt{2r}}S_{\epsilon,B}\right\|_{\rm op}
\le {\beta\over\sqrt2}\left(
 \sqrt{1-1/r}+{\|B\|_{\rm op}\over\sqrt r}
 \right)
\longrightarrow {3\beta\over\sqrt2}
\tag{2.6}
```

almost surely, by the Bai--Yin norm limit for an iid Rademacher matrix.  Hence all three sequences
are eventually below `1/2` under (1.7).

## 3. Spectral laws and deterministic delocalization

### 3.1 Conference children

Since `A_r^2=(r-1)I` and `Tr A_r=0`, the empirical laws of `sA_r` and `tA_r` converge to

```math
\mu_\beta={1\over2}(\delta_{-\beta}+\delta_\beta),
\qquad
\mu_{\beta/\sqrt2}.
\tag{3.1}
```

Condition (2.1) is immediate: even powers are scalar matrices, while every odd power is a scalar
multiple of `A_r`, whose scaled off-diagonal entries are `O(r^{-1/2})` and whose diagonal is zero.

### 3.2 Random-bridge parents

The following records separately the global spectral statement and the stronger entrywise
statement needed by Assumption 2.9.

**Lemma 3.1 (Rademacher conference-bridge delocalization).**  Couple independent iid
Rademacher matrices `B_r` over the Paley orders.  For either fixed orientation, almost surely:

1. the empirical spectral law of
   `X_r=beta S_(epsilon,B_r)/sqrt(2r)` converges to (3.2);
2. for every fixed `j>=1` and `eta>0`,

   ```math
   \max_i\left|(X_r^j)_{ii}-{1\over2r}\operatorname{Tr}X_r^j\right|
   +\max_{i\ne l}|(X_r^j)_{il}|
   =O_{j,\eta}(r^{-1/2+\eta});
   \tag{3.2a}
   ```

3. `limsup ||X_r||_op <= 3 beta/sqrt(2)`.

The lemma is pathwise, so its probability-one event turns the random parent sequence into a
deterministic sequence to which Corollary 2.10(c) applies.

For either fixed orientation and almost every infinite Rademacher bridge sequence,

```math
{\beta\over\sqrt{2r}}S_{\epsilon,B}
```

satisfies (2.1) and has limiting law

```math
\boxed{
\Lambda_{\rm par}
=\mu_{\beta/\sqrt2}
 \boxplus {\rm SC}(\beta^2/2),}
\tag{3.2}
```

where `SC(sigma^2)` is the centered semicircle law of variance `sigma^2`.

**Proof of Lemma 3.1.**  Write `X_r=C_r+Y_r`, where

```math
C_r={\beta\over\sqrt{2r}}\operatorname{diag}(A_r,\epsilon A_r),
\qquad
Y_r={\beta\over\sqrt{2r}}
\begin{pmatrix}0&B_r\\B_r^{\mathsf T}&0\end{pmatrix}.
\tag{3.2b}
```

For the empirical law, first replace `B_r` by an iid Gaussian matrix.  Its two singular-vector
matrices are independent Haar orthogonals, so the rectangular asymptotic-freeness theorem
(Proposition 2.5 in Mai--Speicher--Vargas,
[arXiv:1110.1237](https://arxiv.org/abs/1110.1237), based on the rectangular theorem of
Benaych-Georges) gives the operator-valued Pastur equation for `C_r+Y_r`.  The two block traces
have identical deterministic spectral input `mu_(beta/sqrt(2))`; uniqueness therefore reduces
the operator-valued equation to the scalar equation (3.5), proving (3.2) for Gaussian bridges.

For Rademacher bridges, expand the normalized `k`th trace as a sum over closed walks.  Every
leading term pairs occurrences of the same cross-block edge and hence uses only the common mean
zero and variance one; every block of size at least three in the induced edge partition loses at
least one free vertex and is `O_k(r^{-1})` after normalization.  Thus every fixed trace moment has
the same limit as in the Gaussian model.  Its variance is `O_k(r^{-2})` by the same connected-pair
count; applying the argument to sufficiently high even moments gives summable tails.  Hence the
trace moments converge almost surely.  The norm bound below supplies common compact support,
so moment convergence is equivalent to almost-sure weak convergence of the ESD.  This is the
bounded-entry moment proof of item 1; it also makes explicit why a merely formal Gaussian
replacement is not being used.

For item 2, expand a fixed entry of the `j`th power as a sum over open length-`j` walks.  In its
expectation, every random cross edge must occur an even number of times.  Contract paired bridge
edges successively.  A closed contraction gives a scalar trace, while every remaining open
within-block chain is a power of `A_r`; the identity `A_r^2=(r-1)I` reduces it to either a scalar
multiple of `I` or a scalar multiple of `A_r`.  Terms containing a bridge edge four or more times
lose a free index.  After the `r^{-j/2}` normalization this gives, uniformly in the endpoints,

```math
\max_i\left|\mathbb E(X_r^j)_{ii}
 -{1\over2r}\mathbb E\operatorname{Tr}X_r^j\right|=O_j(r^{-1}),
\qquad
\max_{i\ne l}|\mathbb E(X_r^j)_{il}|=O_j(r^{-1/2}).
\tag{3.3}
```

For the centered second moment, join two open-walk expansions.  Every nonzero connected diagram
uses at least one cross-edge identification between the two walks and therefore loses one free
vertex relative to the product of expectations.  This gives

```math
\sup_{a,b}\left\|(X_r^j)_{ab}-\mathbb E(X_r^j)_{ab}\right\|_2
\le C_{j,\beta}r^{-1/2}.
\tag{3.3a}
```

The centered entry is a Rademacher polynomial of degree at most `j`.  Bonami
hypercontractivity, `||P||_p <= (p-1)^(j/2)||P||_2`, optimized at
`p` of order `r^(2 eta/j)`, now gives, for every `eta>0`,

```math
\Pr\left\{
 |(X_r^j)_{il}-\mathbb E(X_r^j)_{il}|>r^{-1/2+\eta}
 \right\}
\le 2\exp(-c_{j,\eta}r^{2\eta/j}).
\tag{3.4}
```

A union bound over the `4r^2` entries makes the right side summable in `r`; Borel--Cantelli,
(3.3), and the fact that the random trace is the average of the diagonal entries prove (3.2a).

For completeness, the scalar reduction of the operator-valued Pastur equation used above is as
follows.  If `g_1,g_2` are the two block traces, the two child spectral laws are identical, so
uniqueness forces `g_1=g_2=g` and

```math
g(z)=G_{\mu_{\beta/\sqrt2}}
 \left(z-{\beta^2\over2}g(z)\right).
\tag{3.5}
```

This is precisely the Pastur equation for (3.2).

Finally, the triangle inequality and `A_r^2=(r-1)I` give

```math
\|X_r\|_{\rm op}
\le {\beta\over\sqrt2}\left(
 \sqrt{1-1/r}+{\|B_r\|_{\rm op}\over\sqrt r}
 \right).
\tag{3.5a}
```

The rectangular Bai--Yin theorem gives `||B_r||_op/sqrt(r) -> 2` almost surely; see, for
example, Theorem 2.1 in Vershynin's
[survey on extreme singular values](https://www.math.uci.edu/~rvershyn/papers/rv-ICM2010.pdf).
This proves item 3 and compact support.  This completes the lemma.  The global-law portion is also the
two-equal-block special case of the matrix Dyson equation for a rectangular random block plus a
deterministic deformation; see Ajanki--Erdos--Krueger,
[“Stability of the matrix Dyson equation and random matrices with correlations”](https://arxiv.org/abs/1604.08188).

Although the parents are random, this creates no mismatch with the deterministic theorem: apply
Corollary 2.10(c) pathwise on the probability-one event just described.  Finally, convergence also
passes to `E_B log Zbar`.  Indeed,

```math
0\le {1\over2r}\log\overline Z_{2r}(S_{\epsilon,B},t)
\le {1\over2}\left\|{\beta\over\sqrt{2r}}S_{\epsilon,B}\right\|_{\rm op},
\tag{3.6}
```

and the normalized Rademacher operator norms are uniformly integrable (indeed they have uniformly
bounded moments of every fixed order by the standard subgaussian norm tail).  This justifies the
expectation interchange separately for each orientation; averaging the two orientations is then
finite and harmless.

## 4. R-transform calculation

For the symmetric Bernoulli law `mu_c`, inversion of

```math
G_{\mu_c}(z)={z\over z^2-c^2}
```

gives

```math
R_{\mu_c}(u)
={\sqrt{1+4c^2u^2}-1\over2u}.
\tag{4.1}
```

Consequently

```math
{1\over2}\int_0^1R_{\mu_c}(u)\,du=\psi(c).
\tag{4.2}
```

Free convolution adds R-transforms, while

```math
R_{{\rm SC}(\beta^2/2)}(u)={\beta^2\over2}u.
\tag{4.3}
```

Equations (2.5), (3.1), and (3.2) therefore yield

```math
\begin{aligned}
{1\over r}\log\overline Z_r(A_r,s)&\longrightarrow\psi(\beta),\\
{1\over r}\log\overline Z_r(A_r,t)&\longrightarrow\psi(\beta/\sqrt2),\\
{1\over2r}\mathbb E_B\log\overline Z_{2r}(S_{\epsilon,B},t)
&\longrightarrow\psi(\beta/\sqrt2)+{\beta^2\over8}.
\end{aligned}
\tag{4.4}
```

The last limit is the same for both orientations.

The exact output identity now gives

```math
\begin{aligned}
{D_r(\beta)\over r}
&={2\log\overline Z_r(A_r,t)
 +r^2\log\cosh t
 -\mathbb E_{\epsilon,B}\log\overline Z_{2r}(S_{\epsilon,B},t)
 \over r}
 \longrightarrow0,\\
{G_r(\beta)\over r}
&\longrightarrow {\beta^2\over4}
 -2\{\psi(\beta)-\psi(\beta/\sqrt2)\}
 =\gamma(\beta).
\end{aligned}
\tag{4.5}
```

This proves the two limits in (1.8).

To verify the strict sign without a Taylor argument, put

```math
L(\beta)=2\psi(\beta)-2\psi(\beta/\sqrt2)-{\beta^2\over4}
=-\gamma(\beta).
```

Since `psi'(c)=c/(1+sqrt(1+4c^2))`, if

```math
u=\sqrt{1+4\beta^2},
\qquad v=\sqrt{1+2\beta^2},
```

then

```math
{L'(\beta)\over\beta}
={2\over1+u}-{1\over1+v}-{1\over2}
={u-v-\beta^2\over2\beta^2}<0,
\tag{4.6}
```

because `u-v=2 beta^2/(u+v)<beta^2`.  As `L(0)=0`, this proves
`gamma(beta)>0` for every `beta>0`.  Near zero,

```math
\gamma(\beta)
={3\over16}\beta^4-{7\over24}\beta^6+O(\beta^8),
\tag{4.7}
```

recovering and making uniform the earlier linear fourth-order shortfall.

## 5. Lower-quantile consequence

Fix an orientation and write

```math
Y_B=\log\overline Z_{2r}(S_{\epsilon,B},t).
```

Flipping one bridge bit changes `Y_B` by at most `2t`.  Hence McDiarmid's inequality gives

```math
\Pr\{Y_B-\mathbb EY_B\le-a r\}
\le\exp\!\left(-{a^2r\over\beta^2}\right).
\tag{5.1}
```

By (4.4)--(4.5), the mean exceeds the target `2 log Zbar_r(A_r,s)` by
`gamma(beta)r+o(r)`.  Taking `a=gamma(beta)/2` proves (1.9), for example with any
`c_beta<gamma(beta)^2/(4 beta^2)` eventually.

This does not rule out a bridge in an exponentially small basin: the output space has
`2^(r^2+1)` points.  It does prove that a fixed quantile, polynomial sample, or ordinary
uniform-output average cannot repair the interface.

## 6. Consequence for compressed joint states

The conclusion is stronger than a no-go theorem for one overlap ansatz.  On these conference
children the exact full joint quantity satisfies

```math
D_{\rm KL}(U\Vert\Pi)=o(r),
```

whereas compensation requires an additional `gamma(beta)r+o(r)`.  Therefore every compressed
overlap, replica, susceptibility, chain-rule, or strong-data-processing state whose output is a
valid lower bound on this same reverse KL is automatically `o(r)` here.  No refinement of such a
lower bound can close the conference interface; even the uncompressed divergence does not.

What remains logically open is different in kind:

1. finite-temperature minimizers may have a different spectral law and a linear reverse KL;
2. a seed-dependent/nonlocal rule may locate exponentially rare bridges and thereby bypass the
   uniform logarithmic average;
3. the high-temperature theorem applies only while the scaled operator norm is below `1/2`, so it
   cannot by itself establish all fixed-`beta` limits needed for the zero-temperature squeeze.

Thus the new theorem supplies a decisive scalable falsifier for the conference reverse-KL route,
not a convergence theorem and not a falsifier of the minimizer-optimized route.
