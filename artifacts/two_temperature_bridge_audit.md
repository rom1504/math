# Entropy-tilted bridge selection

Date: 2026-08-15.

Status: **exact Rényi composition identity, exact probability-weighted reveal
martingale, a linearly-entropic basin criterion, and a proved fixed-tilt
conference no-go interval**.  The construction does not yet prove a
thermodynamic limit or ground-state convergence.

This note evaluates the proposed second thermodynamic level over bridge
signs.  It is genuinely intermediate between uniform bridge averaging and
exact bridge minimization.  Its surviving obligation is precise: produce a
diffuse bridge law with only linear relative entropy and a linear pressure
gain.

## 1. Exact tilted object

Fix child signings `A,D` of orders `m,n`, let `N=m+n`, and put

```math
t={\beta\over\sqrt N}.
\tag{1.1}
```

Let `Omega={+-1}^{mn+1}` contain a relative child orientation `epsilon` and
all bridge signs `B`.  For `o=(epsilon,B)`, write

```math
z_o=\overline Z_N(A,\epsilon D,B;t),
\qquad
L_o=\log z_o,
\tag{1.2}
```

where `Zbar` is the normalized `cosh` partition function used in the
finite-temperature ledger sections.  If `U` is uniform on `Omega`, the exact
annealed calculation gives

```math
a:=\mathbb E_Uz_o
=(\cosh t)^{mn}\overline Z_m(A,t)\overline Z_n(D,t).
\tag{1.3}
```

Define the size-biased output law and bridge soft minimum

```math
\Pi(o)={z_o\over|\Omega|a},
\qquad
\mathcal R_\lambda
=-{1\over\lambda}\log\mathbb E_Ue^{-\lambda L_o}
=-{1\over\lambda}\log\mathbb E_Uz_o^{-\lambda}.
\tag{1.4}
```

## 2. Exact Rényi identity and rare-basin cost

With

```math
D_\alpha(P\Vert Q)
={1\over\alpha-1}\log\sum_oP(o)^\alpha Q(o)^{1-\alpha},
\tag{2.1}
```

substitution into (1.4) proves

```math
\boxed{
\mathcal R_\lambda
=\log a-D_{1+\lambda}(U\Vert\Pi).}
\tag{2.2}
```

The endpoints are

```math
\lim_{\lambda\downarrow0}\mathcal R_\lambda
=\mathbb E_UL
=\log a-D(U\Vert\Pi),
\qquad
\lim_{\lambda\to\infty}\mathcal R_\lambda=\min_oL_o.
\tag{2.3}
```

If an event `G` of uniform probability `p` satisfies `L_o<=ell` on `G`,
then

```math
\boxed{
\mathcal R_\lambda
\le\ell+{1\over\lambda}\log{1\over p}.}
\tag{2.4}
```

Thus an exponentially rare branch costs `log(1/p)`, not `1/p`.  A single
bridge still costs

```math
0\le\mathcal R_\lambda-\min_oL_o
\le{(mn+1)\log2\over\lambda}.
\tag{2.5}
```

This already identifies the limitation: a basin of probability
`exp(-Theta(N^2))` needs `lambda>>N` to have sublinear cost, the same scale
on which (2.5) resolves the exact minimum.

## 3. Probability-weighted reveal martingale

Reveal the `mn+1` output bits in any order and let `F_j` be the first `j`
bits.  Define

```math
V_j(h)=-{1\over\lambda}\log
\mathbb E_U[e^{-\lambda L}\mid\mathcal F_j=h].
\tag{3.1}
```

Then `V_0=R_lambda`, `V_(mn+1)=L`, and

```math
e^{-\lambda V_{j-1}(h)}
={e^{-\lambda V_j(h,+)}+e^{-\lambda V_j(h,-)}\over2}.
\tag{3.2}
```

The tilted one-bit transition is

```math
q_j(b\mid h)
={e^{-\lambda V_j(h,b)}
 \over e^{-\lambda V_j(h,+)}+e^{-\lambda V_j(h,-)}}.
\tag{3.3}
```

A direct two-point calculation gives the exact identity

```math
\boxed{
V_{j-1}(h)
=\mathbb E_{q_j}V_j(h,b)
+{1\over\lambda}D(q_j(\cdot\mid h)\Vert U_1).}
\tag{3.4}
```

Iterating produces both the Gibbs variational formula and the entropy chain
rule:

```math
\boxed{
\begin{aligned}
\mathcal R_\lambda
 &=\min_q\left\{\mathbb E_qL
      +{1\over\lambda}D(q\Vert U)\right\},\\
 &=\mathbb E_{q_\lambda}L
  +{1\over\lambda}\sum_j\mathbb E_{q_\lambda}
    D(q_j(\cdot\mid\mathcal F_{j-1})\Vert U_1).
\end{aligned}}
\tag{3.5}
```

This is the sought probability-weighted mechanism.  It is not yet an
algorithm: evaluating the two futures in (3.2) is full backward dynamic
programming unless another theorem closes them in a compressed state.

## 4. Same-temperature composition criterion

Let

```math
P_n(\beta)=\min_A\log\overline Z_n(A,\beta/\sqrt n).
\tag{4.1}
```

Choose child minimizers at their contracted scaled temperatures.  Equations
(1.3) and (2.2) give

```math
\boxed{
\begin{aligned}
P_N(\beta)\le{}
 &P_m\!\left(\beta\sqrt{m/N}\right)
 +P_n\!\left(\beta\sqrt{n/N}\right)\\
 &+mn\log\cosh(\beta/\sqrt N)
 -D_{1+\lambda}(U\Vert\Pi).
\end{aligned}}
\tag{4.2}
```

Put

```math
\mathcal T_{m,n}(\beta)
=P_m(\beta)-P_m(\beta\sqrt{m/N})
+P_n(\beta)-P_n(\beta\sqrt{n/N}).
\tag{4.3}
```

The exact sufficient inequality is

```math
\boxed{
\mathcal T_{m,n}(\beta)
+D_{1+\lambda_N}(U\Vert\Pi)
\ge mn\log\cosh(\beta/\sqrt N)-C_\beta N^{1-\delta}.}
\tag{4.4}
```

If (4.4) holds uniformly for every sufficiently large comparable split,
then

```math
P_{m+n}(\beta)
\le P_m(\beta)+P_n(\beta)+C_\beta N^{1-\delta}.
\tag{4.5}
```

Together with the elementary uniform bound `0<=P_n(beta)=O_beta(n)`, the
existing balanced-tree argument gives convergence of `P_n(beta)/n`.

There is a concrete basin version.  Suppose every comparable split has a set
`G_(m,n)` satisfying

```math
U(G_{m,n})\ge e^{-C_1N},
\qquad
\sup_{o\in G_{m,n}}L_o
\le P_m(\beta)+P_n(\beta)+C_2N^{1-\eta}.
\tag{4.6}
```

Taking `lambda_N=N^alpha` in (2.4) proves (4.5) with
`delta=min(alpha,eta)`.  A linearly rare basin is therefore enough.  This is
strictly stronger than finding one isolated bridge and weaker than requiring
a positive uniform probability.

## 5. Conference scaling and the fixed-tilt wall

For two balanced conference children of order `r`, let

```math
s={\beta\over\sqrt r},
\qquad
t={\beta\over\sqrt{2r}},
\qquad
T_r=2\log\overline Z_r(A_r,s).
\tag{5.1}
```

The conference pressure theorem already proved in the repository gives

```math
\log a_r-T_r=\gamma(\beta)r+o(r),
\tag{5.2}
```

where

```math
\gamma(\beta)
={\beta^2\over4}-2\psi(\beta)+2\psi(\beta/\sqrt2)>0.
\tag{5.3}
```

Hence the new proposal faces the exact test

```math
\boxed{
\mathcal R_{\lambda,r}-T_r
=\gamma(\beta)r
-D_{1+\lambda}(U\Vert\Pi_r)+o(r).}
\tag{5.4}
```

Flipping one bridge bit changes a fixed-orientation `L` by at most `2t`.
Hoeffding's martingale lemma yields

```math
\log\mathbb E_Ue^{\lambda(\mathbb E_UL-L)}
\le{\lambda^2\beta^2r\over4}.
\tag{5.5}
```

Consequently,

```math
\mathcal R_{\lambda,r}-T_r
\ge\left(\gamma(\beta)-{\lambda\beta^2\over4}\right)r+o(r).
\tag{5.6}
```

Every fixed tilt below `4 gamma(beta)/beta^2` therefore fails.  Equivalently,
for a law `q` on bridges at a fixed orientation,

```math
\mathbb E_UL-\mathbb E_qL
\le\beta\sqrt{rD(q\Vert U)}.
\tag{5.7}
```

For a joint law on the orientation and bridge, the exact right side has the
additional term

```math
{1\over2}\left|
\mathbb E_Bf_+(B)-\mathbb E_Bf_-(B)
\right|,
\tag{5.8}
```

which is `o(r)` along the conference sequence.  Thus a law gaining
`gamma(beta)r+o(r)` must pay

```math
D(q\Vert U)
\ge{\gamma(\beta)^2\over\beta^2}r+o(r).
\tag{5.9}
```

Linear entropy is the minimum possible scale.

## 6. A nonzero interval of fixed tilts does nothing

> **Theorem 6.1 (small disorder-temperature no-go).**  Fix
>
> ```math
> 0<\beta<{\sqrt2\over6}.
> \tag{6.1}
> ```
>
> There is `lambda_0(beta)>0` such that, along the Paley conference
> sequence, every fixed `0<lambda<lambda_0(beta)` satisfies
>
> ```math
> \boxed{
> {1\over r}D_{1+\lambda}(U\Vert\Pi_r)\longrightarrow0,
> \qquad
> {\mathcal R_{\lambda,r}-T_r\over r}
> \longrightarrow\gamma(\beta)>0.}
> \tag{6.2}
> ```

**Proof.**  Choose `delta>0` and `kappa<1/2` so that

```math
{\beta(3+\delta)\over\sqrt2}<\kappa.
\tag{6.3}
```

This is possible exactly in the range (6.1).  For either orientation, set
`f_epsilon(B)=log Zbar_(2r)(S_(epsilon,B),t)`.  It is convex in the real
bridge entries.  On the convex operator-norm set

```math
\mathcal K_r
=\{B:\|tS_{\epsilon,B}\|_{\rm op}\le\kappa\},
\tag{6.4}
```

the dimension-free high-temperature covariance bound used in
`high_temperature_frobenius_pressure_stability.md` gives

```math
\|\nabla_Bf_\epsilon(B)\|_2\le\beta K_\kappa.
\tag{6.5}
```

The supremum of the supporting affine functions over `K_r` is a convex,
`beta K_kappa`-Lipschitz extension `g_epsilon` agreeing with `f_epsilon` on
`K_r`.  Talagrand convex-Lipschitz concentration on the Rademacher cube gives

```math
\log\mathbb E_Be^{-\lambda(g_\epsilon-\mathbb Eg_\epsilon)}
=O_{\beta,\kappa}(\lambda^2).
\tag{6.6}
```

The conference diagonal blocks cost `beta/sqrt(2)` in operator norm, while
the standard rectangular Rademacher norm tail gives

```math
\Pr\{\|B\|_{\rm op}>(2+\delta)\sqrt r\}
\le e^{-c_\beta r}.
\tag{6.7}
```

Equations (6.3) and (6.7) imply
`Pr(B notin K_r)<=exp(-c_beta r)`.  On the cube, `f_epsilon=O(r^(3/2))`, its
extension is `O(r)`, and `f_epsilon>=0`.  Therefore the exceptional event
changes the mean by `o(1)` and contributes at most `exp(-c_beta r)` to the
negative exponential moment.  The conference pressure theorem supplies,
for both orientations,

```math
{1\over r}\mathbb E_Bf_\epsilon(B)
\longrightarrow
h_\beta:=2\psi(\beta/\sqrt2)+{\beta^2\over4}>0.
\tag{6.8}
```

For every fixed `lambda` with `lambda h_beta<c_beta`, Jensen gives the lower
bound and (6.6)--(6.8) give the upper bound in

```math
\mathbb E_Be^{-\lambda f_\epsilon(B)}
=\exp\{-\lambda h_\beta r+o(r)\}.
\tag{6.9}
```

One may take `lambda_0=c_beta/(2h_beta)`.  Averaging the two orientations
does not alter the exponential rate.  Hence
`R_(lambda,r)/r->h_beta`; the exact annealed identity has
`log(a_r)/r->h_beta`, so (2.2) proves the first limit in (6.2), and (5.2)
proves the second. `square`

The theorem rules out an interval, not every fixed tilt.  A successful
fixed tilt would have to cross a genuine positive threshold.

## 7. The known algebraic lift has quadratic entropy

For a conference child `A_r`, the repository's universal double has bridge
block `A_r+I` and reaches the same-temperature target up to `o(r)` pressure.
For fixed child representatives this is one output among
`2^(r^2+1)`.  A point mass costs

```math
D(q\Vert U)=(r^2+1)\log2.
\tag{7.1}
```

Any algebraic orbit of only `exp(o(r^2))` outputs has the same quadratic
entropy rate.  Formula (3.5) then needs `lambda_r>>r` merely to reduce the
entropy payment to `o(r)`, exactly the resolution at which (2.5) becomes
ordinary exact minimization.

The only genuine remaining implementation is therefore a diffuse law
satisfying, for some admissible fixed or slowly growing `lambda`,

```math
\boxed{
D(q_r\Vert U)=O(r),
\qquad
\mathbb E_{q_r}L
+{1\over\lambda}D(q_r\Vert U)
\le T_r+O(r^{1-\delta}).}
\tag{7.2}
```

Equivalently, a target-reaching bridge basin of probability `exp(-O(r))`
would suffice when combined with a growing tilt.  The uniform bridge law,
polynomial sampling, the small fixed-tilt interval, and the currently known
isolated algebraic lift do not provide it.

Even success at one high-temperature `beta` would not prove the ground-state
limit.  The soft-max reduction needs thermodynamic limits for arbitrary
fixed `beta` before sending `beta` to infinity.
