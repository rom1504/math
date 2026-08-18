# A uniform conference-pressure theorem under arbitrary linear-entropy conditioning

**Status.** Task-local theorem and independent verification; no canonical
files are edited.  This note gives explicit quantifiers and constants for the
observation that, at sufficiently small temperature, *every* bridge event of
probability at least `exp(-Cr)` retains the usual conference pressure.  It
also identifies the finite-disorder-tilt range actually supplied by this
argument and explains why it does not prove the all-fixed-tilt conjecture.

## 1. Setup and frozen dependencies

Let `A_r` run through a sequence of symmetric conference signings,

```math
A_r^2=(r-1)I_r,
```

fix an orientation `epsilon in {+-1}`, and let `U_r` denote the uniform law
on sign bridges `B in {+-1}^{r by r}`.  Put

```math
t_r={\beta\over\sqrt{2r}},
\qquad
S_{\epsilon,B}=
\begin{pmatrix}A_r&B\\B^T&\epsilon A_r\end{pmatrix},
\tag{LC.1}
```

and

```math
f_{\epsilon,r}(B)
=\log\left[2^{-2r}\sum_{x,y}
\cosh\{t_r(H_A(x)+\epsilon H_A(y)+x^TBy)\}\right].
\tag{LC.2}
```

The audited uniform-bridge conference theorem and its regular-sector convex
extension give the following input.  Whenever

```math
0<\beta<\sqrt2/6,
\qquad
{\beta(3+\delta)\over\sqrt2}<\kappa<{1\over2}
\tag{LC.3}
```

for some `delta>0`, define

```math
\mathcal K_{\epsilon,r}(\kappa)
=\{B:\|t_rS_{\epsilon,B}\|_{op}\le\kappa\}.
\tag{LC.4}
```

For every fixed `eta>0`, there is
`c=c_(beta,kappa,eta)>0` such that, for both orientations and all large `r`,

```math
U_r\left(
 \{|f_{\epsilon,r}-h_\beta r|>\eta r\}
 \cap\mathcal K_{\epsilon,r}(\kappa)
 \right)
\le2e^{-cr^2},
\tag{LC.5}
```

where

```math
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4}.
\tag{LC.6}
```

Equation (LC.5) is the two-sided version of the archived regular-sector
theorem.  It follows by extending `f` from the convex set (LC.4) as a
dimension-free convex Frobenius-Lipschitz function, using the audited center
`E g=h_beta r+o(r)`, and applying two-sided convex concentration.  The lower
bound on `kappa` in (LC.3) is essential for this centering: an arbitrary
smaller operator ball need not contain the typical bridge sector.

The exact frozen dependencies used here are:

```text
artifacts/conference_fixed_beta_independent_audit.md
sha256 3711678271772953e1ee91455f5296b7d631545dc3a02dcaf6c89164ca5f9671

artifacts/high_temperature_frobenius_pressure_stability.md
sha256 ae92e40cc56c96c8137f174a3c45b8e460e60dce30c7447e4e359f69a4cc50b5

extremal_information/drafts/conference_regular_conditioned_all_tilts.md
sha256 65e93b956549af59cbf5e41585e691747cfa67d844b9b29dbe0c084ceb0c886a

extremal_information/drafts/conference_regular_and_sublevel_gradient_independent_audit.md
sha256 3016e4f619710bf8a16633238bae4255fecd2263b999c5724e95cdf9cbf7235f

extremal_information/drafts/conference_row_halfcube_pressure.md
sha256 d398b66308c63cdc2c4b00850c91ca769aceb00b66b29b5aec6d9a7a5b172c78

extremal_information/audits/conference_row_halfcube_pressure_independent_audit.md
sha256 90a7f06b2cb3cfc2c92c70be81ee99b8d777c69a6c3e079961b744383cc5f7e3
```

The last two sources record and audit the two-sided form (LC.5).  No
large-deviation or random-matrix theorem beyond these inputs is assumed
below; the bridge-norm estimate is proved directly.

## 2. An explicit bridge-norm entropy exponent

Let `N_r` be a `1/4`-net of the unit sphere with `|N_r|<=9^r`.  For fixed
unit vectors `u,v`, Hoeffding gives

```math
U_r\{|u^TBv|>z\}\le2e^{-z^2/2}.
\tag{LC.7}
```

The standard bilinear-net comparison is

```math
\|B\|_{op}\le2\max_{u,v\in N_r}|u^TBv|.
\tag{LC.8}
```

Consequently, for every fixed `L>0`,

```math
\boxed{
U_r\{\|B\|_{op}>L\sqrt r\}
\le2\exp\left\{-\left({L^2\over8}-2\log9\right)r\right\}.}
\tag{LC.9}
```

Define the crude but explicit regularity exponent

```math
C_*(\beta):=I_{net}(\beta)
={1\over8}\left({1\over\sqrt2\beta}-1\right)^2-2\log9.
\tag{LC.10}
```

This is the supremum of the exponent in (LC.9) over constants `L` for which
block triangle inequality can still put the parent strictly below operator
temperature `1/2`.  Indeed,

```math
\|t_rS_{\epsilon,B}\|_{op}
\le{\beta\over\sqrt2}
\left(\sqrt{1-1/r}+{\|B\|_{op}\over\sqrt r}\right).
\tag{LC.11}
```

Thus, if

```math
L<{1\over\sqrt2\beta}-1,
\tag{LC.12}
```

one may choose

```math
{\beta(1+L)\over\sqrt2}<\kappa<{1\over2}.
\tag{LC.13}
```

Since the useful choices below have `L>2`, (LC.13) also implies (LC.3) for
some `delta>0`, so the center in (LC.5) is valid.

## 3. Arbitrary event theorem

### Theorem LC.1 (linear-entropy conditioning is invisible below `I_net`)

Fix `C>=0` and a temperature `beta>0` satisfying

```math
\boxed{C<I_{net}(\beta).}
\tag{LC.14}
```

For either orientation, let `F_r` be an arbitrary sequence of bridge events,
allowed to depend on `A_r`, `epsilon`, and the complete pressure landscape,
such that

```math
U_r(F_r)\ge e^{-Cr}.
\tag{LC.15}
```

Then, uniformly over all such events,

```math
\boxed{
{f_{\epsilon,r}(B)\over r}\longrightarrow h_\beta
\quad\text{in probability and in }L^1
\quad\text{under }U_r(\,\cdot\mid F_r).}
\tag{LC.16}
```

More quantitatively, choose any `L` satisfying

```math
\sqrt{8(C+2\log9)}<L<{1\over\sqrt2\beta}-1
\tag{LC.17}
```

and put

```math
a={L^2\over8}-2\log9-C>0.
\tag{LC.18}
```

For every fixed `eta>0`, there is `c_(beta,L,eta)>0` such that

```math
\boxed{
\sup_{F_r:\,U_r(F_r)\ge e^{-Cr}}
U_r\left\{
 |f_{\epsilon,r}-h_\beta r|>\eta r\mid F_r
\right\}
\le2e^{-ar}+2e^{-c_{\beta,L,\eta}r^2+Cr}}
\tag{LC.19}
```

for all large `r`.

Equivalently, an explicit sufficient temperature depending only on `C` is

```math
0<\beta<\beta_C^{explicit}
:={1\over\sqrt2\{1+\sqrt{8(C+2\log9+1)}\}}.
\tag{LC.20}
```

With the choice `L=sqrt(8(C+2log9+1))`, (LC.19) then has `a=1`.

**Proof.**  Let

```math
G_r(L)=\{\|B\|_{op}\le L\sqrt r\}.
```

Equations (LC.9), (LC.15), and (LC.18) give

```math
U_r(G_r(L)^c\mid F_r)\le2e^{-ar}.
\tag{LC.21}
```

Choose `kappa` as in (LC.13).  By (LC.11), for all large `r`,
`G_r(L)` is contained in `K_(epsilon,r)(kappa)`.  Dividing (LC.5) by
`U_r(F_r)>=e^(-Cr)` therefore gives

```math
U_r\{|f_{\epsilon,r}-h_\beta r|>\eta r,
      G_r(L)\mid F_r\}
\le2e^{-c_{\beta,L,\eta}r^2+Cr}.
\tag{LC.22}
```

Together (LC.21)--(LC.22) prove (LC.19) and convergence in probability.

It remains to check that the exponentially small irregular part cannot
destroy `L^1`.  On `G_r(L)`, the scaled parent has operator norm at most
`kappa`, and hence

```math
0\le f_{\epsilon,r}(B)
\le {1\over2}\kappa\|(x,y)\|_2^2
=\kappa r.
\tag{LC.23}
```

Globally, the elementary complete-sign cap gives

```math
0\le f_{\epsilon,r}(B)
\le t_r(2r^2-r)\le\sqrt2\beta r^{3/2}.
\tag{LC.24}
```

For fixed `eta>0`, split the conditional expectation of
`|f/r-h_beta|` into: error at most `eta` on the regular nondeviation event;
error at most `kappa+h_beta` on the regular deviation event; and error at
most `sqrt(2)beta sqrt(r)+h_beta` on `G_r(L)^c`.  Equations
(LC.21)--(LC.22) make the latter two contributions vanish uniformly in
`F_r`.  Sending `eta` to zero proves uniform `L^1` convergence. `square`

The proof also allows `U_r(F_r)>=exp(-Cr-o(r))`: choose `L` with a fixed
positive margin in (LC.17), and absorb the sublinear term into that margin.

### Corollary LC.2 (max-density laws)

The same conclusions hold for every bridge law `q_r` satisfying

```math
\left\|{dq_r\over dU_r}\right\|_\infty\le e^{Cr}.
\tag{LC.25}
```

Indeed, `q_r(E)<=e^(Cr)U_r(E)` for every event `E`, so (LC.21)--(LC.24)
apply verbatim.  Event conditioning in LC.1 is the special case
`dq/dU=1_F/U(F)`.

## 4. What the theorem says about lower-tail rates

For every fixed `eta>0`, (LC.5), (LC.9), and (LC.11) imply

```math
\boxed{
\limsup_{r\to\infty}{1\over r}
\log U_r\{|f_{\epsilon,r}/r-h_\beta|>\eta\}
\le-I_{net}(\beta).}
\tag{LC.26}
```

To justify the endpoint in (LC.26), first use any `L` strictly below the
right endpoint of (LC.12), obtaining the corresponding exponent from
(LC.9), and then let `L` increase to that endpoint.  In particular the same
bound holds for every fixed lower sublevel
`f<=(h_beta-eta)r`.

This has two equivalent operational consequences.

1. If `C<I_net(beta)`, no event of mass at least `e^(-Cr)` can be contained
   in a fixed lower-pressure sublevel.
2. Under conditioning by any such event, the remaining fixed pressure
   deviations have exponential rate at least `I_net(beta)-C`.

The exponent is only a certified lower bound on the true rate and is crude
because of the `1/4` net.  Most importantly, it is finite.  The argument
does **not** prove

```math
{1\over r}\log U_r\{f\le(h_\beta-\eta)r\}\to-\infty,
\tag{LC.27}
```

which is the superexponential statement needed to anneal every fixed
negative disorder tilt.  An operator-irregular low-pressure sector of
speed `r` remains possible.

## 5. Exact finite-disorder-tilt consequence

For an admissible event `F_r`, define its negative-moment soft minimum

```math
\mathcal R^F_{\lambda,r}
=-{1\over\lambda}
\log E_{U_r(\cdot\mid F_r)}e^{-\lambda f_{\epsilon,r}}.
\tag{LC.28}
```

### Corollary LC.3 (the certified finite-tilt interval)

If

```math
\boxed{C+\lambda h_\beta<I_{net}(\beta),}
\tag{LC.29}
```

then, uniformly over the events in LC.1,

```math
\boxed{{\mathcal R^F_{\lambda,r}\over r}\longrightarrow h_\beta.}
\tag{LC.30}
```

**Proof.**  Jensen and LC.16 give

```math
\mathcal R^F_{\lambda,r}
\le E[f_{\epsilon,r}\mid F_r]=h_\beta r+o(r).
\tag{LC.31}
```

Choose `d` with

```math
\lambda h_\beta<d<I_{net}(\beta)-C
```

and then choose `eta>0` small enough that
`lambda(h_beta-eta)<d`.  The optimized version of (LC.19) gives

```math
U_r\{f_{\epsilon,r}<(h_\beta-\eta)r\mid F_r\}
\le e^{-dr+o(r)}.
\tag{LC.32}
```

Since `f>=0`, splitting the exponential moment at this sublevel yields

```math
E[e^{-\lambda f}\mid F_r]
\le e^{-\lambda(h_\beta-\eta)r}+e^{-dr+o(r)}.
\tag{LC.33}
```

The first term has the smaller exponential rate, so (LC.33) gives
`liminf R^F_(lambda,r)/r>=h_beta-eta`.  Let `eta` tend to zero and combine
with (LC.31). `square`

Thus, for every prescribed finite `C` and finite tilt ceiling `Lambda`, all
events of mass at least `e^(-Cr)` are annealed for all fixed
`0<lambda<=Lambda` once beta is small enough that

```math
C+\Lambda h_\beta<I_{net}(\beta).
\tag{LC.34}
```

The conclusion also holds for a joint uniform orientation variable because
it holds separately, with the same rate, for both orientations.

At a fixed `beta`, however, (LC.29) only supplies a finite tilt interval.
It gives no conclusion when `lambda h_beta>=I_net(beta)-C`, and it cannot be
iterated to cover all fixed `lambda`.  The archived elementary
Laplace/tail equivalence shows why: annealing every fixed tilt is equivalent
to the superexponential lower tail (LC.27), whereas this proof leaves a
finite speed-`r` operator-norm tail.

## 6. A one-sided KL consequence

An ordinary relative-entropy budget does not imply the two-sided typicality
of LC.1: a law can place fixed mass on a much rarer high-pressure event at
only linear average information cost.  Nevertheless, positivity of the
pressure and negative-MGF transport on the regular sector give a useful
one-sided theorem.  It quantifies how much KL is necessary to lower the
mean pressure.

### Theorem LC.4 (regular mass plus transport forces a KL pressure floor)

Assume `I_net(beta)>0`.  Let `q_r` be arbitrary bridge laws and suppose

```math
\limsup_{r\to\infty}{D(q_r\|U_r)\over r}\le C.
\tag{LC.35}
```

Then

```math
\boxed{
\liminf_{r\to\infty}{E_{q_r}f_{\epsilon,r}\over r}
\ge\left(1-{C\over I_{net}(\beta)}\right)_+h_\beta.}
\tag{LC.36}
```

More precisely, for every fixed `a<I_net(beta)`, the right side may be
replaced by `(1-C/a)_+h_beta`; (LC.36) follows by sending `a` upward to
`I_net(beta)`.

**Proof.**  Fix `a<I_net(beta)`.  Choose `L` below the endpoint in (LC.12)
so that

```math
{L^2\over8}-2\log9>a,
```

and choose `kappa` as in (LC.13).  Let

```math
K_r=\mathcal K_{\epsilon,r}(\kappa),
\qquad p_r=U_r(K_r^c).
```

Since `G_r(L)` is contained in `K_r`, (LC.9) gives

```math
p_r\le2e^{-ar}.
\tag{LC.37}
```

Put `theta_r=q_r(K_r^c)`.  Data processing of relative entropy through the
indicator of `K_r` gives

```math
D(q_r\|U_r)\ge
d_{bin}(\theta_r\|p_r)
\ge\theta_r\log(1/p_r)-\log2.
\tag{LC.38}
```

Here the last inequality follows by expanding binary KL and bounding binary
entropy by `log 2`.  Equations (LC.35), (LC.37)--(LC.38) yield

```math
\limsup\theta_r\le {C\over a}.
\tag{LC.39}
```

Suppose first that `C<a`, so a positive regular mass remains.  Let
`mu_r=U_r(.|K_r)` and `nu_r=q_r(.|K_r)`.  The KL chain rule gives

```math
(1-\theta_r)D(\nu_r\|\mu_r)
\le D(q_r\|U_r)=O(r),
\tag{LC.40}
```

and (LC.39) keeps `1-theta_r` bounded below.  Hence
`D(nu_r||mu_r)=O(r)`.

The archived convex extension supplies the dimension-free negative-MGF
bound

```math
\log E_{\mu_r}
 \exp\{-s(f-E_{\mu_r}f)\}
\le C_{\beta,\kappa}s^2+o(s)+o(1)
\qquad(s\ge0),
\tag{LC.41}
```

and `E_(mu_r)f=h_beta r+o(r)`.  To see that conditioning causes no hidden
loss, write `f=g` on `K_r`, use the unconditional inequality
`log E exp[-s(g-Eg)]<=C_(beta,kappa)s^2`, and divide by
`U_r(K_r)=1-e^{-Omega(r)}`; also `E_(mu_r)f-Eg=o(1)` because `g=O(r)` on
the cube and (LC.37) holds.

Entropy duality applied to (LC.41), with `s` of order
`sqrt(D(nu_r||mu_r)+1)`, gives

```math
E_{\nu_r}f
\ge E_{\mu_r}f
-O_{\beta,\kappa}(\sqrt{D(\nu_r\|\mu_r)+1})-o(r)
=h_\beta r-o(r).
\tag{LC.42}
```

Because `f>=0` also on `K_r^c`, (LC.39) and (LC.42) imply

```math
\liminf {E_{q_r}f\over r}
\ge\left(1-{C\over a}\right)h_\beta.
\tag{LC.43}
```

If `C>=a`, the claimed positive-part bound is merely `E_q f>=0`.  Finally,
send `a` upward to `I_net(beta)` to obtain (LC.36). `square`

Define the conference same-temperature shortfall

```math
\gamma(\beta)
={\beta^2\over4}-2\psi(\beta)+2\psi(\beta/\sqrt2)>0,
\qquad
\tau_\beta=h_\beta-\gamma(\beta)=2\psi(\beta).
\tag{LC.44}
```

### Corollary LC.5 (KL cost of reaching the child target)

If a sequence of laws satisfies

```math
\limsup_{r\to\infty}{E_{q_r}f_{\epsilon,r}\over r}
\le\tau_\beta,
\tag{LC.45}
```

then

```math
\boxed{
\liminf_{r\to\infty}{D(q_r\|U_r)\over r}
\ge I_{net}(\beta){\gamma(\beta)\over h_\beta}.}
\tag{LC.46}
```

Indeed, otherwise a subsequence with KL density below the right side, fed
into LC.4 with a slightly larger budget `C`, would force its mean pressure
strictly above `h_beta-gamma(beta)`, contradicting (LC.45).

The small-temperature limit of this certified information price is finite
and positive.  The exact expansions are

```math
h_\beta={\beta^2\over2}-{\beta^4\over16}+O(\beta^6),
\qquad
\gamma(\beta)={3\beta^4\over16}+O(\beta^6),
\tag{LC.47}
```

and

```math
I_{net}(\beta)
={1\over16\beta^2}-{\sqrt2\over8\beta}+{1\over8}-2\log9.
\tag{LC.48}
```

Therefore

```math
\boxed{
\lim_{\beta\downarrow0}
I_{net}(\beta){\gamma(\beta)\over h_\beta}
={3\over128}.}
\tag{LC.49}
```

This is only a necessary KL cost for lowering the *mean*.  It does not give
conditional typicality under a KL budget, nor does it prove the existence
of a law attaining the target at that cost.

## 7. Independent verification and scope

Every normalization in the proof can be checked without the conference
asymptotics:

- the two sphere nets contribute `2r log 9`;
- the net comparison changes `L sqrt(r)` to `L sqrt(r)/2`, producing
  `L^2r/8` in (LC.9);
- the conference diagonal block costs `sqrt(r-1)`, yielding the `1+L` in
  (LC.11)--(LC.13);
- the regular parent has `2r` spins, so its quadratic cap is `kappa r`;
- the unrestricted parent has `2r^2-r` signed edges, yielding (LC.24);
- conditioning costs at most `Cr`, not `Cr^2`, and therefore cannot disturb
  the regular-sector `r^2` tail.

The theorem is strictly stronger, in its smaller temperature interval, than
closing particular affine, halfcube, template, or row-magnitude fibres: the
event may have arbitrary nonproduct correlations and may be selected using
the pressure itself.  It does not identify the true norm-tail exponent, say
anything uniform over all finite entropy constants at one positive
temperature, exclude a speed-`r` favorable sector outside the certified
entropy range, or prove the all-fixed-tilt conference statement.

The result is a pressure theorem rather than a covariance or operator-norm
claim.  The operator norm is used only to enter a sector where the archived
pressure stability theorem is already valid; its exponentially small
complement is paid directly against the conditioning entropy.
