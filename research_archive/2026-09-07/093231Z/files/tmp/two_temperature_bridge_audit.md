# Two-temperature bridge audit

Status: **exact Renyi composition theorem, exact probability-weighted reveal
martingale, a fixed-tilt conference no-go interval, and a quadratic-entropy
obstruction for the known algebraic conference lift**.  The tilted object is genuinely intermediate between the
uniform logarithmic bridge average and exact bridge minimization.  It does not,
by itself, cancel the conference defect.  A successful implementation now has
a precise new obligation: construct a diffuse bridge law that gains linear
pressure at only linear relative-entropy cost.

This note uses no files outside `/home/math/quadra`; the two numerical scratch
outputs cited at the end are exploratory and live in `/home/math/quadra/tmp`.

## 1. Setup

Fix child signings `A,D` of orders `m,n`, set `N=m+n`, `K=mn`, and put

```math
t={\beta\over\sqrt N}.
```

Let `Omega={-1,1}^{K+1}` consist of a relative orientation `epsilon` and all
bridge bits `B`.  For `o=(epsilon,B)`, write

```math
z_o=\overline Z_N(A,\epsilon D,B;t),
\qquad L_o=\log z_o.
```

If `U` is uniform on `Omega`, the exact annealed identity in the existing
reverse-KL note is

```math
a:=\mathbb E_U z_o
=(\cosh t)^K\overline Z_m(A,t)\overline Z_n(D,t).
\tag{1.1}
```

The associated size-biased output law is

```math
\Pi(o)={z_o\over |\Omega|a}.
\tag{1.2}
```

For `lambda>0`, define the proposed bridge soft minimum

```math
\mathcal R_\lambda
=-{1\over\lambda}\log\mathbb E_U e^{-\lambda L_o}
=-{1\over\lambda}\log\mathbb E_U z_o^{-\lambda}.
\tag{1.3}
```

One may instead first pair the orientations,

```math
\widehat z_B={z_{+,B}+z_{-,B}\over2},
```

and apply everything below on the `K` bridge bits.  Its annealed mean is still
`a`.  Pairing is convenient for bounded-difference arguments; retaining the
orientation bit can only make the terminal soft selection more flexible.

## 2. Exact Renyi identity

Let `D_alpha` denote Renyi divergence in the convention

```math
D_\alpha(P\Vert Q)
={1\over\alpha-1}\log\sum_o P(o)^\alpha Q(o)^{1-\alpha}.
```

Substituting (1.2) into (1.3) gives the exact identity

```math
\boxed{
\mathcal R_\lambda
=\log a-D_{1+\lambda}(U\Vert\Pi).}
\tag{2.1}
```

Indeed,

```math
\mathbb E_Uz_o^{-\lambda}
=a^{-\lambda}|\Omega|^{-1-\lambda}
  \sum_o\Pi(o)^{-\lambda}.
```

Thus the new proposal does not continue the old reverse-KL bound.  It replaces
`D(U||Pi)` by the strictly stronger higher-order divergence
`D_(1+lambda)(U||Pi)`.  Its endpoints are

```math
\begin{aligned}
\lim_{\lambda\downarrow0}\mathcal R_\lambda
 &=\mathbb E_U L_o
  =\log a-D(U\Vert\Pi),\\
\lim_{\lambda\to\infty}\mathcal R_\lambda
 &=\min_o L_o,\\
\min_oL_o\le\mathcal R_\lambda
 &\le\mathbb E_UL_o.
\end{aligned}
\tag{2.2}
```

For every event `G subset Omega` of uniform probability `p` on which
`L_o<=ell`,

```math
\boxed{
\mathcal R_\lambda\le \ell+{1\over\lambda}\log{1\over p}.}
\tag{2.3}
```

The universal one-point version is

```math
0\le\mathcal R_\lambda-\min_oL_o
\le{\log|\Omega|\over\lambda}
={(mn+1)\log2\over\lambda}.
\tag{2.4}
```

Equation (2.3) is exactly the proposed `log(1/p)` payment.  Notice an
important sign issue: if the rare bridges merely reach the desired target
`ell=T+o(N)`, a basin of probability `exp(-cN)` leaves a positive linear
term `cN/lambda` at fixed `lambda`.  Fixed `lambda` succeeds only if the low
bridges lie linearly *below* the target and can pay that entropy, or if their
basin is subexponential.  Alternatively `lambda` must grow with the order.

## 3. Exact probability-weighted reveal martingale

Order the `mn+1` output bits and let `F_j` be the first `j` revealed bits.
For a partial history `h`, define the conditional value function

```math
V_j(h)=-{1\over\lambda}\log
 \mathbb E_U[e^{-\lambda L_o}\mid F_j=h].
\tag{3.1}
```

Then `V_0=R_lambda`, `V_(mn+1)=L`, and

```math
e^{-\lambda V_{j-1}(h)}
={e^{-\lambda V_j(h,+)}+e^{-\lambda V_j(h,-)}\over2}.
\tag{3.2}
```

Define the tilted reveal kernel

```math
q_j(b\mid h)
={e^{-\lambda V_j(h,b)}
 \over e^{-\lambda V_j(h,+)}+e^{-\lambda V_j(h,-)}}.
\tag{3.3}
```

Direct calculation gives the one-step probability-weighted identity

```math
\boxed{
V_{j-1}(h)
=\mathbb E_{q_j}[V_j(h,b)]
 +{1\over\lambda}D(q_j(\cdot\mid h)\Vert U_1).}
\tag{3.4}
```

The full path law is

```math
q_\lambda(o)
={e^{-\lambda L_o}
 \over |\Omega|\mathbb E_Ue^{-\lambda L_o}},
\tag{3.5}
```

and iteration yields the Gibbs variational formula and its entropy chain rule:

```math
\boxed{
\begin{aligned}
\mathcal R_\lambda
 &=\min_q\left\{\mathbb E_qL
       +{1\over\lambda}D(q\Vert U)\right\}\\
 &=\mathbb E_{q_\lambda}L
   +{1\over\lambda}\sum_j
      \mathbb E_{q_\lambda}
       D(q_j(\cdot\mid F_{j-1})\Vert U_1).
\end{aligned}}
\tag{3.6}
```

Equivalently,

```math
V_j(F_j)+{1\over\lambda}
 \sum_{i\le j}D(q_i(\cdot\mid F_{i-1})\Vert U_1)
\tag{3.7}
```

is a martingale under `q_lambda`.  This is the exact counterpart of the
probability-weighted rare-branch mechanism: conditional entropy is charged
with the probability of the branch under the tilted law, and no `1/p`
factor occurs.

There is also an exact limitation.  The transition (3.3) needs the two
future conditional soft minima `V_j(h,+/-)`.  Without a theorem closing them
in a compressed state, (3.2) is full backwards dynamic programming over all
remaining bridge bits.  The martingale is therefore an identity, not yet a
local bridge-selection rule.

## 4. Same-temperature composition criterion

Choose `A,D` among minimizers at the contracted scaled temperatures, so that

```math
\log\overline Z_m(A,t)=P_m(\beta\sqrt{m/N}),
\qquad
\log\overline Z_n(D,t)=P_n(\beta\sqrt{n/N}).
```

Because the parent minimum is no larger than the minimum over this block
family, and the latter is no larger than its soft minimum, (2.1) gives

```math
\boxed{
\begin{aligned}
P_N(\beta)\le{}
 &P_m(\beta\sqrt{m/N})+P_n(\beta\sqrt{n/N})\\
 &+mn\log\cosh(\beta/\sqrt N)
 -D_{1+\lambda}(U\Vert\Pi_{A,D,t}).
\end{aligned}}
\tag{4.1}
```

Let

```math
\mathcal T_{m,n}(\beta)
=P_m(\beta)-P_m(\beta\sqrt{m/N})
 +P_n(\beta)-P_n(\beta\sqrt{n/N}).
```

Then the exact sufficient statement is

```math
\boxed{
\mathcal T_{m,n}(\beta)
 +D_{1+\lambda_N}(U\Vert\Pi_{A,D,t})
\ge mn\log\cosh(\beta/\sqrt N)-C_\beta N^{1-\delta}.}
\tag{4.2}
```

On comparable splits, (4.2) gives

```math
P_{m+n}(\beta)
\le P_m(\beta)+P_n(\beta)+C_\beta N^{1-\delta},
\tag{4.3}
```

and hence the fixed-`beta` thermodynamic limit by the existing balanced-tree
argument.  This is genuinely weaker than exact bridge minimization whenever
`lambda_N` is finite, and genuinely stronger than the old reverse-KL
criterion because Renyi divergence increases with its order.

A particularly concrete entropy-basin corollary follows from (2.3).
Suppose that for every comparable split there is a set `G_(m,n)` with

```math
U(G_{m,n})\ge e^{-C_1N},
\qquad
\sup_{o\in G_{m,n}}L_o
\le P_m(\beta)+P_n(\beta)+C_2N^{1-\eta}.
\tag{4.4}
```

Taking `lambda_N=N^alpha` gives (4.3) with
`delta=min(alpha,eta)`.  Thus a linearly rare *basin* is sufficient and no
individual bridge must be extracted.  This is a clean, falsifiable theorem
that realizes the proposed rare-event principle.

## 5. Exact conference scaling

For balanced conference children of order `r`, use

```math
s={\beta\over\sqrt r},
\qquad t={\beta\over\sqrt{2r}},
\qquad T_r=2\log\overline Z_r(A_r,s).
```

The fixed-temperature conference theorem in the repository proves

```math
\log a_r-T_r=\gamma(\beta)r+o(r),
\qquad
\gamma(\beta)>0,
\tag{5.1}
```

where

```math
\gamma(\beta)={\beta^2\over4}-2\psi(\beta)
 +2\psi(\beta/\sqrt2).
```

Combining (2.1) and (5.1) gives the exact test for the new proposal:

```math
\boxed{
\mathcal R_{\lambda,r}-T_r
=\gamma(\beta)r
 -D_{1+\lambda}(U\Vert\Pi_r)+o(r).}
\tag{5.2}
```

The old theorem says only

```math
D(U\Vert\Pi_r)=o(r).
```

That does **not** imply the same for `D_(1+lambda)`: higher-order Renyi
divergence can detect exponentially rare lower tails that reverse KL does
not see.  Conversely it supplies no evidence that the right side of (5.2)
is nonpositive.  The two-temperature question is exactly whether

```math
D_{1+\lambda_N}(U\Vert\Pi_r)
\ge\gamma(\beta)r-O(r^{1-\delta}).
\tag{5.3}
```

Thus the entropy rate does not automatically cancel `gamma(beta)`; (5.3) is
the new missing theorem.

There is a useful necessary scale calculation.  For either fixed orientation,
flipping one bridge bit changes `L` by at most `2t`.  Hoeffding's martingale
lemma gives

```math
\log\mathbb E_Ue^{\lambda(\mathbb E_UL-L)}
\le {\lambda^2\over8}r^2(2t)^2
={\lambda^2\beta^2r\over4}.
\tag{5.4}
```

The two orientations have the same conference-scale mean, so (5.4), applied
separately and then combined, implies

```math
\mathcal R_{\lambda,r}-T_r
\ge\left(\gamma(\beta)-{\lambda\beta^2\over4}\right)r+o(r).
\tag{5.5}
```

In particular a fixed tilt below

```math
\lambda<{4\gamma(\beta)\over\beta^2}
\tag{5.6}
```

provably cannot cancel the conference defect.  Equivalently, for every bridge
law `q`, entropy transportation gives

```math
\mathbb E_UL-\mathbb E_qL
\le\beta\sqrt{rD(q\Vert U)}.
\tag{5.7}
```

Any law gaining `gamma(beta)r+o(r)` must therefore pay at least

```math
D(q\Vert U)
\ge {\gamma(\beta)^2\over\beta^2}r+o(r).
\tag{5.8}
```

Linear entropy is the correct minimum scale; `o(r)` entropy cannot work.

### 5.1 A fixed interval of disorder temperatures provably does nothing

The preceding Hoeffding bound can be strengthened from a lower bound to an
asymptotic identity for a nonzero interval of fixed `lambda` values.

**Theorem 5.1 (small-disorder-temperature no-go).**  Fix

```math
0<\beta<{\sqrt2\over6}.
```

There is `lambda_0(beta)>0` such that, along the Paley conference sequence,
for every fixed `0<lambda<lambda_0(beta)`, jointly over the bridge and the two
orientations,

```math
\boxed{
{1\over r}D_{1+\lambda}(U\Vert\Pi_r)\longrightarrow0,
\qquad
{\mathcal R_{\lambda,r}-T_r\over r}
\longrightarrow\gamma(\beta)>0.}
\tag{5.9}
```

Thus moving a fixed positive distance away from reverse KL does not
immediately uncover the rare conference bridge.  A successful fixed tilt, if
one exists, must pass a genuine positive threshold.

**Proof.**  Treat one orientation first and regard

```math
f_\epsilon(B)=\log\overline Z_{2r}(S_{\epsilon,B},t)
```

as a convex function of the `r^2` real bridge entries.  Choose `delta>0` and
`kappa<1/2` such that

```math
{\beta\over\sqrt2}(3+\delta)<\kappa<{1\over2}.
\tag{5.10}
```

Concretely, take any
`0<delta<1/(sqrt(2) beta)-3` and then any `kappa` in the open interval in
(5.10).  This is possible precisely because `beta<sqrt(2)/6`.  On the convex set

```math
\mathcal C_r=\{B:\|tS_{\epsilon,B}\|_{\rm op}\le\kappa\},
```

the dimension-free high-temperature covariance theorem used in
`artifacts/high_temperature_frobenius_pressure_stability.md` gives
`||C_J||_op<=K_kappa` for both Ising orientations `J` and `-J`.
If `w_+` and `w_-` are the two normalized partition-function weights, direct
differentiation gives

```math
\nabla_Bf_\epsilon
=t\{w_+(C_J)_{LR}-w_-(C_{-J})_{LR}\}.
\tag{5.10a}
```

Using
`||C_J||_F<=sqrt(2r)K_kappa` therefore gives

```math
\|\nabla_B f_\epsilon(B)\|_2
\le t\sqrt{2r}K_\kappa
=\beta K_\kappa.
\tag{5.11}
```

So `f_epsilon` is dimension-free Euclidean-Lipschitz on `C_r`.  The supremum
of its supporting affine functions on that convex set is a convex,
`beta K_kappa`-Lipschitz extension `g_epsilon` to all real bridges, agreeing
with `f_epsilon` on `C_r`.  Talagrand's convex-Lipschitz concentration theorem
on the Rademacher cube consequently yields, for fixed `lambda`,

```math
\log\mathbb E_Be^{-\lambda(g_\epsilon-\mathbb Eg_\epsilon)}
=O_{\beta,\kappa}(\lambda^2).
\tag{5.12}
```

The two-sided form used here is the Rademacher convex-extension theorem of
Johnson and Schechtman,
[*Remarks on Talagrand's deviation inequality for Rademacher
functions*](https://arxiv.org/abs/math/9201208): up to universal constants,
`Pr{|g-med(g)|>u}<=4 exp(-u^2/(C Lip(g)^2))`.

Here is the norm-tail step with all scale dependencies exposed.  Let `c_T>0`
be a universal constant in Talagrand's upper-tail inequality for convex
one-Lipschitz functions on a Rademacher cube.  Bai--Yin convergence in
probability gives

```math
\operatorname{med}\|B\|_{\rm op}
\le(2+\delta/2)\sqrt r
\tag{5.13a}
```

for all sufficiently large `r`.  The map `B mapsto ||B||_op` is convex and
one-Lipschitz in Frobenius norm.  Hence

```math
\Pr\{\|B\|_{\rm op}>(2+\delta)\sqrt r\}
\le2\exp(-c_T\delta^2r/4).
\tag{5.13b}
```

Off this event, the deterministic block triangle inequality gives

```math
\|tS_{\epsilon,B}\|_{\rm op}
\le {\beta\over\sqrt2}
 \left(\sqrt{1-1/r}+2+\delta\right)<\kappa.
\tag{5.13c}
```

After increasing the finite starting order, set

```math
c_\beta={c_T\delta^2\over8}>0.
\tag{5.13d}
```

Then (5.13b)--(5.13d) prove

```math
\Pr\{B\notin\mathcal C_r\}\le e^{-c_\beta r}.
\tag{5.13}
```

The comparison of `f_epsilon` and its extension on the exceptional set is
also quantitative.  For every sign bridge,

```math
0\le f_\epsilon(B)
\le t{2r\choose2}=O_\beta(r^{3/2}).
\tag{5.13e}
```

The zero bridge belongs to `C_r`.  On `C_r`, the interaction norm is at most
`kappa`, so `f_epsilon=O_kappa(r)`; every point of `C_r` also has Frobenius
norm `O_(beta,kappa)(r)`.  The supporting-plane construction and (5.11)
therefore give

```math
\sup_{B\in\{-1,1\}^{r^2}}|g_\epsilon(B)|
=O_{\beta,\kappa}(r).
\tag{5.13f}
```

Since `f_epsilon=g_epsilon` on `C_r`, equations (5.13), (5.13e), and
(5.13f) imply

```math
|\mathbb Ef_\epsilon-\mathbb Eg_\epsilon|=o(1).
\tag{5.13g}
```

The fixed-temperature conference theorem supplies

```math
{1\over r}\mathbb E_Bf_\epsilon(B)
\longrightarrow
h_\beta:=2\psi(\beta/\sqrt2)+{\beta^2\over4}>0
\tag{5.14}
```

for both orientations.  Talagrand's two-sided convex concentration, integrated
in the lower-tail direction, makes (5.12) quantitative: its right side is at
most

```math
C_T'(1+\lambda^2\beta^2K_\kappa^2)
\tag{5.14a}
```

for a universal finite `C_T'`.  Split the negative moment over `C_r` and its
complement.  On the first piece use `f_epsilon=g_epsilon`, (5.12), and
(5.13g).  On the complement use `f_epsilon>=0`, hence
`e^(-lambda f_epsilon)<=1`, followed by (5.13).  Whenever
`lambda h_beta<c_beta`, this proves

```math
\mathbb E_Be^{-\lambda f_\epsilon(B)}
=\exp\{-\lambda h_\beta r+o(r)\}.
\tag{5.15}
```

The lower bound here is Jensen; for the upper bound use (5.12) on
`C_r` and bound the exceptional contribution by (5.13).  Take, for example,
`lambda_0=c_beta/(2h_beta)`.  Averaging the two orientations leaves (5.15)
unchanged at exponential scale, so
`R_(lambda,r)/r->h_beta`.  The exact annealed identity also has
`log(a_r)/r->h_beta`; (2.1) proves the Renyi limit, and (5.1) gives the second
limit in (5.9).  This completes the proof.

**Proof audit of the exceptional-set step.**

1. A generic nonasymptotic subgaussian norm constant is not enough near
   `beta=sqrt(2)/6`.  Equations (5.13a)--(5.13b) deliberately use the sharp
   Bai--Yin center `2 sqrt(r)` and only then apply dimension-free convex
   concentration.
2. The implication from the rectangular norm event to the full parent event
   is the deterministic block triangle inequality (5.13c); no asymptotic
   freeness or entrywise regularity enters.
3. The extension is legitimate because `C_r` is a convex operator-norm
   sublevel set and the pressure is globally convex.  Every supporting plane
   from `C_r` lies below the pressure globally, while choosing the support
   point equal to the evaluation point proves equality on `C_r`.
4. The bad-set contribution to a negative moment is especially benign:
   cosh pressure is nonnegative, so it is bounded by the bad-event probability
   without an `exp(O(r^(3/2)))` multiplier.  Polynomial bounds are used only
   for the mean comparison (5.13g), where exponential rarity suffices.
5. The strict inequality `lambda h_beta<c_beta` is essential.  Above this
   threshold the exceptional norm set can dominate the crude split, and
   Theorem 5.1 makes no claim.  It is not a no-go theorem for every fixed
   `lambda`.

## 6. Basin-speed dichotomy

Let `Delta_r=(L-T_r)/r`.  If its lower tail obeyed an LDP at speed `r` with
rate `I`, Varadhan's formula would turn (1.3) into

```math
\lim {\mathcal R_{\lambda,r}-T_r\over r}
=\inf_z\left\{z+{I(z)\over\lambda}\right\}.
\tag{6.1}
```

This makes the mechanism exact:

* a bridge basin of probability `exp(-Cr)` can matter at fixed `lambda`, but
  only when its pressure depth pays `C/lambda`;
* if an order-one normalized pressure improvement has probability
  `exp(-Theta(r^2))`, every fixed `lambda` leaves the typical value unchanged;
* if a target-reaching basin has probability at least `exp(-Cr)`, taking
  `lambda_r=r^alpha` converts its entropy payment into `O(r^(1-alpha))`, a
  summable composition defect;
* if only a quadratic-entropy basin is known, one needs `lambda_r` much larger
  than `r` to obtain `o(r)` error.  By (2.4), that is also the scale at which
  the soft minimum uniformly approximates exact minimization, so no
  complexity has been saved.

The existing McDiarmid theorem gives only the upper bound
`Pr{L<=T_r}<=exp(-c_beta r)`.  It gives neither the needed lower basin count
nor an `r^2`-speed no-go.  Therefore the current conference theorem alone
cannot decide (5.3).

## 7. The known algebraic conference lift has quadratic entropy

The repository's universal conference double supplies an exact test.  If

```math
A_r^2=(r-1)I,
```

then the valid signing

```math
P_r=
\begin{pmatrix}
A_r&A_r+I\\
A_r+I&-A_r
\end{pmatrix}
\tag{7.1}
```

has the same Boolean quadratic values as

```math
H_r=P_r+\operatorname{diag}(-I,I),
\qquad H_r^2=2rI.
\tag{7.2}
```

The diagonal correction in (7.2) has identically zero quadratic value.
Consequently, in the strict high-temperature range, the conference theorem
gives

```math
\log\overline Z_{2r}(P_r,\beta/\sqrt{2r})
=T_r+o(r).
\tag{7.3}
```

This is an explicit target-reaching bridge, but for fixed child
representatives it is one point of a space of size `2^(r^2+1)`.  A law
supported on this bridge has

```math
D(q\Vert U)=(r^2+1)\log2.
\tag{7.4}
```

More generally, every algebraic orbit containing only `exp(o(r^2))` outputs
still pays `(log2-o(1))r^2`.  Substitution in (3.6) shows that this known lift
requires `lambda_r>>r` merely to make the entropy cost `o(r)`.  The universal
soft-min error (2.4) has exactly the same requirement.  Hence the proposed
tilt, when fed the currently known conference construction, simply repackages
exact bridge selection.

The only genuine escape is a new bridge law `q_r` satisfying

```math
\boxed{
D(q_r\Vert U)=O(r),
\qquad
\mathbb E_{q_r}L
 +{1\over\lambda}D(q_r\Vert U)
\le T_r+O(r^{1-\delta}),}
\tag{7.5}
```

for some admissible fixed or slowly growing `lambda`.  It must retain almost
all `r^2` bridge bits entropically while producing a linear collective
pressure gain.  Neither the uniform law nor the deterministic conference
double does this.

## 8. Exact order-four falsifier

At child order `r=2`, the children are conference matrices.  The bridge and
orientation outputs contain a switching representative of every signing of
`K_4`.  Direct classification of their energy histograms shows that the minimum normalized partition function
at parent raw temperature `t=beta/2` is

```math
\min_o z_o={\cosh\beta(1+\cosh\beta)\over2}.
\tag{8.1}
```

The same-temperature child target is

```math
e^{T_2}=\cosh^2(\beta/\sqrt2).
\tag{8.2}
```

For every `beta>0`, (8.1) is strictly larger than (8.2).  Indeed their
difference is one quarter of

```math
g(\beta)=\cosh(2\beta)+2\cosh\beta-1
 -2\cosh(\sqrt2\beta),
```

and, after the zero constant term, the coefficient of `beta^(2k)` is

```math
{4^k+2-2^{k+1}\over(2k)!}>0
\qquad(k\ge1).
```

Thus even `lambda=infinity` cannot meet the same-temperature target at this
finite order.  This does not obstruct an asymptotic recurrence, but it is an
exact immediate test showing that tilting cannot manufacture cancellation
when the underlying block family has none.

## 9. Exploratory finite diagnostics

`tmp/audit_tilted_bridge_small.py` exhaustively enumerates all bridges,
orientations, and switching-normalized minimizing child pairs through total
order seven.  Its floating-point output is
`tmp/tilted_bridge_small.json`.  At `beta=0.25,0.5,1`, every balanced case
through order seven has positive defect even at the best bridge; finite
`lambda` therefore cannot help.  At `beta=2`, the order-six best bridge has
negative defect, and the tilt becomes negative only once `lambda` is large
(`lambda=16` in the recorded grid).  This illustrates rather than proves the
energy-entropy tradeoff.

`tmp/sample_conference_tilt.py` samples 20,000 bridges for each of three
temperatures at conference child order six.  Its seeded Monte Carlo output is
`tmp/conference_tilt_mc.json`.  At `beta=0.2`, all 40,000 orientation-resolved
outputs lie above target and changing `lambda` from `0.25` to `64` moves the
sampled soft minimum by only about `2.6e-6`.  This is typical-basin evidence
only; it cannot see an exponentially rare bridge.

## 10. Research judgment

The quantum rare-branch lesson maps correctly to the bridge problem at the
level of (3.4)--(3.7).  It produces a valid new composition interface,
namely (4.2), and the linearly-entropic basin corollary (4.4).  It does not
remove the hard mathematical step.

Theorem 5.1 goes beyond the earlier reverse-KL obstruction: an entire nonzero
interval of higher Renyi orders still has `o(r)` divergence and retains the
full `gamma(beta)r` defect.  This does not cover arbitrarily large fixed
`lambda`, but it prevents an infinitesimal or modest unquantified tilt from
being counted as an escape.

For conference children the hard step is now sharply different from the
closed uniform reverse-KL route:

> Construct a nonproduct, sequentially revealable bridge law that changes
> the pressure by `gamma(beta)r` while spending only `O(r)` entropy, or prove
> that every law achieving that pressure gain has entropy `omega(r)` (ideally
> `Theta(r^2)`).

The latter would close every fixed-disorder-temperature implementation.  The
former would be primary progress and could yield a summable recurrence by
letting `lambda_r` grow mildly.  Merely inserting the known algebraic
conference bridge is not progress: its quadratic entropy forces the same
resolution as full bridge minimization.
