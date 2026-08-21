# Blank-slate direct attack: four architectures and their first obstructions

Date: 2026-08-21.

Status: director report for the campaign started at commit `b5ec773`.  Every
statement labelled **Theorem** below is proved here or reduced to exact finite
arithmetic in `computations/verify_blank_slate_direct_attack.py`.  Asymptotic
saddle calculations not needed by a theorem are labelled **Diagnostic**.

## 1. Problem and normalization

For a hollow symmetric sign matrix `A`, put

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,
\qquad
M_n=\min_A Q(A).
```

The rigorous frontier remains

```math
0.336493364431\ldots
\le \liminf_n {M_n\over n^{3/2}}
\le \limsup_n {M_n\over n^{3/2}}
\le {1\over2}.                                           \tag{1.1}
```

Four proof architectures were frozen before consulting the archive.

| Architecture | First theorem that would move convergence | Outcome |
|---|---|---|
| radial cut-code moment geometry | convergence of the minimum central moment polytope at degrees `d=alpha n`, uniformly as `alpha` grows | exact proxy and sharp Chebyshev recovery proved; the remaining limit is not a strict reduction |
| microcanonical cut-code coverage | a nonperturbative theorem deciding whether the covering multiplicity ever vanishes | raw Bonferroni requires exponentially many replicas; pair data has a transition but does not control holes |
| rectangular/Banach projection | a direct lower bound whose improvement survives on spectrally flat near-minimizers | a new variance correction is exact, but it provably vanishes on every competitive signing |
| exact-minimizer stationarity | local edge-flip conditions force an `O(n^{3/2})` cap or a simpler global structure | scalable one-edge-local and exact radius-two counterexamples kill this implication |

These routes do not use cross-order composition, finite-temperature pressure,
restriction, bridge transport, posterior states, or sparse repair.

## 2. Radial cut-code moment geometry

Let `n>=3` and `N=binom(n,2)`.  Encode a signing by a word `b` in the edge cube.  Let
`C_n^+` be the augmented cut code, consisting of all cuts and their edgewise
complements.  It has `2^n` elements.  For a uniformly random `c in C_n^+`, set

```math
T=d_H(b,c),\qquad q=N-2T.
```

The law `mu_A` of `q` is the antipodal energy law of `A`; in particular

```math
\operatorname{supp}\mu_A\subseteq[-Q(A),Q(A)],
\qquad
\mu_A\{\pm Q(A)\}\ge 2^{1-n}.                           \tag{2.1}
```

For an integer `1<=d<=floor(N/2)`, define `L_d(A)` to be the least `L` for which there is
a probability law `nu`, supported on the correct parity lattice in `[-L,L]`,
such that

```math
\mathbb E_\nu q^{2j}=\mathbb E_{\mu_A}q^{2j}
\quad(0\le j\le d).                                     \tag{2.2}
```

Equivalently, if `K_j^N` is the binary Krawtchouk polynomial, `L_d(A)` is the
least central cap whose Krawtchouk curve has

```math
(s_2(A),\ldots,s_{2d}(A))
```

in its convex hull, where the exact MacWilliams identity gives

```math
s_{2j}(A)
=\mathbb E_{\mu_A}K_{2j}^N(T)
=\sum_{\substack{F\text{ Eulerian}\\|F|=2j}}(-1)^{b\cdot F}. \tag{2.3}
```

Thus `L_d` is a finite linear program in one-dimensional radial moment data.

### 2.1 Recovery from a linear number of moments

**Theorem 2.1 (endpoint Chebyshev recovery).**  If
`w_A=mu_A({-Q(A),Q(A)})`, then

```math
Q(A)\operatorname{sech}\!\left(
 {\operatorname{arcosh}(w_A^{-1/2})\over d}
\right)
\le L_d(A)\le Q(A).                                     \tag{2.4}
```

Consequently, with

```math
\kappa_{n,d}=\operatorname{sech}\!\left(
 {\operatorname{arcosh}(2^{(n-1)/2})\over d}
\right),
\qquad
\Lambda_{n,d}=\min_A L_d(A),
```

one has

```math
\kappa_{n,d}M_n\le\Lambda_{n,d}\le M_n.                \tag{2.5}
```

For `d=floor(alpha n)`,

```math
\kappa_{n,d}\longrightarrow
\operatorname{sech}\!\left({\log2\over2\alpha}\right). \tag{2.6}
```

**Proof.**  The upper bound uses `nu=mu_A`.  For the lower bound suppose that
`nu` in (2.2) is supported in `[-L,L]`.  With `T_d` the Chebyshev polynomial,

```math
P(q)={T_d(q/L)\over T_d(Q(A)/L)}
```

has `P(q)^2` polynomial in `q^2` of degree `d`.  Moment matching gives
`E_nu P^2=E_mu P^2`.  The left side is at most
`T_d(Q(A)/L)^{-2}`, while the endpoint atoms make the right side at least
`w_A`.  Since `T_d(z)=cosh(d arcosh z)` for `z>=1`, solving this inequality
gives (2.4).  Equation (2.1) gives the uniform factor.  QED.

This is stronger than retaining only the top moment.  If

```math
H_d(A)=\left(\mathbb E_{\mu_A}|q|^{2d}\right)^{1/(2d)},
```

then merely

```math
H_d(A)\le L_d(A)\le Q(A)
\le2^{(n-1)/(2d)}H_d(A).                                \tag{2.7}
```

The central moment polytope improves the `1/alpha` recovery error in (2.7) to
`1/alpha^2` in (2.6).

### 2.2 Exact separation certificates

The LP dual is concrete.  For a proposed cap `t`, `L_d(A)>t` exactly when a
polynomial `h(q^2)` of degree at most `d` is nonnegative on the central parity
lattice `|q|<=t` but has `E_mu h<0`.

At order eight, the two gauge-fixed masks `6875` and `6887` have the same
second and fourth moments but caps `14` and `12`.  The exact moments are

```math
\begin{array}{c|rrrr|c}
 &\mathbb E q^2&\mathbb E q^4&\mathbb E q^6&\mathbb E q^8&Q\\ \hline
6875&28&2152&263008&37459072&14\\
6887&28&2152&228448&27782272&12.
\end{array}                                              \tag{2.8}
```

For both masks, `h(z)=z(64-z)` has expectation `-360`, certifying cap bigger
than eight, and `h(z)=z^2(100-z)` has expectations `-47808` and `-13248`,
certifying cap bigger than ten.  On the even energy lattice,

```math
h(z)=-z(z-4)(z-16)(z-144)
```

is nonnegative for `|q|<=12`; its expectations are `-403200` and `3605760`.
It therefore separates the cap-14 mask but correctly does not separate the
cap-12 mask.

### 2.3 Why this is not yet a strict reduction

The exact generating identity

```math
\sum_{j\ge0}s_{2j}(A)z^{2j}
=(1-z^2)^{N/2}\,
\mathbb E_x\cosh\!\left(\operatorname{atanh}(z)H_A(x)\right) \tag{2.9}
```

shows what the linear-degree state contains: signed Eulerian cancellation at
orders proportional to `n`, or equivalently a low-temperature/high-moment
energy transform.  From (2.5), convergence of
`Lambda_{n,floor(alpha n)}/n^{3/2}` for an unbounded sequence
`alpha -> infinity` would imply convergence of `M_n/n^{3/2}`.  More exactly,

```math
{M_n\over n^{3/2}}\text{ converges}
\quad\Longleftrightarrow\quad
\lim_{\alpha\to\infty}
\left(
 \limsup_n{\Lambda_{n,\lfloor\alpha n\rfloor}\over n^{3/2}}
-\liminf_n{\Lambda_{n,\lfloor\alpha n\rfloor}\over n^{3/2}}
\right)=0.                                              \tag{2.10}
```

The reverse implication uses (2.5); the forward implication uses the same
uniform sandwich.  Thus (2.10) is an exact interface, but its present missing
limit statement is quantitatively equivalent to the desired convergence.
The archive's linear-moment/Eulerian calculations already identify the same
signed-cancellation barrier.  The Chebyshev theorem is new; claiming the
remaining limit as a simpler lemma would not be justified.

## 3. Direct microcanonical coverage and its replica barrier

For a radius `r`, define the covering multiplicity

```math
Z_r(b)=|C_n^+\cap B(b,r)|,
\qquad
S_\ell=\sum_b\binom{Z_r(b)}\ell.
```

The order-`K` inclusion--exclusion sum is

```math
B_K=\sum_{\ell=1}^K(-1)^{\ell+1}S_\ell.
```

### 3.1 Exact obstruction to every subexponential truncation

**Theorem 3.1 (Bonferroni multiplicity obstruction).**  Fix `c>0` and an
integer radius sequence of the form below.  Pointwise binomial inversion gives

```math
B_K=|\{b:Z_r(b)>0\}|
-(-1)^K\sum_{Z_r(b)>0}\binom{Z_r(b)-1}{K}.              \tag{3.1}
```

Let

```math
r={N\over2}-{c\over2}n^{3/2}+O(1),
\qquad
\mu=2^{-N}\sum_b Z_r(b).
```

Then

```math
\log\mu=(\log2-c^2)n+o(n).                              \tag{3.2}
```

If the fixed constant satisfies `c<sqrt(log 2)`, an odd Bonferroni upper bound
can be strictly smaller than `2^N` only at order

```math
K\ge\exp((\log2-c^2)n+o(n)).                            \tag{3.3}
```

An even Bonferroni lower bound can equal `2^N` and certify coverage only if

```math
K\ge\max_b Z_r(b)
\ge Z_r(0)
=2^{\,n-(c/\log2)\sqrt n+o(\sqrt n)}.                  \tag{3.4}
```

**Proof.**  Equation (3.1) is the identity

```math
\sum_{\ell=1}^K(-1)^{\ell+1}\binom z\ell
=\mathbf1_{z>0}-(-1)^K\mathbf1_{z>0}\binom{z-1}K.
```

For odd `K`, the pointwise truncated sum satisfies

```math
\mathbf1_{z>0}\left(1+\binom{z-1}K\right)\ge {z\over K}.
```

(For `1<=z<=K` this is immediate; afterward the left-hand discrete slope is
at least one.)  Hence `B_K/2^N>=mu/K`.  If `K<=mu`, the upper bound cannot be
strictly smaller than the universe; this gives (3.3), using the binomial
moderate-deviation estimate (3.2).  For even `K`, equality `B_K=2^N` requires
both full coverage and a zero remainder in (3.1), which is possible only if
`K>=max_b Z_r(b)`.  At `b=0`, ordinary cut words have weight `k(n-k)`.  No
augmented complement contributes for all large `n`, because

```math
\min_k\{N-k(n-k)\}-r={c\over2}n^{3/2}-O(n)>0.
```

The condition `k(n-k)<=r` is

```math
|k-n/2|\ge(\sqrt{c/2}+o(1))n^{3/4}.
```

The corresponding binomial moderate-deviation count is (3.4).  QED.

For a fixed `K`-center intersection, translation by one center reduces its
description to the `2^{K-1}` vertex types relative to the other centers.
More generally, affine rank `s` requires `2^s` types and
`s>=ceil(log_2 K)`.  The necessary order in (3.3) therefore forces
`s=Theta(n)`; (3.4) forces `s>=n-O(sqrt n)`.  Linearity saves the literal
double exponential only by exposing essentially the full cut-code state.
Moreover,

```math
S_\ell=\sum_z h_z\binom z\ell,
\qquad h_z=|\{b:Z_r(b)=z\}|,                             \tag{3.5}
```

is an invertible binomial transform.  Full inclusion--exclusion is exactly
the full coset-multiplicity histogram.

### 3.2 Pair clustering is diagnostic, not coverage

For two code centers with vertex overlap `u`, the leading same-antipodal-class
excess exponent over independent Hamming tails at score
`c n^{3/2}` is

```math
G_c(u)=-I(u)+{2c^2u^2\over1+u^2},
\qquad
I(u)={1\over2}\bigl((1+u)\log(1+u)+(1-u)\log(1-u)\bigr). \tag{3.6}
```

The opposite-antipodal branch has excess

```math
G_c^-(u)=-I(u)-{2c^2u^2\over1-u^2}<0.
```

Since `I(u)>=u^2/2`, the independent saddle is uniquely dominant for
`c<=1/2`; it loses quadratic stability at `c=1/2`, where
`G_{1/2}(u)=-7u^4/12+O(u^6)`, and ceases to be locally maximal above that
value.  This coincidence with the current upper-bound constant is notable,
but second moments control typical multiplicity, not the existence of an
uncovered word.  Theorem 3.1 proves only that no subexponential raw
Bonferroni truncation can bridge that gap.

A concrete next theorem would instead be nonperturbative: control the signed
void generating function

```math
F_n(z)=\mathbb E_b(1-z)^{Z_r(b)},
\qquad F_n(1)=\Pr_b\{Z_r(b)=0\},                         \tag{3.7}
```

at `z=1` without expanding it through subexponentially many factorial moments.
For example, proving `F_n(1)>0` for one fixed `c<1/2` would improve the
all-order upper constant.  No such theorem was found, and this example target
is not presently known to be simpler than the original optimization.

## 4. Rectangular Banach projection and the flat-spectrum barrier

For a sign matrix `B in {+-1}^{m times k}`, set

```math
R_{m,k}=\min_B\max_{u\in\{\pm1\}^m,v\in\{\pm1\}^k}|u^{\mathsf T}Bv|,
\qquad
\mu_k=\mathbb E\left|\sum_{j=1}^k\varepsilon_j\right|.
```

**Theorem 4.1 (exact finite-width rectangular asymptotic).**

```math
m\mu_k\le R_{m,k}\le m\mu_k+k2^k.                       \tag{4.1}
```

If `2^k` divides `m`, equality holds.

**Proof.**  For the rows `r_i` of `B`, maximizing first over `u` gives
`max_v sum_i |r_i dot v|`.  Its average over `v` is `m mu_k`.  Conversely,
use complete equal batches of all `2^k` row types and bound the fewer than
`2^k` residual rows by `k` each.  QED.

For every bipartition `(U,V)` of a symmetric signing `A`, ordinary
polarization of the two full spin configurations `(u,v)` and `(u,-v)` gives

```math
\max_{u,v}|u^{\mathsf T}A_{U,V}v|\le Q(A).              \tag{4.2}
```

Optimizing disjoint side proportions therefore yields only

```math
{Q(A)\over n^{3/2}}
\ge {2\sqrt2\over3\sqrt{3\pi}}-o(1)
=0.3071059\ldots-o(1),                                 \tag{4.3}
```

below the existing lower frontier.

There is an exact second-order gain.  Put

```math
T_{pq}=(B^{\mathsf T}B)_{pq},
\qquad
\eta_k=2^{-(k-2)}\binom{k-2}{\lfloor(k-2)/2\rfloor}.
```

**Theorem 4.2 (degree-two rectangular gain).**  For `m>=1` and `k>=2`,

```math
\|B\|_{\infty\to1}
\ge m\mu_k+{\eta_k^2\over m\mu_k}
\sum_{p<q}T_{pq}^2.                                    \tag{4.4}
```

**Proof.**  For `F(v)=sum_i|r_i dot v|`, the constant Fourier coefficient is
`m mu_k` and its `{p,q}` coefficient is `eta_k T_pq`.  Parseval gives a lower
bound on `E F^2` from these coefficients.  Since `0<=F<=max F`,
`E F^2<=(max F)E F`; rearrange.  QED.

For fixed `alpha,beta>0` with `alpha+beta<=1`, averaging over disjoint row and
column parts of asymptotic proportions `beta` and `alpha` shows that some
partition satisfies, uniformly in `A`,

```math
{Q(A)\over n^{3/2}}
\ge \sqrt{2/\pi}\,\beta\sqrt\alpha
\left(1+\Theta_n(A)\right)-o(1),
\quad
\Theta_n(A)={1\over n^4}\sum_{p<q}(A^2_{pq})^2.         \tag{4.5}
```

The correction cannot move the known frontier.  If `Au=lambda u` with
`||u||_2=1`, independent signs having means
`u_i/||u||_infinity` give the elementary product-rounding bound

```math
Q(A)\ge {|\lambda|\over2\|u\|_\infty^2}
\ge {|\lambda|^3\over2(n-1)}.                          \tag{4.6}
```

To justify (4.5), choose fixed-size disjoint `U,V` uniformly.  Expanding the
squared column correlations and averaging gives, uniformly in `A`,

```math
\mathbb E_{U,V}\sum_{p<q\in V}
\left(\sum_{i\in U}a_{ip}a_{iq}\right)^2
=\alpha^2\beta^2\sum_{p<q}(A^2_{pq})^2+O(n^3).
```

Together with `eta_k^2~2/(pi k)` and
`mu_k~sqrt(2k/pi)`, (4.4) gives (4.5) for at least one partition.

Hence `Q(A)=O(n^{3/2})` implies `||A||_op=O(n^{5/6})`, and

```math
\sum_{p<q}(A^2_{pq})^2
={\operatorname{tr}A^4-n(n-1)^2\over2}
\le {1\over2}\operatorname{tr}A^4
\le {1\over2}\|A\|_{op}^2\operatorname{tr}A^2
=O(n^{11/3}).                                           \tag{4.7}
```

Thus `Theta_n(A)=O(n^{-1/3})` on every competitive signing.  The combination
of product rounding and the new gain disposes of the spectrally concentrated
branch; the gain itself vanishes on the only branch relevant to minimization.
Recovering a leading correction would require information beyond the scalar
spectral defect controlled here; archived cospectral masks already show that
spectrum alone is insufficient.

## 5. Exact-minimizer stationarity and its no-go examples

Switch an exact minimizer so that the all-one spin has energy `q=M_n`.  For a
vertex cut `U`, let `c(U)` be the signed sum of crossed edges.  Then

```math
H_A(\mathbf1-2\mathbf1_U)=q-2c(U),
\qquad 0\le c(U)\le q.                                  \tag{5.1}
```

Flipping a positive edge of an exact minimizer must therefore be witnessed
either by a cut crossing that edge with `c(U)<=1`, or by a cut not crossing
it with `c(U)>=q-1`.  The algebraic roles reverse for a negative edge, but
that statement is vacuous: `U=emptyset` already witnesses every negative edge,
because its flip raises the all-one energy from `q` to `q+2`.  For odd `n`,
`c(U)` is even and `q` has the parity of `N`.  If `n=1 mod 4`, every positive
edge witness therefore has `c=0` or `c=q`; if `n=3 mod 4`, a heavy witness
necessarily has `c=q-1`.

This is a genuine exact-minimizer identity, but local optimality is far too
weak.

**Theorem 5.1 (scalable one-edge local obstruction).**  Let `n=2m`, split the
vertices into equal sets `L,R`, put sign `-1` inside `R`, and `+1` elsewhere.
Then

```math
H_A(x)={X^2+2XY-Y^2\over2},
\qquad X=\sum_{i\in L}x_i,\quad Y=\sum_{j\in R}x_j,
```

The quadratic numerator lies in `[-2m^2,2m^2]`: its maximum as a function of
`X in [-m,m]` occurs at an endpoint, while its minimum occurs at `X=-Y` or an
endpoint.  Hence `Q(A)=m^2=n^2/4`.  Nevertheless, flipping any single edge
does not lower the cap.  A negative `R` edge is protected by the all-one spin, a positive
cross edge by a singleton cut of signed weight one, and a positive `L` edge
by the heavy cut `R` of signed weight `q`.

Even radius two does not recover globality.  At order eight, take negative
edges

```text
16, 26, 27, 35, 37, 45, 48, 68.
```

Its cap is `12`, while the repository's exact value is `M_8=10`.  All 28
single and 378 double edge flips have cap distribution

```text
cap 12: 102,   cap 14: 197,   cap 16: 107.
```

Flipping the triple `13,25,37` reaches cap `10`.  Therefore neither edge
stationarity nor radius-two stability forces even finite global optimality;
the scalable family shows that one-edge stationarity does not force the right
asymptotic scale.

## 6. Archive comparison and director verdict

The archive was consulted only after the four architectures and their first
targets were frozen.

- Linear moments, the Eulerian MacWilliams identity, and their signed
  cancellation barrier already occur in ledger Sections 9.2--9.4.  The
  central moment convex hull and endpoint Chebyshev factor (2.4) do not, but
  their remaining asymptotic obligation is not simpler than convergence.
- Fixed-replica energy/overlap hierarchies were already recognized in
  `entropic_franz_parisi_bernoulli.md`, `eulerian_free_energy_identity.md`,
  and `good_signing_entropy_threshold.md`.  The exact multiplicity identity
  (3.1) and exponential lower bounds (3.3)--(3.4) are new and substantially
  sharpen that obstruction.
- The one-sided switching/cut reformulation and nonnegative Fourier cone were
  already audited in `switching_minimal_graph_extremal.md` and
  `nonnegative_quadratic_fourier_cone.md`.  The exact local counterexamples
  above prevent reviving it through bounded edge stationarity.
- Pure spectral information is already falsified by cospectral order-eight
  masks with different caps.  The degree-two rectangular inequality is new,
  but (4.6)--(4.7) show why it cannot improve the current constant.

The campaign produced several unconditional theorems about `Q(A)` and a
scalable no-go for raw low-replica inclusion--exclusion.  It did not improve
(1.1), produce a strict reduction of convergence, or identify an architecture
whose next lemma is demonstrably weaker than the original rare-extreme
problem.  Under the requested success criterion the campaign is a **STRIKE**,
not a reset.

The best surviving idea is the direct augmented-cut-code coverage problem, but
only in a nonperturbative form that bypasses raw moment/Bonferroni truncation.
No strictly simpler convergence lemma survived.  The exact remaining statement
suggested by the pair transition is:

> For every fixed `epsilon>0`, prove
> `F_{n,1/2-epsilon}(1)=Pr_b{Z_r(b)=0}=0` for all sufficiently large `n` by a
> nonperturbative covering/isoperimetric argument.

Together with the known upper bound this would prove convergence to `1/2`, but
it is currently an exact reformulation, not a strict reduction.  A smaller
discriminating milestone is the same uniform coverage assertion for one fixed
`c>0.336493364431...`, which would improve the rigorous lower frontier.  A
future direct campaign must import a new geometric theorem at this point or
stop.

## 7. Reproduction

Run

```bash
.venv/bin/python computations/verify_blank_slate_direct_attack.py \
  --output computations/results/blank_slate_direct_attack.json
```

The script uses integer arithmetic to verify (2.8), all three radial
separation expectations, the order-eight radius-two obstruction, finite
members of Theorem 5.1, exhaustive tiny instances of Theorem 4.2, and the
complete order-five Bonferroni multiplicity tables at radii two and three.
