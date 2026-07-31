# Reveal-overlap obstruction and the exact theorem a positive route needs

Status: proved high-temperature obstruction plus an exact sufficient overlap
statement. The obstruction applies at sufficiently small fixed scaled
temperature and does not decide the growing-temperature ground-state regime.

Use the notation of `soft_cap_composition_audit.md`. During a bridge reveal,
the exact logarithmic cost of the next edge is

```math
I_\gamma(r)=\log\cosh\gamma
+\frac12\log(1-r^2\tanh^2\gamma),                  \tag{1}
```

where `r=<tau x_i y_j>` in the current extended Gibbs law. Uniformly for
`|r|<=1` and `gamma` tending to zero,

```math
I_\gamma(r)=\frac{\gamma^2}{2}(1-r^2)+O(\gamma^4). \tag{2}
```

## 1. Entropy forces quadratic nonpolarization at small fixed `t`

The following statement is enough for the obstruction.

> **High-temperature reveal theorem.** Fix an aspect-ratio interval
> `eta_0<=m/(m+n)<=1-eta_0`. There are constants `t_0,c>0` such that, for
> exact-minimizer internal blocks of orders `m,n`, a uniformly random order
> and uniformly random signs on the bridge satisfy, for every fixed
> `0<t<=t_0`,
> ```math
> \mathbb E\sum_{e\text{ among the first }c mn\text{ reveals}}
>       (1-r_e^2)\ge c mn.                          \tag{3}
> ```
> Consequently their expected quenched bridge cost is at least
> `c_t(m+n)` at `gamma=t/sqrt(m+n)`.

Here constants can be changed from line to line. The proof has three parts.

### 1.1 Every partial random bridge has the correct cap scale

For a fixed prefix containing `k` bridge edges and a fixed Boolean pair
`(x,y)`, its cross energy is a sum of `k` independent Rademacher signs.
Hoeffding's inequality, followed by a union bound over at most `2^(m+n)`
states and at most `mn` prefixes, gives with probability tending to one,
simultaneously for every prefix,

```math
\max_{x,y}|x^{\mathsf T}C_k y|
\le \sqrt{2k\{(m+n)\log2+O(\log(m+n))\}}
=O((m+n)^{3/2}).                                    \tag{4}
```

The two exact-minimizer internal energies have the same upper scale by the
standard all-order construction bound. Hence every Hamiltonian along such a
good reveal path has absolute cap at most `C N^(3/2)`, where `N=m+n`.

### 1.2 The Gibbs law retains linear entropy

Let `mu_k` be the extended law on `(x,y,tau)` at a prefix and let `U` be
uniform on its `2^(N+1)` states. Since its normalized partition function is
at least one,

```math
D(\mu_k\Vert U)
=\gamma\mathbb E_{\mu_k}[\tau H_k]-\log Z_k
\le\gamma\max|H_k|
\le CtN.                                             \tag{5}
```

For a sufficiently small fixed `t_0`, (5) implies

```math
H(\mu_k)\ge c_0N                                     \tag{6}
```

at every prefix of every good path.

### 1.3 Rank-one rate-distortion turns linear entropy into `mn` variance

Map a Gibbs state to the complete matrix of bridge observables

```math
O_{ij}=\tau x_i y_j.                                 \tag{7}
```

Every `O` is a rank-one sign matrix. The map `(x,y,tau) -> O` has fibers of
size four, so (6) gives `H(O)>=c_0N-log4`.

We use the elementary rank-one rate-distortion lemma:

> If `m,n` have bounded aspect ratio and a law on rank-one `m` by `n` sign
> matrices has entropy at least `eta(m+n)`, then, for a constant
> `c(eta)>0`,
> ```math
> \sum_{i,j}\{1-(\mathbb E O_{ij})^2\}\ge c(eta)mn.  \tag{8}
> ```

To prove the lemma, take independent `O,O'`. The left side of (8) is twice
their expected Hamming distance. If this distance were `o(mn)`, almost all
mass would lie in a Hamming ball of radius `o(mn)` around some rank-one
matrix. A ratio of two rank-one matrices is obtained by flipping a set of
rows and a set of columns. A ball of relative radius `rho=o(1)` therefore
has size

```math
\exp\{O(mh(O(\rho))+nh(O(\rho)))\}=\exp(o(m+n)),     \tag{9}
```

up to the simultaneous-complement redundancy. Markov's inequality and the
entropy decomposition into the ball and its complement would give
`H(O)=o(m+n)`, contradicting the hypothesis. Compactness makes the
qualitative constant in (8) uniform.

At each early prefix for which (4) holds, (8) supplies `cmn` total variance
among all bridge observables. Fewer than `(c/2)mn` edges have yet been
revealed, and each can account for at most one unit. Thus the unrevealed
edges still carry at least `(c/2)mn` variance. Conditional on the current
prefix, the next edge in a random reveal order is uniform among them.
The bad-prefix probability is `o(1)` uniformly over these steps, by the same
Hoeffding union bound, and every increment is nonnegative. Summing the
conditional expected variance for the first `(c/2)mn` steps proves (3).
Finally, (2), with `gamma=t/sqrt N`, gives

```math
\mathbb E\sum_e I_\gamma(r_e)
\ge c\gamma^2mn-O(mn\gamma^4)
=\Theta_t(N).                                       \tag{10}
```

Thus entropy does answer the bounded question: at sufficiently small fixed
`t`, a positive fraction of bridge observables stays bounded away from full
polarization on average, and the uncentered bridge cost in (10) remains
extensive. The argument does not lower-bound `sum r_e^2`, so by itself it
does not quantify the negative correction relative to the fully annealed
cost. In the genuinely perturbative regime one instead expects most `r_e`
to be small, making that centered correction subextensive and leaving nearly
the full annealed cost. In either case there is no entropy-only cancellation
to a summable uncentered defect.

### 1.4 The extensive cost is typical, but this does not control the minimum

Fix either orientation of the second internal block and let

```math
Y(C)=\log Z(A,B,C;\gamma)-\log Z(A,B,0;\gamma).
```

Flipping one bridge sign changes every exponent by at most `2 gamma`, hence

```math
|Y(C)-Y(C^{(e)})|\le2\gamma.                        \tag{10a}
```

The final value is independent of the order used to reveal the signs. Once
(3)--(10) give `E_C Y(C)>=c_tN`, McDiarmid's bounded-differences inequality
gives

```math
\Pr\{Y(C)\le(c_t/2)N\}
\le \exp\!\left{-{c_t^2N^2\over8\gamma^2mn}\right}
=\exp\{-\Omega_t(N)\}.                             \tag{10b}
```

The random reveal order is used to prove the expectation lower bound; the
concentration statement itself is over the independent final bridge signs,
conditional on a fixed orientation. It follows that almost every random
bridge has extensive small-`t` soft cost.

This still says nothing strong enough about a *designed* bridge. There are
`2^{mn}=\exp(Theta(N^2))` bridges, while (10b) suppresses the exceptional
fraction by only `exp(-Omega_t(N))`. The exceptional set can therefore still
contain `exp(Theta(N^2))` bridges. Excluding the minimizing bridge would need
a structural theorem or an `exp(-Omega(N^2))` lower tail, not ordinary
bounded differences.

## 2. Why this does not extend automatically to the ground-state scale

The entropy estimate loses force when `t` grows: (5) permits the Gibbs
entropy to fall to `o(N)`. This is not a technical accident. A distribution
concentrated near one rank-one bridge-response matrix can have
`r_e^2=1-o(1)` on nearly every edge. Mutual information alone cannot exclude
that frozen regime.

A direct high-temperature series has the same boundary. Before edge `e` is
revealed,

```math
r_e={\mathbb E_U[O_e\sinh(\gamma H)]
          \over\mathbb E_U[\cosh(\gamma H)]}.        \tag{11}
```

The linear coefficient vanishes because `e` is absent. The next terms count
signed odd paths between its endpoints. Bounding them requires uniform
spectral/path cancellation for every partial bridge. Without such an
invariant, the number of paths grows too quickly for a series uniform in the
ground-state regime. Thus expansion does not bypass the overlap problem.

## 3. Exact calibrated overlap theorem sufficient for a `b_n` recurrence

Let

```math
b_j=M_j^{2/3},\qquad
T_{m,n}=(b_m+b_n)^{3/2},\qquad
\Delta_{m,n}=T_{m,n}-M_m-M_n.                       \tag{12}
```

For balanced orders, `Delta_(m,n)` is positive and of order `N^(3/2)`.
Take a growing scaled parameter `t_N` and
`gamma_N=t_N/sqrt(N)`. The exact theorem needed from the reveal process is
not zero bridge cost, but the calibrated upper bound

```math
\boxed{
\sum_{e\text{ bridge}} I_{\gamma_N}(r_e)
\le \gamma_N\Delta_{m,n}
   +O(t_NN^{1-\delta}).}                            \tag{13}
```

It may hold in expectation over a reveal distribution, since some realized
orientation and bridge then attain the same upper bound. Choose
`t_N>=N^\delta`; the soft-to-ground entropy cost is

```math
{N\log2\over\gamma_N}
=O(N^{3/2-\delta}).                                 \tag{14}
```

Equations (13)--(14), together with `log Z_j<=gamma_N M_j`, give

```math
M_N\le T_{m,n}+O(N^{3/2-\delta}),                   \tag{15}
```

and hence, in the known nonzero `M_N=Theta(N^(3/2))` range,

```math
b_N\le b_m+b_n+O(N^{1-\delta}).                    \tag{16}
```

This is the desired Hammersley-summable defect on geometric scales.

For `gamma_N=o(1)`, (2) shows what (13) demands at overlap level. Ignoring
only the explicitly controllable `O(mn gamma_N^4)` Taylor remainder, it is

```math
\sum_e(1-r_e^2)
\le {2\Delta_{m,n}\over\gamma_N}
   +O\!\left({t_NN^{1-\delta}\over\gamma_N^2}\right). \tag{17}
```

Since the first term is `Theta(N^2/t_N)`, a polynomially growing `t_N`
requires average polarization `1-r_e^2=O(1/t_N)` with the leading constant
calibrated by the two child optima. Merely proving `r_e^2` is large, or merely
proving a free-energy limit, is insufficient.

The Taylor remainder itself requires, for example, `t_N^4` to fit inside the
allowed logarithmic error in (13); this restricts how rapidly `t_N` can grow
but leaves a nonempty range of small powers.

## 4. Is the overlap theorem genuinely weaker?

Statement (13) has fewer quantifiers than unrestricted bridge optimization:
it concerns exact minimizer blocks, one designed reveal law, and Gibbs
gradients rather than every Boolean cut. In that formal sense it could be a
real reduction.

However, the data it asks to control are

```math
r_e=\partial_{J_e}\log Z
```

at every prefix of a growing bridge. This is the full adaptive Gibbs response
on `Theta(N^2)` coordinates. Scalar child free energies, spectra, row sums,
or bounded quotient states do not determine it. Without a new theorem that
compresses these gradients to a bounded-complexity order parameter, (13) is
the soft version of full bridge/cap response and is not demonstrably simpler
than the original problem.

The evidence-based stopping point is therefore sharp:

- entropy and rank-one geometry prove an extensive high-temperature wall;
- a positive ground-state route needs the calibrated growing-temperature
  theorem (13), not generic polarization;
- absent a new compressed overlap invariant, attempting to prove (13) by
  tracking all prefixes is circular rather than a new composition method.
