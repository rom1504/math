# Compositional carrier growth: exposed cubes and additive atoms

**Status.** Proof source for the cross-benchmark growth theorem.  Covers are
external unless stated otherwise; logarithms in state counts are base two.

## 1. Stable semantic carriers

A response system has an exact stable (r)-carrier if it has a state
(s(A)\in\mathbb R^r), one-Lipschitz readouts (D_q), and exact updates
(U_C) for every continuation, such that

```math
F_A(q)=D_q(s(A)),\qquad s(CA)=U_C(s(A)).             \tag{1.1}
```

If a test slice of the carrier lies in a translate of ([-R,R]^r), coordinate
quantization gives

```math
\log_2\operatorname {Cov}^{\rm ext}_\varepsilon
\le r\log_2(1+2R/\varepsilon).                     \tag{1.2}
```

Conversely, if actual future queries certify a response embedding

```math
\iota:[-a,a]^k\longrightarrow\mathcal A,
\qquad
d_{\rm rsp}(\iota(u),\iota(v))
\ge\alpha\|u-v\|_\infty,                           \tag{1.3}
```

then grid packing gives

```math
\log_2\operatorname {Cov}^{\rm ext}_\varepsilon
\ge k\log_2\left(1+
 \left\lfloor{\alpha a\over2\varepsilon}\right\rfloor\right)
\quad(0<\varepsilon<\alpha a/2).                   \tag{1.4}
```

Thus robust query exposure lower-bounds semantic carrier rank.  Exact update
in (1.1), not the static covers, is what makes the same state reusable.
Shared-parameter max-affine grammars are a different resource: optimizer
fans let (g) binary parameters generate (exp(O(g^2))) robust response
cells, so semantic rank and grammar rank must not be identified.

## 2. Finite additive response atoms

Let (\phi_1,\ldots,\phi_d\in\ell_\infty(\mathcal Q)).  At mass (n), let

```math
\mathcal H_{n,d}=\{c\in\mathbb N^d:\sum_jc_j=n\},
\qquad
F_c=\sum_jc_j\phi_j.                                \tag{2.1}
```

Histograms compose by addition.  On
(V=\{z\in\mathbb R^d:\sum_jz_j=0\}), define

```math
Tz=\sum_jz_j\phi_j,
\qquad
\alpha=\inf_{z\in V,\ \|z\|_\infty=1}\|Tz\|_\infty,
\qquad
L=\max_j\|\phi_j\|_\infty,                         \tag{2.2}
```

and

```math
\sigma=\inf_{0\ne z\in V\cap\mathbb Z^d}\|Tz\|_\infty. \tag{2.3}
```

The exact arithmetic carrier is

```math
\Gamma_\Phi=
\langle\phi_1-\phi_d,\ldots,\phi_{d-1}-\phi_d\rangle_\mathbb Z,
\qquad r_\mathbb Z=\operatorname {rank}_\mathbb Z\Gamma_\Phi. \tag{2.4}
```

For projective responses these objects and all norms are formed intrinsically
modulo constant query functions.  Anchoring at one query changes norms by up
to a factor two and is not silently treated as an isometry.

### Theorem 1 (additive-generator response law)

For `c,c' in H_(n,d)`, contextual equivalence at fixed mass is exactly

```math
c\sim c'\quad\Longleftrightarrow\quad T(c-c')=0,    \tag{2.5}
```

it is a congruence for every future histogram addition, and the number
(N_n) of exact contextual states obeys

```math
\boxed{N_n=\Theta_\Phi(n^{r_\mathbb Z}).}            \tag{2.6}
```

If (d\ge2) and (\alpha>0), the full histogram is the coarsest exact
state, (r_\mathbb Z=d-1), and, with

```math
s_\varepsilon=1+\left\lfloor{2\varepsilon\over\alpha}\right\rfloor,
```

```math
\left(1+\left\lfloor{n\over(d-1)s_\varepsilon}\right\rfloor\right)^{d-1}
\le\operatorname {Cov}^{\rm ext}_\varepsilon\{F_c\}
\le
\left(2+\left\lceil{L(d-1)n\over\varepsilon}\right\rceil\right)^{d-1}.
                                                               \tag{2.7}
```

If (\sigma>0), then at every (0<\varepsilon<\sigma/2),

```math
\operatorname {Cov}^{\rm ext}_\varepsilon\{F_c:c\in\mathcal H_{n,d}\}
=\binom{n+d-1}{d-1}.                                \tag{2.8}
```

#### Proof

At fixed mass,

```math
F_c=n\phi_d+\sum_{j<d}c_j(\phi_j-\phi_d).           \tag{2.9}
```

This proves (2.5), and adding the same histogram preserves it.  Identify the
finitely generated torsion-free group (\Gamma_\Phi) with
(\mathbb Z^{r_\mathbb Z}).  Its atom generators have bounded integer
coordinates, so all mass-(n) sums lie in an (O_\Phi(n)) box, proving the
upper half of (2.6).  Choose (r_\mathbb Z) displayed atom differences
independent over (\mathbb Q), vary each count from zero to
(\lfloor n/r_\mathbb Z\rfloor), and put the unused mass in type (d).
Their sums are distinct, proving the lower half.

When (\alpha>0), (T|_V) is injective, so the exact state is the
histogram.  For (2.7), vary the first (d-1) counts in multiples of
(s_\varepsilon); two resulting responses are more than
(2\varepsilon) apart.  For the upper bound, externally round the first
(d-1) counts on mesh (\varepsilon/[L(d-1)]) and adjust the last one to
preserve total mass.  Its response error is at most (\varepsilon).
Finally, distinct integer histograms are (\sigma)-separated, proving
(2.8).  `square`

Exact arithmetic rank and robust conditioning are different.  A one-real-
parameter query can assign rationally independent atom responses and have
large (r_\mathbb Z), while (\alpha=0); exact state growth then reflects
unlimited precision rather than robust exposed directions.

## 3. Exact finite-grid mean-field rate

Take equally spaced fields

```math
\gamma_j=-B+(j-1)\Delta,
\qquad \Delta={2B\over d-1},
```

and uniform chemical-potential queries (\lambda\in[-B,B]).  This section
uses the literal anchored sup norm.  One site of
type (j) has response atom

```math
\phi_j(\lambda)=(\gamma_j+\lambda)_+.
```

For (z\in V), at the knots (\lambda=-\gamma_k),

```math
{Tz(-\gamma_k)\over\Delta}
=S_k:=\sum_{j>k}z_j(j-k),
\qquad
\sum_{j\ge k+1}z_j=S_k-S_{k+1}.                    \tag{3.1}
```

Hence (\|z\|_\infty\le4\max_k|S_k|), so
(\alpha\ge\Delta/4).  For nonzero integer (z), some (S_k) is a
nonzero integer; moving one site between adjacent bins attains the resulting
bound.  Therefore

```math
\boxed{\sigma=\Delta.}                              \tag{3.2}
```

In the projective half-oscillation norm, the safe corresponding constants
are `alpha>=Delta/8` and `sigma=Delta/2`.

Theorem 1 now gives the exact microscopic response rate

```math
\operatorname {Cov}^{\rm ext}_\varepsilon
=\binom{n+d-1}{d-1},
\qquad0<\varepsilon<\Delta/2,                       \tag{3.3}
```

and polynomial two-sided covers at all larger scales from (2.7).  This is a
quantitative consequence of the response framework, not an assumption that
histograms are the standard mean-field state.

## 4. Benchmark interpretation

The combined carrier/exposure and atom laws predict the observed growth
regimes.

| benchmark | realizable exposed carrier | congruence | growth |
|---|---|---|---|
| width-(w) pure Max-Cut | arbitrary projective boundary profile | gluing/elimination | exponential in (2^{w-1}) coordinates |
| fixed-width Ising | finite projective boundary profile | max-plus transfer | constant in chain length |
| (d)-type heterogeneous mean field | additive atom histogram | count addition | (\Theta(n^{r_\mathbb Z})) exact states |
| (r)-class lumpable automaton | reachable tropical aggregate | quotient transition | constant in word length |

The law has two independent premises: the response image must be small at
the required scale, and its state must be a congruence for every declared
future.  Max-Cut's universal compiler defeats the first; transition tolls
defeat the second.  Mean-field histograms and tropical lumpability satisfy
both for different algebraic reasons.
