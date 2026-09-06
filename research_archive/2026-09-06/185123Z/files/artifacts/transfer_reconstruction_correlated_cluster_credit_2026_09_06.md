# A mandatory correlated-cluster credit with fixed one-spin means

Date: 2026-09-06. The construction below adds correlations while preserving
every one-spin mean exactly. Its energy gain therefore contains no aligned
local-field error. An elementary fourth-moment estimate controls clipping;
no cumulant expansion or imported high-temperature theorem is used.
The homogeneous result improves a specified product lower bound. Section 4
separately proves a uniform gap above the **optimal** product variational
value for bounded-norm hollow sign matrices at fixed inverse temperature.
No claim of novelty relative to the Ising literature is made.

## 1. Normalization and an exact correlated law

Let `A` be a real symmetric hollow sign matrix of order `N>=2`, with
`||A||op<=C sqrt N`. Write

```math
H_A(x)=\frac12x^TAx,\quad b=\frac\beta{\sqrt N},\quad
Z_A(b)=\sum_{x\in\{-1,1\}^N}e^{bH_A(x)},\quad \beta>0.
```

Fix arbitrary product means `m_i in [-1,1]` and variances `v_i=1-m_i^2`.
For a partition of the coordinates into `U,V`, first sample the spins in
`U` independently with these means. For each `i in V`, put

```math
f_i=\frac1{\sqrt N}\sum_{j\in U}A_{ij}(X_j-m_j),\qquad
s_i=\mathbb Ef_i^2=\frac1N\sum_{j\in U}v_j,
\quad s_* = \max_{i\in V}s_i.
```

The last equality uses the sign entries and the disjointness of `U,V`.
Take `kappa>0`, `L=1/(8kappa)`, and let `c_L(z)=max(-L,min(z,L))`.
Given the spins in `U`, sample the spins in `V` conditionally independently
with means

```math
m_i+t_i,\qquad
t_i=\kappa v_i\{c_L(f_i)-\mathbb Ec_L(f_i)\}.             \tag{1}
```

Then `E t_i=0`, so the resulting law has the original marginal means.
Also

```math
|t_i|\le v_i/4\le(1-|m_i|)/2,
\qquad \mathbb Et_i^2\le\kappa^2v_i^2s_i.                \tag{2}
```

Thus all conditional probabilities are valid, including the deterministic
case `v_i=0`, where `t_i=0`.

Let

```math
C_V=\frac1{\sqrt N}
\left\|\operatorname{diag}(\sqrt{v_V})A_{VV}
                 \operatorname{diag}(\sqrt{v_V})\right\|_{op},
\qquad S=\sum_{i\in V}v_i s_i
=\frac1N\sum_{i\in V,j\in U}v_i v_j.                     \tag{3}
```

In particular `C_V<=C max_{i in V}v_i<=C`.

## 2. Finite nonuniform credit

Assume

```math
64\kappa^2(3s_*+4/N)\le\frac12.                          \tag{4}
```

Then

```math
\log Z_A(b)\ge
\sum_i h\!\left(\frac{1-m_i}{2}\right)+\frac b2m^TAm
+\left[\frac{\beta\kappa}{2}
       -\kappa^2\left(1+\frac{\beta C_V}{2}\right)\right]S.
                                                                  \tag{5}
```

Here `h` is binary entropy with natural logarithms. If the bracket is
positive and both sides of the partition have positive total variance,
this strictly improves the specified product variational value.

Proof of clipping estimate. For a centered sign of mean `m_j`,

```math
\mathbb E(X_j-m_j)^4=v_j(1+3m_j^2)\le4v_j.
```

Independence in `U` therefore gives

```math
\mathbb Ef_i^4\le3s_i^2+4s_i/N.
```

Since `f_i c_L(f_i)>=f_i^2-f_i^2 1_{|f_i|>L}`,

```math
\mathbb E f_i c_L(f_i)
\ge s_i-\mathbb Ef_i^4/L^2
\ge s_i\{1-64\kappa^2(3s_i+4/N)\}\ge s_i/2.            \tag{6}
```

For `s_i=0` this remains valid. Centering the clipped variable does not
change its covariance with `f_i`, because `E f_i=0`.

Proof of entropy estimate. The entropy as a function of the mean is
`e(z)=h((1-z)/2)`, with `e''(z)=-1/(1-z^2)`. By (2), throughout the
segment between `m_i` and `m_i+t_i`,

```math
1-z^2\ge v_i/2.
```

Taylor's theorem, `E t_i=0`, and conditional independence yield

```math
\mathcal H(X)\ge\sum_i e(m_i)-\sum_{i\in V}\frac{\mathbb Et_i^2}{v_i}
\ge\sum_i e(m_i)-\kappa^2 S.                            \tag{7}
```

Zero-variance summands are interpreted as zero. There is no entropy
independence assertion about the marginal law on `V`: it is conditional
entropy that is being calculated.

Proof of energy estimate. All means are preserved. Distinct spins in `U`
remain independent. For a cross-edge, its extra covariance is
`E[(X_j-m_j)t_i]`. For distinct `i,l in V`, it is `E[t_i t_l]`.
Hollowness removes all diagonal terms. Hence the exact energy change is

```math
b\left(\mathbb EH_A(X)-\frac12m^TAm\right)
=\beta\sum_{i\in V}\mathbb E f_i t_i
 +\frac\beta{2\sqrt N}\mathbb E t^TA_{VV}t.
```

Equations (2), (3), and (6) give a lower bound
`(beta kappa/2)S-(beta C_V/2)kappa^2S`. Combine with (7) and the finite
Jensen inequality `log Z>=H(X)+b E H_A(X)` to obtain (5).
No stationarity condition, local-field bound, or maximizing-spin
assumption enters this proof.

## 3. Explicit small-flip corollary

Choose any spin `x` and objective sign `sigma`, gauge-transform by
`sigma diag(x) A diag(x)`, and set all means to `q=1-2delta`.
Write `v=1-q^2=4delta(1-delta)`. If

```math
0<\delta\le1/3,\qquad \beta^2v\le1/48,
\qquad \beta Cv\le1,
\qquad N\ge64\beta^2,                                  \tag{8}
```

then

```math
\log\{Z_A(b)+Z_A(-b)\}
\ge Nh(\delta)+bq^2|H_A(x)|
       +\frac{\beta^2v^2(N^2-1)}{128N}.                 \tag{9}
```

For (9), choose the objective sign that makes the center energy positive
and use a balanced coordinate partition. Take `kappa=beta/4`. In the
homogeneous case one may use either the general clipping radius in (1)
or the slightly larger radius `delta/(2kappa v)`; both preserve valid
conditional probabilities and give the same stated estimate.
Indeed (8) implies

```math
64\kappa^2(3s_*+4/N)
\le12\beta^2v+16\beta^2/N\le1/2,
\qquad C_V\le Cv.
```

The bracket in (5) is then at least `beta^2/32`, and
`S=v^2 |U||V|/N>=v^2(N^2-1)/(4N)`. Take an absolute maximizer for
`|H_A(x)|=Q(A)` if extracting a cap.

For the current weave, `beta=2t/sqrt p+o(1)`, `C<=1/sqrt p+o(1)`.
At `p=31/32`, `t=4`, `delta=2^{-24}`, both strict inequalities in (8)
hold with large slack, and the order requirement holds eventually.
The mandatory correlation credit per `N` is asymptotically
`beta^2v^2/128`. Thus an annealed upper pressure `A_*m^2+o(m^2)` gives

```math
c\le\frac{A_*-p h(\delta)-t^2v^2/32}
              {2t\sqrt p(1-2\delta)^2}.                 \tag{10}
```

This is a mechanism improvement, although the additional numerical gain
at the deliberately tiny existing flip probability is very small. It
does not, by itself, compare with the best nonuniform product law.

## 4. A uniform gap above the optimal product variational value

For fixed `beta,C>0`, define

```math
P_A(\beta)=\max_{m\in[-1,1]^N}
 \left\{\sum_i e(m_i)+\frac\beta{2\sqrt N}m^TAm\right\},
\quad
\kappa_0=\min\left\{\frac1{32},
                    \frac\beta{4(1+\beta C/2)}\right\},
\quad v_0=\operatorname{sech}^2(\sqrt2\,\beta C).
```

For every `N>=4`, every hollow sign matrix satisfying the norm bound obeys

```math
\log Z_A(\beta/\sqrt N)
\ge P_A(\beta)+\frac{\beta\kappa_0v_0^2}{128}N.          \tag{11}
```

Proof. A maximizing product vector exists by compactness and is interior:
moving a boundary coordinate inward gains entropy of order
`epsilon log(1/epsilon)`, while its energy change is only `O(epsilon)`
at this fixed finite matrix. Stationarity therefore gives

```math
m_i=\tanh(\beta u_i),\qquad
u_i=(Am)_i/\sqrt N,
\qquad \frac1N\sum_i u_i^2\le C^2.
```

At least `N/2` coordinates have `|u_i|<=sqrt2 C`, hence `v_i>=v_0`.
Split these good coordinates as evenly as possible between `U,V`, placing
all other coordinates arbitrarily. Their contribution gives

```math
S\ge\frac{v_0^2(N^2-4)}{16N}\ge\frac{v_0^2N}{32}
\quad(N\ge4).
```

For any product means, `s_*<=1`, `C_V<=C`; and `kappa_0<=1/32` implies
`64kappa_0^2(3s_*+4/N)<=5/16<1/2` for `N>=2`. Its other defining bound
makes the bracket in (5) at least `beta kappa_0/4`. Applying (5) at the
optimal product vector proves (11).

Thus the positive correction can be made genuinely irreducible to
product-mean retuning. The universal constant in (11) is intentionally
coarse, and depends on both `beta` and `C`; no bound uniform as temperature
decreases to zero is asserted. The stronger variance-weighted finite
estimate (5) should be used when the actual product means are known.

## 5. Reproducible bounded checks

Run

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python computations/transfer_reconstruction_correlated_cluster_exact_checks_2026_09_06.py
```

The checker enumerates the complete correlated law on four spins for every
hollow sign matrix, verifies unchanged marginals and the exact conditional
energy decomposition, and certifies both homogeneous and nonuniform
entropy-plus-energy gains using outward rational logarithms. It also checks
the smallness conditions at the current asymptotic weave parameters.
