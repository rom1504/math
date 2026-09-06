# Gaussian-smoothed optimized pressure: exact bridge and unresolved obstruction

Status: verified elementary reductions; proposed thermodynamic theorem OPEN.
This note records why generic spin-glass results do not currently finish it.

## Exact sufficient theorem

Define

```math
P_n(t)=\min_{A\text{ signing}}n^{-3/2}\mathbb E_G M(A+tG).
```

For standard independent Gaussian edges, a union bound/mgf calculation gives
`E M(G)/n^(3/2) <= sqrt((1+1/n) log 2)`. Jensen and the triangle inequality give

```math
m_n\le P_n(t)\le m_n+t\sqrt{(1+1/n)\log2}.                 \tag{1}
```

Thus convergence of P_n(t) for every t in ANY sequence decreasing to zero
would imply convergence of m_n. This alone is an equivalent regularization
criterion and is not a reduction of the missing mathematical difficulty.

For beta>0 one can instead use

```math
P_n(t,\beta)=\min_A {1\over\beta n}\mathbb E_G\log\sum_{s,x}
\exp\left({\beta s\over\sqrt n}(q_A(x)+tq_G(x))\right),
```

with uniform error `0<=P_n(t,beta)-P_n(t)<=2 log2/beta`.
The finite Gaussian-deformation theorem in the companion artifact adds
`P_n(t)>=sqrt(1+t^2)(m_n-10n^(-1/6))`.

## Why the familiar Gaussian limit proof does not apply

Guerra--Toninelli prove a thermodynamic limit by interpolating Gaussian
covariances between a full system and independent subsystems:
[primary preprint](https://arxiv.org/abs/cond-mat/0204280).
Here the deterministic interaction q_A remains arbitrary, and optimizing over
A does not remove the deterministic interpolation term. Scaling a block from
1/sqrt(n) to 1/sqrt(N), and filling its cross edges with signs, is exactly the
original missing completion estimate. Gaussian smoothing does not make those
deterministic interpolants signings.

Results for additional conventional order parameters, e.g.
[Chen's primary paper](https://arxiv.org/abs/2401.10223), are plausible sources
only if the arbitrary deterministic signing can be represented with controlled
finite-dimensional conventional order. No such approximation at o(n^(3/2))
Boolean error is known here. Bounded operator norm alone is not finite rank.

## Precise generic-perturbation quantifier obstruction

One might add a small generic mixed Gaussian perturbation to a deterministic
near-minimizer, invoke replica identities/synchronization, and hope to obtain a
common variational formula. Even if this succeeds for a fixed signing sequence,
the deterministic interaction/cavity observable is not determined by the spin
overlap distribution alone.

There is a second, concrete quantifier issue. If A minimizes the EXPECTED
perturbed pressure, a single-edge flip gives an inequality for
`E_G <s x_i x_j>_G`, not for `<s x_i x_j>_G` in a fixed disorder sample. Squaring
the former describes independent-disorder replicas and cannot be substituted
for the same-disorder overlap used by Gaussian integration by parts.
If instead one chooses the perturbation first and optimizes A afterwards, A
depends on that perturbation; standard fixed-base independence hypotheses fail.
Uniformity over all `2^(n choose 2)` bases is not supplied by the usual
subextensive-perturbation concentration argument.

Related primary literature on stochastic stability is
[Contucci--Giardina--Giberti](https://arxiv.org/abs/1101.2858). Its equilibrium
stability statement is not an extremal deterministic-disorder limit theorem.

## Precise next theorem that would be substantive

Find a fixed t>0 SMALL-NOISE comparison which proves approximate block
subadditivity of `n P_n(t,beta)` with o(n) error uniformly in beta on a sequence
going to infinity, and repeat at arbitrarily small t. Alternatively identify
a compact order parameter with a variational formula for P_n(t,beta) that
includes the deterministic base interaction and whose admissible set is
independent of n with all-order realization.

Neither is proved. Merely adding a Gaussian term, stating a Parisi-shaped
functional, or proving identities for a fixed-base subsequence does not meet
this requirement.
