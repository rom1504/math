# An exact non-Gaussian Bellman value: symmetric two-point sources

Date: 2026-09-06. Status: verified consequence of the reconstructed
conditional-copy lower bound and conditional-variance supersolution.
This computes the Gaussian-boundary Bellman certificate for one source family,
not the original signing minimax limit.

Let `nu_a=(delta_(-a)+delta_a)/2`, with a>=0. Let H_t be the common limit
of the decreasing self-transport terminal iterates and the increasing
Gaussian-boundary iterates of the original binary Bellman operator.
Then

```math
H_t(\nu_a)=T_t(\nu_a)=E_t(\nu_a)
=\max_{0\le s\le1}\left\{
h\left({1+s\over2}\right)+g_t(a^2(1-s^2))\right\}-\log2.
```

For a=0 the displayed expression also gives zero. For a>0 a posterior
law is determined by its signed bias s in [-1,1]. Its conditional entropy
is h((1+s)/2) and variance a²(1-s²). The classwise objective is the
average of their sum, minus log2. It cannot exceed the displayed maximum.
For each s>=0, take the equally weighted pair of reflected posteriors
having biases s and -s. Their barycentre is exactly nu_a, so this is an
admissible channel attaining that scalar value.

Crucially these two posteriors have the SAME conditional variance.
Therefore this same channel has exactly the same value in the older
mean-variance envelope E_t. The universal inequalities E_t<=H_t<=T_t
force equality throughout. Monotonicity then also gives

```math
(\mathcal B T_t)(\nu_a)=T_t(\nu_a),
```

because B T<=T, while B T>=B H=H at this source.

This disproves the tentative hypothesis that a non-Gaussian source with
T_t>g_t(Var) must have strict one-step supersolution drift. Such examples
exist: at a=1,t=4, g_t(1)<-log2 while T_t>=-log2. In fact finite t gives
T_t>-log2, since a posterior with error probability epsilon has entropy
`epsilon log(1/epsilon)+O(epsilon)` and g_t of its variance is
`-4ta² epsilon+O(epsilon²)`.

The gap disappears here because optimal residual variance can be chosen
homogeneous across labels, not because the posterior distributions are
Gaussian. More generally, for ANY symmetric finite source, if a
T_t-optimal channel has constant conditional variance, the same squeeze
proves E_t=H_t=T_t at that source. Approximate attainment with vanishing
Jensen loss suffices as well.

## Ternary contrast and continuing diagnostic

For the actual ternary root p=31/32,t=4, replaying the archived posterior
hull diagnostic on a 1501 grid gives two occupied posterior-z orbits near
0.9086667 and 0.9973333, not a single homogeneous orbit. Its root exponent
lower witness is approximately -0.0150522. The older E-envelope Gaussian
candidate has root exponent approximately -0.04318795. These floating
observations do not certify a strict E/T gap or optimal posterior structure.

The next bounded test is direct evaluation of B T at this root using the
full reversal/swap-symmetric pair family. It has only two pair parameters,
but its children are five-point laws and require the full variable-precision
envelope. The test is implemented in
`computations/decisive_bridge_supersolution_gap_2026_09_06.py` using a
Gaussian-mixture finite-grid dual. No finite-grid value is being treated as
a continuum Bellman upper bound.
