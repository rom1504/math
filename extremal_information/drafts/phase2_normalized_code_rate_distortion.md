# Normalized rate--distortion for syndrome-rooted code responses

**Status.**  Independently audited theorem draft following
`phase2_coding_rooted.md`.  The theorem below is proved for unrestricted
appended-fragment covering-radius queries; its concise form is promoted as
Theorem 8.3.

## 1. Question and verdict

Let `G=F_2^w`, let a full-rank parity-check fragment `H` have nonzero column
support `S_H`, and write

```math
\mathcal R_H(E)=\rho(\ker[H\ E]).                  \tag{NRD.1}
```

The exact bit-exposing environments from Theorem CR.2 change this response by
only one.  By themselves they cannot yield a lower bound at additive error
`epsilon w`: the constant prediction `3/2` answers every such special query
within `1/2`.

The unrestricted response experiment is stronger.  A direct-sum family of
syndrome supports makes appended-fragment radii implement **all subset-count
queries** on `Theta(w)` latent bits.  Differences across many blocks then add
before the final absolute error.  This gives a normalized lower bound:

> For every fixed `epsilon<1/8`, there is `c_epsilon>0` such that, for
> infinitely many widths `w`, any deterministic summary answering every
> appended-fragment radius query to uniform additive error `epsilon w`
> requires at least `c_epsilon w` bits.

The fragments and the exposing environments may all be chosen full rank and
of length `Theta_epsilon(w)`.  Thus the result does not use the
`Theta(2^w)`-length family from the exact lattice-scale lower bound.

The lower bound is only linear in `w`, not `2^w`.  It proves that normalized
response complexity remains nonzero, while leaving open whether the complete
support state admits a dramatic approximate compression.

## 2. A block family and its exact response algebra

Fix integers `L>=2` and `q>=1`, put `w=Lq`, and decompose

```math
G=V_1\oplus\cdots\oplus V_q,
\qquad V_j\simeq\mathbb F_2^L.                    \tag{NRD.2}
```

Inside each `V_j`, fix a basis `B_j` and let

```math
D_j=V_j\setminus\{0\}.                            \tag{NRD.3}
```

For `a=(a_1,...,a_q) in {0,1}^q`, define

```math
S_a=\bigcup_{j=1}^q S_{a,j},
\qquad
S_{a,j}=
\begin{cases}
B_j,&a_j=0,\\
D_j,&a_j=1.
\end{cases}                                       \tag{NRD.4}
```

Choose a parity-check fragment `H_a` with this support.  All `H_a` can have
the same length

```math
n=q(2^L-1):                                       \tag{NRD.5}
```

in a sparse block, repeat basis columns while retaining every basis type; in
a dense block, use each member of `D_j` once.  Repetition does not change the
coset-leader profile.  Every `H_a` has rank `w`.

For each subset `P subseteq [q]`, let `E_P` be an appended fragment with
support

```math
T_P=
\left(\bigcup_{j\in P}B_j\right)
\cup
\left(\bigcup_{j\notin P}D_j\right).              \tag{NRD.6}
```

Every `E_P` is full rank.  If desired, repeat its basis columns so all
environments also have the common length (NRD.5).

### Theorem NRD.1 (block subset-count response)

For every `a in {0,1}^q` and `P subseteq[q]`,

```math
\boxed{
\mathcal R_{H_a}(E_P)
=q+(L-1)|\{j\in P:a_j=0\}|.}                     \tag{NRD.7}
```

Consequently, for `a,b in {0,1}^q`, the metric induced by this restricted
family of valid appended fragments is exactly

```math
\begin{aligned}
\max_{P\subseteq[q]}
|\mathcal R_{H_a}(E_P)-\mathcal R_{H_b}(E_P)|
&=(L-1)\max\{N_{01}(a,b),N_{10}(a,b)\},\\
N_{01}(a,b)&=|\{j:a_j=0,b_j=1\}|,\\
N_{10}(a,b)&=|\{j:a_j=1,b_j=0\}|.
\end{aligned}                                      \tag{NRD.8}
```

In particular, the full unrestricted response distance satisfies

```math
\|\mathcal R_{H_a}-\mathcal R_{H_b}\|_\infty
\ge {L-1\over2}\,d_H(a,b).                       \tag{NRD.9}
```

#### Proof

If a generator support is a union of supports lying in the direct summands
`V_j`, its coset-leader profile and radius split:

```math
\lambda(s_1,\ldots,s_q)=\sum_{j=1}^q\lambda_j(s_j),
\qquad
\rho=\sum_{j=1}^q\rho_j.                          \tag{NRD.10}
```

Indeed, every representation of a syndrome separates uniquely into its
block contributions, and independently shortest block representations can
be concatenated.  A basis block has radius `L`, while a complete nonzero
support block has radius one:

```math
\rho(B_j)=L,
\qquad
\rho(D_j)=1.                                      \tag{NRD.11}
```

In `S_a union T_P`, block `j` remains a basis block exactly when `j in P`
and `a_j=0`; every other block contains `D_j`.  Summing (NRD.11) proves
(NRD.7).

Let `Z_a={j:a_j=0}`.  Equation (NRD.7) gives

```math
{1\over L-1}
|\mathcal R_{H_a}(E_P)-\mathcal R_{H_b}(E_P)|
=\bigl||P\cap Z_a|-|P\cap Z_b|\bigr|.             \tag{NRD.12}
```

The supremum over `P` is

```math
\max\{|Z_a\setminus Z_b|,|Z_b\setminus Z_a|\},  \tag{NRD.13}
```

attained by choosing one of the two directed differences.  This proves
(NRD.8).  Their sum is `d_H(a,b)`, so their maximum is at least half the
Hamming distance, proving (NRD.9). `square`

## 3. Uniform normalized information lower bound

Call a deterministic summary `Z(H)` an `eta`-accurate response summary on a
fragment family if there are decoded answers `rhat(Z(H),E)` satisfying

```math
\sup_{H,E}
|\widehat r(Z(H),E)-\mathcal R_H(E)|\le\eta.       \tag{NRD.14}
```

The environments in the supremum may be restricted or unrestricted; a
lower bound for the block environments `(E_P)_P` is automatically a lower
bound for the complete experiment.

### Theorem NRD.2 (finite response-packing bound)

Let `1<=d<=q`.  Any summary satisfying (NRD.14) for the family `(H_a)_a`
and all appended fragments, with

```math
2\eta<(L-1)\left\lceil{d\over2}\right\rceil,     \tag{NRD.15}
```

requires at least

```math
q-\log_2\left(\sum_{i=0}^{d-1}\binom qi\right)   \tag{NRD.16}
```

bits in the worst case.

#### Proof

The greedy Hamming packing argument gives a set

```math
\mathcal A\subseteq\{0,1\}^q,
\qquad
d_H(a,b)\ge d\quad(a\ne b),                      \tag{NRD.17}
```

of size

```math
|\mathcal A|
\ge {2^q\over\sum_{i=0}^{d-1}\binom qi}.         \tag{NRD.18}
```

For two distinct members, (NRD.8) and (NRD.17) give response separation at
least

```math
(L-1)\left\lceil{d\over2}\right\rceil.           \tag{NRD.19}
```

If they shared one summary value, their two decoded answers would agree for
every environment.  The triangle inequality and (NRD.14) would then bound
their response distance by `2 eta`, contradicting (NRD.15).  Hence `Z` is
injective on `mathcal A`, and the logarithm of (NRD.18) proves (NRD.16).
`square`

### Corollary NRD.3 (positive normalized rate)

Fix `0<delta<1/2` and

```math
\epsilon<{\delta(L-1)\over4L}.                   \tag{NRD.20}
```

For `w=Lq` tending to infinity, every deterministic summary uniformly
accurate to `eta=epsilon w` for all unrestricted appended-fragment radius
queries requires at least

```math
\left({1-h_2(\delta)\over L}-o(1)\right)w         \tag{NRD.21}
```

bits, even on the fixed-length family (NRD.5).

Here `h_2` is binary entropy.  To prove the claim, take
`d=ceil(delta q)` in Theorem NRD.2 and use

```math
\sum_{i=0}^{d-1}\binom qi
\le 2^{q(h_2(\delta)+o(1))}.                      \tag{NRD.22}
```

The strict margin in (NRD.20) ensures (NRD.15) for all sufficiently large
`q`.

For every fixed `epsilon<1/8`, choose an integer

```math
L>{1\over1-8\epsilon}                             \tag{NRD.23}
```

and then choose

```math
{4\epsilon L\over L-1}<\delta<{1\over2}.         \tag{NRD.24}
```

Equations (NRD.20)--(NRD.21) give a constant `c_epsilon>0` and the
`c_epsilon w`-bit lower bound stated in Section 1.  Since `L` depends only on
`epsilon`, the common fragment length

```math
q(2^L-1)={2^L-1\over L}\,w                       \tag{NRD.25}
```

is linear in `w`.

The latent block vector itself uses `q=w/L` bits and determines the complete
support `S_a`, so it exactly answers **all** unrestricted appended-fragment
responses, not only the block environments.  Thus, on this restricted state
family under the full query class, the normalized response complexity is
`Theta(w)` bits whenever the lower-bound condition holds.

## 4. Shannon version

The packing also yields a stochastic information lower bound without a
coordinatewise bit decoder.  Let `A` be uniform on a packing `mathcal A`
from (NRD.17), and put

```math
\Delta_*=(L-1)\left\lceil{d\over2}\right\rceil.   \tag{NRD.26}
```

Suppose a transcript `Z` determines a reconstructed complete response vector
`Rhat_Z`, and

```math
\mathbb E\|\widehat R_Z-\mathcal R_{H_A}\|_\infty
\le D.                                             \tag{NRD.27}
```

Nearest-neighbor decoding among the packing responses fails only when the
displayed sup error is at least `Delta_*/2`.  Therefore

```math
p_e\le {2D\over\Delta_*}.                          \tag{NRD.28}
```

Whenever `p=2D/Delta_*<=1/2`, Fano's inequality gives

```math
I(A;Z)
\ge \log_2|\mathcal A|-h_2(p)-p\log_2(|\mathcal A|-1). \tag{NRD.29}
```

This is linear in `w` whenever `D` is a sufficiently small fixed fraction of
`Delta_*` and the packing has positive asymptotic rate.  It is a global
response-vector statement; it does not claim a lower bound under an
arbitrary average distribution on individual environments.

## 5. Why the old support-bit exposure does not scale

For the special environments `E_s` in (CR.11),

```math
\mathcal R_H(E_s)\in\{1,2\}
```

for every fragment `H`.  The zero-bit decoder

```math
\widehat r(E_s)={3\over2}                          \tag{NRD.30}
```

has uniform additive error `1/2`.  For every fixed `epsilon>0`, this is at
most `epsilon w` once `w>=1/(2epsilon)`.  Hence the old one-query-per-support-
bit exposure and its isometric Hamming-cube response vector cannot by
themselves prove any positive normalized rate.

The block environments evade this obstruction by aggregating a chosen set
of missing dense blocks into one radius.  Formula (NRD.7) performs the sum
before the final approximation error is paid.  This is the coding analogue
of joint-channel cancellation/exposure: separately visible unit effects are
useless at scale `epsilon w`, whereas a legal composed query turns linearly
many coherent effects into one order-`w` response difference.

## 6. Constant and scope audit

1. **The factor `1/2` in (NRD.9) is necessary for arbitrary binary
   vectors.**  The two directed difference counts sum to the Hamming
   distance; without imposing a coordinatewise order, only their maximum can
   be selected by one query.

2. **The factor `2` in (NRD.15) is the decoder triangle loss.**  Two true
   responses sharing one decoded response can be separated by at most twice
   the uniform approximation error.

3. **The threshold `1/8` follows from positive-rate binary packing.**  In
   this construction, asymptotically positive binary-code rate requires
   `delta<1/2`; combining this with (NRD.20) and allowing `L` to grow through
   fixed constants gives `epsilon<1/8`.  This is a limitation of the block
   construction, not an impossibility theorem above `1/8`.

4. **All ranks and lengths are legitimate.**  Each state and environment
   contains a spanning basis in every direct summand.  Repeated columns pad
   lengths without changing supports or profiles.  For fixed `epsilon`, both
   state and query lengths are `Theta_epsilon(w)`.

5. **The upper bound is source-restricted.**  On the block family the
   sufficient response state is the `q`-bit vector `a`, while the root-
   distance table has `2^n` labeled entries.  If one chooses one canonical
   `H_a` per vector, this is a succinct source coordinate rather than a
   many-to-one quotient.  One may instead include arbitrary duplicate-column
   realizations with the same support: `a` then forgets those multiplicities
   while every appended-radius response remains unchanged.  The strict
   forgetting phenomenon is independently witnessed in Theorem CR.1.

6. **This does not establish exponential normalized complexity.**  The
   exact state for arbitrary supports uses `Theta(2^w)` bits.  The present
   packing proves only `Omega_epsilon(w)` bits and is sharp only on its
   `q`-bit block subfamily.  An `exp(o(2^w))` approximate quotient for the
   unrestricted support family remains possible.

## 7. Director judgment

This is a genuine normalized extension of the exact support-bit theorem, not
a relabeling of its lattice-scale packing.  The new mechanism is a legal
future fragment that performs subset aggregation across independent syndrome
blocks, creating order-`w` response separation while every individual
support-bit probe remains bounded.

It is also deliberately bounded in claim.  Direct-sum additivity makes the
model tractable, and the resulting response experiment is an exact subset-
count sketching problem.  The theorem establishes a positive extremal
information rate for one nontrivial compositional code family, but it does
not determine the rate--distortion function of arbitrary syndrome supports.

The next precise question, if this direction continues, is whether arbitrary
support profiles admit a subexponential-in-`2^w` net under the complete
future-radius response metric at distortion `epsilon w`, or whether a
superlinear packing can be constructed without relying on an explicit direct
sum decomposition.
