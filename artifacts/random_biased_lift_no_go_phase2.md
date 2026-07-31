# Biased random complete lifts: exact scaling and entropy no-go

Date: 2026-07-31.  This is an agent-authored blank-slate derivation.  The
comparison with prior project routes is deferred to Section 5.

## 1. Definition and the scaling that composition requires

Let `A` be an order-`n` signing and write

```math
p=\operatorname{cap}(A)
 =\max_{x\in\{\pm1\}^n}\left|\sum_{i<j}a_{ij}x_ix_j\right|.
```

Fix an integer `k>=2`, a bias `mu in [0,1]`, and put `N=nk`.  Replace
each macro vertex `i` by a fibre `V_i` of `k` clones.  Independently for
each `i<j` and `u in V_i,v in V_j`, choose a sign `C_(uv)` with

```math
\mathbb E C_{uv}=\mu a_{ij}.                            \tag{L.1}
```

Fill every diagonal fibre with any fixed order-`k` signing `D`.  Call the
resulting order-`N` signing `L_(A,k,mu)`.

For a lifted spin `z`, put `s_i=sum_(u in V_i)z_u`.  The inter-fibre
mean is exactly

```math
\mathbb E H_L(z)=\mu\sum_{i<j}a_{ij}s_is_j.             \tag{L.2}
```

The multilinear polynomial on the right attains its maximum absolute
value on the corners of `[-k,k]^n`.  Hence its cap is `mu k^2 p`.  The
within-fibre contribution has cap at most `n cap(D)`.

The ideal `k`-fold composition law in the campaign's
`b_n=M_n^(2/3)` scale would have energy target

```math
(k p^{2/3})^{3/2}=k^{3/2}p.                            \tag{L.3}
```

Thus the bias that retains the base at precisely the compositional scale
is

```math
\mu=k^{-1/2}.                                           \tag{L.4}
```

Constant blocks have `mu=1` and inflate (L.3) by `sqrt(k)`.  Unbiased
blocks have no macro channel at all.

## 2. A rigorous upper bound and its leading defect

Let

```math
E={n\choose2}k^2,\qquad
\sigma^2=1-\mu^2,\qquad
\Lambda=(N+2)\log2.
```

For every `A,k,mu,D`, some realization of the random lift satisfies

```math
\boxed{
\operatorname{cap}(L_{A,k,\mu})
\le \mu k^2p+n\operatorname{cap}(D)
 +\sqrt{2\sigma^2E\Lambda}+{4\over3}\Lambda .}
                                                               \tag{L.5}
```

Indeed, for a fixed lifted spin, subtract (L.2).  The resulting sum has
variance `sigma^2 E` and independent summands bounded by `2`.  Bernstein's
inequality at the displayed threshold has two-sided failure probability
at most `2 exp(-Lambda)`.  A union bound over the `2^N` lifted spins is
strictly below one.  The deterministic within-fibre cap then proves
(L.5).

At the required bias (L.4), `sigma^2=1-1/k`, and

```math
\sqrt{2\sigma^2E\Lambda}=\Theta(N^{3/2}).               \tag{L.6}
```

If `cap(D)=O(k^(3/2))` and `n` tends to infinity, the diagonal-fibre term
is harmless since

```math
{n\operatorname{cap}(D)\over N^{3/2}}=O(n^{-1/2}).     \tag{L.7}
```

The inter-fibre entropy term is not harmless.  From the safe concavity
bound, its contribution to `b_N` is `Theta(N)`, rather than `o(N)` or a
geometrically summable power saving.

More generally, separate mean-plus-fluctuation control has an unavoidable
bias tradeoff.  Keeping the macro term no larger than (L.3), up to a
subleading error, requires `mu<=k^(-1/2)+o(1)` for a fixed `k` and a
base with `p=Theta(n^(3/2))`.  Then `sigma` is bounded away from zero.
Conversely, making the Bernstein term `o(N^(3/2))` by variance reduction
requires `sigma=o(1)`, hence `mu=1-o(1)`; this inflates the macro term by
the leading factor `(1-o(1))sqrt(k)`.

## 3. The centered supremum really is of leading order

The leading scale in (L.6) is not merely an artifact of a crude union
bound.  Let `Xi` be the centered inter-fibre matrix,

```math
\Xi_{uv}=C_{uv}-\mu a_{ij}\qquad(u\in V_i,v\in V_j).
```

Divide the macro fibres into two shores as evenly as possible and let
their clone unions have sizes `r` and `s`.  Set all spins on the first
shore to `1`; for a vertex `v` on the second shore set its spin to the
sign of the centered row sum

```math
R_v=\sum_{u\text{ on the first shore}}\Xi_{uv}.
```

Flipping all second-shore spins reverses the cross term and leaves every
within-shore term fixed.  The larger absolute energy of these two states
is therefore at least

```math
\sum_v |R_v|.                                           \tag{L.8}
```

Each summand in `R_v` is centered, has variance `sigma^2`, and has fourth
moment at most `4sigma^2`.  If `r sigma^2>=1`, the fourth-moment formula
and Paley--Zygmund give

```math
\mathbb E|R_v|\ge {\sigma\sqrt r\over28\sqrt2}.         \tag{L.9}
```

The variables `|R_v|` are independent across `v`, and their variances are
at most `r sigma^2`.  Thus Chebyshev's inequality yields, with probability
at least `1-6272/s`,

```math
\boxed{
\operatorname{cap}(\Xi)
\ge {\sigma s\sqrt r\over56\sqrt2}
=\Omega(\sigma N^{3/2}).}                              \tag{L.10}
```

The constants are intentionally crude.  The point is the exponent: at
the compositional bias (L.4), the centered process itself has a leading
Boolean supremum.  Therefore any proof that controls the macro channel
and the independent residual by separate suprema cannot have a summable
defect.  Equation (L.10) does **not** prove that cancellation between the
two channels can never produce a good individual lift; it proves that
such a result would need a state-dependent absorption/cancellation theorem.

## 4. Exact theorem-level obligation left by the calculation

An independent complete lift can become a convergence mechanism only if
one proves a coupled inequality substantially stronger than (L.5), for
example

```math
\max_z\left|
 \mu\sum_{i<j}a_{ij}s_i(z)s_j(z)+H_\Xi(z)
 \right|
\le k^{3/2}\operatorname{cap}(A)+O(N^{3/2-\delta})       \tag{L.11}
```

at `mu=k^(-1/2)`, uniformly for a structured family of seeds and some
fixed `delta>0`.  The correction must couple the residual to the seed's
soft-spin local fields.  Independence supplies zero conditional mean and
a seed-independent covariance, so (L.11) has no mechanism in the model
defined above.

This is also the falsification criterion for a modified random lift: if
its residual remains conditionally centered with a seed-independent dense
variance floor on the low-slack states, it only renames (L.5).  A viable
modification must exhibit a checkable negative conditional drift, a
low-slack covariance kernel, or another proof-relevant dependence.

## 5. Comparison with the ledger

Only after deriving (L.1)--(L.11), the project history was searched for
random lifts and blow-ups.  The mechanism is not genuinely new:

- `bounded_op_signed_realization.md`, equation (2), proves the exact
  Frobenius residual forced by any sign block whose row sum preserves a
  `sqrt(k)` macro coefficient;
- `regular_microblock_absorption_audit.md` derives exact covariance for
  independent regular-Hadamard blocks and identifies the same need for a
  seed-local-field-dependent correction;
- `correlated_microblock_variance_floor.md` shows that even arbitrary
  cross-edge correlations cannot substantially improve a generic uniform
  variance proxy for exact-sum blocks;
- ledger Section 10.61 uses biased clone lifts successfully for a different
  purpose: transferring a strict finite weighted counterexample.  There the
  macro order is fixed and the signal is `Theta(k^2)`, so microscopic
  `o(k^2)` noise is enough.  That regime does not preserve the
  `n^(3/2)` composition normalization.

The new contribution of this note is a compact exact bias/normalization
calculation, the explicit bound (L.5), and the elementary centered-supremum
lower bound (L.10).  Its strategic classification is a **verified no-go for
independent mean-plus-noise lift proofs**, not primary progress and not a
no-go for dependently corrected microblocks.
