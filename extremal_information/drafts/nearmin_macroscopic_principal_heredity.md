# Near-minimizers are hereditary under every sublinear vertex deletion

Date: 2026-08-17.

Status: proved draft for independent audit.  The result is a near-order
structural theorem, not a fixed-ratio recurrence and not a proof of
convergence.

## 1. Result

For a hollow real symmetric matrix `A`, write

```math
Q(A)=\max_{x\in\{\mathord\pm1\}^{|A|}}|H_A(x)|,
\qquad
H_A(x)=\sum_{i<j}a_{ij}x_ix_j.
```

### Theorem MH.1 (random-bridge composition and principal heredity)

Let `n=m+k`, with `m,k>=1`.  Then

```math
\boxed{
M_n\le M_m+M_k+\sqrt{2(\log2)mkn}.}                \tag{MH.1}
```

Consequently, if `A` is any order-`n` signing satisfying

```math
Q(A)\le M_n+\eta,
```

then every principal `m`-restriction `A[U]` obeys

```math
\boxed{
Q(A[U])-M_m
\le\eta+M_k+\sqrt{2(\log2)mkn}.}                   \tag{MH.2}
```

There is an absolute, fully explicit random-sign upper bound

```math
M_k\le\sqrt{(\log2)(k^3-k)}.                       \tag{MH.3}
```

Thus, uniformly for every `k=o(n)` and every `U` of size `n-k`,

```math
{Q(A[U])-M_{n-k}\over n^{3/2}}
\le {\eta\over n^{3/2}}
 +O\left(\sqrt{k\over n}+left({k\over n}\right)^{3/2}\right).
                                                               \tag{MH.4}
```

More explicitly, if `A_n` is a signing sequence with
`Q(A_n)-M_n=o(n^(3/2))`, if `k_n=o(n)`, and if `U_n` is **any** sequence of
`n-k_n` principal vertex sets, then

```math
Q(A_n[U_n])-M_{n-k_n}=o(n^{3/2}).
```

Thus every vanishing near-minimizer remains a vanishing near-minimizer after
**every** `o(n)`-vertex deletion.  The earlier direct edge-count estimate
only reached `o(sqrt n)` deletions.

### Corollary MH.2 (dense liminf orders suffice)

Suppose there is a strictly increasing sequence `n_j->infinity` such that

```math
{M_{n_j}\over n_j^{3/2}}
\longrightarrow\liminf_n {M_n\over n^{3/2}},
\qquad
n_{j+1}-n_j=o(n_j).                                 \tag{MH.5a}
```

Then `M_n/n^(3/2)` converges.  Indeed, for
`n_j<=n<n_(j+1)`, apply (MH.1) with
`m=n_j` and `k=n-n_j=o(n_j)`.  Equations (MH.1) and (MH.3) make the added
term `o(n_j^(3/2))`, while `n/n_j->1`.  The resulting limsup is at most
the displayed liminf.  This density criterion is an archive rediscovery:
see ledger Sections 1.8 and 10.18 and the all-order action-recovery audit,
Section 7.

## 2. Proof

Take exact minimizers `B` and `C` of orders `m` and `k`.  Let `R` be an
`m`-by-`k` matrix of independent Rademacher signs.  For a fixed Boolean pair
`(x,y)`, the bilinear form `x^TRy` is a sum of `mk` independent signs, so

```math
\Pr\{|x^TRy|>t\}\le2\exp\{-t^2/(2mk)\}.            \tag{MH.5}
```

There are only `2^(m-1)2^(k-1)=2^(n-2)` projective pairs.  At

```math
t=\sqrt{2(\log2)mkn},
```

the union-bound probability is at most `1/2`.  Hence one exact-sign bridge
has Boolean bilinear cap at most `t`.  The block signing

```math
P=\begin{pmatrix}B&R\\R^T&C\end{pmatrix}
```

then satisfies

```math
Q(P)\le M_m+M_k+t,
```

which proves (MH.1).

For any principal set `U`, unbiased independent spins on its complement
give

```math
H_{A[U]}(x_U)=\mathbb E[H_A(x_U,X_{U^c})].
```

Therefore `Q(A[U])<=Q(A)`.  Combine this with the near-minimality
hypothesis and (MH.1) to get (MH.2).

For `k=1`, (MH.3) is the identity `M_1=0`.  For `k>=2`, choose a uniformly
random hollow signing of order `k`.  For each
projective spin its energy is a sum of `E_k=binom(k,2)` independent signs.
A union bound over `2^(k-1)` spins and both tails is strictly below one at

```math
t_k=\sqrt{2E_k(k+1)\log2}
   =\sqrt{(\log2)(k^3-k)}.
```

Some signing therefore has cap at most `t_k`, proving (MH.3).  Substitution
in (MH.2) gives (MH.4).  `square`

## 3. What this changes

1. **Direct near-minimizer contact.**  This is uniform over every exact or
   certified near-minimizer and every principal set.  It is not a claim
   about sampled restrictions or a structured family.
2. **Near-order transfer only.**  Equation (MH.1) says
   `M_n-M_(n-k)=o(n^(3/2))` for `k=o(n)`.  It gives no control across a
   fixed order ratio and therefore cannot by itself force convergence.
   Corollary MH.2 is a useful density criterion, but it was already present
   in the archive and no theorem presently forces liminf orders to have
   relative gaps tending to zero.
3. **No active-set heredity.**  The restricted signing is near-minimal in
   cap; its ground states, shell geometry, and atlas gauges need not be
   inherited from those of `A`.  The order-11 exact example whose every
   one-vertex restriction has cap 17 rather than `M_10=13` remains fully
   consistent with (MH.2).
4. **Mesoscopic two-cap consequence.**  In the collapsed branch of PP.4,
   the exceptional shore has `|S|=o(n)`.  The large principal block
   `A[S^c]` is therefore a rigorously certified vanishing near-minimizer at
   its own order.  What is lost is precisely the active interface response:
   the two high-energy poles agree on `delta(S)` and become opposite after
   deleting it.  Thus deletion removes the leading response channel rather
   than transporting its shell.
5. **Archive comparison.**  The random rectangular bridge upper bound and
   dense-liminf-order criterion are already recorded in the ledger.  The
   only newly isolated item is the elementary synthesis with principal
   heredity, which extends the campaign's prior deterministic `o(sqrt n)`
   deletion window to all `o(n)` for near-minimizers (together with a modest
   finite-constant sharpening).  This does not change the fixed-ratio
   frontier and must not be credited as a new bridge method.

The smallest missing step is now explicit: to turn macroscopic principal
heredity into a reusable state, one must transport some response object
through deletion.  Scalar cap near-minimality alone discards the
mesoscopic interface where PP.4 can store the full leading energy.
