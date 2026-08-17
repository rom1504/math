# Independent audit: one-sided core structure and diffuse counterexample

**Verdict on the repaired frozen source: PASS.**  The cut identities,
moment inequalities, Bollobas--Scott normalization, biased-random
construction, simultaneous row regularity, switching-distance conclusion,
and every ambient exponent in OS.5 are correct.  The construction is an
internal-geometry falsifier, not an exact-minimizer construction.  Two scope
repairs identified during the audit have been made: the harmless domain
condition `k>=2`, and replacement of the overbroad words "moments, spectrum"
by the precise displayed low-moment and spectral-outlier information.

## 1. Frozen source and normalization

```text
extremal_information/drafts/one_sided_core_structure_and_counterexample.md
sha256 82c10f2c62759a1c4293c1bc407563483d370d9d674441835f6d0f67c9c7397a
```

The source uses the half-quadratic normalization

```math
H_B(x)=\sum_{i<j}b_{ij}x_ix_j={1\over2}x^{\mathsf T}Bx,
\qquad K=\binom{k}{2}.
```

The originally submitted hash was
`bcffa83eacc30ee9e729a64bdef75e4f7fcbe0924c106ec0a59a8999cb6e3a27`.
It was changed only by the two scope repairs stated above.  The proof audit
below applies to the repaired frozen hash.

## 2. OS.1: cut window, rows, and switching distance

After switching a positive maximizer to `1`, flipping the vertices in `S`
changes the sign of exactly the crossing edges, hence

```math
H_B(\mathbf1^S)=P-2w(S).
```

The two endpoint inequalities `-N<=H_B(1^S)<=P` give exactly

```math
0\le w(S)\le {P+N\over2}.
```

Every Boolean spin is represented by one such cut (up to global sign), so
the upper endpoint is attained.  A uniform cut crosses each edge with
probability `1/2`; because `sum_(i<j)b_(ij)=P`,

```math
\mathbb E w(S)={P\over2},
\qquad
\max_Sw(S)-\mathbb Ew(S)={N\over2}.
```

For a singleton, `w({i})=r_i`, and summing the rows counts each edge twice.
Thus all factors in OS.4--OS.6 are correct.

For every `u`, an edge agrees with the switched clique `u_i u_j` precisely
when its contribution `b_(ij)u_i u_j` is `+1`.  Therefore

```math
H_B(u)=K-2d_H(B,uu^{\mathsf T}),
```

and maximizing proves OS.7 exactly.  In particular, distance `o(k^2)` from
some switched clique is equivalent to the separate condition
`P=K-o(k^2)`; no one-sided ratio assumption supplies it.

## 3. OS.2: variance and cubic moment

Distinct degree-two Walsh characters are orthogonal.  Hence

```math
\mathbb EH_B=0,
\qquad
\mathbb EH_B^2=K.
```

The pointwise endpoint polynomial

```math
(P-H_B)(H_B+N)\ge0
```

averages to `PN>=K`, with no missing factor two.  In the third moment, the
only surviving triples of edge characters are the six orderings of a base
triangle, so `E H_B^3=6 tau(B)`.

Expanding the second endpoint polynomial gives

```math
0\le\mathbb E[(H_B+N)(H_B-a)^2]
=6\tau(B)+(N-2a)K+Na^2.
```

Optimizing at `a=K/N` gives

```math
6\tau(B)\ge {K^2\over N}-NK.
```

Here `N>0` for `k>=2`, because the nonzero mean-zero quadratic polynomial
cannot be nonnegative everywhere.  The right side is positive exactly when
`N^2<K`, as stated.

As an independent finite check, I exhaustively evaluated all Boolean spins
for 100 random signings at each order `2<=k<=7`.  OS.4--OS.7,
`E H^3=6tau`, and `PN>=K` all passed exactly over the integers.  This is a
check, not an ingredient of the proof.

## 4. OS.3: one-sided discrepancy-product normalization

The archived theorem uses doubled energies

```math
P_d=2P,
\qquad N_d=2N,
\qquad R_d=2(P+N).
```

Its negative-orientation inequality is

```math
N_dR_d\ge
{(1-r_N^2)k^3\over1600},
\qquad
r_N={N_d\over k(k-1)}={N\over K}.
```

Dividing by four gives precisely OS.16:

```math
N(P+N)\ge
{(1-(N/K)^2)k^3\over6400}.
```

The applicability condition also translates exactly to
`(1-(N/K)^2)/4>=1/k`.

If `N/P->0`, then `N/K<=N/P->0`, since `P<=K`.  Thus the theorem applies
eventually and

```math
{N\over P}\left(1+{N\over P}\right){P^2\over k^3}
\ge {1-o(1)\over6400}.
```

This forces `P^2/k^3->infinity`.  The Rayleigh quotient of `1` is `2P/k`,
so `||B||_(2 to 2)/sqrt(k)->infinity`.  OS.15 therefore follows with the
correct normalization and direction.

## 5. OS.4: biased-random construction

Let `mu=k^(-alpha)` with `0<alpha<1/2`.  For a fixed cut of size `s<=k/2`,
the number of crossing edges satisfies `m=s(k-s)>=sk/2`.  Hoeffding gives

```math
\Pr\{w(S)\le0\}\le e^{-\mu^2m/2}
\le e^{-\mu^2sk/4}.
```

Using `binom(k,s)<=(ek/s)^s`, the total failure probability is bounded by

```math
\sum_{s=1}^{k/2}
\exp\left(s\left[\log(ek/s)-{1\over4}k^{1-2\alpha}\right]\right)=o(1).
```

Indeed, the negative term dominates `log(ek/s)` uniformly in `s`.  Thus all
nontrivial cuts are strictly positive with high probability; the empty and
full cuts have weight zero.  The cut flip identity then proves directly
that `1` is a global maximizer.

For `B=mu(J-I)+W`, the centered upper-triangular entries of `W` are
independent and have range length two.  For a fixed unit vector `z`, the
summands in

```math
z^{\mathsf T}Wz=2\sum_{i<j}W_{ij}z_i z_j
```

have squared ranges summing to at most
`16 sum_(i<j)z_i^2z_j^2<=8`.  Hence

```math
\Pr\{|z^{\mathsf T}Wz|>t\}\le2e^{-t^2/4}.
```

A `1/4`-net of size at most `9^k`, followed by the standard symmetric
quadratic-form net comparison, gives `||W||_(2 to 2)<=C sqrt(k)` with
probability `1-o(1)` after increasing the absolute constant `C`.

On this event,

```math
P=H_B(\mathbf1)=\mu K+O(k^{3/2})
=(1+o(1))\mu K,
```

because `mu k^2 >> k^(3/2)`.  For every Boolean `x`,

```math
H_B(x)
={\mu\over2}\bigl((\sum_i x_i)^2-k\bigr)
 +{1\over2}x^{\mathsf T}Wx
\ge-\frac{\mu k}{2}-\frac C2k^{3/2}.
```

Thus `N=O(k^(3/2))` and
`N/P=O(k^(alpha-1/2))=o(1)`.

OS.7 now yields

```math
\min_u d_H(B,uu^{\mathsf T})={K-P\over2}
=\left({1\over2}-o(1)\right)K,
```

which is asymptotically the largest possible nearest-codeword distance,
since averaging over `u` gives mean distance `K/2`.

Finally, each row has mean `mu(k-1)`.  Hoeffding plus a union bound over the
`k` rows gives

```math
\max_i|r_i-\mu(k-1)|=O(\sqrt{k\log k})
=o(\mu k),
```

using `alpha<1/2`.  The cut, operator, and row events each have probability
`1-o(1)`, so their intersection is nonempty for every sufficiently large
order.  No independence between these three events is needed.

The spectral interpretation is also correct.  The mean matrix has a top
eigenvalue `mu(k-1)` and an orthogonal eigenvalue `-mu`; its eigengap is
`mu k >> sqrt(k)`.  Davis--Kahan therefore puts the leading eigenvector of
`B` at angle `o(1)` from `1`.  Meanwhile `||W||_F^2=(1-o(1))k^2`, so the
spectral spike does not imply entrywise rank-one structure.

## 6. OS.5: ambient scaling

For `k=floor(n^gamma)` and

```math
\alpha=2-{3\over2\gamma},
\qquad 3/4<\gamma<1,
```

one has `0<alpha<1/2` and

```math
\gamma(2-\alpha)={3\over2}.
```

Consequently

```math
P(B_k)=(1/2+o(1))k^{2-\alpha}
       =(1/2+o(1))n^{3/2},
```

while `N(B_k)=O(k^(3/2))=O(n^(3gamma/2))=o(n^(3/2))` and `k=o(n)`.
The edge-count floor is exact: a `k`-vertex sign block has
`P<=binom(k,2)`, so carrying a fixed positive multiple of `n^(3/2)`
requires `k=Omega(n^(3/4))`.

For the displayed example, `alpha=1/4` gives `gamma=6/7`; hence

```math
k=n_k^{6/7},
\qquad
k^{2-\alpha}=k^{7/4}=n_k^{3/2},
\qquad
k^{3/2}=n_k^{9/7}.
```

All exponents and the coefficient `1/2` in OS.27--OS.28 are correct.  A
fixed positive multiple of `mu` tunes that coefficient without affecting
the probability estimates for sufficiently large `k`.

## 7. Archive collision and scope

The archive already contains:

- the equivalence between nonnegative signed cuts and an s-maximal graph;
- the exact switched-cut/row identities;
- the Bollobas--Scott one-sided energy-product theorem;
- the exact-minimizer orientation-separation theorem that motivates the
  internal core conditions;
- rank-one switching-minimal examples and warnings that one-sided
  switching minimality alone is weak.

I found no archived result constructing the OS.4 family with all four
simultaneous properties

```text
cut positivity,
N/P -> 0,
uniform row regularity,
nearest switched-clique distance (1/2-o(1))K,
```

nor its tuning to every ambient polynomial core scale
`n^gamma`, `3/4<gamma<1`.  OS.4--OS.5 are therefore a new scalable
falsifier relative to this repository.  The probabilistic mechanism itself
is elementary and should not be presented as a new random-graph technique.

The scope is important.  The family is not embedded in an exact minimizer,
does not constrain the core--complement cross block, and does not falsify
the exact-minimizer balance lemma.  It shows only that the core's internal
cut positivity, row regularity, displayed low moments, spectral outlier,
and edge density cannot by themselves force a near clique.  The repaired
source now states precisely this limited conclusion.

## 8. Disposition

```text
OS.1: PASS.
OS.2: PASS.
OS.3: PASS; doubled-to-half normalization independently checked.
OS.4: PASS; all three high-probability events coexist.
OS.5: PASS; all ambient exponents and the edge floor are correct.
Archive novelty: PASS as a new synthesis/scalable internal falsifier.
Exact-minimizer implication: deliberately not claimed.
```
