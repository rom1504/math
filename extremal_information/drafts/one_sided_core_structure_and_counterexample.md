# One-sided principal cores: exact cut structure and a diffuse counterexample

**Status.** Task-local theorem draft.  All deterministic identities and the
probabilistic construction below are proved here.  This note does **not**
construct such a core inside an exact minimizer; it determines what can and
cannot follow from the core's internal one-sided geometry alone.

## 1. Normalization

Let `B` be a hollow symmetric order-`k` sign matrix, where `k>=2`, and put

```math
H_B(x)=\sum_{i<j}b_{ij}x_ix_j={1\over2}x^TBx,
\qquad
P=P(B)=\max_xH_B(x),
\qquad
N=N(B)=-\min_xH_B(x).                         \tag{OS.1}
```

Write `K=binom(k,2)`.  Switch by a positive maximizer, so throughout
Sections 2--3

```math
H_B(\mathbf 1)=P.                               \tag{OS.2}
```

For `S subseteq [k]`, let

```math
w(S)=\sum_{i\in S,\ j\notin S}b_{ij}.           \tag{OS.3}
```

## 2. The exact structural atlas

### Proposition OS.1 (cut window, rows, and switching distance)

Under (OS.2),

```math
\boxed{0\le w(S)\le {P+N\over2}\quad(S\subseteq[k]),}
                                                               \tag{OS.4}
```

and, for a uniform random subset `S`,

```math
\mathbb E w(S)={P\over2},
\qquad
\max_Sw(S)-\mathbb Ew(S)={N\over2}.             \tag{OS.5}
```

The signed row sums `r_i=sum_(j ne i)b_(ij)` consequently satisfy

```math
0\le r_i\le {P+N\over2},
\qquad
\sum_i r_i=2P.                                  \tag{OS.6}
```

Moreover, the minimum Hamming distance of `B` from a switched all-positive
signing is exactly

```math
\boxed{
\min_{u\in\{\pm1\}^k}
\#\{ij:b_{ij}\ne u_i u_j\}={K-P\over2}.}       \tag{OS.7}
```

#### Proof

If `x^S` is obtained from `1` by flipping `S`, then

```math
H_B(x^S)=P-2w(S).                                \tag{OS.8}
```

Maximality at `1` and the lower endpoint `-N` give (OS.4).  Each edge
crosses a uniform cut with probability `1/2`, proving (OS.5).  Singleton
cuts give (OS.6).  Finally

```math
H_B(u)=K-2\#\{ij:b_{ij}\ne u_i u_j\},           \tag{OS.9}
```

and maximizing the left side proves (OS.7). `square`

Thus a near-clique conclusion is available **if and only if** the separate
edge-count statement `P=K-o(k^2)` is available.  The asymmetry `N=o(P)`
does not itself say this.

In graph language, the graph of positive entries has at least half of the
edges of every cut and is maximum-sized in its Seidel-switching class.  The
extra datum `N=o(P)` says that its largest signed cut exceeds the mean
signed cut by only `o(P)`.

### Proposition OS.2 (variance and triangle constraints)

For a uniform Boolean spin `X`,

```math
\mathbb EH_B(X)=0,
\qquad
\mathbb EH_B(X)^2=K.                             \tag{OS.10}
```

Consequently

```math
\boxed{PN\ge K.}                                \tag{OS.11}
```

If

```math
\tau(B)=\sum_{i<j<ell}b_{ij}b_{jell}b_{ell i},
```

then

```math
\mathbb EH_B(X)^3=6\tau(B)                      \tag{OS.12}
```

and

```math
\boxed{
6\tau(B)\ge {K^2\over N}-NK.}                  \tag{OS.13}
```

#### Proof

Walsh orthogonality gives (OS.10).  The pointwise inequality

```math
(P-H_B)(H_B+N)\ge0
```

gives `H_B^2 <= (P-N)H_B+PN`; averaging proves (OS.11).  In the cubic
moment only ordered triples forming a triangle survive, which proves
(OS.12).  Finally, for every real `a`,

```math
(H_B+N)(H_B-a)^2\ge0.
```

Averaging and using (OS.10) yields

```math
6\tau(B)\ge (2a-N)K-Na^2.
```

Choosing `a=K/N` proves (OS.13). `square`

The triangle constraint becomes genuinely positive only in the much
stronger range `N^2<K`.  It gives no clique stability at the ambient-scale
core parameters relevant below.

### Proposition OS.3 (one-sidedness forces a spectral outlier, but not
rank-one approximation)

Let `B_k` be any sequence with

```math
{N(B_k)\over P(B_k)}\longrightarrow0.            \tag{OS.14}
```

Then

```math
{P(B_k)\over k^{3/2}}\longrightarrow\infty,
\qquad
{\|B_k\|_{2\to2}\over\sqrt k}\longrightarrow\infty.  \tag{OS.15}
```

#### Proof

The archived Bollobas--Scott one-sided discrepancy-product theorem, in the
half-energy normalization (OS.1), gives

```math
N(P+N)\ge
{\bigl(1-(N/K)^2\bigr)k^3\over6400}              \tag{OS.16}
```

whenever `(1-(N/K)^2)/4 >= 1/k`.  Under (OS.14), `N/K -> 0`: otherwise
`P>N` and `P<=K` would make `N/P` bounded away from zero along a
subsequence.  Hence (OS.16) applies with right side `(1-o(1))k^3/6400`.
Because `N=o(P)`, this forces `P^2/k^3 -> infinity`.  Finally

```math
\|B\|_{2\to2}\ge {\mathbf1^TB\mathbf1\over k}={2P\over k},
```

which proves the second assertion. `square`

This is the strongest general spectral consequence found: the gauged
ground state has a super-bulk Rayleigh quotient.  It does not say that the
matrix is close in entries or Frobenius norm to that rank-one direction.

## 3. A scalable diffuse one-sided core

### Theorem OS.4 (cut-positive biased-random counterexample)

Fix `0<alpha<1/2` and put `mu_k=k^(-alpha)`.  For all sufficiently large
`k` there is a hollow symmetric sign matrix `B_k` such that

```math
\begin{aligned}
&H_(B_k)(\mathbf1)=P(B_k)=(1+o(1))\mu_k K,\\
&N(B_k)=O(k^{3/2}),\\
&{N(B_k)\over P(B_k)}=O(k^{alpha-1/2})=o(1),      \tag{OS.17}\\
&\min_u\#\{ij:(B_k)_{ij}\ne u_i u_j\}
  =\left({1\over2}-o(1)\right)K.                 \tag{OS.18}
\end{aligned}
```

Every signed cut is nonnegative.  In fact the rows can simultaneously be
chosen uniformly regular:

```math
r_i=(1+o(1))\mu_k k\quad\hbox{for every }i.       \tag{OS.19}
```

Thus the one-sidedness is not carried by a near-clique, a switching-rank-one
matrix, or a small exceptional set of heavy rows.

#### Proof

Choose the upper-triangular entries independently with

```math
\Pr(b_{ij}=1)={1+\mu_k\over2},
\qquad
\Pr(b_{ij}=-1)={1-\mu_k\over2}.                  \tag{OS.20}
```

For a fixed cut with `m=s(k-s)` crossing edges, Hoeffding's inequality
gives

```math
\Pr\left(\sum_{i\in S,j\notin S}b_{ij}\le0\right)
\le \exp(-\mu_k^2m/2).                           \tag{OS.21}
```

For `s<=k/2`, `m>=sk/2`, and therefore the union bound is at most

```math
\sum_{s=1}^{k/2}
\exp\left[s\log(ek/s)-{1\over4}\mu_k^2sk\right]=o(1),  \tag{OS.22}
```

because `mu_k^2k=k^(1-2alpha) >> log k`.  Hence all cuts are positive with
probability tending to one.

Write

```math
B=\mu_k(J-I)+W,                                  \tag{OS.23}
```

where the upper-triangular entries of `W` are independent, centered, and
uniformly bounded.  A standard epsilon-net proof gives, with probability
`1-o(1)`,

```math
\|W\|_{2\to2}\le C\sqrt k.                      \tag{OS.24}
```

For completeness, for each fixed unit vector `z`, Hoeffding applied to
`z^TWz=2sum_(i<j)W_(ij)z_i z_j` gives
`Pr(|z^TWz|>t)<=2exp(-c t^2)`.  A `1/4`-net of the unit sphere has at most
`9^k` points, and the usual net comparison proves (OS.24) after taking
`t=C sqrt(k)`.

On the intersection of (OS.22) and (OS.24), Proposition OS.1 makes `1` a
global maximizer, while

```math
P=\mu_kK+{1\over2}\mathbf1^TW\mathbf1
  =\mu_kK+O(k^{3/2})
  =(1+o(1))\mu_kK.                               \tag{OS.25}
```

For every Boolean `x`,

```math
H_B(x)
={\mu_k\over2}\left(\left(\sum_i x_i\right)^2-k\right)
 +{1\over2}x^TWx
\ge-{\mu_k k\over2}-{Ck^{3/2}\over2}.           \tag{OS.26}
```

This proves the first three lines of (OS.17).  Equation (OS.18) follows
exactly from (OS.7) and `P=o(K)`.  Finally, a rowwise Hoeffding union bound
gives

```math
\max_i|r_i-\mu_k(k-1)|=O(\sqrt{k\log k})=o(\mu_k k),
```

which proves (OS.19).  The three required high-probability events have a
nonempty intersection for large `k`, establishing existence. `square`

The same construction has a clean spectral interpretation.  The rank-one
mean has top eigenvalue `mu_k(k-1)`, while (OS.24) is `O(sqrt(k))`; hence
there is one coherent spectral outlier and its leading eigenvector is
`o(1)`-close to `1`.  Nevertheless the centered residual contains almost
all Frobenius mass and (OS.18) says that half of all entries still disagree
with every switched rank-one signing.  A spectral spike and a near-clique
are therefore very different conclusions.

### Corollary OS.5 (every polynomial ambient core scale above the edge floor)

Fix any `3/4<gamma<1`, put

```math
\alpha=2-{3\over2\gamma}\in(0,1/2),
\qquad k=\lfloor n^\gamma\rfloor,
```

and use OS.4 at order `k`.  Then

```math
P(B_k)=\left({1\over2}+o(1)\right)n^{3/2},
\qquad
N(B_k)=O(n^{3\gamma/2})=o(n^{3/2}),              \tag{OS.27}
```

while `k=o(n)` and the switching distance remains
`(1/2-o(1))binom(k,2)`.  Thus diffuse counterexamples occur at every
polynomial zero-density core size strictly above the information-theoretic
edge floor `n^(3/4)`.

For example, take `alpha=1/4` and an ambient order `n_k=k^(7/6)` (rounding
to an integer does not matter).  Then

```math
k=n_k^{6/7}=o(n_k),
\qquad
P(B_k)=\left({1\over2}+o(1)\right)n_k^{3/2},
\qquad
N(B_k)=O(n_k^{9/7})=o(n_k^{3/2}),                \tag{OS.28}
```

while remaining at Hamming distance `(1/2-o(1))binom(k,2)` from every
switched clique.  Replacing `mu_k` by a fixed constant multiple tunes the
positive coefficient in (OS.27)--(OS.28).

Thus the orientation-separated core produced by the current rare-tail
reduction can have exactly the required ambient energy and zero density,
yet be internally diffuse.  Any theorem excluding it for exact minimizers
must exploit the core--complement coupling or coefficient minimality, not
just signed-cut positivity, row regularity, the displayed low-moment
constraints, existence of a spectral outlier, or internal edge density.

## 4. Archive comparison

The following ingredients are archived and are not claimed as new:

1. `artifacts/nonnegative_quadratic_fourier_cone.md` identifies (OS.4) with
   an s-maximal graph and notes that fixed graphon/cycle data are blind to
   second-order cut fluctuations.
2. `artifacts/one_sided_energy_product.md` proves the
   Bollobas--Scott inequality used in OS.3, with all normalization factors.
3. `extremal_information/drafts/exact_minimizer_oriented_core_separation.md`
   already proves (OS.4)--(OS.6) for a tail-generated exact-minimizer core.

The new result is OS.4--OS.5: a scalable, uniformly row-regular,
cut-positive family with `N/P -> 0` that is maximally far from every
switching-rank-one signing, including at the precise zero-density ambient
core scale.  The elementary cubic inequality (OS.13) and the general
spectral-outlier consequence (OS.15) sharpen the structural atlas, but do
not exclude this family.

## 5. Research classification

```text
PROVES:
  one-sided core => exact cut window + moment constraints;
  N/P -> 0 => a super-sqrt(k) spectral outlier;
  near-clique iff the separate edge-count condition P=K-o(k^2) holds.

FALSIFIES:
  one-sided core => near clique / entrywise rank one;
  one-sided core + regular rows => near clique;
  existence of a spectral outlier => small switching distance;
  generic internal structure alone excludes an ambient-scale diffuse core.

DOES NOT FALSIFY:
  an exact-minimality theorem coupling the core to its opposite-orientation
  near-minimal complement.
```
