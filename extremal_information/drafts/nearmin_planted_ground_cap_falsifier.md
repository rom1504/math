# A planted-ground perturbation confines genuine near-minimizer shells

**Status:** proved and independently audited.  This is an
adversarial falsifier for structural statements about *all vanishing
near-minimizers*.  It does not construct, or obstruct, an exact minimizer with
the same geometry.

## 1. Result

Let

```math
E={n\choose2},\qquad
H_a(x)=\sum_{i<j}a_{ij}x_ix_j,
```

and write an oriented augmented cut as `z=sigma c(x)`, where
`c(x)_(ij)=x_i x_j`.  Thus

```math
Q(a)=\max_z\langle a,z\rangle .
```

For `Delta>=0`, define the positive oriented shell

```math
\mathcal S_a^+(\Delta)
=\{z:\langle a,z\rangle\ge Q(a)-\Delta\}.       \tag{PG.1}
```

The edge-Hamming distance is `d_E(z,z')=|{e:z_e ne z'_e}|`.  Its
projective version is `min{d_E(z,z'),E-d_E(z,z')}`.

### Theorem PG.1 (planted-ground cap)

Let `a` be a hollow signing with `Q=Q(a)=O(n^(3/2))`, and choose an
oriented ground word `z_0` with `<a,z_0>=Q`.  Put

```math
D=\{e:a_e(z_0)_e=-1\},\qquad N_D=|D|={E-Q\over2}.
```

For every integer `r` such that

```math
n=o(r),\qquad r=o(n^2),                              \tag{PG.2}
```

and, at each finite order, `1<=r<=N_D`,

there is a set `F subseteq D`, `|F|=r`, such that the exact signing `b`
obtained by flipping precisely the edges of `F` has

```math
Q(b)=Q(a)+2r                                             \tag{PG.3}
```

and every `z in \mathcal S_b^+(Delta)` satisfies

```math
{d_E(z,z_0)\over E}
\le {32n\over r}+{\Delta\over r}+{\Delta\over2E}.       \tag{PG.4}
```

The numerical constant `32` is inessential.  In particular, if

```math
\Delta=o(r),\qquad n=o(r),                              \tag{PG.5}
```

then the complete positive `Delta`-shell lies in one edge-Hamming, hence
one projective, cap of radius `o(1)` around `z_0`.

### Corollary PG.2 (Level-5 vanishing-near-minimizer falsifier)

Let `Delta_n=o(n^(3/2))` be any prescribed shell width.  For all sufficiently
large `n`,
start with an exact minimizer `a_n`, choose a ground word, and apply Theorem
PG.1 with, for example,

```math
h_n=\max\{n,\Delta_n\},\qquad
r_n=\left\lceil\sqrt{n^{3/2}h_n}\right\rceil .          \tag{PG.6}
```

Then

```math
Q(b_n)=M_n+2r_n=M_n+o(n^{3/2}),                         \tag{PG.7}
```

so `b_n` is a genuinely `epsilon_n`-near-minimizing sequence with
`epsilon_n=2r_n/n^(3/2)->0`, while

```math
\sup_{z\in\mathcal S_{b_n}^+(\Delta_n)}
{d_E(z,z_{0,n})\over E}\longrightarrow0.               \tag{PG.8}
```

Consequently the shell contains no fixed-scale two-point projective
packing for all sufficiently large `n`.  Its projective cardinality is at
most `exp(o(n))`: (PG.8) forces the underlying vertex cut to differ from
`z_0` on `o(n)` vertices, up to global vertex sign.

Indeed, if the right side of (PG.4) is `eta_n`, any two shell words satisfy

```math
{ |\langle z,z'\rangle|\over E}
\ge 1-4\eta_n,                                         \tag{PG.8a}
```

by the Hamming triangle inequality.  Thus the failure is directly in the
absolute-overlap metric used by the conditional physical compiler.

The same construction gives a fixed-halo statement.  For any fixed
`epsilon>0`, take `r=floor(epsilon n^(3/2)/4)`.  Then
`Q(b)<=M_n+epsilon n^(3/2)` and every shell of width
`Delta=o(n^(3/2))` is confined to one vanishing projective cap.

Equivalently, if `epsilon_n->0` but `epsilon_n sqrt(n)->infinity`, take
`r` proportional to `epsilon_n n^(3/2)`.  The parent lies in
`N_n(epsilon_n)`, and every shell satisfying
`Delta=o(epsilon_n n^(3/2))` has radius

```math
O\left({1\over\epsilon_n\sqrt n}
       +{\Delta\over\epsilon_n n^{3/2}}
       +{\Delta\over n^2}\right)=o(1).             \tag{PG.8b}
```

## 2. Exact gap identity

For an augmented cut `z`, let

```math
T_z=\{e:z_e ne (z_0)_e\},\qquad
t_z=|T_z|,
```

and put

```math
d_a(z)=Q-\langle a,z\rangle\ge0,
\qquad m_z=|D\cap T_z|.
```

Since `a_e(z_0)_e=-1` on `D`, flipping `F subseteq D` gives

```math
\langle b,z\rangle
=\langle a,z\rangle+2\sum_{e\in F}(z_0)_ez_e.
```

The target word gains `2r`.  Subtraction gives the exact identity

```math
\boxed{
Q+2r-\langle b,z\rangle
=d_a(z)+4|F\cap T_z|.}                                \tag{PG.9}
```

Both terms on the right are nonnegative.  Therefore `z_0` is a ground word
of `b` and (PG.3) follows with equality; no upper-bound relaxation or
unknown assertion about the other cuts is being used.

There is a second exact identity.  Writing `s_e=a_e(z_0)_e`,

```math
d_a(z)=2\sum_{e\in T_z}s_e
      =2(t_z-2m_z),
```

and hence

```math
t_z=2m_z+{d_a(z)\over2}.                              \tag{PG.10}
```

Thus a near-top competitor far from `z_0` must disagree with `z_0` on many
of the initially negative ground edges `D`.  A sparse uniform hitting set
inside `D` detects all such competitors simultaneously.

## 3. Uniform hitting lemma

### Lemma PG.3

Let `mathcal G` be at most `2^n` subsets of a set `D` of size `N_D`, and
let `1<=r<=N_D`.  There is an `r`-element set `F subseteq D` such that for
every `G in mathcal G`, either

```math
|G|<{16nN_D\over r},                                  \tag{PG.11}
```

or

```math
|F\cap G|\ge {r|G|\over2N_D}.                         \tag{PG.12}
```

#### Proof

Choose `F` uniformly among the `r`-subsets of `D`.  For fixed `G`, the
hypergeometric variable `X=|F cap G|` has mean
`mu=r|G|/N_D`.  The standard lower-tail bound gives

```math
Pr(X<mu/2)\le e^{-mu/8}.
```

If (PG.11) fails, this is at most `e^(-2n)`.  A union bound over at most
`2^n` sets is strictly below one.  Therefore a simultaneous choice exists.
`square`

Apply the lemma to

```math
\mathcal G=\{D\cap T_z:z\text{ an augmented cut}\}.
```

There are at most `2^n` distinct oriented augmented-cut words.  If
`z in mathcal S_b^+(Delta)`, (PG.9) gives

```math
d_a(z)\le\Delta,\qquad |F\cap T_z|\le\Delta/4.         \tag{PG.13}
```

If (PG.11) holds then `m_z<16nN_D/r`.  Otherwise (PG.12) and (PG.13) give
`m_z<=Delta N_D/(2r)`.  In either case, and harmlessly summing the two
bounds,

```math
m_z\le {16nN_D\over r}+{\Delta N_D\over2r}.           \tag{PG.14}
```

Combining (PG.10), (PG.13), (PG.14), and `N_D<=E` proves (PG.4).
For a bounded-cap base, `N_D=(E-Q)/2=Theta(n^2)`, so (PG.2) also ensures
`r<=N_D` for large `n`.

The proof is actually code-generic.  For a sign code `\mathcal Z subseteq
{+-1}^E`, replace `2^n` in Lemma PG.3 by `|mathcal Z|` and the threshold
`16n` by a constant times `log|mathcal Z|`.  Thus the mechanism is not a
hidden quadratic identity: what makes it effective here is the enormous
ambient edge length `Theta(n^2)` together with only `2^n` augmented cuts.

## 4. Projective cardinality

For completeness, an oriented augmented cut within `eta E` edge-Hamming
distance of `z_0` has, after fixing the common orientation, vertex-Hamming
distance `k` satisfying

```math
k(n-k)\le\eta E.                                     \tag{PG.15}
```

For `eta=o(1)`, the edge orientation must agree with that of `z_0`:
an oppositely oriented augmented cut differs on at least
`E-floor(n^2/4)=(1/2-o(1))E` edges.  Projectivizing the vertex word then
lets us take `k<=n/2`, whence `k=o(n)`.  Therefore the number of possible projective
words is at most

```math
\sum_{k\le o(n)}{n\choose k}=\exp(o(n)).              \tag{PG.16}
```

This is much stronger than failure of the fixed-`gamma` premise in the
absolute-overlap compiler: the entire positive shell has subexponential
projective support.

## 5. Assumption-distance and quantifier audit

1. **Actual quadratic energies.**  The construction is an exact hollow
   `+-1` signing.  Equations (PG.9)--(PG.10) use the quadratic augmented-cut
   code, not an abstract two-word marginal law.
2. **Certified near-minimality.**  Starting from any exact minimizer gives
   the equality `Q(b)=M_n+2r`; this is a mathematical Level-5 certificate,
   not a heuristic upper bound.  It is existential rather than an explicit
   efficiently computable family because the starting minimizer is not
   known in closed form.
3. **What is falsified.**  Any claim that *every* fixed-`epsilon` or
   vanishing-`epsilon_n` near-minimizer has a fixed-scale projective packing
   in its exact active shell is false.  More strongly, for every prescribed
   `Delta_n=o(n^(3/2))`, some vanishing-near-minimizer sequence has its
   complete positive `Delta_n`-shell in one vanishing cap.
4. **What survives.**  The theorem does **not** falsify a statement
   restricted to exact minimizers.  It also does not falsify a halo theorem
   whose admitted shell width is required to dominate the parent's
   near-minimality slack: our confinement needs `Delta=o(r)`, whereas the
   excess cap is `2r`.
   In particular it does not contradict Theorem 36.2, whose balanced shell
   for an `epsilon`-near-minimizer is deliberately chosen at normalized
   width `kappa` with `kappa` larger than the normalized excess.
5. **Exact versus one-step near.**  The perturbation is a rigorous scalable
   version of the finite phenomenon that one-step-near active shells can be
   spiky even when all audited exact-minimizer active shells are diffuse.
   Exact optimality is therefore not a cosmetic strengthening: it is the
   entire remaining possible source of projective-shell rigidity.
6. **Bounded cap.**  The construction has `Q(b)=O(n^(3/2))`.  Thus bounded
   cap alone, spectral-scale normalization, and vanishing relative excess
   do not prevent planted projective concentration.
7. **Compatibility with universal affine shells.**  Universal
   multiscale-partition shells may still exist inside the cap in (PG.8):
   flipping `o(n)` vertices creates exponentially many words while moving
   only `o(E)` edge coordinates.  Large thin-shell cardinality and
   projective diffusion are genuinely different requirements.
8. **Archive collision and new increment.**  The equality (PG.3) and the
   one-sided frozen-edge inequality implicit in (PG.9) are the existing
   geodesic planted-face mechanism from
   `nearmin_deterministic_inequalities.md`.  The new step here is Lemma PG.3
   applied to the *entire augmented-cut family*, which converts frozen
   coordinates into the uniform geometric conclusion (PG.4), and hence the
   subexponential projective-shell falsifier (PG.8).  The result should not
   be credited as a new planting mechanism.

## 6. Research judgment

The projective-shell route now has a sharp boundary:

- for arbitrary vanishing near-minimizers it is rigorously false, even at
  any preassigned subleading shell width;
- for exact minimizers it survives all current finite tests and is not
  touched by this planted perturbation;
- any proof for exact minimizers must use a non-perturbatively stable
  consequence of *exact* global optimality.  A property continuous under
  `o(n^(3/2))` coefficient edits cannot supply it.

This makes the exact-minimizer three-point statement a cleaner target, but
also increases its assumption distance: it sits precisely at a structural
discontinuity that disappears under arbitrarily small normalized cap slack.
