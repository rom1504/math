# Multiscale affine shells from a partition cap budget

Date: 2026-08-17.

Status: proof draft.  This strengthens the low-local-field affine cube from
dimension `Theta(sqrt n)` to dimension `floor(n/q)` with relative shell
defect `O(1/q)`, uniformly for every signing.  It is a one-shot response
theorem, not a physical all-endpoint composition theorem.

## 1. A cap budget across every vertex partition

For a hollow symmetric real matrix `D`, write

```math
Q_+(D)=\max_x H_D(x),
\qquad Q_-(D)=-\min_x H_D(x),
\qquad Q(D)=\max(Q_+(D),Q_-(D)).
```

### Lemma MP.1 (partition cap budget)

For every vertex partition `[n]=J_1 sqcup ... sqcup J_q`,

```math
\sum_(a=1)^q Q_+(D[J_a])\le Q(D),
\qquad
\sum_(a=1)^q Q_-(D[J_a])\le Q(D),                 \tag{MP.1}
```

and consequently

```math
\boxed{\sum_(a=1)^q Q(D[J_a])\le2Q(D).}           \tag{MP.2}
```

#### Proof

In each block choose a spin maximizing the positive internal energy.
Multiply those block spins by independent global block signs.  Internal
energies do not change and every cross-block edge has mean zero.  Hence the
mean full energy over the block signs is the sum of the positive block
maxima; some choice attains at least that mean, proving the first inequality.
Choose negative minimizers instead and take a choice no larger than its
mean to prove the second.  Finally
`max(u,v)<=u+v` for the nonnegative one-sided caps. `square`

### Lemma MP.2 (subset-edge mass costs only the one-sided block cap)

For every `S subseteq J` and every hollow symmetric `D[J]`,

```math
-\sum_(\{i,j\}\subseteq S)d_ij\le Q_-(D[J]),
\qquad
\left|\sum_(\{i,j\}\subseteq S)d_ij\right|
\le Q(D[J]).                                       \tag{MP.3}
```

#### Proof

Fix every spin in `S` to `+1` and choose the spins in `J\setminus S`
independently and uniformly.  Every edge not contained in `S` has mean
zero, so the expected block energy is exactly the displayed subset-edge
sum `P_S`.  Since every realized energy lies in
`[-Q_-(D[J]),Q_+(D[J])]`, one has
`-Q_-<=P_S<=Q_+`, proving both assertions. `square`

## 2. The multiscale cube

### Theorem MP.3 (universal multiscale affine shell algebra)

Let `A` be a hollow signing of order `n`, let `2<=q<=n`, and put
`Q=Q(A)`.  There are

* an absolute ground state `x` and orientation `rho in {+-1}` with
  `rho H_A(x)=Q`;
* a vertex set `I` of size at least `floor(n/q)`;

such that all projectively distinct spins

```math
\mathcal C_I=\{x^S:S\subseteq I\}
```

obey the stronger one-sided estimate

```math
\boxed{\rho H_A(x^S)\ge\left(1-{8\over q}\right)Q.} \tag{MP.4}
```

In particular, when `q>8` the whole cube has one common positive
orientation.  The family is closed under every odd coordinatewise product.

#### Proof

Orient and switch the signing:

```math
D=\rho\,\operatorname {diag}(x)A\operatorname {diag}(x).
```

Then `H_D(1)=Q`.  Its row sums

```math
\ell_i=\sum_jd_ij
```

satisfy `ell_i>=0` and `sum_iell_i=2Q`.  Partition the vertices into `q`
blocks of sizes `floor(n/q)` or `ceil(n/q)`.  Put

```math
L_a=\sum_(i\in J_a)\ell_i,
\qquad Q_a^-=Q_-(D[J_a]).
```

By Lemma MP.1,

```math
\sum_a(2L_a+4Q_a^-)\le4Q+4Q=8Q.
```

Choose a block `I=J_a` with `2L_a+4Q_a^-<=8Q/q`.  For every `S subseteq I`,
the exact flip identity is

```math
H_D(1^S)=Q-2\sum_(i\in S)\ell_i
             +4\sum_(\{i,j\}\subseteq S)d_ij.     \tag{MP.5}
```

Lemma MP.2 gives directly

```math
Q-H_D(1^S)
\le2L_a+4Q_a^-\le {8Q\over q}.
```

Switching back proves (MP.4), including the common orientation when `q>8`.
Since `q>=2`, the complement of `I` is
nonempty; distinct masks stay distinct modulo global spin reversal.  Odd
products XOR the masks, proving closure. `square`

### Corollary MP.4 (near-linear vanishing-width shell entropy)

If `Q(A_n)=O(n^(3/2))` and `q_n->infinity` with `q_n<=n/2`, then every
`A_n` has

```math
2^{\lfloor n/q_n\rfloor}
```

projectively distinct spins in an `O(n^(3/2)/q_n)` absolute shell.  For
`q_n>8` this entire cube has one common energy orientation.  For
example, `q_n=log^2n` gives affine dimension `Theta(n/log^2n)` and vanishing
normalized deficit `O(1/log^2n)`.

After choosing a constant-distance binary subcode of the masks, the shell
contains

```math
\exp(\Omega(n/q_n))
```

oriented cut words at pairwise edge-Hamming distance

```math
\Theta(n^2/q_n)                                             \tag{MP.6}
```

whenever `q_n->infinity`.  This is carrier-word packing with one short
affine grammar, not contextual information complexity.

Equivalently, for any shell budget `Delta` with
`16Q/n<=Delta<=2Q`, choose `q=ceil(8Q/Delta)`.  Then the absolute
`Delta`-shell contains an affine cube of dimension at least

```math
\left\lfloor {n\Delta\over16Q}\right\rfloor.       \tag{MP.6a}
```

For `Q=O(n^(3/2))`, this gives the universal entropy law
`log|S_Delta|=Omega(Delta/sqrt n)` throughout that range.  It strengthens
the earlier cardinality-only deep-hole shell estimate there, but because it
is universal it carries no selective near-minimality information.

## 3. Growing declared response language

Take an even-sized subset `I' subseteq I`, losing at most one coordinate,
and form the star port frame

```math
W=(x,(x^{\{i\}})_(i\in I')).
```

For a real field `g`, write

```math
\mathcal B_A(g)=
\max_(y\in\{+-1\}^n,\sigma\in\{+-1\})
\{\sigma H_A(y)+g\mathbin\cdot y\}.
```

The majority selector of every endpoint and every odd port product belongs
projectively to `mathcal C_I`.  Hence, for every real `m>=0`,

```math
0\le Q+m||W\epsilon||_1-\mathcal B_A(mW\epsilon)
\le {8Q\over q}.                                   \tag{MP.7}
```

More strongly, the **one-sided** response of the oriented child obeys the
same bound, with no endpoint-dependent quadratic sign:

```math
0\le Q+m||W\epsilon||_1
 -\max_y\{H_(\rho A)(y)+m(W\epsilon)\mathbin\cdot y\}
\le {8Q\over q}.                                  \tag{MP.8}
```

Indeed, for `epsilon in {+-1}^{|I'|+1}` put
`t=sum_jepsilon_j`.  The frame field has coordinate `x_it` off `I'` and
`x_i(t-2epsilon_i)` at its indexed exceptional row.  Both integers are odd,
so its sign selector differs projectively from `x` only on a subset of
`I'`.  It therefore has the common oriented energy supplied by MP.4 and
pays the full field norm.  Cap plus Holder gives the matching upper bounds
in (MP.7)--(MP.8).

The already-declared projective frame has a star-histogram presentation of
`O((n/q)log n)` bits.  Choosing `q/log n->infinity` makes this sublinear,
while `q->infinity` makes the response error `o(n^(3/2))` on every
bounded-cap sequence.

This is a strict sub-landscape response quotient on a designed, growing
language.  It is much stronger than scalar shell cardinality.  It still
does not encode the labelled ground-state gauge independently of `A`.

## 4. Physical and compositional ceiling

The theorem sharply separates one-shot response compression from reusable
physical composition.

1. The full exact-sign port block has an aligned endpoint of order
   `n|I|=Theta(n^2/q)`.  This is `o(n^(3/2))` only for `q>>sqrt n`, exactly
   when the interface is `o(sqrt n)`.  The larger affine dimension therefore
   cannot be physicalized by the raw frame at target scale.
2. If `|I'|>=2`, a `t=3` endpoint has field norm `3n-4` and selects the exact
   ground state (for the degenerate empty frame, the single ground port has
   norm `n`).  Scalar microcanonical compilation of this endpoint has residual
   `O(sqrt(n|I|(n+|I|)))=O(n^(3/2)/sqrt q)`, which is subleading for
   `q->infinity` but is the known random-bridge scale.
3. Traversing a fixed multiplicative change of order through steps of
   relative size `1/q` requires `Theta(q)` separately paid steps.  Their
   normalized `O(1/sqrt q)` compiler errors do not sum.  Recomputing the
   favorable block supplies no congruence between steps.

Thus MP.3 proves that very large vanishing-width affine shells and their
one-shot response algebra are **universal bounded-cap phenomena**.  They
cannot serve as near-minimizer-specific structure.  What remains missing is
joint cancellation of physical residual channels or a cross-level
congruence; neither follows from the partition budget.
