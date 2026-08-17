# Independent audit: equivariant matched-roof query bank

**Final verdict: PASS.**

The simultaneous switching refinement is correct: one realization at the
all-positive pole produces every biased fill by diagonal conjugation, and
all estimates used by BR.2 survive verbatim.  The resulting bounded-cap
physical response complexity is indeed `Theta(n)` bits on every fixed
bounded-cap switching orbit (conditional on the public base child and query
architecture).

The first inspected version required three scope/proof repairs:

1. If `P(A)<N(A)`, the proof negates `A` and thereby constructs contexts for
   the orbit of `-A`, not literally the stated orbit of `A`.  Pull the
   construction back by globally negating each complete parent: the bridge
   becomes `-B` and the fill bank becomes `-D_q=S_q(-D_1)S_q`, while every
   absolute cap is unchanged.
2. Prove that the full fill bank has exactly `2^(n-1)` members.  This turns
   “specified using `n-1` bits” into an exact conditional fixed-length
   query-index statement.  For the selected subbank, a public enumeration
   uses `ceil(log_2 |I|)` bits; alternatively declare the full projective
   pole language and avoid charging a codebook per query.
3. Replace the unsupported phrase that the public base “can require
   `Theta(n^2)` bits.”  The proof gives a one-time explicit storage upper
   bound of `binom(n,2)` bits and proves no description-complexity lower
   bound for the chosen realization.  State only that no subquadratic or
   explicit construction is supplied.

The source initially inspected was

```text
extremal_information/drafts/equivariant_matched_roof_query_bank.md
sha256 8e1e74871e8da36a772160b386d528613c5af075ee91b5b1ea42431caa28dc91
```

All three repairs were implemented.  The final source audited and frozen
for this verdict is

```text
extremal_information/drafts/equivariant_matched_roof_query_bank.md
sha256 6d4d8bf8d06917cd976fc8a0f3c0ad5f44eab4819fbf6707eb44ff3ac6a2e44b
```

## 1. Simultaneous fill algebra: PASS

Let `S_q=diag(q)`, `R_1=J-I`, and suppose the single frozen realization
satisfies

```math
D_1=pR_1+E_1,
\qquad \|E_1\|_{2\to2}\le K_0\sqrt n.
```

Then

```math
S_qR_1S_q=qq^T-I=R_q,
```

so

```math
D_q=S_qD_1S_q=pR_q+S_qE_1S_q.
```

Orthogonal conjugation preserves the error norm.  Diagonal switching also
preserves exact signs, hollowness, the complete Boolean energy multiset, and
therefore `Q(D_q)=Q(D_1)`.  One successful all-positive realization thus
works simultaneously for all projective poles; there is no exponential
union bound and no independently chosen error matrix.

In fact the full fill orbit is free modulo the global sign of the pole.  If
`D_q=D_r`, put `s=q\odot r`.  Since every off-diagonal entry of `D_1` is
nonzero,

```math
s_is_j=1\qquad(i\ne j).
```

Hence all coordinates of `s` agree and `q=+-r`.  Conversely `D_q=D_{-q}`.
Thus, for `n>=2`, the full bank has exactly `2^(n-1)` fills.  Conditional on
one shared `D_1`, its exact fixed-length pole index is `n-1` bits.

## 2. Every BR.2 estimate survives: PASS

The simultaneous construction preserves more than BR.2 needs.

For a boundary word `y`, gauge the pole `q` to `1` and write
`S={i:q_iy_i=-1}`, `d=|S|`, and `h=d(n-d)`.  Only crossing entries remain in
the error-energy difference, so

```math
\begin{aligned}
H_{D_q}(q)-H_{D_q}(y)
&=2ph+H_{E_q}(q)-H_{E_q}(y)\\
&\ge 2ph-2K_0\sqrt{nh}.
\end{aligned}
```

This is exactly BR.8.  Also

```math
H_{R_q}(q)+H_{R_q}(y)=n(n-1)-2h
```

and each error energy has modulus at most `K_0n^(3/2)/2`, giving BR.9.
Consequently:

- BR.15 uses the unchanged target equality and `H_{D_q}(q)`;
- BR.16 uses BR.8 and the unchanged bridge Lipschitz estimate;
- BR.17 uses the same far-pole rank-one loss;
- BR.18 uses BR.9 to suppress the opposite absolute channel;
- BR.20 uses exact signs, hollowness, and the uniform cap of `D_q`.

The constants `theta` and `lambda` are chosen once from
`(L,delta,C_P)`.  Hence the same base works for every pole in that one
uniform matched-roof family.  This does not assert a single base for
families with unboundedly varying `L` or vanishing `delta`.

## 3. Pullback after global orientation: PASS after repair

BT.3 supplies the appropriate one-sided thin tail for both signs, but BR.2
is written for a positive matched roof.  Set

```math
sigma=\begin{cases}
+1,&P(A)=Q(A),\\
-1,&N(A)=Q(A),
\end{cases}
\qquad \widetilde A=\sigma A.
```

Apply Theorem 21.8 and BR.2 to `widetilde A`.  If the resulting parent is

```math
\widetilde{\mathcal P}_{s,q}
=\begin{pmatrix}
\sigma A^s&B\\ B^T&D_q
\end{pmatrix},
```

then

```math
\sigma\widetilde{\mathcal P}_{s,q}
=\begin{pmatrix}
A^s&\sigma B\\ \sigma B^T&\sigma D_q
\end{pmatrix}.
```

It is a complete exact-sign context for the original switching child,
`Q(sigma widetilde{mathcal P})=Q(widetilde{mathcal P})`, and

```math
\sigma D_q=S_q(\sigma D_1)S_q.
```

The final source includes exactly this pullback.  Thus the children really
belong to the stated orbit of `A`, including when its negative cap is the
larger one.

## 4. Information bounds and their exact scope

Let the physical gap in EQ.8 be `c n^(3/2)`.  For a deterministic encoder
and one common decoder with pointwise error `epsilon n^(3/2)`, any
`epsilon<c/2` (the source's `c/3` is safe) forces injectivity on the
constructed set `I`: at query `i`, two children sharing a state would have
true responses at distance at most `2epsilon n^(3/2)`, contrary to EQ.8.
Therefore

```math
\log_2 |\operatorname{range}(\mathrm{Enc})|
\ge\log_2|I|
\ge {\gamma\over\log 2}n.
```

On the other hand, the switching action on any complete signing is free
modulo the global switch, by the same nonzero-edge argument used for the
fill orbit.  Relative to the public base `A`, exactly `n-1` bits identify a
child in the full orbit.  A decoder may reconstruct it and answer every
declared query by exhaustive optimization; no efficiency is claimed.  Thus
the proved fixed-accuracy response rate is exactly `Theta(n)` in order,
with bounds

```math
{\gamma\over\log2}n\le R_n(\epsilon)\le n-1.
```

The coefficient-identity problem itself has exactly `n-1` conditional bits.
The physical fixed-scale response theorem does **not** prove a sharp leading
constant or an `n-1`-bit response lower bound.

The same pigeonhole proof applies to any reusable compositional carrier,
congruence, feature algebra, or finite-precision state whose child-owned
range has `exp(o(n))` distinguishable values and from which all queries in
this public language are decoded uniformly.  It therefore rules out more
than a named encoding algorithm: it rules out every `o(n)`-bit reusable
state for this language.

It does not rule out:

- the exact `n-1`-bit switching label, or any other `O(n)`-bit carrier;
- a continuous state with uncharged infinite precision;
- query-dependent re-encoding of the child after the query is known;
- average-case or distributional error without an additional Fano-type
  argument;
- a narrower, gauge-covariant future language that co-switches the context;
- a proof concerned only with the isolated scalar cap;
- near-minimality of the order-`2n` parents, cross-order transfer, or
  convergence.

## 5. Public query description

The equivariant refinement genuinely removes quadratic *per-query* data.
The common bridge `B` and base fill `D_1` are public once per order.  A full
query is then generated from a projective pole in exactly `n-1` bits.  If
only the selected bank is declared and a public enumeration is stored, its
index costs `ceil(log_2|I|)` bits.  More cleanly, one may declare all
projective poles as the query language; the proof merely selects a
separating subbank inside it and no query codebook is needed to construct a
context.

The one-time base has `binom(n,2)` edge signs and can be stored using that
many bits.  The probabilistic existence proof neither supplies an explicit
subquadratic description nor proves that `Theta(n^2)` bits are necessary.

## 6. Archive comparison

| Archived item | Actual relation |
|---|---|
| BR.1--BR.3 / Theorem 36.25 | Supplies the biased selector one pole at a time.  EQ.1 is a strict simultaneous/query-description refinement, not a new scalarization inequality. |
| BT.3 / Theorem 36.26 | Supplies the all-bounded-cap thin tail.  Together with Theorem 21.8 and BR, it already contains the response-packing implication; EQ.3 packages the all-bounded-cap corollary and adds the common-base query bank. |
| BCL.0--BCL.1 | Already proves the `Theta(n)` upper/lower response rate abstractly from a tail, and physically for the regular-Walsh special family.  EQ extends physical scalarization to every bounded-cap complete signing through BT.3, but the information-theoretic pigeonhole and `n-1` upper bound are not new. |
| BCX | Already has an explicit bounded-cap anti-pin for a special regular-Hadamard switching code.  It does not give the present all-bounded-cap matched-roof compiler. |
| OV.2 | Supplies the earlier rectangular biased exact-sign rounding idea.  EQ's genuinely new observation is the square fill's simultaneous diagonal-switch orbit. |
| Theorem 36.22 | Already proves the universal metric switching broadcast and the free projective switch label.  EQ supplies the bounded-cap scalar future language. |

Accordingly the new theorem-level increment is narrow but useful: the
entire matched-roof physical query bank has one shared exact-sign base and
only a linear varying pole label.  The `Theta(n)` child response rate is an
existing implication sharpened in query-description scope, not a new rate
law.

## 7. Disposition

With the three repairs listed at the start implemented, the classification
is:

```text
PASS / PROVES:
  one common switched fill simultaneously preserves every BR.2 inequality;
  the varying physical query label costs exactly n-1 bits for the full bank;
  every bounded-cap switching orbit has Theta(n) bounded-cap response bits.

SCOPE:
  the lower bound excludes exp(o(n))-state reusable carriers for this
  A-dependent public future language, but not O(n)-bit states or arbitrary
  compositional proofs.

NO NEW:
  convergence, cross-order transfer, near-minimality of parents, or a sharp
  response-rate constant.
```

No normalization, simultaneous-existence, information-counting, or archive-
scope defect remains.  The small displayed source phrase `i ne j` in the
stabilizer explanation is typographical prose for `i\ne j` and does not
affect the argument.
