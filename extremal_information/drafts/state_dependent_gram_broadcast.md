# State-dependent Gram broadcast at the total extremal scale

Status: task-local rigorous draft.  This answers the escape question left by
Theorem 21.23 inside a deliberately structured, semantic class of dense
bridges.  It does **not** encode independent bits on independent edges.  A
hidden alternating form is evaluated on one public list of labels, so all
edge phases form a bicharacter code and obey its linear identities.  A fixed
public base signing is then dressed by those phases.

The conclusion is a positive-rate family whose scalar contextual responses
are separated by a fixed multiple of `N^(3/2)`, while every child and every
pairwise contrast remains on that same scale.  Thus state-dependent
cross-block coefficients genuinely escape the state-local ceiling; they do
so by broadcasting each compatibility coordinate into quadratically many
interaction atoms.

## 1. A public alternating-form sampler

Let `V=F_2^r`, and let `Alt(V)` be the vector space of alternating bilinear
forms on `V`.  In characteristic two these forms are symmetric with zero
diagonal, and

```math
h:=\dim Alt(V)={r(r-1)\over2}.                         \tag{GB.1}
```

For a public ordered list `P=(p_1,...,p_k) in V^k`, define the evaluation
word of `B in Alt(V)` by

```math
c_B(i,j)=B(p_i,p_j),\qquad 1\le i<j\le k.             \tag{GB.2}
```

When both coordinate orders are used, extend this symmetrically by
`c_B(j,i)=c_B(i,j)` and put `c_B(i,i)=0`.

### Lemma GB.1 (linear-size labels test every alternating form densely)

For every `r>=2`, with

```math
k=64r^2,                                               \tag{GB.3}
```

there is a list `P in V^k` such that every nonzero `B in Alt(V)` satisfies

```math
#\{i<j:B(p_i,p_j)=1\}
\ge {1\over4}{k\choose2}
={k(k-1)\over8}.                                      \tag{GB.4}
```

In particular `B -> c_B` is an injective linear binary code of dimension
`h`, block length `binom(k,2)`, and relative distance at least `1/4`.

#### Proof

Choose the `p_i` independently and uniformly from `V`.  If `B!=0`, its rank
is even and at least two.  For independent uniform `p,q`, conditioning on
whether `p` lies outside the radical gives

```math
Pr\{B(p,q)=1\}
={1-2^{-rank(B)}\over2}\ge {3\over8}.                 \tag{GB.5}
```

Put `X_B=sum_(i<j)1_{B(p_i,p_j)=1}`.  Changing one `p_i` changes at most
`k-1` summands.  Since

```math
E X_B-{1\over4}{k\choose2}
\ge {1\over8}{k\choose2}={k(k-1)\over16},            \tag{GB.6}
```

McDiarmid's inequality yields

```math
Pr\left\{X_B<{1\over4}{k\choose2}\right\}
\le \exp(-k/128).                                     \tag{GB.7}
```

There are fewer than `2^h` nonzero forms.  For (GB.3), the union-bound
logarithm is at most

```math
h\log2-k/128
\le {\log2\over2}r(r-1)-{r^2\over2}<0.               \tag{GB.8}
```

Hence a list satisfying (GB.4) exists. `square`

The code is strongly constrained.  It is a linear evaluation code, and for
every relation `sum_j lambda_jp_j=0`,

```math
sum_j lambda_jc_B(i,j)=0\quad\hbox{for every }i.       \tag{GB.9}
```

Thus its `h=Theta(k)` information coordinates do not specify the
`Theta(k^2)` edge phases independently.

## 2. Dense coefficient disagreement forces a Boolean response

For a real symmetric hollow `k by k` matrix `D`, write

```math
Q(D)=\max_{s\in\{+-1\}^k}
\left|\sum_{i<j}D_{ij}s_is_j\right|.                  \tag{GB.10}
```

### Lemma GB.2 (support-only discrepancy lower bound)

If `D_ij in {0,+-2}` and at least `k^2/16` unordered entries are nonzero,
then

```math
Q(D)\ge {\sqrt2\over32}k^{3/2}.                       \tag{GB.11}
```

#### Proof

A random vertex bipartition cuts each nonzero edge with probability `1/2`.
Fix a partition `(L,R)` with at least `e>=k^2/32` nonzero cross entries.
For uniform `y in {+-1}^R`, optimize `x_i` separately on `L`.  The sharp
`p=1` Khintchine inequality gives

```math
\begin{aligned}
\max_{x,y}|x^TD_{LR}y|
&\ge \sum_{i\in L}E_y\left|\sum_{j\in R}D_{ij}y_j\right|\\
&\ge \sqrt2\sum_{i\in L}\sqrt{d_i}
 \ge {\sqrt2 e\over\sqrt k}
 \ge {\sqrt2\over32}k^{3/2},                         \tag{GB.12}
\end{aligned}
```

where `d_i` is the number of nonzero entries in row `i`, and
`sqrt(d_i)>=d_i/sqrt(k)`.  For the two full sign vectors `(x,y)` and
`(-x,y)`, the within-part terms agree while the cross term changes sign.
Therefore the larger absolute full quadratic value is at least the cross
value in (GB.12). `square`

The lemma is deliberately sign-blind: no pseudorandomness assumption on the
nonzero coefficients is used.

## 3. A spectrally flat bicharacter coset

For the sampler in Lemma GB.1 put

```math
\chi_B(i,j)=(-1)^{B(p_i,p_j)}.                        \tag{GB.13}
```

These are alternating bicharacters restricted to the public label list:
`chi_(B+C)=chi_B chi_C`, and (GB.9) supplies many parity checks.

### Lemma GB.3 (one base signing makes the whole coset flat)

There is a hollow symmetric signing `A` such that every

```math
A_B(i,j)=A(i,j)\chi_B(i,j),\qquad B\in Alt(V),         \tag{GB.14}
```

satisfies

```math
\|A_B\|_(2\to2)\le8\sqrt k.                          \tag{GB.15}
```

#### Proof

Choose the upper-triangular entries of `A` independently and uniformly from
`{+-1}`.  Fix `B` and a unit vector `z`.  Then

```math
z^TA_Bz=2\sum_{i<j}A_{ij}\chi_B(i,j)z_iz_j
```

is a Rademacher sum whose squared coefficient norm is at most two.  Hence

```math
Pr\{|z^TA_Bz|>t\}\le2e^{-t^2/4}.                     \tag{GB.16}
```

A `1/4`-net of the unit sphere has at most `9^k` points and satisfies
`||M||_op<=2max_z|z^TMz|` for symmetric `M`.  Taking `t=4sqrt(k)` and
union-bounding over the net and all `2^h` forms gives failure probability at
most

```math
2\exp\{k\log9+h\log2-4k\}<1,                         \tag{GB.17}
```

because `h/k<1/128`.  Some `A` therefore obeys (GB.15) simultaneously.
`square`

Thus this is not a family whose response separation comes from a hidden
`Theta(k^2)` ferromagnetic term: every individual state has Boolean cap at
most `4k^(3/2)`.

## 4. Total-scale contextual packing

The hidden child with state `B` is the ordinary quadratic landscape

```math
H_B(s)={1\over2}s^TA_Bs.                              \tag{GB.18}
```

Declare, before the hidden state is chosen, the public same-support additive
overlay language

```math
\mathcal T=\{-H_T:T\in Alt(V)\}.                      \tag{GB.19}
```

A continuation is specified by its own alternating form `T`; it does not
inspect `B`.  Define the scalar absolute response

```math
R_T(B)=\max_s|H_B(s)-H_T(s)|=Q(A_B-A_T).              \tag{GB.20}
```

### Theorem GB.4 (state-dependent Gram broadcast escapes the local ceiling)

For every `r>=2`, there are `k=64r^2`, a public label list, a public base
signing, and `2^h` hidden states, `h=r(r-1)/2`, such that

```math
R_B(B)=0,                                             \tag{GB.21}
```

whereas for every `B!=T`,

```math
{\sqrt2\over32}k^{3/2}
\le R_T(B)\le8k^{3/2}.                                \tag{GB.22}
```

Consequently the contextual response vectors are pairwise separated by
`(sqrt(2)/32)k^(3/2)`.  Any summary answering every continuation to uniform
error below `(sqrt(2)/64)k^(3/2)` needs at least

```math
h={r(r-1)\over2}\ge{k\over256}                       \tag{GB.23}
```

bits.  This is a positive information rate at the total `k^(3/2)` extremal
scale.

#### Proof

For `B!=T`, the difference `L=B+T` is a nonzero alternating form.  By
Lemma GB.1, `A_B-A_T` has entries of magnitude two on at least
`k(k-1)/8>=k^2/16` unordered pairs and is zero elsewhere.  Lemma GB.2 gives
the lower bound in (GB.22).  Lemma GB.3 gives

```math
Q(A_B-A_T)
\le {k\over2}\|A_B-A_T\|_op
\le8k^{3/2},                                          \tag{GB.24}
```

which proves the upper bound.  If two children shared an approximate
summary, their decoded responses at context `-H_T` would differ by less
than twice the allowed error, contradicting (GB.21)--(GB.22).  Finally,
`r>=2` gives `h>=(r^2)/4=k/256`. `square`

The same response theorem lifts to genuine dense block bridges.  Let `W_n` be any
symmetric regular Hadamard matrix of order `n`, with Boolean vector `u`
satisfying `W_nu=sqrt(n)u`, and use the cross block

```math
(A_B)_{ij}W_n.                                        \tag{GB.25}
```

For `N=kn`, the operator upper bound and the Boolean witness `s_i u` turn
(GB.22) into

```math
{\sqrt2\over32}N^{3/2}
\le R_T^{(n)}(B)\le8N^{3/2}.                          \tag{GB.26}
```

Common hollow sign diagonal blocks may be inserted to make every child an
exact hollow sign matrix; they cancel in the contrast continuation.  If an
individual-child `O(N^(3/2))` cap is also claimed while `n` varies, choose
these common diagonal blocks with operator norm `O(sqrt(n))` (a random
hollow signing supplies such a choice).  Taking `n=1` already gives the exact
hollow signings and cap bound in Theorem GB.4.

The hidden information is `h=Theta(k)=Theta(N/n)` bits.  It therefore has a
positive rate per total vertex for fixed `n`, including `n=1`, but rate
`Theta(1/n)` if the inner lift order itself grows.  The normalized response
gap in (GB.26) stays fixed in either regime.

## 5. Why this is semantic rather than arbitrary coefficient storage

The hidden object is the alternating form `B`, not an edge table.  Every
edge phase is forced by the public labels through one law (GB.13), the image
has only `h=Theta(k)` bits inside `Theta(k^2)` edge positions, and the parity
identities (GB.9) hold.  Equivalently, every `B` is the pullback of an
alternating (possibly degenerate, hence presymplectic) Gram form, so
`chi_B(p,q)` is the commutator phase of the
corresponding Weyl/Walsh labels.  This is precisely compatibility data.

This hidden compression is conditional on the shared public pair `(A,P)`.
The probabilistically chosen base signing `A` itself has `Theta(k^2)`
nonuniform bits; the theorem does not compress arbitrary dense signings or
give a uniform explicit constructor for that base.

The child and continuation are also distinct.  A child supplies `A_B`; a
public future supplies its own `-A_T` from the fixed declared language.
The future never reads `B`.  Contextual separation merely chooses, after
two states are compared, one member of a language that existed for every
state.

The child and each negative-clone future are exact signings separately, but
their same-edge overlay has coefficients in `{0,+-2}`.  Thus the response
language is a legitimate predeclared additive contextual language, not a
claim that the composed parent remains in the exact-sign class or that the
future lives on disjoint appended variables.

What makes the escape possible is visible quantitatively.  Flipping any one
nonzero compatibility coordinate changes a nonzero alternating form, hence
at least `k^2/16` bridge atoms.  The `Omega(k^(3/2))` influence threshold in
Theorem 21.23 is therefore exceeded by a factor `Omega(sqrt(k))`.  The
construction does not refute that theorem; it realizes its only remaining
bounded-coefficient escape.

## 6. Scope and falsification value

1. The continuation language has `2^h` members.  The theorem proves a
   response-information lower bound, not a small decoder or an efficient
   optimization algorithm.
2. The result uses absolute quadratic response, matching the motivating
   extremal normalization.  It does not claim the same packing for a single
   one-sided upper response.
3. The construction is an algebraically structured dense bridge family,
   not a statement about arbitrary near-minimizers or the original
   convergence problem.
4. It decisively falsifies any proposed total-scale ceiling based only on
   the number `h=O(k)` of hidden coordinates, or on the bounded cap/spectral
   norm of every child.  A valid ceiling must also bound the incidence with
   which one hidden coordinate can alter interaction atoms, or quotient the
   joint child--continuation bicharacter gauge.
5. The public base is shared nonuniform advice of quadratic size, and the
   context acts by same-support overlay.  Both limitations are part of the
   theorem's declared architecture.

The accompanying exact wind tunnel is
[`../experiments/verify_state_dependent_gram_broadcast.py`](../experiments/verify_state_dependent_gram_broadcast.py).
