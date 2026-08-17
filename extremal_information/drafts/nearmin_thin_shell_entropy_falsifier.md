# Thin-shell entropy: a sharp lower scale, the Walsh boundary, and the halo obstruction

Date: 2026-08-17.

Status: rigorous theorem/counterexample-first report.  The deterministic
statements below use only the exact flip identity.  The Walsh count imports
two stated primary-source theorems and checks the normalization explicitly.

## 1. Question and verdict

Consider the candidate that every genuine near-minimizer has at most

```math
\exp\{O(\sqrt n\,\operatorname {polylog}n)\}
```

oriented cuts within additive `O(n)` of its absolute cap.  None of the
bounded-cap families tested here falsifies this candidate.  There is,
however, a sharp obstruction to strengthening it:

* every bounded-cap signing already has
  `exp(Omega(sqrt(n) log n))` cuts in an `O(n)` shell;
* the regular Walsh/PC.3 tensor family has that many **exact** top cuts;
* present counting results do not rule out exponentially many self-dual
  bent functions, so even the exact Walsh-shell upper bound is open;
* an `o(n^(3/2))` coefficient-edit halo does not preserve an `O(n)` shell.
  Black-box edit stability preserves that thin scale only for `O(n)` edits.

Thus the proposed upper bound is a precise open candidate, essentially
optimal in its `sqrt(n) log n` exponent if true.  It is a shell-cardinality
statement, not a contextual response-entropy theorem.

Throughout,

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad Q(A)=\max_x|H_A(x)|.
```

## 2. Every bounded-cap signing has a large `O(n)` shell

### Theorem TS.1 (universal thin-shell entropy floor)

Let `A` be a hollow signing of order `n` with

```math
Q(A)\le Cn^{3/2}.
```

Fix `c>0` and let `k=floor(c sqrt(n))`.  For all sufficiently large `n`,
there are at least

```math
{1\over2}{\lfloor n/2\rfloor\choose k}
=2^{(c/2+o(1))\sqrt n\log_2n}                    \tag{TS.1}
```

projectively distinct spins with one common energy orientation and absolute
deficit at most

```math
(8cC+2c^2+o(1))n.                                \tag{TS.2}
```

Consequently their correspondingly oriented signed cut words are distinct.

#### Proof

Orient a ground state `x` so `H_A(x)=Q`, replacing `A` by `-A` if needed,
and put

```math
\ell_i=x_i(Ax)_i,
\qquad s_{ij}=a_{ij}x_ix_j.
```

One-spin optimality and double counting give

```math
\ell_i\ge0,
\qquad \sum_i\ell_i=2Q.                          \tag{TS.3}
```

Hence the set

```math
J=\{i:\ell_i\le4Q/n\}                            \tag{TS.4}
```

has cardinality at least `n/2`.  If `x^S` is obtained by flipping `S`,
direct expansion gives

```math
H_A(x^S)
=Q-2\sum_{i\in S}\ell_i
  +4\sum_{\{i,j\}\subseteq S}s_{ij}.             \tag{TS.5}
```

For every `S subset J` of size `k`, therefore,

```math
Q-|H_A(x^S)|
\le Q-H_A(x^S)
\le {8kQ\over n}+2k(k-1),                         \tag{TS.6}
```

which is (TS.2).  The `binom(|J|,k)` spins are projectively distinct because
`k<n/2`.  Pigeonhole on the sign of their energies retains half with one
common orientation.  On a complete graph, two different projective spins
give different cut words. `square`

The extra `log n` relative to the affine-cube theorem LA.1 is important.
LA.1 freezes one set of `Theta(sqrt n)` low-field coordinates to obtain an
odd-product algebra of size `2^(Theta(sqrt n))`.  TS.1 takes **all**
`sqrt(n)`-subsets of a linear-size low-field pool.  It proves cardinality,
but that union has no asserted small composition algebra.

## 3. Walsh matrices attain the same scale at exact cap

Let `d` be even, `N=2^d`, and let `W_d` be the unnormalized Sylvester Walsh
matrix.  Thus

```math
W_d^2=NI,
\qquad \operatorname {tr}W_d=0.
```

Put `A_d=W_d-diag(W_d)`.  This is a hollow signing and, for every Boolean
`z`,

```math
H_{A_d}(z)={1\over2}z^TW_dz.                       \tag{TS.7}
```

### Theorem TS.2 (many exact Walsh maximizers)

The number of projectively distinct positive exact maximizers of `A_d`
obeys

```math
\log_2 |\mathcal G_d^+|
\ge(1/2+o(1))\sqrt N\log_2N,                       \tag{TS.8}
```

and

```math
Q(A_d)={1\over2}N^{3/2}.                           \tag{TS.9}
```

The same conclusion holds for the unflipped hollow child matrices in the
diagonal-conjugate PC.3 tensor tower.  It is not asserted for their
sparse-flipped or physically completed descendants.

#### Proof and imported inputs

A sign vector `F` is a positive `sqrt(N)` eigenvector of `W_d` exactly when
it is the sign of a self-dual bent Boolean function.  The spectral bound
gives (TS.9), and each such vector attains it.

For even `d>=6`, Carlet--Danielsen--Parker--Sole, Theorem 4.9 and the
explicit example after it, give an injective lift from every bent function
`f` in `d-2` variables
to a self-dual bent function in `d` variables whose sign vector is

```math
(F,\widetilde F,\widetilde F,-F).                  \tag{TS.10}
```

The first block makes the lift injective.  Haugland's 2025 bent-square
lower bound states, for even `t`,

```math
\log_2 b_t\ge t2^{t/2}(1+O(1/t)),                  \tag{TS.11}
```

where `b_t` counts bent functions on `t` variables.  Apply this with
`t=d-2` and divide by two for global complementation/projectivization:

```math
\log_2 |\mathcal G_d^+|
\ge(d-2)2^{(d-2)/2}(1+o(1))-1
=(1/2+o(1))\sqrt N\log_2N.
```

The order-16 regular Walsh matrix used by PC.3 is a diagonal conjugate of
`W_4`; its tensor powers are diagonal conjugates of `W_(4j)`, so Boolean
eigenvectors are carried bijectively. `square`

Primary sources:

* C. Carlet, L. E. Danielsen, M. G. Parker, P. Sole,
  [Self-Dual Bent Functions](https://ii.uib.no/~larsed/sdbent.pdf),
  *International Journal of Information and Coding Theory* 1 (2010),
  384--399, especially Theorem 4.9.
* J. K. Haugland,
  [A lower bound on the number of bent squares](https://arxiv.org/abs/2508.14605)
  (2025), main counting theorem.

For comparison, exhaustive enumeration at order `N=16` gives 20 positive
and 20 negative raw exact Walsh eigenvectors, hence 10 positive and 10
negative projective classes.  The four poles used in the PC.3 presentation
are a strict subfamily of this full exact shell.

### What TS.2 does not prove

The best imported lower bound is
`2^(Theta(sqrt(N) log N))`, not `exp(cN)`.  Known upper bounds for all bent
functions remain exponential in `N`, so they do not prove the proposed
thin-shell upper bound even for Walsh.  In particular, the possibility that
the Walsh eigenspace has `exp(cN)` Boolean vertices is not excluded by the
results used here.  Walsh **realizes** the universal lower scale; no matching
Walsh-shell upper bound is claimed.  It is an open counting problem, not a
counterexample.

## 4. Coefficient-edit halos do not preserve the thin scale

For `Delta>=0`, write

```math
\mathcal S_A(\Delta)
=\{x:Q(A)-|H_A(x)|\le\Delta\}.                     \tag{TS.12}
```

### Lemma TS.3 (black-box shell stability under edge edits)

Let `A` be an exact minimizer of order `n`, and let `B` be obtained from it
by flipping `r` edge signs.  Then

```math
\mathcal S_B(\Delta)\subseteq\mathcal S_A(\Delta+2r),
\qquad
\mathcal S_A(\Delta)\subseteq\mathcal S_B(\Delta+4r).       \tag{TS.13}
```

#### Proof

Uniformly in `x`,

```math
|H_A(x)-H_B(x)|\le2r.                               \tag{TS.14}
```

Also `Q(B)>=M_n=Q(A)` by minimality and `Q(B)<=M_n+2r`.  If
`x in S_B(Delta)`, then

```math
|H_A(x)|\ge |H_B(x)|-2r\ge Q(B)-\Delta-2r
\ge M_n-\Delta-2r.
```

This is the first inclusion.  If `x in S_A(Delta)`, then

```math
Q(B)-|H_B(x)|
\le(M_n+2r)-(M_n-\Delta-2r)=\Delta+4r,
```

proving the second. `square`

Therefore an `o(n^(3/2))` edit construction transfers only an
`o(n^(3/2))` shell by coefficientwise Lipschitz control.  To retain an
`O(n)` shell this argument requires `r=O(n)`.  Any larger edit construction
needs a new cap-relative cancellation theorem, not merely proximity in the
near-minimizer halo.

In particular, putting the Walsh signing within `o(N^(3/2))` edits of an
exact minimizer would imply

```math
M_N\ge {1\over2}N^{3/2}-o(N^{3/2}),                 \tag{TS.15}
```

by (TS.14), while the Walsh construction supplies the reverse inequality.
That proposed transfer would already establish the conjectural constant
`1/2` on the Walsh subsequence.  It cannot be treated as a cheap halo move.

## 5. A precise mesoscopic falsifier

The exact flip formula isolates what would decisively kill the proposed
upper bound.  At an oriented ground state define, for `J subset [n]`,

```math
L(J)=\sum_{i\in J}\ell_i,
\qquad
P(J)=\max_{S\subseteq J}
\left|\sum_{\{i,j\}\subseteq S}s_{ij}\right|.      \tag{TS.16}
```

### Proposition TS.4 (low-field mesoscopic block criterion)

For every `S subset J`,

```math
Q-|H_A(x^S)|\le2L(J)+4P(J).                         \tag{TS.17}
```

Consequently, if a sequence of genuine near-minimizers has ground states
and sets `J_n` such that, for a specified exponent `C_0`,

```math
{|J_n|\over\sqrt n(\log n)^{C_0}}\longrightarrow\infty,
\qquad L(J_n)+P(J_n)=O(n),                          \tag{TS.18}
```

then its `O(n)` shell has
`2^{|J_n|-1}` same-orientation projective states and the candidate is false.
Here the known lower bound `Q=Omega(n^(3/2))` for genuine near-minimizers
ensures that the `O(n)` defect cannot reverse the ground orientation for
large `n`; the factor `1/2` is retained to cover projectivization uniformly.
To falsify every unspecified polylogarithmic exponent at once, require the
displayed ratio to diverge for every fixed `C_0`.

This is a serious falsifier rather than vocabulary.  A pseudorandom signing
on `k` vertices naturally has subset discrepancy on the `k^(3/2)` scale;
at `k=n^(2/3)` this is `O(n)`.  Thus the unresolved issue is whether genuine
near-minimality excludes a mesoscopic induced block with total ground-state
local-field mass `O(n)`.  The identity `sum_i ell_i=2Q=Theta(n^(3/2))`
alone does not exclude it.

## 6. Exact shell, thin shell, and contextual entropy are different

1. **Exact-shell cardinality.**  TS.2 gives
   `2^(Omega(sqrt n log n))` exact Walsh maximizers.  It does not give a
   response packing.
2. **`O(n)`-shell entropy.**  TS.1 gives the same exponent for every
   bounded-cap signing.  The states share a concise ground-state/low-field
   description and need not encode that many reusable response states.
3. **Contextual response entropy.**  Theorem 21.8 in the repository starts
   from an entropy deficit in a much thicker, fixed-`n^(3/2)` shell and uses
   a random low-operator bridge to obtain `Omega(n)` response bits.  Its
   conclusion neither follows from nor bounds the `O(n)` shell above.

Conflating these quantities would reverse the interpretation: a large thin
shell can be evidence of many witnesses while still possessing a very small
generative presentation.

## 7. Director judgment

The proposed

```math
#\mathcal S_A(Kn)\le\exp\{O_K(\sqrt n\,\operatorname {polylog}n)\}
```

for genuine near-minimizers survives this audit, but it is not presently a
selected compositional lemma:

* TS.1 shows it could improve the exponent only by polylogarithmic factors;
* Walsh/PC.3 does not falsify it and instead realizes the lower scale;
* proving it even for Walsh requires new self-dual-bent counting input;
* shell cardinality alone has no proved cap-relative composition payoff;
* Proposition TS.4 gives the most discriminating next falsifier: search for
  mesoscopic low-field, low-subset-discrepancy blocks in certified
  near-minimizers.

No construction has been certified inside an `o(n^(3/2))` halo of exact
minimizers while retaining `exp(cn)` states at `O(n)` deficit.  TS.3 explains
why the obvious sparse-edit route does not supply one.
