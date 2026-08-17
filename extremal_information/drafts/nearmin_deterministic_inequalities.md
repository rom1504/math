# Deterministic inequalities for quadratic-cap near-minimizers

Date: 2026-08-17.

Status: independent extremal-combinatorics/spectral report.  The three
candidate statements in Section 3 were frozen before reading
`ACTIVE_STATE.md` or the current contextual-response drafts.  Sections 4--6
record the subsequent archive comparison.  Nothing here is promoted to the
canonical theorem files.

## 1. Normalization and elementary geometry

Let

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
Q(A)=\max_x|H_A(x)|,\qquad
M_n=\min_AQ(A),
\]

where `A` is hollow, symmetric, and has signs off the diagonal.  Write
`N=binom(n,2)` and identify `A` with its edge-sign word `a`.  The augmented
cut code is

\[
\mathcal V_n=\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\}.
\]

For `n>=3`, it has `2^n` distinct words and

\[
Q(a)=\max_{v\in\mathcal V_n}a\mathbin\cdot v
     =N-2d(a,\mathcal V_n).
\tag{1.1}
\]

Thus an exact minimizer is a deepest hole of `\mathcal V_n`, and if `rho_n`
is its covering radius, then

\[
M_n=N-2\rho_n,\qquad
Q(a)-M_n=2\bigl(\rho_n-d(a,\mathcal V_n)\bigr).
\tag{1.2}
\]

Two elementary invariances are important when assessing any proposed
observable.

* Switching `A` by `DAD`, with `D` diagonal and Boolean, just permutes the
  values of `H_A`; hence `Q` is switching invariant.  In edge language this
  flips a cut `delta(S)`.
* If `A` and `B` differ on `r` edges, then

  \[
  |Q(A)-Q(B)|\le 2r.
  \tag{1.3}
  \]

Consequently every exact minimizer has an `o(n^(3/2))` Hamming halo of
asymptotic near-minimizers.  A property not stable under such edits cannot be
forced on the class requested in the prompt.

## 2. Consequences requiring no structural ansatz

### 2.1 The exact near-minimal flip inequality

For `v in \mathcal V_n`, define its cap deficit and its signed correlation
with an edge set `F` by

\[
d_A(v)=Q(A)-a\mathbin\cdot v\ge0,
\qquad
c_F(v)=\sum_{e\in F}a_ev_e.
\tag{2.1}
\]

Let

\[
\eta(A)=Q(A)-M_n.
\]

If `a^F` is obtained by flipping precisely `F`, then

\[
Q(a^F)=\max_v\{Q(A)-d_A(v)-2c_F(v)\}.
\tag{2.2}
\]

Since `Q(a^F)>=M_n=Q(A)-eta(A)`, a maximizer for `a^F` gives the exact
certificate

\[
\boxed{
\text{for every }F\subseteq E(K_n)\text{ there is }v\in\mathcal V_n
\text{ with }d_A(v)+2c_F(v)\le\eta(A).
}
\tag{2.3}
\]

In particular,

\[
c_F(v)\le\frac{\eta(A)}2,
\qquad
d_A(v)\le\eta(A)+2|F|.
\tag{2.4}
\]

If `N_v={e:a_ev_e=-1}`, (2.3) is equivalently

\[
2|F\cap N_v|
\ge |F|+\frac{d_A(v)-\eta(A)}2.
\tag{2.5}
\]

For an exact minimizer, every single edge is opposed by a state at deficit
at most two.  More generally, every `r`-edge set is at least half opposed by
one state at deficit at most `2r`.  For an `eta`-near minimizer this remains
informative only on flip scales appreciably larger than `eta`.

The quantifiers matter.  The witnessing state in (2.3) can depend on `F`.
There is no common probability law or common active face in this argument.

The full family (2.3), over all `F`, is exactly equivalent to the numerical
inequality `Q(A)-M_n<=eta`: every other signing is some `a^F`.  Its
restriction to `|F|<=r` is genuinely weaker and is the only version treated
as structural below.

### 2.2 Heredity over a sub-square-root order window

If `S` has `m` vertices, then

\[
\boxed{Q(A[S])\le Q(A).}
\tag{2.6}
\]

Indeed, fix a spin on `S` and extend it by independent unbiased spins on the
complement.  The expected full energy is exactly the energy on `S`, and its
absolute value is therefore at most `Q(A)`.

Conversely, extending an optimal order-`m` signing by arbitrary signs gives

\[
M_n\le M_m+\binom n2-\binom m2.
\tag{2.7}
\]

Combining (2.6)--(2.7), every principal `m`-restriction of a signing with
`Q(A)<=M_n+eta` obeys

\[
\boxed{
Q(A[S])-M_m
\le \eta+\binom n2-\binom m2.
}
\tag{2.8}
\]

For `m=n-k`, the extra loss is `k(2n-k-1)/2`.  Hence every deletion of
`k=o(sqrt(n))` vertices preserves asymptotic near-minimality at the
`n^(3/2)` scale.  This is uniform over all restrictions, not merely most of
them.  It covers only a vanishing relative order window and does not by
itself yield a cross-scale recurrence.

### 2.3 Local fields at an absolute ground state

Orient an absolute ground state so that `H_A(x)=Q(A)` and put

\[
\ell_i=x_i(Ax)_i.
\]

A one-spin flip and the two-sided cap give

\[
\boxed{0\le\ell_i\le Q(A),\qquad \sum_i\ell_i=2Q(A).}
\tag{2.9}
\]

Thus the average oriented field is `2Q(A)/n=O(sqrt(n))` on every bounded-cap
sequence, and

\[
|\{i:\ell_i>t\sqrt n\}|\le\frac{2Q(A)}{t\sqrt n}.
\tag{2.10}
\]

There is no corresponding square-tail conclusion: a vanishing set of
vertices can carry fields of order `n` while costing only
`o(n^(3/2))` in the cap.

### 2.4 What the cap does force spectrally

Let `u` be a unit eigenvector of `A`, with eigenvalue `lambda`, and put
`r=||u||_infty`.  Choose independent signs `X_i` with
`E X_i=u_i/r`.  Since the diagonal is zero,

\[
\mathbb E H_A(X)=\frac{u^TAu}{2r^2}.
\]

Therefore

\[
\boxed{|\lambda|\le2Q(A)\,||u||_\infty^2.}
\tag{2.11}
\]

In particular, a `K/sqrt(n)`-flat eigenvector of a bounded-cap signing has
eigenvalue `O(K^2 sqrt(n))`.  Large eigenvalues are therefore forced to be
localized, but the operator norm itself need not be `O(sqrt(n))`.

## 3. The three frozen candidate lemmas

### Candidate 1: mesoscopic flip-shell coverage

For `Delta>=0`, put

\[
\mathcal T_\Delta(A)
=\{v\in\mathcal V_n:a\mathbin\cdot v\ge Q(A)-\Delta\}.
\]

The candidate was that near-minimality forces a nontrivial collective
coverage property on every mesoscopic flip scale.  The precise statement is
already a theorem:

> If `Q(A)<=M_n+eta`, then for every `r` and every edge set `F` with
> `|F|<=r`, there is `v in \mathcal T_{eta+2r}(A)` such that
> `sum_(e in F) a_e v_e<=eta/2`.

This is (2.3)--(2.4).  For `eta=o(r)` it says that at least
`(1/2-o(1))|F|` edges of `F` are opposed by one near-cap state.

This statement is demonstrably weaker than computing `Q`.  The archived
exact order-11 signing of cap 19 is stable under every batch of at most four
edge flips, but a five-edge batch lowers its cap to 17, which is `M_11`.
Thus even the `eta=0`, `r=4` form can hold for a nonglobal signing.

Verdict: **proved, but too weak for the current transfer problem.**  The
state depends on `F`; separate witnesses need not admit one balanced law.
The all-`F` strengthening ceases to be a reduction because it is equivalent
to near-minimality itself.

### Candidate 2: thin-shell fractional eutaxy

To express the missing common-law quantifier, define

\[
\tau_\Delta^+(A)
=\min_{\mu\in\mathcal P(\mathcal T_\Delta(A))}
  \max_e\mathbb E_\mu[a_ev_e].
\tag{3.1}
\]

Finite minimax gives the equivalent dual form

\[
\tau_\Delta^+(A)
=\max_{w\in\Delta_N}\min_{v\in\mathcal T_\Delta(A)}
  \sum_e w_ea_ev_e.
\tag{3.2}
\]

Thus this records only `N` first moments of one near-cap law; it neither
records all energies nor determines `Q`.

The frozen strong candidate was:

> For some absolute `C`, every sequence with
> `Q(A_n)-M_n=o(n^(3/2))` satisfies
> `tau^+_(Cn)(A_n)=o(1)`.

It is false by a generic geodesic surgery, not a special algebraic example.
Start from an exact minimizer `a` and an active word `v_0`, so
`a dot v_0=M_n`.  Its negative support

\[
N_0=\{e:a_e(v_0)_e=-1\}
\]

has size `(N-M_n)/2`.  Choose `F subset N_0`, `|F|=r`, and flip `F` to
obtain `b=a^F`.  Then

\[
b\mathbin\cdot v_0=M_n+2r,
\qquad
Q(b)\le M_n+2r,
\]

so in fact

\[
\boxed{Q(b)=M_n+2r.}
\tag{3.3}
\]

The known (and elementary probabilistic) bound `M_n=O(n^(3/2))` makes
`|N_0|=(1/2+o(1))N`, so all choices `r=o(n^2)` used below are available.

For any `v`, let

\[
k_F(v)=|\{e\in F:b_ev_e=-1\}|.
\]

Since the signs on `F` were flipped toward `v_0`, one has

\[
b\mathbin\cdot v
\le Q(b)-4k_F(v).
\tag{3.4}
\]

Consequently every law supported on `T_Delta(b)` satisfies

\[
\sum_{e\in F}\Pr\{b_ev_e=-1\}\le\frac\Delta4.
\]

Some planted edge therefore has aligned bias at least

\[
\boxed{
\tau_\Delta^+(b)\ge1-\frac{\Delta}{2r}
                 =1-\frac{\Delta}{Q(b)-M_n}.
}
\tag{3.5}
\]

Take, for example, `r=floor(n^(4/3))`.  Then `b` is an
`o(n^(3/2))`-near-minimizer, while for every `Delta=O(n)`,
`tau_Delta^+(b)->1`.  More generally, no universal eutaxy theorem can use a
shell `Delta=o(Q(A)-M_n)`.  In coding language, moving `r` steps from a deep
hole toward one nearest codeword produces a certified radius-`r` near-deep
hole whose nearest-codeword face freezes all `r` moved coordinates.

This does not disprove a statement with shell thickness at least the full
optimality excess, for example `Delta>=eta+O(n)`.  That is the precise
remaining boundary, not a promoted fourth candidate.

There is an informative existential counterpart.  For `lambda>0`, define

\[
Z_\lambda(a)=\sum_{v\in\mathcal V_n}e^{\lambda a\cdot v}
\]

and choose a signing `a_lambda` minimizing `Z_lambda`.  Under its Gibbs law,
one-edge comparison gives

\[
\mathbb E[a_ev_e]\le\tanh\lambda\quad\text{for every }e.
\tag{3.6}
\]

Also

\[
Q(a_\lambda)\le M_n+\frac{n\log2}{\lambda},
\qquad
\Pr\{d_{a_\lambda}(v)\ge t\}
\le2^ne^{-\lambda t}.
\tag{3.7}
\]

Taking `lambda=b_n/sqrt(n)`, where `b_n->infinity` and
`b_n=o(sqrt(n))`, and conditioning on a shell of width
`O(n^(3/2)/b_n)`, produces a specially selected
`o(n^(3/2))`-near-minimizer with one common law satisfying

\[
\max_e\mathbb E[a_ev_e]
\le\frac{b_n}{\sqrt n}+o(1).
\tag{3.8}
\]

The universal candidate fails while this soft-selection theorem holds.  The
gap is exactly the quantifier requested by the campaign: near-minimality of
an arbitrary signing does not imply minimization of its finite-temperature
partition function.

Exploratory exact LPs support an exact-minimizer-only remnant.  On every
minimizer orbit through order 8, the values of `tau_2^+` are

\[
\begin{array}{c|cccccc}
n&3&4&5&6&7&8\\ \hline
\tau_2^+&1/3&1/2&2/5&1/5&1/3\text{ to }7/15&2/7.
\end{array}
\]

For saved certified exact representatives of orders 9--14 the corresponding
LP values are approximately

\[
0.500, 0.244, 0.308, 0.250, 0.256, 0.209.
\]

These finite values are evidence only.  Even if an exact-minimizer shell
bound were proved, the archive's cut-compatible rank example shows that
coordinatewise first-moment balance alone does not imply low response rank
or a composable state.

Verdict: **falsified for arbitrary asymptotic near-minimizers in a thin
shell.**  A shell at least as thick as the optimality excess remains open,
but first-moment eutaxy by itself has no proved compositional payoff.

### Candidate 3: unpeeled spectral flatness

The precise candidate was:

> Every sequence with `Q(A_n)-M_n=o(n^(3/2))` has
> `||A_n||_op=O(sqrt(n))`.

It is false by an even simpler sparse implant.  Start from any exact
minimizer `A_n`, choose `S` of size `k=floor(n^(2/3))`, and change the
principal block on `S` to the all-positive signing.  At most `binom(k,2)`
edges change, so the resulting `B_n` satisfies

\[
Q(B_n)\le M_n+k(k-1)=M_n+O(n^{4/3})
                    =M_n+o(n^{3/2}).
\tag{3.9}
\]

On the other hand, the normalized indicator of `S` has Rayleigh quotient
`k-1`, whence

\[
\boxed{||B_n||_{op}\ge k-1\asymp n^{2/3}.}
\tag{3.10}
\]

Thus neither a fixed `O(sqrt(n))` operator bound nor unpeeled spectral
flatness is forced even when the normalized optimality excess tends to zero.
The surviving statement (2.11) correctly localizes the obstruction: the
implanted eigenvector lives on only `n^(2/3)` coordinates.

Verdict: **falsified.**  Any viable spectral theorem must permit peeling or
charge localization quantitatively; it cannot be a uniform equality-case
stability theorem.

## 4. Comparison with the active state and archived counterexamples

The independent pass meets the existing archive at four precise points.

1. Equation (2.2) is the near-minimal extension of the archived exact
   edge-flip/deep-hole certificate.  The archive already develops its
   hypergeometric weighted-layer consequences and shows that local stability
   through four edges need not be global.  Candidate 1 therefore changes no
   frontier arrow.
2. The soft-selection argument (3.6)--(3.8) is the direct absolute-cap
   analogue of the archived common Gibbs cap law for centered width.  That
   archive also supplies the decisive limitation: a common law with small
   coordinate means can still have cut-compatible affine rank merely
   `O(sqrt(n))`.  Balance is not synchronization.
3. The active state records universal-vertex and conference-deletion
   obstructions to unpeeled spectral and square-field rigidity.  The clique
   implant (3.9)--(3.10) is a shorter generic obstruction applying to a
   certified `o(n^(3/2))` halo around every exact minimizer.
4. Equation (2.8) agrees with the principal-deletion program: normalized
   near-minimality is hereditary over `o(sqrt(n))` deleted vertices, but that
   window is far too short to compare fixed-ratio orders.  The exact
   deletion/insertion profile, not scalar heredity, remains the missing
   information.

No external equality/stability theorem was imported.  Conference-matrix
equality theorems require `A^2=(n-1)I`, while (3.9) shows that near-minimality
does not even force the operator scale on the unpeeled matrix.  Lattice
eutaxy theorems concern continuous local maxima with a common subgradient;
the discrete Hamming deep-hole condition supplies only the set-dependent
witness (2.3).  Their hypotheses therefore do not apply to Candidate 2.

## 5. What is actually forced, and what is not

The strongest clean deterministic conclusion from naked near-minimality is
the following hierarchy:

\[
\text{near-minimality}
\Longrightarrow
\begin{cases}
\text{set-dependent near-cap flip witnesses, (2.3);}\\
\text{principal heredity for }o(\sqrt n)\text{ deletions, (2.8); }\\
\text{first-moment local-field control, (2.9); }\\
\text{flat-eigenvector control, (2.11).}
\end{cases}
\]

None gives a common collective response state.  Moreover, the geodesic
planting theorem (3.3)--(3.5) is a Level-5 negative result: it constructs
genuine certified near-minimizers, at every order admitting a minimizer,
whose entire thin near-cap face carries an arbitrarily large frozen edge
set.  Therefore any proposed implication

\[
Q(A)\le M_n+\eta
\Longrightarrow \text{collective shell regularity}
\]

must allow a response window at least of order `eta`, or explicitly quotient
the planted radial motion.  Sparse-edit robustness is not optional.

The one plausible unpromoted boundary question left by this report is:

\[
\text{Does an exact minimizer, or an }\eta\text{-near minimizer on the shell }
\Delta\ge\eta+O(n),
\text{ admit a genuinely collective law stronger than first moments?}
\]

The phrase "stronger than first moments" is essential.  Proving only
`tau_Delta^+=o(1)` would reproduce a known common-balance module without
supplying the rank, overlap linkage, or response-diameter control required
for composition.

## 6. Frontier accounting for this specialist branch

| Item | Classification | Benchmark | Consequence |
|---|---|---:|---|
| Mesoscopic flip-shell coverage | PROVES a weak implication, but archive-rediscovery | 5 | exact quantitative consequence for all near-minimizers; no common law |
| Geodesic planted-face theorem | FALSIFIES thin-shell eutaxy | 5 | shell width must be at least the optimality excess |
| Soft-selected eutaxy (3.6)--(3.8) | PROVES a selected-sequence statement; NO FRONTIER CHANGE | 5 | balanced common law exists for some `o(n^(3/2))`-near-minimizer, not every one |
| Clique spectral implant | FALSIFIES unpeeled spectral rigidity | 5 | peeling/localization is necessary |
| Principal heredity (2.8) | WEAKENS a local-order hypothesis only | 5 | all `o(sqrt(n))` deletions remain near-minimal |

The smallest missing lemma was not reduced to a proved compositional
statement.  It was narrowed negatively: a viable lemma must be both
`eta`-thick and genuinely collective beyond coordinatewise barycentric
balance.  The response loss, state complexity, and orders covered by the
live convergence routes are unchanged.
