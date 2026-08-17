# Spectral/harmonic audit of quadratic-cap near-minimizers

Date: 2026-08-17.

Status: specialist report for the near-minimizer campaign.  The two outputs
below are a proved edit-robust spectral envelope and a proved scalable
falsifier for stronger peeling.  The only imported ingredient in the first
output is the real symmetric Grothendieck factorization, cited precisely
below.  No claim is made that either output characterizes minimizers.  No
canonical theorem or frontier file was edited.

## 0. Verdict

Let

\[
 H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
 Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,
 \qquad M_n=\min_AQ(A),
\]

where `A` is hollow, symmetric, and has signs off the diagonal.  The audit
produced exactly two statements.

| Output | Evidentiary status | Edit robustness | Response/composition payoff |
|---|---|---|---|
| SH.1, Grothendieck spectral envelope | **PROVED**, using a standard cited factorization | normalized trace norm changes by `o(1)` under `o(n^(3/2))` flips; the peeling conclusion re-applies after every such edit | **none at the target scale**: it gives sublinear rank only above a diverging multiple of `sqrt(n)`, whose uniform Boolean tail error is too large |
| SH.2, multi-clique peeling obstruction | **PROVED**, elementary construction around every exact minimizer | it lies inside the mandated `o(n^(3/2))` edit halo | **negative**: an `O(sqrt(n))` core can require `n^(1-o(1))` exceptional spins, and the remaining bounded-operator sign core still need not have a sublinear SVD interface |

Thus naked near-minimality does force a useful weighted/peeled spectral
description, but the description does not supply the campaign's missing
collective response state.  The stronger hope

\[
 \text{delete }O(n^{1-\delta})\text{ vertices and obtain }O(\sqrt n)
 \text{ operator norm}
\]

is false for every fixed `delta>0`, even for
`o(n^(3/2))`-near-minimizers.  The purely qualitative possibility of an
`o(n)` peel followed by an `O(sqrt(n))` core is not disproved.

## 1. Output SH.1: an edit-robust Grothendieck spectral envelope

Write

\[
 \beta(A)=\|A\|_{\infty\to1}
 =\max_{x,y\in\{\pm1\}^n}|x^TAy|.
\]

Let `K_G` denote the real Grothendieck constant.

### Theorem SH.1

For every hollow symmetric real matrix `A` on `n` vertices,

\[
 \boxed{\beta(A)\le 4Q(A).}                                      \tag{SH.1}
\]

There are a nonnegative diagonal matrix `D` and a symmetric matrix `T`
such that

\[
 A=DTD,\qquad \operatorname{tr}(D^2)=1,
 \qquad \|T\|_{2\to2}\le 4K_GQ(A).                              \tag{SH.2}
\]

Consequently:

1. the trace norm satisfies

   \[
   \boxed{\|A\|_*\le4K_GQ(A);}                                  \tag{SH.3}
   \]

2. for every `L>=1`, there is a vertex set `S_L` with

   \[
   \boxed{|S_L|<\frac nL,\qquad
   \|A[V\setminus S_L]\|_{2\to2}
   \le {4K_GLQ(A)\over n};}                                     \tag{SH.4}
   \]

3. for every `t>0`, the above-root spectral tail obeys

   \[
   \boxed{\#\{j:|\lambda_j(A)|>t\sqrt n\}
   \le {4K_GQ(A)\over t\sqrt n}.}                               \tag{SH.5}
   \]

The theorem uses only `Q(A)`, not a maximizing spin, the value of `M_n`, or
the response landscape.

### Proof

First note the hereditary inequality `Q(A[U])<=Q(A)`: fix spins on `U`,
extend them by independent unbiased spins off `U`, and average the full
energy.  For signs `x,y`, put

\[
 z_i=x_i\mathbf1_{\{x_i=y_i\}},\qquad
 w_i=x_i\mathbf1_{\{x_i=-y_i\}}.
\]

Then `x=z+w`, `y=z-w`, the supports of `z,w` are disjoint, and symmetry gives

\[
 x^TAy=z^TAz-w^TAw.
\]

Each term has absolute value at most `2Q(A)` by heredity, proving (SH.1).

The symmetric Grothendieck factorization says that a symmetric matrix `G`
admits `G=DTD`, with `D>=0` diagonal, `tr(D^2)=1`, and

\[
 \|T\|_{2\to2}\le K_G\|G\|_{\infty\to1}.
\]

This is Theorem 5.3 in Tropp's arXiv/technical-report version (Theorem 5.2
in the published SODA version); Section 5.3 of the published version gives
the semidefinite construction.  See Joel A. Tropp,
[*Column Subset Selection, Matrix Factorization, and Eigenvalue
Optimization*](https://arxiv.org/abs/0806.4404), SODA 2009, pp. 978--986.
Together with (SH.1), this proves (SH.2).

Schatten Holder gives

\[
 \|DTD\|_*
 \le\|D\|_F\,\|TD\|_F
 \le\|D\|_F^2\|T\|_{2\to2}
 =\|T\|_{2\to2},
\]

which proves (SH.3).  If `D=diag(d_i)`, take

\[
 S_L=\{i:d_i^2>L/n\}.
\]

Since `sum_i d_i^2=1`, `|S_L|<n/L`.  On the complement,
`||D||_op^2<=L/n`, so compressing (SH.2) proves (SH.4).  Finally, every
eigenvalue counted in (SH.5) contributes more than `t sqrt(n)` to the trace
norm, and (SH.3) finishes the proof.  `square`

For the underlying sign/SDP form of Grothendieck's inequality, see Noga Alon
and Assaf Naor,
[*Approximating the Cut-Norm via Grothendieck's
Inequality*](https://doi.org/10.1137/S0097539704441629), SIAM J. Comput. 35
(2006), 787--803.  This citation is supporting literature, not an additional
unproved hypothesis.

### Explicit near-minimizer corollary

A union bound supplies a completely explicit universal upper bound for
`M_n`.  Choose the off-diagonal signs independently.  For a fixed Boolean
`x`, Hoeffding gives

\[
 \Pr\{|H_A(x)|\ge t\}\le2e^{-t^2/(2N)},\qquad N=\binom n2.
\]

Union over all `2^n` spins, with

\[
 t=\sqrt{2N(n+2)\log2},
\]

shows that

\[
 M_n\le\sqrt{n(n-1)(n+2)\log2}
      =(\sqrt{\log2}+o(1))n^{3/2}.                               \tag{SH.6}
\]

Therefore every `A` satisfying

\[
 Q(A)\le M_n+\varepsilon n^{3/2}
\]

obeys

\[
 \frac{\|A\|_*}{n^{3/2}}
 \le4K_G(\sqrt{\log2}+\varepsilon+o(1)),                         \tag{SH.7}
\]

and, for every `L>=1`, has a set of fewer than `n/L` vertices outside
which

\[
 \|A[V\setminus S_L]\|_{2\to2}
 \le4K_G(\sqrt{\log2}+\varepsilon+o(1))L\sqrt n.                 \tag{SH.8}
\]

For `L=L_n->infinity`, (SH.5) also says that only `o(n)` eigenmodes can lie
above `L_n sqrt(n)`.

### Exact robustness under the edit halo

If `A` and `B` differ on `r` unordered edges, their difference is a sum of
`r` matrices of the form

\[
 \pm2(e_ie_j^T+e_je_i^T),
\]

each of trace norm `4`.  Hence

\[
 \boxed{|\|A\|_*-\|B\|_*|\le4r.}                                \tag{SH.9}
\]

Thus `||A||_*/n^(3/2)` is genuinely invariant under the full
`o(n^(3/2))` Hamming halo.  Also `Q(B)<=Q(A)+2r`, so (SH.4) re-applies to
`B` with its right side enlarged by at most `8K_G Lr/n`.  The exceptional
set need not be the same; the forced property is existential, not a stable
vertex labeling.

### Why SH.1 is weaker than optimizing `Q`

SH.1 holds for every matrix with `Q(A)=O(n^(3/2))`, including ordinary
random sign matrices.  It neither distinguishes exact minimizers from a
constant-factor bounded-cap class nor identifies an active Boolean state.
Conversely, (SH.7) retains only one scalar spectral mass bound, and (SH.8)
retains one weighted factorization/peel; neither determines the `2^n`
Boolean energy values.  It is therefore strictly less informative than the
optimization problem in the operational sense required by the campaign.

### Response/composition assessment: no target-scale payoff

The repository's scale-sensitive spectral bridge theorem gives uniform
Boolean error `n sigma_(r+1)` at balanced size.  SH.1 makes the number of
modes above `L sqrt(n)` sublinear only when `L->infinity`, but discarding the
rest then has certified error `L n^(3/2)`, larger than the target scale.  At
the useful cutoff `epsilon sqrt(n)`, (SH.5) gives no sublinear rank bound.

The peel in (SH.8) is also insufficient.  Taking `L->infinity` stores only
`o(n)` exceptional spins, but leaves a dense sign core with operator norm
`O(L sqrt(n))`, not a low-rank or synchronized core.  Dense bounded-operator
sign matrices have extensive target-visible Frobenius mass, as already
proved in `drafts/bounded_operator_rank_barrier.md`.  Therefore SH.1 supplies
no current response roof, amalgamation law, or `o(n^(3/2))` replacement.

Classification: **PROVES a Level-5 implication, but NO FRONTIER CHANGE.**

## 2. Output SH.2: a scalable multi-clique obstruction to strong peeling

The single `n^(2/3)` all-plus clique from the independent deterministic
report is only the first point on a sharper tradeoff.

### Theorem SH.2 (multi-clique peeling lower bound)

Let `A` be any exact order-`n` minimizer.  Let `S_1,...,S_m` be disjoint
sets of size `k`, where `mk<=n`.  Form `B` by overwriting every off-diagonal
entry inside each `S_j` by `+1`.  Then

\[
 \boxed{Q(B)\le M_n+m k(k-1).}                                  \tag{SH.10}
\]

Moreover, for every vertex set `T` and every `R>=0`,

\[
 \|B[V\setminus T]\|_{2\to2}\le R
 \quad\Longrightarrow\quad
 \boxed{|T|\ge m(k-1-R)_+.}                                    \tag{SH.11}
\]

In particular, let `L=L_n->infinity` with `L=o(n^(1/6))`, and take

\[
 k=\lfloor L\sqrt n\rfloor,
 \qquad m=\left\lfloor{\sqrt n\over L^3}\right\rfloor.
\]

Then

\[
 Q(B)\le M_n+(1+o(1)){n^{3/2}\over L},                           \tag{SH.12}
\]

while, for every fixed `C`, every `T` for which
`||B[V\setminus T]||_op<=C sqrt(n)` satisfies

\[
 \boxed{|T|\ge(1-o(1)){n\over L^2}.}                            \tag{SH.13}
\]

Taking, for example, `L=log log n`, the required peel is `n^(1-o(1))`.
Consequently no fixed `delta>0` admits a universal theorem saying that every
`o(n^(3/2))`-near-minimizer becomes `O(sqrt(n))` in operator norm after
deleting `O(n^(1-delta))` vertices.

### Proof

At most `m binom(k,2)` edges are changed.  Flipping one edge changes every
Boolean energy by at most `2`, hence changes `Q` by at most `2`.  This proves
(SH.10).

For each `j`, put `s_j=|S_j\setminus T|`.  If `s_j>=1`, the principal block
on `S_j\setminus T` is `J_(s_j)-I_(s_j)` and has operator norm at least its
largest eigenvalue `s_j-1`.  A principal compression cannot have larger
operator norm than the full symmetric matrix, so

\[
 s_j-1\le\|B[V\setminus T]\|_{2\to2}\le R.
\]

If `s_j=0`, all `k` vertices of that block were deleted.  Thus
`|T intersect S_j|>=(k-1-R)_+` in either case.
Summing over the disjoint blocks proves (SH.11).

For the displayed choice of `k,m`, the blocks use

\[
 mk=(1+o(1)){n\over L^2}\le n
\]

vertices and edit at most

\[
 m\binom k2=(1+o(1)){n^{3/2}\over2L}=o(n^{3/2})
\]

edges.  Equations (SH.12)--(SH.13) follow from (SH.10)--(SH.11).  Finally,
`n/L^2` exceeds `n^(1-delta)` for every fixed `delta>0` if `L` grows
subpolynomially.  `square`

### Scope and response consequence

SH.2 is a certified construction inside the exact-minimizer edit halo at
every order: it does not assume that the starting minimizer has conference,
Hadamard, random, or bounded-operator structure.  It strengthens the planted
`n^(2/3)` clique warning by showing that the number of vertices needed for a
fixed-root-scale peel can be arbitrarily close to linear.

It does **not** disprove a bare `o(n)`-vertex peeling theorem; its lower bound
`n/L^2` is still `o(n)`.  Nor does it show that every possible structural
state is large.  It specifically falsifies polynomially sublinear
fixed-`O(sqrt(n))` spectral peeling and any response architecture that pays
one exact spin per peeled vertex while claiming `O(n^(1-delta))` state bits.

Classification: **FALSIFIES a natural strengthened Level-5 spectral arrow.**
The negative conclusion is compositional: the hoped-for polynomially
sublinear exceptional-spin carrier cannot be inferred from near-minimality.
It does not itself yield a replacement or convergence theorem.

## 3. Literature boundary: what transfers and what does not

The literature search emphasized primary sources and equality/inverse
results.  None located bridges quadratic-cap near-optimality to an
`O(sqrt(n))` unpeeled or cheaply peeled operator theorem.

1. **Foundational graph discrepancy.**  Erdős and Spencer proved the
   `Theta(n^(3/2))` induced-subgraph discrepancy scale in
   [*Imbalances in k-colorations*](https://doi.org/10.1002/net.3230010407),
   Networks 1 (1971), 379--385.  This establishes the ambient order of
   `M_n`-type quantities but contains no inverse or stability description of
   colorings close to minimum discrepancy.

2. **Grothendieck/cut norm.**  Alon--Naor's cited paper gives the SDP/sign
   comparison, and Tropp gives the exact diagonal factorization used in
   SH.1.  These are the appropriate edit-robust harmonic tools because the
   quadratic cap polarizes to the `infinity-to-1` norm.  Their conclusion is
   a weighted factorization, not ordinary operator rigidity.

3. **Inverse expander mixing.**  Bilu and Linial,
   [*Constructing expander graphs by 2-lifts and discrepancy vs. spectral
   gap*](https://doi.org/10.1007/s00493-006-0029-7), Combinatorica 26
   (2006), 495--519, control spectral radius from a *size-normalized*
   disjoint-set discrepancy parameter.  From the cap one only gets, for
   disjoint `S,T`,

   \[
   {|\mathbf1_S^TA\mathbf1_T|\over\sqrt{|S||T|}}
   \le\min\left\{\sqrt{|S||T|},{4Q(A)\over\sqrt{|S||T|}}\right\}
   \le2\sqrt{Q(A)}=O(n^{3/4}).
   \]

   Their inverse lemma therefore gives at best an `n^(3/4)`-scale input
   (up to logarithms), fully compatible with the clique implants.  The
   needed local size-sensitive discrepancy hypothesis is not forced by
   `Q(A)=O(n^(3/2))`.  Nikiforov's
   [*Cut-norms and spectra of matrices*](https://arxiv.org/abs/0912.0336)
   similarly does not remove localized spikes from a global cap bound.

4. **Conference/Seidel equality.**  Haemers and Parsaei Majd,
   [*Spectral symmetry in conference
   matrices*](https://arxiv.org/abs/2004.05829), record that a symmetric
   conference matrix has spectrum `+-sqrt(n-1)` and analyze principal
   submatrices.  The 2026 preprint of Guterman and Saha,
   [*Extremal problems on the p-Seidel energy of
   graphs*](https://arxiv.org/abs/2606.17828), proves for `p>2` that the
   minimum `p`-Seidel energy is attained exactly when
   `A^2=(n-1)I` (when such a conference matrix exists).  This equality is a
   power-mean consequence of the fixed identity
   `sum lambda_i^2=n(n-1)`.  It is not an equality theorem for `Q`, and SH.2
   shows why importing its stability conclusion would be invalid: cap-near
   matrices may contain many mesoscopic all-plus blocks.

5. **Relevant 2024--2026 discrepancy work.**  Maillard,
   [*Average-case matrix discrepancy: satisfiability
   bounds*](https://doi.org/10.1002/rsa.70033), Random Structures &
   Algorithms (2025), proves phase transitions for choosing signs on
   independent Gaussian matrix summands under an operator-norm constraint.
   The randomness, the variable being signed, and the objective differ from
   the present deterministic Seidel problem, so no stability implication
   transfers.  Christoph, Gishboliner, and Krivelevich,
   [*Subgraph discrepancies in the complete
   graph*](https://arxiv.org/abs/2602.04069) (2026), obtain sharp discrepancy
   guarantees for copies of prescribed guest graphs.  Their quantifier is
   over embeddings of a guest graph, not Boolean switchings of one signing,
   and again gives no near-equality structure for `Q`.

These sources support a negative literature conclusion, not a priority
claim: current spectral equality theorems optimize spectral moments or
assume local normalized mixing, whereas the campaign objective optimizes a
global Boolean quadratic cap and permits mesoscopic localization.

## 4. Frontier accounting

The strongest rigorously forced statement from this branch is now

\[
 \text{near-minimal cap}
 \Longrightarrow
 \begin{cases}
 \text{edit-stable }O(n^{3/2})\text{ trace norm},\\
 \text{after }<n/L\text{ vertices, operator norm }O(L\sqrt n),\\
 o(n)\text{ modes above }L_n\sqrt n\text{ for }L_n\to\infty.
 \end{cases}
\]

This is demonstrably weaker than `Q` optimization, but no proved
response/composition theorem consumes it at `o(n^(3/2))` loss.  SH.2 also
shows that asking for a fixed `O(sqrt(n))` core with a polynomially
sublinear peel is impossible.  Hence the spectral version of the smallest
missing lemma is narrowed to something genuinely collective:

> Either prove an `o(n)`-peel `O(sqrt(n))` theorem together with a response
> law that does not pay one unconstrained bit per exceptional vertex, or
> find a non-SVD synchronization property of the weighted Grothendieck core.

Even this narrowed lemma has no established composition payoff.  The
appropriate campaign label is therefore:

- SH.1: **PROVES AN ARROW only to a weak spectral envelope; NO FRONTIER
  CHANGE for R1/R3**.
- SH.2: **FALSIFIES polynomially sublinear fixed-root spectral peeling**.
- Benchmark: Level 5 for both theorem and falsifier.
- Assumption distance: 4 (uniform over genuine near-minimizers).
- Target loss/state: unchanged; no `o(n^(3/2))` reusable carrier obtained.
