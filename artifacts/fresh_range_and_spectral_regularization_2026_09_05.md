# The stronger range theorem and what it does, and does not, regularize

Date: 2026-09-05. Independent derivation of the parent's range consequence.
The diagonal-factorization argument in Section 4 is classical in substance;
it is included with a proof to separate it from the new range information.

Write

\[
 P(A)=\max_x H_A(x),\quad N(A)=-\min_xH_A(x),\quad
 R(A)=P(A)+N(A),\quad Q(A)=\max(P(A),N(A)).
\]

All spins are Boolean and `H_A(x)=x^TAx/2`. Since its uniform mean is zero,
both `P,N` are nonnegative. In particular `Q≤R≤2Q`.

## 1. The paired construction bounds range, not just absolute optimum

Let `B=A/sqrt(n-1)`. The hierarchical construction gives conditional spin
means

\[
 \mu^+=F(U)+S H(U),\qquad \mu^-=-F(U)+S H(U),
\]

with `|F|+|H|≤1`, followed by coordinatewise independent Boolean rounding.
Because `B` is hollow, rounding preserves each conditional quadratic
expectation. Exactly,

\[
 E H_B(\mu^+)-E H_B(\mu^-)=2 E F(U)^TB[S H(U)].       \tag{1}
\]

Thus every proved hierarchical certificate `J` gives

\[
 \boxed{\qquad R(A_n)\ge(2J-o(1))n^{3/2}.\qquad}     \tag{2}
\]

Initially the tree theorem applies to sequences with `Q(A_n)=O(n^{3/2})`.
This entails no restriction on the universal range conclusion: a putative
sequence violating (2) has bounded range and therefore bounded `Q`, so it
lies within the theorem's hypotheses. A subsequence contradiction also
shows the error in the infimum-over-all-signings formulation is uniform.

The completed exact minorant allows any fixed `J` strictly below

\[
 0.429786450737628041355434677764144259463671564581969931627991.
\]

This is the lower endpoint of an independently audited interval, not a
floating-point optimizer. The claim (2) uses any strictly smaller fixed
constant if one prefers to avoid endpoint notation.

If `Q(A_n)≤(C+o(1))n^{3/2}`, then (2) implies separately

\[
 \boxed{\quad
 P(A_n),N(A_n)\ge(2J-C-o(1))n^{3/2}.
 \quad}                                                    \tag{3}
\]

In particular actual asymptotic near-minimizers, using the all-order upper
bound `C=1/2`, have both positive and negative extrema at least
`(2J-1/2-o(1))n^{3/2}> (0.35957290147-o(1))n^{3/2}`.

## 2. Exact superadditivity under vertex partitions

Partition the vertices into blocks `D,R`, keeping all original bridges.
Choose spins maximizing both internal energies. Reversing all spins of
one block preserves both internal energies and negates the bridge term.
One of the two choices therefore gives at least their sum. Applying the
same argument to minima proves the exact inequalities

\[
 P(A)\ge P(A_D)+P(A_R),\qquad
 N(A)\ge N(A_D)+N(A_R).                              \tag{4}
\]

There is no signed-cancellation assumption. For multiple blocks, one may
average over independent whole-block reversals: every bridge has mean
zero, while all chosen internal energies remain unchanged. Consequently

\[
 P(A)\ge\sum_bP(A_b),\quad N(A)\ge\sum_bN(A_b),
 \quad R(A)\ge\sum_bR(A_b).                          \tag{5}
\]

For any induced block of order `k=o(n)`, (2) applied to its complement
and (4) give

\[
 \boxed{\quad
 R(A_D)\le 2Q(A)-2J n^{3/2}+o(n^{3/2}).
 \quad}                                                    \tag{6}
\]

This is stronger than the corresponding bound on `Q(A_D)`. The error is
uniform over the choice of block, since the complement's range lower
bound is uniform. More generally replace the final `n^{3/2}` in (6) by
`(n-k)^{3/2}` without requiring `k=o(n)`.

At a near-minimizer with `Q≤(1/2+o(1))n^{3/2}`, every sublinear block has
range budget at most

\[
 (1-2J+o(1))n^{3/2}< (0.140427098525+o(1))n^{3/2}.     \tag{7}
\]

The sum of ranges of any disjoint family whose union is sublinear obeys
the same budget.

## 3. Localized spikes are bounded, but not eliminated

A coherent all-positive block of order `k` has
`P=k(k-1)/2` and `N=floor(k/2)`. Hence if `k=gamma n^{3/4}+o(n^{3/4})`,
(7) forces

\[
 \gamma\le\sqrt{2-4J}<0.5300.                        \tag{8}
\]

This improves a constant in the allowable spike budget; it does not force
`gamma=0`, nor exclude unbounded normalized operator norms.

For any symmetric hollow sign matrix one also has

\[
 \|A\|_{op}^2\le\beta(A)\le2R(A),\qquad
 \beta(A)=\max_{x,y\in\{\pm1\}^n}|x^TAy|.            \tag{9}
\]

The first inequality follows by bounding the largest absolute row sum of
`A²` using a row of `A` as a cube vector; here it uses flat sign entries.
For the second, put `u=(x+y)/2`, `v=(x-y)/2`; these are disjoint-support
cube vectors and `x^TAy=2(H_A(u)-H_A(v))`. Completing each partial spin at
random shows `P≥H_A(u)` and `N≥-H_A(v)`, and reversing the roles handles
the opposite sign. Thus (6) also limits the operator norm of an induced
sublinear sign block, but only at scale `n^{3/4}`.

There is no proved replacement argument converting (6) into a vanishing
budget. Replacing a block changes the internal energy, but its optimizing
spins can depend on the bridges. The exact bridge-reversal inequalities
are lower bounds on extrema, not upper bounds that survive arbitrary
replacement. In particular they cannot justify assuming that a favorable
internal energy change lowers the whole matrix's optimum.

## 4. Classical fixed-fraction spectral regularization, proved directly

Every low-cap signing, without any minimizing hypothesis, already admits
bounded-operator regularization by removing an arbitrarily small *fixed*
fraction of its vertices. We give constants and the scope precisely.

Put `c=asinh(1)` and `K=pi/(2c)<1.783`. The elementary real Grothendieck
inequality with this constant is

\[
 \left|\sum_{ij}a_{ij}\langle u_i,v_j\rangle\right|
 \le K\beta(A)                                           \tag{10}
\]

for arbitrary families of unit vectors. To verify it, replace the vectors
by orthogonal direct sums of odd tensor powers with weights from
`sin(c t)`. Assign the alternating coefficient signs to the second family.
The sum of absolute coefficient weights is `sinh(c)=1`, so the new vectors
are unit and their cross inner products equal `sin(c<u_i,v_j>)`. Common
Gaussian hyperplane rounding has expected sign product
`(2/pi)arcsin(sin(c<u_i,v_j>))=(2c/pi)<u_i,v_j>`. The maximum bilinear
Boolean value bounds the absolute expected value, proving (10). The
Gaussian sign-correlation formula follows directly from the uniform angle
of a two-dimensional standard Gaussian.

There exists a nonnegative diagonal matrix `D` such that

\[
 D\succeq A,\qquad D\succeq-A,\qquad
 \operatorname{Tr}D\le K\beta(A).                         \tag{11}
\]

Indeed the dual of minimizing `Tr D` subject to the two matrix inequalities
is

\[
 \max\{\operatorname{Tr}A(X-Y):X,Y\succeq0,
                              \operatorname{diag}(X+Y)=1\}.
\]

Write the two Gram systems as `X_ij=<a_i,a_j>` and
`Y_ij=<b_i,b_j>`. The vectors `u_i=(a_i,b_i)` and `v_i=(a_i,-b_i)` are unit,
and their cross products equal `X_ij-Y_ij`. Equation (10) bounds the dual
by `K beta(A)`. Strict feasibility of a large scalar diagonal gives strong
duality and attainment. Nonnegativity of `D` follows from the hollow
diagonal of `A` and either semidefinite constraint.

Discard coordinates with

\[
 D_{ii}>\frac{K\beta(A)}{\epsilon n}.
\]

There are fewer than `epsilon n` such coordinates. On the remaining
principal submatrix `A_R`, (11) implies

\[
 \boxed{\quad
 \|A_R\|_{op}\le\frac{K\beta(A)}{\epsilon n}
 \le\frac{2K R(A)}{\epsilon n}
 \le\frac{4K Q(A)}{\epsilon n}.
 \quad}                                                     \tag{12}
\]

In particular if `Q(A)≤C n^{3/2}`, then
`||A_R||op≤(4KC/epsilon)sqrt(n)`. Relative to its own order, this is a
fixed operator bound depending only on `C,epsilon`.

No claim of novelty is made for this Grothendieck/Pietsch-type argument.
It yields bounded operator norm after `epsilon n` deletions for each fixed
epsilon. Taking `epsilon→0` makes the proved operator bound diverge; it
does not yield an `o(n)` deletion theorem with a fixed operator bound.
It also imposes no Gram-coherence condition and does not remove twin-row
or fixed local tensor correlations.

## 5. A useful transfer principle for future lower bounds

Suppose a constant `c_*` is a universal asymptotic lower bound for `Q` on
every uniformly bounded-normalized-operator family, with the *same*
constant independent of that family's operator bound. Then it is a
universal lower bound on all signings. To see this, any counterexample has
bounded `Q`; use (12) with fixed epsilon, then induced-submatrix monotonicity
`Q(A)≥Q(A_R)` to obtain

\[
 Q(A)\ge(c_*-o(1))(1-\epsilon)^{3/2}n^{3/2}.
\]

Let epsilon tend to zero after the dimension limit. The same argument
works for range using (4). This can legitimately reduce a future universal
lower-bound proof to bounded operator norm; it cannot add independence,
coherence, or an exact involution identity to the remaining matrix.

## 6. Full-order bounded-operator approximation of the minimum

There is a further variational consequence: the deleted vertices can be
reintroduced with a controlled error in the Boolean optimum. This uses no
new range inequality, but it does give a genuine reduction of the original
convergence problem.

Suppose `Q(A)≤C n^{3/2}`, and choose `0<epsilon≤1/2`. Use Section 4 to
remove a set `D` of order `k<epsilon n`, leaving `R`. Keep the induced
matrix on `R`, but independently resample every undirected edge incident
to `D` as a uniform sign. Call the resulting full-order signing `A_tilde`.
For every spin, its energy is the sum of the old `R` energy, the new bridge
bilinear energy, and the new `D` energy. Consequently, deterministically,

\[
 Q(\widetilde A)\le Q(A_R)+\beta(C_{R,D})+Q(C_{D,D}).   \tag{13}
\]

Here `beta(C_R,D)` is the rectangular maximum over the two independent
spin vectors. The basic induced-submatrix monotonicity `Q(A_R)≤Q(A)` is
exact, by averaging uniformly over the missing coordinates.

Hoeffding bounds and union bounds show, simultaneously with probability
tending to one,

\[
 \begin{split}
 \beta(C_{R,D})
 &\le(\sqrt{2\log2}+o(1))\sqrt{k(n-k)n},\\
 Q(C_{D,D})
 &\le\sqrt{\log2}\,k^{3/2}+o(n^{3/2}).
 \end{split}                                               \tag{14}
\]

The errors can be made uniform in `k≤epsilon n`: use the explicit
thresholds

\[
 \sqrt{2k(n-k)(n\log2+2\log n)},\qquad
 \sqrt{k(k-1)(k\log2+2\log n)}.
\]

The corresponding failure probabilities are at most `2n^{-2}` each.
The conventions for empty edge sets are harmless.

The random added symmetric matrix `E`, zero on `R×R`, also satisfies
`||E||op≤8sqrt(n)` with probability tending to one. For completeness,
for every fixed unit vector `x`,

\[
 E e^{\lambda x^TEx}\le e^{\lambda^2},\qquad
 \Pr\{|x^TEx|\ge t\}\le2e^{-t^2/4}.
\]

A one-quarter net of the Euclidean unit sphere has size at most `9^n`.
The standard quadratic-net inequality is
`||E||op≤2 max_net |x^TEx|`. Put `t=4sqrt(n)` and take a union bound;
the failure probability is at most `2 exp[-(4-log9)n]`.

Thus a deterministic completion exists for which

\[
 \boxed{\begin{split}
 Q(\widetilde A)&\le Q(A)+[f(\epsilon)+o(1)]n^{3/2},\\
 \|\widetilde A\|_{op}&\le
            (4KC/\epsilon+8)\sqrt n,\\
 f(\epsilon)&=\sqrt{2\log2\,\epsilon(1-\epsilon)}
                         +\sqrt{\log2}\,\epsilon^{3/2}.
 \end{split}}                                               \tag{15}
\]

For the first line, the function of `k/n` in (14) is increasing on
`[0,epsilon]` when `epsilon≤1/2`. The second line combines the old
principal-submatrix bound (12) with the random-added-matrix bound. In
particular the error is `O(sqrt(epsilon))`, at the same matrix order.

Define the restricted minima

\[
 M_n(L)=\min\{Q(A): A\text{ is a hollow sign matrix},
                         \ \|A\|_{op}\le L\sqrt n\}.
\]

Use any fixed `C` strictly above the known all-order asymptotic upper
bound on `M_n/n^{3/2}`. For all sufficiently large fixed `L`, take
`epsilon=4KC/(L-8)≤1/2` in (15), starting from an exact minimizer. Then

\[
 \boxed{\quad
 0\le\limsup_n\frac{M_n(L)-M_n}{n^{3/2}}
 \le f\left(\frac{4KC}{L-8}\right)=O(L^{-1/2}).
 \quad}                                                     \tag{16}
\]

The construction ensures these restricted feasible sets are nonempty
for large `n`; no existence assumption has been hidden in the notation.

**Conditional convergence reduction.** If, for every sufficiently large
fixed `L`, the sequence `M_n(L)/n^{3/2}` converges, then the original
sequence `M_n/n^{3/2}` converges. Indeed (16) bounds its limsup-minus-liminf
by `f(4KC/(L-8))`, and this tends to zero. The same conclusion follows if
the oscillations of the restricted sequences themselves tend to zero as
`L→infinity`.

This does not prove convergence for a fixed operator cap. It isolates a
bounded-operator version of the original problem and quantifies that no
additional asymptotic obstruction comes solely from unbounded spectral
spikes. It also does not assert that the modified matrix is close to the
original one in normalized Frobenius norm: the construction deliberately
resamples every edge incident to the discarded vertices.
