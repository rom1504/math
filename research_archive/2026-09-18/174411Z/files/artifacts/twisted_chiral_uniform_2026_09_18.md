# Twisted chiral doubling: joint selection and exact convergence obligations

Campaign: 2026-09-18, 17:09:11--20:09:11 UTC. This is the independent
uniform-theorem track. Convention throughout:

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad Q(A)=\max_x|H_A(x)|,
\qquad M_n=\min_A Q(A),\qquad f(n)=M_n/n^{3/2}.
\]

The older `chiral_scale_preserving_lift.md` and
`dependent_lift_analytic_audit.md` use twice this norm. Their compressed-lift
obstructions do **not** automatically apply to the twisted double studied here.
The two-multiplier power-saving convergence implication was already proved in
`fresh_limit_algebra_2026_09_05.md`; section 5 below sharpens its accounting and
records explicit countermodels to weaker claims.

## 1. Two exact forms of the objective

Let B be switching/permutation-equivalent to A, let d have sign entries, put
C=B+diag(d), and set

\[
D=\begin{pmatrix}A&C\\C&-A\end{pmatrix}.
\]

For arbitrary Boolean x,y,

\[
H_D(x,y)=H_A(x)-H_A(y)+x^TCy.
\]

Flipping all coordinates of x leaves its internal energy unchanged and
negates the bridge. Therefore

\[
\boxed{Q(D)=\max_{x,y}\bigl(|H_A(x)-H_A(y)|+|x^TCy|\bigr).} \tag{1}
\]

Equivalently write y=t*x, I={i:t_i=1}, J=I^c, z=x. Then

\[
H_A(x)-H_A(y)=2z_I^TA_{IJ}z_J,
\]

\[
x^TCy=2H_{B[I]}(z_I)-2H_{B[J]}(z_J)+d(I)-d(J).
\]

This proves the proposed cut/internal-energy identity exactly, with the two
terms evaluated at the **same** partition and spin configuration.

Let D_0 denote the same matrix with d=0. The diagonal matching has n edges,
so for every choice of d,

\[
\boxed{|Q(D)-Q(D_0)|\le n.} \tag{2}
\]

Thus diagonal optimization is important at finite order but cannot remove a
fixed positive leading-scale defect.

For the complex representation, put u=(x+y)/2, v=(x-y)/2, w=u+iv. Every
w_i lies in {1,-1,i,-i}, and

\[
\boxed{H_D(x,y)=\operatorname{Re}\bigl(w^T(C-iA)w\bigr).} \tag{3}
\]

The real linear operator D represents the antilinear map
z -> (A+iC) conjugate(z). Consequently

\[
\|D\|_{op}=\|A+iC\|_{op}.
\]

This identity does not replace the four-phase Boolean constraint by an
arbitrary Euclidean vector.

## 2. A genuine quantitative joint-selection lemma

The following sufficient condition involves only two **marginal** profiles of
the seed. It neither assumes a favorable twist nor defines a new name for the
desired minimum.

For 0<=k<=n and nonnegative integers a,b, define

\[
K_k(a)=\#\{(x,y):x,y\in\{\pm1\}^n,
\ d_H(x,y)=k,\ |H_A(x)-H_A(y)|=a\},
\]

\[
J_k(b)=\#\{(u,v):u,v\in\{\pm1\}^n,
\ d_H(u,v)=k,\ |u^TAv|=b\}.
\]

Write J_k(>r)=sum_{b>r}J_k(b). For L>=2Q(A), set

\[
\mathcal E_A(L)=
\sum_{k=0}^n\sum_a
\frac{K_k(a)J_k(>L-a)}{2^n\binom nk}. \tag{4}
\]

**Selection lemma.** If E_A(L)<8, a signed permutation g exists such that,
with B=g^TAg,

\[
\boxed{Q\begin{pmatrix}A&B\\B&-A\end{pmatrix}\le L.} \tag{5}
\]

Every diagonal sign completion then has cap at most L+n. A uniform random
signed permutation succeeds with probability at least 1-E_A(L)/8.

**Proof.** The signed permutation group acts transitively on the ordered
Boolean pairs at a fixed Hamming distance. The orbit at distance k has
2^n binom(n,k) elements. For fixed (x,y), a uniform g therefore makes
(gx,gy) uniform on that orbit, and

\[
\Pr_g\{|H_A(x)-H_A(y)|+|x^Tg^TAg y|>L\}
=\frac{J_k(>L-|H_A(x)-H_A(y)|)}{2^n\binom nk}.
\]

Summing proves that (4) is the exact expected number of violating ordered
pairs, without any independence assumption between different pairs. No pair
with y=+/-x violates when L>=2Q(A). On all other pairs, independent global
reversals of x,y and their interchange form an eight-element symmetry orbit
on which the tested expression is constant. Thus the number of violating
orbits is an integer, with expectation E_A(L)/8. If that expectation is below
one, some g has no violation; Markov gives the probability claim. QED.

The finite group permits conditional-expectation derandomization, although
this statement makes no polynomial-time complexity claim. Permutations enter
materially: random switching alone is transitive only after the entire
coordinatewise product x*y is fixed, rather than just its Hamming weight.

## 3. Exact tests and their limitation

`computations/twisted_chiral_uniform_2026_09_18.py` computes (4) in rational
arithmetic from exact integer profiles. Its separate order-four check averages
the actual violation count over all 384 signed permutations, at 13 thresholds,
and exactly matches (4). Output and all seed matrices are in
`computations/results/twisted_chiral_uniform_2026_09_18.json`.

Selected values:

| seed | Q(A) | smallest L certified by E_A(L)<8 | 2 sqrt(2) Q(A) |
|---|---:|---:|---:|
| stored order 5 minimizer | 4 | 16 | 11.314 |
| stored order 6 minimizer | 5 | 22 | 14.142 |
| stored order 10 minimizer | 13 | 54 | 36.770 |
| dependent chiral order 12 witness | 20 | 76 | 56.569 |

These are **sufficient certificate values**, not minima over twists. For the
order-12 witness the exact expected count at L=56 is

\[
\mathcal E_A(56)=72777149/16128>8.
\]

Thus this first-moment certificate does not establish the target there. It
does not imply that all twists fail: bad configurations may occur in large,
strongly dependent clusters, and favorable twists can exist despite a large
average number of violations. A uniform theorem would require a substantially
stronger profile inequality or a different selection argument.

## 4. What a doubling theorem actually supplies

Suppose an admissible selection rule is available at every step, and

\[
Q(A_{j+1})\le 2^{3/2}Q(A_j)+C n_j,
\qquad n_{j+1}=2n_j.
\]

Then direct summation gives, for every depth r,

\[
\frac{Q(A_r)}{(2^r n)^{3/2}}
\le\frac{Q(A_0)}{n^{3/2}}
+\frac{C(1+\sqrt2)}{2\sqrt n}. \tag{6}
\]

No compression identity is used here. Conversely, if the selection theorem
only holds for exact global minimizers, it gives a recurrence for M_n by
reselecting an actual minimizer at each size; one may not claim it literally
iterates the particular constructed children.

More generally let normalized one-step defects be bounded by epsilon(n), let
e(n)=sup_{m>=n}epsilon(m), and suppose

\[
S(n):=\sum_{j=0}^{\infty}e(2^j n)<\infty,
\qquad S(n)\longrightarrow0. \tag{7}
\]

Every increasing multiplier word whose letters are at least two has total
normalized defect at most S(n). In particular a power saving suffices. A bare
epsilon(n)->0 assertion does not imply (7).

## 5. Coverage: sufficient supplements and explicit countermodels

### Two multiplicatively independent integer multipliers

If, for k=2 and k=3 and all sufficiently large n,

\[
f(kn)\le f(n)+\epsilon(n), \tag{8}
\]

with (7), then f(n) converges. Indeed every descendant n 2^a 3^b has
normalized cap at most f(n)+S(n). The semigroup {2^a 3^b:a,b>=0} has relative
gaps tending to one: a finite set of residues b log 3 modulo log 2 forms an
arbitrarily fine net, and the remaining sufficiently large logarithm is
covered by a nonnegative multiple of log 2. Principal deletion gives
M_m<=M_N for m<=N. Choosing the least descendant N>=m therefore yields

\[
\limsup_{m\to\infty} f(m)\le f(n)+S(n).
\]

Taking n along a liminf subsequence proves the claim. The same proof uses any
two multiplicatively independent integer multipliers, with the corresponding
minimum-generator Dini tail.

### One doubling rule plus a selected proportional-thinning rule

Instead of a second multiplier, it is enough that the constructed dyadic
descendants A_N admit, uniformly for N/2<=m<=N, a principal set S of size m
with

\[
Q(A_N[S])\le (m/N)^{3/2}Q(A_N)+\eta(N)N^{3/2},
\qquad \eta(N)\longrightarrow0. \tag{9}
\]

Only **one** thinning step is used: for each target m choose the first dyadic
descendant N>=m. The normalized thinning error is at most 2^{3/2}eta(N), so
(6), or (7), again implies limsup f<=liminf f. No summability of eta is needed.
Condition (9) is needed only on the chosen descendants, not on all signings
and not on a uniform average of principal subsets. Earlier counterexamples to
average proportional contraction must not be overstated as refuting this
selected-subset statement.

### Why doubling alone is insufficient, even with an O(1) error

Set alpha=3/2 and, for t>=1,

\[
F(t)=t^\alpha\{c+\rho\sin(2\pi\log_2 t)\},
\]

where c>rho>0 and

\[
\rho\sqrt{\alpha^2+(2\pi/\log2)^2}<\alpha c.
\]

Then F is strictly increasing, F(2t)=2^alpha F(t), and F(n)/n^alpha does not
converge. Thus even exact doubling plus monotonicity leaves multiplicative
phase oscillation. One may take c=.46 and rho=.005, keeping this countermodel
inside the currently reported asymptotic interval. Rounding F(n) upward to an
integer with the parity of binom(n,2) changes it by less than two. For all
sufficiently large n the rounded sequence is increasing, has increments
O(sqrt(n))<=n, and satisfies the doubling upper inequality with additive error
two. Hence parity and the elementary one-vertex extension bound do not repair
the coverage gap. This is a countermodel to a proposed implication, **not** a
claim that the actual M_n behaves this way.

### Why two multipliers with only o(n^{3/2}) error are insufficient

For large t let

\[
G(t)=t^\alpha\{c+\rho\sin(\log\log(t+e))\},\qquad 0<\rho<c.
\]

This is eventually increasing and its normalized values do not converge.
For each fixed k, the mean-value theorem gives

\[
G(kn)/(kn)^\alpha-G(n)/n^\alpha=O(1/\log n).
\]

Consequently it satisfies normalized o(1) upper defects simultaneously for
k=2 and k=3. The corresponding dyadic Dini tail diverges. This countermodel
separates the two independent obligations: multiplicatively dense coverage
and summable accumulated error.

## Current conclusion

The exact objective, a computable joint-selection lemma, and the convergence
accounting are proved. No uniform low-cap/minimizer selection theorem has been
proved. In particular a successful O(n)-defect twisted doubling estimate would
be a meaningful step but would still need an order-coverage supplement such as
(8) for a second multiplier or the selected thinning statement (9).
