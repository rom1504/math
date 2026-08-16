# Obstruction atlas for compressed combinatorial energy landscapes

Date: 2026-08-16.  Scope: obstruction report only.  Nothing below proves or
attempts to prove convergence or nonconvergence of `M_n/n^(3/2)`.

## 1. Conventions and evidence labels

For a hollow symmetric signing `A` of order `n`, write

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|.
```

The following labels are used in the atlas.

* **Proved** means that the stated separation is a theorem, possibly using a
  probabilistic existence argument or a displayed exact finite certificate.
* **Computational** means that the claim is an exact finite exhaustive
  enumeration, but no scalable theorem is asserted.
* **Schematic** means that the construction is an abstract energy landscape
  or scalar-pressure model and is not asserted to be realizable by complete
  quadratic signings.

A **true pair** consists of two actual signings (or two actual bridges with
the same fixed children) whose summaries agree exactly or asymptotically and
whose target extrema or composition responses separate.  A **one-sided
blindness construction** exhibits an exceptional Boolean direction invisible
to a summary but does not supply a second signing with the same summary and a
separated cap.  Some finite examples are true matrix pairs but only one-sided
as proposed *ground-state* summaries; this distinction is stated explicitly.

## 2. A common separation template

Let `S_n(A)` be any proposed summary, with metric `d_S`, and let
`F_n(A)` be a normalized extremum or a normalized composition response.  If a
decoder `g_n` has modulus of continuity `omega`, then for every pair `A,B`,

```math
\max\{|F_n(A)-g_n(S_n(A))|,|F_n(B)-g_n(S_n(B))|\}
\ge {1\over2}\left(
 |F_n(A)-F_n(B)|-\omega(d_S(S_n(A),S_n(B)))
\right).                                                   \tag{2.1}
```

This is just the triangle inequality, but it is the precise test needed
here.  An exact collision has `d_S=0`; an asymptotic collision has
`d_S -> 0`.  A fixed positive normalized gap therefore rules out every
uniformly continuous sharp estimator based only on that summary.  The same
inequality applies with

```math
F_n(A)=Q(\mathcal L(A))/|\mathcal L(A)|^{3/2}
```

for a fixed lift or composition operation `\mathcal L`.

The repository's strongest examples are all manifestations of one elementary
zero-entropy-spike identity.  If `B` is obtained from `A` by reversing an
edge set `F`, then

```math
\Delta_F(x):=H_B(x)-H_A(x)
=-2\sum_{ij\in F}a_{ij}x_ix_j,                           \tag{2.2}
```

and for a uniform Boolean spin `X`, orthogonality of distinct edge
characters gives

```math
\mathbb E\Delta_F(X)^2=4|F|.                            \tag{2.3}
```

If the reversed edges are chosen against one selected spin `s`, however,

```math
\Delta_F(s)=2|F|.                                      \tag{2.4}
```

At `|F|=Theta(n^(3/2))`, the normalized endpoint effect in (2.4) is
`Theta(1)`, while the normalized `L^2` effect in (2.3) is
`O(n^(-3/4))`.  Fixed replica laws, bounded moments, fixed signed-subgraph
densities, normalized trace data, and high-temperature pressure can all see
the latter scale while the maximum sees the former.  The Walsh planting in
Section 7 additionally arranges `||A||_op,||B||_op=O(\sqrt n)`, so this is not
merely a large-spectral-spike phenomenon.

Two standard entropy inequalities mark the repair boundary.  For a function
on a set of `L` states,

```math
\|f\|_\infty\le L^{1/p}\|f\|_p,                         \tag{2.5}
```

so a generic moment-to-maximum argument on the projective Boolean cube pays
`2^((n-1)/p)` and needs `p=Omega(n)` for constant-factor control.  For the
normalized cosh pressure

```math
p_{n,\beta}(A)
:={1\over n}\log\left[
 2^{-n}\sum_x\cosh\left({\beta H_A(x)\over\sqrt n}\right)
\right],
```

the antipodal maximizing pair and the trivial upper bound give

```math
{p_{n,\beta}(A)\over\beta}
\le {Q(A)\over n^{3/2}}
\le {p_{n,\beta}(A)\over\beta}+{\log2\over\beta}.        \tag{2.6}
```

Thus bounded inverse temperature leaves a nonvanishing entropy window.  The
actual signing pair in Section 7 shows that this window is not only an
abstract concern, at least on a whole strict-high-temperature interval.

## 3. Atlas at a glance

| ID | Summary declared equal or asymptotically equal | Separated quantity | Evidence | Pair type |
|---|---|---|---|---|
| S | Complete eigenvalue multiset at order 8 | `Q=14` versus `Q=12` | **Computational** exact enumeration | true finite pair |
| L | Every oriented switching/permutation restriction profile through 6 vertices; hence the associated bounded Eulerian, spectral, and energy moments | Hadamard-lift cap gap at least `N^(3/2)/(8*10^(3/2))` | **Proved** from an exact finite seed/certificate | true scalable pair |
| H | Complete projective energy histogram, hence also the one-coset augmented-cut distance enumerator, for two order-8 minimizers | fixed universal-double caps `40` versus `32` | **Computational** exact enumeration | true finite pair with different composition response |
| R | Full signed energy extrema on every magnetization slice, and therefore both radial support functions `F_A^\pm(t\mathbf1)` for all `t` | best one-vertex insertion `12` versus `10` | **Computational** exact enumeration | true finite pair |
| G | Full bridge Gram pair, or even energy-shell-conditioned second moments | parent caps `9` versus `11` (and separate `9` versus `11` collisions) | **Computational** exhaustive bridge enumeration | true finite bridge pairs |
| A | Every fixed action profile at the available `\infty -> 1` endpoint (`d_M -> 0`) | normalized caps separated by a fixed positive constant | **Proved** existential construction | true scalable pair; no common `2 -> 2` bound |
| T | Every fixed trace polynomial, fixed signed-graph density, and every fixed-replica sampled overlap-energy law | `limsup Q(A_n)/n^(3/2)<=1/2`, `liminf Q(B_n)/n^(3/2)>=2/3` | **Proved** Walsh planting | true scalable pair with bounded normalized operator norm |
| P | Normalized cosh pressures agree on an entire interval `|beta|<=beta_0` | the same cap gap of at least `1/6-o(1)` | **Proved**, derived from row T and high-temperature Frobenius stability | true scalable pair |
| O | One selected labelled row-field law `(5,3,3,1,1,1)` at a strict one-flip maximum | caps `11` versus `7` | **Computational** exact enumeration | true matrix pair, but one-sided as a *ground-state-law* claim |
| D | Vanishing fourth-moment/two-walk defect `Delta_4 -> 0` | a Boolean direction with order-one normalized `A^2` defect | **Proved** conference deletion identity | one-sided blindness, not a cap-separated pair |
| X | All pressures on any prescribed bounded beta interval for an abstract landscape | normalized maxima `0` versus `delta` | **Schematic** | abstract pair, not asserted quadratic-sign realizable |

Rows L, A, T, and P are the theorem-level separations.  Rows S, H, R, G,
and O are exact finite falsifiers: they reject universal finite-state claims,
but by themselves say nothing asymptotic.  Row D is useful only as a
one-sided stress test.

## 4. Exact finite collision cards

### 4.1 Same full spectrum, different Boolean extrema (row S)

The root-gauged order-8 signings with internal-edge masks `6875` and `6887`
have identical power traces

```text
(tr A^k)_(k=1)^8 = (0,56,0,680,0,10328,0,169544)
```

and hence the same characteristic polynomial

```math
\lambda^8-28\lambda^6+222\lambda^4-620\lambda^2+425.     \tag{4.1}
```

Their exact energy ranges are respectively `[-14,12]` and `[-12,12]`, so

```math
Q(A)=14,\qquad Q(B)=12.                                  \tag{4.2}
```

This is equality of the complete eigenvalue multiset, not merely the first
few moments.  The claim is finite and computational: all `2^7` projective
spins are evaluated exactly, and the characteristic coefficients and power
traces are integer checks.  Reproduction is in
`computations/audit_constructive_family_obstructions.py`; the exact matrices
and hashes are in
`computations/results/constructive_family_phase2_audit.json`.

The same pair also gives a composition warning.  Its zero-diagonal
polynomial fractional doubles have exact caps `40` and `44`, respectively,
despite having the same child spectrum.  This second statement concerns a
fractional bridge with cross diagonal zeros, whereas (4.2) is already a
separation between genuine signings.

### 4.2 Same full energy histogram, different lift response (row H)

The exhaustive order-8 minimizer classification has two switching/
permutation/global-sign classes.  Both have `Q=10` and exactly the same
projective energy histogram

```text
-10:4, -8:10, -6:12, -4:16, -2:16, 0:12,
  2:16,  4:16,  6:12,  8:10, 10:4.                     (4.3)
```

For any hollow signing `S`, define the fixed universal double

```math
\mathcal D(S)=
\begin{pmatrix}
S&S+I\\
S+I&-S
\end{pmatrix}.                                          \tag{4.4}
```

It is a hollow order-16 signing.  Applying (4.4) to the two classes gives

```math
Q(\mathcal D(S_0))=40,
\qquad
Q(\mathcal D(S_1))=32.                                  \tag{4.5}
```

Thus the complete unlabeled one-body landscape determines the child cap but
does not determine even this one fixed composition response.  This is a true
pair of exact minimizer classes, not a comparison with a nonrealizable
summary.  The source matrices have hashes
`ce853c77df6700ec17d43c5cba7cef88bd7aa2b91b470283854fa368661a42a9` and
`ccfdf94157d2f36afc5d9909aa85526110e5057126a1882144945df0c1cb6f96`;
the doubles have hashes
`200bb2906a127df52a560bf945d03b389611c24344a0d27f037996c57e4c88de` and
`2e47c3408a9fb11b3206d8f21272500dd54b75e8ea9b5b499c5f55fe592a111a`.
The exhaustive check is
`computations/compare_m8_class_composition.py`, with committed output
`computations/results/m8_class_composition_comparison.json`.

### 4.3 Same radial one-profile, different optimal insertion (row R)

For a signing `A`, set

```math
U_A(m)=\max_{\sum_i x_i=m}H_A(x),
\qquad
L_A(m)=\min_{\sum_i x_i=m}H_A(x).                       \tag{4.6}
```

The two explicit order-7 matrices `A_1,A_2` in
`artifacts/scale_transfer_profile_no_go.md` have the identical complete
table

```text
m                 -7  -5  -3  -1   1   3   5   7
L_A(m)              3  -9  -9  -9  -9  -9  -9   3
U_A(m)              3   7   7   7   7   7   7   3.       (4.7)
```

In particular both have `Q=9`, and the table makes both external-field
support functions in the single radial direction `h=t 1` identical for every
real `t`.  Nevertheless, if

```math
E(A)=\min_{b\in\{\pm1\}^7}\max_x
\left(|H_A(x)|+|b\mathbin\cdot x|\right),                \tag{4.8}
```

then exact enumeration gives

```math
E(A_1)=12,\qquad E(A_2)=10.                              \tag{4.9}
```

More precisely, the histograms of the row response
`R_A(b)=max_x(|H_A(x)|+|b\mathbin\cdot x|)` are

```text
          10  12  14  16
A_1        0  42  72  14
A_2        2  62  56   8.                               (4.10)
```

For `A_2`, `b=(-1,1,-1,-1,1,1,-1)` attains `10`.  The exact repair is the
full support function

```math
F_A^\sigma(h)=\max_x(\sigma H_A(x)+h\mathbin\cdot x),     \tag{4.11}
```

but it is injective because

```math
H_A(x_0)=\lim_{t\to\infty}
\bigl(F_A^+(t x_0)-tn\bigr).                             \tag{4.12}
```

So the closed exact state restores the entire labelled energy word.

### 4.4 Same Gram/second-moment states, different parent caps (row G)

Fix both order-3 children to be the positive triangle `J_3-I_3`.  The two
bridges

```math
C_{13}=\begin{pmatrix} 1&-1& 1\\ 1&-1&-1\\-1&-1&-1\end{pmatrix},
\qquad
C_{66}=\begin{pmatrix}-1& 1&-1\\-1&-1&-1\\ 1&-1&-1\end{pmatrix}             \tag{4.13}
```

have the identical full Gram pair

```math
CC^{\mathsf T}=
\begin{pmatrix}3&1&-1\\1&3&1\\-1&1&3\end{pmatrix},
\qquad
C^{\mathsf T}C=
\begin{pmatrix}3&-1&1\\-1&3&1\\1&1&3\end{pmatrix}.       \tag{4.14}
```

Yet the corresponding order-6 parent caps are `9` and `11`.  Separate
exhaustive collisions show:

* bridge codes `10,11` have the same energy-conditioned left and right
  marginal variance profiles but parent caps `11,9`;
* bridge codes `78,85` have the same conditional second moment in every
  internal-energy-shell pair but parent caps `11,9`.

These are exact finite bridge pairs.  The repair again displays the
information boundary.  If

```math
R_C(e,f)=\max_{H_A(x)=e,\ H_B(y)=f}|x^{\mathsf T}Cy|,      \tag{4.15}
```

then global sign symmetry gives the identity

```math
Q\!\left(\begin{pmatrix}A&C\\C^{\mathsf T}&B\end{pmatrix}\right)
=\max_{e,f}\bigl(|e+f|+R_C(e,f)\bigr).                    \tag{4.16}
```

But constructing `C` under all inequalities
`R_C(e,f)<=T-|e+f|` is exactly the original family of parent spin
constraints.  The exhaustive certificates are in
`computations/results/phase2j_bridge_gram_response_collision.json` and are
summarized in `artifacts/phase2j_augmented_cut_gram_response_audit.md`.

### 4.5 One selected row-field law does not certify globality (row O)

The explicit order-6 signings in
`artifacts/ar_near_minimizer_rigidity_audit.md`, equation (4.3), have the same
labelled row-sum vector

```math
(5,3,3,1,1,1).                                          \tag{4.17}
```

Hence at the all-one spin they have exactly the same selected one-profile
law of `(x_i,(Ax)_i/\sqrt6)`.  Every row sum is positive, so this spin is a
strict one-flip local maximum for both, of energy `7`.  Exact enumeration
nevertheless gives

```math
Q(A)=11,\qquad Q(B)=7.                                  \tag{4.18}
```

This is a true pair of matrices, but it is only a one-sided falsifier to a
*ground-state-law* realization theorem: the selected state is a ground state
for `B` and is not one for `A`.  It proves that the law plus one-flip
stability does not certify globality.  Adding the condition
`|H_A(y)|<=H_A(1)` for every `y` simply reinstates the full Boolean maximum.

## 5. The scalable local-profile theorem (row L)

Define `phi_6(A)` to count every switching/permutation/global-negation class
among all principal restrictions of orders 4, 5, and 6.  The two root-gauged
order-10 signings with codes

```text
5850642905,  28771662001                              (5.1)
```

have the common profile

```text
(45,165; 3,64,127,58; 0,2,11,2,11,78,39,8,55,4),       (5.2)
```

but energy ranges `[-17,19]` and `[-17,21]`.  Thus their caps are `19` and
`21`.  In fact, before quotienting global negation, their complete oriented
switching/permutation restriction-class histograms agree at every order
`r<=6`.

This is also an exact bounded-moment collision.  On 4, 5, and 6 vertices,
the class histograms are related by invertible integer character transforms
to all even-edge Eulerian signed orbit sums.  Consequently (5.2) fixes, among
other quantities,

```math
\operatorname{tr}(A^4),\quad \operatorname{tr}(A^6),
\quad \mathbb EH_A(X)^4,\quad \mathbb EH_A(X)^6,          \tag{5.3}
```

with the second moments fixed trivially.  The exact formulas and transform
determinants are in `artifacts/phase2_phi6_moment_theory_report.md`.

The finite collision becomes a theorem-level, correct-scale separation.
For `k=4^r`, let `H_k` be the symmetric Sylvester Hadamard matrix and
`D_k=diag(H_k)`.  Put

```math
S_A(k)=(A+I_{10})\otimes H_k-I_{10}\otimes D_k,           \tag{5.4}
```

and define `S_B(k)` similarly.  These are hollow signings of order
`N=10k`.  Every oriented restriction profile through six vertices remains
equal: a small restriction is determined by its macro support, its selected
micro-coordinates/occupancies, and the oriented base class, and the lift is a
fixed transformation on that data.

For the first base, the rational certificate

```math
y={1\over200}(445,490,661,668,436,645,405,427,485,513)
```

satisfies

```math
\operatorname{Diag}(y)\pm{A+I\over2}\succ0,
\qquad \sum_i y_i={207\over8}.                           \tag{5.5}
```

Exact leading-principal-minor or rational `LDL^T` checks prove

```math
Q(S_A(k))\le {207\over8}k^{3/2}.                         \tag{5.6}
```

For the second base, a Boolean spin of positive energy `21` has
`s^T(B+I)s=52`; tensoring it with the Boolean `+\sqrt k` eigenvector of
`H_k` gives

```math
Q(S_B(k))\ge26k^{3/2}.                                  \tag{5.7}
```

Therefore

```math
Q(S_B(k))-Q(S_A(k))
\ge {1\over8}k^{3/2}
={1\over8\,10^{3/2}}N^{3/2}.                            \tag{5.8}
```

Both sides have `O(N^(3/2))` cap: for the second family this follows directly
from
`||S_B(k)||_op <= ||B+I||_op sqrt(k)+1=O(sqrt(k))` and the spectral cap
bound.  This is a true scalable pair at the
project scale, although its constants are not near the unknown optimum.  It
refutes universal control by any fixed local profile through six vertices,
and hence by the corresponding finite set of moments, to `o(N^(3/2))`.
The base equality and certificates were independently rechecked by
`computations/phase2b_verify_phi6_collision.py`; the lift proof is in
`artifacts/phase2b_phi6_collision_report.md`.

## 6. One-profile topology without uniform integrability (row A)

Let `s=floor(n^(3/4))`, `t=n-s`.  Choose a competitive `t`-vertex signing
`D`, a competitive `s`-vertex signing `E`, and an `s` by `t` sign bridge `C`
with

```math
\|C\|_{\infty\to1}=O(n^{11/8}).                          \tag{6.1}
```

Such a bridge exists by Hoeffding plus a union bound.  Define

```math
A_n=\begin{pmatrix}E&C\\C^{\mathsf T}&D\end{pmatrix},
\qquad
B_n=\begin{pmatrix}J_s-I_s&C\\C^{\mathsf T}&D\end{pmatrix}.              \tag{6.2}
```

The triangle inequality gives

```math
Q(A_n)\le(1/2+o(1))n^{3/2}.                             \tag{6.3}
```

Let `P(D)` be the positive one-sided extremum and use the universal lower
bound `P(D)>=c_*t^(3/2)` after a global sign choice.  On an all-one clique
spin and a positive endpoint spin for `D`, flipping the entire clique block
makes the cross term nonnegative while preserving both internal energies.
Hence

```math
Q(B_n)\ge(1/2+c_*-o(1))n^{3/2}.                         \tag{6.4}
```

For `T_A=A/\sqrt n` on the uniform vertex probability space, the two
operators differ only on the `s` clique coordinates.  Coupling the same
bounded inputs shows, for every fixed action profile level `k`,

```math
d_H(\mathcal S_k(T_{A_n}),\mathcal S_k(T_{B_n}))
\le {s\over n}\longrightarrow0,                         \tag{6.5}
```

and hence `d_M(T_A,T_B)->0`, even though (6.3)--(6.4) have a fixed gap.
This is a true pair and is stronger than equality of one selected profile.

The caveat is exact and essential: the clique creates output of size
`s/\sqrt n=n^(1/4)` on a set of mass `s/n=n^(-1/4)`.  Weak profile distance
forgets the set, while its `L^1` energy contribution is order one.  Thus the
pair has no common normalized `2 -> 2` bound and fails uniform
integrability.  Under a common `2 -> 2` bound the repository's action
continuity inequality does control the objective, so row A must not be used
against that stronger topology.  Source:
`artifacts/action_convergence_boolean_spikes.md`.

## 7. Bounded-operator traffic, low replicas, and bounded-beta pressure (rows T and P)

### 7.1 The true Walsh/planted pair

Let `n=2^d`, index coordinates by `F_2^d`, and let

```math
W_{uv}=(-1)^{u\mathbin\cdot v}.
```

Delete the diagonal to obtain the hollow signing `A_n`.  Since `W^2=nI`,

```math
\|A_n\|_{op}\le\sqrt n+1,
\qquad
{Q(A_n)\over n^{3/2}}\le {1\over2}+o(1).                \tag{7.1}
```

The negative-edge graph of `A_n` has `n(n-2)/4` edges.  Retain each such
edge with probability `4/(3\sqrt n)` and choose a deterministic realization
`F_n` supplied by concentration and matrix Bernstein, so that

```math
|F_n|=(1/3+o(1))n^{3/2},
\qquad
\|\operatorname{Adj}(F_n)\|_{op}=O(\sqrt n).             \tag{7.2}
```

Reverse those negative edges to obtain `B_n`.  Then

```math
\left\|{A_n\over\sqrt n}\right\|_{op}
+\left\|{B_n\over\sqrt n}\right\|_{op}=O(1),            \tag{7.3}
```

whereas the all-one spin gives

```math
\liminf_n {Q(B_n)\over n^{3/2}}\ge {2\over3}.            \tag{7.4}
```

At the same time,

```math
\left\|{A_n-B_n\over\sqrt n}\right\|_{2,\tau}^2
={8|F_n|\over n^2}=O(n^{-1/2})\longrightarrow0.          \tag{7.5}
```

Telescoping under (7.3) makes every fixed normalized trace polynomial agree
asymptotically.  A fixed signed test-graph density changes by at most
`O_F(|F_n|/n^2)=O_F(n^(-1/2))`.  Finally, for a uniform spin `X`,

```math
\mathbb E\left|
{H_{B_n}(X)-H_{A_n}(X)\over n^{3/2}}
\right|^2
={4|F_n|\over n^3}=O(n^{-3/2}).                          \tag{7.6}
```

Coupling any fixed number of uniform replicas with the same spins therefore
makes their joint normalized energy laws converge in Wasserstein distance.
If the summary also records the finite overlap array
`((x^a\mathbin\cdot x^b)/n)^2`, that array agrees *exactly* under this coupling.
Thus every fixed sampled overlap-energy law agrees as well.

This last assertion is about sampled laws, not support sets.  The planted
all-one direction is an exponentially sparse atom.  Retaining the support of
the one-replica law would already retain the maximum and is not contradicted
by this pair.  The theorem and construction are Proposition 4.1 of
`artifacts/traffic_laplace_principle.md`.

### 7.2 A bounded-temperature pressure corollary

Define

```math
\overline Z_n(X)=2^{-n}\sum_x
\cosh\left({1\over2}x^{\mathsf T}Xx\right).
```

The proved strict-high-temperature stability theorem says that for every
`kappa<1/2`, if `||X||_op,||Y||_op<=kappa`, then

```math
|\log\overline Z_n(X)-\log\overline Z_n(Y)|
\le {K_\kappa\over2}\|X-Y\|_*.                          \tag{7.7}
```

Let `C` bound both normalized operator norms in (7.3), and choose
`beta_0>0` and `kappa<1/2` with `beta_0 C<kappa`.  Apply (7.7) to

```math
X={\beta A_n\over\sqrt n},
\qquad
Y={\beta B_n\over\sqrt n}.
```

Nuclear/Frobenius comparison and (7.5) give, uniformly for
`|beta|<=beta_0`,

```math
\begin{aligned}
{1\over n}|\log\overline Z_n(X)-\log\overline Z_n(Y)|
&\le {K_\kappa\over2}\|X-Y\|_{2,\tau}\\
&\le {K_\kappa\beta_0\over2}
\left\|{A_n-B_n\over\sqrt n}\right\|_{2,\tau}
\longrightarrow0.                                      \tag{7.8}
\end{aligned}
```

Equations (7.1), (7.4), and (7.8) give an actual signing separation:
normalized cosh pressure agrees on a nontrivial compact beta interval, while
the normalized extrema differ by at least `1/6-o(1)`.  This corollary is a
new combination of two proved repository results, not an additional
computation.  The pressure input is
`artifacts/high_temperature_frobenius_pressure_stability.md`.

The interval cannot be enlarged from these inputs alone.  The covariance
bound behind (7.7) requires the strict `1/2` operator-temperature margin, and
at sufficiently low temperature a planted state may become thermodynamically
visible.  Row P therefore says “bounded high-temperature interval,” not
“all fixed beta.”

## 8. One-sided blindness and the schematic pressure boundary

### 8.1 Vanishing fourth defect can hide a complete Boolean spike (row D)

Split one vertex from a symmetric conference signing of order `n+1`:

```math
C=\begin{pmatrix}0&b^{\mathsf T}\\b&A\end{pmatrix},
\qquad C^2=nI.
```

Exact multiplication gives

```math
Ab=0,
\qquad
A^2-(n-1)I=I-bb^{\mathsf T}.                            \tag{8.1}
```

Therefore

```math
\Delta_4(A):={\|A^2-(n-1)I\|_F^2\over n^3}
={n(n-1)\over n^3}\longrightarrow0,                     \tag{8.2}
```

while the Boolean coloring `b` obeys

```math
\left\|{A^2-(n-1)I\over n}b\right\|_{L^2([n])}
={n-1\over n}\longrightarrow1.                          \tag{8.3}
```

This proves that empirical fourth-moment/two-walk flatness is not a universal
Boolean-profile statement.  It is deliberately classified one-sided:
(8.1)--(8.3) do not exhibit another signing with the same `Delta_4` and a
separated cap.  Source: `artifacts/ar_near_minimizer_rigidity_audit.md`,
Section 2.2.

The same source's universal-positive-vertex extension is likewise
one-sided.  It preserves normalized near-optimality and an `O(\sqrt n)`
operator bound while putting a field of size `n` on one vanishing-mass
coordinate.  It falsifies square-field uniform-integrability forcing, but it
is not a same-summary cap-separated pair.

### 8.2 The abstract entropy-spike landscape (row X)

Let one landscape have energy zero on all `2^n` spins.  Let a second,
spin-flip-symmetric landscape have energy `delta n^(3/2)` on one antipodal
pair and zero elsewhere.  Its normalized maximum is `delta`, but its pressure
is

```math
{1\over n}\log\left[
1-2^{1-n}+2^{1-n}\cosh(\beta\delta n)
\right].                                                \tag{8.4}
```

For every finite `B` and every `delta<log(2)/B`, (8.4) tends to zero
uniformly on `0<=beta<=B`, exactly as for the zero landscape.  This is the
sharp schematic explanation of the entropy window (2.6).  It is not asserted
that either arbitrary table is the energy word of a complete signing.  The
actual row-P pair is the quadratic-sign realization of the same blindness on
some nontrivial high-temperature interval.

The separate abstract scalar-pressure theorem in
`artifacts/finite_temperature_scalar_no_go.md` has a different role.  It
constructs analytic convex pressures satisfying centered subadditivity,
monotonicity, the correct second derivative, entropy squeezes, and uniform
beta-Lipschitz bounds while their diagonal values and zero-temperature slopes
oscillate.  It is a proved theorem about abstract functions but remains
**schematic** as evidence about signings.

## 9. Coding/coset reinterpretation

Let `E=binom(n,2)` and identify a signing with a bit word `a` on the edges.
For the augmented cut code

```math
\mathcal C_n^+
=\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\},
```

the distance enumerator of the coset `a+\mathcal C_n^+` is

```math
W_a(z)
=\sum_{\sigma=\pm1}\sum_{[x]}
 z^{(E-\sigma H_A(x))/2},                               \tag{9.1}
```

where `[x]` is a spin modulo global sign.  In particular,

```math
d(a,\mathcal C_n^+)={E-Q(A)\over2}.                      \tag{9.2}
```

Thus `Q(A)` is the deficit associated with the coset-leader weight, and the
complete absolute energy histogram is precisely the one-coset weight data.

Rows L and H become two distinct coding obstructions.

1. **Bounded puncturing decks do not determine global coset depth.**  For a
   vertex set `U`, restricting `A` to `U` punctures the edge word to the
   clique-coordinate set `E(K_U)`.  Quotienting by switching, global
   negation, and vertex permutation records the local augmented-cut coset
   type.  The order-10 pair (5.1) has the same multiset of these punctured
   types for every `|U|<=6`, but (9.2) gives leader weights

   ```math
   d(a,\mathcal C_{10}^+)=13,
   \qquad
   d(b,\mathcal C_{10}^+)=12.                            \tag{9.3}
   ```

   More strongly, the Hadamard lifts have the same bounded puncturing deck
   and

   ```math
   d(S_A(k),\mathcal C_{10k}^+)
   -d(S_B(k),\mathcal C_{10k}^+)
   ={Q(S_B(k))-Q(S_A(k))\over2}
   \ge {1\over16}k^{3/2}.                                \tag{9.4}
   ```

   This is a scalable separation between local coset data and global
   coset-leader depth.

2. **A complete source-coset enumerator is not functorial under a fixed
   lift.**  The two order-8 minimizer cosets have the same enumerator (9.1)
   and common leader weight `(28-10)/2=9`.  Their universal doubles live at
   length `binom(16,2)=120`, but (4.5) gives leader weights

   ```math
   d(\mathcal D(S_0),\mathcal C_{16}^+)=40,
   \qquad
   d(\mathcal D(S_1),\mathcal C_{16}^+)=44.              \tag{9.5}
   ```

   Hence no map from the source coset weight enumerator alone can predict
   even this fixed lifted leader weight.  The missing invariant is the
   alignment of source energy layers with the lift's new coordinates, not
   another scalar moment of the same enumerator.

This coding view also explains why a growing Krawtchouk or dual-degree repair
eventually ceases to be a compression: the full signed Eulerian/Krawtchouk
transform is invertibly equivalent to (9.1).  Bounded dual degree sees only
bounded moments; degree high enough to isolate the extreme tail reconstructs
the coset histogram/support that contains the original optimization.

## 10. Conclusions for summary design

The examples support four precise, limited conclusions.

1. **Typical laws are not support sets.**  Fixed moments, fixed replicas,
   graph/traffic limits, and bounded-beta pressure average over states.  A
   planted antipodal pair can carry a leading maximum and zero entropy.  Row T
   proves this even under a common normalized operator bound.
2. **Local decks are not global coset depth.**  Row L is the strongest
   theorem-level finite-profile obstruction: exact equality of every local
   oriented type through six vertices survives a low-scale lift while the
   caps remain separated at order `N^(3/2)`.
3. **Composition needs alignment, not only marginals.**  Rows H, R, and G
   successively defeat the complete child energy histogram, a full radial
   extremal profile, bridge Gram data, and shell-conditioned second moments.
   Their exact repairs, (4.11) or (4.15), restore the full labelled response.
4. **Uniform integrability is a real boundary.**  Row A defeats weak action
   profiles only by a vanishing-mass/unbounded-output spike.  It does not
   contradict quantitative action continuity under a common `2 -> 2` bound.
   Conversely, row T shows that bounded operator norm alone does not make
   trace, replica, or high-temperature-law summaries extremum-continuous.

Accordingly, none of spectrum, bounded moments, bounded local profiles, an
unlabeled energy histogram, one selected profile, finitely sampled overlap
laws, or bounded high-temperature pressure is a lossless summary of the
Boolean extremum or its composition response.  The proved repairs all retain
one of two expensive objects: the extremal support over Boolean directions,
or the complete response envelope against future fields/bridges.
