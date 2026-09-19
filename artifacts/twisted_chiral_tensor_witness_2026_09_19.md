# Explicit tensor witnesses and optimal diagonal Gram defect

Status: rigorous director derivation, independently reconstructed. The even
self-tensor classification is ALREADY in the archive; this gives a much
shorter explicit witness, removes its Gaussian-rounding dependency, and
quantifies every diagonal completion of the new chiral witnesses.
It is not an obstruction to arbitrary twists, restrictions, or convergence.

## 1. A direct Boolean witness, with no rounding or polarization

Let H be a symmetric full sign matrix of order q>=2, including its diagonal.
Put

```math
\rho=\frac{\operatorname{tr}H^4}{q^3},\qquad
\tau=\frac{(\operatorname{tr}H)^2}{q^3}.
```

The vector z=vec(H) is an actual sign vector and

```math
(H\otimes H)\operatorname{vec}(H)=\operatorname{vec}(H^3),
\qquad z^T(H\otimes H)z=\operatorname{tr}H^4.       \tag{1}
```

For each integer k>=1, the Boolean vector z^tensor k therefore proves

```math
\frac{Q(\operatorname{hollow}(H^{\otimes2k}))}{q^{3k}}
\ge\tfrac12(\rho^k-\tau^k).                       \tag{2}
```

Hollowing subtracts the spin-independent trace (tr H)^(2k). There is no
bilinear-to-quadratic loss: (1) is already a same-spin witness.
Since tr H²=q², Cauchy--Schwarz on the squared eigenvalues gives
tr H⁴>=q³, with equality exactly when H²=qI. Also tau<=1/q.
Thus non-Hadamard fixed seeds give divergence in(2); symmetric Hadamard
seeds have normalized even-power cap tending to1/2, by the matching
spectral upper bound. This is the classification previously proved in
`transfer_seed_tensor_power_instability_2026_09_06.md`, now without its
vector relaxation, Gaussian sign identity, or scalar polarization.

For a non-Hadamard seed, odd powers diverge too. Choose a fixed sign v with
w=v^T H v nonzero. The all-positive v works unless its value is zero; in
that case q is even, and one coordinate flip changes the value by -4 times
an odd off-diagonal row sum, so supplies such v. Tensoring v with z^tensor k
gives

```math
\frac{Q(\operatorname{hollow}(H^{\otimes(2k+1)}))}{q^{3k+3/2}}
\ge\frac{|w|}{2q^{3/2}}\rho^k
 -\frac{|\operatorname{tr}H|}{2q^{3/2}}\tau^k.       \tag{3}
```

This does NOT identify the limiting odd-power cap for a general Hadamard
seed. The non-Hadamard divergence extension is immediate from the displayed
witness, not a classification of arbitrary varying-seed constructions.

The witness also applies to heterogeneous P=tensor_j(H_j tensor H_j),
of total order N=product_j q_j²:

```math
Q(\operatorname{hollow}P)/N^{3/2}
\ge\tfrac12\prod_j\frac{\operatorname{tr}H_j^4}{q_j^3}
       -\frac1{2\sqrt N}.                          \tag{4}
```

Thus literal paired self-tensors never retain a fixed asymptotic discount
below1/2. Every off-diagonal Gram defect contributes multiplicatively.
The statement allows growing seed orders, but requires the exact paired
tensor structure. Selecting a principal restriction can evade that structure,
as the project's strict upper constructions do.

## 2. Exact optimization over all sign diagonals

Let D be any symmetric hollow full signing of order q. Define

```math
c_i=(D^3)_{ii},\qquad
\Omega_0=\sum_{i\ne j}(D^2)_{ij}^2,
\qquad H=D+\operatorname{diag}h,\quad h_i\in\{\pm1\}.
```

Its full Gram defect has the exact expression

```math
\Omega(H)=\|H^2-qI\|_F^2
=\Omega_0+2q(q-2)+2\left(\sum_i h_i\right)^2+4\sum_i c_i h_i.
                                                               \tag{5}
```

Indeed the diagonal of H² is q and its off-diagonal entries are
(D²)_ij+D_ij(h_i+h_j). Expanding their squared sum gives(5), using
sum_j D_ij(D²)_ij=(D³)_ii. In particular tr H⁴=q³+Omega(H).

For exactly k positive h_i, minimize the linear term by choosing the k
smallest c_i. If c_(1)<=...<=c_(q), the minimum over all sign diagonals is

```math
\min_{0\le k\le q}\left[
\Omega_0+2q(q-2)+2(2k-q)^2-4\sum_i c_i+8\sum_{i=1}^k c_{(i)}
\right].                                                       \tag{6}
```

After computing D² and the triangle row sums c, this needs only sorting and
a prefix sweep, O(q log q) operations; that is NOT the total preprocessing
cost. It covers all2^q diagonal completions without searching them.

## 3. Stronger exact consequence of chirality

For D=[[A,C],[C,-A]], the signed quarter-turn J satisfies JDJ^T=-D.
Hence c_i and c_(i+q/2) are opposites. The linear term in(5) is minimized
by h_i=-sign(c_i), with opposite h values in each zero pair. This choice
also has sum h=0, so minimizes the quadratic term simultaneously. Therefore

```math
\Omega_*:=\min_h\Omega(D+\operatorname{diag}h)
=\Omega_0+2q(q-2)-4\sum_i|c_i|.                    \tag{7}
```

For every other diagonal, writing s=sum h, (5) implies
Omega(H)>=Omega_*+2s². Thus, for every k>=1,

```math
(\operatorname{tr}H^4)^k-s^{2k}
\ge(q^3+\Omega_*+2s^2)^k-s^{2k}
\ge(q^3+\Omega_*)^k.
```

The same-spin tensor witness gives the uniform, all-diagonal bound

```math
\boxed{\quad
Q(\operatorname{hollow}(H^{\otimes2k}))/q^{3k}
\ge\tfrac12(1+\Omega_*/q^3)^k.
\quad}                                                       \tag{8}
```

There is no diagonal correction in(8). If Omega_*=0, every completion's
even self-tensors have normalized cap at least1/2, and a Hadamard completion
attains equality. If Omega_*>0, every completion diverges uniformly for this
fixed D. This is a statement about an explicitly declared tensor family,
not about the original minimum over all signings.

## 4. The two cap40 witnesses give different tensor failures

The exact verifier
`computations/twisted_chiral_tensor_gram_2026_09_19.py` checks6,096 sign
diagonals of84 random small full signings, checks the vectorization identity
directly, and computes(7) for both stored order20 parents.

The Hadamard-completable witness has Omega_*=0. Its best completion's
self-tensor square has cap exactly4,000 at order400, normalized1/2.
The non-completable witness has Omega_*=960. For EVERY sign diagonal its
self-tensor square has cap at least4,480, normalized0.56, and(8) gives
the growth lower bound(1/2)(28/25)^k for even powers2k.

These are exact integer certificates, not spectral numerics or heuristic
tensor optimization. They explain why neither finite40 witness by itself
supplies an asymptotic improvement via literal self-tensoring. They do not
exclude a genuinely different recursive twist, a restriction, or a global
replacement operation with controlled cap loss.
