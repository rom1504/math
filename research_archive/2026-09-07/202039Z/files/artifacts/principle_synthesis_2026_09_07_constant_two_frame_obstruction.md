# Constant two-spin weave: an actual independent-frame obstruction

2026-09-07. **Proved, with independent reconstruction of the conditioning
and repair martingale by the construction reviewer. Full director audit
requested.** This is a new actual-parent cap statement at the constant
multiplier `N=2n`. It is not a fixed-temperature pressure obstruction.

## 1. The exact seed-preserving operator and the theorem

Let `S` be any prescribed hollow symmetric signing of order `n`. For
directed signs `eta_ij`, define hollow matrices

```
C_ij=S_ij eta_ji,       D_ij=S_ij eta_ij eta_ji,
W=[[S,C],[C^T,D]].                                      (1)
```

Fill the `n` missing matching edges, the diagonal entries of `C`, with
arbitrary fixed signs. This makes `W` a full hollow signing of order `2n`
and changes its cap by at most `n`. Each off-diagonal macro tile is the
rank-one matrix

```
S_ij [1;eta_ij] [1,eta_ji].
```

In particular the first-half principal child is EXACTLY `S`, without
old-edge edits, random final column gauges, or a relabeling of seed data.

Put

```
c_two=(2/pi+sqrt(2/pi))/(2sqrt(2))
     =0.5071738708131547... .                            (2)
```

There are two statements.

* If all directed off-diagonal `eta_ij` are independent fair signs, then
  for EVERY fixed seed `S`, uniformly in that seed,
  `Q(W)/(2n)^(3/2)>=c_two-o_P(1)`.
* If `n` is even, independently choose each full length-`n` row `eta_i`
  uniformly among the exactly balanced sign rows, and use its
  off-diagonal coordinates in (1), then the same lower bound holds
  uniformly for every fixed `S` with `Q(S)<=K n^(3/2)`, for fixed `K`.

The second law consists of genuinely orthogonal TWO-ROW PHYSICAL FRAMES
`H_i=[1;eta_i]`, since `H_i H_i^T=n I_2`. Orthogonality is within each
frame, not between the incidence rows belonging to different fibres.
At Hadamard orders these frames are exactly the first two rows of a
dephased Hadamard after independent ordinary column permutations.

Thus the independent-frame constant-multiplier revival fails even when
one first selects any desired exact child minimizer. The quantifier is
EVERY fixed seed, followed by random frames, not a single event holding
simultaneously over adaptively chosen seeds. Correlated frames, favorable
rare frame choices, and global edge rewrites are not excluded.

## 2. Exact conditional bridge response in the iid law

The matrix `D` in (1) is an iid hollow symmetric signing. Conditional on
`D`, each unordered edge has one independent fair orientation variable
`C_ij`, and

```
C_ji=S_ij D_ij C_ij.                                    (3)
```

Choose `y` maximizing the POSITIVE energy `H_D(y)`, and set

```
g_i=sum_(j!=i) C_ij y_j,      x_i=sign(g_i),
```

with independent unbiased tie coins. This definition uses hollow `C`;
the already filled matching edges are paid separately by at most `n`.
Conditional on `D,y`, each `g_i` has the distribution of a sum of `n-1`
independent fair signs. Therefore

```
E[x^T C y |D,y]=n E|epsilon_1+...+epsilon_(n-1)|.         (4)
```

The two fields `g_i,g_j` share only their edge variable. Their other
summands are independent. If `R` is a sum of `n-2` fair signs, put

```
lambda_n=P(R=0)  if n is even,
lambda_n=P(R=1)  if n is odd,
kappa_n=lambda_n^2 ~ 2/(pi n).
```

Equivalently `lambda_n=2^(-(n-2))*binom(n-2,floor((n-2)/2))`.
The unbiased-tie convention gives
`E sign(R+a)=a lambda_n` for `a=+-1`. Using (3) on the shared edge gives
the exact identity

```
E[x_i x_j|D,y]=kappa_n S_ij D_ij y_i y_j,
E[H_S(x)|D,y]=kappa_n H_D(y).                            (5)
```

The seed contribution therefore has a FAVORABLE expectation; it is not
being bounded below by `-Q(S)`. Combining (4),(5), and the matching cost,

```
E Q(W) >= (1+kappa_n) E P(D)
             +n E|sum_(1)^(n-1) epsilon_j|-n,            (6)
P(D)=max_y H_D(y).
```

No assumption about the child cap, polarity balance, or operator norm is
needed in this iid step.

## 3. An elementary `2/pi` random-child witness

For completeness, the random-child estimate needed in (6) is

```
liminf E P(D)/n^(3/2) >= 2/pi.                           (7)
```

It follows from fixed-degree Gaussian spectral rounding and the elementary
closed-walk semicircle moments, without the Parisi formula or an SK
ground-state limit. Here are the uniform-integrability details.

Write `B=D/sqrt(n)`. Fix a real polynomial `p` and `epsilon>0`, put
`R=p(B)^2+epsilon I`, and Gaussian-round a centered vector with covariance
`R`. Let `q=p^2`, `v=int q dmu_sc+epsilon`, and
`Lambda=diag(R_ii^(-1/2))`. For this FIXED polynomial the usual expansion
into finitely many walks gives

```
E sum_i (R_ii-v)^2=O(1),
E tr(R^2)=O(n),
E sum_(i!=j) R_ij^4=O(1),
E tr(B R)=n int x q(x) dmu_sc(x)+O(1).                  (8)
```

To check the only less familiar estimate in (8), expand four walks from
`i` to `j`, with `i!=j`. A nonzero expectation requires every edge to
occur an even number of times. If their union contains a cycle, its
vertex count is at most its edge count, at most half the total length.
If it is a tree, every edge on the path from `i` to `j` occurs at least
four times across the four walks. The same vertex-count bound follows.
After summing labels, every contributing shape is consequently `O(1)`.
There are finitely many shapes. The diagonal variance follows by the
same connected-pair walk count after subtracting the disconnected
expectations; the remaining Catalan tree terms give the other two limits.

Since every `R_ii>=epsilon`, the first estimate implies
`E||Lambda-v^(-1/2)I||_F^2=O(1)`. The deterministic identity
`sum_i B_ij^2=(n-1)/n` then gives, by Frobenius Cauchy--Schwarz,

```
E |tr[B Lambda R Lambda]-v^(-1)tr(BR)|=O(sqrt(n)).       (9)
```

For example, expand the difference as
`Delta R Lambda+v^(-1/2)R Delta`, and pair `B Delta` with `R Lambda`.
This avoids any implicit operator-norm uniform integrability assumption.

The Gaussian sign identity is
`E[y_i y_j|D]=(2/pi)arcsin((Lambda R Lambda)_ij)`.
The universal inequality `|arcsin z-z|<=C|z|^3`, together with (8),
implies that the sum of the off-diagonal cubic errors has expectation
`O(sqrt(n))`. Equations (8),(9) therefore give

```
liminf E P(D)/n^(3/2)
 >= int x p(x)^2 dmu_sc(x) /
       [pi (int p(x)^2 dmu_sc(x)+epsilon)].              (10)
```

Take `p(x)=(x+2)^r`, first send `epsilon` to zero, and then let the fixed
integer `r` increase. A beta integral gives the exact quotient

```
int x(x+2)^(2r) dmu_sc / int (x+2)^(2r) dmu_sc
 =2-6/(2r+3) -> 2.
```

This proves (7). All order limits precede increasing the polynomial
degree. The elementary fair-sign absolute-sum asymptotic in (6) is
`sqrt(2(n-1)/pi)+o(sqrt(n))`, so (6) yields the constant (2).

Finally changing one directed `eta_ij` changes only two physical parent
edges, hence changes `Q(W)` by at most `4`. Bounded differences over
`n(n-1)` independent signs turns the expectation bound into a
high-probability lower bound, with failure at most `exp(-c_delta n)`
for any fixed gap `delta` below (2). This proves the first theorem.

## 4. Coupling iid rows to exactly balanced physical frames

Let `n` be even and extend each original iid incidence row to length `n`
by an unused fair diagonal port. In row `i` write its signed imbalance
as `b_i`, let `r_i=|b_i|/2`, and let `M_i=(n+|b_i|)/2` be its majority
count. Uniformly select `r_i` majority coordinates to flip. Do this
independently at every row, obtaining `eta'_i`.

The output is exactly balanced. Coordinate-permutation invariance of
the coupling makes its marginal uniform on all balanced rows, and the
rows remain independent. Thus this constructs the second law exactly,
not through conditioning an exponentially rare event.

Conditional on the original rows put

```
p_i=r_i/M_i,       a_i=sign(b_i)p_i,
mu_ij=E[eta'_ij|eta]=(1-p_i)eta_ij-a_i.                  (11)
```

If `r_i=0`, take `p_i=a_i=0`. Write `X_i=eta'_i-mu_i`. Sampling without
replacement gives the deterministic norm and covariance estimates

```
||X_i||_2^2=4r_i(1-p_i)<=4r_i,
Cov(X_i|eta)
 =4p_i(1-p_i)[M_i/(M_i-1) P_i
                    -(1/(M_i-1))1_i 1_i^T]
 <=8p_i I,                                              (12)
```

where `P_i` projects onto majority coordinates and `1_i` is their
indicator. The zero case is interpreted separately. Restricting to
off-diagonal coordinates preserves these bounds.

With probability tending to one, uniformly over rows,

```
r_max=O(sqrt(n log n)),
p_max, |a|_max=O(sqrt(log n/n)),
||C||op+||D||op=O(sqrt(n)).                              (13)
```

The first two follow from scalar Hoeffding and a union bound. The last
follows from the elementary sphere-net subgaussian estimate: `C` is iid
off-diagonal and `D` is symmetric iid off-diagonal. Their dependence on
each other is immaterial for taking a union bound.

## 5. The conditional means cost only a subleading operator perturbation

Let `C',D'` be formed from `eta'`, and use diagonal matrices `P=diag(p_i)`,
`A_0=diag(a_i)`, `F=I-P`. Equations (11) and row independence give EXACTLY

```
E[C'|eta]=C F-S A_0,
E[D'|eta]=F D F-A_0 C F-F C^T A_0+A_0 S A_0.             (14)
```

All displayed matrices are hollow, so no diagonal correction is hidden
in (14). For the seed class `Q(S)<=K n^(3/2)`, the already established
bounded-entry fourth-moment inequality gives `||S||op=O_K(n^(3/4))`.
It is preserved explicitly in
`principle_synthesis_2026_09_07_logarithmic_local_law.md`, equation (3).
Combining this with (13),(14) shows

```
||E[C'|eta]-C||op=O_K(n^(1/4)sqrt(log n)),
||E[D'|eta]-D||op=O_K(sqrt(log n)+n^(-1/4)log n).         (15)
```

The term `S A_0` is why a seed condition is needed for this coupling.
No bounded operator-norm assumption stronger than the actual cap bound
has been imposed.

## 6. Row-exposure fluctuations are also subleading

Conditional on the original rows, the repaired rows are independent.
For `C'-E[C'|eta]`, expose one repaired row at a time. The corresponding
column increment has norm at most `2sqrt(r_i)`. The two predictable
quadratic variations are bounded by
`O(r_max+sum_i p_i)I`, by (12).

For `D'`, its Doob martingale has a particularly useful exact form.
When row `i` is exposed, the increment is

```
Y_i=e_i z_i^T+z_i e_i^T,
z_i(j)=S_ij X_i(j) * eta'_ji     if j<i,
z_i(j)=S_ij X_i(j) * mu_ji       if j>i,
z_i(i)=0.                                               (16)
```

The multipliers in (16) are predictable and have magnitude at most one.
Thus `||Y_i||=||z_i||<=2sqrt(r_i)` and
`E[z_i z_i^T|past]<=8p_i I`. Since

```
Y_i^2=||z_i||^2 e_i e_i^T+z_i z_i^T,
```

the total predictable quadratic variation is at most
`(16r_max+8sum_i p_i)I`. In particular it is
`O(sqrt(n log n))I` on (13). There is no omitted quadratic interaction
between two newly exposed rows: the Doob formula already includes the
earlier row as an actual coefficient and the later row as its mean.

Apply the self-adjoint and rectangular matrix Freedman inequalities,
with increment bound `O((n log n)^(1/4))` and the preceding variance
bound. The exact primary source is [Tropp, Theorem 1.2 and Corollary 1.3](https://tropp.caltech.edu/papers/Tro11-Freedmans-Inequality.pdf).
They give, with conditional probability tending to one,

```
||C'-E[C'|eta]||op+||D'-E[D'|eta]||op
 =O(n^(1/4)(log n)^(5/4))=o(sqrt(n)).                    (17)
```

For example a sufficiently large constant times the displayed threshold
has a polynomially small tail after the matrix dimension factor; all
parameters in this application are deterministic conditional on the
original rows satisfying (13).

Combining (15),(17), and leaving the matching edges unchanged, gives an
actual coupling of the two full parents with

```
||W'-W||op=O_K(n^(1/4)(log n)^(5/4)),
|Q(W')-Q(W)|<=n ||W'-W||op=o(n^(3/2)).                  (18)
```

Therefore the iid-law lower bound transfers to independent exactly
orthogonal two-row physical frames. This proves the second theorem.

## 7. Exact original-problem consequence and remaining escape

For every fixed bounded-cap seed, including any exact child optimizer
selected before the frame randomness, this concrete `N=2n` construction
has normalized cap above `0.50717-o(1)` with high probability. It cannot
land the current original upper bound below `0.493608094`, let alone
give a lossless recurrence from a hypothetical smaller limiting value.

Unlike the old `k~n` compiler, it genuinely retains the seed in an actual
principal child and does not suppress its contribution by a parent-scale
factor. Its failure is instead caused by a random opposite child plus
a bridge field that can be exploited jointly, even after exact local
orthogonality is restored. The proof does not forbid correlated frames
or existence of specially selected completions. In particular it is
consistent with the separately proved correlated anti-invariant
subhalf family and does not establish an asymptotic obstruction for
all selectable exact-child lifts.

## 8. Bounded exact checks

The fraction-arithmetic script
`computations/principle_synthesis_2026_09_07_constant_two_frame_check.py`
checks the shared-edge identity at orders `3,...,9`, and exhausts every
input row and every permitted majority repair at lengths `4,6,8`.
It verifies the conditional means, every covariance entry in (12), the
deterministic centered norm, and the exact uniform balanced-row
pushforward. All tests PASS in about one second. The output is preserved
in `computations/results/principle_synthesis_2026_09_07_constant_two_frame_check.json`.
These checks audit finite algebra only; the asymptotic theorem follows
from the proof above, not from extrapolating these examples.

## 9. A positive repair theorem for correlated inputs

The repair itself is not confined to the unfavorable random-frame law.
Suppose a DETERMINISTIC full parent of the rank-one form (1) has
`Q(W)<=K n^(3/2)`, and each length-`n` incidence row has imbalance at most
`2r_*`, with `n` even. An unused diagonal port may be appended first,
changing this bound by at most one. Apply the same independent majority
repairs conditional on these given rows. Then, with high probability,

```
|Q(W')-Q(W)|
 <= O_K(r_* n^(3/4)+n sqrt(r_*) log n),                 (19)
```

while the first child remains EXACTLY `S` and every two-row physical
frame is exactly orthogonal. The input incidence rows may have arbitrary
cross-row dependence or have been selected by a global optimization.

Indeed the cap bound and the fourth-moment theorem give
`||S||op+||C||op+||D||op=O_K(n^(3/4))`, by compression of the full parent.
Now `p_max<=2r_*/n`, so the deterministic conditional-mean change in (14)
is `O_K(r_* n^(-1/4))` in operator norm. Equations (12),(16) give
predictable variance `O(r_*)` and increment norm `O(sqrt(r_*))` for the
repair noise; Freedman bounds it by `O(sqrt(r_*)log n)`. Multiplication
by `n` converts operator norm to the cap difference, proving (19).

In particular if `r_*=o(n^(3/4))`, this is an actual
`o(n^(3/2))`-cost orthogonalization theorem for the chosen parent, even
for a correlated or globally selected input. It removes a local-frame
constraint in that quantified near-balanced class; it does not construct
the required low-cap parent from its child. No iid-law obstruction is
being extrapolated to these correlated inputs.
