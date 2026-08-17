# Independent audit: Walsh-prefix nonconvergence

Audited files:

- [`walsh_prefix_nonconvergence.md`](walsh_prefix_nonconvergence.md)
- [`verify_walsh_prefix_nonconvergence.py`](../experiments/verify_walsh_prefix_nonconvergence.py)

**Verdict: PASS.**  The draft defines one coherent infinite sign matrix in
the standard Kronecker order, every `C_n` is a principal prefix of that same
matrix, and both displayed phases have the stated normalization.  The
48-coordinate certificate is in the correct `B tensor H_2` order and its
right-tensor propagation lands in the actual `B tensor H_(2+s)` prefix
without a coordinate permutation.  The resulting theorem is genuine
nonconvergence of the explicitly constructed sequence `(q_n)`, and the
draft correctly makes no nonconvergence claim for the minimizing values
`M_n`.

No canonical theorem file, draft, or verifier was changed in this audit.

## 1. Direct-limit and Kronecker order

Use the usual block convention for the Kronecker product and associate the
tensor power as

```math
H_{r+1}=H\otimes H_r.
```

The first `4^r` row and column indices of `H_(r+1)` have outer coordinate
zero.  Their principal block is therefore

```math
H_{00}H_r=H_r,
```

because `H_(00)=1`.  Thus the finite matrices really form a nested family.
Every pair of nonnegative integer indices belongs to some such finite
block, so their union defines a unique infinite symmetric sign matrix `S`.

Equivalently, write sufficiently long base-four expansions

```math
i=(i_1,\ldots,i_r),\qquad j=(j_1,\ldots,j_r).
```

Then

```math
S_{ij}=\prod_{t=1}^r H_{i_tj_t}.                       \tag{A.WP.1}
```

Padding both expansions on the left adds the factor `H_(00)=1`, so
(A.WP.1) is independent of the chosen length.  This digit description also
shows directly that every finite `C_n=S[:n,:n]` is a prefix of the same
matrix, not a separately chosen signing at order `n`.  Since all entries of
`S` are signs, `A_n=C_n^circ` is symmetric, hollow, and sign-valued at every
off-diagonal position.

For the other phase, an index below `3*4^r` has an `(r+1)`-digit expansion

```math
(a,alpha),\qquad a\in\{0,1,2\},\quad 0\le alpha<4^r.
```

Consequently

```math
C_{3\cdot4^r}=(H_{ab})_{0\le a,b<3}\otimes H_r
=B\otimes H_r.                                       \tag{A.WP.2}
```

This proves (WP.8) in exactly the coordinate order used by the certificate.
I also checked the literal arrays through `r=4`:

```text
H_(r+1)[:4^r,:4^r]       = H_r,
H_(r+1)[:3*4^r,:3*4^r]   = B tensor H_r.
```

There is a potentially confusing but harmless distinction here.  Passing
from one principal prefix to a larger one pads integer labels by a new
leading base-four digit, whereas the propagated witness below appends
copies of `u` on the right.  The witnesses are not claimed to be literal
prefix-extensions of one another.  What is needed, and what holds with no
permutation, is

```math
(B\otimes H_2)\otimes H_s
=B\otimes H_{2+s}.                                    \tag{A.WP.3}
```

All inner factors are the same `H`, so the tensor-power convention in the
verifier and the association in (A.WP.3) produce the same numerical matrix.

## 2. Exact `4^r` phase and hollowing

The displayed `H` is symmetric, `H^2=4I`, and `Hu=2u`.  Hence, for
`n=4^r`,

```math
H_r^2=nI,\qquad \lVert H_r\rVert_{2\to2}=\sqrt n,
\qquad H_r u^{\otimes r}=\sqrt n\,u^{\otimes r}.      \tag{A.WP.4}
```

Every Boolean vector has squared Euclidean norm `n`, so the Rayleigh bound
and the displayed Boolean eigenvector give matching inequalities:

```math
|x^TH_rx|\le n\sqrt n,
\qquad
(u^{\otimes r})^TH_r(u^{\otimes r})=n\sqrt n.         \tag{A.WP.5}
```

Thus (WP.7) is exact, including the absolute value.

For any symmetric matrix `C` and Boolean `x`, hollowing satisfies the
pointwise identity

```math
x^TC^circ x=x^TCx-\sum_i C_{ii}x_i^2=x^TCx-tr(C).     \tag{A.WP.6}
```

Here `tr(H)=0`, and hence

```math
tr(H_r)=tr(H)^r=0\qquad(r\ge1).                        \tag{A.WP.7}
```

Hollowing therefore leaves every Boolean quadratic value unchanged, not
merely the maximum.  With the project's convention

```math
Q(A)=\frac12\max_x|x^TAx|,
```

(A.WP.5)--(A.WP.7) give

```math
q_{4^r}=\frac{n^{3/2}}{2n^{3/2}}=\frac12.
```

There is no missing factor of two: `x^TAx` counts the two orientations of
each off-diagonal edge, and the definition of `Q` supplies the factor
`1/2` exactly once.

## 3. The 48-coordinate certificate

Split the displayed `z` into its three length-16 blocks `z_1,z_2,z_3` and
form their exact `H_2` cross-correlation matrix.  Direct integer
multiplication gives

```math
K=(z_i^TH_2z_j)_{i,j=1}^3
=\begin{pmatrix}
52&28&28\\
28&-44&52\\
28&52&-44
\end{pmatrix}.                                        \tag{A.WP.8}
```

All 48 displayed coordinates are in `{+-1}`.  In the standard block order,

```math
z^T(B\otimes H_2)z
=\sum_{i,j=1}^3B_{ij}K_{ij}=356.                      \tag{A.WP.9}
```

The sum over ordered pairs in (A.WP.9) is the full quadratic-form
convention used in (WP.2); no extra off-diagonal factor belongs in the
certificate.  Also

```math
tr(B\otimes H_2)=tr(B)tr(H_2)=(-1)\cdot0=0,            \tag{A.WP.10}
```

so the same value `356` is obtained after hollowing.

For propagation, put `v_s=u^(tensor s)`, with `v_0=(1)`.  Then

```math
H_sv_s=2^sv_s,
\qquad v_s^Tv_s=4^s,
\qquad v_s^TH_sv_s=8^s.                               \tag{A.WP.11}
```

Using (A.WP.3), the Boolean vector `z tensor v_s` therefore obeys

```math
\begin{aligned}
&(z\otimes v_s)^T(B\otimes H_{2+s})(z\otimes v_s)\\
&\qquad=(z^T(B\otimes H_2)z)(v_s^TH_sv_s)
=356\,8^s.                                            \tag{A.WP.12}
\end{aligned}
```

I independently evaluated (A.WP.12) for `s=0,1,2,3`; the exact values were

```text
356, 2848, 22784, 182272.
```

The trace remains zero for every claimed level because
`tr(B tensor H_(2+s))=tr(B)tr(H_(2+s))=0`.

## 4. Normalization, separation, and logical conclusion

At `r=2+s`, the order is

```math
n=3\cdot4^{2+s}=48\cdot4^s,
\qquad n^{3/2}=48^{3/2}8^s.                           \tag{A.WP.13}
```

Combining the positive witness (A.WP.12), trace removal, and the factor
`1/2` in (WP.2) gives

```math
q_{3\cdot4^{2+s}}
\ge {356\,8^s\over2\,48^{3/2}8^s}
={89\over96\sqrt3}.                                  \tag{A.WP.14}
```

The simplification is exact because `48^(3/2)=192sqrt(3)`.  The strict gap
over `1/2` is also exact:

```math
{89\over96\sqrt3}>{1\over2}
\iff 89>48\sqrt3,
\qquad 89^2-3\cdot48^2=1009>0.                        \tag{A.WP.15}
```

The indices `4^r` and `3*4^r` both tend to infinity.  The first subsequence
is identically `1/2`, while every term of the second subsequence (from
`r=2` onward) is at least the fixed larger constant in (A.WP.14).  Hence
the full explicitly defined sequence `(q_n)` cannot converge, and the
liminf/limsup inequalities in (WP.6) follow.  Exact optimization on the
high phase is neither asserted nor needed.

This logic does not transfer to `M_n`.  If `M_n` is the minimum over all
order-`n` signings, a large value for this one signing at `3*4^r` is not a
universal lower bound on that minimum.  Proving a separated high phase for
`M_n` would require such a universal lower bound; the draft explicitly
disclaims it in the title status, theorem discussion, and final scope
paragraph.

The broader fixed-phase sentence in Section 3 is also consistent with the
same ordering.  For fixed `d`, choose `k` with `d<=4^k`.  If `C_d` is the
leading `d` block of `H_k`, then

```math
C_{d4^r}=C_d\otimes H_r,
```

as the leading `d` outer blocks of `H_k tensor H_r`.  Regular-Hadamard
amplification therefore applies to each fixed phase separately.

## 5. Verifier audit

Running the supplied command produced

```text
Walsh prefix nonconvergence checks passed: 22;
high=0.535251812061216; gap=0.035251812061216
```

The substantive certificate operations use exact `int64` arithmetic, and
all tested magnitudes are far below overflow.  The script correctly checks
`H^2=4I`, `Hu=2u`, finite nesting, Boolean pole eigenvectors, zero traces,
the value `356`, two further tensor levels, and the exact integer inequality
equivalent to the radical gap.

There are three optional hardening improvements, none of which is a proof
defect in the current draft:

1. The verifier could assert `abs(Z)==1` entrywise and explicitly compare
   `powers[r+1][:3*4^r,:3*4^r]` with `kron(B,powers[r])`.  At present both
   facts are evident from the literals/proof but are not dedicated asserts.
2. It could compare each propagated matrix with the corresponding coherent
   principal prefix.  I performed those exact comparisons independently;
   they pass.
3. The pole check uses the floating expression `(4**r)**1.5`.  It is exact
   for the small powers of two tested, but `8**r` would express the intended
   integer identity more robustly.  The theorem's decisive gap already has
   a separate exact-integer assertion.

Checking only two propagated levels is appropriate for regression; the
all-level claim follows from the algebraic tensor identity (A.WP.12), not
from extrapolating finite numerical tests.

## 6. Final scope verdict

No factor-of-two, trace, sign, tensor-association, prefix, or index-order
counterexample was found.  The construction gives one all-order sequence
of dense hollow signings whose normalized energies have two uniformly
separated geometric subsequences.  Its correct promotion scope is an
explicit Walsh direct-limit nonconvergence benchmark, not a theorem about
the optimizing sequence `M_n`.

## Addendum: Theorem WP.2, continuous phase profile

This addendum audits Theorem WP.2, which was added to the draft after the
audit above began.

**Addendum verdict: PASS.**  The base-four rational prefix has the stated
outer/inner Kronecker order, the principal-deletion estimate has the right
operator norms and half-quadratic factors, floors contribute only a uniform
`O(4^(-r))` increment to the scale variable, and dense-set convergence plus
the displayed asymptotic modulus does imply uniform convergence to a
continuous function.  The final varying-mantissa statement follows from
that uniform convergence, not merely from pointwise convergence.

### A. Exact prefix at every base-four rational

Let

```math
t={p\over4^k}\in[1,4],\qquad r\ge k,\qquad s=r-k.
```

Necessarily `4^k<=p<=4^(k+1)`.  With the same standard Kronecker order used
in the first audit,

```math
H_{r+1}=H_{k+1}\otimes H_s.                           \tag{A.WP.16}
```

Each outer coordinate of `H_(k+1)` indexes a consecutive block of `4^s`
coordinates.  Therefore the leading

```math
n=p4^s=t4^r
```

rows and columns consist of the first `p` complete outer blocks, and

```math
C_n=H_{k+1}[:p,:p]\otimes H_s
=B_{p,k}\otimes H_{r-k}.                              \tag{A.WP.17}
```

Thus (WP.16) has neither an off-by-one tensor depth nor a hidden coordinate
permutation.  The endpoints also work literally:

- at `t=1`, taking `p=4^k` makes `B_(p,k)` the nested `H_k` prefix and
  (A.WP.17) becomes `H_r`;
- at `t=4`, taking `p=4^(k+1)` makes `B_(p,k)=H_(k+1)` and
  (A.WP.17) becomes `H_(r+1)`.

Nonreduced representations of a base-four rational give the same prefix
because the finite powers are nested.  I checked (A.WP.17) as an exact
array identity for 198 small triples `(p,k,r)`, covering every admissible
`p` for `k=0,1,2` and both endpoints.

For completeness, put `G_s=B_(p,k) tensor H_s` and `N_s=p4^s`.  The
regular-Hadamard theorem gives convergence of

```math
{1\over2N_s^{3/2}}\max_x|x^TG_sx|.                   \tag{A.WP.18}
```

Hollowing changes each full quadratic by the constant `tr(G_s)`, so

```math
\left|Q(G_s^circ)-{1\over2}\max_x|x^TG_sx|\right|
\le {1\over2}|tr(G_s)|\le {N_s\over2}.               \tag{A.WP.19}
```

After division by `N_s^(3/2)` this is at most
`1/(2sqrt(N_s))`.  In the present Walsh construction it is even exactly
zero for every `s>=1`, because `tr(H_s)=0`.  Finally,

```math
F_{k+s}(t)
=t^{3/2}{Q(G_s^circ)\over N_s^{3/2}},                 \tag{A.WP.20}
```

so the regular-amplification limit proves convergence of `F_r(t)` for
every base-four rational `t`, with the normalization claimed in the draft.

### B. Principal deletion and the operator-norm bound

Write the project energy as

```math
\mathcal H_A(x)={1\over2}x^TAx,
\qquad Q(A)=\max_x|\mathcal H_A(x)|.
```

For the block decomposition (WP.17),

```math
\mathcal H_{A_m}(x,y)
=\mathcal H_{A_n}(x)+x^TRy+{1\over2}y^TDy.           \tag{A.WP.21}
```

The coefficient of the cross term is one, not one half, because the two
off-diagonal blocks occur twice in the full symmetric quadratic.

Since `m<=4R_r=4^(r+1)`, all blocks come from the single ambient matrix
`H_(r+1)`.  Its norm is

```math
\lVert H_{r+1}\rVert=2^{r+1}=2\sqrt{R_r}.             \tag{A.WP.22}
```

A rectangular submatrix is `P H_(r+1) Q` for coordinate-selection
contractions `P,Q`, so its norm cannot exceed (A.WP.22).  If `E` is the
unhollowed principal block underlying `D`, then

```math
\lVert R\rVert\le2\sqrt{R_r},
\qquad
\lVert D\rVert
\le\lVert E\rVert+\lVert diag(E)\rVert
\le2\sqrt{R_r}+1,                                    \tag{A.WP.23}
```

because `diag(E)` has sign entries and hence norm one.  Boolean vectors of
lengths `n,h` have norms `sqrt(n),sqrt(h)`, giving exactly

```math
\left|x^TRy+{1\over2}y^TDy\right|
\le2\sqrt{R_rnh}+{h\over2}(2\sqrt{R_r}+1).            \tag{A.WP.24}
```

This confirms (WP.18)--(WP.19), including all factors of `1/2`.

Principal deletion is indeed lossless for the absolute Boolean energy.
Fix `x` with

```math
sigma\mathcal H_{A_n}(x)=Q(A_n),\qquad sigma\in\{+-1\}.
```

For an independent uniform Boolean vector `Y`, hollowing of `D` gives

```math
\mathbb E_Y[\sigma\mathcal H_{A_m}(x,Y)]
=\sigma\mathcal H_{A_n}(x)=Q(A_n).                   \tag{A.WP.25}
```

Some realization is at least its expectation, so `Q(A_m)>=Q(A_n)`.
Conversely, (A.WP.21)--(A.WP.24) and the triangle inequality give

```math
Q(A_m)\le Q(A_n)
+2\sqrt{R_rnh}+{h\over2}(2\sqrt{R_r}+1).              \tag{A.WP.26}
```

Dividing by `R_r^(3/2)` yields precisely (WP.20).

As a finite numerical regression, I exhaustively inspected every
consecutive split `R_r<=n<=m<=4R_r` inside `H_3` (`r=2`).  The ambient norm
was `8`; every cross-block norm was at most `8`, and every hollow new-block
norm was at most `9`, as (A.WP.23) requires.

### C. Floors and a uniform modulus

Let `delta=s-t`, `R=R_r`, and

```math
a={n\over R},\qquad b={h\over R}.
```

The floor inequalities give

```math
0\le a\le4,
\qquad
0\le b
={\lfloor sR\rfloor-\lfloor tR\rfloor\over R}
\le\delta+{1\over R}.                                \tag{A.WP.27}
```

Consequently (WP.20) implies the explicit uniform modulus

```math
|F_r(s)-F_r(t)|
\le4\sqrt{\delta+R^{-1}}
+(\delta+R^{-1})
 +{\delta+R^{-1}\over2\sqrt R}                       \tag{A.WP.28}
```

for `t<=s`; symmetry gives the same statement in terms of `|s-t|`.
For every positive epsilon, first choose `delta` small and then `r` large;
the right side becomes smaller than epsilon.  This is exactly the
asymptotic equicontinuity needed in the proof.  In particular, a one-point
floor jump has size `O(R^(-1/2))`, so the step discontinuities of `F_r`
vanish uniformly.

Uniform boundedness is immediate from deletion monotonicity and the full
Walsh endpoint:

```math
0\le F_r(t)\le F_r(4)
={Q(A_{4R})\over R^{3/2}}
={\tfrac12(4R)^{3/2}\over R^{3/2}}=4.                \tag{A.WP.29}
```

No floor or endpoint error occurs at `t=4`.

### D. Dense convergence implies uniform convergence here

The proof's finite-net step is valid even though each `F_r` is a step
function rather than a continuous function.  Given an error tolerance,
use (A.WP.28) to choose a common small mesh for all sufficiently large
`r`, and choose the finitely many mesh points from the dense set of
base-four rationals.  Pointwise convergence at those finitely many points
makes the values uniformly Cauchy on the mesh.  For an arbitrary `t`, a
nearby mesh point and (A.WP.28) on each of two large levels give

```math
|F_r(t)-F_s(t)|
\le |F_r(t)-F_r(d)|+|F_r(d)-F_s(d)|+|F_s(d)-F_s(t)|.
```

All three terms are uniformly small.  Thus `(F_r)` is uniformly Cauchy in
the sup norm and converges uniformly to some `F`.

Continuity does not follow merely from uniform convergence of these step
functions, but it does follow from the retained modulus: taking a limit in
(A.WP.28) gives, for example,

```math
|F(s)-F(t)|\le4\sqrt{|s-t|}+|s-t|.                   \tag{A.WP.30}
```

Hence `F` is continuous.  Each `F_r` is nondecreasing by (A.WP.25), so its
uniform limit is nondecreasing as well.  This supplies every detail needed
for the theorem's finite-net sentence.

### E. Denominator floors and the varying mantissa

For fixed `t`, set

```math
rho_r(t)={\lfloor tR_r\rfloor\over R_r}.
```

Uniformly on `[1,4]`,

```math
|rho_r(t)-t|\le R_r^{-1},\qquad rho_r(t)\ge1.         \tag{A.WP.31}
```

Therefore

```math
{Q(A_{\lfloor tR_r\rfloor})
 \over\lfloor tR_r\rfloor^{3/2}}
={F_r(t)\over rho_r(t)^{3/2}}
\longrightarrow {F(t)\over t^{3/2}}=L(t).            \tag{A.WP.32}
```

In fact, (A.WP.31) and the uniform convergence of `F_r` make (A.WP.32)
uniform in `t`, which is stronger than the pointwise assertion (WP.13).
Continuity of `L` follows because `t` is bounded away from zero.

At `t=1` and `t=4`, (WP.4) gives `L(1)=L(4)=1/2`.  At `t=3`, the uniform
limit exists and (WP.5) gives

```math
L(3)\ge {89\over96\sqrt3}>{1\over2}.
```

Thus the endpoint identification is compatible with a continuous
base-four mantissa circle, while the interior phase is nonconstant.

Finally, for an arbitrary integer sequence `n->infinity`, put

```math
r=\lfloor\log_4n\rfloor,
\qquad R=4^r,
\qquad t_n={n\over R}\in[1,4).
```

Here there is no floor error at all: `floor(t_n R)=n`.  Hence

```math
q_n={F_r(t_n)\over t_n^{3/2}},
\qquad
|q_n-L(t_n)|
\le\lVert F_r-F\rVert_\infty\longrightarrow0.        \tag{A.WP.33}
```

This verifies the draft's final phase-profile conclusion for a varying
mantissa.  It is precisely the uniform, rather than merely pointwise,
convergence proved above that licenses (A.WP.33).

### F. Verifier coverage

The current `verify_walsh_prefix_nonconvergence.py` is unchanged from the
WP.1 version.  It still passes its 22 exact checks, but it contains no
dedicated assertions for WP.2's rational-prefix identity, floor modulus,
or deletion bounds.  The reproducibility paragraph accurately lists only
the WP.1 checks and does not claim otherwise.  This is a computational
coverage gap, not a flaw in WP.2's self-contained argument.  The independent
exact prefix and submatrix-norm regressions reported in this addendum passed.

No edit to the theorem draft or verifier was made during this second audit.
