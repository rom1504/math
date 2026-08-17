# Approximate Boolean recovery forces amplification-carrier limits

Status: rigorous task-local theorem draft.  This extends the mechanism in
`regular_hadamard_amplification_limit.md`; it is not a statement about the
minimum defining `M_n` and does not modify the canonical files.

## 1. The metric statement that is not yet a recovery theorem

For nonempty compact sets in a metric space write

```math
e(K,L)=\sup_{x\in K}\operatorname{dist}(x,L)
```

for directed Hausdorff excess.  Thus
`d_H(K,L)=max(e(K,L),e(L,K))`.

If

```math
d_H(K_r,K_{r+1})\leq\delta_r,
\qquad \sum_r\delta_r<\infty,                         \tag{AR.1}
```

then the triangle inequality makes `(K_r)` Hausdorff Cauchy, and completeness
of the compact-set hyperspace gives a limit with

```math
d_H(K_r,K_\infty)\leq\sum_{s\geq r}\delta_s.          \tag{AR.2}
```

This is useful bookkeeping but is tautological: it assumes the two-sided
response comparison that an amplification theorem is supposed to produce.
Exact Hadamard amplification does not normally satisfy
`K_{r+1} subset K_r`; new Boolean response points are allowed to appear.

The one-sided version is less immediate and is the compactness engine used
below.

### Lemma AR.1 (summable directed recovery)

Let `X` be a compact metric space and let `K_r` be nonempty compact subsets of
`X`.  If

```math
e(K_r,K_{r+1})\leq\epsilon_r,
\qquad \sum_r\epsilon_r<\infty,                       \tag{AR.3}
```

then `K_r` converges in Hausdorff distance to a nonempty compact set
`K_infinity`.

Only the old-to-limit direction has a defect-only rate:

```math
e(K_r,K_\infty)\leq\sum_{s\geq r}\epsilon_s.          \tag{AR.4}
```

There is no universal bound on `e(K_infinity,K_r)` that tends to zero as a
function only of the displayed tail.  The tail may even vanish identically
while the reverse excess remains positive, because new points can be
introduced at arbitrarily late levels.

#### Proof

Starting from any `x_r in K_r`, choose successively `x_{s+1} in K_{s+1}`
with `d(x_s,x_{s+1})<=epsilon_s`.  The resulting sequence is Cauchy and has a
limit in `X`.  Let `E` be the set of all limits of all such forward chains,
started at every level and every point, and put `K_infinity=closure(E)`.
Every `x_r` has a forward-chain limit within the tail in (AR.4), proving that
inequality.

It remains to prove the other directed excess tends to zero.  Given `eta>0`,
take a finite `eta/3`-net `z_1,...,z_m` of `K_infinity` with each `z_i in E`;
this is possible because `E` is dense in the compact set `K_infinity`.
Let `s_i` be a starting level of a chain ending at `z_i`.  Once
`r>=max_i s_i`, that chain has a point of `K_r` within
`sum_(s>=r) epsilon_s` of `z_i`.  For all sufficiently large `r` this tail is
at most `eta/3`, so every point of `K_infinity` is within `2eta/3` of `K_r`.
Thus `e(K_infinity,K_r)->0`. `square`

Convexity is not used in this lemma.  Compactness, rather than mere
boundedness in a varying-dimensional space, is essential.

## 2. An intrinsic Boolean lifting and retraction theorem

Let

```math
Omega_r={+-1}^{n_r},
\qquad
kappa_r(x,y)={x^TA_ry\over a_r},                      \tag{AR.5}
```

where `A_r` is a real symmetric `n_r by n_r` matrix and `a_r>0`.  For fixed
outer dimension `d`, define the raw response matrix and its convex carrier by

```math
\begin{aligned}
Phi_r^{(d)}(x_1,...,x_d)
  &=\bigl(kappa_r(x_i,x_j)\bigr)_{1\leq i\leq j\leq d},\\
K_r^{(d)}&=\operatorname{conv} Phi_r^{(d)}(Omega_r^d).
\end{aligned}                                         \tag{AR.6}
```

The ambient norm below is entrywise `l_infinity`.  Assume the intrinsic
kernel bound

```math
\sup_r\sup_{x,y\in Omega_r}|kappa_r(x,y)|\leq M<\infty. \tag{AR.7}
```

A **Boolean forward recovery** is any map
`L_r:Omega_r -> Omega_(r+1)`.  Its pair distortion is

```math
alpha_r=\sup_{x,y\in Omega_r}
 |kappa_{r+1}(L_rx,L_ry)-kappa_r(x,y)|.                \tag{AR.8}
```

A **Boolean backward recovery** is a map
`P_r:Omega_(r+1) -> Omega_r`, with distortion

```math
beta_r=\sup_{z,w\in Omega_{r+1}}
 |kappa_r(P_rz,P_rw)-kappa_{r+1}(z,w)|.                \tag{AR.9}
```

These notions depend only on the finite Boolean response kernels, not on a
chosen matrix presentation or on `d`.

### Theorem AR.2 (Boolean shadowing and quantitative retraction)

Under (AR.7), the following hold simultaneously for every fixed `d`.

1. Every forward recovery gives

   ```math
   e_infinity(K_r^{(d)},K_{r+1}^{(d)})\leq alpha_r.    \tag{AR.10}
   ```

   Consequently, if `sum_r alpha_r<infinity`, then `K_r^{(d)}` has a
   compact convex Hausdorff limit.

2. If backward recoveries also exist, then

   ```math
   d_H^infinity(K_r^{(d)},K_{r+1}^{(d)})
   \leq\max(alpha_r,beta_r).                           \tag{AR.11}
   ```

   If the right sides are summable, the limit has the quantitative bound

   ```math
   d_H^infinity(K_r^{(d)},K_infinity^{(d)})
   \leq\sum_{s\geq r}\max(alpha_s,beta_s).            \tag{AR.12}
   ```

3. For every linear query `theta` on the symmetric response coordinates,
   the support values converge.  The convergence is uniform on sets with
   bounded dual `l_1` norm.  The same holds for absolute support values.

If additionally `P_r L_r` is the identity, the pair is literally a
lift/retraction.  That identity is not needed for (AR.10)--(AR.12); the two
kernel-distortion inequalities are the exact response obligations.

#### Proof

Given a raw response point from `(x_1,...,x_d)`, apply the *same* map `L_r`
to every block.  Every coordinate changes by at most `alpha_r`.  Given a
convex combination of raw points, use the corresponding convex combination
of their lifted raw points.  This proves (AR.10).  Bound (AR.7) places every
carrier in the fixed compact cube `[-M,M]^{d(d+1)/2}`.  Lemma AR.1 now proves
the first conclusion, and convexity passes to Hausdorff limits.

Applying `P_r` blockwise proves the reverse directed excess bound by the same
argument, hence (AR.11).  The Hausdorff-Cauchy argument (AR.1)--(AR.2) gives
(AR.12).  Finally,

```math
|h_K(theta)-h_L(theta)|
\leq ||theta||_1 d_H^infinity(K,L),                   \tag{AR.13}
```

and `max |theta dot x|=max(h_K(theta),h_K(-theta))`.
This proves the support claims. `square`

The common-map and pairwise requirements matter.  Separate witness-dependent
lifts do not preserve a full response matrix, and preservation only of
`kappa_r(x,x)` does not control cross-block queries.

## 3. A checkable matrix-level certificate

Call a matrix

```math
T_r in {0,+-1}^{n_{r+1} by n_r}                      \tag{AR.14}
```

a signed replication matrix if every row has exactly one nonzero entry.
Then `T_rx` is Boolean whenever `x` is.  Put `C_r=A_r/a_r` and

```math
D_r^up=T_r^TC_{r+1}T_r-C_r.                           \tag{AR.15}
```

For a real matrix `D`, use

```math
||D||_(infinity->1)
=\max_{x,y\in{+-1}^n}|x^TDy|.                        \tag{AR.16}
```

### Theorem AR.3 (compressed-kernel certificate)

The lift `L_r x=T_rx` has exactly

```math
alpha_r=||D_r^up||_(infinity->1).                     \tag{AR.17}
```

In particular, either computable bound

```math
alpha_r
\leq ||D_r^up||_(entrywise 1)
\quad\hbox{or}\quad
alpha_r\leq n_r||D_r^up||_(2->2)                     \tag{AR.18}
```

is a sufficient certificate.  If (AR.7) holds and the chosen upper bounds
are summable, all fixed-`d` carriers converge.

Likewise, let

```math
S_r in {0,+-1}^{n_r by n_{r+1}}
```

have one nonzero per row, and put

```math
D_r^down=S_r^TC_rS_r-C_{r+1}.                         \tag{AR.19}
```

Then `P_rz=S_rz` has

```math
beta_r=||D_r^down||_(infinity->1),                    \tag{AR.20}
```

with the analogous entrywise and operator-norm bounds (using dimension
`n_(r+1)`).  If `S_rT_r=I`, these maps are an actual linear Boolean
lift/retraction.

#### Proof

For Boolean `x,y`, the error after lifting is exactly `x^TD_r^up y`.
Maximizing gives (AR.17).  Also

```math
|x^TDy|\leq\sum_{ij}|D_{ij}|,
\qquad
|x^TDy|\leq ||x||_2||D||_(2->2)||y||_2
=n||D||_(2->2),                                      \tag{AR.21}
```

which proves (AR.18).  The backward statement is identical. `square`

The exact `infinity->1` norm is generally a hard cut-norm computation; the
two bounds in (AR.18) are deliberately conservative, directly checkable
certificates.  Compression is important: perturbations `E` satisfying
`T_r^TET_r=0` cost no forward-recovery error even if their full
operator norm is large.

## 4. Non-exact Hadamard/tensor amplification

Now take the natural quadratic scale `a_r=n_r^{3/2}`.  Suppose

```math
n_{r+1}=h_r n_r,
\qquad
A_{r+1}=A_r\mathbin{\mathop\otimes}H_r+E_r,            \tag{AR.22}
```

where `H_r` and `E_r` are symmetric, and choose
`u_r in {+-1}^{h_r}`.  Define

```math
rho_r={u_r^TH_ru_r\over h_r^{3/2}},
\qquad
M_r={||A_r||_(2->2)\over\sqrt{n_r}},
\qquad
e_r={||E_r||_(2->2)\over\sqrt{n_{r+1}}}.              \tag{AR.23}
```

### Corollary AR.4 (approximately regular tensor factors)

If

```math
\sup_r M_r<\infty,
\qquad
\sum_r\bigl(|rho_r-1|+e_r\bigr)<\infty,               \tag{AR.24}
```

then every fixed-dimensional Boolean quadratic carrier `K_r^(d)` converges,
and hence so does every fixed outer quadratic support query.

Explicitly, for a fixed symmetric outer matrix `B`, let `theta_B` have
diagonal coordinates `B_(ii)` and off-diagonal coordinates `2B_(ij)`.  Then

```math
{1\over2(dn_r)^{3/2}}\max_{x\in Omega_r^d}
x^T(B\mathbin{\mathop\otimes}A_r)x
={1\over2d^{3/2}}h_{K_r^{(d)}}(theta_B),              \tag{AR.24a}
```

so both the upper and absolute normalized Boolean quadratic maxima converge.

A vector-residual certificate for the first summand is

```math
|rho_r-1|
\leq {||H_ru_r-\sqrt{h_r}u_r||_2\over h_r}.           \tag{AR.25}
```

The uniform bound on `M_r` is itself local-checkable.  For example, it follows
if

```math
sigma_r\geq0,
\qquad ||H_r||_(2->2)\leq(1+sigma_r)\sqrt{h_r},
\qquad \sum_r sigma_r<\infty,
\qquad \sum_r e_r<\infty.                             \tag{AR.26}
```

#### Proof

Use the signed replication `T_r=I_(n_r) tensor u_r`.  Direct compression
gives

```math
D_r^up
=(rho_r-1){A_r\over n_r^{3/2}}
+{T_r^TE_rT_r\over n_{r+1}^{3/2}}.                   \tag{AR.27}
```

The first term has `infinity->1` norm at most `M_r|rho_r-1|`.  Since
`||T_rx||_2^2=n_{r+1}` for every Boolean `x`, the second has Boolean bilinear
norm at most `e_r`.  Hence

```math
alpha_r\leq M_r|rho_r-1|+e_r,                         \tag{AR.28}
```

which is summable under (AR.24), while the same `M_r` bound proves (AR.7).
Theorem AR.2 applies.  Cauchy--Schwarz proves (AR.25).  Finally,

```math
M_{r+1}\leq(1+sigma_r)M_r+e_r,                        \tag{AR.29}
```

so the standard finite-product bound proves the claim following (AR.26).
`square`

Exact regular Hadamard amplification has `rho_r=1`, `E_r=0`, and
`M_r=M_0`, so this recovers the earlier nesting theorem.  Positivity is not
cosmetic: `rho_r` close to `-1` reverses the signed carrier.  Absolute scalar
support values may survive such a reversal, but convergence of the signed
response carrier need not.

### Corollary AR.5 (explicit scalar quasi-monotonicity)

Under the hypotheses and notation of Corollary AR.4, put

```math
\begin{aligned}
p_r^+&={1\over2n_r^{3/2}}\max_{x\in Omega_r}x^TA_rx,\\
p_r^abs&={1\over2n_r^{3/2}}\max_{x\in Omega_r}|x^TA_rx|,\\
gamma_r&=M_r|rho_r-1|+e_r.
\end{aligned}                                         \tag{AR.29a}
```

Then

```math
p_{r+1}^+\geq p_r^+-{gamma_r\over2},
\qquad
p_{r+1}^abs\geq p_r^abs-{gamma_r\over2},              \tag{AR.29b}
```

and both scalar sequences converge.

Indeed, lift an upper optimizer, or an absolute optimizer followed by the
reverse triangle inequality, and use (AR.28) with `x=y`.  Both sequences are
bounded above by `M_r/2`.  In general, if
`z_(r+1)>=z_r-epsilon_r`, the sequence is bounded above, and
`sum epsilon_r<infinity`, then

```math
z_r-\sum_{j\geq r}\epsilon_j                          \tag{AR.29c}
```

is nondecreasing and bounded above.  It converges, and adding back the
vanishing tail proves convergence of `z_r`.

## 5. A genuinely non-tensor dense-sign hierarchy

The operator certificate is realized by constrained dense sign matrices, not
only by arbitrary real perturbations.

### Corollary AR.6 (perfect-matching flips)

Fix a symmetric regular Hadamard matrix `H` of nontrivial order `h>1`, with
`Hu=sqrt(h)u` for a Boolean `u`.  (In particular, `h` is even.)  Start from
any symmetric full `+-1` matrix `C_0` of order `N_0`.  Having constructed
`C_r`, form `C_r tensor H` and reverse the two symmetric off-diagonal signs
on every edge of a fixed-point-free perfect matching.  Call the resulting
full sign matrix `C_(r+1)`, where `N_(r+1)=hN_r`.

Then the fixed-port carriers of `C_r` converge.  If
`A_r=C_r^circ` is obtained by zeroing the diagonal, it is a dense hollow sign
matrix and

```math
{Q(A_r)\over N_r^{3/2}}                               \tag{AR.29d}
```

converges.  Each level differs from exact tensor amplification on
`N_(r+1)/2=Theta(N_r)` undirected edges.

#### Proof

The perturbation `E_r=C_(r+1)-C_r tensor H` is, after a permutation, a direct
sum of blocks

```math
\begin{pmatrix}0&\pm2\\\pm2&0\end{pmatrix},           \tag{AR.29e}
```

with the same sign in the two off-diagonal positions of each block.

Therefore `||E_r||_(2->2)=2` and

```math
e_r={2\over\sqrt{N_{r+1}}}.                           \tag{AR.29f}
```

This is summable because `N_r=N_0h^r`.  Exact regularity gives `rho_r=1`,
and (AR.29) gives the required uniform spectral bound.  Corollary AR.4 proves
carrier convergence.

For the hollowing statement, every Boolean `x` satisfies

```math
x^TC_rx=x^TA_rx+\operatorname{tr}(C_r).               \tag{AR.29g}
```

Since `Q(A_r)=(1/2)max_x|x^TA_rx|` and
`|tr(C_r)|<=N_r`,

```math
\left|{Q(A_r)\over N_r^{3/2}}-p_r^abs(C_r)\right|
\leq {1\over2\sqrt{N_r}}\longrightarrow0.            \tag{AR.29h}
```

Corollary AR.5 finishes the proof. `square`

The condition `h>1` is required here: it both supplies even growing orders
for the perfect matchings and makes the hollowing correction vanish.

## 6. Sharp boundaries and small falsifiers

### 6.1 Vanishing, nonsummable lift error is insufficient

Let `n_r=2^r`, `c_r=sin(log(r+2))`, and

```math
A_r=c_r\sqrt{n_r}I_{n_r},
\qquad a_r=n_r^{3/2}.                                 \tag{AR.30}
```

Duplicate every coordinate, `L_rx=x tensor (1,1)`.  Then

```math
alpha_r=|c_{r+1}-c_r|\longrightarrow0,                \tag{AR.31}
```

but `K_r^(1)={c_r}` does not converge.  Rounding
`exp(pi/2+2pi k)` and `exp(3pi/2+2pi k)` gives subsequences tending to `1`
and `-1`.  Moreover, the Lipschitz bound for sine gives
`alpha_r<=log((r+3)/(r+2))<=1/(r+2)`, so the errors are square summable.
Thus `alpha_r->0`, square summability, or informal
"asymptotic regularity" cannot replace a finite total drift without an
additional cancellation mechanism.

### 6.2 Preserving self-quadratics does not preserve the carrier

At order two and scale one, take

```math
A=0,
\qquad
\widetilde A=\operatorname{diag}(1,-1).               \tag{AR.32}
```

The identity map preserves `x^TAx=x^T\widetilde Ax=0` for every Boolean
`x`.  But for `x=(1,1)` and `y=(1,-1)`,

```math
x^TAy=0,
\qquad x^T\widetilde Ay=2.                            \tag{AR.33}
```

Thus diagonal-only optimizer recovery controls the scalar `d=1` query but
not the query-complete `d>=2` response carrier.  Pair distortion in (AR.8)
is the necessary hypothesis for the present proof.

### 6.3 Absolute regularity does not orient a signed carrier

With `n_r=1`, `A_r=(-1)^r`, and the identity Boolean map, successive kernels
are exact negatives.  Their absolute scalar maxima agree, but
`K_r^(1)={(-1)^r}` oscillates.  A hypothesis `|rho_r|=1` therefore cannot
replace `rho_r` close to `+1` in a signed-carrier theorem.

### 6.4 One-sided recovery has no defect-only convergence rate

Let positive `b_j` have finite sum `B`, set

```math
C_j=b_j\begin{pmatrix}-1&1\\1&-1\end{pmatrix},
\qquad
A_r=(0)\mathbin{\mathop\oplus}C_1\mathbin{\mathop\oplus}\cdots
       \mathbin{\mathop\oplus}C_r,                   \tag{AR.34}
```

and use scale one.  The lift appends `(1,1)`.  Since
`(1,1)C_j(1,1)^T=0`, it preserves every pair response exactly, so
`alpha_r=0`.  Nevertheless

```math
K_r^(1)=\left[-4\sum_{j\leq r}b_j,0\right],
\qquad
d_H(K_r^(1),K_\infty^(1))=4\sum_{j>r}b_j.             \tag{AR.35}
```

The tail can decay arbitrarily slowly.  A backward recovery, or another
quantitative bound on late innovations, is genuinely needed for (AR.12).
Demanding an exact backward recovery in the regular-Hadamard theorem would,
however, be too strong: it would force equality of consecutive carriers,
contradicting known strict finite amplification examples.

## 7. Scope and next falsifier

The theorem proves a whole-sequence limit for each fixed outer dimension
from a common, all-pairs Boolean lift.  It tolerates non-tensor and
non-Hadamard perturbations visible through the compressed kernel, and it
does not require optimizer consistency.  It does not compute the limit,
give an innovation-side convergence rate from a forward lift alone, transfer
minimizers between unrelated orders, show that the matching-flip hierarchy
is extremally competitive, or tolerate nonsummable spectral drift without an
additional cancellation hypothesis.

The next useful falsifier should target any proposed weakening of finite
total drift.  The diagonal drift family (AR.30) already rules out bare
`alpha_r->0` and even `sum alpha_r^2<infinity`.  To go weaker than `l_1`, one
must state a real cancellation or martingale hypothesis and test it against
the same slowly rotating scalar mode embedded in a growing Boolean kernel.
For matrix sharpness in the other direction, test perturbations with
`T_r^TE_rT_r=0` and large full operator norm: they leave forward recovery
exact and show that the convenient `||E_r||/sqrt(n_(r+1))` condition is
sufficient, not necessary.

The accompanying verifier is
`experiments/verify_approximate_boolean_recovery_amplification_limit.py`.
