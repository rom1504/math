# Global non-diagonal quadratic forms: exact mapping and scale audit

Status: **Verified exact formulation and exact fixed-order Schatten-moment
obstruction; no improvement to the `0.336493364431...` lower bound and no
composition recurrence.**  The non-diagonal form below retains the same Boolean switch on
both sides, but its first and second matrix moments live at scale `n`, whereas
the project asks for an extreme eigenvalue at scale `n^(3/2)`.  The recent
zeta-function rank argument is genuine, but its decisive rank-versus-moment
step does not bridge that scale gap.

Throughout, let

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad M(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|.
\tag{1}
```

Thus the repository's `M_n` is `min_A M(A)`.

## 1. What the recent zeta argument actually supplies

The primary sources are Anthropic's
[research note](https://www.anthropic.com/research/riemann-zeta), the
[full paper](https://www-cdn.anthropic.com/564f962e60643842f5fcb4a17c9dbc8f608f1c37.pdf),
and the accompanying
[five-page expert summary](https://www-cdn.anthropic.com/23455459f8832d06bb175cc0f88d019aed962ef8.pdf).
The mathematical mechanism is more specific than the phrase "use a
non-diagonal quadratic form."

The paper compresses Weil's Hermitian form to a finite Gabor-type function
space.  In that compression:

- a zero on the critical line contributes a positive rank-one square;
- a functional-equation pair off the line contributes a `2 x 2` hyperbolic
  block of signature `(1,1)`; and
- an explicit-formula calculation on the prime side supplies the trace and
  Hilbert--Schmidt norm of the complete, non-diagonal matrix.

Its Lemma 3.2 is a finite-dimensional statement.  If `P` is positive
semidefinite with rank at most `r`, and a Hermitian `Q` has at most `b`
positive eigenvalues, then, for `c>0`,

```math
\|P+Q\|_F^2
\ge c\,\operatorname{tr}P-{c^2r\over4}
   +2c\,\operatorname{tr}Q-c^2b.                 \tag{2}
```

The source then uses (2) to lower-bound the rank contributed by the desired
zeros.  The important inputs are not non-diagonality alone.  They are the
canonical hyperbolic decomposition and an independent dual calculation of
both low matrix moments at the same bulk scale as the rank being counted.

The reported agentic/Astra-style results are useful methodological evidence
for searching such a form, but no theorem from those reports maps to (1).
Only the proved rank lemma (2) is imported below.

## 2. An exact same-switch non-diagonal form

Let `Omega={+1,-1}^n` with uniform measure and let
`chi_S(x)=prod_(i in S)x_i` be its Walsh characters.  Multiplication by the
Hamiltonian is the self-adjoint operator

```math
(\mathcal M_Af)(x)=H_A(x)f(x)                    \tag{3}
```

on `L^2(Omega)`.  Its spectrum is exactly the multiset of Boolean energies,
so

```math
\|\mathcal M_A\|=M(A).                            \tag{4}
```

For a set of levels `L subset {0,...,n}`, put

```math
V_L=\operatorname{span}\{\chi_S:|S|\in L\},
\qquad T_L=P_{V_L}\mathcal M_AP_{V_L}.            \tag{5}
```

In the Walsh basis this is the explicit real symmetric matrix

```math
(T_L)_{S,T}=
\begin{cases}
a_{ij},&S\mathbin\triangle T=\{i,j\},\\
0,&\text{otherwise}.
\end{cases}                                      \tag{6}
```

Thus (6) is not a spectral norm or a bilinear row/column relaxation.  Each
matrix entry is a genuine same-switch edge interaction.

More generally, for any orthonormal test functions `f_1,...,f_d`, with
`v(x)=(f_1(x),...,f_d(x))`, the associated compression has the exact
decomposition

```math
T=\mathbb E_x H_A(x)v(x)v(x)^{\mathsf T},
\qquad
T-\theta I
=\mathbb E_x(H_A(x)-\theta)v(x)v(x)^{\mathsf T}.  \tag{7}
```

Consequently

```math
n_+(T-\theta I)
\le \#\{x:H_A(x)>\theta\}.                        \tag{8}
```

This is the exact Boolean counterpart of the zeta positive-index map.  The
difference is that configurations below the threshold contribute negative
squares, not canonical hyperbolic blocks carrying separate arithmetic
meaning.  Establishing a positive eigenvalue above `theta` is already a
weighted high-energy witness.

### 2.1 Exact layer comparison theorem

For one level `L={k}`, write `T_k=T_{\{k\}}` and `d_k=binom(n,k)`.  For a
Boolean `x`, set `v_x(S)=chi_S(x)` for `|S|=k`.  Direct counting in (6) gives

```math
{v_x^{\mathsf T}T_kv_x\over\|v_x\|_2^2}
=\rho_{n,k}H_A(x),
\qquad
\rho_{n,k}={2k(n-k)\over n(n-1)}.                 \tag{9}
```

Combining (4), compression, and (9) proves the exact sandwich

```math
\boxed{\rho_{n,k}M(A)\le\|T_k\|\le M(A).}        \tag{10}
```

This retains the right Boolean orbit but, at a central level,
`rho_(n,k)=1/2+o(1)`.  It therefore has a constant landing loss if used by
itself.

There is an exact nearly lossless version.  Here `w` means the **number of
levels**, not their integer span.  Let
`L_w={k_0,k_0+2,...,k_1}` contain `w` consecutive levels of one parity,
centered at `n/2`, so `k_1=k_0+2(w-1)`.  If
`d_w=sum_(k in L_w)binom(n,k)`, then

```math
\rho_{L_w}
=1-{{n-2\choose k_0-2}+{n-2\choose k_1}\over d_w},
\qquad
\rho_{L_w}M(A)\le\|T_{L_w}\|\le M(A).            \tag{11}
```

Here and below a binomial coefficient is zero when its lower argument is
outside its natural range.  To see the boundary term directly, fix an edge
`{i,j}`.  Among all subsets in the band, toggling that edge remains in the
band except when a bottom-level subset contains both endpoints or a
top-level subset contains neither.  Those two classes have sizes
`binom(n-2,k_0-2)` and `binom(n-2,k_1)`, respectively.  This proves (11)
without an asymptotic approximation.

For `w -> infinity` and `w=o(sqrt(n))`, central-binomial comparison gives
`rho_(L_w)=1-O(1/w)`.  Hence this banded non-diagonal form recovers `M(A)` to
relative `o(1)`.  Its dimension, however, is

```math
d_w=(1+o(1))w{n\choose\lfloor n/2\rfloor}
=2^{n-o(n)}.                                      \tag{12}
```

So (11) is an exact mapping, not a lower-complexity reduction.

More quantitatively, on the relevant class `M(A)=O(n^(3/2))`, (11) implies

```math
0\le M(A)-\|T_{L_w}\|
\le (1-\rho_{L_w})M(A)=O(n^{3/2}/w)=o(n^{3/2}).   \tag{12a}
```

Therefore an additive `o(n^(3/2))` evaluation or composition theorem for the
band norm would give one for the original cap, and conversely the original
cap bounds the band norm.  The near-isometry has not discarded the
project-scale obligation: it has embedded a test vector `v_x` for every
Boolean state into an exponentially dimensional matrix.  A reduction would
still need an additional polynomial-complexity description of its extreme
eigenvalue or bridge update.

## 3. The moment scale obstruction

The matrix in (6) has zero trace.  Counting its nonzero entries gives another
exact identity:

```math
\operatorname{tr}T_L^2
=d_L{n\choose2}\rho_L.                            \tag{13}
```

Here `rho_L` is the fraction of subsets in `L` for which toggling one fixed
edge remains in `L`; permutation symmetry makes it independent of the edge.
It is (9)'s `rho_(n,k)` for one level and (11)'s `rho_(L_w)` for the parity
band.

Therefore trace and Hilbert--Schmidt norm alone certify only

```math
\|T_L\|\ge
\left({\operatorname{tr}T_L^2\over d_L}\right)^{1/2}
=\sqrt{{n\choose2}\rho_L}=\Theta(n),             \tag{14}
```

a factor `sqrt(n)` below the needed scale.  In particular, for every
`theta>0`, `tr(T_L-theta I)=-theta d_L`; the thresholded
Cauchy--Schwarz/rank step used in the zeta proof cannot even start from these
first two moments at `theta=Theta(n^(3/2))`.

### 3.1 Exact fourth moment

The fourth moment does retain genuine global signing information, but still
at the wrong scale.  Define the signed four-cycle sum

```math
Z_4(A)=\sum_{C\cong C_4}\prod_{e\in C}a_e,
\qquad
\operatorname{tr}A^4=n(n-1)(2n-3)+8Z_4(A).       \tag{15}
```

An exact closed-walk enumeration gives

```math
\operatorname{tr}T_k^4=B_{n,k}+c_{n,k}Z_4(A),    \tag{16}
```

where

```math
c_{n,k}=32{n-4\choose k-2}
 +8{n-4\choose k-1}+8{n-4\choose k-3}.           \tag{17}
```

Here `B_(n,k)` is signing-independent.  Explicitly, with

```math
\theta_j=(k-j)(n-k-j)-j,
\qquad m_j={n\choose j}-{n\choose j-1},           \tag{18}
```

one may take

```math
B_{n,k}=\sum_{j=0}^{\min(k,n-k)}m_j\theta_j^4
       -3c_{n,k}{n\choose4}.                      \tag{19}
```

Equations (16)--(19) follow by comparing with the all-positive Johnson graph.
For one fixed four-cycle, the valid ordered edge uses number `8,32,8`
according as the starting subset contains `1,2,3` of its vertices, which is
exactly (17).

The coefficient (17) and the boundary formula in (11) were independently
verified by exact integer construction for four random signings at every
order `4<=n<=9`: 180 `(n,k,A)` trace-four cases and 448 parity-band cases,
with three independent Boolean-spin Rayleigh checks per band.  The verifier
is [`computations/verify_global_layer_formulas.py`](../computations/verify_global_layer_formulas.py),
and its compact certificate is
[`computations/results/global_layer_formulas_verification.json`](../computations/results/global_layer_formulas_verification.json).

Since `|Z_4(A)|<=3 binom(n,4)`, (16) yields

```math
\left({\operatorname{tr}T_k^4\over d_k}\right)^{1/4}
=O(n).                                            \tag{20}
```

This is not peculiar to order four.  A closed walk of length `2q` in (6)
uses a multigraph in which every active vertex has even degree, hence at most
`2q` active vertices.  For every fixed `q`, uniformly in the signing,

```math
0\le {\operatorname{tr}T_k^{2q}\over d_k}
=O_q(n^{2q}),
\qquad
\left({\operatorname{tr}T_k^{2q}\over d_k}\right)^{1/(2q)}
=O_q(n).                                          \tag{21}
```

The precise scope of this ceiling is important.  With normalized trace
`tau_k(X)=tr(X)/d_k`, (21) says

```math
\|T_k\|_{L^{2q}(\tau_k)}
:=\tau_k(|T_k|^{2q})^{1/(2q)}=O_q(n)              \tag{21a}
```

for each **fixed** `q`, uniformly in `A` and for `k` in any range.  Hence the
direct Schatten certificate
`||T_k|| >= ||T_k||_(L^(2q)(tau_k))`, moment-ratio certificates made from a
fixed list of these normalized even traces, and the zeta-style
trace/Hilbert--Schmidt step applied to this unshifted homogeneous compression
all remain on the `O_q(n)` bulk scale.  At a positive threshold
`theta=Theta(n^(3/2))`, the first moment of `T_k-theta I` is negative, so the
thresholded rank--trace argument has no positive trace surplus.

This does **not** prove that every imaginable inference using inertia plus
extra signing-specific structure is impossible, nor does it cover
`q=q(n)->infinity`.  It proves a ceiling for the direct fixed-order
trace-polynomial/Schatten transplant.  Inertia by itself is separately
falsified as a magnitude detector in Section 3.2.  A growing-order tail
theorem or another structural statistic could escape the ceiling, but that
would be the genuinely new ingredient rather than an application of the
zeta rank lemma alone.

### 3.2 Self-complementarity kills inertia as a magnitude detector

Suppose a permutation `P` satisfies `PAP^T=-A`.  Its induced permutation
`U_k` on `k`-subsets obeys

```math
U_kT_k(A)U_k^{\mathsf T}=-T_k(A).                 \tag{22}
```

Hence the nonzero spectrum of `T_k` is exactly symmetric and its positive
and negative indices are equal.  This applies both to exact
self-complementary minimizers and to the Wigner-scale self-complementary
ensemble already constructed in the ledger.  Inertia therefore cannot
distinguish a conference-like matrix from a large trace-four defect within
that class.  Even in the full multiplication space, the paired eigenvalues
`(+H,-H)` form perfect hyperbolic planes, but their number contains no
information about the magnitude of `H`.

This pinpoints the mismatch with the zeta proof: there a bulk fraction of
objects is the target and first/second moments see that fraction.  Here the
target is an exponentially thin extreme, `sqrt(n)` standard deviations above
the bulk, while the hyperbolic/inertia data see only the sign pairing.

## 4. Finite stress test

The central layer `k=floor(n/2)` was built exactly from (6).  Caps were
enumerated exactly; sparse symmetric eigensolves were used only for the
displayed operator norms.  The second- and fourth-moment columns are exact
integer traces before taking roots.

| signing | `n` | `M(A)` | `dim T_k` | `||T_k||` | `(tr T_k^2/d)^(1/2)` | `(tr T_k^4/d)^(1/4)` |
|---|---:|---:|---:|---:|---:|---:|
| exact `M_6` conference | 6 | 5 | 20 | 4.2361 | 3 | 3.3701 |
| exact `M_8`, class 0 | 8 | 10 | 70 | 7.1126 | 4 | 4.7987 |
| exact `M_8`, class 1 | 8 | 10 | 70 | 7.1542 | 4 | 4.7987 |
| exact nonconference `M_10` minimizer | 10 | 13 | 252 | 10.1911 | 5 | 6.1712 |
| `GF(9)` conference | 10 | 15 | 252 | 10.0623 | 5 | 6.0170 |
| exact `M_12` witness | 12 | 18 | 924 | 12.8773 | 6 | 7.3995 |
| exact `M_14` conference witness | 14 | 21 | 3432 | 15.8215 | 7 | 8.6538 |
| self-complementary Wigner sample | 8 | 16 | 70 | 9.3808 | 4 | 5.4487 |
| self-complementary Wigner sample | 12 | 28 | 924 | 17.8054 | 6 | 7.6842 |

The two order-eight minimizer classes already have identical second and
fourth layer moments but different layer norms.  For the aligned
self-complementary samples the positive/kernel/negative inertia counts are
`28/14/28` at order eight and `442/40/442` at order twelve, exactly as (22)
predicts.  Their normalized spectral fourth moments are respectively
`tr(A^4)/n^3=2.453125` and `1.645833...`, versus
`(n-1)^2/n^2=0.81` for the order-ten conference matrix.  The larger global
defect changes the constant in (20), but not its order.

The exploratory finite-table reproducer is
`/home/math/quadra/tmp/audit_global_layer_form.py` (SHA-256
`23fa981676a56c243669a744e7914fcd24bb1d0a41e875e5cbb7b4595e8c02ad`);
its output is `/home/math/quadra/tmp/audit_global_layer_form.log` (SHA-256
`7f04d0f0099006e2da05c0143a59c53c8c10ae9faf1c9b142c4b6945a7ef8770`).
These temporary files deliberately remain outside the repository history.

## 5. SDP and composition checks

The standard global correlation SDP is

```math
S(A)=\max\{|\langle A,G\rangle|:
                 G\succeq0,\ \operatorname{diag}G=\mathbf1\}.   \tag{23}
```

The correlations `G=(A+tI)^2/(n-1+t^2)` give the exact universal bound
`S(A)>=n sqrt(n-1)`.  This is the previously audited global Gaussian form.
Ordinary Grothendieck followed by polarization yields only

```math
M(A)\ge {S(A)\over4K_G},                           \tag{24}
```

far below `0.336493 n^(3/2)`, even with the new numerical improvement to
`K_G`.  Avoiding (24)'s loss is precisely a symmetric same-switch
Grothendieck/Krivine theorem; labeling the SDP as an indefinite quadratic
form does not make it a separate route.  The already proved diagonal SDP
duals `D_+-A >= 0`, `D_-+A >= 0` remove positive/negative spectral
cancellation but likewise do not select one common Boolean vector.

There is also an exact cross-order description.  For

```math
G=\begin{pmatrix}A&B\\B^{\mathsf T}&D\end{pmatrix},
```

decompose the `k`-subset space by the number `r` of vertices in the first
block.  The diagonal block of `T_k(G)` is

```math
T_r(A)\otimes I+I\otimes T_{k-r}(D),               \tag{25}
```

and the bridge is a block-tridiagonal hard-core lift of `B`, joining
`r` to `r+1` and `r-1`.  This is a genuine composition identity, but the
triangle bound on its bridge is on the leading `N^(3/2)` scale, and a
central-band version approximating `M(G)` retains `2^{N-o(N)}` dimensions.
Controlling (25) with a summable defect therefore still requires the same
energy-dependent anti-alignment as the original bridge problem.

## 6. Research judgment

The global-form analogy yields one clean exact framework:

```text
Boolean energy maximum
  = norm of a full multiplication form
  ~= norm of a central-band non-diagonal compression
  -> block-structured hard-core composition form.
```

It does not yet yield mathematical progress on the frontier.  The verified
obstruction is sharper than "the analogy is vague":

1. the exact first two moments of every permutation-invariant Walsh-level
   compression give only a `Theta(n)` normalized Schatten certificate;
2. every fixed higher normalized Schatten moment has the same `O_q(n)` root
   scale (without claiming a no-go for growing moments or extra structure);
3. self-complementary conference-like and Wigner-like systems have identical
   positive/negative inertia symmetry; and
4. making the compression nearly lossless requires exponentially many
   Walsh states, after which its norm is an equivalent version of the
   original maximum.

A viable successor would need one genuinely new ingredient: a growing-order
moment/tail theorem for these hard-core forms, or a polynomial-complexity
invariant controlling their extreme eigenvalue under block composition.
Without that ingredient, the zeta rank--trace method loses exactly the
`sqrt(n)` scale needed to beat the current lower bound or produce a summable
composition defect.
