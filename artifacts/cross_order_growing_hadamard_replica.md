# Growing Hadamard replicas: exact core and annealed sign completion

Status: **proved exact all-order replica identity and a direct cross-order
inequality whose sign-completion defect is only `O_beta(k)`**.  In
particular, the `r` missing matching edges in the universal double cost
`O_beta(1)`, not the previously used pointwise `O_beta(sqrt(r))`.  The
remaining core term is an exact rotated-cube pressure; no sublinear bound
for that term is proved here.

Let `A` be a hollow symmetric signing of order `r`, and write

```math
H_A(x)={1\over2}x^{\mathsf T}Ax,
\qquad
\phi_A(s)=\log\mathbb E_x\cosh(sH_A(x)).             \tag{1}
```

Let `T` be a symmetric Hadamard matrix of order `k`:

```math
T\in\{\pm1\}^{k\times k},
\qquad T=T^{\mathsf T},
\qquad T^2=kI_k.                                     \tag{2}
```

The statements below apply only at orders for which such a `T` is
available.  In particular, the Sylvester matrices
`T=H_2^{\otimes d}` give the infinite valid family `k=2^d`.  No existence
claim is made for arbitrary `k`.

## 1. Exact diagonalized core identity

Index the `kr` parent vertices by `(a,i) in [k] times [r]`.  The Kronecker
core

```math
K_0=T\otimes A                                      \tag{3}
```

has sign entries except on the pairs `((a,i),(b,i))`, `a<b`, where its
entry is zero.  There are exactly

```math
L={k\choose2}r                                      \tag{4}
```

such unordered matching-coordinate pairs.

Represent a Boolean parent configuration by a `k` by `r` matrix `X`, whose
rows are the replica spins.  If

```math
T=U\,\sqrt{k}\operatorname{diag}(\sigma_1,\ldots,\sigma_k)U^{\mathsf T},
\qquad \sigma_p\in\{\pm1\},                         \tag{5}
```

and `Y=U^T X` has rows `y_p`, then

```math
\boxed{
H_{K_0}(X)
=\sqrt{k}\sum_{p=1}^k\sigma_p H_A(y_p).}            \tag{6}
```

This is just

```math
{1\over2}\operatorname{vec}(X)^{\mathsf T}
 (T\otimes A)\operatorname{vec}(X)
={\sqrt{k}\over2}\sum_p\sigma_p y_p^{\mathsf T}Ay_p.
```

Put

```math
s={\beta\over\sqrt r},
\qquad t={\beta\over\sqrt{kr}}.                    \tag{7}
```

For each coordinate `i`, the column `Y_{\bullet i}` has the rotated-cube
law

```math
\nu_U=(U^{\mathsf T})_\#\operatorname{Unif}\{\pm1\}^k,
```

and these `r` columns are independent.  Hence (6) gives the exact pressure
identity

```math
\boxed{
\phi_{K_0}(t)
=\log\mathbb E_{Y\sim\nu_U^{\otimes r}}
 \cosh\!\left(s\sum_{p=1}^k\sigma_pH_A(y_p)\right).} \tag{8}
```

Thus diagonalization does not produce `k` independent Boolean children:
it replaces their one-site product law by the joint rotated-cube law
`nu_U`.  Formula (8), rather than a scalar child pressure, is the exact
core obligation.

Two normalizations are preserved exactly:

```math
t\|T\otimes A\|_{\rm op}=s\|A\|_{\rm op},           \tag{9}
```

and, for every integer `j>=1`,

```math
t^{2j}\operatorname{Tr}(T\otimes A)^{2j}
=k s^{2j}\operatorname{Tr}A^{2j}.                  \tag{10}
```

Indeed, `T^2=kI`.  Consequently increasing the replica number does not by
itself enter a stricter operator-temperature regime; it preserves the
child's normalized spectral scale.

## 2. Annealed completion of all missing signs

For each of the `L` pairs in (4), choose an independent sign `b_e`.  Let
`K_b` be the resulting hollow symmetric signing, and let
`chi_e(X)` be the Boolean edge monomial at that pair.  For every fixed
configuration,

```math
\begin{aligned}
\mathbb E_b\cosh\!\left(t\left[H_{K_0}(X)
 +\sum_eb_e\chi_e(X)\right]\right)
&=\cosh(tH_{K_0}(X))\prod_e\cosh(t\chi_e(X))\\
&=\cosh(tH_{K_0}(X))(\cosh t)^L.
\end{aligned}                                       \tag{11}
```

Averaging also over `X` and then taking the best filling proves

```math
\boxed{
\min_b\phi_{K_b}(t)
\le\phi_{K_0}(t)+L\log\cosh t
\le\phi_{K_0}(t)+{\beta^2(k-1)\over4}.}             \tag{12}
```

The last inequality uses `log cosh t<=t^2/2` and the exact value of `L`.
This annealed completion is deterministic by the probabilistic method and
does not require `k=o(r)`.

## 3. Direct pressure-defect implication

Let

```math
P_r(\beta)=\min_A\phi_A(\beta/\sqrt r),             \tag{13}
```

and let `A` be an actual order-`r` minimizer.  Define the exact rotated-core
defect

```math
\mathcal D_{k,r}(T,A;\beta)
=\phi_{T\otimes A}(\beta/\sqrt{kr})-kP_r(\beta).    \tag{14}
```

Since every `K_b` in (12) is a valid order-`kr` signing, (12) gives the
required immediate arrow to the permanent cross-order defect:

```math
\boxed{
P_{kr}(\beta)-kP_r(\beta)
\le\mathcal D_{k,r}(T,A;\beta)
 +{\beta^2(k-1)\over4}.}                            \tag{15}
```

In particular, for any valid Hadamard orders `k=k(r)`,

```math
\boxed{
\mathcal D_{k,r}(T,A;\beta)=O_\beta((kr)^{1-\delta})
\Longrightarrow
P_{kr}(\beta)-kP_r(\beta)
=O_\beta((kr)^{1-\delta}+k).}                       \tag{16}
```

Because `k/(kr)=1/r`, the completion term is `o(kr)` along every sequence
with `r->infinity`, even if `k` grows rapidly.  All substantive content is
therefore in the rotated-cube term (14), not in integrality of the parent.

### Universal-double correction

For `k=2`, take the Sylvester matrix

```math
T=\begin{pmatrix}1&1\\1&-1\end{pmatrix}.
```

The archived induced-clique orbit identity identifies
`phi_(T tensor A)(beta/sqrt(2r))` exactly.  Formula (15) sharpens its
matching-edge payment to

```math
\boxed{
P_{2r}(\beta)-2P_r(\beta)
\le\mathcal D_{2,r}(T,A;\beta)
 +r\log\cosh\!\left({\beta\over\sqrt{2r}}\right)
\le\mathcal D_{2,r}(T,A;\beta)+{\beta^2\over4}.}    \tag{17}
```

Thus an orbit estimate `mathcal D_(2,r)=O(r^(1-delta))` gives the same
exponent, without the former `O_beta(sqrt r)` floor.

## 4. Zero-temperature completion

The same core also has a uniformly harmless sign completion for the cap.
For random `b`, Hoeffding's inequality and a union bound over all `2^(kr)`
configurations show that some filling satisfies

```math
\max_X\left|\sum_eb_e\chi_e(X)\right|
\le\sqrt{2L(kr+2)\log2}.
```

Therefore

```math
\boxed{
M(K_b)\le M(T\otimes A)
 +\sqrt{k(k-1)r(kr+2)\log2}.}                       \tag{18}
```

After division by `(kr)^(3/2)`, the added term is
`O(r^(-1/2))`, uniformly in `k`.  Again, the growing replica number is not
the completion obstruction.

## 5. Verdict

### A genuine growing-`k` core obstruction

Increasing `k` does not, by itself, make the exact core defect in (14)
small.  The following amplification makes that statement quantitative.

Let `S` be a symmetric Hadamard matrix of order `h`, and suppose a signing
`A` of order `r` satisfies

```math
M(S\otimes A)=(1+\eta)h^{3/2}M(A)                  \tag{19}
```

for some `eta>0`.  For `q=4^d`, put

```math
R_q=(J_4-2I_4)^{\otimes d}.
```

Then `R_q` is a symmetric regular Hadamard matrix, and
`R_q 1=sqrt(q) 1`.  Therefore `T=R_q tensor S` is a symmetric Hadamard
matrix of order `k=qh`.  If `z` maximizes the core in (19), the Boolean
configuration `1 tensor z` proves

```math
M(T\otimes A)
\ge q^{3/2}M(S\otimes A)
=(1+\eta)k^{3/2}M(A).                              \tag{20}
```

The two extremal configurations `w,-w` and
`cosh u>=e^|u|/2` give

```math
\phi_{T\otimes A}(\beta/\sqrt{kr})
\ge {\beta\over\sqrt{kr}}M(T\otimes A)-kr\log2.
```

On the other hand `phi_A(beta/sqrt(r))<=beta M(A)/sqrt(r)`.
Consequently, if `A` is an actual order-`r` minimizer, then its exact core
defect obeys

```math
\boxed{
\mathcal D_{k,r}(T,A;\beta)
\ge k\left({\beta\eta M(A)\over\sqrt r}-r\log2\right).} \tag{21}
```

Thus one fixed positive tensor-cap excess creates a pressure-core floor
linear in `kr` for every

```math
\beta>{r^{3/2}\log2\over\eta M(A)}.                 \tag{22}
```

This is not merely hypothetical.  Take the all-temperature exact order-four
pressure minimizer

```math
A=\begin{pmatrix}
0&1&-1&-1\\
1&0&-1&-1\\
-1&-1&0&-1\\
-1&-1&-1&0
\end{pmatrix},
```

whose thermal minimality is proved in the universal-double audit.  It has
`M(A)=4`, while for

```math
S=\begin{pmatrix}1&1\\1&-1\end{pmatrix}
```

the Boolean vector
`(-1,-1,-1,-1,1,1,-1,-1)` gives `M(S tensor A)=12`.
Hence `eta=3/(2sqrt(2))-1`, and for `k=2*4^d`, (21) becomes

```math
\boxed{
\mathcal D_{k,4}(T,A;\beta)
\ge k\left[\left({3\over\sqrt2}-2\right)\beta
 -4\log2\right].}                                  \tag{23}
```

The coefficient is positive for

```math
\beta>{4\log2\over3/\sqrt2-2}
=22.8534526105817\ldots.                            \tag{24}
```

Equations (23)--(24) rigorously rule out the claim that the rotated-core
term becomes sublinear merely because the Hadamard replica count tends to
infinity, even when the seed is an actual pressure minimizer.  Their scope
is deliberately precise: the child order here is fixed.  They do not rule
out a theorem for a selected sequence of actual minimizers with
`r->infinity`, nor do they lower-bound the best completed parent's pressure;
they obstruct the specific core certificate (15).

### Exact finite wind tunnel: the core defect is not monotone in `k`

Exact Gray-code enumeration gives a complementary falsifier to any simple
monotonicity claim.  The following entries are
`mathcal D_(k,r)/(kr)`; every child shown is an exact pressure minimizer.

| child | `k` | `beta=1` | `beta=2` | `beta=4` | `beta=8` |
|---|---:|---:|---:|---:|---:|
| order-3 triangle | 2 | -0.007264 | -0.067598 | -0.338923 | -1.006338 |
| order-3 triangle | 4 | -0.007152 | -0.005947 | +0.048199 | +0.057665 |
| order-3 triangle | 8 | -0.008550 | -0.055698 | -0.187731 | -0.449769 |
| order-4 minimizer above | 2 | -0.003715 | -0.010917 | +0.031532 | +0.155851 |
| order-4 minimizer above | 4 | -0.005224 | -0.020207 | -0.037513 | -0.043196 |

These signs are exact consequences of the recorded integer energy
histograms followed by high-precision evaluation, not heuristic sampling.
They show both phenomena:

1. the order-four positive universal-double core defect is removed, rather
   than amplified, by the Sylvester order-four outer matrix;
2. the order-three defect is negative at `k=2`, positive at `k=4` for
   `beta=4,8`, and negative again at `k=8`.

The low-temperature limits in the equal-cap cases are also exact from the
top energy multiplicities:

```math
\lim_{\beta\to\infty}\mathcal D_{4,3}=+\log2,
\qquad
\lim_{\beta\to\infty}\mathcal D_{4,4}=-\log2.       \tag{25}
```

Thus neither `mathcal D_(k,r)` nor its normalized version has a universal
one-step monotonicity in the replica order, even on exact minimizers.
Growing `k` can repair one seed and spoil another.  Any asymptotic theorem
must use a structural property of the selected large-order minimizers and
cannot follow from replica count alone.

The exact histograms and all displayed values are generated by
[`audit_growing_hadamard_core.cpp`](../computations/audit_growing_hadamard_core.cpp)
and stored in
[`growing_hadamard_core_exact.json`](../computations/results/growing_hadamard_core_exact.json).
Reproduce them without using a system temporary directory via

```text
g++ -O3 -std=c++17 computations/audit_growing_hadamard_core.cpp \
  -o /home/math/quadra/tmp/audit_growing_hadamard_core
/home/math/quadra/tmp/audit_growing_hadamard_core \
  computations/results/growing_hadamard_core_exact.json
```

## 6. Verdict

Growing Hadamard replication removes the sign-completion loss completely
at leading order, and (15) is a genuine all-order quantitative composition
inequality for every available replica order.  It does **not** yet improve
the exponent for actual optimizing children: the exact core defect (14)
is uncontrolled, and (9)--(10) show that increasing `k` supplies no small
operator-temperature or spectral-moment parameter.  Any successful use of
(15) must prove a non-spectral rotated-cube universality estimate for the
actual minimizing child, while a scalable obstruction must exhibit a
linear positive value of (14) along actual large-order minimizers.
