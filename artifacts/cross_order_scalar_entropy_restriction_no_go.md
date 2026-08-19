# Scalar entropy--restriction data cannot force a cross-order recurrence

**Status.**  Rigorous proof-class countermodel.  It does not falsify the
quadratic signing recurrence.  It proves that the explicitly listed scalar
heredity, extension, symmetry, Gibbs-variational, and Shearer-entropy
properties cannot by themselves prove even a qualitative fixed-temperature
canonical-disorder recurrence.  Any positive entropy-compensated
restriction theorem must use information absent from those axioms, such as
the special cut-query/restriction incidence at leading scale.

## 1. A log-periodic hereditary cap profile

Put

```math
L_n={n\choose2},\qquad
\omega={\pi\over\log2},qquad
c(x)=0.40+0.01\sin(\omega\log x),
\qquad \tau(x)=x^{3/2}c(x).                             \tag{SC.1}
```

Direct differentiation gives

```math
\tau'(x)=\sqrt x\left\{
 {3\over2}c(x)+0.01\omega\cos(\omega\log x)
 \right\}.                                            \tag{SC.2}
```

The expression in braces lies in `[0.5396,0.6604]`.  Consequently,
after changing finitely many initial terms and rounding to the nearest
integer congruent to `L_n` modulo two, there is a sequence `T_n` satisfying

```math
0<T_{n+1}-T_n<n,qquad
{T_n\over n^{3/2}}=c(n)+O(n^{-3/2}).                   \tag{SC.3}
```

The parity rounding changes `tau(n)` by at most one.  Since its consecutive
increment is between positive constant multiples of `sqrt(n)`, both
inequalities in (SC.3) survive for all sufficiently large `n`.

On the natural edge-signing cube `A_n={+-1}^{L_n}`, define the abstract cap

```math
\widehat Q_n(A)=T_n\qquad(A\in\mathcal A_n).           \tag{SC.4}
```

This model has all of the following properties.

1. It is invariant under vertex permutations, switching, and global sign.
2. Its edge Lipschitz constant is zero.
3. Principal restriction is hereditary, because `T_m<=T_n` for `m<n`.
4. Every labelled order-`m` signing has exactly
   `2^(L_n-L_m)` extensions for a fixed principal restriction.
5. Its one-vertex increment obeys the sharp universal range
   `0<=T_(n+1)-T_n<=n`.
6. Its normalized values remain in `[0.39+o(1),0.41+o(1)]`, strictly inside
   the rigorous asymptotic interval for the true signing problem.

It is even a max-linear Boolean response if no restriction is placed on the
query family:

```math
\widehat Q_n(A)
={T_n\over L_n}\max_{z\in\{\pm1\}^{L_n}}|\langle A,z\rangle|.
                                                               \tag{SC.5}
```

The difference from the real problem is decisive: (SC.5) uses
`2^(Theta(n^2))` arbitrary queries, whereas the quadratic cap uses only the
`2^(n-1)` multiplicatively closed cut queries `z_(ij)=x_ix_j`.

The rounding and phase asymptotics can be reproduced numerically with

```bash
.venv/bin/python computations/verify_scalar_ecr_countermodel.py \
  --min-n 100 --max-n 100000 --min-j 4 --max-j 10 --beta 1
```

## 2. Exact macroscopic failure of ECR

For fixed `beta>0`, the canonical-disorder partition and law are

```math
\widehat Z_n(\beta)
=2^{L_n}e^{-\beta\sqrt nT_n},
\qquad \widehat\mu_{n,\beta}=U_n.                    \tag{SC.6}
```

Every restriction marginal is uniform.  Therefore the Shearer slack and
the marginal canonical mismatch both vanish exactly:

```math
\widehat{\mathcal S}_{N,m}=0,
\qquad
D((\widehat\mu_N)_S\Vert\widehat\mu_m)=0.             \tag{SC.7}
```

Choose `N_j` to be the nearest integer to `4^j/sqrt(2)` and put
`m_j=floor(N_j/2)`.  Then

```math
\omega\log N_j=-{\pi\over2}+o(1)\pmod {2\pi}.          \tag{SC.8}
```

Since `omega log(1/2)=-pi`,

```math
{T_{N_j}\over N_j^{3/2}}=0.39+o(1),
\qquad
{T_{m_j}\over m_j^{3/2}}=0.41+o(1).                   \tag{SC.9}
```

With `q=L_m/L_N`, the restriction-energy excess from the canonical identity
is

```math
\begin{aligned}
\widehat{\mathcal D}_{N_j,m_j}
&={\sqrt {m_j}\over q}T_{m_j}-\sqrt {N_j}T_{N_j}\\
&=N_j(N_j-1){m_j\over m_j-1}
  {T_{m_j}\over m_j^{3/2}}
  -N_j^2{T_{N_j}\over N_j^{3/2}}\\
&=(0.02+o(1))N_j^2.                                    \tag{SC.10}
\end{aligned}
```

Combining (SC.7)--(SC.10) gives the quantitative falsifier

```math
\boxed{
\beta\widehat{\mathcal D}_{N_j,m_j}
-\widehat{\mathcal S}_{N_j,m_j}
=(0.02\beta+o(1))N_j^2.}                               \tag{SC.11}
```

Thus neither `O(N^(2-delta))` nor merely `o(N^2)` ECR follows from the
listed scalar resources.  Equivalently, for

```math
\widehat\psi_n=L_n^{-1}\log\widehat Z_n,
```

one has

```math
\widehat\psi_{N_j}(\beta)-\widehat\psi_{m_j}(\beta)
=0.04\beta+o(1).                                       \tag{SC.12}
```

## 3. The degeneracy is not essential

For a nonconstant version, let

```math
h_n(A)=2\,\mathbf1\left\{
 \left|\sum_{i<j<k}a_{ij}a_{jk}a_{ki}\right|
 \ge {1\over2}{n\choose3}\right\},
\qquad
\widehat Q_n^{\rm nc}(A)=T_n+h_n(A).                  \tag{SC.13}
```

The triangle statistic is invariant under switching, permutations, and
global edge-sign reversal.  The cap remains parity-correct, one edge flip
changes it by at most two, and (SC.2) makes restriction hereditary and the
one-vertex increment at most `n` for all sufficiently large orders.
Moreover, if `R_n` denotes the triangle sum in (SC.13), orthogonality of
distinct triangle monomials under uniform disorder gives

```math
\mathbb E_U R_n^2={n\choose3}.
```

Thus some signing has
`|R_n|<=sqrt(binomial(n,3))<binomial(n,3)/2`, whereas the all-positive
signing activates the layer.  The model is genuinely nonconstant for all
large `n`.

The extra canonical tilt has log-density range at most
`2 beta sqrt(n)`.  Hence

```math
D(\widehat\mu_n^{\rm nc}\Vert U_n)=O_\beta(\sqrt n),
\qquad
0\le\widehat{\mathcal S}_{N,m}^{\rm nc}
=O_\beta(\sqrt N),                                     \tag{SC.14}
```

while the contribution of `h` to `mathcal D_(N,m)` is only `O(sqrt N)` on
comparable splits.  The quadratic term in (SC.11) and the fixed pressure
gap in (SC.12) therefore survive unchanged.

### Proposition SC.1 (stability under subleading scalar structure)

More generally, let `R_n` be arbitrary functions on the signing cubes with

```math
r_n:=\sup_A|R_n(A)|=o(n^{3/2}),                        \tag{SC.15}
```

and set `Q_n^R=T_n+R_n`.  Whenever the resulting cap satisfies the desired
structural axioms, its canonical law has density relative to uniform whose
logarithmic range is at most

```math
2\beta\sqrt n\,r_n.
```

Consequently, on every comparable restriction,

```math
\begin{aligned}
D(\mu_{n,\beta}^R\Vert U_n)&\le2\beta\sqrt n\,r_n,\\
0\le \mathcal S_{N,m}^R&\le2\beta\sqrt N\,r_N,\\
\mathcal D_{N,m}^R
&=\mathcal D_{N,m}^{\rm base}
  +o(N^2).                                              \tag{SC.16}
\end{aligned}
```

Quantitatively, if `r_n<=Cn^gamma` for some `gamma<3/2`, all three error
terms in (SC.16) are `O_(beta,C)(N^(gamma+1/2))` on comparable splits.

In particular, along the pairs in (SC.8),

```math
\boxed{
\beta\mathcal D_{N_j,m_j}^R-\mathcal S_{N_j,m_j}^R
=(0.02\beta+o(1))N_j^2.}                              \tag{SC.17}
```

The associated normalized pressure also satisfies

```math
\psi_{N_j}^R(\beta)-\psi_{m_j}^R(\beta)
=0.04\beta+o(1).                                      \tag{SC.17a}
```

**Proof.**  A probability density obtained by tilting a reference law by a
function with range `Delta` has relative entropy at most `Delta`.  Apply
this with `Delta<=2 beta sqrt(n) r_n`.  The exact information-loss
identity

```math
\mathcal S_{N,m}^R
=D(\mu_{N,\beta}^R\Vert U_N)
-q^{-1}\mathbb E_SD((\mu_{N,\beta}^R)_S\Vert U_m)
```

then bounds the nonnegative slack by the first relative entropy.  Finally,
on a comparable split `sqrt(m)/q=O(sqrt(N))`, so the two `R` terms in
`mathcal D` are `o(N^2)`.  The correction to normalized log partition is
at most `beta sqrt(n)r_n/L_n=o(1)`.  Equations (SC.17) and (SC.17a)
follow.  `square`

The exponent threshold is sharp for arbitrary perturbations: a
parity-rounded monotone profile
`bar(T)_n=0.40n^(3/2)+O(1)` and the deterministic choice
`R_n=bar(T)_n-T_n=O(n^(3/2))` erase the oscillation completely.

The stability statement can retain genuine information about the real cut
cap, rather than only an artificial local statistic.  Let `Q_n^{cut}` be
the actual quadratic Boolean cap and put

```math
g_n(A)=2\left\lfloor\sqrt{Q_n^{cut}(A)}\right\rfloor,
\qquad
\widehat Q_n^{cut,sub}(A)=T_n+g_n(A).                  \tag{SC.18}
```

This perturbation is parity-correct and invariant under switching,
permutations, and global sign.  The elementary restriction inequality
`Q_m^{cut}(A[S])<=Q_N^{cut}(A)` makes it hereditary.  An edge flip changes
`Q^{cut}` by at most two and hence changes `g` by at most two.  Also

```math
Q_{n+1}^{cut}(A)\le Q_n^{cut}(A[n])+n,
\qquad Q_n^{cut}(A)\ge\sqrt{L_n},                      \tag{SC.19}
```

so `g_(n+1)-g_n=O(sqrt(n))`; together with (SC.2), the one-vertex increment
of (SC.18) lies in `[0,n]` for all large `n`.  Since `g_n=O(n)`, Proposition
SC.1 applies with `gamma=1`.  Thus even a subleading, parity-correct scalar
record of the genuine cut cap does not repair ECR.  This does **not** test
the exact cut-query algebra at leading scale: that is precisely the
structure which remains outside the countermodel.

### Proposition SC.2 (few cut queries do not suffice with a vanishing coefficient)

There is also an exact response representation using only the genuine
augmented cut group.  Define the real-valued cap

```math
F_n(A)=T_n+{1\over n}Q_n^{cut}(A)
=\max_{x\in\{\pm1\}^n,\,\sigma\in\{\pm1\}}
 \left\{T_n+{\sigma\over n}
   \sum_{i<j}a_{ij}x_ix_j\right\}.                    \tag{SC.20}
```

Thus (SC.20) has only `2^n` distinct affine queries, whose directions are
exactly the multiplicatively closed augmented cut characters.  It is
symmetric and hereditary for all large orders.  Indeed, if `B=A[S]`, then

```math
{Q_m^{cut}(B)\over m}-{Q_N^{cut}(A)\over N}
\le {N-m\over2},                                      \tag{SC.21}
```

whereas `T_N-T_m` is at least a positive constant times
`sqrt(m)(N-m)` away from finitely many initial orders.  For a one-vertex
extension,

```math
{Q_{n+1}^{cut}\over n+1}-{Q_n^{cut}\over n}\le1,
```

so the extension increment is again in `[0,n]` eventually.  The edge
Lipschitz constant of the added term is `2/n`.

The deterministic `T_n` drops out of the canonical law, while the remaining
tilt is `beta Q_n^{cut}/sqrt(n)`.  Its log-density range is `O_beta(n^(3/2))`.
Proposition SC.1 therefore applies with `gamma=1` and proves the same
quadratic ECR failure.  This shows that merely having `exp(O(n))` genuine
cut directions and their multiplicative closure is still insufficient if
they appear with a vanishing coefficient beside an order-dependent affine
offset.

That qualification cannot be removed.  If `w>0`,
`mathcal C subseteq {+-1}^L`, and a homogeneous response

```math
G(a)=w\max_{z\in\mathcal C}|\langle a,z\rangle|
```

is identically `wL` on the `L`-dimensional sign cube, then
`|<a,z>|=L` holds exactly when `z=+-a`.  Thus, for every `a`, the query set
must contain one of `a,-a`.
Hence `|C|>=2^(L-1)`.  The constant construction (SC.5) cannot be made into
an `exp(O(n))` homogeneous flat response.  What remains untested, and is
therefore the only identified escape for a positive ECR theorem, is some
leading-scale property absent from these models.  The simultaneous package
of zero affine offset, unit leading edge coefficient, and joint restriction
incidence of the genuine cut characters is the canonical such candidate.

## 4. Consequence for the outward campaign

This theorem is a proof-class no-go, not a theorem about actual minimizers.
It rules out deriving entropy-compensated restriction from the explicitly
listed scalar properties, even after adding arbitrary `o(n^(3/2))`
structure and, in particular, a subleading function of the genuine cut
cap.  Proposition SC.2 additionally rules out treating the mere number and
multiplicative closure of cut directions as sufficient while their
coefficient is subleading.  A positive theorem must use a leading property
absent from the countermodels, such as unit-normalized cut-query incidence
or a genuinely quadratic restriction-incidence inequality.  Merely
invoking unrestricted max-linearity, canonical optimality, cap parity,
extension counts, or entropy contraction cannot improve the cross-order
exponent.
