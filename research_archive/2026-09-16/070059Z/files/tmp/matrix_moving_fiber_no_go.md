# Matrix-valued moving fibers: exact conditions and a scoped no-go

Date: 2026-08-15.

Status: theorem-first derivation.  The operator coset inequality and the
constant-root compression theorem below are exact.  The final support/rank
obstruction uses the already verified asymptotic cube-subgraph theorem of
Bollobas--Lee--Letzter.  No bound on `M_n` is improved.

The purpose is to test whether a genuinely matrix-valued, multi-representative
fiber can evade the scalar partial-transversal obstruction without paying its
channels separately.  It cannot do so in the constant-root cone below: the
complete Loewner inequality has a scalar compression with exactly the same
strength.  At the required `n^(-1/2)` scale that compression needs
`exp(Omega(n log n))` Fourier support.  A natural termwise-orthogonal
multi-transversal implementation consequently needs hidden dimension of the
same order.

## 1. Operator-valued joint-Gram cone

Write the edge cube additively as

```math
G=\mathbb F_2^E,
\qquad E=\binom n2,
```

let `C=C_n^+` be the augmented cut code, and put

```math
D=C^\perp,
\qquad
\tau(z)={1\over E}\sum_{e=1}^E(-1)^{z_e}.
```

Thus `|C|=2^n`, `|G/D|=2^n`, and for an edge signing `a`,

```math
\mu(a)=\max_{c\in C}\tau(a+c)={Q(a)\over E}.
```

Let `H` be a finite-dimensional real or complex Hilbert space and let

```math
K:G\longrightarrow \operatorname{Herm}(H).
```

Use normalized Fourier coefficients

```math
Q_R=\widehat K(R)=2^{-E}\sum_{z\in G}\chi_R(z)K(z).
```

Assume the following three *joint* conditions:

```math
K(z)\succeq0 \quad(z\in G),                                      \tag{1.1}
```

```math
Q_R\succeq0 \quad(R\in G),                                      \tag{1.2}
```

and, for some `lambda>=0`,

```math
B_R:=(A_EQ)_R-\lambda Q_R\succeq0 \quad(R\in G),                 \tag{1.3}
```

where normalized cube adjacency acts on the Fourier index:

```math
(A_EQ)_R={1\over E}\sum_{e=1}^E Q_{R+e}.                         \tag{1.4}
```

Condition (1.2) says that `K` is operator-valued positive definite.
Condition (1.3) says that the complete matrix-valued remainder

```math
F(z)=(\tau(z)-\lambda)K(z)                                      \tag{1.5}
```

is operator-valued positive definite.  No scalar matrix entry, eigenchannel,
or left/right response is bounded separately in these assumptions.

## 2. Exact Loewner coset inequality

Define

```math
T_a=\sum_{c\in C}K(a+c),
\qquad
S_a=\sum_{c\in C}F(a+c),
\qquad
J=\sum_{c\in C}F(c).                                            \tag{2.1}
```

> **Theorem 2.1 (operator moving-fiber inequality).**  Under
> (1.1)--(1.3), for every `a in G`,
>
> ```math
> \boxed{(\lambda-\mu(a))_+T_a\preceq J.}                        \tag{2.2}
> ```
>
> On the support of `T_a`, this implies
>
> ```math
> \boxed{
> \mu(a)\ge
> \lambda-
> \lambda_{\min}
> \left(T_a^{-1/2}JT_a^{-1/2}\right).}                           \tag{2.3}
> ```

**Proof.**  Subgroup Fourier orthogonality gives

```math
S_a=|C|\sum_{d\in D}\chi_d(a)B_d,
\qquad
J=|C|\sum_{d\in D}B_d.                                         \tag{2.4}
```

Every `B_d` is positive semidefinite, so the *complete* sums satisfy

```math
-J\preceq S_a\preceq J.                                         \tag{2.5}
```

Pointwise positivity of `K` and `tau(a+c)<=mu(a)` also give

```math
S_a
=\sum_{c\in C}(\tau(a+c)-\lambda)K(a+c)
\preceq(\mu(a)-\lambda)T_a.                                    \tag{2.6}
```

If `mu(a)<lambda`, combine the lower inequality in (2.5) with (2.6).
If `mu(a)>=lambda`, (2.2) is automatic.  Taking the least generalized
Rayleigh quotient of the pair `(J,T_a)` proves (2.3). `square`

This proof preserves cancellation in the entire fiber until after (2.2).
In particular, it is not a scalar-channel proof written in matrix notation.

## 3. Exact matrix-square realization conditions

A natural multi-representative realization starts with maps

```math
V_R:H\longrightarrow W
```

of finite Fourier support `U={R:V_R!=0}` and forms

```math
G_V(z)=\sum_R\chi_R(z)V_R,
\qquad
K(z)=G_V(z)^*G_V(z).                                             \tag{3.1}
```

Pointwise positivity (1.1) is then automatic.  The exact Fourier
coefficients are

```math
\boxed{Q_T=\sum_R V_R^*V_{R+T}.}                                 \tag{3.2}
```

They are Hermitian because the cube has exponent two, but they need not be
positive semidefinite.  If

```math
(A_EV)_R={1\over E}\sum_eV_{R+e},
```

then the exact coefficients of the complete remainder are

```math
\boxed{
B_T
=\sum_RV_R^*\bigl((A_EV)_{R+T}-\lambda V_{R+T}\bigr)
=(A_EQ)_T-\lambda Q_T.}                                         \tag{3.3}
```

Consequently, within the square ansatz (3.1), the two exact Gram
requirements are

```math
\sum_RV_R^*V_{R+T}\succeq0,                                    \tag{3.4}
```

```math
\sum_RV_R^*\bigl((A_EV)_{R+T}-\lambda V_{R+T}\bigr)\succeq0     \tag{3.5}
```

for every `T`.  Merely requiring the individual `V_R` to be positive, or
diagonalizing all of them, does not establish (3.4)--(3.5).  Simultaneous
diagonalization would instead return to separately positive scalar channels,
the case that this attempt is supposed to avoid.

## 4. Constant-root condition and exact quotient leakage

The clean optimizer-independent target is

```math
Q_d=0\qquad(d\in D\setminus\{0\}).                              \tag{4.1}
```

It makes the root mass constant:

```math
T_a=|C|Q_0=|C|P,
\qquad
P:=Q_0=\sum_RV_R^*V_R,                                          \tag{4.2}
```

for every `a`.

There is an exact quotient formula for the remaining leakage.  For each
coset `x in G/D`, put

```math
Y_x=\sum_{R\in x}V_R,                                           \tag{4.3}
```

and let quotient adjacency be

```math
(\overline A Y)_x={1\over E}\sum_eY_{x+\overline e}.            \tag{4.4}
```

Then

```math
\sum_xY_x^*Y_x=\sum_{d\in D}Q_d=P                              \tag{4.5}
```

and

```math
\mathcal R
:=\sum_xY_x^*(\overline A Y)_x
=\sum_{d\in D}(A_EQ)_d.                                       \tag{4.6}
```

Thus

```math
\boxed{{J\over|C|}=\mathcal R-\lambda P\succeq0.}               \tag{4.7}
```

Equations (4.2) and (4.7) turn Theorem 2.1 into

```math
(2\lambda-\mu(a))P\preceq\mathcal R.                           \tag{4.8}
```

Restricting to `ran P`, define

```math
\rho_{\min}
=\lambda_{\min}(P^{-1/2}\mathcal RP^{-1/2})\ge\lambda.         \tag{4.9}
```

The resulting uniform lower bound is

```math
\boxed{\mu(a)\ge2\lambda-\rho_{\min}.}                         \tag{4.10}
```

Equivalently, if the least quotient leakage is
`ell_min=rho_min-lambda`, then (4.10) is

```math
\mu(a)\ge\lambda-\ell_{\min}.                                 \tag{4.11}
```

This is the exact matrix-valued analogue of the scalar quotient leakage,
with all representatives interacting before any generalized eigenvalue is
taken.

## 5. Constant-root matrix fibers compress without loss

The hoped-for genuinely joint advantage fails in this cone.

> **Theorem 5.1 (lossless scalar compression).**  Every nonzero
> matrix-valued constant-root certificate satisfying (1.1)--(1.3) and
> (4.1) contains a scalar kernel satisfying all the same positivity and
> constant-root conditions and attaining exactly the bound (4.10).

**Proof.**  Choose a generalized eigenvector `u in ran P` such that

```math
u^*\mathcal Ru=\rho_{\min}u^*Pu.                                \tag{5.1}
```

Compress only after the complete operator theorem has been formed:

```math
\kappa_u(z)=u^*K(z)u,
\qquad
q_R^{(u)}=u^*Q_Ru,
\qquad
b_R^{(u)}=u^*B_Ru.                                               \tag{5.2}
```

Then `kappa_u(z)>=0`, `q_R^(u)>=0`, and `b_R^(u)>=0`.  Moreover,

```math
b^{(u)}=A_Eq^{(u)}-\lambda q^{(u)},
\qquad
q_d^{(u)}=0\quad(d\in D\setminus\{0\}).                       \tag{5.3}
```

Its root mass is `u^*Pu`, while its quotient Rayleigh quotient is exactly
`rho_min` by (5.1).  The scalar version of (4.10) is therefore the same
inequality, not a weakened trace average. `square`

The proof of Theorem 2.1 never paid scalar channels.  Theorem 5.1 is instead
an obstruction: once the root mass is made constant and the numerator is a
positive Loewner operator, a worst generalized direction necessarily
witnesses the full bound.  Matrix rank by itself cannot improve it.

## 6. Correct-scale support is necessarily superexponential

Let `u` be any nonzero scalar compression, not necessarily the minimizing
one, and put

```math
q_R=u^*Q_Ru,
\qquad
S=\operatorname{supp}q.
```

From (1.2)--(1.3),

```math
q\ge0,
\qquad
A_Eq\ge\lambda q.                                               \tag{6.1}
```

After restricting to `S`, Collatz--Wielandt gives

```math
\rho\bigl(Q_E[S]\bigr)\ge\lambda E,                            \tag{6.2}
```

where the left side uses unnormalized cube adjacency.  Theorem 4 of
Bollobas--Lee--Letzter, *Eigenvalues of subgraphs of the cube*, together
with their Hamming-ball spectral asymptotics, then implies the following
already verified inversion:

```math
\lambda\ge {c\over\sqrt n}
\quad\Longrightarrow\quad
\boxed{|S|\ge\exp(\Omega_c(n\log n))}.                           \tag{6.3}
```

Indeed, a Hamming ball capable of unnormalized spectral radius
`cE/sqrt(n)` must have radius `Omega_c(n)` in dimension
`E=Theta(n^2)`, and such a ball has
`exp(Omega_c(n log n))` vertices.  This is precisely the inverse form of
the scalar support estimate already used for partial transversals.

For a square realization (3.1),

```math
\operatorname{supp}q^{(u)}\subseteq U+U,                         \tag{6.4}
```

because (3.2) can be nonzero only when two elements of `U` differ by the
given Fourier index.  Hence

```math
\boxed{|U|\ge\exp(\Omega_c(n\log n))}                            \tag{6.5}
```

whenever `lambda>=c/sqrt(n)`.

Now consider the natural *multi-representative transversal* implementation
of (4.1), in which root orthogonality is imposed term by term:

```math
V_R^*V_S=0
\quad\text{for distinct }R,S\text{ in the same }D\text{-coset}. \tag{6.6}
```

The nonzero ranges of the maps in any one `D`-coset are then mutually
orthogonal subspaces of `W`.  If `r=dim W`, that coset contains at most `r`
nonzero representatives.  Since `|G/D|=2^n`,

```math
|U|\le r2^n.                                                     \tag{6.7}
```

Combining (6.5)--(6.7) proves the scalable rank obstruction

```math
\boxed{r\ge\exp(\Omega_c(n\log n)).}                            \tag{6.8}
```

Thus polynomial rank, `exp(O(n))` rank, and more generally
`exp(o(n log n))` rank cannot realize a correct-scale orthogonal
multi-transversal moving square.

## 7. The partial-matching orbit: exact reduction

The proposed high-rank orbit is the first natural support large enough to
meet (6.3).  Let `M_l` denote the `l`-edge partial matchings of `K_n`.  Then

```math
|\mathcal M_l|
={n!\over(n-2l)!2^l l!}
=\exp(\Theta(n\log n))                                          \tag{7.1}
```

when `l=alpha n` with `0<alpha<1/2`.  Normalized cube adjacency restricted
to the matching layers has radial off-diagonal coefficient

```math
c_l
={\sqrt{(l+1)\binom{n-2l}{2}}\over E}.                          \tag{7.2}
```

Thus `c_l=Theta(n^(-1/2))` on linear layers.  This orbit has exactly the
support growth demanded by (6.3); unlike bounded-rank or `exp(O(n))`
states, it is not excluded for size reasons.

### 7.1 Exact cut twirl and signed hafnian fibers

For a matching `M`, let

```math
\partial M=\{i:\deg_M(i)=1\}
```

be its covered vertex set.  Two partial matchings satisfy

```math
M+N\in D
\quad\Longleftrightarrow\quad
\partial M=\partial N.                                         \tag{7.3}
```

Indeed, equality of the boundaries is equivalent to even degree in the
symmetric difference; that difference is a disjoint union of alternating
even cycles and automatically has even size.  (The converse follows from
even degree.)

Suppose the square amplitude (3.1) is supported on partial matchings and
define the matrix-valued signed matching sums

```math
H_U(a)
=\sum_{M:\partial M=U}\chi_M(a)V_M,
\qquad
Y_U=H_U(0)=\sum_{M:\partial M=U}V_M.                             \tag{7.4}
```

Character orthogonality and (7.3) give the exact root mass

```math
\boxed{
{T_a\over|C|}
=\sum_{U\subseteq[n],\ |U|\ {\mathrm{even}}}H_U(a)^*H_U(a).}     \tag{7.5}
```

For scalar layer-constant amplitudes this is

```math
{T_a\over|C|}
=\sum_l v_l^2\sum_{|U|=2l}
\left|\sum_{M\text{ perfect on }U}\chi_M(a)\right|^2,          \tag{7.6}
```

the sum of squares of signed principal hafnians.  This formula exhibits the
remaining nonconstant-root obligation exactly; it does not query a maximizing
spin.

The quotient `G/D` is labelled by vertex boundary together with edge-count
parity.  Only one parity is occupied by matchings for a fixed boundary.  If
`Abar` toggles the two boundary vertices and the parity of one edge, then

```math
\boxed{
{J\over|C|}
=\sum_UY_U^*(\overline A Y)_U
-\lambda\sum_UY_U^*Y_U.}                                      \tag{7.7}
```

Combining (7.5) and (7.7) with Theorem 2.1 is an exact, unsplit
matching-fiber Loewner certificate:

```math
(\lambda-\mu(a))_+
\sum_UH_U(a)^*H_U(a)
\preceq
\sum_UY_U^*(\overline A Y)_U
-\lambda\sum_UY_U^*Y_U.                                       \tag{7.8}
```

This is the minimal theorem that any nonconstant-root matching-scheme
attempt must strengthen: its right side is algebraic and independent of
`a`, while its left root mass is the matrix signed-hafnian sum.

### 7.2 Why an irreducible matrix fiber does not avoid compression

If the construction is `S_n`-equivariant and the input fiber `H` is an
irreducible `S_n` representation, then the constant-root operators `P` and
`R` commute with that representation.  Schur's lemma makes both scalar on
`H`.  Consequently every unit vector has the same quotient leakage, and
the normalized trace kernel

```math
\kappa_{\rm tr}(z)={1\over\dim H}\operatorname{tr}K(z)           \tag{7.9}
```

is a scalar certificate with exactly the same `lambda`, root mass, leakage,
and bound.  For a reducible fiber, Theorem 5.1 selects the least generalized
block/direction and reaches the same conclusion.  Hence the matching orbit's
large size evades the support obstruction, but an irreducible matrix fiber
does **not** evade the lossless scalar-compression obstruction when its root
is constant.

There is also a simple positivity obstruction to trying to cancel every
boundary sum.  From (3.2) and (7.4),

```math
\sum_{d\in D}Q_d=\sum_UY_U^*Y_U.                                \tag{7.10}
```

Since every `Q_d` is positive semidefinite and `Q_0=P`, the condition
`Y_U=0` for every `U` forces `P=0`.  Thus a nontrivial irreducible orbit
cannot kill its hafnian aggregates by representation-theoretic
cancellation while retaining operator Fourier positivity.

### 7.3 Scalar matching amplitudes lose the full leading scale

The exact hafnian twirl gives a scalable obstruction to the first scalar
square implementation.  Suppose `V_M=w_l>=0` for every
`M in M_l`.  Put

```math
p_l=(2l-1)!!,
\qquad
N_l=\binom n{2l}p_l,
\qquad
\alpha_l=\sqrt{N_l}\,w_l.                                      \tag{7.11}
```

The radial matching super-eigenvector condition is

```math
c_{l-1}\alpha_{l-1}+c_l\alpha_{l+1}\ge\lambda\alpha_l,          \tag{7.12}
```

with `c_l` from (7.2).  At the all-positive root, (7.6) is

```math
{T_{\mathbf1}\over|C|}
=D:=\sum_lp_l\alpha_l^2.                                       \tag{7.13}
```

The quotient aggregate gains the square root of the matching-fiber
multiplicity.  In normalized radial coordinates
`beta_l=sqrt(p_l) alpha_l`, its off-diagonal coefficient is

```math
\overline c_l
={\sqrt{\binom{n-2l}{2}\binom{2l+2}{2}}\over E}
=\sqrt{2l+1}\,c_l.                                              \tag{7.14}
```

Since `p_{l+1}=(2l+1)p_l`, the quotient Rayleigh numerator is

```math
\rho_{\rm quot}D
=2\sum_lp_{l+1}c_l\alpha_l\alpha_{l+1},                         \tag{7.15}
```

and `J/T_1=rho_quot-lambda`.  Multiplying (7.12) by
`p_l alpha_l` and summing gives

```math
\lambda D
\le\sum_l(p_l+p_{l+1})c_l\alpha_l\alpha_{l+1}.                  \tag{7.16}
```

Therefore

```math
2\lambda-\rho_{\rm quot}
\le {2\sum_lp_lc_l\alpha_l\alpha_{l+1}\over D}.                 \tag{7.17}
```

The right side is the Rayleigh quotient of a tridiagonal path in the
`beta` coordinates, with off-diagonal entries

```math
d_l={c_l\over\sqrt{2l+1}}\le{1\over\sqrt E}.                    \tag{7.18}
```

Its spectral radius is at most `2 max_l d_l`.  Finally,
all scalar autocorrelation coefficients are nonnegative, so
`T_a<=T_1`; also `J>=0`.  The exact rooted bound is consequently

```math
\boxed{
\lambda-{J\over T_a}
\le\lambda-{J\over T_{\mathbf1}}
\le {2\over\sqrt E}=O(n^{-1}).}                                \tag{7.19}
```

Thus the scalar layer-constant matching square has enough hidden support and
the right internal eigenvalue, but cut-code fiber multiplicity destroys the
entire `n^(-1/2)` leading scale.  This is the matching-scheme version
of the independently paid scalar-channel obstruction.

### 7.4 Direct matching Fourier support is too small at every rank

There is a stronger operator obstruction before radialization.  Suppose

```math
\mathbf K(a)=\sum_{M\ {\rm partial\ matching}}\chi_M(a)Q_M,
\qquad Q_M\succeq0,                                             \tag{7.20a}
```

is pointwise positive semidefinite, and suppose the empty-index coefficient
of the complete remainder is positive:

```math
B_\varnothing
={1\over E}\sum_eQ_e-\lambda Q_\varnothing\succeq0.             \tag{7.20b}
```

Fix a vertex `i` and average `K` over all cube variables outside its
star.  Pointwise positivity is preserved.  A partial matching contained in
one star has size at most one, so the resulting matrix-valued marginal is

```math
Q_\varnothing+\sum_{e\ni i}x_eQ_e\succeq0
\qquad(x_e\in\{\pm1\}).                                        \tag{7.20c}
```

Set every `x_e=-1` and sum (7.20c) over the `n` vertices.  Since
each edge occurs twice,

```math
2\sum_eQ_e\preceq nQ_\varnothing.                               \tag{7.20d}
```

Combining (7.20b) and (7.20d), on the support of
`Q_emptyset`,

```math
\lambda Q_\varnothing
\preceq {1\over E}\sum_eQ_e
\preceq {n\over2E}Q_\varnothing
={1\over n-1}Q_\varnothing.                                    \tag{7.20e}
```

Hence

```math
\boxed{\lambda\le{1\over n-1}.}                                 \tag{7.20f}
```

This closes *arbitrary-rank*, noncommuting kernels whose Fourier
coefficients themselves are supported on partial matchings.  No
`S_n` symmetry or scalar-channel decomposition is used.

For completeness, lossless compression gives the following scalar radial
specialization.  Choose numbers

```math
q_0,q_1,\ldots,q_{\lfloor n/2\rfloor}\ge0
```

and define the `S_n`-radial matching Fourier polynomial

```math
K_q(a)
=\sum_{l=0}^{\lfloor n/2\rfloor}
q_l\sum_{M\in\mathcal M_l}\chi_M(a).                            \tag{7.20}
```

It is positive definite because all Fourier coefficients are nonnegative.
It is pointwise nonnegative exactly when

```math
K_q(a)\ge0\qquad\text{for every edge signing }a.                \tag{7.21}
```

The complete remainder `(tau-lambda)K_q` is positive definite exactly when

```math
lq_{l-1}+\binom{n-2l}{2}q_{l+1}
\ge\lambda E q_l
\quad(0\le l\le\lfloor n/2\rfloor),                             \tag{7.22}
```

with missing terms set to zero.  At nonmatching Fourier indices the
coefficient inequality is automatic.

No nonempty matching lies in `D`, so this kernel has constant root.  More
strongly, if `d in D\setminus{0}`, no one-edge toggle of `d` is a matching:
if `d+e=M` were a matching, `partial d=0` would force
`partial M=partial e`, hence `M=e` and `d=0`.  Therefore

```math
\sum_{d\in D}(A_Eq)_d=(A_Eq)_0=q_1.                             \tag{7.23}
```

The exact bound is consequently

```math
\boxed{
\mu(a)\ge2\lambda-{q_1\over q_0}.}                              \tag{7.24}
```

In particular, the following finite statement would reach the sharp leading
scale:

```math
K_q(a)\ge0\ \forall a,
\qquad
(7.22),
\qquad
{q_1\over q_0}=\lambda={1-o(1)\over\sqrt n}.                    \tag{7.25}
```

This state is genuinely smaller than full parent maximization: it consists
of `O(n)` radial coefficients and one universal pointwise positivity theorem
for a matching polynomial.  However, one star marginal disproves (7.25).

Normalize `K_q` to a probability density on the edge cube.  Fix a vertex and
marginalize to its `n-1` incident edge variables `x_1,...,x_{n-1}`.  The only
Fourier sets contained in this star that are matchings are the empty set and
the singletons.  Therefore the marginal density is exactly

```math
2^{-(n-1)}
\left(1+{q_1\over q_0}\sum_{j=1}^{n-1}x_j\right).               \tag{7.26}
```

Evaluating this nonnegative density at `x_1=...=x_{n-1}=-1` gives the
exact ceiling

```math
\boxed{{q_1\over q_0}\le{1\over n-1}.}                           \tag{7.27}
```

But the `l=0` instance of (7.22) is

```math
{q_1\over q_0}\ge\lambda.                                       \tag{7.28}
```

Consequently every pointwise-nonnegative kernel whose Fourier support is
contained in the partial matchings satisfies

```math
\boxed{\lambda\le{1\over n-1}=o(n^{-1/2}).}                      \tag{7.29}
```

This is a scalable no-go, not finite evidence.  It closes the direct radial
Fourier-matching kernel even though its support has the correct
`exp(Theta(n log n))` size.  Thus support size is necessary but far from
sufficient: pointwise root positivity detects the large intersecting family
of star edges.

## 8. Scope and research judgment

The proved obstruction has five levels.

1. **No matrix advantage at constant root.**  This is unconditional inside
   the operator positive-definite cone (1.1)--(1.3): Theorem 5.1 preserves
   the exact bound under one scalar compression.
2. **Minimum Fourier support.**  Any nonzero compressed certificate with
   `lambda=Theta(n^(-1/2))` needs `exp(Omega(n log n))` Fourier support.
3. **Minimum hidden rank for the natural multi-transversal square.**  Under
   the explicit termwise orthogonality condition (6.6), hidden dimension
   must also be `exp(Omega(n log n))`.
4. **Scalar matching-square leakage.**  Even though linear matching layers
   have the correct internal eigenvalue and orbit size, their exact cut-fiber
   multiplicity reduces the rooted certificate to `O(1/n)`, by (7.19).
5. **Direct matching Fourier support at every rank.**  The operator star
   marginal (7.20a)--(7.20f) forces `lambda<=1/(n-1)` without
   radialization, commutativity, or a scalar-channel payment.

The last conclusion must not be overstated.  The algebraic condition
`Q_d=0` allows cancellation among the terms in (3.2), so (6.6) is sufficient
but not necessary.  There is no general rank lower bound here for such
cancelling implementations.  Nor does the theorem exclude an implicitly
represented `exp(Theta(n log n))` harmonic or `S_n` hierarchy, a nonconstant
root mass `T_a`, or a nonabelian certificate outside the operator-Fourier
cone.

It closes both obvious ways to use the partial-matching scheme: scalar
matching amplitudes are killed by quotient leakage, while Fourier
coefficients directly supported on matchings are killed by one star
marginal.  The only matching-based door left by these theorems is a
genuinely noncommuting amplitude whose autocorrelation coefficients occupy
broader symmetric differences of matchings and whose boundary-fiber
cancellation is not equivalent to scalar positive channels.

A viable continuation must therefore use at least one genuinely new
ingredient:

* nonconstant root mass controlled by an algebraic inequality rather than
  made constant;
* cancellations implementing (4.1) that are not termwise orthogonality and
  whose huge difference-of-matching Fourier support has an implicit
  `S_n` description;
* an indefinite intermediate object whose *complete* final Gram remainder
  is positive but whose root kernel is outside (1.2); or
* a nonabelian moving-stabilizer theorem that does not reduce to the
  constant-root operator cone above.

No generic numerical search is justified before one of these alternatives
is stated as a precise theorem-level construction.
