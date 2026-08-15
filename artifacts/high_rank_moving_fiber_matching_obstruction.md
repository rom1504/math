# High-rank moving fibers and the partial-matching obstruction

Date: 2026-08-15.

Status: **proved operator extension, exact matching-fiber reduction, and
scalable no-go theorems for the two direct partial-matching realizations**.
No bound on `M_n` is improved.

The purpose of this note is to test the suggestion that a growing
`S_n`-representation might compress the augmented cut-code covering problem.
The calculation keeps the complete matrix-valued Gram remainder intact.  It
does not bound scalar channels separately.  The outcome is useful but
negative: partial matchings have the required hidden entropy and the correct
order of transition eigenvalue, yet their direct implementations cannot reach
the coefficient required for convergence.

## 1. Operator-valued moving-Gram theorem

Let

```math
G=\mathbb F_2^E,qquad E=\binom n2,qquad
C=\mathcal C_n^+,qquad D=C^\perp,
```

and write

```math
\tau(z)={1\over E}\sum_{e=1}^E(-1)^{z_e},qquad
\mu(a)=\max_{c\in C}\tau(a+c)={Q(a)\over E}.
```

Let `K:G -> Herm(H)` be a finite-dimensional operator kernel.  With
normalized Fourier coefficients `Q_R=K_hat(R)`, assume jointly that

```math
K(z)\succeq0,qquad Q_R\succeq0,qquad
B_R:=(A_EQ)_R-\lambda Q_R\succeq0,                 \tag{1.1}
```

where

```math
(A_EQ)_R={1\over E}\sum_eQ_{R+e}.
```

Thus `K` is pointwise positive, while both `K` and the complete remainder

```math
F(z)=(\tau(z)-\lambda)K(z)
```

are operator-valued positive-definite kernels.  Define

```math
T_a=\sum_{c\in C}K(a+c),qquad
S_a=\sum_{c\in C}F(a+c),qquad
J=\sum_{c\in C}F(c).
```

> **Theorem 1.1 (operator moving-fiber inequality).** For every `a in G`,
>
> ```math
> \boxed{(\lambda-\mu(a))_+T_a\preceq J.}          \tag{1.2}
> ```
>
> On the support of `T_a`,
>
> ```math
> \mu(a)\ge\lambda-
> \lambda_{\min}(T_a^{-1/2}JT_a^{-1/2}).           \tag{1.3}
> ```

**Proof.** Subgroup Fourier orthogonality and (1.1) give

```math
S_a=|C|\sum_{d\in D}\chi_d(a)B_d,qquad
J=|C|\sum_{d\in D}B_d,qquad
-J\preceq S_a\preceq J.                            \tag{1.4}
```

Pointwise positivity and `tau(a+c)<=mu(a)` give

```math
S_a\preceq(\mu(a)-\lambda)T_a.                    \tag{1.5}
```

Combine the lower half of (1.4) with (1.5) when `mu(a)<lambda`.
The other case is automatic. `square`

This is a genuinely joint Loewner inequality: scalarization occurs nowhere
in its proof.

## 2. Constant roots nevertheless scalarize exactly

For a matrix-square realization, write

```math
G_V(z)=\sum_R\chi_R(z)V_R,qquad K(z)=G_V(z)^*G_V(z).
```

Its exact Fourier conditions are

```math
Q_T=\sum_RV_R^*V_{R+T}\succeq0,                  \tag{2.1}
```

```math
B_T=\sum_RV_R^*((A_EV)_{R+T}-\lambda V_{R+T})
=(A_EQ)_T-\lambda Q_T\succeq0.                  \tag{2.2}
```

Suppose the cut root is made constant by

```math
Q_d=0qquad(d\in D\setminus\{0\}).                \tag{2.3}
```

Then `T_a=|C|P` for every `a`, where `P=Q_0`.  Aggregate amplitudes over
quotient fibers,

```math
Y_x=\sum_{R\in x}V_R,qquad x\in G/D,
```

and let `Abar` be normalized quotient adjacency.  Direct summation gives

```math
{J\over|C|}=\mathcal R-\lambda P,qquad
\mathcal R=\sum_xY_x^*(\overline A Y)_x.           \tag{2.4}
```

Consequently Theorem 1.1 gives

```math
\mu(a)\ge2\lambda-\rho_{\min},qquad
\rho_{\min}=\lambda_{\min}
(P^{-1/2}\mathcal RP^{-1/2}).                     \tag{2.5}
```

Choose a generalized eigenvector attaining `rho_min` and compress only now:

```math
\kappa_u(z)=u^*K(z)u.
```

It is a scalar pointwise-positive, positive-definite kernel; its complete
remainder is positive definite; it obeys (2.3); and it attains exactly (2.5).
Thus matrix rank supplies no intrinsic advantage in the constant-root cone.
For an irreducible `S_n` fiber, the same conclusion follows immediately from
Schur's lemma.  A high-rank representation can still be valuable as an
implicit construction of the resulting scalar kernel, but it must provide
new algebraic control rather than strength from rank alone.

Every nonzero scalar compression has Fourier weights `q>=0` satisfying

```math
A_Eq\ge\lambda q.                                  \tag{2.6}
```

The verified Bollobas--Lee--Letzter cube-subgraph theorem therefore implies

```math
\lambda\ge {c\over\sqrt n}
\quad\Longrightarrow\quad
|\operatorname{supp}q|\ge\exp(\Omega_c(n\log n)). \tag{2.7}
```

This is the quantitative reason a growing hidden representation is necessary.
For a termwise-orthogonal multi-transversal square, it also forces hidden rank
`exp(Omega(n log n))`.  That rank conclusion uses termwise orthogonality;
algebraic cancellation in (2.1) is not covered.

## 3. Exact partial-matching fiber

Let `M_l` be the `l`-edge partial matchings of `K_n`.  For `l=alpha n`,
`0<alpha<1/2`,

```math
|\mathcal M_l|
={n!\over(n-2l)!2^ll!}=\exp(\Theta(n\log n)).       \tag{3.1}
```

The normalized add/delete adjacency between layers has coefficient

```math
c_l={\sqrt{(l+1)\binom{n-2l}{2}}\over E}
=\Theta(n^{-1/2}).                                 \tag{3.2}
```

Thus this orbit passes the support-size and order-of-magnitude tests.

For a matching `M`, its cut-code syndrome is its covered vertex set together
with the parity of `|M|`.  Hence

```math
M+N\in D\quad\Longleftrightarrow\quad
\partial M=\partial N;                             \tag{3.3}
```

the symmetric difference is then a union of alternating even cycles.  For
matrix amplitudes supported on matchings, define

```math
H_U(a)=\sum_{\partial M=U}\chi_M(a)V_M,qquad
Y_U=H_U(0).
```

The exact cut twirl and numerator are

```math
\boxed{{T_a\over|C|}=\sum_UH_U(a)^*H_U(a),}         \tag{3.4}
```

```math
\boxed{{J\over|C|}
=\sum_UY_U^*(\overline A Y)_U
-\lambda\sum_UY_U^*Y_U.}                         \tag{3.5}
```

Equations (3.4)--(3.5) are the exact nonconstant-root matching interface.
For scalar layer weights, (3.4) is a weighted sum of squared signed principal
hafnians.

## 4. The radial matching square loses the leading scale

Take a nonnegative scalar amplitude `V_M=w_l` on every `M in M_l`.  Put

```math
p_l=(2l-1)!!,qquad
N_l=\binom n{2l}p_l,qquad
\alpha_l=\sqrt{N_l}w_l.
```

The matching super-eigenvector condition is

```math
c_{l-1}\alpha_{l-1}+c_l\alpha_{l+1}
\ge\lambda\alpha_l.                                \tag{4.1}
```

At the all-positive signing,

```math
{T_{\mathbf1}\over|C|}=D_0:=\sum_lp_l\alpha_l^2.   \tag{4.2}
```

Cut-fiber aggregation changes the quotient path coefficient to

```math
\overline c_l=\sqrt{2l+1}\,c_l.                    \tag{4.3}
```

Using `p_{l+1}=(2l+1)p_l`, multiply (4.1) by
`p_l alpha_l` and sum.  If `rho_quot` is the resulting quotient Rayleigh
quotient, this gives

```math
2\lambda-\rho_{\rm quot}
\le {2\sum_lp_lc_l\alpha_l\alpha_{l+1}\over D_0}. \tag{4.4}
```

In the coordinates `beta_l=sqrt(p_l)alpha_l`, the right side is the Rayleigh
quotient of a path with coefficients

```math
d_l={c_l\over\sqrt{2l+1}}\le {1\over\sqrt E}.       \tag{4.5}
```

Moreover `T_a<=T_1`, because every autocorrelation coefficient is
nonnegative.  Therefore, whenever `T_a>0`,

```math
\boxed{
\lambda-{J\over T_a}
\le\lambda-{J\over T_{\mathbf1}}
\le {2\over\sqrt E}=O(n^{-1}).}                   \tag{4.6}
```

When `T_a=0`, Theorem 1.1 is vacuous.  Thus the scalar matching square has
the correct internal transition scale but loses it completely to boundary
multiplicity.

## 5. Direct matching Fourier support fails at every rank

The stronger obstruction needs neither radiality nor commutativity.  Suppose

```math
K(a)=\sum_{M\ {m matching}}\chi_M(a)Q_M,qquad
Q_M\succeq0,                                       \tag{5.1}
```

is pointwise positive semidefinite, and suppose the empty Fourier coefficient
of the complete remainder is positive:

```math
{1\over E}\sum_eQ_e-\lambda Q_\varnothing\succeq0. \tag{5.2}
```

Fix a vertex `i` and average `K` over every edge outside its star.  A matching
contained in a star has size at most one, so pointwise positivity gives

```math
Q_\varnothing+\sum_{e\ni i}x_eQ_e\succeq0
\qquad(x_e\in\{\pm1\}).                            \tag{5.3}
```

Set every `x_e=-1` and sum over vertices.  Each edge occurs twice:

```math
2\sum_eQ_e\preceq nQ_\varnothing.                  \tag{5.4}
```

Combining (5.2) and (5.4) yields, on the nonzero support of
`Q_emptyset`,

```math
\boxed{\lambda\le {n\over2E}={1\over n-1}.}        \tag{5.5}
```

This closes arbitrary-rank operator kernels whose Fourier coefficients are
supported on partial matchings.  No `S_n` symmetry, simultaneous
diagonalization, or scalar-channel payment is used.

## 6. A representation-graph ceiling

The incidence operator from `M_l` to `M_{l+1}` has largest singular value

```math
\sqrt{(l+1)\binom{n-2l}{2}}.
```

Every common-irreducible or multiplicity block is a compression of this
operator.  Consequently every add/delete partial-matching representation
graph has

```math
\lambda
\le2\max_l c_l
=\left({4\over3\sqrt3}+o(1)\right){1\over\sqrt n}. \tag{6.1}
```

The maximum occurs at `l/n -> 1/6`.  Since

```math
{4\over3\sqrt3}=0.769800\ldots<1,
```

even a perfect root estimate in this architecture cannot reach the
conference coefficient `1/sqrt(n)` needed by the convergence criterion.

As a concrete positive check, the common standard module `1^perp` embeds
isometrically in every matching layer and gives a valid unsplit
moving-projection Gram remainder with `lambda=Theta(n^(-1/2))`.  Its cut root
is exactly `(n-1)^2` times the scalar signed-hafnian statistic.  Oddness of
every principal hafnian gives a uniform floor, but the resulting finite tests
through order eight are vacuous beyond the first layer.  The representation
therefore adds no new root channel.

## 7. Scope and literature boundary

The perfect-matching association scheme is the multiplicity-free Gelfand-pair
scheme described by Srinivasan,
[*The perfect matching association scheme*](https://arxiv.org/abs/1807.00481).
Its orbital relations are indexed by the alternating-cycle structure of two
matchings.  This precisely matches (3.3), but ordinary orbital
diagonalization is unrooted: an arbitrary signing inserts a diagonal twist in
(3.4), which the scheme eigenvalues do not control.

The following classes are now closed:

1. constant-root operator rank as a source of strength by itself;
2. scalar `S_n`-radial matching amplitude squares;
3. arbitrary-rank kernels with Fourier support directly on matchings; and
4. every pure add/delete matching representation as a route to coefficient
   `1` in the convergence theorem.

The proof does **not** close a noncommuting amplitude whose autocorrelation
occupies broader symmetric differences of matchings, including alternating
paths and cycles.  It also does not close a nonconstant-root theorem that
uniformly controls the matrix signed-hafnian sum (3.4), or a moving
representation built from graph orbits other than partial matchings.

Those surviving formulations did not yield a polynomially closed root-mass
bound in this attempt.  They remain open, but are not presently demonstrated
to be simpler than the original cut-code tail.
