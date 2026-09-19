# Wave 48B: row-truncated ground-lift coarea

The main result below is a **Verified exact sharpening** of the proposed
ground-lift criterion.  It applies after an arbitrary center truncation,
including the project-row truncation.  It removes the need to prove a strict
`(1-eta)` coarea improvement: an elementary argument makes a non-strict
normalized boundary bound extract degree `Theta(1/n)`, and the balanced-slice
FKN theorem upgrades this to constant degree whenever the first two Johnson
eigenvalues are separated by a constant fraction of the spectral gap.  Both
are much stronger than the required `exp(-O(L_0))` degree.

The project-row *nonzero-mass* clause remains separate.  A scalable genuine
signing (not a minimizer) shows that it cannot follow from ground incidence
and row truncation alone.

All finite and symbolic checks are in
`tmp/row_truncated_coarea_r48_check.py`; its output is
`tmp/row_truncated_coarea_r48.out`.

## 1. A sharp first-level bound for a small slice family

Let `U_m` be uniform on the `m`-slice of `[n]`, put `p=m/n`, and let
`f=1_F`, `a=E f`.  If `P_j` denotes the Johnson harmonic projection and
`W_j=||P_j f||_2^2`, then

```math
\boxed{W_1\le (n-1)a^2.}
\tag{R48B.1}
```

Indeed, for `X_i=1{i in S}` and

```math
\mu_i=\Pr(i\in S\mid F),\qquad
v_i=\mathbb E[f(X_i-p)]=a(\mu_i-p),
```

we have `sum_i v_i=0`.  The slice covariance on the sum-zero coordinate
space is scalar with eigenvalue `p(1-p)n/(n-1)`.  Therefore the level-one
projection has the exact norm

```math
W_1=
\frac{(n-1)a^2}{p(1-p)n}
\sum_{i=1}^n(\mu_i-p)^2.
\tag{R48B.2}
```

The vector `mu` lies in `[0,1]^n` and has sum `m=np`.  Convexity gives

```math
\sum_i(\mu_i-p)^2\le np(1-p),
```

with equality at a zero-one vector of weight `m`, proving (R48B.1).  This
elementary estimate is especially strong when `a=o(1/n)`: an exponentially
small Boolean family has negligible level-one mass.

## 2. Exact Johnson sharpening of the coarea ratio

For the down-up kernel `K_ell`, write

```math
\lambda_1=\frac{\ell(n-m)}{m(n-\ell)},\qquad
\delta=1-\lambda_1,
```

and

```math
\lambda_2=
\frac{(\ell)_2(n-m)_2}{(m)_2(n-\ell)_2},qquad
g=\lambda_1-\lambda_2,qquad
\kappa=\frac g\delta.
```

The convention is `lambda_2=0` if a falling factorial vanishes.  For

```math
B=\langle f,(I-K_\ell)f\rangle,
```

the exact spectrum gives

```math
B=\delta(a-a^2)+
\sum_{j\ge2}(\lambda_1-\lambda_j)W_j.
```

Since the eigenvalues decrease with `j`, and since
`sum_{j>=2} W_j=a-a^2-W_1`, (R48B.1) yields the **Verified pointwise bound**

```math
\boxed{
\frac{B}{\delta a}
\ge 1+\kappa-(1+\kappa n)a.
}
\tag{R48B.3}
```

No minimizer property is used here.

Now let `C` be any selector-independent class of centers (in particular,
`C={z:R_2(z)<=R_*}`), and for each center let `f_z,a_z,B_z` be its hard
favorable family, degree, and boundary.  Assume only

```math
D_C:=\mathbb E[1_Ca_z^2]>0.
```

Under the truncated ground-lift law `nu(z) proportional 1_C(z)a_z`, the
normalized boundary is

```math
\mathcal R_C=
\frac{\mathbb E[1_Ca_zB_z]}{\delta\,\mathbb E[1_Ca_z^2]}.
```

Average (R48B.3) under the double-incidence law
`pi_C(z) proportional 1_C(z)a_z^2`.  This gives

```math
\boxed{
\mathcal R_C
\ge 1+\kappa-(1+\kappa n)\mathbb E_{\pi_C}a_z.
}
\tag{R48B.4}
```

Consequently, the formerly proposed strict criterion has the stronger exact
consequence

```math
\boxed{
\mathcal R_C\le1-\eta
\quad\Longrightarrow\quad
\max_{z\in C}a_z
\ge\mathbb E_{\pi_C}a_z
\ge\frac{\kappa+\eta}{1+\kappa n}.
}
\tag{R48B.5}
```

In particular, **the strict factor is unnecessary**:

```math
\boxed{
\mathbb E[1_Ca_zB_z]
\le(1-\lambda_1)\mathbb E[1_Ca_z^2],\quad D_C>0
\ \Longrightarrow\ 
\max_{z\in C}a_z\ge\frac{\kappa}{1+\kappa n}.
}
\tag{R48B.6}
```

If `ell=m-s`, direct cancellation gives

```math
\boxed{
\kappa=
\frac{(m-s)(n-m)(n-2)}
{n(m-1)(n-m+s-1)}.
}
\tag{R48B.7}
```

Thus for fixed `m/n -> p in (0,1)` and `s=o(n)` (including
`s=Theta(n/L_0)`), `kappa=1-o(1)` and (R48B.6) extracts degree
`(1+o(1))/n`.  More generally any scale with `kappa` bounded below extracts
`Theta(1/n)`.  Since `1/n >> exp(-O(L_0))`, (R48B.6), plus nonzero
project-row mass, is a strictly weaker sufficient coarea lemma than
(10.1213) with a prescribed positive `eta`.

The converse obstruction is equally concrete.  If every project-row center
has `a_z=o(1/n)` and `kappa` is bounded below, then (R48B.4) forces
`mathcal R_C >= 1+kappa-o(1)`, not merely a tiny failure above one.  Hence a
successful truncated coarea theorem must create polynomial-degree,
first-harmonic/dictator-scale structure; generic exponentially small
families cannot satisfy it.

There is a further **Verified literature-based upgrade** on a balanced
fixed-density slice.  Since `mathcal R_C` is the `pi_C`-average of
`B_z/(delta a_z)`, `mathcal R_C<=1` gives some positive-degree center in `C`
with

```math
\frac{B_z}{\delta a_z}\le1.
```

For that center, the exact spectral decomposition above implies

```math
\sum_{j\ge2}W_j
\le\frac{a_z^2}{\kappa}.
\tag{R48B.7a}
```

Indeed, the nonconstant baseline contributes `delta(a_z-a_z^2)`, so the
remaining spectral term is at most `delta a_z^2`; it is at least
`(lambda_1-lambda_2) sum_{j>=2}W_j=delta kappa sum_{j>=2}W_j`.

[Filmus's balanced-multislice FKN theorem, Theorem 1](https://arxiv.org/abs/1809.03089)
states, in the two-color case, that if a Boolean `F` has
`epsilon=||F^{>1}||_2^2`, then it differs from a Boolean dictator on at most

```math
4\epsilon+O_\rho(\epsilon^2)
```

of the slice, provided both color densities are at least fixed `rho>0` and
`n>=N(rho)`.  The harmonic normalization is exactly the one used above.
Boolean degree-one functions on the two-color slice are `0`, `1`, `X_i`,
and `1-X_i`, with means `0`, `1`, `p`, and `1-p`.

Fix a ratio window `p in [p_0,p_1] subset (0,1)` and `kappa>=kappa_0>0`.
If `a_z` tended to zero, (R48B.7a) would give
`epsilon<=a_z^2/kappa_0`.  Approximation by the zero dictator would require

```math
a_z\le4a_z^2/\kappa_0+O_{p_0,p_1}(a_z^4),
```

which is impossible for positive sufficiently small `a_z`.  Approximation
by any nonzero dictator is also impossible, since its mean is at least
`min(p_0,1-p_1)>0`, whereas the approximation error is `O(a_z^2)`.
Therefore there is a constant `c_0=c_0(p_0,p_1,kappa_0)>0` such that

```math
\boxed{
D_C>0,\quad \mathcal R_C\le1,\quad\kappa\ge\kappa_0
\quad\Longrightarrow\quad
\max_{z\in C}a_z\ge c_0.
}
\tag{R48B.7b}
```

The quantifiers matter: the constant conclusion uses a density window away
from zero and one, a uniform lower bound on `kappa`, and sufficiently large
`n`.  It does not apply at `ell=1`, where `kappa=Theta(1/n)`, nor at a
degenerating restriction ratio.  The elementary bound (R48B.5) remains
valid without those hypotheses.

## 3. Project-row mass is genuinely independent

The denominator condition `D_C>0` cannot be discarded.  Here is a scalable
actual-signing obstruction to any law-free mass theorem.

Take the all-positive signing

```math
A^+=J-I.
```

For `m>=3`, `Q(A^+[S])=m(m-1)`, and a full spin `z` is a ground lift on
`S` exactly when `z` is constant on `S`.  If `k` coordinates of `z` are
positive, then

```math
a_z=
\frac{\binom{k}{m}+\binom{n-k}{m}}{\binom nm}.
\tag{R48B.8}
```

Also `(A^+)^2=(n-2)J+I`, so

```math
R_2(z)=n+(n-2)(2k-n)^2.
\tag{R48B.9}
```

Fix `m/n -> p>1/2`.  If `a_z>0`, then
`max(k,n-k)>=m`, hence

```math
|2k-n|\ge2m-n=(2p-1)n+O(1)
```

and therefore `R_2(z)=Omega(n^3)`.  For every fixed `c>0`, all such centers
eventually lie outside a project cap `O(n^(9/4-c))`.  Thus

```math
\mathbb E[1_Ca_z^2]=0
```

for the project-row class, despite every selector having ground lifts.
The checker evaluates this separation at `p=3/4,c=1/8` and verifies the
exact formulas.

This is **not an exact-minimizer obstruction**: `A^+` has quadratic norm
`n(n-1)`.  Its scoped conclusion is that ground incidence, cube-average row
`E R_2=n(n-1)`, coarea algebra, and truncation alone cannot prove project
mass.  Any positive theorem must use exact-minimizer geometry to establish
`D_C>0` (or construct a different low-row center law).

## 4. Exact completion criterion for nonzero project mass

There is a strongest straightforward extension criterion which combines the
uniform-completion identity with full-energy Parseval.  Fix a selector `S`,
put `T=[n]\S`, `k=n-m`, and take an oriented child ground `(sigma,y)`.
Replace `A` by `sigma A` in this calculation.  This leaves `A^2`, all row
squares, and every absolute quadratic norm unchanged, and lets us write

```math
y^{\mathsf T}A[S]y=Q_S:=Q(A[S])\ge0.
```

Write

```math
u=A[S]y,\qquad t=A[T,S]y.
```

For a uniform outside word `w in {+-1}^T` and `z=(y,w)`, direct expansion
gives the **Verified exact completion mean**

```math
\boxed{
\mathbb E_w R_2(z)
=\lVert A[:,S]y\rVert_2^2+k(n-1)
=\lVert u\rVert_2^2+\lVert t\rVert_2^2+k(n-1).
}
\tag{R48B.10}
```

Consequently any child ground satisfying

```math
\lVert A[:,S]y\rVert_2^2+k(n-1)\le R_*
\tag{R48B.11}
```

has at least one completion in the row class `C={R_2<=R_*}`, and hence
forces `D_C>0`.

There is an exact shortfall-sensitive sufficient form of (R48B.11).  The
full energy on the outside cube is

```math
E(w)=Q_S+2t^{\mathsf T}w+w^{\mathsf T}A[T]w.
```

Its constant, linear, and quadratic Walsh parts are orthogonal.  Since
`|E(w)|<=q_n=Q(A)` for every completion,

```math
\boxed{
Q_S^2+4\lVert t\rVert_2^2+2k(k-1)
=\mathbb E_wE(w)^2\le q_n^2.
}
\tag{R48B.12}
```

Thus the following is an exact sufficient condition for project mass:

```math
\boxed{
\lVert A[S]y\rVert_2^2
+\frac{q_n^2-Q_S^2-2k(k-1)}4
+k(n-1)\le R_*.
}
\tag{R48B.13}
```

If `Delta=q_n-Q_S`, then
`q_n^2-Q_S^2<=2q_n Delta`.  Therefore a cleaner, slightly stronger-to-assume
package is

```math
\boxed{
\lVert A[S]y\rVert_2^2=O(R_*),
\qquad
\Delta=O(R_*/q_n).
}
\tag{R48B.14}
```

At the project scale `R_*=Theta(n^(9/4-c))` and
`q_n=Theta(n^(3/2))`, the second clause is the concrete near-parent condition

```math
q_n-Q(A[S])=O(n^{3/4-c}).
\tag{R48B.15}
```

The first clause is independent margin regularity.  If
`ell_i=y_i(A[S]y)_i` are the oriented child local margins, child one-spin
stability gives exactly

```math
\ell_i\ge0,\qquad
\sum_{i\in S}\ell_i=Q_S,
\qquad
\lVert A[S]y\rVert_2^2=\sum_i\ell_i^2
\le(m-1)Q_S.
\tag{R48B.16}
```

This does **not** prove the first clause of (R48B.14).  For a fixed-density
child with `Q_S=Theta(n^(3/2))`, (R48B.16) is only `O(n^(5/2))`.

The exact-minimizer operator estimate reaches the same wall:

```math
\lVert A[:,S]y\rVert_2^2
\le m\lVert A\rVert_{op}^2
\le2m q_n=O(n^{5/2}).
\tag{R48B.17}
```

Combining conditional expectation with (R48B.17) guarantees a completion
only at row `O(n^(5/2))`.  Compared with `n^(9/4-c)`, the unresolved factor
is

```math
n^{(5/2)-(9/4-c)}=n^{1/4+c}.
\tag{R48B.18}
```

Completion Parseval and one-spin stability do not repair this exponent.
Using only the scalar child cap, their upper bound has the shape

```math
(m-1)Q_S+\frac{q_n^2-Q_S^2}{4}+O(n^2).
\tag{R48B.19}
```

As a function of `Q_S in [0,q_n]` the displayed main part is concave, so
its best endpoint scale is still `m q_n=Theta(n^(5/2))` (the other endpoint
is `q_n^2=Theta(n^3)`).  This is a **proof-method exponent wall**, not a
lower bound on the true best completion row.  It says that Parseval,
one-spin stability, and the operator estimate, without either the separate
local-margin regularity in (R48B.14) or a new cancellation theorem, cannot
force `D_C>0` at project scale.

The checker independently verifies (R48B.10), (R48B.12), and (R48B.16) on
a nontrivial exact finite completion cube.

## 5. Research judgment

The row-truncated coarea route has a cleaner live target than (10.1213):

1. prove **nonzero project-row double-incidence mass** `D_C>0`; and
2. at one scale with `kappa` bounded below, prove only the non-strict
   inequality (R48B.6).

Those two claims already yield a constant-mass hard exceptional center at a
balanced scale with `kappa` bounded below (and polynomial mass by the
elementary argument without FKN).
The adjacent-port identity (10.1211) would therefore only need to pay the
unimproved spectral-gap constant; no exponentially small strict margin is
needed.  However, neither claim follows formally from minimality, and the
all-positive family shows that the first is a substantive, minimizer-specific
obligation.  No convergence proof or exact-minimizer counterexample is
obtained.
