# Row-magnitude fibres cannot lower conference pressure in the full conference high-temperature range

**Status.** Task-local theorem report.  This strengthens the small-temperature
row-magnitude theorem to the full range `0<beta<sqrt(2)/6`, for the one-sided
conclusion relevant to favorable basins.  It proves that the conditioned
pressure cannot have a smaller leading rate.  It does not claim that a large
population spike cannot raise the pressure.

## 1. Setup and statement

Let `A_r` be a symmetric conference signing, let
`v_r in {+-1}^r`, and let

```math
E_r=\{R\in\{+-1\}^r:|\langle R,v_r\rangle|\in I_r\},
\qquad
p_r=2^{-r}|E_r|\ge p_0>0.                         \tag{MF.1}
```

Let `B_r` have independent rows, each uniform on `E_r`.  For
`epsilon in {+-1}` put

```math
f_(epsilon,r)(B)=
\log\left[2^{-2r}\sum_(x,y)
\cosh\left\{{\beta\over\sqrt{2r}}
\big(H_A(x)+\epsilon H_A(y)+x^TBy\big)\right\}\right]. \tag{MF.2}
```

Recall

```math
h_\beta=2\psi(\beta/\sqrt2)+{\beta^2\over4}.
```

### Theorem MF.1 (full-range no-gain for constant-mass magnitude fibres)

For every `p_0>0`, every fixed

```math
0<\beta<{\sqrt2\over6},                              \tag{MF.3}
```

and both orientations, uniformly over (MF.1),

```math
\boxed{
\mathbb E\left[\left(h_\beta-{f_(epsilon,r)(B_r)\over r}\right)_+\right]
\longrightarrow0.}                                   \tag{MF.4}
```

In particular, for every fixed `eta>0`,

```math
\boxed{
\Pr\{f_(epsilon,r)(B_r)\le(h_\beta-\eta)r\}
\longrightarrow0.}                                   \tag{MF.5}
```

Thus no such exact speed-`r` fibre reaches the smaller conference child
target `tau_beta=h_beta-gamma(beta)`.  The statement is deliberately
one-sided: a sufficiently strong tail constraint can create a finite-rank
mean-field instability and raise the pressure.

## 2. Project away the exceptional population direction

Column switching permits taking `v_r=1`: the same change of variables
conjugates the second child by a diagonal sign matrix, which remains a
conference signing and has the same uniform-bridge pressure theorem.  Write

```math
P={11^T\over r},\qquad B^\circ=B(I-P),\qquad B^\parallel=BP.
                                                               \tag{MF.6}
```

Every row of `B^circ` has sum zero, while `B^parallel` has rank at most one.
For a row `R` from (MF.1), with `S=<R,1>`,

```math
\mathbb E[S^2\mid E_r]\le {r\over p_0}.               \tag{MF.7}
```

Consequently

```math
\mathbb E\|B_r1\|_2
\le\left(r\,\mathbb E[S^2\mid E_r]\right)^{1/2}
\le {r\over\sqrt{p_0}}.                              \tag{MF.8}
```

The issue is therefore not the size of the rank-one part, whose *scaled
nuclear norm* is bounded in mean, but whether `B^circ` has the sharp bulk
operator edge.  The next lemma proves exactly what is needed by coupling,
without importing a dependent-coordinate Bai--Yin theorem.

## 3. The projected layer coupling has subcritical operator cost

Use the row-layer coupling from the small-temperature magnitude theorem.
Thus `W_r` is a fully iid Rademacher bridge, `B_r` has the conditioned law,
the rows are coupled independently, and the plus sets in each row are
nested.  Put

```math
D=B_r-W_r,\qquad D^\circ=D(I-P).                       \tag{MF.9}
```

### Lemma MF.2 (the projected coupling error is `o(sqrt r)` in operator norm)

Under this coupling,

```math
\boxed{\|D^\circ\|_(op)=o_\Pr(\sqrt r),}              \tag{MF.10}
```

and the exceptional event can be chosen to have probability `O(r^(-10))`
plus `exp(-c_(p_0)sqrt(r))`.  Moreover,

```math
\mathbb E\|D^\circ\|_F=O_(p_0)(r^(3/4))=o(r).         \tag{MF.11}
```

**Proof.**  Let `K_i,K_i'` be the two plus counts and
`d_i=|K_i-K_i'|`.  Ordinary binomial concentration, followed by division by
`p_r>=p_0` for the conditioned count, gives

```math
\Pr\{\max_i d_i>2r^(3/4)\}
\le e^{-c_(p_0)\sqrt r}                               \tag{MF.12}
```

for all large `r` (a harmless polynomial factor is absorbed).

Conditional on the counts, the changed coordinates in row `i` form a
uniform `d_i`-subset `T_i`.  Up to an irrelevant common sign,

```math
D_i^\circ=2\left(1_(T_i)-{d_i\over r}1\right),
\qquad
\|D_i^\circ\|_2^2=4d_i(1-d_i/r).                     \tag{MF.13}
```

Permutation symmetry gives the exact covariance identity

```math
\mathbb E[(D_i^\circ)^TD_i^\circ\mid d_i]
={4d_i(r-d_i)\over r(r-1)}(I-P).                     \tag{MF.14}
```

On (MF.12)'s good event, the independent positive-semidefinite summands
`(D_i^circ)^T D_i^circ` have norm at most `8r^(3/4)`, while the norm of their
conditional expected sum is at most `8r^(3/4)`.  Matrix Chernoff (or matrix
Bernstein after subtracting the conditional expectation), at threshold
`C r^(3/4)log r`, yields conditional failure probability `O(r^(-10))` after
increasing `C`: the summand and expected-sum bounds are both
`O(r^(3/4))`, so the Bernstein exponent is a constant multiple of `log r`.
Hence

```math
\|D^\circ\|_(op)
\le C r^(3/8)\sqrt{\log r}=o(\sqrt r),                \tag{MF.15}
```

which proves (MF.10).  Finally projection is a Frobenius contraction and
the layer coupling gives

```math
\mathbb E\|D^\circ\|_F
\le\mathbb E\|D\|_F
\le\left(4r\,\mathbb E d_i\right)^{1/2}
=O_(p_0)(r^(3/4)).
```

This proves (MF.11). `square`

Since `||W_r(I-P)||op<=||W_r||op` and the rectangular Rademacher edge is
`(2+o(1))sqrt(r)`, Lemma MF.2 gives

```math
\boxed{
\|W_r^\circ\|_(op),\ \|B_r^\circ\|_(op)
\le(2+o_\Pr(1))\sqrt r.}                             \tag{MF.16}
```

## 4. The projected bridge retains the uniform pressure

For a real bridge `C`, write

```math
X(C)={\beta\over\sqrt{2r}}
\begin{pmatrix}A_r&C\\C^T&\epsilon A_r\end{pmatrix}.
                                                               \tag{MF.17}
```

Fix `beta` as in (MF.3).  Choose `delta>0` and `kappa<1/2` such that

```math
{\beta(3+\delta)\over\sqrt2}<\kappa.                 \tag{MF.18}
```

Equations (MF.10), (MF.16), and `||A_r||op=sqrt(r-1)` put both
`X(W_r^circ)` and `X(B_r^circ)` in this common strict operator ball with
probability tending to one, with a complement negligible even after
multiplication by the crude pressure bound.

The archived high-temperature pressure-stability theorem then gives

```math
|f(B_r^\circ)-f(W_r^\circ)|
\le {K_\kappa\beta\over\sqrt2}\|D^\circ\|_F=o_(L^1)(r).  \tag{MF.19}
```

The iid bridge and its projection also have the same leading pressure.  On
their common regular event, pressure stability and

```math
\|W_rP\|_*={\|W_r1\|_2\over\sqrt r},
\qquad \mathbb E\|W_r1\|_2\le r                     \tag{MF.20}
```

give

```math
\mathbb E|f(W_r)-f(W_r^\circ)|=O_(\beta,\kappa)(1)+o(r). \tag{MF.21}
```

The uniform conference theorem and (MF.19)--(MF.21) therefore imply

```math
\boxed{
{f_(epsilon,r)(B_r^\circ)\over r}\longrightarrow h_\beta
\quad\hbox{in probability and in }L^1.}              \tag{MF.22}
```

## 5. Convexity restores an arbitrarily large rank-one spike one-sidedly

Let

```math
Y={\beta\over\sqrt{2r}}
\begin{pmatrix}0&B_rP\\PB_r^T&0\end{pmatrix},
\qquad g(s)=\log\overline Z_(2r)(X(B_r^\circ)+sY).    \tag{MF.23}
```

The log-cosh partition pressure is convex in its interaction matrix, hence
`g` is convex on the whole real line.  Only its base point needs to be in
the strict high-temperature ball.  At that point the archived covariance
bound and nuclear/operator duality give

```math
g'(0)\ge-{K_\kappa\over2}\|Y\|_*
=-K_\kappa{\beta\over\sqrt{2r}}
  {\|B_r1\|_2\over\sqrt r}.                           \tag{MF.24}
```

The global supporting-line inequality for a convex function now yields

```math
f(B_r)=g(1)
\ge f(B_r^\circ)
-{K_\kappa\beta\over\sqrt2\,r}\|B_r1\|_2.          \tag{MF.25}
```

Crucially, no operator bound on the full endpoint `X(B_r)` is used.  By
(MF.8), the expected loss in (MF.25) is at most
`K_kappa beta/sqrt(2p_0)=O(1)`.  The negligible failure of the base regular
event contributes `o(r)` by the crude `O_beta(r^(3/2))` pressure bound.
Combining this one-sided estimate with (MF.22) proves (MF.4), and Markov
proves (MF.5). `square`

## 6. Why a raw sharp edge is the wrong full theorem

The conditioned row covariance has one eigenvalue

```math
\alpha_r={1\over r}\mathbb E[S^2\mid E_r]
```

and bulk eigenvalue `(r-alpha_r)/(r-1)`.  A fixed Gaussian-tail event can
make `alpha_r` arbitrarily large while retaining constant one-row mass.  In
fact

```math
\|B_r\|_(op)\ge {\|B_rv_r\|_2\over\sqrt r},
```

whose typical leading coefficient is `sqrt(alpha_r)`.  Thus no uniform
`(2+o(1))sqrt(r)` edge is possible for the *unprojected* class.  The theorem
works because the exceptional direction is removed before applying the
strict-high-temperature comparison, and convexity makes its restoration
one-sided and finite-cost at the pressure scale.

## 7. Scope

This closes the declared constant-mass, single-magnitude, row-product class
as a lower-pressure conference basin throughout the full range in which the
uniform conference parent theorem is available.  It does not show equality
of the conditioned pressure rate: a strong population spike may raise it.
It also does not cover non-product correlations between rows, several
distinguished directions whose span grows with `r`, or row events of
vanishing one-row mass.
