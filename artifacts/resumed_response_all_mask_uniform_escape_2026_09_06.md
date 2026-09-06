# Uniform strict escape from every nonnegative marked creation mask

Date: 2026-09-06. The finite nonlinear center theorem has been reconstructed
independently by the director and the response track; see
`resumed_response_center_theorem_independent_audit_2026_09_06.md`.
This note supplies the complete infinite-mask approximation and supremum
argument. No iterability or convergence assertion is made.

## 1. Universal variational statement

Let `U` be the marked-tree Gaussian creation isometry, and let

```
C_mark=sup {J(H): H jointly even, 0<=H<=1},
J(H)=E |UH|(1-H).
```

This is the whole nonnegative marked-mask variational closure, not only
central masks or masks produced by a contraction. Set

```
c_inf=liminf_(n->infinity) M_n/n^(3/2),
s0=1/40000000000,
x0=sqrt(19/6),
g_s(x)=E|x+sN|-x=2[s phi(x/s)-x Phi(-x/s)],
delta0=(6/25) g_s0(x0)>0.
```

Then

```
c_inf >= C_mark+delta0 > C_mark.                           (1)
```

The displayed `delta0` is a rigorously specified positive real constant,
not a numerically evaluated decimal improvement. It is extremely small
and is not claimed sharp. The original evaluated lower endpoint remains
`.4306581794055286` unless a separate new numerical certificate is supplied.

## 2. Extending the finite center theorem to arbitrary masks

Fix any jointly even measurable mask `0<=H<=1`. Write

```
W=UH, F=sign(W)(1-H), K=U*F,
s=sqrt(sum_(r=3,5,...,51) |E F h_r(G0)|^2), G0=U1.
```

If `s>0`, the finite center theorem and the approximation below prove

```
c_inf >= Phi(H,F):=E H E_N|K+sN|.                         (2)
```

Choose increasing finite ancestor-closed coordinate sigma-fields containing
`G0`. First condition **both** `H` and `F` on the same finite field. Their
pointwise constraint is preserved:

```
|E[F|A_j]|+E[H|A_j] <= E[|F|+H|A_j] <=1.
```

Apply the same finite-dimensional Ornstein--Uhlenbeck smoothing to both
conditioned functions, with smoothing time tending to zero. Positivity,
conditional Jensen, and preservation of constants show that the resulting
finite smooth pair `(H_j,F_j)` still obeys

```
0<=H_j<=1, |F_j|+H_j<=1,
H_j even, F_j odd,
H_j->H and F_j->F in Gaussian L2.
```

Its finite coordinates are fixed before any matrix-size limit. Its odd
response need not equal `sign(UH_j)(1-H_j)`; the finite center theorem
requires feasibility, not that special optimization identity.

Define `K_j=U*F_j` and let `s_j` be its degree-51 nonlinear edge mass.
The adjoint and the finite coefficient projection are contractions, so

```
||K_j-K||_2<=||F_j-F||_2,
|s_j-s|<=||F_j-F||_2.
```

In particular `s_j>0` eventually. Apply the already unrestricted finite
center theorem separately to each such fixed pair. It gives

```
c_inf >= E H_j E_N|K_j+s_jN|.
```

These right sides converge to (2). For example, with `N` independent,
`||K+sN||_2^2=||K||_2^2+s^2<=||F||_2^2<=1`, since the two selected
Gaussian chaos components are orthogonal. Hence the absolute difference
of the right sides is at most

```
||H_j-H||_2 +(1+E|N|)||F_j-F||_2,
```

which tends to zero. This is an approximation of fixed Gaussian targets
after valid fixed-construction matrix inequalities; it does not run an
infinite algorithm at growing matrix order.

## 3. One common positive increment for every `J(H)>=.43`

The independently checked rational Gram/monomial argument in
`resumed_response_center_bestresponse_2026_09_06.md`, Section 4, gives

```
J(H)>=.43 => E G0 F>=.141,
s>97/(100000*2^25)>s0.                                   (3)
```

The exact replay is
`computations/resumed_response_degree51_certificate_2026_09_06.py`.
Its selected set has only the 25 odd degrees from 3 through 51, fixed
independently of the mask, its input dimension, and matrix size.

Put `mu=EH`, `p=EH^2`. The creation isometry and pointwise mask bounds give

```
J=E H K,
||K||_2^2<=||F||_2^2=1-2mu+p<=1-mu,
J^2<=p(1-2mu+p)<=mu(1-mu).
```

Therefore `J>=.43` forces `mu>6/25`: for `mu<=6/25`, monotonicity of
`mu(1-mu)` on `[0,1/2]` would give
`J^2<=114/625=.1824<.1849=.43^2`.

For `x>=0`, `g_s(x)` is decreasing and convex in `x`, and increasing
in `s>0`. Weighted Jensen and Cauchy--Schwarz imply

```
Phi(H,F)-J
 >=E H g_s(|K|)
 >=mu g_s(E[H|K|]/mu)
 >=mu g_s(sqrt((1-mu)/mu))
 >=(6/25) g_s0(sqrt(19/6))=delta0.                        (4)
```

The last expression increases with `mu` because both its positive
prefactor and its decreasing-argument factor increase. Together with
(2)--(3), this proves `c_inf>=J(H)+delta0` for every high-value mask.

For explicit positivity without subtracting almost equal Gaussian-tail
terms, set `a=x0/s0` and `u=1/(a+1)`. Then

```
g_s0(x0)=2s0 integral_a^infinity (z-a)phi(z) dz
        >=2s0 u^2 phi(a+2u)>0.
```

Thus the positive increment is not based on floating-point subtraction
or a conjectured numerical tail size.

## 4. Supremum and operator-limit order

The banked 21-anchor mask has value above `.4306`, so `C_mark>.43`.
For every sufficiently small `epsilon>0`, choose one fixed mask with
`J(H)>C_mark-epsilon` and `J(H)>=.43`. Equation (4) gives

```
c_inf>=C_mark-epsilon+delta0.
```

Letting `epsilon` decrease to zero proves (1). There is no compactness
assumption on the set of masks and no claim that a maximizing mask exists.

The finite theorem itself removes an operator cap safely: fix a deletion
fraction, obtain its fixed normalized operator bound on a principal
submatrix, apply the finite response theorem there, and use exact principal
monotonicity. Its Gaussian right side is independent of that fixed cap,
so the deleted fraction can then tend to zero. The above infinite-mask
approximation is applied to this already valid finite inequality.

## 5. Why the zero-strip counterexample and this escape coexist

The separately verified mask in
`resumed_response_high_value_zero_strip_2026_09_06.md` has `J>.4301875`
and `H=1` throughout `|UH|<=.001`. It falsifies the old proposed general
local-slack mechanism. The center response uses that available marked
amplitude rather than reducing it, so no zero-strip slack is needed.

Nor does this theorem iterate automatically. The updated center is
`H sign(S K+sZ)`, whose new spin depends on the old spin and unmarked
noise. It is not an independent own-spin mark, and the old creation
identity cannot be reapplied after renaming it. A closed enlarged
transport theory would be required to justify further rounds. In
particular (1) supplies neither a universal sharp bound nor an
all-order upper construction nor convergence/nonconvergence of the
original normalized minima.
