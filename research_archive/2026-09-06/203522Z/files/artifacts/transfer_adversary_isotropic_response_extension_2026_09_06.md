# Isotropic absolute response is stable under sublinear vertex extensions

Date: 2026-09-06. Root-proposed theorem, independently reconstructed and
audited here: PASS. No independence assumption on an optimizing law is
needed for the upper bound. This is a statement about actual Boolean
laws, not replacement of their value by a spectral relaxation.

For a hollow symmetric signing A of order n define

```math
I(A)=\max_{\mu:\ E_\mu xx^T=I_n} E_\mu |H_A(x)|,
\qquad H_A(x)=x^T A x/2.
```

The feasible laws form a nonempty compact finite-dimensional polytope,
so the maximum exists. In the earlier finite isotropy screen this is
exactly `Q(A)` minus minimum isotropic mean slack.

## Finite extension theorem

Let n,r>=1 and let

```math
B=\begin{pmatrix}A&C\\C^T&D\end{pmatrix}
```

be ANY hollow signing extension by r vertices. Then

```math
\boxed{0\le I(B)-I(A)
\le r\sqrt n+\frac r2\sqrt{r-1}.}                      (1)
```

For the upper bound take any jointly isotropic law on `(x,y)`.
Its x and y marginals are isotropic. Since `||y||_2=sqrt(r)`,
Cauchy--Schwarz and `E xx^T=I_n` give

```math
E|x^TCy|\le\sqrt r\,E\|C^Tx\|_2
\le\sqrt r\,(\operatorname{tr}CC^T)^{1/2}=r\sqrt n.
```

This uses only the x marginal and the deterministic length of y;
x and y need NOT be independent. For the new principal block, spectral
absolute value gives `|y^TDy|<=y^T|D|y`. Thus

```math
E|H_D(y)|\le\tfrac12\operatorname{tr}|D|
\le\tfrac12\sqrt r\|D\|_F
=\tfrac r2\sqrt{r-1}.
```

The first term satisfies `E|H_A(x)|<=I(A)`. Apply the triangle
inequality to `H_B=H_A+x^TCy+H_D` and maximize over joint laws.

For the lower bound use a maximizing isotropic law for x and add r
independent unbiased signs y, independent of x. The full raw second
moment is identity, even if x itself has nonzero mean. Conditional on
x, the bridge and new-block energies have mean zero. Jensen yields
`E[|H_B(x,y)| | x]>=|H_A(x)|`, proving monotonicity.

## Uniform near-order continuity and the response-gap consequence

For every signing A, the same spectral-absolute-value argument gives

```math
0\le I(A)\le\tfrac12\operatorname{tr}|A|
\le\tfrac n2\sqrt{n-1}.                               (2)
```

This is an upper bound for I, NOT for Q. From (1)--(2), if
`t=r/n`, then

```math
\left|\frac{I(B)}{(n+r)^{3/2}}-\frac{I(A)}{n^{3/2}}\right|
\le t+\tfrac12t^{3/2}.                                (3)
```

For the upward difference, use the increment in (1). For the downward
difference, monotonicity and (2) bound it by
`(1/2)[1-(1+t)^(-3/2)]<=(3/4)t`. Hence every extension by `r=o(n)`
has `I(B)-I(A)=o(n^(3/2))`, uniformly in all signings and bridges.
No bound on the original cap Q(A) is required.

Principal-cap monotonicity also gives the exact consequence

```math
Q(B)-I(B)\ge Q(A)-I(A)-r\sqrt n-\tfrac r2\sqrt{r-1}.   (4)
```

Consequently, IF a given parent sequence has a macroscopic isotropic
response gap, no sublinear-vertex extension can erase that gap at leading
order. Making its oriented extrema balanced would not change (4).
This is conditional: no macroscopic response gap for actual exact
minimizers is asserted here, and neither (1) nor (4) settles convergence
of the minimum normalized Boolean cap.

## Exact finite replay

`computations/transfer_adversary_isotropic_response_extension_2026_09_06.py`
checks every root-normalized signing of total order two through five,
and every nonempty proper choice of old vertices: 2046 extension
comparisons. Its 76 distinct response LPs have exact rational primal
isotropic laws and exact pointwise quadratic dual majorants with equal
objectives. Monotonicity and the radical upper bound are checked using
rational arithmetic, not a floating-point tolerance. It also verifies
(2) in squared rational form on every matrix.

Full certificates are preserved in
`computations/results/transfer_adversary_isotropic_response_extension_2026_09_06.json`.
The checker reuses rational reconstruction helpers from the earlier
`transfer_adversary_minimizer_isotropy_2026_09_06.py`; it has no external
download dependency. Final run: PASS all 2046 comparisons and 76 exact
primal/dual values. The finite replay is a regression check; the proof
of (1)--(4) applies at all orders.
