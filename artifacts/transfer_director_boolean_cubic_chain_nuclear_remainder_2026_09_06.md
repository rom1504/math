# A cubic Boolean chain remainder at nuclear scale

2026-09-06. Proved elementary auxiliary theorem. This controls a precise
chain-rule error, NOT the remaining Stein brackets or the full rich-feedback
comparison. Low influence is an explicit hypothesis, not a property asserted
for arbitrary near-minimizing signings or their coherent returns.

## 1. Exact conditional cancellation

On the uniform Boolean cube write
`Delta_a U=(U(S_a=1)-U(S_a=-1))/2`. This derivative is independent of S_a.
For a vector query V write `V=V_0+S_a d`, with `d=Delta_a V`. For any C3
function f and any scalar Z, conditioning on all coordinates except a gives

```math
\mathbb E_a\left[\Delta_a Z\left\{\Delta_a f(V)
                  -\nabla f(V)\cdot d\right\}\right]
=-{\Delta_a Z\over4}\int_{-1}^1(1-t^2)
     D^3f(V_0+td)[d,d,d],dt.                         (1)
```

Indeed, with h(t)=f(V0+td), the expression in braces after averaging is
`[h(1)-h(-1)-h'(1)-h'(-1)]/2`. Two integrations by parts give (1).
Since the integral of `1-t^2` is 4/3,

```math
\left|\mathbb E\Delta_a Z\left\{\Delta_a f(V)
                     -\nabla f(V)\cdot\Delta_a V\right\}\right|
\le {\|D^3f\|_\infty\over3}
       \mathbb E|\Delta_a Z|\,\|\Delta_a V\|_2^3.      (2)
```

The operator norm of the trilinear third derivative is used. The factor
1/3 is exact: f(t)=t^3 with V0=0 attains equality. Keeping an unaveraged
quadratic Taylor remainder would discard the cancellation responsible for
the stronger influence power in (2).

## 2. A self-contained fourth-moment estimate

For a multilinear Boolean polynomial P of degree at most d,

```math
\|P\|_4^2\le\sum_S3^{|S|}\widehat P(S)^2
                   \le3^d\|P\|_2^2.                 (3)
```

For completeness, split P=P0+s P1 at one coordinate. Cauchy--Schwarz gives
`E P^4 <= a^2+6ab+b^2 <= (a+3b)^2`, with
`a=||P0||4^2`, `b=||P1||4^2`. Induction on the number of coordinates proves
the first inequality; the degree bound proves the second. Thus no imported
probabilistic approximation theorem is a dependency of this note.

For a q-dimensional polynomial vector V of degree at most d, putting
`I_a(V)=sum_l E(Delta_a V_l)^2`, Minkowski in L2 and (3) give

```math
\big(\mathbb E\|\Delta_a V\|_2^4\big)^{1/4}
\le3^{(d-1)/2}\sqrt{I_a(V)}.                         (4)
```

## 3. Nuclear-small chain defect

Consider n probe polynomials Z_j, each a pure homogeneous Boolean Fourier
chaos of a fixed positive degree K, and n vector queries V_i. Their base
cube may have any finite number of coordinates. Assume:

- each V_i has degree at most d, where d and K do not grow with n;
- `max_(j,a) E(Delta_a Z_j)^2 <= c/n`;
- `sup_i sum_a I_a(V_i) <= I_0`;
- `max_(i,a) I_a(V_i) <= kappa_n`;
- each f_i is C3 and `sup_i ||D3 f_i||_infty <= M`.

The vector dimension can be fixed or varying: the stated aggregate
influences and trilinear derivative norms already account for it. Define
the ACTUAL random brackets

```math
\Gamma_{j,i,l}=\sum_a\Delta_a Z_j\Delta_a V_{i,l},
\quad
R_{ij}=\mathbb E[Z_j f_i(V_i)]
 -{1\over K}\sum_l\mathbb E[
            \partial_l f_i(V_i)\Gamma_{j,i,l}].       (5)
```

Then, writing `||.||_*` for the matrix nuclear norm,

```math
\max_{i,j}|R_{ij}|
 \le C\sqrt{\kappa_n/n},\qquad
 {\|R\|_*\over n}\le C\sqrt{\kappa_n},
 \quad
C={M I_0\sqrt c\over3K}
           3^{(K+3d-4)/2}.                          (6)
```

Proof: Fourier orthogonality gives the exact integration-by-parts identity
`E[Z_j U]=(1/K) sum_a E[Delta_a Z_j Delta_a U]` for arbitrary cube U.
Apply (2), Hölder with four factors, (3) to Delta Z, and (4) to Delta V.
The resulting sum is bounded by

```math
{M\over3K}3^{(K+3d-4)/2}\sqrt{c/n}
          \sum_a I_a(V_i)^{3/2}
\le C\sqrt{\kappa_n/n}.
```

Finally `||R||_* <= sqrt(n)||R||_F <= n^(3/2) max|R_ij|` proves (6).
In particular vanishing aggregate coordinate influence gives nuclear o(n)
chain error. If kappa_n=O(1/n), the normalized nuclear error is O(n^-1/2).
Neither a finite moment catalog nor moment determinacy of the query law is
assumed.

## 4. What has and has not been removed

The exact brackets in (5) have NOT been replaced by covariance constants.
Controlling their fluctuation against the derivatives of f is a separate
obligation. The theorem therefore supplies one nuclear-small remainder in
a possible multiroot Stein proof; it is not a Gaussian comparison theorem.

This distinction matters for signing feedback: diffuse probe coordinates
alone do not imply kappa_n -> 0 for the coherent query list. Terms such as
QS may retain a coordinate of order one. In that case (6) only gives an
O(n) nuclear bound, and cannot authorize bounded-response closure.

Conversely, finite-degree low-influence queries may have non-Gaussian,
non-moment-determinate limiting laws. For example the multilinear cube
representation of `(n^(-1/2) sum_a S_a)^3` has degree three, bounded total
influence and maximum influence O(1/n), but converges to X^3 for Gaussian X.
The separate continuous-moment counterexample shows why polynomial density
in that OUTPUT cannot replace the direct bounded-smooth calculation here.
This theorem controls the chain remainder for that example without claiming
to turn the X^3 query into a Gaussian variable.
