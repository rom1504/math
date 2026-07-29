# Finite-type conference fibres

## Status

This note analyzes the explicit obstruction
\[
A_k=C_k\otimes R+I_k\otimes D,\qquad
R=J_3-2I_3,\quad D=J_3-I_3,
\tag{1}
\]
where \(C_k\) is a symmetric conference matrix.  Three conclusions
are separated carefully:

1. the four-state reduction below is exact;
2. the finite-channel cavity functional is the natural fixed-depth
   limit, and its diagram proof is the conference proof with finite
   channel tensors inserted;
3. its numerical value for the conference two-step rule is computed,
   not yet interval-certified.

## 1. Exact four-state reduction

Choose representatives
\[
v_0=(1,1,1),\quad v_1=(-1,1,1),\quad
v_2=(1,-1,1),\quad v_3=(1,1,-1).
\]
Every \(x\in\{\pm1\}^{3k}\) is uniquely described fibrewise as
\[
x_u=\sigma_uv_{t_u},\qquad
\sigma_u\in\{\pm1\},\quad t_u\in\{0,1,2,3\}.
\]
The two finite kernels are
\[
K=(v_a^\top Rv_b)_{a,b=0}^3
=
\begin{pmatrix}
3&1&1&1\\
1&-5&3&3\\
1&3&-5&3\\
1&3&3&-5
\end{pmatrix},
\tag{2}
\]
and
\[
L_a=v_a^\top Dv_a=(6,-2,-2,-2)_a.
\tag{3}
\]
Therefore
\[
\boxed{
Q(A_k)=
\max_{\sigma,t}
\left|
\sum_{u,v}(C_k)_{uv}\sigma_u\sigma_vK_{t_u,t_v}
+\sum_uL_{t_u}
\right|.
}
\tag{4}
\]
This is an exact reduction from \(3k\) Boolean spins to a signed
four-state model on the \(k\)-vertex conference graph.  The last sum
is \(O(k)\), hence negligible on the \(k^{3/2}\) scale.

The equivalent three-layer form, with
\(x,y,z\in\{\pm1\}^k\), is
\[
x^\top(C_k\otimes R)x
=-x^\top C_kx-y^\top C_ky-z^\top C_kz
+2x^\top C_ky+2x^\top C_kz+2y^\top C_kz.
\tag{5}
\]

## 2. Numerical ground-state scale

Random-restart coordinate ascent on (1) gives:

| \(k\) | \(n=3k\) | best \(Q(A_k)/(n\sqrt{n-1})\) |
|---:|---:|---:|
| 42 | 126 | 1.06621 |
| 74 | 222 | 1.05628 |
| 110 | 330 | 1.03948 |
| 150 | 450 | 1.04013 |
| 194 | 582 | 1.03775 |

These are heuristic lower bounds, not exact maxima.  They show that
the fibre obstruction itself is far above
\(c_2=0.7833875\ldots\), even though a proof must not rely on the
heuristic search.

## 3. Finite-channel cavity law

Normalize the base and fibre matrices as
\[
U=C_k/\sqrt{k-1},\qquad T=R/\sqrt3.
\]
Then
\[
\frac{A_k}{\sqrt{3k-1}}=U\otimes T+O(k^{-1/2})
\tag{6}
\]
in the fixed-fibre cavity calculation.  At one base vertex let
\[
S\sim\operatorname{Unif}\{\pm1\}^3,\qquad
G\sim N(0,T^2)
\]
be independent.  For a threshold \(t\), define coordinatewise
\[
F_a=\operatorname{sign}(G_a+tS_a),
\]
\[
a_t=2\phi(t),\qquad b_t=2\Phi(t)-1,\qquad
R_0=F-a_tG-b_tS,
\tag{7}
\]
and put \(K_0=\mathbb E R_0R_0^\top\).  The matrix regression in
(7) is exact: \(\mathbb EFS^\top=b_tI\) and
\(\mathbb EFG^\top=a_tT^2\).

The fresh residual field is
\[
W\sim N(0,TK_0T),
\tag{8}
\]
independent of \((S,G,R_0)\).  For a Boolean response
\(Y=Y(S,G,R_0,W)\in\{\pm1\}^3\), define
\[
\mathcal A=\mathbb EYS^\top,\qquad
\mathcal H=\mathbb EYG^\top(T^2)^{-1},
\]
\[
\mathcal M=\mathbb EYR_0^\top,\qquad
\mathcal J=\mathbb EYW^\top(TK_0T)^{-1}.
\tag{9}
\]
The paired finite-channel energy is
\[
\boxed{
\mathcal E_T(Y)=\frac23\left[
\operatorname{tr}(T\mathcal A T\mathcal H^\top)
+\operatorname{tr}(T\mathcal M T\mathcal J^\top)
\right].
}
\tag{10}
\]
For a one-dimensional fibre \(T=1\), (10) is exactly
\[
2(\alpha\beta+\gamma\delta),
\]
the scalar two-step functional.

To derive (10), remove a base edge \(U_{uv}=\varepsilon\).  Its two
first-order cavity perturbations are
\[
G_u\mapsto G_u+\varepsilon TS_v,\qquad
W_u\mapsto W_u+\varepsilon TR_{0,v},
\tag{11}
\]
and symmetrically with \(u,v\) exchanged.  Gaussian integration by
parts on the receiving fields gives (10); summing uses
\(\sum_vU_{uv}^2=1\).  At fixed Hermite degree, the base-conference
contraction estimate is unchanged up to constants depending on the
fixed channel dimension three.  Dithering and Hermite truncation then
give the same order of limits as in the scalar conference theorem.

## 4. Evaluation of the inherited conference rule

At
\[
t=0.8414699114,
\]
the residual covariance has
\[
(K_0)_{aa}=0.3265071075\ldots,\qquad
(K_0)_{ab}=-0.0001741121\ldots\quad(a\ne b).
\tag{12}
\]
The off-diagonal entry is small but nonzero.  It is obtained exactly
from a shifted bivariate Gaussian sign integral at correlation
\(-1/3\):
\[
(K_0)_{ab}
=\mathbb E\!\left[
\operatorname{sign}(G_a+tS_a)
\operatorname{sign}(G_b+tS_b)\right]
-a_t^2(-1/3).
\tag{13}
\]

### 4.1 Exact variance rigidity

The small number in (12) has a useful general interpretation.  For
an arbitrary fixed sign fibre, write
\[
C=T^2,\qquad \operatorname{diag}C=1.
\]
Decompose the scalar residual in the hybrid Hermite basis as
\[
r(S,G)
=\sum_{\substack{\ell\ge3\\\ell\ {\rm odd}}}
u_\ell h_\ell(G)
+S\sum_{\substack{\ell\ge2\\\ell\ {\rm even}}}
v_\ell h_\ell(G),
\tag{14a}
\]
where the \(h_\ell\) are orthonormal Gaussian Hermite polynomials.
For two distinct fibre coordinates, the spins are independent, so
the entire \(S v_\ell\) branch disappears.  Consequently
\[
(K_0)_{ij}
=\kappa(C_{ij}),\qquad
\kappa(q)=
\sum_{\substack{\ell\ge3\\\ell\ {\rm odd}}}
u_\ell^2q^\ell
\quad(i\ne j).
\tag{14b}
\]
Therefore the average variance of the fresh field satisfies the exact
identity
\[
\boxed{
\frac1s\operatorname{tr}(TK_0T)
=s_t^2+
\frac1s\sum_{i\ne j}
\sum_{\substack{\ell\ge3\\\ell\ {\rm odd}}}
u_\ell^2C_{ij}^{\ell+1}
\ge s_t^2.
}
\tag{14c}
\]
All powers in the correction are even.  Thus the conference fibre
\(C=I\) minimizes the average residual-field variance.  At the
threshold used here,
\[
u_3
=\frac{2(t^2-1)\phi(t)}{\sqrt6}
=-0.06673981045\ldots,
\]
so quantitatively
\[
\boxed{
\frac1s\operatorname{tr}(TK_0T)-s_t^2
\ge
0.00445420230\ldots\,
\frac1s\sum_{i\ne j}C_{ij}^4.
}
\tag{14d}
\]
Equality holds only when \(C_{ij}=0\) for every \(i\ne j\), hence
\(T^2=I\).  This proves the desired matrix-level equality rigidity
for the first nonlinear residual channel.

Use the scalar conference response coordinatewise, with the
unnormalized residual divided by
\[
s_t=\sqrt{1-a_t^2-b_t^2}=0.5714080044\ldots:
\]
\[
Y_a=\operatorname{sign}\left(
0.5859761744S_a+0.6179560304G_a
+\frac{0.2396817825}{s_t}(R_0)_a
+\frac{0.4661704739}{s_t}W_a
\right).
\tag{14}
\]
Scrambled-Sobol evaluations of the six-dimensional Gaussian integral
in (10), through \(2^{21}\) points per scramble, give
\[
\boxed{
\mathcal E_T(Y)\approx0.78344,
}
\tag{15}
\]
about \(5\times10^{-5}\) above the scalar value
\(c_2=0.783387533648\ldots\).  This margin is not yet
interval-certified.  Direct finite-\(k\) simulations of the
Onsager-corrected rule increase toward this value:

| \(k\) | \(n=3k\) | mean normalized energy |
|---:|---:|---:|
| 602 | 1806 | 0.76652 |
| 1010 | 3030 | 0.77201 |
| 1202 | 3606 | 0.77290 |
| 1602 | 4806 | 0.77378 |

The slow approach is consistent with the fixed-depth cavity limit,
but again is not a proof of (15).

## 5. General finite-type target

For a fixed fibre size \(s\), a symmetric sign matrix \(R\), and a
zero-diagonal fibre filling \(D\), the same reduction uses
\[
T=R/\sqrt s,\qquad \operatorname{diag}(T^2)=1.
\]
It produces an \(s\)-channel version of (10).  A sufficient finite-type
theorem for the universal depth-two bound is
\[
\boxed{
\sup_{t,Y}\mathcal E_T(Y)\ge c_2
\quad\text{for every fixed symmetric sign fibre }T.
}
\tag{16}
\]
The three-fibre obstruction passes this test numerically, almost at
equality for the inherited scalar rule.  Proving (16), preferably by
showing that \(T^2=I\) minimizes the optimized functional, is the next
finite-dimensional lemma.  It would handle all bounded correlation
components, though a separate approximation theorem is still needed
to reduce arbitrary correlation structures to finite types.

Equation (14c) proves the variance half of that lemma.  The exact
remaining implication can now be stated as a robust variational
inequality:
\[
\boxed{
\sup_Y\mathcal E_T(Y)
\ge c_2+
c_0\left[
\frac1s\operatorname{tr}(TK_0T)-s_t^2
\right]
}
\tag{17}
\]
for some universal \(c_0>0\), or even with a lower-order error that is
uniform over fixed fibres.  The obstacle is that the coordinates of
\(W\) are correlated and can have unequal variances; average variance
monotonicity alone does not automatically imply (17).  This is now
the single finite-dimensional gap rather than an unspecified
``non-conference correlation'' problem.

### 5.1 Whitening removes the apparent noise loss

There is a second exact structural fact.  On the supports of the
covariance matrices, whiten the receiving fields by
\[
\widetilde G=(T^2)^{-1/2}G,\qquad
\widetilde R=K_0^{-1/2}R_0,\qquad
\widetilde W=(TK_0T)^{-1/2}W.
\]
The normalized direct-coupling operators are
\[
O_1=T(T^2)^{-1/2}=\operatorname{sgn}(T)
\tag{18}
\]
and
\[
O_2=K_0^{1/2}T(TK_0T)^{-1/2}.
\tag{19}
\]
Both are orthogonal:
\[
O_1O_1^\top=I,\qquad
O_2O_2^\top
=K_0^{1/2}T(TK_0T)^{-1}TK_0^{1/2}=I.
\tag{20}
\]
For singular fibres the same statement holds on the nonzero supports,
using Moore--Penrose inverses.

Thus the increased and unequal covariance of \(W\) is not an
information-theoretic noise loss.  After whitening, every canonical
correlation in both paired channels is still one; only an orthogonal
rotation in the finite channel space remains.  The finite-type lemma
can equivalently be phrased as:

> the Boolean response variational problem cannot fall below its
> scalar value when its two perfectly paired channels are subjected
> to the compatible orthogonal rotations (18)--(19).

This identifies the cube/rotation mismatch as the sole obstruction.
A generic Grothendieck rounding of the rotations loses a constant and
is too weak; a successful proof must exploit that the response may
depend jointly on all \(s\) channel coordinates.

### 5.2 An explicit conditional-score conjecture

The full matrix first variation makes the remaining claim concrete.
Let the scalar optimal moments for normalized residuals be
\[
(\alpha,\beta,\gamma,\delta)
=(0.5618430437,\ 0.5327670921,\
0.4238402492,\ 0.2179176762).
\]
Set \(C=T^2\), \(K=K_0\), and \(D=TK T\).  Insert the scalar target
moment matrices
\[
\mathcal A_0=\alpha I,\quad
\mathcal H_0=\beta I,\quad
\mathcal M_0=\gamma s_t I,\quad
\mathcal J_0=(\delta/s_t)I
\]
into the gradient of (10).  The resulting explicit score is
\[
\boxed{
L_T
=\beta CS+\alpha G
+\frac{\delta}{s_t}CR_0
+\gamma s_t\,C D^\dagger W,
\qquad
Y_T=\operatorname{sign}(L_T),
}
\tag{21}
\]
where \(D^\dagger\) is the Moore--Penrose inverse.  For \(C=I\), (21)
is exactly the scalar stationary response.

The finite-type problem is reduced to the following checkable
inequality:
\[
\boxed{
\mathcal E_T(Y_T)\ge c_2
\quad\text{for every symmetric sign fibre }T=R/\sqrt s.
}
\tag{22}
\]
This is stronger than merely asserting
\(\sup_Y\mathcal E_T(Y)\ge c_2\).

The identity behind (21) is useful for a future proof.  If
\((\mathcal A,\mathcal H,\mathcal M,\mathcal J)\) are the moments of a
response \(Y\), the first-variation score is
\[
L=
T\mathcal HTS+
T\mathcal ATC^\dagger G+
T\mathcal JTR_0+
T\mathcal MTD^\dagger W.
\tag{23}
\]
At a stationary point \(Y=\operatorname{sign}L\), homogeneity gives
\[
\boxed{
\mathcal E_T(Y)=\frac1s\mathbb E\|L\|_1.
}
\tag{24}
\]
Thus (22) is an \(\ell_1\)-support-function inequality for the
conditional score, not an opaque AMP assertion.

Computed evidence is strong:

- all switching/permutation classes of symmetric sign fibres for
  \(s\le4\) were tested;
- nonorthogonal \(s=2\) fibres give about \(1.0111\);
- the two nonorthogonal \(s=3\) classes give at least \(0.820\);
- every nonorthogonal \(s=4\) class tested is above the conference
  value, while the orthogonal classes return \(c_2\) within quadrature
  error;
- random fibres of sizes \(5,6,7\) also returned values above \(c_2\)
  after higher-precision rechecks.

These are computations, not an interval proof.  Attempts to deduce
(22) from average variance alone fail because some diagonal entries of
\(D\) can be slightly smaller than \(s_t^2\).  Anderson's inequality
does not directly apply: the non-Gaussian source pair
\((S,R_0)\) is rotated along with the Gaussian receiving fields.  A
black-box Borell/Grothendieck comparison again loses a fixed constant.
No counterexample to (22) was found.

## 6. Passage beyond finite fibres

A bounded-component decomposition is too restrictive.  A correlation
matrix can have a giant bounded-degree graph of constant correlations,
so deleting \(o(n)\) vertices need not reduce it to finite components.
The appropriate compact object is a rooted weighted correlation
graphing.

A sufficient quantitative approximation lemma for depth two is:

> **Mesoscopic correlation lemma.**  For every competing signing
> sequence and every \(\varepsilon>0\), either its Boolean quadratic
> witness already exceeds \(c_2-\varepsilon\), or, after removing
> \(o(n)\) vertices, its weighted row-correlation graph can be
> truncated to uniformly bounded degree so that the discarded
> correlations have
> \[
> \frac1n\sum_{\text{discarded }(i,j)}q_{ij}^4<\varepsilon.
> \tag{25}
> \]

Bounded-degree rooted neighborhoods then have subsequential local weak
limits, and every fixed-depth dithered AMP observable is continuous in
that topology.  Finite rooted types approximate the limiting graphing;
one does not need finite connected components.

The Schatten--Grothendieck bound proves only
\[
\sum_{i\ne j}q_{ij}^2=O(n^{5/4})
\]
at the \(c_2\) target.  It does not yet prove (25).  A high
correlation-degree vertex can carry enough fourth-moment mass to evade
bounded-degree truncation.  Converting that alternative into a direct
Boolean witness is the precise mesoscopic gap between the finite-type
conjecture (22) and a universal theorem.

## Reproducibility

- `_finite_fibre_amp.py`: exact construction, local search, and direct
  finite-\(k\) rounding.
- `_finite_type_population.py`: Sobol evaluation of (10).
