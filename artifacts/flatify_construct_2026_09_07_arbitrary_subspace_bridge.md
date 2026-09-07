# Actual-sign Hadamard bridge adapted to arbitrary prescribed subspaces

Root proposed the operation below. This is a positive construction theorem,
not a sufficient condition for the original flatification target.

## Theorem

Let an order n sign Hadamard H exist, and let X,Y be arbitrary real
subspaces of R^n of dimensions rx,ry with rx+ry<=n. There is an actual
n-by-n full-sign bridge C such that, uniformly for all Boolean x,y,

 |x^T C y|<=sqrt n ||(I-P_X)x|| ||(I-P_Y)y||
             +O(n^(5/4)sqrt(rx+ry)+n).                (1)

The constants are absolute. Thus rx+ry<=O(n^(1/2-2delta)) gives an
O(n^(3/2-delta)) remainder, for any fixed0<delta<=1/4. The zero-rank case
can use H itself with zero remainder. No coordinate-pattern or finite
partition assumption on the prescribed subspaces is made.

## 1. An orthogonal modification with localized exceptional terms

Normalize H0=H/sqrt n. Choose ANY coordinate subspace E of dimension ry.
There is an orthogonal U mapping H0Y onto E with rank(U-I)<=2ry: extend
the isometry within the span of these two spaces and use identity on its
orthogonal complement. Put W1=U H0.

Let Z=P_(Y-perp) W1^T X, of dimension z<=rx. Choose a coordinate set S
of size rx+ry. Since dim(Y-perp intersect R^S)>=rx, choose a z-dimensional
subspace F there. An orthogonal V fixing Y pointwise can map Z onto F,
with rank(V-I)<=2z<=2rx. Put

 W=W1 V^T.

Then W is orthogonal, WY=E, and
P_(Y-perp)W^T X is contained in F. Also

 rank(W-H0)<=2(rx+ry),
 ||sqrt n W-H||_*<=4(rx+ry)sqrt n.                   (2)

The latter uses ||W-H0||op<=2. All choices can be made explicitly by
orthonormal basis extension and nullspace computations.

There is a simplification of the root's initial proof: E need not be
chosen to have small overlap with X. Set x_perp=(I-P_X)x and
y_perp=(I-P_Y)y. The exact expansion is

 x^T W y=x_perp^T W y_perp+x^T W P_Y y
                              +(P_Xx)^T W y_perp.

The second term has magnitude at most sqrt(ry n), because WP_Yy lies in
the coordinate space E. In the third term only
f=P_(Y-perp)W^T P_Xx contributes. This f belongs to F, so f is supported
on S and orthogonal to Y. Therefore

 |<f,y_perp>|=|<f,y>|<=||f|| sqrt(|S|)<=sqrt(n(rx+ry)).

Multiplying by sqrt n proves the deterministic real-matrix estimate

 |x^T(sqrt n W)y|
 <=sqrt n||x_perp||||y_perp||
    +n[sqrt(ry)+sqrt(rx+ry)].                         (3)

In particular the deterministic exceptional error is O(n sqrt(rx+ry)),
slightly better than the initially proposed O(nr).

### One-shot alternative with exact feature-feature cancellation

The independent/adversarial agents supplied a cleaner simultaneous
construction. Choose coordinate sets of size rx+ry and subspaces
V_R contained in X-perp of dimension ry, and V_L contained in Y-perp
of dimension rx, each supported on its respective coordinate set.
These exist by the dimension bound for intersecting a coordinate space
with an orthogonal complement.

Prescribe an isometry from the orthogonal sum H0Y plus H0V_L onto the
orthogonal sum V_R plus X, mapping the first summand onto V_R and the
second onto X. Extend it to an orthogonal U, equal to identity outside
the span of its domain and target. Then rank(U-I)<=2(rx+ry). With W=UH0,

 WY=V_R, W^T X=V_L, and P_X W P_Y=0 EXACTLY.

The two exceptional terms in the orthogonal decomposition can each be
paired with the full Boolean vector, because V_R is orthogonal to X and
V_L to Y. Their coordinate supports give deterministic error at most
2n sqrt(rx+ry) after scaling by sqrt n. Equations(2) and all subsequent
rounding bounds remain unchanged. This avoids two sequential rotations
and has the additional exact feature-feature cancellation. The earlier
two-rotation proof is retained as an independently valid alternative.

## 2. Rectangular target contraction and uniform discrepancy rounding

The following rounding lemma is stated separately to make its stronger
discrepancy conclusion clear. Let A be any n-by-n full-sign matrix,
T=A-Delta with ||T||op<=sqrt n, and L=||Delta||_*. From a singular value
decomposition of Delta, define

 ell_i=(sqrt(Delta Delta^T))_ii,
 r_j=(sqrt(Delta^T Delta))_jj.

Then sum ell_i=sum r_j=L and |Delta_ij|<=sqrt(ell_i r_j). Put
s_i=(1+ell_i)^(-1/2), t_j=(1+r_j)^(-1/2), and M_ij=s_i T_ij t_j.
The inequality

 (1+ell_i)(1+r_j)>=(1+sqrt(ell_i r_j))²

shows that every M_ij lies in[-1,1]. Independently round each entry to a
sign C_ij with expectation M_ij. The total variance is at most4nL:

 sum(1-M_ij²)
 <=sum[1-s_i²t_j²]+2sum|Delta_ij|
 <=2nL+2nL.

Here sum|Delta_ij|<=nL follows from the rank-one magnitude envelope and
Cauchy--Schwarz. For a=(2n+2)log2, Bernstein's inequality and a union bound
over all2^(2n) Boolean pairs show, with probability at least1/2,

 beta(C-M)<=sqrt(8nLa)+4a/3.                          (4)

The important additional observation is that contraction did not move the
target far in bilinear cube norm. Since
sum(1-s_i)²<=L and sum(1-t_j)²<=L,

 beta(M-T)<=2n sqrt L.                              (5)

Indeed expand M-T=(diag(s)-I)T diag(t)+T(diag(t)-I), and use
||T||op<=sqrt n, ||x||=||y||=sqrt n, and the displayed diagonal-square
bounds. Thus the rounding produces the UNIFORM discrepancy estimate

 beta(C-T)<=2n sqrt L+sqrt(8nLa)+4a/3
           =O(n sqrt L+n).                          (6)

This is stronger than merely preserving the cap of T. It relies on the
operator norm bound on T; that hypothesis must not be dropped.

## 3. Completing the actual construction and its scope

Apply the lemma with A=H and T=sqrt n W. Equation(2) gives
L<=4(rx+ry)sqrt n. Combining(3) and(6) proves(1). The sampling procedure
has a positive, explicitly bounded success probability; the resulting
matrix has exactly sign entries, not weighted or approximately flat ones.

This supplies a genuinely broader bridge operation than localization of
a bounded number of coordinate-pattern classes. It works for arbitrary
prescribed subspaces with total dimension o(sqrt n), with explicit power
savings when their dimension has a power margin below sqrt n.

It does not by itself prove the requested favorable flatification: no
theorem here chooses X,Y so that the old optimal-child energies plus the
residual term in(1) fit the target cap. Both child reversals remain part of
the parent maximum. Old-edge changes are authorized in the main problem,
but this theorem only constructs a bridge for fixed chosen subspaces.

## 4. An explicit child-energy threshold certificate and its limitation

For equal-order children write p=||P_Xx||²/n and q=||P_Yy||²/n. Suppose
each child has cap at most c n^(3/2), and EVERY spin with feature fraction
at most tau has absolute energy at most a n^(3/2). Equation(1) then gives
the following fully paid parent upper coefficient, apart from its rounding
remainder:

 max{2a+1, c+a+sqrt(1-tau), 2c+1-tau}.                (7)

These are respectively the two-low, one-high, and two-high feature cases.
Thus sufficient numerical requirements for the equal-child target are

 2a+1<=2sqrt2 c,
 c+a+sqrt(1-tau)<=2sqrt2 c,
 2c+1-tau<=2sqrt2 c,

with a fixed positive margin if one wants to absorb the rounding error.
In particular a<=sqrt2 c-1/2 and
tau>=1-2(sqrt2-1)c are necessary for this particular certificate. No such
actual-optimizer energy threshold is assumed to hold.

### The fixed-rank threshold is contradicted by the width theorem

For ANY prescribed subspace X of fixed rank r, including subspaces that
vary with n, and any fixed tau>0, the slice ||P_Xx||²/n<=tau has energy
oscillation at least(2L0-o(1))n^(3/2), where
L0=2/(3sqrt(3pi)). This follows from the finite-partition balanced-width
theorem by the following explicit approximation.

Let U be an n-by-r orthonormal basis matrix. Classify its rows by the
scaled vectors sqrt n U_i. Rows of norm exceeding T number at most
rn/T². Put these in one exceptional class. Quantize each coordinate of
the remaining rows to mesh epsilon inside[-T,T], producing a fixed finite
partition (r,T,epsilon fixed). Spins balanced in every class annihilate
the quantized basis, up to O(1/sqrt n) from odd class sizes. The low-row
quantization error has Frobenius norm at most sqrt r epsilon; the high-row
contribution has norm at most the square root of the number of high rows,
since coordinate restriction of U has operator norm at most1. Hence

 ||U^T x||/sqrt n<=sqrt r epsilon+sqrt r/T+o(1).

Choose T large and epsilon small to make this below sqrt(tau). The width
theorem supplies two such nearly feature-orthogonal spins whose energy
difference is at least(2L0-o(1))n^(3/2). Therefore any threshold a in(7)
must satisfy a>=L0-o(1). Already the first term of(7) is at least
1+2L0-o(1)=1.4343133439...-o(1), which exceeds2sqrt2 c throughout c<=1/2.

Thus this elementary threshold application cannot work with fixed-rank
features, even though the actual bridge construction itself is valid.
This is not a no-go theorem for joint bridge/child control, for polynomially
growing feature rank, or for authorized global changes to child edges.

## 5. Stronger Gaussian audit: residual energy survives the whole rank range

Root subsequently supplied a stronger argument, independently checked
here. It extends the scope obstruction to arbitrary rank r=o(sqrt n),
the entire vanishing-error range of the bridge construction. It does not
invalidate that construction or establish a large actual parent cap.

Let A be any low-cap full-sign matrix, X any rank-r subspace, and split
vertices deterministically into I,J of sizes s=n/3+O(1), t=2n/3+O(1).
For independent standard Gaussians g indexed by I, set

 x_I=sign(g), x_J=sign(A_JI g),
 x_plus=(x_I,x_J), x_minus=(x_I,-x_J).

Their cross energy L=x_I^T A_IJ x_J has the EXACT expectation

 E L=st(2/pi)arcsin(1/sqrt s)
     =[4/(3pi sqrt3)+o(1)]n^(3/2).                   (8)

This uses the Gaussian sign-correlation identity coordinate by coordinate;
each correlation is A_ij/sqrt s. There is no central-limit approximation.
Moreover |L|<=Q(A), because the two full energies are T+L and T-L.

The underlying standardized Gaussian vector has Gram matrix R=B B^T,
where B stacks I_s above A_JI/sqrt s. Thus diag R=1 and
||R||op<=1+||A||op²/s. Its sign covariance is K=(2/pi)arcsin(R), entrywise.
The arcsine series has positive coefficients whose normalized sum is1.
Every Schur power of the correlation matrix R has operator norm at most
||R||op: Schur multiplication by a positive semidefinite diagonal-one
matrix is positive, unital, and contractive on self-adjoint operator norm.
It follows that ||K||op<=||R||op.

For clarity, the sharper cap-to-operator bound needed here is
||A||op²<=8Q(A). Let beta be the real bilinear cube norm. The complex
infinity-to-one norm is at most2beta: rotate the total phase and split
the real part into two real bilinear terms. Interpolating this bound with
the complex1-to-infinity norm1 gives ||A||op²<=2beta; polarization gives
beta<=4Q. Consequently ||K||op=O(sqrt n).

The sign-reversed witness has covariance D_J K D_J of the same norm, so

 E[||P_X x_plus||²+||P_X x_minus||²]=O(r sqrt n)=o(n).

Markov's inequality makes both feature fractions o(1) with probability
1-o(1). Removing the bad event from (8) costs only o(n^(3/2)), since
|L|<=Q(A)=O(n^(3/2)). Some good outcome therefore has L at least the
right side of(8), up to o(n^(3/2)), while BOTH witnesses have feature
fraction o(1). Their energy difference is2L.

Thus the low-feature slice has half-width at least

 G0=4/(3pi sqrt3)=0.245035064631907...,

uniformly for every prescribed rank r=o(sqrt n). Any threshold a in(7)
must be at least G0-o(1), making its first term at least
1+2G0-o(1)=1.49007012926...-o(1). More generally, separately optimized
SIGNED child profiles have the same obstruction by adding their interval
widths, as in balanced_overlap_obstruction.md Section3. Their critical
coefficient is(.5+G0)/sqrt2=0.5268193464..., above the entire current regime.

Accordingly, the positive bridge theorem is not a completed flatification:
one must exploit joint child/bridge dependence beyond just the residual
norms, or perform additional authorized old-edge changes. Merely increasing
the prescribed feature rank while keeping it o(sqrt n) cannot repair this
particular scalar-profile certificate.
