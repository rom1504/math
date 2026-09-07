# Independent audit: full-sign approximation of a low-rank orthogonal target

2026-09-07. Root's proposed rectangular operation passes reconstruction.
It controls the bilinear norm of the DIFFERENCE from the real target,
not merely the cap of a diagonally contracted surrogate.

Let H be an n-by-n sign Hadamard and T a real matrix satisfying
TT^T=T^T T=nI. Put Delta=H-T, L=||Delta||_*,
ell_i=[(Delta Delta^T)^(1/2)]_ii and
kappa_j=[(Delta^T Delta)^(1/2)]_jj. Singular-value factorization gives
|Delta_ij|<=sqrt(ell_i*kappa_j), and both leverage sums equal L.

Define D_ii=(1+ell_i)^(-1/2), E_jj=(1+kappa_j)^(-1/2), B=DTE.
Since H_ij is a sign,

    |T_ij|<=1+sqrt(ell_i*kappa_j)
       <=sqrt((1+ell_i)(1+kappa_j)).

Thus every B_ij lies in [-1,1]. Independently round it to a sign C_ij
with mean B_ij. This is always a full actual sign matrix.

## Uniform contraction error to T itself

For every Boolean x,y,
||(I-D)x||^2<=sum ell_i=L, and similarly for E. Using ||T||op=sqrt(n)
and ||D x||,||E y||<=sqrt(n),

    beta(T-DTE)<=2n*sqrt(L).

This step is necessary: a cap upper bound for the contracted matrix alone
would not establish approximation of the original weighted target.

## Sharper variance from exact orthogonality

The independent rounding variance is exactly

    V=sum_ij(1-B_ij^2)=n^2-||DTE||F^2.

The left contraction loses n*sum_i(1-D_ii^2), since each row of T has
squared norm n. The additional right contraction loses at most
n*sum_j(1-E_jj^2), since each column of DT has squared norm at most n.
Consequently V<=2nL. This argument retains the full weighted orthogonal
target; no entrywise absolute-sum surrogate or false independent leverage
assumption is used.

For u=(2n+2)log2, Bernstein for a fixed Boolean x,y gives deviation at
most sqrt(2Vu)+(4/3)u except with probability at most 2exp(-u). Union over
all 2^(2n) pairs has probability at most 1/2. Hence some rounded full sign
C obeys

    beta(C-T)<=2n sqrt(L)+2sqrt(nL*u)+(4/3)u.              (1)

If Delta has rank at most R, then ||Delta||op<=2sqrt(n), so
L<=2R sqrt(n). Equation (1) is O(n^(5/4)sqrt(R)+n), uniformly over target
orientations and leverage coherence. In particular R=o(sqrt(n)) gives
o(n^(3/2)) approximation in the ACTUAL full bilinear norm.

This does not assert that every relevant dense flatification difference
has low nuclear norm. It does validate the proposed low-rank orthogonal-
rotation localization mechanism whenever its rank bound is established.

## Direct simultaneous localization with one finite-rank correction

Let E_L,E_R be arbitrary subspaces of dimensions rL,rR, and d=rL+rR<=n.
Choose coordinate spaces S_L,S_R of dimension d. Dimension counting
permits V_R subset S_L intersect E_L-perp of dimension rR and V_L subset
S_R intersect E_R-perp of dimension rL. The domain E_R direct-sum V_L and
range V_R direct-sum E_L are orthogonal decompositions. Choose an isometry
between them sending E_R to V_R and V_L to E_L.

Starting from O=H/sqrt(n), realize this isometry by a left orthogonal
correction U which fixes the orthogonal complement of the span of the
old and desired images. That span has dimension at most 2d, so
rank(U-I)<=2d. The target T=sqrt(n)*UO satisfies

    T E_R subset S_L, T^T E_L subset S_R,
    P_L T P_R=0, rank(T-H)<=2d.

For arbitrary Boolean x,y, decompose using the orthogonal projections
P_L,P_R. The cross-projected term vanishes exactly. The two remaining
terms with one projected factor are bounded by n sqrt(d), using the
coordinate support on the opposite side and the ORIGINAL Boolean vector.
No assertion that an arbitrary orthogonal projection preserves the cube
is made. Thus

    |x^T T y|<=sqrt(n)*sqrt(n-||P_Lx||^2)*sqrt(n-||P_Ry||^2)
                   +2n sqrt(d).

Combining with (1) gives an actual full sign bridge with the same residual
bound plus O(n^(5/4)sqrt(d)+n). For d=o(sqrt(n)), the error vanishes after
n^(3/2) normalization. This accommodates arbitrary subspaces, not only
sign-pattern partitions or mutually orthogonal Boolean centers.

The all-order version, including a primary-source-checked polynomial
Hadamard-order gap and the sharper zero-embedding argument, is proved in
`flatify_adversary_2026_09_07_power_saving_all_order_localization.md`.
Its error remains O(n^(5/4)sqrt(d)) for d>=1; it does not pay a needless
square-root deletion error from Boolean padding.
