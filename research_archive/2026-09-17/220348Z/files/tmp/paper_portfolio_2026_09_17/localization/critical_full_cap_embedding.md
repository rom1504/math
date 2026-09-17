# Critical block size: unrestricted full-sign cap embedding

2026-09-17. New deduction by the localization track from the director's
covariance-matched row law. This upgrades the restricted-query sharpness
example to an unrestricted cap, but NOT to a bounded-normalized-cap family.

Let B be any n by m full sign matrix, m=3n, n>=3. Form the full hollow
sign matrix W with left child -offdiag(J_n), right child +offdiag(J_m),
and bridge B. Put N=n+m=4n, r=B 1, and beta=binom(m,2).

For fixed x, replacing a nonconstant right word y by its majority sign
gains 2k(m-k) in right-child energy, where 1<=k<=m/2 is its minority
size. The bridge can lose at most 2nk. Since m-k>=m/2>n, the total
positive energy strictly increases. Thus the positive maximum has a
constant right word, which simultaneous reversal reduces to y=1.

Also, for every x,y,

    -H_W(x,y) <= binom(n,2)+m/2+nm = (7/2)n^2+n.

The positive maximum is at least beta=(9/2)n^2-(3/2)n: average H_W(x,1)
over uniform x. For n>=3, beta exceeds the preceding negative bound.
Consequently Q(W) is always the positive maximum and EXACTLY

    Q(W) = beta+n/2+max_x [r dot x-(sum_i x_i)^2/2].          (A)

Suppose now B has independent globally symmetric rows with the same row
law. A fair randomized tie-break makes sign(r_i) independent fair signs.
Evaluating the maximum at those signs, and bounding it above by sum|r_i|,
gives the two deterministic expectation inequalities

    beta+n E|r_1| <= E Q(W) <= beta+n E|r_1|+n/2.           (B)

Apply this to the director's tilted row law in dimension m and to iid
fair sign rows. Both laws have EXACT mean zero and covariance I_m,
and a uniform Euclidean subGaussian constant. Their row-sum absolute
expectations differ by (eta+o(1))*sqrt(m), where eta>0 (the limiting
integral evaluates to about 0.02230342349). The row proof alone gives
eta>26/11025, using only elementary Gaussian tail inequalities.

Subtracting (B),

    E Q(W_iid)-E Q(W_tilted)
       = (eta+o(1))*n*sqrt(m)+O(n)
       = (eta*sqrt(3)/8+o(1))*N^(3/2).                     (C)

All entries of BOTH compared matrices are actual signs. All randomness
is in n independent blocks of size m=(3/4)N; deterministic child terms
are common arbitrary offsets in the block theorem. Thus an all-offset
unrestricted-full-sign cap comparison cannot have o(N^(3/2)) error at
linear block size using only matching block covariances and a uniform
subGaussian constant. Comparing both laws to a common Gaussian law
would also contradict (C), so an explicit Gaussian pinning argument is
unnecessary.

Important scope: Q(W) itself is Theta(N^2), from the deliberately strong
right ferromagnet. This example proves sharpness of the GENERAL
all-offset theorem. It does not refute a stronger theorem restricted
to actual minimizing children or Q(W)=O(N^(3/2)).

The row-law exact rational diagnostic is
computations/paper_localization_2026_09_17_critical_block.py and its
output critical_block_audit.json in this same directory. No numerical
calculation is a premise of (A)--(C).
