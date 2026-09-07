# Delocalized star-input Gram: a first frame, not the full hierarchy

2026-09-07. **Proved fixed-polynomial Gram lemma.**

Let B be real symmetric and hollow, with row square sums one,
||B||op<=L fixed, and b=max|B_ij|->0. Put Q=B^2 and P=B circ B.
Let h,g be fixed even polynomials. Constants may depend on h,g,L.

## 1. Gaussian inputs: exact identity

For independent standard Gaussian N, let G=BN and
V_i=N_i h(G_i), U_i=N_i g(G_i). Since G_i excludes N_i, the diagonal
covariance is E[h(Z)g(Z)]. For i!=j, integration by parts in N_i and
then N_j gives exactly

    E[V_i U_j]=B_ij^2 E[h'(G_i)g'(G_j)].              (1)

There are no additional terms: B_ii=B_jj=0, and i!=j removes the
derivative of the other explicit root factor. In the self-covariance
case define F_h(q)=E[h'(Z_1)h'(Z_2)] at correlation q. Because h' is
odd, its Hermite expansion gives

    F_h(q)=sum_(r>=1 odd) a_r^2 q^r.

Thus F_h(Q) is PSD and

    ||F_h(Q)||op <= E[h'(Z)^2] ||Q||op.

Every row of P has squared Euclidean norm at most b^2. Cauchy--Schwarz
therefore gives

    sum_j |P_ij F_h(Q)_ij| <= b ||F_h(Q)||op.

Since P has zero diagonal, this proves

    ||Cov(V)-E[h(Z)^2] I||op=O_h,L(b).               (2)

Cross covariances follow by polarization. Also Cov(N,V)=E[h(Z)] I
exactly. Hence, writing p=E h(Z), v=E h(Z)^2,

    Cov(BV)=v Q+O_op(b),       Cov(G,BV)=p Q.         (3)

This establishes every fixed finite family of star inputs, not merely
the quadratic Hermite star.

## 2. Rademachers: a weighted error rather than n times a local error

Let S be independent signs, G=BS, V_i=S_i h(G_i). At i!=j write
G_i=X+a S_j and G_j=Y+a S_i, where a=B_ij. The pair (X,Y) is independent
of the two endpoint spins. Exact averaging gives

    E[V_i V_j]=E[delta_a h(X) delta_a h(Y)],
    delta_a h(x)=(h(x+a)-h(x-a))/2.

For fixed polynomial h, delta_a h(x)=a h'(x)+a^3 r_a(x), where all
fixed moments of r_a(X) are uniformly bounded for |a|<=1. The pair
(X,Y) has variances 1-a^2 and covariance Q_ij. Its fixed polynomial
moments differ from the corresponding Gaussian moments by O_h(b^2):
every non-pairing even block has size at least four, whose coefficient
sum is at most b^2 times a fixed row-square sum. Restoring variance
one costs O_h(a^2). Consequently

    E[V_i V_j]=P_ij F_h(Q_ij)+O_h(P_ij b^2).          (4)

Since P has row sum one, this error has absolute row sum O_h(b^2),
not O(n b^2). The diagonal differs from v by O_h(b^2). Thus (2) also
holds on Rademachers. Cross covariances follow by polarization.
Cov(S,V) is exactly diagonal, with entries E h(G_i)=p+O_h(b^2), so
the second statement in (3) acquires only an O_op(b^2) error.

## 3. Associated star Gaussianization

The Gaussian response B[N h(BN)] is a finite sum of star chaoses.
Each nonconstant Hermite piece has fixed-root proper cuts O_h,L(b),
by `principle_invent_2026_09_07_dimension_free_tree_cuts.md`. Combining
these cuts with (3) and the fixed-chaos Gaussian theorem gives

    (G_i, [B(N h(G))]_i)
       -> (Z, p Z+sqrt(v-p^2) Z'),

uniformly in i; Z,Z' are independent standard Gaussians. The constant
Hermite component is already linear Gaussian.

Deleting repeated marked slots costs O_h(b) in row Hilbert norm.
The squarefree star has maximal influence O_h,L(b^2), so fixed-degree
coordinate replacement transfers its law to signs. The raw Rademacher
Hermite input differs from its squarefree star by O_h(b) in each row
L2: this is the same pairing-versus-large-block moment expansion as
in Section 2. Bounded root transport controls the averaged error.
Thus the displayed limiting pair holds for actual Rademacher
Y=B[S h(G)] in the root-averaged sense needed for bounded energy tests.

For an even H in [0,1], fixed-polynomial L2 approximation and fixed
||B||op give the same formula, first for smooth masks and then by
Gaussian L2 approximation. Degrees and masks are fixed before the
dimension limit.

## 4. Literal width witness and failed numerical target

For G=BS and Y=B[S H(G)], both vectors

    mu_+=S H(G)+(1-H(G)) sign(Y),
    mu_-=S H(G)-(1-H(G)) sign(Y)

lie in the cube. Their quadratic difference is exactly

    H_B(mu_+)-H_B(mu_-)=2 sum_i (1-H(G_i)) |Y_i|.

Cube rounding and the local limit therefore give

    liminf width(B)/n
      >= E[(1-H(Z)) |p Z+sqrt(v-p^2) Z'|].            (5)

Rank-one variance completion removes the fixed operator restriction
for this particular delocalized conclusion, with ordered limits.

An exploratory floating integral for H=1_(|Z|<=alpha) has maximum
approximately .3510653361304389 at alpha approximately .7108917249.
At alpha=.7 it is .3510412817089524. These are NOT interval certificates
or original signing improvements. The hoped value above .4 from this
one-stage mask was falsified; the recursive frame is necessary for
that proposed route.

The unproved next extension is the Gram for even responses depending
on recursively generated non-star fields. Identity (1) does not
automatically persist in that setting.
