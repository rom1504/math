# Hidden-partition count concentration for the ORIGINAL absolute cap

Date: 2026-09-07. Status: proved cleanup theorem; the endpoint drift inequality
needed for convergence remains unproved. Root and the adversary independently
noticed the sharper hypergeometric variance that improves the initial 3/4
pressure exponent to 2/3 below.

**Superseding improvement:** Section 7 proves the stronger bound
O_delta(beta^2 sqrt(N)), using vertex-star contraction. Its direct absolute-cap
consequence is O_delta(N^(-1/4)) in normalized expected error at beta=N^(1/4).
Sections 1--6 retain the valid earlier proof and isolate exactly what the
stronger monotonicity improves. The endpoint inequality remains unproved.

## 1. Definitions and conclusion

Fix 0<delta<=1/2 and N=m+n with delta N<=m,n. All auxiliary profiles use
the ambient denominator sqrt(N-1), including after vertex deletion. Define

```math
F(v;beta)=min_A log[(Z_A^+(v;beta)+Z_A^-(v;beta))/2],
Z_A^sigma=sum_x exp[sigma beta sum_{i<j} A_ij sqrt(v_ij)x_i x_j/sqrt(N-1)].
```

The minimum is over independent signs on the nonzero edges. The joint
orientation variable sigma is essential: this is absolute pressure, NOT
half-width pressure. Define also

```math
Q(v)=min_A max_x |sum_{i<j} A_ij sqrt(v_ij)x_i x_j|.
```

Select a uniform m-set S, reveal its membership in vertex order, and let K_k
be the number of positive labels among the first k vertices. Write u=N-k,
R=m-r, a=(N-1)/(m-1), b=(N-1)/(n-1). The posterior mean profile v(k,r)
has known-positive pairs a, known-negative pairs b, known-opposite pairs 0,
positive-known/unknown pairs aR/u, negative-known/unknown pairs b(u-R)/u,
and unknown pairs

```math
c(R)=[aR(R-1)+b(u-R)(u-R-1)]/[u(u-1)].
```

Absent pair types need no value. Every row sum is exactly N-1. Let r_0(k)
be an integer nearest mk/N in the feasible count interval; distance is at
most 1/2. Optimized pressures and caps depend on counts, not which labels
were revealed.

For all sufficiently large N depending only on delta, uniformly in k and
ALL beta>0,

```math
E |F(v(k,K_k);beta)-F(v(k,r_0(k));beta)|
    <= C_delta beta^2 N^(2/3).                         (1)
```

Consequently

```math
sup_k E |Q(v(k,K_k))-Q(v(k,r_0(k)))| / N^(3/2)
 <= C_delta beta N^(-1/3) + 2 log(2)/beta.            (2)
```

Taking beta=N^(1/6) makes (2) O_delta(N^(-1/6)). This is a genuine
original weighted absolute-cap concentration statement. It does not compare
the initial and terminal profiles.

## 2. Comparison lemmas, including the one-sided caveat

For each fixed A the absolute partition function is sum_x cosh(beta H_A).
It is nondecreasing in beta>=0. Therefore F(cv;beta)<=F(v;beta) for
0<=c<=1. This is GLOBAL radial monotonicity, not edgewise monotonicity.

At a minimizing signing, put lambda_e=beta sqrt(v_e)/sqrt(N-1) and
d_e=<sigma A_e x_i x_j> under the joint sigma/spin Gibbs law. Flipping
edge e changes the partition function by the ratio

```math
cosh(2lambda_e)-d_e sinh(2lambda_e) >= 1.
```

Thus d_e<=tanh(lambda_e), and on each smooth active branch

```math
partial_e F = beta d_e/[2 sqrt(N-1) sqrt(v_e)]
            <= beta^2/[2(N-1)].
```

The finite minimum is locally Lipschitz away from zero coordinates;
integrating there and taking continuous limits at zero proves, for w>=v,

```math
F(w)<=F(v)+ beta^2/[2(N-1)] sum_e(w_e-v_e).           (3)
```

Only this upper derivative bound is asserted. In particular F need NOT
increase coordinatewise, and no two-sided coordinate Lipschitz estimate is
used. Combine (3) with radial monotonicity: if w>=(1-epsilon)v and both
profiles have equal edge mass T, then

```math
F(w)<=F(v)+ beta^2 epsilon T/[2(N-1)].               (4)
```

A reverse bound requires and will receive reverse scaled domination.

Vertex insertion lemma: suppose the inserted row has variance sum <=N-1.
For any fixed signing, summing the new spin inserts a factor 2 cosh of
its row field. It follows that the optimized pressure increment is at
least log 2. For the upper bound fix an optimizing remainder and average
independent row signs. The averaged new absolute partition function is
exactly 2 times the old one times product_e cosh(lambda_e). Some row signing
is no larger than that average. Hence

```math
0 <= F_full-F_rest-log 2 <= sum_e log cosh(lambda_e)
                           <= beta^2/2.             (5)
```

Repeated deletion uses the same ambient denominator and can only decrease
remaining row sums. Finally every row-regular N-vertex profile obeys

```math
N log 2 <= F(v;beta) <= N log 2 + beta^2 N/4.        (6)
```

The lower bound follows from cosh>=1. The upper bound follows by averaging
all independent edge signs, yielding 2^N product_e cosh(lambda_e).

## 3. Consecutive-count comparison in the central band

Let u>=24/delta and suppose all counts between r and r+1 have remaining
positive fraction in [delta/2,1-delta/2]. Choose canonical known-label
arrangements for these two counts that differ at exactly one vertex w.
Delete w. Both remainder profiles have identical zero pattern, and their
total edge mass is exactly (N-2)(N-1)/2: both full profiles are row regular,
and deleting w removes N-1 edge mass.

Known-known weights agree. The known/unknown nonzero weights have relative
changes at most 2/(delta u) in either direction. For unknown pairs, a,b>=1
and u>=4 imply

```math
c(R) >= [u^2/2-u]/[u(u-1)] >= 1/3.
```

For N large depending on delta, a,b<=2/delta. Also

```math
|c(R-1)-c(R)|
 = 2|b(u-R)-a(R-1)|/[u(u-1)] <= 4/(delta u).
```

Thus every nonzero remainder coordinate changes relatively by at most
epsilon=12/(delta u)<=1/2, in BOTH directions. Apply (4) twice and (5)
for the deleted vertex. Writing F_{k,r}=F(v(k,r);beta), we obtain

```math
|F_{k,r+1}-F_{k,r}| <= 3 beta^2 N/(delta u)+beta^2/2
                     <= (4/delta) beta^2 N/u.       (7)
```

Deletion and reinsertion do not cost beta^2 instead of beta^2/2: the two
increments after subtracting log 2 each lie in the SAME interval
[0,beta^2/2], so their difference has magnitude at most beta^2/2.

## 4. Hypergeometric bulk and endpoint-window cleanup

Put p=m/N, mu=pk, D=|K_k-r_0(k)|. The exact variance is

```math
Var(K_k)=p(1-p) k u/(N-1) <= u/2.
```

Thus ED<=sqrt(u/2)+1/2<=2 sqrt u. On the event D<=delta u/4, every
intermediate count between K_k and r_0 has remaining positive fraction in
[delta/2,1-delta/2], provided u>=2/delta. Indeed the centered remaining
count differs from pu by at most 1/2. Consequently (7) telescopes to give
good-event expected error at most (8/delta) beta^2 N/sqrt u.

For u>=4/delta, D>delta u/4 implies |K_k-mu|>delta u/8. Chebyshev gives
probability at most 32/(delta^2 u). Applying (6) on this bad event adds
at most (8/delta^2) beta^2 N/u. No exponential concentration is required.
In particular the total bulk bound is

```math
E |F_{k,K_k}-F_{k,r_0}| <= (8/delta+8/delta^2) beta^2 N/sqrt u. (8)
```

For the endpoint window choose ANY terminal partition consistent with a
given revealed count. Its profile and the posterior profile agree on the
known-known submatrix. Delete the u unknown vertices from both. By (5)
both pressures belong to the same interval

```math
[F_known+u log2, F_known+u log2+u beta^2/2].
```

The optimized terminal pressure T_{N,m}(beta) is independent of the
particular terminal partition by vertex relabeling. Therefore for EVERY
feasible r,

```math
|F_{k,r}-T_{N,m}(beta)| <= u beta^2/2.
```

Applying this separately at K_k and r_0 gives the deterministic tail bound
|F_{k,K_k}-F_{k,r_0}|<=u beta^2. Split at u=N^(2/3). Equation (8) in the
bulk and this deletion bound in the tail prove (1), for example with
C_delta=1+8/delta+8/delta^2 and N^(2/3)>=24/delta.

## 5. Passage to absolute caps with growing beta

For any N-vertex profile, an orientation/spin attaining the absolute cap
contributes at least half its exponential to the absolute partition sum,
whereas that sum has total normalized counting mass 2^N. Minimizing gives

```math
beta Q(v)/sqrt(N-1)-log2 <= F(v;beta)
 <= beta Q(v)/sqrt(N-1)+N log2.                      (9)
```

Apply (9) in both directions to two profiles, then use (1) and
sqrt(N-1)<=sqrt N, N+1<=2N. This proves (2). Every error used was explicitly
proportional to beta^2, so beta=N^(1/6) is permitted without a hidden
fixed-temperature remainder.

## 6. The exact remaining inequality and its consequence

Write Psi_j(beta) for uniform-profile optimized ABSOLUTE pressure with
denominator sqrt(j-1), and E_j=M_j/sqrt(j-1). The endpoint normalization
above matches children exactly. Moreover

```math
T_{N,m}(beta) <= Psi_m(beta)+Psi_n(beta).             (10)
```

To prove this, select minimizing child signings, whose positive/negative
partition sums are (a_+,a_-) and (b_+,b_-). Reversing every sign in the
second child interchanges its two sums. The smaller of
a_+b_++a_-b_- and a_+b_-+a_-b_+ is at most
(a_++a_-)(b_++b_-)/2. Dividing by 2 for the parent orientation average
gives (10). Thus no width-to-cap bridge is necessary here.

One sufficient STILL UNPROVED deterministic inequality is

```math
F(v(0,0);beta) <= F(v(N,m);beta)+C_delta beta^2 N^(2/3) (11)
```

uniformly for comparable splits and, at least, beta=N^(1/6). One may seek
(11) along the centered-count three-species path v(k,r_0(k)), but (1)
does NOT prove its cumulative signed drift. In particular one must NOT sum
the O(beta^2 N^(2/3)) one-time cleanup error over all N reveal steps.
The two endpoints themselves have no count fluctuation, so concentration
alone gives no endpoint inequality whatsoever.

If (11) were proved, (9)--(10) would give directly

```math
E_N <= E_m+E_n+C_delta beta N^(2/3)+(N+1)log2/beta.
```

At beta=N^(1/6) this is

```math
E_N <= E_m+E_n+O_delta(N^(5/6)).                     (12)
```

Whether a suitable balanced-split almost-subadditivity principle plus the
already available local size regularity closes original convergence is
separate from the cleanup proved here. Neither (11) nor (12) is claimed.
The finite reveal drift obstructions and optimizer-switching signs in the
companion generator note remain relevant.

## 7. Stronger vertex-star contraction: square-root pressure cleanup

Here is the improvement not used in Sections 1--6. For a fixed signing A,
scale every coefficient incident to a single vertex i by t in [0,1].
Writing the Hamiltonian as H_rest+x_i L, summing x_i in the absolute
partition function gives

```math
sum_{x_i=+-1} cosh[beta(H_rest+t x_i L)/sqrt(N-1)]
 = 2 cosh[beta H_rest/sqrt(N-1)] cosh[beta t L/sqrt(N-1)].
```

Every summand is nondecreasing in t>=0. Minimizing over signings preserves
this monotonicity. Iterating over vertices proves

```math
F(D v D;beta) <= F(v;beta)                           (13)
```

for every diagonal matrix D with entries in [0,1]. Here D acts on VARIANCE
profiles, so a star step for D_ii=d_i scales its coefficients by sqrt(d_i).
This is NOT monotonicity under arbitrary individual edge contraction.

Now revisit the consecutive-count remainders of Section 3. All differing
coordinates are incident to one of the u unknown vertices; known-known
coordinates agree exactly. They still have equal edge mass. Choose D=1 on
known vertices and D=1-epsilon on unknown vertices, where
epsilon=12/(delta u)<=1/2. The relative-coordinate estimate already proved
implies D v D<=w: a known/unknown coordinate is multiplied by 1-epsilon,
and an unknown/unknown coordinate by (1-epsilon)^2<=1-epsilon. Similarly
D w D<=v. Moreover

```math
sum_e [v_e-(D v D)_e]
 = epsilon sum_KU v_e+(2epsilon-epsilon^2) sum_UU v_e
 <= epsilon sum_{i unknown} sum_{j neq i} v_ij
 <= epsilon u(N-1).
```

Use the one-sided derivative comparison (3), equal edge masses, and (13).
It gives F(w)-F(v)<=beta^2 epsilon u/2=6 beta^2/delta. Apply it in reverse,
then reinsert the single deleted known vertex using (5). The result is

```math
|F_{k,r+1}-F_{k,r}| <= (6/delta+1/2) beta^2
                     <= (7/delta) beta^2.           (14)
```

Thus on the same good event as Section 4, the expected error is at most
(14/delta) beta^2 sqrt(u). The identical Chebyshev and global-range estimate
give, for u>=24/delta,

```math
E |F_{k,K_k}-F_{k,r_0}|
 <= (14/delta) beta^2 sqrt(u)
    +(8/delta^2) beta^2 N/u.                        (15)
```

Split now at u=sqrt(N), using the SAME common-terminal-value deletion
argument for the tail. For sqrt(N)>=24/delta this proves, uniformly in k
and all beta>0,

```math
sup_k E |F_{k,K_k}-F_{k,r_0}|
 <= C'_delta beta^2 sqrt(N),
C'_delta=1+14/delta+8/delta^2.                      (16)
```

Using (9) gives the sharpened direct original-cap consequence

```math
sup_k E |Q(v(k,K_k))-Q(v(k,r_0(k)))|/N^(3/2)
 <= C'_delta beta/sqrt(N)+2log(2)/beta
 = O_delta(N^(-1/4))  at beta=N^(1/4).              (17)
```

Again the beta dependence is exact; no fixed-temperature asymptotic is
being substituted at growing beta.

For clarity, a correspondingly strong sufficient endpoint statement would
be the STILL UNPROVED inequality

```math
Psi_N(beta)<=T_{N,m}(beta)+C_delta beta^2 sqrt(N)
```

at beta=N^(1/4). Together with (9)--(10) it would imply

```math
E_N <= E_m+E_n+O_delta(N^(3/4)).
```

Equation (16) does not establish that endpoint statement. Its contribution
is removal of the random count fluctuations from one-time optimized
absolute-pressure/cap observables, with a quantitative vanishing normalized
error. The signed drift and any min-envelope switching loss are untouched.
