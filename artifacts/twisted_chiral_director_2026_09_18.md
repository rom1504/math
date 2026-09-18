# Director derivations: twisted chiral doubles

Status: work in progress, 2026-09-18. Original normalization throughout:
Q(A)=max_x|sum_(i<j)a_ij x_i x_j|. Numerical, exact finite, and uniform
statements are explicitly separated. No convergence result is claimed.

## 1. Two exact forms of the interaction

For C=B+diag(d), direct expansion and interchange of x,y give

```math
Q\!\begin{pmatrix}A&C\\C&-A\end{pmatrix}
=\max_{x,y}\left\{|H_A(x)-H_A(y)|+|x^TCy|\right\}.
```

Indeed interchanging x,y reverses the first summand and fixes the second;
max(|u+v|,|-u+v|)=|u|+|v|. This formula is valid for arbitrary real
symmetric A,C with A hollow. It is not cancellation between two absolute
maxima. The two summands must be evaluated at the SAME pair x,y.

Putting y=t*x and I={t=1} gives the independently audited cut identity in
twisted_chiral_adversary_2026_09_18.md. Chiral symmetry centers the energy
range but does not bound its size.

## 2. A smaller, necessary invariant with exact diagonal optimization

Write beta(A)=max_y sum_i|(Ay)_i| and
gamma(A)=min_(d in {+-1}^n) beta(A+diag(d)). Every allowed twisted double
satisfies Q(D)>=gamma(A)>=beta(A)-n. The first inequality follows from
Section1 and conjugating B back to A: signed permutation conjugacy maps
diagonal signs bijectively to diagonal signs. It is an invariant necessary
bound, NOT a sufficient parent cap bound.

For a hollow sign matrix let h_i(y)=y_i(Ay)_i, l(y)=sum_i|h_i(y)|,
z(y)=#{i:h_i(y)=0}, sigma_i(y)=sign(h_i(y)), with sign(0)=0. Integrality gives

```math
\beta(A+\operatorname{diag}d)
=\max_y\{l(y)+z(y)+\sigma(y)\cdot d\}.
```

This is exact because |h+d|=|h|+1_(h=0)+sign(h)d for every integer h and
d=+-1. Thus diagonal optimization is a finite linear discrepancy problem
in n bits, with 2^(n-1) distinct rows, rather than a 2^(2n)-pair search.
It does not resolve the larger switch/permutation optimization.

The director's integer program and the adversarial program independently
agree on the following data for the stored representatives:

| n | child Q | beta | gamma |
|---|---:|---:|---:|
|3|3|6|5|
|4|4|8|8|
|5|4|8|11|
|6|5|12|14|
|7, class0|9|18|17|
|7, classes1,2|9|18|19|
|8, class0|10|24|20|
|8, class1|10|24|22|
|9, stored representative|12|28|27|
|10, stored representative|13|40|30|

These are exact finite necessary bounds. In particular beta(A)/Q(A)>2sqrt2
at order10 does not itself obstruct the finite target: diagonal optimization
can remove ten units. Conversely gamma does not distinguish the successful
order7 class1 from class2, whose complete twist family has a larger minimum.

Code/output: computations/twisted_chiral_director_invariants_2026_09_18.py
and computations/results/twisted_chiral_director_invariants_2026_09_18.json.

## 3. Inverse reconstruction from previously stored parent witnesses

The director independently searches for a signed permutation S with
S^2=-I and S^T D S=-D. Rooted signed triangle graphs turn the switching
isomorphism condition into ordinary finite graph isomorphism. An exact
check of every recovered map precedes use. For each tested transversal
of the paired coordinates, write D in canonical blocks A,C,-A and test
whether C-diag(C) is in A's signed-permutation orbit.

This reconstructs explicit twists for stored parent orders4,6,8,10,12,16,
with exact full-cube Walsh-transform caps4,5,10,13,18,30 respectively.
The first recovered seeds at orders10,12,16 have caps6,7,14, so they are
not optimal children. This is a property of those recovered presentations,
not evidence that optimal children cannot attain the same parent cap.
The forward search independently finds optimal-child presentations.

Code: computations/twisted_chiral_director_inverse_2026_09_18.py.
Initial raw output is preserved at tmp/twisted_chiral_2026_09_18/director_inverse.jsonl.
This output is to be archived at a checkpoint; the program reads only tracked
input matrices. Bounded unsuccessful inverse searches are never called UNSAT.

## 4. Independent larger witness check

For the stored order9 minimizer, take p=(0,1,2,3,5,4,6,7,8) and
s=d=(1,1,1,1,1,-1,-1,-1,-1), with B_ij=s_i s_j A_(p_i,p_j).
Independent integer Walsh transformation over all262144 parent spins gives
exact minimum -33 and maximum33. This is an order18 upper witness, not a
global lower certificate. It is not conference: D^2-17I has entries -4,0,4
with counts28,252,44, and max|AC-CA|=4.

## 5. Unsuccessful stronger shortcuts

The tempting estimate min_twist Q(D)<=beta(A)+O(n) cannot hold uniformly
on all conference seeds: beta(A)<=n sqrt(n-1), whereas the preserved
universal lower bound at order2n has coefficient greater than1 in child
n^(3/2) units. This observation depends on the reported asymptotic lower
chain; it is a warning against extrapolating small finite beta gaps, not
a newly independent proof of that chain.

A pure operator-norm upper bound Q(D)<=n||D|| cannot certify normalized
caps below1/2: every order2n full signing has ||D||>=sqrt(2n-1) by its
Frobenius norm. A useful construction below that scale must retain a
Boolean-versus-spectral gap rather than pay the operator bound alone.
This restricts that certificate, not the underlying sign matrices.
