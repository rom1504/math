# Switching-only lifts are partial reversals of the seed

Status: exact algebra, with independent matrix and full-spin checks.
This distinguishes the two freedoms in the requested family; it is not
an asymptotic cap bound.

For a hollow signing A, a subset J, and diagonal signs d, let A^J be A
with exactly the edges whose BOTH endpoints lie in J negated. Let
s_i=-1 on J and +1 outside, and S=diag(s). Then

```math
 Q\begin{pmatrix}A&SAS+\operatorname{diag}d\\
 SAS+\operatorname{diag}d&-A\end{pmatrix}
 =Q\begin{pmatrix}A^J&A^J+\operatorname{diag}(s d)\\
 A^J+\operatorname{diag}(s d)&-A^J\end{pmatrix}.    \tag{1}
```

Proof: for i in J substitute x_i=y'_i and y_i=-x'_i, leaving the other
coordinate pairs unchanged. This is an actual signed permutation of the
2n Boolean coordinates. On a cross-pair edge i,j the energy is

```math
 a_{ij}(x_i x_j-y_i y_j)
 +(s_i s_j a_{ij})(x_i y_j+y_i x_j).
```

If neither endpoint is in J it stays unchanged. If one is in J it becomes
a_ij times BOTH new expressions, and if both are in J it becomes -a_ij
times both. The matching d_i x_i y_i changes to s_i d_i x'_i y'_i.
This proves(1), including exact matching signs and not merely the zero-diagonal
core.

Consequently minimizing over switching-only lifts is exactly the same as
minimizing untwisted lifts over these partial reversals of the seed and all
matchings. The transformed seed need NOT retain its original cap or be an
optimizer. Allowing arbitrary seeds, switching-only gives no larger set of
parent switching classes than untwisted doubling. For a FIXED child its
freedom is real, because partial reversal changes that child.

## Why permutation freedom is different

For two hollow full signings A,B put r_ij=a_ij b_ij. The following are
equivalent:

1. B=SAS for some diagonal signs S.
2. Every triangle obeys r_ij r_jk r_ki=+1.

The forward direction cancels the three vertex signs. Conversely fix a root0,
set s_0=1 and s_i=r_0i; its triangles force r_ij=s_i s_j.
For n<3 the same conclusion follows directly without triangles.

A local quarter-rotation of one chiral pair replaces its incident edge
parameters (a,b) by either (b,-a) or(-b,a). Hence it flips r on precisely
the incident edges. Rotating both ends leaves r unchanged. The triangle
products of r therefore survive ALL such local pair rotations, simultaneous
pair switches, and pair relabeling (up to relabeling the triangles).

For B a signed-permuted copy of A, triangle signs satisfy

```math
 \tau_B(i,j,k)=\tau_A(p(i),p(j),p(k)),
 \qquad \tau_A(i,j,k)=a_{ij}a_{jk}a_{ki}.
```

Thus a permutation can create a nontrivial relative triangle pattern
that switching alone cannot produce. These are familiar switching/two-graph
identities, not a newly claimed theory or an efficient test of all possible
permutations. Arbitrary pair rotations preserve full signs and chirality,
but need NOT preserve B being switching/permutation-equivalent to A.

## Exact finite discrimination

The companion switching-only enumeration checks all switches and ALL
matchings, with no permutation search. Permutations strictly improve the
best fixed-child cap for both successful order7 classes (25 to21), both
order8 classes (32 to30), order9 classes0 and2 (37 to33), and order9
class5 (39 to33). The other recorded classes have equal minima with or
without permutation. Class labels are those of the frozen census inputs.

All1024 partial reversals of the conference order10 seed have cap at least15:
82 have15, 762 have17, and180 have19. Thus (1) cannot turn its cap40 parent
into a parent of an optimal cap13 child. This agrees independently with the
exhaustive fixed-optimal-child cap44 result. These facts do not exclude a
different chiral decomposition of the parent or an asymptotic transfer theorem.
