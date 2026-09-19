# Independent audit: near-order repair and the sharp diamond bound

September19, uniform-selection track. This audit independently reconstructs
the proofs in `twisted_chiral_near_order_2026_09_19.md` and the sharp-norm
section of `twisted_chiral_uniform_2026_09_18.md`. All energies are
`H_A(x)=x^T A x/2`; `beta(A)=max_(x,y signs)|x^T A y|`.

The accompanying solver-free script
`computations/twisted_chiral_uniform_final_audit_2026_09_19.py` independently
checks10000 scalar contraction instances,1700 cycle repairs,800 insertion
coefficient expansions, and all4096 bilinear spin pairs of the weighted
sharpness matrix using plain scalar arithmetic. All pass; the results
JSON has the same basename under `computations/results/`. These finite
checks supplement, rather than replace, the universal proofs below.

## 1. Near-order theorem: PASS

The stated `300[(r/n)beta(A)+n sqrt(r)]+r` simultaneous error is valid.
In fact the displayed estimates suffice with288 in place of300. The
larger stated constant is retained.

### Factor-one contraction

For arbitrary offsets h_t and a 1-Lipschitz scalar function phi,

```math
 E_e\sup_t(h_t+e\phi(z_t))
 ={1\over2}\sup_{s,t}(h_s+h_t+\phi(z_s)-\phi(z_t))
 \le {1\over2}\sup_{s,t}(h_s+h_t+|z_s-z_t|)
 ={1\over2}\sup_{s,t}(h_s+h_t+z_s-z_t).
```

The last equality follows by interchanging s and t. Iterating this
one-coordinate inequality proves contraction with constant one for a
supremum **without** an additional outside absolute value. It applies
conditionally to `w_i=|Z_i-Z'_i|`, since the nonzero differences have
independent uniform signs. Cauchy--Schwarz then gives exactly

```math
 E L_R(A) <= p beta(A)+n sqrt(2p(1-p)(n-1)).
```

The initial symmetrization is also valid: first separate
`p sum_i |(Ay)_i|<=p beta(A)`, then apply Jensen to the independent
ghost Bernoulli variables inside the remaining supremum. No hidden
factor two is needed at either step.

### Simultaneous reservoir and cycle repair

For p=2r/n and `W=(r/n)beta(A)+n sqrt(r)`, the expected sum of the
two row norms is at most4W. Paley--Zygmund gives probability at least1/6
of a reservoir with at least r vertices; Markov excludes row-norm sum
larger than48W with probability at most1/12. Hence the two desired
events intersect with probability at least1/12. Restricting to any
r-point subset preserves the bounds by row-norm monotonicity.

The permutation containment was checked independently, including cycles
with zero or one surviving vertex. If i outside R changes image after
cycle skipping, its original successor lies in R; its new successor is
the first surviving point after a nonempty run of R vertices and
therefore belongs to p0(R). If i lies in R, its repaired image is i.
Thus the changed-column set obeys

```math
 T subset R union p0^(-1)(R),
 p0(T) union p1(T) subset U:=R union p0(R).
```

Retaining signs on all unchanged columns is important. Under that
convention the difference of the two conjugated bridges is supported
on rows or columns in T, and
`L_T(g^T A g)=L_(p0(T))(A)` is exact, not merely an upper bound.
The three-term incident-row decomposition then gives bridge norm at
most6L_U(A).

The reverse extension also preserves the required family: extend the
restricted signed permutation by the identity on R. The incident
diagonal-block error costs at most3L_R, its conjugated bridge costs
at most3L_R, and the new matching costs r. This proves the reverse
error6L_R+r. Together with L_U<=48W, both directions are bounded by
288W+r. The Q and beta errors are smaller, at most72W and144W.

### Exact coverage scope

If Q(A)<=K n^(3/2), put delta=r/n<=1/4. The theorem and the elementary
bound `F(A)<=2Q(A)+beta(A)+n<=6Q(A)+n` imply

```math
 |F(A[C])/(n-r)^(3/2)-F(A)/n^(3/2)|
 <= (1-delta)^(-3/2)
       [300(4K delta+sqrt(delta))+delta/sqrt(n)]
    +(6K+n^(-1/2))[(1-delta)^(-3/2)-1].
```

Thus all three normalized quantities have modulus `O_K(sqrt(delta))`.
For an already available good sequence of child orders N_j with
N_(j+1)/N_j->1, choose the next N_j above an arbitrary target m and
apply the theorem with r=N_j-m. This fills all child orders while
retaining one common restricted child and its admissible twisted parent.
The resulting parents have all even orders; ordinary principal deletion
of one parent vertex supplies odd-order original signings with no
asymptotic normalization loss. Odd orders are not themselves asserted
to be twisted doubles.

A single dyadic sequence has no such density. Nor does summing small
deletions give a vanishing macroscopic error: the square-root costs do
not telescope. Therefore this audit finds no implicit convergence or
order-coverage claim beyond the theorem's stated relative-o(1) scope.

The lower row-response formula
`L_R >= (n-r)E|S_r|+rE|S_(r-1)|` also checks directly by averaging
the dual row signs. It proves sharpness of the row-norm scale, not a
lower bound on the actual change in F.

## 2. Sharp diamond theorem: PASS with the hollow hypothesis retained

Define Delta(A) as the maximum absolute value of
`H_A(a)-H_A(c)+b^T A d` under the two coordinatewise diamond constraints.
For hollow symmetric A, this objective is affine separately in each
pair (a_i,b_i) and each pair (c_i,d_i). Starting from a global maximizer,
successively replacing each pair by a maximizing diamond vertex gives
a vertex maximizer. This is the exact reason hollowness is required;
it is not a general assertion that a nonhollow quadratic maximizes at
vertices of a polytope.

At vertices, simultaneous coordinate switching by a_i+b_i puts every
coordinate into one of the eight listed types. The certificate includes
all64 type-pair entries, including equal-type entries. The latter matter
when a type has more than one physical coordinate, even though A itself
has zero physical diagonal. The exact identity

```math
 sum_r w_r(X_r Y_r^T+Y_r X_r^T)=24K,
 sum_r w_r=32,
```

gives `F=sum_r(w_r/24) X_r^T A Y_r`, because A is symmetric and
`F=Tr(AK)/2`. Hence Delta(A)<=4beta(A)/3 with the stated normalization.

The weighted six-type sharpness matrix M is NONHOLLOW. Its use is
nevertheless sound: the proof does not apply the vertex theorem to M
or assert that its displayed diamond witness is a global maximizer.
It only needs the exact finite facts beta(M)=18 and a diamond value24.
For the hollow blowup

```math
 M_r=M tensor J_r-diag(M tensor J_r),
```

the replicated witness has value24r^2-10r, while
`|beta(M_r)-18r^2|<=14r`. The loss10r is the explicit diagonal
contribution of the witness;14r is the diagonal absolute mass.
Thus hollow matrices approach ratio4/3 from below, and the already
proved hollow upper bound finishes sharpness.

### Full-sign realization and its scale

Scale M by1/3 and realize each off-diagonal block mean using independent
full-sign edges, with symmetry imposed by sampling only unordered
pairs. For N=6r, the centered error E obeys beta(E)=O(N^(3/2)) with
positive probability: for fixed x,y,
`x^T E y=sum_(i<j)E_ij(x_i y_j+x_j y_i)`, and Hoeffding plus a union
bound over4^N sign pairs gives this estimate. For example threshold
3N^(3/2) already makes the failure probability exponentially small.

The same fixed diamond witness changes by at most2beta(E), using the
real-cube bilinear norm for its four bounded vectors. Consequently
the full-sign realization has

```math
 beta(A)=6r^2+O(r^(3/2)),
 Delta(A)>=8r^2-O(r^(3/2)),
 Q(A)=3r^2+O(r^(3/2)).
```

The matching upper bounds follow from norm perturbation and the hollow
diamond theorem. This is valid asymptotic full-sign sharpness at
**quadratic cap scale**, not at the competitive n^(3/2) scale.

### Nonhollow block averages: finite correction confirmed

For a scalar diagonal entry the diamond expression has absolute maximum
one. To see this without incorrectly assuming vertex reduction, write
a=|a|,c=|c| in [0,1] and bound the positive side by
`.5a^2-.5c^2+(1-a)(1-c)`. It is convex in a, so a=0 or1 maximizes it;
both endpoints are at most one. Interchanging the two pairs and negating
one bridge variable treats the negative side.

Therefore `Delta(diag t)=sum|t_i|`. Hollowing a general symmetric P gives

```math
 Delta(P)<=4beta(P)/3+(7/3)sum_i|P_ii|.
```

This validates the corrected coarse native bound
`F(A)<=4beta(A)/3+(40/3)epsilon n^2+4kn+(10/3)n`.
The previous last term n omitted the hollowing cost. The correction
does not affect any hollow theorem, the clone corollary, or the coarse
O(n^2/sqrt(log n)) asymptotic remainder.

No native4beta/3+O(n) theorem or convergence conclusion is inferred.

## 3. Same-child-orbit insertion and optimized-family modulus: PASS

This independently audits Section16 of
`twisted_chiral_symmetry_followup_2026_09_19.md`.

Fix a child A and a particular allowed parent D using signed permutation
G. Add r vertices with independent fair child-edge signs and use G'=G
direct-sum I_r. This retains both A and the **chosen** D as exact induced
submatrices. It also retains the defining relation between the new child
and its bridge. Sampling independent parent edges instead would not.

At an old--new edge (i,j), put a=x_j,b=y_j. Direct expansion from the
parent energy gives coefficient

```math
 c_ij=a[x_i+(Gy)_i]+b[(Gx)_i-y_i].
```

Hence `c_(.j)=T(x,y)` for `T=[aI+bG,aG-bI]`. Orthogonality of G gives
`TT^T=4I`: the two cross terms `ab(G+G^T)` cancel exactly. Thus its
squared norm is at most8n. The new--new coefficient is
`x_i x_j-y_i y_j+x_i y_j+x_j y_i`, always +/-2. The total squared
coefficient bound is therefore precisely

```math
 V=8nr+4 binom(r,2).
```

The claimed sharpness of this coefficient bound also checks: with even
n, take G skew orthogonal, y=-Gx, and all new x_j=y_j=1. Then the
old--new vector is `2(I+G)x` and has squared norm8n; all new--new
coefficients simultaneously have square4. This is sharpness of the
uniform coefficient argument, not a parent-cap lower theorem.

Let `N=n+r` and `L=nr+binom(r,2)`. The three union bounds are valid:

| increment | squared-coefficient bound | number of one-sided tests |
| --- | ---: | ---: |
| child Q(E) | L | 2^N |
| child beta(E) | 4L | 2^(2N-1) |
| matching-free parent core | V | 2^(2N-1) |

For the first row, the test count includes both tails after quotienting
by overall spin reversal. For beta, its absolute maximum equals its
positive maximum by reversing one argument; simultaneous reversal
preserves each value. For the parent core, the chiral rotation
`(x,y)->(y,-x)` reverses the energy, while overall reversal preserves it.
Thus neither of the latter rows needs an extra factor two.

Substituting the displayed thresholds in (35) gives failure probability
at most1/4 in each row. Their union has probability at most3/4, so one
and the same completion satisfies all bounds. The arbitrary prescribed
new matching costs only r. There is no independence assumption between
the three events or between different parent configurations.

For the optimized family `T_n=min_A F(A)`, the crude identity-twist
random construction has variance4 binom(n,2), yielding exactly

```math
 T_n<=n+sqrt(4n(n-1)(2n+1)log2)
    <=[1+sqrt(8log2)] n^(3/2).
```

For any chosen bridge, take a positive child ground state and a negative
one, and reverse an entire parent half to choose the favorable bridge
sign. This proves parent cap at least U+V>=Q(A), so any child optimizing
T_n has Q(A)<=T_n and beta(A)<=4T_n. The cap hypothesis needed for the
near-order deletion theorem is therefore automatic and uniform.

Insertion from an optimizing order-n child gives
`T_(n+r)<=T_n+O(n sqrt(r))`. Deleting r coordinates from an optimizing
order-(n+r) child gives the reverse bound. The restrictions
r<=n/4 imply both r<=(n+r)/4 and n+r<=5n/4, so all constants remain
absolute. Consequently

```math
 |T_(n+r)-T_n|<=C n sqrt(r),
 |T_(n+r)/(n+r)^(3/2)-T_n/n^(3/2)|=O(sqrt(r/n)).
```

The same-child family now has genuine constructive operations in both
near-order directions. These estimates still allow slow oscillation and
do not supply a fixed-ratio transfer or convergence theorem.

## 4. Exact order20 response diagnostics

After completing the requested proof audits, direct integer enumeration
of all2^19 projective spins independently checked the two retained cap40
parents in `twisted_chiral_2026_09_19_joint_target38.json`:

| parent source SHA prefix | Q | beta | projective Q extrema per sign | projective beta-max responses |
| --- | ---: | ---: | ---: | ---: |
| d327680f (Hadamard-completable) |40|84|96|3600|
| 247d887a (non-completable) |40|88|60|180|

All row fields are odd and therefore nonzero, so each maximizing response
spin has a unique aligned bilinear partner. The beta-maximizer coordinate
product has positive-support histogram
`{4:120,8:720,10:1920,12:720,16:120}` for the first parent and
`{6:30,8:60,12:60,14:30}` for the second.

The first parent satisfies beta/Q=2.1<3sqrt(2)/2, and hence its **diamond
resource** is bounded by112<2sqrt(2)*40. The second is not certified by
that sufficient condition. Neither assertion supplies a native twist.

The direct-enumeration program and output are
`computations/twisted_chiral_parent20_beta_2026_09_19.py` and
`computations/results/twisted_chiral_parent20_beta_2026_09_19.json`.
The JSON contains both matrices, all counts, and exact bilinear witnesses
independently rechecked by plain scalar summation.
