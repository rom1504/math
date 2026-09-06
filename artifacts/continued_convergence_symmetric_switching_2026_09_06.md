# Symmetric switching as an order-transfer candidate (2026-09-06)

Status: the local move and witness-transport identities below are proved and
verified exactly. No order-transfer estimate, asymptotic cap improvement, or
convergence theorem follows from the finite experiments.

## 1. Why this candidate was considered

The original target is convergence of
\(M_n/n^{3/2}\), where \(M_n\) is the minimum absolute Boolean quadratic
cap over symmetric hollow sign matrices. An upper recurrence at two
multiplicatively independent integer multipliers, with summable normalized
defects, would suffice. The exact sufficient criterion and its proof already
appear in `fresh_limit_algebra_2026_09_05.md`; they are not new here.

Ordinary tensor lifts expose extra vector-valued spin choices. The candidate
here was to change those choices by non-monomial, sign-preserving orthogonal
conjugations of the lift, retaining its spectrum but not its Boolean cap.
This is an attempted refinement of the tensor-transfer architecture, not a
fourth convergence architecture or a claimed theorem about minimizers.

## 2. Exact local move

Let \(H\) be a full symmetric sign matrix, let \(q\) be four coordinates,
and let \(v\in\{-1,1\}^4\), extended by zero outside \(q\). Set
\[
 R=I-\frac12vv^T.
\]
Since \(v^Tv=4\), \(R^T=R\) and \(R^2=I\). Suppose every coordinate of
\(v^TH[q,:]\) belongs to \(\{-4,0,4\}\). Then the off-block entries of
\(RHR\) are signs: a four-coordinate sign column orthogonal to \(v\)
is fixed, and a column equal to \(\pm v\) is negated.

The remaining obligation is the internal block. With \(K=H[q,q]\) and
\(V=vv^T\), it is precisely
\[
 4K'=4K-2VK-2KV+VKV\in\{-4,4\}^{4\times4}.                 \tag{1}
\]
When (1) holds, \(H'=RHR\) is a full symmetric sign matrix and has the
same spectrum and trace as \(H\). In particular, a symmetric Hadamard
matrix remains a symmetric Hadamard matrix. Hollowing both matrices changes
each quadratic form by the common trace; it must not be silently omitted
for full matrices with nonzero trace.

The four-row closed-quadruple condition (their entrywise product is constant)
supplies the first hypothesis when \(v\) is any one of the column patterns,
but does **not** automatically supply (1). For example, a full order-four
Sylvester block fails (1). The verifier checks both conditions and the full
integer identity \((2R)H(2R)=4H'\).

## 3. Boolean eigenvectors have births as well as survivors

For a sign vector \(x\), the restriction \((Rx)_q\) is Boolean exactly
when
\[
 \prod_{i\in q}v_i x_i=1.
\]
These are the patterns with \(v^Tx\in\{0,\pm4\}\). In the other half
of the patterns, \(v^Tx=\pm2\), and \((Rx)_q\) has three zero entries
and one entry of magnitude two. Its squared Euclidean norm is still four.

Write \(E_\lambda(H)\) for an eigenspace. Orthogonal similarity gives the
exact identity
\[
 E_\lambda(H')\cap\{-1,1\}^n
 =R\left(E_\lambda(H)\cap R\{-1,1\}^n\right).             \tag{2}
\]
Thus the new Boolean eigenset consists of:

- old Boolean eigenvectors obeying the even-parity constraint, transported
  by \(R\);
- old eigenvectors Boolean outside \(q\), with the specified zero/zero/
  zero/\(\pm2\) pattern inside \(q\), transported by \(R\).

Consequently, an argument counting only surviving old Boolean witnesses is
invalid. With overlapping switches, the preimage of the cube develops a
larger rational-coordinate reservoir, and (2), not just a list of parity
constraints on the original cube, is the exact condition to control.

## 4. Finite tests and actual scope

The reproducible verifier is
`computations/continued_convergence_symmetric_switching_2026_09_06.py`.
Starting with the indicated Sylvester ordering, the numbers of valid
nontrivial moves at orders 4, 8, and 16 are respectively 0, 24, and 228.

At order 16, the initial cap is 32, exactly the spectral half bound.
There are 20 Boolean eigenvectors at each of the eigenvalues \(+4,-4\),
counting antipodes separately. The recorded eight-step walk changes those
counts to 12 each but leaves the cap 32. The separately seeded 200-step walk
also has cap 32 at every step. Its full matrices and moves are recorded in
`computations/results/continued_convergence_symmetric_switching_walk200_2026_09_06.json`.

This is evidence that the move changes the Boolean geometry, not evidence
of a smaller asymptotic cap. Even a finite matrix without Boolean spectral
eigenvectors would not establish a uniform asymptotic gap or transport a
growing exact minimizer.

The first sufficient theorem still missing is an optimized-orbit estimate
for growing minimizer seeds, at two independent multipliers or a dense
all-order mesh, with a power-saving normalized loss. Neither the local
identity nor the finite walk supplies any estimate of that kind.

## 5. Primary literature checked, with quantifiers retained

- [Orrick, *Switching operations for Hadamard matrices*](https://arxiv.org/abs/math/0507515):
  the usual closed-quadruple switch preserves the Hadamard property. It is
  a one-sided operation, so it does not establish symmetry of our conjugate.
- [Crnković, Egan and Švob, *A universal theory of switching for combinatorial
  objects, and applications to complex Hadamard matrices*](https://arxiv.org/html/2511.07020v1),
  Section 4.1 and Theorems 6.1 and 6.3: orthogonal field/block conditions
  justify Hadamard-preserving switches, but do not remove the internal
  block check in (1).
- [*Switching graphs and Hadamard matrices*](https://arxiv.org/html/2410.10638v1),
  Theorems 2.5 and 3.1: the graph correspondence uses McKay's doubled
  row/column representation, not symmetry of the original Hadamard matrix.
- [Shi et al., *Self-dual Hadamard bent sequences*](https://arxiv.org/pdf/2203.16439),
  Conjecture 1: the conjecture asserts existence for **some** Hadamard matrix
  at each even-square order. It is not a conjecture that every Hadamard
  matrix, or every symmetric one, has a Boolean eigenvector. Their
  low-dimensional eigenspace searches for nonsymmetric representatives
  must not be confused with the half-dimensional eigenspaces here.

The intended next useful step is a scalable orbit description or a precise
asymptotic obstruction. More finite walks alone are not the missing theorem.

## 6. Uniform arbitrary-seed rectangle gates

The following observation was supplied independently by the root researcher:
for a tensor \(B\otimes H\), the four rows indexed by
\((a,u),(a,v),(b,u),(b,v)\) are always closed. This requires no special
row-product relation in the seed \(B\).

There is a complete local classification. Let \(S=B[\{a,b\},\{a,b\}]\)
and \(T=H[\{u,v\},\{u,v\}]\). For either sign \(\sigma,\tau\), take
\(w=(1,\sigma)\otimes(1,\tau)\). The four-row dot products factor into
two factors in \(\{0,\pm2\}\), so the off-block condition holds.
The internal conjugate of \(S\otimes T\) is sign-valued if and only if
it is **not** the case that both \(S\) and \(T\) have opposite diagonal
entries.

Here is a proof, not just the enumeration check. Diagonal conjugations
normalize \(w\) to the all-one vector without changing either diagonal.
If both diagonals are opposite, the sums of the entries of \(S,T\) are
both \(\pm2\). In (1), division by four therefore makes every transformed
entry even, so it is not a sign. Otherwise suppose \(S\) has constant
diagonal \(d\). In the normalized basis, either \(S\mathbf1=0\), in
which case the internal block is unchanged, or \(S=dJ_2\). In the latter
case the reflector acts on the range of \(S\) by the signed permutation
\(I_2-\mathbf1\mathbf1^T\); the internal conjugate is
\(S\otimes (PTP^T)\), again sign-valued.

Thus completing an arbitrary hollow seed by \(B=A+I\) makes every such
cross-seed rectangle initially admissible, for every symmetric full sign
outer matrix \(H\). Exhaustion of all 256 local choices gives exactly
192 valid choices, in agreement with this classification.

## 7. Simultaneous disjoint-pair compatibility

The root researcher integrated the preceding rectangles into a stronger
exact gate: in the two seed channels \((1,1)/\sqrt2,(1,-1)/\sqrt2\),
apply independent signed outer-coordinate permutations \(O_+,O_-\).
For one seed pair and \(B=A+I\), the transformed matrix is always full
sign and symmetric. This follows because the seed pair's internal block
has rank one, and every outside seed link lies in a single channel.

For two disjoint seed pairs, normalize \(O_+=I\) and write the relative
permutations as \(P,Q\). A rank-one seed crossblock remains sign-valued
automatically. For a rank-two seed crossblock the exact compatibility
criterion is
\[
 H_{ij}(HQ^T)_{ij}(PH)_{ij}(PHQ^T)_{ij}=1
 \quad\hbox{for every }i,j.                              \tag{3}
\]
Indeed a rank-two 2-by-2 sign matrix remains a sign Hadamard matrix under
the two-channel basis change. Reconstructing each physical entry amounts
to applying a normalized order-four Hadamard transform to four signs.
That transform is Boolean precisely for the odd-parity four-sign inputs.
The seed channel coefficients already have odd parity, leaving (3).

For a Walsh outer matrix \(H_{xy}=(-1)^{x\cdot y}\), take translations
\(P=P_t,Q=P_u\). Then the left side of (3) is the constant
\((-1)^{t\cdot u}\). Thus mutually orthogonal translation directions
assigned to disjoint seed pairs give simultaneous admissible gates for
**arbitrary** seeds. An outer binary dimension at least the number of
seed pairs supplies such directions.

These gates retain all scalar seed witnesses. In fact every pure Boolean
tensor \(z\otimes f\) is taken to a Boolean vector: the two coordinates
of \(z\) in each seed pair occupy exactly one of the two channels, on
which the action is a signed permutation of \(f\). The possible useful
effect is therefore removal of additional vector-valued tensor witnesses,
not removal of the scalar seed objective itself.

This is a concrete larger admissible family. It does not yet bound its
optimized Boolean cap, and an arbitrary second pair gate after an arbitrary
first one is not automatically admissible; condition (3) is essential.

## 8. A stabilized landing obstruction for all pure-tensor-preserving gates

There is a stronger obstruction than retention of individual scalar seed
witnesses. Suppose \(R\) is orthogonal on \(\mathbb R^n\otimes\mathbb R^s\)
and maps every \(u\otimes v\), for Boolean \(u,v\), to a Boolean vector.
For a fixed Boolean \(v\), define the linear map
\[
 V_vu=R(u\otimes v).
\]
It maps the seed cube into the parent cube. In fact each row of \(V_v\)
is a signed coordinate functional: a linear form taking only the values
\(\pm1\) on an entire sign cube has exactly one nonzero coefficient,
of magnitude one. Orthogonality also gives \(V_v^TV_v=sI_n\).

For an arbitrary symmetric seed \(B\), put
\(L=R(B\otimes H_s)R^T\). Then exactly
\[
 V_v^T L V_v=(v^TH_s v)B.                                 \tag{4}
\]
For every further symmetric outer matrix \(K\), tensoring \(V_v\) with
the identity maps every Boolean witness for \(B\otimes K\) to a Boolean
witness for \(L\otimes K\). Consequently
\[
 Q(L\otimes K)\ge |v^TH_s v|\,Q(B\otimes K).             \tag{5}
\]
This statement uses the full quadratic form for matrices with diagonals.
Taking the common tensor-regularization supremum in the archived functional
\(\mathcal R\) gives
\[
 \mathcal R(L)\ge |v^TH_s v|\,\mathcal R(B).              \tag{6}
\]
In particular, if \(H_s\) has a Boolean spectral eigenvector, then
\[
 \frac{\mathcal R(L)}{s^{3/2}}\ge\mathcal R(B).             \tag{7}
\]

Thus a single integrated gate, or a compatible disjoint collection of such
gates, cannot lower the seed's **stabilized** tensor coefficient when the
initial outer matrix is Boolean-saturating. This rules out the proposed
strategy of gating once and then using arbitrary dense regular outer
padding to prove a smaller landing coefficient. Finite core-cap reductions
do not contradict (7).

The scope is exact: (7) does not exclude a direct order-transfer program
with new interleaved gates at each size, gates that abandon the full
pure-tensor invariant, nonsaturating initial outer factors, or a proof that
the ungated regularization already agrees with the original liminf. Those
are additional obligations, not consequences of the local switch.
