# A failed quotient-rigidity implication and a canonical Walsh obstruction

Date: 2026-09-06. These are exact structural results about a candidate
order-transfer construction. They do not prove convergence or an asymptotic
upper bound for the original minimization problem.

## 1. Universal sign preservation does not propagate every Boolean quotient

Let \(H\) be a symmetric Hadamard matrix of order \(s\), and let \(O\)
be orthogonal such that \(OH\) and \(OHO^T\) are full sign matrices.
On one selected pair of seed fibres define
\[
 R=P_+\otimes I_s+P_-\otimes O,
 \qquad P_\pm=\frac12\begin{pmatrix}1&\pm1\\\pm1&1\end{pmatrix},
                                                               \tag{1}
\]
and make \(R\) the identity on every other seed fibre.

For **every** symmetric full sign seed \(B\) with diagonal entries one,
\(R(B\otimes H)R^T\) is a full symmetric sign matrix. The selected
principal pair block is \(2P_b\otimes H\), where \(b\) is its seed edge;
it becomes \(2P_b\otimes H\) or \(2P_b\otimes OHO^T\). Each outside
seed link has its two signs in one of the two channels and becomes a
signed copy of \(H\) or \(OH\). All remaining blocks are unchanged.

The all-one seed specialization is unchanged:
\[
 R(J_n\otimes H)R^T=J_n\otimes H.                           \tag{2}
\]
If \(Hv=\sqrt s\,v\) with \(v\) Boolean, then \(\mathbf1\otimes v\)
is therefore a Boolean spectral vector of this specialization. But a seed
spin with opposite signs on the selected pair is sent to vectors involving
\(Ov\). Whenever \(Ov\) is non-Boolean, that spectral vector has no
Boolean orbit of the asserted form over all seed spins.

This refutes the proposed implication

> universal sign preservation, together with a Boolean spectral vector of
> the all-one specialization, forces its common Boolean seed orbit.

The example satisfies the universal disjoint-support and Jordan-algebra
constraints, because those constraints follow from the universal sign
preservation that was just proved. They cannot establish the false
implication.

### Explicit integer example

Use the standard Walsh matrix of order 16, indexed by integers 0 through 15.
Take the four-coordinate reflection
\[
 O=I-ww^T/2,\qquad
 \operatorname{supp}w=(0,1,6,7),\quad w|_{\operatorname{supp}w}=(1,-1,-1,1).
\]
Both \(OH\) and \(OHO^T\) are full sign matrices. The vector
\[
 v=(1,-1,-1,-1,-1,1,1,1,-1,1,1,1,-1,1,1,1)
\]
satisfies \(Hv=4v\), but
\[
 Ov=(0,0,-1,-1,-1,1,2,0,-1,1,1,1,-1,1,1,1)
\]
is not Boolean. The combined gate (1) is itself the eight-coordinate
reflection \(I-zz^T/4\), with \(z=(1,-1)\otimes w\) on the selected
two fibres.

The script `computations/continued_convergence_quotient_counterexample_2026_09_06.py`
checks all identities using integer arithmetic, and checks all 64 full
constant-diagonal seeds at seed order four. The exact output is in the
corresponding `computations/results/` JSON. The symbolic argument above
proves universality at every seed order at least two.

Other common quotients remain in this particular example. The assertion
refuted is propagation of **every specified** Boolean spectral vector,
not existence of at least one surviving common quotient.

## 2. All initially admissible Walsh quadruple reflections have common bent witnesses

One proposed next construction was to place different operators \(O_j\)
from Section 1 on disjoint seed pairs, hoping to make the common preserved
Boolean spectral-vector intersection empty. If every \(O_j\) is a single
admissible symmetric closed-quadruple reflection of the **initial Walsh
matrix**, this is impossible even before imposing mutual compatibility of
the different gates.

Let the Walsh matrix be indexed by \(\mathbb F_2^d\), with \(d\) even,
and define
\[
 q_0(x)=\sum_{i<j}x_ix_j,
 \quad J=I+\mathbf1\mathbf1^T,
 \quad c=d/2+1\pmod2.
\]
For every linear coefficient \(\ell\) of parity \(c\), the Boolean vector
\[
 f_\ell(x)=(-1)^{q_0(x)+\ell\cdot x}                       \tag{3}
\]
is a positive or negative Walsh spectral eigenvector. Moreover **every**
initially admissible symmetric closed-quadruple reflection takes it to a
Boolean vector.

### Spectral assertion

The polar form of \(q_0\) is
\[
 \beta(a,b)=a\cdot b+(a\cdot\mathbf1)(b\cdot\mathbf1)=a^TJb.
                                                               \tag{4}
\]
Because \(d\) is even, \(J^2=I\), so this quadratic form is nonsingular
and its Walsh transform has magnitude \(2^{d/2}\). Completing the square
shows that its dual has quadratic part \(q_0(Jy)\), and
\[
 q_0(Jy)=q_0(y)+c\,\mathbf1\cdot y.
\]
For the linear perturbation (3), the dual linear coefficient is therefore
\(c\mathbf1+J\ell\). If \(\mathbf1\cdot\ell=c\), this equals
\(\ell\). Hence the transform is \(\pm2^{d/2}f_\ell\), as claimed.

### Admissibility assertion

A closed quadruple of distinct Walsh rows is an affine plane
\(a+\operatorname{span}\{r,t\}\). Its principal four-by-four block,
up to sign conjugation, is the character table of the restricted bilinear
form. Its real rank is \(2^k\), where \(k\) is the binary rank of
\[
 \begin{pmatrix}r\cdot r&r\cdot t\\r\cdot t&t\cdot t\end{pmatrix}.
\]
A four-coordinate reflection is internally sign-preserving exactly when
\(k\le1\). At rank zero it fixes the internal rank-one block. At rank
one its nontrivial action on the two repeated row classes is a signed
permutation. At rank two the block is an order-four Hadamard matrix and
the reflected internal entries are even, hence are not signs.

The condition \(k\le1\) is
\[
 (r\cdot r)(t\cdot t)+(r\cdot t)^2
 =\beta(r,t)=0.                                           \tag{5}
\]
But the product of (3) around the affine plane is precisely
\((-1)^{\beta(r,t)}\). Every admissible plane thus has even parity for
every vector in (3). Its reflection preserves Booleanity by the exact
four-coordinate parity criterion.

There are \(2^{d-1}\) such common projective witnesses (each (3) has first
coordinate one). At Walsh order 16 this gives eight. Direct enumeration of
the fifteen admissible two-dimensional direction spaces confirms that all
eight survive every one of them.

### Scope of this obstruction

For controlled gates built from these single initial reflections, the
common vector from (3) supplies the signed-coordinate embedding and the
stabilized lower bound in Section 8 of
`continued_convergence_symmetric_switching_2026_09_06.md`. Pairwise
compatibility of the operators cannot remove a vector already preserved by
all individual operators.

The proof does not apply to arbitrary products whose intermediate outer
matrices change, arbitrary non-quadruple orthogonal operators, or
minimizer-dependent gates. It should not be stated as saturation of every
non-Cayley symmetric Hadamard family.
