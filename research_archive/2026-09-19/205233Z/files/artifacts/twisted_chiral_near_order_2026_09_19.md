# Near-order restriction while retaining the twisted-double family

Status: director proof, September19; full independent reconstruction PASS,
with105 exact finite rerouting checks and exact Bernoulli expectation checks
in computations/results/twisted_chiral_near_order_audit_2026_09_19.json. This is a
uniform positive structural theorem, not the desired doubling inequality.
All norms use H_A(x)=sum_(i<j)a_ij x_i x_j and Q(A)=max|H_A|.

## Statement

For an order-n hollow full signing define

```math
 F(A)=\min_{g,d}Q\begin{pmatrix}
 A&g^TAg+\operatorname{diag}d\\
 g^TAg+\operatorname{diag}d&-A
 \end{pmatrix},\qquad
 \beta(A)=\max_{x,y\in\{\pm1\}^n}|x^TAy|.
```

Here g ranges over signed permutations and d over all sign diagonals.
For every integer 1<=r<=n/4 there is ONE principal set C, |C|=n-r, such
that all three differences

```math
 |Q(A[C])-Q(A)|,\quad |\beta(A[C])-\beta(A)|,
 \quad |F(A[C])-F(A)|
```

are at most

```math
 300\left[\frac rn\beta(A)+n\sqrt r\right]+r.       \tag{1}
```

No spectral assumption, optimality, or prescribed limiting constant is
needed. In particular, for Q(A)<=K n^(3/2), beta(A)<=4Q(A) makes (1)
O_K(n sqrt(r)). If r=O(n^(1-eta)), the defect is
O_K(n^(3/2-eta/2)). The child and the parent remain exact full signings,
and the recovered parent remains a twisted double of the restricted child.
The statement is existential; no polynomial-time optimizer is asserted.

## 1. Random row restriction

For R a coordinate set write

```math
 L_R(A)=\max_{y\in\{\pm1\}^n}\sum_{i\in R}|(Ay)_i|
       =\|P_R A\|_{\infty\to1}.
```

Let independent Bernoulli(p) coordinates choose R. Then

```math
 \mathbb E L_R(A)
 \le p\beta(A)+n\sqrt{2p(1-p)(n-1)}.               \tag{2}
```

Proof: put f_i(y)=|(Ay)_i|, and use an independent Bernoulli copy Z'.
Jensen and the supremum inequality give

```math
 \mathbb E\sup_y\sum_iZ_i f_i(y)
 \le p\beta(A)+\mathbb E\sup_y\sum_i(Z_i-Z'_i)f_i(y).
```

Conditionally on w_i=|Z_i-Z'_i|, the nonzero differences are independent
uniform signs. The scalar Rademacher contraction inequality for the
1-Lipschitz map u->|u| bounds the second term by

```math
 \mathbb E\sup_y\sum_i\epsilon_iw_i(Ay)_i
 =\mathbb E\|A^T(\epsilon w)\|_1
 \le\sum_j\left(2p(1-p)\sum_i a_{ij}^2\right)^{1/2}.
```

The contraction used here needs no factor2 because the supremum is not
made absolute. An elementary coordinate proof is: for arbitrary offsets
h_t, average sup_t(h_t+epsilon phi(t_i)) over one sign; the resulting
half-supremum over s,t contains phi(s_i)-phi(t_i), bounded by |s_i-t_i|.
Interchanging s,t makes its supremum equal to that with s_i-t_i. Iterate
over the coordinates. Thus (2) is self-contained.

## 2. Selecting a small set simultaneously for a prescribed permutation

Fix a permutation p0, and use a Bernoulli reservoir S of density p=2r/n.
Both S and p0(S) have the same distribution. Equation(2) gives

```math
 \mathbb E[L_S(A)+L_{p_0(S)}(A)]
 \le4\left[\frac rn\beta(A)+n\sqrt r\right]=:4W.
```

Paley--Zygmund gives Pr{|S|>=r}>=1/6 since E|S|=2r and
E|S|^2<=(2r)^2+2r. Markov gives probability at least11/12 that the
displayed sum is <=48W. The events therefore intersect. Select any r
points R inside such a reservoir. Monotonicity of L yields

```math
 L_{R\cup p_0(R)}(A)\le L_R(A)+L_{p_0(R)}(A)\le48W. \tag{3}
```

For a hollow symmetric matrix, removing all entries touching R produces
the error E_R=P_R A+A P_R-P_R A P_R. Each of its three bilinear terms
has norm <=L_R(A), hence

```math
 \beta(E_R)\le3L_R(A),\qquad Q(E_R)\le\tfrac32L_R(A). \tag{4}
```

These give the asserted Q and beta continuity bounds immediately.

## 3. Repairing the twist rather than abandoning its constraint

Choose a minimizer g,d defining F(A), and let p0 be its underlying
permutation, with g e_i=s_i e_(p0(i)). Select R by (3), and put C=R^c.
On C replace each cycle of p0 by the cycle obtained by deleting its R
vertices; on R use the identity. Call the resulting permutation p1.
Retain the old signs on unchanged columns and choose arbitrary signs on
changed columns, giving g1. Then g1(C)=C and g1(R)=R.

The changed-column set T is contained in R union p0^(-1)(R), and BOTH
p0(T) and p1(T) are contained in U=R union p0(R). For B0=g^TAg and
B1=g1^TAg1, B0-B1 is supported on rows or columns in T. Its row-restricted
bilinear norm is at most

```math
 L_T(B_0)+L_T(B_1)
 =L_{p_0(T)}(A)+L_{p_1(T)}(A)\le2L_U(A).
```

Apply the three-term incident-row identity (4) to the difference to get
beta(B0-B1)<=6L_U(A). Changing the bridge therefore costs at most this
amount in Q. Principal deletion of both copies of R cannot increase Q,
and its remaining bridge is exactly g1_C^T A[C] g1_C plus diag(d_C).
Consequently

```math
 F(A[C])\le F(A)+6L_U(A).                          \tag{5}
```

In the reverse direction, extend an optimal twist of A[C] to fix R.
Adding back A's incident edges in the two diagonal blocks costs at most
2Q(E_R)<=3L_R(A), the incident bridge costs at most3L_R(A), and its r
matching edges cost at most r. Therefore

```math
 F(A)\le F(A[C])+6L_R(A)+r.                        \tag{6}
```

Equations(3)--(6) prove (1), with room in its numerical constant.

## 4. What this removes, and what it does not

Ordinary principal restriction preserves the signing constraint but need
not preserve the relation B=g^TAg. The cycle repair above removes THAT
obligation with an explicit subleading error when r=o(n). It also shows
that discarded vertices cannot be charged at the crude nr scale here.

This is NOT a proof of F(A)<=2sqrt2 Q(A)+o(n^(3/2)). It provides no
transfer across a fixed proportional gap. Summing many small deletions
does not turn the square-root modulus into a vanishing macroscopic loss.
Thus it can fill relative-o(1) gaps in an already dense good sequence,
but cannot by itself turn a single doubling subsequence into convergence.
The separate multiplier/thinning obligations in
`twisted_chiral_uniform_2026_09_18.md` remain.

The random-row estimate is a classical symmetrization/contraction mechanism;
the application retaining a same-seed twisted orbit uses the explicit
permutation-cycle repair. External novelty has not been established.

## 5. The discarded-row norm has the sharp square-root scale

The n sqrt(r) term is not merely an artifact of the upper-bound proof for
L_R. For EVERY full signing A and EVERY r-set R, random signs epsilon on R
give

```math
 L_R(A)\ge\mathbb E_\epsilon\|A^T P_R\epsilon\|_1
       =(n-r)a_r+r a_{r-1},\qquad
 a_k=\mathbb E|\epsilon_1+\cdots+\epsilon_k|.        \tag{7}
```

There are exactly r nonzero unit coefficients in each column outside R
and r-1 inside R; changing their fixed signs does not affect the law.
In particular a_k~sqrt(2k/pi), so if r->infinity and r=o(n), the right
side is (sqrt(2/pi)+o(1))n sqrt(r). The elementary exact formulas are
a_(2j)=2j binom(2j,j)/4^j and
a_(2j+1)=(2j+1)binom(2j,j)/4^j, with a_0=0.

Thus an argument paying the full discarded-row response separately cannot
replace this scale by o(n sqrt(r)). This is a lower bound on the response
norm, NOT a lower bound on |F(A[C])-F(A)|, which may be much smaller.

## 6. An audited unsuccessful connection to earlier regularization

The earlier cloned-block preparation theorem supplies actual near-minimizers
whose near-ground window o(n^(3/2)) has exp(o(n)) physical states. That does
not immediately make the signed-permutation selection lemma here effective:
a dangerous pair in |H_A(x)-H_A(y)|+|x^TBy| can have an internal-energy
difference a fixed amount below 2Q(A). It need not belong to the prepared
shrinking near-ground window. Paying the whole bridge by its bilinear norm
does not remove these outside states. No uniform control of that missing
macroscopic energy band has been established. This attempted combination
does not supply a doubling defect, and is not counted as a new reduction.

## 7. Complementary insertion theorem

Section16 of `twisted_chiral_symmetry_followup_2026_09_19.md` constructs an
actual extension in the other direction, retaining a prescribed child AND
its prescribed twisted parent as induced submatrices. Extending the signed
permutation by the identity and sampling only new child edges gives parent
error at most

```math
r+\sqrt{2[8nr+4\binom r2]\,[2(n+r)+1]\log2}.
```

The same completion controls the child's quadratic and bilinear caps.
Combined with the deletion theorem above and a crude uniform cap bound,
this proves `|T_(n+r)-T_n|<=C n sqrt(r)` for `T_n=min_A F(A)` and
`1<=r<=n/4`. The proof and its independent audit retain all dependencies
of the parent entries. Ordinary independent parent-edge completion would
not preserve the defining twisted-family relation.
