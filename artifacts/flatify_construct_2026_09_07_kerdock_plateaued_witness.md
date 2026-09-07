# An exact plateaued pentagon witness for quadratic-phase MUB lifts

Date: 2026-09-07. Status: **Pending independent audit**. This is a growing-
dimension Boolean obstruction to a concrete construction class, not a lower
bound on all full signings or a counterexample to the selected-child theorem.

## 1. An explicit maximal real MUB family

Let m be odd, K=F_(2^m), W=K direct-sum F2, and t=|W|=2^(m+1). For a in K,
put, with absolute field trace Tr,

```math
 B_a((x,s),(y,r))=
   \operatorname{Tr}(a^2xy)+\operatorname{Tr}(ax)\operatorname{Tr}(ay)
   +s\operatorname{Tr}(ay)+r\operatorname{Tr}(ax).          (1)
```

This is an alternating binary bilinear form: on a repeated vector the first
two terms cancel because Tr(a²x²)=Tr(ax), and the last two cancel.
For a!=b, let c=a+b, alpha=a/c, u=cx, v=cy. The difference becomes

```math
 \operatorname{Tr}(uv)+\operatorname{Tr}(\alpha u)\operatorname{Tr}(v)
 +\operatorname{Tr}(u)\operatorname{Tr}(\alpha v)
 +\operatorname{Tr}(u)\operatorname{Tr}(v)
 +s\operatorname{Tr}(v)+r\operatorname{Tr}(u).              (2)
```

If (u,s) is in its radical, the coefficient of r gives Tr(u)=0, and the
coefficient of v gives u+[Tr(alpha u)+s]1=0. Since m is odd, Tr(1)=1;
taking trace gives u=0 and then s=0. Thus B_a+B_b is nonsingular.

Choose a quadratic Q_a polarizing B_a, by fixing coordinates on W and
using the upper-triangular coefficients of B_a. Let H be the Walsh matrix
H(z,w)=(-1)^(z dot w), and set

```math
 U_a=t^{-1/2}H\operatorname{diag}((-1)^{Q_a}).              (3)
```

Together with Id, these are t/2+1 mutually unbiased real bases. Indeed the
Walsh transform of a quadratic with nonsingular polar has every magnitude
sqrt(t). One elementary proof squares the transform and changes variables
by the difference of the two inputs; orthogonality leaves only difference
zero. Therefore U_a U_b^T has every entry of magnitude t^(-1/2).

The result below applies more broadly to ANY five quadratic-phase MUBs of
the form (3), without requiring this particular field construction.

## 2. A four-coordinate real vector witness for the actual optimal A5

Use the following switching of the pentagon signing, exactly matching the
archived order-five witness:

```math
 A=\begin{pmatrix}
 0&1&1&1&1\\1&0&1&1&-1\\1&1&0&-1&-1\\
 1&1&-1&0&1\\1&-1&-1&1&0
 \end{pmatrix},\qquad Q(A)=4.
```

The five unit vectors

```math
 v_0=(-1,1,1,-1)/2,\quad
 v_1=v_2=(-1,1,-1,-1)/2,\quad
 v_3=v_4=(0,0,1,0)
```

have the exact seed Gram energy

```math
 \sum_{i<j}a_{ij}\langle v_i,v_j\rangle=5.                 (4)
```

The first flat vector has an affine Boolean sign phase on F2², while the
second has a quadratic sign phase with polar bit1. This parity difference
is the mechanism exploited below; merely using vectors on the sphere
would not be enough.

## 3. Realization on a common two-dimensional spectral subspace

Let the five MUB phases have polars B_0,...,B_4. For any u!=0, the two
linear functionals (B_0+B_1)(u,.) and (B_0+B_2)(u,.) are nonzero and distinct:
each of these forms and their difference B_1+B_2 is nonsingular. Hence
the functionals are linearly independent. Choose v such that both have
value1 at v. Necessarily u,v are independent. On L=span(u,v), put
c=B_1(u,v)=B_2(u,v). Then B_0(u,v)=c+1.

Embed the five vectors of section2 at the four coordinates 0,u,v,u+v.
Multiply ALL of them by the same sign phase (-1)^(c alpha beta), where
w=alpha u+beta v in L. This common diagonal multiplication preserves
every Gram entry in (4). Call the resulting unit vectors w_i in R^t.

For i=0,1,2, the phase of diag((-1)^Q_i)w_i on L has polar bit1: the
bits are respectively (c+1)+0+c, c+1+c, and c+1+c. A two-dimensional
quadratic phase with polar bit1 has Walsh sums of magnitude2. Consequently

```math
 x_i=\sqrt t\,U_iw_i=H\operatorname{diag}((-1)^{Q_i})w_i
          \quad\text{is Boolean for every }i.             (5)
```

For i=3,4 this is immediate because w_i is a signed coordinate vector.
This is an exact actual-Boolean construction, not Gaussian rounding.

Form the full cross-fibre sign blocks

```math
 C_{ij}=\sqrt t\,a_{ij}U_iU_j^T\quad(i\ne j).
```

Equation (5) and (4) give the EXACT cross energy

```math
 H_C(x)=5t^{3/2}.                                        (6)
```

Thus even generic high-rank sums of the Kerdock phases cannot force the
lossless bound t^(3/2)Q(A)=4t^(3/2): a fixed two-dimensional adapted spectral
subspace already supplies a larger Boolean response.

## 4. Arbitrary within-fibre completion is cancelled exactly

For p,l in W, apply the SAME real Heisenberg signed permutation to every
spectral vector:

```math
 (T_{p,l}w)(z)=(-1)^{l\cdot z}w(z+p).                     (7)
```

It preserves all pair Gram products. It translates L to an affine subspace
and adds only a linear phase there. The restricted quadratic polar bits
in section3 are unchanged. Therefore every
x_i(p,l)=sqrt(t)U_i T_(p,l)w_i remains Boolean and has cross energy (6).

For each unit vector w, averaging (7) first over l kills all off-diagonal
entries of its rank-one covariance, and then averaging p makes the diagonal
uniform. Hence

```math
 \mathbb E_{p,l}[(T_{p,l}w)(T_{p,l}w)^T]=\operatorname{Id}/t,
 \qquad \mathbb E_{p,l}[x_i(p,l)x_i(p,l)^T]=\operatorname{Id}. (8)
```

It follows that ANY fixed hollow matrix D_i inserted in diagonal fibre i
has mean energy zero under this finite law. Thus every full sign completion F
of these cross blocks has an actual Boolean assignment with

```math
 Q(F)\ge\max_xH_F(x)\ge5t^{3/2}.                         (9)
```

This cancellation is justified by a proved covariance identity; it is not
unpaid cancellation between an unchanged child and a bridge.

## 5. Subsets containing the identity basis

Suppose a selected five-basis subset comes from a quadratic Kerdock family
and at least one quadratic basis U_b is unused. Right-multiply every selected
basis by U_b^T. This leaves every cross block U_iU_j^T exactly unchanged.

For a selected U_a with a!=b, the relative matrix is

```math
 U_aU_b^T=t^{-1}H\operatorname{diag}((-1)^{Q_a+Q_b})H.
```

The quadratic Q_a+Q_b has nonsingular polar. Completing its quadratic square
in the Walsh sum writes this Hadamard matrix as a row-sign/permutation times
t^(-1/2)H diag((-1)^R), where R is another quadratic phase. More explicitly,
its entries are a constant sign times
t^(-1/2)(-1)^(Q^*(z+w)); expanding the dual quadratic Q^* gives a row phase,
an invertible bilinear pairing (absorbed by a row permutation), and a column
quadratic phase. If Id itself was selected, Id U_b^T=D_b H/sqrt(t), which
becomes H/sqrt(t) after a row switching.

The row signs and permutations are physical cube automorphisms. They only
conjugate the arbitrary D_i, so (9) remains valid. The new column quadratic
phases have nonsingular pair differences because their bases remain mutually
unbiased. Section3 applies. Therefore (9) holds for EVERY five-basis subset
of the explicit maximal real family whenever an unused quadratic reference
is available, in particular for t>=16.

At N=5t the forced coefficient is 1/sqrt(5), versus the finite child's
row-normalized coefficient Q(A)/(5sqrt(4))=2/5. The conclusion is a precise
construction-class obstruction, not a claim about arbitrary minimizing
children whose orders tend to infinity.

## 6. Exact replay and failed numerical evidence

`computations/flatify_construct_2026_09_07_kerdock_diagnostic.py` constructs
the nine bases at t=16 using F8 and verifies every integer Hadamard identity.
Its original 40-restart searches found caps316--324 over126 five-basis subsets;
the lower316 was a heuristic miss, not a smaller cap.

`computations/flatify_construct_2026_09_07_kerdock_exact_witness.py` now
constructs an exact energy320 witness for ALL126 subsets. For subset
(0,2,4,6,8), which had the failed316 search, it explicitly checks all256
Heisenberg translates/modulations, their Boolean entries, preservation of
every Gram entry, and covariance sum256 Id in every fibre. The generated
JSON preserves all126 witnesses and unused references under computations/results.

The general proof above uses no external theorem. The finite code is a replay
of the exact algebra, not a proof of numerical global optimality.
