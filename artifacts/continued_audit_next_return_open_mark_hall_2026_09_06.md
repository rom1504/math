# Partial next-return theorem: open-mark Hall bounds

Date: 2026-09-06. Fixed-polynomial-stage independent reconstruction.
This closes the full-contraction condition for the part of the new
return having at least two old-noise branches. It does not close the
entire second-return criterion or establish arbitrary-depth closure.

Use the first marked-history notation

```math
W=(S,G,Y,QS,QD),\quad D=S h_2(G),\quad Y=BD,
\quad Z=Br(G,Y),
\quad L=Bc_0+BD_aB r(G,Y).
```

The next remainder is `eta=B(C-c0-D_a Z)`. At a fixed polynomial stage,
its squarefree main separates into a centered-coefficient one-noise
part and terms with at least two Z branches. All local collisions
touching a Z branch are removed first by the already audited quantitative
averaged-L2 surgery. Internal coherent Boolean contractions remain exact.
Write epsilon=n^-1/2 and allow fixed powers of log n in the bounds below.

## 1. The proved partial statement

Let X be a fixed main term having at least two homogeneous Z branches,
possibly multiplied by a fixed polynomial coherent coefficient in W.
Let K_i,p be an exact degree-p Walsh kernel of BX. For either

```math
L^{(0)}=Bc_0(W),\qquad L^{(1)}=BD_aB r(G,Y),
```

and every fixed q>p, its exact degree-q Walsh kernel J_i,q satisfies

```math
\frac1n\sum_i\|K_{i,p}\star_p J_{i,q}\|_F^2\longrightarrow0.
                                                               (1)
```

The contraction includes exact distinct-label restrictions and the
q-p uncontracted right labels. These open labels are the point of this
note; they are not replaced by a scalar energy test.

## 2. The source contraction is small with its marks left open

Write the right source as a fixed polynomial P_b in W_b (for L^(0)),
or in G_b,Y_b (for L^(1)). Expand its exact Walsh kernel into products
of primitive exact Walsh kernels. The primitives have original degree
one or three. All repeated internal coherent seed labels are represented
by their exact finite equality patterns; no coherent pre-Walsh tensor
is declared close in L2 to its squarefree projection.

For each source pair a,b form the full-left contraction with open marks

```math
H_{a,b,U}=X_{a,p}\star_p P_{b,q},\qquad |U|=q-p.
```

Each left Z branch has at least three distinct slots. Every such slot
is contracted into the right kernel. The local surgery ensures that
different left noise branches and the left coherent factor do not share
these slots. A right primitive has at most three slots. Thus any set
of k left noise branches has at least k neighboring right primitives:
there are at least 3k distinct matched labels, and each right primitive
can supply at most three of them. This remains true if a label appears
in several right primitives through a coherent Boolean hyperedge.
Hall's theorem supplies distinct matched primitives for two selected
left noise branches.

Each matched merge costs O(epsilon), with three cases:

1. A proper portion of the noise is summed against its matched primitive.
   Its small proper cut supplies epsilon.
2. Some shared label remains open because it is a right free mark or
   also occurs at another right primitive. Fixing that label supplies
   epsilon through the noise's small fixed-slot influence. The remaining
   shared labels can be summed without losing this bound.
3. A degree-three noise and a degree-three coherent primitive are fully
   contracted, with no retained label. This is the exceptional comparison
   of transported h3(G) with Y or QD. Its covariance is O(epsilon)
   uniformly in the two roots: the exact source h3(G)/D covariance has
   bounded absolute row and column sums O(epsilon), and multiplication
   by the bounded-op B,Q preserves its operator bound.

No larger-degree noise can be consumed entirely by a single degree-at-
most-three primitive. The two chosen merges are vertex-disjoint, so
their gains multiply. Every subsequent merge is Hilbert-contractive,
including retained Boolean hyperedge labels: blockwise multiplication
and summation satisfy the usual Frobenius Cauchy--Schwarz bound. The
remaining primitive row norms and fixed polynomial moments are bounded.

Exact distinctness is a finite inclusion-exclusion over label equalities.
A repeated label within a primitive vanishes. Other retained equalities
are covered by the same merge argument. In particular a right free mark
does not become an uncharged scalar summation. It remains a tensor index
through the final Frobenius norm. Consequently

```math
\sup_{a,b}\|H_{a,b,\cdot}\|_F
       =O(\epsilon^2\operatorname{polylog}n),
\qquad \|H\|_F=O(\operatorname{polylog}n).             (2)
```

This is a two-root Frobenius estimate, not merely an entrywise scalar
covariance estimate.

## 3. Outer transport and the order of norms

Let M=B for L^(0), and M=BD_aB for L^(1). The output-root contraction
is exactly

```math
(K\star_p J)_{i,U}=\sum_{a,b}B_{ia}M_{ib}H_{a,b,U}.
```

Apply B and M separately on the two root axes of H, then restrict their
two output roots to the diagonal i=i. The restriction is a coordinate
projection of norm one. Thus

```math
\|K\star_p J\|_F
 \le\|B\|_{op}\|M\|_{op}\|H\|_F
 =O(\operatorname{polylog}n).
```

Dividing its square by n proves (1). An entrywise row-l1 estimate would
lose this conclusion; it is essential to retain both outer root axes
and apply the bounded operators before the diagonal projection.

The removed left local-collision errors have quantitative averaged L2
size n^-1 times a fixed logarithmic factor. The right fixed-polynomial
kernel has at most polylogarithmic local norm by its global-cut bounds.
The full contraction is bounded by the product of these norms, so these
errors remain negligible in (1). This passage is only at fixed degree.

## 4. The remaining centered-coefficient obstruction

For X=(A-EA)Z there is only one left noise branch. Hall alone supplies
one epsilon, which is insufficient after summing two unrestricted roots.
The old scalar-energy proof separated split and unsplit right G/Y
forests. In an unsplit open-mark extension, the centered coherent bridge
has global Frobenius size at most sqrt(n) times a logarithmic factor,
using its root map and the other kernel's all-global cut. A uniformly
epsilon noise factor would then suffice. Pins to free marks can be
charged by fixed-slot influence rather than discarded.

However, for the right source c0(W), a QD primitive may split between
the centered coherent factor and the noise. Unlike old Y, QD does not
have small proper local cuts. The old two-gain proof therefore cannot
be imported verbatim. The explicit marked structure of QD may still
give a global bound, but that estimate is not proved in this note.

A stronger sufficient route would be to show that every original Walsh
degree above three in Bc0(W) has small proper local cuts. That would
make the full left contraction a proper cut of the higher-degree right
kernel. Existing all-global-cut bounds alone do not prove this: exact
internal coherent contractions must be handled rather than suppressed.
No such stronger theorem is claimed here.

Accordingly, equal-degree covariance remains retained, the actual L is
unchanged, and the full second-return conditional Gaussian comparison
is still open at the centered-coefficient/coherent-return subcase.

## 5. A proved contracted-two-factor transport lemma

The previously established uncontracted transport lemma has this exact
extension. Suppose rooted squarefree tensors A_j and C_j have all global
cut norms at most C_A,C_C. Contract any k seed slots between them and
leave at least one free seed slot in EACH factor. For any row weights w_j,
the tensor

```math
T=\sum_j w_j(A_j\star_k C_j)
```

has every proper local cut bounded by
`max_j|w_j| C_A C_C`, up to the fixed-degree factors for symmetrization
and distinctness. Thus flat B transport supplies epsilon even when
the two factors have been partially contracted with each other.

For completeness, take a proper output cut L|R. As both factors retain
free slots, possibly swapping their names or transposing the cut gives
nonempty L_A from A and nonempty R_C from C. Write their other free
slots as R_A and L_C, and their contracted slots as K. Regard C as the
global map

```math
R_C\longrightarrow(j,L_C,K),
```

and A as the transposed global map

```math
(j,R_A,K)\longrightarrow L_A.
```

Start with R_A,R_C, apply C with R_A as a spectator, apply the diagonal
weight w_j, and then apply A with L_C as a spectator. The resulting
map is exactly the required flattening: both j and K are contracted
once, not traced over a surviving repeated tensor axis. Its norm is
bounded by the product of the two global norms and max|w|. If the first
orientation is unavailable, the swapped one is available unless one
factor has no free slot, precisely the excluded case.

The contracted K labels meeting a free label within either source give
zero by individual squarefreeness. Remaining cross-free exclusions are
rectangular pinching complements and cost only degree-dependent factors.
There is no uncontrolled partially traced identity in this proof.

In a primitive-times-polynomial expansion of a coherent source, this
settles every product-formula term for which neither factor is completely
consumed. The remaining terms completely consume a primitive of degree
at most three into the other exact polynomial kernel. A full high-degree
diagram still has at least two original primitive vertices carrying free
slots, but selected diagram subkernels have not been shown to inherit
the all-global-cut bounds of the unselected coherent polynomial. That
additional graph-level issue remains open here.
