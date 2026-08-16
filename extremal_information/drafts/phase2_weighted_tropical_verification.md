# Independent audit: query-weighted tropical exposure

**Scope.**  This report audits
`phase2_weighted_tropical_exposure.md` only.  It checks the min-plus sign,
the four-error and weighted Cauchy constants, the matching disjointness step,
the zero-mass convention, the `D_r` calculation, and the max-plus remark.
It does not edit a surface theory file.

## 1. Verdict

**Theorem WE.1 is false as stated.**  An ordinary matching in the graph on
anchor indices need not give disjoint four-cell witness sets.  Repeated
anchor rows or columns can make rectangles belonging to endpoint-disjoint
edges overlap.  The proof then double-counts the squared error on their
intersection.

There is an exact `3 by 3` counterexample with four distinct anchors, all
cell masses positive, and a rank-one min-plus approximant.  Thus this is not
an edge case caused by duplicate anchors or zero query mass.

The theorem has two clean repairs:

1. require the selected witness rectangles themselves to be pairwise
   cell-disjoint; or
2. replace ordinary matchings by fractional packings of the four-cell
   witness sets.

Alternatively, the original matching formulation is valid if at least one
of the maps `i -> x_i` and `i -> y_i` is globally injective.  The `D_r`
application has this property and survives unchanged.

After that repair, the result is a useful and rigorous average-error
extension of the uniform tropical fooling-set theorem.  It is not a
characterization of average-error rank, however, and the weighted matching
parameter should be presented as a sufficient certificate rather than
“the exact additional datum.”

## 2. The local four-cell calculation is correct

Write

```math
\widetilde M(x,y)=\min_{1\le t\le k}
\{u_t(x)+v_t(y)\}.
```

Every factor term majorizes `Mtilde`.  If term `t` is tight at both anchors
`(x_i,y_i)` and `(x_j,y_j)`, separability and majorization give

```math
\widetilde M(x_i,y_j)+\widetilde M(x_j,y_i)
\le
\widetilde M(x_i,y_i)+\widetilde M(x_j,y_j).       \tag{VWE.1}
```

For `e= Mtilde-M`, subtracting the target four-cell contrast gives

```math
G_{ij}
\le
-e(x_i,y_j)-e(x_j,y_i)
+e(x_i,y_i)+e(x_j,y_j)
\le \sum_{z\in C_{ij}}|e(z)|.                    \tag{VWE.2}
```

Thus the direction of min-plus majorization and the four absolute-error
inequality in the draft are correct.  Positivity of `G_ij` also forces
`x_i ne x_j` and `y_i ne y_j`, so the four cells within one witness
rectangle are distinct.

If every cell in `C_ij` has positive mass, weighted Cauchy--Schwarz gives

```math
\left(\sum_{z\in C_{ij}}|e(z)|\right)^2
\le
\left(\sum_{z\in C_{ij}}\mu(z)e(z)^2\right)
\left(\sum_{z\in C_{ij}}{1\over\mu(z)}\right).
                                                               \tag{VWE.3}
```

Consequently the local bound

```math
\sum_{z\in C_{ij}}\mu(z)e(z)^2\ge w_{ij}          \tag{VWE.4}
```

has exactly the stated constant.  Cauchy--Schwarz is sharp as an inequality
under only the four-cell `l_1` constraint.

The convention `w_ij=0` when one witness cell has zero mass is also safe.
Equivalently, interpret `1/0=+infinity` in the denominator.  The four-cell
inequality alone allows all required error to sit on an uncharged cell, so
no positive universal local MSE certificate is available in that case.

## 3. The disjointness assertion is false

The draft claims that if `ij` and `kl` are disjoint as edges on anchor
indices, then `C_ij` and `C_kl` are disjoint as sets of matrix cells.  This
does not follow.  For example, anchors with indices in disjoint edges may
reuse a row:

```math
(x_1,y_1)=(1,1),\quad (x_2,y_2)=(2,2),
\qquad
(x_3,y_3)=(1,3),\quad (x_4,y_4)=(3,2).
```

For the disjoint graph edges `12` and `34`,

```math
C_{12}=\{1,2\}\times\{1,2\},
\qquad
C_{34}=\{1,3\}\times\{2,3\},
```

and both contain `(1,2)`.  Endpoint disjointness is therefore insufficient
for summing (VWE.4).

### Exact `3 by 3` counterexample

Let the anchors and graph be as above, with

```math
E=\{12,34\},
```

and put the uniform law on all nine cells.  Let

```math
M={1\over5}
\begin{pmatrix}
-1&2&-1\\
 1&-1&0\\
 0&-1&1
\end{pmatrix},
\qquad
\widetilde M=0.                                  \tag{VWE.5}
```

The zero matrix has min-plus factorization rank one.  Both declared gaps
are exactly one:

```math
M_{12}+M_{21}-M_{11}-M_{22}=1,
\qquad
M_{12}+M_{33}-M_{13}-M_{32}=1.                  \tag{VWE.6}
```

Hence each exposure weight is

```math
w_{12}=w_{34}={1\over 4\cdot9}={1\over36}.
```

For `k=1`, the unique coloring makes both edges monochromatic, and `E`
itself is an ordinary matching.  The parameter in the draft is therefore

```math
\mathfrak m_1={1\over18}.                         \tag{VWE.7}
```

On the other hand,

```math
\|M\|_F^2={2\over5},
\qquad
\mathbb E_{\rm unif}(\widetilde M-M)^2
={2\over45}<{1\over18}.                          \tag{VWE.8}
```

This contradicts (WE.4).  The shared cell `(1,2)` has the same positive
coefficient in both gap constraints, so its squared error is precisely what
the invalid summation counts twice.

There is no analogous nondegenerate `2 by 2` construction with four
distinct anchors split between two endpoint-disjoint edges: the only two
anchor pairs that span the whole rectangle are its two diagonals, and their
crossing gaps have opposite signs.  The `3 by 3` example is already enough
to falsify the theorem under its stated hypotheses.

## 4. Exact repairs

### 4.1 Integral witness packing

In (WE.3), replace “`F` is a matching in `E`” by

```math
C_e\cap C_{e'}=\varnothing
\quad\text{for all distinct }e,e'\in F.          \tag{VWE.9}
```

Then the written proof is valid verbatim.  This is a matching in the
four-uniform witness hypergraph on matrix cells, not generally a matching
in the anchor graph.

The original graph-matching definition is also valid under the additional
hypothesis that `i -> x_i` is injective, or that `i -> y_i` is injective.
For endpoint-disjoint graph edges, the corresponding row sets, respectively
column sets, are then disjoint, forcing the cell rectangles to be disjoint.

### 4.2 Fractional witness packing

A stronger and more natural repair retains overlapping witnesses.  For a
fixed tight-term coloring `c`, assign numbers `alpha_e>=0` to monochromatic
edges, subject to

```math
\sum_{e:z\in C_e}\alpha_e\le1
\quad\text{for every matrix cell }z.              \tag{VWE.10}
```

Multiplying (VWE.4) by `alpha_e` and summing gives

```math
\mathbb E_\mu e^2
\ge
\max_{\alpha\text{ satisfying }(\mathrm{VWE.10})}
\sum_{e\text{ monochromatic}}\alpha_e w_e.       \tag{VWE.11}
```

Taking the minimum over all `k`-colorings yields a correct weighted
fractional-packing theorem.  Integral cell-disjoint families are the special
case `alpha_e in {0,1}`.  This formulation records witness congestion
exactly at the level used by the proof.

## 5. The `D_r` calculation survives

For `D_r`, use diagonal anchors `(i,i)` and `E=K_r`.  Anchor rows and columns
are globally injective, so disjoint graph edges do have disjoint witness
rectangles.  Under the uniform law on `r^2` cells,

```math
G_{ij}=2,
\qquad
w_{ij}={4\over4r^2}={1\over r^2}.                 \tag{VWE.12}
```

If the color classes have sizes `s_1,...,s_q`, `q<=k`, the maximum
monochromatic matching has size

```math
\sum_{a=1}^q\left\lfloor{s_a\over2}\right\rfloor.
```

Minimizing this expression over colorings gives exactly

```math
\left\lceil{(r-k)_+\over2}\right\rceil:          \tag{VWE.13}
```

use singletons when `r<=k`; when `r>k`, `k-1` singleton classes and one
class of size `r-k+1` attain the displayed value, while at most one vertex
per nonempty class can remain unmatched.  Therefore (WE.9) is correct after
either repair.  For `k=1`, its order is `1/(2r)`, while the all-one
rank-one approximant has MSE `1/r`, exactly as claimed up to a factor two.

## 6. Max-plus needs the reversed contrast

Negating a max-plus factorization of `Mtilde` produces a min-plus
factorization of `-Mtilde`, but it also reverses the four-cell gap.  Thus the
max-plus counterpart requires

```math
M(x_i,y_i)+M(x_j,y_j)
-M(x_i,y_j)-M(x_j,y_i)>0,                         \tag{VWE.14}
```

with that reversed quantity used in the exposure weight.  The phrase “the
same statement holds ... after negating both matrices” is defensible only
if this change of hypothesis is understood.  It should be made explicit;
the positive crossed-minus-anchor gap in (WE.1) cannot simply be reused for
max-plus tight terms.

## 7. Research judgment

The corrected statement is not a tautology.  Its nontrivial link is that a
`k`-term min-plus representation induces a `k`-coloring of exposed anchors,
after which four-cell contrasts and a query law force MSE through a packing
certificate.  It answers the precise failure exhibited by the unweighted
`D_r` example and supplies a quantitative bridge from uniform tropical rank
to average response error.

It is nevertheless a **certificate reduction, not yet a generative model
theorem**:

- the witness system `(I,E)` is chosen externally;
- the packing value is defined by the same worst coloring that the proof
  must defeat;
- it need not characterize the true best rank-`k` average error;
- it may vanish even when a different witness geometry gives a strong
  lower bound.

Accordingly, the draft's phrase “the exact additional datum” is too strong.
The repaired theorem is substantive Level-2 machinery relative to the
repository's prior uniform-error result.  It becomes Level-3/generative only
if a natural model supplies a polynomially describable witness hypergraph
whose fractional packing value can be bounded at the target scale without
enumerating the full response table.

**Promotion recommendation:** do not promote WE.1 as written.  Replace its
matching parameter by the cell-disjoint or fractional witness-packing
version, state the reversed max-plus contrast explicitly, retain the `D_r`
application, and describe the result as a sufficient average-error
certificate rather than an exact characterization.

## 8. Addendum after the live draft was strengthened

The live theorem draft was edited independently before this audit was
delivered.  Its Section 2 now requires that the `x_i` are pairwise distinct
and that the `y_i` are pairwise distinct.  **Under that strengthened
hypothesis, Theorem WE.1 is correct.**  This addendum supersedes the opening
verdict as an assessment of the current live draft; Sections 3 and 4 above
remain an audit of why the injectivity hypothesis matters and of what would
be needed without it.

Indeed, if graph edges `ij` and `kl` have disjoint endpoint sets, pairwise
distinctness of the `x_i` gives

```math
\{x_i,x_j\}\cap\{x_k,x_l\}=\varnothing.
```

Since

```math
C_{ij}=\{x_i,x_j\}\times\{y_i,y_j\},
```

this already implies `C_ij cap C_kl` is empty.  Pairwise distinctness of
the `y_i` gives the same conclusion independently.  Thus either coordinate
injectivity would suffice; assuming both is stronger but harmless.  The
`3 by 3` counterexample in Section 3 repeats the row `x_1=x_3`, so it is
explicitly excluded.  With disjoint witness rectangles, summing (WE.6) is
valid, and the rest of the min-plus proof and all constants survive.

One wording issue remains in the max-plus sentence.  “After negating both
matrices” is mathematically correct if it means applying the entire theorem
to

```math
N=-M,\qquad \widetilde N=-\widetilde M,
```

and recomputing the gaps and `mathfrak m_k` for `N`.  Written directly in
terms of `M`, this requires the reversed contrast (VWE.14).  It does **not**
give a max-plus theorem with the original positive crossed-minus-anchor gap
or the unchanged quantity `mathfrak m_k(M,E,mu)`.  An explicit sentence to
that effect would remove the only remaining ambiguity.

**Corrected final verdict for the current draft:** the min-plus theorem is
rigorous as stated under the new pairwise-distinct anchor-coordinate
hypothesis; the `D_r` corollary is correct; no constant or zero-mass defect
remains.  The max-plus claim is correct under full negation and recomputation
but should spell out its reversed contrast.  The Level-2-versus-generative
research judgment in Section 7 is unchanged.
