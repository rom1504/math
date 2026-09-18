# Wave 41: global Johnson-cycle averaging

## Verdict

Global Johnson averaging does not manufacture the missing anchored-conflict
modulus.  A fixed selector-independent tie-break does improve the Wave 40
pairwise conclusion: distinct canonical responses have a **strict quantized**
gap, at least `8`.  Johnson Poincare then gives a valid averaged conflict
bound.  Its coefficient loses one factor for selector diffusion and another
for the possible Hamming size of one response jump, so it is far outside the
project scale.

Elementary four-cycles add no symmetric curvature.  After adding the reverse
orientation, a square is exactly the sum of its four adjacent pair gaps; after
summing all squares, this is only a fixed multiplicity of the adjacent-gap
sum.  The unsymmetrized part is a signed circulation in `A` and the response
labels.

A scalable common-penalty canonical map has corrected conflict `Omega(n)`,
parent row square `n`, bounded adjacent gaps, and bounded four-cycle regrets.
It is a complete-signing/common-penalty wall, not an exact-minimizer incidence
example.  Thus a theorem using exact-minimizer complement incidence remains
open, but it must add genuine quantitative normal-fan curvature or a sharper
signing-specific upper bound on the cyclic gap.

## 1. Adjacent swaps and strict canonicality

Fix an anchor `v`.  Put `N=n-1`, `k=m-1`, and `r=n-m`; the anchored slice is
`J(N,k)`.  For `T=S-a+b`, where `a in S\{v}`, `b notin S`, and
`C=S\{a}=T\{b}`, write the off-diagonal oriented-cut entries as `d^S_ij`.
The full pair gap is exactly

```math
\boxed{
g_{S,T}:=\langle B_S-B_T,d_S-d_T\rangle
=4\sum_{i\in C}\left[
a_{ai}(d^S_{ai}-d^T_{ai})
-a_{bi}(d^S_{bi}-d^T_{bi})
\right].}
\tag{R41.1}
```

Suppose `d_S` maximizes `<B_S,d>-P(d)` and ties are resolved by one fixed
total order on responses.  If `g_(S,T)=0`, both pairwise optimality
inequalities are equalities.  Hence `d_T` is also optimal at `S` and `d_S`
is also optimal at `T`.  A fixed total order cannot then select distinct
responses.  Therefore

```math
d_S\ne d_T\quad\Longrightarrow\quad g_{S,T}>0.
```

Because both `B_S-B_T` and `d_S-d_T` have even integral off-diagonal entries
and the matrix inner product counts both orientations,

```math
\boxed{d_S\ne d_T\quad\Longrightarrow\quad g_{S,T}\in8\mathbb Z_{>0}.}
\tag{R41.2}
```

This strictness uses the common tie-break; it does not apply to arbitrary
preselected grounds such as the flat Wave 40 pair.

## 2. What global adjacent averaging actually gives

For a uniform directed adjacent edge let

```math
h_{S,T}=\sum_{i\in(S\cap T)\setminus\{v\}}
\mathbf1\{x_i^S\ne x_i^T\}.
```

For each fixed `i ne v`, apply the exact spectral-gap Poincare inequality on
the conditional Johnson graph

```math
J(n-2,m-2),\qquad
\lambda_i=\frac{n-2}{(m-2)(n-m)}.
```

Converting its directed-edge normalization back to `J(n-1,m-1)` gives

```math
\boxed{
\mathcal C_v
\le\frac{(m-1)(n-m)}{n-2}
\mathbb E_{S\to T}h_{S,T}.}
\tag{R41.3}
```

Since an adjacent overlap has `m-2` nonanchor vertices, (R41.2) yields only

```math
\boxed{
\mathcal C_v
\le
\frac{(m-1)(m-2)(n-m)}{8(n-2)}
\mathbb E_{S\to T}g_{S,T}.}
\tag{R41.4}
```

At fixed selector density the prefactor in (R41.4) is `Theta(n^2)`.

There is also an exact signed formula for the mean gap.  For a response `d`
define

```math
e_S(d)=\sum_{i<j\in S}a_{ij}d_{ij},\qquad
\rho_j(d)=\sum_{i\ne j}a_{ij}d_{ij},\qquad
D_S(d)=\sum_{j\in S}\rho_j(d),
```

```math
I_{v,S}(d)=\sum_{i\in S\setminus\{v\}}a_{vi}d_{vi},
\qquad
O_{v,S}(d)=\sum_{b\notin S}a_{vb}d_{vb}.
```

Uniformly averaging the swap shore in (R41.1), and then using reversibility,
gives

```math
\boxed{
\mathbb E_{S\to T}g_{S,T}
=8\,\mathbb E_S
\frac{2(n-2)e_S(d_S)-(k-1)D_S(d_S)
-rI_{v,S}(d_S)-O_{v,S}(d_S)}{kr}.}
\tag{R41.5}
```

This is linear and signed, not a square.  If
`E R_2(d_S)<=R_*`, principal monotonicity gives
`|2e_S(d_S)|<=q_n`, while Cauchy--Schwarz bounds the degree sum.  Hence

```math
\boxed{
\mathbb E g_{S,T}
\le8\left[
\frac{(n-2)q_n}{kr}
+\frac{(k-1)\sqrt{mR_*}}{kr}
+1+\frac1k
\right].}
\tag{R41.6}
```

At fixed density and `R_*=O(n^(9/4-c))`, this is only
`O(n^(5/8-c/2)+n^(1/2))`.  Substitution in (R41.4) is much worse than the
trivial `C_v=O(n)` and therefore cannot prove the project conflict target.

## 3. Exact elementary four-cycle form

Let `C_0` have size `m-2` and contain the anchor.  For four distinct vertices
`a,b,c,d` outside `C_0`, use

```math
S_{00}=C_0\cup\{a,c\},\quad
S_{10}=C_0\cup\{b,c\},\quad
S_{11}=C_0\cup\{b,d\},\quad
S_{01}=C_0\cup\{a,d\}.
```

Write `z^pq_ij=d^(S_pq)_ij` and let `c_S` use the usual ordered quadratic
normalization.  The forward cyclic expression in (10.1104) is exactly twice

```math
\begin{aligned}
&\sum_{i\in C_0}\bigl[
a_{ai}(z^{01}_{ai}-z^{10}_{ai})
+a_{bi}(z^{10}_{bi}-z^{01}_{bi})\\
&\hspace{34mm}
+a_{ci}(z^{00}_{ci}-z^{11}_{ci})
+a_{di}(z^{11}_{di}-z^{00}_{di})\bigr]\\
&+a_{ac}(z^{00}_{ac}-z^{10}_{ac})
+a_{bc}(z^{10}_{bc}-z^{11}_{bc})\\
&+a_{bd}(z^{11}_{bd}-z^{01}_{bd})
+a_{ad}(z^{01}_{ad}-z^{00}_{ad}).
\end{aligned}
\tag{R41.7}
```

Every term is signed.  There is no `|H|`, squared label difference, or
unsigned cross-overlap energy.

If `F_square` is the full `B`-energy cyclic sum and `F_square^rev` its reverse,
then exactly

```math
\boxed{
F_\square+F_\square^{\rm rev}
=\sum_{\{S,T\}\in E(\square)}g_{S,T}.}
\tag{R41.8}
```

Every anchored Johnson edge belongs to

```math
(k-1)(r-1)=(m-2)(n-m-1)
```

elementary squares.  Consequently, summing (R41.8) over all squares is just
that multiplicity times the adjacent-gap sum.  The antisymmetric part is the
signed circulation (R41.7); the symmetric part contains no new curvature.

## 4. Exact finite canonical audit

For the exact order-five minimizer `A_5` in Wave 40, take `m=4`, anchor `0`,
`P=0`, and the same lexicographic order on all oriented projective cuts for
every selector.  The resulting canonical exact-complement-ground map has

```text
corrected anchored conflict = 1,
six undirected pair gaps     = (32,32,16,16,16,32),
overlap disagreements       = (0,0,2,0,2,2).
```

Thus even a genuinely common canonical rule on an exact minimizer need not
be coherent.  Exhaustion of every exact minimizer at orders `3` and `4` finds
no such conflict for this `m=n-1`, `P=0`, fixed-lex rule.  Order five is
minimal only in this precisely stated audit.  Unlike the arbitrary Wave 40
flat pair, all distinct canonical outputs here have the strict positive gap
required by (R41.2).

## 5. Scalable canonical wall

For every `h=1 mod 4`, let `n=2h`, `m=n-2`, split the vertices into two
`h`-sets `K,L`, and put the anchor in `L`.  Choose the `K x L` signs so that
every bipartite row and column sum is `+1` or `-1`, with the two signs nearly
balanced.  Inside `K` and `L`, use balanced circulant signings with every row
sum zero.

Let `x=1` and let `y` flip all of `K`.  Then exactly

```math
R_2(x)=R_2(y)=n.
```

For `delta_ij=a_ij(d^x_ij-d^y_ij)`, every signed `delta`-degree is `+2` or
`-2` and `U=sum_(i<j)delta_ij=2`.  Index an anchored selector by its two
omitted nonanchor vertices `O={o,p}`.  Use the fixed common penalty

```text
P(d_y)=0,  P(d_x)=4,  P(other response)=4n(n-1)+10,
```

and one fixed tie-break preferring `d_y`.  The exact objective difference is

```math
\boxed{
[H_{B_S}(d_x)-P(d_x)]-[H_{B_S}(d_y)-P(d_y)]
=-4(D_o+D_p)+4\delta_{op},}
\tag{R41.9}
```

where `D_o,D_p in {+2,-2}`.  Its possible values are
`{-24,-16,-8,0,8,16,24}`, so this defines a canonical map using only `x,y`.

There are `Theta(h^2)` selectors of each response type.  Independently drawn
selectors of opposite types disagree on at least `h-4` common vertices, so

```math
\mathcal C_v=\Omega(h)=\Omega(n).
```

Nevertheless every adjacent cyclic gap is at most `48`, and every directed
elementary-four-cycle regret is at most `96`.  Thus their normalized global
averages are `O(1)`.  At `h=5` the exact checker obtains

```text
n=10, row=10,
corrected conflict=391/252,
mean adjacent gap=44/7,
mean adjacent overlap disagreement=74/63.
```

This wall uses one complete signing and one selector-independent penalty; it
is not an arbitrary assignment of grounds.  It has extremely low parent row,
but the signing is not asserted to be an exact minimizer and the chosen
responses do not satisfy the complement-incidence threshold.  It therefore
rules out a theorem from canonicality, row, and Johnson-cycle algebra alone,
not an exact-minimizer incidence theorem.

## 6. Sharp surviving hypothesis

Let `G_*` be an available upper budget for the mean adjacent Bregman gap.  An
exact sufficient normal-fan hypothesis is

```math
\mathbb E_{S\to T}g_{S,T}
\ge\kappa_n\mathbb E_{S\to T}h_{S,T},
\qquad
\frac{(m-1)(n-m)}{n-2}\frac{G_*}{\kappa_n}
\le C_{\rm conflict}(n).
\tag{R41.10}
```

Together with `E g<=G_*`, (R41.3) proves precisely the conflict clause of
(10.1089); complement incidence and the parent-row/information clauses can
then be supplied independently.  One may take the generic row upper budget
from (R41.6), but at project parameters it is not affordable.  Therefore the
useful successor must prove either:

1. a much stronger minimizer-specific cross-overlap curvature and a matching
   small gap budget; or
2. direct edge coherence
   `E h <= [(n-2)/((m-1)(n-m))] C_conflict(n)`
   for the incidence-preserving low-row canonical map.

Ordinary strict `8Z` canonicality, all elementary cycle inequalities, and the
known row budget do not approach (R41.10).

## Verification

`tmp/johnson_cycle_r41_check.py` verifies (R41.1), strict gap quantization,
the exact order-five canonical map and its minimal scoped audit, the common
penalty wall (including exhaustive response maximization at `n=10`), all
directed elementary-square identities, and the displayed rational values.
It ends with

```text
PASS johnson_cycle_r41_check
```
