# Nine-type obstruction: limiting allocation/path minimax

All numbers below use the integer matrix `M=672W` from (10.448).  Scaling
all capacities and obligations by `1/672` (or by the blow-up factor
`kappa L^2/672`) leaves `K_*` unchanged.

## 1. Discrete macro endpoint tree

Exact enumeration gives

```text
root: P=2056, N=632, R=2688, d=1424; 1 positive and 1 negative endpoint
S:    P= 648, N=216,        d= 432; 1 positive and 3 negative endpoints
T:    P=1384, N=392,        d= 992; 1 positive and 6 negative endpoints
```

All endpoint choices give the same LP signature.  With children unordered,
the full signature is

```text
root d=1424
  T edge caps (0,24), child d=992
    edge caps (116,888), child F d=220
      edge caps (276,276), leaf
      edge caps (56,276), child d=0
        two edge caps (110,110), leaves
    edge caps (668,888), child d=0
      two edge caps (110,110), leaves
  S edge caps (48,1344), child d=432
    two edge caps (216,432), child d=0
      two edge caps (108,108), leaves
```

Here `F` is the three-type child with weights `83,-55,-55`.

The exact optimum for root obligation `W=R=2688` is

```math
K_*=\frac{175}{99}=1+\frac{76}{99}.
```

An exact primal certificate is as follows.  At the root, allocate `608/33`
to the capacity-24 bucket on the `T` edge and pass `29488/33` to `T`;
allocate `1344` to the capacity-1344 bucket on the `S` edge and pass `432`
to `S`.  At `T`, use the branch leading to `F`: allocate `22496/33` into
its capacity-888 bucket and pass `6992/33` to `F`.  At `F`, allocate
`6992/33` into a capacity-276 bucket.  At `S`, allocate all `432` into
one capacity-432 bucket.  All other variables vanish.  Every active `T`
load is `76/99`, every active `S` load is `1`, and conservation and all
imbalance caps hold exactly.

For the matching dual, use the following edge-antichain weights:

```text
q(root,T)=2/99,       q(root,S)=7/11,
q(T,F)=74/99,         q(T,other)=74/99,
q(F,child)=23/99      on both F children,
q(S,child)=4/11       on both S children,
all deeper q=0.
```

Put `y=1/1188` on the root, `T`, `F`, the directly used leaf below `F`,
and `S`, and put all other `y` equal to zero.  Put bucket slack

```math
z_{\mathrm{root},S,\mathrm{high}}=\frac7{19008},
```

and zero on every other bucket.  Put obligation slack `1/1188` on the
four zero-imbalance intermediate edge nodes (their objective cost is zero),
and zero elsewhere.  Root-to-leaf `q` sums are at most one, and the active
bucket inequalities are equalities.  The dual objective is

```math
2688\frac1{1188}-1344\frac7{19008}
=\frac{224}{99}-\frac{49}{99}=\frac{175}{99}.
```

The scratch implementation independently reconstructs and exactly checks
both certificates.

## 2. Flat endpoint faces in dense blow-ups

The six negative endpoints of `T` are not isolated in the graphon cube.
They induce a six-cycle of cube edges: vertices `7,8` have opposite signs,
two of the three `B={4,5,6}` signs are opposite, and the remaining `B`
coordinate can be any value in `[-1,1]`.  These are all continuous endpoint
ties.  Indeed,
the multilinear extension is a convex combination of vertex values; an
extremizer can have fractional coordinates only if every vertex in the
corresponding cube face is a discrete extremizer.  The six discrete minima
contain exactly these six edges and no higher-dimensional face.

Let `t` be the fraction of the split clone class sent to one child and set
`u=min(t,1-t) in [0,1/2]`.  The two children of `T` are three-type systems
`F_s`, for `s=u,1-u`, with

```text
d(F_s)=220s,
T-to-F_s caps=(668-552s,888),
F_s child caps=(110+166s,110+166s) and (110-54s,110+166s),
last two-type caps=(110s,110s).
```

Consequently the maximum obligation that `F_s` can accept with congestion
at most `z` is

```math
H_s(z)=\min\{220s,(110+166s)z\}.
```

After the `S` side is filled at its more efficient marginal rate, it serves
exactly `1776` units at cost one.  The remaining `912` units use `T`.  If
`theta` is the `T` contribution, exact elimination of the allocation LP gives

```math
888\theta+
\max_{z_1+z_2=\theta}
\bigl(H_u(z_1)+H_{1-u}(z_2)\bigr)
=912-24\theta,
\qquad K(t)=1+\theta.
```

The larger slope is used first.  If

```math
u_0=\frac{18311-\sqrt{313271161}}{9130}=0.0669805\ldots,
```

then

```math
\theta(u)=
\begin{cases}
\displaystyle\frac{912}{1188-166u},&0\le u\le u_0,\\[6pt]
\displaystyle\frac{215192-41832u-73040u^2}
{(276-166u)(1022+166u)},&u_0\le u\le\frac12.
\end{cases}
```

Both pieces are increasing.  For the second piece, the numerator of the
derivative after clearing the positive squared denominator is

```math
14848880608-29345416256u+7892258848u^2>0
\quad(0\le u\le1/2).
```

Therefore every possible flat-face endpoint selection satisfies the exact
range

```math
\boxed{\frac{175}{99}\le K(t)\le\frac{2017}{1105}}
\qquad
(1.76768\ldots\le K(t)\le1.82534\ldots).
```

Exact rational primal and dual certificates were checked at both endpoints.
Thus neither macro endpoint ties nor a fractional clone-class tie approaches
the previous `10/3` lower bound, let alone four.

## 3. Alternate strictification

For the alternate `49/50` strictification, the correct denominator-6000
weight on the internal `C-D` edges is `245`, not `251`; the `1/1000`
addition applies only to the parent cross block.  With this correction,
exact primal and dual certificates give

```math
K_*^{\rm vertex}=\frac{1861}{1053}=1.767331\ldots.
```

For clarity, the exact dual is as follows.  On the root edges to the flat
five-type child and the other four-type child, respectively, put
`q=8/351,661/1053`.  Put `q=784/1053` on both edges below the flat child,
`q=245/1053` on both edges below its active three-type child, and
`q=392/1053` on both edges below the four-type child.  All deeper `q` vanish.
Put `y=1/10530` on the root, the flat child, its active three-type child and
directly used leaf, and the four-type child.  The sole positive-cost bucket
slack is

```math
z_{\mathrm{root},\,\mathrm{four\mbox{-}type},\,\mathrm{high}}
=\frac{539}{12636000}.
```

Put obligation slack `1/10530` on the four intervening zero-imbalance edge
nodes.  The longest path sums are
`8/351+784/1053+245/1053=1` and
`661/1053+392/1053=1`.  The dual objective is exactly

```math
24000\frac1{10530}-12000\frac{539}{12636000}
=\frac{2400-539}{1053}=\frac{1861}{1053}.
```

The generic checker verifies this certificate when its rationalization
denominator is raised from the default `10^7` to at least `12,636,000`.

Its analogous flat-face parameter is again monotone and has maximum

```math
K_*^{\rm half\ split}=\frac{3575}{1959}=1.824911\ldots.
```

This half-split value also has fully exact certificates.  The standalone
checker is `tmp/verify_alt_halfsplit_kmin_r10.py`; it hardcodes the complete
rational tree, checks every conservation, capacity, path-load, antichain,
bucket-dual, and obligation-dual inequality using `fractions.Fraction`, and
checks equality of the two objectives.  Its node labels are: `0` root, `1`
the flat five-type child, `2,7` its two half-split three-type children, and
`12` the other four-type child.  Nodes `3,8` are directly used leaves and
`4,9,13,16` are intervening zero-imbalance edge nodes.  Branch and orientation
indices below agree with the `TREE` dictionary in the checker.

All nonzero primal variables are

```text
a(0,0,1)=258560/653       a(0,1,1)=24000
a(1,0,1)=8960             a(1,1,1)=7786240/1959
a(2,0,0)=1960             a(7,0,0)=1703240/1959
a(12,0,1)=7840

w(1)=10293920/653         w(2)=1960
w(7)=1703240/1959         w(12)=7840

theta(0)=3575/1959        theta(1)=1616/1959
theta(2)=4/7              theta(7)=3476/13713
theta(12)=1

path-z(0,0)=1616/1959     path-z(0,1)=1
path-z(1,0)=4/7           path-z(1,1)=3476/13713
path-z(2,0)=4/7           path-z(7,0)=3476/13713
path-z(12,0)=1.
```

All nonzero dual variables are

```text
q(0,0)=16/653             q(0,1)=1175/1959
q(1,0)=q(1,1)=1568/1959
q(2,0)=q(2,1)=343/1959    q(7,0)=q(7,1)=343/1959
q(12,0)=q(12,1)=784/1959

y(0)=y(1)=y(2)=y(3)=y(7)=y(8)=y(12)=1/19590

bucket-slack(0,1,1)=49/1880640
s(4)=s(9)=s(13)=s(16)=1/19590.
```

All omitted variables are zero.  The longest dual paths have weights

```math
\frac{16}{653}+\frac{1568}{1959}+\frac{343}{1959}=1,
\qquad
\frac{1175}{1959}+\frac{784}{1959}=1.
```

The sole positive-cost dual slack gives objective

```math
48000\frac1{19590}-24000\frac{49}{1880640}
=\frac{4800-1225}{1959}=\frac{3575}{1959},
```

matching the primal root value exactly.

## 4. Recursive substitution (exploratory only)

A formal leaf-graft test replaces each of the nine zero-scale macro leaves
by a copy scaled by `1/81`, while retaining the dominant outer range.  The
optimized values at depths one, two, and three are

```text
175/99 = 1.767677,   401/228 = 1.758772,   401/228 = 1.758772.
```

A more detailed block-uniform hierarchy, including the induced changes in
child imbalances and capacities, gives approximately

```text
1.767677, 1.592514, 1.590581.
```

The second experiment assumes the hierarchical endpoints remain block
uniform; that is automatic only after choosing the inner scale sufficiently
small compared with the outer strict endpoint margins.  It is therefore an
exploratory diagnostic, not a realizability theorem.  In either model,
recursion decreases rather than amplifies congestion.  The rigorous result
from this memo is the one-level exact certificate and the full flat-face
range above.
