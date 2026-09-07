# A hard-constraint compiler for fixed symmetric finite alphabets

Date: 2026-09-07. Status: **Pending independent audit**. This is a positive
construction proposed after the director's signed Euler transition count.
It uses that count, hence inherits its ordinary Eulerian-orientation lower
bound dependency. It is not an original-problem convergence theorem.

## 1. Target theorem

Fix a symmetric finite real alphabet
`E={0,+-a_1,...,+-a_q}`, with distinct positive a_r; the zero atom may be
omitted. Fix a symmetric law nu with positive masses on its support. Let
m be odd, d=m-1, and assume every `d nu(a)` is integral. At each vertex,
independently uniformly permute the d-symbol multiset with this type. For
any symmetric seed signing S, put

```math
Z_S(t)=\mathbb E\prod_{i<j}
\exp[-t(u_{ij}-S_{ij}u_{ji})^2],\qquad e=md/2.
```

The proposed uniform construction proves

```math
\log Z_S(t)\ge-eH(\nu)
-C_{E,\nu}m^{3/2}(\log m)^{3/2}
-4t\sum_{r=1}^q a_r^2.
\tag{1}
```

In particular, for any `t_m -> infinity` with `t_m=o(m²)`,

```math
\sup_S\left|m^{-2}\log Z_S(t_m)+\tfrac12H(\nu)\right|\longrightarrow0.
\tag{2}
```

Thus fixed symmetric finite alphabets retain no leading seed-pressure gain
at the critical t~sqrt(m) scale. The construction stays inside magnitude
matching; its energetic payment is at most ONE parity defect per nonzero
magnitude, not the number of row-count corrections.

## 2. A dense multicolor degree compiler

Let the magnitude colors have probabilities
`p_r=nu(+a_r)+nu(-a_r)` for r>0 and `p_0=nu(0)` if present. Every positive
p_r is fixed. The desired color-r degree is `d_r=d p_r`. For r>0 this
equals `2d nu(a_r)` and is even. Since d is even, the zero-color degree is
also even. Choose any positive-probability color as a buffer b.

Initially color each undirected edge independently with law p. With
probability tending to one, the following simultaneous properties hold:

1. Every color degree differs from its target by at most
   `D=C_p sqrt(m log m)`.
2. For every two distinct vertices u,v and every ordered pair of positive-
   probability colors r,s, there are at least c_p m vertices w with
   `color(uw)=r` and `color(vw)=s`.

Both claims follow from binomial Hoeffding/Chernoff bounds and a union
bound over O(m²) fixed-color tests. They include r=s, which later ensures
connectivity. Fix a deterministic repair algorithm as follows.

For each nonbuffer color r in turn, let delta_v be its target degree minus
its current degree. Earlier stages change only their own color and b, so
the initial errors for r are still the original errors, bounded by D.

- If delta_u>0 and delta_v<0, choose w with uw colored b and vw colored r.
  Exchange those two colors. Only u and v change their r-degrees, by +1
  and -1 respectively.
- If all nonzero deltas are positive, choose two deficit units at vertices
  u,v, allowing u=v if its deficit is at least two. Find an alternating
  path u-a-b'-v with colors b,r,b and interchange the two colors along
  the path. This increases the r-degrees at u,v and preserves the two
  internal degrees. If u=v, use the corresponding three-edge triangle.
- If all nonzero deltas are negative, use the reversed pattern r,b,r.

The total degree discrepancy has even sum, so the same-sign step always
has two units available. Every step decreases `sum_v |delta_v|` by two,
and uses at most three edge recolorings. Thus each color uses at most
`(3/2)mD` recolorings, and all stages use

```math
R=O_p(m^{3/2}\sqrt{\log m}).
\tag{3}
```

## 3. Why the alternating paths remain available

Require every internal path vertex to have accumulated at most
`m^(3/4)/2` internal incidence changes. The total number of incidence
changes is O(R), so at any time only
`O_p(m^(3/4)sqrt(log m))=o(m)` vertices are ineligible. The endpoint changes
at a fixed vertex are at most O_p(D), because each color's discrepancy
is reduced monotonically and internal steps preserve that color's degree.
Consequently the total number of changed edges incident to any vertex is
at most `m^(3/4)` for sufficiently large m.

The initial two-color common-neighbor count loses at most twice that
number, and hence remains at least `(c_p/2)m`. Two-edge transfers therefore
have an eligible mediator. For a three-edge path, choose an eligible a in
the required color neighborhood of u; this neighborhood has linear size.
Then choose b' in the required two-color common neighborhood of a and v.
Remove the o(m) ineligible vertices and the bounded list of forbidden
repetitions. A choice remains. This also handles the triangle case u=v.

The repair is therefore well-defined. Once every nonbuffer color has its
target degree, the buffer has its target automatically. Same-color common
neighborhoods remain linear, so every positive-probability color graph is
connected. Each graph is regular of its prescribed even degree.

## 4. Counting exact colorings

Let P_p be the independent edge-color law. Restrict to the preceding good
event, which has probability at least one half for all sufficiently large
m. Repair changes at most R edges. Put
`rho_p=min_r p_r/max_r p_r>0`. If u repairs to v, then
`P_p(v)>=rho_p^R P_p(u)`. Each v has at most

```math
\sum_{j\le R}\binom ej(q+1)^j
=\exp[O_p(m^{3/2}(\log m)^{3/2})]
```

preimages, where q+1 can be replaced by the actual number of colors.
Thus the probability of repaired exact colorings is at least the inverse
of that subexponential factor. Every exact coloring has e p_r edges of
color r, so every such coloring has probability exactly exp[-eH(p)].
There are consequently at least

```math
\exp[eH(p)-O_p(m^{3/2}(\log m)^{3/2})]
\tag{4}
```

connected regular color decompositions with the desired degrees.

## 5. Signed Euler assignments on every nonzero color

For a connected color-r graph of even regular degree 2k_r, the director's
transition-system argument and the ordinary Eulerian-orientation lower
bound give, whenever its signed global parity is compatible,

```math
\#\{\text{balanced signed assignments on }G_r\}
\ge\left[\frac{\binom{2k_r}{k_r}}{2^{k_r}k_r}\right]^m
=2^{e_r}\exp[-O_p(m\log m)].
\tag{5}
```

If the parity fails, change the required seed sign on one selected edge
of G_r. The same count then produces assignments satisfying the original
seed everywhere except that one edge. Assign values +-a_r from these
incidence signs. Exactly half the color-r incidences at every vertex have
each sign, as required by nu. Zero-color incidences receive zero and impose
no sign constraint.

Colorings are recovered from the magnitudes of the resulting incidences,
so counts from distinct magnitude colorings do not overlap. Multiplying
(4) and (5) over r>0 gives at least

```math
\exp[eH(p)+e(1-p_0)\log2
-O_{E,\nu}(m^{3/2}(\log m)^{3/2})]
=\exp[eH(\nu)-O_{E,\nu}(m^{3/2}(\log m)^{3/2})]
```

actual row-type assignments. They satisfy every hard Gaussian equality
except at most one edge per nonzero magnitude. Their total defect is at
most `4 sum_r a_r²`. The number of all independent row-type assignments
has logarithm `md H(nu)+O_E(m log m)=2eH(nu)+O_E(m log m)`.
Dividing the weighted count by this number proves (1).

## 6. Upper comparison for the growing-temperature statement

Let Delta>0 be the minimum gap between distinct symbols of E. For a
coupling (X,Y) with both marginals nu, set delta=Pr(X!=Y). Its transport
objective satisfies

```math
\mathbb E\log K_t(X,Y)-D(\pi\Vert\nu\otimes\nu)
\le-H(\nu)+h(\delta)+\delta\log(|E|-1)-t\Delta^2\delta.
```

Optimizing the final three terms yields

```math
\mathcal V_{K_t}(\nu)
\le-H(\nu)+\log(1+(|E|-1)e^{-t\Delta^2}).
```

The diagonal coupling gives the opposite lower bound `V>=-H(nu)`.
The exact symmetric-type canonical upper identity is valid at every t;
it pays only `-log P_0(exact types)=O_E(m log m)` in the upper direction.
Together with (1), this proves (2).

## 7. Audit obligations

The main new obligation is the load-balanced alternating-path degree
compiler in Sections 2--3. Every operation must preserve internal degrees
and all previously fixed colors. The signed Euler lower bound on arbitrary
connected even regular color graphs is a dependency, not reproved here.
All constants may depend on the fixed alphabet law; no uniform growing-
alphabet theorem is claimed.
