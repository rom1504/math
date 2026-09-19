# Exact weighted witness and a `\{\pm1\}` blow-up

## 1. The finite LP

For a complete weighted graph on `S\sqcup T`, write `\delta_R` for the
cut-incidence vector, `\tau=\delta_T`, `e_X` for the incidence vector of the
internal edges of `X`, and `r_X` for a prescribed child-minimum cut.  The
weighted same-positive relaxation, normalized by `c=1`, is

```math
\begin{aligned}
\max\quad &(e_S+e_T-r_S-r_T)\cdot w,\\
\text{s.t.}\quad
&\tau\cdot w=1,\\
&0\le \delta_R\cdot w\le1 &&(R\subseteq S\sqcup T),\\
&r_X\cdot w\le\delta_W\cdot w &&(W\subseteq X),\\
&\delta_W\cdot w\le(e_X-r_X)\cdot w &&(W\subseteq X),\\
&(e_S+e_T)\cdot w\ge0.
\end{aligned}
```

The fourth line says that the prescribed cut is a child minimum.  The fifth
is exactly positive dominance `M_X\le a_X-m_X`.

## 2. An exact `5+4` primal witness

Partition

```math
S=U\sqcup V,\quad |U|=3,\ |V|=2,
\qquad
T=C\sqcup D,\quad |C|=3,\ |D|=1.
```

Give every edge between a fixed pair of cells the following weight (the
diagonal entries mean distinct vertices in that cell):

| pair | weight | pair | weight |
|---|---:|---|---:|
| `UU` | `1/8` | `UV` | `-1/12` |
| `VV` | `1/6` | `CC` | `1/8` |
| `CD` | `1/24` | `UC` | `1/24` |
| `UD` | `1/24` | `VC` | `1/12` |
| `VD` | `0` | | |

If a cut selects `u,v,r,d` vertices from `U,V,C,D`, respectively, its weight
is

```math
\begin{aligned}
f_0(u,v,r,d)
={}&\frac{u(3-u)+r(3-r)}8+\frac{v(2-v)}6\\
&-\frac{u(2-v)+(3-u)v}{12}
+\frac{r(1-d)+(3-r)d}{24}\\
&+\frac{u(4-r-d)+(3-u)(r+d)}{24}
+\frac{v(3-r)+(2-v)r}{12}.
\end{aligned}
```

Checking the `4\cdot3\cdot4\cdot2=96` count tuples gives `0\le f_0\le1`.
Equivalently, the projective cut-value histogram is

```text
value:  0  1/6  1/4  1/3  5/12  7/12  2/3  3/4  5/6  11/12  1
count:  4    2   14    6     8     18   22   48   30      24  80.
```

The endpoint cut is `S|T`, and its weight is

```math
c=12\frac1{24}+6\frac1{12}=1.
```

On `S`, take the cut `U|V`; on `T`, take the empty cut.  Exact enumeration of
the `16` and `8` projective child cuts gives

```math
(a_S,m_S,M_S)=\left(\frac1{24},-\frac12,\frac16\right),
\qquad
(a_T,m_T,M_T)=\left(\frac12,0,\frac13\right).
```

Hence

```math
M_S\le a_S-m_S=\frac{13}{24},
\qquad
M_T\le a_T-m_T=\frac12,
```

so both children are positive-dominant, but

```math
a_S+a_T-m_S-m_T=\frac{25}{24}>c.
```

Both `a_S` and `a_T` are positive, so this also violates the weaker stopping
inequality (10.435), not only the sufficient inequality (10.437).

## 3. Strictification

Let `G` have weight `1/20` on every one of the twenty `S-T` edges and zero
internally, and put

```math
\bar W=\frac{49}{50}W_0+\frac1{50}G.
```

Both summands have all cuts in `[0,1]`.  For `G`, only the empty/full cuts
have value zero and only `S,T` have value one.  Thus `\bar W` has projectively
unique parent endpoints.  Direct exact enumeration gives lower and upper
non-endpoint gaps `1/125` and `1/200`.  Its child data are

```math
(a_S,m_S,M_S)=
\left(\frac{49}{1200},-\frac{49}{100},\frac{49}{300}\right),
```

```math
(a_T,m_T,M_T)=
\left(\frac{49}{100},0,\frac{49}{150}\right).
```

The child-minimum gaps are `49/120` and `49/400`, and the positive-dominance
gaps are

```math
a_S-m_S-M_S=\frac{147}{400},
\qquad
a_T-m_T-M_T=\frac{49}{300}.
```

The violation remains strict:

```math
a_S+a_T-m_S-m_T-c=\frac1{48}.
```

For later use, the distinct strictified edge weights are

```math
\begin{array}{c|ccccccccc}
\text{pair}&UU&UV&VV&CC&CD&UC&UD&VC&VD\\ \hline
\bar w&49/400&-49/600&49/300&49/400&49/1200&
251/6000&251/6000&31/375&1/1000.
\end{array}
```

In particular,

```math
\max_{i<j}|\bar w_{ij}|=\frac{49}{300}<1.
```

No rescaling is needed before using these numbers as sign biases.

In the ledger normalization,

```math
h_S=\frac{49}{600},\quad h_T=\frac{49}{50},\quad
H=\frac{637}{600},
```

```math
Q_S=\frac{49}{24},\quad Q_T=\frac{49}{50},\quad
b_S=\frac1{25},\quad b_T=2.
```

Consequently

```math
\mathcal B=\frac{51}{25}
<\frac{637}{300}=I,
\qquad
I-\mathcal B=\frac1{12}.
```

## 4. Lifting to complete sign matrices

For each of the nine weighted vertices make a clone class of size `L`.  For
clones in two different classes `i,j`, independently choose an edge sign
`A_{uv}\in\{\pm1\}` with

```math
\mathbb E A_{uv}=\bar w_{ij}.
```

This is the Bernoulli law

```math
\mathbb P(A_{uv}=1)=\frac{1+\bar w_{ij}}2,
\qquad
\mathbb P(A_{uv}=-1)=\frac{1-\bar w_{ij}}2,
```

which is admissible by `\max|\bar w_{ij}|=49/300<1`.  Within a clone class
use independent unbiased signs.  Every unordered pair of distinct clones is
assigned exactly one sign, so this is a complete symmetric zero-diagonal
`\{\pm1\}` matrix, not a multigraph or a graph with missing edges.

For a clone cut `R`, let `\alpha_i` be the fraction selected from class `i`.
If `F` is the multilinear extension of the base cut function, then exactly

```math
\mathbb E\,d_A(R)=L^2F(\alpha),
\qquad
F(\alpha)=\mathbb E\,f(B_\alpha),
```

where `B_\alpha` independently selects base vertex `i` with probability
`\alpha_i`.  Since the only zero vertices of `f` are empty/full and the only
vertices with value `c` are `S,T`, compactness and the positive one-coordinate
derivatives at those vertices give constants `\gamma_0,\gamma_1>0` such that

```math
F(\alpha)\ge\gamma_0\,d_1(\alpha,\{\mathbf0,\mathbf1\}),
```

```math
c-F(\alpha)\ge
\gamma_1\,d_1(\alpha,\{\mathbf1_S,\mathbf1_T\}).
```

For completeness, the compactness assertion is quantitative.  Away from the
two zero vertices, `F` is strictly positive: as an expectation of nonnegative
Boolean cut values it could vanish only if the independent Bernoulli support
were contained in the two complementary zero vertices, and no positive-
dimensional product subcube has that property.  Hence the ratio of `F` to
their `\ell^1` distance has a positive minimum on every compact set away from
them.  Near either zero vertex `e`, multilinearity gives

```math
F(e+t)=\sum_i\bigl(f(e\mathbin\triangle\{i\})-f(e)\bigr)|t_i|
+O(\lVert t\rVert_1^2),
```

and every displayed one-coordinate coefficient is positive because the
endpoint is projectively unique.  This proves the claimed linear lower bound.
The same argument applies to `c-F` at the two maximum endpoints.

The child cut functions have projectively unique prescribed minima: on `S`
the minimum is the clone union corresponding to `U|V`, and on `T` it is the
empty/full pair.  Their next discrete cut gaps are, respectively, `49/120`
and `49/400`.  Thus their multilinear extensions satisfy analogous linear
lower bounds in the distance from the prescribed minimum/complement pair.

Suppose a clone cut is at Hamming distance `k` from the relevant endpoint.
Its expected lower or upper gap is at least `\gamma kL`.  The difference from
the endpoint cut uses at most `k(9L-k)` independent signed edges.  Hoeffding's
inequality therefore gives a failure probability

```math
\mathbb P(\text{failure for this cut})
\le \exp\left(-\frac{\gamma^2kL}{18}\right).
```

There are at most `2\binom{9L}{k}` cuts at distance `k` from either member of
an endpoint pair.  Thus

```math
\sum_{k=1}^{\lfloor9L/2\rfloor}
2\binom{9L}{k}\exp\left(-\frac{\gamma^2kL}{18}\right)=o(1),
```

because `\binom{9L}{k}\le(e9L/k)^k` while the negative exponent is linear in
`kL`.  This is a projective shell count: distance is taken from either member
of the complementary endpoint pair, so every cut is counted.

With probability tending to one, every parent cut is therefore in
`[0,c_L]`, and both prescribed child cuts are actual child minima.

The last assertion is exact, not merely an approximation statement.  For a
child cut `R` and prescribed child cut `R_*`, put
`k=\min\{|R\triangle R_*|,|R\triangle R_*^c|\}`.  The random difference
`d_A(R)-d_A(R_*)` uses only the at most `k(n_XL-k)` edges whose cut status
changes, and its expectation is at least `\gamma_XkL`.  The same shell union
bound proves simultaneously that it is positive for every
`R\notin\{R_*,R_*^c\}`.  Hence the actual minimum equals
`m_{X,L}=d_A(R_*)` projectively uniquely.  No comparison of expectations is
being substituted for the assertion about the random extrema.

For child positive dominance, the base slack is uniformly at least
`\eta_S=147/400` or `\eta_T=49/300` at every Boolean cut.  Its multilinear
extension has the same lower bound.  Each fixed clone cut violates with
probability `\exp(-\Omega(L^2))`; a union bound over all child cuts is still
`o(1)`.  Explicitly, this controls the random inequalities

```math
d_A(R)\le a_{X,L}-d_A(R_*)=a_{X,L}-m_{X,L}
\qquad (R\subseteq X_L)
```

simultaneously, including the random internal total `a_{X,L}`.  Hence both
actual children are positive-dominant.

Finally, the random variable

```math
a_{S,L}+a_{T,L}-m_{S,L}-m_{T,L}-c_L
```

has mean `L^2/48` and standard deviation `O(L)`.  The two internal totals are
also positive with probability tending to one.  The fixed parent endpoint
energies have expectations `2L^2(a+c)>0` and `2L^2(a-c)<0`; their signs are
likewise stable with probability tending to one.

On the simultaneous event just proved, put `a=a_{S,L}+a_{T,L}` and
`m=m_{S,L}+m_{T,L}`.  Because each `a_{X,L}>0` and each child is
positive-dominant,

```math
Q(X)=2(a_{X,L}-2m_{X,L}),
\qquad |h_X|=2a_{X,L}.
```

Using the exact exposure formula (10.433),

```math
b_X=2c_L-Q(X)+|h_X|=2c_L+4m_{X,L}.
```

Therefore

```math
\mathcal B=4c_L+4m,
\qquad I=4a,
\qquad I-\mathcal B=4(a-m-c_L)>0.
```

Thus, for every sufficiently large `L`, there exists an actual complete
`\{\pm1\}` matrix of order `9L` whose endpoint pair and positive-dominant
children violate (10.440) and (10.435) exactly.

The finite weighted witness lacks equal edge magnitudes, but this is not a
property that can save the proposed universal inequality: biased random
blow-up restores exact `\{\pm1\}` entries while preserving the strict
violation.
