# Exact weighted obstruction and a probabilistic signing blowup

This is a scratch report for Wave 9.  It uses the ledger's undoubled edge
cut normalization.  Thus, for symmetric off-diagonal edge weights `w_ij`,

```math
g(R)=\sum_{i\in R,\,j\notin R}w_{ij},
```

and changing a spin from `\mathbf 1` across `R` lowers its quadratic energy
by `4g(R)`.  The eventual random blowup has honest off-diagonal signs in
`\{\pm1\}`.

## 1. An exact strict `4+5` weighted certificate

Let

```math
S=\{0,1,2,3\},\qquad T=\{4,5,6,7,8\},
\qquad B=\{4,5,6\}.
```

Take `w_{ij}=M_{ij}/672`, where `M` is symmetric with zero diagonal and its
upper-triangular entries are specified by

```math
\begin{array}{c|c}
\text{edge class}&M_{ij}\\ \hline
S\text{--}S&54\\
S\text{--}B&28\\
S\text{--}\{7\}&55\\
S\text{--}\{8\}&29\\
B\text{--}B&83\\
B\text{--}\{7,8\}&-55\\
\{7,8\}&113.
\end{array}
```

Equivalently, the full scaled matrix is

```math
M=\begin{pmatrix}
0&54&54&54&28&28&28&55&29\\
54&0&54&54&28&28&28&55&29\\
54&54&0&54&28&28&28&55&29\\
54&54&54&0&28&28&28&55&29\\
28&28&28&28&0&83&83&-55&-55\\
28&28&28&28&83&0&83&-55&-55\\
28&28&28&28&83&83&0&-55&-55\\
55&55&55&55&-55&-55&-55&0&113\\
29&29&29&29&-55&-55&-55&113&0
\end{pmatrix}.
```

Exact enumeration of all `2^9` cuts gives

```math
g(\varnothing)=g([9])=0,
\qquad g(S)=g(T)=1,
```

and, for every other cut,

```math
\frac1{112}=\frac6{672}
\le g(R)\le
\frac{666}{672}=\frac{111}{112}.
\tag{A}
```

There are no other zero or unit cuts.  Thus `p=\mathbf1` is the unique
projective positive endpoint and the `S|T` spin is the unique projective
negative endpoint.  The cross cut has `c=1`.

Write `a_X` for total internal edge weight and `m_X,M_X` for the minimum and
maximum internal cut weights.  Exact child enumeration gives

```math
\begin{array}{c|ccc}
X&a_X&m_X&M_X\\ \hline
S&27/56&0&9/28\\
T&1/21&-55/112&19/112.
\end{array}
\tag{B}
```

The minimizing cut in `T` is `B|\{7,8\}`.  Both children are strictly
positive-dominant because

```math
a_S-m_S-M_S=\frac9{56}>0,
\qquad
a_T-m_T-M_T=\frac{31}{84}>0.
\tag{C}
```

Also `a_S+a_T=89/168>0`.  Consequently the hard same-positive inequality
(10.440) fails strictly:

```math
a_S+a_T-(m_S+m_T)
=\frac{49}{48}>c=1.
\tag{D}
```

For completeness, in the doubled quadratic-energy normalization,

```math
h_S=2a_S=\frac{27}{28},
\qquad h_T=2a_T=\frac2{21},
\qquad H=h_S+h_T=\frac{89}{84}.
```

The parent data are

```math
P=2c+H=\frac{257}{84},
\qquad N=2c-H=\frac{79}{84},
\qquad I=|P-N|=\frac{89}{42}.
```

From (B),

```math
Q(S)=2a_S-4m_S=\frac{27}{28},
\qquad
Q(T)=2a_T-4m_T=\frac{173}{84}.
```

Hence the stronger candidate (10.437) has exact deficit

```math
\bigl(Q(S)+Q(T)\bigr)-\bigl(4c-|H|\bigr)
=\frac1{12}>0.
\tag{E}
```

Using the exact endpoint capacity identity (10.433),

```math
b_S=2c-Q(S)+|h_S|=2,
\qquad
b_T=2c-Q(T)+|h_T|=\frac1{28}.
```

Therefore

```math
\mathcal B-I
=\frac{57}{28}-\frac{89}{42}
=-\frac1{12}<0.
\tag{F}
```

So the weighted model violates both (10.435) and (10.437), each by `1/12`.
All values above are rational identities; the exhaustive audit is implemented
in `tmp/strict_weighted_cut_r9.py`.

## 2. Honest `\{\pm1\}` blowups

The strict weighted obstruction implies the existence of finite, unweighted
complete signed graphs with the same failure.

Fix `\kappa=5`.  Since

```math
\max_{i<j}|\kappa w_{ij}|=\frac{565}{672}<1,
```

the following probabilities are valid.  Replace each macro vertex `i` by a
cluster `V_i` of `L` vertices.  Independently assign every edge between
`V_i` and `V_j` a sign with mean

```math
\mathbb E A_{uv}=\kappa w_{ij},
\qquad
\Pr(A_{uv}=1)=\frac{1+\kappa w_{ij}}2.
\tag{G}
```

Assign edges inside each `V_i` independent unbiased signs.  This is an honest
signing of `K_{9L}`.

Let `E=V_4\cup\cdots\cup V_8` be the intended endpoint cut.  For a micro cut
`R`, put

```math
x_i=\frac{|R\cap V_i|}{L}.
```

If `Z_i` are independent Bernoulli variables with parameters `x_i`, then
the expected cut is exactly

```math
\mathbb E C_L(R)=\kappa L^2G(x),
\qquad
G(x)=\mathbb E\,g(\{i:Z_i=1\}).
\tag{H}
```

Thus `G` is the multilinear extension of the exact macro cut function.

### 2.1 Endpoint stability on every Hamming shell

Let `\delta=1/112`.  For a cut `R`, choose its complement if necessary and
put

```math
d=\min\{|R|,9L-|R|\},
\qquad s=\sum_i x_i=\frac dL\le\frac92.
```

The expected number of disagreeing pairs among the nine Bernoulli variables
is

```math
\begin{aligned}
D(x)
&=\sum_{i<j}(x_i+x_j-2x_ix_j)\\
&=8s-s^2+\sum_i x_i^2\\
&\ge8s-\frac89s^2
\ge4s.
\end{aligned}
\tag{I}
```

A nonconstant nine-bit state has at most `\lfloor9^2/4\rfloor=20`
disagreeing pairs.  Consequently

```math
\Pr(Z\text{ is nonconstant})\ge\frac{D(x)}{20}
\ge\frac{s}{5}=\frac{d}{5L}.
\tag{J}
```

By (A), every nonconstant macro cut has value at least `\delta`.  Equations
(H)--(J) imply

```math
\mathbb E C_L(R)\ge\frac{\kappa\delta}{5}Ld.
\tag{K}
```

For the upper endpoint inequality, XOR `R` with `E`.  If

```math
d_E=\min\{|R\mathbin\triangle E|,
9L-|R\mathbin\triangle E|\},
```

then (A), applied to `1-g` and to the endpoint pair in place of the trivial
pair, gives identically

```math
\mathbb E\bigl[C_L(E)-C_L(R)\bigr]
\ge\frac{\kappa\delta}{5}Ld_E.
\tag{L}
```

For (K), only the `d(9L-d)\le9Ld` crossing edges occur.  For (L), the
coefficient of an edge is the difference of its two cut indicators; it is
nonzero exactly on the edge boundary of `R\mathbin\triangle E`.  Hence again
at most `9Ld_E` independent signs occur, each with coefficient `\pm1` and
range length two.

Hoeffding's inequality, asking for at most half the expected gap, therefore
gives for either boundary

```math
\Pr(\text{failure at a cut of shell }d)
\le
\exp\!\left(-\frac{\kappa^2\delta^2}{1800}Ld\right).
\tag{M}
```

There are at most two copies of each shell around either projective boundary.
A union bound for the lower and upper endpoint conditions is thus

```math
4\sum_{d=1}^{\lfloor9L/2\rfloor}
{9L\choose d}
\exp\!\left(-\frac{\kappa^2\delta^2}{1800}Ld\right).
\tag{N}
```

This tends to zero: use
`{9L\choose d}\le(9eL/d)^d`; the negative exponent is linear in `Ld`,
whereas the shell entropy is at most linear in `d\log L`.  In the present
constants the exponential coefficient is

```math
\frac{\kappa^2\delta^2}{1800}=\frac1{903168}.
```

Thus, with probability tending to one, every nontrivial, non-endpoint cut
satisfies

```math
0<C_L(R)<C_L(E).
\tag{O}
```

It follows that `\mathbf1` and the `S|T` cluster spin are the exact unique
projective positive and negative endpoints of the finite signing.

### 2.2 Uniform child extrema

For any fixed `\varepsilon>0`, Hoeffding applied to a child micro cut gives

```math
\Pr\bigl(|C_L^X(R)-\mathbb EC_L^X(R)|
>\varepsilon L^2\bigr)
\le2\exp(-c\varepsilon^2L^2)
```

for an absolute constant `c>0`.  There are at most `2^{5L}` such cuts in
either child.  Hence a union bound gives, with probability tending to one,
uniform `o(L^2)` convergence of all child cut values to their multilinear
weighted expectations.  A multilinear function on a cube attains its
minimum and maximum at vertices.  Therefore

```math
\frac{a_X^{(L)}}{\kappa L^2}\to a_X,
\qquad
\frac{m_X^{(L)}}{\kappa L^2}\to m_X,
\qquad
\frac{M_X^{(L)}}{\kappa L^2}\to M_X.
\tag{P}
```

The internal totals in (P) follow directly from ordinary Hoeffding
concentration; the unbiased intra-cluster edges contribute only `o(L^2)`.
The exact identities

```math
P_X=2a_X-4m_X,
\qquad N_X=4M_X-2a_X
```

then give uniform convergence of `P_X,N_X,Q_X`.  The endpoint cross total
satisfies `c_L/(\kappa L^2)\to1` as well.

The endpoint-stability event (O) and all convergence events in (P) have
probability tending to one, so their intersection does too.  On that
intersection, the exact endpoint identity (10.433) applies.  Dividing by
`\kappa L^2` and using (F) gives

```math
\frac{\mathcal B_L-I_L}{\kappa L^2}
\longrightarrow-\frac1{12}.
\tag{Q}
```

Thus, for every sufficiently large `L`, at least one finite complete
`\{\pm1\}` signing violates (10.435).  Equation (E) and the same convergence
show that it also violates (10.437).

The strict dominance gaps in (C) persist under (P), and

```math
\frac{a_S^{(L)}+a_T^{(L)}-(m_S^{(L)}+m_T^{(L)})}{c_L}
\longrightarrow\frac{49}{48}>1.
```

Hence these honest signings also violate (10.440) in its stated hard
same-positive regime.  Since any choice satisfying the proposed sufficient
condition (10.442) would imply (10.440), the universal existence form of
(10.442) fails on the same eventual blowups.

## 3. Scope

This disproves the local stopping targets (10.435), (10.437), (10.440), and
the universal proposed selection statement (10.442).  It invalidates the
specific depth-two argument that used (10.435) to obtain `K\le4`.

It **does not** disprove a universal optimized congestion bound `K\le4`.
A different allocation, a different endpoint choice, or a deeper stopping
construction might still achieve that bound.  Nor does this construction by
itself give `K_{\min}>4` or solve the original quadratic signing limit
problem.
