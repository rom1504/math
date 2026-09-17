# Joint Bernoulli increments and all-window full-sign regularization

Date: 2026-09-17. Primary mechanism: information-controlled adaptive
selection, applied to a permutation-symmetric collection of random
incident edges. This is a finite theorem for actual full signings, not
a Gaussian surrogate or a separate payment over new-spin patterns.

The canonical campaign theorem is the director's
[exchangeable-sign regularization](paper_director_exchangeable_sign_regularization_2026_09_17.md).
This file is its independent reconstruction, the earlier column
argument, and a source/replay record, not a competing canonical version.
The director proposed the all-incident-edge upgrade in Section 4 below.
The Bernoulli track reconstructed its fractional information bound and
the complete physical compiler independently. The earlier independent
column argument is retained in Section 7 to explain the route and its
strictly narrower scope. No external novelty claim is made.

## 1. Notation and finite theorem

For a full symmetric hollow signing A on n vertices, write

```
H_A(x)=sum_{i<j} A_ij x_i x_j,   Q(A)=max_x |H_A(x)|.
```

For a nonempty Boolean code C on N coordinates define its Bernoulli
and Gaussian widths by

```
b(C)=E_epsilon max_{x in C} epsilon dot x,
w(C)=E_g max_{x in C} g dot x.
```

Start with an arbitrary fixed A on n vertices. Add q new vertices and
choose EVERY edge incident to a new vertex independently and uniformly
from {-1,1}. Let W be the resulting order-N signing, N=n+q, and put

```
ell=nq+q(q-1)/2,
U=sqrt(2 ell [(N+1)log 2+log 4]),
mu_N=E|epsilon_1+...+epsilon_N|,
L=mu_N+2 sqrt(N(N+1)log 2/(q+1)).                     (1)
```

With probability at least 1/2, simultaneously

```
Q(W)<=Q(A)+U,                                         (2)
b(E_W(T))<=min{N,T+4L}    for EVERY T>=0,             (3)
```

where E_W(T)={z: Q(W)-|H_W(z)|<=T} is the FULL absolute near-level
code on all N vertices. In particular one realization satisfies both
claims. There is no restriction q^2<=n or q^3<=n. The theorem remains
true with n=0 if the empty initial form is interpreted as zero.

For B=min{N,T+4L}, this same realization obeys

```
log |E_W(T)| <= B log(eN/B),                           (4)
w(E_W(T)) <= B sqrt(2log(eN/B)) +sqrt(2/pi) B/e.       (5)
```

All logarithms are natural. The zero-width singleton convention in
Sections 5--6 resolves B=0, although here L>0 whenever N>=1.

## 2. Two elementary information lemmas

### Adaptive subGaussian selection

Suppose R is random, U is a finite-valued selected index, and for every
deterministic u,

```
E f(R,u)=mu,
log E exp(lambda[f(R,u)-mu]) <=lambda^2 v/2.
```

Then, without any independence assumption between U and R,

```
E f(R,U)<=mu+sqrt(2v I(U;R)).                         (6)
```

Indeed, compare the joint law of (U,R) with its product marginal law.
The variational inequality for relative entropy gives
lambda(Ef-mu)<=I(U;R)+lambda^2 v/2; minimize over lambda>0.
This is the standard information-theoretic selection-bias inequality,
not a new probabilistic principle. The decisive primary proof is
Russo--Zou, AISTATS 2016, Proposition 1, including the short proof in
their [supplement](https://proceedings.mlr.press/v51/russo16-supp.pdf).
The [main paper](https://proceedings.mlr.press/v51/russo16.pdf) states
the result in adaptive-data-analysis notation.

### Fractional packing of information over independent coordinates

Let (X_i) be independent, and let S_1,...,S_r be coordinate sets, with
each coordinate appearing in at most a of the sets. For arbitrary U,
including U with additional independent tie-breaking randomness,

```
sum_j I(U;X_{S_j}) <= a I(U;X_all) <= a H(U).         (7)
```

To see this directly, order the coordinates. For i in S, independence
and conditioning reduce entropy to give

```
I(U;X_i | X_{S intersect {indices<i}})
 <= I(U;X_i | X_{indices<i}).
```

Sum first over i in S using the chain rule, then over the sets, and
use their maximum multiplicity. Thus no unproved Shearer variant or
independence conditional on U is being assumed.

For independent fair signs R=(R_1,...,R_m), every fixed z in {-1,1}^m
has response f(R,z)=|R dot z| with common mean mu_m and centered
subGaussian proxy m. Changing a single sign changes f by at most 2;
the bounded-difference MGF bound is lambda^2 sum_i 2^2/8=m lambda^2/2.

## 3. Conditional add-one increment controls all near-level windows

For ANY fixed full signing W on N vertices, append one new vertex with
an independent fair incident row h. If W+h denotes that signing, then

```
Q(W+h)=max_z [|H_W(z)|+|h dot z|].                    (8)
```

This follows from max_{t=+-1}|a+t b|=|a|+|b|. Define

```
Delta(W)=E_h Q(W+h)-Q(W) >=0.                        (9)
```

Restricting the maximum in (8) to E_W(T) gives

```
Delta(W)>= E_h max_{z in E_W(T)} |h dot z|-T
            >= b(E_W(T))-T.                         (10)
```

In fact E_W(T) is antipodal, so the last inequality is equality for
the widths, but equality is not needed. Crucially, ONE number Delta(W)
controls every window T at once. There is no union bound over windows
or over the 2^q new-spin choices.

## 4. Expected increment for the random-incident-edge model

Consider the model with q+1 new vertices, of total order N+1. From all
its absolute groundstates choose Z uniformly, using an independent
tie-breaking random variable. The law of this choice is equivariant
under permutations of the new vertices; a lexicographic selector on
all vertices would not automatically have that property.

For each new vertex j let R_j be its complete incident random row.
Every random edge belongs to at most two such rows: old--new edges
belong to one and new--new edges to two. Formula (7) therefore gives

```
sum_{j=1}^{q+1} I(Z;R_j) <=2H(Z)<=2(N+1)log 2.
```

Permutation symmetry makes all terms equal. By (6), with proxy N and
the deterministic index z consisting of Z's remaining coordinates,

```
E |R_j dot Z_{-j}|
 <=mu_N+sqrt(2N I(Z;R_j))<=L.                        (11)
```

Using the full vector Z instead of its restriction only increases the
information and is harmless. All N entries of the marginal R_j are
independent fair signs, including its edges to the other new vertices.

Deleting vertex j from the chosen maximizing spin configuration gives

```
Q(W_{q+1})-Q(W_{q+1} minus j)
 <= |R_j dot Z_{-j}|.                                (12)
```

This is just |a+b|-|a|<=|b|; it does not fix an energy polarity.
The deleted signing has exactly the q-new-vertex distribution. Thus

```
E_W Delta(W)=E Q(W_{q+1})-E Q(W_q)<=L.               (13)
```

Markov's inequality gives P(Delta(W)>4L)<=1/4. On its complement,
(10) proves (3) for all T. Separately, for each fixed full spin z the
random part of H_W(z) is a sum of ell independent fair signs. Hence

```
P(max_z |H_random(z)|>u) <=2^(N+1)exp[-u^2/(2ell)]
```

and the value U in (1) gives failure at most 1/4. Intersect the two
events to obtain probability at least 1/2 and (2). When ell=0 the
random-part bound is interpreted deterministically as U=0.

## 5. Boolean entropy directly from Bernoulli width

For EVERY nonempty C subset {-1,1}^N,

```
log|C| <= b(C) log(eN/b(C)).                          (14)
```

Here the right side is zero when b(C)=0. Let d be the VC dimension of
C viewed as a family of subsets of its N-coordinate domain. If a set
S of d coordinates is shattered, choose, for each epsilon_S, one fixed
word x(epsilon_S) realizing that pattern. This choice depends only on
epsilon_S. Independence of all remaining random signs yields

```
b(C)>=E epsilon dot x(epsilon_S)=d.                  (15)
```

The Sauer bound is |C|<=sum_{i=0}^d binom(N,i). For completeness, split
the family at its last coordinate. The union of the two sections has
dimension at most d, and their intersection has dimension at most d-1;
cardinality is the sum of union and intersection sizes. Induction with
the Pascal recursion proves the displayed bound, including the empty
section and d=0 cases. For 1<=d<=N, put t=d/N<=1 and use

```
sum_{i<=d}binom(N,i)<=t^(-d)(1+t)^N<=(eN/d)^d.
```

The function u log(eN/u) increases on [0,N], and d<=b(C)<=N, giving
(14). If d=0 there is only one word and both entropy and width vanish.
Classical attribution: N. Sauer, *On the density of families of sets*,
JCTA 13 (1972), 145--147,
[DOI](https://doi.org/10.1016/0097-3165(72)90019-2).
The entire combinatorial argument needed here is supplied above.

The logarithmic loss cannot generally be removed: the code consisting
of the all-ones word and its N single-coordinate flips has width
2(1-2^(-N)) but entropy log(N+1).

## 6. Optional Gaussian width and same-order consequence

For any B0>0, write each independent Gaussian as its symmetric clipping
to [-B0,B0] plus its tail. A centered variable in this interval is
dominated in convex order by B0 times a fair sign. Successively apply
this one-dimensional fact to the coordinatewise convex function
max_{x in C} z dot x. Consequently

```
w(C)<=B0 b(C)+N E(|G|-B0)_+
     <=B0 b(C)+sqrt(2/pi)N exp(-B0^2/2).
```

Choose B0=sqrt(2log(eN/b(C))) to obtain (5), first with b(C) itself
and then with its larger bound B; the resulting expression increases
on (0,N]. If b(C)=0, let B0 tend to infinity. This clipping is the
elementary bounded-coordinate/tail decomposition behind the comparison;
the full Bednorz--Latala theorem is not required for this application.

Given ANY signing A_N and any 1<=q<N, keep its principal block on
n=N-q vertices and refill every edge incident to the other q vertices
as above. Principal restriction cannot increase Q: average over the
deleted signs, keeping the surviving signs fixed. The resulting
same-order signing W differs only on at most Nq edges and satisfies

```
Q(W)<=Q(A_N)+O(N sqrt(q)),
b(E_W(N/sqrt(q)))=O(N/sqrt(q)),
w(E_W(N/sqrt(q)))=O((N/sqrt(q))sqrt(log(eq))),
log|E_W(N/sqrt(q))|=O((N/sqrt(q))log(eq)).             (16)
```

These are uniform bounds with universal constants. If q tends to
infinity while q=o(N), the cap error is o(N^(3/2)) and the displayed
full near-level code is subexponential. In particular q=floor(N/log N)
gives cap cost O(N^(3/2)/sqrt(log N)), window and Bernoulli width
O(sqrt(N log N)), Gaussian width O(sqrt(N)log N), and entropy
O(sqrt(N)(log N)^(3/2)). The change is incident to o(N) vertices.

This applies to exact minimizers at EVERY order; it gives an equivalent
asymptotic variational class with controlled shrinking-window geometry.
It does NOT prove convergence of M_N/N^(3/2), improve either sharp
constant bound, or give low response under an isotropic sign law.
The window is o(N^(3/2)); universal fixed-positive-deficit entropy and
the Hadamard all-law response obstruction are unaffected.

## 7. Earlier independent-column version, with its original scope

Let H0 be any deterministic function on {-1,1}^n and take independent
fair columns h_1,...,h_q. Define

```
F_q=max_x [H0(x)+sum_{j<=q}|h_j dot x|/sqrt(n)].
```

At q+1 columns select the lexicographically first maximizing x, using
a fixed lexicographic order on x only. The objective is symmetric in
the column labels, so this selector is invariant under column order.
Independence and (7) with multiplicity one give

```
I(X;h_j)<=n log2/(q+1).
```

The normalized absolute response has mean mu_n/sqrt(n)<=1 and proxy
one. Removing a column therefore yields

```
E(F_{q+1}-F_q)<=mu_n/sqrt(n)+sqrt(2nlog2/(q+1)).      (17)
```

Conditioning on the first q columns, their conditional increment
controls every objective near-level Bernoulli width by
sqrt(n)[margin+increment]. Taking H0=|H_A|/sqrt(n) produces the exact
old-spin joint maximum after optimizing all new signs. If the new
child D is fixed separately, however, the actual full-parent window
inflates by 2Q(D), and the new-coordinate width costs q. With
Q(D)=O(q^(3/2)) this version only delivers the scale n/sqrt(q) when
q^2=O(n). Section 4 genuinely removes this compiler loss by treating
the whole random child together with the bridge.

## 8. Verification and attribution boundary

The director and localization track independently reconstructed the
column increment, the adaptive information bound, and the direct
VC/Sauer entropy conversion. The director supplied the decisive
all-incident-edge extension, including overlap multiplicity two.
The Bernoulli track's finite replay
`computations/paper_bernoulli_2026_09_17_joint_increment.py` passed:
4,840 exact random-graph atoms across six information examples,
442 conditional near-level windows, and all 273 nonempty Boolean
codes of orders 1--3. Uniform groundstate tie probabilities, expected
caps, conditional increments, widths, and VC checks used exact rational
arithmetic. Mutual informations and the logarithmic bounds were
evaluated numerically. The add-one mean identity was exact in every
case; all row informations agreed to numerical precision. Output is
preserved at `tmp/paper_portfolio_2026_09_17/bernoulli/joint_increment_audit.json`.
The script also passed Python compilation. No numerical calculation
substitutes for the finite proof.

The results are a composition of classical selection-bias information,
independent-coordinate entropy, and Boolean shattering with the exact
full-sign add-one identity. They are not being attributed to a new
Bernoulli-process theorem, and no claim of literature novelty is made.
