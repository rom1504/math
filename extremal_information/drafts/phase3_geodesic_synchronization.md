# Geodesic cycle contraction and sharp deterministic synchronization

**Status.**  The statements in Sections 2--9 are proved.  The finite claims
and the sharp example at `(D,k)=(6,4)` are checked by
[`verify_phase3_geodesic_synchronization.py`](../experiments/verify_phase3_geodesic_synchronization.py).

This note extracts the global compatibility that was not used in the
fibre-by-fibre count of `phase3_geodesic_fibre_hard_core.md`.  A shortest
binary word is equivalent to a contraction inequality on *every* additive
cycle in the quotient.  When the quotient support is complete, that
inequality forces every transversal to lie at bounded **uniform** Hamming
distance from a linear section.  This in turn gives an all-raw-future
response quotient: replacing the entire support by the linear-section
support changes every appended-fragment word metric by at most eleven.  The
average bound is a vector-valued BLR theorem, and its constant `3` is
asymptotically sharp even for supports whose chosen word realizes the full
Cayley diameter.

Thus dense geodesic supports supply a strict all-context response quotient
with quadratic-size description, despite potentially exponential raw
support complexity.  For incomplete projection, the same mechanism has an
error controlled by the quotient Cayley diameter.  Exact summary-only
composition is not claimed.

## 1. Setup

Let `G=F_2^w`, let `S subseteq G\{0}` span `G`, and let

```math
B=\{b_1,\ldots,b_D\}\subseteq S
```

be linearly independent.  Put

```math
W=\operatorname{span}B,\qquad
t=b_1+\cdots+b_D,qquad
Q=G/W,qquad
\pi:G\longrightarrow Q.                         \tag{GS.1}
```

For `u in W`, write `|u|_B` for the Hamming weight of its unique
`B`-coordinate vector.  For a subset `R`, `sum R` always means its binary
sum.  Repetitions need not be considered: two copies of the same generator
cancel.

## 2. The exact cycle-contraction criterion

### Theorem GS.1 (geodesic iff quotient cycles contract)

The representation `t=sum B` is shortest among representations of `t` by
subsets of `S` if and only if

```math
\left|\sum_{s\in R}s\right|_B\le |R|             \tag{GS.2}
```

for every `R subseteq S\B` satisfying

```math
\sum_{s\in R}\pi(s)=0.                           \tag{GS.3}
```

Thus, if `B` is a shortest representation of a diametral point, (GS.2)
holds.  Conversely, (GS.2) certifies that `B` is a length-`D` geodesic; by
itself it does **not** assert that no other point has distance greater than
`D`.

#### Proof

Suppose first that `B` is shortest.  For a projected-zero `R`, put

```math
I=\operatorname{supp}_B\left(\sum_{s\in R}s\right).
```

Then

```math
t=\sum_{s\in R}s+\sum_{i\notin I}b_i.
```

After cancelling any repeated generators this is a representation of `t`
of length at most `|R|+D-|I|`.  Shortness gives
`D<=|R|+D-|I|`, hence `|I|<=|R|`.

Conversely, let `A subseteq S` represent `t`, write
`I=A cap B` and `R=A\B`.  Projection gives (GS.3), while

```math
\sum_{s\in R}s=t+\sum_{b_i\in I}b_i
               =\sum_{b_i\notin I}b_i.
```

Consequently (GS.2) says `D-|I|<=|R|`, and therefore
`|A|=|I|+|R|>=D`.  `square`

This is an exact reformulation of shortest-path minimality, not yet a new
compression theorem.  Its value is that it turns minimality into a family
of finite compatibility inequalities on quotient cycles.

## 3. A pointwise linear chart

Because `S` spans `G`, the set `pi(S\B)` spans `Q`.  Choose a basis
`q_1,...,q_k` of `Q`, where `k=w-D`, and representatives
`c_j in S\B` with `pi(c_j)=q_j`.  There is a unique linear section

```math
L_0:Q\longrightarrow G,\qquad L_0(q_j)=c_j.     \tag{GS.4}
```

Let `|q|_q` denote Hamming weight in this quotient basis.

### Corollary GS.2 (triangularization along a geodesic)

If (GS.2) holds, then every `s in S` obeys

```math
|s+L_0(\pi s)|_B\le |\pi s|_q+1.                \tag{GS.5}
```

Here the left side is well defined because the argument lies in `W`.

#### Proof

The assertion for `s=b_i` is `1<=1`.  For `s notin B`, take the symmetric
difference of `{s}` with the selected representatives `c_j` occurring in
the quotient-coordinate expansion of `pi(s)`.  Its projection is zero, its
cardinality is at most `1+|pi(s)|_q`, and its sum is
`s+L_0(pi(s))`.  Apply (GS.2). `square`

This already reduces a general generator to a quotient label plus at most
`k+1` bits of geodesic error.  Complete quotient support improves `k+1` to
a uniform constant.

## 4. Vector-valued BLR synchronization

We use the following elementary aggregate form of the BLR stability lemma.

### Lemma GS.3 (joint BLR bound with a uniform upgrade)

Let `Q=F_2^k`, let `f:Q to F_2^D` satisfy `f(0)=0`, and set

```math
\Delta_f(x,y)=f(x)+f(y)+f(x+y).
```

Then there is a linear map `A:Q to F_2^D` such that

```math
2^{-k}\sum_{x\in Q}|f(x)+A(x)|
\le
\mathbb E_{x,y}|\Delta_f(x,y)|.                 \tag{GS.6}
```

If in addition

```math
|\Delta_f(x,y)|\le\delta\qquad\text{for every }x,y, \tag{GS.7}
```

then the same `A` may be chosen so that

```math
\max_x |f(x)+A(x)|\le3\delta.                  \tag{GS.8}
```

#### Proof

Write `f=(f_1,...,f_D)` and

```math
\delta_j=\Pr_{x,y}[f_j(x)+f_j(y)+f_j(x+y)=1].
```

For `g_j=(-1)^{f_j}`, Fourier expansion on `Q` gives

```math
1-2\delta_j
=\mathbb E_{x,y}g_j(x)g_j(y)g_j(x+y)
=\sum_\alpha \widehat g_j(\alpha)^3.
```

If `M_j=max_alpha \widehat g_j(alpha)`, Parseval implies

```math
\sum_\alpha\widehat g_j(\alpha)^3
\le M_j\sum_\alpha\widehat g_j(\alpha)^2=M_j.
```

Choose a character attaining `M_j`; its corresponding linear Boolean
function `a_j` differs from `f_j` on a fraction
`p_j=(1-M_j)/2<=delta_j` of `Q`.  Taking
`A=(a_1,...,a_D)` and summing over `j` proves (GS.6).

For the uniform conclusion, put `e_j=f_j+a_j` and define the local rejection
probability

```math
r_j(x)=\Pr_y[e_j(x)+e_j(y)+e_j(x+y)=1].
```

If `e_j(x)=1`, then rejection occurs exactly when
`e_j(y)=e_j(x+y)`.  Translation invariance and the union bound give

```math
r_j(x)\ge1-\Pr_y[e_j(y)\ne e_j(x+y)]\ge1-2p_j.
```

Therefore, for every fixed `x`,

```math
\begin{aligned}
|f(x)+A(x)|
&\le\sum_j r_j(x)+2\sum_jp_j\\
&\le\mathbb E_y|\Delta_f(x,y)|
  +2\,\mathbb E_{u,v}|\Delta_f(u,v)|\\
&\le3\delta.
\end{aligned}
```

Linearity of `A` was used in the equality of the defects.  This proves
(GS.8). `square`

### Theorem GS.4 (dense geodesic synchronization)

Assume (GS.2) and that every nonzero `q in Q` occurs in `pi(S\B)`.  Choose
an arbitrary representative

```math
s_q\in S cap \pi^{-1}(q)\quad(q\ne0),qquad s_0=0.
```

Then some linear section `L:Q to G` satisfies

```math
\sum_{q\in Q}|s_q+L(q)|_B
\le3\,2^k                                      \tag{GS.9}
```

More precisely the right side may be replaced by

```math
3\,(2^k-1)(2^k-2)/2^k.                         \tag{GS.10}
```

The same section has the uniform bound

```math
\max_{q\in Q}|s_q+L(q)|_B\le9.                 \tag{GS.11}
```

Moreover, for the complete fibres
`F_q=S cap pi^{-1}(q)` one has

```math
\max_{q\ne0}\max_{s\in F_q}|s+L(q)|_B\le11.   \tag{GS.12}
```

#### Proof

Start from any linear section `L_0` and identify `W` with `F_2^D` in the
`B` coordinates.  Put `f(q)=s_q+L_0(q)`.  If `x,y` are nonzero and
distinct, then `x,y,x+y` are three distinct nonzero quotient vectors, so
(GS.2), applied to their three representatives, gives

```math
|\Delta_f(x,y)|\le3.                           \tag{GS.13}
```

In the remaining pairs `Delta_f(x,y)=0`.  Hence its expectation is bounded
by (GS.10) divided by `2^k`, while its pointwise norm is at most three.
Lemma GS.3 gives a linear `A:Q to W`; setting `L=L_0+A` proves
(GS.9)--(GS.11).

Two generators in one quotient fibre form a projected-zero two-set, so
(GS.2) also gives

```math
|s+s_q|_B\le2\qquad(s\in F_q).
```

Thus the maximum in (GS.12) is at most
`2+|s_q+L(q)|_B`.  Use (GS.11). `square`

The theorem is deterministic: no disorder average, exchangeability, or
assumed overlap linkage is present.  It also keeps the `D` output channels
joint until after the cycle inequality; bounding the channels first would
lose the crucial uniform total defect bound.

## 5. All-raw-future response compression

For a spanning support `U`, let `ell_U(x)` be its binary word metric and
`rho(U)=max_x ell_U(x)`.  Given the section in Theorem GS.4, define its
linear support

```math
S_L=B\cup\{L(q):q\in Q\setminus\{0\}\}.         \tag{GS.14}
```

### Theorem GS.5 (dense geodesic raw-context quotient)

Under the hypotheses of Theorem GS.4, for every arbitrary appended support
`T subseteq G\setminus\{0\}` and every `x in G`,

```math
\ell_{S_L\cup T}(x)\le\ell_{S\cup T}(x)+11,
\qquad
\ell_{S\cup T}(x)\le\ell_{S_L\cup T}(x)+9.      \tag{GS.15}
```

Consequently

```math
\sup_T\sup_x
|\ell_{S\cup T}(x)-\ell_{S_L\cup T}(x)|\le11,  \tag{GS.16}
```

and in particular every appended-fragment covering-radius response differs
by at most eleven.

#### Proof

Take a representation using `S union T` and partition its generators into
those retained from `T`, those in `B`, and a subset `R subseteq S\B`.
Elements belonging to both source and future may be assigned arbitrarily;
any cancellations in the replacements below only shorten the word.  Put

```math
q=\sum_{s\in R}\pi(s).
```

If `q=0`, cycle contraction replaces all of `R` by at most `|R|` members of
`B`, with no overhead.  If `q` is nonzero, compare `R` with the selected
representative `s_q`.  The symmetric difference
`R` symmetric-difference `{s_q}` is projected-zero, so

```math
\left|\sum_{s\in R}s+s_q\right|_B
\le |R\mathbin\triangle\{s_q\}|\le |R|+1.
```

Thus `R` can be replaced by `s_q` and members of `B` at overhead at most
two.  By (GS.11), `s_q` can then be replaced by `L(q)` and at most nine
members of `B`.  This proves the first inequality in (GS.15).

Conversely, group all generators from `S_L\B` in a representation using
`S_L union T`.  Their sum is exactly `L(q)`, where `q` is the sum of their
labels.  If `q=0`, delete them; if `q` is nonzero, replace them by the single
generator `L(q)`.  Replacing that generator by `s_q` and at most nine
members of `B` gives a representation using `S union T` at overhead at
most nine.  This proves the reverse inequality, hence (GS.16). `square`

The summary `(B,L)` uses at most `wD+wk=w^2` raw matrix bits.  Choosing, for
example, the lexicographically first valid `(B,L)` makes it a deterministic
summary map.  The theorem is stronger than a one-shot estimate of
`rho(S)`: it answers every raw future support and every rooted word-length
query with constant error.  It is not asserted here that `(B,L)` carries an
exact homomorphic product for combining two already-compressed summaries;
the proved statement is the all-raw-context response property (GS.16).

### Corollary GS.6 (full-class response entropy and bounded composition)

Fix the chart `B subseteq G` and quotient `Q`.  In the all-context metric
on dense geodesic supports defined by the left side of (GS.16), the class
has an `11`-cover with at most

```math
2^{Dk}                                           \tag{GS.17}
```

centres.  If the chart may also vary, `O(w^2)` bits still suffice.  More
generally, if `S_1,...,S_m` each satisfy the dense hypothesis in their own
charts and `S_{L_i}` are their summaries, then for every raw future `T`,

```math
\sup_x\left|
\ell_{T\cup S_1\cup\cdots\cup S_m}(x)
-\ell_{T\cup S_{L_1}\cup\cdots\cup S_{L_m}}(x)
\right|\le11m.                                  \tag{GS.18}
```

#### Proof

For fixed `B`, two linear sections differ by an element of
`Hom(Q,W)`, so there are exactly `2^(Dk)` possible sections.  Theorem GS.5
maps every source to one of their supports at radius eleven.  A varying
ordered basis and section can be stored by at most `wD+wk=w^2` bits.
For (GS.18), replace the sources one at a time.  At each replacement, the
union of `T` with every other raw or already summarized source is an
arbitrary future support allowed in Theorem GS.5.  The pointwise triangle
inequality gives `11m`. `square`

This is bounded-composition control, not an exact summary-only product.  It
is nevertheless a strict quotient of the full source landscape whenever
the fibres contain exponentially many raw membership choices.

## 6. Quotient diameter controls the price of fibre information

The dense case has quotient diameter one.  The following stripping law
shows what replaces the constant two when quotient labels are incomplete.
It is naturally a theorem about group extensions, not specifically about
binary codes.

### Theorem GS.7 (group-extension fibre stripping)

Let

```math
0\longrightarrow K\longrightarrow G
\overset{\pi}{\longrightarrow}Q\longrightarrow0
```

be an extension of finite abelian groups.  Let `B=-B` generate `K`, let
`A=-A` be a source generator set in `G\setminus K`, and assume the
cycle-contraction
property

```math
\ell_B(s_1+\cdots+s_r)\le r                    \tag{GS.19}
```

for every source word `s_1,...,s_r in A` whose projected sum is zero.
Assume `U=pi(A)` generates `Q`, let its Cayley diameter be `h`, and suppose we have
chosen lifts `s_u in A` satisfying `s_{-u}=-s_u` for every `u in U`.  Write
`A^{tr}={s_u:u in U}`.  Then, for every future symmetric generator set
`T` and every `x in G`,

```math
0\le
\ell_{A^{tr}\cup B\cup T}(x)-\ell_{A\cup B\cup T}(x)
\le2h.                                          \tag{GS.20}
```

#### Proof

The first inequality is monotonicity because `A^{tr} subseteq A`.  In a
word over `A union B union T`, group the `r` letters from `A`.  If their
quotient sum is `q`, choose a quotient word of length `m<=h` for `q` and
replace its labels by the selected lifts.  Concatenating the original word
with the inverses of those lifts gives a projected-zero source word of
length `r+m`.  By (GS.19), its sum is represented by at most `r+m` kernel
letters.  The original `r` source letters can therefore be replaced by
`m` selected lifts and at most `r+m` kernel letters, for total length at
most `r+2h`.  The future and pre-existing kernel letters are left
unchanged; cancellations can only help. `square`

The sign-compatible lift is an explicit hypothesis in a general extension;
it is automatic in an elementary abelian `2`-group.  In the binary setting
of Theorem GS.1, take `K=W`, use the independent basis `B` as the kernel
generators, take `A=S\B`, and put `U=pi(A)`.  Applying (GS.2) to a singleton
in `S cap W` shows that it has `B`-weight at most one, hence `S cap W=B`;
therefore `A` is disjoint from `W` as required.  Repetitions cancel, so the
theorem applies verbatim to subset words.

### Proposition GS.8 (the factor two is sharp)

For every `h>=1`, there is a binary example satisfying (GS.2) for which the
upper gap in (GS.20) is exactly `2h`, already with no future support.

#### Proof

Take `Q=F_2^h` with basis `q_1,...,q_h` and let `W=F_2^{2h}` have basis
`B={e_1,...,e_{2h}}`.  In the fibre over `q_i`, put

```math
c_i=(0,q_i),\qquad
c_i'=c_i+e_{2i-1}+e_{2i},                       \tag{GS.21}
```

and set `S=B union {c_i,c_i':1<=i<=h}`, selecting `c_i` for the transversal.
A projected-zero subset of `S\B` contains either neither or both members of
each fibre.  Its `B`-weight is exactly its cardinality, so (GS.2) holds.

The point `x=sum_i c_i'` has length `h` in `S`.  In the selected-transversal
support, quotient independence forces all `h` generators `c_i`, after which
all `2h` kernel coordinates remain.  Its length is exactly `3h`.  The gap
is `2h`. `square`

There is also a full selector cube behind this pair.  For `J subseteq[h]`,
let `S_J` contain every `c_i` and contain `c_i'` exactly when `i in J`.
For the rooted queries `x_P=sum_{i in P}c_i'`, disjoint coordinates give

```math
\ell_{S_J}(x_P)=|P|+2|P\setminus J|,
\qquad
\sup_P|\ell_{S_J}(x_P)-\ell_{S_K}(x_P)|
=2\max\{|J\setminus K|,|K\setminus J|\}.       \tag{GS.22}
```

Ordinary Hamming packings therefore give `Omega(h)` response bits at fixed
linear distortion.  This recovers only the established selector-cube rate;
its role here is to show that the sharp `2h` gap is an information family,
not an isolated pair.

This gives a precise information-growth law: microscopic choices inside a
fibre cost at most twice the quotient diameter in every future word query,
and they can attain that cost.  Hence they are uniformly negligible at
scale `w` when `h=o(w)`, but can amplify to a macroscopic response when
`h=Theta(w)`.

## 7. Partial-quotient synchronization

The quotient-diameter parameter also extends the linear-section summary
beyond complete quotient projection.

### Theorem GS.9 (small quotient diameter gives a linear response summary)

Return to the binary hypotheses of Theorem GS.1 and assume that `Q` is
nontrivial.  Put

```math
P=\pi(S\setminus B),\qquad
h=\max_{q\in Q}\ell_P(q).                       \tag{GS.23}
```

Choose one source representative for every label in `P`.  For every `q`,
choose a shortest word of at most `h` selected representatives with
quotient sum `q`, and let `g(q)` be its sum in `G`.  Then there is a linear
section `L:Q to G` such that

```math
\max_{q\in Q}|g(q)+L(q)|_B\le9h.               \tag{GS.24}
```

For its complete linear support `S_L` from (GS.14), every raw future `T`
and every root `x` satisfy

```math
\ell_{S_L\cup T}(x)\le\ell_{S\cup T}(x)+10h+1,
\qquad
\ell_{S\cup T}(x)\le\ell_{S_L\cup T}(x)+10h-1. \tag{GS.25}
```

Thus `(B,L)` is an `O(w^2)`-bit all-context summary with response distortion
at most `10h+1`.  In particular, `h=o(w)` gives `o(w)` distortion.

If `Q=0`, cycle contraction forces `S=B`, so the exact summary is immediate;
the separate assumption only avoids writing the vacuous reverse bound
`10h-1=-1` when `h=0`.

#### Proof

For `x,y in Q`, take the symmetric difference of the three chosen source
words for `x,y,x+y`.  It is a projected-zero subset of at most `3h`
generators.  The cycle criterion gives

```math
|g(x)+g(y)+g(x+y)|_B\le3h.
```

Choose any linear section `L_0` and put `f(q)=g(q)+L_0(q) in W`.  The
linear terms cancel from its additive defect, so Lemma GS.3, with
`delta=3h`, gives `A:Q to W`; setting `L=L_0+A` proves (GS.24).

For the first inequality in (GS.25), group a subset `R subseteq S\B` used
by a source word and let its quotient sum be `q`.  If `q=0`, replace `R` by
at most `|R|` basis elements.  Otherwise compare it with the chosen word
`C_q` defining `g(q)`.  Cycle contraction and (GS.24) give

```math
\left|\sum_{s\in R}s+L(q)\right|_B
\le |R|+h+9h.
```

Using one generator `L(q)` therefore costs at most `10h+1` beyond `|R|`.
Conversely, group all linear-section generators in a word over `S_L`.
Their sum is `L(q)`.  For nonzero `q`, replace it by the at most `h` generators
in `C_q` and at most `9h` basis elements; relative to the nonempty original
group the overhead is at most `10h-1`.  For `q=0`, delete the group.  This
proves (GS.25). `square`

Theorem GS.5 is the `h=1` instance with the same constants.  Proposition GS.8 shows that
dependence on `h` cannot be removed from a theorem based only on fibre
stripping and cycle contraction, although the optimal constants in
(GS.24)--(GS.25) remain open.

## 8. An information consequence for transversals

For fixed `B`, quotient coordinates, and full quotient support, an arbitrary
transversal has `D(2^k-1)` raw bits.  Theorem GS.4 gives a smaller covering
description.

### Corollary GS.10 (metric entropy of geodesic transversals)

For `D>=8` and `k>=2`, the number of transversal maps satisfying (GS.2) is
at most

```math
2^{Dk}
\sum_{j\le3\,2^k}{D(2^k-1)\choose j}.           \tag{GS.26}
```

Consequently its base-two logarithm is

```math
O\!\left(Dk+2^k\log D\right),                  \tag{GS.27}
```

instead of the unrestricted `D(2^k-1)` bits.

#### Proof

Encode the linear map in (GS.9), then the set of at most `3*2^k` erroneous
entries in its `D` by `(2^k-1)` truth table.  This proves (GS.26).
For `D>=8`, the error fraction is at most `1/2`; the standard Hamming-ball
entropy estimate gives (GS.27). `square`

This is sub-landscape complexity when `D` grows.  It is a statement about
one representative per quotient fibre, not about the number of full
supports and not yet about their response-metric entropy.

## 9. The constant three is asymptotically sharp

The BLR constant in Theorem GS.4 cannot be improved under the stated
geodesic hypothesis, even after imposing that `t` is diametral.

### Theorem GS.11 (three coupled bent pairs)

For every `D>=6` and every positive multiple `k` of four, there are
`G=W direct-sum Q`, `dim W=D`, `dim Q=k`, a spanning support `S`, and a
diametral geodesic basis `B` of `W` such that `pi(S\B)=Q\{0}` and the unique
representative `s_q` in each nonzero quotient fibre satisfies

```math
\min_{L\text{ linear section}}
\sum_{q\in Q}|s_q+L(q)|_B
=3\left(2^k-2^{k/2}-1\right).                  \tag{GS.28}
```

In particular, the optimal universal leading constant in (GS.9) is exactly
three.

#### Construction and proof

On `F_2^4`, take the alternating matrices

```math
A_1=\begin{pmatrix}
0&0&0&1\\0&0&1&0\\0&1&0&0\\1&0&0&0
\end{pmatrix},\quad
A_2=\begin{pmatrix}
0&0&1&1\\0&0&0&1\\1&0&0&0\\1&1&0&0
\end{pmatrix},\quad
A_3=\begin{pmatrix}
0&0&1&0\\0&0&1&1\\1&1&0&0\\0&1&0&0
\end{pmatrix}.                                  \tag{GS.29}
```

Each is nonsingular, and `A_3=A_1+A_2`.  For `k` divisible by four, take
block-diagonal copies.  Let `beta_i` be the associated alternating forms and
let `phi_i:Q to F_2` be quadratic forms with polar forms `beta_i`.  Put

```math
c(q)=1_{q\ne0},
```

and, in the first six `B` coordinates of `W`, define

```math
f(q)=\bigl(
\phi_1(q),c(q)+\phi_1(q),
\phi_2(q),c(q)+\phi_2(q),
\phi_3(q),c(q)+\phi_3(q)
\bigr),                                         \tag{GS.30}
```

with all remaining coordinates zero.  Identify `G=W direct-sum Q` and use
the canonical section `L_0(q)=(0,q)`; set

```math
S=B\cup\{L_0(q)+f(q):q\ne0\}.                  \tag{GS.31}
```

We first verify cycle contraction.  Let `R subseteq Q\{0}` have binary sum
zero and write `r=|R|`.  In the `i`th coordinate pair, the two output
parities are

```math
p_i=\sum_{q\in R}\phi_i(q),qquad
p_i+(r\bmod2).
```

If `r` is odd, exactly one coordinate in each pair is one, so the total
weight is three; every nonempty odd zero-sum set of distinct nonzero vectors
has `r>=3`.  If `r` is even and `r>=6`, the total weight is at most six.  The
only remaining case is `r=4`.  Such a set is an affine two-flat
`{a,a+u,a+v,a+u+v}` with independent `u,v`, and quadratic polarization gives
`p_i=beta_i(u,v)`.  Since `beta_3=beta_1+beta_2`, not all three `p_i` can be
one.  The output weight is therefore at most four.  Theorem GS.1 now says
that `B` is a shortest representation of `t=sum B`.

It is also diametral.  A point `(z,0)` uses at most `D` members of `B`.  For
nonzero `q`, the one-generator representation has length
`1+|z+f(q)|_B` and is at most `D` unless `z+f(q)=t`.  In that exceptional
case choose `q_1 notin {0,q}` and put `q_2=q+q_1`.  For each coordinate pair
in (GS.30), the triple defect has exactly one one because `c` has triple
defect one.  Hence

```math
|f(q)+f(q_1)+f(q_2)|_B=3.
```

The two quotient generators together with the `D-3` required members of
`B` represent `(z,q)` in `D-1` steps.  Thus the diameter is at most `D`, and
the point `t` shows equality.

Finally, nonsingularity of `beta_i` makes `phi_i` bent.  A direct Fourier
calculation gives

```math
\left|\sum_{q\in Q}(-1)^{\phi_i(q)+a(q)}\right|=2^{k/2}
```

for every linear Boolean `a`.  Thus the distance of `phi_i` from the nearest
linear function is `2^{k-1}-2^{k/2-1}`.  On `Q\{0}`, adding `c` complements
the truth table, so the corresponding nearest distance is
`2^{k-1}-2^{k/2-1}-1`.  Minimization of a vector-valued linear map separates
over its six coordinates.  Summing the three pairs proves (GS.28).
`square`

Since a supremum dominates an average, the same construction also gives

```math
\min_L\max_q|s_q+L(q)|_B
\ge3\left(1-2^{-k/2}-2^{-k}\right).
```

Thus geodesicity cannot force a vanishing uniform synchronization radius;
the present universal bounds leave only the constant-factor gap from three
to nine.

This sharp example is important conceptually.  The nonlinear information is
not carried by one bad scalar channel: six channels are coupled so that
every additive circuit obeys the joint contraction inequality.  A theorem
which tries to improve synchronization by charging channels separately
cannot work.

## 10. Director assessment: what is and is not new

1. **The cycle criterion itself is a reformulation.**  Theorem GS.1 is the
   exact optimality certificate for the chosen geodesic.  Calling it a new
   invariant would merely rename shortest-path duality.

2. **The synchronization consequence is substantive but uses known
   machinery.**  Theorem GS.4 was not present in the fibre count: a natural
   finite extremal hypothesis forces every point of an entire dense
   transversal to lie within nine of one homomorphism, independent of `D`
   and `k`.  Its proof is a genuinely joint application of Fourier/BLR
   stability, not a new replacement for BLR.

3. **Sharpness identifies the surviving nonlinear information.**  Three coupled bent
   pairs retain `3-o(1)` nonlinear bits per quotient point while satisfying
   every cycle inequality and preserving a diameter-`D` geodesic.  Thus
   geodesicity synchronizes to uniformly bounded distortion, but cannot
   force vanishing average distortion.  This also rules out a scalar-channel
   strengthening of the BLR step.

4. **Full fibre membership has a sharp worst-case bound at the quotient
   scale.**  Theorem GS.7 strips all but one lift per quotient generator at
   all-context cost at most `2h`, and Proposition GS.8 attains `2h` on a
   family of examples.  This is not an iff assertion for each individual
   source.
   Theorem GS.9 then compresses the retained quotient geometry to a linear
   section at cost `O(h)`.  This is a general group-extension response law,
   not merely an enumeration of syndrome supports.

5. **This reaches generative-theory level for a nontrivial class.**  The
   natural finite hypothesis “there is a geodesic whose projected generator
   set has diameter `h`” produces an `O(w^2)`-bit summary answering every
   future rooted word query to `O(h)` error.  For `h=o(w)` this is a strict
   macroscopic quotient.  The sharp `2h` example explains when microscopic
   information becomes macroscopic under composition.

6. **Exact closure remains open.**  The summaries can be substituted into
   any raw context and into a bounded union by telescoping, but the present
   map is not an exact monoid homomorphism.  Replacing `m` separately
   summarized sources incurs an `O(mh)` bound.  A closed product whose error
   depends only on the quotient diameter of the final union would be
   strictly stronger.

The strongest next statement suggested by this note is therefore:

> Construct, or rule out, an approximate summary-only composition law for
> geodesic charts whose error is controlled by the quotient diameter of the
> final union rather than by the number of summarized fragments.  A converse
> should determine whether large quotient diameter necessarily forces large
> all-context response complexity, beyond the sharp fibre-stripping example.
