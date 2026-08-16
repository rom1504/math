# Coding geometry audit: pair data, rooted holes, and a parity hierarchy

Date: 2026-08-16

Status: independent specialist report. This file does not modify a final
deliverable or any main state file.

## 1. Conclusions

Let \(Q_n=\mathbb F_2^n\), let \(d\) be Hamming distance, and put
\(q(x,y)=n-2d(x,y)\). For a code \(A\subseteq Q_n\), use the membership
landscape

\[
H_A(x)={\bf 1}_{\{x\in A\}}.
\]

For

\[
C=\{0,3,5,6\},\qquad D=\{0,3,5,9\}\subset Q_4,
\]

the following statements are exact.

1. Both codes are four-point equilateral sets: every two distinct codewords
   are at distance \(2\). Their ordered internal distance polynomial is
   \(I_C(z)=I_D(z)=4+12z^2\).

2. Equality is stronger than the repository experiment explicitly checks.
   The **entire ambient ordered census** of
   \((H_A(x),H_A(y),q(x,y))\), with \(x,y\in Q_4\), is identical for \(C\)
   and \(D\). The same is true in every Cartesian power, both for hard
   product membership and for additive block-membership energy.

3. Their covering radii are
   \(\rho(C)=2\) and \(\rho(D)=3\). Consequently
   \(\rho(C^m)=2m\) and \(\rho(D^m)=3m\), leaving normalized gap \(1/4\).

4. A labeled strong field separates the membership landscapes immediately.
   A code-forcing reward followed by a rooted field extracts the nearest-code
   distance and hence the covering-radius gap. These are two complementary
   query regimes, made precise below.

5. The smallest repair is query-dependent. For a fixed root \(z\), nearest
   code response needs exactly \(\delta_A(z)=d(z,A)\), not a new unrooted
   moment. For every radial linear query at \(z\), the minimal state is the
   upper convex hull of
   \(\{(q(z,x),H_A(x)):x\in Q_n\}\). The full rooted distance enumerator is a
   convenient sufficient statistic but is generally not minimal.

6. Full unrooted triple data do **not** repair the failure. In fact the full
   membership-and-distance distributions of \(C\) and \(D\) agree through
   four points and first differ at five points. An infinite parity-coset
   hierarchy shows that, for every fixed \(k\), global data through order
   \(k\) can miss a different covering radius.

The obstruction is rooted and extrinsic. It cannot be repaired universally
by merely moving from pairs to triples.

## 2. Exact pair theorem, including the full ambient energy census

For \(A\subseteq Q_n\), define the ordered inner distribution

\[
I_A(j)=\#\{(a,b)\in A^2:d(a,b)=j\}.
\]

For \(\epsilon,\eta\in\{0,1\}\), also define

\[
N_A^{\epsilon\eta}(j)
=\#\{(x,y)\in Q_n^2:H_A(x)=\epsilon,\ H_A(y)=\eta,\ d(x,y)=j\}.
\]

### Theorem 2.1 (inner distribution determines the full two-point membership census)

If \(|A|=s\), then

\[
\begin{aligned}
N_A^{11}(j)&=I_A(j),\\
N_A^{10}(j)=N_A^{01}(j)&=s\binom nj-I_A(j),\\
N_A^{00}(j)&=(2^n-2s)\binom nj+I_A(j).
\end{aligned}                                                   \tag{2.1}
\]

Consequently, equal size and equal ordered inner distribution imply equality
of the complete ambient histogram of
\((H(x),H(y),d(x,y))\), equivalently of
\((H(x),H(y),q(x,y))\).

**Proof.** From each fixed \(a\in A\), exactly \(\binom nj\) ambient words
are at distance \(j\). Removing the codewords among them gives the second
line. There are \(2^n\binom nj\) ambient ordered pairs at distance \(j\);
subtracting the other three membership patterns gives the final line.
No linearity assumption is used. \(\square\)

For \(C,D\), the complete common table is

| \(j\) | \(N^{11}\) | \(N^{10}=N^{01}\) | \(N^{00}\) |
|---:|---:|---:|---:|
| 0 | 4 | 0 | 12 |
| 1 | 0 | 16 | 32 |
| 2 | 12 | 12 | 60 |
| 3 | 0 | 16 | 32 |
| 4 | 0 | 4 | 8 |

Thus the collision really is a complete energy--energy--overlap collision
for the membership landscapes, not merely equality of a code-only distance
enumerator.

### Theorem 2.2 (Cartesian powers)

For the hard product code \(A^m\subseteq Q_{nm}\),

\[
I_{A^m}(z)=I_A(z)^m,\qquad I_A(z):=\sum_j I_A(j)z^j.     \tag{2.2}
\]

Hence

\[
I_{C^m}(z)=I_{D^m}(z)
=(4+12z^2)^m
=4^m\sum_{r=0}^m\binom mr3^r z^{2r}.                   \tag{2.3}
\]

Applying (2.1) with \(n=4m\) and \(s=4^m\) proves equality of the **full**
ambient pair membership/overlap census in every power.

There is a second common convention. If

\[
\bar H_A^{(m)}(x_1,\ldots,x_m)=\sum_{b=1}^m H_A(x_b),
\qquad q_m=\sum_{b=1}^m q(x_b,y_b),
\]

then the joint generating function of
\((\bar H(x),\bar H(y),q_m(x,y))\) is the \(m\)-th power of the base joint
generating function. The additive-energy convention therefore also gives
exact equality in every Cartesian power.

The coefficients in (2.3) agree exactly with the experiment JSON: for
example \(m=2\) gives \(16,96,144\), and \(m=5\) gives
\(1024,15360,92160,276480,414720,248832\).

## 3. Covering radius and exact field separation

Write a word as a sign vector

\[
\sigma(x)=((-1)^{x_1},\ldots,(-1)^{x_n}),\qquad
\langle\sigma(x),\sigma(z)\rangle=n-2d(x,z).            \tag{3.1}
\]

### 3.1 The base radii and tensor law

In the displayed bit order,

\[
C=\{0\}\times\{u\in Q_3:\operatorname{parity}(u)=0\}.
\]

The even-parity code in \(Q_3\) has radius \(1\), while the fixed first
coordinate contributes at most one more. Both contributions can occur, so
\(\rho(C)=2\).

For \(D\), the word \(z_*=1110_2=14\) is at distance \(3\) from every
codeword:

\[
(d(z_*,a):a\in D)=(3,3,3,3).
\]

No word can be farther than \(3\): a word of weight at most \(3\) is within
three of \(0\in D\), and \(1111\) is at distance \(2\) from every nonzero
word of \(D\). Thus \(\rho(D)=3\).

For products, distance separates blockwise:

\[
d((x_1,\ldots,x_m),A^m)=\sum_{b=1}^m d(x_b,A).
\]

Maximizing each block independently proves

\[
\rho(A^m)=m\rho(A).                                    \tag{3.2}
\]

This proves the reported \(2m\) versus \(3m\) radii rather than merely
extrapolating the base computation.

### 3.2 A genuinely strong external field

For \(h>1/2\), set

\[
V_A(h;z)=\max_{x\in Q_n}
\{H_A(x)+h\langle\sigma(x),\sigma(z)\rangle\}.
\]

Every \(x\ne z\) loses at least \(2h>1\) in field energy and can gain at
most one unit of membership energy. Therefore

\[
V_A(h;z)=hn+H_A(z).                                    \tag{3.3}
\]

At the common labeled root \(z=6\), \(6\in C\setminus D\), and hence

\[
V_C(h;6)=4h+1,\qquad V_D(h;6)=4h.
\]

Thus the identical global pair census does not determine this elementary
labeled strong-field response. Equation (3.3) also makes the information
cost transparent: responses to all sufficiently strong labeled fields recover
the complete membership function.

### 3.3 A code-forcing reward extracts the holes

To interrogate covering rather than membership at the root, fix \(h>0\) and
introduce a hard code constraint, or its finite version

\[
V_A^{B,h}(z)=\max_x
\{B H_A(x)+h\langle\sigma(x),\sigma(z)\rangle\},
\qquad B>2hn.                                           \tag{3.4}
\]

The bound on \(B\) forces every maximizer into \(A\), so

\[
V_A^{B,h}(z)=B+h(n-2\delta_A(z)),
\qquad \delta_A(z):=d(z,A).                             \tag{3.5}
\]

At \(z_*=14\),

\[
\delta_C(z_*)=1,\qquad \delta_D(z_*)=3,
\]

and (3.5) gives \(B+2h\) versus \(B-2h\). At
\(z_*^{(m)}=(14,\ldots,14)\), the gap is \(4hm\). Taking the worst root gives

\[
\min_z V_A^{B,h}(z)=B+h(n-2\rho(A)),                    \tag{3.6}
\]

so for \(C^m,D^m\) the worst-root values are \(B\) and \(B-2hm\).

## 4. The smallest rooted repair

There are three useful levels, and they should not be conflated.

### 4.1 Natural sufficient table

For a retained root label \(z\), define

\[
R_H(z;e,j)=\#\{x:H(x)=e,\ d(x,z)=j\}.                  \tag{4.1}
\]

This rooted energy--distance enumerator determines every radial response

\[
V_H(h;z)=\max_{e,j:R_H(z;e,j)>0}\{e+h(n-2j)\}.          \tag{4.2}
\]

For a hard linear-code constraint, \(R_C(z;1,j)\) counts the weight-\(j\)
words of the coset \(z+C\).

### 4.2 Minimal statistic for all radial linear responses

Counts in (4.1) are irrelevant to a zero-temperature maximum. The exact
minimal quotient for all \(h\in\mathbb R\) is

\[
K_H^+(z)=\operatorname{upper\,hull}\operatorname{conv}
\{(n-2d(x,z),H(x)):x\in Q_n\}.                          \tag{4.3}
\]

Its support function is (4.2), and convex duality recovers the upper hull
from all responses. This is the root-specialized upper response roof from
query_response_body.md.

### 4.3 Minimal statistic for this covering failure

Under a hard code constraint, or in regime (3.4), (4.3) collapses for
\(h>0\) to its rightmost code point. The minimal per-root statistic is exactly

\[
\delta_A(z)=\min\{j:R_A(z;1,j)>0\}.                    \tag{4.4}
\]

It is necessary as well as sufficient, since (3.5) recovers it from any one
fixed \(h>0\). If only the covering radius is queried and root labels are not,
the minimal scalar is \(\max_z\delta_A(z)=\rho(A)\). If all labeled
nearest-code queries are allowed, the labeled map \(z\mapsto\delta_A(z)\)
must be retained. That map is lossless for the code itself, since
\(A=\{z:\delta_A(z)=0\}\); the declared labeled interface leaves no universal
compression here.

The complete rooted profiles make the collision visible:

\[
\begin{array}{c|c|c}
&\text{sorted distances to the four codewords}&
\text{number of roots}\\ \hline
C&(0,2,2,2),(1,1,1,3),(1,3,3,3),(2,2,2,4)&4,4,4,4\\
D&(0,2,2,2),(1,1,1,1),(1,1,3,3),(2,2,2,4),(3,3,3,3)
&4,1,6,4,1
\end{array}
\]

In particular, the distance-layer counts at distances \(0,1,2,3\) are

\[
C:(4,8,4,0),\qquad D:(4,7,4,1).                        \tag{4.5}
\]

Pair averaging erased exactly this root-to-code information.

## 5. Triples do not suffice

Define the full unrooted \(t\)-point census

\[
T_t(A)=\operatorname{hist}\left(
(H_A(x_i))_{i=1}^t,(d(x_i,x_j))_{1\le i<j\le t}
\right),\quad (x_1,\ldots,x_t)\in Q_n^t.               \tag{5.1}
\]

This is stronger than recording only \(t\) codewords: it records every
membership pattern and the complete ambient distance matrix.

### Proposition 5.1 (the base pair agrees through order four)

\[
T_t(C)=T_t(D)\qquad(1\le t\le4),                       \tag{5.2}
\]

but \(T_5(C)\ne T_5(D)\).

**Proof.** A count with some tuple positions constrained to be codewords is
controlled by the embedded configuration of the distinct constrained
codewords. One codeword is handled by cube transitivity; two by the
intersection numbers of the Hamming scheme. Every three distinct codewords
in either code form a binary equilateral triangle of side \(2\), and all such
triangles are cube-isometric. If four positions are all distinct codewords,
there is no remaining ambient point and their recorded distance matrix is the
same equilateral matrix. Exact zero/one membership patterns follow from
inclusion--exclusion. This proves (5.2).

For separation at five points, let the first four entries list the four
distinct codewords and let the fifth be a noncodeword. In \(D\), root \(14\)
is at distance \(3\) from all four codewords. Hence there are \(4!=24\)
ordered tuples with energy pattern \((1,1,1,1,0)\), all six internal
distances \(2\), and all four root distances \(3\). There are none for \(C\),
because \(\rho(C)=2\). \(\square\)

Thus complete triple data do not suffice; complete four-point data do not
suffice either. What succeeds at order five is effectively “all four code
points plus a retained external root,” exactly the missing interface.

### Theorem 5.2 (an infinite fixed-order hierarchy)

For every fixed \(k\), there are binary codes \(A,B\) whose full censuses
\(T_t\) agree for every \(t\le k\), but whose covering radii differ.

**Construction and proof.** Let \(r=2s+1\) be odd with \(s\ge1\), put \(M=r+1\) and
\(N=2^{r-1}\). For \(\epsilon\in\{0,1\}\), index the \(N\) coordinates by

\[
P_\epsilon=\{v\in\mathbb F_2^r:|v|\equiv\epsilon\pmod2\}.
\]

Define \(A_\epsilon\subset Q_N\) to consist of the zero word and the \(r\)
coordinate-function rows

\[
a_i(v)=v_i,\qquad 1\le i\le r.                         \tag{5.3}
\]

Each pair of distinct words is at distance \(N/2\), for both parity choices.
More strongly, every proper selection of at most \(r\) distinct codewords has
the same embedded column-pattern multiset for \(\epsilon=0\) and \(1\). To
see this, translate one selected word to zero. The remaining patterns are the
values of a space \(W\) of linear forms on \(v\). The parity form
\({\bf1}\cdot v\) is not in \(W\): if the zero label was selected, a basis
label is missing; if it was not, \(W\) consists of even-weight combinations,
whereas \(r\) is odd. Therefore every fiber of the selected-pattern map meets
the two parity hyperplanes equally.

It follows, by the same ambient-completion and inclusion--exclusion argument,
that

\[
T_t(A_0)=T_t(A_1)\qquad(1\le t\le M).                  \tag{5.4}
\]

It remains to prove that the radii differ. Define

\[
S_\epsilon=
\sum_{\substack{0\le w\le r\\w\equiv\epsilon\ (2)}}
\binom rw\max\{w,M-w\}.                                \tag{5.5}
\]

For any root \(z\), averaging its distances over the \(M\) codewords gives

\[
\min_{a\in A_\epsilon}d(z,a)
\le {1\over M}\sum_{a\in A_\epsilon}d(z,a)
\le {S_\epsilon\over M}.                               \tag{5.6}
\]

Let \(p=s\bmod2\). On coordinates of weight \(w\equiv p\pmod2\), choose
the root bit opposite the strict majority of the \(M\) codeword bits. There
is no tie, since the tie weight \(s+1\) has parity \(1-p\). The resulting
sets of codeword labels disagreeing with the root are exactly all parity-\(p\)
subsets of \([M]\) of size \(>M/2\), one from each complementary pair. This
family is invariant under every permutation of the \(M\) labels, so all \(M\)
distances are equal to \(S_p/M\). Thus

\[
\rho(A_p)=S_p/M.                                       \tag{5.7}
\]

The elementary alternating-binomial identity

\[
S_p-S_{1-p}=\binom{2s}{s}                              \tag{5.8}
\]

follows by pairing terms \(w\) and \(2s+1-w\) in
\(S_0-S_1\): each pair contributes
\((-1)^w\binom{2s+1}{w}\), and
\(\sum_{w=0}^s(-1)^w\binom{2s+1}{w}
=(-1)^s\binom{2s}{s}\).
Together with (5.6), it gives

\[
\rho(A_{1-p})
\le\left\lfloor {S_{1-p}\over M}\right\rfloor
< {S_p\over M}=\rho(A_p).                              \tag{5.9}
\]

Choosing odd \(r\) with \(r+1\ge k\) proves the theorem. Moreover, the first
possible distinction is at order \(M+1\): list all \(M\) codewords and append
a deepest ambient root. The different radii force different signatures.

For \(r=3\), the even columns
\(000,011,101,110\) produce \(C\), while the odd columns
\(001,010,100,111\) produce \(D\). The original four-bit example is the
first member of this hierarchy, not an isolated coincidence.

## 6. Relation to primary coding-theory literature

- **Distance and weight enumerators.** MacWilliams' original theorem shows
  how the weight spectrum of a linear code determines that of its dual
  ([MacWilliams 1963](https://doi.org/10.1002/j.1538-7305.1963.tb04003.x)).
  Delsarte placed the inner distribution of arbitrary codes in the Hamming
  association scheme and related its transform to external parameters
  ([Delsarte 1973 thesis](https://dial.uclouvain.be/pr/boreal/object/boreal%3A205698),
  [Delsarte 1973 paper](https://doi.org/10.1016/S0019-9958(73)80007-5)).
  Equation (2.1) is the elementary two-color specialization of that viewpoint.

- **Joint and higher enumerators.** MacWilliams, Mallows, and Sloane
  introduced the joint weight enumerator and its MacWilliams transform
  ([1972](https://doi.org/10.1109/TIT.1972.1054898)). Wei's generalized
  Hamming weights minimize the union support of an \(r\)-dimensional linear
  subcode ([1991](https://doi.org/10.1109/18.133259)); Kløve developed support
  weight distributions ([1992](https://doi.org/10.1016/0012-365X(92)90559-X));
  Britz related higher support enumerators to matroid polynomials
  ([2002](https://doi.org/10.37236/1636)). These refinements retain coordinate
  support information that an ordinary radial distance enumerator discards.
  Concretely, \(C\) is a degenerate linear \([4,2,2]\) code with generalized
  weights \(d_1=2,d_2=3\), while \(D\) is nonlinear. A tuple-support
  enumerator can see that the union of supports of the three nonzero words is
  \(3\) for \(C\) and \(4\) for \(D\). This is useful coordinate-resolved
  information, but it is not outer/coset data: even all intrinsic metric
  \(t\)-point distributions of the equilateral sets \(C,D\) agree.

- **Covering radius and cosets.** Delsarte's four-parameter paper gives the
  external-distance bound on covering radius. Goethals and van Tilborg's
  original uniformly packed-code paper studies controlled outer distributions
  ([1975](https://research.tue.nl/en/publications/uniformly-packed-codes/)).
  For a linear code, the least weight in coset \(z+C\) is \(\delta_C(z)\),
  and the maximum coset-leader weight is the covering radius; coset weight
  enumerators are the classical rooted object. See Schatz's original note
  ([1980](https://doi.org/10.1080/00029890.1980.11995087)).

- **Completely regular codes.** Neumaier's definition makes the distance
  partition equitable
  ([1992](https://doi.org/10.1016/0012-365X(92)90565-W)). Equivalently in a
  distance-regular graph, the rooted weight distribution depends only on the
  root's distance layer. This is a positive condition under which the rooted
  table compresses. Neither present code is completely regular: the
  distance-one layer in the profile table above contains multiple profiles.

- **Homometric sets.** Patterson introduced homometric structures in the
  diffraction setting ([1939](https://www.nature.com/articles/143939b0));
  Rosenblatt and Seymour give an algebraic structure theorem in terms of
  Patterson functions/difference multisets
  ([1982](https://doi.org/10.1137/0603035)). Terminology matters here:
  \(C,D\) are homometric only after **radializing** differences to Hamming
  weight. They are not homometric as subsets of \(\mathbb F_2^4\). The exact
  XOR autocorrelation of \(C\) is \(4\) on \(C\) and zero elsewhere, while
  that of \(D\) is \(4\) at zero and \(2\) on
  \(3,5,9,6,10,12\).

The literature points to the same dividing line as the direct proof: inner
enumerators describe code-code averages; covering is an outer, rooted/coset
extreme. Complete regularity can make the latter compressible, but it does
not make unrooted inner data universally sufficient.

## 7. Comparison with the repository artifacts

This comparison was made only after the direct derivations and literature
check above.

1. entropy_overlap_lab.py computes the ordered internal enumerator correctly
   and tensors it by polynomial convolution. It computes the base radii
   exactly. Its product radii are obtained by multiplying the base radii
   rather than enumerating product cubes; equation (3.2) supplies the proof.

2. The experiment does **not** explicitly form the ambient
   \((H(x),H(y),q(x,y))\) histogram for the membership landscape. Theorem 2.1
   proves that stronger equality and removes any ambiguity in “complete pair
   data.”

3. entropy_overlap_results.json agrees coefficient-for-coefficient with
   (2.3) through all five reported powers and with radii \(2m,3m\).

4. Using the repository virtual environment, importing code_collision and
   rerunning it through power five reproduced the saved enumerators and radii.
   A separate standard-library enumeration reproduced the full ambient base
   table, rooted profiles, equality of \(T_t\) for \(t\le4\), and the explicit
   order-five separation.

5. The rooted upper hull (4.3) is consistent with the response-roof duality
   in query_response_body.md. The coding example sharpens its interpretation:
   declaring the root interface converts an averaged pair invariant into a
   coset/outer-distribution problem, and for nearest-code response the roof
   reduces exactly to the distance-to-code map.

## 8. Recommended theorem-level formulation

> The complete global two-replica membership-energy/overlap law, even under
> every Cartesian power, is not sufficient for labeled linear-field response
> or for a rooted covering extreme. The missing minimal datum for the
> code-forced radial query is the labeled distance-to-code function; for all
> radial fields it is the labeled upper rooted response hull. No fixed order
> of unrooted replica data repairs this universally.

The last sentence is supported by Theorem 5.2. The four-bit pair gives a
particularly transparent instance: pair, triple, and four-point global data all
collide, the rooted statistic separates, and the first unrooted witness
appears only when all four codewords and one external root can occur together.
