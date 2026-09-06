# Conditioning selectors to exclude every sparse spin spectrum

Date: 2026-09-06. Director's exact-span observation, proved here with a
quantitative approximate version and its precise permanent limitation.
No unknown Boolean cap or general spectrum-type counting estimate is used.
The full proof has passed the independent audit agent's reconstruction.

Let O be ANY real orthogonal m-by-m matrix (in the weave take
`O=H^T/sqrt(m)`). Let T be uniform among k-subsets, p=k/m. A spin row is
the ternary vector `xi` with support EXACTLY T and nonzero entries ±1.

## 1. Exact sparsity: a finite counting theorem

For every integer s in [0,m],

```math
\Pr_T\{\exists\xi:\ \operatorname{supp}\xi=T,
\ \xi_T\in\{\pm1\}^T,
\ |\operatorname{supp}(O\xi)|\le s\}
\le\frac{\binom ms3^s}{\binom mk}.                        (1)
```

Proof. For each s-element frequency set U, the vectors with transformed
support in U lie in an s-dimensional real linear subspace E. Choose s
coordinate evaluations whose restriction to E is injective: such a set
exists by selecting linearly independent rows of a basis matrix for E.
Each ternary vector in E is determined by its values on these coordinates,
so E contains at most 3^s ternary vectors. Every vector of Fourier support
at most s lies in some E_U with |U|=s; no extra sum over smaller sizes is
needed. Each vector determines one selector. Divide this union count by
the number of selectors. This also includes the all-zero vector in the
upper count, harmlessly, although k>0 excludes it from the event.

Thus, if s/m->q and

```math
h(q)+q\log3<h(p),                                       (2)
```

the probability is exponentially small. At p15/16 the positive root
of equality in (2) is approximately 0.04520065110394. Any fixed smaller
q excludes, with high probability, sparse spectra for ALL 2^k spin rows
simultaneously. This is stronger than a statement about a typical spin.

## 2. Approximate sparsity: nets and a Hamming repair

Fix rho,delta>0, put `r=4(rho+delta)^2`, and assume r<2/3. Then

```math
\Pr_T\{\exists\xi,U:\ |U|\le s,
\ \operatorname{supp}\xi=T,\ \xi_T\in\{\pm1\}^T,
\ \|P_{U^c}O\xi\|_2\le\rho\sqrt m\}
\le\frac{\binom ms(1+2\sqrt p/\delta)^s}
             {\binom mk}
       \sum_{j\le\lfloor rm\rfloor}\binom mj2^j.        (3)
```

To prove it, fix U and project xi onto E=O^T R^U. Its projection has
norm at most sqrt(k). The radius-sqrt(k) ball of E has a delta-sqrt(m)
net of size at most `(1+2sqrt(p)/delta)^s`, by the elementary disjoint
ball volume argument. A net center is within `(rho+delta)sqrt(m)` of xi.
Round that center coordinatewise to a nearest member of {-1,0,1}, with
a deterministic tie rule. Every disagreement with xi costs at least
1/4 in squared distance from the UNROUNDED center. The rounded ternary
word and xi therefore differ in at most rm positions. A ternary Hamming
ball of that radius has exactly the upper count in (3): choose changed
positions and one of two new symbols at each. Union over nets and U,
then divide by binomial(m,k).

In particular, with s/m->q the logarithmic bad-probability rate is at
most

```math
h(q)+q\log(1+2\sqrt p/\delta)+h(r)+r\log2-h(p).          (4)
```

The restriction r<2/3 ensures the ternary Hamming-ball exponent is
`h(r)+r log2`; below this maximum the summands are exponentially increasing.

An explicit example at p15/16 is

```text
q=.02, rho=.01, delta=.02, r=.0036;
first four terms of (4) = .2160480307418113;
h(p)                    = .2337916587064593;
positive probability gap = .0177436279646480.
```

Thus a selector can be required, with probability tending exponentially
to one, to have the following uniform property:

```math
\|P_{U^c}O\xi\|_2^2>10^{-4}m
\quad\hbox{for EVERY supported spin xi and every }|U|\le.02m.       (5)
```

Integer floors and small adjustments in q do not affect the strict gap.
This excludes the finite-spike selectors in
`continued_feedback_weave_rare_selector_obstruction_2026_09_06.md` and
also a fixed neighborhood of such spectral concentration. It does not
exclude every low-entropy or approximately plateaued profile.

## 3. Conditioning is legitimate for the actual weave

Let G_m be any selector event supplied by (1) or (3), with
`P(G_m)=1-exp(-Omega(m))`. Choose the m fibre selectors independently
from the UNIFORM LAW CONDITIONED ON G_m. The original exact signing and
soft graph-contraction proof are unchanged; their averaged bound becomes

```math
2e^{t\gamma m^2}
 \left\{\frac{\mathbb E_T[Z_T(t)1_{G_m}]}{\Pr(G_m)}\right\}^{m}.
                                                               (6)
```

Since `m^{-1}log P(G_m)->0`, any strict exponential one-row upper bound
on the truncated expectation transfers at the SAME exponential rate.
This is a legitimate change of selector law for an existence construction,
not a proof that the remaining expectation satisfies that bound. In
particular the unconditional rare-selector obstruction is not an
obstruction to (6) once its selectors have been removed.

## 4. What the property alone gives for the permanent

There is a deterministic, but quantitatively weak, direct consequence.
Normalize the magnitudes by sqrt(k), writing a_j=|(H^T xi)_j|/sqrt(k);
then sum a_j²=m. Property (3)'s complement implies for every |U|<=s

```math
\sum_{j\notin U}a_j^2>\epsilon m,\qquad\epsilon=\rho^2/p.        (7)
```

For the normalized soft kernel put
`phi_t(x)=-log[(1+exp(-4tx))/2]/2`. This is increasing and concave on
[0,infinity), with phi_t(0)=0. The PSD kernel satisfies
`K(a,b)<=sqrt(K(a,a)K(b,b))`; hence every permutation product is at most
the product of the diagonal entries. Its square-root permanent is at
most `exp(-sum phi_t(a_j²))`.

After ANY allowed diagonal deletion, form U from that deleted coordinate
and the s-1 largest remaining magnitudes. Outside U, every a_j² is at
most `m/(s-1)`. Concavity gives

```math
\phi_t(x)\ge\frac{s-1}{m}\phi_t(m/(s-1))x
\quad(0\le x\le m/(s-1)).
```

Using (7), taking the maximum over deletions, and then s/m->q yields

```math
\frac1m\log L_t(a)
\le-\epsilon q\phi_t(1/q)+o(1).                         (8)
```

For the explicit parameters above and a moderate positive t, its decay
constant is only about `7.39357e-7`. Paying all row choices gives
`p log2-epsilon q phi_t(1/q)`, still very positive before the weave tilt.

There is a basic limitation to treating (7) as the only information:
the flat profile a_j=1 obeys it whenever epsilon<1-q, but has exactly

```math
L_t(a)=\left[(1+e^{-4t})/2\right]^{(m-1)/2}.
```

Its decay rate is at most `(log2)/2`, smaller than p log2 when p>1/2.
Thus an upper argument that merely maximizes L over admissible profiles
and then pays 2^k cannot succeed from anti-sparsity alone. One must also
count how MANY actual supported spin rows have each dangerous profile,
or use a stronger averaged inequality. No assertion is made here that
the flat profile is actually realized on every good selector.

The theorem removes a rigorous exceptional class at negligible selector
probability cost. Its remaining limitation is now quantitative and
explicit, not an implicit appeal to a general half-entropy count.

## 5. Which lower obstructions survive conditioning

The typical-Gaussian profile event has probability tending to one under
the joint uniform selector/spin law. Intersecting with ANY selector event
of probability tending to one preserves this fact. Consequently the
typical-profile lower bound, and its p<=0.9225232669 fixed-tilt obstruction,
survive the conditioning in Section 3 unchanged.

More generally, suppose a selector event A_m has probability
`exp(-I m+o(m))`, whereas the removed event has probability at most
`exp(-c m+o(m))`, with I<c. Then
`P(A_m intersect G_m)=P(A_m)(1-o(1))`. Any uniform lower bound on Z_T
over A_m therefore survives at the same exponential rate. Conversely,
there is no such protection when I>c: the finite-spike class has rate
h(p), larger than the explicit removal-rate guarantee .0177436, and is
indeed completely excluded by (5).
