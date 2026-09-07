# Global monotone balancing with operator-norm-compatible sign recovery

2026-09-07. New synthesis-track actual-sign construction; complete proof below,
pending independent audit. This strengthens the older same-order clique
repair in two respects: the final polarity gap is at most 2, and a fixed
normalized operator bound is preserved up to an additive constant.

Write `P(A)=max_x H_A(x)`, `R(A)=max_x -H_A(x)`, `Q(A)=max(P,R)`, with
`H_A(x)=x^T A x/2` for a hollow symmetric full signing. No efficient search for
the exact extrema is claimed. The baseline is the original ABSOLUTE cap Q,
not the half-width `(P+R)/2`.

## 1. Finite theorem

Orient A by overall negation if necessary so P=Q>=R and put Delta=P-R. If
Delta<=2, leave A unchanged. Otherwise choose a Boolean negative maximizer y,
so `H_A(y)=-R`, and let

```math
d={n\choose2},\qquad
\mathcal E=\{\{i,j\}: A_{ij}=y_iy_j\},\qquad
L=|\mathcal E|=\frac{d-R}{2}.                            \tag{1}
```

Suppose an integer T with `1<=T<=L` satisfies the following explicit test.
Define

```math
\begin{split}
u&=(n+1)\log2+\log\bigl(8(T+1)(L+1)\bigr),\\
E&=\sqrt{8Tu}+\frac43u,\\
z&=\log\bigl(16n(T+1)(L+1)\bigr),\\
F&=\sqrt{\frac{8T(n-1)z}{L}}+\frac43z,
\end{split}                                               \tag{2}
```

and assume

```math
2T>\Delta+E.                                             \tag{3}
```

Then there exists an actual hollow full signing B of the SAME order such that

```math
\boxed{\begin{split}
|P(B)-R(B)|&\le2,\\
Q(B)&\le Q(A)+E+2,\\
\#\{e:B_e\ne A_e\}&\le\min\left(T,\frac{\Delta+E+2}{2}\right),\\
\|B\|_{\rm op}&\le\|A\|_{\rm op}
                 +\frac{(n-1)T}{L}+F.
\end{split}}                                             \tag{4}
```

If the initial orientation was reversed, reverse the final B back as well.
The cap, norm, gap and the stated number of changed edges are unaffected.
Only initially eligible edges are changed; they may be distributed over the
entire old graph. No principal block or bounded-rank input is assumed.

## 2. Exact monotonicity of the negative extremum

Let S be ANY subset of eligible edges and let A^S be obtained by flipping
those edges. Then for every spin x,

```math
-H_{A^S}(x)=-H_A(x)+2\sum_{\{i,j\}\in S}y_iy_jx_ix_j
\le R+2|S|.
```

Equality is attained by x=y. Therefore

```math
R(A^S)=R+2|S|                                             \tag{5}
```

EXACTLY. The original negative maximizer y remains a negative maximizer
throughout every eligible-flip path, regardless of the ordering.

A single edge flip changes P by at most 2. All energy values have the same
parity d, so each change in P belongs to {-2,0,2}. Along a path flipping one
eligible edge at a time, (5) implies

```math
(P_{t+1}-R_{t+1})-(P_t-R_t)\in\{-4,-2,0\}.              \tag{6}
```

The gap is monotone. At its first nonpositive value, its value is at least
-2, because the preceding strictly positive gap was an even integer >=2.
Thus the first crossing has absolute polarity gap at most 2. This integrality
argument is why a one-shot path improves the older O(n)-gap repair.

## 3. Uniform random prefixes have favorable means

Choose a uniform random ordering of the L eligible edges. Let A_t be the
matrix after its first t edges are flipped, for `0<=t<=T`. Every prefix set
is a uniform t-element subset. With `p=t/L` and
`Y=offdiag(yy^T)=yy^T-I`, its mean is exactly

```math
\mu_t=(1-p)A-pY.                                        \tag{7}
```

For any Boolean x, `H_Y(x)=((y^T x)^2-n)/2>=-n/2`. Also

```math
P=Q(A)\ge\sqrt{\mathbb E_xH_A(x)^2}
=\sqrt d\ge n/2\qquad(n\ge2).
```

Consequently, for EVERY prefix size t,

```math
P(\mu_t)\le(1-p)P+pn/2\le P,\qquad
\|\mu_t\|_{\rm op}\le\|A\|_{\rm op}+p(n-1).             \tag{8}
```

The negative extreme of the mean is `(1-p)R+pd=R+2t`, with equality at y.
The operator bound uses the convex coefficient (1-p) on A; bounding the raw
matrix difference by an unpriced multiple of `||A||op` is unnecessary.

## 4. Bernoulli conditioning pays all prefix dependence

For a fixed t, generate independent Bernoulli(p) selectors on the eligible
edges, with p=t/L. Conditional on selecting exactly t edges their set is
uniform, hence their signing has precisely the law of A_t. The conditioning
event has probability at least `1/(L+1)`: t is a mode of Binomial(L,t/L), and
its L+1 probabilities sum to one. The endpoint cases t=0,L have probability
one and are harmless.

This comparison is made for each DETERMINISTIC prefix size. No assertion of
independence after conditioning is made. A union bound over t then controls
the whole random permutation path, even though different prefixes are highly
dependent.

For a fixed spin x in the unconditioned model, the centered energy is a sum
of independent centered summands of magnitude at most 2 and total variance
at most `4Lp(1-p)<=4t<=4T`. Scalar Bernstein gives

```math
\Pr\{|H_{A_t-\mu_t}(x)|>E\}
\le 2e^{-u}
```

before the conditioning factor. Using all `2^n` spins, the factor L+1, and
all T+1 prefix sizes, the probability that

```math
Q(A_t-\mu_t)\le E\quad\text{for EVERY }0\le t\le T       \tag{9}
```

fails is at most

```math
2^{n+1}(T+1)(L+1)e^{-u}=1/8.                             \tag{10}
```

For spectral control, each unconditioned matrix summand is
`-2A_ij(Bernoulli(p)-p)(e_i e_j^T+e_j e_i^T)`. It is independent, centered,
self-adjoint, and has operator norm at most 2. The sum of its expected
squares is diagonal, with largest entry at most

```math
4p(1-p)\max_i\deg_{\mathcal E}(i)
\le\frac{4T(n-1)}{L}.                                  \tag{11}
```

The self-adjoint matrix Bernstein inequality, applied to both signs, gives

```math
\Pr\{\|A_t-\mu_t\|_{\rm op}>F\}\le2n e^{-z}
```

before conditioning. Thus the simultaneous spectral event

```math
\|A_t-\mu_t\|_{\rm op}\le F\quad(0\le t\le T)            \tag{12}
```

fails with probability at most

```math
2n(T+1)(L+1)e^{-z}=1/8.                                  \tag{13}
```

There is therefore a random edge ordering satisfying BOTH (9) and (12).
Every bound is uniform before the later adaptive choice of a crossing time.

The only non-elementary concentration input is [Tropp, *User-Friendly Tail
Bounds for Sums of Random Matrices*, Theorem 6.1(ii)](https://tropp.caltech.edu/papers/Tro11-User-Friendly-preprint.pdf).
Its hypotheses and variance parameter were checked directly in the primary
author PDF, printed pages 23--24. Here dimension is n, summand bound is 2,
and (11) is the full matrix variance. The displayed square-root-plus-linear
threshold follows by substituting into the theorem's Bernstein denominator.
The scalar use has the identical elementary Bernstein substitution.

## 5. Stop at the first crossing

Fix an ordering on both good events. Equation (8) and (9) give
`P(A_t)<=P+E` for every t<=T. But (3) and (5) give

```math
R(A_T)=R+2T>P+E\ge P(A_T).
```

The first crossing t_* thus occurs by T. Its polarity gap is at most 2 by
(6). Immediately before crossing the positive extreme is the absolute cap,
so `Q(A_(t_*-1))<=P+E`; the crossing flip changes cap by at most 2. Therefore
`Q(A_(t_*))<=P+E+2`. Since the negative extreme at crossing equals its cap,
(5) also gives `R+2t_*<=P+E+2`, proving the edit bound in (4).

Finally apply (8) and (12) at this SAME selected t_* to obtain the operator
bound. Setting B=A_(t_*) proves every part of the finite theorem.

## 6. Uniform bounded-cap consequence

Fix C>0 and suppose `Q(A)<=C n^(3/2)`. Put

```math
w=\sqrt{n(\Delta+n)},\qquad
T=\left\lceil\frac\Delta2+64w\right\rceil.              \tag{14}
```

For all sufficiently large n depending only on C:

- `L=(d-R)/2>=n^2/8`, since R<=C n^(3/2);
- `T<=L`, since `T=O_C(n^(3/2))=o(n^2)`;
- the u in (2) is at most 2n, since T,L<=n^2 and all extra logarithms are
  O(log n);
- `nT<=66w^2`, using `w>=n` and `n Delta<=w^2`, so
  `E<=4sqrt(nT)+(8/3)n<40w`;
- `2T>=Delta+128w>Delta+E`, verifying the finite hypothesis.

Thus the actual B satisfies

```math
\boxed{\begin{split}
|P(B)-R(B)|&\le2,\\
Q(B)&\le Q(A)+O\!\left(\sqrt{n\Delta}+n\right),\\
\#\{e:B_e\ne A_e\}
 &\le\Delta/2+O\!\left(\sqrt{n\Delta}+n\right),\\
\|B\|_{\rm op}&\le\|A\|_{\rm op}+O_C(\sqrt n).
\end{split}}                                             \tag{15}
```

The cap/edit constants may be absolute once n is sufficiently large for the
displayed input class. For the norm, (2) and `L>=n^2/8` give more explicitly

```math
\|B\|_{\rm op}-\|A\|_{\rm op}
\le \frac{8T}{n}
 +8\sqrt{\frac{Tz}{n}}+\frac43z,
\qquad z=O(\log n).                                    \tag{16}
```

Since T=O_C(n^(3/2)), the first term is O_C(sqrt n) and the other two are
`O_C(n^(1/4)sqrt(log n)+log n)=o_C(sqrt n)`.

In particular, at every order an exact minimizer admits a same-order
selectable near-minimizer with gap <=2 and cap excess O(n^(5/4)). The norm
increment is controlled even when the original normalized norm diverges;
fixed bounded normalized norm is preserved if it was initially present.

## 7. A simultaneous balanced, spectrally regular near-minimizer normal form

The previously audited regular-core/refill theorem gives, from an exact
minimizer at order n and for each fixed epsilon in (0,1/2], an actual same-order
signing A_epsilon with

```math
Q(A_\epsilon)\le M_n+2\sqrt\epsilon\,n^{3/2},\qquad
\|A_\epsilon\|_{\rm op}
 \le (4K_GC/\epsilon+8)\sqrt n,
```

where C is any fixed eventual upper bound for M_n/n^(3/2). Apply (15) to this
already regularized input. Its cap bound is uniform in epsilon<=1/2, so the
new norm increment and n^(5/4) cost are uniform over this epsilon interval.
The resulting B_epsilon,n satisfies simultaneously

```math
\boxed{\begin{split}
|P(B_{\epsilon,n})-R(B_{\epsilon,n})|&\le2,\\
Q(B_{\epsilon,n})&\le M_n+2\sqrt\epsilon\,n^{3/2}+O_C(n^{5/4}),\\
\|B_{\epsilon,n}\|_{\rm op}
 &\le (4K_GC/\epsilon+O_C(1))\sqrt n.
\end{split}}                                             \tag{17}
```

This compatibility was not supplied by the older clique repair: its planted
block could have operator size n^(3/4). Equation (17) says that polarity
balance can be imposed for free at the leading scale inside fixed-accuracy
spectrally regular models. It does NOT say that one fixed operator constant
captures an o(1)-near-minimizing sequence as epsilon decreases to zero.

## 8. Microcanonical counting consequence

The transformation in (15) changes O_C(n^(3/2)) edges for every input of cap
at most Cn^(3/2). Choose one admissible output deterministically for each
input. Any fixed output has at most

```math
\sum_{j\le h_n}{d\choose j}
=\exp\bigl(O_C(n^{3/2}\log n)\bigr),
\qquad h_n=O_C(n^{3/2}),
```

possible preimages. Therefore the number of signings with cap at most Cn^(3/2),
and optionally with normalized operator norm at most L_0, is at most
`exp(O_C(n^(3/2)log n))` times the number with cap at most
`Cn^(3/2)+O_C(n^(5/4))`, gap<=2, and normalized operator norm at most
`L_0+O_C(1)`. The leading n^2-speed disorder counting entropy is unchanged
by this balanced repair. This is a deterministic preimage count, with no
unproved independence of selected witnesses or masks.

## 9. Finite checks

`computations/principle_synthesis_2026_09_07_global_balancing_check.py` checks
90 complete eligible-edge paths at orders 3 through 12, including 68 with
nonzero initial polarity gap. Seed 20260907 fixes all random matrices and
edge orders. Integer arithmetic verifies (1), (5), (6), every prefix's exact
mean energy and negative extremum, the favorable positive mean, and the
gap<=2 first-crossing conclusion. The realized prefix discrepancy is kept as
an exact Fraction and verifies the pathwise cap bound.

The saved full record is
`computations/results/principle_synthesis_2026_09_07_global_balancing_check.json`.
Nine large-n parameter substitutions additionally check the chosen finite
hypotheses numerically; they are explicitly NOT signing realizations or
asymptotic certificates. The theorem rests on the preceding analytical
inequalities, not on those finite tests. The checker deliberately follows
some initial gap-2 paths even though the theorem may leave them unchanged.

## 10. Exact scope and distinction from width recovery

The construction is valid for EVERY bounded-cap input, but only guarantees a
SELECTABLE modified near-minimizer when starting from an exact minimizer.
It does not prove that every original exact minimizer was balanced, nor that
the output remains exactly minimizing. All old edges are available for changes;
the proof chooses a sparse-in-density, globally distributed subset.

Most importantly, it raises the smaller polarity until it meets the larger
one, while controlling the larger. Its reference is P=Q(A), NOT
`W(A)=(P+R)/2`. Thus this does not convert width minimizers to original-cap
minimizers, does not prove `M_n-W_n=o(n^(3/2))`, and does not settle the old
width-to-absolute gap. Nor does it construct the desired weighted-child
replacement at a new order. The new content is an actual stronger normal-form
operation that is compatible with the spectral reductions already proved.
