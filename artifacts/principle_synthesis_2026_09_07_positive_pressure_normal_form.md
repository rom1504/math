# Positive-pressure compatibility of global balancing; exact-minimizer boundary

2026-09-07. Continuation of the new global balancing theorem. The first part
is a proved actual-sign refinement. The final part records an exact optimizer
consequence and why it does not currently close the original limit.

Let A be hollow full signing, oriented so `P(A)=Q(A)>=R(A)`. Use the eligible
negative-ground edge set, random-prefix path, parameters T,L,E,F,u,z, and good
event from `principle_synthesis_2026_09_07_global_balancing.md`. Assume its
finite condition `2T>Delta+E`. The good event has probability at least 3/4.
For beta>=0 put

```math
Z_A^+(\beta)=\sum_x\exp\bigl(\beta H_A(x)/\sqrt n\bigr).
```

## 1. Finite simultaneous positive-pressure refinement

For any finite list `beta_1,...,beta_r` in `[0,sqrt(n)/2]`, the output B can be
chosen to satisfy ALL the cap, gap, edit and spectral guarantees of the global
balancing theorem and, simultaneously for j=1,...,r,

```math
\boxed{\quad
\log Z_B^+(\beta_j)\le\log Z_A^+(\beta_j)
 +\frac{\beta_j T\sqrt n}{2L}
 +\frac{3\beta_j^2T}{n}
 +\log\!\left(\frac43r(T+1)(L+1)\right).
\quad}                                                     \tag{1}
```

When `Q(A)<=Cn^(3/2)`, the chosen `T=O_C(n^(3/2))` and `L>=n^2/8` give

```math
\frac{\log Z_B^+(\beta_j)-\log Z_A^+(\beta_j)}n
\le O_C(\beta_j^2n^{-1/2})
      +O_C(\beta_j/n)+O(\log(rn)/n).                    \tag{2}
```

In particular every fixed finite temperature list can be retained with a
vanishing pressure-per-vertex error by one and the same spectrally compatible,
gap-at-most-two actual output. This does not assert preservation of the
negative-pressure profile.

### Proof of (1)

The unconditioned Bernoulli model at deterministic prefix size t has mean
`mu_t=(1-p)A-pY`, where `p=t/L` and `Y=offdiag(yy^T)`. For a fixed spin, each
centered eligible-edge contribution xi has |xi|<=2. The standard scalar
Bernstein moment bound gives, for `lambda=beta/sqrt(n)<=1/2`,

```math
\log\mathbb E\exp\!\left(\lambda\sum_e\xi_e\right)
\le\frac{\lambda^2\sum_e\mathbb E\xi_e^2}
             {2(1-2\lambda/3)}
\le3\beta^2t/n,                                        \tag{3}
```

since the variance sum is at most 4t. The signs of the spin features do not
affect that variance or magnitude bound.

Pointwise `H_(mu_t)(x)<=(1-p)H_A(x)+pn/2`. Moreover the function
`s -> log Z_A^+(s)` is convex, has derivative zero at s=0 because H_A has
uniform mean zero, and hence is nondecreasing for s>=0. Therefore

```math
Z_{\mu_t}^+(\beta)
\le e^{\beta p\sqrt n/2}Z_A^+((1-p)\beta)
\le e^{\beta p\sqrt n/2}Z_A^+(\beta).                   \tag{4}
```

Sum (3) over all spins and use (4). Conditioning on the number of Bernoulli
selectors being exactly t costs at most L+1, as proved in the global theorem.
Thus under the actual random-prefix law,

```math
\mathbb E Z_{A_t}^+(\beta)
\le (L+1)
 \exp\!\left(\frac{\beta T\sqrt n}{2L}
             +\frac{3\beta^2T}{n}\right)Z_A^+(\beta).
                                                               \tag{5}
```

For an arbitrary first-crossing stopping time tau on the good event, positivity
alone gives

```math
Z_{A_\tau}^+(\beta)\mathbf1_{\rm good}
\le\sum_{t=0}^T Z_{A_t}^+(\beta).
```

Divide each expression by its corresponding upper bound in (5), and sum over
the r declared temperatures. Conditional expectation on the good event is at
most `(4/3)r(T+1)`. Some ordering in that event therefore has this bound on the
sum of all r nonnegative ratios, hence on each individual ratio. Taking logs
proves (1). No optional-stopping identity or conditional independence at tau
has been asserted.

## 2. An actual balanced one-sided soft-max formulation of the original target

Define the constrained one-sided free energy

```math
F_n^{\rm bal}(\beta)=
\min_{B:\,|P(B)-R(B)|\le2}
 \frac1{\beta n}\log Z_B^+(\beta),\qquad\beta>0.
```

The new global balancing theorem makes the constraint feasible at every
sufficiently large order with only O(n^(5/4)) extra cap. The elementary
soft-max inequalities give

```math
\boxed{\quad
\frac{M_n}{n^{3/2}}-\frac2{n^{3/2}}
\le F_n^{\rm bal}(\beta)
\le\frac{M_n}{n^{3/2}}+O(n^{-1/4})+\frac{\log2}{\beta}.
\quad}                                                     \tag{6}
```

Thus all-order convergence of these BALANCED one-sided pressures for every
fixed beta would imply original convergence, by sending beta to infinity
after the dimension limit. This is a sufficient formulation, not a claim
that it has been proved or made intrinsically easier. The unrestricted
one-sided minimum is a different and essentially trivial quantity; the
balance constraint in (6) cannot be dropped.

The operator-compatible theorem also permits a fixed-accuracy spectral cutoff
before this balanced reduction, with the explicit `2sqrt(epsilon)` cap
payment already recorded in the global theorem. This does not prove closure
of the balance constraint under a Gaussian block interpolation.

## 3. Exact-minimizer high-state hitting condition

Suppose now that A itself is EXACT minimizing, with `P=M_n` and
`Delta=P-R>2`. For every integer `0<=t<Delta/2` and EVERY eligible t-subset S,
equation `R(A^S)=R+2t<P` and original minimality force

```math
P(A^S)\ge P.
```

Equivalently, for every such S there is a positive-side spin x satisfying

```math
H_A(x)-2\sum_{e=\{i,j\}\in S}y_iy_jx_ix_j\ge P.          \tag{7}
```

This is an exact all-subsets statement on an actual optimizer, not merely
stationarity under one edge. The corresponding general near-minimizer version
with excess eta requires `R+2t<M_n` and replaces the right side by M_n.

The statement has a clear counting consequence but no leading contradiction
at present. If `a_x=P-H_A(x)` and p=t/L, then the centered deviation required
in the Bernoulli comparison is at least

```math
(1-p)a_x+p(P-n/2).                                      \tag{8}
```

Even for ground states this is only of order `t/sqrt(n)` when
P=Theta(n^(3/2)); its scalar concentration exponent is only of order t/n.
At the largest possible t=O(n^(3/2)), this is order sqrt(n), whereas the
Boolean spin entropy is order n. Thus a direct all-spin union bound does not
force `Delta=o(n^(3/2))`, let alone Delta<=2 for exact minimizers. This is a
scale diagnosis of this implementation, not an assertion that every abstract
profile achieving those counts is realizable.

For the NEW selectable repaired near-minimizers Delta<=2 already, so the
nonempty high-state test t<Delta/2 is no longer available. Consequently (7)
cannot be recycled as a stronger ground-face theorem for those repaired
inputs. The added near-minimality error O(n^(5/4)) also exceeds the order-n
radial mean benefit available in (8) at t=O(n^(3/2)).

## 4. What remains outside these results

Neither (1) nor (6) asserts `Z_B^-(beta)<=Z_A^-(beta)e^{o(n)}`. The repair
raises the smaller negative endpoint, and its negative finite-temperature
entropy can change. An absolute-pressure or width-to-absolute transformation
would need additional joint landscape control. The present paid statements
preserve original cap to a power-saving error, provide near-exact polarity
balance, preserve normalized operator regularity, and can simultaneously avoid
increasing declared POSITIVE pressure profiles at leading order.
