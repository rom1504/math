# Wave 19 memo: the restriction game and its adaptivity gap

## Status

The game/LP identities, aggregate certificate, selector variance formula, and
the order-nine rational certificate below are **proved and independently
machine checked** by `tmp/restriction_dual_r19_verify.py`.  They do not prove a
power-saving restriction theorem.  They isolate its missing content as a
pure-versus-hidden-mixed adaptivity gap.  The order-nine calculation is a
finite obstruction only, not an asymptotic counterexample.

Throughout, use the ordered matrix pairing

```math
\langle B,C\rangle=\operatorname{tr}(B^{\mathsf T}C).
```

Thus `x^T A x` counts every unordered edge twice.  Let
`\mathcal D_n={\sigma xx^T:\sigma\in{\pm1},x\in{\pm1}^n/{\pm1}}`, with
diagonals ignored, and put

```math
c_A(S,d)=\langle A[S],d[S]\rangle.
```

Then `max_{d\in\mathcal D_n}c_A(S,d)=Q(A[S])`.

## 1. Exact pure and hidden-mixed games

Let `\mathcal S_m=\binom{[n]}m`.  Define

```math
V_{\rm ad}(A,m)=\min_{S\in\mathcal S_m}\max_{d\in\mathcal D_n}c_A(S,d)
=\min_{|S|=m}Q(A[S])
```

and

```math
V_{\rm hid}(A,m)
=\min_{\pi\in\Delta(\mathcal S_m)}
  \max_{d\in\mathcal D_n}\mathbb E_{S\sim\pi}c_A(S,d).
```

This is a finite zero-sum game, so finite-dimensional LP duality (not an
unjustified swap) gives the exact dual

```math
\boxed{
V_{\rm hid}(A,m)
=\max_{\lambda\in\Delta(\mathcal D_n)}
  \min_{S\in\mathcal S_m}\mathbb E_{d\sim\lambda}c_A(S,d).}
```

The primal LP minimizes `t` subject to

```math
\pi_S\ge0,\qquad \sum_S\pi_S=1,\qquad
\sum_S\pi_Sc_A(S,d)\le t\quad(d\in\mathcal D_n),
```

and the displayed dual chooses `\lambda`, with every row expectation at
least its objective value.

For uniform `S`, every ordered off-diagonal pair survives with probability

```math
p_2=\frac{(m)_2}{(n)_2}.
```

Consequently, for every fixed oriented cut `d`,

```math
\mathbb E_Sc_A(S,d)=p_2\langle A,d\rangle,
```

and hence, when `Q(A)=q_n`,

```math
\boxed{V_{\rm hid}(A,m)\le p_2q_n.}
```

The inequality may be strict because the optimal hidden mixture need not be
uniform.  If

```math
I(A,m)=V_{\rm ad}(A,m)-V_{\rm hid}(A,m),
```

then the following genuinely sufficient (but stronger than necessary)
rounding estimate would prove the steering lemma at this pair:

```math
I(A,m)\le
\left[\left(\frac mn\right)^{3/2}-p_2\right]q_n
+O(n^{3/2-c}).
```

Indeed `V_ad=V_hid+I` and `V_hid\le p_2q_n`.  This preserves the correct
quantifiers: for each requested `m`, one may choose a different exact
order-`n` minimizer.  It does not require one nested chain.

## 2. Aggregate rank-one certificate

For each `S`, choose an oriented local ground `d_S` and embed `d_S[S]` by
zero outside `S`.  For any distribution `\pi` on `m`-sets put

```math
H_\pi=\mathbb E_{S\sim\pi}\bigl[1_S1_S^{\mathsf T}\odot d_S\bigr].
```

Then exactly

```math
\boxed{\langle A,H_\pi\rangle
=\mathbb E_{S\sim\pi}Q(A[S]).}
```

Each zero-extended local rank-one direction lies in
`\operatorname{conv}\mathcal D_n`: keep its signs on `S` and choose the
outside signs independently and uniformly.  Thus

```math
H_\pi\in\operatorname{conv}\mathcal D_n,
\qquad \langle A,H_\pi\rangle\le Q(A).
```

For uniform `\pi`, also `|(H_\pi)_{ij}|\le p_2`.  Therefore the assertion
that every `m`-subset has excess above `L` has the exact aggregate
certificate

```math
H\in\operatorname{conv}\mathcal D_n,\quad
|H_{ij}|\le p_2,\quad \langle A,H\rangle>L.
```

These convex facts alone give only the parent cap `L<Q(A)`, not the desired
`(m/n)^{3/2}` cap.  The missing fact is correlation: the maximizing column
is allowed to depend on the revealed row `S`.

There is an exact way to see this.  Randomly extend each `d_S` outside `S`
as above and call the resulting full oriented cut `D`.  If `S'` is an
independent uniform `m`-set, then

```math
\mathbb E c_A(S,D)=\mathbb E_SQ(A[S]),\qquad
\mathbb E c_A(S',D)=p_2\mathbb E_SQ(A[S]).
```

Thus this canonical coupling's revealed-selector gain is exactly
`(1-p_2)\mathbb E_SQ(A[S])`.  Merely averaging the selected subgradients
cannot remove the adaptive correlation.

## 3. Exact fixed-state selector variance

Fix `x\in{\pm1}^n`, write `w_{ij}=a_{ij}x_ix_j` for unordered edges, and
let

```math
Y_S=\sum_{i<j}w_{ij}1_{\{i,j\subset S\}},\quad
W=\sum_{i<j}w_{ij},\quad
r_i=\sum_{j\ne i}w_{ij},\quad N=\binom n2,
```

with `p_k=(m)_k/(n)_k`.  Partitioning pairs of edges into equal, adjacent,
and disjoint pairs gives

```math
\boxed{
\operatorname{Var}Y_S
=(p_3-p_4)\sum_i r_i^2
+N(p_2-2p_3+p_4)
+(p_4-p_2^2)W^2.}
```

Here `\sum_i r_i^2=x^TA^2x`, `2W=x^TAx`, and the matrix payoff
`X_S=x_S^TA[S]x_S=2Y_S` has four times this variance.  The verifier checks
the identity by exact enumeration.

Global minimality supplies only a subcritical fixed-state estimate.  By
polarization, `||A||_{\infty\to1}\le2Q(A)`.  Therefore

```math
x^TA^2x=\|Ax\|_2^2
\le\|Ax\|_\infty\|Ax\|_1
\le2(n-1)Q(A).
```

For a competitive minimizer `Q(A)=O(n^{3/2})`, the variance bound is at
most `O(n^{5/2})`, so a fixed state's standard deviation can be
`O(n^{5/4})`.  Variance alone supplies no subgaussian tail.  Even after a
hypothetical tail upgrade, a naive union bound must pay for `2^n` oriented
projective states and returns a leading-scale (or worse) selector error.
This is a scoped no-go for second-moment-plus-raw-union-bound arguments,
not a theorem excluding correlated chaining or minimizer-specific ground
compression.

A rigorous conditional (but currently far too strong) criterion is obtained
directly from one-sided Chebyshev.  For `L>p_2Q(A)`, let `\mu_d=p_2\langle
A,d\rangle` and let `v_d` be the displayed variance (times four).  If

```math
\sum_{d\in\mathcal D_n}
\frac{v_d}{v_d+(L-\mu_d)^2}<1,
```

then some `m`-set has `Q(A[S])\le L`.  A useful new theorem would have to
replace this raw exponential state sum by a minimizer-specific effective
family or correlated selector-process bound.  This is a concrete structural
parameter, not an endpoint telescope, but present estimates do not control
it.

## 4. Exact finite audit: `A_9`, `m=8`

For the exact minimizer `A_9` in (10.550), all nine rows have maximum 24, so

```math
V_{\rm ad}(A_9,8)=24.
```

Exact rational primal and dual certificates give

```math
\boxed{V_{\rm hid}(A_9,8)=\frac{464}{25}=18.56.}
```

The primal distribution over the deleted vertex is

```math
\pi=\frac1{25}(4,2,4,0,4,4,2,5,0).
```

The verifier exhausts all 512 oriented projective columns and checks that
their payoff is at most `464/25`.  A dual distribution on eight displayed
columns in the verifier has numerator weights

```math
(388,64,149,173,160,340,51,75)/1400.
```

Its nine row means are `464/25`, except the fourth (zero-based deleted
vertex 3), which is `488/25`.  Hence weak duality proves equality without
floating-point optimization.  In particular,

```math
I(A_9,8)=24-\frac{464}{25}=\frac{136}{25},
```

while uniform mixing has value `p_2q_9=(7/9)24=56/3` and uniform-based
rounding loss `24-56/3=16/3`.  The target is
`24(8/9)^{3/2}=20.113\ldots`.  Thus exact global minimality does not force a
zero-error or finite sharp rounding of the hidden selector.  This example
does not amplify and does not falsify a power-saving asymptotic lemma.

## 5. Updated route-level conclusion

The dual does not turn global minimality into a power saving by itself.
Uniform mixing already has the desired scalar coefficient; all difficulty is
rounding a hidden selector when the ground state is chosen after `S` is
revealed.  The aggregate rank-one certificate remains inside the full cut
polytope and therefore sees only `Q(A)`.  The exact variance identifies the
extra data as the full family of ground-state fields `x^TA^2x` together with
the correlated metric/entropy of those states; the scalar cap and its
standard operator-norm consequence do not control that family sharply
enough.

The most precise next positive lemma is therefore a minimizer-specific
rounding theorem for this finite game: bound its revealed-selector advantage
by the coefficient gap `[(m/n)^{3/2}-p_2]q_n` with a power-saving remainder,
possibly after choosing a target-specific exact minimizer.  A legitimate
negative result would be an asymptotic exact-minimizer family whose game gap
exceeds that budget by `\Omega(n^{3/2})`.  `A_9` is only the required finite
wall that any proposed proof must survive.
