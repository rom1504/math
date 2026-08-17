# Composition audit for the universal low-field affine interface

Date: 2026-08-17.

Status: proof draft.  The positive statement below is a deterministic
one-step near-order insertion theorem.  The two lower bounds show why this
particular affine algebra does not improve the random-bridge scale or give a
reusable cross-order recurrence.

## 1. The one bounded endpoint hidden in the canonical frame

Let `A` be a hollow symmetric signing of order `n`, choose an absolute ground
state `x` with `|H_A(x)|=Q(A)=Q`, and let

```math
W=(x,x^{\{i_1\}},\ldots,x^{\{i_k\}})
```

be the canonical frame of LA.2.  Put `s=k+1` and assume that `k>=4` is even,
so `s` is odd.  Choose an endpoint `epsilon_*` with `epsilon_0=1`, with
`k/2+1` positive single-flip coordinates and `k/2-1` negative ones.  Then

```math
t=\sum_(j=0)^k\epsilon_j=3.
```

Consequently, for `g=W epsilon_*`,

```math
\operatorname {sgn}(g)=x,
\qquad |g_i|\in\{1,3,5\},
\qquad \boxed{\|g\|_1=3n-4}.                       \tag{LC.1}
```

Indeed, outside `I` the field is `3x_i`; at the `k/2+1` positive port
indices it is `x_i`, and at the other `k/2-1` port indices it is `5x_i`.
This endpoint is better than the balanced `t=1` endpoint: its field is still
only `O(n)`, but its selector is the exact ground state rather than merely an
`O(kQ/n+k^2)` near-ground state.  In particular, for every real `a`,

```math
\boxed{\mathcal B_A(ag)=Q+|a|\|g\|_1.}             \tag{LC.2}
```

The upper bound is cap plus Holder.  For the lower bound use `x` when
`a>=0` and `-x` when `a<0`, and choose the independent quadratic sign in
the absolute trust response to orient `H_A(x)`.

## 2. A sharp-order microcanonical compiler

The rowwise microcanonical construction admits a stronger global estimate
than the rowwise-absolute concentration bound in Theorem 21.66 when
`s>>sqrt n`.  At the canonical `s=Theta(sqrt n)` scale the two bounds have
the same `n^(5/4)` order.

### Lemma LC.1 (bilinear microcanonical compilation)

Let `g in Z^n`, let `s>=2`, and suppose

```math
|g_i|<=s,
\qquad g_i=s\pmod2.
```

There are `B in {+-1}^{n times s}` and `eta_* in {+-1}^s` such that

```math
B eta_*=g
```

and, for an absolute constant `C`,

```math
\boxed{
 \max_(eta in {+-1}^s)
 \left\|B eta-{\langle eta,eta_*\rangle\over s}g\right\|_1
 \le C\sqrt{ns(n+s)}.}                           \tag{LC.3}
```

#### Proof

Switch columns so that `eta_*=1`.  Independently in every row, sample `b_i`
uniformly from the Boolean slice of sum `g_i`.  For fixed `eta`, write

```math
a={\langle eta,1\rangle\over s},
\qquad X_i=b_i^Teta-a g_i.
```

Then `E X_i=0`.  Moreover `X_i` is twice a centred sample-without-replacement
sum from a population in `[-1,1]`.  Hoeffding's comparison theorem for
sampling without replacement, followed by Hoeffding's lemma, gives

```math
E exp(lambda X_i)<=exp(C_0lambda^2s)              \tag{LC.4}
```

with an absolute `C_0`, uniformly in the row slice and in `eta`.  The rows
are independent.  Hence, for each fixed `(z,eta) in {+-1}^n times
{+-1}^s`,

```math
\sum_i z_iX_i
```

is subgaussian with variance proxy `C_0ns`.  A union bound over the
`2^(n+s)` pairs gives, with positive probability,

```math
\max_(z,eta)\left|\sum_i z_iX_i\right|
 <=C\sqrt{ns(n+s)}.                              \tag{LC.5}
```

For fixed `eta`, maximizing the left side over `z` is exactly
`sum_i|X_i|`.  Thus one sampled matrix obeys (LC.3). `square`

The scale in (LC.3) is the same as the ordinary random bipartite-signing
bound.  Unlike an iid bridge, however, the compiler realizes the selected
field `g` exactly.

## 3. The strongest one-step composition consequence

### Theorem LC.2 (anchored near-order insertion)

Let `A`, `x`, `W`, `s`, and `g` be as in Section 1, with `s<=n`.  There is a
hollow exact signing

```math
P=\begin{pmatrix}A&B\\B^T&C\end{pmatrix}          \tag{LC.6}
```

of order `n+s` such that

```math
\boxed{
 \left|Q(P)-\{Q(A)+3n-4\}\right|
 \le C\sqrt{ns(n+s)}+C_1s^{3/2}.}                \tag{LC.7}
```

Here `C,C_1` are absolute.  In particular, if `s=o(n)`, then

```math
Q(P)=Q(A)+o(n^{3/2}).                             \tag{LC.8}
```

#### Proof

The field in (LC.1) satisfies the parity and magnitude hypotheses of Lemma
LC.1, so choose `B` with target `eta_*` and error `R` bounded by (LC.3).
A standard random-sign union bound supplies a hollow order-`s` signing `C`
with `Q(C)<=C_1s^(3/2)`.  For arbitrary `eta`, put

```math
a={\langle eta,eta_*\rangle\over s},
\qquad e_eta=B eta-a g.
```

Trust-response Lipschitzness and (LC.2) give

```math
\mathcal B_A(B eta)
 <=\mathcal B_A(ag)+\|e_eta\|_1
 <=Q+\|g\|_1+R.                                  \tag{LC.9}
```

The block trust identity, or just introduction of the global absolute sign,
therefore gives

```math
Q(P)<=Q+\|g\|_1+R+Q(C).                          \tag{LC.10}
```

At the target endpoint `eta_*`, the child trust response is exactly
`Q+\|g\|_1`, while the internal shore can subtract at most `Q(C)`.  Hence

```math
Q(P)>=Q+\|g\|_1-Q(C).                             \tag{LC.11}
```

Equations (LC.1), (LC.3), and (LC.10)--(LC.11) prove (LC.7). `square`

Thus the affine frame gives a deterministic *anchored* insertion: one
declared shore word sees the exact old ground state and adds the known
linear baseline `3n-4`.  This is stronger information than an unconditioned
random bridge supplies.  As a cap upper bound, however, its error has the
same `Theta(n sqrt s)` order when `s=o(n)`.

At the canonical interface width `s=Theta(sqrt n)`, the statement is exactly

```math
\boxed{Q(P)=Q(A)+3n+O(n^{5/4}),\qquad |P|=n+Theta(\sqrt n).} \tag{LC.11a}
```

Thus both the known random bridge and the affine compiler give a subleading
single-step increment, but neither gives a derivative-scale `O(n)` main
increment together with an `o(n)` defect that could plausibly survive
`Theta(sqrt n)` successive insertions.

## 4. Two ceilings

The failure to improve the cap scale is structural, not just a weak proof of
Lemma LC.1.

### Proposition LC.3 (scalar-affine balancing has an unavoidable residual)

Let `B in {+-1}^{n times s}`, `eta_* in {+-1}^s`, and `g=B eta_*`, and
suppose `|g_i|<=G`.  Then

```math
\boxed{
 \max_eta\left\|B eta-{\langle eta,eta_*\rangle\over s}g\right\|_1
 \ge n\left(\sqrt{s/2}-{G\over\sqrt s}\right).}  \tag{LC.12}
```

In particular, for the field (LC.1), the right side is
`Omega(n sqrt s)` when `s->infinity`.

#### Proof

Average over uniform `eta`.  For every Boolean row `b_i`, the sharp
Khinchine inequality at exponent one gives

```math
E_eta|b_i^Teta|>=\sqrt{s/2}.
```

Also

```math
E_eta{|\langle eta,eta_*\rangle|\over s}
 <=s^{-1/2}.
```

Reverse triangle inequality, summed over rows, proves that the average of
the norm in (LC.12) is at least its displayed right side.  Its maximum is
therefore at least the average. `square`

This proposition is a ceiling for every argument that compresses arbitrary
shore endpoints to the one scalar `a` and pays the omitted channel by
uniform `l_1` response Lipschitzness.  It does **not** rule out a genuinely
joint theorem in which the residual cancels against the child quadratic
energy before absolute values.

There is a complementary obstruction to preserving the full affine
selector orbit exactly.

### Proposition LC.4 (all-endpoint affine physicalization is unbalanced)

Let `s` be odd, let `I subset[n]` have size `k`, and let
`B in {+-1}^{n times s}`.  Suppose that, for every endpoint `eta`,

```math
\operatorname {sgn}(B eta)
```

belongs projectively to the affine coset
`{x^S:S subseteq I}`.  Then all switched rows `x_iB_(i,.)`, `i notin I`,
are identical.  Consequently there is an endpoint `r in {+-1}^s` such that

```math
x^TBr>=(n-2k)s.                                  \tag{LC.13}
```

For every hollow exact-sign shore `C`, the resulting parent obeys

```math
\boxed{
Q\begin{pmatrix}A&B\\B^T&C\end{pmatrix}
 >=Q(A)+(n-2k)s-{s(s-1)\over2}.}                 \tag{LC.14}
```

For the canonical size `k=s-1` and `s=o(n)`, this is
`Q(A)+(1-o(1))ns`.  In particular `s=Theta(sqrt n)` creates a fixed leading
`Theta(n^(3/2))` increment.

#### Proof

Because `s` is odd, no row field has a tie.  For every two rows `i,j` outside
`I`, the hypothesis says

```math
\operatorname {sgn}((x_iB_i)^Teta)
=\operatorname {sgn}((x_jB_j)^Teta)
```

for all Boolean `eta`.  Boolean odd-majority threshold functions with
Boolean weights are injective: if `u!=v`, switch so `u=1`, put all signs on
the disagreement set equal, and choose the sum on its complement to have
minimum absolute value; then `u^Teta` and `v^Teta` have opposite signs.
Thus `x_iB_i=x_jB_j=:r`.

At endpoint `r`, every row outside `I` contributes `s` to `x^TBr`; every
row in `I` contributes at least `-s`, proving (LC.13).  Evaluate the parent
at `(x,r)` and use `H_C(r)>=-Q(C)>=-s(s-1)/2` to obtain (LC.14). `square`

Propositions LC.3 and LC.4 give an exact dichotomy for the two obvious uses
of the affine algebra:

1. preserving all affine selectors forces an `Omega(ns)` endpoint;
2. balancing those endpoints through a scalar affine compiler forces an
   `Omega(n sqrt s)` omitted channel if it is paid separately.

At `s=Theta(sqrt n)`, the first is leading and the second is the familiar
`Theta(n^(5/4))` subleading random-bridge scale.

## 5. Repeated recomputation does not close the scale gap

Suppose (LC.7) is applied repeatedly with relative increments
`theta_j=s_j/n_j ->0`.  Its normalized one-step error is

```math
O(\sqrt{theta_j}+theta_j^{3/2}).                  \tag{LC.15}
```

Traversing a fixed positive logarithmic change of order requires
`sum_j theta_j=Theta(1)`.  If `max_j theta_j<=delta`, then

```math
\sum_j\sqrt{theta_j}
 >= {\sum_jtheta_j\over\sqrt delta}
 =Omega(delta^{-1/2}).                            \tag{LC.16}
```

Thus additive reuse of the one-step certificate has no summable geometric
defect: vanishing relative steps make the accumulated certified loss worse,
not better.  Recomputing the low-field frame does not create a congruence
between successive interfaces, so the anchored target word also carries no
depth-controlled state through the iteration.

This is a ceiling for recurrence by triangle accumulation.  It does not
exclude a future global construction that correlates all levels and cancels
their errors jointly.

## 6. Comparison with the existing frontier (after the derivation)

- Theorem 21.50 identifies the exact histogram carrier and its fixed-scale
  coreset.  The canonical low-field frame has a tiny histogram presentation,
  but Proposition LC.4 shows that exact dynamic physicalization of all of
  its selectors is macroscopically unbalanced.
- Theorems 21.52--21.55 explain when product-closed Boolean ports can be
  reused.  LA.1 supplies exact *near-top* odd-product closure but no merge
  law between frames recomputed from different children.  Equation (LC.16)
  quantifies the resulting missing dynamic congruence.
- Corollary 21.65 is the PC.3-specific super-target endpoint obstruction.
  Proposition LC.4 is its universal affine-coset analogue: it uses only the
  fact that all bulk selector signs must agree.
- Lemma LC.1 sharpens the rowwise concentration scale in Theorem 21.66 to
  the natural bilinear scale when `s>>sqrt n` (and matches its order at
  `s=Theta(sqrt n)`).  Theorem LC.2 then gives an exact selected-field
  anchor, but Proposition LC.3 proves that scalar affine compilation cannot
  improve the order of this error.
- Theorem 21.67 uses special coherent information to preserve a *difference*
  between two parents.  Nothing in LA.1 relates its generic residual fields
  to the child energy, so that joint mechanism is unavailable here.
- The general random insertion inequality (ledger (10.1684)) already gives

  ```math
  Q(A_{n+s})<=Q(A_n)+O(\sqrt{ns(n+s)}).
  ```

  Therefore LC.2 strengthens the response information at the selected
  endpoint but does not strengthen known near-order cap transfer.

## 7. Verdict

The universal low-field affine shell algebra has one honest composition
consequence: an exact-sign, deterministic, selected-field-anchored
near-order insertion theorem.  Its cap error is nevertheless exactly at the
known random-bridge order.  The interface cannot be made into a better
cross-order carrier by either direct exact physicalization, scalar
microcanonical balancing, or repeated recomputation.

To move beyond this ceiling one would need a theorem that pays the
`Theta(n sqrt s)` balanced residual **jointly** against the child quadratic
channel, or a cross-level construction whose errors cancel without being
bounded one step at a time.  The low-field affine algebra alone supplies
neither ingredient.
