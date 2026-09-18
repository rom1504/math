# Genuine nonconvergence diagnostic

## Verdict

No genuine nonconvergence construction or theorem was found.  More strongly,
none of the arithmetic/design candidates tested survives both obligations
needed for separated subsequences of the **actual optima**:

1. a strict upper construction on one infinite class; and
2. a larger lower bound for **every signing** on another infinite class.

Conference, Paley, Hadamard, bent/semi-bent, maximal-excess, and covering-radius
integrality phenomena naturally address a selected algebraic matrix or a
bounded lattice effect.  They do not provide obligation 2.  Fixed arithmetic
classes also fail a rigorous density test: the normalized actual optima have a
relative-order modulus, so every ratio-dense class has the same complete
cluster set as the full sequence.

The only logically surviving nonconvergence shape is an alternating sequence
of multiplicatively separated, macroscopic log-scale epochs.  Realizing it
would require a new order-specific universal lower theorem whose hypotheses
hold on one sparse scale hierarchy but not another.  No candidate for such a
theorem emerged from the ledger or the literature search.

This is an explicit account of why the direction lacks leverage, not a proof
of convergence.

## 1. A structural theorem for the actual optima

Write

\[
 a_n=\frac{M_n}{n^{3/2}}.
\]

Throughout, `H_A(x)=sum_(i<j) a_ij x_i x_j` counts each edge once.  If
`Q(A)=max_x |x^T A x|` is used instead, then `Q(A)=2 max_x|H_A(x)|` and
`min_A Q(A)=2M_n`.  Every numerical constant below uses the one-copy `M_n`
normalization.

The following statements apply to `M_n` itself, not to a chosen construction.

### 1.1 Padding and deletion inequalities

For all `n`,

\[
 M_n\le M_{n+1}\le M_n+n. \tag{1}
\]

For the first inequality, restrict any order-`n+1` signing to `n` vertices.
For a fixed old spin `x`, its two extensions have energies
`H_A(x) +/- b dot x`, whose maximum absolute value is
`|H_A(x)|+|b dot x| >= |H_A(x)|`.  The second inequality follows by adjoining
one arbitrary signed row to an order-`n` optimizer.

For arbitrary `n,h`, independently sign the `n by h` cross block between
optimal children.  Hoeffding plus a union bound over the `2^(n+h)` spin pairs
gives

\[
 M_{n+h}\le M_n+M_h+
 \sqrt{2nh(n+h+2)\log 2}. \tag{2}
\]

The same argument on a random signing of `K_k` gives

\[
 M_k\le R_k:=\sqrt{k(k-1)(k+2)\log2}. \tag{3}
\]

The exact maximum of `R_k/k^(3/2)` over integers `k>=2` is at `k=4`, so

\[
 0\le a_k\le C_0,
 \qquad C_0=\sqrt{\frac{9\log2}{8}}
 =0.883057516887\ldots . \tag{4}
\]

### 1.2 Explicit finite relative-order modulus

Let `1 <= n <= m`, put `h=m-n` and

\[
 \beta=\frac{m-n}{m}.
\]

Equations (1)--(4) imply the two one-sided bounds

\[
 \boxed{
 a_m-a_n\le
 C_0\beta^{3/2}+
 \sqrt{2\log2\,\beta(1-\beta)(1+2/m)} }
 \tag{5}
\]

and

\[
 \boxed{
 a_n-a_m\le
 C_0\bigl[1-(1-\beta)^{3/2}\bigr]. }
 \tag{6}
\]

For (5), write

\[
 a_m-a_n
 \le\frac{M_m-M_n}{m^{3/2}}
 \le\frac{R_h+\sqrt{2nh(m+2)\log2}}{m^{3/2}}.
\]

For (6), monotonicity gives

\[
 a_n-a_m
 \le M_n\left(n^{-3/2}-m^{-3/2}\right).
\]

In particular, uniformly as `beta -> 0`,

\[
 |a_m-a_n|\le \sqrt{2\log2}\,\sqrt\beta+O(\beta). \tag{7}
\]

Thus the sequence is asymptotically Holder-`1/2` in relative order (and in
`log n` locally).

One completely explicit version is useful below.  For `epsilon>0`, define

\[
 \delta(\epsilon)=\min\left\{
 \frac12,
 \left(\frac{\epsilon}{2C_0}\right)^{2/3},
 \frac{\epsilon^2}{16\log2},
 \frac{2\epsilon}{3C_0}
 \right\}. \tag{8}
\]

If `beta <= delta(epsilon)`, then (5)--(6) give

\[
 |a_m-a_n|\le\epsilon. \tag{9}
\]

The constants are deliberately conservative.  The checker verifies (5)--(9)
numerically without using any conjectural value of the limit.

### 1.3 The cluster set is a closed interval

From (5)--(6) with `m=n+1`,

\[
 |a_{n+1}-a_n|\longrightarrow0. \tag{10}
\]

Every bounded real sequence whose successive increments tend to zero has a
connected cluster set.  Indeed, if `L=liminf a_n < c < U=limsup a_n`, choose
arbitrarily late indices below `c` followed by indices above `c`.  At the first
crossing, one term lies within the largest intervening one-step increment of
`c`; that error tends to zero.  Hence `c` is a cluster point.  The endpoints
are cluster points by boundedness.  Therefore

\[
 \boxed{\operatorname{Clust}(a_n)=[L,U].} \tag{11}
\]

This is only a necessary structural condition.  It does **not** show `L=U`.

### 1.4 Necessary log-epoch consequence if nonconvergence is real

Assume hypothetically that `g=U-L>0`, and put `eta=g/8` and
`delta=delta(eta)`.  There are arbitrarily large low centers `n` with
`a_n <= L+eta` and high centers `m` with `a_m >= U-eta`.  For a center `v`,
define the explicit integer interval

\[
 I_v=\left\{k\in\mathbb N:
 \left\lceil(1-\delta)v\right\rceil\le k\le
 \left\lfloor\frac{v}{1-\delta}\right\rfloor\right\}.
\]

If `k<=v`, the relative parameter for `(k,v)` is `(v-k)/v<=delta`; if
`k>=v`, it is `(k-v)/k<=delta`.  Equation (9) therefore gives
`a_k<=L+2eta` throughout `I_v` for a low center and `a_k>=U-2eta`
throughout `I_v` for a high center.  These bands are separated by `g/2`.
On the logarithmic axis, `I_v` has asymptotic width
`-2 log(1-delta)>0`, independent of `v`.

Low and high centers can be selected alternately.  Every transition between
opposite types then obeys the fixed multiplicative separation

\[
 \frac{\max(m,n)}{\min(m,n)}>\frac1{1-\delta}>1. \tag{12}
\]

Thus genuine nonconvergence would not consist of isolated exceptional orders.
It would create alternating macroscopic intervals in order, or epochs of
positive width on the `log n` axis.  This is **only a necessary consequence
conditional on nonconvergence**, not evidence against nonconvergence and not
a falsifier of it: a hierarchy such as powers of two could still have enough
multiplicative separation.

### 1.5 Ratio-dense order classes are rigorously excluded as two phases

If an infinite set `S={s_j}` has `s_(j+1)/s_j -> 1`, then every large `n` lies
between two members of `S` whose ratios to `n` tend to one.  Equations (5)--(7)
show that the restriction `(a_s)_(s in S)` and the full sequence have exactly
the same cluster set.

Consequently, none of the following can carry a different constant from the
full sequence:

- a fixed residue or fixed `2`-adic class;
- primes in a fixed admissible progression;
- all Paley conference orders `q+1` with `q=1 mod 4` prime;
- any finite prescription of Legendre symbols, which is again a union of
  fixed progressions when compatible.

A sparse sub-subsequence with growing arithmetic conditions is not excluded
by this theorem, but it still faces the universal-lower-bound quantifier in
Section 3.

## 2. A rigorous falsifier of conference-nonexistence stability

Exact conference nonexistence at selected orders cannot produce a leading
spectral gap.  This remains true for substantially more spectral data than
the operator norm alone.

For any `n`, choose a symmetric Paley conference matrix `C` of order
`N=n+r` with `r=o(n)`, and write a principal decomposition

\[
 C=\begin{pmatrix}A&B\\B^T&D\end{pmatrix},
 \qquad |A|=n.
\]

Since `C^2=(N-1)I`,

\[
 \|A\|_{op}\le\sqrt{N-1}=(1+o(1))\sqrt n, \tag{13}
\]

and hence

\[
 \max_x|H_A(x)|\le\frac12n\sqrt{N-1}
 =\left(\frac12+o(1)\right)n^{3/2}. \tag{14}
\]

More sharply,

\[
 A^2-(n-1)I=rI-BB^T. \tag{15}
\]

Using `||B||_F=sqrt(nr)` and `||B||_op<=sqrt(N-1)`,

\[
 \frac{\|A^2-(n-1)I\|_F}{n^{3/2}}
 \le \frac rn+\frac{\sqrt{nr(N-1)}}{n^{3/2}}
 =O\!\left(\frac rn+\sqrt{\frac rn}\right)=o(1). \tag{16}
\]

It follows that the empirical spectrum of `A/sqrt(n-1)` converges to
`(delta_{-1}+delta_{+1})/2`: (16) puts all squared eigenvalues at `1` in mean
square, `tr A=0` balances the two signs, and (13) gives uniform boundedness.
Thus all fixed spectral moments also approach the exact-conference values.

Here `N>=n` with `N-n=o(n)` is available from Paley orders by the prime number
theorem in the progression `1 mod 4`.  Therefore exact conference
nonexistence cannot force a positive gap in
operator norm, fourth-moment defect, or any fixed collection of normalized
spectral moments: asymptotically conference-flat sign matrices exist at every
order.  More precisely, this falsifies any proposed lower correction that is
a positive continuous function of those normalized defects.  Exact
determinants and other discontinuous arithmetic invariants can still
distinguish the matrices, but (13)--(16) show that nonexistence alone supplies
no macroscopic spectral defect from which a Boolean lower gap follows.  A
surviving theorem would have to see Boolean-cube geometry not encoded in
these spectral statistics.  Equations (13)--(16) do not rule out such a
Boolean theorem and hence do not prove convergence.

## 3. Candidate audit with the required quantifiers

Let `P(A)=max_x |H_A(x)|`.  A separated-subsequence proof needs, for some
`c` and fixed `epsilon>0`, both

\[
 \exists A_n\ (n\in S):\quad P(A_n)\le(c+o(1))n^{3/2}, \tag{17}
\]

and

\[
 \forall A\ (\text{order }m\in T):\quad
 P(A)\ge(c+\epsilon-o(1))m^{3/2}. \tag{18}
\]

The candidates fail as follows.

### 3.1 Conference existence versus gaps

An exact conference matrix supplies only the upper bound `1/2`.  It is not a
strict low construction.  Its nonexistence supplies no lower bound at all,
and (13)--(16) show that no positive spectral defect can be forced.  Paley
orders are ratio-dense, so exact existence gaps cannot define a distinct
actual-optimum phase.

### 3.2 Paley resonance versus nonresonance

The ledger's Paley square-wave theorem proves that a selected Paley matrix has
energy `(1/2-o(1))p^(3/2)` along resonant primes.  This says that this one
upper construction is poor; it does **not** imply `M_p` is large, because a
different signing may be better.  The inequality points the wrong way for
(18).

A rigorous all-Boolean nonresonance theorem could conceivably give (17) with
`c<1/2` on a sparse prime sequence.  Even granting that currently open
theorem, no candidate supplies (18) on a second sequence.  Fixed-character
resonant classes are ratio-dense and cannot be a uniformly high phase of the
actual optima if the full sequence has a lower cluster point.  Growing
character prescriptions evade density but make the missing universal claim
harder, not easier.

### 3.3 Hadamard, bent/semi-bent, and maximal excess

Sylvester powers and some other design orders pass the multiplicative-sparsity
test.  But flat spectrum again gives only the `1/2` upper construction.
Bent/semi-bent divisibility constrains selected Cayley/Walsh constructions,
not arbitrary edge signings.  “Maximum excess” theorems explicitly construct
a large Boolean witness for a selected Hadamard or conference matrix; they
certify the opposite of a strict low construction and still say nothing about
the minimum over all signings.

The augmented cut code has bounded, not growing, divisibility.  Its nonzero
codeword-weight gcd is

\[
 \gcd\{w(c):c\in\mathcal C_n,\ c\ne0\}=
 \begin{cases}2,&n=1\pmod4,\\1,&\text{otherwise},\end{cases} \tag{19}
\]

and `M_n=binom(n,2)-2 rho(C_n)`.  Hence the universal radius/energy lattice
has spacing `O(1)`, which is `o(n^(3/2))`.  Fixed parity and residue classes
are independently eliminated by Section 1.5.

### 3.4 Blowups, products, and deletions

Ordinary clone blowup scales a macro edge by the square of the fibre size and
therefore inflates, rather than preserves, the `n^(3/2)` constant.  Hadamard
block lifts introduce a leading fibre-orthogonal spin channel; the ledger has
a deterministic lower floor for that channel.  Tensor submultiplicativity is
already false.  Principal deletion transfers information at ratio `1+o(1)`
but no verified proportional deletion theorem preserves a hypothetical
strict constant below `1/2`.

Thus no known lift turns a sparse low order into a controlled scale hierarchy.
More importantly, even a successful low lift would address (17), not the
universal high requirement (18).

## 4. Why current quantitative bounds cannot be combined into separation

The verified all-order lower bound is

\[
 a_n\ge0.336493364431\ldots-o(1),
\]

and the verified all-order upper bound is `a_n<=1/2+o(1)`.  An explicit low
construction must lie at or above the universal lower constant, so the
existing lower theorem can never sit strictly above it on a second class.
Conversely, conference/Paley/Hadamard high energy is a property of one matrix,
not a lower bound on `M_n`.  A genuine separation therefore requires a new
order-sensitive universal lower theorem before any design construction can
close the argument.

Known universal lower mechanisms (Gaussian rounding, field-plus-spin,
spectral moments, SDP/Grothendieck arguments) depend smoothly on `n` and have
no arithmetic hypothesis that could distinguish a sparse hierarchy.  Exact
conference stability cannot supply one by Section 2.  The coding literature
located in this diagnostic treats the ordinary cocycle-code covering radius
or general dual-distance bounds; neither gives an order-sensitive theorem for
the augmented cut code at the required `n^(3/2)` deficit scale.

The primary sources checked were Sole--Zaslavsky's coding formulation of the
ordinary cocycle-code radius
([paper](https://people.math.binghamton.edu/zaslav/Tpapers/cas.sidma1994.pdf)),
the conference maximum-excess constructions of Momihara--Suda
([arXiv](https://arxiv.org/abs/1611.01305)), the Hadamard maximum-excess
families of Hirasaka--Momihara--Suda
([arXiv](https://arxiv.org/abs/1712.08984)), and de Launey's survey of sparse
Hadamard existence results
([arXiv](https://arxiv.org/abs/1003.4001)).  They supply useful analogies and
selected designs, but no theorem with the universal quantifier (18).

## 5. Final assessment

The arithmetic/design nonconvergence program should be retired unless a new
proposal begins by stating an order-sensitive theorem of the form (18) for
**all** signings.  Searching for further resonant examples, exact design
orders, residue effects, or bad behavior of a selected construction cannot
meet the requested standard.

The sparse log-epoch scenario remains logically possible, so this diagnostic
does not support declaring convergence.  It does show that nonconvergence is
currently the lower-leverage direction: its missing ingredient is not a
technical bridge around an existing construction, but an entirely absent
universal high-order mechanism.

Reproducibility: run

```bash
.venv/bin/python tmp/diagnostic_nonconvergence_check.py
```

The checker verifies the random-bound constant, the explicit modulus width,
the augmented-cut-code gcd through order 100, and Paley conference/minor
identities and defect bounds on a sample of orders.
