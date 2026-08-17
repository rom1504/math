# Independent audit: bounded-cap near-top tails

**Verdict:** **PASS.**  The fixed-rate two-sided thin-tail theorem is
correct, its constants are conservative, and its strongest quantifier
really is uniform over every complete signing with
`Q(A)<=C n^(3/2)`.  Exact minimality is unnecessary.  The proof has no
hidden use of convergence, a tail hypothesis, or optimizer heredity.

## 1. Frozen source

```text
extremal_information/drafts/bounded_cap_near_top_tail.md
sha256 6d01de4d8dbe07bf87758560a452b95d9b65e62bfce1a28bb37cf6eac4c4c6ef
```

All calculations below use the half-quadratic normalization

```math
H_A(x)={1\over2}x^TAx,
\quad P(A)=\max H_A,
\quad N(A)=-\min H_A,
\quad Q(A)=\max(P,N).
```

## 2. BT.1 and principal monotonicity

For positive block maximizers `u,v`, the two spins `(u,v)` and `(-u,v)`
retain the same two internal energies and reverse the complete cross term.
One cross orientation is nonnegative, proving

```math
P(A)>=P(A[T])+P(A[R]).
```

Applying this to `-A` proves the `N` statement.  There is no factor two in
the unordered-edge normalization.  This is the archived one-sided block
superadditivity law.

The other deterministic input used later is also exact: for any fixed spin
on a principal set, averaging independent unbiased spins on its complement
recovers the principal energy.  Hence

```math
Q(A[R])<=Q(A).
```

## 3. BT.2: one-sided product normalization and constant

Let `B=A[R]` have order `m`, with half-energy one-sided caps `p,q`.  In the
archived doubled normalization,

```math
P_2=2p,\qquad N_2=2q,\qquad R_2=2(p+q).
```

The Bollobas--Scott translation gives

```math
P_2R_2>={ (1-r_B^2)m^3\over1600},
\qquad r_B={P_2\over m(m-1)}.
```

If `M=Q(A)`, then `Q(B)<=M`, so

```math
P_2R_2=4p(p+q)<=8pM.
```

Therefore

```math
p>={ (1-r_B^2)m^3\over12800M}.
```

This verifies the source's factor `12800`.  Under `M<=C n^(3/2)` and
`m>=3n/4`, one has `r_B=O_C(n^(-1/2))` uniformly.  Thus the graph-density
hypothesis is valid and `1-r_B^2>=1/2` for all sufficiently large `n`.
It follows that

```math
p>={27\over1638400C_0}n^(3/2)
  >{1\over100000C_0}n^(3/2),
\qquad C_0=\max(C,1).
```

The numerical comparison is correct (`27>16.384`).  Because the argument
also applies to `-B`, both one-sided caps of every such linear principal
block have a uniform positive lower bound, although BT.3 only needs the
orientation matching the tail under consideration.

## 4. BT.3: conditional concentration and exact constants

At `epsilon=1/4`, the repaired PC.1 factorization gives

```math
|T|<n/4,
\quad ||A[R]||_(2 to2)<=32K_GC_0\sqrt n,
\quad ||A_(R,T)x_T||_2<=8\sqrt2K_GC_0n.
```

BT.2 and BT.1 imply

```math
P(A[T])<=P(A)-gamma_C n^(3/2),
\qquad gamma_C={1\over100000C_0}.
```

Thus the event

```math
P(A)-H_A(x)<d_Cn^(3/2),
\qquad d_C=gamma_C/2,
```

forces, conditionally on every `x_T`,

```math
h^TX_R+H_(A[R])(X_R)>{gamma_C\over2}n^(3/2).
```

This direct inclusion is important: the near-top threshold may depend on
`A`, but the probabilistic argument depends only on the fixed gap
`gamma_C/2`.  No unproved uniformity in PC.1's displayed threshold is
being invoked.

For the linear half-event, the standard Rademacher bound
`Pr(h^TX>s)<=exp[-s^2/(2||h||_2^2)]`, with
`s=(gamma_C/4)n^(3/2)`, gives exactly

```math
\exp\left[-{gamma_C^2\over4096K_G^2C_0^2}n\right].
```

For the quadratic half-event, Hanson--Wright is applied to
`X^TA[R]X` at threshold `(gamma_C/2)n^(3/2)`.  Since
`||A[R]||_F^2<=n^2`, its two branches are

```math
{gamma_C^2\over4}n,
\qquad
{gamma_C\over64K_GC_0}n.
```

These are BT.15.  Both bounds are uniform in `x_T`; averaging and absorbing
the prefactor three yields the fixed `kappa_C>0` in BT.8.

Applying the same proof to `-A` gives the lower-end thin tail.  Finally,
for any spin with `Q(A)-|H_A(x)|<d_Cn^(3/2)`, either `H_A(x)>=0`, in which
case

```math
P(A)-H_A(x)<=Q(A)-|H_A(x)|,
```

or `H_A(x)<0`, in which case the analogous inequality holds with `N`.
Thus the absolute event is contained in the union of the two one-sided
events, and another constant-factor absorption proves BT.8a.

The theorem does not need to orient `A` or assume `P(A)=Q(A)`.  It supplies
a fixed-rate tail at **both** one-sided endpoints.  Orientation is needed
only later, when the matched roof used for scalarization must coincide with
the absolute cap.

## 5. BT.4 and absence of circularity

The archived random-sign upper bound gives

```math
M_n<=\sqrt{(\log2)(n^3-n)}<n^(3/2),
```

so exact minimizers lie in BT.3 with `C=1`.  The advertised

```math
d_0=1/200000
```

is exactly `gamma_1/2`.  After orienting an exact minimizer so that its
larger endpoint is positive, BT.8 is precisely `L_tail`.

The dependency graph is acyclic:

```text
archived one-sided discrepancy product
  + archived block superadditivity
  + new PC.1 conditional concentration
  => BT.3
  => L_tail.
```

The proof never assumes an upper-tail deficit, exact-minimizer structure,
near-order density, or convergence of `M_n/n^(3/2)`.  The bound
`Q(A[R])<=Q(A)` is elementary averaging, not a minimizer heredity claim.
There is therefore no hidden circularity.

## 6. BT.5: exact-minimizer orientation anatomy

For `m=n-k` with `k=o(n)`, Theorem 36.15 gives

```math
0<=Delta_n=M_n-M_m
 <=M_k+\sqrt{2(\log2)mkn}=o(n^(3/2)).
```

For an oriented exact minimizer and a core with
`P(A[T])>=(t-o(1))n^(3/2)`, positive block superadditivity makes

```math
P(A[R])<=M_n-(t-o(1))n^(3/2).
```

But `Q(A[R])>=M_m=M_n-Delta_n`; hence, for large `n`, its negative cap
must equal its absolute cap.  Negative superadditivity then gives

```math
N(A[T])<=M_n-N(A[R])<=Delta_n.
```

All statements in BT.18--BT.21 follow, including
`N(A)=M_n-o(n^(3/2))`.  Applying the already checked one-sided lower bound
to the `n-o(n)` complement gives

```math
P(A[R])>=\left({1\over12800c_n}-o(1)\right)n^(3/2),
\qquad c_n={M_n\over n^(3/2)},
```

and BT.23 follows by positive superadditivity.  Since
`limsup c_n<=1/2`, the coefficient is at least `1/6400-o(1)`.  The
trivial core edge count yields BT.27 with the correct factor `sqrt(2t)`.

## 7. Archive comparison and classification

The archive already contains:

- BT.1 as one-sided partition superadditivity;
- the doubled-energy Bollobas--Scott product theorem and its half-energy
  `12800` consequence;
- principal averaging/monotonicity;
- the random-bridge near-order estimate used only in BT.5;
- Hanson--Wright and Grothendieck--Pietsch separately.

Searches for combinations of "Pietsch", "one-sided product", "upper
tail", and "principal core" found no archived feedback theorem of BT.3's
form.  The new increment is exact and substantive: apply the one-sided
product theorem to the linear complement of the Pietsch heavy set, feed its
positive cap through one-sided superadditivity, and then use the resulting
fixed conditional gap.  This proves a previously open uniform fixed-rate
tail law.

BT.5's qualitative opposite-polarity conclusion is a parallel duplicate of
OC.2 in `exact_minimizer_oriented_core_separation.md`, developed
concurrently.  Its explicit retained-complement estimates BT.22--BT.23 are
the sharpening.  The final source records this collision, so it does not
affect BT.3's novelty classification.

The proper classification is:

```text
BT.3: PROVED Level-5 theorem for every bounded-cap complete signing.
L_tail: PROVED, no longer conditional.
BT.5: PROVED structural refinement for exact minimizers.
Scalar physical packing: follows only after the separate BR.2 audit.
Cross-order recurrence/convergence: not implied.
```

This result warrants a frontier reset for the contextual incompressibility
route, but it does not alter the numerical asymptotic interval or reopen a
direct convergence claim.

## 8. Physical contextual corollary

The downstream implication was checked against the independently audited
selector files

```text
bounded_cap_boundary_roof_selector.md
sha256 631d5ddcc79fc868c0086a8f9bb469d201980df71455461e7d0b5f3675251e87

bounded_cap_boundary_roof_selector_independent_audit.md
sha256 b5662236e945fddd2f361c56fff9b33afdc93c7681de28885e63ac4319d775f0
```

Orient an exact minimizer so `P(A)=Q(A)=M_n`.  BT.4 gives precisely the
one-sided shell-count hypothesis (21.38).  The proof of Theorem 21.8 then
produces one exact-sign bridge of operator norm `O(sqrt n)`, exponentially
many switched children (all still of absolute cap `M_n`), matched-roof
equalities at their named poles, and both directed cross deficits of order
`n^(3/2)`.  These are exactly BR.2--BR.3's hypotheses.  For each public
pole, the audited biased exact-sign fill is common to every child, and the
resulting complete order-`2n` parents have cap `O(n^(3/2))` while the named
scalar caps differ by `Omega(n^(3/2))`.

Consequently the chain

```text
BT.3/BT.4
  + Theorem 21.8
  + audited BR.2--BR.3
=> Omega(n) all-spins-free scalar physical contextual bits
```

is rigorous and unconditional for exact minimizers.  It charges no public
query-description complexity: the query bank can contain exponentially
many separately frozen exact fills of quadratic bit complexity.  It proves
information heaviness under the project's contextual model, not an
efficient compiler, a near-minimal order-`2n` construction, or convergence.

## 9. Frozen disposition

```text
BT.1: PASS (archived input).
BT.2: PASS, including 12800 and gamma_C.
BT.3: PASS for all bounded-cap complete signings, both one-sided tails.
BT.4: PASS; L_tail is proved.
BT.5: PASS.
Archive novelty: PASS as a new feedback synthesis, not new ingredients.
Physical contextual corollary: PASS with uncharged-query scope.
Overall verdict: PASS.
```

I rechecked the archive-classification repair in the final source with hash
`6d01de4d8dbe07bf87758560a452b95d9b65e62bfce1a28bb37cf6eac4c4c6ef`;
the frozen verdict remains **PASS**.
