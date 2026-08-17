# Second independent audit: bounded-cap near-top tails

**Verdict:** **PASS.**  The requested archive-classification repair has
now been incorporated into the frozen source.  The main
theorem BT.3 is correct: every hollow complete signing with
`Q(A)<=C n^(3/2)` has a fixed exponential deficit in each one-sided
near-endpoint tail, uniformly over the whole bounded-cap class.  Corollary
BT.4 therefore proves the exact-minimizer hypothesis `L_tail`.  BT.5 is also
mathematically correct, but its qualitative opposite-polarity conclusion is
already Theorem OC.2 in
`exact_minimizer_oriented_core_separation.md`; only its explicit retained-cap
bounds BT.22--BT.23 sharpen that archived statement.

## 1. Frozen source

```text
extremal_information/drafts/bounded_cap_near_top_tail.md
sha256 6d01de4d8dbe07bf87758560a452b95d9b65e62bfce1a28bb37cf6eac4c4c6ef
```

This audit uses the half-energy normalization

```math
H_A(x)={1\over2}x^TAx,\qquad
P(A)=\max H_A,\qquad N(A)=-\min H_A,\qquad Q(A)=\max(P,N).
```

## 2. Deterministic inputs

For a principal split `T sqcup R`, positive block maximizers `u,v` give
full spins `(u,v)` and `(-u,v)` with identical internal energy and opposite
cross energy.  Hence

```math
P(A)\ge P(A[T])+P(A[R]).
```

Applying this to `-A` proves the corresponding `N` inequality.  The same
argument with a fixed principal spin and random complementary spins also
checks `Q(A[R])<=Q(A)`.  Thus BT.1 and every use of principal monotonicity
are exact and have no hidden factor two.

## 3. Bollobas--Scott normalization

Let `B=A[R]` have order `m` and half-energy caps `p=P(B)`, `q=N(B)`.
The archived theorem is in doubled energy:

```math
P_2=2p,\qquad N_2=2q,\qquad R_2=2(p+q),
```

```math
P_2R_2\ge{(1-r_B^2)m^3\over1600},\qquad
r_B={P_2\over m(m-1)}.
```

If `M=Q(A)`, principal monotonicity gives `p+q<=2M`; consequently

```math
4p(p+q)=P_2R_2\le8pM,
```

and therefore

```math
\boxed{p\ge{(1-r_B^2)m^3\over12800M}.}
```

This confirms BT.5's denominator `12800`.  If `m>=3n/4` and
`M<=C_0 n^(3/2)`, then `r_B=O_C(n^(-1/2))`; hence the graph-density
hypothesis `p_G(1-p_G)>=1/m` holds and `1-r_B^2>=1/2` eventually.  Thus

```math
p\ge {27\over1638400C_0}n^{3/2}
   > {1\over100000C_0}n^{3/2}.
```

The generic first sentence of BT.2 (arbitrary fixed `epsilon`) follows by
replacing `(3/4)^3` by `(1-epsilon)^3`; the displayed proof specializes to
the only value later used, `epsilon=1/4`.

## 4. Conditional tail calculation

At Pietsch parameter `epsilon=1/4`, PC.1 supplies `|T|<n/4` and, uniformly
for every frozen `x_T`,

```math
\|A[R]\|_{2\to2}\le32K_GC_0\sqrt n,
\qquad
\|A_{R,T}x_T\|_2\le8\sqrt2K_GC_0n.
```

BT.2 and BT.1 imply

```math
P(A[T])\le P(A)-\gamma_Cn^{3/2},\qquad
\gamma_C={1\over100000C_0}.
```

Therefore, with `d_C=gamma_C/2`, the event

```math
P(A)-H_A(x_T,X_R)<d_Cn^{3/2}
```

forces

```math
h^TX_R+H_{A[R]}(X_R)>{\gamma_C\over2}n^{3/2}.
```

This direct inclusion is the key quantifier check: although the endpoint
`P(A)` moves with `A`, the conditional excess is the same fixed gap for
every signing and every `x_T`.

Splitting at `gamma_C n^(3/2)/4`, the Rademacher linear tail gives

```math
\exp\left[-{\gamma_C^2\over4096K_G^2C_0^2}n\right].
```

For the quadratic term, `H_{A[R]}=X_R^TA[R]X_R/2`, hollowness centers the
form, and `\|A[R]\|_F^2<=n^2`.  Hanson--Wright therefore gives

```math
2\exp\left[-c_{HW}n\min\left{
 {\gamma_C^2\over4},{\gamma_C\over64K_GC_0}
\right}\right].
```

These are exactly BT.14--BT.15.  They are uniform in the conditioned spin,
so averaging over `x_T`, multiplying by `2^n`, and absorbing the fixed
prefactors proves BT.8.

Applying the same argument to `-A` proves the negative endpoint.  Finally,
if `Q(A)-|H_A(x)|<d_Cn^(3/2)`, then according to the sign of `H_A(x)` one
has either

```math
P(A)-H_A(x)\le Q(A)-|H_A(x)|
```

or the analogous inequality with `N(A)+H_A(x)`.  Thus the absolute event
is contained in the union of the two one-sided events.  BT.8a has the
correct event direction and loses only a harmless factor two.

## 5. Exact-minimizer implication

The random-sign bound

```math
M_n\le\sqrt{(\log2)(n^3-n)}<n^{3/2}
```

puts every exact minimizer in BT.3 with `C=1`.  After globally orienting it
so that `P(A)=Q(A)`, BT.8 is exactly Theorem 21.8's hypothesis (21.38),
with `d_0=1/200000` and a fixed positive entropy deficit.  There is no use
of exact minimality, convergence, optimizer heredity, or `L_tail` inside
the proof; hence the implication is not circular.

## 6. BT.5

For `m=n-k`, `k=o(n)`, the archived random-bridge inequality gives

```math
0\le\Delta_n=M_n-M_m
\le M_k+\sqrt{2(\log2)mkn}=o(n^{3/2}).
```

If an oriented exact minimizer has
`P(A[T])>=(t-o(1))n^(3/2)`, then

```math
P(A[R])\le M_n-P(A[T])<M_n-\Delta_n\le Q(A[R])
```

eventually.  Hence `Q(A[R])=N(A[R])`, and negative block
superadditivity gives

```math
N(A[T])\le M_n-N(A[R])\le\Delta_n.
```

This verifies BT.18--BT.21.  Reapplying the checked one-sided product bound
to the linear complement yields

```math
P(A[R])\ge\left({1\over12800c_n}-o(1)\right)n^{3/2},
\qquad c_n={M_n\over n^{3/2}},
```

and hence BT.23.  The known `limsup c_n<=1/2` gives the stated
`1/6400-o(1)` lower coefficient.  The edge-count bound
`P(A[T])<=|T|(|T|-1)/2` gives BT.27 with factor `sqrt(2t)`.

## 7. Archive collision and final classification

The following are archived inputs rather than new results:

- BT.1: ledger (10.13);
- BT.2's product estimate: ledger (10.148)--(10.151) and
  `artifacts/one_sided_energy_product.md`;
- BT.5's near-order estimate: Theorem 36.15;
- BT.5's qualitative opposite-polarity core/complement conclusion:
  Theorem OC.2 in `exact_minimizer_oriented_core_separation.md`.

The repaired source now records the last collision explicitly.  BT.22--BT.23
are a useful explicit strengthening, while BT.5 is no longer classified as
an entirely new structural theorem.

No archived statement found by searches for the joint terms `Pietsch`,
`one-sided product`, `linear complement`, `near-top tail`, and `principal
core` contains BT.3's feedback step.  That step is the substantive new
result:

```text
Pietsch linear complement
  + Bollobas--Scott retained one-sided cap
  + oriented block superadditivity
  + conditional Hanson--Wright
  => uniform fixed-rate endpoint tails for every bounded-cap signing.
```

Accordingly:

```text
BT.3: PASS, genuinely new Level-5 bounded-cap theorem.
BT.4 / L_tail: PASS, unconditionally proved.
BT.5: PASS mathematically; mostly an archived theorem plus explicit bounds.
Convergence or a cross-order recurrence: NOT implied.
```

## 8. Frozen post-repair check

The only substantive change from the originally audited source is the
archive-classification repair in item 2 of Section 4: it now identifies
OC.2 as the prior qualitative opposite-polarity theorem and restricts the
new BT.5 increment to BT.22--BT.23.  No definition, hypothesis, estimate,
constant, proof step, or downstream implication changed.  I rechecked the
repaired paragraph against the archive and retain the unconditional
**PASS** verdict for source SHA-256

```text
6d01de4d8dbe07bf87758560a452b95d9b65e62bfce1a28bb37cf6eac4c4c6ef
```
