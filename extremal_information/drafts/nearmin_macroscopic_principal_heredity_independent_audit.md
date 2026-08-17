# Independent audit: macroscopic principal heredity

Date: 2026-08-17.

Audited file:
`extremal_information/drafts/nearmin_macroscopic_principal_heredity.md`.

## Verdict

**PASS WITH REQUIRED SCOPE REPAIRS.**  The mathematical implications
(MH.1)--(MH.4), including their constants and their uniformity over every
principal set, are correct.  No counterexample was found.  The dense-liminf
corollary is also correct.

The draft must, however, repair its archive/novelty accounting.  Near-order
random padding and the ratio-dense-subsequence corollary were already proved
in the ledger (Sections 1.8 and 10.18) and restated in
`artifacts/minimal_all_order_action_recovery.md`, (AR.12).  The present
argument slightly sharpens finite constants and isolates a useful
near-minimizer principal-restriction corollary, but it is not a new
cross-order mechanism and does not change the rigorous frontier.

One proof cross-reference should also be corrected: Corollary MH.2 uses
(MH.1) and (MH.3), not principal-heredity estimate (MH.4).

## 1. MH.1: random rectangular bridge

Put `N=m+k` and `K=mk`.  For a fixed pair `(x,y)`,

```math
S_{x,y}=x^TRy
```

is a sum of `K` independent Rademacher variables, so

```math
P(|S_{x,y}|>t)\le 2\exp(-t^2/(2K)).                 \tag{A.1}
```

Because the event in (A.1) is invariant under independently replacing
`x` or `y` by its negative, it is enough to union-bound over

```math
2^{m-1}2^{k-1}=2^{N-2}
```

projective pairs.  With

```math
t^2=2K N\log 2,
```

the bad probability is at most

```math
2^{N-2}\,2e^{-N\log2}=\frac12.                     \tag{A.2}
```

Thus a bridge with `|x^TRy|<=t` for every Boolean pair exists.  For exact
minimizers `B,C`, the block energy has the normalization

```math
H_P(x,y)=H_B(x)+H_C(y)+x^TRy,
```

with no missing factor of two, since `Q` sums over unordered edges.  The
triangle inequality gives exactly

```math
M_N\le M_m+M_k+\sqrt{2(\log2)mkN}.
```

This checks the count and constant in (MH.1), including `m=1` or `k=1`.
The finite constant is slightly sharper than the archived
`sqrt(2mk(N+2)log2)` baseline because the present proof uses both projective
symmetries and only the bridge fluctuation.  That is a finite constant
cleanup, not an asymptotic mechanism.

## 2. MH.2 and the every-deletion quantifier

For any fixed spin `x_U` on a principal set `U`, extend it by independent
unbiased spins `X` on `U^c`.  Every cross monomial has one mean-zero factor,
and every monomial internal to `U^c` has expectation zero.  Therefore

```math
H_{A[U]}(x_U)=E_X H_A(x_U,X),
```

and hence

```math
Q(A[U])\le Q(A).                                    \tag{A.3}
```

If `|U|=m`, `k=N-m`, and `Q(A)<=M_N+eta`, then (A.3) and (MH.1) give

```math
0\le Q(A[U])-M_m
\le \eta+M_k+\sqrt{2(\log2)mkN}.                   \tag{A.4}
```

The first inequality uses that `A[U]` is itself a signing.  Crucially, the
right side of (A.4) is independent of `U`; no union bound over subsets is
needed.  Thus the theorem really does hold simultaneously for **every**
deletion set of the declared size.

Using (MH.3), `m<=N`, and `k<=N`, one gets the explicit normalized estimate

```math
{Q(A[U])-M_{N-k}\over N^{3/2}}
\le {\eta\over N^{3/2}}
 +\sqrt{\log2}\left({k\over N}\right)^{3/2}
 +\sqrt{2\log2\,{m\over N}{k\over N}}.            \tag{A.5}
```

This is (MH.4).  Consequently, for any sequence of signings with
`eta_N=o(N^(3/2))`, any sequence `k_N=o(N)`, and **any choices** of sets
`U_N` of size `N-k_N`,

```math
Q(A_N[U_N])-M_{N-k_N}=o((N-k_N)^{3/2}).             \tag{A.6}
```

The draft's every-`o(n)`-deletion claim is therefore valid.  To prevent
misreading, “vanishing near-minimizer” should be defined explicitly as
`eta_N/n^(3/2)->0`, preferably by displaying the quantified sequence form
(A.6).  No ground-state, shell, response, or gauge heredity follows from
(A.3); the draft correctly disclaims all of those stronger conclusions.

## 3. MH.3: standalone random-sign cap

For order `k>=2`, a fixed projective spin has energy equal in law to a sum
of

```math
E_k=\binom{k}{2}
```

independent signs.  Union-bounding the two-tail estimate over `2^(k-1)`
projective spins at

```math
t_k^2=2E_k(k+1)\log2=(\log2)(k^3-k)
```

gives

```math
2^{k-1}\,2e^{-(k+1)\log2}=\frac12.                 \tag{A.7}
```

Thus some signing has cap at most `t_k`.  For `k=1`, `M_1=0` and the formula
also has right side zero, although the probabilistic proof divides by
`E_1=0` and must remain separated as in the draft.  There is no small-order
exception: for example the bound is `sqrt(6 log2)>1=M_2`.  Direct comparison
with the recorded exact values through order 14 also satisfies (MH.3).

The archived generic random-sign estimate used a slightly looser `k+2`
factor.  The present `k+1` is justified by the projective count in (A.7).

## 4. Dense-liminf corollary

Assume explicitly that `n_j` is strictly increasing, tends to infinity,
and

```math
n_{j+1}-n_j=o(n_j).
```

For `n_j<n<n_{j+1}`, put `k=n-n_j`.  Uniformly in that interval,
`k=o(n_j)` and `n/n_j->1`.  Equations (MH.1) and (MH.3) give

```math
M_n\le M_{n_j}+O(k^{3/2})+O(n_j\sqrt{k})
     =M_{n_j}+o(n_j^{3/2}).                         \tag{A.8}
```

If `M_{n_j}/n_j^(3/2)->L=liminf M_n/n^(3/2)`, then (A.8) yields
`limsup M_n/n^(3/2)<=L`; the reverse inequality is the definition of `L`.
The endpoint `n=n_j` needs no padding (equivalently, take the trivial
`k=0` case outside the stated `m,k>=1` theorem).

This proof is correct.  The draft should replace “Equations
(MH.3)--(MH.4) make the added term ...” by “Equations (MH.1) and (MH.3) ...”.
It should also say that the criterion is an archive rediscovery:

- ledger Section 1.8 already states the random-padding inequality and the
  ratio-dense subsequence convergence criterion;
- ledger Section 10.18 states the stronger cluster-set consequence;
- `minimal_all_order_action_recovery.md`, Section 7, says near-order transfer
  is already proved in both directions and restates the ratio-dense version.

## 5. Counterexample and overclaim search

The following possible failure modes were checked and do not invalidate the
mathematics.

1. **Independent projective signs.**  Although the full block energy is not
   invariant when only one block spin is negated, the event being union-bounded
   is `|x^TRy|<=t`, which is invariant.  The count `2^(N-2)` is valid.
2. **Factor-two normalization.**  The cross term is `x^TRy`, not
   `2x^TRy`, under the unordered-edge definition of `H`.
3. **Restriction expectation.**  Internal complement edges also average to
   zero because their two distinct spins are independent; no constant term
   remains.
4. **Arbitrary versus selected deletion.**  Estimate (A.4) contains no
   set-dependent term, so it proves every deletion, not merely existence or
   a random-subset assertion.
5. **Finite nonheredity.**  The archived order-11 exact minimizer all of whose
   one-vertex deletions have cap 17 while `M_10=13` is compatible with (A.4):
   the theorem is normalized asymptotic heredity, not exact optimality or a
   bounded absolute gap.
6. **PP.4 use.**  In PP.4 the parent is an exact minimizer and the exceptional
   shore has `|S|=o(n)`.  Applying (A.6) to `A[S^c]` is legitimate.  It says
   nothing about transporting the two poles, and the draft correctly notes
   that deleting `S` destroys their interface response.

No scalable counterexample to (MH.1)--(MH.4) exists because their proofs are
elementary and uniform.  The only substantive overclaim found is historical:
the dense-order transfer is not new, and macroscopic cap heredity is a simple
newly isolated corollary of two archived modules rather than a new structural
or compositional state.

## 6. Required repairs before promotion

1. Mark Corollary MH.2 as archive rediscovery and cite ledger Sections 1.8,
   10.18 and `minimal_all_order_action_recovery.md`, Section 7.
2. Reclassify the “new increment” in archive comparison as an elementary
   synthesis: the only new item is the explicit all-`o(n)` near-minimizer
   restriction corollary (plus harmless finite constant sharpening).  State
   explicitly that it causes no frontier change.
3. Correct the Corollary MH.2 proof reference from `(MH.3)--(MH.4)` to
   `(MH.1) and (MH.3)`.
4. Define “vanishing near-minimizer” with the sequence quantifiers used in
   (A.6), and say `n_j` is strictly increasing in Corollary MH.2.

After these repairs, the draft is mathematically and evidentially sound.
