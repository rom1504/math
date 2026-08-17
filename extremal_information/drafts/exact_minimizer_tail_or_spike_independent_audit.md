# Independent audit: exact-minimizer tail or spike

**Verdict: PASS.**

The Hanson--Wright theorem, its converse, the cap-only bootstrap, the
exact-minimizer implication and localization, and the spectral-spike
construction are correct.  The final two-clique complete-sign construction
also handles the absolute-cap orientation correctly.  The earlier halo
quantifier, clique count, sharpness witness, and eigenvalue-sign defects have
all been repaired in the frozen source.

## 1. Frozen source

```text
extremal_information/drafts/exact_minimizer_tail_or_spike.md
sha256 7a9f4369b6d1c3b5efcc3961a3dd48e2a24dfd9733e5e36fa73424dacb9f3c51
```

This audit used the edge-sum normalization

```math
H_A(x)=\frac12x^TAx,\qquad \|A\|_F^2=n(n-1).
```

## 2. TS.1 constants and quantifiers: PASS

Substituting `s=2tn^(3/2)` into Hanson--Wright gives

```math
\frac{s^2}{\|A\|_F^2}
=n\frac{4t^2n}{n-1},
\qquad
\frac{s}{\|A\|_{2\to2}}
=n\frac{2t}{\Lambda(A)}.
```

Thus (TS.4) has the correct factor two and both correct branches.  With
`t,L` fixed, the leading prefactor two is absorbed for all sufficiently
large `n` by halving the exponent, as in (TS.6).  Equations (TS.7)--(TS.8)
strictly exclude the Frobenius branch, and rearranging the operator branch
gives (TS.9) with the displayed direction and constant.  No maximizing spin
or minimizer hypothesis is hidden in this argument.

The use in TS.2 is also correct.  If the entropy deficits tend to zero and
the orders tend to infinity, (TS.9) forces the normalized operator norms to
diverge.  For a single fixed failed deficit it gives only a fixed lower
bound, exactly as the draft says.

## 3. Cap-only and localization consequences: PASS

For (TS.11a), let `v` be a unit eigenvector with
`|lambda|=||A||_(2 to2)` and `mu=||v||_infinity`.  With
`s=sign(v)` and independent Boolean `Y_i` of means `v_i/mu`,

```math
|\mathbb E s^TAY|=|\lambda|\,\|v\|_1/\mu.
```

At a coordinate attaining `mu`, the eigenvector equation gives
`|lambda|mu<=||v||_1`.  Thus the Boolean bilinear norm is at least
`||A||_(2 to2)^2`; ordinary cube polarization bounds it by `4Q(A)` in the
half-quadratic normalization.  This proves (TS.11a).  Equations
(TS.11b)--(TS.11c) then have the correct constants: inserting
`Lambda<=2sqrt(C)n^(1/4)` into the operator branch of TS.1 gives exponent
`c_HW t n^(3/4)/sqrt(C)`.

The normalization in (TS.16a) is also correct.  For independent signs with
means `v_i/mu`, hollowness gives

```math
\mathbb E H_A(X)=\lambda/(2\mu^2),
\qquad |\lambda|\le2Q(A)\mu^2.
```

Consequently `|lambda|/sqrt(n)->infinity` and `Q(A)=O(n^(3/2))` imply
`1/mu^2=o(n)`.  The frozen source correctly uses
`|lambda|=||A||_(2 to2)`, so this also covers a negative extremal
eigenvalue.

The Grothendieck--Pietsch archive comparison is accurate at the stated
order scale.  Applying common-support removal at threshold `t sqrt n` to
`Q(A)=O(n^(3/2))` deletes `O(n/t)` vertices and leaves the complementary
principal block with operator norm `O(t sqrt n)`.  This localizes all larger
modes on a common exceptional support but does not eliminate them from the
original matrix or control the cross fields; the draft preserves that
distinction.

## 4. TS.3 and the indefinite principal block: PASS

Overwriting the unordered edges of a `k`-vertex block changes every Boolean
energy by at most `2*binom(k,2)=k(k-1)`.  Hence (TS.17) is correct.  The
principal block is `J-I`, so Cauchy interlacing gives

```math
\lambda_{\max}(B_n)\ge k-1
```

even if the full matrix is indefinite.  Therefore its operator norm is at
least `k-1`.  With `k=L_n sqrt n`, the edit loss is
`O(L_n^2n)=o(n^(3/2))` precisely under `L_n=o(n^(1/4))`.

The final quantifier is now exact.  The construction supplies a sequence in
the unrestricted `o(n^(3/2))`-near-minimizer class.  For a prescribed halo
`epsilon_n n^(3/2)`, it supplies a diverging normalized spike when
`epsilon_n sqrt n -> infinity`; it explicitly does not claim arbitrarily
thin prescribed halos.  The cap-only spectral range
`Lambda=O(n^(1/4))` from (TS.11b) is also correct, and the complete-sign
construction reaches that endpoint scale.

## 5. The two-clique complete-sign construction: PASS

Let `K=binom(k,2)`.  The positive and negative planted cliques contribute
at most `K` each in absolute value.  Together with the complement cap and
the `o(n^(3/2))` cross cap, this gives

```math
Q(C_n)\le (1/2+a^2+o(1))n^{3/2}.
```

Fixing the positive-clique spins creates energy `K+Z_+`; fixing the
negative-clique spins creates `-K+Z_-`.  Orthogonality of the remaining
Boolean characters gives

```math
\mathbb E Z_\pm=0,
\qquad
\mathbb E Z_\pm^2=O(nk^2+n^2)=o(n^3).
```

Whichever global orientation realizes the absolute cap, the corresponding
diffuse family has oriented energy at least `K-o(n^(3/2))`.  Its deficit is
therefore at most

```math
(1/2+a^2-a^2/2+o(1))n^{3/2}
=(1/2+a^2/2+o(1))n^{3/2}.
```

This proves the claimed wide-window warning for every fixed `d_0>1/2`.
The independent-random-cross-block union bound has the stated
`O(n sqrt k+k sqrt n)=o(n^(3/2))` scale.  The construction makes no claim
at the narrow exact-minimizer window.

## 6. Clique count and hypercontractive sharpness: PASS

### (TS.20)

For a fixed fractional threshold, the probability among the `2^k` clique
assignments is `exp(-Theta_eta(k))`.  Hence their **number** is
`2^k exp(-Theta_eta(k))`, and the full count is

```math
2^{n-k}\,[2^k\exp(-\Theta_\eta(k))]
=2^n\exp(-\Theta_\eta(k)).
```

The frozen (TS.20) now includes this `2^k` factor, so both displayed forms
are correct.

### (TS.26) sharpness

The source now uses the correct degree-only weighted witness: a clique on
`k=floor(sqrt n)` vertices with every edge weight `sqrt n` and zero
elsewhere.  Then

```math
\|H\|_2=\Theta(n),\qquad
\max H=\Theta(n^{3/2}),
```

and a fixed-fraction top tail has probability
`exp(-Theta(sqrt n))`.  This witness is weighted, not a complete signing,
which is exactly the scope of its degree-only sharpness statement.  The
separate unit-weight `k=n^(3/4)` clique witnesses the `n^(3/4)` cap-only
Hanson--Wright wall, as claimed.

## 7. Disposition

The theorem-level conclusion survives the audit:

```text
uniform O(sqrt n) operator norm for exact minimizers => L_tail;
failure with vanishing entropy deficit => an unbounded normalized spike.
```

It remains a conditional dichotomy, not yet a strict reduction, because no
proved exact-minimality law excludes the localized/common-support spike
branch.  The archive-collision wording is appropriately conservative: the
cap-only bootstrap and both individual and common-support localization are
identified as pre-existing, not claimed as new.  No further repair is
required.
