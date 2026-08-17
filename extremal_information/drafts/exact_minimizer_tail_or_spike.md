# Exact-minimizer upper tails: the Hanson--Wright spike alternative

**Status.** Task-local theorem and obstruction audit.  The tail theorem
below is rigorous for every hollow signing.  It reduces `L_tail` to exclusion
of one explicit spectral-spike branch.  No known consequence of exact
minimization excludes that branch.  In particular, this note does **not**
prove `L_tail`.

## 1. Normalization

For a hollow symmetric signing `A`, put

```math
H_A(x)={1\over2}x^TAx=\sum_{i<j}A_(ij)x_ix_j,
\qquad P(A)=\max_xH_A(x),
\qquad Q(A)=\max_x|H_A(x)|.                    \tag{TS.1}
```

Let

```math
\Lambda(A)={\|A\|_(2 to2)\over\sqrt n},
\qquad
N_A(t)=#\{x:H_A(x)\ge t n^{3/2}\}.              \tag{TS.2}
```

The full symmetric Frobenius normalization is

```math
\|A\|_F^2=n(n-1).                                \tag{TS.3}
```

## 2. The exact Hanson--Wright dichotomy

### Theorem TS.1 (tail versus spectral spike)

There is an absolute constant `c_HW>0` such that every hollow symmetric
signing, every `t>0`, and every `n>=2` obey

```math
\boxed{
N_A(t)
\le2\exp\left[
n\log2-c_(HW)n\min\left\{
 {4t^2n\over n-1},{2t\over\Lambda(A)}
\right\}\right].}                              \tag{TS.4}
```

Consequently, if `Lambda(A)<=L`, then for all sufficiently large `n`,

```math
N_A(t)\le
\exp\{(\log2-\kappa(t,L))n\},                  \tag{TS.5}
```

where one may take

```math
\kappa(t,L)={c_(HW)\over2}
             \min\{4t^2,2t/L\}.                \tag{TS.6}
```

Conversely, suppose

```math
N_A(t)>\exp\{(\log2-\kappa)n\}.                \tag{TS.7}
```

If

```math
{\kappa+\log2/n\over c_(HW)}
 <{4t^2n\over n-1},                             \tag{TS.8}
```

then necessarily

```math
\boxed{
\Lambda(A)>
{2c_(HW)t\over\kappa+\log2/n}.}                \tag{TS.9}
```

Thus failure of a prescribed fixed entropy deficit is witnessed by a
quantitatively large spectral spike.  The quantifiers are uniform over the
signing; neither a maximizing spin nor exact minimization is used.

#### Proof

Let `X` be uniform on the Boolean cube.  The Rademacher Hanson--Wright
inequality gives, since `tr(A)=0`,

```math
\Pr\{|X^TAX|\ge s\}
\le2\exp\left[-c_(HW)
 \min\left\{{s^2\over\|A\|_F^2},
             {s\over\|A\|_(2 to2)}\right\}\right].          \tag{TS.10}
```

The event in (TS.2) implies `X^TAX>=2tn^(3/2)`.  Substitute this value of
`s`, use (TS.3), and multiply the probability by `2^n`; this is (TS.4).
Equation (TS.5) absorbs the prefactor two into half of the exponent.
Finally, (TS.7) and (TS.4) imply

```math
c_(HW)\min\left\{{4t^2n\over n-1},
                         {2t\over\Lambda(A)}\right\}
<\kappa+{\log2\over n}.                         \tag{TS.11}
```

Under (TS.8), the second branch must be smaller, which rearranges to
(TS.9). `square`

The sole imported result is Theorem 1.1 of Rudelson and Vershynin,
[*Hanson--Wright inequality and sub-gaussian concentration*](https://doi.org/10.1214/ECP.v18-2865),
Electronic Communications in Probability 18 (2013).  This is the same
primary theorem already used in Theorems 21.8 and 21.8a; no new literature
hypothesis is being introduced.

### Corollary TS.1a (the sharp consequence of cap control alone)

The archived asymmetric-rounding inequality, translated to the present
half-quadratic normalization, is

```math
\boxed{Q(A)\ge {1\over4}\|A\|_(2 to2)^2.}        \tag{TS.11a}
```

For completeness, if `v` is a unit eigenvector of extremal eigenvalue
magnitude and
`mu=||v||_infinity`, choose one Boolean vector as `sign(v)` and a second
with independent coordinates of means `v_i/mu`.  Their expected bilinear
form has magnitude
`||A||_(2 to2)||v||_1/mu>=||A||_(2 to2)^2`, by the eigenvector equation.
Ordinary polarization bounds the Boolean bilinear norm by `4Q(A)`, proving
(TS.11a).

Consequently, if `Q(A)<=Cn^(3/2)`, then

```math
\|A\|_(2 to2)\le2\sqrt C n^{3/4},
\qquad \Lambda(A)\le2\sqrt C n^{1/4}.            \tag{TS.11b}
```

Substitution in TS.1 gives, for every fixed `t,C>0` and all sufficiently
large `n`,

```math
N_A(t)\le2\exp\left\{n\log2-
 {c_(HW)t\over\sqrt C}n^{3/4}\right\}.           \tag{TS.11c}
```

Thus cap control alone gives a universal `Theta(n^(3/4))` entropy deficit,
not the fixed-rate `Theta(n)` deficit required by `L_tail`.  The spike
examples in Section 5 show that this exponent is scale-sharp even after
requiring complete signs (for a sufficiently wide upper window).

## 3. Exact implication for `L_tail`

Fix any `0<c_-<0.336493364431...`.  The established asymptotic lower bound
implies, for all sufficiently large `n`,

```math
M_n\ge c_-n^{3/2}.                              \tag{TS.12}
```

Let `A_n` be an exact minimizer and globally negate it when needed so that
`P(A_n)=M_n`.  Choose any fixed

```math
0<d_0<c_-,\qquad t_0=c_--d_0>0.                \tag{TS.13}
```

Its `d_0n^(3/2)` upper level set is contained in

```math
\{x:H_(A_n)(x)\ge t_0n^{3/2}\}.                 \tag{TS.14}
```

Therefore TS.1 proves the following sufficient lemma and obstruction
alternative.  It is not an equivalence: a spectral spike need not create a
large Boolean upper tail.

### Corollary TS.2 (spectral-root sufficient lemma)

If there is a fixed `L<infinity` such that every sufficiently large exact
minimizer satisfies

```math
\|A_n\|_(2 to2)\le L\sqrt n,                    \tag{TS.15}
```

then `L_tail` holds with the fixed `d_0` from (TS.13) and
`kappa=kappa(t_0,L)>0` from (TS.6).

Conversely, if `L_tail` fails at that `d_0` along exact minimizers with
entropy deficits `kappa_j downarrow0`, then TS.9 forces

```math
{\|A_(n_j)\|_(2 to2)\over\sqrt {n_j}}
\longrightarrow\infty.                          \tag{TS.16}
```

For a **fixed** failed deficit `kappa>0`, TS.9 gives a fixed large constant
lower bound, not divergence.  This distinction is important: TS.1 is a
dichotomy, not an inverse theorem identifying a unique spike.

The archive already localizes the spike forced by (TS.16).  If `v_j` is a
unit eigenvector for `|lambda_j|=||A_(n_j)||_(2 to2)`, the biased-sign
inequality (2.11) in `nearmin_deterministic_inequalities.md` says

```math
|\lambda_j|\le2Q(A_(n_j))\|v_j\|_infinity^2.
```

Since `Q(A_(n_j))=O(n_j^(3/2))`, (TS.16) implies

```math
{1\over\|v_j\|_infinity^2}
\le {2Q(A_(n_j))\over|\lambda_j|}=o(n_j).        \tag{TS.16a}
```

Thus any vanishing-deficit failure of `L_tail` is forced into an already
identified **localized** spectral-spike branch.  This localization is
pre-existing, not a new consequence claimed by this note.

## 4. Why exact minimality has not closed the spike branch

The currently proved exact-minimizer laws are:

1. `Q(A)=M_n=Theta(n^(3/2))`;
2. non-improvement under every coefficient flip;
3. the fractional near-top shell equilibrium of Theorem 36.2;
4. the asymmetric-rounding envelope (TS.11a), giving only
   `||A||_(2 to2)=O(n^(3/4))`, and the trace-norm/factorization envelope
   `||A||_*<=4K_GQ(A)=O(n^(3/2))`;
5. the individual localization inequality used in (TS.16a), and the
   stronger multi-mode Grothendieck--Pietsch common-support removal theorem;
6. the deep-hole local covering identities.

None implies (TS.15).  The strongest cap-only estimate (TS.11a) still
permits a top eigenvalue `L_nsqrt n` throughout the range
`1<<L_n<=O(n^(1/4))`, and Section 5 attains the endpoint scale.
Fractional shell balance is a law on exposed cut coordinates and supplies
no bound on a real eigenvector.
Edge-flip minimality says that some active response blocks every proposed
coefficient descent; it does not make that response correlate with the top
eigenvector.  Turning any of these facts into (TS.15) would require a new
rooted synchronization theorem.

The common-support theorem does prove that, for every `t->infinity`, all
modes above order `t sqrt n` can be removed on `O(n/t)=o(n)` vertices.  It
does not eliminate those modes.  Nor does it give a fixed-rate tail bound:
conditioning on the exceptional spins creates large cross fields, and the
available Hanson--Wright exponent and exceptional-state entropy are both on
the same non-linear scale.  This is exactly the archived selector half left
open by Grothendieck--Pietsch removal.

There is also a decisive robustness obstruction.

### Proposition TS.3 (the spectral root is not robust over the unrestricted vanishing near-minimizer class)

Let `L_n->infinity` with `L_n=o(n^(1/4))`, put
`k_n=floor(L_nsqrt n)`, and start from any exact minimizer `A_n`.  Overwrite
one `k_n`-vertex principal block by a positive clique, obtaining `B_n`.
Then

```math
Q(B_n)\le M_n+k_n(k_n-1)=M_n+o(n^{3/2}),          \tag{TS.17}
```

but

```math
\|B_n\|_(2 to2)\ge k_n-1,
\qquad
{\|B_n\|_(2 to2)\over\sqrt n}\longrightarrow\infty.       \tag{TS.18}
```

The cap bound is edit Lipschitzness.  The clique principal submatrix is
`J-I`, whose top eigenvalue is `k_n-1`; interlacing proves (TS.18).

This is the single-clique case of the archived multi-clique peeling
obstruction SH.2, not a new construction.  Its role here is exact: there
exist `o(n^(3/2))`-near-minimizer sequences on which the Hanson--Wright
sufficient condition fails.  More quantitatively, a prescribed additive
halo `epsilon_n n^(3/2)` permits this construction whenever
`epsilon_n sqrt n->infinity`, by taking
`1<<L_n<<min(n^(1/4),sqrt(epsilon_n sqrt n))`.  The proposition does **not**
cover every prescribed vanishing halo and does not falsify `L_tail` for
exact minimizers; exactness may be genuinely discontinuous.

## 5. Sharpness and bounded-cap falsifiers

### 5.1 Weighted spike: exact sharpness of the missing branch

Let `k=floor(an^(3/4))` and let a hollow weighted matrix be a positive
clique on `k` vertices and zero elsewhere.  Its quadratic cap is

```math
Q={k(k-1)\over2}=Theta(n^{3/2}),
\qquad \|A\|_(2 to2)=k-1=Theta(n^{3/4}).          \tag{TS.19}
```

For every fixed `0<eta<1`, the set above `(1-eta)Q` has size

```math
2^{n-k}\left[2^k\exp(-Theta_eta(k))\right]
=\exp\{n\log2-Theta(n^{3/4})\}.                  \tag{TS.20}
```

Indeed the energy is `((sum_(i<=k)x_i)^2-k)/2`, and Cramer/binomial large
deviations on those `k` spins give (TS.20); the other `n-k` spins are free.
Thus the operator branch in TS.1 is scale-sharp for quadratic Boolean
landscapes.  Any proof of `L_tail` must use the complete-sign and/or exact-
minimizer hypotheses to rule out this effective sublinear support.

### 5.2 Complete signs: a two-sided wide-window bounded-cap falsifier

The same phenomenon survives exact signs, although the construction does
not reach the narrow window allowed in `L_tail`.  Fix `d_0>1/2` and choose
`a>0` so small that

```math
{1\over2}+{a^2\over2}<d_0.                       \tag{TS.21}
```

Take two disjoint sets `S_+,S_-`, each of size
`k=floor(an^(3/4))`.  Put a positive clique on `S_+`, a negative clique on
`S_-`, and on the linear-sized complement put any signing of cap at most
`(1/2+o(1))n^(3/2)`.  Choose every remaining cross block with total Boolean
cap `o(n^(3/2))`.  Such a choice exists by independent signs and a union
bound: the maximum of the cross sum over all Boolean spins is
`O(n\sqrt k+k\sqrt n)=O(n^(11/8))`.

Writing `K=k(k-1)/2`, the resulting complete signing `C_n` obeys

```math
Q(C_n)\le
\left(a^2+{1\over2}+o(1)\right)n^{3/2},
\qquad
\|C_n\|_(2 to2)\ge k-1=(a+o(1))n^{3/4}.         \tag{TS.22}
```

Fix all spins of `S_+` positive and choose the other `n-k` spins uniformly.
The energy is `K+Z_+`, where orthogonality of distinct free Boolean
monomials gives

```math
\mathbb E Z_+=0,
\qquad
\mathbb E Z_+^2=O(nk^2+n^2)=o(n^3).             \tag{TS.23}
```

The `nk^2` term allows the worst possible coalescence of the `k` fixed-to-
free cross edges into each free linear monomial.  Hence Chebyshev shows that
a `1-o(1)` fraction of these extensions have energy at least `K-o(n^(3/2))`.
Likewise, after fixing `S_-` positive, a `1-o(1)` fraction of extensions
have energy at most `-K+o(n^(3/2))`.

Now orient `C_n` by whichever global sign makes its positive maximum equal
to `Q(C_n)`.  The appropriate one of these two diffuse families has oriented
energy at least `K-o(n^(3/2))`.  By (TS.22), its deficit from the oriented
maximum is at most

```math
\left({1\over2}+{a^2\over2}+o(1)\right)n^{3/2}
<d_0n^{3/2}.                                    \tag{TS.24}
```

It contains `2^{n-k}(1-o(1))=exp\{n log2-o(n)\}` spins.  Thus every fixed
`d_0>1/2` admits bounded-cap complete signings whose correctly oriented
upper level set has no fixed entropy deficit.

The constant `1/2` is the cap paid by the unavoidable linear-sized
complete-sign core, using the repository's all-order upper construction.
This is not a falsifier for `L_tail`, which may choose `d_0<c_-<1/2`, and
it is not a near-minimizer.  It is the strongest elementary complete-sign
warning obtained here: bounded cap alone does not give a tail deficit for
arbitrary fixed windows.

## 6. Why existing moment and shell inequalities stop short

For every signing,

```math
\|H_A\|_2=\sqrt{{n\choose2}}=Theta(n).
```

Degree-two hypercontractivity gives

```math
\|H_A\|_p\le(p-1)\|H_A\|_2.                     \tag{TS.25}
```

Optimizing Markov's inequality at a threshold `Theta(n^(3/2))` yields only

```math
\Pr\{|H_A|\ge tn^{3/2}\}
\le\exp(-Theta_t(\sqrt n)).                      \tag{TS.26}
```

This `sqrt n` exponent is sharp if one retains only degree and `L_2`
normalization: put edge weight `sqrt n` on a positive clique of
`floor(sqrt n)` vertices and zero elsewhere.  Then `||H||_2=Theta(n)`, its
maximum is `Theta(n^(3/2))`, and a fixed-fraction upper level has probability
`exp(-Theta(sqrt n))`.  For complete signings with cap `O(n^(3/2))`, the
stronger asymmetric-rounding input improves (TS.26) to the
`exp(-Theta(n^(3/4)))` deficit in (TS.11c), and Section 5 shows that exponent
is scale-sharp for wide windows.  Neither exponent is linear in `n`.

A fixed sixth-moment improvement gives merely a polynomial tail bound.
Hence the finite sixth-moment suppression observed in the structural audit
cannot imply `L_tail` without a growing-moment or spike-exclusion theorem.

Theorem 36.2's first-marginal shell equilibrium and Theorem 36.7's affine
near-top cube are likewise compatible with upper tails of either
`exp(o(n))` or positive exponential rate; neither supplies an upper count.
In fact Theorem 36.7 gives a universal **lower** bound of
`exp(Omega(sqrt n))` on an `O(n)`-deficit shell, which is far below the
`exp((log2-kappa)n)` upper ceiling sought in `L_tail`.

## 7. Archive comparison and verdict

- Theorem 21.8a already applies Hanson--Wright to a regular Walsh child with
  `||A||op=sqrt n`; TS.1 extracts the exact general dichotomy and converse.
- The `O(n^(3/4))` cap-only consequence (TS.11a)--(TS.11c) is an archive
  collision with the asymmetric Boolean-rounding bootstrap and the existing
  `n^(3/4)` Hanson--Wright wall.  It is included to make the sharp baseline
  explicit, not claimed as new.
- `nearmin_deterministic_inequalities.md` (2.11) already turns a large
  eigenvalue into the individual localization (TS.16a).
- `artifacts/orientation_even_grothendieck_localization.md` and
  `artifacts/multispike_grothendieck_spectral_removal.md` are stronger than
  any single-mode localization: they aggregate mixed-sign spikes and give a
  common `o(n)` exceptional support.  They localize, but do not rule out, the
  only branch left by TS.1.
- Theorem 36.4 / SH.1 gives trace-norm control and spectral peeling, but not
  a fixed operator root.  TS.3 is the archived SH.2 obstruction specialized
  to the present quantifiers.
- The near-minimizer sixth-moment audit is finite evidence only, and even a
  proved fixed-moment gap would not imply the exponential deficit.
- The universal affine-shell theorem is lower-tail-population information,
  not an upper-tail deficit.

The research judgment is therefore:

```text
L_tail: OPEN.

PROVED conditional implication:
  L_tail follows from a uniform O(sqrt n) operator norm for exact
  minimizers; failure with vanishing entropy deficit forces an unbounded
  normalized spectral spike.

NOT YET A STRICT REDUCTION:
  no existing exact-minimality theorem makes the uniform operator bound
  demonstrably simpler than L_tail, and the bound is false on the
  unrestricted vanishing near-minimizer class.

OBSTRUCTION:
  the operator-root statement is false on the unrestricted class of
  vanishing near-minimizers (and in every prescribed halo satisfying
  `epsilon_n sqrt n->infinity`), while no existing exact-minimality law
  controls the spike branch.

NEXT NON-EQUIVALENT LEMMA:
  localization itself is already archived.  Prove that exact coefficient
  optimality excludes the localized/common-support spike, or prove L_tail
  directly conditional on that spike.  Reapplying Hanson--Wright, fixed
  moments, individual localization, Grothendieck--Pietsch removal, or
  fractional shell balance cannot advance the arrow.
```
