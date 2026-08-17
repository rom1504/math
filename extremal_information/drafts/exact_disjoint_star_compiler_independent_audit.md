# Independent audit: exact disjoint star compiler barrier

**Verdict:** PASS, with two scope clarifications and no mathematical repair
required.  I checked Propositions/Theorems SC.1--SC.6, including the later
full-support Fourier inversion, interacting-selector, exposed-set, and
bounded-cap additions.  I also ran
`experiments/verify_exact_disjoint_star_compiler.py`; all finite checks pass.

This audit was performed after freezing an independent compiler hierarchy:
a sparse exact-edge identity; a complete-sign contextual compiler modulo a
query-only baseline; and a rate-preserving compiler, which must use
`N=O(k)` vertices if a `Theta(k^(3/2))` response is to remain at positive
`N^(3/2)` scale.  Thus the note's distinction between sparse compilation and
complete-sign closure is substantive rather than terminological.

## 1. Sparse identity and normalization

For one edge,

```math
\max_{y=\pm1}y(x_i-T_{ij}x_j)
=|x_i-T_{ij}x_j|=1-T_{ij}x_ix_j.
```

Summation proves (SC.2), and adding the old energy proves (SC.3).  No factor
of two is missing under the convention
`H_T=\sum_{i<j}T_{ij}x_ix_j`.  The construction is genuinely disjoint and
uses exact unit coefficients on its present edges, but it is sparse and
one-sided.  With `N=k+binom(k,2)`, a `k^(3/2)` contrast is indeed
`Theta(N^(3/4))`, hence has normalized size `Theta(N^(-3/4))`.

## 2. Pair Fourier coefficient and Gram barrier

Conditioning on the other `d-2` spins gives (SC.8).  When `d` is even it is
the indicator that their sum is zero; when `d` is odd it is one half of the
indicator that the sum has absolute value one.  Both cases equal

```math
\gamma_d=2^{-(d-2)}{d-2\choose\lfloor(d-1)/2\rfloor}.
```

The signs contribute the factor `sigma_i sigma_j`.  The quoted central
binomial estimate gives
`gamma_d<=1/sqrt(d-1)`, and hence
`gamma_d binom(d,2)<=d^(3/2)/2`; the constants in (SC.7)--(SC.10) follow.

Exact equality of functions in (SC.4) forces equality of their pair Fourier
coefficients.  Therefore

```math
K=\sum_a\gamma_{d_a}v_av_a^T
```

is positive semidefinite, has rank at most `m`, and has the required
off-diagonal entries `-T_ij`.  For `T_ij=1`, every completion is
`diag(p)-11^T` with `p_i>0`; its kernel is at most one-dimensional, so
`rank K>=k-1`.  This verifies the PSD completion and all scale consequences.

## 3. Every even Fourier level and full-support stars

For `k=2n`, the cube Laplacian convention with level-`s` eigenvalue `2s`
does satisfy

```math
L|\sum_i x_i|=2|\sum_i x_i|-4n1_{\{\sum_i x_i=0\}}.
```

The middle-slice Krawtchouk coefficient in (SC.15c) is the standard identity

```math
2^{-2n}{2n\choose n}(-1)^r{n\choose r}/{2n\choose2r},
```

so comparison at level `2r` gives (SC.15d).  Averaging the odd-dimensional
function over its last coordinate yields `f_(2n)+I_0`, which gives
(SC.15e).  The displayed factors never vanish in the stated ranges.

After antipodal symmetrization, Fourier inversion of the empirical star law
has only its constant and pair levels left.  It is exactly

```math
\mu(\sigma)=2^{-k}\left(1-{H_T(\sigma)\over m\gamma_k}\right).
```

Nonnegativity gives (SC.15f).  Applying this separately to `T` and `-T`
gives `m gamma_k>=Q(T)`.  Lemma SC.6 supplies the uniform
`Q(T)=Omega(k^(3/2))` lower bound and `gamma_k=Theta(k^(-1/2))`, hence the
quadratic conclusion (SC.15g).  Allowing the symmetrized law to have
half-integral empirical masses causes no issue: only positivity and its
moments are used.

## 4. Interacting selectors

At a balanced `x`, an affine selector active for `F` lies below the envelope
at both antipodes.  Each comparison is at least
`a=k^2/2-2eta`; adding and dividing by two gives
`x dot b(y)>=a`.  This verifies the factor in (SC.23).  Hoeffding with
`||b(y)||_2^2<=D^2`, followed by a union bound over at most `2^m` selectors
and `binom(k,k/2)>=2^k/(k+1)`, gives (SC.20).

For a complete bipartite block, `D^2=km^2`; solving the resulting cubic
gives exactly

```math
m/k>=(8\log2)^{-1/3}-o(1)
```

and solving it for `eta/k^2` gives (SC.22).  The exposed-set version SC.5a
uses the same two-antipode comparison and is correct:

```math
\log K\ge\log|X|-k\log2+{a^2\over2D^2}.
```

For slopes `By`, `D<=||B||_(2->2)sqrt(m)` is immediate.

The oscillation claim SC.5b is also correct.  Uniformly over two old-spin
inputs, each affine slope changes by at most `2||By||_1`, so
`osc F<=2||B||_(infinity->1)`.  Flipping the entire auxiliary shore changes
the sign of `x^TBy` while preserving both internal quadratic energies;
therefore `||B||_(infinity->1)<=Q(P)`.  At `N=O(k)` and bounded normalized
cap this is only `O(k^(3/2))`, whereas `osc C_+=k^2/2`.

## 5. Complete-sign cap

With `|U|=N/3+O(1)`, Khintchine gives a cross value at least

```math
|V|\sqrt{|U|/2}
=\left({2\over3\sqrt6}+o(1)\right)N^{3/2}.
```

Flipping all spins on `V` reverses that cross term and preserves both shore
energies, so one of the two complete quadratic values has at least this
absolute magnitude.  This verifies (SC.25) and the perturbative-completion
scale claim.  The draft's proof now references the correct label (SC.25).

## 6. Scope clarifications

1. SC.3 and SC.4b classify **independent-star** selectors.  They do not
   cover mixed-support interacting auxiliary spins.  SC.5 covers arbitrary
   interactions only for the high-oscillation all-positive cut shell, and
   SC.5a requires a large antipodally exposed set to transfer that argument
   to another target.
2. SC.5b rules out uniformly approximating the cut-shell future by a
   bounded-cap linear-overhead parent.  It does not rule out a contextual
   response difference riding on a large common background, nor does it
   yet apply to the flat alternating-form family whose relevant oscillation
   is only `Theta(k^(3/2))`.

Subject to those already explicit qualifications, the results are rigorous
and theorem-level.  They decisively show why the sparse identity is not yet
an exact-complete, rate-preserving compositional closure theorem.
