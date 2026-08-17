# Independent audit: witness transversals and orbit-query large deviations

**Files audited**

- `drafts/extremal_witness_transversals.md`;
- `drafts/orbit_query_large_deviations.md`;
- `experiments/verify_extremal_witness_transversals.py`;
- `experiments/verify_orbit_query_large_deviations.py`.

**Verdict:** PASS after two minor finite-statement repairs and one explicit
hypothesis repair.  The group-cover directions, tail exponent, flat
exact-sign realization, coding equivalence, Cramer/type argument, and
heterogeneous composition formula are mathematically sound.  No claimed
exponential rate changes.

## 1. The group-cover theorem

For `G={+-1}^k`, condition WT.2 says that for every `s` there are
`x in X` and `w in W_alpha` with `sx=w`.  Since every group element is its
own inverse,

```math
s=wx,
```

so the condition is exactly

```math
G=W_alpha X.                                       \tag{A.1}
```

The draft therefore has the product-set direction right:
`|G|<=|W_alpha||X|`.  Conversely, a uniform point lies in any fixed
translate with probability `p_alpha`; `q` independent samples miss it with
probability at most `exp(-p_alpha q)`.  Union over `2^k` translates proves
the displayed ceiling in WT.4.

The additive version in LD.4 similarly gives `G_n=W_n-X_n`; commutativity
is more than enough, and cardinality gives the same lower bound.

### Minor repair A

The upper bound LD.6 should contain a ceiling:

```math
L_n(a)\le
\left\lceil{n\log q+1\over p_n(a)}\right\rceil.    \tag{A.2}
```

Random sampling proves (A.2), not literally the unrounded real-valued upper
bound currently printed.  This changes `log L_n` by `O(1)` and has no effect
on LD.7.

## 2. The flat block example and exact-sign realization

For `k=s^2`, direct expansion gives

```math
H_D(x)=\sum_{b=1}^s(M_b^2-s),
\quad
\max H_D=s(s^2-s)=s^3-s^2,
\quad
\min H_D=-s^2                                      \tag{A.3}
```

when `s` is even.  Thus `Q=s^3-s^2` for every even `s>=2` (equality of the
two absolute endpoints occurs at `s=2`).

The moment calculation is correct.  The Gaussian identity gives, for
`0<lambda<1/2`,

```math
\mathbb E\exp(\lambda M_b^2/s)
=\mathbb E_Z\prod_{i=1}^s
 \cosh(\sqrt{2\lambda/s}\,Z)
\le\mathbb E_Ze^{\lambda Z^2}
=(1-2\lambda)^{-1/2}.                              \tag{A.4}
```

The event `H_D>=alpha(s^3-s^2)` implies

```math
\sum_b M_b^2/s
\ge\alpha s^2+(1-\alpha)s.                        \tag{A.5}
```

At `lambda=1/4`, (A.4)--(A.5) give

```math
\Pr\{H_D\ge\alpha Q\}
\le\exp\{-\alpha s^2/4+O_alpha(s)\}.             \tag{A.6}
```

The negative absolute tail is empty once `alpha(s^3-s^2)>s^2`.  This proves
WT.7 with exactly the claimed leading constant.

The exact-sign realization in WT.9 is also valid.  Let `K` be block
diagonal with `J_s-I_s` on every block, and let `R` be symmetric, zero on
those diagonal blocks, and independent Rademacher between blocks.  Then

```math
A=R+K,
\qquad A'=R-K,
\qquad A-A'=2K=D.                                  \tag{A.7}
```

Both `A,A'` are hollow complete sign matrices.  Also `||K||=s-1`.  For a
fixed unit vector `v`, `v^TRv` is subgaussian with an absolute variance
proxy; an `epsilon`-net union bound over the sphere therefore gives
`||R||=O(sqrt(k))=O(s)` with positive probability.  The same one `R` occurs
in both matrices, so

```math
||A||+||A'||\le2||R||+2||K||=O(s).                 \tag{A.8}
```

Thus there is no hidden use of two independently favorable completions.
Switching is orthogonal conjugation and preserves all stated norms and the
contrast.

## 3. Coding corollary

Take the base distance landscape `f(x)=d_H(x,1)`.  Then
`f_s(x)=f(sx)=d_H(x,s)`.  Since

```math
d_H(-x,s)=k-d_H(x,s),                              \tag{A.9}
```

the requirement `d_H(x,s)>=alpha k` is exactly

```math
d_H(-x,s)\le k-\lceil\alpha k\rceil
=\lfloor(1-\alpha)k\rfloor.                       \tag{A.10}
```

Hence `-X` is a covering code with precisely the radius claimed.  The
sphere-cover lower bound is correct, as is the random-cover exponent.

### Minor repair B

The exact upper inequality WT.14 should be printed with a ceiling (or an
additive `+1`):

```math
|X|\le
\left\lceil
{(k\log2+1)2^k\over
 \sum_{j\le\lfloor(1-\alpha)k\rfloor}{k\choose j}}
\right\rceil.                                      \tag{A.11}
```

Again this has no effect on WT.15.  For fixed `alpha>1/2`, the standard
binomial-volume estimate gives the stated `O(log k)` remainder.

## 4. Cramer/type proof

Let `mu` be the uniform pushforward law of `f`.  Chernoff proves

```math
\liminf_n-{1\over n}\log p_n(a)\ge I_f(a).         \tag{A.12}
```

For the other direction, a type `nu_n` has probability
`exp(-nD(nu_n||mu)+O(log n))`, with the `O(log n)` uniform because the
alphabet is fixed.  The constrained entropy minimum

```math
\inf\{D(\nu||\mu):\mathbb E_\nu f\ge a\}           \tag{A.13}
```

equals `sup_(theta>=0){theta a-Lambda_f(theta)}`.  Because
`E f<a<max f`, a rational type can approximate the minimizer from the
feasible side (move `O(1/n)` mass toward a maximizer if rounding initially
falls below `a`).  One such type is contained in the upper tail, proving

```math
\limsup_n-{1\over n}\log p_n(a)\le I_f(a).         \tag{A.14}
```

Thus LD.10 and the orbit-query exponent are correct.  The prose after LD.8
should ideally call Chernoff the **lower bound on the rate**, not “the upper
rate”; the formulas themselves have the right direction.

## 5. Heterogeneous composition

The formula is valid, but its theorem-level hypotheses should be made
explicit.  Let `f` and `g` live on finite abelian groups (possibly of
different orders), let `N=m+n`, and suppose

```math
{m\over N}\longrightarrow\lambda.                 \tag{A.15}
```

For a nondegenerate upper-tail threshold require

```math
\lambda\mathbb Ef+(1-\lambda)\mathbb Eg
<a<
\lambda\max f+(1-\lambda)\max g.                  \tag{A.16}
```

The exact normalized log moment is

```math
{m\over N}\Lambda_f(\theta)+{n\over N}\Lambda_g(\theta),
                                                               \tag{A.17}
```

which converges to LD.11.  Equivalently, the two-type method minimizes

```math
\lambda D(\nu_f||\mu_f)
 +(1-\lambda)D(\nu_g||\mu_g)                       \tag{A.18}
```

under the weighted energy constraint.  Convex duality gives exactly LD.12.
The group-transversal factor is at most linear in `N` before taking its
logarithm, even when the two group orders differ.  Hence

```math
{1\over N}\log L_{m,n}(a)\longrightarrow I_\lambda(a).       \tag{A.19}
```

### Hypothesis repair C

Section 2 should say that `f,g` are finite **group landscapes**, add
(A.16), and state (A.19).  Without (A.16), LD.12 still defines a convex
quantity, but it need not describe a nonempty exponentially small upper
tail.  With these hypotheses the heterogeneous composition claim passes.

## 6. Verifier audit

Both supplied programs run successfully in the project virtual environment.

The witness-transversal verifier correctly checks:

- the product-set lower-bound direction for every nonempty `W,X` through
  cube dimension three;
- exact extrema and tail cardinalities of the block family for `s=2,4`.

It does not attempt to certify the probabilistic covering upper bound or the
asymptotic random-matrix norm statement.  That is acceptable because both
are analytic existence proofs, but the verifier description should not be
read as testing WT.4 or WT.9.

The large-deviation verifier exactly counts binary tails at
`n=20,40,80,160` and confirms convergence to the KL rate for three
thresholds.  It tests the Rademacher benchmark, not the general type proof or
the heterogeneous law.  Again the scope is appropriate, though an optional
heterogeneous finite test would improve regression coverage.

## 7. Final classification

| Item | Verdict | Required action |
|---|---|---|
| WT.1 lower/product direction | PASS | none |
| WT.1 random-cover direction | PASS | none |
| WT.2 tail constant | PASS | none |
| WT.9 flat exact signs | PASS | optionally expand the net argument |
| coding antipode/radius | PASS | add ceiling to exact upper bound |
| LD.1 Cramer exponent | PASS | clarify “lower bound on rate” wording |
| LD.6 finite upper bound | PASS after rounding | add ceiling |
| heterogeneous composition | PASS after hypotheses | state group domains, threshold range, and limit (A.19) |

The mathematical contribution survives intact: reciprocal extremal mass is
the exact orbit-query rate up to subexponential factors, and log-moment
functions form a genuine composition algebra for heterogeneous independent
products.  The result remains deliberately scoped to coordinatewise group
translations; no claim about interacting composition is smuggled into the
proof.
