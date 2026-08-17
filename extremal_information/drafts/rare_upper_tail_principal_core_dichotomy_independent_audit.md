# Independent audit: rare upper tails and principal cores

**Verdict on the frozen source:** **REPAIR, THEN PASS.**  The three main
implications are mathematically sound for hollow matrices with
`|a_(ij)|<=1` (in particular, complete signings), and PC.3 really does
follow with the stated diagonal choice.  The source needs one hypothesis
repair, one false robustness sentence removed, and a more conservative
frontier classification.  In particular, this is a useful inverse
concentration theorem, but not a strict reduction of `L_tail` and not by
itself Level-5 progress on exact minimizers.

## 1. Frozen source

```text
extremal_information/drafts/rare_upper_tail_principal_core_dichotomy.md
sha256 7f429ee66041b4c9fceb66037f23e37e753c6680d5762bfb5f250945509183ec
```

I used throughout the source normalization

```math
H_A(x)={1\over2}x^TAx,\qquad
Q(A)=\max_x|H_A(x)|.
```

## 2. PC.1: factor and normalization reconstruction

Put `B=||A||_(infinity to 1)`.  For Boolean `x,y`, set
`u=(x+y)/2` and `v=(x-y)/2`.  The two vectors are in the unit cube and
have disjoint supports, whence symmetry gives

```math
x^TAy=u^TAu-v^TAv.
```

Hollowness makes `z^TAz` multiaffine, so its absolute maximum on the cube
is attained at a Boolean corner and equals `2Q(A)`.  Therefore

```math
\boxed{B\le4Q(A).}
```

This is the correct factor in the half-quadratic normalization.  It agrees
with the archived `B<=2Q_full` when `Q_full=max|x^TAx|=2Q(A)`.

The Grothendieck--Pietsch factorization used in the source gives

```math
|u^TAv|\le K_GB
 (\sum_i\mu_i u_i^2)^{1/2}
 (\sum_j\nu_j v_j^2)^{1/2}.
```

Each of the two heavy sets in PC.7 has cardinality strictly less than
`epsilon*n/2`; hence their union has size less than `epsilon*n`.  On the
complement `R`, both weights are at most `2/(epsilon*n)`, and so

```math
\|A[R]\|_(2 to 2)
\le {2K_GB\over\epsilon n}
\le {8K_GC\over\epsilon}\sqrt n.
```

For `h=A_(R,T)x_T`, dualizing the Euclidean norm over vectors supported on
`R` uses the light `mu` factor on `R` and only the trivial bound
`sum_(j in T)nu_j<=1` on the fixed shore.  Thus

```math
\|h\|_2\le K_GB\sqrt{2\over\epsilon n}
\le {4\sqrt2K_GC\over\sqrt\epsilon}\,n.
```

Both displayed constants in PC.8--PC.9 are therefore correct.

Conditionally on `x_T`, the exact decomposition is

```math
H_A(x_T,X_R)=H_(A[T])(x_T)+h^TX_R+H_(A[R])(X_R).
```

If `P(A[T])<=(t-eta)n^(3/2)`, an upper-tail event forces one of the latter
two terms to exceed `(eta/2)n^(3/2)`.  The Rademacher linear tail and the
last norm estimate give

```math
Pr\{h^TX_R\ge(eta/2)n^(3/2)\}
\le \exp[-c eta^2\epsilon n/C^2].
```

For a hollow sign/bounded-entry matrix,
`||A[R]||_F^2<=n^2`.  Hanson--Wright, applied to `X_R^TA[R]X_R` at the
equivalent threshold `eta*n^(3/2)`, then gives

```math
2\exp[-c n\min\{eta^2,eta\epsilon/C\}].
```

No independence between the linear and quadratic pieces is needed; a
union bound suffices.  The estimates are uniform in the frozen spin, so
averaging over `x_T` gives PC.5 with the prefactor three.

### Required hypothesis repair

The setup currently says only "real symmetric hollow matrix", but the
proof uses `||A[R]||_F^2<=n^2`.  Add

```text
|a_(ij)|<=1
```

to PC.1 (or state it only for hollow signings).  Without an entry/Frobenius
hypothesis that displayed step is unjustified.  A general-real variant is
possible, but has a weaker Frobenius branch: PC.8 and
`||A[R]||_F^2<=n||A[R]||_op^2` give a term of order
`eta^2 epsilon^2/C^2`, not the source's `eta^2`.  The signing application
does not need that variant.

## 3. PC.2: variance and counting reconstruction

Let `k=|T|` and freeze a maximizing core spin.  On the free shore write

```math
Z=h^TX_R+H_(A[R])(X_R).
```

Distinct degree-one and degree-two Walsh characters are orthogonal, so

```math
E Z=0,
\qquad
E Z^2=\sum_(j in R)h_j^2+\binom{|R|}{2}
\le nk^2+n^2.
```

The `nk^2` term is the correct worst-case cross-field variance: as many as
`k` fixed-to-free edges can coalesce into one linear coefficient.  If
`k=o(n)`, the variance is `o(n^3)`.  Chebyshev therefore says that a
`1-o(1)` fraction of the `2^(n-k)` extensions lose less than the fixed
margin `eta*n^(3/2)`.  Dividing by `2^n` proves PC.15.  The factor
`2^(-|T|)` and all normalizations are correct.

## 4. PC.3: the diagonal argument really works

Assume `r_n=-(1/n)log p_(A_n)(t)->0` and put

```math
delta_n=max(r_n,1/n),\qquad
epsilon_n=eta_n=delta_n^(1/4).
```

For small `delta_n`, the minimum in PC.5 is bounded below, with constants
depending only on fixed `C`, by

```math
c_C delta_n^(3/4).
```

Indeed the three powers are `delta^(1/2)`, `delta^(1/2)`, and
`delta^(3/4)`.  Moreover

```math
{r_n\over delta_n^(3/4)}\le delta_n^(1/4)->0,
\qquad
n delta_n^(3/4)\ge n^(1/4)->infinity.
```

The second relation absorbs the prefactor three.  Consequently PC.5 is
incompatible with the observed probability for all large `n`; PC.4 must
fail for the factorization heavy set.  That set has
`|T_n|<delta_n^(1/4)n=o(n)`, and failure of PC.4 gives

```math
P(A_n[T_n])>(t-delta_n^(1/4))n^(3/2).
```

This is PC.18.  Finally `P(A[T])<=binom(k,2)` for a `k`-vertex signing, so

```math
k\ge(\sqrt{2t}-o(1))n^(3/4),
```

which verifies PC.19.  Thus there is no hidden diagonal-choice or
quantifier failure in PC.3.

## 5. Repairs to the exact-minimizer consequence

The phrase "every `T=o(n)`" is not a pointwise finite-order quantifier.
State `L_core` sequentially:

> there is no sequence of suitably oriented exact minimizers `A_(n_j)` and
> sets `T_j` with `|T_j|/n_j->0` for which
> `P(A_(n_j)[T_j])>=(t_0-o(1))n_j^(3/2)`.

This is the weakest conclusion directly matched to PC.3.  If the stronger
fixed-gap version is retained, require explicitly
`0<eta_0<t_0` and formulate its uniform/sequential quantifiers.  The
displayed fixed-gap condition at `(t_0-eta_0)` is stronger than needed, not
the "exact" minimal obstruction.

There is also a false robustness sentence in Section 5.  The archived
single- and multi-clique overwrites in an `o(n^(3/2))` near-minimizer halo
create operator spikes whose total planted principal energy is
`o(n^(3/2))`.  They do **not** implant a fixed-positive-level PC.18 core.
A clique carrying `t n^(3/2)` for fixed `t>0` has size
`Theta(n^(3/4))` and costs `Theta(n^(3/2))` edge edits, not
`o(n^(3/2))`.  Therefore remove the claim that the archive proves the
fixed-threshold `L_core` quantifier "necessarily discontinuous" over every
vanishing halo.  The archive proves discontinuity for uniform spectral
regularity, not for this principal-energy lemma.

## 6. Archive comparison and frontier classification

The following ingredients are already archived:

- `B<=4Q` in the present normalization;
- the two-measure Grothendieck--Pietsch heavy-coordinate deletion and its
  `O(sqrt(n)/epsilon)` complementary operator norm;
- Hanson--Wright upper tails and the tail-versus-spectral-spike theorem;
- sublinear planted clique/operator-bubble examples;
- endpoint and operator concentration--compactness language.

I found no archived statement proving the conditional cross-field estimate
and the resulting zero-rate-tail `=>` global-scale principal-core theorem.
PC.1--PC.3 are therefore a new assembly and a useful inverse theorem, not
an exact duplicate.  Their genuinely new point is that conditioning on the
Pietsch heavy set controls the entire cross field strongly enough at a
diagonal scale, so a zero-rate upper tail must leave leading energy inside
the heavy principal block.

The classification in the source is nevertheless too strong in two ways.

1. PC.2--PC.3 themselves show that, up to a threshold slack, zero-rate
   tails and sublinear principal cores are converse descriptions of the
   same phenomenon.  Thus `L_core` is a sharper structural
   **reformulation/inverse characterization**, not yet a demonstrated
   strict reduction of `L_tail` in difficulty.
2. No theorem here proves the no-core property for exact minimizers, and
   the archived vanishing-halo examples do not falsify it at fixed energy
   level.  It is safe to call PC.1--PC.3 a rigorous generic bounded-cap
   structural dichotomy (or a conditional Level-5 diagnostic under the
   project's loose benchmark vocabulary), but it should not reset the
   Level-5 near-minimizer frontier.  The exact-minimizer branch remains
   open, and the scalar all-spins-free selector remains separately open.

## 7. Disposition

After adding the bounded-entry hypothesis, repairing the `L_core`
quantifiers, and deleting the unsupported halo-discontinuity claim:

```text
PC.1: PASS.
PC.2: PASS.
PC.3: PASS.
Archive novelty: PASS as a new corollary/structural inverse theorem.
Strict-reduction claim: REJECT.
Level-5 frontier-reset claim: REJECT.
```

## 8. Frozen repaired-source recheck

The author added the bounded-entry hypothesis, replaced the finite-order
`T=o(n)` wording by the minimal sequential `L_core` formulation, explicitly
disclaimed a proved strict reduction/frontier reset, and corrected the
planted-halo discussion to distinguish subleading principal energy from a
fixed positive global-scale core.  I rechecked the repaired file in full.

```text
extremal_information/drafts/rare_upper_tail_principal_core_dichotomy.md
sha256 43bf9ac6fe40bd99a3c38bf791a355d0f30d2a79d5ff1453643356238d806318
final verdict PASS
```
