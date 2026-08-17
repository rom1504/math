# Independent audit of the Level-5 near-minimizer theorems

Date: 2026-08-17.

## 0. Scope and verdicts

I independently checked

- `nearmin_pinned_response_packing.md` (Theorem PR.1),
- `nearmin_spectral_harmonic_report.md` (Theorems SH.1 and SH.2), and
- the current post-audit version of `near_minimizer_fractional_balance.md`
  (Theorems FB.1, FB.3, and Corollary FB.2).

The checks covered the matrix-versus-edge normalization, all asymptotic and
uniform quantifiers, the random-code and perturbation union bounds, the
independence used in pairwise concentration, the amplitude-`n` pinning
identity, the one-hot information inference, the exact symmetric
Grothendieck factorization and its Schatten consequences, and the
principal-compression step in the multi-clique obstruction.

| Result | Verdict | Required action |
|---|---|---|
| PR.1 | **PASS** | No mathematical repair. Choose representatives of the projective code explicitly. |
| SH.1 | **PASS** | No mathematical repair. Make the Tropp theorem number version-specific. |
| SH.2 | **PASS** | No mathematical repair. Treat an empty surviving clique block separately in one proof sentence. |
| FB.1 | **PASS** | The prior audit's normalization and strict-inequality repairs are present and consistent. |
| FB.2 | **PASS** | No repair. |
| FB.3 | **PASS** | The prior audit's quantifier, integrality, rounding, and interpretation repairs are present. |

Thus none of the audited Level-5 claims needs weakening, and none is a
FAIL.  The suggested edits below are copy-ready statement/proof hygiene;
they do not alter a displayed theorem bound.

## 1. Theorem PR.1: PASS

### 1.1 Code rate and quantifiers

For independent uniform `u,v in {+-1}^n`, `u dot v` is a sum of `n`
independent Rademachers, so Hoeffding gives

```math
 Pr\{|u\mathbin\cdot v|>n/2\}\le 2e^{-n/8}.
```

Sampling `M=floor(e^(c_0n))` vectors and union-bounding over fewer than
`M^2/2` pairs succeeds whenever `2c_0<1/8`; the draft's
`c_0<1/32` is safely inside this range.  A repeated or antipodal pair has
absolute inner product `n`, so the same good event already excludes
projective repetitions.  Consequently the passage to the quotient loses
no exponential rate.

The concentration failure probability later is
`exp(-c kappa n^(3/2))`.  At the lower endpoint
`kappa sqrt(n)>=C`, this is `exp(-c Cn)`.  After fixing any sufficiently
small absolute code rate `c<c_0`, one can choose the absolute constant `C`
so that union bounds over `e^(cn)` vertices and `e^(2cn)` ordered pairs
close.  The upper restriction `kappa<=1/10` merely ensures a valid sparse
flip probability.  If the interval in (PR.2) is nonempty, then
`n>=(10C)^2`, so choosing `C` large also absorbs the initial “all large
`n`” in the code lemma.  The theorem therefore has the stated all-order
quantifier without an unmentioned exceptional finite range.

The quotient notation does require a representative choice because `u`
and the query `nu` are not functions of a projective class as literal
vectors.  (By contrast, `z_u` is invariant under `u -> -u`.)  This is only
a notation issue: choose one representative of every class once and use it
throughout.

### 1.2 Cancellation and independence

Write `g_u(e)=a_eu_iu_j`.  Since `F_u subset N_u={g_u=-1}`,

```math
 H_(b^u)(u)=H_a(u)+2|F_u|,
 \qquad
 H_(b^v)(u)=H_a(u)-2\sum_(e\in F_v)g_u(e),
```

which proves (PR.8) with the displayed signs.  Also

```math
 2q\sum_e1_(N_v)(e)g_u(e)
 =q\left(\sum_eg_u(e)-\sum_eg_u(e)g_v(e)\right)
 =q\left(H_a(u)-\langle z_u,z_v\rangle\right).
```

Thus the base energy cancels exactly.  The edge-vector normalization is

```math
 \langle z_u,z_v\rangle
 ={(u\mathbin\cdot v)^2-n\over2},
```

so

```math
 q(E-\langle z_u,z_v\rangle)
 ={q\over2}(n^2-(u\mathbin\cdot v)^2)
 \ge {3\over8}\kappa n^{3/2}.
```

There is no hidden dependence between the two channels used in this
calculation.  The construction samples the entire families `F_u` and
`F_v` independently when `u ne v`; within either family its edge
indicators are independent.  Even when the same unordered edge belongs to
both `N_u` and `N_v`, its two indicators are different independent random
variables.  Accordingly the centered version of (PR.8) is a sum of
independent variables of magnitude at most two and variance at most

```math
 4q(|N_u|+|N_v|)\le8qE\le4\kappa n^{3/2}.
```

Bernstein at a downward displacement of at least
`(3/16)kappa n^(3/2)` gives exponent
`Omega(kappa n^(3/2))`, as claimed.  The size bound also has the advertised
exponent: `|F_u|` is dominated by `Bin(E,q)`, whose mean is less than
`kappa n^(3/2)/2`, while the threshold in (PR.13) is
`kappa n^(3/2)`.  A standard multiplicative Chernoff bound therefore gives
`exp(-Omega(kappa n^(3/2)))`.  The simultaneous ordered-pair and vertex
union bound is valid.

Each changed unordered edge changes every `H` value by at most two.  Hence

```math
 Q(b^u)\le Q(a)+2|F_u|
 \le M_n+2\kappa n^{3/2},
```

with no missing factor of two or matrix double-counting.

### 1.3 Exact pinning at amplitude `L=n`

If `d=d_H(x,u)`, exactly `d(n-d)` unordered pairs have
`x_ix_j ne u_iu_j`.  Therefore

```math
 H_b(x)-H_b(u)\le2d(n-d).
```

Meanwhile

```math
 (nu)\mathbin\cdot x-(nu)\mathbin\cdot u=-2nd.
```

The total objective difference is at most

```math
 2d(n-d)-2nd=-2d^2\le0.
```

Thus `u` really is a global maximizer at the exact boundary value `L=n`,
including the cases `d=0` and `d=n`, and

```math
 R_b(nu)=n^2+H_b(u).
```

No strict-amplitude margin is needed.

### 1.4 One-hot separation really gives `Omega(n)` bits

The conclusion does not require independently writable coordinates.  If
two different landscapes `b^u,b^v` shared the same summary state, then at
the single common decoder query `nu` they would receive the same decoded
number.  Uniform error `<Delta/2`, where
`Delta=delta kappa n^(3/2)`, would imply by the triangle inequality that
their true responses differ by less than `Delta`.  This contradicts the
oriented inequality

```math
 R_(b^u)(nu)-R_(b^v)(nu)\ge Delta.
```

Therefore the `|U|>=e^(cn)` landscapes require distinct states and at least
`log_2|U|=Omega(n)` bits.  The fact that the separating query depends on the
ordered pair through its first state is harmless: the decoder is assumed
to answer every query in the common bank `{nu:u in U}`.  The draft's
classification as a one-hot packing, rather than a Boolean cube, is exact.

### Copy-ready notation repair for PR.1

Replace the first sentence of Section 2 by:

> Fix `rho=1/2`.  For all sufficiently large `n`, the standard random-code
> argument gives a set of projective classes `U subset {+-1}^n/{+-1}` with
> `|U|>=e^(c_0n)` and `|u dot v|<=rho n` for distinct classes.  Choose and
> fix one representative `u in {+-1}^n` of each class; all subsequent
> occurrences of `u`, `nu`, and `b^u` refer to these representatives.

## 2. Theorem SH.1: PASS

### 2.1 Polarization and the constant four

For a principal set `U` and Boolean `s` on `U`, extend `s` by independent
unbiased spins off `U`.  Conditional expectation of the full quadratic
energy is `H_(A[U])(s)`, so

```math
 Q(A[U])\le Q(A).
```

For Boolean `x,y`, the draft's disjoint-support decomposition
`x=z+w`, `y=z-w` gives, using symmetry,

```math
 x^TAy=z^TAz-w^TAw.
```

If `U=supp(z)`, then `z^TAz=2H_(A[U])(z_U)`, and similarly for `w`.
Each term therefore has modulus at most `2Q(A)`.  This proves

```math
 \|A\|_(infinity->1)\le4Q(A)
```

under exactly the draft's convention
`H_A(x)=sum_(i<j)a_ijx_ix_j`; no factor of two is missing.  The definition
using Boolean `x,y` is the usual `infinity->1` norm because both cube
optimizations attain their extrema at sign vectors.

### 2.2 Symmetric Grothendieck factorization and trace norm

The cited Grothendieck factorization supplies, for real symmetric `A`, a
nonnegative diagonal `D` with `tr(D^2)=1` and a common-side factorization

```math
 A=DTD,
 \qquad \|T\|_(2->2)\le K_G\|A\|_(infinity->1).
```

The factor `T` may be taken symmetric.  On the positive support of `D`, it
is forced to equal `D^(-1)AD^(-1)`, which is symmetric; zero rows and
columns can be set to zero.  Hence the draft's symmetric formulation is
legitimate, not an extra uncited strengthening.  Combining this with the
previous paragraph gives the constant `4K_GQ(A)` in (SH.2).

The Schatten calculation is also exact:

```math
 \|DTD\|_*
 \le\|D\|_F\|TD\|_F
 \le\|D\|_F^2\|T\|_(2->2)
 =\|T\|_(2->2).
```

Thus (SH.3) follows with constant `4K_G`; it does not require a hidden
dimension factor.

There is one bibliographic ambiguity, but no theorem gap.  Tropp's longer
arXiv/technical-report version labels the Grothendieck factorization
Theorem 5.3.  The published SODA version labels the factorization Theorem
5.2 and uses Theorem 5.3 for the semidefinite characterization.  Since the
draft cites both the arXiv link and the SODA pagination, it should identify
the version explicitly.

### 2.3 Peeling and spectral tail

For

```math
 S_L={i:d_i^2>L/n},
```

the identity `sum_i d_i^2=1` gives `|S_L|<n/L`.  On its complement `C`,

```math
 A[C]=D_C T[C]D_C,
 \qquad
 \|A[C]\|_(2->2)
 \le {L\over n}\|T[C]\|_(2->2)
 \le {4K_G LQ(A)\over n}.
```

Principal compression cannot increase operator norm (and here `T` is
symmetric in any case).  This proves (SH.4) for every real `L>=1`, including
`L>n`.  Since the trace norm of a symmetric matrix is
`sum_j|lambda_j|`, every eigenvalue above `t sqrt(n)` consumes more than
that amount of trace norm, proving (SH.5) exactly.

The random-sign upper bound on `M_n` also has the right normalization.  A
fixed `H_A(x)` is a sum of `N=binom(n,2)` Rademachers, so

```math
 Pr\{|H_A(x)|\ge t\}\le2e^{-t^2/(2N)}.
```

At `t^2=2N(n+2)log2`, union over `2^n` spins has total failure probability
at most `1/2`.  This gives (SH.6), then (SH.7)--(SH.8) directly.

Finally, one changed unordered sign edge contributes
`+-2(e_ie_j^T+e_je_i^T)`, whose two nonzero singular values are both two,
so its trace norm is four.  Hence (SH.9) and the `8K_G Lr/n` enlargement of
the peeling bound are both correct.  This trace-norm robustness does not
assert stability of the particular exceptional set, and the draft
correctly says so.

The near-minimizer classification is appropriately weak: SH.1 applies to
the entire `Q=O(n^(3/2))` class and yields neither an active-state
classification nor target-scale low rank.  At cutoff `L_n sqrt(n)` the
tail rank becomes `o(n)` only for `L_n->infinity`, while the associated
uniform Boolean error scale is `L_n n^(3/2)`.  No response or convergence
claim is being smuggled into the theorem.

### Copy-ready citation repair for SH.1

Replace the citation sentence by:

> This is Theorem 5.3 in Tropp's arXiv/technical-report version (Theorem
> 5.2 in the published SODA version); Section 5.3 of the published version
> gives the semidefinite construction.

## 3. Theorem SH.2: PASS

At most `m binom(k,2)` unordered edges are overwritten, and an actual sign
change alters each Boolean quadratic energy by two.  Thus

```math
 Q(B)\le Q(A)+2m\binom k2=M_n+mk(k-1),
```

so (SH.10) has the correct edge rather than symmetric-matrix
normalization.

Let `C=V\setminus T` and `s_j=|S_j\setminus T|`.  When `s_j>=1`, the
principal compression `B[S_j\setminus T]=J_(s_j)-I_(s_j)` has largest
eigenvalue `s_j-1`; hence

```math
 s_j-1\le\|B[S_j\setminus T]\|_(2->2)
 \le\|B[C]\|_(2->2)\le R.
```

This inference is valid even though there are arbitrary edges from the
clique to the rest of `C`: a principal compression of a symmetric matrix
cannot have larger operator norm than the full matrix.  If `s_j=0`, then
all `k` block vertices were deleted and the desired deletion bound is
automatic.  In both cases

```math
 |T\cap S_j|\ge(k-1-R)_+.
```

Summing over the disjoint blocks proves (SH.11).  Thus there is no flaw in
the principal-compression inference; only the proof sentence saying the
empty block has norm `s_j-1=-1` should be avoided.

For `k=floor(Lsqrt(n))` and
`m=floor(sqrt(n)/L^3)`, the condition `L=o(n^(1/6))` makes both floors
asymptotically negligible and gives

```math
 mk=(1+o(1)){n\over L^2},
 \qquad
 m\binom k2=(1+o(1)){n^{3/2}\over2L}.
```

Consequently (SH.12) is an `o(n^(3/2))` additive halo, and for fixed `C`,

```math
 m(k-1-C\sqrt n)_+=(1-o(1)){n\over L^2},
```

which proves (SH.13).  Choosing a subpolynomial diverging `L`, such as
`log log n`, makes this lower bound `n^(1-o(1))`, larger than
`O(n^(1-delta))` for every fixed `delta>0`.  The theorem therefore really
falsifies every polynomially sublinear fixed-`O(sqrt(n))` peeling claim.
Because `n/L^2=o(n)`, it does not falsify a qualitative `o(n)` peel, exactly
as the draft states.

### Copy-ready empty-block repair for SH.2

Replace the three sentences beginning “For each `j`” by:

> For each `j`, put `s_j=|S_j setminus T|`.  If `s_j>=1`, the principal
> block on `S_j setminus T` is `J_(s_j)-I_(s_j)` and has operator norm at
> least its largest eigenvalue `s_j-1`.  Principal compression cannot
> increase operator norm, so `s_j-1<=||B[V setminus T]||_(2->2)<=R`.  If
> `s_j=0`, all `k` vertices of the block were deleted.  In either case
> `|T intersect S_j|>=(k-1-R)_+`; summing over the disjoint blocks proves
> (SH.11).

## 4. Fractional-balance consistency check: PASS

This section checks that the current draft consistently incorporates the
earlier independent audit rather than redoing that audit from scratch.

### 4.1 FB.1

The minimax order and positive-part sign in (FB.7) remain correct.  For the
maximizing `w`, every shell member satisfies
`sum_e w_ea_ez_e>=Ev(S)`.  The random perturbation has fixed-query variance
at most `4qE=O(kappa n^(3/2))` and summands bounded by two.  Paying
`log|Z_n^+|=O(n)` in Bernstein gives exactly

```math
 C(\sqrt\kappa n^{5/4}+n)=eta_n(\kappa)n^{3/2}.
```

The normalization

```math
 2qE=\kappa(1-1/n)n^{3/2}
```

is exact.  For a non-shell cut, the actual expectation bound is the
slightly stronger

```math
 Q(a)-\kappa(1+1/n)n^{3/2},
```

so the draft's (FB.11) is safe.  The hypothesis
`epsilon+eta_n(kappa)<kappa/2` implies the needed
`kappa-eta_n(kappa)>epsilon`.  If `v(S)` strictly exceeded the right side
of (FB.5), every realized energy would be strictly below
`Q(a)-epsilon n^(3/2)<=M_n`, a valid strict contradiction.  The same
minimizing measure is used for (FB.5) and (FB.6), and
`Q(a)>2kappa n^(3/2)` gives `P-N>0`, hence `N<=P` and the displayed factor
two.

The parameter consequences are consistent.  For fixed sufficiently small
`epsilon`, `kappa=sqrt(epsilon)` eventually satisfies the theorem condition
and, using the positive lower bound on `M_n/n^(3/2)`, the hypothesis of
(FB.6).  For an exact minimizer and `kappa=n^(-1/6)`,
`eta_n(kappa)=O(n^(-1/3))`, so `eta_n/kappa=O(n^(-1/6))`.

### 4.2 FB.2

For each edge the empirical mean has expected absolute deviation at most
`K^(-1/2)` by Cauchy--Schwarz and the variance bound.  Averaging those
expectations over edges proves existence of one sample satisfying (FB.13).
No simultaneous edgewise event and hence no union bound over `E` is needed.

### 4.3 FB.3

The present statement now assumes a positive `c_0`, makes `r=r_n` an
integer sequence, requires `n<=r`, and states
`r/n^(3/2)->0`.  These are the quantifier and sampling repairs requested by
the prior audit.  On the deficit-`2r` shell, the population mean is
`Omega(n^(-1/2))`; Hoeffding--Serfling therefore gives
`exp(-Omega(r/n))`.  Since `r/n>=1`, choosing the theorem's `c_1` strictly
below the tail constant makes the shell union bound strictly less than one.

Minimality is used with the correct orientation.  A maximizing signed cut
`z_*` for the flipped signing satisfies

```math
 0\ge d_a(z_*)+2\sum_(e\in F)a_e(z_*)_e.
```

The sum is at least `-r`, so `d_a(z_*)<=2r`, placing `z_*` back in the
shell, while the same inequality makes its sampled sum nonpositive.  This
contradicts the chosen `F`.

The corollary now uses `r=floor(n^(3/2)/L_n)` and explicitly allows a
constant decrease to absorb rounding.  Its assumptions
`L_n->infinity`, `L_n<=sqrt(n)` ensure both `r/n^(3/2)->0` and eventually
`r>=n`.  Finally, the interpretation has been repaired to say that shell
cardinality rules out a bounded explicit enumeration, not a bounded
generative or algebraic description.  This is consistent with what FB.3
actually proves.

## 5. Final classification

The audited statements establish three different, noninterchangeable
facts:

1. PR.1 gives an exponential one-hot response packing for unrestricted
   amplitude-`n` fields and therefore a genuine linear state-bit lower
   bound.  It does not construct a low-cap physical continuation.
2. SH.1 gives an edit-stable weighted spectral envelope but no target-scale
   response carrier; SH.2 rules out polynomially sublinear fixed-root
   peeling but leaves a qualitative `o(n)` peel open.
3. FB.1--FB.3 force fractional balance and a large near-top shell, but
   neither first marginals nor shell cardinality alone yield contextual
   response separation.

All of those scope qualifications in the current drafts are accurate.  No
near-minimizer classification is overstated after the small copy edits
identified above.
