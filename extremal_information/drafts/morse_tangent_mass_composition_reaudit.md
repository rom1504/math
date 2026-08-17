# Re-audit addendum: repaired Morse tangent-mass composition

Scope: second independent audit of
`drafts/morse_tangent_mass_composition.md` and
`experiments/verify_morse_tangent_mass.py`, using the first audit as a
checklist.  This addendum focuses on the newly added localization lemma,
Gaussian finite-parameter semigroup, determinant/amplitude algebra,
uniformity, integer recovery, numbering, and the precise scope of the
finite-dimensional claim.  I made no edits to the audited draft or verifier.

## Verdict

**REPAIR (local, not structural).**  The repaired draft now contains a
genuine strict finite-dimensional composable class: the full-lattice
Gaussian family in TM.5.  Its parameter law, tangent amplitude, determinant
identity, associativity, and fixed-compact-set uniformity are correct.  The
localization lemma also repairs the earlier binomial/multinomial boundary-
summand gap.

Two points should be fixed before promotion.

1. The proof of TM.6 says that the total flooring error is “at most
   polynomially many units.”  That is false for a convolution: replacing
   `X_k,Y_k` by their floors can lose as much as `X_k+Y_k` in one product.
   The corollary is nevertheless true, by splitting into a uniform saddle
   neighbourhood (where **both** factors are exponentially large and hence
   flooring is relatively negligible) and its exponentially suppressed
   complement.  The proof needs this argument.
2. The last paragraph calls TM.6 “the minimal compositional repair.”  No
   minimality theorem is proved.  “A finite-dimensional compositional repair
   on a natural model class” is supported; “minimal” is not.

There are also two editorial ambiguities worth repairing: “boundary types”
in TM.4 should read “boundary summands” unless lower-dimensional boundary
target types are separately intended, and theorem labels and equation tags
currently reuse `TM.1`--`TM.6` while the equation tags jump from 17 back to
9.  The theorem numbering itself is sequential and nonduplicated.

## 1. Lemma TM.3: localization

The localization mechanism is correct.  Write the fixed saddle
neighbourhood as `U`.  On `U`, the local two-sided estimates give the same
`Theta(n^(alpha+beta+d/2)e^(nh))` contribution as TM.1, and the uniform local
relative expansions give TM.2.  If the global envelopes carry at most
`n^D` for a fixed `D`, then the at most `O(n^d)` terms outside `U` contribute

```math
O(n^(D+d)e^(n(h-eta)))=O(e^(n(h-eta/2))),
```

for all sufficiently large `n`.  This is exponentially negligible relative
to every fixed polynomial multiple of `e^(nh)`.  The uniform version is valid
when `D`, the envelope constants, the gap, the saddle neighbourhood, and the
local remainders are uniform.

For canonical wording, “a global upper bound of the form polynomial times”
should explicitly mean a polynomial of uniformly bounded degree with a
uniform coefficient.  The present prose conveys that intent in the last
sentence but does not put it into the numbered hypotheses.

This lemma is exactly what was missing from the first draft.  Uniform
Stirling estimates are needed only near an interior entropy saddle; the
method-of-types envelope and strict concavity dispose of summands near the
simplex boundary.  It does **not** imply the same `d/2` formula when the
*target type itself* lies on a simplex face, since then tangent dimension can
drop.  Thus the sentence “Boundary types are handled” is correct only if it
means boundary summands in an interior-target convolution.

## 2. TM.5: Gaussian semigroup algebra

Let `P,Q` be the two precision matrices and put

```math
Sigma=P^(-1), \qquad Tau=Q^(-1), \qquad
R=(Sigma+Tau)^(-1).
```

Completing the square in the convolution gives output centre `mu+nu`,
precision `R`, leading constant `c+c'`, and saddle precision `P+Q`.  The
standard shifted Gaussian lattice sum is

```math
(1+o(1))n^(d/2)(2pi)^(d/2)/sqrt(det(P+Q)),
```

uniformly in the lattice shift when the eigenvalues remain in a fixed compact
subset of `(0,infinity)`.  Therefore

```math
alpha_out=alpha+alpha'+d/2,
qquad
a_out={(2pi)^(d/2)aa'\over sqrt(det(P+Q))}.
```

The determinant identity used in the draft is valid without assuming that
`P` and `Q` commute:

```math
P^(-1)(P+Q)Q^(-1)=P^(-1)+Q^(-1),
```

so

```math
det(P+Q)=det(P)det(Q)det(P^(-1)+Q^(-1)).
```

After passing to total-mass amplitude

```math
m=a(2pi)^(d/2)/sqrt(det P),
```

one obtains exactly `m_out=mm'`.  Hence

```math
(c,mu,Sigma,alpha,m) star
(c',mu',Sigma',alpha',m')
=(c+c',mu+mu',Sigma+Sigma',alpha+alpha'+d/2,mm').
```

This operation is genuinely associative and commutative coordinatewise.  It
also agrees with the amplitude produced by either bracketing of a
three-factor convolution.  There is no omitted `2pi`, determinant, or power
of `n`.

The parameter count

```math
d(d+1)/2+d+3
```

is correct.  On any compact coordinate set its covering number satisfies
`log N(epsilon)=O_d(log(1/epsilon))`, with the implicit constant depending on
the compact set and the chosen coordinate metric.

### Uniformity and depth

The one-step relative error is uniform for compact parameter families and
compact `ell/n` query sets.  In fact the shifted full-lattice Gaussian sum
admits a shift-uniform estimate, so every **fixed finite** iterated
convolution is represented by the same parameter law.  The statement does
not establish an error bound uniform in a composition depth that grows with
`n`; the draft does not claim such a bound.  Likewise a bounded initial
parameter set is not invariant under arbitrary depth because `c`, `alpha`,
and `Sigma` add.  Neither fact invalidates the finite-dimensional semigroup
claim.

The phrase “strict finite-dimensional composable response state” is now
justified **for this Gaussian family**, understood as an asymptotic response
state under convolution.  It is not a minimal quotient, and it supplies no
closure theorem for arbitrary Morse fields.  The draft otherwise draws this
scope boundary correctly in Sections 3 and 4.

## 3. TM.6: finite integer recovery

The proposed truncation-and-floor construction can be made rigorous on the
declared compact parameter and query sets.  Compactness lets one choose `R`
so every convolution saddle lies a fixed scaled distance inside both
truncation boxes.  One can then choose `C` and `delta>0` so that, throughout a
fixed common saddle neighbourhood,

```math
X_k=e^(nC)G_n^theta(k)>=e^(delta n),
\qquad
Y_k=e^(nC)G_n^phi(ell-k)>=e^(delta n).
```

There,

```math
0<=X_kY_k-floor(X_k)floor(Y_k)
 <=X_k+Y_k,
```

and division by `X_kY_k` bounds the relative loss by
`2e^(-delta n)`.  Outside that neighbourhood, the Gaussian quadratic gap
shows that the **entire unrounded product sum** is exponentially smaller than
the main saddle mass; the rounded product is bounded above by that same sum.
The parts cut off by the two boxes are treated identically.  Consequently,
uniformly on the compact query set,

```math
(widehat G_n^theta*widehat G_n^phi)(ell)
=(1+o(1))G_n^((theta+C) star (phi+C))(ell),
```

where `(theta+C)` means replacing `c` by `c+C`.

This proves the intended all-order finite-integer realization.  The current
one-sentence proof instead calls the rounding error “polynomially many
units,” which is not a valid estimate for products and should be replaced by
the saddle/complement split above.  It would also help to display the last
formula explicitly, since “the same uniform convolution law” leaves its
normalization implicit.

The recovery is abstract: supports have polynomially many sites and weights
may have linearly many bits in `n`.  It is not a realization in a code, graph,
or sign-matrix class, and the draft correctly says so.  No state-minimality
conclusion follows.

## 4. Verifier

The repaired verifier runs successfully and reports 20 checks.  Relative to
the first version, it now usefully checks:

- an off-centre two-dimensional Gaussian convolution;
- the full matrix determinant/amplitude prediction;
- the `n^(d/2)=n` factor for `d=2`; and
- the covolume-two correction on the sublattice `2Z`.

The reported Gaussian ratios agree with one to floating-point precision, as
expected for the rapidly convergent theta-sum correction.  The Vandermonde
and quartic checks remain correct finite smoke tests.

The script does not test TM.6 flooring/truncation, three-factor associativity,
parameter-uniform suprema, or a multinomial example.  This is not a defect in
the proof, but it means the script should not be cited as independent
verification of those claims.  A small flooring test would be useful after
the proof wording is repaired.

## 5. Numbering, scope, and required final edits

The logical theorem sequence is clean:

```text
TM.1 theorem -> TM.2 theorem -> TM.3 lemma -> TM.4 corollary
             -> TM.5 theorem -> TM.6 corollary.
```

However, displayed equations also use `(TM.1)` through `(TM.6)`, so phrases
such as “TM.6” can denote either the two-speed valuation equation or the
integer-recovery corollary.  Equation tags additionally appear in the order
1--8, 13--17, 9--12, 18--25.  Use a separate equation namespace or renumber
them monotonically before canonical integration.

The remaining mathematical scope is sound after these edits:

- TM.1--TM.2 are classical local discrete-Laplace results;
- TM.3 supplies the required localization;
- TM.5 is a genuine finite-dimensional Gaussian response semigroup;
- TM.6 supplies abstract finite-integer recovery once its rounding proof is
  repaired;
- no arbitrary-Morse finite closure is claimed; and
- no minimality theorem is available.

Accordingly the status line is defensible, but the final literature sentence
must replace “the minimal compositional repair” by “a finite-dimensional
compositional repair.”

## 6. Audit classification

**REPAIR**, for the TM.6 flooring proof and unsupported last use of
“minimal.”  These are local repairs.  The core theorem TM.5 passes, and no
counterexample or structural obstruction was found to the repaired program.

## 7. Final recheck after repair

**PASS mathematically; two editorial substitutions remain.**  The revised
TM-A--TM-F sequence implements the substantive repairs above.  In particular,
TM-F now states the recovered convolution explicitly and proves flooring by
the correct two-region argument: both factors are uniformly exponentially
large near the saddle, while the complete unrounded product mass is
exponentially suppressed off the saddle and outside the truncation boxes.
The two `c -> c+C` shifts give exactly the normalization on its right-hand
side.  TM-E's determinant, amplitude, covariance-addition law,
associativity, compact-family uniformity, and finite parameter count remain
correct.  The draft no longer claims finite closure for arbitrary Morse
profiles or minimality of the Gaussian state.

Before canonical promotion, make these two notation-only edits:

1. define `theta+C` in TM-F to mean “replace only the `c` coordinate of
   `theta` by `c+C`”; and
2. replace the final literature cross-reference `(TM.6)` by
   `Corollary TM-F`, because `(TM.6)` now denotes the two-speed valuation
   equation.

Neither point changes a theorem.  Subject to those two substitutions, the
final audit verdict is **PASS**.
