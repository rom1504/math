# Adversarial verification of the syndrome-rooted coding report

**Scope.**  This report independently audits
`phase2_coding_rooted.md` and
`verify_phase2_code_syndrome_profiles.py`.  It does not promote any result to
the surface theory files.

**Verdict.**  The mathematical cores of Theorems CR.1 and CR.2, the strict
quotient pair, and the outer-spectrum falsifier are correct.  The response
quotient is noncircular: the exposing environments are fixed from the labeled
syndrome interface and do not depend on the unknown fragment.  Four scope
clarifications should be made before promotion:

1. “Coarsest” is relative **only** to exact covering-radius queries obtained
   by appending arbitrary parity-check fragments over the same labeled
   syndrome group.  It is false for named-root, puncturing, weighted, or
   finite-temperature queries.
2. The present outer-spectrum witness uses a rank-deficient environment.  It
   is valid under the model as defined, but an equally small full-rank repair
   is available and removes an ambiguity in the phrase “full-rank
   fragments.”
3. The `Theta(2^w)` bit lower bound is a worst-case, exact/sub-half-unit
   statement.  Its fixed-length packing uses length `Theta(2^w)`; it does not
   establish the same rate for fragments of length polynomial in `w` or for
   distortion growing with `w`.
4. The state is strictly smaller than the full **labeled root-distance table
   or code**, not than its declared complete future-response vector.  Theorem
   CR.1 proves that the latter and the syndrome profile are equivalent.

With those qualifications, the report is suitable as a positive known-model
validation of a query-relative feature algebra.

## 1. Independent reconstruction of CR.1

Let `G=F_2^w`, let `H` be a `w` by `n` matrix of row rank `w`, and assume its
columns are nonzero.  Write

```math
C_H=\ker H,
\qquad
\lambda_H(s)=\min\{|e|:He=s\}.
```

### Root distance and radius

For every root `x`,

```math
x+e\in C_H
\quad\Longleftrightarrow\quad
H(x+e)=0
\quad\Longleftrightarrow\quad
He=Hx.
```

Taking the minimum weight proves

```math
d(x,C_H)=\lambda_H(Hx).
```

Full row rank is used here to make `H` surjective.  Therefore every
`s in G` equals `Hx` for some root and

```math
\rho(C_H)=\max_{s\in G}\lambda_H(s).
```

If `H` is not full row rank, the first formula remains true but the maximum
must be restricted to `im H`; the profile is infinite off that image.  Thus
the full-rank hypothesis is substantive exactly at this point.

### Multiplicities disappear

Suppose a correction of minimum weight uses two coordinates whose columns
are the same vector `g`.  Turning off both coordinates preserves its syndrome
because `g+g=0` and lowers its weight by two.  Hence no minimum correction
uses a column type twice.  It follows that `lambda_H` depends only on

```math
S_H=\{h_i:1\le i\le n\}\subseteq G\setminus\{0\}.
```

Conversely, for nonzero `s`,

```math
\lambda_H(s)=1
\quad\Longleftrightarrow\quad
s\in S_H.
```

Thus `S_H` and `lambda_H` determine each other.  Zero columns are irrelevant
to these radius queries: their coordinates can always be matched by a
codeword and may indeed be deleted.  They would cease to be irrelevant for
queries that inspect the code itself.

### Concatenation and the min-plus law

For fragments `H_1,H_2`, split a correction to syndrome `s` as `(e_1,e_2)`
and set `u=H_1e_1`.  Its second syndrome is necessarily `s+u`.  Minimizing
first inside each block and then over `u` gives

```math
\lambda_{[H_1\ H_2]}(s)
=\min_{u\in G}
 \{\lambda_{H_1}(u)+\lambda_{H_2}(s+u)\}.
```

Both inequalities are attained: every combined correction supplies a
candidate on the right, and minimizers for any candidate `u` concatenate to
a correction on the left.  This is therefore an equality, not a relaxation.
The support version is simply

```math
S_{[H_1\ H_2]}=S_{H_1}\cup S_{H_2}.
```

Associativity follows either from concatenation or from associativity of
group min-convolution.  For finite functions, comparing every candidate in
the two minima gives

```math
\|f\star g-f'\star g'\|_\infty
\le \|f-f'\|_\infty+\|g-g'\|_\infty.
```

There is no min/max reversal here: min-plus computes coset-leader distances,
and covering radius subsequently takes the outer maximum over syndromes.

### Exact future-environment sufficiency

For a future fragment `E`, the same calculation gives

```math
\mathcal R_H(E)
=\rho(\ker[H\ E])
=\max_{s\in G}(\lambda_H\star\lambda_E)(s).
```

Thus `lambda_H`, together with the presented query `E`, answers every query
in the declared class.  The decoder does not require `M_n`, an optimizer at
the target length, or any hidden information about `H`.

### Exact future-environment minimality

For each nonzero `s`, let `E_s` contain one copy of every nonzero syndrome
except `s`.  This environment depends only on `(G,s)`, not on `H`.  It is
full row rank for every `w>=2`.  If `s in S_H`, the composite support is all
of `G\setminus\{0\}` and its radius is one.  If `s notin S_H`, its support is
exactly all nonzero elements except `s`.  Choose
`u notin {0,s}`; then both `u` and `s+u` occur and

```math
s=u+(s+u),
```

so the missing syndrome has distance exactly two.  Hence

```math
\mathcal R_H(E_s)=1_{\{s\in S_H\}}+2_{\{s\notin S_H\}}.
```

The complete response vector therefore recovers every support bit.  Two
fragments have the same responses to all future fragments if and only if
they have the same `S_H`, equivalently the same `lambda_H`.  This proves the
claimed coarsest deterministic quotient, up to injective recoding.

The quantifiers are important.  The environments `E_s` have `2^w-2`
columns.  The theorem establishes minimality for unrestricted future
fragments over a fixed labeled `G`; it does not establish minimality for
environments whose length is bounded by `poly(w)`.  This is a scope boundary,
not circularity.

For `w=1`, every full-rank fragment has support `{1}` and the quotient is
trivial.  Excluding this case from the exposing argument is correct.

## 2. Independent reconstruction of CR.2

Fix a basis `B` of `G` and put

```math
N=2^w-1-w,
\qquad t=\lfloor N/2\rfloor.
```

For every `t`-subset `U` of the nonzero nonbasis vectors, take a fragment
with one column of each type in `B union U`.  These fragments

- all have the same length `w+t`;
- all have full row rank because they contain `B`; and
- have distinct response vectors, since an `E_s` query exposes any support
  bit on which two sets differ.

The response separation is at least one in sup norm.  In fact it is exactly
one for this middle-layer family: its supports have at least half the group
and hence diameter at most two, while all nontrivial composite radii are one
or two.  Only the lower separation is needed.

If a deterministic message were shared by two fragments and a common
decoder approximated every response with error strictly below `1/2`, then at
the separating query the triangle inequality would give a true distance
strictly below one, a contradiction.  Consequently the range has at least
`binom(N,t)` elements.  Since

```math
\binom N{\lfloor N/2\rfloor}\ge {2^N\over N+1},
```

the message length is at least

```math
N-\log_2(N+1)
=2^w-1-w-\log_2(2^w-w),
```

which is exactly (CR.13).  The support indicator uses `2^w-1` bits.  Thus
the worst-case exact response complexity is `Theta(2^w)` bits, with an
`O(w)` gap between the displayed bounds.

This is a genuine bit-count theorem, not a count of tropical factors.  Its
scope must remain explicit:

- the packed fragment length `w+t` is `Theta(2^w)`;
- the response error is unnormalized and less than half the integer lattice
  spacing; and
- no lower bound at additive error `epsilon w` follows from it.

### Shannon corollary normalization

Let the `N` optional support bits be independent and uniform, keep the basis
fixed, and use the `N` special environments with uniform query measure.  The
response at coordinate `s` is one or two according to the corresponding
bit.  Therefore

```math
\|R_U-R_V\|_{L^2(\mathrm{Unif}[N])}^2
={d_H(U,V)\over N},
```

so the inverse-Hamming modulus is `kappa=1/N`.  If `Delta` denotes **mean**
squared error over those environments, the general posterior-width theorem
indeed gives

```math
I(U;Z)\ge N[1-g(\min\{4\Delta,1\})].
```

Thus (CR.16) has the correct constants.  Calling the raw response vector an
“isometric copy” and then using averaged error is harmless but potentially
confusing; the report should state the normalized identity above.

## 3. Finite witnesses

### Strict quotient pair

For

```math
H_A=(1,1,1,2,3),
\qquad H_B=(1,1,2,2,3)
```

over `F_2^2`, both supports contain all three nonzero syndromes, hence both
profiles are `(0,1,1,1)` and all future responses agree.  Direct kernel
enumeration gives

```math
W_{C_{H_A}}(z)=1+3z^2+3z^3+z^5,
```

```math
W_{C_{H_B}}(z)=1+2z^2+4z^3+z^4.
```

Coordinate isometries preserve weight enumerators, so the codes are not
isometric.  This proves that the response quotient does not reconstruct the
code.  It also does not reconstruct the labeled root-distance table, whose
zero set is the code itself.

### Same outer spectrum, different future response

For the two rank-three fragments

```math
H_A=(1,2,3,4),
\qquad H_B=(1,2,4,7),
```

the independently recomputed profiles are

```math
(0,1,1,1,1,2,2,2),
\qquad
(0,1,1,2,1,2,2,1).
```

Both histograms are `(1,4,3)`.  Since both kernels have size two,

```math
O_{C_{H_A}}(z)=O_{C_{H_B}}(z)=2+8z+6z^2.
```

The draft appends `E=(3,5,6)`.  Its columns span only the even-weight plane,
so `rank(E)=2`; nevertheless `[H_A E]` and `[H_B E]` remain full rank because
the left fragments are full rank, and the stated radii `(2,1)` are correct.

If “future fragment” is intended to be full rank on its own, replace it by

```math
E'=(1,3,5,6).
```

This support spans `F_2^3`.  The first union is
`{1,2,3,4,5,6}` and has radius two, while the second union is all seven
nonzero types and has radius one.  Thus the counterexample survives the
stronger convention with no conceptual change.

The environment was selected to distinguish the displayed pair, but it is
the same environment for both fragments and is a legitimate member of the
predeclared universal query class.  This is an ordinary separating query,
not an environment that receives hidden access to the unknown input.

### Reproducibility check

Running

```text
.venv/bin/python \
  extremal_information/experiments/verify_phase2_code_syndrome_profiles.py
```

reproduced the saved JSON byte-for-byte.  The script exhaustively checks
the four spanning supports over `F_2^2`, the 92 spanning supports over
`F_2^3`, all ordered pair compositions, and every special exposing
environment.  An independent direct check of the full-rank replacement
`E'=(1,3,5,6)` returned radii `(2,1)`.

The script deliberately deduplicates columns when computing syndrome
profiles and retains duplicates when computing code weight enumerators.  That
is exactly the distinction required by the proofs rather than an accidental
loss of data.

## 4. Comparison with the Sheshadri door identity

Sheshadri's theorem concerns a fixed code `C` and a fixed coordinate cut
`L sqcup R`.  It writes the complete conditional distance table as

```math
W(x_L,x_R)
=\min_{\tau\in P_R(C)/C_R}
 \{D(x_L,\tau)+d(x_R,\tau)\}
```

and proves that the `2^s` door states are also necessary among arbitrary
min-plus factorizations, where

```math
s=\dim C-\dim C_L-\dim C_R.
```

For `C=ker[H\ E]` with `H` full row rank and `rank(E)=r`, a direct dimension
count gives

```math
\dim C=n+m-w,
\quad \dim C_L=n-w,
\quad \dim C_R=m-r,
\quad s=r.
```

In particular, a full-rank future fragment has `2^w` door states.  Reindexing
the door state by its syndrome gives the same min-plus decomposition behind
the convolution formula.  Therefore the following parts of the coding draft
are classical or direct specializations, not new consequences of extremal
information theory:

- syndrome/coset-leader decoding;
- the min-plus split over a syndrome boundary;
- associative concatenation; and
- the use of `2^w` syndrome states in the full-rank worst case.

What the coding report genuinely adds relative to the cited door theorem is
orthogonal to its tropical-rank result:

1. It varies the left fragment and declares the complete family of appended-
   fragment **covering-radius** queries, then proves that response
   equivalence is exactly equality of the labeled syndrome profile.
2. Binary Hamming cancellation collapses that profile further to the support
   set of distinct column types, so multiplicities and the code itself are
   forgotten under repeated composition.
3. The special environments give a deterministic message-bit lower bound
   for this response experiment.  Tropical factor count alone does not imply
   a bit or mutual-information lower bound.
4. The explicit outer-spectrum collision identifies the missing datum as the
   alignment of profile levels with the additive syndrome group.

Conversely, CR.1--CR.2 do not improve Sheshadri's exact rank theorem.  They do
not prove that every min-plus representation of a fixed conditional root
table needs `2^w` terms; that is precisely the independent content of the
door lower bound.  The two results should remain presented as complementary.

## 5. Final classification

| Claim | Audit result | Required qualification |
|---|---|---|
| Root distance and radius identities | Verified | Full row rank is needed for the all-`G` maximum. |
| Support/profile equivalence | Verified | Binary unweighted Hamming setting only. |
| Min-plus convolution and union algebra | Verified | Common labeled syndrome group is essential. |
| Exact future-query sufficiency | Verified | Appended-fragment covering-radius queries only. |
| Coarsest exact quotient | Verified | Unrestricted future fragments; `w>=2` exposing proof. |
| `Theta(2^w)` response bits | Verified | Worst case, length `Theta(2^w)`, error `<1/2`. |
| Shannon response bound | Verified | `Delta` must mean mean squared error over uniform special queries. |
| Strict quotient pair | Verified | Shows forgetting of code/root table, not response vector. |
| Outer-spectrum falsifier | Verified | Current `E` has rank two; use `E'` for full-rank convention. |
| Novelty beyond door identity | Partially new | Quotient/minimality, bit packing, and outer falsifier are additions; convolution itself is classical. |

The strongest next question in the draft is correctly identified as an
approximate one: determine response-metric entropy at distortion proportional
to `w`.  Neither CR.2 nor the exact tropical-rank theorem controls that scale.
