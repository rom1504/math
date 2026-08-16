# Second audit: contextual responses, Rees collapse, and closed summaries

**Status.** Independent verification.  This audit covers CRL.1--CRL.3,
HRC.1--HRC.5, CSC.1--CSC.2, and the prime-cycle experiment.  All central
claims survive.  Two bookkeeping corrections and one useful strengthening
are recorded below.

## 1. Contextual response law

### CRL.1: verified

For a commutative monoid,

```math
d_F(a*u,b*u)
=\sup_c|F(a*u*c)-F(b*u*c)|
\le d_F(a,b).
```

Applying this once to each factor and using the triangle inequality gives
CRL.3.  The zero-distance relation is therefore a congruence.  Conversely,
two objects assigned the same exact closed state must have identical answers
after every continuation, so the syntactic quotient is coarsest.

The context-ideal variant is also correct provided, as stated, the context
set is closed under multiplication by arbitrary monoid elements.

### CRL.2: verified

The power-set construction is antitone in every possible enlargement:
values can move from one to an entry of `R`, or from an entry of `R` to
zero, but never upward.  Since the range is contained in `[0,1]`, distinct
sets have response difference at most their symmetric-difference distance.
Thus arbitrary finite response tables really do occur inside monotone
idempotent union systems.

### CRL.3: verified, including the `L^2` constants

For `a_j=lambda_j(z_j-z'_j)` and independent
`B_j=1_{j in P}`,

```math
E\left(\sum_j a_jB_j\right)^2
={1\over4}\sum_j a_j^2
 +{1\over4}\left(\sum_j a_j\right)^2.
```

This is exactly CRL.10.  With equal weights, division by Hamming distance
gives at least `lambda^2/4`; two oppositely directed changes attain equality
when `q>=2`.  Substitution into Theorem 7.1 gives
`4 Delta/(kappa q)=16 Delta/(lambda^2 q)`, so CRL.12 has the correct
constant.  The syndrome and matching-Max-Cut mappings are exact.

## 2. Hard-core/Rees law

### HRC.1: verified, with a needed saturation sentence

The easy set `I_eta` is an ideal by translation contraction.  The proposed
Rees quotient is well defined because `I_eta` is also saturated under exact
response equivalence:

```math
d_F(q,q')=0
\quad\Longrightarrow\quad
d_F(q,0)=d_F(q',0).
```

Thus if one choice of representatives enters the ideal, every equivalent
choice does.  This sentence should be added to the proof before promotion.
Terminal decoding and the factor-two converse are then immediate and sharp.

### HRC.2: verified

For a nonnegative antitone deficit with absorbing zero, every continuation
lies in `[0,F(q)]` and the identity context attains `F(q)`.  Hence
`d_F(q,0)=F(q)`.  The ideal `{F<=2 eta}`, decoded by its midpoint `eta`,
has error at most `eta` after every continuation.  Any cell containing zero
must lie inside it, because two response profiles sharing a decoder are at
distance at most `2 eta`.  This proves genuine maximality, not merely a
sufficient collapse.

### HRC.3--HRC.4: verified

For matroid flats `X,Y`, the contexts `Y,X` attain the two directed join
increments.  For an arbitrary context `T`, monotonicity and diminishing
returns give

```math
r(X\cup T)-r(Y\cup T)
\le r(X\vee Y)-r(Y)
```

whenever the left side is nonnegative, and symmetrically otherwise.  This
proves HRC.10 and makes equality of closure exactly equivalent to equality of
all future responses.  Truncating at codimension `k` is precisely HRC.2
with target error `k/2`.

### HRC.5: valid as stated and stronger up to `epsilon<1/4`

The printed `epsilon<1/8` packing is correct.  More generally, for
`w=2d`, choose

```math
2\epsilon<\gamma<1/2
```

and pack middle-dimensional subspaces at injection distance at least
`gamma*w`.  A ball has exponent at most
`(gamma-gamma^2+o(1))w^2`, versus total exponent `w^2/4+o(w^2)`;
the greedy packing therefore has logarithmic size

```math
((1/2-\gamma)^2-o(1))w^2.
```

This yields `Theta(w^2)` response complexity for every fixed
`epsilon<1/4`.  In odd width, take a fixed line complementary to a fixed
even-dimensional hyperplane and adjoin that line to every packed subspace;
injection distances are unchanged.

### Bookkeeping corrections

The current HRC draft uses the tag `HRC.8` twice: once for the hard-core
state count and once for the matroid response definition.  The sentence
“Equation (HRC.7) is a criterion” should refer to the hard-core count, not to
the midpoint ideal `J_eta`.  Renumbering is required before surface
promotion.

## 3. Minimum-index closed summaries

CSC.1 is exact under its declared scalar deterministic model.  A closed
summary map is a monoid homomorphism after restricting the target to its
reachable image, and its kernel is a congruence.  One decoder value within
`epsilon` of every value in a class forces class oscillation at most
`2epsilon`.  Conversely, a real interval of diameter `2epsilon` has a
midpoint of radius `epsilon`.  Therefore

```math
K_{\rm closed}(\epsilon)
=\min_\theta
\{|M/\theta|:\operatorname{osc}_F(C)\le2\epsilon
\text{ for every class }C\}
```

is not only a lower bound but the exact minimum.

This relies on scalar real output.  For a general normed vector output,
diameter `2epsilon` need not imply an `epsilon`-radius common center.

## 4. Prime-cycle separation and experiment

For `M=Z_p` and `F(x)=cos(2*pi*x/p)`, every congruence is induced by a
subgroup.  Prime order leaves only equality and the universal congruence.
The universal class has optimal one-state error

```math
e_1(p)={1+\cos(\pi/p)\over2},
```

so every fixed `epsilon<1` forces `p` exact closed states for sufficiently
large primes.

In contrast, translated response profiles lie on the phase circle.  With
`k=ceil(2*pi/epsilon)`, choosing the nearest residue to each regular
`k`-mesh point gives angular covering error at most

```math
{\pi\over k}+{\pi\over p}<\epsilon
```

when `p>k`; if `p<=k`, use every residue.  Cosine is one-Lipschitz, so
the actual all-context response net has `O(1/epsilon)` states.

The script
`experiments/verify_phase3_prime_cycle_summary.py` correctly evaluates the
full context supremum, checks the mesh error for `p=5,11,31,101`, and
matches the one-state formula above.  Its output
`experiments/phase3_prime_cycle_summary_results.json` passes.  The script's
center-count assertion has a harmless extra `+1`; the proof gives at most
`k` distinct centers.

The apparent two-dimensional compression

```math
\cos(2\pi(x+c)/p)
=\cos(2\pi x/p)\cos(2\pi c/p)
 -\sin(2\pi x/p)\sin(2\pi c/p)
```

does not contradict the lower bound.  Exact complex multiplication retains
all `p` distinct phases (and hence `log_2 p` bits); quantizing the phase
gives a small response net but destroys homomorphic closure.  This is the
precise distinction between a low-dimensional response roof, a small
one-shot net, and a small exact feature algebra.

## 5. Judgment

The audits support one general law and one hard boundary:

1. exact reusable scalar feature algebras are exactly low-index congruences
   with classwise `F`-oscillation at most `2epsilon`;
2. one-shot future-response metric entropy can be bounded independently and
   can be asymptotically much smaller.

The prime cycle makes the separation maximal.  Thus the congruence law is
mathematically exact but classical in mechanism; its generative value is to
prevent one-shot nets from being misclassified as reusable composition
states.
