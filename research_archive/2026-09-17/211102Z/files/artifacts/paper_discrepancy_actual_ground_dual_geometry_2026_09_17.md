# Actual finite ground codes: nonradial closure and unavoidable covariance

2026-09-17. Exact finite mechanism tests, not asymptotic structure claims.
The matrices' global-minimum labels are imported from the archive; the
caps, complete ground codes, and all identities below are independently
enumerated. The underlying 42 response-game certificates are frozen and
are not changed by this diagnostic.

## 1. Two different isotropy questions

For the complete absolute ground code C of A, write

```
beta(C)=min_nu max_(x in C) E_nu |h.x|
       =max_mu supported on C min_(h in {+-1}^n) E_mu |h.x|.
```

Column isotropy constrains the MINIMIZING law nu by E_nu hh^T=I.
Query isotropy instead constrains the MAXIMIZING law mu by E_mu xx^T=I.
These are different restrictions. The existence of an isotropic primal
does not certify an isotropic dual, or conversely.

The [full-column games](paper_discrepancy_full_column_games_2026_09_17.md)
give exact primal and dual values against every physical sign column.
An additional feasibility computation, followed by rational all-column
verification, finds isotropic optimal query laws at n=6,12,14. For n=12
and n=14 the uniform distribution on the entire projective ground code
already works; the same is true at n=6. Global fair signs can always be
added to center these laws without changing any quantity here.

The finite absence of an isotropic ground law in the n4, both n8, and
n10 examples was already established by exact separating certificates in
[the earlier isotropy audit](transfer_adversary_exact_minimizer_isotropy_finite_2026_09_06.md).
It is not a new asymptotic obstruction. The projection argument below
strengthens the n8 cases to their exact best possible covariance norm.

## 2. A concrete nonradial mechanism at n=12

The stored n12 signing has Q=18 and twenty projective ground words,
ten in each energy sector. Let mu be uniform on these twenty words,
sigma(x)=sign H_A(x), and put

```
Sigma=E_mu xx^T,       K=E_mu sigma(x)xx^T.
```

Exact integer calculations give

```
Sigma=I,       E_mu sigma=0,       K^2=I,
Kx=sigma(x)x for every ground word x.                 (1)
```

Thus P_+=(I+K)/2 and P_-=(I-K)/2 are orthogonal rank-six projections.
The ten words in each sector form a tight frame in the corresponding
six-dimensional subspace:

```
sum_(sigma(x)=+1) xx^T=20 P_+,
sum_(sigma(x)=-1) xx^T=20 P_-.
```

Every opposite-sector overlap is zero. For each fixed ground word x,
the multiset of absolute overlaps with all twenty words is

```
13 copies of 0,   6 copies of 4,   1 copy of 12.
```

Consequently the same uniform ground law, used as a physical column
law, is isotropic and has response exactly (6*4+12)/20=9/5 at every
ground query. Conversely exhaustive exact evaluation gives

```
E_mu |h.x| >= 9/5 for every one of the 2^11
projective physical columns h.                       (2)
```

Equations (1)--(2) therefore supply a particularly transparent matching
primal/dual certificate: beta=beta_iso=9/5.

Crucially K is NOT proportional to A. Its off-diagonal entries include
0, +-1/5 and +-2/5. The best balanced law constrained to the two maximal
radial covariance endpoints instead has response exactly 9/4, by the
separate [radial certificates](paper_discrepancy_radial_correlation_2026_09_17.md).
That is a genuine 25 percent cost of the radial subclass on this actual
finite example. The gain is explained by a nonradial separation of the
two signed ground sectors, not by replacing a Gaussian variance bound
with an assumed absolute-moment formula.

The general elementary implication used here is worth recording:
if a law on signed vectors has Sigma=I and K^2=I, then its positive
and negative covariance measures are (I+K)/2 and (I-K)/2. Positivity
forces every supported vector into its respective eigenspace; the
opposite sectors are orthogonal. This implication does not assume a
special matrix A. Whether comparable closure can be found for large
actual minimizers is unresolved.

## 3. Different exact closure at n=14

The stored n14 matrix obeys A^2=13I and has Q=21. Its complete
projective ground code has 156 words. Uniform mu satisfies

```
Sigma=I,       K=(3/13)A,       K^2=(9/13)I.
```

Every row of its absolute Gram matrix has the histogram

| Absolute overlap | 0 | 2 | 4 | 6 | 10 | 14 |
|---|---:|---:|---:|---:|---:|---:|
| Multiplicity | 57 | 49 | 21 | 21 | 7 | 1 |

The row sum is 392, so its mean is 98/39. Exact enumeration against
all 8192 projective physical columns shows that 98/39 is the minimum
and that precisely the 156 ground columns attain it. The same uniform
law is therefore an isotropic primal and unrestricted dual. Here the
radial subclass ALSO attains the exact optimum, unlike n12.

At n6, uniformity similarly gives Sigma=I and a constant own-code
response 5/3, which is the all-column minimum. Its absolute-overlap
row histogram is 0 five times, 2 five times, 4 once, and 6 once.
These three examples demonstrate finite self-dual closure, not an
asserted common asymptotic design or an implication of low entropy.

## 4. Sharp covariance obstruction on both n8 ground codes

Let Sigma_0 be the covariance of the uniform ground law. In the first
n8 class the exact matrix P=(Sigma_0-I/2)/2 is a rank-two orthogonal
projection, and EVERY ground word obeys x^T P x=5. Hence for ANY law
mu supported on the full ground code,

```
||E_mu xx^T||_op >= tr(P E_mu xx^T)/rank(P)=5/2.       (3)
```

Uniform mu attains equality because Sigma_0=I/2+2P. It is also an
optimal query law for the unrestricted response game, with value 3/2.
Its own-code response is 5/2, however: the minimizing physical columns
are all off the ground code. Thus an optimal query law need not itself
be an optimal column law.

In the second n8 class, P=Sigma_0-I/2 is a rank-four orthogonal
projection and every ground word satisfies x^T P x=6. The same argument
gives the sharp lower bound 3/2 on every supported covariance norm,
again attained by uniform mu. This law has own-code and all-column
minimum response 3/2, so it is now also an unrestricted primal.
The two ground codes have the same size and cap but distinct geometry.

These bounds exclude isotropic query laws without any numerical SDP.
They do NOT disprove an O(1) covariance bound for response-optimal
query laws of asymptotically large EXACT minimizers.

For asymptotically near-minimizing signings, a different obstruction is
already rigorous: the archived planted-clique construction forces
EVERY near-ground-supported law to have covariance norm of order
N^(2/3), while the cap exceeds M_N by only O(N^(4/3)). The exact
scope and archive collision are recorded in
[the covariance-spike audit](paper_discrepancy_nearminimizer_covariance_spike_2026_09_17.md).
That construction is not an exact minimizer, and its coherent spike
can be deleted at negligible normalized response cost. It therefore
does not settle the structural question left open above.

## 5. Reproducibility and limitations

Run the solver-free exact verification:

```
.venv/bin/python computations/paper_discrepancy_2026_09_17_ground_covariance_geometry.py
```

It enumerates all physical words, recomputes the caps, checks the
projection identities with integers/rationals, verifies the signed
sector identities, and evaluates every physical column. No floating
eigenvalue is used for any claim in Sections 2--4.

The separate scripts `paper_discrepancy_2026_09_17_column_dual_moments.py`
and `paper_discrepancy_2026_09_17_optimal_dual_isotropy.py` inspect the
frozen response duals and construct exact positive feasibility witnesses.
Their floating eigenvalue lists and infeasibility solver statuses remain
diagnostics only. The latter are not used in place of the exact archived
separators or the explicit projection proof (3).

No statement here gives a uniform subGaussian column law, an actual
macroscopic extension with favorable value, control of all energy shells,
or a dimension-free covariance theorem for actual minimizers. The n12
ground response already deteriorates from 9/5 to 1323/536 under the
isotropic column constraint at deficit two. The entire nearcode, not
just these exact ground frames, remains essential to the original task.
