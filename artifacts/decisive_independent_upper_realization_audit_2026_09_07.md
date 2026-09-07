# Independent audit: full-spin weave to the all-order .494515125 bound

Date: 2026-09-07. Status: **PASS at the scope stated here.**

I reconstructed the orbital/type realization and full-spin construction
from their formulas, rather than using earlier positive audit verdicts.
The input analytic statement is that, for p=31/32 and t=4, every strict
margin below

```math
a=-[p\log2+g_4(1)+4(1-\sqrt p)]>0
```

is attained by some finite-depth upper certificate f_r=B^r Phi_4.
The direct-E stopping source proves precisely this from BE<=E and the
certified Gaussian phase E_4(nu_p)=g_4(1); H=E is unnecessary. The initial
audit checked the subsequent original-signing implication and final rational
conversion. The additional independent phase reconstruction and complete
directed-interval replay are recorded in Section 8 below.

Sources read completely for the relevant chain:

- `transfer_reconstruction_standalone_2026_09_06.md`, Sections 1--7;
- `decisive_audit_standalone_direct_E_upper_2026_09_07.md`;
- `decisive_bridge_improved_all_order_cap_audit_2026_09_06.md`;
- `decisive_director_gaussian_phase_all_order_upper_2026_09_07.md`.

## 1. The orbital estimate really is uniform over terminal bases

For P_t(v)=E_g exp(-t||v-gv||^2), L_t=sqrt(P_t), the Gaussian Fock feature
orbit covariance has norm P_t(v). Its finite normalized Gram matrix has
nonnegative entries and constant row sum P_t(v); its top eigenvalue is
that row sum. At degree at most D, the invariant rank is at most

```math
\sum_{j\le D/2}p(j)
\le \inf_{s>0}\exp(sD/2+\pi^2/(6s))
=\exp(\pi\sqrt{D/3}).
```

The feature tail is a Poisson tail of mean at most 2tCm. With
D=ceil(e^2(2tC+1)m), the Chernoff estimate is at most
exp(-(4tC+1)m), while P_t(v)>=exp(-4tCm). An orthogonal transformation
preserves degree. Projection onto the low-degree invariant space and
Cauchy--Schwarz therefore give

```math
\mathbb E_g L_t(Ugv)\le\exp(O_{t,C}(\sqrt m))L_t(v)
```

for every orthogonal U, with no entrywise assumption on v or U.
Consequently a terminal Hadamard containing arbitrarily many H_2/H_12
tensor factors is covered by the same estimate. Such terminal complexity
is not a growing binary Bellman depth.

The full signed-permutation invariant space is contained in the two-block
invariant space, so Gaussian Fock factorization gives the split inequality
L(v_+,v_-)<=L(v_+)L(v_-). Permanent terms fixing coordinate i give
P(v)>=K_t(v_i,v_i)P(v without i)/m, with K_t(a,a)>=1/2. Hence the maximum
deleted-coordinate L is bounded by sqrt(2m)L(v), including repetitions
and large deleted coordinates.

## 2. Fixed-depth counting and amplitudes

For a root absolute type, the number of signed words is
s!2^(nonzero count)/product c_a!, whereas an ordered pair table has
(s/2)!/product n_ab! realizations. Although its two signed marginals may
be unequal, their average absolute marginal is fixed. Since the reference
source nu is symmetric, this is enough for the cross-entropy identity

```math
\log\Pr(\pi)=-{s\over2}D(\pi\Vert\nu\otimes\nu)+O(\log s).
```

Simultaneous reversal and swapping make both signed marginals nu, leave
both child absolute laws unchanged, and decrease entropy. Independent
reversal of only one input is not required or used.

At fixed depth r every reachable alphabet is finite. In the application
its entries are finite linear combinations of the common root amplitude
1/sqrt(k/m), which tends to 1/sqrt(p); no new small or large amplitude
class appears as m grows. The number of types is polynomial and the sum
of terminal orbital errors is O_r(sqrt(m)). A terminal moment may be as
large as 2^r times the root moment, but r is fixed. Thus the normalized
logarithmic error is o(1), uniformly in all chosen terminal Hadamards.

## 3. Symmetry, diagonal, and the exact doubled-edge kernel

For arbitrary full Hadamards H_i (they need not be symmetric), define

```math
W_{(i,a),(j,b)}=S_{ij}H_i(a,j)H_j(b,i),\qquad S_{ij}=S_{ji}\in\{-1,1\}.
```

This is symmetric and every entry is a sign. Its diagonal entry in fibre
i is exactly S_ii, independent of a. Retaining k rows per fibre gives
N=mk. With h_i=H_i[T_i,:]^T x_i and E=x^T W_T x, row orthogonality gives

```math
\sum_i\|h_i\|^2=m^2k,\qquad
D_\sigma=\sum_{i,j}(h_i(j)-\sigma S_{ij}h_j(i))^2
          =2(m^2k-\sigma E).
```

For an undirected edge, its two squared terms are identical. In
exp[-tD_sigma/(2k)], averaging the fair edge sign therefore gives

```math
\frac{e^{-t(a-b)^2}+e^{-t(a+b)^2}}2,
\qquad a=h_i(j)/\sqrt k,\quad b=h_j(i)/\sqrt k.
```

There is no extra factor two. Diagonal defects are nonnegative and may
be discarded for an upper bound on this exponential moment, regardless
of the chosen S_ii or sigma.

Conditional on each independently permuted row's selected absolute
diagonal entry, the remaining multiset is uniformly permuted. A Gram
factorization of the folded kernel gives a vertex tensor equal to its
average over these permutations. Its squared norm is exactly
per K_t[b]/(m-1)!, since a pair of independent permutations has a uniform
relative permutation. Repeated values do not alter the quotient. Graph
Cauchy--Schwarz bounds the contraction by the product of these norms.

## 4. All spin choices and the expectation power

Let Z_i be the sum over all 2^k fibre spins of the maximum deleted-row L.
Conditional on the fibre bases, exponential Markov and the complete union
over spin vectors and the two energy orientations give

```math
\Pr\{\exists x,\sigma:D_\sigma(x)\le2\gamma m^2k\}
\le 2e^{t\gamma m^2}\prod_i Z_i.
```

Use physical H_i=sqrt(m)(U_m^(i))^T, with independent recursive U_m^(i).
The normalized spectrum is exactly
U_m^(i)(1_T x)/sqrt(k/m). Fresh signed input permutations make its law
the same for every fixed support and row spin. Linearity of expectation,
not independence among different spins, gives the factor 2^k in E Z_i.
The deletion inequality then shows

```math
\limsup_m m^{-1}\log\mathbb E Z_i
\le p\log2+f_r(\nu_p).
```

The fibre bases themselves are independent. Hence averaging their product
gives product_i E Z_i, or (E Z)^m when the selectors agree. It is never
E[Z^m]. Fixed different selectors of the same size have the same upper
bound and cause no change. The selectors are not selected adaptively.

## 5. Strict margin, cap normalization, and hollowing

Fix 0<a'<a and then a finite r with
p log2+f_r(nu_p)+4(1-sqrt(p))<=-a'. Fix 0<eta<a'/4. Set
gamma_m=1-sqrt(k/m)+eta. The failure exponent divided by m^2 has limsup
at most -a'+4eta<0. Thus for every sufficiently large admissible m an
actual full W_T exists with

```math
\max_x|x^T W_Tx|\le m^2k(\sqrt{k/m}-\eta).
```

The diagonal trace is exactly k sum_i S_ii and has magnitude at most N.
Its removal therefore costs at most N/2 in the hollow Hamiltonian cap.
Since N=mk, the resulting signing satisfies

```math
{Q(A)\over N^{3/2}}
\le {1\over2}-{\eta\over2\sqrt{k/m}}+{1\over2\sqrt N}.
```

The normalization is valid for both energy signs, and does not assume
the full matrix diagonal vanishes or that its terminal bases are symmetric.

## 6. Explicit all-order quantifiers

The skew Paley matrix over F_11 gives a nonsymmetric H_12=I+R with
R^T=-R and RR^T=11I. Together with H_2 it supplies all terminal orders
s=2^a12^b, a,b nonnegative. Irrationality of log12/log2 follows from
unique factorization.

For any epsilon>0 there is a finite set of nonnegative b whose residues
b log12 modulo log2 have maximum circular gap at most epsilon. For every
sufficiently large target log X, choose the next such residue and add a
nonnegative multiple a log2. This puts an available log s in
[log X,log X+epsilon]. The nonnegativity of a follows because the finite
b set was chosen before X tends to infinity. The available terminal
orders therefore have consecutive ratios tending to one.

Multiplication by the fixed 2^r preserves that property. For
N=m floor(pm), one has N/(pm^2)->1; thus the constructed orders also
have consecutive ratios tending to one. For every n take the least
constructed N>=n. A principal restriction has no larger hollow cap,
because averaging the omitted unbiased spins recovers any prescribed
restricted energy. The normalized loss is only (N/n)^(3/2)->1.

The complete quantifier order is therefore: choose a strict margin; fix
finite r; pass through all sufficiently large admissible m and then all n;
only then let the strict margin tend to zero. No uniform growing-depth
statement is necessary. This yields

```math
\limsup_n{M_n\over n^{3/2}}
\le {4+p\log2+g_4(1)\over8\sqrt p}.
```

## 7. Independent exact replay

The canonical cap conversion and standalone exact checker were rerun and
passed. In addition, the new independent program

`computations/decisive_independent_upper_realization_audit_2026_09_07.py`

uses its own rational logarithm series and 10^18-scale integer square-root
brackets. It gives a strictly tighter enclosure inside the published cap
and verifies the exact positive decimal gap

```math
{494515125\over10^9}-{7787631971809\over15748015748016}
={2161550371\over7874007874008000000}>0.
```

The numerator of the cap expression is checked positive; hence the upper
cap uses the lower sqrt(p) endpoint. No floating value enters an asserted
inequality. The approximate value is .49451512472517395.

The new integer checker also verifies the explicitly nonsymmetric H_12,
a recursive order-24 basis with nonsymmetric terminals, all 4096 spins
of an order-12 restricted weave, and 128 additional order-48 spin tests.
Both orientations, exact diagonal traces, the physical transpose, and the
doubled-edge defect identity are checked. The finite tests supplement,
but do not replace, the preceding all-spin and all-order proof.

No mathematical correction was found in the audited realization chain.
Its conclusion is an original all-order upper bound, not convergence or
an identification of the construction with optimal signings.

## 8. Additional independent phase and minimal-closure check

After completing the realization audit I read the precision supersolution,
the exact ternary Gaussian-phase proof, and its entire rectangle verifier.
The following were reconstructed directly.

For a fixed precision, the rate-distortion representation follows by
comparing each channel with its Gibbs kernel. At a maximizing reproduction
measure, the KKT witness q(y) is at most one and equals one on its support,
which makes the Gibbs channel's output marginal self-consistent. Symmetry
can be imposed by concavity. Its witness is
exp(-lambda y^2)[a+b cosh(2lambda y/sqrt(p))], with a,b strictly positive.

Writing z=2lambda y/sqrt(p) and r=a/b, its positive critical points solve
h(z)=p/(2lambda), where h=sinh(z)/[z(r+cosh(z))]. Direct differentiation
shows sign(h')=sign(r-R(z)). The numerator and denominator of R have
coefficient ratios 4^j/(2j), strictly increasing for j>=1. Pairing terms
in the quotient derivative proves R strictly increasing from 2 to infinity.
Thus h decreases, or first increases and then decreases. The derivative
of q has the sign of h-p/(2lambda); q has at most one positive local
maximum. At y=1/sqrt(p), its zero-source and negative-source Gaussian
terms have strictly negative derivatives, while its positive-source term
has derivative zero. This endpoint cannot be a maximum. Outside the source
hull all three terms decrease. Consequently the reproduction support
really is contained in {0,+a0,-a0}, rather than being an assumed ansatz.

The transformation from its central weight w to v is a bijection of
[0,1], reversing its endpoints. Substitution gives
G=p log(1+vB)-log(1+vA), with derivative numerator
pB-A-(1-p)ABv. This proves the exact clipped optimizer and all three
branches used in the program, including A=0 and B=0.

The normalized ternary source is unit subgaussian because its even moments
are coefficientwise dominated by Gaussian moments for p>=1/3. Entropy
duality gives I>=E(posterior mean)^2/2, so J_lambda=lambda for lambda<=1/2.
The maximizing precision of c_4(lambda)-lambda is
4(1-(sqrt(257)-1)/16)<1/2. At equality, the strict coefficient
1/2-lambda forces zero posterior-mean variance, and then zero information.

On lambda>=1/2, c_4(lambda)-lambda is decreasing: its derivative at 1/2
is 7/15-1<0 and it is strictly concave. G increases with B and decreases
with A for every v, hence also after maximizing v. The rectangular upper
bound is therefore in the correct direction. I inspected the verifier's
exact rational comparisons of the dyadic A/B endpoints; in its interior
branch these guarantee B>A before evaluating log(B-A). Float arithmetic
chooses only the subdivision and proposed prune. Every accepted rectangle
is rechecked by directed intervals; a wrong proposed prune fails an
assertion. Even rounded dyadic midpoints preserve the covering partition.

The installed mpmath 1.3.0 source was inspected: mpi_exp and mpi_log call
their scalar routines with round_floor at the lower endpoint and
round_ceiling at the upper endpoint. The complete fresh run of

`computations/decisive_bridge_gaussian_phase_rectangles_2026_09_07.py --verify`

finished in 93.54 seconds with 425009 tree nodes and 212505 accepted
rectangles, status `directed interval exclusion`. Every accepted box is
strictly below -777671/1000000, while the separately exact Gaussian lower
bound is -622136276211/800000000000, strictly larger.

Finally, the binary precision proof's weighted square and information
identities check directly. The arithmetic/harmonic parent precisions remain
admissible, preserve the product, and reduce logarithmic spread, so the
concavity of c_t(exp(s)) proves BE<=E without equal parent precisions.
In the direct stopping argument the true budget is Phi-G. Finite-depth
near-optimality bounds expected accumulated drift by Phi-G+zeta; the
second moment is an exact martingale along the uniform branch. The three
stopping contributions are bounded by epsilon, E0 sup_{s>C}K(s)/s, and
K(C)(Phi-G)/(r kappa). This only needs E as a supersolution, its Gaussian
boundary value, and moment-ball continuity. No lower comparison E<=f_r
or H=E identification is needed. No new gap was found in this extension.
