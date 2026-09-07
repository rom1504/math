# A finite joint-profile certificate for the rank-two sign construction

Date: 2026-09-07. Status: proved finite theorem, independently reconstructed
by both the construction and adversarial researchers. See
`flatify_adversary_2026_09_07_typed_and_regular_support_audits.md`.
This is a rigorous finite upper-certificate mechanism. Its asymptotic
hypothesis is NOT yet verified and no new bound on M_n is claimed here.

## Actual ensemble and normalization

Let k=2m be a Hadamard order and N=mk. Fix any order-k Hadamard F.
Independently in every fibre use F_i=F P_i D_i, where P_i uniformly permutes
columns and D_i supplies independent column signs. Use the rank-two cross
tiles from `flatify_construct_2026_09_07_rank_two_weave.md`, and set all
within-fibre entries temporarily to zero. If x_i is a fixed Boolean row,
put a_i=|F^T x_i|/sqrt(k), with empirical magnitude law rho_i. The random
feature vector is a uniform permutation of this word with fresh fair signs.

Write Z=x^T C_cross x/sqrt(N). One cross edge between fibres with feature
magnitudes (a,b) and (c,d) contributes

```math
sqrt(2)[ac epsilon_1 eta_1+ad epsilon_1 eta_2
       +bc epsilon_2 eta_1-bd epsilon_2 eta_2].
```

Its exact exponential kernel K_t is the average of the exponential of t
times this quantity over the four independent fair signs. In particular

```math
K_t(a,b;c,d)=prod_{w in {ac,ad,bc,bd}} cosh(sqrt(2)t w)
            -prod_{w in {ac,ad,bc,bd}} sinh(sqrt(2)t w).       (1)
```

The whole correction in (1) is retained. K_t is positive, but no
positive-semidefinite kernel hypothesis is used below. Different unordered
cross edges use disjoint column signs. The two loop coordinates per fibre
are unused in C_cross.

## Exact finite inequality

For a realizable magnitude law rho with multiplicities n_a summing to k,
let

```math
p_k(rho) = k! prod_a rho(a)^(n_a) / prod_a n_a!.
```

This is the probability that k iid rho samples have exactly that type.
For any real potential phi_i on the finite support of rho_i, define

```math
L_ij = log E_{rho_i^2 tensor rho_j^2}
       [K_t(a,b;c,d) exp(phi_i(a)+phi_i(b)+phi_j(c)+phi_j(d))],
U_t(rho_1,...,rho_m;phi_1,...,phi_m)
 = sum_{i<j} L_ij
   +sum_i [2 log E_{rho_i} exp(phi_i)
           -k E_{rho_i}phi_i-log p_k(rho_i)].                (2)
```

Then, for every fixed full spin x having these row types,

```math
log E exp(t Z) <= U_t.                                    (3)
```

Proof: replace each uniform fixed-type row word by iid rho_i coordinates
conditioned on its exact type. On the conditioning event its sum of phi_i
is the deterministic k E rho_i phi_i. Insert this exponential factor,
then drop the conditioning indicator in the nonnegative expectation and
divide by prod_i p_k(rho_i). Under the unconditioned product distribution,
all cross-edge four-tuples and all unused loop coordinates are independent.
The former give L_ij; the latter give the two single-coordinate factors.
This proves (2)--(3). Potentials may be optimized separately for each tuple
of row types; no minimax interchange or attainment is required.

## From row counts to an actual full-sign upper bound

Let R_k be the finite set of magnitude types occurring among F^T x/sqrt(k)
for x in the k-cube, and let c_k(rho) count their Boolean preimages. If,
for some t>0 and q>=0,

```math
2 sum_{(rho_1,...,rho_m) in R_k^m}
 exp[sum_i log c_k(rho_i)+inf_phi U_t-2t q N] < 1,           (4)
```

then some actual choice of signed column permutations has
Q(C_cross)<q N^(3/2). Indeed the cross law is invariant under reversing
all edge contributions, by reversing both column signs at one chosen end
of each edge. Apply Markov in both polarities and union over the counted
Boolean rows. For a nonattained infimum choose potentials sufficiently
close to it using the strict margin in (4).

Fill each fibre with any actual O(k^(3/2))-cap signing. Triangle inequality
adds O(m k^(3/2))=O(N^(5/4)) for k=2m. Thus (4) is a genuine actual-sign
certificate with power-saving within-fibre completion, not a statement
about weighted entries alone. An O(k^(3/2)) signing exists, for example,
by the elementary random-sign union bound; the current stronger upper
construction is not required for this completion step.

## Uniformly subleading type overhead

For Hadamard F the unnormalized magnitudes |F^T x| are integers and their
squared sum is k^2. If D distinct positive magnitudes occur, then
1^2+...+D^2<=k^2, so D=O(k^(2/3)). The standard multinomial type lower
bound p_k(rho)>=(k+1)^(-|support rho|) therefore gives

```math
-sum_i log p_k(rho_i)=O(m k^(2/3) log k)=o(N),             (5)
```

uniformly in all row types. Also |R_k|<=exp(O(k^(2/3))). For completeness,
ignore the constraint on the number of nonzero entries and bound the number
of solutions to sum_{j>=1} j^2 h_j=k^2 by its generating product:

```math
log prod_{j>=1}(1-exp(-s j^2))^(-1)
 = sum_{r>=1}(1/r)sum_{j>=1}exp(-s r j^2)
 <= (sqrt(pi)/2) zeta(3/2) s^(-1/2).
```

Multiplying by exp(s k^2) and taking s=k^(-4/3) proves the assertion.
The zero multiplicity is then determined by the total word length.
Consequently log |R_k^m|=O(m k^(2/3))=o(N). In (4) the maximum type-tuple
exponent, plus this explicit subleading term, can replace the sum.

## Exact unresolved obligation and scope

The k<=16 experiments use genuine c_k but only a small set of homogeneous
and two-type trials. They are not a verification of (4) for growing k.
Neither Gaussian moments nor an arbitrary source entropy replaces c_k.
Diluted perfectly aligned relaxed sources already falsify doing that.
Heterogeneous row types and energy in unused loop coordinates are handled
in (2), rather than omitted: the latter have their exact exponential factors.

This certificate needs row counts (at most 2^k inputs) and four-coordinate
transport kernels, not the 2^(mk) parent landscape. It is nevertheless only
an upper-construction certificate. Independent column signs erase an outer
seed signing; proving (4) alone would not transfer actual liminf children.
