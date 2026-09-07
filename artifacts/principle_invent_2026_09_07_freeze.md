# Independent principle freeze — 2026-09-07

Status: authored research hypotheses, not proved claims. This note was frozen
from the naked problem and README before consulting the ledger, active theory
notes, or upper-bound construction. The assigned task supplies only
`0.4333221116640807 <= liminf M_n/n^(3/2)` and the reported strict upper endpoint
`0.493608094`, together with the child reversal identity.

Write `F(A)=max_{x in {+-1}^n}|sum_{i<j} A_ij x_i x_j|` and
`M_n=min_A F(A)` over symmetric hollow sign matrices.

## Three independent architectures

### 1. Renormalization through a nonclassical signing ensemble

A bounded `F(A)/n^(3/2)` forces the usual dense graphon to vanish, so ordinary
fixed-order subgraph statistics cannot identify the optimum. The relevant
state must retain microscopic spectral/eigenvector information. A plausible
positive mechanism is an explicit ensemble closed under arbitrary-order
realization, rich enough to approximate *selected* exact minimizers.

Exact proposed theorem target: there is a class `P` of dimension-independent
randomized signing generators, each producing a signing `A_N(P)` at every
sufficiently large integer `N`, such that:

1. For every `P`, `F(A_N(P))/N^(3/2)` converges in probability to `c(P)`.
2. For every sequence of exact minimizers `A_n`, some subsequence and some
   `P in P` satisfy `c(P) <= liminf F(A_n)/n^(3/2)`.

This would prove convergence to `inf_P c(P)`. It is not enough to define a
generator by taking a favorable subsequence and leaving its missing sizes
unspecified. The generator must be a concrete order-free probability law,
e.g. a matrix functional of independent high-dimensional rotations or finite
local gadgets. A first derivation target is to identify a strict-subhalf
generator and the retained statistic absent from the ordinary graphon.

Concrete falsifier: the proposed state identifies two realizable ensembles
with different Boolean ground-state constants, or the generator's
finite-dimensional re-embedding has a vector-spin gap bounded away from zero.

### 2. Extremal free energy with a lower-tail variational principle

Let `A` be a uniform random signing and define

```math
p_n(c)=\Pr\{F(A)\le c n^{3/2}\}.
```

The existence question is the location of the first nonempty lower-tail
level, not the ordinary quenched SK ground state. The natural speed of this
rare event is `n^2`, because the outer disorder has `n(n-1)/2` bits, while the
inner spin entropy is only order `n`.

Exact proposed theorem target: there is a lower semicontinuous rate function
`I:[0,infinity)->[0,log(2)/2] union {infinity}` such that, for every continuity
point `c`,

```math
-n^{-2}\log p_n(c)\longrightarrow I(c),
```

and the effective domain is an interval with left endpoint `c_*`, with
`p_n(c)>0` eventually for every `c>c_*` and `p_n(c)=0` eventually for every
`c<c_*`. The final support statement is indispensable: an LDP alone may not
see superexponentially rare isolated constructions.

Concrete mechanism to investigate: a variational law for
`sum_A (sum_x exp(beta Q_A(x)/sqrt(n)) + exp(-beta Q_A(x)/sqrt(n)))^(-gamma n)`.
The disorder sum has negative replicas and a second extensive scale. A proof
must justify an interpolation or variational principle at that two-scale
normalization; merely invoking the SK thermodynamic limit is invalid.

Concrete falsifier: competing dimension-dependent microcanonical classes
whose support edges oscillate while all positive-probability lower-tail
rates agree.

### 3. A finite-depth signing calculus, not tensor amplification

An auxiliary Hadamard/orthogonal lift of a fixed signing tends to enlarge the
Boolean spin set to vector spins. The classical `infinity -> 1` norm already
fails to multiply under tensor product: `H_2` has norm 2, whereas
`H_2 tensor H_2` has norm 8, exceeding 4. Repeated tensor amplification is
therefore not a neutral operation on a favorable Boolean constant.

Exact proposed theorem target: construct a finite list of signing operations
`T_1,...,T_r`, using signed inter-block patterns and optionally changing all
internal edges, with an explicit finite-dimensional response state `s(A)` and
maps `s(T_j(A))=R_j(s(A))+o(1)`. The associated normalized cost `C(s)` must be
nonincreasing along a dense-size realization path, and every selected exact
minimizer must admit such a state at cost `F(A)/n^(3/2)+o(1)`.

The first concrete derivation target is weaker but positive: identify the
exact state and operation behind a scalable strict-subhalf construction,
then determine whether its operation extends beyond that special carrier.
The state must retain the dependence of positive and negative internal
energies on bridge response. In particular it must respect

```math
\max_{\epsilon=\pm1}|a+b+\epsilon c|=|a+b|+|c|.
```

Concrete falsifier: a pair of signings with identical proposed states has
different amplified cost by a fixed fraction of `n^(3/2)`, or the operation
necessarily returns to the vector-spin/spectral half barrier.

## Preliminary ranking

Architecture 3 is the first working route because it has an immediately
inspectable positive instance: the known strict-subhalf construction. Its
objective is to extract a genuine reusable operation, not rename the
extension inequality. Architecture 1 is the ambitious global theorem that
could organize that operation. Architecture 2 is a more independent
variational fallback, but the support-edge obligation is a major issue.

No convergence result is claimed by this freeze.
