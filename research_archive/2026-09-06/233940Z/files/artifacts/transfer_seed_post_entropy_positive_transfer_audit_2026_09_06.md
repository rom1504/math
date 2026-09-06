# Post-entropy positive seed-transfer audit and a weighted decoder boundary

Date: 2026-09-06. Working research record, not a convergence theorem.
The useful flat threshold-entropy theorem is banked separately. This
note preserves the subsequent attempts to turn it into a positive
all-order CAP transfer, including the exact point where they stop.

## 1. The direction of the new entropy theorem

The Gaussian smoothing/resampling proof constructs a probability law
on Boolean spins of the SAME signing `A`. Gibbs variational converts
its entropy and correlation into a LOWER bound for that matrix's
partition function. This improves an upper construction only because
the already-built weave ensemble separately supplies an annealed
UPPER partition bound.

The only change of order inside the entropy proof is conditioning a
Gaussian likelihood on a coordinate subset. Its norm is a determinant
of a Gaussian principal covariance, not the cap of a newly constructed
full signing. The determinant comparison uses an operator norm and
flat squared coefficients. There is no inequality there transporting
an arbitrary seed's Boolean cap to a larger full signing.

Thus inserting `Q(seed)` as an unexplained substitute for the seed's
Gaussian or Gram quantity would simply reintroduce the previously
identified Boolean/vector gap. No such substitution was made.

## 2. A cap-only block-decoder idea and its exact obstacle

One possible positive architecture is to split a large target into
`d` blocks of `ell` physical spins. For each seed edge `ij`, choose a
flat cross matrix `H_ij` and put the seed sign on that entire block.
For physical block spins `x_i`, the normalized cross quantities are

```math
C_{ij}(x)=\ell^{-3/2}x_i^T H_{ij}x_j.
```

Hadamard blocks ensure each separate `C_ij` lies in `[-1,1]`.
That alone does not let one apply the original seed cap to their
simultaneous collection. A sufficient universal decoder would express
all these correlations as a mixture of seed-spin products, with an
additional common energy sign allowed for the absolute cap. Its
values would have to lie in the convex hull of the signed cut vectors.
Gram constructions produce an elliptope-type condition instead, which
is weaker. Introducing that universal decoder is therefore not itself
progress unless its required Boolean compatibility can be proved.

There is a scalable obstruction to the stronger requirement of a
cap-contractive block lift for EVERY WEIGHTED real seed. Let the seed
be the four-edge frustrated cycle, written bipartitely as

```math
C=\begin{pmatrix}1&1\\1&-1\end{pmatrix}.
```

Its same-spin hollow quadratic cap is exactly its bipartite norm,
`Q(C)=2`. After replacing each of its four edges by ANY full sign
`ell`-by-`ell` cross matrix, the lifted form is the bilinear form of
an arbitrary full sign matrix `D` of order `2ell` on each side.
For independent uniform left spins `z`, every coordinate of `D^T z`
has the law of a sum `S_(2ell)` of `2ell` independent signs. Hence

```math
Q(\mathrm{lift})=\beta(D)
\ge\mathbb E\|D^Tz\|_1
=2\ell\,\mathbb E|S_{2\ell}|
=\left(\frac4{\sqrt\pi}+o(1)\right)\ell^{3/2}.           (1)
```

The proposed unit-loss transfer would instead give
`2 ell^(3/2)`. The limiting multiplicative gap is `2/sqrt(pi)>1`.
This proof does not require the four cross blocks to be equal,
Hadamard, independent, or selected by a linear rule.

The scope is essential: the four-vertex seed here has zero missing
edges. It is NOT a full signing of the original complete graph.
Therefore (1) rules out an overstrong weighted cap-isometry or
universal positive decoder; it does not falsify a nonlinear construction
adapted only to a chosen sequence of full near-minimizing sign seeds.

The exact finite expectation used above is
`E|S_(2ell)|=2ell binom(2ell,ell)/2^(2ell)`.
Already at `ell=2`, the averaged lower bound is `6>2*2^(3/2)`.

## 3. A low-rank seed-memory covariance, and why it is not a construction

A separate attempt was to retain actual seed cap information in a
replicated covariance. If `A` is a fixed hollow seed of order `d`,
then on `d` blocks of size `ell`,

```math
R=I+\frac{s}{\sqrt d}\left(A\otimes\frac{J_\ell}{\ell}\right)
```

has unit diagonal and is positive when `|s| ||A||op/sqrt d<1`.
For every physical Boolean spin `x`, its seed term obeys the exact
bound

```math
|x^T(R-I)x|\le\frac{2|s|\ell}{\sqrt d}Q(A).
```

Indeed the block sums divided by `ell` lie in `[-1,1]^d`, and a
multilinear quadratic form attains its absolute maximum on the cube's
vertices. Thus this covariance remembers the ACTUAL seed cap, not
its spectral norm, in that particular energy calculation.

But this is only a covariance. No square flat sign basis with this
Gram matrix and a controlled Boolean weave cap was constructed.
Taking many independent sign columns would approximate its Gram only
at a much wider rectangular aspect ratio, and would leave the actual
large-signing cap obligation untouched. In addition its seed-dependent
perturbation has rank at most `d` when `ell` grows. No assertion was
made that this changes the extensive generic row-pressure baseline.
This attempt is preserved as an uncompleted calculation, not promoted
to a seed-transfer state.

## 4. Current genuine gap

The requested positive comparison still needs an ACTUAL full signing
at arbitrarily large, relatively dense orders with normalized cap at
most that of one liminf-realizing seed plus vanishing seed-order loss.
The threshold theorem supplies a useful same-order pressure witness;
the block-decoder attempt requires too strong a weighted compatibility;
the low-rank covariance does not supply a flat signing realization.

Sharp vertex insertion and ordinary degree-two Sidon tensor machinery
were screened against the archive and not rerun: their precise
obligations and thick-cap/vector-gap failures are already documented.
The fresh convergence track is independently investigating a lower
reverse-Fekete inequality for `M_n^(2/3)`, a distinct order-aggregation
architecture. No positive upper transfer or convergence implication
has been obtained in this post-entropy segment.

## 5. A finite below-half aggregation diagnostic, not an asymptotic falsifier

The committed signing in
`computations/results/heuristic_m14_from_conference.json` has exact cap
21. Among its 1,716 unordered partitions into two sets of size seven,
exactly 624 have both principal child caps equal to nine; the other
1,092 have both caps equal to eleven. One low-cap partition is
`(0,8,9,10,11,12,13)` and `(1,2,3,4,5,6,7)` (zero-based indices).
These facts are exhaustively verified by
`computations/transfer_seed_subhalf_partition_diagnostic_2026_09_06.py`.

Both low child caps are strictly below the half benchmark since
`4*9^2 < 7^3`. Exact zero-error powered aggregation fails since
`21^2 < 8*9^2`, equivalently `21^(2/3) < 2*9^(2/3)`.
This only shows that the below-half qualification does not rescue a
FINITE zero-error pointwise inequality. There is no amplification here
preserving these constants, so it says nothing against the proposed
asymptotic `o(N)` pointwise version, or the actual-minimum reverse-Fekete
target with an `O(sqrt(N))` defect. Self-tensoring this finite signing is
not a valid amplification, by the banked even-power classification.

## 6. Additional screened upper architectures

An ordinary weighted blow-up `A tensor J_ell/sqrt(ell)` has exactly the
desired homogeneous cap: its block means lie in the Boolean cube, so
its cap is `ell^(3/2) Q(A)`. But its entries have magnitude
`1/sqrt(ell)`, not one. Independently rounding its entries to signs
introduces an order-`N^(3/2)` fluctuation term, not a vanishing term.
The seed bound therefore does not survive by an elementary rounding
estimate. A dependent rounding theorem with the required cap control
would still be the principal missing construction; it is not supplied
by Gaussian covariance rounding.

A possible enforcement gadget was also screened. Divide `N=d ell`
vertices into d groups and attempt to enforce a small set of Boolean
group codewords using within-group interactions. Its worst possible
internal absolute energy is `O(d ell^2)`, whose normalized scale is
`O(sqrt(ell/d))`. That overhead can be small when `ell<<d`.
However a basic ferromagnetic penalty then has only an `O(ell)` local
enforcement field, whereas a generic flat bounded-norm bridge has
natural local-field scale `sqrt(N)=sqrt(d ell)`. The naive sufficient
enforcement regime `ell>>d` conflicts with the small-overhead regime.
This scaling calculation does not rule out sophisticated gadgets:
neither scale is asserted to be a necessary bound for every design.
It records why the simple proposed construction was not a proof.
The two-sided absolute cap also prevents silently discarding a large
negative penalty as one could in a one-sided ground-state reduction.

No theorem of rounding, coding, or enforcement is imported in this
section, and no new obstruction for all original near-minimizers is
claimed. These calculations close the bounded post-entropy architecture
screen rather than supplying a positive convergence step.
