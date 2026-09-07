# Independent reconstruction: first synthesis-to-principle results

2026-09-07. Root independently read the complete proofs below and reconstructed
their mathematical steps. This audit does not establish convergence. No earlier
positive verdict was used as a proof premise.

## 1. Logarithmic local law — PASS

Source: `principle_synthesis_2026_09_07_logarithmic_local_law.md`.

For the with-replacement sampler, the next row's odd Fourier coefficients
vanish by its fresh vertex switch. Averaging the square of an even coefficient
of order j over the old independent vertex labels gives the average j-th power
of the row correlation r_vw. Summing subsets gives exactly the even binomial
expansion in the source. Since |r_vw|<=1, every even j>=2 term is bounded by
r_vw^2. The entropy chain rule is used in the correct direction: conditioning
on the full latent past can only increase the divergence from a uniform fresh
row. Old labels need not remain independent conditional on observed edges.
They are independent under the unconditional average which the bound uses.

Summing t=1,...,k-1 gives 2^(k-1)-k, including the zero k=1,2 cases. Pinsker
and the separate repeated-label coupling give precisely the displayed TV
bound. The Schatten-four implication uses B=A+I and beta(B)<=4Q(A)+n;
there is no diagonal normalization error. The support-count converse has
the correct threshold above twice log_2 n and concerns one fixed orbit.

The actual cap is preserved on every orbit point, including exact minimizers.
Therefore small local laws cannot identify that cap. The theorem does NOT
show that the optimum requires a local-law description, or exclude selecting
exceptional principal restrictions. The iid-greedy consequence uses a
separate probabilistic lower bound, not an upper estimate for typical caps.

## 2. Global monotone balancing — PASS

Source: `principle_synthesis_2026_09_07_global_balancing.md`, Sections 1--6.

The exact equality R_t=R_0+2t follows for every eligible prefix by evaluating
the same negative maximizer and the matching universal upper bound. The gap
decreases by 0,2,4; first crossing is consequently in {-2,0}. Its mean is
(1-p)A-p(yy^T-I). Since the latter rank-one form has quadratic minimum -n/2,
the positive mean maximum never increases when P=Q>=sqrt(binomial(n,2)).

For fixed t, conditioning independent Bernoulli(t/L) choices on their count
costs at most L+1; the exact count is a binomial mode. The proof pays this
factor and ALL deterministic prefixes before choosing an adaptive crossing.
Scalar variance is <=4T and summand bound 2. Matrix variance is at most
4T(n-1)/L and summand norm 2. Substitution into Bernstein with
sqrt(2 variance*u)+(4/3)u gives the advertised thresholds. The two failure
probabilities are each <=1/8. Root checked the matrix theorem against
[Tropp's primary paper, Theorems 1.4 and 6.1](https://tropp.caltech.edu/papers/Tro11-User-Friendly-preprint.pdf).

With w=sqrt(n(Delta+n)), T=ceil(Delta/2+64w), eventual L>=n^2/8 and u<=2n
give E<40w, so the crossing condition is strict. The matrix error is O_C(sqrt
n), while the cap error is O_C(n^(5/4)). The latter is measured from Q,
NOT from half the energy range. The operation does not preserve exact
optimality; it constructs balanced additive near-minimizers. Its joint use
with spectral-core recovery retains that theorem's fixed-accuracy limit order.

## 3. Finite-type canonical seed universality — PASS

Source: `principle_invent_2026_09_07_finite_type_seed_universality.md`.

The entropy maximization has independent edge variables and aggregate local
type constraints. Positive product feasibility and positive kernels put its
unique optimizer in the relative interior. The dual multiplier normalization
gives the displayed infimum of Psi; the lower multiplier-ratio estimate is
uniform because each kernel entry ratio lies in [kappa,1/kappa]. Thus each
edge atom is >=alpha^2 kappa^3, independently of the number of rows.

On the exact type event the likelihood ratio is constant apart from the edge
kernel product, with the sign in the partition identity as written. Repairing
O(m^(3/2)) incidences costs at most that many edges; the atom ratio and Hamming
preimage count give exp[-O(m^(3/2)log m)], in the needed lower-probability
direction. Zero type masses mean restricting the relevant local support.

For the pressure comparison, log(r^T K s) is smooth on a neighborhood of
the product simplex because K is strictly positive there. A smooth compact
extension has absolutely summable Fourier coefficients. The separated series
therefore bounds its pairing with a matrix by a constant times the bilinear
Boolean norm; four real/imaginary terms suffice. This is uniform over every
choice of the local simplex vectors, so it passes through the infimum.
All constants depend on alphabet/kernel; no growing-condition-number claim
follows. Annealed row-type pressure is not the actual parent maximum.

## 4. Signed Eulerian entropy — independent second reconstruction available

Root's proof is `principle_director_signed_eulerian_entropy_2026_09_07.md`.
The invention agent separately reconstructed its transition count and
last-exit injection in `principle_invent_2026_09_07_transition_count_audit.md`.
Both reconstructions import only the ordinary Eulerian-orientation lower
bound, whose exact even-degree hypotheses root checked in
[Borbényi--Csikvári, Theorem 1.1](https://arxiv.org/html/1905.06215).
The general graph is connected and loopless. The required object is one
Euler trail, not a simple-cycle decomposition. This distinction avoids a
false signed-cycle decomposition premise.

## Remaining work

The fixed-alphabet hard compiler is a later extension and needs its own
audit. The results above neither preserve a seed's normalized cap at new
orders nor prove a value characterization. They justify further mathematical
work, not a solution claim or a campaign stop.
