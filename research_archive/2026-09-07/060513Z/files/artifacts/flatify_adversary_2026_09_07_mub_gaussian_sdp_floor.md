# A scalable Gaussian SDP lower test for MUB fibre constructions

Status: elementary theorem for the specified construction. No actual selected-
minimizer obstruction has been found from it. Developed while independently
checking the constructive agent's mutually-unbiased-basis diagnostic.

Let A be a fixed hollow symmetric full signing of order k. For each d suppose
U_1,...,U_k are real orthogonal matrices of order d with every entry of
U_i U_j^T equal to +-1/sqrt(d) for i!=j. Build a full parent signing C with
off-fibre blocks

```math
C_{ij}=\sqrt d\,a_{ij}U_iU_j^T,
```

and **arbitrary** hollow full-sign diagonal blocks. Define

```math
\operatorname{SDP}_{abs}(A)
=\max_{R\succeq0,\ R_{ii}=1}
 \left|\sum_{i<j}a_{ij}R_{ij}\right|.
```

For every d>=2,

```math
{Q(C)\over d^{3/2}}
\ge {2\over\pi}\operatorname{SDP}_{abs}(A)
   -{2\over3\pi d}\binom{k}{2}.                    (1)
```

## Proof

Fix any correlation matrix R. Construct jointly Gaussian d-vectors g_i with
`Cov(g_i,g_j)=R_ij I_d`, and choose Boolean fibre spins
`x_i=sign(U_i g_i)`. Each individual U_i g_i has independent standard normal
coordinates. Consequently the mean energy of **every** diagonal full-sign
block is zero, without knowing or choosing its entries.

For a cross-fibre coordinate pair, the correlation is
`R_ij (U_i U_j^T)_{ab}=+-R_ij/sqrt(d)`. The elementary Gaussian sign identity
(from rotational symmetry of a two-dimensional Gaussian) is
`E sign(Z)sign(W)=(2/pi)arcsin(corr(Z,W))`. Multiplying by the corresponding
full-sign entry and summing the d² coordinate pairs gives exactly

```math
\mathbb E[x_i^T C_{ij}x_j]
={2\over\pi}d^2 a_{ij}\arcsin(R_{ij}/\sqrt d).
```

For |z|<=1/sqrt(2),
`|arcsin(z)-z|<=|z|³/3`: integrate
`1/sqrt(1-s²)-1<=s²`, whose rationalized denominator is at least one in
this interval. Since |R_ij|<=1, summing the errors gives

```math
\left|\mathbb E H_C(x)
-{2\over\pi}d^{3/2}\sum_{i<j}a_{ij}R_{ij}\right|
\le {2\over3\pi}\sqrt d\binom{k}{2}.
```

Use Q(C)>=|E H_C(x)| and optimize R. This proves (1).

## Meaning for seed transfer

Along any unbounded-dimensional family with k fixed, a coefficient-one seed
transfer `Q(C)<=d^(3/2)Q(A)+o(d^(3/2))` requires

```math
\operatorname{SDP}_{abs}(A)\le(\pi/2)Q(A).
```

Thus a seed violating this inequality would give a scalable obstruction to
this MUB construction, even allowing arbitrary within-fibre completion. It
would not disprove arbitrary global flatification. To challenge selected
actual minimizers, the seed must additionally have the required genuine
optimality status; no such seed is asserted here.

The affine-plane/Hadamard MUB family used by the constructive diagnostic has
d=q² and q+1 bases, so for fixed k the standard power-of-two construction
can supply arbitrarily large dimensions. The theorem itself only requires
the explicitly stated orthogonality and cross-entry identities.

## Numerical screen and failed obstruction

`tmp/flatify_adversary_2026_09_07_seed_sdp_probe.py` checks both SDP polarities
for the stored full-sign witnesses at orders 3 through 10. Their Boolean caps
are exhaustively replayed. CLARABEL's SDP outputs are numerical, not exact
certificates; all matrices and solver statuses are preserved in its JSON.

The SDP/cap ratios range from 1.00 to approximately 1.39755, below pi/2.
In particular, the pentagon's value is approximately 5sqrt(5)/8=1.397542486.
This bounded test finds no actual-minimizer obstruction. The small exact
minima were independently enumerated only through order seven in this task;
the higher witnesses' global minimality remains their archived provenance.

An initially considered spectral-spike shortcut does not help: changing
o(k^(3/2)) total edge l1 mass changes both Q and SDP_abs by at most that
amount. The known planted-clique/universal-row near-minimizer perturbations
can create large spectral norms but cannot by themselves manufacture a
fixed SDP/cap gap. That attempted inference was rejected before use.

The Python environment emitted unrelated GLOP/PDLP import-version warnings;
the explicitly selected CLARABEL solves all completed with status optimal.

The supplementary `tmp/flatify_adversary_2026_09_07_all_n7_sdp_screen.py`
exhausts all 3,240 first-row-gauged exact order-seven minimizers and bins them
by integer traces of powers one through seven. The six bins have counts
360,840,420,420,840,360. Solving one numerical SDP per bin gives values
approximately 12.228857226,10.362371463,10.5 and their polarity copies;
all remain below (pi/2)*9. Spectral binning here is only a candidate-selection
device: this calculation does not assert isomorphism completeness or certify
an upper SDP bound for each omitted candidate.

For literature scope, the primary abstract of Friedland--Lim,
[Symmetric Grothendieck inequality](https://arxiv.org/abs/2003.07345), was
checked. It explicitly separates unrestricted symmetric matrices from the
positive-semidefinite cone, where the Nesterov pi/2 theorem applies. No PSD
pi/2 bound is silently imported for our hollow indefinite sign seeds. This
literature check is not a dependency of the elementary Gaussian proof above.
