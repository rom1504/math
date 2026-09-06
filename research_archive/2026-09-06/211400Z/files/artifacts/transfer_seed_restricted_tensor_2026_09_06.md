# Random restriction of a literal seed tensor: inherited witnesses and barriers

Date: 2026-09-06. This route retains the actual seed matrix in
`H_s otimes B` and is distinct from the scalar PSD-kernel weave. The
inherited-witness theorem and the scoped upper-certificate barriers below
are proved. Numerical probes are lower witnesses only. There is no
seed-transfer upper theorem here.

## 1. Exact inherited-witness concentration

Fix a full symmetric sign seed `B` of order `d`. Let `H_s` be any
symmetric Hadamard matrix of order `s`, and put
`C_s=H_s otimes B`, of order `N=ds`. Let `K_s` be its hollowing.
For a fixed Boolean vector `x`, and independent Bernoulli-`p` vertex
selectors `eta_i`, define

```math
Z_x=\sum_{i<j}(K_s)_{ij}x_ix_j\eta_i\eta_j.
```

Writing `a_ij=(K_s)_ij x_i x_j`, direct expansion of the independent
selectors gives the EXACT variance identity

```math
\mathbb E Z_x=p^2P_{K_s}(x),
\operatorname{Var}(Z_x)
=p^3(1-p)\|K_sx\|^2
 +p^2(1-p)^2\sum_{i<j}(K_s)_{ij}^2.                         (1)
```

Indeed disjoint edge pairs have zero covariance, whereas two edges
sharing one endpoint have covariance `p^3(1-p)`. The identity
`sum_i (sum_(j!=i) a_ij)^2=2sum_edges a_ij^2+2sum_adjacent-pairs a_e a_f`
then gives (1).

At fixed seed, `||K_s||op<=sqrt(s)||B||op+1`, so (1) is `O_B(N^2)`,
uniformly in the chosen witness. Chebyshev therefore gives fluctuations
`o(N^(3/2))` with probability tending to one. Also `|T|/N -> p`.
Apply this to either objective sign of a maximizing full witness, paying
`O(N)` for hollowing. It follows that

```math
\frac{Q((K_s)_T)}{|T|^{3/2}}
\ge \sqrt p\,\frac{Q(C_s)}{(ds)^{3/2}}-o_{\mathbb P}(1).    (2)
```

The assertion also holds for a uniform selector of exact size `k=floor(pN)`.
For a direct proof, put `pi_j=(k)_j/(N)_j`. Separating identical, adjacent,
and disjoint edge pairs gives

```math
\operatorname{Var}(Z_x)
=(\pi_2-2\pi_3+\pi_4)\sum_e a_e^2
 +(\pi_3-\pi_4)\|K_sx\|^2
 +(\pi_4-\pi_2^2)P_{K_s}(x)^2
\le\sum_e a_e^2+\|K_sx\|^2.
```

The first two coefficients are probabilities of specified inclusion and
exclusion events and lie in `[0,1]`; the last is nonpositive by negative
dependence of two disjoint inclusion pairs. The mean is `pi_2 P_K(x)`.
Thus the same concentration conclusion follows without conditioning or a
local limit estimate.

## 2. Regularized seed witnesses survive with factor sqrt(p)

Use the existing cofinal regular outer family
`H_s=H_4^(otimes a) otimes H_144^(otimes b)`. Its normalized caps increase
cofinally to the actual regularized norm `R(B)`. Equation (2) gives

```math
\liminf_{s\to\infty}
\frac{Q((H_s\otimes B)_T\text{ hollowed})}{|T|^{3/2}}
\ge\sqrt p\,\frac{\mathcal R(B)}{d^{3/2}}
\quad\text{in probability}.                               (3)
```

Thus a typical random restriction cannot achieve a seed-cap landing
`Q(B)/d^(3/2)+epsilon_d` unless

```math
\sqrt p\,\mathcal R(B)/d^{3/2}
\le Q(B)/d^{3/2}+\epsilon_d.                               (4)
```

In particular, retention tending to one cannot remove a persistent
regularized-seed gap. If that gap is fixed, a nonvanishing fraction of
vertices must be removed. This is necessary, not sufficient: new
restriction-adapted witnesses are absent from (3).

The quantifier is important. Equation (3) concerns typical random
selectors. It does not rule out exceptional selected restrictions, and
it does not assume `R(B)=Q(B)` or `R(B)=T(B)`.

## 3. Two obvious upper bounds do not remove the gap

Put `r=R(B)/d^(3/2)` and `L=||B||op/sqrt(d)`, so `r<=L/2`.
The spectral upper bound on a restricted tensor is `L/(2sqrt(p))+o(1)`,
which is never less than `r`.

An expansion retaining the original cap is also insufficient. Write
`P_T=pI+Delta`; the squared norm `||Delta x||^2` is independent of
the Boolean `x` and is `p(1-p)N+o(N)`. Expanding the quadratic form and
bounding the two remainder terms spectrally gives

```math
\frac{Q((K_s)_T)}{|T|^{3/2}}
\le \sqrt p\,r
 +L\sqrt{1-p}+\frac{L(1-p)}{2\sqrt p}+o(1).                (5)
```

Since `L>=2r`, the right side is at least
`r[1/sqrt(p)+2sqrt(1-p)]>=r`. Therefore this elementary perturbation
bound cannot certify any reduction of the stabilized coefficient,
even though the inherited witnesses themselves decrease by `sqrt(p)`.

There is a stronger scope check for absolute PSD-majorant upper bounds.
If `D>=+B,-B`, then `sqrt(s) I_s otimes D` majorizes both signs of
`H_s otimes B`, and its principal restriction remains a majorant.
However every full symmetric sign matrix of order `n` has absolute
PSD-majorant value at least `n^(3/2)/2`. Hence ANY valid upper bound
obtained solely through such a majorant of the restricted tensor is at
least one half after normalization. Since the original all-order upper
bound is now strictly below one half, that entire upper mechanism is
insufficient for lossless transfer of asymptotically minimizing seeds.
This last statement concerns PSD-majorant certification, not the actual
restricted matrix.

## 4. Bounded finite lower-witness diagnostic

The companion C++ program constructs actual hollow principal restrictions
of `H_4^(otimes r) otimes B` and searches BOTH objective signs. The order-six
seed is the exact hollow conference minimizer of normalized code220,
completed by the diagonal `(+,+,+,-,-,-)`; its trace is zero and its full
Boolean cap remains exactly five. The order-two seed is `H_2`, of cap one.
All reported energies are recomputed by a separate integer quadratic
evaluation inside the program. No upper estimate is produced.

At RNG seed20260906 and30 restarts per objective sign, the following
lower witnesses were found:

| Full seed | Outer order | Retention | Restricted order | Witness energy | Normalized lower witness |
|---|---:|---:|---:|---:|---:|
| `H_2` |1024|.95|1945|40602|.4733343455|
| `H_2` |1024|.80|1638|33439|.5044085450|
| `H_2` |1024|.60|1228|23906|.5555323640|
| `H_2` |1024|.40|819|14397|.6142513958|
| completed order-six minimizer |256|.95|1459|33271|.5970118019|
| completed order-six minimizer |256|.80|1228|25840|.6004750391|
| completed order-six minimizer |256|.60|921|17396|.6223861987|
| completed order-six minimizer |256|.40|614|9959|.6545806690|

These finite seeds have substantial fixed-order completion/regularization
effects, so their values are not extrapolated to growing near-optimal
seeds. In particular, the sub-half LOWER witness in the first row is not
an upper construction below one half. The diagnostic suggests that simply
reducing retention is not a reliable way to remove the seed's entangled
witnesses: at moderate retentions, other configurations already give
substantially larger values.

Reproduce, for example:

```text
g++ -O3 -std=c++17 computations/transfer_seed_restricted_tensor_probe_2026_09_06.cpp -o /home/math/quadra/tmp/transfer_seed_restricted_tensor_probe_2026_09_06
/home/math/quadra/tmp/transfer_seed_restricted_tensor_probe_2026_09_06 2 5 0.8 30
```

## 5. External hypothesis check and current remaining task

The nearby rigorous orthogonally invariant spin-glass literature does
not currently supply the needed transfer theorem. Its coupling law is
Haar orthogonally invariant, and the cited results are high-temperature
free-energy/TAP theorems, not uniform ground-state theorems for exact
restricted sign tensors. See
[Fan and Wu](https://arxiv.org/abs/2105.02797) and
[Fan, Li and Sen](https://arxiv.org/abs/2202.09325).
No theorem from those papers is imported here.

The next useful statement would be an UPPER bound on restriction-adapted
Boolean witnesses at the nonvanishing deletion required by (4), with
leading coefficient `Q(B)` rather than `R(B)` or a spectral/PSD norm.
Neither the concentration theorem nor the finite lower searches provide
it. Further vector-profile notation without such an inequality would not
advance the original convergence gap.
