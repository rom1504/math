# Finite actual-child inverse-escort low-degree falsifier

**Status.**  Rigorous parity identity plus complete-cube numerical evidence on
exhaustively selected thermal-minimizing children.  This note does **not**
prove an asymptotic lower bound against degree `o(N)`.

## 1. Declared approximation problem

Fix the physical raw temperature and inverse escort

```math
t={4\over\sqrt N},\qquad \lambda=1,
\qquad {dq\over dU}(B)\propto P(B)^{-1}.
```

For a bridge edge `e`, let `r_e(B_{-e})` be the exact deleted-edge cavity
response.  For `K=1,3`, the experiment computes

```math
\mathcal E_{e,K}
=\min_{\deg f\le K,\ f=f(B_{-e})}
 E_q\{r_e-f\}^2.                                  \tag{LF.1}
```

This is the true weighted least-squares optimum, not the error of a chosen
Walsh truncation.  If `chi_S` are the Walsh characters not containing `e`,
the normal equations are

```math
G_{S,T}=E_q\chi_{S\mathbin\triangle T},
\qquad b_S=E_qr_e\chi_S.                          \tag{LF.2}
```

The escort has full support, so `G` is positive definite.

## 2. Exact parity orthogonality

**Lemma LF.1.**  For every antipodally symmetric actual-child rank-one law,

```math
q(B)=q(-B),\qquad r_e(-B_{-e})=-r_e(B_{-e}).       \tag{LF.3}
```

Consequently

```math
\boxed{E_qr_e\chi_S=0\quad\text{whenever }|S|\text{ is even}.} \tag{LF.4}
```

In particular the optimal degree-`2k` error equals the optimal
degree-`(2k-1)` error.

*Proof.*  Antipodal symmetry makes the cosh likelihood even under the global
bridge flip, hence so is every power escort.  In the deleted-edge formula,
replace the latent word `Q` by `-Q` after globally flipping the retained
bridge bits.  The denominator is unchanged and the numerator changes sign.
Under the same flip `chi_S` acquires the sign `(-1)^|S|`; (LF.4) follows by
pairing bridge points. `square`

The stored numerical parity errors are exactly zero at floating precision.
Thus the even-degree plateau in the experiment is a theorem, rather than an
empirical pattern.

## 3. Actual child selection and exhaustive calculation

For the splits

```text
N=6:  3+3,       N=8:  4+4,       N=10: 3+7,
```

the child signing spaces are enumerated exhaustively in root switching
gauge.  Exact absolute-energy histograms are compared at `t=4/sqrt(N)` with
80-decimal outward-rounded interval arithmetic.  Every competing histogram
has a strictly positive partition-sum gap.  There is one signed-permutation/
global-sign minimizing class in every child problem.  The JSON records the
exact signing counts, histogram counts, representative hashes, interval
certificate, and gap to the next histogram whenever one exists.

For each of both relative orientations the program then enumerates every
bridge: `2^9`, `2^16`, and `2^21` points.  Signed-automorphism vertex orbits
reduce only the repeated edge calculations, not the bridge cube.  The
reported edge-orbit sizes cover every bridge edge.

The degree-three results, pooled over both orientations and every bridge-edge
orbit, are:

| `N` | split | range of `E_(e,3)` | unexplained fraction of `E_q r_e^2` |
|---:|:---:|:---:|:---:|
| 6 | `3+3` | `0.0177978`--`0.0179984` | `2.71%`--`3.43%` |
| 8 | `4+4` | `0.0901517`--`0.0920519` | `18.93%`--`20.44%` |
| 10 | `3+7` | `0.0901513`--`0.0973884` | `21.92%`--`23.52%` |

For comparison, the degree-one residual ranges are

| `N` | range of `E_(e,1)` | unexplained fraction |
|---:|:---:|:---:|
| 6 | `0.151846`--`0.301360` | `23.16%`--`57.38%` |
| 8 | `0.277153`--`0.292111` | `58.19%`--`64.85%` |
| 10 | `0.298724`--`0.330753` | `69.20%`--`83.78%` |

The largest condition number of a displayed degree-three Gram matrix is
`132.20`; every maximum normal-equation residual is below `9e-16`.  Hence
the order-eight and order-ten separation is far larger than the observed
floating error.

This is a genuine finite falsifier of the proposal that the universal
rank-one rectangle response is already captured by degree three on the
actual optimizing-child escort.  It also shows that the positive raw overlap
at these orders is not merely the energy of the best degree-three cavity
component.

## 4. Why this does not scale yet

The data do not falsify degree `o(N)`.  Three finite orders cannot distinguish
a fixed positive residual from a crossover at a growing odd degree.
More importantly, none of the presently proved actual-child identities gives
a lower bound on (LF.1):

1. Mandatory Eulerian likelihood coefficients give the exact
   `p_{-e}^2`-weighted collision identity.  They do not survive replacement
   by the negative weight `P^{-lambda}`.
2. The escort may spend `Theta(N)` relative entropy to avoid the
   `Theta(N)`-rare positive-likelihood peaks responsible for the ordinary
   Walsh obstruction.
3. Child edge-flip minimality constrains prior Gibbs moments.  It supplies no
   inequality for the best projection of `r_e` after inverse retuning.
4. The finite calculation has not isolated one high odd-degree character or
   response packing whose `q`-mass remains bounded below with `N`.

Thus the precise missing falsifier is a uniform actual-minimizer statement

```math
\inf_{\deg f\le k_N}
 {1\over mn}\sum_e E_{q_\lambda}(r_e-f_e)^2\ge c>0
 \quad\text{for every }k_N=o(N),                 \tag{LF.5}
```

or a weaker version for one explicit `k_N -> infinity`.  Proving (LF.5)
requires a negative-tail optimizer identity absent from the current archive.
The complete finite evidence identifies the gap but does not bridge it.

## 5. Reproduction

From the repository root:

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_escort_low_degree_falsifier.py
```

Output:
[`../../computations/results/actual_child_escort_low_degree_falsifier.json`](../../computations/results/actual_child_escort_low_degree_falsifier.json).
