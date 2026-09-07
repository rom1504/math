# Hybrid exact-support bridge: arbitrary paired spectra are paid, not discarded

2026-09-07. An actual full-sign construction step. The main agent derived
the hybrid potential; independent researchers reconstructed the conditional
integration and pin accounting. The theorem controls a declared sector,
not the whole parent cube or the original convergence problem.

## 1. Exact physical construction and target sector

Use the SAME center-adapted Walsh frames, active/inactive frequency pairs,
four reserved G fibres, and nonreciprocal center-pin layout as
`flatify_independent_2026_09_07_mixed_paired_bridge_sector.md`, Sections1--2.
In particular k=2^d=2m, n=mk, alpha=4/5, s=3/5. Each physical row has a
fixed nonzero translation a_i and a frequency b_i with a_i dot b_i=1.
Pair j with j+b_i, listing a_i dot j=0 first. The zero-frequency pair
is pinned to a reserved G row at the opposite side. No edge is pinned
at both ends. All other pair permutations and column signs are chosen
as in that actual ensemble. C is a full n-by-n sign matrix,
||C||op=sqrt(n), and x0^T C y0=0 EXACTLY.

For a target spin, a row may be declared:

* G: its pinned pair and the remaining ordered absolute-pair histogram
  obey the Gaussian residual profile condition (4) of the cited theorem.
* P: its relative physical word eta=x0*x is invariant under the FIXED
  translation a_i. Require only

      d_i^2=(u0-alpha sqrt(k))^2/k,

  where u0 is the actual zero-frequency coefficient. There is NO Gaussian,
  moment, empirical-spectrum, or distributional condition on the remaining
  P coefficients. Their inactive coordinates vanish exactly, including
  the second coordinate of the pinned pair.

Require average d_i^2<=tau_k^2 on each side, tau_k->0. The four reserved
rows must be G. Every other row may be G or P, and the SAME bridge will
control the union over all such declarations. The translations and
pairings are fixed before the target spin is chosen.

## 2. Statement

For all sufficiently large compatible orders there exists one output C
of this ensemble such that, simultaneously for all target pairs in the
above union of sectors,

    |x^T C y| <= (18/25)n^(3/2),       x0^T C y0=0.      (1)

This strengthens the mixed Gaussian/paired-Gaussian theorem by removing
the ENTIRE residual-profile hypothesis from P rows. It includes highly
non-Gaussian paired words, arbitrary additional periodicities, and sparse
or multilevel Walsh spectra, provided the fixed pairing constraint and
mean condition hold. It does not choose a new translation for each target.

## 3. Couple only G rows; keep P words exact

For G rows use the exact-law ordered-pair matching/conditioning argument
from the cited proof. Their comparison residual coordinates are independent
N(0,s^2); the first pinned coordinate is alpha sqrt(k). The additional
pinned partner and empirical approximation cost are o(sqrt(k)) per row
in squared-average scale. Dropping conditioning costs at most 2 per G row.

For P rows retain their ENTIRE actual coefficient word, with its actual
pair permutation and signs, and only replace u0 by alpha sqrt(k). This
costs exactly k d_i^2 in squared Euclidean error. There is no iid surrogate,
no empirical conditioning, and no independence assumption on P coordinates.

Consequently, on each side the physical and comparison feature arrays
differ by at most eta sqrt(n), eta=tau_k+o(1). Their squared norms differ
by at most e n, e=eta(2+eta), and the bilinear form differs by at most e n.
The type-weighted squared-norm potential has error at most 2lambda_G e n
over both sides. These estimates are uniform over every allowed declaration.

## 4. A type-dependent potential pays arbitrary P residuals

Set

    t=10/3,    lambda_G=35/18,    lambda_P=5/4,
    v=s^2=9/25,    a=1+2lambda_G v=12/5.

Use exp[t g^T K h - sum_rows lambda_type(||g_row||^2-k)],
with the analogous sum for h. Condition first on ALL actual P residuals,
their pair permutations and signs. Gaussian G coordinates are independent
over edge slots, so the following edge integrals multiply. No product law
for P coordinates is needed.

For ordinary GG edges the exact factor remains

    M_GG=(a^2-t^2 v^2)^(-1)=25/108.

On a GP edge, write the P vector as (u,0), with arbitrary real u. Rotate
the isotropic G pair by H2/sqrt(2) and integrate it exactly. The factor is

    a^(-1) exp[-(lambda_P-t^2 v/(2a))u^2]
       =(5/12)exp[-(5/12)u^2] <= 5/12.                (2)

On a PP edge the exponent is

    -(5/4)(u^2+w^2)+(t/sqrt(2))u w <= 0,              (3)

because 2lambda_P=5/2>10/(3sqrt(2))=t/sqrt(2).
Thus its factor is at most one, POINTWISE for arbitrary actual residuals.
This is the key reason no P spectrum needs to be reconstructed.

## 5. Center pins and exact norm compensation

Every pinned edge meets a reserved G row. For a pinned mean
mu=alpha sqrt(k) originating in a row of type T, its exponential factor
is exp[(t^2 v/(2a)-lambda_T)mu^2]. The coefficient t^2 v/(2a)=5/6
does not depend on T. A P-origin pin has determinant factor exactly 5/12,
the ordinary GP upper bound in (2). A G-origin pin differs from the
ordinary GG factor by sqrt(9/5). All pin determinant corrections cost
only exp[O(m)]. No pinned center-center term exists.

Let eps_L,eps_R be the fractions of P rows on the two sides. Combining
all norm compensations with all pinned means gives, per n,

    37/15 - (1/4)(eps_L+eps_R).                        (4)

Indeed changing one row's lambda from 35/18 to5/4 changes this sum by
(lambda_P-lambda_G)(1-alpha^2)k= -k/4. Replacing a P reference by an
arbitrary word did not change this accounting: its norm is fixed up to
the explicitly paid pinned-coordinate replacement error.

## 6. Entropy and the same-ensemble union

A G row with flip fraction delta_i has at most exp[k h(delta_i)] spins.
A P row has at most exp[(k/2)h(delta_i)], since its k/2 physical pairs
must agree. The pin error forces average |delta_i-1/10|<=tau_k/2.
Binary-entropy continuity and summing over row overlap values therefore
give the SAME upper count as in the mixed proof:

    log count <= n[2h(1/10)-(eps_L+eps_R)h(1/10)/2
                    +2h(tau_k/2)] +2m log(k+1).       (5)

From (2)--(4) the leading log moment per n is at most

    F0 + a1(eps_L+eps_R)+a2 eps_L eps_R,
    F0=37/15-(1/2)log(108/25),
    a1=(1/2)log(9/5)-1/4,
    a2=(1/2)log(4/3).                                 (6)

Since eps_L eps_R<=(eps_L+eps_R)/2, the Chernoff rate minus (5) at bridge
level18/25 is at least, before vanishing errors,

    D0 + D1(eps_L+eps_R),
    D0=(1/2)log(108/25)-2h(1/10)-1/15 >1/75,
    D1=h(1/10)/2+1/4-(1/2)log(9/5)-(1/4)log(4/3)>0.   (7)

Both are strict numerical inequalities to be accepted only via the
separate directed rational/interval certificate, not their decimal values.
The positive margin pays norm/transport error (t+2lambda_G)e n,
entropy continuity, O(m log k) overlap enumeration, O(m) pin and
conditioning factors, and the at most 2^(2m-4) G/P declarations.
All those costs are o(n). Both bridge polarities have the same estimates.
The resulting failure probability is below one, proving (1).

## 7. Quantitative implication and unpaid original obligation

For actual children of normalized cap at most c, in energy windows
|H_A(x)|,|H_D(y)|<=(16c/25+o(1))n^(3/2), (1) gives

    |H_A(x)+H_D(y)|+|x^T C y|
       <=[18/25+32c/25+o(1)]n^(3/2).                  (8)

At c>=.47 this is strictly below 2sqrt(2)c n^(3/2). This pays the
translation-paired exceptional spectra by an actual common bridge,
instead of falsely counting them as negligible or Gaussian. It also
explains why their large entropy is not itself a construction obstruction.

It leaves arbitrary non-paired exceptions, the four reserved G rows,
other energy windows, and other noise levels uncontrolled. No global
bound on Q(parent), power-saving flatification, or original convergence
is inferred. The exact remaining original obligation has not changed.
