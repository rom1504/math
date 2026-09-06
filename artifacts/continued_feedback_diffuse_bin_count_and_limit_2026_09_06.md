# Diffuse finite-bin Walsh counting and a quantitative limit of this bound

Date: 2026-09-06. A scoped extension of the certified eight-block count.
The final obstruction concerns its explicit RELAXATION, not the actual
partition sum or the existence of actual rows with the trial profile.

Use p15/16, the uniform k-selector law, the Walsh basis, and the constant

```math
A=-65899554490727779/10485760000000000,
\qquad C_0=p\log2+A/8=-.13575842548\ldots.                (1)
```

These are from `continued_feedback_eight_block_stable_flat_certificate_2026_09_06.md`.
The Walsh basis is essential: this coding theorem cannot be combined
unchanged with a different randomly recursive Hadamard ensemble.

## 1. Finite-bin theorem, with its diffuse hypothesis explicit

Condition T on the high-probability event `||1_T-p||_{U²}=o(1)`. Declare
a class of supported spins xi satisfying `||xi||_{U²}=o(1)` uniformly.
For example a fixed upper bound on every normalized spectral magnitude
implies this condition, since the fourth Fourier moment is O(1/m).
Macroscopic Walsh spikes are NOT covered by this diffuse hypothesis.

Choose a fixed finite magnitude quantizer with centers c_j>=0, including
a distinguished zero center c_0=0, and deterministic bin assignments.
Let nu_j be a declared global bin profile. Suppose every row in the
class has normalized squared quantization error at most d²m:

```math
\sum_u(a_u-c_{\operatorname{bin}(a_u)})^2\le d^2m,
\qquad a_u=|(H\xi)_u|/\sqrt k.                          (2)
```

Put e=64pd² and assume e<17/18. If N_T(nu,d) counts such diffuse rows,
then

```math
\limsup_m\frac1m\log\mathbb E_T[1_{G_m}N_T(\nu,d)]
\le\min\left\{p\log2,
 C_0+\tfrac12[H(\nu)+(1-\nu_0)\log2]
       +\tfrac12[h(e)+e\log17]\right\}.                 (3)
```

The interpretation is literal: the zero bin costs no coefficient sign
bit, even when its actual coefficients are small and nonzero. Their
replacement by zero is already charged in (2). The bin profile need not
be a limiting continuous density.

### Proof

Use the same affine eight-block projection onto the four nonplane Walsh
frequencies. Diffuseness and selector quasirandomness give the same
four-wise product-P marginal law used for the fixed rational dual.
Thus the averaged conditional recovery cost is still
`m[H3+A/8+o(1)]`, independent of the declared bin profile.

Encode the selected m/2 coefficients by their bin and, outside the zero
bin, their sign. Averaging over the three block directions includes
each nonzero global frequency with probability 1/2+O(1/m); the exceptional
DC coordinate contributes only O(log m) to this finite-alphabet count.
Concavity of the selected-word entropy gives averaged coding cost at most

```math
\tfrac m2[H(\nu)+(1-\nu_0)\log2]+o(m).                 (4)
```

Equivalently, use the product symbol code with probabilities nu0 for the
zero bin and nu_j/2 for each signed nonzero bin. Its mean log length is
(4); zero-count bins are omitted. Because the number of bins is fixed,
the smallest positive empirical probability is at least 1/m, controlling
the exceptional DC payment. The direction is selected to minimize the
SUM of the coefficient code and the conditional recovery code, not each
separately; their averages bound that sum directly.

Replacing magnitudes by their centers gives total block-sum squared
error at most8d²k. Round the m/2 integer block sums and correct at most
32d²k of them. This is exactly the correction fraction e and entropy
payment in (3). The conditional recovery probabilities sum to one for
every fixed exact block-output word, without any assumption of actual
independent input blocks. Divide the total supported-vector count by
binomial(m,k), and cap it by the trivial 2^k count. This proves (3).

## 2. An explicit one-row profile upper expression

For actual finite bins I_j form the entrywise kernel upper matrix

```math
\overline K_{ij}=\sup_{a\in I_i,b\in I_j}K_t(a,b).
```

The director's finite table theorem bounds the square-root permanent
at rate `Psi_Kbar(nu)/2`, with only subexponential factors for a fixed
number of bins. A removed coordinate changes the profile by O(1/m),
and the finite positive-kernel pressure is continuous in the margins.
Therefore (3), plus the profile pressure and the weave tilt, is an
ACTUAL sufficient upper expression for each declared diffuse class:

```math
\min\{p\log2,
 C_0+\tfrac12[H(\nu)+(1-\nu_0)\log2]
          +\tfrac12[h(e)+e\log17]\}
 +\tfrac12\Psi_{\overline K}(\nu)+t(1-\sqrt p).          (5)
```

Summing over the polynomially many empirical profiles in fixed bins is
legitimate. Letting the bin mesh tend to zero is NOT free: entropy can
diverge, and the repair payment competes with finite-tilt kernel overlap.

## 3. A two-atom profile rigorously defeats this relaxation

Even with zero quantization error, the expression (5) cannot be made
uniformly negative over all diffuse profiles. Consider the formal trial

```math
\nu=\tfrac12\delta_a+\tfrac12\delta_b,
\qquad a=\tfrac12,\quad b=\tfrac{\sqrt7}2.               (6)
```

It has second moment one, no zero bin, and bounded magnitudes. Its code
rate from (3) is `C0+log2=.55738875508...`, below the trivial bound.
For equal masses on two types, the pressure is exactly

```math
\Psi_K(\nu)=\log\frac{\sqrt{K(a,a)K(b,b)}+K(a,b)}2.     (7)
```

To verify (7), a feasible coupling has diagonal masses r/2,r/2 and
off-diagonal masses (1-r)/2,(1-r)/2. Its variational objective is
`r log sqrt(Kaa Kbb)+(1-r)log Kab+h(r)-log2`; optimize r.

Put `gamma=(a-b)^2=2-sqrt(7)/2` and `Delta=1-sqrt(p)`. Since the two
diagonal entries are at least1/2 and
`K(a,b)>=exp(-gamma t)/2`, (7) implies that (5) at this profile is at least

```math
C_0+\tfrac12\log(1+e^{-\gamma t})+\Delta t.              (8)
```

This lower bound on the RELAXED UPPER EXPRESSION has minimum

```math
C_0+\tfrac12 h(2\Delta/\gamma).                         (9)
```

Indeed its derivative vanishes where the logistic probability is
2Delta/gamma; this lies between zero and1/2, so the minimizing t is
positive. Here `2Delta/gamma=(4-sqrt(15))/(4-sqrt(7))>.09` and is below
1/2. Also C0>-.136 and h(.09)/2>.151. These elementary inequalities can
be enclosed using the same rational logarithm series as the fixed dual.
Thus (8) is **greater than .015 for EVERY t>0**.

Direct non-certified minimization of the EXACT two-type expression
gives the stronger value about .02660434 at t3.65733; that number is not
needed for the rigorous conclusion.

No assertion is made that the actual ternary inverse-Walsh constraint
realizes the trial profile at the count permitted by (3). The point is
precise: the certified eight-block code plus its profile entropy and the
exact permanent variational bound admit this positive relaxation value,
so these inequalities alone do not prove the complete one-row bound.
A stronger profile-sensitive recovery count, or an additional actual
feasibility constraint, is still required. This quantitative limitation
survives arbitrarily fine binning because the trial has only two atoms.
