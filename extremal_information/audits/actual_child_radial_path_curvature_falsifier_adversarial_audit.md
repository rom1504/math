# Adversarial audit: actual-child radial path curvature falsifier

## Verdict

**PASS, with an essential scope qualification.**  I independently checked
the low-channel expansion RP.3, the exact ground-sector norms RP.10, and the
tropical normalization and rational slopes RP.13.  The saved result is
reproduced by the committed script, and the two order-eight inputs are actual
pressure-minimizing children in the temperature range used here, not
conference or Paley surrogates.

The qualification is that all diagonal statements in this note are
**fixed-order statements for the split `2+8`**.  Here “physical diagonal”
means that the internal raw temperature and bridge-channel amplitude are set
equal, `u=t`.  The limit in RP.13 is `t=u -> infinity` with parent order
`N=10` fixed, hence `beta=t sqrt(10) -> infinity`.  It is not the contracted
large-order regime `t=beta/sqrt(N)` at fixed `beta`.  Accordingly, the exact
collision disproves universal exact determination by radial child data, but
does not by itself disprove an asymptotic radial theorem that permits finite
exceptions or errors vanishing with order.

## Inputs and reproducibility

I audited:

- `drafts/actual_child_radial_path_curvature_falsifier.md`;
- `experiments/actual_child_radial_path_curvature_falsifier.py`;
- `computations/results/actual_child_radial_path_curvature_falsifier.json`;
- the cited canonical-path identities and the actual-child overlap
  decomposition.

I reran the experiment into
`/home/math/quadra/tmp/rp_audit.json`.  It reproduced the committed result,
including

```text
class=0 ground_trace=14 J/N=0.443573173 TCshare=0.816813379
class=1 ground_trace=10 J/N=0.523848938 TCshare=0.474757829
MATCH
```

The temporary output was then removed.  The ground-sector and tropical
parts use integer/rational arithmetic.  The `t=u=3` path integrals are
floating-point quadrature, not interval certificates.

## 1. Check of RP.3

### 1.1 Fourier order and coefficient

Let `U` be the fair bridge law and `rho=tanh(u)`.  The exact channel
likelihood is

```math
p_u(B)=E_{\mu_\epsilon}\prod_e(1+\rho B_eQ_e).
```

The central symmetry of the augmented child law kills every odd Walsh
level.  Thus, at fixed finite child orders,

```math
p_u(B)=1+\rho^2\sum_{e<f}E(Q_eQ_f)B_eB_f+O(\rho^4).
```

Subtracting the logarithms of the exact row marginals deletes precisely
the degree-two characters whose two edges lie in the same bridge row.  It
leaves

```math
h_u(B)=\rho^2 H_2(B)+O(\rho^4),
\qquad
H_2(B)=\sum_{i<k}\sum_{j,\ell}
\Gamma_{ik;j\ell}B_{ij}B_{k\ell}.
```

There is no missing factor of two: each cross-row unordered edge pair has a
unique representation with `i<k`, while `j,l` are freely ordered row
coordinates.  Distinct quadruples `(i,k,j,l)` give distinct Walsh
characters.

Since `r_u=U+O(rho^2)` on a fixed finite cube and the surviving Walsh
characters are orthonormal under `U`,

```math
\operatorname{Var}_{r_u}(h_u)
=\rho^4\sum_{i<k,j,\ell}\Gamma_{ik;j\ell}^2+O(\rho^6).
```

The canonical centered-cumulant formula integrates this variance with
weight of total mass `lambda^2/2`; hence

```math
\mathcal J_u=
{\lambda^2\rho^4\over2}
\sum_{i<k,j,\ell}\Gamma_{ik;j\ell}^2+O(\rho^6).
```

The coefficient, power of `rho`, and factor `1/2` in RP.3 are therefore
correct.

### 1.2 Total correlation versus row retuning

Every monomial of `H_2` uses one bit in each of two distinct rows.  Each
canonical reference row is central.  Conditional expectation over all
rows except any fixed row therefore annihilates `H_2`.  The first
perturbation of a one-row marginal is consequently `O(s rho^4)`, and its KL
cost is `O(s^2 rho^8)`.  After the RP.2 division by `s^2` and integration,

```math
R_u=O(\rho^8).
```

Using the exact ES decomposition `J_u=T_u+R_u` yields the asserted leading
term for `T_u`.  This portion of RP.3 is correct.

### 1.3 Uniformity caveat

All `O(.)` bounds above are on a fixed finite bridge cube, with child
orders, internal raw temperature `t`, and `lambda` fixed.  Their constants
may depend on those parameters.  RP.3 therefore supplies no uniform
large-order or contracted-temperature estimate.  The draft states the
fixed-order hypothesis correctly; later uses must retain it.

Also, the algebra of RP.3 itself applies to any central rank-one binary
channel.  Its actual-minimizer content enters through the certified pair in
Section 2, rather than through the local expansion alone.

## 2. Check of RP.9 and RP.10

For the order-two left child, direct sector conditioning gives

```math
E[X_1X_2\mid\tau=s]=s\tanh(t).
```

Substitution into the augmented sector law gives exactly

```math
\Gamma_{12;j\ell}^{(\epsilon)}
=\epsilon\tanh(t)S_D(t)_{j\ell},
```

so RP.9 has the correct sign and normalization.  The sign disappears from
the squared curvature coefficient.

For each order-eight class, there are eight active projective ground states
with `|H|=10`.  If `G_D` denotes the stored signed integer sum

```math
G_D=\sum_{|H_D(y)|=10}\operatorname{sign}(H_D(y))yy^{\mathsf T},
```

then the zero-temperature tangent is `S_D(infinity)=G_D/8`.  The committed
integer matrices satisfy

```math
\|G_{A_0}\|_F^2=896,
\qquad
\|G_{A_1}\|_F^2=640.
```

They are symmetric, so `Tr(S_D^2)=||S_D||_F^2`.  Division by `8^2=64`
gives respectively `14` and `10`, exactly as in RP.10.  Combining with
RP.3 and `tanh(t)->1` gives the two low-channel coefficients
`7 lambda^2` and `5 lambda^2`.

The equal signed energy histograms were also reproduced.  Together with the
cited exhaustive order-eight classification and pressure comparison, this
pair consists of actual contracted-temperature pressure minimizers for the
stated raw-temperature range `t>=3`.  It is not a surrogate pair.

## 3. Check of RP.13

### 3.1 Correct tropical parent rate

For orientation `epsilon=-1`, the exact projective parent summand is

```math
\cosh(t(H_A(x)-H_D(y)))\cosh(t x^{\mathsf T}By).
```

Its exponential rate is therefore

```math
|H_A(x)-H_D(y)|+|x^{\mathsf T}By|,
```

not `|H_A-H_D+x^TBy|`.  Maximizing this expression over projective child
spins gives the `kappa_D(B)` used in RP.12 and in the script.

### 3.2 Leading coefficients and escort limits

For every bridge in this particular certificate, the script checks that
every maximizing summand has both cosh arguments nonzero.  Consequently
each active summand contributes the same factor `1/4` at leading order, and
the remaining leading coefficient is proportional to the exact active
multiplicity.  All averaging and `1/4` constants common to a table cancel
in the normalized escort laws.

For a fixed row word `b_i`, summing over the other bridge row is dominated
by

```math
\kappa_i(b_i)=\max_{B_{-i}}\kappa_D(B).
```

The negative escort at `lambda=1` first selects row words minimizing this
row-minimax rate and then weights them by the inverse of the summed active
coefficient among its maximizing completions.  The endpoint escort instead
selects global minimizers of `kappa_D`.  Thus

```math
\lim_{t\to\infty}{\mathcal J_t\over t}
=E_{r_\infty}\kappa_D-min_B\kappa_D(B).
```

This normalization is correct for `lambda=1`, the specialization explicitly
made in RP.13.  For general `lambda`, the escort exponents and resulting
weights would change; the draft does not claim the displayed fractions for
general `lambda`.

### 3.3 Exact rational values

The `2 x 8` bridge cube has `2^16` elements and is completely enumerated.
The exact integer tables give:

```math
\begin{array}{c|c|c|c|c}
 D&\min\kappa_D&\min\kappa_i&|\operatorname{supp}r_{\infty,i}|
 &E_{r_\infty}\kappa_D\\ \hline
 A_0&15&21&24&205/12\\
 A_1&15&21&8&84/5.
\end{array}
```

Therefore

```math
205/12-15=25/12,
\qquad
84/5-15=9/5.
```

RP.13 is exact.  I also checked the bridge-bit reshape convention: the low
eight bits form one row and the high eight bits the other.  The two axis
labels are used consistently in the row-support calculation, and in any
event the final product expectation is symmetric under exchanging the two
rows.

## 4. Scope of the falsifier

What is rigorously falsified is an **exact universal determination** of the
canonical interaction path by complete radial pressure/entropy data:

- at fixed finite order, the low-channel curvature differs exactly;
- at fixed order `N=10`, the same-amplitude zero-temperature diagonal slope
  differs exactly;
- at `t=u=3`, complete finite enumeration gives a large numerical
  separation in both `J/N` and its TC/retuning allocation.

What is **not** established is:

- a separation along `N->infinity` with `t=u=beta/sqrt(N)` and fixed
  `beta`;
- a positive asymptotic density of either interaction quantity;
- failure of a radial statistic supplemented by a genuinely nonradial,
  compressed overlap state;
- failure of an approximate radial theorem whose allowed error absorbs
  fixed-order collisions.

Accordingly, phrases such as “physical-temperature ray” should be read, or
ideally qualified, as **the fixed-order same-amplitude/zero-temperature
ray**.  The finite `t=u=3` computation is likewise at `N=10` and
`beta=3 sqrt(10)`, not in the contracted asymptotic regime.

## Final audit judgment

The theorem and certificate are mathematically sound within their stated
finite-order scope.  They give a legitimate actual-minimizer no-go for
radial-only summaries, with no hidden conference surrogate, scalar-channel
factor, or tropical normalization error.  They should not be cited as an
asymptotic obstruction to contracted-temperature child control without a
new growing-order construction or a uniform version of RP.3/RP.13.
