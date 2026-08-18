# Convex/superconcentration audit for the actual child collision free energy

**Question.**  Can standard convex concentration or superconcentration
machinery control

```math
h(B)=\log G_{A,D,u}(B)
```

under the hybrid actual-child escorts `q_s`, using the fact that `A,D` are
pressure minimizers?

**Verdict.**  Not from the optimizer information currently proved.  There
is a precise one-row curvature theorem to which Poincare, logarithmic
Sobolev, or superconcentration estimates could be applied.  The universal
actual-law inputs prove only an `O(N)` bound there, and that scale is sharp
for both the rank-one and block-parity falsifiers.  The missing
optimizer-specific input is an external-field cavity superconcentration
theorem; AC.32--AC.33 are zero-field, one-sided lower-moment constraints and
do not provide it.

The falsifiers below are not claimed to be optimizing children.  They rule
out every attempt that stops before using optimizer-specific structure.

## 1. Exact local-curvature reduction

Fix an actual child pair and the notation of IC.1--IC.4.  For a row `i`, a
configuration `B_-i`, and `0<=theta<=s`, define the conditional local tilt

```math
{d\mu_{i,\theta}^{B_{-i}}\over dr_i}(b)
={e^{-\theta h(b,B_{-i})}
  \over E_{r_i}e^{-\theta h(R_i,B_{-i})}}.           \tag{CSA.1}
```

At `theta=s` this is exactly `q_s^{i|-i}`.  Conditional exponential-family
calculus gives

```math
\boxed{
D(q_s^{i\mid-i}\Vert r_i)
=\int_0^s\theta\,
 \operatorname{Var}_{\mu_{i,\theta}^{B_{-i}}}
       h(R_i,B_{-i})\,d\theta.}                     \tag{CSA.2}
```

Indeed the derivative of the conditional KL at `theta` is `theta` times
the displayed variance, and the KL vanishes at zero.

Define

```math
\Xi_N
=\sup_{0<s\le\lambda}\sum_{i=1}^m
 E_{(q_s)_{-i}}
 \sup_{0\le\theta\le s}
 \operatorname{Var}_{\mu_{i,\theta}^{B_{-i}}}
       h(R_i,B_{-i}).                                \tag{CSA.3}
```

**Theorem CSA.1 (one-row superconcentration closure).**  For every finite
actual-child channel,

```math
\boxed{\mathcal J\le {\lambda^2\over2}\Xi_N.}       \tag{CSA.4}
```

Consequently `Xi_N=o(N)` closes the canonical no-gain branch, while
`J>=eta N` forces `Xi_N>=2 eta N/lambda^2`.

*Proof.*  Substitute (CSA.2) into IC.21.  The inner integral is at most
`s^2/2` times its supremum, hence `E_s<=Xi_N/2`.  IC.23 gives
`J<=lambda int_0^lambda E_s ds<=lambda^2 Xi_N/2`.
`square`

This is the minimal place where ordinary one-row concentration machinery
can enter without paying the full and row-erased channels separately.

## 2. What a functional inequality would actually prove

For a row bit `e=(i,j)`, let `B^e` denote its flip and put

```math
\nabla_e h(B)=h(B)-h(B^e).
```

Suppose every conditional tilt in (CSA.1) satisfies, with a common constant
`C_P`, the flip-form Poincare inequality

```math
\operatorname{Var}_{\mu}(f)
\le C_P\sum_{j=1}^nE_\mu
       [f(B)-f(B^{ij})]^2.                           \tag{CSA.5}
```

Define the actual conditional cavity gradient energy

```math
\Gamma_N
=\sup_{0<s\le\lambda}\sum_iE_{(q_s)_{-i}}
 \sup_{0\le\theta\le s}\sum_j
 E_{\mu_{i,\theta}^{B_{-i}}}(\nabla_{ij}h)^2.       \tag{CSA.6}
```

Then CSA.1 immediately gives the rigorous sufficient theorem

```math
\boxed{\mathcal J\le {\lambda^2C_P\over2}\Gamma_N.} \tag{CSA.7}
```

Thus a uniform functional inequality together with `Gamma_N=o(N)` would
close the phase.

The available deterministic estimates stop exactly one order too early.
Flipping one bridge bit changes `log p` and its row marginal `log p_i` by at
most `2u` each, so

```math
|\nabla_{ij}h|\le4u,
\qquad
\Gamma_N\le16u^2mn.                                 \tag{CSA.8}
```

At `u=beta/sqrt(N)` and a comparable split, (CSA.7)--(CSA.8) give only
`J=O(N)`.  Even a dimension-free Poincare or logarithmic-Sobolev constant
does not produce the required `o(N)` result.  A genuine superconcentration
gain must prove a power saving in `Gamma_N`, or improve (CSA.5) on the
specific cavity functions by a factor tending to zero.

This is not an artifact of constants.  Under the uniform law on `n` bits,
`f(B)=n^{-1/2}sum_jB_j` has bit oscillation `2/sqrt(n)`, zero density
complexity, and variance one.  Small coordinate oscillation and bounded
conditional Renyi complexity therefore permit order-one variance per row.

## 3. Convex concentration does not apply jointly

Embed bridge signs in real fields `b` and write

```math
F(b)=\log E_Qe^{u\langle b,Q\rangle},
\qquad
F_i(b_i)=\log E_{Q_i}e^{u\langle b_i,Q_i\rangle}.
```

Up to constants that cancel,

```math
h(b)=F(b)-\sum_iF_i(b_i).                            \tag{CSA.9}
```

Each term is convex, but their difference need not be.  Its Hessian is

```math
\nabla^2h
=u^2\left[
 \operatorname{Cov}(Q\mid b)
 -\bigoplus_i\operatorname{Cov}(Q_i\mid b_i)
 \right],                                           \tag{CSA.10}
```

which has no fixed sign.  Hence convex Lipschitz concentration,
Brascamp--Lieb, and related variance inequalities cannot be applied to `h`
as a single cancellation-preserving observable.  Applying them separately
to `F` and `sum_iF_i` discards the exact cancellation and returns the same
leading `O(N)` scale.

The rank-one falsifier makes the sign failure explicit.  With the notation
of CC.19,

```math
h=g\left(\sum_iV_i\right)-\sum_ig(V_i),
\qquad g(v)=\log\cosh v.                            \tag{CSA.11}
```

For one row and fixed `S=sum_(k ne i)V_k`, the second derivative is

```math
g''(S+v)-g''(v).
```

For `S>0` it is negative at `v=0` and positive at `v=-S`; the cavity
function is neither convex nor concave.

## 4. The falsifiers hit the local-curvature criterion sharply

### 4.1 Rank-one channel

Under the canonical product at `s=0`, the row summaries `V_i` in CC.20 are
iid, nondegenerate, and converge to a fixed tilted Gaussian law.  Given
`S=sum_(k ne i)V_k`, the row cavity function is, up to a constant,

```math
v\longmapsto g(S+v)-g(v).                            \tag{CSA.12}
```

By the central limit theorem, `|S|` tends to infinity in probability.  On
bounded `v`-sets, (CSA.12), after subtracting a constant, converges according
to the sign of `S` to `v-g(v)` or `-v-g(v)`.  Both have strictly positive
variance under the nondegenerate limiting row law.  Truncation and the
subgaussian row tails therefore give

```math
\sum_iE_{r_{-i}}\operatorname{Var}_{r_i}(h\mid B_{-i})
=\Theta(m)=\Theta(N).                               \tag{CSA.13}
```

Thus `Xi_N=Omega(N)` already at the zero-tilt limit.  This agrees with the
linear `J` in CC.22--CC.23 and shows that no generic convex or Poincare
estimate can improve the scale.

### 4.2 Block-parity channel

For CC.25--CC.29,

```math
h=\sum_{I}\log\left(1+\delta\prod_{i\in I}a_n(B_i)\right),
                                                               \tag{CSA.14}
```

with independent fixed-size blocks under `r`.  Conditional on the other
rows in its block, one row sees

```math
b\longmapsto\log(1+c_Ia_n(b)).                       \tag{CSA.15}
```

The tilted row CLT makes `a_n(B_i)` converge to a nondegenerate variable;
`c_I` is nonzero almost surely in the limit.  Hence the averaged conditional
variance in (CSA.15) has a strictly positive limit.  Summing over rows gives

```math
\sum_iE_{r_{-i}}\operatorname{Var}_{r_i}(h\mid B_{-i})
=\Theta(m)=\Theta(N).                               \tag{CSA.16}
```

This example additionally has all left overlaps through any prescribed
fixed order equal to the iid values.  Finite-order overlap control cannot
supply the missing power saving in (CSA.6).

## 5. What pressure minimality supplies, and why it is the wrong sign

For an actual minimizing child, the genuine extra input is AC.32:

```math
E_{\nu_A}\exp\left{-2t\tau\sum_{e\in S}a_ex_e\right}\ge1
\quad\text{for every internal edge set }S.          \tag{CSA.17}
```

These are **lower** bounds on exponential moments at zero external field,
the opposite direction from a subgaussian or superconcentration estimate.
Viewed only through this lower-moment form, variance and higher fluctuations
make the inequality easier to satisfy; extracting the upper conditional
variances in (CSA.3) or (CSA.6) would itself require a new theorem using the
joint Gibbs/signing structure.  The one-edge consequence AC.33 controls only
the signed first moment

```math
a_eE_{\nu_A}(\tau x_e)\le\tanh t.                   \tag{CSA.18}
```

It gives no square or covariance upper bound.

There is also a mismatch of environments.  The cavity formula CC.8 queries
the child at external fields

```math
u s B^{\mathsf T}X
```

and then averages under inverse bridge tilts `q_s`.  AC.32 concerns the
unperturbed child Gibbs law and internal sign flips.  A zero-field pressure
minimizer need not remain minimizing, stable, or weakly correlated under
these external fields.  No proved implication transfers (CSA.17) to a
Poincare constant, convexity, derivative intermittency, or an `o(N)` bound
on the cavity gradient energy.

The full inhomogeneous contraction FC.8 does not repair this gap.  It is
equivalent to all sign-flip comparisons, hence to the complete discrete
minimization statement; using all its boundary data without an additional
compression reconstructs the child optimization.  Its homogeneous/radial
part retains only the absolute-energy histogram, already proved
insufficient for extension response by the actual order-eight minimizer
collision in Section 5 of the flip-averaging note.

## 6. Sharp missing optimizer-specific theorem

Known machinery would apply if one proved the following new statement.

> **Actual-minimizer external-field cavity superconcentration.**  For
> contracted-temperature minimizing children and fixed `beta,lambda`, the
> local tilts (CSA.1), under the actual hybrid outer laws, satisfy
> `Xi_N=o(N)` in (CSA.3).  A more structural sufficient version is a uniform
> row functional inequality (CSA.5) together with
> `Gamma_N=o(N)` in (CSA.6).

This is strictly about the actual child-induced law and never invokes a
target-order optimizer.  It is also exactly the optimizer-specific input
missing from current convex/superconcentration tools: stability under the
continuum of cavity fields, with a power saving over the universal
`O(N)` gradient budget.

The theorem is not implied by weak coordinates, bounded conditional row
Renyi complexity, rank-one support, any fixed overlap hierarchy, or
convexity of the constituent log partitions.  The two scalable falsifiers
prove those ceilings.  Conversely, the falsifiers do not disprove the
theorem for actual pressure minimizers.  They show that a proof must extract
new external-field stability from minimization, rather than invoke a known
concentration inequality on the inputs already available.

**Director conclusion.**  There is presently no applicable off-the-shelf
convex/superconcentration closure.  The best rigorous output is CSA.1 and
the ceiling (CSA.8): existing machinery localizes the issue to one-row
cavity curvature but cannot improve its leading scale.  The next step would
need a genuinely new optimizer-to-external-field stability lemma; absent
that, further generic concentration work is a strike rather than a reset.
