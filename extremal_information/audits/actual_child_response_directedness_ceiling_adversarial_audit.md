# Adversarial audit: response-directedness ceiling

**Object audited:**
[`drafts/actual_child_response_directedness_ceiling.md`](../drafts/actual_child_response_directedness_ceiling.md)

**Verdict:** **PASS after making DC.9--DC.11 quantitative.**  The conditional
bit and conditional row-Renyi bounds, the reverse product-projection bound,
and the `Omega(r)` cross-row erasure-information conclusion are correct.
I patched the informal recovery argument into an explicit Pinsker/MMSE bound
and replaced a convergence statement with the exact martingale-variance
identity.  The construction remains a generic inference ceiling, not an
actual-child counterexample.

## 1. Conditional bit and row complexity

Given any observations outside a retained bit, the only shared uncertainty
is a posterior law on `W`.  Hence

```math
E[B_e\mid\text{rest}]=a_rE[W\mid\text{rest}],
```

whose magnitude is at most `a_r`.  This proves item 1 with no typical-event
qualification.

After arbitrary conditioning away from a retained set `S`, its law is still
a mixture of the three product components with means `0,+a_r,-a_r`.
Squared Hilbert norm is convex, and the component likelihoods satisfy

```math
\|f_0\|_2^2=1,
\qquad
\|f_\pm\|_2^2=(1+a_r^2)^{|S|}.
```

Therefore (DC.6), and hence
`D_2<=|S|log(1+a_r^2)<=c^2`, is valid even after conditioning changes the
three mixture weights.  The same proof also permits conditioning on some
coordinates of the row and marginalizing the unretained ones.

## 2. Reverse projection direction

The mixture contains `(1-epsilon)U_d` pointwise, so

```math
{dU_d\over dq_r}\le {1\over1-\epsilon}.
```

Thus `D(U_d||q_r)<=-log(1-epsilon)`.  Since `U_d` is a row-product law, it is
a feasible point for the reverse projection and (DC.4) has the correct KL
direction.  No comparison with forward total correlation is being used.

## 3. One-bit erasure information

Conditionally on `A,C`, the omitted sign has mean `a_rM_(AC)`; conditionally
only on `A`, it has mean `a_rM_A`.  For two sign laws, total variation is
half the difference of their means.  Pinsker therefore gives exactly

```math
D(P_{B_e|AC}\Vert P_{B_e|A})
\ge {a_r^2\over2}(M_{AC}-M_A)^2,
```

and averaging proves (DC.8).  The factor `1/2` is correct.

## 4. DC.9--DC.11 and the extensive sum

For any reference law `U_A`, the information-radius identity gives

```math
I(W;A)\le\sum_wP(w)D(P_{A|w}\Vert U_A).
```

Taking `U_A` fair and using conditional independence proves (DC.9).  Since
the one-bit divergence is `O(a_r^2)`, the result is `O(epsilon c^2)`
uniformly in `r`.

The original text invoked Pinsker/Fano informally.  The patch now uses

```math
E\operatorname{Var}(W|A)=\epsilon-E(E[W|A])^2
```

and bounds the second term by `sqrt(2I(W;A))`.  A sufficiently small fixed
`c` therefore leaves at least `epsilon/2` posterior variance.

There are `k=r(r-1)` observations in `C`.  Their conditional means are
`0,+a_rk,-a_rk`, and Hoeffding thresholding has error at most
`2exp(-a_r^2k/8)=exp(-Omega_c(r))`.  Hence the MMSE after `A,C` is `o(1)`.
The exact martingale identity is

```math
E(M_{AC}-M_A)^2
=E\operatorname{Var}(W|A)-E\operatorname{Var}(W|A,C),
```

so it is at least `epsilon/2-o(1)`.  Summing (DC.8) over the `r^2`
exchangeable bits gives

```math
r^2{c^2\over2r}(\epsilon/2-o(1))
=\left({c^2\epsilon\over4}-o(1)\right)r.
```

Thus one may take any fixed
`kappa_(c,epsilon)<c^2epsilon/4` for all sufficiently large `r`.  The claimed
`Omega(r)` bound is rigorous.

## 5. Scope

All four properties coexist: weak conditional bits, bounded conditional row
`D_2`, extensive cross-row erasure information, and `O(1)` reverse distance
to a row product.  Therefore no direction-free functional inequality can
turn the RR cross-row erasure sum into the reverse projection.  The example
does not have the optimized-child rank-one channel structure, so it only
forces any successful theorem to use an additional actual-channel
hypothesis.
