# Adversarial audit: the identical row-erasure reduction CR.0

**Object audited:** Proposition CR.0 in
[`../drafts/actual_child_cross_row_response_decomposition.md`](../drafts/actual_child_cross_row_response_decomposition.md).

**Verdict:** **PASS, including unequal child orders.**  The row-erased
forward likelihood is the same one-child mixture for every left row, and
the canonical inverse-escort product is therefore
`r_(row,u)^(tensor m)`.  I added the explicit conditioned sector law to make
the normalization and the unequal-order dependence checkable.

## Sector normalization

At zero bridge, condition on the relative orientation `epsilon` and write
`tau_1=s`, `tau_2=epsilon s`.  The augmented child law is

```math
 \Pr(s,x,y\mid\epsilon)
 \propto e^{stH_A(x)}e^{\epsilon stH_D(y)}.
```

Using the normalized one-sided partition functions from (CR.4a), its sector
probability is exactly

```math
 \pi_s^{(\epsilon)}
 ={Z_A^s(t)Z_D^{\epsilon s}(t)
   \over\sum_{a=\pm1}Z_A^a(t)Z_D^{\epsilon a}(t)}.
```

The factors `2^m` and `2^n` cancel between numerator and denominator.  No
extra factor two survives the conditioning on `epsilon`.

## Erased-row prior

Given `s`, the two children are independent.  Global-flip invariance of the
quadratic Hamiltonian makes `X_i` fair for every `i`, even when rows of `A`
are inequivalent.  Thus `sX_i` is a fair sign independent of `Y`.  The law
`mu_(D,epsilon s)` is itself globally flip invariant, so

```math
 Q_{i,*}=sX_iY\ \stackrel{d}{=}\ Y
 \quad\text{under }\mu_{D,\epsilon s}.
```

After mixing the sectors, every planted row has the same law
`sum_s pi_s^(epsilon) mu_(D,epsilon s)`.  Passing it through the independent
binary channel gives exactly (CR.4c).  Hence all row marginal likelihoods,
their inverse-power escorts, and their reverse-Renyi works agree.

## Unequal orders and scope

The proof uses only global flips and conditional product structure, not a
row permutation symmetry and not `m=n`.  For arbitrary `m,n`, the row word
has length `n`; `m` affects only `Z_A^s` and thus the two scalar mixture
weights.  The column statement follows analogously with the children
interchanged.  This is an exact finite identity for the actual optimized
children; it does not claim that the full bridge output or its inverse
escort is row-product.
