# Main-agent Wave 40 note: cyclic monotonicity of canonical complement grounds

This is an independent derivation for later comparison with the assigned
agent's work.  It uses the ordered matrix pairing and the notation of
(10.1043)--(10.1045).

For an `m`-selector `S`, let `B_S=A^{F_S}` be obtained by flipping every edge
outside `E(S)`.  For an oriented cut word `d`,

```math
H_{B_S}(d)=2c_S(d)-E_d.                                      (M1)
```

Choose `d_S` to maximize

```math
H_{B_S}(d)-P(d)                                               (M2)
```

over all oriented cut words, where `P` is **any selector-independent
penalty** (row price, Hamming price to a fixed center, or another common
regularizer).  For two selectors, optimality gives

```math
H_{B_S}(d_S)-P(d_S) >= H_{B_S}(d_T)-P(d_T),
H_{B_T}(d_T)-P(d_T) >= H_{B_T}(d_S)-P(d_S).
```

The penalty and the parent-energy terms cancel after addition.  Hence

```math
\boxed{
c_S(d_S)-c_S(d_T)+c_T(d_T)-c_T(d_S)\ge0.
}                                                               (M3)
```

Equivalently,

```math
\langle B_S-B_T,d_S-d_T\rangle\ge0.                            (M4)
```

More generally, for a selector cycle `S_1,...,S_k,S_(k+1)=S_1`, comparison
of `d_(S_j)` with `d_(S_(j+1))` and summation gives

```math
\boxed{
\sum_{j=1}^k c_{S_j}(d_{S_j}-d_{S_{j+1}})\ge0.
}                                                               (M5)
```

Thus the canonical complement-ground map is cyclically monotone even after
an arbitrary common penalty.  If `d_S` has optimization error `zeta_S` in
(M2), the right side of (M3) becomes at worst
`-(zeta_S+zeta_T)/2`.

This does not yet control anchored conflict.  When two projective labels have
the same auxiliary orientation and differ by a vertex shore `U`, (M4) is only
a **signed** sum over the intersection of `cut(U)` with
`E(S) triangle E(T)`.  It can vanish for macroscopic `U`; no modulus of
monotonicity follows from (M3).  Low row gives an upper bound on such a shore,
not the lower bound in distance that conflict control would require.  The
Wave 39 complete-signing block construction exhibits exactly the relevant
zero-shore degeneracy for arbitrary grounds, though not yet for canonical
complement grounds of an actual minimizer.

The precise positive successor would be a minimizer-specific **normal-fan
curvature** statement: for a common penalty that preserves the certificate
threshold and row budget, show that the Bregman gap in (M2) dominates a
project-scale function of anchored Hamming disagreement after averaging over
an affordable selector family.  Without such a gap, cyclic monotonicity is a
new exact organization of the exchange information but not a conflict bound.
