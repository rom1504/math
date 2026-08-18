# Independent audit: rowwise antipodal halfcube pressure

**Frozen source:**
`extremal_information/drafts/conference_row_halfcube_pressure.md`

**SHA-256:**
`d398b66308c63cdc2c4b00850c91ca769aceb00b66b29b5aec6d9a7a5b172c78`

**Verdict:** **PASS.**  The exact count and pushforward, uniform switched-
child concentration, orbit union, two-sided and `L^1` conclusions, and
scope are all correct.  No repair is required.

## Exact family and pushforward

The odd tie rule selects exactly one row from every antipodal pair, including
when `r` is even and zero inner products occur.  Hence every row halfcube
has `2^(r-1)` elements and the product family has exactly `2^(r^2-r)`
elements.

For every input row `W_i` there is a unique sign `s_i(W)` for which
`s_i(W)W_i` is its selected representative.  Conversely, every selected
bridge has exactly `2^r` preimages, one for each independent input-row sign.
Thus `D_(s(W))W` is exactly uniform on the halfcube, not merely asymptotically
so.

## Uniform switching concentration

The spin substitution `x -> D_sx` proves (RH.12).  Equivalently,
conjugating by `diag(D_s,I)` transforms the fixed-child parent with bridge
`D_sW` into the switched-left-child parent with bridge `W`; operator norms
are identical.

If `||W||_op<=(2+delta)sqrt(r)`, block-triangle inequality gives

```math
\left\|{\beta\over\sqrt{2r}}
\begin{pmatrix}D_sAD_s&W\\W^T&\epsilon A\end{pmatrix}\right\|_{op}
\le {\beta\over\sqrt2}
   (\sqrt{1-1/r}+2+\delta)<\kappa
```

simultaneously for every `s`.  This verifies (RH.15).

For fixed `s`, the map `W -> D_sW` preserves both the uniform bridge law and
the operator-norm event.  Via (RH.12), the pressure therefore has exactly the
same law and the same mean error as the unswitched conference pressure.
There is no nonuniform `o(r)` term hidden among the `2^r` switched children.

On the regular set, the dimension-free convex Frobenius-Lipschitz extension
has two-sided Talagrand concentration.  Its mean differs from the conference
center by `o(r)`, uniformly as just observed, so (RH.16) follows for each
fixed `s`.  Union over `2^r` switches changes `exp(-c r^2)` only by
`exp(O(r))`; adding the rectangular norm tail proves (RH.17) with the stated
`exp(-c_0r)+exp(-c_1r^2)` form.

## Adaptive selector, two-sided convergence, and `L^1`

The event in (RH.17) is simultaneous in `s`, so it remains valid for the
fully bridge-dependent selector `s(W)`.  The exact pushforward then proves
the two-sided quantitative halfcube estimate (RH.19) and convergence in
probability.

Every sign parent satisfies `0<=f<=C_beta r^(3/2)`.  After division by `r`,
the exceptional contribution is at most
`O(sqrt(r))(exp(-c_0r)+exp(-c_1r^2))=o(1)`.  On the good event the normalized
error is arbitrarily small.  This proves uniform integrability and the
claimed `L^1` convergence.  The fixed positive gap `h_beta-tau_beta` then
implies the vanishing target-reaching fraction in (RH.20).

## Scope

The theorem closes the declared antipodal selector family, including
arbitrary deterministic `u,v` and ties.  It does not assert that the
operator-irregular exceptional subfamily is absent or superexponentially
small relative to the full bridge cube.  The source explicitly preserves
this limitation and correctly distinguishes the switching mechanism from a
small-Hamming repair.

## Corrections

None required.
