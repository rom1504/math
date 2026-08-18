# Repair verification: speed-`r` conference basin report

**Current frozen source:**
`extremal_information/drafts/conference_speed_r_basin_candidate_audit.md`

**SHA-256:**
`c0d1725e22c160bdd5e034e53c6c022072e2dda932f06338e067ddc181ece247`

**Verdict:** **PASS.**

Relative to the independently audited source at SHA-256
`688bcfd5303932728f6b8fc6ec1675dbb1e3c00ceece034fad2b6512e033d13a`,
the substantive repair is exactly the requested clarification at
(CB.24)--(CB.25): with the original mean profile `M` fixed, the source now
freezes

```math
W_M=(\sqrt{1-m_e^2}\,V_e)_e
```

and defines the translated convolution

```math
G_M(Z)=\mathbb E_V f(Z+W_M).
```

This function is convex in `Z`, and symmetry of `W_M` together with
`f(-B)=f(B)` makes it even.  Therefore

```math
\mathbb E f(M+W_M)=G_M(M)\ge G_M(0)=\mathbb E f(W_M),
```

which is exactly the comparison required in CB.3.  The variance profile no
longer appears to vary with the translation argument.  The subsequent
notation `W=(sqrt(1-m_e^2)V_e)` refers to the same frozen random vector, so
(CB.23) and (CB.26) remain consistent.

No further repair is required; the current hash is fully **PASS**.
