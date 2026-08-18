# Independent audit: expanded speed-`r` conference basin report

**Frozen source:**
`extremal_information/drafts/conference_speed_r_basin_candidate_audit.md`

**SHA-256:**
`688bcfd5303932728f6b8fc6ec1675dbb1e3c00ceece034fad2b6512e033d13a`

**Verdict:** **PASS WITH THE SAME NOTATIONAL REPAIR IN CB.3.**  The source
at this hash contains the previously audited CB.1--CB.4 plus the new CB.5
adaptive-gauge theorem.  All theorem statements and quantitative conclusions
are correct.  The only repair remains: freeze
`W_M=(sqrt(1-m_e^2)V_e)` and define `G_M(Z)=E f(Z+W_M)` before applying
convex evenness in (CB.24)--(CB.25).

The detailed audit of CB.1--CB.4 is in
`conference_speed_r_basin_candidate_independent_audit.md`, explicitly tied
to the earlier hash
`fba0a166c67337776aa7c2ecac42a00a3acc3830000992957ce26ab6fb8a7f47`.
Those sections are unchanged in the expanded source.  This artifact checks
the added theorem and the resulting expanded scope.

## CB.5 adaptive-gauge audit

1. **Group size and freeness.**  If `D_sBD_u=B`, then every nonzero bridge
   entry gives `s_i u_j=1`; hence all `s_i` and all `u_j` equal one common
   sign.  After quotienting the simultaneous global sign, the action is
   free and has `2^(2r-1)` elements.  The asserted transversal sizes follow.
   Row-only switching is free with group size `2^r`.

2. **Uniform norm event for every gauge.**  Row/column switching is left and
   right multiplication by orthogonal diagonal matrices, so it preserves
   `||B||_op`.  On
   `||B||_op<=(2+delta)sqrt(r)`, the triangle bound gives, simultaneously
   for all gauges,

   ```math
   \|tS_{\epsilon,D_sBD_u}\|_{op}
   \le {\beta\over\sqrt2}
      (\sqrt{1-1/r}+2+\delta)<\kappa
   ```

   for all large `r`, by the stated choice of `delta,kappa`.

3. **Fixed-gauge law.**  For each fixed `(s,u)`, `D_sBD_u` is again a
   uniform bridge and has the same operator norm as `B`.  Therefore the
   already proved operator-regular lower-tail theorem applies with exactly
   the same constants.  Equation (CB.42), obtained by the corresponding
   spin change of variables, is also correct, although distributional
   invariance of the bridge alone is enough for this step; no theorem for a
   new child pair is being silently assumed.

4. **Adaptive union.**  If an arbitrary bridge-dependent selector produces
   a low regular output, at least one of the at most `2^(2r-1)` fixed gauges
   produces it.  Thus the fixed-gauge `exp(-c r^2)` bound survives the union
   as `exp(-c_1 r^2)`.  Adding the standard rectangular norm-tail probability
   `exp(-c_0r)` proves (CB.40).

5. **Cross-section consequence.**  A deterministic orbit selector has one
   selected point per free orbit and the same number of preimages for each
   selected point.  Its pushforward of the uniform cube is therefore uniform
   on the cross-section.  Applying (CB.40) proves the claimed vanishing
   lower-deviation fraction for row, column, and row--column transversals.

## Scope

CB.5 closes adaptive gauge selection only on the operator-regular sector;
its surviving `exp(-O(r))` term is precisely the operator-irregular norm
tail.  The expanded report correctly does not interpret this term as a
favorable basin and does not claim a superexponential full lower-pressure
tail.  The final research judgment and table preserve that limitation.

## Required repair

Only the frozen-profile notation in CB.3 described above.  CB.5 requires no
repair.
