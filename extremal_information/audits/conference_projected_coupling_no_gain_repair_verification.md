# Repair verification: projected-coupling no-gain criterion

**Repaired source:**
`extremal_information/drafts/conference_projected_coupling_no_gain.md`

**SHA-256:**
`5df7d16d7a144c255f0f4ab520d732ebadf04206bc6a38b1081850c054e1d597`

**Verdict:** **PASS.**  The repair implements the stronger direct one-sided
proof correctly.  PC.8 now requires only `Pr(G_r^c)=o(1)`, PC.4--PC.7, and
the archived iid conference pressure theorem.  No exceptional-event uniform
integrability or sign constraint is needed.

On `G_r`, PC.9--PC.11 give

```math
f_r(B_r)\ge f_r(W_r)
-{K_\kappa\beta\over\sqrt2}
 \|(B_r-W_r)(I-P_r)\|_F
-{K_\kappa\beta\over\sqrt{2r}}
 (\|W_rP_r\|_*+\|B_rP_r\|_*).
```

After division by `r`, PC.6 makes the first error `o(1)` in mean and PC.7
makes the other two errors `o(1)` in mean.  On `G_r^c`, nonnegativity of the
cosh pressure bounds the positive normalized shortfall by `h_beta`.  This
is exactly sufficient for PC.8 and leaves no hidden integrability step.

The nuclear constants and the rank corollary are unchanged and remain
correct.

## Editorial debris from the repair

Three nonmathematical cleanup edits are advisable:

1. “Assume the following four conditions” should now say “three
   conditions.”
2. The setup sentence saying exact signs are needed for a crude
   exceptional-event pressure bound is stale; exact signs are no longer
   needed anywhere in PC.1.
3. The proof ends with two consecutive `` `square` `` markers; delete one.

These do not affect the theorem or proof.

## Final editorial verification

The final source with SHA-256
`d02a34238d5339437d46cdc0fbed838f335cf7660167c1450d179504cb3cc11e`
removes all three items above.  I rechecked the resulting statement and
direct proof.  **Final verdict: PASS.**
