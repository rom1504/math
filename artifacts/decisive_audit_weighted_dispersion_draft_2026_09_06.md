# Frozen weighted-dispersion screen

2026-09-06. Negative finite diagnostic; this branch is frozen at the
director's request in favor of original-limit compatibility.

For positive integer weights define
`Q_w(B)=max_x |sum_(i<j) B_ij w_i w_j x_i x_j|`.
The tempting finite inequality

`min_B Q_w(B)/(sum_i w_i sqrt(sum_i w_i²)) >= M_d/d^(3/2)`

is false. At `d=7`, `w=(2,2,2,2,2,1,1)`, exhaustive switching-gauged
enumeration gives minimum 21, achieved by free-edge code 956 in the
checker convention. Its normalized value is
`21/(12 sqrt(22))=.3731012536223183`, below
`M_7/7^(3/2)=9/7^(3/2)=.485954322440435`.

`computations/decisive_audit_weighted_dispersion_screen_2026_09_06.py`
enumerates all switching-gauged signings through `d=7` and all patterns
of weights in `{1,2}`. It prints the complete exact-cap table and witness
codes. It imports no research data and writes no files.

This does not falsify an asymptotic inequality involving the global
liminf constant, and does not classify partially repeated rows of
near-minimizers. For uniform fixed multiplicity `r`, the elementary
blowup identity gives `Q(B tensor J_r)=r²Q(B)+O(n)` after bounded
within-class diagonal choices, so such uniformly cloned families
cannot be liminf-minimizing: their normalized cap is at least
`sqrt(r) liminf M_d/d^(3/2)+o(1)`. Extending this to nonuniform clone
weights would need a genuinely new weighted-cap argument.
