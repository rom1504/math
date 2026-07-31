# A non-Kronecker correction kills old witnesses but worsens fresh ones

Date: 2026-07-31. This is an agent-authored research report. It tests a
bounded-description cross-fiber correction after the complete diagonal audit.

## 1. Exact correction model

For each of the `2^14=16384` diagonal completions `D`, let `z_D` be the saved
order-56 witness from the all-diagonal certificate, with positive energy
`E_D` between 210 and 240. Index a vertex as `(i,a)`, where `i` is one of 14
macro vertices and `a` is one of four micro positions.

For a fixed ordered micro position `(a,b)`, define one **bundle** to toggle
the 91 edges

```math
\{\{(i,a),(j,b)\}:0\le i<j<14\}.                    \tag{MC1}
```

This changes the cross-macro micro template itself. It is not a diagonal
term of the form `D tensor Delta` and is therefore a materially different
representation from the diagonal-completion audit.

Let `y_(a,b)` be the Boolean decision to toggle bundle `(a,b)`, and put

```math
c_{D,a,b}=\sum_{i<j}A_{ij}(H_4)_{ab}
              (z_D)_{i,a}(z_D)_{j,b}.                \tag{MC2}
```

The saved witness's corrected energy is exactly

```math
E_D'=E_D-2\sum_{a,b}c_{D,a,b}y_{a,b}.                \tag{MC3}
```

The finite certificate-hitting problem is the 16-variable binary program

```math
\min\sum_{a,b}y_{a,b}
\quad\text{subject to}\quad
-208\le E_D'\le208\quad\text{for every }D.           \tag{MC4}
```

Both sides of the absolute-value constraint are explicitly imposed.

## 2. Exact optimum for the saved witness set

The exact optimum of (MC4) is two. One optimum is

```math
y_{0,1}=y_{1,0}=1,\qquad y_{a,b}=0\text{ otherwise}. \tag{MC5}
```

It modifies 182 of the 1,456 cross-fiber edges, exactly `1/8=12.5%`.
For all 16,384 old witnesses, direct exact evaluation after (MC5) gives

```math
122\le E_D'\le176.                                   \tag{MC6}
```

The lower bound does not rely only on the optimizer's status. The independent
verifier exhausts the empty choice and all 16 singleton bundles; each leaves
some witness with absolute energy above 208. It then checks (MC5) against all
16,384 rows. Hence the optimum two is an independently checkable finite
theorem.

The corrected-energy vector hash is

```text
f05b6e72c00ac8efcb8a40bca7bc68d7c71bc758ad92cab8c7d65b6012f10ac9
```

The model and independent verifier are:

- `computations/phase2f_cross_fiber_certificate_hitting.py`;
- `computations/phase2f_verify_cross_fiber_micro_hitting.py`;
- `computations/results/phase2f_cross_fiber_micro.json`;
- `computations/results/phase2f_cross_fiber_micro_verify.log`.

## 3. Fresh-witness falsification

Solving (MC4) controls one previously saved witness per signing. It does
**not** upper-bound any corrected signing's cap. To test whether (MC5) was a
real structural repair, a fresh two-sided coordinate-ascent search was run
independently on every one of the 16,384 corrected signings. There were no
holdouts: every corrected signing produced a new explicit witness of absolute
energy at least 230.

All new witnesses were then independently checked by direct integer matrix
multiplication. The exact finite-family statement is

```math
\boxed{
 \min_D\;\max_{\text{saved fresh }w_D}
 |H_{S_D^{\mathrm{corr}}}(w_D)|=230.}                \tag{MC7}
```

Here the minimum means the minimum over the 16,384 explicitly saved lower
bounds, not an exact computation of the individual caps. The distribution is

```text
230: 2, 234: 74, 236: 10, 238: 5646,
240: 204, 242: 10448.
```

Thus the correction lowers all old certificates into `[122,176]` yet every
corrected signing immediately exposes a new certificate in `[230,242]`. The
minimum certified cap lower bound actually rises from 210 to 230.

The fresh certificate and verifier log are:

- `computations/results/phase2f_corrected_micro_pair_all_diagonal_audit.json`;
- `computations/results/phase2f_corrected_micro_pair_all_diagonal_audit_verify.log`.

Their canonical record hash is

```text
f24a16404257c693c1a82f70bf290fee00e276b2d95d8557e3253c1eeb607b89
```

## 4. Scope and resulting invariant requirement

The positive part is concrete: a bounded non-Kronecker micro-template change
can simultaneously destroy all 16,384 witnesses selected under the old
template, and two bundles are exactly necessary within this representation.

The negative part is more important for route selection. A state consisting
of one active witness per diagonal is not stable under correction. The
optimizer can move that finite active face far below threshold while another
face becomes uniformly worse. A useful correction state must therefore
control an **envelope** of near-ground states, or a dual quantity that is
stable when the micro template changes. Individual witness separation is not
a cap-control invariant.

This report proves no tensor persistence for the corrected template and no
cap upper bound. It also does not rule out a richer scale-dependent micro
algebra. It supplies a falsifiable next question: whether alternating exact
template correction and fresh-witness generation stabilizes after a bounded
number of rounds, or instead keeps exposing new faces until the active state
has unbounded complexity.

## 5. Bounded alternating-game diagnostic

One bounded continuation was performed, and then stopped rather than turning
into an indefinite separator loop. For a pooled set `W` of saved witnesses,
the exact 16-bit robust-template problem was

```math
\min_{T\in\{0,1\}^{16}}
\max_{(D,w)\in W}|H_{S_{T,D}}(w)|.                   \tag{MC8}
```

This is deliberately a **robust common-template diagnostic**. It is not the
landing-family quantifier `min_(T,D) cap(S_(T,D))`, and therefore an objective
above threshold would still require careful interpretation. Here the
objectives are far below threshold, so the quantifier issue cannot create a
false obstruction.

Pooling the original and first fresh families gives 32,768 witnesses. Exact
CP-SAT constraint generation proves that (MC8) has optimum 72. One optimum
toggles the eight positions

```text
(0,2), (1,3), (2,0), (2,2), (2,3), (3,1), (3,2), (3,3).
```

All pooled energies then lie in `[-72,52]`. A third all-diagonal fresh search
at this exact minimax template again found a threshold witness for every
`D`, with no unresolved cases. The 16,384 new exact witnesses have minimum
absolute energy 210 and record hash

```text
33cec568dafd3d6d504a0f837c18c844d93e90c6f65a91c21347cb630718f8a0
```

Finally, pooling all three rounds gives 49,152 witnesses. The exact robust
common-template optimum is still only 90. One optimal ten-bundle template
puts the entire three-round pool in `[-88,90]`.

The reproducible program and records are:

- `computations/phase2g_micro_template_minimax.py`;
- `computations/results/phase2g_two_round_micro_template_minimax.json`;
- `computations/results/phase2g_three_round_micro_template_minimax.json`;
- `computations/results/phase2g_minimax_template_all_diagonal_audit.json`;
- `computations/results/phase2g_minimax_template_all_diagonal_audit_verify.log`.

This is a compressed negative conclusion from the alternating game. Three
large witness families do not approximate a template-stable cap envelope:
one common bounded template simultaneously suppresses every active witness
to less than half the target threshold. Fresh optimization then exposes a
nearly disjoint face family. Continuing to append individual witness cuts is
therefore not justified. A viable state must encode a face envelope, a
switching-invariant dual measure, or an exact/relaxed cap functional that
transports with the template; raw active witnesses have demonstrated severe
search-path instability.
