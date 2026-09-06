# Wave 16: exact finite response-to-temporal audit

Status: **Verified finite enumeration and exact algebra.**  All exhaustive
claims below are checked by `tmp/check_finite_bridge_r16.py`; no floating-point
comparison is used.  This memo does not claim an asymptotic theorem.

Run

```bash
/home/math/quadra/.venv/bin/python tmp/check_finite_bridge_r16.py
/home/math/quadra/.venv/bin/python tmp/check_finite_bridge_r16.py --full
```

The optional `--dump-active` flag prints all 904 oriented active spins and
their energy/deficit records.

## 1. The `A_9` capture wall

Use `A_9` from (10.550) and

```math
V=(0,1,2,4,5,6),\qquad V^c=(3,7,8).
```

Write `D=A_9[V]`, and let `C` be `A_9` with the `V\times V` block erased.
Exact Boolean enumeration reconstructs (10.579):

```math
Q(A_9)=24,qquad Q(C)=22,qquad Q(D)=22,qquad q_6=10.
```

Thus the captured excess is `X=12`.  For every order-six minimizer `G`, put
`H_G=C\oplus_VG`.  There are 384 labelled `G` and exactly 40 satisfy
`Q(H_G)=24`.

Here is a compact exact enumeration of the 40 matrices.  Reindex `V` as
`0,...,5`; in lexicographic upper-triangle edge order, bit `k` is one iff the
corresponding edge of `G` is `+1`.  The 15-bit codes are

```text
069a 0a4e 0c8e 129c 1a07 1b9f 1c0b 1e48
1e84 222e 2607 279f 2cdf 2efa 2fae 308e
340d 3628 3682 38bf 3afc 3bce 3eb3 3ed5
3f2b 3f4d 460b 4adf 520d 5ed9 5f8d 62bf
6e6b 6ea7 76b9 778b 7a6d 7ac7 7cad 7ccb
```

The checker asserts this list, and `--full` gives the endpoint counts for
each code.

The common-mosaic responses are

```math
\Phi_D=Q(A_9)-Q(D)=2,qquad
\Phi_G=Q(H_G)-q_6=14,qquad
\Delta_{\rm resp}=12.
```

## 2. Every active signed state

For an oriented active state `\omega=(\sigma,z)` of `H_G`, define

```math
c=\sigma z^{\mathsf T}Cz,quad
g=\sigma z_V^{\mathsf T}Gz_V,quad
d=\sigma z_V^{\mathsf T}Dz_V,
```

and deficits `\delta_G=10-g`, `\delta_D=22-d`.  Always `c+g=24`.
Across the 40 completions there are 904 completion--active-state incidences
(46 distinct oriented states):
496 of orientation `+` and 408 of orientation `-`.  Twenty-four completions
have 21 active states `(12+,9-)`; sixteen have 25 `(13+,12-)`.

The following table is the complete orientation-suppressed energy profile.
The column `A` is the oriented energy `c+d` in the original `A_9`, and
`gain=\delta_D-\delta_G` is the pointwise response-score increase.

| count | `c` | `g` | `delta_G` | `d` | `delta_D` | `A` | gain |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 36 | 14 | 10 | 0 | -14 | 36 | 0 | 36 |
| 156 | 14 | 10 | 0 | -6 | 28 | 8 | 28 |
| 216 | 14 | 10 | 0 | 2 | 20 | 16 | 20 |
| 168 | 14 | 10 | 0 | 10 | 12 | 24 | 12 |
| 64 | 18 | 6 | 4 | -10 | 32 | 8 | 28 |
| 120 | 18 | 6 | 4 | -2 | 24 | 16 | 20 |
| 144 | 18 | 6 | 4 | 6 | 16 | 24 | 12 |

In particular, every active state has pointwise gain at least 12, but the
old response score ranges from `-22` to `2`.  A max-of-affine-functions
comparison therefore does not select one canonical active witness.

There are 312 **persistent** completion--state incidences (14 distinct
oriented states), active for both `A_9` and `H_G`.  Their only profiles are

```math
(c,d,g)=(14,10,10)\quad(168\text{ states}),qquad
(18,6,6)\quad(144\text{ states}).
```

Every completion has between 4 and 13 persistent grounds.  Consequently, for

```math
A_t=C+(1-t)D+tG,qquad 0\le t\le1,
```

one persistent affine energy line is identically 24.  Convexity gives the
opposite inequality `Q(A_t)\le24`, so exactly

```math
\boxed{
Q(A_t)=24,\qquad
\Phi_t:=Q(A_t)-[(1-t)Q(D)+tq_6]=2+12t.
}
```

Thus the response rise is wholly the falling local baseline; it is not work
done by the global norm along the replacement segment.

## 3. Endpoint and induced-deletion audit

Pair every positive ground `p` of `H_G` with every negative ground `n`.
All `H_G` are balanced, so (10.556) makes the total root residual

```math
\mathfrak r_{H_G}(p,n)=|p^{\mathsf T}H_Gn|.
```

There are 5,088 pairs over the 40 completions, with exact histogram

| residual | 0 | 4 | 8 | 16 |
|---:|---:|---:|---:|---:|
| pairs | 1,944 | 1,152 | 1,896 | 96 |

Every completion has two or three residual-16 pairs.  For each pair, retain
each of its two endpoint shores `X`.  Put `h_X=p_X^TH_G[X]p_X`; the total
oriented residual on that shore is `|h_X|`.  With parent normalization
`q_9/9^{3/2}=8/9`, the sign column below is the exact sign of
`Q(H_G[X])-(8/9)|X|^{3/2}`.

| count | `|p^T H_G n|` | `|X|` | `Q(H_G[X])` | `|h_X|` | demand sign |
|---:|---:|---:|---:|---:|:---:|
| 1,824 | 0 | 4 | 8 | 0 | `+` |
| 120 | 0 | 4 | 12 | 0 | `+` |
| 1,032 | 0 | 5 | 8 | 0 | `-` |
| 912 | 0 | 5 | 16 | 0 | `+` |
| 1,152 | 4 | 3 | 6 | 2 | `+` |
| 288 | 4 | 6 | 14 | 2 | `+` |
| 864 | 4 | 6 | 18 | 2 | `+` |
| 1,368 | 8 | 4 | 8 | 4 | `+` |
| 528 | 8 | 4 | 12 | 4 | `+` |
| 1,896 | 8 | 5 | 12 | 4 | `+` |
| 96 | 16 | 4 | 8 | 8 | `+` |
| 96 | 16 | 5 | 8 | 8 | `-` |

Hence every completion has many neutral positive-demand shores: 55 such
directed shores for each 21-state completion and 96 for each 25-state
completion.  Also, neither `V` nor `V^c` is an endpoint shore of any of the
40 hybrids—or of the original `A_9`.  The response block therefore fails
exact retained-set endpoint compatibility; broader fractional compatibility
is not excluded.

The direct non-endpoint deletion retaining `V` has

```math
x_V=22-\frac{16}{3}\sqrt6>0,qquad
\varepsilon_V=Q(D)-q_6=12,qquad
x_V-\varepsilon_V=10-\frac{16}{3}\sqrt6<0.
```

Thus the same 12 units that constitute the response rise occur as terminal
excess and overpay this root demand.  Reclassifying those 12 units as new
endpoint service would double-count them.

## 4. An explicit candidate and exact falsification

For a balanced global minimizer `A`, define its best original root endpoint
resource

```math
R_{\rm end}(A)=
\max_{p\in\mathcal G_+(A),\ n\in\mathcal G_-(A)}
|p^{\mathsf T}An|.
```

The weakest direct scalar response-to-endpoint proposal is

```math
\boxed{\Delta_{\rm resp}\le R_{\rm end}(A).}
\tag{B0}
```

It is non-tautological: its right side is actual residual in the original
temporal matrix, not in an auxiliary replacement.

**Falsified exactly by `A_8`.**  Use (10.445) and the partition

```math
(0,3,4,5)\mid(1,2,6,7).
```

Both original block norms are 12, while `q_4=8` and `Q(C)=16`.  Exhausting
all `48^2` labelled minimizing replacements gives best hybrid norm 20, with
16 attaining pairs.  Therefore

```math
\Phi_D=20-24=-4,qquad
\Phi_G=20-16=4,qquad
\Delta_{\rm resp}=8.
```

But `A_8` has four positive and four negative projective grounds and every
one of their 16 cross-Gram entries is zero.  Hence

```math
\boxed{8=\Delta_{\rm resp}>R_{\rm end}(A_8)=0.}
```

The quantity that repairs the scalar arithmetic is the already-counted
captured internal excess:

```math
\sum_i(Q(D_i)-q_4)=4+4=8.
```

Adding it repairs this example only tautologically, because here
`\Delta_{\rm resp}=X+E_G=X`.  It is not an established temporal terminal
credit: the response blocks are not endpoint shores, and the actual endpoint
shores have zero terminal excess.  Each best hybrid does acquire endpoint residual
8, but that residual belongs to the changed matrix, not the temporal `A_8`.
Transporting it back requires an orientation/replacement theorem.

There is a second, witnesswise failure on `A_9`.  For an active state
`\omega=(\sigma,z)` of `H_G`, let

```math
r_*(\omega)=
\max_{w\text{ an opposite-orientation ground of }H_G}
|z^{\mathsf T}H_Gw|.
```

This is a generous total-root-residual upper bound; it ignores shore
allocation and temporal compatibility.  Even the natural statewise variant
`\Delta_{\rm resp}\le r_*(\omega)` is false:
712 of the 904 active states have `r_*=8<12`; this includes 248 of the 312
persistent grounds.  Every completion has a bad persistent ground, and in
eight completions every persistent ground is bad.  Existential active-state
selection
variant survives (`\max_\omega r_*=16` for every completion), but it changes
the active witness.  A genuine proof would need an active-state mixing/Hall
theorem, not a pointwise charge.

## 5. Comparison with the timing walls

The exact checks reproduce the earlier obstructions:

| matrix | temporal fact | endpoint resource | obstruction |
|---|---|---:|---|
| `A_5` | singleton deletion demand `8-64\sqrt5/25>0` | 0 on the matching pair/tree | positive prefix cannot be charged pointwise |
| `A_7` | `x_1<0`, then `x_2=(108\sqrt{42}-90\sqrt{35})/49-2>0`, while `x_1+x_2<0` | 0 in both orientations on the `6\to5` suffix | ancestor credit cannot serve a descendant-only Hall cut |
| `A_8` | every endpoint shore has demand `8-5\sqrt2>0` and terminal excess 0 | every root endpoint residual is 0 | complete root-neutrality wall |

The `A_9` hybrid audit contains the same phenomenon locally: 2,856 neutral
directed shores have positive raw centered demand, and 1,944 remain positive
after their terminal excess is paid.  Therefore a scalar or outcome-by-outcome endpoint charge
cannot replace the laminar causal cuts of (10.543).

## 6. Finite conclusion

The audit falsifies a direct response-to-original-endpoint bridge and a
statewise hybrid-active repair.  Any surviving theorem must simultaneously:

1. count the internal excess only once, and prove separately when it becomes
   a legitimate temporal terminal credit rather than fresh response service;
2. use resources of the original temporal matrix, or prove a transfer from
   hybrid orientations;
3. mix active states with explicit compatibility, since most individual
   `A_9` witnesses see only 8 units although exceptional pairs see 16; and
4. respect causal timing/antichain cuts, as forced by `A_5`, `A_7`, and `A_8`.

All four requirements are genuine; deleting any one is defeated by an exact
matrix above.
