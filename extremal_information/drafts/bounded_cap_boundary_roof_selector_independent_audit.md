# Independent audit: bounded-cap boundary roof selector

**Verdict: PASS.**

The selector theorem itself is correct.  In particular, the biased
exact-sign fill has the claimed operator norm and cap, the crossing-block
estimate has the right factor two, the near/far split closes with the stated
constants, and (BR.9) suppresses the negative outer absolute-value channel.
The proof of Theorem 21.8 also supplies the two *directed* matched-roof
deficits asserted in (BR.24), not merely their projective average.

The final source explicitly orients the source signing so that its positive
maximum `P` equals its absolute cap `Q(A)`.  This supplies the common-roof
hypothesis for all switched children.  It also parameterizes the cap bound
by an explicit coefficient `C_P`, so the parent-cap constant has the stated
uniform dependence.

The result is query-owned and noncircular under the repository's declared
contextual model.  It does not compute a child optimizer or assume the
order-`2n` optimum.  Its public exact query fills can, however, have
quadratic description length; query information is uncharged here.  This
scope should remain explicit.

## 1. Frozen source

```text
extremal_information/drafts/bounded_cap_boundary_roof_selector.md
sha256 631d5ddcc79fc868c0086a8f9bb469d201980df71455461e7d0b5f3675251e87
```

The audit uses the repository's edge-sum normalization

```math
H_A(x)=\frac12x^TAx=\sum_{i<j}A_{ij}x_ix_j.
```

## 2. Profiles, roofs, and projective Lipschitz constants: PASS

For the block parent, the cross energy is exactly `x^TBy`; there is no
missing factor two.  Since `Q(C)<=P`, for either `sigma=+-1`,

```math
f_C^\sigma(y)
\le P+\max_x x^TBy
=P+\|By\|_1=U(y).
```

Both `f_C^sigma` and `U` are even in `y`: in the first expression one also
changes `x` to `-x`.  Thus one may choose the representative of `y` at
ordinary Hamming distance `d=d_P(y,q)<=n/2` from `q`.  Then

```math
\begin{aligned}
|f_C^\sigma(y)-f_C^\sigma(q)|
&\le\max_x|x^TB(y-q)|,\\
|U(y)-U(q)|
&\le\|B(y-q)\|_1,\\
\max_x|x^TB(y-q)|
&\le\sqrt n\,\|B\|_{2\to2}\,2\sqrt d
\le2Ln\sqrt d.
\end{aligned}
```

This proves both inequalities in (BR.3) with the displayed constant.

## 3. The biased exact-sign fill: PASS

For `p=lambda/sqrt(n)`, the prescribed upper-triangular distribution has

```math
\mathbb E(D_q)_{ij}=p q_iq_j.
```

Hence `E_q=D_q-pR_q` is hollow, symmetric, centred, and has independent
bounded upper-triangular entries.  For a fixed unit vector, its quadratic
form is a sum with squared coefficient mass controlled by
`sum_(i<j)v_i^2v_j^2<=1/2`.  Hoeffding plus a constant sphere net therefore
gives one realization with

```math
\|E_q\|_{2\to2}\le K_0\sqrt n
```

for an absolute `K_0`.  This remains uniform in every fixed `lambda`, since
`p<=1` eventually and the centred entries stay uniformly bounded.

The hollow rank-one matrix is `R_q=qq^T-I`, so its eigenvalues are `n-1`
and `-1`; therefore `||R_q||=n-1`.  Consequently

```math
Q(D_q)
\le\frac n2\|D_q\|
\le\frac n2\{p(n-1)+K_0\sqrt n\}
=\frac{\lambda+K_0+o(1)}2n^{3/2}.
```

Thus (BR.6)--(BR.7), exact hollowness, and the complete-sign constraint all
check.

## 4. Crossing block and two-sided rank-one bounds: PASS

Gauge `q` to the all-positive word and let `S` be the `d` coordinates where
`y` differs.  With `h=d(n-d)`, the rank-one energy drop is exactly

```math
H_{pR_q}(q)-H_{pR_q}(y)=2ph.
```

Only edges crossing `(S,S^c)` survive in the error difference, and

```math
|H_{E_q}(q)-H_{E_q}(y)|
=2|\mathbf1_S^TE_q\mathbf1_{S^c}|
\le2K_0\sqrt{nh}.
```

This proves (BR.8), including its factor two.

Also

```math
H_{R_q}(q)+H_{R_q}(y)=n(n-1)-2h.
```

Since `h<=n^2/4` and each error energy has absolute value at most
`K_0n^(3/2)/2`,

```math
H_{D_q}(q)+H_{D_q}(y)
\ge(\lambda/2-K_0-o(1))n^{3/2}.
```

This is (BR.9).  The important asymmetry is real: the rank-one mean is
`+(lambda/2+o(1))n^(3/2)` at `q`, whereas its minimum over Boolean words is
only `-O(lambda sqrt(n))`.

## 5. BR.2 near/far and absolute-channel audit: PASS

Let `T=U(q)+H_{D_q}(q)`.  The matched equality produces an actual target
parent configuration of value `T`, proving (BR.15).

### Near positive channel

For `d<=theta n`, the directed deficit at `q`, profile Lipschitzness, and
the crossing estimate give

```math
f_C^+(y)+H_{D_q}(y)
\le T-\delta n^{3/2}
 +2Ln\sqrt d+2K_0\sqrt{nd(n-d)}.
```

Both errors together are at most
`2(L+K_0)sqrt(theta)n^(3/2)`, hence at most
`delta n^(3/2)/4` by (BR.13).  This proves (BR.16).

### Far positive channel

On `theta n<=d<=n/2`,

```math
h=d(n-d)\ge\theta(1-\theta)n^2.
```

The coarse bounds

```math
U(y)-U(q)\le2Ln^{3/2},
\qquad
2K_0\sqrt{nh}\le2K_0n^{3/2}
```

are deliberately loose but valid.  Subtracting the rank-one drop
`2lambda theta(1-theta)n^(3/2)` and using the first condition in (BR.14)
proves (BR.17).

### Negative absolute channel

The exact negative parent channel is

```math
f_C^-(y)-H_{D_q}(y).
```

Using its roof, the coarse change in `U`, and (BR.9) gives

```math
f_C^-(y)-H_{D_q}(y)
\le T-(\lambda/2-K_0-2L-o(1))n^{3/2}.
```

The second condition in (BR.14) makes this at most
`T-delta n^(3/2)` for all sufficiently large `n`.  Thus the proof really
controls the outer negative absolute-cap channel rather than silently
dropping it.

The projective near/far cases cover every boundary word because all three
relevant functions are even.  Taking the maximum of both signs yields
(BR.19), and comparison with (BR.15) actually gives a
`3delta/4` gap; the theorem's `delta/2` is safe.

Finally,

```math
\|B\|_{\infty\to1}
\le\sqrt n\,\|B\|_{2\to2}\sqrt n
\le Ln^{3/2}.
```

Therefore the triangle bound (BR.20) proves the complete-parent cap.  All
three blocks have exact signs and the diagonal blocks are hollow, so the
parent is indeed a complete hollow signing of order `2n`.

## 6. What Theorem 21.8 supplies: directed deficits PASS

For `s_y=u_* odot sign(By)`, Lemma 21.7 gives exactly

```math
f_{A^{s_y}}^+(y)=P+\|By\|_1.
```

For a distinct codeword `z`, changing variables under switching gives the
stronger identity

```math
f_{A^{s_z}}^+(y)
=P+\|By\|_1-\Delta_A(s_z\mathbin\odot By).
```

The probabilistic event constructed in Theorem 21.8 is union-bounded over
**ordered** code pairs.  For every near-top spin it supplies at least
`delta_0 n` mismatch rows of magnitude at least `a sqrt(n)` in the first
query field; outside the near-top set the energy deficit is at least
`d_0n^(3/2)`.  With that theorem's `d=2a delta_0<d_0`, this proves

```math
\Delta_A(s_z\mathbin\odot By)\ge d n^{3/2}.
```

Repeating the ordered event with `(z,y)` gives the second field deficit.
Thus both inequalities in (BR.24) are genuinely present in the proof of
Theorem 21.8.  Equation (21.37) records only half their oscillation, but the
individual directed statements are stronger and valid.

### Orientation hypothesis

BR.2 uses

```math
Q(C_y)=Q(A)\le P,
```

whereas Theorem 21.8 in isolation defines only `P=max H_A`.  A general
hollow signing can have `Q(A)>P`.  The frozen source now correctly assumes
the global orientation

```text
P=max H_A=Q(A)
```

and qualifies its conclusion by this ambient assumption.  Switching
preserves `Q`, so every child meets BR.10 and the application is rigorous.

## 7. Query ownership and circularity: PASS with a scope warning

For each public pole `q`, BR.1 chooses a fill using only `q`, the public
constants, and an operator-norm event.  The same frozen `D_q` is attached to
every child.  The exponentially many fills need not satisfy one joint
random event: each nonempty finite existence problem can be solved and then
the resulting public query bank frozen.

The construction does not evaluate the child landscape, encode its full
energy table, or use `M_(2n)`.  The source maximizer `u_*` is used upstream
to define the switching orbit in Theorem 21.8, but no target-order optimum
is assumed.  The scalar response packing is therefore noncircular.

There is nevertheless no query-complexity bound.  Although its structured
mean is described by `(q,lambda)`, one exact error realization `E_q` may
require `Theta(n^2)` public bits, and the query bank may contain
`exp(Theta(n))` such matrices.  The result lower-bounds reusable **child
state** when arbitrary declared public contexts are uncharged.  It should
not be presented as an efficient or low-information query compiler without
a separate derandomization/description theorem.  Item 6 of the frozen
source states precisely this limitation.

## 8. Archive comparison

| Archived result | Actual relation |
|---|---|
| UP.1 | No collision. UP.1 assumes one old coordinate is pinned against every complete child. BR.2 only localizes a named boundary pole for one spectrally Lipschitz matched-roof family and allows old optimizers to change. |
| Theorem 21.29 | Its exact coordinate lock pays a quadratic calibration. BR.2 obtains only projective localization and pays `Theta(n^(3/2))`, so it is a different, weaker selector premise with the desired cap scale. |
| AO.2 | AO.2 exposes fixed-scale augmented-cut geometry through a repeated rank-one cross shore of width `Theta(sqrt n)`. It does not prove the present matched-roof/two-absolute-channel scalarization. |
| BCX.1--BCX.3 | BCX already supplies an all-spins-free bounded-cap scalar packing for the special regular-Hadamard switching code. Hence BR.2 is not new on that special family. It extends the physical compiler to arbitrary directed matched roofs output by Theorem 21.8, including the conditional exact-minimizer application. |
| OV.1 | BR.2 respects the orientation ceiling: its boundary fill itself spends `Theta(n^(3/2))` internal cap, enough for target-scale absolute-channel visibility. |
| OV.2 | BR.1 is the square-symmetric analogue of OV.2's biased exact-sign rounding with a centred spectral error. The frozen source identifies this archive reuse; the new content is its crossing estimate and combination with matched-roof stability and (BR.9). |
| WS.1--WS.3 | Compatible. Those theorems charge reusable old-witness dictionaries; BR uses a query-indexed public boundary selector and proves a child-state packing. |

The earlier independent audit of the switching-broadcast route explicitly
left a low-cap boundary selector/restricted anti-pin open.  With the frozen
orientation hypothesis, BR.2 fills exactly that conditional gap; it does
not prove `L_tail`, convergence, or a cross-order recurrence.

## 9. Disposition

The frozen source implements all first-pass repairs:

1. Section 4 globally orients the source signing so `P=Q(A)` and states the
   conclusion only under that ambient assumption.
2. BR.2 assumes `P<=C_Pn^(3/2)` and returns an explicit cap coefficient
   depending on `(L,delta,C_P)`.
3. The archive comparison records the OV.2 reuse, the pre-existing BCX
   special case, and the uncharged query-description complexity.

The result is therefore a rigorous conditional implication:

```text
globally oriented bounded-cap sign quadratic
+ upper-tail hypothesis (21.38)
+ Theorem 21.8 directed matched roofs
=> an exp(Omega(n)) all-spins-free scalar physical contextual packing
   with exact order-2n parents of cap O(n^(3/2)).
```

It proves no upper-tail estimate, convergence, recurrence, or target-order
optimality statement.  No further repair is required.
