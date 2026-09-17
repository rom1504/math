# Wave 36: simultaneous replacement reduces to width-weighted agreement

## Status and verdict

The orientation-doubling, completion-width, conflict, field-dual, and game
upper-bound identities below are **proved**.  The `A_6,A_8,A_9` calculations
are **exact finite calculations**, checked by
`tmp/simultaneous_replacement_r36_check.py`.  The asymptotic width-weighted
agreement statement is an **open target**.

The principal conclusion is two-sided.

1. There is a safe weaker version of (10.1007) in which every favorable
   projective label may use either auxiliary response orientation.  Pairing
   those orientations removes every internal fractional-block moment and
   reduces the desired lower bound to one explicit nonnegative cost:
   completion-width deficit plus anchored conflict.
2. Conversely, the simultaneous game cannot manufacture agreement from
   independent one-block response bounds.  A general upper bound shows that
   (10.1007) already requires more than half of its selected family to lie in
   one radius-`d` center degree, with an additional margin paying the
   one-block response gaps.  At the natural price, the full singleton-
   anchored slices of each of `A_6,A_8,A_9` fail the target for every anchor.
   This is finite/full-family evidence only: the target permits thinning to a
   family of relative density `exp{-O(TL_0)}`.

## 1. Doubling the response orientation is safe

Fix a singleton anchor `v_*`, canonically orient every projective label to
have value `+1` there, and let `F_S` be the favorable **projective** list from
(10.925).  A response state has two logically distinct pieces: its
projective label `[y]` and the auxiliary absolute-value orientation
`sigma in {+-1}`.  Define the enlarged state set

```math
\widetilde\Omega_S^F
=\{(\sigma,y):[y]\in F_S,\ y_{v_*}=1,\ \sigma\in\{\pm1\}\}.
\tag{R36.S1}
```

This enlargement preserves favorable support.  Membership `[y] in F_S`
means that **some** child orientation of this same projective word obeys the
favorable deficit cutoff.  The decoder and the eventual coset hit use only
`[y]`; after decoding one simply reinstates that favorable child orientation.
The `sigma` used in the response game is auxiliary and need not be that
orientation.

It also preserves the nonnegative accounting exactly.  For every `sigma`
and every `y`, favorable or not,

```math
R^S_{A[S]}(\sigma,y)
=\max_{x_{S^c}}\sigma x^{\mathsf T}Ax
\le Q(A)=q_n,
\qquad
\Delta_S(\sigma,y):=q_n-R^S_{A[S]}(\sigma,y)\ge0.
\tag{R36.S2}
```

Thus finite minimax and (10.1005)--(10.1006) apply verbatim to the enlarged
game.  If its row-qualified value satisfies (10.1007), the same accounting
still supplies a low-conflict law of favorable projective labels, hence the
same decoded center.  This is a genuinely weaker sufficient game because its
maximizer has more states.

## 2. Opposite response orientations give an exact width certificate

For one canonically oriented projective label define the outside-completion
polynomial

```math
f_{S,y}(x_{S^c})
=x_{S^c}^{\mathsf T}A[S^c]x_{S^c}
 +2y^{\mathsf T}A[S,S^c]x_{S^c}.
```

Then

```math
Z_{S,+}(y)=\max f_{S,y},
\qquad
Z_{S,-}(y)=-\min f_{S,y},
```

and define the **completion-width deficit**

```math
\boxed{
b_S(y)
=q_n-\frac{Z_{S,+}(y)+Z_{S,-}(y)}2
=q_n-\frac{\max f_{S,y}-\min f_{S,y}}2
\ge0.
}
\tag{R36.S3}
```

The last inequality follows either by adding the two parent bounds in
(R36.S2), whose internal signed energies cancel, or directly from the fact
that every fixed-fiber interval lies inside the full parent interval of
width at most `2q_n`.

Now fix any assignment of favorable projective labels
`mathbf y=(y^S)_(S in G)`.  In the dual law, independently choose the two
response orientations of every fixed label with probability `1/2`.  The
edge words are opposite, so for every selector and every internal edge

```math
m_e^S=\mathbb E[\sigma y_i^Sy_j^S]=0.
```

The label assignment, its anchored conflict, and all of its Hamming medians
are unchanged.  Hence, for the enlarged game (with the row restriction when
the displayed assignment is row-qualified),

```math
\boxed{
\widetilde W^F_{\mathcal G,C}(\kappa)
\ge q_n-
\left\{\mathbb E_{S\in\mathcal G}b_S(y^S)
       +\kappa\mathcal C_{v_*}(\mathbf y)\right\}.
}
\tag{R36.S4}
```

Therefore the completely explicit **width-weighted agreement** condition

```math
\boxed{
\min_{\substack{y^S\in F_S\ (S\in\mathcal G)\\
                 R_{\rm med}(\mathbf y)\le2n(n-1)}}
\left\{\mathbb E_S b_S(y^S)+\kappa\mathcal C_{v_*}(\mathbf y)\right\}
\le\frac{\kappa(d+1)}2
}
\tag{R36.S5}
```

implies the enlarged version of (10.1007), and therefore (10.967).  This
certificate survives the `A_6/A_9` witness-migration wall: it never asks the
unrestricted response optimizer to enter a child-favorable signed state.
It fixes a favorable projective label first, then uses both auxiliary
response orientations of that same label.

## 3. Exact row-free external-field form

For completeness, the conflict penalty has an exact one-field
linearization.  Put

```math
p_i=\Pr_{S\in\mathcal G}\{i\in S\},
\qquad
r_i=\frac{\mathbb E_S[\mathbf1_{\{i\in S\}}y_i^S]}{p_i}.
```

The anchor is omitted below.  Directly from (10.1001),

```math
\mathcal C_{v_*}(\mathbf y)
=\frac12\sum_i p_i(1-r_i^2),
```

so Legendre completion of the square gives

```math
\boxed{
-\kappa\mathcal C_{v_*}(\mathbf y)
=\max_{h\in[-1,1]^{n-1}}
\left\{
\kappa\mathbb E_S\sum_{i\in S\setminus\{v_*\}}h_i y_i^S
-\frac\kappa2\sum_i p_i(1+h_i^2)
\right\}.
}
\tag{R36.S6}
```

Consequently the **row-free** enlarged game has the exact primal form

```math
\widetilde W^F_{\mathcal G}(\kappa)
=\min_{(\xi^S)}\max_h
\left\{-\frac\kappa2\sum_i p_i(1+h_i^2)
+\mathbb E_S\max_{(\sigma,y)\in\widetilde\Omega_S^F}
\left[R^S_{\xi^S}(\sigma,y)
+\kappa\sum_{i\in S\setminus\{v_*\}}h_i y_i\right]\right\}.
\tag{R36.S7}
```

Holding `h` fixed before minimizing gives the valid lower bound

```math
\widetilde W^F_{\mathcal G}(\kappa)
\ge\max_h\left\{-\frac\kappa2\sum_i p_i(1+h_i^2)
+\mathbb E_S V_S^F(h)\right\},
\tag{R36.S8}
```

where `V_S^F(h)` is the corresponding one-block field-tilted minimax.
There is no asserted minimax equality in (R36.S8): the maximum of the local
quadratic branches need not be concave in `h`.  The row restriction also
destroys the selectorwise decomposition and remains a separate input.

## 4. A necessary center-degree bound for the game itself

The reverse direction shows why independent one-block replacement cannot by
itself prove (10.1007).  Let

```math
v_S=\min_{\xi^S}\max_{u\in\widetilde\Omega_S^F}R^S_{\xi^S}(u)
\le q_n,
\qquad
\bar g=q_n-\mathbb E_Sv_S\ge0,
\tag{R36.S9}
```

and define the radius-`d` favorable center degree inside `G` by

```math
\delta_{\mathcal G,d}
=\max_z\frac1{|\mathcal G|}
\left|\{S\in\mathcal G:a_z(S)\le d\}\right|.
\tag{R36.S10}
```

For any favorable label assignment choose a coordinatewise Hamming median
`z`.  Its mean local error is at most the corrected conflict by (10.1002).
On the other hand, every selector outside the radius-`d` center degree has
integer error at least `d+1`.  Therefore, for every assignment,

```math
\boxed{
\mathcal C_{v_*}(\mathbf y)
\ge\mathbb E_Sd_H(z_S,y^S)
\ge(d+1)(1-\delta_{\mathcal G,d}).
}
\tag{R36.S11}
```

Choose independently, for every `S`, a fractional block attaining `v_S`.
Every assignment has response at most `E_Sv_S`; (R36.S11) then gives the
exact game upper bound

```math
\boxed{
\widetilde W^F_{\mathcal G,C}(\kappa)
\le q_n-\bar g
-\kappa(d+1)(1-\delta_{\mathcal G,d}).
}
\tag{R36.S12}
```

The row restriction can only decrease the left side.  In particular,
(10.1007), even for the enlarged game, necessarily requires

```math
\boxed{
\delta_{\mathcal G,d}
\ge\frac12+\frac{\bar g}{\kappa(d+1)}.
}
\tag{R36.S13}
```

Thus the response target already contains a majority center-degree theorem;
it does not create that coherence from averaged block minimality.  This is
not circular as a sufficient lemma, but it sharply limits what an attempted
proof from independent one-block games can accomplish.

## 5. Exact finite full-slice falsifiers and the role of slack

Use the integer favorable allowance `floor(B_(n,m))`, take the entire
singleton-anchored slice, `d=0`, and `kappa=4(n-1)`.  For each matrix, all
full words satisfy the row cap, so the row-free and row-qualified domains
coincide.  Exact enumeration and rational LP certificates give:

| matrix, `m` | `q_n` | allowance | best anchored `delta_(G,0)` | anchored mean `bar g` | best upper bound on `W` | target |
|:--|--:|--:|--:|--:|--:|--:|
| `A_6,3` | 10 | 1 | `1/2` | 2 | -2 | 0 |
| `A_8,5` | 20 | 2 | `5/7` | `58/7` | `26/7` | 6 |
| `A_9,6` | 24 | 2 | `27/56` | anchor-dependent | `-47/7` | 8 |

The displayed `A_9` number is the largest (least obstructive) upper bound
over its nine anchors; every anchor fails.  Hence the full anchored slice
fails (10.1007) on all three exact minimizers even after enlarging the game
to both response orientations.

For `A_6,m=3`, every anchored selector has a unique favorable projective
label and its doubled one-block value is exactly `8`.  On the full ten-set
slice the conflict is exactly `1`, so

```math
\widetilde W^F_{\mathcal G}(20)=8-20=-12
```

exactly.  Exhausting all `2^10-1` nonempty subfamilies shows that every
subfamily of size at least six fails the target, while two size-five
subfamilies pass.  This makes the scope crucial: the finite result rules out
the **full** anchored family and any density above one half in this example;
it does not rule out the existential thinning allowed in (10.1007).

The wall is also genuinely sensitive to favorable slack.  Enlarging the
child-deficit allowance monotonically raises the one-block values and center
degrees.  Explicit conflict-zero coherent assignments certify the doubled
target at allowance `4` for `A_6`, allowance `8` for `A_8`, and allowance
`16` for `A_9`; their worst mean width deficits are respectively

```text
2/5 < 10,       202/35 < 14,       111/14 < 16,
```

against the natural target budgets.  Hence neither the old zero-slack
witness migration nor the new full-slice obstruction is a slack-uniform
finite counterexample.

## 6. Best asymptotic successor

For the replacement implementation, the clean next target is no longer the
opaque value (10.1007).  Seek an anchored family `G` with

```math
\log\beta^{-1}=O(TL_0)
```

and row-qualified favorable projective labels satisfying, for
`kappa=4(n-1)`,

```math
\boxed{
\mathcal C_{v_*}(\mathbf y)
+\frac{\mathbb E_S b_S(y^S)}{4(n-1)}
\le\frac{d+1}{2}.
}
\tag{R36.S14}
```

At the project scales, the new quantitative response requirement is

```math
\mathbb E_S b_S(y^S)=O(n(d+1))
=O\!\left(\frac{Tn^{3/2-2c_0}}{(\log n)^2}+n\right).
\tag{R36.S15}
```

This target explicitly separates the two missing mechanisms:

- **affordable thinning/coherence:** find a subfamily of density
  `exp{-O(TL_0)}` whose favorable labels have project-scale conflict;
- **fiber completion width:** on the same selected labels, show that the
  outside-completion interval misses the parent width by only `O(n(d+1))`
  on average.

The finite evidence says full-slice averaging cannot replace the first item,
while slack can materially improve both.  An asymptotic falsifier for this
replacement implementation would require an unbounded exact-minimizer family
for which every affordable anchored thinning either has center degree below
the necessary threshold (R36.S13), has average completion-width deficit
`omega(n(d+1))`, or has no row-good median.  The finite `A_6,A_8,A_9`
calculations establish none of those unbounded statements.

