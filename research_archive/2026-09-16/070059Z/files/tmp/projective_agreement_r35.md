# Wave 35 B: anchored projective agreement and the exact replacement interface

## Verdict

There is a clean elementary decoder on the complete uniform slice, but its
input is substantially stronger than pairwise existence of compatible labels.
After fixing a common anchor, one must choose **one favorable label per
selector simultaneously** so that a degree-corrected two-selector conflict is
small.  A coordinatewise plurality then gives one global projective center and
loses exactly the anchor probability and the density of the selector
subfamily.  At the Wave 35 scales these losses are affordable.

The decoder does **not** control the row square of the plurality.  Perfect
agreement may decode an arbitrary high-row word.  Thus a row-qualified version
of the multi-block replacement inequality below, or an independent row
theorem, is still necessary for (10.967).

The two cited agreement papers do not supply the required statement in the
present regime.  Their inputs are prescribed local functions rather than the
large favorable lists, and their quantitative conclusions do not reach
distance `d=floor(D/s)`: here the soundness can be
`exp{-O(T L_0)}` and `d=o(T L_0)`, whereas the primary small-soundness bounds
give much coarser local error.  This is a parameter mismatch, not a
counterexample to a possible sharper complete-slice theorem.

The decoder, its conflict identity, and the multi-block minimax identity below
are **Verified**.  The needed replacement lower bound and its row-qualified
form are **Open targets**.  The finite tables are **Exact finite evidence**
only.

## 1. Scales and favorable projective lists

Retain the Wave 34 scales

```math
L_0=n^{3/4-c_0},\qquad
r=\left\lceil\frac{n\log(2k_0)}{L_0}\right\rceil,\qquad
s=\lceil r/T\rceil,\qquad
D=\Theta(k_0),\qquad k_0=\Theta(L_0/\log n),
```

where `1<=T<=n^eta` and `eta<c_0`.  Then

```math
d:=\lfloor D/s\rfloor
=O\!\left(\frac{Tn^{1/2-2c_0}}{(\log n)^2}\right),
\qquad d=o(TL_0).
\tag{R35.B1}
```

For every `m`-set `S`, let `F_S` be the nonempty projective list from
(10.925): its elements are spin classes `[y] in {+-1}^S/{+-1}` satisfying
the favorable child-deficit cutoff inherited from (10.896).  Thus

```math
a_z(S)=\min_{[y]\in F_S}d_{\rm pr}(z_S,y).
```

No useful bound on `|F_S|` is known; it may be exponential.  In particular,
the statement that every pair of lists contains a compatible pair does not
choose labels coherently around cycles and is not an agreement theorem.

## 2. Exact anchored plurality decoder

The anchor deals with the projective orientation exactly.  Fix a nonempty
`A subset [n]`, `|A|=a`, and a projective anchor word `[alpha]`.  Let

```math
\Omega_A=\{S\in\tbinom{[n]}m:A\subseteq S\},
\qquad \mathcal G\subseteq\Omega_A,
\qquad
\beta=\frac{|\mathcal G|}{|\Omega_A|}>0.
```

For every `S in G`, suppose one label `[y^S] in F_S` has been selected with
`[y^S|_A]=[alpha]`.  Choose the unique representative whose restriction to
`A` equals the fixed representative `alpha`.  (For a singleton anchor, every
projective label admits this canonical orientation, so there is no restriction
or list-size loss.)

Let `S` be uniform on `G` and, for `i notin A`, put

```math
p_i=\Pr_{S\sim U(\mathcal G)}\{i\in S\}.
```

Coordinates with `p_i=0` never occur and are omitted from the following sum.
For independent `S,T` uniform on `G`, define the **anchored,
degree-corrected conflict**

```math
\boxed{
\mathcal C_A(\mathbf y)=
\mathbb E_{S,T}
\sum_{i\in(S\cap T)\setminus A}
\frac{\mathbf1\{y_i^S\ne y_i^T\}}{p_i}.
}
\tag{R35.B2}
```

This is ordinary Hamming disagreement between the anchor-oriented
representatives, not projective distance on `S cap T`.  That distinction is
necessary: two labels can agree on the anchor, be opposite on almost every
other overlap coordinate, and nevertheless have small *projective* overlap
distance after flipping one of them.  Such an orientation flip is exactly the
finite shield seen in (10.976).  If one starts instead with a projective
pair-error bound `e`, an anchor of more than `e` exactly matching coordinates
forces the minimizing overlap orientation to be the anchored one; the direct
condition (R35.B2) avoids needing this extra reduction.

For each `i notin A`, write

```math
u_i=\Pr\{i\in S,y_i^S=+1\},\qquad
v_i=\Pr\{i\in S,y_i^S=-1\},\qquad u_i+v_i=p_i.
```

Let `z_A=alpha` and choose `z_i` by plurality between `u_i,v_i` (arbitrarily
if `p_i=0` or in a tie).  The contribution of coordinate `i` to the mean
local error is `min(u_i,v_i)`, while its contribution to (R35.B2) is
`2u_iv_i/p_i`.  Since `max(u_i,v_i)>=p_i/2`,

```math
\boxed{
\mathbb E_{S\sim U(\mathcal G)}d_H(z_S,y^S)
=\sum_{i\notin A}\min(u_i,v_i)
\le \mathcal C_A(\mathbf y).
}
\tag{R35.B3}
```

This proof also covers all zero-degree and tie cases.  Because the local error
is integer, Markov gives the exact center-degree conclusion

```math
\boxed{
U_m\{S:a_z(S)\le d\}
\ge
\beta\frac{\binom{n-a}{m-a}}{\binom nm}
\left(1-\frac{\mathcal C_A(\mathbf y)}{d+1}\right)_+
=\beta\frac{(m)_a}{(n)_a}
\left(1-\frac{\mathcal C_A(\mathbf y)}{d+1}\right)_+ .
}
\tag{R35.B4}
```

Ordinary distance to the selected representative upper-bounds projective
distance, so no further orientation loss occurs in (R35.B4).

If `m/n>=rho>0`, `a<=m/2`, and
`C_A(y)<=(d+1)/2`, then

```math
\frac{(m)_a}{(n)_a}
=\prod_{j=0}^{a-1}\frac{m-j}{n-j}
\ge(\rho/2)^a,
```

and hence

```math
a+\log\beta^{-1}=O(TL_0)
\quad\Longrightarrow\quad
U_m\{S:a_z(S)\le d\}\ge e^{-O(TL_0)}.
\tag{R35.B5}
```

In particular a singleton anchor costs exactly `m/n>=rho`; only the relative
density `beta` then enters the project exponent.

For the full anchored slice `G=Omega_A`, every outside coordinate has
`p_i=(m-a)/(n-a)=:theta`.  The condition becomes the simpler raw pair bound

```math
\mathbb E_{S,T}
d_H\!\left(y^S_{(S\cap T)\setminus A},
           y^T_{(S\cap T)\setminus A}\right)
\le \frac{\theta(d+1)}2.
\tag{R35.B6}
```

The intersection has linear expected size, whereas the permitted number of
conflicts is only `O(d)=o(n)`.  Thus a constant-accuracy agreement theorem is
far too weak.

### Exact sufficient agreement condition

Combining (R35.B4)--(R35.B5), a precise sufficient interface for (10.967) is:

> There exist an anchor `A`, an anchored selector family `G` with
> `a+log beta^{-1}=O(TL_0)`, and a simultaneous transversal
> `[y^S] in F_S` such that `C_A(y)<=(d+1)/2`, whose plurality `z` also obeys
> `R_2(z)<=2n(n-1)`.

Everything except the last row clause is decoded by (R35.B2)--(R35.B5).  The
last clause is genuinely separate; see Section 5.

## 3. Why the cited agreement theorems do not land

### Dinur--Steurer, *Direct Product Testing*, ECCC TR13-179

The primary source is
[ECCC TR13-179](https://eccc.weizmann.ac.il/report/2013/179/).  Its low-
soundness result is Lemma 1.2 for a tuple map
`f:[N]^k->[M]^k`.  With intersection parameter `t<=k/2`, it writes

```math
\delta=(1-\eta)^{t/2},\qquad
\eta=1-\delta^{2/t},
```

and produces *local* direct-product structure with average coordinate error
`O(eta)`.  It is not the high-accuracy global set-ensemble decoder required
here.  At `delta=e^{-O(TL_0)}` and the largest possible
`t=Theta(n)`, its displayed parameter gives

```math
\eta=O(TL_0/n),\qquad k\eta=O(TL_0),
```

so even the most favorable reading of its stated error scale is
`O(TL_0)` local bit errors, while (R35.B1) needs `d=o(TL_0)`.  Taking smaller
`t` only worsens `eta`.

Appendix A discusses the symmetrized set version.  For `k<<sqrt N` it reduces
sets to tuples, and for larger `k<=N/2` it sketches the high-acceptance
reduction.  It does not state the shrinking-soundness, linear-density,
high-accuracy theorem needed here; our ratio window may also have `m>n/2`.

### Dikstein--Dinur, arXiv:2308.09582

The primary source is
[arXiv:2308.09582v2](https://arxiv.org/abs/2308.09582).  Its introductory,
informal main-theorem formulation assumes acceptance
`epsilon>Omega(1/log k)` and describes agreement on a `poly(epsilon)`
fraction with a `poly(1/epsilon)`-sheet cover/list.  Its formal theorem has
additional suitability and base-test-soundness hypotheses and accuracy
`gamma=exp(poly(1/epsilon_0))eta`.  The complete complex is cover-free, so
its topological obstruction is not the problem here.  The parameter regime
is: our acceptance may be
`e^{-O(TL_0)}<<1/log m`, and our permitted error fraction is

```math
d/m=O\!\left(\frac{Tn^{-1/2-2c_0}}{(\log n)^2}\right)=o(1),
```

not a fixed one percent.

Both papers start with one prescribed local function per set.  Here `F_S` is
a potentially exponential list.  If one first proved compressed sublists of
size `L<=e^{O(TL_0)}`, then choosing uniformly would lose at most `L^2` in a
two-query success probability, still affordable at the project exponent.
No such compressed-list theorem or coherent selector is presently available.
Pairwise *existence* of compatible members can use different witnesses on
every pair and supplies neither a transversal nor a cocycle.

Finally, projective labels could be encoded by pair products
`y_iy_j`, but approximate decoding would then need a separate rank-one repair
and changes the error normalization from vertices to edges.  Restricting to
sets containing one anchor coordinate and canonically orienting there is an
exact reduction to the ordinary binary alphabet and is cleaner.

## 4. Exact multi-block replacement game that would suffice

This section isolates, without asserting it follows from global minimality,
the response inequality needed to manufacture the simultaneous transversal
in Section 2.

Specialize the anchor throughout this section to a singleton
`H={i_0}`; every projective label is then compatible with the anchor and has
a unique anchor-oriented representative.  Let `A` denote the exact
minimizing signing, `q=Q(A)=q_n`.  For each selector `S in G`, let
`Omega_S^F` be the anchor-oriented local states whose underlying projective
spin lies in `F_S`.  As in (10.971), write

```math
R^S_\xi(u)=Z_S(u)+2\sum_{e\subset S}\xi_e d_e(u),
\qquad \xi\in[-1,1]^{E(S)},
```

and, for the actual internal block `a[S]`,

```math
\Delta_S(u)=q-R^S_{a[S]}(u)\ge0.
```

An assignment `mathbf u=(u_S)_(S in G)` canonically orients its spin labels
at the anchor and has conflict `C_A(mathbf u)` from (R35.B2).  For `kappa>0`
define the hard-favorable simultaneous game

```math
\boxed{
\begin{aligned}
W_{\mathcal G}^{F}(\kappa)
=\min_{(\xi^S)_{S\in\mathcal G}}
\max_{\mathbf u\in\prod_{S\in\mathcal G}\Omega_S^F}
\bigg\{&\mathbb E_{S\sim U(\mathcal G)}R^S_{\xi^S}(u_S)
-\kappa\mathcal C_A(\mathbf u)\bigg\}.
\end{aligned}
}
\tag{R35.B7}
```

Finite minimax, with a law `mu` on entire assignments, gives

```math
\boxed{
W_{\mathcal G}^{F}(\kappa)
=\max_\mu\left\{
\mathbb E_{\mu,S}Z_S(u_S)
-2\mathbb E_S\sum_{e\subset S}|m_e^S|
-\kappa\mathbb E_\mu\mathcal C_A(\mathbf u)
\right\},
\quad m_e^S=\mathbb E_\mu d_e(u_S).
}
\tag{R35.B8}
```

For every optimizing law, comparison with the actual blocks is the exact
multi-selector analogue of (10.972):

```math
\boxed{
q-W_{\mathcal G}^{F}(\kappa)
=\mathbb E_{\mu,S}\Delta_S(u_S)
+4\mathbb E_S\sum_{e\subset S}(a_em_e^S)_+
+\kappa\mathbb E_\mu\mathcal C_A(\mathbf u).
}
\tag{R35.B9}
```

All terms are nonnegative.  Consequently the precise row-free replacement
inequality

```math
\boxed{
W_{\mathcal G}^{F}(\kappa)
\ge q-\frac{\kappa(d+1)}2
}
\tag{R35.B10}
```

forces `E_mu C_A<=(d+1)/2`; one assignment in the support then has this
conflict, and (R35.B4) decodes a global word with the desired center degree.
It need not have the required row square.

For an exact row-qualified statement, let `Med(mathbf u)` be the set of
coordinatewise plurality words minimizing
`E_S d_H(z_S,y^{S})` (all tie choices are retained), and set

```math
R_{\rm med}(\mathbf u)=\min_{z\in\operatorname{Med}(\mathbf u)}R_2(z).
```

Restrict the maximum in (R35.B7) to assignments satisfying
`R_med(mathbf u)<=C`; denote the resulting value by
`W_(mathcal G,C)^F(kappa)`.  The same proof gives (R35.B8)--(R35.B9), now with
`mu` supported on row-qualified assignments.  Therefore

```math
\boxed{
W_{\mathcal G,\,2n(n-1)}^{F}(\kappa)
\ge q-\frac{\kappa(d+1)}2
}
\tag{R35.B11}
```

together with `a+log beta^{-1}=O(TL_0)` proves (10.967).  This is a precise
multi-block replacement lemma; (R35.B11), not generic global minimality, is
the missing input.

A natural conflict price is `kappa=4(n-1)`, the maximum full quadratic-energy
change per spin flip.  Then (R35.B11) asks for loss at most
`2(n-1)(d+1)`.  By (R35.B1), this is

```math
O\!\left(\frac{Tn^{3/2-2c_0}}{(\log n)^2}+n\right),
```

within the desired power-saving tolerance.  The known independent
sign-rounding loss for one `m=Theta(n)` block is instead

```math
\eta_m=2m\sqrt{(m-1)\log2}=\Theta(n^{3/2}),
```

larger by a power (and logarithms).  More decisively, it applies only to the
unrestricted unpenalized game.  It supplies neither the hard favorable
restriction, the cross-selector penalty, nor the row-qualified support.

The zero-slack `A_9` certificate (10.976) audits the first failure already for
one block: at the exhibited fractional block every favorable state scores at
most `16` while `q_9=24`, and all full-game optimizers are bad.  Thus one may
not infer even the `kappa=0` hard-favorable version of (R35.B10) from exact
global minimality.  Project-scale slack or genuinely simultaneous structure
would have to overcome this witness migration.

## 5. The row square is logically independent

Equation (R35.B3) minimizes a Hamming objective coordinate by coordinate.
It contains no information about

```math
R_2(z)=\|Az\|_2^2.
```

An abstract local system with every list equal to the restrictions of one
fixed high-row word has `C_A=0` and decodes that high-row word.  Hence no
agreement theorem based only on overlap consistency can imply the row clause
in (R35.B11).  Favorability of the restrictions controls child quadratic
energy, not the squared full row fields; no proved inequality converts one to
the other at the target exponent.

There are three honest ways forward:

1. prove the row-qualified replacement inequality (R35.B11) directly;
2. prove (R35.B10) and a separate minimizer-specific theorem that one
   row-good plurality can be chosen among the resulting assignments/ties; or
3. bypass agreement centers and prove the bare one-coset/arbitrary-cut tail
   (10.964)/(10.795).

Merely decoding a global projective word is not enough.

## 6. Exact finite audit

`tmp/projective_agreement_r35_check.py` verifies (R35.B3) with rational
arithmetic, including `p_i=0` and plurality ties, and enumerates the displayed
`A_6,A_8,A_9`.  At `t=d=0`, using the integer favorable allowance
`floor(B_(n,m))`, the best full-center degrees are

| signing, `m` | `floor(B_(n,m))` | best degree | row at selected maximizer |
|:--|--:|--:|--:|
| `A_6,3` | 1 | `1/4` | 30 |
| `A_8,5` | 2 | `15/28` | 64 |
| `A_9,6` | 2 | `29/84` | 128 |

For the anchor `A={0}`, choosing on each anchored selector a favorable label
nearest the displayed maximizing center gives respectively

```text
                 mean local error     corrected conflict
A6,m=3                 1                       1
A8,m=5                24/35                  27/25
A9,m=6                81/56                  86/49
```

in exact agreement with (R35.B3).  All full spins in these three small
examples already obey `R_2<=2n(n-1)` (row ranges `30`, `8..72`, and
`16..128`), so they cannot test the asymptotic row obstruction.  The table is
only a consistency audit; it proves no uniform agreement or replacement
bound.

## 7. Updated target from this attack

The viable agreement route is not “apply a known list-agreement theorem.”
Its exact remaining statement is the anchored, row-qualified simultaneous
replacement inequality (R35.B11), for some selector family whose anchor and
relative-density cost is `O(TL_0)`.  A weaker first milestone is (R35.B10),
which would prove that global minimality creates a coherent favorable
projective transversal, while leaving the row square explicitly open.

Conversely, a scalable family in which every affordable anchored subfamily
and favorable transversal has either
`C_A>(d+1)/2` or only high-row Hamming medians would falsify this particular
agreement implementation.  Pair-specific compatible witnesses, constant-
accuracy decoding, and the existing finite zero-slack walls do not meet that
falsification standard.
