# Independent audit: radial shell covering and local affine composition

Date: 2026-08-17.

Audited drafts:

* `nearmin_radial_shell_covering.md`;
* `local_affine_interface_composition.md`.

## 0. Executive verdict

| Item | Verdict | Main finding |
|---|---|---|
| RS.1--RS.5 identities and constants | **PASS** | The shell cutoff is inclusive, so the possible equality case is retained and the restricted maximum is exactly equal to the unrestricted one. |
| RS theorem novelty | **ARCHIVE REDISCOVERY / USEFUL REPACKAGING** | Equations (2.2)--(2.5) of `nearmin_deterministic_inequalities.md` already prove the same optimizer, deficit, and half-covering statements.  The common-shell response wording is an immediate corollary, not a new structural theorem. |
| RS information claim | **REPAIR SCOPE** | A proper thin subset of the augmented cut code need not have smaller exponential description complexity.  No strict response quotient follows without a covering/entropy bound. |
| LC.1 endpoint and `3n-4` norm | **PASS AFTER ORIENTATION REPAIR** | The row calculation is exact.  The opening claim that one can orient a quadratic ground state by replacing `x` with `-x` is false; introduce `rho H_A(x)=Q(A)` and carry `rho` through the lower witnesses. |
| LC.1 microcanonical compiler | **PASS** | The bilinear union bound proves `C sqrt(ns(n+s))`.  It genuinely improves Theorem 21.66 for `sqrt(n) << s <= n`, but has the same `n sqrt(s)` order at the canonical `s=Theta(sqrt(n))` scale. |
| LC.2 anchored insertion | **PASS AFTER ORIENTATION REPAIR** | The two-sided interval about `Q(A)+3n-4` is correct.  Its cap conclusion is exactly the order of the archived random near-order transfer (10.1684), though the selected-field anchor is additional response information. |
| LC.3 scalar residual lower bound | **PASS WITH ASYMPTOTIC QUALIFIER** | The displayed inequality is correct.  For the `G=5` target it is `Omega(n sqrt(s))` only as `s -> infinity` (the displayed lower bound is negative for the smallest allowed `s=5`). |
| LC.4 all-endpoint ceiling | **PASS AFTER ORIENTATION REPAIR** | Odd Boolean majority functions with Boolean weights are injective, and the `ns` lower bound follows.  If `[n]\I` is empty, choose an arbitrary row word separately; the canonical regime has a nonempty bulk. |
| LC.15--LC.16 accumulation | **PASS AS A CEILING ON THIS CERTIFICATE** | The parameter sum cannot be made summable by refining relative increments.  This proves failure of triangle accumulation of the stated upper bounds, not a lower bound on the actual errors of every correlated multilevel construction. |

The required repairs are local.  After them the composition draft is a valid
one-step response theorem plus two useful no-go statements.  It is not an
improvement of near-order transfer and does not yet use the low-field
near-top property in its positive insertion theorem.

## 1. Radial shell theorem

For every augmented cut word `z`, flipping `F` gives

```math
\langle a^F,z\rangle
=\langle a,z\rangle-2\sum_{e\in F}a_ez_e
=Q(a)-d_a(z)-2s_F(z).
```

Maximization proves (RS.1) exactly.  If `z_F` maximizes the flipped signing,
then

```math
Q(a)-d_a(z_F)-2s_F(z_F)
=Q(a^F)\ge M_n\ge Q(a)-\eta,
```

so

```math
d_a(z_F)+2s_F(z_F)\le\eta.                    \tag{A.1}
```

Because `s_F(z_F)>=-|F|`, (RS.2) and (RS.3) follow with exactly the stated
constants.  Also, if `|F|<=r` and `d_a(z)>eta+2r`, then

```math
-d_a(z)-2s_F(z)
\le-d_a(z)+2|F|<-\eta.
```

The full maximum is at least `-eta`, so no word outside the common shell can
maximize.  Crucially, a word with

```math
d_a(z)=\eta+2r,
\qquad s_F(z)=-r
```

has value exactly `-eta`; it is included because (RS.4) uses `<=`, not `<`.
Thus the common-shell restriction remains exact at equality.

Finally,

```math
s_F(z)=|F\setminus D_z|-|F\cap D_z|
      =|F|-2|F\cap D_z|,
```

which rearranges (A.1) to (RS.5).  The displayed RS.5 quantifier is for
sets of exactly `r` edges.  The natural ball version is the same formula
with `r` replaced on its right-hand side by `|F|`.

### Archive comparison and information content

The draft `nearmin_deterministic_inequalities.md` already records exactly

```math
Q(a^F)=\max_z\{Q(a)-d_a(z)-2s_F(z)\}
```

and (A.1), the deficit bound `d_a(z)<=eta+2|F|`, and the disagreement form
(its equations (2.2)--(2.5)).  It even freezes the `|F|<=r` statement as
Candidate 1 and labels it an archive rediscovery.  RS.1's useful new
presentation is that the *same shell* answers all queries in the declared
Hamming ball, but this is an immediate quantifier extraction from those
inequalities.

The sentence saying that the retained information is visibly smaller than
the full landscape needs qualification.  Even when
`eta+2r=o(n^(3/2))`, the shell can still contain exponentially many words
and require the same `Theta(n)` response bits as the full projective code.
At most, once the known lower bound `Q(a)=Omega(n^(3/2))` is invoked, a thin
shell cannot contain both members of an antipodal pair; that saves no
macroscopic information rate.  A strict reduction requires a proved small
presentation, metric-entropy bound, or congruence.  The final two paragraphs
of the draft correctly acknowledge precisely this missing step.

## 2. The orientation repair needed throughout the composition draft

For a quadratic Hamiltonian, `H_A(-x)=H_A(x)`.  Thus one cannot orient an
absolute ground state onto the positive side by replacing `x` by `-x`.
Choose instead

```math
\rho\in\{\mathord\pm1\},
\qquad \rho H_A(x)=Q(A).
```

The low-field set used to construct the canonical frame must be defined
with `rho A`, as in the independent audit of
`local_field_affine_shell_algebra.md`.  The frame vectors themselves remain
`x,x^{\{i_1\}},...,x^{\{i_k\}}`, so all endpoint and compiler calculations
are unchanged.

For (LC.2), the exact trust identity becomes

```math
\mathcal B_A(ag)=Q(A)+|a|\|g\|_1
```

by choosing quadratic orientation `sigma=rho` at `x` when `a>=0` and at
`-x` when `a<0`.  For LC.4, evaluate the physical parent at endpoint
`rho r` rather than necessarily at `r`; then

```math
H_A(x)+x^TB(\rho r)+H_C(\rho r)
=\rho\{Q(A)+x^TBr\}+H_C(r),
```

whose absolute value is at least `Q(A)+x^TBr-Q(C)`.  These changes preserve
all displayed constants and preserve the original top-left block `A`.

## 3. Endpoint calculation (LC.1--LC.2)

Let `t=sum_j epsilon_j=3`.  Outside `I`,

```math
x_i g_i=3.
```

At port `i_j`,

```math
x_{i_j}g_{i_j}=3-2\epsilon_j,
```

which is `1` at the `k/2+1` positive single-flip coordinates and `5` at
the `k/2-1` negative ones.  Hence

```math
\|g\|_1
=3(n-k)+(k/2+1)+5(k/2-1)
=3n-4.
```

All switched coordinates are positive, so `sgn(g)=x`.  This verifies every
claim in (LC.1), including parity and the bound `|g_i|<=s` for `s=k+1>=5`.
The trust identity then follows as above.  It is the *selected response*
that has the exact increment `3n-4`; the completed parent's cap lies in the
interval (LC.7), not at that value exactly.

## 4. Bilinear microcanonical compiler (LC.3)

After a column switch take `eta_*=1`.  A row with sum `g_i` is obtained by
choosing a subset `P_i` of size `(s+g_i)/2` uniformly and setting its entries
to `+1`.  For fixed `eta`,

```math
X_i=b_i^Teta-{\langle eta,1\rangle\over s}g_i
=2\left(\sum_{j\in P_i}\eta_j
 -{|P_i|\over s}\sum_j\eta_j\right).
```

This is a centred sample-without-replacement sum.  Hoeffding comparison and
Hoeffding's lemma give a two-sided subgaussian MGF with variance proxy
`C_0s`.  Independence of rows gives proxy `C_0ns` for
`sum_i z_iX_i`.  A union bound over `2^(n+s)` pairs `(z,eta)` yields

```math
\max_{z,eta}|\sum_i z_iX_i|
\le C\sqrt{ns(n+s)}.
```

For fixed `eta`, maximizing over `z` is exactly `sum_i|X_i|`, proving
(LC.3).  There is no missing factor of `s` or logarithm.

Compared with Theorem 21.66, whose error is

```math
O(n\sqrt s+s^{3/2}\sqrt n),
```

the new bound is genuinely sharper for `sqrt(n)<<s<=n`; it removes the
rowwise-absolute concentration penalty.  At `s=Theta(sqrt n)`, both bounds
are `Theta(n^(5/4))`.  The label “sharp-order” is justified for bounded
`g_i` and `s<=n` by LC.3, but not by the supplied lower bound when `s>>n`.

## 5. Anchored insertion (LC.7)

For each shore endpoint, write

```math
B\eta=a g+e_\eta,
\qquad a={\langle\eta,\eta_*\rangle\over s},
\qquad \|e_\eta\|_1\le R.
```

The trust response is `l_1`-Lipschitz, so

```math
\mathcal B_A(B\eta)
\le Q(A)+|a|\|g\|_1+R
\le Q(A)+\|g\|_1+R.
```

Introducing the global absolute sign in the block Hamiltonian gives the
upper bound `Q(P)<=Q(A)+||g||_1+R+Q(C)`.  At `eta_*`, the declared old spin
and the appropriate quadratic orientation give at least
`Q(A)+||g||_1-Q(C)`.  A standard random hollow shore has
`Q(C)<=C_1s^(3/2)`, so (LC.7) follows exactly.

For `s=o(n)`, both `sqrt(ns(n+s))` and `s^(3/2)` are `o(n^(3/2))`, proving
(LC.8).  For `s=Theta(sqrt n)`, the error is `O(n^(5/4))` and (LC.11a) is
correct.

The cap conclusion is not a stronger near-order transfer theorem.  Ledger
(10.1684), with `h=s`, already gives for every child

```math
Q(P)\le Q(A)+O(\sqrt{ns(n+s)})
```

after random completion.  LC.2 adds a deterministic exact selected field
and a matching lower anchor, but its cap error has the same order.  Also,
none of LC.1--LC.2 uses the low-field near-top estimate `Delta_k`: for the
positive theorem, any absolute ground state and any `k` coordinates supply
the same anchor.  The low-field algebra becomes relevant only if one tries
to preserve many selectors, which LC.4 obstructs in the direct
physicalization.

## 6. The two ceilings

For LC.3, averaging and reverse triangle give row by row

```math
\mathbb E_\eta|b_i^T\eta-a g_i|
\ge \mathbb E|b_i^T\eta|-|g_i|\mathbb E|a|
\ge \sqrt{s/2}-{G\over\sqrt s}.
```

Here the first constant is the sharp `p=1` Khinchine lower constant and
`E|<eta,eta_*>|<=sqrt(s)`.  Summation and `max>=average` prove (LC.12).
For `G=5`, this displayed bound becomes positive only for sufficiently
large odd `s` (in fact from `s=9` onward), so the `Omega(n sqrt s)` sentence
must be read asymptotically, not uniformly from the draft's smallest
`s=5`.

For LC.4, switch two bulk row words so one is `u=1`.  If the other Boolean
word `v` differs on a nonempty set `D`, take all endpoint signs on `D`
equal and choose the complementary sign sum to have absolute value at most
one.  Then `u^Teta` and `v^Teta` have opposite signs.  Hence the odd
threshold map `u -> (sign(u^Teta))_eta` is injective, and all bulk switched
rows equal one word `r`.  At endpoint `r`,

```math
x^TBr\ge(n-k)s-ks=(n-2k)s.
```

After the `rho` repair above, the parent lower bound (LC.14) follows.  If
there are no bulk rows, one may choose any `r`; the bound is then only the
trivial negative lower bound.  In the intended canonical regime
`k=s-1<=n-1`, a bulk row exists.

## 7. Repeated accumulation

With `theta_j=s_j/n_j` and `s_j<=n_j`, the remainder in (LC.7), divided by
`n_j^(3/2)`, is

```math
O(\sqrt{\theta_j(1+\theta_j)}+\theta_j^{3/2})
=O(\sqrt{\theta_j}+\theta_j^{3/2}).
```

Across a fixed logarithmic order change with vanishing maximum step,
`sum_j theta_j=Theta(1)`.  If `theta_j<=delta`, then

```math
\sqrt{\theta_j}\ge {\theta_j\over\sqrt\delta},
```

which proves (LC.16).  Therefore the *available triangle-inequality error
certificate* worsens under finer subdivision and cannot supply a summable
recurrence.  This does not prove that the realized errors themselves have
that sign or magnitude, and it does not exclude a correlated construction
that cancels errors jointly across levels.  The draft's final qualification
states this distinction correctly; the earlier prose should consistently
say “certified/triangle-accumulated loss.”

## 8. Recommended canonical disposition

1. Keep RS.1 as a concise radial-response corollary or declared-query
   formulation of the already proved flip theorem; do not number it as an
   independent new frontier theorem.
2. Repair `rho` in every LC statement and proof before promotion.
3. Promote the bilinear compiler improvement and the selected-field
   anchored insertion as response-level results, while explicitly noting
   that (10.1684) already supplies the same cap-transfer order.
4. Promote LC.3--LC.4 as scoped no-go results.  Neither rules out the joint
   residual/child cancellation mechanism named in the draft.
5. Treat LC.15--LC.16 only as failure of repeated *separately paid* error
   bounds, not as an impossibility theorem for multilevel composition.
