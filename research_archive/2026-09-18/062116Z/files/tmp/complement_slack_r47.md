# Wave 47A: complement slack cannot follow from minimality near the diagonal

## Status

- **Verified:** if `T=S^c`, the complement signing is a switching of the
  signing obtained by flipping only the internal edges of `T`.  Consequently
  every active complement slack is at most `2|T|(|T|-1)`.
- **Verified scalable obstruction:** for a fixed noncentral switching
  multiplier `0<theta<1`, the high-slack event in (10.1191) is empty whenever
  the complement size `r=n-m=o(n^(3/4))`.  At `m=n-1` every active slack is
  exactly zero for every exact minimizer and every scalar optimizer, although
  every nonempty column has mass at least `1/n`.
- **Verified:** the exact fractional replacement LP and its dual expose the
  quantifier gap.  Low-slack incidences control selected responses, whereas a
  contradiction requires a single integral replacement to beat every
  response after it is revealed.  Convexifying the replacement changes this
  game and can create a spurious positive margin as large as `q_n/4`.
- **Falsified:** a uniform theorem of the form "small conditional high-slack
  fraction implies a strict improving signing replacement" cannot follow
  from exact minimality and scalar optimality alone.  The universal
  `m=n-1` orbit already has zero slack and no improving signing.  The stored
  `A_8,A_9` examples additionally exhibit explicit adaptive witness migration
  after the natural low-slack replacements are combined.
- **Open:** none of this falsifies the fixed-density arithmetic route.  When
  `m/n` stays below one by a fixed amount, the new slack cap is vacuous at the
  `q_n=Theta(n^(3/2))` scale.  The fixed-density mass and high-slack claims,
  and hence the conditional row conclusion, remain open.

Throughout this note `r=n-m` is the complement size.  It is unrelated to the
power-saving/exponent parameter used elsewhere in the ledger.

## 1. Exact small-complement factorization and slack cap

Let `A` be an exact order-`n` minimizer, `Q(A)=q=q_n`, and put

```math
B_S=-A+2P_SAP_S,
\qquad T=S^c,
\qquad r=|T|.
```

Write `A^{E(T)}` for the signing obtained by flipping precisely the internal
edges of `T`, and let `D_T` be diagonal with entries `-1` on `T` and `+1` on
`S`.  Inspection of the three edge classes (inside `S`, crossing `S,T`, and
inside `T`) gives the exact identity

```math
\boxed{B_S=D_TA^{E(T)}D_T.}
\tag{R47.A1}
```

Thus switching invariance gives

```math
Q(B_S)=Q(A^{E(T)}).
```

For every Boolean word `x`, the triangle inequality gives

```math
\begin{aligned}
|x^{\mathsf T}A^{E(T)}x|
&=|x^{\mathsf T}Ax-2x_T^{\mathsf T}A[T]x_T|\\
&\le q+2Q(A[T])
\le q+2r(r-1).
\end{aligned}
```

Therefore, if `(d,S)` is active and

```math
L_S(d)=\langle B_S,d\rangle=q+\Gamma_S,
```

then

```math
\boxed{
0\le\Gamma_S
\le Q(B_S)-q
\le2Q(A[T])
\le2r(r-1).}
\tag{R47.A2}
```

This is a deterministic pointwise statement; conditioning on `I_d` cannot
improve it.  The verified lower bound in ledger section 1.7 implies that, for
some universal `c_0>0` and all sufficiently large `n`,
`q_n>=c_0n^(3/2)`.  Hence, for any fixed `theta in (0,1)`,

```math
r=o(n^{3/4})
\quad\Longrightarrow\quad
2r(r-1)<\frac{1-\theta}{\theta}q_n
```

eventually.  In this regime

```math
\boxed{
\Pr_{S\mid I_d}
\{\theta\Gamma_S>(1-\theta)q_n\}=0
\quad\text{for every nonempty column }d.}
\tag{R47.A3}
```

In particular, this applies when the block-flip ratio `k/n` tends to a fixed
nonzero noncentral value, so that the multiplier in (10.1185) tends to a
constant `theta=(1-2k/n)^2<1`.  It does not apply if that switching ratio is
itself sent to zero and `theta` tends to one; in that limit the row-contraction
factor `1-theta` simultaneously degenerates.

The endpoint `m=n-1` is even sharper.  Then `E(T)` is empty, so `B_S` is
switching-equivalent to `A` and `Q(B_S)=q`.  Every active incidence therefore
has

```math
\boxed{\Gamma_S=0.}
\tag{R47.A4}
```

This does not reflect a mass failure.  There are `n` selectors, so every
nonempty incidence column has `alpha_d>=1/n` and
`-log alpha_d<=log n`.  Equivalently, flipping one vertex of any parent
ground gives a column incident to the selector omitting that vertex.  Thus a
scalar optimizer over nonempty columns automatically has polynomial mass but
zero high-slack fraction.

The row remains separate.  The generic exact-minimizer estimate gives only

```math
R_2(d)=x^{\mathsf T}A^2x
\le n\lVert A\rVert_{\rm op}^2
\le2nq_n=O(n^{5/2}),
```

which misses the required `O(n^(9/4-c))` scale.  Equations (R47.A3)--(R47.A4)
say that the arithmetic retention proof cannot supply the row estimate here;
they do not assert that every such scalar optimizer actually has large row.

### Fixed-density qualification

The established convergence bridge needs the restriction recurrence on one
fixed active high-ratio window, followed by exact landing.  It may take a
compact window

```math
[p_0,p_1]\subset(1/\sqrt2,1),
\qquad p=m/n,
```

with `p_1<1`.  Then `r=Theta(n)` and (R47.A2) is only `O(n^2)`, hence gives
no information at the `q_n=Theta(n^(3/2))` scale.  The scalable obstruction
therefore rules out a uniform near-identity implementation; it does **not**
retire the fixed-density complement arithmetic route.

## 2. Exact replacement game and the adaptive-witness LP wall

For any signing `C`, let `Omega` be the oriented projective states and write

```math
s_e^C(d)=c_e\,\sigma x_ix_j,
\qquad E_C(d)=2\sum_es_e^C(d).
```

Let `H` be an allowed edge block.  Flipping `F subset H` gives

```math
E_{C^F}(d)=E_C(d)-4\sum_{e\in F}s_e^C(d).
```

The exact integral all-response margin at the global minimum `q=q_n` is

```math
\begin{aligned}
\eta_{\rm int}(C,H)
&=\max_{F\subseteq H}\min_{d\in\Omega}
\left\{\sum_{e\in F}s_e^C(d)+\frac{q-E_C(d)}4\right\}\\
&=\frac14\left[q-\min_{F\subseteq H}Q(C^F)\right]
\le0.
\end{aligned}
\tag{R47.A5}
```

Thus a strict contradiction to signing minimality is exactly the statement
`eta_int>0`.  Lowering the currently selected incidence responses is not
enough: the minimizing state in (R47.A5) is chosen after the replacement.

Now fractionally flip with parameters `p_e in [0,1]`, replacing `c_e` by
`(1-2p_e)c_e` on `H`.  Finite minimax gives the exact LP dual

```math
\boxed{
\begin{aligned}
\eta_{\rm frac}(C,H)
&=\max_{p\in[0,1]^H}\min_d
\left\{\sum_{e\in H}p_es_e^C(d)+\frac{q-E_C(d)}4\right\}\\
&=\frac14\left[q-\min_{p\in[0,1]^H}Q(C(p))\right]\\
&=\min_{\mu\in\Delta(\Omega)}
\left\{
\frac14\mathbb E_\mu[q-E_C(d)]
+\sum_{e\in H}(\mathbb E_\mu s_e^C(d))_+
\right\}.
\end{aligned}}
\tag{R47.A6}
```

The last equality follows because maximizing separately over `0<=p_e<=1`
produces the positive part of each mean signed edge.  This dual is rigorous,
but it is the game in which one randomizes/fractionalizes the replacement
before a common response is chosen.  It is not the adaptive integral game
(R47.A5).

The gap can be maximal.  If `H` is the full edge set, `p_e=1/2` makes the
zero matrix, so

```math
\eta_{\rm frac}(C,H)=q/4,
```

while (R47.A5) is nonpositive because every integral corner is a signing of
cap at least `q`.  For `C=B_S`, the full block contains the reverse flip which
returns `A`, so `eta_int(B_S,H)=0` exactly.  Any slack-to-replacement proof
must therefore control integral rounding and every migrated response; an LP
using only the low-slack incidences cannot do so.

There is also an exact interpretation of the selected responses.  If `d^T`
denotes vertex-switching the physical word of `d` on `T`, then (R47.A1)
gives

```math
\boxed{L_S^A(d)=E_{A^{E(T)}}(d^T).}
\tag{R47.A7}
```

Thus low complement slack merely says that the particular response `d^T` is
near the threshold for the small internal replacement `A^{E(T)}`.  It says
nothing about the other states which determine `Q(A^{E(T)})`, or about the
new maximizers after several such replacements are combined.

## 3. Exact `A_6,A_8,A_9` audit

Every positive-price scalar optimizer is Pareto-undominated in
`(R_2,alpha)`, so auditing all Pareto columns is a superset of auditing all
scalar optimizers.  Exhaustion gives:

| signing, selector size | Pareto `(R_2,alpha)` types | every active `(Gamma, Q(B_S)-q, Q(B_S)-L_S)` |
|:--|:--|:--|
| `A_6,m=5` | `(30,1/2)` | `(0,0,0)` |
| `A_8,m=6` | `(40,1/14)`, `(72,3/14)` | `(0,4,4)` |
| `A_9,m=7` | `(56,1/36)`, `(72,1/18)`, `(88,1/12)` | `(4,4,0)` |

Here `r=1` for `A_6` and `r=2` for `A_8,A_9`.  The latter exactly match the
universal bound `Gamma<=2r(r-1)=4`.  The two examples isolate opposite
failure modes:

- On `A_8`, every selected incidence is threshold-tight but lies four below
  the actual complement cap.  Small `Gamma` does not even make the selected
  state an approximate complement maximizer at the integer scale.
- On `A_9`, every selected incidence is an exact complement ground and has
  the maximal allowed slack four.  Reversing its one missing-pair edge lowers
  the cap from `28` to `q=24`, but not below `q`; equality is shielded by the
  parent ground face.

For a still more direct combination test, associate to each active
`m=n-2` incidence its missing-pair edge and flip the union of all those edges
in `A`.  Exact exhaustive results are:

| signing and Pareto type | selected energies after union flip | cap after union flip | oriented shield states with energy `>=q` |
|:--|:--|--:|--:|
| `A_8,(40,1/14)` | `16,16` | `28` | `10` |
| `A_8,(72,3/14)` | six copies of `8` | `24` | `14` |
| `A_9,(72,1/18)` | `24,24` | `32` | `22`--`27` |
| `A_9,(88,1/12)` | `20,20,20` | `36` | `10`--`12` |

Thus the natural common replacement can push every selected response to or
strictly below `q` while the true cap rises well above `q`.  New response
states, absent from the low-slack incidence list, shield the signing.  For
every audited Pareto type, exhaustive enumeration of all subsets of these
missing-pair edges has minimum cap exactly `q`, never below it.

The `A_6` star case and an additional audit of `A_8,A_9` at `m=n-1` verify
(R47.A4) on every active incidence.  This is finite confirmation of the
universal proof, not the source of it.

## Verification

`tmp/complement_slack_r47_check.py`:

1. verifies (R47.A1) exactly for every active incidence of every Pareto
   column in the displayed `A_6,A_8,A_9` cases;
2. exhausts complement caps and the three-part slack profiles;
3. exhausts the natural union replacements, all their sub-replacements, and
   the migrated shield states; and
4. verifies at `m=n-1` that every active incidence has zero slack.

It terminates with

```text
PASS complement_slack_r47_check
```

## Research judgment

The proposed minimality-based dichotomy is sharply false without an extra
full-face or integral-rounding hypothesis.  Near the diagonal it is
universally false at asymptotic scale, not merely defeated by small examples.
At fixed density, the exact surviving complement target remains unchanged:
one must separately prove project-scale column mass and a positive
conditional fraction of genuinely `Theta(q_n)`-slack incidences.  If such a
fraction is small, minimality supplies only adaptive replacement witnesses;
turning those into one strict all-response improvement is essentially the
missing integral congestion/rounding theorem, not a consequence of the
near-tight incidences themselves.
