# Greedy augmented recoupling as a weighted-set potential

## Status

**Verified conditional reduction; no bound improvement is claimed.**  This
note gives an exact potential and terminal-margin identity for augmented
row-sign recoupling.  It then uses only the project-scale cap hypothesis and
the universal Rademacher tails of switched row fields to prove a diffuse-core
dichotomy:

- sublinear unmatched terminal cores imply an `o(n^(3/2))` expected greedy
  defect;
- any scalable counterexample must create simultaneous linear-sized
  unmatched cores on both selected shores with nonvanishing probability.

This is theorem-only work.  No computation was used.

All energies use doubled normalization

```math
Q(A)=\max_{z\in\{\pm1\}^n}|z^TAz|.               \tag{1}
```

## 1. Switched row-sign variables

Let `X` be uniform Boolean and switch

```math
B=\operatorname{Diag}(X)A\operatorname{Diag}(X),
\qquad \ell=B\mathbf1.                            \tag{2}
```

Use the zero-field convention `sign(0)=+1`, and put

```math
I=\{i:\ell_i\ge0\},\qquad J=\{j:\ell_j<0\}.
```

Write

```math
C=B[I],\quad D=B[J],\quad K=B[I,J],
\quad k=K\mathbf1_J,\quad h=K^T\mathbf1_I,       \tag{3}
```

and

```math
P=\mathbf1_I^TC\mathbf1_I,
\qquad R=\mathbf1_J^TD\mathbf1_J,
\qquad H=\mathbf1_I^TK\mathbf1_J.               \tag{4}
```

The two shore masses are

```math
A_I=\sum_{i\in I}\ell_i=P+H,
\qquad
A_J=-\sum_{j\in J}\ell_j=-(R+H).                \tag{5}
```

Both are nonnegative, and the row-sign response is

```math
X^TA\operatorname{sign}(AX)=P-R=A_I+A_J.         \tag{6}
```

If `PR>=0`, the existing two-shore theorem recouples (6) without loss.  If
`PR<0`, (5) forces

```math
\boxed{P>0>R.}                                   \tag{7}
```

Indeed, `P<0` would force `H>0` from `P+H>=0`, and then `R+H<0`
would force `R<0`; the case `R>0` is symmetric.  Thus there is only one hard
orientation.

## 2. Exact set-function form

Consider the branch anchored on `I`.  For the free shore `J`, set

```math
a_j=-\ell_j=-(D\mathbf1+h)_j>0.                  \tag{8}
```

For `S subset J`, define

```math
\boxed{
\Phi_J(S)=a(S)+2e_D(S),
\quad
a(S)=\sum_{j\in S}a_j,
\quad
e_D(S)=\sum_{\{i,j\}\subset S}d_{ij}.}          \tag{9}
```

Let `r^S=1-2 1_S`.  A direct expansion gives

```math
(r^S)^TDr^S+2h^Tr^S
=R+2H+4\Phi_J(S).                                \tag{10}
```

The collapsed-coordinate spin is redundant: replacing `(r,t)` by `tr`
leaves `r^TDr` unchanged and turns `2t h^Tr` into `2h^T(tr)`.  Consequently
the exact aligned augmented cap is the maximum of (10).  Its shortfall from
the opposite witnessed shore energy `-R` is

```math
\boxed{
\delta_J^{\rm exact}
=4\left[{A_J\over2}-\max_{S\subset J}\Phi_J(S)\right]_+.} \tag{11}
```

The other branch has exactly the same form after replacing

```math
(D,h,R,a,J)\quad\hbox{by}\quad
(-C,-k,-P,\ell_I,I).                              \tag{12}
```

Thus augmented recoupling is two instances of one nonnegative-vertex-weight
quadratic set problem (the `J`-branch weights are strictly positive; the
symmetric `I` branch can contain zero row fields).  This is not a renaming of
the full cap: each instance has only one shore free, but it retains that
shore's complete signed edge set.

## 3. Terminal-margin identity

After eliminating the collapsed spin, the implemented augmented
best-improvement ascent becomes monotone ascent for `Phi_J` using singleton
toggles when a free coordinate flips and the whole-set complementation
`S -> J setminus S` when the collapsed coordinate flips.  Its terminal set
is nevertheless singleton-stable.  More generally, let `S` be any
singleton-stable terminal set, and put `T=J setminus S`.  Define its removal
and addition margins by

```math
\begin{aligned}
p_i&=a_i+2\sum_{u\in S\setminus\{i\}}d_{iu}
&& (i\in S),\\
q_j&=-a_j-2\sum_{u\in S}d_{ju}
&& (j\in T).                                     \tag{13}
\end{aligned}
```

One-flip stability is exactly

```math
p_i\ge0\quad(i\in S),
\qquad q_j\ge0\quad(j\in T).                    \tag{14}
```

Put `p(S)=sum_(i in S)p_i`.  Since

```math
p(S)=a(S)+4e_D(S),                               \tag{15}
```

the greedy certificate shortfall has the exact form

```math
\begin{aligned}
\delta_J^{\rm gr}
&=4\left[{A_J\over2}-\Phi_J(S)\right]_+\\
&=\boxed{2\,[a(T)-p(S)]_+.}                     \tag{16}
\end{aligned}
```

This is the desired potential invariant.  It is path-independent: the
complete history of flips disappears, leaving only original row-field mass
on the terminal outside and removal-stability margin on the terminal inside.

It also has a local-field interpretation.  If

```math
L=Dr^S+h,
```

then (8) and (13) give

```math
L_i=-p_i\quad(i\in S),
\qquad L_j=q_j\quad(j\in T).                     \tag{17}
```

Thus `p(S)` is exactly the terminal absolute local-field mass carried by the
selected side.  A failure is not caused merely by small terminal energy; it
requires the unselected *original* negative-row mass to exceed this selected
terminal stability mass.

The symmetric construction (12) gives `delta_I^gr`.  The combined greedy
recoupling defect is

```math
\Delta_{\rm gr}=0\quad(PR\ge0),
\qquad
\Delta_{\rm gr}=\min\{\delta_I^{\rm gr},
                       \delta_J^{\rm gr}\}
\quad(PR<0).                                     \tag{18}
```

## 4. Unmatched cores

When `a(T)>p(S)`, regard the scalar margin `p(S)` as payment against the
weights `{a_j:j in T}`.  Define the unmatched-core number `kappa_J` to be the
smallest integer `k` for which the `k` largest weights in `T` have total at
least

```math
u_J=a(T)-p(S)={\delta_J^{\rm gr}\over2}.          \tag{19}
```

Set `kappa_J=0` when `u_J<=0`.  Equivalently, after distributing total
payment `p(S)` arbitrarily among the outside weights, `kappa_J` is the
smallest possible support size of the residual unpaid mass.  Define
`kappa_I` by (12), and put

```math
\kappa_* = \min\{\kappa_I,\kappa_J\}.            \tag{20}
```

This definition applies in the hard branch; set `kappa_*=0` when `PR>=0`.

The use of the minimum is important: the augmented theorem pays only the
better of the two shore shortfalls.

## 5. Universal order-statistic bound for row fields

For every fixed signing and every coordinate `i`, the switched field

```math
\ell_i=X_i(AX)_i
=\sum_{j\ne i}a_{ij}X_iX_j                       \tag{21}
```

is distributed as a sum of `n-1` independent Rademachers.  Hence

```math
\Pr\{|\ell_i|\ge t\}
\le2\exp\left(-{t^2\over2(n-1)}\right).          \tag{22}
```

No independence between different rows is asserted or needed.  Let `L_k`
be the sum of the `k` largest values among `|ell_1|,...,|ell_n|`.  Integrating
(22) and summing over coordinates proves, for a universal constant `C`,

```math
\boxed{
\mathbb E L_k
\le C k\sqrt{n\log(2en/k)}}
\qquad(1\le k\le n).                             \tag{23}
```

For completeness, for any threshold `t_0`,

```math
L_k\le kt_0+\sum_i(|\ell_i|-t_0)_+.
```

Take expectations, use (22), and choose
`t_0=sqrt(2(n-1)log(2n/k))`; the Gaussian tail integral gives (23), after
enlarging `C` to cover `k` near `n`.

## 6. Diffuse-core theorem under project-scale cap

Assume a sequence of signings satisfies

```math
Q(A_n)\le K n^{3/2}                              \tag{24}
```

with fixed `K`.  Initialize each augmented branch by the field-aligned spin
and run its monotone best-improvement ascent (singleton moves plus the
possible whole-set complementation described above).  Then, for every
`0<alpha<=1`,

```math
\boxed{
{\mathbb E\Delta_{\rm gr}\over n^{3/2}}
\le
C'\alpha\sqrt{\log(2e/\alpha)}
+K\Pr\{\kappa_*>\alpha n\}.}                   \tag{25}
```

#### Proof

On `{kappa_*<=alpha n}`, choose the shore realizing `kappa_*`.  By
(16), (19), and the definition of `L_k` (with `L_0=0`),

```math
\Delta_{\rm gr}\le2L_{\lfloor\alpha n\rfloor},
```

and (23) supplies the first term of (25).

For the complementary event, each witnessed principal energy has magnitude
at most `Q(A_n)`: randomly completing a fixed principal spin makes all other
terms mean zero, so the full cap dominates its absolute expected energy.
The augmented ascent terminates at a coordinate-stable point, whose target
energy is a sum of absolute local fields and is therefore nonnegative.
Consequently each branch defect is at most its opposite witnessed principal
energy, and

```math
\delta_I^{\rm gr},\delta_J^{\rm gr}\le Q(A_n).
```

This proves the second term in (25).

### Corollary 6.1 (conditional greedy theorem)

If

```math
{\kappa_*\over n}\longrightarrow0
\quad\hbox{in probability},                     \tag{26}
```

uniformly over the signings in (24), then

```math
\boxed{\mathbb E\Delta_{\rm gr}=o(n^{3/2}).}    \tag{27}
```

Indeed, first fix `alpha`, use (26) in (25), and then send `alpha` to zero.

### Corollary 6.2 (necessary scalable obstruction)

Conversely, suppose along a subsequence

```math
\mathbb E\Delta_{\rm gr}\ge\varepsilon n^{3/2}. \tag{28}
```

Choose `alpha=alpha(epsilon,K)>0` so that the first term of (25) is at most
`epsilon/2`.  Then

```math
\boxed{
\Pr\{\kappa_I>\alpha n\ \hbox{and}\
       \kappa_J>\alpha n\}
\ge {\varepsilon\over2K}}                       \tag{29}
```

along a further subsequence (with harmless adjustment of the universal
constant).

Thus isolated large fields, a sublinear exceptional shore, or one bad branch
cannot falsify greedy augmented recoupling.  A true counterexample must
produce, with nonvanishing probability, **two simultaneous diffuse stable
cores**: on both shores a linear number of original row fields remain unpaid
by the terminal removal margins.

## 7. Research consequence

The desired theorem has been reduced to a concrete dynamical statement that
is weaker than evaluating either restricted cap:

> Prove that field-initialized best-improvement ascent makes the minimum
> unmatched-core number `kappa_*` sublinear in probability, or directly
> prove that its unmatched mass is `o(n^(3/2))` in expectation.

The cap hypothesis is used only to control the rare diffuse-core event; the
small-core contribution is handled uniformly by the exact Rademacher law.
If this conditional statement is false, (29) is the required form of a
scalable obstruction.  It must be diffuse on both shores and survive the
full row-sign law, rather than being a single bad principal block or a rare
large-field configuration.
