# Polynomial fractional bridge audit

This report answers a bounded constructive question: can a direct polynomial
formula in the child matrices give a low-variance fractional bridge whose cap
is certified without enumerating all parent states?

## 1. A positive theorem for Boolean-saturating conference children

Let `S` be a symmetric signing of order `n` satisfying

```math
S^2=(n-1)I.                                           \tag{P1}
```

Take the second child to be `-S` and use the fractional bridge `C_0=S`.
This bridge has its `n` diagonal entries equal to zero and all other entries
in `{+1,-1}`, so

```math
V(C_0)=\sum_{i,j}(1-(C_0)_{ij}^2)=n.                 \tag{P2}
```

The fractional parent is

```math
P_0=\begin{pmatrix}S&S\\S&-S\end{pmatrix},
\qquad P_0^2=2(n-1)I.                                \tag{P3}
```

Therefore its Boolean cap is at most

```math
\frac{2n}{2}\lVert P_0\rVert_{op}
=n\sqrt{2(n-1)}.                                     \tag{P4}
```

Suppose in addition that `S` has a Boolean eigenvector and saturates its
spectral cap bound:

```math
p=\operatorname{cap}(S)=\frac n2\sqrt{n-1}.          \tag{P5}
```

Then the right side of (P4) is exactly the ideal equal-child target

```math
n\sqrt{2(n-1)}=2\sqrt2\,p
=(p^{2/3}+p^{2/3})^{3/2}.                            \tag{P6}
```

Thus `C_0` is an explicitly supplied ideal-cap fractional bridge with only
`n`, rather than `n^2`, units of rounding variance. Independently round its
diagonal zeros, or simply choose any diagonal sign matrix `D`. The integral
parent with bridge `S+D` obeys, by operator-norm perturbation,

```math
\operatorname{cap}\begin{pmatrix}S&S+D\\S+D&-S\end{pmatrix}
\le 2\sqrt2\,p+n.                                    \tag{P7}
```

The corresponding `2/3`-power defect is `O(sqrt(n))`, hence geometrically
summable. Equations (P1)--(P7) are a proved uniform composition theorem for
every Boolean-saturating symmetric conference signing. Square-field Paley
conferences from the ledger satisfy (P5), by their explicit Boolean
eigenvectors.

The retained state is genuinely simpler than full parent minimization:
(P1), a Boolean eigenvector certificate, and the displayed polynomial bridge
certify everything by matrix multiplication and one norm bound.

### Fixed-block symmetric-Hadamard extension

The mechanism is not limited to two blocks. Let `K` be a symmetric Hadamard
matrix of fixed order `t`, so `K` has sign entries and `K^2=tI`. The
fractional order-`tn` signing

```math
P_0=K\mathbin\otimes S                               \tag{P7a}
```

has diagonal child blocks `+S` or `-S` and satisfies

```math
P_0^2=t(n-1)I.                                       \tag{P7b}
```

Its only invalid off-diagonal entries are the zeros joining the same vertex
coordinate in two different blocks: there are exactly `binom(t,2)n` such
edges. Fill them arbitrarily by signs. After reordering vertices by their
within-block coordinate, the perturbation is a direct sum of `n` symmetric
zero-diagonal order-`t` sign matrices, so its operator norm is at most `t-1`.
The resulting integral signing `P` obeys

```math
\operatorname{cap}(P)
\le\frac{tn}{2}\left(\sqrt{t(n-1)}+t-1\right)
=t^{3/2}p+\frac{tn(t-1)}2.                           \tag{P7c}
```

The first term is exactly the ideal target for `t` equal children,

```math
t^{3/2}p=(t p^{2/3})^{3/2}.                          \tag{P7d}
```

For every fixed symmetric-Hadamard order `t`, (P7c) has `O(n)` energy defect
and `O(sqrt(n))` defect after taking the `2/3` power. The fractional variance
is only `binom(t,2)n=O(n)`. Sylvester matrices supply every fixed dyadic `t`.

This is a uniform multi-block composition theorem, but all children are
copies of one Boolean-saturating conference signing. It controls geometric
scaling inside that structured family; it still supplies no landing theorem
for the true minima.

## 2. Exact relation to the universal integral double

If `U` is a signed permutation with `T=-U^T S U`, then

```math
C=(I+S)U                                             \tag{P8}
```

is an integral sign bridge (the diagonal of `I+S` is `+1`). Direct
multiplication gives

```math
SC+CT=0,
\quad CC^T=(I+S)^2,
\quad C^TC=U^T(I+S)^2U.                              \tag{P9}
```

For a conference child, `(I+S)^2=nI+2S`; these are exactly the order-12
`6+6` identities after switching and permuting the right child. Hence the
apparently exceptional bridge in that witness is the universal double.

The fractional formula replaces the diagonal `I` in (P8) by zero. It gains
variance `V=n`, but its ideal cap conclusion requires the saturation
condition (P5).

## 3. Saturation is essential, not a spectral technicality

There is first a sharp obstruction that does not assume a polynomial formula.
Let `A,B` be any two order-`n` sign children of common cap at most `p`, and
let `C_0` be any fractional bridge with variance

```math
V=n^2-\lVert C_0\rVert_F^2.
```

The fractional parent `P` has

```math
\lVert P\rVert_F^2
=2n(n-1)+2(n^2-V)=4n^2-2n-2V.
```

Since `P` has order `2n`,

```math
\lVert P\rVert_{op}^2
\ge\frac{\lVert P\rVert_F^2}{2n}
=2n-1-\frac Vn.                                     \tag{P10*}
```

Certifying the ideal target `T=2\sqrt{2}p` by the ordinary operator-norm bound
`cap(P)<=n||P||op` therefore requires

```math
V\ge2n^2-n-\frac{8p^2}{n}.                           \tag{P10**}
```

If `p=(c+o(1))n^(3/2)`, any such certificate requires

```math
\frac{V}{n^2}\ge2-8c^2-o(1).                        \tag{P10***}
```

Thus **subquadratic variance plus an ideal spectral-norm certificate forces
`c -> 1/2`**. If `c<1/(2\sqrt{2})`, the right side exceeds one and even the
zero bridge cannot be certified at the ideal target by this norm. For a
Boolean-saturating conference child, (P10**) gives exactly `V>=n`, attained
by the zero-diagonal bridge in Section 1.

This proves that the positive theorem is extremal for ordinary spectral-norm
certification. A fractional bridge intended to land near a family with a
constant strictly below `1/2` needs a genuinely Boolean, state-weighted norm
or inequality; no choice of polynomial bridge can fix the Frobenius barrier.

For an arbitrary symmetric child `S`, (P3) still gives

```math
\lVert P_0\rVert_{op}=\sqrt2\lVert S\rVert_{op}.
```

The norm certificate is at most the ideal target precisely when

```math
\frac n2\lVert S\rVert_{op}\le\operatorname{cap}(S). \tag{P10}
```

The reverse inequality always holds by the spectral cap upper bound, so
(P10) demands equality. Thus this polynomial/norm proof applies exactly to
Boolean saturation and cannot cover generic low-cap children with spectral
slack.

At order six, the conference child in the exact `M_12` witness has cap `5`,
whereas its spectral upper bound is `3sqrt(5)=6.708...`. Its zero-diagonal
fractional double has exact cap `18`, above the ideal target
`10sqrt(2)=14.142...`; varying only the diagonal coefficient does not lower
that cap.

The order-8 cospectral collision from the constructive-family report makes
the loss sharper. The two children have the same characteristic polynomial
but caps `14` and `12`. Their zero-diagonal polynomial doubles have exact caps
`40` and `44`, respectively. Hence the full spectrum neither determines the
polynomial double's Boolean cap nor compensates for spectral slack.

These finite caps are exhaustive arithmetic, reproduced by
`computations/audit_constructive_family_obstructions.py`.

The same obstruction is quantitative for a scaled polynomial bridge
`C_0=alpha S` (with zero diagonal). Put

```math
rho=\frac{2\operatorname{cap}(S)}{n\lVert S\rVert_{op}}\le1.
```

The parent norm is `sqrt(1+alpha^2)||S||`. Certifying the ideal target by
this norm requires

```math
\alpha^2\le2\rho^2-1,\qquad
V(C_0)=n(n-1)(1-\alpha^2)+n.                         \tag{P10a}
```

Thus any fixed spectral-saturation gap `rho<=1-epsilon` forces
`V(C_0)=Omega(n^2)` in this entire scaled-polynomial family. If
`rho<1/sqrt(2)`, even the zero-bridge parent cannot be certified at the ideal
target by its operator norm, although its exact Boolean cap is trivially at
most `2 cap(S)`. This pinpoints why replacing Boolean margins by a norm is
not harmless for generic children.

## 4. Research judgment

The formula is a genuine bounded-state, summable-defect composition theorem
inside the Boolean-saturating conference family. It does not solve the
campaign's central landing problem: the same square-field family has
normalized cap tending to `1/2`, and the constructive-family audit shows that
all `o(n)` deletions retain that limit. The theorem should be retained as the
cleanest positive composition mechanism, but it gives no evidence that this
family lies within `o(n)` of `b_n`.

For non-saturating children, every norm-only use of `P_0` loses exactly the
spectral-slack gap in (P10). Any extension must use Boolean state information
beyond spectrum; the order-8 collision proves that this is not optional.

## 5. Relative switching/permutation does not randomize the universal double

Let `U` be any signed permutation, set `T=-U^T S U`, and take the covariant
bridge

```math
C=(S+D)U,                                             \tag{P11}
```

where `D` is diagonal (integral or fractional). Changing Boolean variables
from `y` to `z=Uy` removes `U` completely. If
`J={i:x_i=-z_i}`, direct expansion gives the exact energy formula

```math
H_S(x)-H_S(z)+x^T(S+D)z
=2H_S(x)-4H_{S[J]}(x_J)+\sum_i d_i x_i z_i.          \tag{P12}
```

Thus every signed permutation produces exactly the same maximization problem.
Randomizing `U` cannot turn the dangerous principal-restriction term into a
concentration problem: it only renames the `2^n` right-child states. This is
an exact invariance obstruction, not a failed finite experiment.

One can decouple the right child and bridge permutations, but then the
polynomial square identity is lost. To see the precise norm cost, keep the
zero-diagonal bridge `S`, put `R=W^T S W` for the independently conjugated
right child, and assume `S^2=R^2=r^2I`. For

```math
P_W=\begin{pmatrix}S&S\\S&-R\end{pmatrix},
\qquad J=S/r,\quad K=R/r,
```

exact multiplication gives

```math
\frac1{r^2}P_W^2=
\begin{pmatrix}2I&I-JK\\I-KJ&2I\end{pmatrix}.       \tag{P13}
```

Since `J,K` are symmetric orthogonal matrices, the off-diagonal blocks are
transposes. Therefore

```math
\lVert P_W\rVert_{op}
=r\sqrt{2+\lVert I-JK\rVert_{op}}.                  \tag{P14}
```

The ideal spectral constant `sqrt(2)r` is retained with a little-oh relative
error only if

```math
\lVert I-JK\rVert_{op}=o(1),
\quad\text{equivalently}\quad
\lVert S-R\rVert_{op}=o(r).                          \tag{P15}
```

A generic relative conjugation has no such guarantee; the trivial bound is
`||I-JK||<=2`, which raises the norm certificate from `sqrt(2)r` to as much
as `2r`, a leading-order loss. Adding a diagonal `D` changes the parent norm
by at most `1` and cannot repair a fixed loss in (P14).

Hence the two natural choices exhaust the immediate permutation idea:

1. make the bridge covariant, in which case (P12) proves exact invariance and
   there is no random concentration;
2. decouple it, in which case a norm proof needs the near-automorphism
   condition (P15), while a raw union bound still ranges over all `2^(2n)`
   parent states.

This does not rule out a new Boolean inequality for carefully chosen
non-automorphic `W`, but such an inequality must control the full joint state
geometry. Neither switching/permutation randomness nor the polynomial norm
identity supplies it.
