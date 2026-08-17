# Independent audit of the multiscale partition composition note

Date: 2026-08-17.

Audited draft: `multiscale_partition_composition_audit.md` (the version with
the one-sided constant `8`).

## 0. Verdict

**REPAIR.**  MC.2--MC.4 are correct, including the cell-slice
subgaussian scale, the Frobenius/projection lower bound, the Boolean-anchor
corollary, and the fixed-ratio `n^(3/2)` normalization.  MC.1 is correct
after one necessary hypothesis is restored: the fixed partition must be
**equitable**, with every cell of size `floor(n/q)` or `ceil(n/q)`.  For an
arbitrary partition, the averaging proof does not control the size of the
selected cell, and the stated theorem is false even for one exact signing.

There are also two scope repairs worth making explicit.

1. In MC.3, the retained channel has dimension **at most** `d`; it can be
   smaller when the prescribed vectors `u_(.,c)` are dependent or zero.
2. The Boolean-basis conclusion MC.19 does not automatically apply to an
   arbitrary block-averaging projection from MC.3.  It applies when the
   retained subspace actually has an orthonormal Boolean basis whose exact
   bridge responses obey MC.18 (for example, suitable equal cells with a
   Hadamard basis).  MC.17 remains valid for every orthogonal projection.

| Claim | Verdict | Audit finding |
|---|---|---|
| MC.1 simultaneous support | **REPAIR** | The cost averaging and constant `8` pass, but the size claim requires an equitable partition. |
| MC.2 many-child obstruction | **PASS** | The clique calculation has cap and defect `Theta(n^(3/2))`; its weighted/non-near-minimizer scope is stated correctly. |
| MC.3 joint compiler | **PASS WITH WORDING REPAIR** | Slice conditioning gives proxy `O(s-d)` and the Boolean union bound gives the displayed scale. Say “at most `d`-dimensional.” |
| MC.4 general projection bound | **PASS** | Sharp `p=1` Khintchine plus row-norm truncation gives MC.17 exactly. |
| MC.4 Boolean-anchor corollary | **PASS WITH SCOPE** | MC.19--MC.21 and all powers of `n,s` are correct under the additional Boolean-basis hypothesis. |
| Fixed-ratio interpretation | **PASS** | At `s/n -> theta`, MC.20 is a fixed positive multiple of the parent scale `(n+s)^(3/2)`. |
| Novelty | **SCOPED** | MC.4 is fundamentally a Frobenius/Khintchine rank argument. Its useful addition is the balanced-anchor-to-large-Boolean-residual specialization, not a new general rank principle or an information lower bound. |

## 1. MC.1: the averaging is correct, but equitability is necessary

For child `r`, after orientation and switching, let

```math
L_(r,a)=\sum_(i\in J_a)\ell_i^(r),
\qquad R_(r,a)=Q_-(D^(r)[J_a]).
```

The ground-state row sums are nonnegative and sum to `2Q_r`.  The one-sided
partition budget gives

```math
\sum_aR_(r,a)\le Q_r.
```

Consequently

```math
\sum_a\{2L_(r,a)+4R_(r,a)\}
\le4Q_r+4Q_r=8Q_r.                              \tag{A.1}
```

Summing (A.1) over the children and selecting a minimum-cost cell gives

```math
\sum_r\{2L_(r,a)+4R_(r,a)\}
\le {8\sum_rQ_r\over q}.                       \tag{A.2}
```

Every summand is nonnegative.  The exact flip identity and
`-sum_(i,j in S)d_ij <= R_(r,a)` therefore give, for every child and every
`S subseteq J_a`, precisely

```math
\rho_rH_(A^(r))((x^(r))^S)
\ge Q_r-{8\sum_uQ_u\over q}.                   \tag{A.3}
```

Thus the one-sided conclusion and constant `8` are correct.  If every
partition cell has size `floor(n/q)` or `ceil(n/q)`, the selected cell also
has the claimed size.  Nothing in the averaging argument gives this for an
arbitrary partition.

Here is an exact-sign counterexample to MC.1 as presently worded.  Let `A`
be the all-positive signing of `K_n`, let `q=floor(sqrt n)`, and partition
the vertices into `q-1` singletons and one cell of size `n-q+1`.  The only
cell of size at least `floor(n/q)` is the large cell.  Its ground state is
`1`, and flipping a subset of about `n/2` vertices inside that cell leaves
total magnetization `O(1)`, so

```math
H_A(1^S)=O(n),
\qquad Q(A)-H_A(1^S)=Theta(n^2).                \tag{A.4}
```

The asserted allowance is only

```math
{8Q(A)\over q}=O(n^(3/2)),                     \tag{A.5}
```

so the large cell fails.  All low-cost singleton cells are too small.
Hence “fix one partition” must be replaced by “fix one equitable
partition,” or the size conclusion must be dropped.

With this repair, `Q_r<=Cn^(3/2)` and `t=o(q)` give normalized defect at
most `8Ct/q=o(1)`, exactly as claimed.  MC.2 correctly shows that a
`Theta(q)`-sized family cannot be synchronized by this budget argument
alone.  It is not a signing-specific obstruction because its matrices have
zero entries, a limitation the draft already records.

## 2. MC.3: independent derivation of the cell-slice scale

Fix one row and cell of size `h`, write

```math
a_j=\eta_j-\bar\eta,
\qquad \sum_ja_j=0,
```

and let `K=(h+u)/2`.  A uniform Boolean vector of sum `u` is obtained by
choosing a uniform `K`-subset and putting `+1` there.  Its residual cell
contribution is

```math
X_c=\sum_jb_ja_j=2\sum_(j\text{ selected})a_j. \tag{A.6}
```

It is centred.  Hoeffding comparison for sampling without replacement,
followed by Hoeffding's lemma for the population `a_j` (whose range has
length at most two), yields

```math
E\exp(\lambda X_c)\le\exp(C\lambda^2h).         \tag{A.7}
```

The estimate is uniform in the feasible slice, including the deterministic
extreme slices.  Singleton cells contribute zero.  Independence across
cells makes the full row residual subgaussian with proxy

```math
C\sum_(c:h_c>1)h_c
\le2C\sum_c(h_c-1)=2C(s-d).                    \tag{A.8}
```

Rows are independent.  Hence, for fixed Boolean `(z,eta)`,
`sum_i z_iX_i` has proxy `O(n(s-d))`.  A union bound over the
`2^(n+s)` pairs proves

```math
\max_(z,eta)|\sum_i z_iX_i|
\le C\sqrt{n(s-d)(n+s)}.                       \tag{A.9}
```

For fixed `eta`, maximizing in `z` is exactly the row `l_1` norm.  This
proves MC.12.  If `d=s`, every cell is a singleton and the residual is
identically zero; this is the harmless degenerate case of the formula.

Every block-constant Boolean endpoint lies in `range(P)` and is therefore
preserved exactly.  The response `BP eta` depends on `d` scalar cell
averages, but its actual linear dimension can be below `d`; “genuinely
`d`-dimensional” should be weakened to “at most `d`-dimensional” unless a
rank hypothesis on the prescribed `u_(.,c)` is added.  The construction is
joint in the relevant sense: cell residuals are summed before the rowwise
absolute value and before the global union bound.  It is not `d` separate
Lipschitz payments.

MC.3 is a real multichannel extension of the scalar rowwise
microcanonical compiler (LC.1/Theorem 21.66), but it uses the same
sampling-without-replacement mechanism rather than importing a new
concentration principle.

## 3. MC.4: projection and Boolean-anchor calculations

Let `r_i=b_i(I-P)`.  Orthogonal projection gives
`||r_i||_2<=||b_i||_2=sqrt s`.  For uniform Boolean `eta`, sharp real
Khintchine at exponent one and `u>=u^2/sqrt s` on `[0,sqrt s]` give

```math
\begin{aligned}
E_\eta||B(I-P)\eta||_1
&=\sum_iE|r_i^T\eta|\\
&\ge {1\over\sqrt2}\sum_i||r_i||_2\\
&\ge {1\over\sqrt{2s}}\sum_i||r_i||_2^2\\
&={||B(I-P)||_F^2\over\sqrt{2s}}.
\end{aligned}                                   \tag{A.10}
```

Since `P` and `I-P` are orthogonal right projections,

```math
||B(I-P)||_F^2=||B||_F^2-||BP||_F^2
=ns-||BP||_F^2.                                 \tag{A.11}
```

The maximum is at least the mean, proving MC.17 with exactly its displayed
constant.

For a Boolean orthonormal basis `u_a=eta^(a)/sqrt s` of `range(P)`, put
`g^(a)=B eta^(a)`.  Then

```math
||BP||_F^2
=\sum_a||Bu_a||_2^2
={1\over s}\sum_a||g^(a)||_2^2
\le {1\over s}\sum_a||g^(a)||_1^2
\le {dL^2n^2\over s}.                           \tag{A.12}
```

Substitution in (A.10)--(A.11) is exactly MC.19.  If `s/n -> theta>0`,
`L=O(1)`, and `d=o(n)`, its right side is

```math
(1-o(1))\sqrt{\theta/2}\,n^(3/2).              \tag{A.13}
```

In terms of the full parent order `N=n+s`, this is

```math
\left({\sqrt{\theta/2}\over(1+\theta)^(3/2)}-o(1)\right)N^(3/2),
                                                               \tag{A.14}
```

a fixed positive target-scale loss.  Driving this displayed obstruction to
`o(n^(3/2))` requires

```math
d\ge(1-o(1)){s^2\over L^2n},                   \tag{A.15}
```

so MC.21 is correct.  If the right side exceeds `s`, the conclusion is that
no full collection satisfying the balanced-anchor hypothesis exists; there
is no contradiction with `P=I`.

The additional Boolean-basis assumption is substantive.  A general
block-averaging subspace has a basis of normalized cell indicators, not
automatically a Boolean orthonormal basis.  For equal cells and an available
Hadamard system on the cell labels, Boolean block-constant basis vectors do
exist.  Otherwise MC.19 must not be cited as a direct corollary for the
specific MC.3 projection.  MC.17, or a separate estimate using the actual
vectors `u_(.,c)`, is the universally applicable statement.

## 4. Novelty and exact no-go scope

MC.4 is not merely the earlier bounded-operator rank theorem in different
notation, but it is built from the same classical resource accounting.

* Theorem 18.8/BR.1 uses the exact sign Frobenius mass plus an operator-norm
  ceiling to force linearly many scale-visible singular values.
* MC.4 uses exact sign Frobenius mass plus Khintchine and balanced responses
  on a Boolean basis to force a target-scale `infinity -> 1` residual after
  projection.

Neither statement implies the other without extra hypotheses: MC.4 assumes
no operator-norm bound, while BR.1 assumes no balanced Boolean anchor
responses.  The genuinely useful project-level addition in MC.4 is the
bridge

```text
few exact balanced Boolean anchors
    => little captured Frobenius mass
    => a large separately paid Boolean residual.
```

As a general matrix theorem, however, MC.4 is an elementary
Frobenius/Khintchine projection bound and should be described that way.  It
does not prove contextual information complexity: linear rank can have a
short algebraic presentation.  It does not lower-bound the actual parent
cap, since the endpoint exposing the residual may obtain cancellation from
the child quadratic energy and the shore energy.  It does not obstruct a
nonlinear quotient, an algebraically closed linear-size channel, or a joint
same-switch estimate evaluated before absolute values.  The scope paragraph
in the draft correctly preserves all of these escape routes.

After repairing MC.1's partition hypothesis and the two wording/scope
points above, the final implication graph is accurate: the note proves
finite-family support compatibility and a genuine multichannel one-shot
compiler, while establishing only a scoped no-go for sublinear balanced
linear projections with separately paid residuals.  It does not prove a
fixed-ratio recurrence or a signing-specific synchronization theorem.
