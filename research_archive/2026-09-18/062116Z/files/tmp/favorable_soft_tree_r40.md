# Wave 40: favorable distance transforms on completion trees

## 1. Outcome

There is an exact favorable soft potential which makes every descendant flip
admissible and removes the ancestor-multiplicity loss:

```math
\pi_S([y])=d_{\rm pr}([y],F_S),
\tag{R40.1}
```

where `F_S` is the actual favorable projective fiber.  When constrained
terminals are leaves (unconstrained full-spin Steiner nodes are allowed), the
soft problem

```math
\min_{\mathbf y}\left[D_T(z;\mathbf y)+\sum_v\pi_{S_v}(y^v)\right]
\tag{R40.2}
```

is **exactly equal** to the hard favorable completion problem.  Nearest-fiber
projection charges each terminal once.  The min-plus dynamic program is thus
well defined on a flip-closed domain, and its minimizer can be chosen
favorable.

This is a positive exact continuation of (10.1099), but not a new estimate.
Comparing (R40.2) with the common-root scaffold gives exactly the center-star
total cost (10.926); minimizing over arbitrary scaffolds gives exactly the
original terminal min-cut (10.923).  A finite-temperature log partition also
retains, rather than eliminates, the unresolved cross-fiber boundary
compatibility.  The smallest obstruction occurs already for an exact
order-four minimizer and two actual `t=0` favorable fibers.

## 2. Leaf Lipschitz lemma

Let `T` contain a fixed full projective root `[z]`, any number of
unconstrained full-projective Steiner nodes, and constrained terminal leaves
`v`.  Leaf `v` has selector `S_v` and prescribed projective word `[y^v]`.
Let `D_T(z;\mathbf y)` be the minimum total projective-Hamming edge length of
full completions.

If only leaf `v` changes from `[y]` to `[f]`, then

```math
\left|D_T(z;\mathbf y^{v\leftarrow f})-D_T(z;\mathbf y)\right|
\le d_{\rm pr}([y],[f]).
\tag{R40.3}
```

Indeed, take optimal full completions for `y`.  Orient `f` relative to the
old leaf restriction so that exactly `d_pr(y,f)` specified coordinates
change, flip only those coordinates in the full leaf completion, and leave
all outside coordinates and all other nodes fixed.  Since the leaf has one
incident edge, its length rises by at most that Hamming distance.  Reverse
the roles of `y,f` for the other direction.

Terminals can be made leaves by allowing unconstrained Steiner completion
nodes.  If the final forest formalism insists on a tree whose vertices are
terminals, a projective-metric terminal MST has cost at most twice the Steiner
cost, so this changes only the absolute constant.

## 3. Exact favorable distance-transform theorem

For nonempty favorable fibers `F_v=F_(S_v)`, define (R40.1).  Then

```math
\boxed{
\min_{\mathbf y}
\left[D_T(z;\mathbf y)+\sum_v d_{\rm pr}(y^v,F_v)\right]
=
\min_{\mathbf f\in\prod_vF_v}D_T(z;\mathbf f).
}
\tag{R40.4}
```

All labels on the left range over the full projective cubes, so the domain is
closed under every descendant coordinate flip.

**Proof with both directions.**  Denote the left and right sides by `L` and
`H`.  Taking `y^v=f^v in F_v` in the soft problem makes every potential zero,
so `L<=H`.  Conversely, fix arbitrary `y`.  For every leaf choose a nearest
`f^v in F_v`.  Change the leaves one at a time.  Equation (R40.3) gives

```math
D_T(z;\mathbf f)
\le D_T(z;\mathbf y)+\sum_vd_{\rm pr}(y^v,F_v).
\tag{R40.5}
```

The left side is at least `H`.  Taking the infimum over `y` gives `H<=L`,
proving equality.  Moreover, projecting any soft minimizer by this argument
produces a favorable soft minimizer.  Thus canonical joint selection, not an
arbitrary preselection of grounds, is built into (R40.4).

Choose such a favorable minimizer and an optimal completion with edge
disagreement sets `U_e`.  The exact descendant-flip Euler inequality becomes

```math
\boxed{
|U_e|\le
\sum_{v\in H_e}
d_{\rm pr}\!\left((y^v)^{U_e\cap S_v},F_v\right),
}
\tag{R40.6}
```

because `pi_v(y^v)=0`.  This is (10.1099) for the explicit potential (R40.1).
Summing (R40.6) would still reuse leaves along ancestors.  The point of
(R40.4)--(R40.5) is that global nearest-fiber projection supplies the correct
one-charge argument instead.

## 4. Edge-conditional min-plus normalization

The theorem has an exact descendant dynamic program.  For an unconstrained
node `u`, let `X_u` be its full projective spin; if `u` is a terminal leaf,
put `pi_u(X_u)=d_pr((X_u)_(S_u),F_u)`, and otherwise put `pi_u=0`.  Given a
parent state `a`, define recursively

```math
M_u(a)=\min_{x\in\{\pm1\}^n/\{\pm1\}}
\left\{d_{\rm pr}(a,x)+\pi_u(x)
+\sum_{w:\,w\text{ child of }u}M_w(x)\right\}.
\tag{R40.7}
```

At the fixed root, sum the child messages.  Every terminal potential occurs
in exactly one leaf message.  Standard tree elimination shows that the root
value is (R40.2), hence by (R40.4) the exact hard favorable forest cost.  This
is the requested edge/descendant conditional normalization; it does not hide
an ancestor multiplicity.

## 5. Scaffold bound and exact project-scale hypothesis

Equation (R40.4) gives, for **every** unconstrained scaffold `bar y`,

```math
\boxed{
\min_{\mathbf f\in\prod_vF_v}D_T(z;\mathbf f)
\le D_T(z;\bar{\mathbf y})
+\sum_vd_{\rm pr}(\bar y^v,F_v).
}
\tag{R40.8}
```

If `bar y^v=z|_(S_v)`, its completion cost is zero and (R40.8) is

```math
\min_{\mathbf f\in\prod_vF_v}D_T(z;\mathbf f)
\le\sum_va_z(S_v).
\tag{R40.9}
```

This is precisely the total-cost center-star criterion (10.926).  The
equal-radius exceptional-center statement (10.967) is a convenient stronger
uniform-law way to force (R40.9).  Allowing a nonconstant scaffold can be
strictly better through free-coordinate screening, but minimizing the right
side of (R40.8) over all scaffolds returns equality (R40.4), hence exactly the
terminal min-cut target (10.923), not a new bound.

At `L_0=n^(3/4-c_0)`, `k_0=Theta(L_0/log n)`, and group size
`s=Theta(n^(1/4+c_0-eta)log n)`, the exact additional theorem needed is:

> Against every selector law and on an event of probability
> `exp{-O(rL_0)}`, partition the batch into at most `n^eta` groups and, in
> every group, find a low-row root, a leaf-Steiner tree, and an unconstrained
> scaffold satisfying
> ```math
> D_T(z;\bar{\mathbf y})+
> \sum_vd_{\rm pr}(\bar y^v,F_(S_v))=O(k_0).
> \tag{R40.10}
> ```

Then (R40.8) gives favorable forest cost `O(k_0)` and closes the established
restriction chain.  For (R40.10) to be more than a restatement, the scaffold
must come from a canonical soft measure or another actual-minimizer theorem.
No consequence currently recorded in the ledger bounds its distance term.

## 6. Finite-temperature audit

Replace each minimum in (R40.7) by the normalized soft minimum

```math
-\tau\log\left(|\mathcal X|^{-1}
\sum_{x\in\mathcal X}e^{-E(x)/\tau}\right).
\tag{R40.11}
```

This preserves full support and descendant-flip admissibility.  But eliminating
one full projective spin has deterministic approximation error at most
`tau(n-1)log 2`; over `O(s)` states the error is `O(tau n s)`.  Retaining an
`O(k_0)` conclusion without an additional entropy theorem therefore requires

```math
\tau=O\!\left(\frac{k_0}{ns}\right)
=O\!\left(
\frac{n^{-1/2-2c_0+\eta}}{(\log n)^2}
\right)
\tag{R40.12}
```

at the smallest group size.  Thus generic finite-temperature smoothing is
forced almost to the min-plus limit.  It can improve the route only if exact
minimality supplies a much smaller **conditional** entropy or a direct lower
bound on the normalized compatibility partition.  Local terminal
normalizations alone do not do this.

The boundary obstruction is already present at the smallest possible order.
Consider the exact order-four minimizer

```math
A=\begin{pmatrix}
0&-1&-1&-1\\
-1&0&-1&-1\\
-1&-1&0&1\\
-1&-1&1&0
\end{pmatrix},
\qquad q_4=8.
\tag{R40.13}
```

For `S={0,1,2}` and `T={0,2,3}`, both principal norms are `6`.  Their unique
absolute-ground projective words are respectively

```math
[+,+,+],\qquad[+,-,-].
\tag{R40.14}
```

Since

```math
B_{4,3}=3\sqrt3-4\in(0,4)
\tag{R40.15}
```

and the only child deficits are `0,4`, these are the actual `t=0` favorable
fibers.  They disagree projectively on their two-coordinate overlap and their
completion cylinders have distance one.  The low-row root `z=(+,+,+,+)` has
row-square `20<=24`, hits the first fiber, and is one flip from the second.
Thus the hard cost and the distance-transform cost both equal one.

If each singleton fiber carries its normalized point mass, both local log
normalizers are zero, but the global compatibility partition is

```math
Z_\lambda=e^{-\lambda}.
\tag{R40.16}
```

Therefore (R40.13)--(R40.16) falsify only a claim that local/descendant
normalizations cancel the boundary compatibility cost.  They do **not**
falsify a correct finite-temperature DP: that DP retains the nonzero root
message, which is precisely the unresolved cross-selector quantity.  Orders
at most three cannot give this example: every proper two-vertex principal
signing has all projective words as absolute grounds, so order four is
minimal.

## 7. Verification and recommendation

`tmp/favorable_soft_tree_r40_check.py` exhaustively verifies `q_4=8`, the two
actual favorable fibers, their cylinder distance one, the exact hard/soft
identity, the low-row star scaffold, and the retained log-partition boundary
cost.

The favorable distance transform is the clean canonical continuation of
(10.1099).  It should be recorded as an exact equivalence and used to prevent
future ancestor overcounting.  It does not by itself raise terminal min-cut
above its current fallback rank: progress now requires a genuinely
minimizer-specific proof of the soft-scaffold event (R40.10), a conditional
entropy bound strong enough to beat (R40.12), or the already identified
exceptional-center estimate.  Merely replacing min-plus by log-sum-exp does
not supply that input.
