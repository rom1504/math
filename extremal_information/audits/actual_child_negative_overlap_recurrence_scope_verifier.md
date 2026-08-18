# Adversarial verification of the negative-overlap recurrence scope

Status: **PASS after one quantifier correction**.  The algebra, inequality
directions, power-saving propagation, and the main conclusion are correct.
Equation (RC.9), however, is not by itself a sufficient all-order criterion
unless `epsilon_N` already denotes a monotone/dyadic-shell envelope.

## 1. Algebra and directions

With

```math
G=E_UL-V_\lambda^{\rm row},\qquad
I=\lambda(V_\lambda^{\rm row}-V_\lambda),\qquad
\Delta=V_\lambda^{\rm row}-T,
```

direct substitution gives

```math
\Delta=(V_\lambda-T)+I/\lambda=(E_UL-T)-G.
```

All signs are correct.  The fair law is an admissible row-product law, so
`G>=0`; the full Gibbs variational class contains the row-product class, so
`I>=0`.  Theorem 37.50 therefore yields exactly

```math
G\le C_{\rm LS}\lambda t^2mn\widehat\rho_N^-(\lambda),
\qquad
I\le C_{\rm LS}\lambda^2t^2mn\widehat\rho_N^-(\lambda).
```

Under `V_lambda<=T+E`, division of the second inequality by `lambda` gives

```math
\Delta\le E+C_{\rm LS}\lambda t^2mn\widehat\rho_N^-(\lambda).
```

Since `t^2=beta^2/N`, this is precisely (RC.7).  On comparable splits,
`t^2mn=Theta_beta(N)`, so overlap decay `O(N^(-alpha))` and target error
`O(N^(1-gamma))` give `O(N^(1-min(alpha,gamma)))`.  No power or factor of
`lambda` is missing.

The basin statement needs the standard harmless qualification from Theorem
37.19: choose the slack `a>Delta` with `h+a>0`.  Target reach supplies this at
the same rate because `V_lambda>=L(0)` implies
`(-h)_+=(L(0)-T)_+<=E`.  Thus a slack controlling `Delta_+`, `(-h)_+`, and one
subleading positive term preserves the stated power accuracy and the
`exp[-O(N)]` mass scale.

## 2. Recurrence check

The power-saving assertion is correct provided the estimate is uniform over
the comparable splits used by the archived balanced-tree theorem.  A defect
`C N^(1-delta)` contributes `O((2^j k)^(-delta))` per vertex at level `j`,
whose sum is `O(k^(-delta))`.

The more general display (RC.9) needs a quantifier correction.  For an
arbitrary nonregular sequence `epsilon_N`, the point samples

```math
\sum_{j\ge0}\epsilon_{2^j k}
```

do not control the nearby orders that occur when leaves have sizes `k` and
`k+1` or when a general target order is merged.  For example, set
`epsilon_(2^r+1)=1/r` and set all other `epsilon_N` to zero.  Every exact
dyadic ray meets at most one nonzero odd index, so the displayed sum tends to
zero as `k` tends to infinity.  But every sufficiently large dyadic shell
contains such an index and the shell suprema have a harmonic, hence
nonsummable, tail.

A sufficient all-order formulation is

```math
\lim_{k\to\infty}\sum_{j\ge0}
 \sup_{2^jk\le N\le 2^{j+1}k}\epsilon_N=0,
```

up to harmless fixed changes of the shell endpoints dictated by the balanced
window.  The original (RC.9) is valid if `epsilon_N` is nonincreasing, or if
it is explicitly defined to be this shell envelope.  This issue does not
affect the power-law consequence (RC.8).

## 3. Scalar-recurrence conclusion

For positive `lambda`, the soft minimum satisfies

```math
\min_B L(B)\le-\lambda^{-1}\log E_Ue^{-\lambda L(B)}=V_\lambda.
```

Parent minimization gives `P_N(beta)<=min_B L(B)`.  Hence, when `T_N` is the
declared child recurrence right-hand side, target reach already gives

```math
P_N(\beta)\le T_N+E_N.
```

Thus the scope note's main verdict is correct: overlap controls product
approximation and basin abundance only *after* the scalar recurrence has
already been supplied.  Also, because both overlap statements are upper
bounds, positive overlap cannot imply positive `G` or `I`; the row-product
rank-one example with `I=0` is a decisive witness.

**Disposition.**  Correct (RC.9) as above, or state its missing regularity
hypothesis.  With that qualification, the source note passes adversarial
verification and its Level-5/Level-6 conclusion is unchanged.
