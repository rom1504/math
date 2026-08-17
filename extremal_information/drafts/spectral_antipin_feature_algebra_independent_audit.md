# Independent audit: spectral anti-pin feature algebra

**Audited files**

- `drafts/spectral_antipin_feature_algebra.md`;
- `experiments/verify_spectral_antipin_feature_algebra.py`.

**Verdict:** PASS WITH ONE REQUIRED WORDING/INEQUALITY REPAIR.  SA.1--SA.3
are correct with the stated constants and scope.  SA.4 gives an exact
Boolean cap for the repeated port, and a rigorous global Boolean upper bound
for the orthogonal-port cap.  It therefore proves that the cap gap is **at
least** `(2-sqrt(2))rn`; it does not prove that this is the exact gap.

## 1. SA.1 normalization and resolvent positivity

Let `u=s odot x`.  Because deleting the diagonal subtracts the Boolean
constant `tr(H)/2`, the trace-zero hypothesis gives exactly

```math
H_(A_s)(x)={1\over2}u^THu.                          \tag{A.1}
```

The spectral hypothesis gives the upper bound `|u^THu|/2<=rn/2`, while
`u=1` attains `rn/2` because `H1=r1`.  Thus

```math
Q(A_s)=rn/2                                         \tag{A.2}
```

with no missing calibration or factor of two.

Since `H` is symmetric and `||H||<=r`, every eigenvalue lies in `[-r,r]`.
For either sign,

```math
K_\sigma=2mI-\sigma H\succeq(2m-r)I\succ0.         \tag{A.3}
```

Thus both resolvents in SA.4 exist and the square completion is legitimate.
On the Boolean sphere,

```math
{\sigma\over2}u^THu+mw^Tu
=mn-{1\over2}u^TK_\sigma u+mw^Tu
\le mn+{m^2\over2}w^TK_\sigma^{-1}w.              \tag{A.4}
```

The normalization in SA.4 is precisely chosen so that the last term is
`mn Psi_sigma(w)`.  Condition SA.5 therefore yields

```math
R_\sigma(w)
\le mn+{rn\over2}-\delta mn.                       \tag{A.5}
```

This verifies SA.8.

On the diagonal, `x=s,y=1` simultaneously attains child energy `rn/2`,
bridge energy `mn`, and `H_C(1)=E_C`.  The separate absolute bounds give the
matching upper bound, so SA.7 is exact.  The two diagonal queries give
opposite gaps `delta mn`, and cap nonexpansiveness gives `d_C<=d_0`.
Finally `d_0<=Q(A_s)+Q(A_t)=rn`; hence

```math
d_C\ge\delta mn\ge{\delta m\over r}d_0.            \tag{A.6}
```

All constants in SA.9 pass.

The nonzero-trace paragraph is also correctly scoped.  Hollowing changes
every Boolean energy by the same scalar `tr(H)/2`, of magnitude at most
`n/2`.  It does not preserve the exact diagonal identity, but it perturbs
all cap inequalities by only `O(n)`, lower order when `r,m=Theta(sqrt n)`.

## 2. SA.2 constants

Under `H^2=r^2I` and `m=r`,

```math
(2rI\mp H)^{-1}={2rI\pm H\over3r^2}.               \tag{A.7}
```

For Boolean `w`, `||w||^2=n` and `w^THw=rnrho(w)`.  Substitution gives

```math
\Psi_+(w)={2+\rho(w)\over6},
\qquad
\Psi_-(w)={2-\rho(w)\over6}.                      \tag{A.8}
```

Since `r/(2m)=1/2`, the code condition `|rho|<=theta` gives

```math
\max_\sigma\Psi_\sigma
\le{2+\theta\over6}
={1\over2}-{1-\theta\over6}.                      \tag{A.9}
```

Thus `delta=(1-theta)/6` is exact.  The Hanson--Wright existence statement
has the right scales: a sign matrix has `||H||_F=n`, while the involution has
`||H||=r=sqrt(n)`.

## 3. SA.3 spherical relaxation and Gram formula

After independently optimizing the endpoint sign of each uncompleted port
and the outer absolute channel, the exact Boolean cap is the maximum of
SA.17 over `sigma` and `epsilon`.  No channel is omitted.

For fixed `(sigma,epsilon)`, put
`v=sum_a epsilon_a w_a`.  For every `a_0>r/2`,

```math
{\sigma\over2}u^THu+mv^Tu
\le a_0n+{m^2\over2}v^T(2a_0I-\sigma H)^{-1}v.    \tag{A.10}
```

This is a valid upper bound on every Boolean `u`.  It is also the exact
Lagrange dual of the maximization over the Euclidean sphere `||u||^2=n`.
Because `H^2=r^2I` and `tr H=0`, both eigenvalues `+-r` occur; the dual
domain really is `a_0>r/2`, with a boundary limit covering the hard case.

Set `a_0=alpha r`.  The inverse is

```math
(2\alpha rI-\sigma H)^{-1}
={2\alpha rI+\sigma H\over r^2(4\alpha^2-1)}.      \tag{A.11}
```

Using

```math
||v||^2=n\epsilon^TG\epsilon,
\qquad
v^THv=rn\epsilon^TR\epsilon                       \tag{A.12}
```

gives exactly SA.19.  Therefore SA.19 is simultaneously:

1. the exact spherical optimum for that channel; and
2. a genuine upper bound on the Boolean multi-port cap.

The theorem does not claim equality with the Boolean cap, and its explicit
warning about a fixed leading integrality gap is correct.

The state count also passes.  The two symmetric `l by l` Gram matrices have
`l(l+1)` scalar entries total; each normalized integer numerator needs
`O(log n)` bits.  Filling all missing edges among the `lm` auxiliary
vertices changes the Hamiltonian uniformly by at most

```math
{lm\choose2}=O_l(m^2)=O_l(n)                       \tag{A.13}
```

when `l` is fixed and `m=Theta(sqrt n)`.  This completion claim is only an
additive cap perturbation.  It does not say that endpoint optimizers or the
exact Gram formula survive the completion, and the draft scopes it that way.

## 4. SA.4 tensor construction

The vector `v_0` in the draft is balanced and satisfies
`H_16v_0=4v_0`; the verifier checks both facts exactly.  Tensoring gives

```math
H\mathbf1=r\mathbf1,
\qquad Hv=rv,
\qquad \mathbf1^Tv=0.                              \tag{A.14}
```

Thus both one-port objects have `G=R=(1)`.  For either top Boolean
eigenvector `w`, the exact one-port cap before completion is

```math
\max_u\{|H_A(u)|+r|w^Tu|\}={3\over2}rn,            \tag{A.15}
```

because the spectral/Cauchy upper bound is attained at `u=w`.

For two repeated ports, the exact cap is

```math
\mathcal Q(\mathbf1,\mathbf1)={5\over2}rn.         \tag{A.16}
```

For the orthogonal pair, each choice of field signs gives
`z=+-1+-v` with `||z||=sqrt(2n)`.  Hence, for every Boolean `u`,

```math
|H_A(u)|+r|z^Tu|
\le\left({1\over2}+\sqrt2\right)rn.               \tag{A.17}
```

This is a direct bound on the exact global Boolean cap, not merely an
evaluation of the spherical carrier SA.19.  Therefore

```math
\mathcal Q(\mathbf1,\mathbf1)
-\mathcal Q(\mathbf1,v)
\ge(2-\sqrt2)rn.                                   \tag{A.18}
```

### Required repair

The sentence before SA.28 currently says the gap “is”
`(2-sqrt(2))rn`, and SA.28 is displayed as an equality-like bare
expression.  SA.26 is only an upper bound, so the proved statement is
**at least** (A.18).  Replace that sentence and SA.28 by

```math
\mathcal Q(\mathbf1,\mathbf1)
-\mathcal Q(\mathbf1,v)
\ge(2-\sqrt2)rn=\Theta(n^{3/2}).                   \tag{A.19}
```

The distinction is real: at `n=16`, the verifier finds

```text
same-pair cap = 160,
orthogonal-pair cap = 96,
actual gap = 64,
```

whereas `(2-sqrt(2))rn` is approximately `37.49`.  The theorem's scalable
separation is stronger than the displayed lower bound, not equal to it.

An arbitrary public exact-sign completion costs only `O(n)` by (A.13), so it
cannot erase this leading gap.  No optimizer identity or spherical
attainment is needed for that conclusion.

## 5. Scope and verifier

The compositional conclusion is accurate after the inequality repair:
separate one-port self-Gram states omit cross entries, and those entries can
change the exact two-port Boolean cap at leading scale.  SA.4 does not prove
an `Omega(p^2)` lower bound for every representation, and the draft explicitly
declines that stronger claim.

The verifier runs successfully.  It checks:

- the SA.2 resolvent formulas on three vectors at `n=16`;
- the self and cross `(G,R)` matrices;
- exact one-port and two-port Boolean caps at `n=16`;
- the scalable tensor eigenvector identities at the next tensor order.

It does not test the general resolvent theorem SA.1 or numerically minimize
the trust-region formula SA.19.  Those are proved algebraically above.  Its
inequality assertion for SA.4 is correctly written as `same-cross >= ...`,
which is further evidence that only the prose/display in the draft needs
repair.

## 6. Final classification

| Item | Verdict | Repair |
|---|---|---|
| SA.1 normalization | PASS | none |
| resolvent positivity | PASS | none |
| SA.1 metric constants | PASS | none |
| SA.2 involution constants | PASS | none |
| SA.3 spherical formula | PASS | none |
| Boolean multi-port upper bound | PASS | none |
| `O_l(n)` exact-sign completion | PASS | retain fixed-`l` scope |
| SA.4 exact repeated cap | PASS | none |
| SA.4 orthogonal global cap bound | PASS | none |
| SA.4 claimed gap | PASS after repair | replace equality wording by `>=` |

The finite-port Gram algebra is mathematically sound and genuinely
generative.  Its failure of congruence is established at the exact Boolean
cap level, not only in a relaxation.
