# Exact conditional-entropy upper scheme, and a parity-relaxation falsifier

Date: 2026-09-06. Scope: the full asymmetric recursive-Hadamard Bellman
operator in `continued_convergence_recursive_orbit_bound_2026_09_06.md`.
The algebraic reductions below have been independently reconstructed by
the audit agent. No negative full-row upper certificate is claimed.

## 1. Exact entropy telescope

Write `f_r=B^r Phi_t`, and use only the safe global reversal and input-swap
averages. At a node, the signed pair law is `pi`, its orthogonal image is
`rho=(U,V)`, and both signed child marginals are already symmetric. Thus

`H(pi)=H(U)+H(V)-I(U;V)`.

Let `C_t=2H-F_t`, as in the convergence artifact. Summing the exact
finite-tree objective gives

`f_r(nu_root) = max { -H(nu_root)
  + 2^(-r-1) sum_leaves C_t(nu_leaf)
  - (1/2) sum_internal_nodes 2^(-depth) I(U;V) }`.                 (1)

The maximization is over the original linear pair/child-flow constraints.
There is no independent-relative-sign restriction at any internal node.

In integer coordinates, `U=i+j` and `V=i-j` always have the same parity.
Consequently `I(U;V)>=h(q_odd)`. On any declared interval `[a,b]` for
`q_odd`, the chord of the concave binary entropy is a lower bound. Using
these chords in (1) gives a concave upper relaxation: expose each leaf
`C_t` through its self-coupling entropy maximization, retain all linear
tree constraints, and maximize a concave function. This is a correct
box upper scheme, but it discards too much dependence.

## 2. An all-tilt obstruction to the entire parity-only relaxation

This obstruction concerns the relaxation obtained by replacing each full
mutual information in (1) by its EXACT shared-parity entropy. It therefore
also applies to every globally covering parity-chord relaxation.

Put `p=15/16`. Choose independent input coordinates at every node. At
depth `d`, every signed law is that of

`X_d = 2^(-d/2) (X_1+...+X_(2^d))`,

where the summands have the original symmetric ternary law `nu_p`.
Their scale has no effect on entropy. Let `H_d` denote this finite-law
entropy, and let `b_d` be the probability that its unscaled integer sum
is odd. The policy is admissible at every node. The identity signed
self-coupling gives `F_t(mu)<=H(mu)` and hence `C_t(mu)>=H(mu)` at
every `t>=0`. Therefore the parity-relaxed candidate exponent obeys

`p log 2 + relaxed_f_r + t(1-sqrt(p))`
` >= p log 2 - H_0 + H_r/2 - (1/2) sum_(d=1)^r h(b_d)
     + t(1-sqrt(p))`.                                        (2)

All laws in (2) are exact convolutions of the integer weights `(15,2,15)`
with denominator `32`. Exact outward rational logarithm enclosures give

| depth | tilt-independent lower constant in (2) |
| --- | --- |
| 1 | greater than 0.265080457775 |
| 2 | greater than 0.271143586529 |
| 3 | greater than 0.195876879873 |
| 4 | greater than 0.053895806507 |

In particular **shared parity alone cannot establish the desired
negative exponent at depth 3 or 4, for any positive tilt**. This is an
exact limitation of the relaxation, not a lower bound of this size for
the actual Bellman value. The omitted within-parity output dependence
is essential.

The checker uses the positive atanh logarithm series after dyadic range
reduction, including its explicit geometric remainder, and rounds
outward to rational denominator `10^12` before adding terms. No optimizer
or floating logarithm enters the asserted inequalities. Run

```
OPENBLAS_NUM_THREADS=1 .venv/bin/python \
  computations/continued_feedback_recursive_entropy_upper_2026_09_06.py \
  --parity-falsifier \
  --output computations/results/continued_feedback_recursive_parity_falsifier_2026_09_06.json
```

## 3. Preserve all dependence: exact concavity of the final recursive level

Let `w` be an absolute child law and `K_t` the folded Gaussian kernel.
Its leaf potential is exactly

`Phi_t(w) = (1/2) max_(eta with both marginals w)
              [H(eta)+<log K_t,eta>] - H(w)`.

At one parent let `q` be its signed pair law and `w_+,w_-` its two
absolute child laws. Set `C_1=f_1+H_signed(input)`. Exposing the two
transport plans gives

`C_1 = max (1/4) sum_(s=+,-) {
   [H_signed(q)-H(w_s)] + [H(eta_s)-H(w_s)]
   + <log K_t,eta_s> }`.                                     (3)

The first entropy difference is the conditional entropy of the signed
pair given its absolute output. The second is the conditional entropy
of the transport destination given its source. Both are concave jointly
in the displayed probability variables, since their marginals are fixed
linear maps. All constraints are linear. Thus (3) is a genuine exposed
concave maximization, including every asymmetric pair policy.

For orbit class masses `q_c`, class multiplicities `m_c`, and deterministic
absolute-output map `o(c)`, the first conditional entropy is explicitly

`-sum_c q_c log(q_c / w_(o(c))) + sum_c q_c log m_c`.

This is the expression used in the code; the zero class and folded
kernel normalization are retained.

Group the bottom level of the full depth-r tree using (3). The ONLY
remaining nonconcave terms are `-H(w)` for absolute input profiles at
depths `1,...,r-1`. The difference between signed and absolute entropy
is the linear function `(1-w_0) log 2`. The exact diamond support,
parity, multiplicities, and parent/child constraints are all still present.

## 4. A finite globally certifiable upper scheme and its error

For each component `x` of an intermediate profile, replace the convex
function `x log x` by its secant on a declared interval `[a,b]`. Retain
the interval constraints. The resulting objective is concave, so each
box has a standard entropy/transport dual upper certificate.

If `b-a<=h`, the pointwise secant error is at most `h/e`. To see this,
write `x=a+theta h`. The error as a function of `a>=0` decreases, since
its derivative is

`(1-theta) log a + theta log(a+h) - log(a+theta h) <= 0`.

At `a=0` the error is `-h theta log theta<=h/e`. This includes intervals
touching zero and requires no lower bound on any probability.

There are `2^d` intermediate nodes at depth `d`, each with weight
`2^(-d)` and alphabet size `2^d+1`. Therefore every box relaxation
overstates the exact objective at any feasible point by at most

`h/e sum_(d=1)^(r-1)(2^d+1)
 = h(2^r+r-3)/e`.                                           (4)

Consequently, after covering all feasible profiles by such boxes, the
largest exact box upper value lies between the true Bellman value and
that value plus (4). At depth 3 the bound is `8h/e`. This supplies an
explicit finite certificate scheme rather than a local-optimizer claim.
It does not prove that a small number of boxes suffices.

## 5. Bounded numerical attempt and stopping criterion

The companion program implements both relaxations with explicit conic
entropy variables. At `r=3,p=15/16,t=5`, the unboxed parity relaxation
has numerical candidate exponent `+1.51816`. Restricting all seven parity
masses to intervals of radius `.001` around that optimizer still gives
`+0.70273`. These numbers are diagnostic; Section 2 proves the stronger
all-tilt impossibility of this relaxation.

For the exact conditional-entropy/secant scheme, the unbranched numerical
upper relaxation is `+1.82810`. A deliberately bounded 100-split adaptive
probe reduces the largest numerical relaxed value only to `+0.84637`.
This is much too loose to justify extending the enumeration during this
campaign. The probe values are NOT interval-certified global bounds;
the script labels solver failures and rounding limitations explicitly.

Saved results:

- `computations/results/continued_feedback_recursive_parity_falsifier_2026_09_06.json`
- `computations/results/continued_feedback_recursive_entropy_upper_depth3_t5_2026_09_06.json`
- `computations/results/continued_feedback_recursive_exact_boxes_depth3_t5_2026_09_06.json`
- `computations/results/continued_feedback_recursive_exact_boxes_branch100_depth3_t5_2026_09_06.json`

The independent convergence track has additionally found rigorously
positive feasible-policy exponents at depth 3 for `t=4` and `t=6`.
Those are genuine Bellman obstructions, unlike the relaxation obstruction
in Section 2. Negative local-policy values at intermediate tilts remain
lower bounds and do not establish the required full-row upper estimate.
