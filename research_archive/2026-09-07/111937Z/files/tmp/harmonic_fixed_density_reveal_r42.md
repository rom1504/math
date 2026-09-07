# Wave 42: fixed-density reveal identity and an exact-minimizer pointwise wall

## 1. Normalizations

For every retained set `S` define

```math
G_S(d)=2^{-|S^c|}K_{\beta,S}(d[S]),
\qquad F_S=\log G_S,
\qquad b_S=e^{-F_S}.
```

At size `m`, use the unnormalized sum and uniform average

```math
Z_m(d)=\sum_{|S|=m}b_S(d),
\qquad
U_m(d)=\binom nm^{-1}Z_m(d).
```

The posterior selector law is `pi_m(S|d)=b_S(d)/Z_m(d)`.  For a vertex
edge `d^+,d^-` which flips `i`, put

```math
N_{m,i}(e)=\sum_{\substack{|S|=m\\i\notin S}}b_S(e[S]).
```

Every summand is unchanged across the edge, so

```math
r_i(d^\pm)=\frac{N_{m,i}(e)}{Z_m(d^\pm)},
\qquad
\chi_i=\log\frac{U_m(d^-)}{U_m(d^+)}
=\log\frac{r_i(d^+)}{r_i(d^-)}.
```

The fixed-size identity remains `sum_i r_i=n-m`.

## 2. Exact outside-spin reveal

Take `S` not containing `i`, put `R=S union {i}`, and let `d^+,d^-` be
the two restrictions to `R`.  The composition identity (10.916) is

```math
G_S
=\frac12\{e^{\beta a}G_R(d^+)+e^{-\beta a}G_R(d^-)\}.
```

Define the coarse half-log-odds

```math
B_{S,i}
=\beta a+\frac12\{F_R(d^+)-F_R(d^-)\}.
```

Factoring the geometric mean gives the exact reciprocal identity

```math
\boxed{
b_S(e)=\sqrt{b_R(d^+)b_R(d^-)}\,\operatorname{sech}B_{S,i}(e).
}
```

Under the external Gibbs completion conditioned on the retained boundary,
the revealed spin has probabilities

```math
Q_{S,i}(\pm\mid e)=\frac{e^{\pm B_{S,i}}}{2\cosh B_{S,i}}.
```

Consequently

```math
D(U_{\{\pm\}}\Vert Q_{S,i})=\log\cosh B_{S,i},
\qquad
\operatorname{sech}B_{S,i}
=e^{-D(U\Vert Q_{S,i})}.
```

Iterating along any order of the outside spins proves the exact reverse-KL
reveal chain

```math
F_S(y)=\sum_{\ell=1}^{|S^c|}
\mathbb E_{\text{uniform earlier reveals}}
\log\cosh B_\ell.
```

Writing `m_ell=tanh B_ell`, the elementary inequality
`log cosh B >= (tanh B)^2/2` gives the reverse-KL square-function bound

```math
\sum_\ell\mathbb E_Um_\ell^2\le2F_S.
```

The expectation here is under **uniform** earlier reveals.  Conditional
magnetization is a martingale under the Gibbs completion law, not under this
uniform law, so ordinary martingale orthogonality cannot be inserted without
a change-of-measure cost.

## 3. Exact affinity/cross-level factorization

At size `m+1`, condition the selector to contain `i` and set

```math
Z_{m+1,i}^\pm
=\sum_{\substack{|R|=m+1\\i\in R}}b_R(d^\pm),
\qquad
\alpha_i^\pm(R)=\frac{b_R(d^\pm)}{Z_{m+1,i}^\pm}.
```

Their Hellinger affinity and geometric-mean selector law are

```math
H_i=\sum_{R\ni i}\sqrt{\alpha_i^+(R)\alpha_i^-(R)},
\qquad
\rho_i(R)=\frac{\sqrt{\alpha_i^+(R)\alpha_i^-(R)}}{H_i}.
```

Summing the reveal identity over `R=S union {i}` gives

```math
\boxed{
N_{m,i}
=\sqrt{Z_{m+1,i}^+Z_{m+1,i}^-}\,
H_i\,\mathbb E_{\rho_i}\operatorname{sech}B_{R\setminus i,i}.
}
```

Therefore the exact fixed-density analogue of the one-deletion shared-floor
formula is

```math
\boxed{
\sqrt{r_i(d^+)r_i(d^-)}
=\Gamma_{m,i}(d^+,d^-)
H_i(d^+,d^-)
\mathbb E_{\rho_i}\operatorname{sech}B,
}
```

where

```math
\Gamma_{m,i}
=\sqrt{\frac{Z_{m+1,i}^+Z_{m+1,i}^-}{Z_m(d^+)Z_m(d^-)}}
=\frac{n-m}{n}
\sqrt{
\frac{U_{m+1\mid i}(d^+)U_{m+1\mid i}(d^-)}
     {U_m(d^+)U_m(d^-)}
}.
```

Here `U_(m+1|i)` is the uniform average over `(m+1)`-sets containing `i`;
the binomial factor is
`binom(n-1,m)/binom(n,m)=(n-m)/n`.

Equivalently, define

```math
C_{\rm lev}=-\log\Gamma_{m,i},
\quad C_{\rm aff}=-\log H_i,
\quad C_{\rm rev}=-\log\mathbb E_{\rho_i}e^{-D(U\Vert Q_{S,i})}.
```

Then

```math
\boxed{
-\frac12\log\{r_i(d^+)r_i(d^-)\}
=C_{\rm lev}+C_{\rm aff}+C_{\rm rev}.
}
```

The affinity and reverse-KL costs are nonnegative, and Jensen gives
`C_rev <= E_rho D(U||Q)`.  The cross-level cost has no sign.

If the endpoint response reverses base odds `omega_i`, (10.981) gives
`min(r_i(d^+),r_i(d^-))<=e^{-|omega_i|}`.  Hence every crossing obeys the
exact trichotomy

```math
\boxed{
C_{\rm lev}+C_{\rm aff}+C_{\rm rev}
\ge\frac12|\omega_i|.
}
```

This is a genuine fixed-density anti-resonance reduction, but not a bound:
large odds may be paid by a reveal KL, collapse of state-edge selector
affinity, or a cross-level incidence imbalance.

The `H_i` here is **not** the adjacent-selector Hellinger quantity required
elsewhere in the harmonic route.  It compares posterior laws on
`(m+1)`-selectors containing `i` at two neighboring **states**.  The Johnson
Hellinger input compares component state laws for two neighboring
**selectors**.  They are transposed axes of the joint law, and no direct
inequality between them is currently proved.  Adjacent-selector control only
suggests a possible Bayes/transport argument; it cannot presently be cited as
control of `H_i`.

## 4. Exact-minimizer wall to pointwise fixed-density anti-resonance

The displayed exact order-nine minimizer `A_9` has `Q(A_9)=q_9=24`.  Take
`m=4`, orientation `sigma=+1`, and the two projective spins

```text
x = (1,-1,-1, 1,-1,-1, 1,-1,-1),
y = (1,-1,-1, 1, 1,-1, 1,-1,-1).
```

They differ only at vertex `4` and both have parent energy `-8`.  Therefore
the base edge odds are exactly

```math
\omega_4=0
```

for every `beta`.

For a state `d`, define the tropical completion exponent

```math
a(d)=\min_{|S|=4}\max_{z:z[S]=d[S]}
\{\langle A_9,z\rangle-c_S(d[S])\}.
```

Exact enumeration gives

```math
a(x)=16,
\qquad a(y)=8.
```

If `N_S` is the multiplicity of the maximum external energy, the sums of
leading reciprocal multiplicities over selectors attaining `a(d)` are

```math
\sum_{S:a_S(x)=16}\frac1{N_S}=\frac{1318}{105},
\qquad
\sum_{S:a_S(y)=8}\frac1{N_S}=\frac18.
```

The common factors `2^5/binom(9,4)` cancel, so the actual matched score jump
has the exact low-temperature asymptotic

```math
\boxed{
\chi_4(x,y)
=8\beta+\log\frac{105}{10544}+o(1).
}
```

At `beta=8`, the exact finite sum gives `59.39064809341531`, while the
displayed leading expression is `59.39064809341524`.  Thus an actual
quadratic completion, common Gibbs base, and exact minimizer have a
zero-base-odds fixed-density edge with unbounded score jump.  No uniform-in-
temperature pointwise analogue of (10.1129) is possible.

The affinity factorization explains the failure.  On this edge, at increasing
`beta`, the three costs have slopes

```text
C_rev  ~  4 beta,
C_aff  ~  4 beta,
C_lev  ~ -4 beta,
total  ~  4 beta.
```

For example, at `beta=2`,

```text
chi = 11.3924395965,
r_4(x)=0.5509719858,
r_4(y)=0.00000621521,
(C_rev,C_aff,C_lev,total)
=(6.60333354,6.76796723,-7.07900966,6.29229111).
```

Thus reveal composition alone permits both exponentially collapsing
state-edge posterior affinity and exponentially favorable cross-level
incidence.

## 5. Why this does not falsify signed migration

The bad edge is Gibbs-rare.  Exact enumeration of `a(d)` over all 512
oriented states shows that every energy-`24` ground has `a(d)=12` and that

```math
\max_d\{\langle A_9,d\rangle-s a(d)\}=24-12s,
\qquad 0\le s\le1.
```

The two bad endpoints have tropical lines

```math
-8-16s,
\qquad -8-8s.
```

Their gap below the dominant line is therefore at least `28`, uniformly in
time.  Finite-sum Laplace bounds give

```math
\mu_s\{x,y\}=O(e^{-28\beta})
```

uniformly for `0<=s<=1`.  Since `|chi|=O(beta)`, the edge curvature and its
entire weighted-energy or signed-migration contribution are

```math
O(\beta^2e^{-28\beta}).
```

Hence this is a sharp pointwise/composition mechanism wall, not an
asymptotic project-temperature counterexample and not a falsifier of the
integrated signed-migration target (10.1083).  A viable theorem must weight
the reveal/affinity trichotomy by actual context mass and prove that costly
edges are correspondingly rare.  Bounding the three factors pointwise is
impossible.
