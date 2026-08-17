# Spectral anti-pins and the finite-port Gram algebra

**Status.** Rigorous task-local theorem and scalable counterexample.  This
abstracts BCX from the Walsh formula.  The one-port compiler needs only two
resolvent quadratic features; for a symmetric two-eigenvalue matrix these
collapse to one Rayleigh coordinate.  A fixed number of ports has an
`O(l^2)` Gram certificate, but that certificate is not a reusable
congruence: composing independently summarized port families creates new
cross-Gram information at leading scale.

## 1. The resolvent criterion is the weakest input used by the spherical proof

Let `H in {+-1}^{n times n}` be symmetric and satisfy

```math
H\mathbf1=r\mathbf1,
\qquad ||H||_(2->2)\le r,
\qquad \operatorname{tr}H=0,                         \tag{SA.1}
```

where `r>0`.  Put

```math
A=H-\operatorname{diag}(H),
\qquad A_s=D_sAD_s.                                  \tag{SA.2}
```

Then `A_s` is a hollow complete signing and, on Boolean vectors,

```math
H_(A_s)(x)={1\over2}(s\odot x)^TH(s\odot x),
\qquad Q(A_s)={rn\over2}.                            \tag{SA.3}
```

Fix an integer `m` with `2m>r`.  For a Boolean pair `s,t`, write
`w=s odot t` and define the two normalized resolvent features

```math
\Psi_\sigma(w)
={m\over2n}w^T(2mI-\sigma H)^{-1}w,
\qquad \sigma\in\{+1,-1\}.                          \tag{SA.4}
```

### Theorem SA.1 (general one-port anti-pin criterion)

Let `mathcal S subset {+-1}^n` and suppose there is `delta>0` such that
for all distinct `s,t in mathcal S`,

```math
\max_{\sigma=+-1}\Psi_\sigma(s\odot t)
\le {r\over2m}-\delta.                               \tag{SA.5}
```

For query `t`, append `m` spins, use the query-owned bridge
`t 1_m^T`, and put any fixed public complete signing `C` on the new shore.
Let

```math
F_s(t)=Q\left(
 H_(A_s)(x)+(t\mathbin\cdot x)(\mathbf1_m\mathbin\cdot y)
 +H_C(y)\right).                                    \tag{SA.6}
```

If `E_C=Q(C)` and `H_C(\mathbf1)=E_C`, then

```math
F_s(s)={rn\over2}+mn+E_C                             \tag{SA.7}
```

provided `C` has a Boolean maximizer at `mathbf1` (in particular for the
positive clique), while for `s!=t`,

```math
F_s(t)\le {rn\over2}+mn-\delta mn+E_C.               \tag{SA.8}
```

Consequently the projective metric of the cap-response profiles has gap at
least `delta mn`.  Moreover, if

```math
d_0(s,t)=Q(A_s-A_t),
```

then

```math
d_C(s,t)\le d_0(s,t)\le rn,
\qquad
d_C(s,t)\ge {\delta m\over r}d_0(s,t).              \tag{SA.9}
```

Thus `m=Theta(r)=Theta(sqrt(n))` and fixed `delta` give a constant-gain,
linear-order, bounded-cap exact-sign metric compiler whenever the pairwise
resolvent condition (SA.5) holds.

#### Proof

After switching `x` by `s`, omission of the auxiliary internal energy gives
the two absolute channels

```math
R_\sigma(w)=\max_{u\in\{\pm1\}^n}
 \left\{{\sigma\over2}u^THu+m w^Tu\right\}.         \tag{SA.10}
```

Because `||u||_2^2=n`, putting `K_sigma=2mI-sigma H` and completing the
square gives

```math
R_\sigma(w)
\le mn+{m^2\over2}w^TK_\sigma^{-1}w
=mn(1+\Psi_\sigma(w)).                               \tag{SA.11}
```

Both `K_sigma` are positive definite by (SA.1).  Condition (SA.5) turns
(SA.11) into the cap without `C` bounded by

```math
mn+{rn\over2}-\delta mn.
```

Adding `C` changes a cap by at most `E_C`, proving (SA.8).

On the diagonal, `x=s` and `y=mathbf1` attain
`rn/2+mn+E_C`.  The separate child and bridge bounds give the matching
upper bound before `C`, and `E_C` gives it afterward.  This proves (SA.7).
The two diagonal queries for a pair give opposite response gaps, hence
`d_C>=delta mn`.  Changing only the old child block gives
`d_C<=d_0`, while (SA.3) gives `d_0<=rn`.  These inequalities imply
(SA.9). `square`

The trace-zero condition is only an exact normalization convenience.
For a sign matrix, deleting a nonzero diagonal changes every Boolean energy
by the scalar `tr(H)/2`, whose magnitude is at most `n/2=o(n^(3/2))` in the
spectral regime above.  The same argument then has an additive `O(n)` loss.

The criterion contains strictly less information than a Boolean response
landscape: it stores two real numbers per child--query pair.  It also makes
no assertion about the old or auxiliary optimizer off the diagonal.  It is
the weakest condition used by this resolvent certificate, not a claimed
necessary condition for the exact Boolean response.

## 2. Two eigenvalues collapse the resolvents to one Rayleigh bit

### Corollary SA.2 (the abstract Hadamard compiler)

Suppose, in addition,

```math
H^2=r^2I,
\qquad m=r\in\mathbb N.                              \tag{SA.12}
```

For a Boolean `w`, put

```math
\rho(w)={w^THw\over rn}.                             \tag{SA.13}
```

Then

```math
\Psi_+(w)={2+\rho(w)\over6},
\qquad
\Psi_-(w)={2-\rho(w)\over6}.                        \tag{SA.14}
```

Therefore every code satisfying

```math
|\rho(s\odot t)|\le\theta<1\qquad(s\ne t)          \tag{SA.15}
```

has the compiler of Theorem SA.1 with

```math
\delta={1-\theta\over6}.                             \tag{SA.16}
```

#### Proof

The involution identity gives

```math
(2rI\mp H)^{-1}={2rI\pm H\over3r^2}.
```

Substitution in (SA.4) proves (SA.14), and

```math
\max(\Psi_+,\Psi_-)={2+|\rho|\over6}
\le{1\over2}-{1-\theta\over6}.
```

Apply Theorem SA.1. `square`

Thus Walsh coordinates are one way to produce `(H,mathcal S)`, not part of
the compiler mechanism.  Any trace-zero regular symmetric Hadamard matrix,
together with a two-sided Rayleigh code, works.  Hanson--Wright produces an
`exp(Omega(n))` code for every such `H` because
`||H||_F=n` and `||H||=sqrt(n)`.

## 3. A fixed number of ports has a quadratic Gram certificate

Now let a query have `l` repeated ports `w_1,...,w_l`, each with `m` copies.
Ignore the public auxiliary completion for the moment.  Endpoint
optimization in each port and the absolute child channel reduce the cap to
the maximum, over `sigma in {+-1}` and `epsilon in {+-1}^l`, of

```math
\max_{u\in\{\pm1\}^n}
\left\{{\sigma\over2}u^THu
 +m\left(\sum_a\epsilon_aw_a\right)^Tu\right\}.     \tag{SA.17}
```

Assume `H^2=r^2I`, retain `2m>r`, but do not require `m=r`, and define the
two port Gram matrices

```math
G_(ab)={w_a^Tw_b\over n},
\qquad
R_(ab)={w_a^THw_b\over rn}.                          \tag{SA.18}
```

### Theorem SA.3 (finite-port spherical carrier)

For `epsilon in {+-1}^l`, put

```math
g_epsilon=\epsilon^TG\epsilon,
\qquad h_epsilon=\epsilon^TR\epsilon.
```

The spherical relaxation of (SA.17) is bounded by, and in the usual
trust-region dual formulation equals,

```math
rn\inf_{\alpha>1/2}
\left\{
\alpha+{(m/r)^2(2\alpha g_\epsilon
                         +\sigma h_\epsilon)
                  \over2(4\alpha^2-1)}
\right\}.                                           \tag{SA.19}
```

Consequently the complete one-layer spherical certificate for `l` ports is
a function only of the `l(l+1)` entries of `(G,R)` and the finite
`2^(l+1)` channel list.  At fixed `l` it has `O(l^2 log n)` exact
description bits.  A public exact-sign completion on the `lm` auxiliary
vertices perturbs the cap by at most

```math
{lm\choose2}=O_l(n)                                  \tag{SA.20}
```

when `m=Theta(sqrt(n))`.

#### Proof

Let `v=sum_a epsilon_a w_a`.  For any `a_0>r/2`, completing the square with
`2a_0I-sigma H` gives

```math
{\sigma\over2}u^THu+mv^Tu
\le a_0n+{m^2\over2}v^T(2a_0I-\sigma H)^{-1}v.
                                                                  \tag{SA.21}
```

The trust-region Lagrange dual is exact for the Euclidean sphere; taking the
infimum over `a_0` therefore gives its value (with the boundary understood
by a limit in the hard case).  Write `a_0=alpha r` and use

```math
(2\alpha rI-\sigma H)^{-1}
={2\alpha rI+\sigma H\over r^2(4\alpha^2-1)}.
```

Since `||v||^2=ng_epsilon` and
`v^THv=rn h_epsilon`, (SA.19) follows.  The Gram count is immediate.
Finally uniform cap Lipschitzness and the trivial cap of a complete signing
on `lm` vertices give (SA.20). `square`

Theorem SA.3 is useful as a constant-scale **certificate**.  It is not an
`o(n^(3/2))` reconstruction theorem for the Boolean cap: the spherical
relaxation can retain a fixed leading integrality gap.

## 4. Why the fixed-port carrier is not a reusable congruence

The obstruction already occurs when two individually identical one-port
states are combined.

Let `H_16` be the regularized Walsh matrix (BCX.0) of order `16`.  It has a
balanced Boolean `+4` eigenvector

```text
v_0 = (-,-,-,+;  -,-,+,-;  +,-,+,+;  -,+,+,+).
```

For `j>=1`, set

```math
H=H_16^{\otimes j},
\qquad n=16^j,
\qquad r=4^j=\sqrt n,
\qquad v=v_0\otimes\mathbf1_{16^{j-1}}.              \tag{SA.22}
```

Then

```math
H\mathbf1=r\mathbf1,
\quad Hv=rv,
\quad \mathbf1^Tv=0.                                 \tag{SA.23}
```

### Theorem SA.4 (composition creates a leading cross-Gram variable)

The one-port objects `w=mathbf1` and `w=v` have identical self state

```math
G=(1),\qquad R=(1),                                  \tag{SA.24}
```

and identical one-port absolute response `3rn/2`.  Nevertheless, combining
two copies gives

```math
\mathcal Q(\mathbf1,\mathbf1)={5\over2}rn,           \tag{SA.25}
```

whereas

```math
\mathcal Q(\mathbf1,v)
\le\left({1\over2}+\sqrt2\right)rn.                 \tag{SA.26}
```

Here

```math
\mathcal Q(w_1,w_2)
=\max_u\left\{|H_A(u)|
 +r|w_1^Tu|+r|w_2^Tu|\right\}                       \tag{SA.27}
```

is the two-port cap before an `O(n)` public completion.  The cap gap obeys

```math
\mathcal Q(\mathbf1,\mathbf1)-\mathcal Q(\mathbf1,v)
\ge(2-\sqrt2)rn=Theta(n^{3/2}).                     \tag{SA.28}
```

#### Proof

For either Boolean `+r` eigenvector `w`, the one-port response is at most
`rn/2+rn` by the spectral and Cauchy--Schwarz bounds, and `u=w` attains it.
This proves (SA.24) and the one-port assertion.

For two copies of `mathbf1`, `u=mathbf1` attains (SA.25), again matching
the separate upper bounds.  For `(mathbf1,v)`, choose the signs of the two
absolute fields and put `z=+-mathbf1+-v`.  Orthogonality in (SA.23) gives
`||z||_2=sqrt(2n)` for every choice.  Therefore

```math
|H_A(u)|\le rn/2,
\qquad r|z^Tu|\le r||z||_2||u||_2=\sqrt2rn,
```

which proves (SA.26). `square`

The full two-port Gram carrier distinguishes the examples: the cross entries
of `(G,R)` are both one in (SA.25) and both zero in (SA.26).  But those cross
entries are absent from the two separate one-port states.  More generally,
combining `p` ports creates all cross matrices

```math
{w_a^Tz_b\over n},
\qquad {w_a^THz_b\over rn}.                          \tag{SA.29}
```

They are not functions of the two internal Gram pairs.  Hence:

- for fixed `l`, `(G,R)` is a polynomial-size one-layer certificate;
- under repeated composition, this Gram presentation must adjoin new
  cross-Gram data and has an `O(p^2)` table after `p` ports; SA.4 proves that
  the product of the separate states is insufficient, but does not claim a
  universal `Omega(p^2)` lower bound for every possible representation;
- a bounded reusable quotient needs an additional synchronization law which
  makes the cross entries functions of the retained state.

This is the precise boundary between BCX's one-layer metric embedding and a
compositional congruence.  The anti-pin theorem is positive without supplying
the missing synchronization.

The finite identities and the scalable tensor witness are checked in
[`../experiments/verify_spectral_antipin_feature_algebra.py`](../experiments/verify_spectral_antipin_feature_algebra.py).
