# Latent child entropy does not orient the rank-one overlap floor

Status: **rigorous exact identities and optimizer-specific finite no-go**.
This audit starts from the actual zero-bridge law of two contracted-
temperature children.  It identifies the exact posterior-information budget
behind the rank-one support obstruction and separates sector retuning, child-
factor retuning, and induced cross-child dependence.  The separation is
valid under every raw negative disorder tilt.

The identities do not decide whether the overlap floor is row-additive or
irreducible reverse row dependence.  In fact, the required relative-entropy
term is evaluated under the already retuned latent law, so controlling it
from scalar child Gibbs entropy would restate the retuning obligation.  The
two certified equal-entropy order-eight optimizing children provide a finite
actual-minimizer witness that entropy does not determine even the first
orientation-independent cross-row response coefficient.

## 1. Exact sector normalization

Fix child orders `m,n`, put `d=mn`, and use the common raw temperature `t`.
For `a in {+-1}` define

```math
 Z_A^a=2^{-m}\sum_x e^{taH_A(x)},
 \qquad
 \mu_{A,a}(x)={e^{taH_A(x)}\over 2^m Z_A^a},       \tag{LE.1}
```

and define `Z_D^a,mu_(D,a)` analogously.  Start with the two independent
augmented Gibbs laws in equations (2.1)--(2.2) of
[`artifacts/finite_temperature_reverse_kl_interface.md`](../../artifacts/finite_temperature_reverse_kl_interface.md).
Condition on the relative orientation

```math
 \epsilon=\tau_1\tau_2.
```

Writing `S=tau_1`, the resulting zero-bridge latent law is exactly

```math
 \boxed{
 \nu_\epsilon(s,x,y)
 =\pi_s^{(\epsilon)}\mu_{A,s}(x)\mu_{D,\epsilon s}(y),
 \qquad
 \pi_s^{(\epsilon)}
 ={Z_A^sZ_D^{\epsilon s}\over
   \sum_{r=\pm1}Z_A^rZ_D^{\epsilon r}}.}           \tag{LE.2}
```

This is the sector law in (EO.5) of
[`actual_child_entropy_overlap_orientation_no_go.md`](../drafts/actual_child_entropy_overlap_orientation_no_go.md).
The induced rank-one bridge word is

```math
 Q=Q(S,X,Y),
 \qquad Q_{ij}=SX_iY_j.                             \tag{LE.3}
```

The binary-channel density conditional on the latent state is

```math
 k_t(B\mid Q)
 ={e^{t\langle B,Q\rangle}\over(\cosh t)^d}        \tag{LE.4}
```

relative to the fair bridge law `U`.  Hence the conditional forward output
density is

```math
 p_\epsilon(B)
 =\mathbb E_{\nu_\epsilon}k_t(B\mid Q)
 ={\mathbb E_{\nu_\epsilon}e^{t\langle B,Q\rangle}
   \over(\cosh t)^d}.                               \tag{LE.5}
```

The exact joint output identity (2.4) of the cited reverse-KL interface
shows that, for a constant `c_epsilon` independent of `B`,

```math
 L_\epsilon(B)=c_\epsilon+\log p_\epsilon(B).       \tag{LE.6}
```

Consequently, for every raw disorder exponent `a`, in particular every
`a in [-lambda,0]`, the actual bridge tilt is

```math
 {dq_a\over dU}(B)
 ={p_\epsilon(B)^a\over\mathbb E_Up_\epsilon^a}.   \tag{LE.7}
```

No conference or other surrogate law enters (LE.1)--(LE.7).

## 2. Exact posterior-information budget

Suppress `epsilon` in the notation.  Let the ordinary **forward-channel**
Bayes posterior be

```math
 \nu_B(dz)
 ={e^{t\langle B,Q(z)\rangle}\nu(dz)
   \over\mathbb E_\nu e^{t\langle B,Q\rangle}},    \tag{LE.8}
```

and put

```math
 m(B)=\mathbb E_{\nu_B}Q.
```

For any bridge law `q` define

```math
 \eta(dB,dz)=q(dB)\nu_B(dz),
 \qquad
 \bar\nu=\eta^Z.                                  \tag{LE.9}
```

Thus the following theorem applies in particular to every actual `q_a` in
(LE.7).

**Theorem LE.1 (posterior budget and factor split).**  One has

```math
 \boxed{
 \mathbb E_qD(\nu_B\Vert\nu)
 =I_\eta(B;Z)+D(\bar\nu\Vert\nu)
 =I_\eta(B;Q)+D(\bar\nu\Vert\nu).}                \tag{LE.10}
```

Moreover,

```math
 \boxed{
 \begin{aligned}
 D(\bar\nu\Vert\nu)
 &=D(\bar\pi\Vert\pi)\\
 &\quad+\sum_s\bar\pi_s\bigl\{
   D(\bar\nu_{X\mid s}\Vert\mu_{A,s})
  +D(\bar\nu_{Y\mid s}\Vert\mu_{D,\epsilon s})
  +I_{\bar\nu}(X;Y\mid S=s)\bigr\}.
 \end{aligned}}                                    \tag{LE.11}
```

Finally the posterior energy satisfies

```math
 \boxed{
 t\,\mathbb E_q\langle B,m(B)\rangle
 =\mathbb E_q\log p(B)+d\log\cosh t
  +I_\eta(B;Q)+D(\bar\nu\Vert\nu).}               \tag{LE.12}
```

*Proof.*  The first equality in (LE.10) is the relative-entropy chain rule:

```math
 \mathbb E_qD(\nu_B\Vert\nu)
 =D(\eta\Vert q\otimes\nu)
 =I_\eta(B;Z)+D(\bar\nu\Vert\nu).                 \tag{LE.13}
```

The likelihood in (LE.8) depends on `z` only through `Q(z)`.  Therefore
both `nu_B(z|Q)` and `bar nu(z|Q)` equal `nu(z|Q)`.  Conditional relative
entropy vanishes on the fibres of `Q`, and `B-Q-Z` is a Markov chain under
`eta`.  This proves the second equality in (LE.10).

Apply the KL chain rule first to `S`, then use the product reference
`mu_(A,s) tensor mu_(D,epsilon s)` in each sector.  The elementary identity

```math
 D(r_{XY}\Vert p_X\otimes p_Y)
 =D(r_X\Vert p_X)+D(r_Y\Vert p_Y)+I_r(X;Y)         \tag{LE.14}
```

gives (LE.11).  Pointwise in `B`, (LE.8) gives

```math
 D(\nu_B\Vert\nu)
 =t\langle B,m(B)\rangle
  -\log\mathbb E_\nu e^{t\langle B,Q\rangle}
 =t\langle B,m(B)\rangle-\log p(B)-d\log\cosh t. \tag{LE.15}
```

Average (LE.15) and use (LE.10) to obtain (LE.12). `square`

Equation (LE.11) is the strongest exact entropy allocation furnished by
the latent representation.  Its first term is sector retuning, its next two
terms are coherent retuning of the two actual child Gibbs factors, and its
last term is induced cross-child dependence.

## 3. What the child entropy scalar knows

For fixed `(s,Q)`, the equation `Q=sxy^T` has exactly the two preimages
`(x,y)` and `(-x,-y)`.  They have equal conditional Gibbs weight because
both child Hamiltonians are quadratic.  Therefore

```math
 H(X,Y\mid S,Q)=\log2.                              \tag{LE.16}
```

Since `X,Y` are independent conditional on `S`, this gives the exact latent
entropy identity

```math
 \boxed{
 H(Q)
 =\sum_s\pi_s\{H(\mu_{A,s})+H(\mu_{D,\epsilon s})\}
  -\log2+I(S;Q),
 \qquad 0\le I(S;Q)\le\log2.}                     \tag{LE.17}
```

At every positive finite temperature all spin configurations have positive
weight.  Thus, independently of the optimizing children,

```math
 \boxed{|\operatorname {supp}Q|=2^{m+n-1}.}        \tag{LE.18}
```

This is the exact support count used in Theorem 37.52 and in
[`actual_child_rank_one_support_overlap_obstruction.md`](actual_child_rank_one_support_overlap_obstruction.md).
Hence the separate conditional child entropies determine the Shannon
entropy of `Q` only up to the single sector bit, while support entropy is
the same for every actual child pair.

Neither fact controls the retuning term in (LE.10).  Indeed,

```math
 D(\bar\nu\Vert\nu)
 =-H(\bar\nu)-\mathbb E_{\bar\nu}\log\nu,          \tag{LE.19}
```

whereas child Gibbs entropy contains `-E_nu log nu`, under the **prior**.
The missing cross-expectation in (LE.19) is exactly the response of the
child weights to the negative-tilt posterior.  Proving that it is `o(N)`
rules out the factor-retuning terms in (LE.11); proving that it is linear
selects a retuning/dependence branch.  Thus such a bound would be a new
optimizer-specific retuning theorem, not a consequence of the scalar
entropy identity.

There is a second exact mismatch.  Under `q_a`, Bayes' formula gives

```math
 \eta_a(dB\mid Q)
 \ \propto\ U(dB)k_t(B\mid Q)p(B)^{a-1}.           \tag{LE.20}
```

For the forward law `a=1`, the last factor disappears and the bridge bits
are conditionally independent.  On the negative path `a<=0`, the global
factor `p^(a-1)` remains.  Therefore the mutual-information term
`I_eta(B;Q)` in (LE.10) is not a common-latent product representation of
the negative escort and is not, by data processing alone, a lower bound on
the reverse row-product projection.  Extracting such a bound would require
the missing reverse-tensorization or no-product-background theorem.

## 4. Exact optimizing-child entropy collision

The lack of geometric information is witnessed inside the actual optimizer
class.  The two order-eight signings `A_0,A_1` displayed in (FC.22)--(FC.24)
of
[`actual_child_flip_averaging_ceiling.md`](../drafts/actual_child_flip_averaging_ceiling.md)
have the same exact projective absolute-energy histogram

```text
|H|       0   2   4   6   8  10
count    12  32  32  24  20   8.
```

The certified classification
[`computations/results/m8_minimizer_orbits.json`](../../computations/results/m8_minimizer_orbits.json)
and the comparison (FC.22) show that both are exact thermal-pressure
minimizers for every `t>=3`.  The common histogram gives identical pressure
and Gibbs entropy profiles at every temperature.  The finite calculations
are reproduced by
[`actual_child_radial_ceiling_witness.py`](../experiments/actual_child_radial_ceiling_witness.py).

Pair either class with the same order-two minimizing child.  The exact
oriented-overlap calculation in Theorem EO.4 of
[`actual_child_entropy_overlap_orientation_no_go.md`](../drafts/actual_child_entropy_overlap_orientation_no_go.md)
is independent of orientation and gives

```math
 K(A_0,D;\infty)=20,
 \qquad K(A_1,D;\infty)=12.                        \tag{LE.21}
```

The strict separation persists at every sufficiently large finite `t`.
By Theorem EO.2, the corresponding integrated cross-row response tangents
obey

```math
 \min_\epsilon\mathsf T_u(A_0,D,\epsilon)
 -\min_\epsilon\mathsf T_u(A_1,D,\epsilon)
 ={\lambda^2u^4\over2}\{K(A_0,D;t)-K(A_1,D;t)\}
  +O(u^6)>0                                        \tag{LE.22}
```

for all sufficiently small positive channel amplitudes `u`.  Thus even the
complete scalar child pressure/entropy profiles do not determine the first
orientation-independent cross-row response coefficient of exact optimizing
children.

There is a complementary physical-amplitude fixed-sector witness.  The two
order-three minimizing pairs in Theorem EO.3 have identical separate-child
pressure and entropy data but satisfy

```math
 \mathsf T_+(q)-\mathsf T_-(q)
 =54q^{10}+O(q^{12})>0                              \tag{LE.23}
```

at the physical channel amplitude.  Reversing one child exchanges the two
orientation labels, so (LE.23) is deliberately not claimed as an
orientation-minimized obstruction; (LE.21)--(LE.22) supply the
orientation-invariant tangent obstruction.

## 5. No-go scope

The exact posterior budget (LE.10)--(LE.12) is useful bookkeeping, but it
does not narrow the current directional lemma:

1. child Gibbs entropy controls the prior entropy, not the retuned
   cross-entropy in (LE.19);
2. the factor KL terms in (LE.11) are the coherent-retuning obligation one
   hoped entropy would decide;
3. the remaining latent mutual information does not productize the negative
   escort because of (LE.20);
4. equal-entropy actual minimizers have different cross-row response even
   after the orientation minimum.

Therefore a scalar child entropy/effective-support criterion cannot, by the
present exact identities, distinguish the support-induced overlap floor
into row-additive retuning versus irreducible reverse row dependence.  A
successful theorem must add nonradial posterior geometry--for example a
genuine reverse-tensorization, synchronization, or negative-tail theorem.
This audit proves no asymptotic lower bound on the reverse row-product KL
and no Level-6 recurrence.
