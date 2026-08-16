# Adversarial verification of the posterior-width information inequality

**Scope.** This report audits only
`phase2_general_information_inequality.md`.  I reconstructed the argument
without using its proof order, checked the constants against finite exhaustive
tables, and compared it with the previous hard Walsh-decoding bound.

## Verdict

**Overall: ACCEPT after the explicit minor corrections listed below.**

The binary entropy curve, nonlinear posterior-width identity, weighted and
homogeneous forms of Theorem 5.1, and the sharpness constant are correct.  The
Ising and nearest-code applications are correct.  The optional Max-Cut
application is correct only after its omitted probabilistic, distortion, and
exposure hypotheses are restored.  No tested value of `kappa` is too large.

| Claim | Verdict | Reason |
|---|---|---|
| Lemma 3.1 | **ACCEPT** | The identity and both Jensen steps are exact; `g` is increasing and concave. |
| Lemma 4.1 and (4.4) | **ACCEPT** | The two-copy variance calculation gives exactly `1/4`, and a two-point posterior attains it. |
| Weighted Theorem 5.1 | **ACCEPT** | Conditional-mean projection plus the weighted posterior-width inequality gives exactly the displayed feasible region. |
| Homogeneous bound (5.2) | **ACCEPT** | Concavity makes equal posterior coordinate variances extremal; the cutoff at one is necessary. |
| Sharpness | **ACCEPT** | Independent binary symmetric channels attain equality.  In fact the weighted statement is sharp as well. |
| Ising application | **ACCEPT**, clarify | The centering constant cannot cancel a Walsh component.  One may project it away, and switching gives exact `kappa=4a^2`. |
| Rooted-code application | **ACCEPT**, clarify | Exact `kappa=1/K`; qualify the informal “below `1/4`” statement asymptotically. |
| Max-Cut application | **CORRECT** | Add uniform prior, decoder/distortion definition, and `M>(n-1)/2`; without a uniform prior the information conclusion is false. |

## 1. Independent reconstruction of Lemma 3.1

For a sign `X`, write `p=P(X=1)` and `w=2p-1`.  Then

```math
v=\operatorname{Var}(X)=1-w^2=4p(1-p).
```

The smaller of `p,1-p` is

```math
p_*(v)={1-\sqrt{1-v}\over2},
```

so symmetry of binary entropy gives

```math
H(X)=h_2(p)=h_2(p_*(v))=g(v).
```

With `t=sqrt(1-v)`, direct differentiation in bits gives

```math
g'(v)={\operatorname{arctanh}t\over2t\log2}>0
```

and

```math
g''(v)=
-{t/(1-t^2)-\operatorname{arctanh}t\over4t^3\log2}<0.
```

The numerator in the last expression starts at zero and has derivative
`2t^2/(1-t^2)^2`.  Endpoint continuity handles `v=0,1`.  Thus the claimed
monotonicity and concavity are correct.

For every posterior value `z`, entropy subadditivity gives

```math
H(A\mid Z=z)\le\sum_e g(v_e(z)).
```

Jensen first over `z`, then over coordinates, gives exactly (3.3).  No
posterior independence is assumed.  Uniformity and independence are used
only in the prior identity `H(A)=N` later.

## 2. Independent reconstruction of Lemma 4.1 and the infimum identity

Let `A,A'` be independent posterior draws.  For every Hilbert-valued response,

```math
\operatorname{Var}_\pi(R)
={1\over2}\mathbb E\|R_A-R_{A'}\|^2.
```

If `gamma in Gamma(R)`, pairwise separation implies

```math
\operatorname{Var}_\pi(R)
\ge {1\over2}\sum_e\gamma_e P(A_e\ne A'_e).
```

For a sign coordinate of posterior mean `w_e`, the disagreement probability
is `(1-w_e^2)/2=Var(A_e)/2`.  Hence

```math
\operatorname{Var}_\pi(R)
\ge {1\over4}\sum_e\gamma_e\operatorname{Var}_\pi(A_e).
```

This verifies both the factor and direction in (4.3).

For the exact infimum (4.4), choose a minimizing pair `a,b` in the finite
cube and put probability `1/2` on each.  Its response variance is
`||R_a-R_b||^2/4`, while every changed sign coordinate has variance one, so
`V_pi(A)=d_H(a,b)`.  The ratio is `kappa/4`.  The preceding lower bound proves
that no posterior has a smaller ratio.  This also covers `kappa=0`: a
colliding pair attains zero.  Thus (4.4) has no missing convexity or injectivity
hypothesis.

## 3. Theorem 5.1, including weights and sharpness

Conditionally on `Z=z`, the posterior mean response is the Hilbert least
squares estimator.  Therefore an arbitrary decoder satisfies

```math
\Delta\ge
\mathbb E_Z\operatorname{Var}(R_A\mid Z)
\ge {1\over4}\sum_e\gamma_e\bar v_e,
```

where `bar v_e=E Var(A_e|Z)`.  Lemma 3.1 independently gives

```math
I(A;Z)\ge N-\sum_e g(\bar v_e).
```

The actual posterior-variance vector is consequently feasible in (5.1), so
maximizing over the larger displayed feasible set has the correct inequality
direction.  For constant weight `kappa`, its mean coordinate variance is at
most `min{4Delta/(kappa N),1}`.  Increasing concavity of `g` gives (5.2).
Random decoder coins may be appended to the transcript and do not alter the
argument.

The claimed homogeneous sharpness is exact.  For `R_a=s a`,
`kappa=4s^2`.  Passing every bit through a BSC of crossover `p` and decoding
by the posterior mean gives

```math
I=N[1-h_2(p)],\qquad
\Delta=4s^2Np(1-p),\qquad
g(4p(1-p))=h_2(p).
```

There is also a useful omitted weighted sharpness statement.  For any
nonnegative weights, take

```math
R_a=\left({\sqrt{\gamma_e}\over2}a_e\right)_{e=1}^N.
```

Then pairwise squared distance is exactly the sum of the changed `gamma_e`.
Independent BSCs with posterior variances `v_e` give

```math
\Delta={1\over4}\sum_e\gamma_ev_e,
\qquad
I=N-\sum_e g(v_e).
```

Thus the weighted entropy-allocation envelope is also the best possible one
under only the `Gamma(R)` hypothesis.  Adding this observation would make the
scope of the sharpness claim clearer, but it is not required for correctness.

## 4. Exposed charts and projection

Equations (6.2)--(6.6) are correct.  In an affine chart,
`(a-b)/2` ranges over every nonzero ternary vector, and

```math
\|R_a-R_b\|^2=4t^TGt,
\qquad d_H(a,b)=|\operatorname{supp}t|.
```

The restricted-eigenvalue expression therefore has the right factor four.
It evaluates all channels jointly and does not separately charge their
absolute values.

Projection needs a wording clarification, not a proof repair.  If `P` is a
contraction, apply it to both the true response and decoder.  The projected
distortion is at most the declared distortion, so a lower bound for the
projected task is a valid lower bound for the original task.  This is valid
even for a parameter-dependent component, provided the text explicitly says
that it is proving a lower bound for a contracted task; it would not be valid
to claim that the contracted response is sufficient for the original query
experiment.  The last two sentences of Section 6 currently blur these two
uses.

## 5. Application I: shifted dense Ising

### Exposure

If `x` differs from `u` on a set of size `k`, the field loses `2Mk`, while
the interaction can gain at most `2ak(n-k)`.  Thus `M>a(n-1)` exposes `u`.
The rank-one query loses `2Lk(n-k)` and the latent energy gains at most
`2ak(n-k)`, so `L>a` exposes exactly `{u,-u}`.  Both yield
`R_A(u)=q_A(u)-c_A` after subtracting the declared query constant.

### The nuisance shift cannot hurt

The parameter-dependent centering term is a degree-zero Walsh coefficient,
orthogonal to every degree-two character.  Hence

```math
\|R_A-R_B\|_2^2
=4a^2d_H(A,B)+(c_A-c_B)^2.
```

There is no cross term and therefore no possible cancellation.  Equivalently,
projecting onto the degree-two Walsh subspace removes the nuisance constant
and leaves squared distance exactly `4a^2 d_H(A,B)`.  This contraction is
legitimate for a lower bound even though the full query response does observe
`c_A`.

In fact the full response map has

```math
\boxed{\kappa=4a^2}\qquad(n\ge2),
```

not merely `kappa>=4a^2`.  For the reverse inequality, take any signing and
switch one vertex.  Exactly `n-1` edge signs change, the two quadratic forms
are related by a spin relabeling and hence have the same `c_A`, and their
response squared distance is `4a^2(n-1)`.  This is a clean recommended
addition to Section 7.

Substitution in (5.2) gives (7.5).  For
`D=Delta/(a^2N)`, the old result is `N[1-h_2(D)]` on `D<=1/2`, while the new
one is

```math
N\left[1-h_2\left({1-\sqrt{1-D}\over2}\right)\right].
```

For every `0<D<=1/2`, the new binary argument is strictly smaller than `D`,
so the information lower bound is strictly larger.  It remains positive for
every `D<1`, as claimed.

## 6. Application II: rooted nearest-code response

The exposure proof is correct: `-d(x,C)` can gain at most one per changed
spin, while `M<u,x>` loses `2M` per changed spin, so `M>1/2` uniquely exposes
`u`.  Thus the complete response is `Mm-d(u,C)`.

If two anchored codes differ in `r` membership coordinates, then at each
root in the symmetric difference one profile is zero and the other is an
integer at least one.  Under the uniform measure on `K=2^m` roots,

```math
\|R_a-R_b\|_2^2\ge r/K.
```

Taking the full cube and deleting one nonanchor word gives equality, so
`kappa=1/K` exactly.  With `N=K-1`, (8.6) follows with no missing factor.

The only needed correction is to the informal threshold wording.  A fixed
`Delta<1/4` yields a positive `Theta(K)` information lower bound **for all
sufficiently large `m`**.  At finite `m`, the exact condition from (8.6) is

```math
\Delta<{K-1\over4K}.
```

For example, merely saying `Delta<1/4` at `K=4` would be too weak.  The
claims for `Delta=o(1)` and uniform error below `1/2` are correct.

## 7. Application III: counterfactual Max-Cut

The algebraic modulus is correct:

```math
C_B(u)-C_D(u)
={1\over2}\sum_e(B_e-D_e)
-{1\over2}\sum_e(B_e-D_e)\chi_e(u),
```

and orthogonality gives (9.1).  Thus `kappa>=1/4`.  It is exact for `n>=3`:
take two one-edge graphs on distinct edges.  Their Hamming distance is two,
their constant terms cancel, and their squared response distance is `1/2`.
For `n=2`, the exact constant is `1/2`, so the stated lower certificate
`1/4` remains valid but is not sharp.

As currently written, however, (9.2) lacks hypotheses and is false for, say,
a deterministic `B`, when `I(B;Z)=0`.  Replace the opening of Section 9 by an
explicit statement of the following form:

> Let `B` be uniform on `{0,1}^N`.  Let `M>(n-1)/2`, so the vertex-prize
> query exposes `u`, and let a decoder have
> `Delta=E_{B,Z,U}[Rhat_Z(U)-R_B(U)]^2`.  Then (9.2) holds.

The exposure threshold follows because flipping `k` vertices can gain at
most `k(n-k)` cut edges and loses `2Mk` in the vertex prize.  With those
hypotheses restored, (9.2) follows exactly from `4Delta/(kappa N)` with the
valid certificate `kappa=1/4`.

## 8. Exhaustive falsification results

I independently enumerated the finite tables without writing an artifact.
All query averages used the uniform measure specified in the draft.

- For all anchored codes through `m=3`, the minimum ratios were respectively
  `1/2`, `1/4`, and `1/8`, exactly `1/K`.
- For every dense Ising signing pair through `n=5`, the minimum
  `||R_A-R_B||_2^2/(a^2 d_H(A,B))` was exactly `4`.
- For all Max-Cut graph pairs through `n=4`, the minimum was `1/2` at `n=2`
  and `1/4` at `n=3,4`.

These checks search the exact pair falsifier (10.1), rather than sampling
posteriors.  The analytic proofs above cover all orders.

## 9. Is the theorem genuinely stronger, or only packing?

It is genuinely stronger than the repository's old hard-decoding bound in a
precise quantitative sense.  Hard Walsh sign decoding converts squared error
to an edge-error probability and yields `h_2(D)` only while `D<=1/2`.  The
posterior-width proof retains the posterior magnitude and yields the sharp
binary squared-error curve `g(D)`, which is strictly smaller on the old range
and nontrivial until `D=1`.  A minimum-separation packing argument alone has
no such continuous sharp curve and typically becomes silent once exact
decoding fails.

The information-theoretic half is nevertheless classical binary
rate--distortion geometry, not a new entropy principle.  For the affine Ising
and Max-Cut frames, the theorem is a clean abstraction of that classical
converse after orthogonal projection.  Its genuinely generative content is
the nonlinear response-width reduction: the same statement applies to the
nearest-code distance profile, and a model can now be accepted or falsified
by proving one inverse-Hamming modulus.  The draft's director checkpoint
states this distinction accurately and should retain it.

## Exact edits required before promotion

1. In the one-paragraph summary, replace “average squared error below `1/4`”
   by “any fixed average squared error below `1/4`, for all sufficiently
   large code dimension,” or state the exact threshold `(K-1)/(4K)`.
2. In Section 6, distinguish explicit contraction for a lower bound from a
   claim of sufficiency.  A parameter-dependent constant may be projected
   away for the former, but not silently discarded for the latter.
3. In Section 7, explicitly state that the constant component is orthogonal
   and optionally strengthen `kappa>=4a^2` to the exact identity
   `kappa=4a^2` using a one-vertex switching witness.
4. In Section 9, add `B` uniform, `M>(n-1)/2`, and the exact definition of
   `Delta` before asserting (9.2).  Without the uniform-prior hypothesis the
   displayed mutual-information lower bound is false.
5. Optional but useful: record weighted sharpness using the coordinate-scaled
   cube construction above, and record exact Max-Cut `kappa=1/4` for `n>=3`.

No change is required to Lemmas 3.1 or 4.1, Theorem 5.1, equations
(7.4)--(7.5), or equations (8.4)--(8.6).

## 10. Follow-up: universal binary boundary kernels

**Verdict: ACCEPT with the normalization and scope stated explicitly.**

Let `A in {0,1}^{B times B}` be uniform, with `|B|=Q`, and regard its
`N=Q^2` entries as independent latent bits (equivalently replace each entry
by the sign `2A_ab-1` when invoking Theorem 5.1).  For the one-state-per-fibre
component in Theorem FG.4, an endpoint-pinning context with penalties of
magnitude greater than one exposes a prescribed fibre `(a,b)`.  After
removing the known context constant its response is exactly

```math
R_A(a,b)=A_{ab}.
```

Give the `Q^2` endpoint pairs the uniform probability measure and define the
**per-query** mean-square error

```math
\Delta=
\mathbb E_{A,Z,(a,b)\sim\operatorname{Unif}(B^2)}
 [\widehat R_Z(a,b)-A_{ab}]^2.                       \tag{V.1}
```

If two kernels differ in `r` entries, then

```math
\|R_A-R_D\|_{L^2(B^2)}^2={r\over Q^2}.
```

Consequently the exact inverse-Hamming modulus is

```math
\boxed{\kappa={1\over Q^2}}.                         \tag{V.2}
```

Since `kappa N=1`, Theorem 5.1 gives

```math
\boxed{
I(A;Z)\ge
Q^2\left[1-g\left(\min\{4\Delta,1\}\right)\right].
}                                                     \tag{V.3}
```

The factor four would be different if `Delta` denoted the unnormalized sum
over entries: for total squared error `Delta_tot=Q^2 Delta`, the argument of
`g` is `4Delta_tot/Q^2`.  The corollary should therefore say “uniform-query
average MSE” every time it displays (V.3).

The curve is sharp for `0<=Delta<=1/4`.  Pass the `Q^2` independent bits
through independent binary symmetric channels of crossover `p<=1/2` and
decode each bit by its posterior mean.  A Bernoulli bit then contributes
MSE `p(1-p)`, while

```math
I(A;Z)=Q^2[1-h_2(p)],
\qquad
g(4p(1-p))=h_2(p).
```

Thus equality holds in (V.3).  At `Delta=1/4`, the constant decoder `1/2`
uses no information and attains the zero-rate endpoint.  The cutoff for
larger declared distortions is therefore also correct.

This materially **extends**, but does not pointwise dominate, Theorem FG.4.
FG.4 is a deterministic worst-case packing theorem: uniform error below
`1/2` for every kernel and context forces the full `Q^2` bits.  Equation
(V.3) instead permits randomized encoders and decoders and controls Bayesian
average MSE with a sharp continuous rate curve.  In particular, every fixed
`Delta<1/4` still costs a positive constant fraction of `Q^2` bits even when
exact kernel recovery is not required.  Under FG.4's stronger worst-case
uniform-error hypothesis, FG.4 remains the sharper lossless conclusion.

Finally, the same scope restriction as FG.4 is essential.  The source ranges
over **all** binary boundary kernels, realized by an arbitrary factor on the
endpoint pair.  No `Q^2` lower bound follows for kernels generated by a fixed
local language or restricted transfer semigroup until that realizable class
contains a comparably separated high-entropy source.
