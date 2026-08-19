# The centered channel form of the exact cross-order defect

Status: **exact identity and proved channel-level no-go**.  The identity
uses actual thermal optimizers.  The countermodel is an analogous binary
code with actual optimizers, not the complete-graph cut code; it rules out a
generic Blackwell/data-processing closure rather than an optimizer-specific
theorem for the original problem.

## 1. Exact cancellation of the quadratic bridge cost

For a hollow sign matrix `A` of order `k`, write

```math
\phi_A(u)=\log\overline Z_A(u),
\qquad
\overline Z_A(u)=2^{-k}\sum_x\cosh(uH_A(x)),
```

and put

```math
K_k=\binom k2,
\qquad
\ell(u)=\log\cosh u,
\qquad
R_A(u)=\phi_A(u)-K_k\ell(u).
```

The minimized pressure and its exact cross-order difference are

```math
P_k(\beta)=\min_A\phi_A(\beta/\sqrt k),
\qquad
E_{m,n}(\beta)=P_{m+n}(\beta)-P_m(\beta)-P_n(\beta).    \tag{1.0}
```

Fix `N=m+n`, `t=beta/sqrt(N)`, `s_m=beta/sqrt(m)`, and
`s_n=beta/sqrt(n)`.  Let `A,D` be exact minimizers defining
`P_m(beta),P_n(beta)` at `s_m,s_n`.  The joint bridge/orientation output
identity gives

```math
P_N(\beta)\le P_m(\beta)+P_n(\beta)+G_{m,n},             \tag{1.1}
```

where

```math
G_{m,n}=mn\ell(t)-\Delta_A-\Delta_D-
        D_{\rm KL}(U\Vert\Pi),                           \tag{1.2}
```

`Pi` is the exact positive bridge/orientation output law and `U` is
uniform on its `2^(mn+1)` bridge/orientation outputs.  When the child
minimizers are not unique, `G`, `Pi`, and the two `Gamma` terms below refer
to the selected pair.  Identity (1.1) holds for every selected minimizing
pair; a bound for one such pair is enough for the upper recurrence.  Here

```math
\Delta_A=\phi_A(s_m)-\phi_A(t),
\qquad
\Delta_D=\phi_D(s_n)-\phi_D(t).
```

Define the centered radial losses

```math
\Gamma_A=R_A(t)-R_A(s_m),
\qquad
\Gamma_D=R_D(t)-R_D(s_n).                               \tag{1.3}
```

They are nonnegative.  Indeed, independently flip the edges of `A` by
signs of mean

```math
\alpha={\tanh t\over\tanh s_m}.
```

The binary-channel identity is

```math
\mathbb E_\xi e^{R_{A\odot\xi}(s_m)}=e^{R_A(t)}.         \tag{1.4}
```

Since `A` minimizes `phi_A(s_m)`, and the subtracted baseline is independent
of `A`, every term inside the average is at least `e^{R_A(s_m)}`.  This
proves `Gamma_A>=0`; the other child is identical.

Substitution in (1.2) gives the exact decomposition

```math
\boxed{
G_{m,n}=C^0_{m,n}(\beta)+\Gamma_A+\Gamma_D
        -D_{\rm KL}(U\Vert\Pi),}                         \tag{1.5}
```

with

```math
C^0_{m,n}(\beta)
=mn\ell(t)
 -K_m[\ell(s_m)-\ell(t)]
 -K_n[\ell(s_n)-\ell(t)].                               \tag{1.6}
```

All order-`N` quadratic terms in (1.6) cancel.  More precisely, let

```math
q(u)={u^2\over2}-\log\cosh u.
```

Then

```math
C^0_{m,n}(\beta)
={\beta^2\over4}+K_mq(s_m)+K_nq(s_n)-K_Nq(t).            \tag{1.7}
```

For `u>=0`, the elementary inequalities

```math
0\le q(u)\le {u^4\over12}                               \tag{1.8}
```

follow by integrating `0<=u-tanh(u)<=u^3/3`; evenness handles negative
`u`.  Consequently, for every
split and every positive order,

```math
\boxed{
{\beta^2\over4}-{\beta^4\over24}
\le C^0_{m,n}(\beta)
\le {\beta^2\over4}+{\beta^4\over12}.}                  \tag{1.9}
```

Thus, on comparable splits, the annealed `Theta_beta(N)` bridge term was
not itself the true obstruction: after exact centering and the two thermal
payments its entire quadratic part is `O_beta(1)`.  The only possible
linear term is the joint quantity

```math
\Gamma_A+\Gamma_D-D_{\rm KL}(U\Vert\Pi).                \tag{1.10}
```

This statement has an immediate defect implication.  For any fixed
`delta>0`,

```math
\boxed{
\Gamma_A+\Gamma_D-D_{\rm KL}(U\Vert\Pi)
\le C_\beta N^{1-\delta}
\quad\Longrightarrow\quad
E_{m,n}\le C_\beta N^{1-\delta}
 +{\beta^2\over4}+{\beta^4\over12}.}                    \tag{1.11}
```

Equation (1.11), rather than a bound on any one of its three terms, is the
weakest direct target exposed by this identity.  It is not presently proved
for the original optimizing children.

There is an exact strictly stronger certificate which uses the parent's
freedom to choose the relative child orientation.  For a bridge `B`, let

```math
L_\epsilon(B)=\log\overline Z_N
 \left(\begin{pmatrix}A&B\\B^{\mathsf T}&\epsilon D\end{pmatrix},t\right).
```

Since the parent may choose both `epsilon` and `B`, while (1.5) is their
joint uniform average,

```math
\begin{aligned}
E_{m,n}(\beta)
&\le \mathbb E_B\min_{\epsilon\in\{\pm1\}}L_\epsilon(B)
      -P_m(\beta)-P_n(\beta)\\
&=C^0_{m,n}(\beta)+\Gamma_A+\Gamma_D-D_{\rm KL}(U\Vert\Pi)
  -{1\over2}\mathbb E_B|L_+(B)-L_-(B)|.               \tag{1.12}
\end{aligned}
```

Thus any estimate

```math
\Gamma_A+\Gamma_D-D_{\rm KL}(U\Vert\Pi)
-{1\over2}\mathbb E_B|L_+-L_-|
\le C_\beta N^{1-\delta}                              \tag{1.13}
```

immediately gives the same power-saving defect as (1.11).  Equation (1.12)
is not a proof that orientation contrast is large; it prevents the joint
uniform output identity from charging a relative-orientation bit which the
actual parent can optimize for free.

The contrast admits a second exact decomposition which sharply limits this
possibility.  Put

```math
Z_A=\mathbb E_x\cosh(tH_A(x)),
\quad S_A=\mathbb E_x\sinh(tH_A(x)),
\quad u_A={S_A\over Z_A},
```

and define `Z_D,S_D,u_D` analogously.  With

```math
a_0=(\cosh t)^{mn}Z_AZ_D,
\qquad r_\epsilon(B)={e^{L_\epsilon(B)}\over a_0},
```

fixed-orientation bridge averaging gives

```math
\mathbb E_Br_\epsilon(B)=1+\epsilon u_Au_D.          \tag{1.14}
```

If

```math
V={1\over2}\mathbb E_B|r_+(B)-r_-(B)|,
\qquad
\mathcal C_{\rm or}={1\over2}\mathbb E_B|L_+(B)-L_-(B)|,
```

then direct algebra gives

```math
\boxed{
D_{\rm KL}(U\Vert\Pi)+\mathcal C_{\rm or}
=-\mathbb E_B\log\min\{r_+(B),r_-(B)\}.}             \tag{1.15}
```

Since `E min(r_+,r_-)=1-V` and `V>=|u_Au_D|`, Jensen yields

```math
D+\mathcal C_{\rm or}
\ge-\log(1-V)
\ge-\log(1-|u_Au_D|).                                \tag{1.16}
```

Consequently the exact direct implications are

```math
\boxed{
E_{m,n}\le C^0+\Gamma_A+\Gamma_D+\log(1-V)
\le C^0+\Gamma_A+\Gamma_D+\log(1-|u_Au_D|).}         \tag{1.17}
```

This is the sharp use of the two scalar means in (1.14).  If
`Gamma_A+Gamma_D=lambda N+o(N)`, the scalar signed imbalance can make
(1.17) sublinear only when
`1-|u_Au_D|<=exp{-lambda N+o(N)}`: both children must be exponentially
one-sided.

There is an even more useful obstruction.  Write

```math
h_B={L_+(B)-L_-(B)\over2},
\qquad w_B={r_+(B)+r_-(B)\over2}.
```

The bridge marginal of `Pi` has density `w_B` relative to `U_B`, and

```math
D(U\Vert\Pi)=D(U_B\Vert\Pi_B)+\mathbb E_B\log\cosh h_B,
\qquad
\mathcal C_{\rm or}=\mathbb E_B|h_B|.                \tag{1.18}
```

Because `0<=|h|-log cosh(h)<=log2`, one always has

```math
\boxed{\mathcal C_{\rm or}\le D(U\Vert\Pi)+\log2.}   \tag{1.19}
```

Thus whenever the joint reverse divergence is `o(N)`, orientation selection
can save only `o(N)`.  In particular, under the conditional regular
actual-minimizer assumptions of Section 4, (4.3)--(4.5) and (1.19) give

```math
\mathbb E_B\min_\epsilon L_\epsilon(B)-2P_r(\beta)
=\gamma(\beta)r+o(r).                                \tag{1.20}
```

The free orientation therefore does not repair the linear certificate on
that branch.  This is not a lower bound on the globally minimized parent
pressure; it is a scoped no-go for the orientation-aware composition of the
selected children.

Finally, scalar imbalance is not forced by exact minimizer stationarity.
At order four, the 48 signings with energy histogram
`{-4:2,-2:4,0:4,2:4,4:2}` are exact pressure minimizers at every positive
temperature.  Their partition function is
`Z_good=cosh(2u)cosh(u)^2`; the other 16 signings satisfy

```math
16(Z_{\rm bad}-Z_{\rm good})
=8(\cosh(2u)-1)^2(\cosh(2u)+1)>0.                    \tag{1.21}
```

These minimizers have `S_A=0` identically, yet

```math
R_A(u)=\log\cosh(2u)-4\log\cosh u,
\qquad
R_A'(u)={-4\tanh^3u\over1+\tanh^2u}<0,               \tag{1.22}
```

so their radial loss `Gamma_A` is strictly positive.  Signed partition
imbalance cannot universally pay the centered radial loss even for exact
all-temperature minimizers.

## 2. Bayesian interpretation

Let `Q=(tau x_i x_j)_(i<j)` be the antipodal cut word and pass every
coordinate independently through the binary symmetric channel with
`rho=tanh(u)`.  If `p_(k,u)` is the output law and `U_k` is uniform on edge
signings, then

```math
{dp_{k,u}\over dU_k}(A)
={\overline Z_A(u)\over(\cosh u)^{K_k}}
=e^{R_A(u)}.                                           \tag{2.1}
```

Thus `A` is a least-likely channel output, `Gamma_A` is the increase of its
centered log likelihood under BSC degradation from `s_m` to `t`, and `Pi`
is the conditional bridge/orientation output.  The desired estimate (1.10)
would say that bridge reverse information pays for the two degradation
losses, up to a power saving.

That statement is **not** a generic consequence of Blackwell ordering,
data processing, low latent rate, or rank-one joining.

The pointwise Blackwell identity already has the wrong direction.  If
`p_t=T p_s`, where `T` is the doubly stochastic degrading BSC kernel, and
`a` minimizes `p_s`, set

```math
k_a(z)=T(a\mid z),
\qquad
\widehat k_a(z)={k_a(z)p_s(z)\over p_t(a)}.
```

Then direct expansion of relative entropy gives

```math
\log{p_t(a)\over p_s(a)}
=D(k_a\Vert\widehat k_a)
 +\mathbb E_{k_a}\log{p_s(z)\over p_s(a)}.              \tag{2.2}
```

The second term is nonnegative by minimality, so the backward-kernel
divergence is **at most** the radial loss.  The usual reverse-KL chain rule
also averages (2.2) over uniform outputs, whereas (1.10) selects deepest
outputs.  Standard data processing therefore proves neither the direction
nor the pointwise estimate required by (1.11).

## 3. A rank-one binary-code countermodel with actual minimizers

Take `r>=5` with `r=1 mod 4`.  Partition `K_r=binom(r,2)` coordinates into `r` blocks
`I_1,...,I_r`, each of the even size `(r-1)/2`.  Let
`tau,X_1,...,X_r` be independent fair signs and define the antipodal latent
encoder

```math
q_e=\tau X_i\qquad(e\in I_i).                           \tag{3.1}
```

For an output signing `a`, its uncentered channel partition is

```math
Z_a(u)=\mathbb E_{\tau,X}e^{u\langle a,q(\tau,X)\rangle}
=\prod_{i=1}^r\cosh\!\left(u\sum_{e\in I_i}a_e\right).
                                                                    \tag{3.2}
```

For the product equality, absorb the common fair sign `tau` into the
independent fair vector `X`.

Every factor is at least one.  Hence every signing `a_*` which is balanced
on each `I_i` has

```math
Z_{a_*}(u)=1=\min_aZ_a(u)                               \tag{3.3}
```

at **every** temperature.  These are exact optimizing children, not
surrogates.

Join two such order-`r` systems using one common antipodal sign `tau` and
`r^2` bridge coordinates.  The three latent blocks are

```math
(q^A_e,q^D_f,q^{\rm br}_{ij})
=(\tau X_i,\tau Y_j,\tau X_iY_j),                       \tag{3.4}
```

where `X,Y` are independent.  The total coordinate count is

```math
2K_r+r^2=K_{2r},                                        \tag{3.5}
```

so this has the same quadratic coordinate scale as the motivating model.
For balanced children all internal fields vanish identically, and every
sign bridge `B` has parent partition

```math
Z_B(t)=\mathbb E_{\tau,X,Y}e^{t\tau X^{\mathsf T}BY}
=\mathbb E_{X,Y}e^{tX^{\mathsf T}BY}
=\mathbb E_X\prod_{j=1}^r
 \cosh\!\left(t(B^{\mathsf T}X)_j\right).               \tag{3.6}
```

Jensen's inequality and the fact that every `(B^T X)_j` has the law of a
sum `S_r` of `r` independent signs give

```math
\log Z_B(t)\ge r\,\mathbb E\log\cosh(tS_r).             \tag{3.7}
```

Set `N=2r` and `t=beta/sqrt(N)`.  Since

```math
\mathbb ES_r^2=r,
\qquad
\mathbb ES_r^4=3r^2-2r,
```

Paley--Zygmund yields

```math
\Pr\{|S_r|\ge\sqrt{r/2}\}\ge {1\over12}.              \tag{3.8}
```

On this event `t|S_r|>=beta/2`.  Therefore, uniformly over **every** bridge,

```math
\boxed{
\log Z_B(\beta/\sqrt{2r})
\ge {r\over12}\log\cosh(\beta/2)
={N\over24}\log\cosh(\beta/2).}                       \tag{3.9}
```

Both exact child optima in (3.3) are zero.  Thus, with these exact
optimizing internal children held fixed, this rank-one, latent-rate
`Theta(N)` code has a compulsory linear bridge-completion defect:

```math
\min_B\log Z_B(\beta/\sqrt{2r})\ge c_\beta N,
\qquad
c_\beta={1\over24}\log\cosh(\beta/2)>0.                \tag{3.10}
```

This is not a lower bound on the globally minimized order-`2r` pressure:
that minimization may replace the two internal output signings.  It is the
fixed-optimizing-child defect relevant to a generic channel theorem of the
form (1.11).

In the notation of (1.5), (3.9) equivalently forces

```math
D_{\rm KL}(U\Vert\Pi)
\le\Gamma_A+\Gamma_D-c_\beta N+O_\beta(1).              \tag{3.11}
```

Hence the reverse inequality needed in (1.11) fails by `Theta(N)` in an
explicit family with exact all-temperature minimizing children, a
linearly-sized latent space, BSC degradation, and a rank-one bridge.

The missing ingredient for the original problem must therefore use the
triangle/cut algebra or another optimizer-specific property.  No abstract
data-processing theorem at the level of (2.1) can supply it.

## 4. The regular actual-minimizer branch has a linear defect (conditional)

There is also a sharp conditional conclusion inside the original signing
problem.  It uses the fixed-temperature replica-symmetric theorem and the
Bernoulli spectral extremizer already verified in
[`conference_reverse_kl_fixed_temperature_obstruction.md`](conference_reverse_kl_fixed_temperature_obstruction.md)
and
[`symmetrized_r_transform_bernoulli_extremizer.md`](symmetrized_r_transform_bernoulli_extremizer.md).

Fix `0<beta<sqrt(2)/6`, take `N=2r`, and suppose exact order-`r` pressure
minimizers `A_r,D_r` have all of the following properties.

1. Their normalized matrices have compact symmetric variance-one limiting
   spectral laws.
2. They satisfy the fixed-power diagonal/off-diagonal delocalization
   hypotheses of the cited replica-symmetric theorem.
3. Along a coupled iid-bridge sequence, the parents are almost surely
   eventually inside the strict high-temperature operator-norm regime with
   a uniform margin, and their pressures are uniformly integrable.
4. Along that sequence the iid bridge is almost surely asymptotically free
   from the child blocks in the precise fixed-power sense used in the
   archived conference proof, with the uniform integrability needed to
   pass to expected pressure.

Then exact minimality, the all-order conference upper bound, and the strict
Bernoulli spectral extremizer force the two child limiting laws to be
`(delta_(-1)+delta_1)/2`.  Hence

```math
{1\over r}\phi_{A_r}(\beta/\sqrt r),
{1\over r}\phi_{D_r}(\beta/\sqrt r)
\longrightarrow\psi(\beta),                           \tag{4.1}
```

and their contracted pressures converge to `psi(beta/sqrt(2))`, where

```math
\psi(c)={1\over4}\left[
 \sqrt{1+4c^2}-1
 -\log\!\left({1+\sqrt{1+4c^2}\over2}\right)
\right].                                               \tag{4.2}
```

The verified bridge free-convolution calculation then gives

```math
\boxed{
{D_{\rm KL}(U\Vert\Pi)\over r}\longrightarrow0,
\qquad
{G_{r,r}(\beta)\over r}\longrightarrow\gamma(\beta)>0,}
                                                               \tag{4.3}
```

with

```math
\gamma(\beta)={\beta^2\over4}
-2\psi(\beta)+2\psi(\beta/\sqrt2).                    \tag{4.4}
```

Strict positivity follows by differentiating (4.4), and near zero
`gamma(beta)=3 beta^4/16+O(beta^6)`.  Since (1.9) is bounded,

```math
\boxed{
\Gamma_{A_r}+\Gamma_{D_r}-D_{\rm KL}(U\Vert\Pi)
=\gamma(\beta)r+o(r).}                                 \tag{4.5}
```

Thus the sufficient estimate (1.11) is false by a linear amount on the
regular/delocalized branch **even when the children are the actual exact
minimizers**.  The conclusion is conditional only because the project has
not proved that exact minimizers obey assumptions 1--4.

This creates a genuine architecture fork.  Uniform-output reverse-KL
composition can have a sublinear defect only by exploiting an
irregular/localized optimizer phase or by selecting exponentially rare
bridges outside the uniform-output comparison.  Proving spectral
regularity of the minimizers would falsify this composition route rather
than complete it.
