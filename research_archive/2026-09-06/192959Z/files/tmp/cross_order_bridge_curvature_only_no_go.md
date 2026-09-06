# Bridge curvature gives an exact recurrence interface but not a physical Taylor bound

Status: **proved exact actual-child curvature identity and a scalable
covariance-only method obstruction**.  The first part applies to the actual
own-scale pressure-minimizing children at every order.  It shows exactly what
a curvature estimate would have to prove in order to improve the cross-order
defect.  The second part constructs general centrally symmetric spin laws,
not quadratic-signing Gibbs laws.  It proves that covariance/curvature,
even supplemented by vanishing KL from the fair cube and uniform conditional
min-entropy on every coordinate block, cannot supply the required physical
Taylor remainder.  Thus an optimizer-specific higher-order input is
necessary.

## 1. Exact bridge-amplitude path

For a hollow signing `A` of order `m`, write

```math
\phi_A(u)=\log\mathbb E_x\cosh\{uH_A(x)\}.
```

Put

```math
N=m+n,\qquad t={\beta\over\sqrt N},\qquad
s_m={\beta\over\sqrt m},\qquad s_n={\beta\over\sqrt n},
```

and let `A,D` be exact minimizers defining `P_m(beta),P_n(beta)` at
`s_m,s_n`.  For a relative orientation `epsilon in {+-1}`, a sign bridge
`B`, and a real bridge amplitude `u`, define

```math
L_{\epsilon,B}(u)
=\log\mathbb E_{x,y}
 \cosh\{tH_A(x)+\epsilon tH_D(y)+u x^{\mathsf T}By\}.
\tag{1.1}
```

At `u=t`, this is the pressure of the complete block signing

```math
\begin{pmatrix}A&B\\B^{\mathsf T}&\epsilon D\end{pmatrix}.
```

Put

```math
u_A(t)={\mathbb E_x\sinh(tH_A(x))
              \over\mathbb E_x\cosh(tH_A(x))},
\qquad u_D(t)\ \hbox{ analogously}.
```

Direct multiplication of the two child sector partitions gives

```math
\boxed{
L_{\epsilon,B}(0)
=\phi_A(t)+\phi_D(t)+\log\{1+\epsilon u_A(t)u_D(t)\}.}
\tag{1.2}
```

Let `nu_(epsilon,0)` be the probability law on `(tau,x,y)` proportional to

```math
2^{-m-n-1}\exp\{\tau tH_A(x)+\epsilon\tau tH_D(y)\}.
\tag{1.3}
```

Conditional on `tau=a`, its two spin vectors are independent.  Denote their
correlation matrices by

```math
C_A^a=\mathbb E[xx^{\mathsf T}\mid\tau=a],
\qquad C_D^{\epsilon a}
 =\mathbb E[yy^{\mathsf T}\mid\tau=a],
```

and the sector probabilities by `pi_a^epsilon`.  Every conditional one-spin
mean is zero, because a quadratic Hamiltonian is invariant under global spin
flip.  Differentiating (1.1) at zero therefore gives

```math
L_{\epsilon,B}'(0)=0
```

and

```math
\boxed{
K_\epsilon(B):=L_{\epsilon,B}''(0)
=\sum_{a=\pm1}\pi_a^\epsilon
 \operatorname {Tr}
 \left(B C_D^{\epsilon a}B^{\mathsf T}C_A^a\right).}
\tag{1.4}
```

This is the exact sector-resolved version of the proposed
`Tr(B C_D B^T C_A)` curvature.  Since every correlation matrix has unit
diagonal, averaging a uniform sign bridge kills every off-diagonal monomial
and gives

```math
\boxed{\mathbb E_{B\sim U}K_\epsilon(B)=mn,
\qquad \min_BK_\epsilon(B)\le mn.}                 \tag{1.5}
```

No approximation has entered (1.2)--(1.5).

## 2. The direct cross-order implication

Define the exact physical Taylor remainder

```math
\mathcal T_\epsilon(B)
=L_{\epsilon,B}(t)-L_{\epsilon,B}(0)
 -{t^2\over2}K_\epsilon(B).                         \tag{2.1}
```

Also use the contracted radial payments

```math
\Delta_A=\phi_A(s_m)-\phi_A(t),\qquad
\Delta_D=\phi_D(s_n)-\phi_D(t).                    \tag{2.2}
```

The parent is free to choose `epsilon,B`, so (1.2) and (2.1) give the exact
requested arrow

```math
\boxed{
E_{m,n}(\beta)
\le {t^2\over2}K_\epsilon(B)+\mathcal T_\epsilon(B)
 -\Delta_A-\Delta_D
 +\log\{1+\epsilon u_A(t)u_D(t)\}.}                \tag{2.3}
```

In particular, any hypothesis `P` proving that the right-hand side of
(2.3) is at most `C_beta N^(1-delta)` proves

```math
\boxed{P\quad\Longrightarrow\quad
E_{m,n}(\beta)\le C_\beta N^{1-\delta}.}            \tag{2.4}
```

There is a useful centered form.  Put

```math
R_A(u)=\phi_A(u)-{m\choose2}\log\cosh u,
\qquad
\Gamma_A=R_A(t)-R_A(s_m),
```

and analogously for `D`.  The binary-channel/minimality identity proves
`Gamma_A,Gamma_D>=0`.  Substitution in (2.3) yields

```math
\boxed{
\begin{aligned}
E_{m,n}(\beta)\le{}&C^{(2)}_{m,n}(\beta)
 +\Gamma_A+\Gamma_D
 +{t^2\over2}\{K_\epsilon(B)-mn\}\\
&+\mathcal T_\epsilon(B)
 +\log\{1+\epsilon u_A(t)u_D(t)\},
\end{aligned}}                                      \tag{2.5}
```

where

```math
C^{(2)}_{m,n}(\beta)
={t^2mn\over2}
 -{m\choose2}\{\log\cosh s_m-\log\cosh t\}
 -{n\choose2}\{\log\cosh s_n-\log\cosh t\}.
\tag{2.6}
```

The entire scalar quadratic term is tight:

```math
\boxed{C^{(2)}_{m,n}(\beta)
\le {\beta^2\over4}+{5\beta^4\over48}.}            \tag{2.7}
```

Indeed, the centered-channel constant `C^0` obtained by replacing
`t^2mn/2` with `mn log cosh(t)` is at most
`beta^2/4+beta^4/12`, while

```math
0\le mn\{t^2/2-\log\cosh t\}
\le {mn t^4\over12}\le{\beta^4\over48}.
```

Consequently the weakest curvature-form power-saving target exposed here
is

```math
\boxed{
\begin{aligned}
&\Gamma_A+\Gamma_D
 +{t^2\over2}\{K_\epsilon(B)-mn\}
 +\mathcal T_\epsilon(B)
 +\log\{1+\epsilon u_Au_D\}
 \le C_\beta N^{1-\delta}\\
&\hspace{35mm}\Longrightarrow\quad
E_{m,n}(\beta)
\le C_\beta N^{1-\delta}
 +{\beta^2\over4}+{5\beta^4\over48}.
\end{aligned}}                                      \tag{2.8}
```

Thus optimizing curvature is useful only together with a physical
higher-cumulant estimate.  Merely observing `min_B K<=mn` leaves the two
nonnegative radial losses and the Taylor remainder uncontrolled.

## 3. A covariance-only physical-remainder obstruction

The missing remainder cannot be bounded from covariance, bulk entropy, or
conditional spread alone.  The following theorem is deliberately stated
for general spin laws so that its quantifiers are clear.

**Theorem 3.1 (identity covariance with a linear physical remainder).**
Fix

```math
0<\beta<{4\over\sqrt\pi}.
```

Choose positive constants `c,d` satisfying

```math
c+d<{\beta\over2\sqrt\pi}-{\beta^2\over8},          \tag{3.1}
```

and put `p=1-e^(-d)`.  For every sequence of sign matrices
`B_r in {+-1}^(r by r)`, along `r=2^k-1` there are centrally symmetric,
strictly positive laws `mu_r,nu_r` on `{+-1}^r` such that:

1. both laws have exact identity covariance;
2. `D(mu_r||U_r)+D(nu_r||U_r)=o(1)`;
3. for every coordinate set `S`, every outside value of positive
   probability, and either law `lambda_r`,

   ```math
   \max_z\lambda_r(X_S=z\mid X_{S^c})
   \le(1-p)^{|S|}=e^{-d|S|};                       \tag{3.2}
   ```

4. with independent `X~mu_r,Y~nu_r`, `N=2r`, and
   `t=beta/sqrt(N)`, one has

   ```math
   \boxed{
   \liminf_{r\to\infty}{1\over N}
   \left[
    \log\mathbb E e^{tX^{\mathsf T}B_rY}
    -{t^2\over2}
      \operatorname {Tr}(B_rB_r^{\mathsf T})
   \right]
   \ge {\beta\over2\sqrt\pi}-c-d-{\beta^2\over8}>0.}
   \tag{3.3}
   ```

In particular, every bridge has curvature `r^2`, so minimizing
`Tr(B C_D B^T C_A)` gives no preference at all, while the physical Taylor
remainder can be positive linear.  Notice that the largest possible
right-hand side before subtracting `c+d` is `1/(2 pi)`; hence (3.1)
automatically gives `d<log 2` and therefore `p<1/2`.

### Proof

Let

```math
\kappa_r=\mathbb E\left|\sum_{i=1}^r\varepsilon_i\right|
```

for independent fair signs.  For uniform `x`, every column of `B_r^T x`
has the law of the displayed sum.  Hence

```math
\mathbb E_x\|B_r^{\mathsf T}x\|_1=r\kappa_r.
```

There is therefore a pair `x_*,y_*` with

```math
x_*^{\mathsf T}B_ry_*\ge r\kappa_r.                 \tag{3.4}
```

Here `y_*` is the coordinatewise sign of `B_r^T x_*`.  The elementary
central-limit asymptotic is

```math
{\kappa_r\over\sqrt r}\longrightarrow\sqrt{2/\pi}. \tag{3.5}
```

Index the `r=2^k-1` coordinates by the nonzero vectors of `F_2^k`.  Let
`a in F_2^k` and `b in F_2` be uniform and define the affine-simplex word

```math
X_v^0=x_{*,v}(-1)^{\langle a,v\rangle+b}.            \tag{3.6}
```

It is centrally symmetric, has zero coordinate means, and satisfies

```math
\mathbb E X_v^0X_w^0=0\quad(v\ne w).                \tag{3.7}
```

Pass this word through independent binary noise: multiply every coordinate
by a sign which equals `+1` with probability `1-p`.  Call the resulting law
`sigma_(x_*,p)`.  Noise preserves (3.7), so its covariance is still the
identity.  It is strictly positive, and

```math
\sigma_{x_*,p}(x_*)
\ge {1\over2(r+1)}(1-p)^r.                          \tag{3.8}
```

Now let `q_r=e^(-cr)` and set

```math
\mu_r=(1-q_r)U_r+q_r\sigma_{x_*,p}.                 \tag{3.9}
```

Construct `nu_r` in the same way from `y_*`.  Mixtures preserve central
symmetry, positivity, and identity covariance.  Convexity of relative
entropy in its first argument and the trivial bound
`D(sigma||U_r)<=r log 2` give

```math
D(\mu_r\Vert U_r)\le q_rr\log2=o(1),                \tag{3.10}
```

and similarly for `nu_r`.

For (3.2), condition first on which mixture component and on the affine
word in (3.6).  Inside `S`, the binary noises remain independent after any
observation outside `S`, so every point has conditional probability at most
`(1-p)^|S|`.  The uniform component has point masses `2^(-|S|)`, no larger
than this.  Posterior mixing of these conditional laws preserves the same
upper bound.

Finally, (3.8)--(3.9) imply

```math
\mathbb P\{X=x_*,Y=y_*\}
\ge {e^{-2cr}(1-p)^{2r}\over4(r+1)^2}.              \tag{3.11}
```

Use this one term in the exponential moment and then (3.4):

```math
\log\mathbb E e^{tX^{\mathsf T}B_rY}
\ge tr\kappa_r-2cr+2r\log(1-p)-2\log\{2(r+1)\}.
\tag{3.12}
```

Because `N=2r`, (3.5) gives

```math
{tr\kappa_r\over N}\longrightarrow{\beta\over2\sqrt\pi}.
```

Identity covariance gives

```math
\mathbb E(X^{\mathsf T}B_rY)^2
=\operatorname {Tr}(B_rB_r^{\mathsf T})=r^2,
\qquad {t^2r^2\over2N}={\beta^2\over8}.             \tag{3.13}
```

Since `log(1-p)=-d`, division of (3.12) by `N` and subtraction of
(3.13) prove (3.3).  `square`

## 4. Exact method boundary

Consider any deterministic rule which chooses a bridge using only the two
child sector probabilities and correlation matrices.  Install each law in
Theorem 3.1 identically in two fair sectors.  On this identity-covariance
input, the rule's tie-breaking selects some `B_r`; the theorem then
constructs laws with exactly the same covariance data for which its
quadratic Taylor prediction has a positive linear error.  Therefore

```math
\boxed{
\begin{gathered}
\text{sector covariance data}
+\text{curvature minimization}
+\text{bulk KL/conditional-spread bounds}\\
\not\Longrightarrow\quad
\mathcal T_\epsilon(B)=o(N).
\end{gathered}}                                      \tag{4.1}
```

The quantifier excludes only a rule whose bridge choice and remainder
estimate use no information beyond those data.  A tie-breaker which inspects
higher child correlations or the bridge MGF has left the covariance-only
class.  Likewise, Theorem 3.1 is not an actual-minimizer counterexample:
the constructed laws need not be Gibbs laws of quadratic signings.

For the original recurrence, (2.8) shows the only surviving curvature route
precisely.  It must prove, for the actual optimizing-child laws and for one
jointly chosen `(epsilon,B)`, an optimizer-specific estimate on

```math
\Gamma_A+\Gamma_D+\mathcal T_\epsilon(B)
```

which cancels the curvature deviation and orientation term to `o(N)`.
Neither exact covariance, near-uniform entropy, nor macroscopic conditional
spread can provide that estimate.  This is consistent with the archived
Wishart--Wigner interpolation: its second-moment match yields a direct
power-saving implication, but minimizer optimality controls the endpoint
regrets in the wrong direction.  The present theorem additionally rules out
repairing that mismatch by a covariance-only Taylor expansion at zero
bridge amplitude.
