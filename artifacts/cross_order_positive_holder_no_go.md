# Positive Hölder separation has an unavoidable linear cross-order floor

Status: **proved scalable no-go for a broad direct-recurrence class on
actual own-scale pressure minimizers**.  Any positive generalized-Hölder
argument which separates the two child Hamiltonians and the bridge into
independently paid factors has a positive linear certificate at every fixed
`beta>2.554944595...` on all sufficiently large balanced splits.  The
result includes arbitrary choices of the bridge and arbitrary Hölder
weights, including weights depending on the order and temperature.

This is a lower bound on the numerical value of that upper-bound
certificate, not on the true cross-order defect.  It does not apply to a
signed interpolation or to a genuinely joint inequality which preserves
cancellation before separating the three channels.

## 1. The complete positive-Hölder family

For a hollow signing `A`, write

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
\phi_A(u)=\log\mathbb E_x\cosh(uH_A(x)).             \tag{1.1}
```

Fix `N=m+n`, `t=beta/sqrt(N)`, two child signings `A,D`, an `m` by `n`
sign bridge `B`, and a relative orientation `epsilon`.  The parent pressure
is

```math
L_\epsilon(B)
=\log\mathbb E_{\tau,x,y}
 e^{t\tau\{H_A(x)+\epsilon H_D(y)+x^TBy\}}.          \tag{1.2}
```

Let `alpha,delta,gamma>0` with

```math
\alpha+\delta+\gamma=1.                              \tag{1.3}
```

Generalized Hölder, with exponents `1/alpha,1/delta,1/gamma`, gives

```math
\boxed{
L_\epsilon(B)
\le \alpha\phi_A(t/\alpha)
   +\delta\phi_D(t/\delta)
   +\gamma\psi_B(t/\gamma),}                        \tag{1.4}
```

where

```math
\psi_B(u)=\log\mathbb E_{x,y}e^{u x^TBy}.            \tag{1.5}
```

The orientation disappears from the second factor because `cosh` is even;
the bridge law is symmetric, so its augmented and ordinary exponential
moments agree.  Formula (1.4) also covers zero weights by the corresponding
finite-support Laplace limits.

Any finite iteration of ordinary positive Hölder which ultimately pays
the three channels separately reduces to (1.4): multiply the conditional
weights along the Hölder tree.  The three resulting leaf weights are
nonnegative and sum to one.  Thus (1.4), rather than any particular choice
of exponents, is the full class under audit.

This remains true if an argument first splits one Hamiltonian into several
positive copies.  For every convex `f` with `f(0)=0`, its perspective

```math
g(a,z)=a f(z/a),\qquad a>0,
```

is convex and positively homogeneous, hence subadditive.  Therefore

```math
\sum_j a_j f(z_j/a_j)
\ge\left(\sum_j a_j\right)
 f\left({\sum_jz_j\over\sum_ja_j}\right).           \tag{1.4a}
```

Consolidating all positive copies of each of the three channels can only
decrease the resulting certificate.  Since their total coefficients in the
original Hamiltonian are all `t`, the consolidated lower envelope is
exactly (1.4).  Thus repeated temperature splitting, or a deeper positive
Hölder tree, cannot evade the floor proved below.

Now take exact own-scale pressure minimizers, so

```math
P_m(\beta)=\phi_A(\beta/\sqrt m),
\qquad
P_n(\beta)=\phi_D(\beta/\sqrt n).
```

Define the Hölder certificate

```math
\mathcal H_{\alpha,\delta,\gamma,B}
=\alpha\phi_A(t/\alpha)+\delta\phi_D(t/\delta)
 +\gamma\psi_B(t/\gamma)-P_m(\beta)-P_n(\beta).
                                                               \tag{1.6}
```

Parent minimization and (1.4) give the required immediate defect arrow:

```math
\boxed{
E_{m,n}(\beta)
\le\mathcal H_{\alpha,\delta,\gamma,B}.}             \tag{1.7}
```

In particular, a proof that some member of this family is
`O(N^(1-eta))` would prove the desired cross-order estimate with exactly
that exponent.

## 2. Three finite-order extremal lower bounds

Let

```math
M(A)=\max_x|H_A(x)|.
```

The two configurations `x,-x` realizing `M(A)` and
`cosh z>=e^|z|/2` give, for every `u>=0`,

```math
\phi_A(u)\ge uM(A)-m\log2.                           \tag{2.1}
```

Consequently, for every positive Hölder weight,

```math
\boxed{
\alpha\phi_A(t/\alpha)
\ge tM(A)-\alpha m\log2.}                           \tag{2.2}
```

The same bound holds for `D`.

For a balanced `r` by `r` bridge let

```math
S_r=\varepsilon_1+\cdots+\varepsilon_r,
\qquad \mu_r=\mathbb E|S_r|.
```

Integrating one shore, applying Jensen, and using
`log cosh z>=|z|-log2` gives, for every bridge and every `u>=0`,

```math
\psi_B(u)\ge ur\mu_r-r\log2.                         \tag{2.3}
```

Hence

```math
\boxed{
\gamma\psi_B(t/\gamma)
\ge tr\mu_r-\gamma r\log2.}                        \tag{2.4}
```

The leading terms in (2.2) and (2.4) are independent of the Hölder
weights.  Because the three weights sum to one, their entropy losses add to
exactly one `r log2` on a balanced split.

## 3. Uniform linear floor on actual balanced children

Take `m=n=r`, `N=2r`, and let `A,D` be any exact order-`r` own-scale
pressure minimizers.  Put

```math
M_r=\min_C\max_x|H_C(x)|,
\qquad U_r={M_r\over r^{3/2}}.                       \tag{3.1}
```

Every selected child has `M(A),M(D)>=M_r`.  On the other hand, evaluating
the pressure minimum at a ground-cap minimizer gives

```math
P_r(\beta)\le {\beta M_r\over\sqrt r}=\beta U_r r. \tag{3.2}
```

Substitute (2.2), (2.4), and (3.2) in (1.6), with
`t=beta/sqrt(2r)`.  Uniformly over every bridge and every admissible
Hölder triple,

```math
\boxed{
{\mathcal H_{\alpha,\delta,\gamma,B}\over r}
\ge
\beta\left\{(\sqrt2-2)U_r+{\mu_r\over\sqrt{2r}}\right\}
-\log2.}                                             \tag{3.3}
```

This is an exact finite-order certificate floor.

The rigorous all-order construction gives

```math
\limsup_{r\to\infty}U_r\le{1\over2},                \tag{3.4}
```

while the central limit theorem gives

```math
{\mu_r\over\sqrt r}\longrightarrow\sqrt{2/\pi}.   \tag{3.5}
```

Since `sqrt(2)-2<0`, (3.3)--(3.5) imply

```math
\boxed{
\liminf_{r\to\infty}
 {\mathcal H_{\alpha_r,\delta_r,\gamma_r,B_r}\over r}
\ge c_H\beta-\log2}                                \tag{3.6}
```

for every sequence of weights and bridges, where

```math
c_H={1\over\sqrt\pi}-{2-\sqrt2\over2}
   ={1\over\sqrt\pi}-1+{1\over\sqrt2}
   =0.27129636473430385\ldots .                      \tag{3.7}
```

Thus the critical temperature is

```math
\boxed{
\beta_H={\log2\over c_H}
=2.5549445944061366\ldots .}                         \tag{3.8}
```

For every fixed `beta>beta_H`, there are `c_beta>0` and `r_beta` such
that, for all `r>=r_beta`, every bridge and every positive Hölder
separation satisfies

```math
\boxed{
\mathcal H_{\alpha,\delta,\gamma,B}
\ge c_\beta(2r).}                                    \tag{3.9}
```

Combining (1.7) and (3.9), this entire direct-recurrence class retains a
linear upper certificate and cannot establish an `o(N)` defect in the
low-temperature range needed for the zero-temperature limit.

## 4. Why the obstruction is structural

There is a simple exponent-budget interpretation.  To evaluate the left
child factor at no hotter than its own scale requires

```math
\alpha\ge {t\over\beta/\sqrt m}=\sqrt{m/N},          \tag{4.1}
```

and similarly `delta>=sqrt(n/N)`.  But on every nontrivial split,

```math
\sqrt{m/N}+\sqrt{n/N}>1,                             \tag{4.2}
```

whereas positive Hölder has only the unit budget
`alpha+delta+gamma=1`.  Thus it must overheat at least one child before it
has paid anything for the bridge.  Equations (2.1)--(3.9) quantify this
incompatibility and show that exponent tuning cannot hide it.

The theorem is broader than ordinary polarization or a particular
three-factor estimate: every positive iterated separation tree has the
same final exponent budget.  It leaves open precisely the architectures
which do not pay all channels separately—for example, a signed interpolation,
a joint same-switch inequality, or a direct exceptional bridge exploiting
cancellation inside the common logarithm.
