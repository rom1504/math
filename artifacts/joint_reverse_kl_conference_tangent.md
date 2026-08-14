# Conference tangent obstruction to joint reverse-KL compensation

Status: **exact output-law theorem, exact conference-implementation Taylor
obstruction, and reproducible finite diagnostics**. The Taylor coefficient is
uniform algebraically, but its remainder is not uniform in the order.
Conference signings are also not known to minimize the finite-temperature
child pressure uniformly. Thus this note obstructs a structured
implementation and perturbative proof; it does not falsify the
minimizer-optimized criterion at one fixed positive temperature.

## 1. The joint output law

For an order-$m$ signing $A$, put

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
\overline Z_m(A,t)=2^{-m}\sum_x\cosh(tH_A(x)).
\tag{1.1}
\]

Let $A,D$ have orders $m,n$, let $N=m+n$, and set
$t=\beta/\sqrt N$. Independently draw

\[
\mu_{A,t}(\tau,x)
=\frac{2^{-m-1}e^{t\tau H_A(x)}}{\overline Z_m(A,t)},
\qquad
\mu_{D,t}(\sigma,y)
=\frac{2^{-n-1}e^{t\sigma H_D(y)}}{\overline Z_n(D,t)},
\tag{1.2}
\]

and iid signs $\eta_{ij}$ with

\[
\Pr(\eta_{ij}=a)=\frac{e^{ta}}{2\cosh t}.
\tag{1.3}
\]

Output the relative orientation and bridge

\[
\epsilon=\tau\sigma,\qquad
B_{ij}=\tau x_i y_j\eta_{ij}.
\tag{1.4}
\]

Call this output law $\Pi_{A,D,t}$. If $U$ is uniform on the same
$2^{mn+1}$ outputs and

\[
S_{\epsilon,B}=
\begin{pmatrix}A&B\\B^{\mathsf T}&\epsilon D\end{pmatrix},
\]

direct substitution gives

\[
\boxed{
\Pi_{A,D,t}(\epsilon,B)
=\frac{\overline Z_N(S_{\epsilon,B},t)}
 {2^{mn+1}(\cosh t)^{mn}
  \overline Z_m(A,t)\overline Z_n(D,t)}.}
\tag{1.5}
\]

Taking the uniform average of its logarithm proves the quenched identity

\[
\boxed{
\begin{aligned}
\mathbb E_U\log\overline Z_N(S_{\epsilon,B},t)
={}&\log\overline Z_m(A,t)+\log\overline Z_n(D,t)\\
&+mn\log\cosh t-D_{\rm KL}(U\Vert\Pi_{A,D,t}).
\end{aligned}}
\tag{1.6}
\]

This is genuinely joint: it does not pay left and right scalar channels
separately, and cancellation occurs before the logarithm.

Choose $A,D$ to minimize their pressures at
$s_m=\beta/\sqrt m$ and $s_n=\beta/\sqrt n$, and define

\[
\Delta_A=\log\overline Z_m(A,s_m)-\log\overline Z_m(A,t),
\quad
\Delta_D=\log\overline Z_n(D,s_n)-\log\overline Z_n(D,t).
\tag{1.7}
\]

Equation (1.6) implies

\[
F_N(t)\le F_m(s_m)+F_n(s_n)+mn\log\cosh t
-\Delta_A-\Delta_D-D_{\rm KL}(U\Vert\Pi).
\tag{1.8}
\]

Hence same-scaled-temperature approximate subadditivity follows from

\[
\boxed{
\Delta_A+\Delta_D+D_{\rm KL}(U\Vert\Pi)
\ge mn\log\cosh t-\omega(N).}
\tag{1.9}
\]

A balanced-tree argument proves convergence if, for example,
$\omega(N)=O(N^{1-\delta})$; more generally it is enough that
$\sum_k\omega(2^k)/2^k<\infty$.

## 2. The exact information-loss obligation

Let $V$ be uniform on all latent variables above, let $Q_t$ be their product
tilted law, and let $T$ be the deterministic output map. Then
$T_\#V=U$ and $T_\#Q_t=\Pi$. Reverse data processing supplies only

\[
D_{\rm KL}(U\Vert\Pi)
\le D_{\rm KL}(V\Vert Q_t)
=\log\overline Z_m(A,t)+\log\overline Z_n(D,t)
 +mn\log\cosh t.
\tag{2.1}
\]

This has the wrong direction for (1.9). The reverse chain rule is exact. Put

\[
\Lambda_t=
\mathbb E_{\omega\sim U}
D_{\rm KL}\bigl(V(\cdot\mid\omega)\Vert
Q_t(\cdot\mid\omega)\bigr).
\tag{2.2}
\]

Then

\[
\Lambda_t=
\log\overline Z_m(A,t)+\log\overline Z_n(D,t)
+mn\log\cosh t-D_{\rm KL}(U\Vert\Pi),
\tag{2.3}
\]

so (1.9) is equivalent to

\[
\boxed{
\Lambda_t\le
\log\overline Z_m(A,s_m)+\log\overline Z_n(D,s_n)+\omega(N).}
\tag{2.4}
\]

Thus the missing statement is an upper bound on the conditional reverse
information lost by a many-to-one noisy rank-one output map. Generic
information inequalities do not provide it.

## 3. Exact conference Taylor theorem

Take two copies of a symmetric conference signing $A$ of order $r$:

\[
A^2=(r-1)I,\qquad a_{ij}\in\{\pm1\}\quad(i\ne j).
\tag{3.1}
\]

The calculation applies at every order satisfying (3.1); Paley matrices
supply an infinite family of such orders. Put $N=2r$,
$t=\beta/\sqrt{2r}$, $s=\beta/\sqrt r$, and define

\[
\mathcal M_r(\beta)=
2\log\overline Z_r(A,s)
-\mathbb E_{\epsilon,B}\log\overline Z_{2r}
\left(
\begin{pmatrix}A&B\\B^{\mathsf T}&\epsilon A\end{pmatrix},t
\right).
\tag{3.2}
\]

This is the conference implementation of (1.9):
$\mathcal M_r(\beta)\ge-\omega(2r)$. It is the minimizer-optimized
criterion only at orders and temperatures where the conference child is an
appropriate pressure minimizer.

**Theorem 3.1 (linear fourth-order shortfall).** At every order satisfying
(3.1),

\[
\boxed{
\mathcal M_r(\beta)
=-\frac{\beta^2}{4}
+\left[-\frac{(r-1)(3r-5)}{16r}+\frac1{48}\right]\beta^4
+O_r(\beta^6).}
\tag{3.3}
\]

The coefficient of $\beta^4$ is

\[
-\frac{9r^2-25r+15}{48r}
=-\frac{3r}{16}+O(1)
=-\frac{3N}{32}+O(1).
\tag{3.4}
\]

The joint reverse KL is invisible at this order. More precisely, with
$\rho=\tanh t$,

\[
\boxed{
D_{\rm KL}(U\Vert\Pi_{A,A,t})
=\frac32{r\choose2}^{\!2}\rho^8+O_r(\rho^{10})
=\frac{3(r-1)^2}{128r^2}\beta^8+O_r(\beta^{10}).}
\tag{3.5}
\]

Its first nonzero scaled coefficient is bounded and tends to $3/128$,
whereas the adverse fourth-order coefficient in (3.3) is linear in $N$.

### Proof

Write

\[
W_A(\rho)=(\cosh t)^{-\binom r2}\overline Z_r(A,t).
\]

The first nonconstant term of the even-Eulerian expansion is
$T_4(A)\rho^4$, where

\[
T_4(A)=\sum_{C_4}\prod_{e\in C_4}a_e
=\frac{\|A^2\|_F^2-r(r-1)(2r-3)}8.
\]

Using (3.1),

\[
T_4(A)=-\frac{r(r-1)(r-2)}8,
\tag{3.6}
\]

and therefore

\[
\log\overline Z_r(A,u)
=\frac{r(r-1)}4u^2
-\frac{r(r-1)(3r-5)}{24}u^4+O_r(u^6).
\tag{3.7}
\]

Identity (1.6) rewrites (3.2) as

\[
\mathcal M_r(\beta)
=2\{\log\overline Z_r(A,s)-\log\overline Z_r(A,t)\}
-r^2\log\cosh t+D_{\rm KL}(U\Vert\Pi).
\tag{3.8}
\]

Substituting (3.7) and
$\log\cosh t=t^2/2-t^4/12+O(t^6)$ gives (3.3), once the KL term is
known to be $O_r(\beta^8)$.

For that term, (1.5) and the Eulerian expansion give

\[
\frac{\Pi(\epsilon,B)}{U(\epsilon,B)}
=1+X_{\epsilon,B}\rho^4+O_r(\rho^6),
\qquad \mathbb E_UX=0.
\tag{3.9}
\]

Here $X$ is the sum of mixed four-cycles. Cycles using three vertices in one
block have, for a fixed pair of bridge edges, coefficient

\[
\sum_{k\ne i,j}a_{ik}a_{kj}=(A^2)_{ij}=0.
\]

For each $i<j$ and $\alpha<\gamma$, the remaining $2+2$ cycles give three
distinct Walsh characters:

\[
B_{i\alpha}B_{j\alpha}B_{j\gamma}B_{i\gamma},\quad
\epsilon a_{ij}a_{\alpha\gamma}B_{i\alpha}B_{j\gamma},\quad
\epsilon a_{ij}a_{\alpha\gamma}B_{i\gamma}B_{j\alpha}.
\tag{3.10}
\]

They are mutually orthogonal, so Parseval gives

\[
\mathbb E_UX^2=3{r\choose2}^2.
\tag{3.11}
\]

The likelihood ratio has uniform mean one at every $\rho$. Expanding
$-\mathbb E_U\log(\Pi/U)$ therefore gives first nonzero term
$\frac12\mathbb E X^2\rho^8$, proving (3.5).

## 4. Reproducible conference diagnostics

The program
[audit_joint_reverse_kl_conference.py](../computations/audit_joint_reverse_kl_conference.py)
exactly sums every child and parent spin state for each bridge. It exhausts
all $16$ bridges and both orientations at $r=2$, and uses a seeded uniform
bridge sample at $r=6,10$. The canonical result is
[joint_reverse_kl_conference.json](../computations/results/joint_reverse_kl_conference.json).

The exact run was:

~~~bash
.venv/bin/python computations/audit_joint_reverse_kl_conference.py \
  --orders 2 6 10 --samples-10 32 \
  --output computations/results/joint_reverse_kl_conference.json
~~~

The default supplies 8,192 order-six bridge draws; the explicit option
supplies 32 order-ten draws; and the program fixes the seed.
The canonical JSON SHA-256 is
$51c53f873359b6d141dd5b0a8ec0edcf3d7145de89f7f0066147d310600761e5$.

The table reports $\mathcal M_r(\beta)$. Parentheses contain the Monte Carlo
standard error of the mean parent log pressure; the $r=2$ values are
exhaustive.

| child $r$ | bridge draws | $\beta=0.25$ | $\beta=0.5$ | $\beta=1$ | $\beta=2$ |
|---:|---:|---:|---:|---:|---:|
| 2 | 16 | -0.015665 | -0.063080 | -0.254574 | -0.912557 |
| 6 | 8,192 | -0.018072 (0.000004) | -0.097017 (0.000070) | -0.618799 (0.000968) | -3.154822 (0.006453) |
| 10 | 32 | -0.020572 (0.000054) | -0.128728 (0.000804) | -0.871811 (0.009909) | -3.904860 (0.078696) |

All sampled margins have the adverse sign. At $\beta=1$, the margins per
parent vertex are $-0.06364,-0.05157,-0.04359$. These are finite numerical
evidence, not an asymptotic certificate.

## 5. Scope and judgment

The output-law identity (1.6) escapes the scalar-pressure no-go theorem, but
generic information contraction does not pay the bridge: data processing
has the wrong direction and (2.4) exposes a full-fibre obligation.
Conference children have a linear adverse fourth-order tangent while their
joint reverse KL first responds at eighth order with bounded coefficient.

Two limitations are essential:

1. Conference children are not proved to minimize the scaled
   finite-temperature pressure at all fixed $\beta$ along an infinite
   family. The theorem obstructs this structured implementation, not the
   minimizer-optimized criterion (1.9).
2. The remainder in (3.3) is $O_r(\beta^6)$, not uniform in $r$. Fixed-order
   Taylor coefficients do not establish failure at one fixed $\beta>0$ as
   $r\to\infty$.

Thus the rigorous result rules out coefficientwise or uniformly
perturbative compensation with sublinear defect; the fixed-$\beta$ finite
data give evidence in the same direction. A scalable falsifier would need a
uniform high-temperature theorem showing, for some fixed small $\beta$,

\[
D_{\rm KL}(U\Vert\Pi)=o(r),
\qquad
\mathcal M_r(\beta)\le-c_\beta r
\]

along a conference family.

Finally, (1.9) is not demonstrably simpler than parent composition. It
avoids choosing a parent minimizer, but it asks for an average parent log
bound, which is stronger than existence of one good bridge, and its state is
a full $2^{mn+1}$-point output law. Without a bounded-complexity overlap
compression theorem, it is a precise reformulation rather than a verified
reduction.
