# Joint finite-fibre action: an exact SDP bound and its conference floor

## Status

This note starts from the scalable obstruction to scalar/rank-one
microkernel decompositions and keeps the finite-fibre action as one object.
It proves an exact operator-space inequality for every symmetric microkernel.
No scalar atom expansion, separately paid left/right channel, polarization,
or parent-spin enumeration is used.

The state is genuinely compressed: one polynomial-size signed-elliptope SDP
on the child order.  The result is nevertheless negative for landing.  The
SDP has a signing-independent conference floor, so it can never certify a
normalized constant below `1/2`.  A skew-Clifford variant has the same floor
by a separate Frobenius argument.  These are scalable no-go theorems for the
two joint norm interfaces, not for all nonlocal finite-fibre constructions.

Use the one-copy normalization

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
\operatorname{cap}(A)=\max_x|H_A(x)|.                \tag{1}
```

## 1. The signed-elliptope state

For a symmetric zero-diagonal order-`n` matrix `A`, define

```math
\begin{aligned}
\Gamma(A)=\max\ &\left|
 \sum_{i<j}a_{ij}(P_{ij}-Q_{ij})\right|,\\
\text{subject to }&P\succeq0,\quad Q\succeq0,
\quad \operatorname{diag}(P+Q)=\mathbf1.
\end{aligned}                                        \tag{2}
```

This is not the full Boolean or microtype maximum.  It has two `n` by `n`
positive-semidefinite variables and `n` affine constraints, hence polynomial
description size.  Every Boolean cut is feasible by taking `P=xx^T,Q=0`, so

```math
\Gamma(A)\ge\operatorname{cap}(A).                   \tag{3}
```

The two matrices in (2) share one pointwise norm budget.  They are not
bounded separately and then added.

This state is distinct from the ordinary degree-two endpoints used in
`nonlocal_lift_joint_sos_inheritance.md`.  In that notation,

```math
\begin{aligned}
U_1(A)&=\max_{G\succeq0,\,\operatorname{diag}G=1}
             \sum_{i<j}a_{ij}G_{ij},\\
L_1(A)&=\min_{G\succeq0,\,\operatorname{diag}G=1}
             \sum_{i<j}a_{ij}G_{ij},\\
R_1(A)&=\max\{U_1(A),-L_1(A)\},\qquad
W_1(A)={U_1(A)-L_1(A)\over2}.                        \tag{2a}
\end{aligned}
```

Taking `P=G,Q=0` in (2) shows `Gamma(A)>=R_1(A)`, but equality need not
hold.  `Gamma` optimizes one *indefinite* contraction Gram matrix `P-Q`,
whereas `R_1` uses one positive correlation matrix.  For the saved first
order-eight minimizer class, for example, numerical evaluation gives
`R_1=11.313708...` but `Gamma=14.422205...`.  Thus the lift mapping below is
not a relabeling of `R_1` or `W_1`, although all three share a conference-scale
floor.

Strong SDP duality gives the useful equivalent form

```math
\boxed{
\Gamma(A)=\min_{y\in\mathbb R^n}\sum_i y_i
\quad\text{subject to}\quad
\operatorname{diag}(y)-{A\over2}\succeq0,
\quad
\operatorname{diag}(y)+{A\over2}\succeq0.}          \tag{4}
```

Indeed, for the positive sign in (2), the Lagrangian coefficients of `P`
and `Q` are respectively `diag(y)-A/2` and `diag(y)+A/2`.  The negative sign
interchanges the same two constraints.  Slater feasibility follows from a
sufficiently large constant `y`.

## 2. Exact joint lifting theorem

Let `R` be any symmetric real `k` by `k` matrix and let
`rho=||R||_op`.  Let `D_0` be a zero-diagonal order-`k` signing.  If `R` is a
full-sign matrix, then

```math
B=A\otimes R+I_n\otimes D_0                         \tag{5}
```

is an order-`nk` signing.  The following theorem applies more generally even
when `R` is not full sign.

> **Verified joint-action theorem.** For every child signing `A`,
>
> ```math
> \boxed{
> \operatorname{cap}(B)
> \le k\rho\,\Gamma(A)+n\operatorname{cap}(D_0).}    \tag{6}
> ```

To prove it, write a parent spin as Boolean microvectors
`x_i in {+-1}^k` and put `v_i=x_i/sqrt(k)`.  Set `O=R/rho`; it is a
self-adjoint contraction.  Write its positive and negative parts as

```math
O=O_+-O_-,\qquad O_+,O_-\succeq0,qquad O_++O_-=|O|\preceq I. \tag{7}
```

Let `P,Q` be the Gram matrices of `O_+^(1/2)v_i` and
`O_-^(1/2)v_i`.  Then

```math
P_{ij}-Q_{ij}=\langle v_i,Ov_j\rangle,
\qquad \operatorname{diag}(P+Q)\le\mathbf1.          \tag{8}
```

Add mutually orthogonal private coordinates to `P` to fill the diagonal
slack in (8).  This changes no off-diagonal entry and makes `(P,Q)` feasible
for (2).  The complete inter-fibre action is therefore bounded at once:

```math
\begin{aligned}
\left|\sum_{i<j}a_{ij}x_i^{\mathsf T}Rx_j\right|
&=k\rho\left|\sum_{i<j}a_{ij}(P_{ij}-Q_{ij})\right|\\
&\le k\rho\,\Gamma(A).                               \tag{9}
\end{aligned}
```

The internal term in (5) is `sum_i H_(D_0)(x_i)` and has absolute value at
most `n cap(D_0)`, proving (6).

If `R` is full sign, then `||R||_F=k`, so

```math
\rho\ge\sqrt k.                                      \tag{10}
```

Equality in (10) means that `R` is Hadamard.  In particular, for a symmetric
Hadamard microkernel the best coefficient supplied by (6) is exactly

```math
\operatorname{cap}(B)
\le k^{3/2}\Gamma(A)+n\operatorname{cap}(D_0).       \tag{11}
```

This is the requested correct-scale finite-fibre inequality with a state
strictly smaller than full parent optimization.

## 3. The exact conference floor

Assume in this section that `A` is a full zero-diagonal signing, so
`a_ij in {+-1}` for every `i ne j`.  The compressed state in (4) then has
an unavoidable leading gap.  Let
`Y=diag(y)` be any feasible dual point and put `S=sum_i y_i`.  The constraints
force every `y_i>0` when `n>=2`.  Define

```math
C={1\over2}Y^{-1/2}AY^{-1/2}.                        \tag{12}
```

Congruence in (4) gives `-I preceq C preceq I`, hence

```math
\operatorname{tr}(C^2)le n,
\qquad
\sum_{i\ne j}{1\over y_i y_j}\le4n.                \tag{13}
```

On the other hand, Cauchy--Schwarz over the `n(n-1)` ordered pairs gives

```math
[n(n-1)]^2
\le
\left(\sum_{i\ne j}y_i y_j\right)
\left(\sum_{i\ne j}{1\over y_i y_j}\right).         \tag{14}
```

For fixed `S`,

```math
\sum_{i\ne j}y_i y_j
=S^2-\sum_i y_i^2
\le {n-1\over n}S^2.                                \tag{15}
```

Combining (13)--(15) proves the universal theorem

```math
\boxed{\Gamma(A)\ge {n\sqrt{n-1}\over2}.}           \tag{16}
```

All equality conditions are rigid.  Equality in (15) makes `Y` constant;
equality in (13) makes `C^2=I`.  Thus equality in (16) occurs exactly when

```math
A^2=(n-1)I,                                          \tag{17}
```

that is, when `A` is a symmetric conference matrix.  Conversely, for such
an `A`, the constant choice `Y=sqrt(n-1) I/2` is feasible in (4), so equality
does hold.

Combining (10), (11), and (16), the joint certificate can never return a
leading normalized upper constant below `1/2`:

```math
k\rho\,\Gamma(A)
\ge k^{3/2}{n\sqrt{n-1}\over2}.                     \tag{18}
```

This is a limitation of the **certificate value**, not a lower bound on the
actual lifted Boolean cap.  It says that the signed elliptope has already
forgotten exactly the Boolean information that could place an optimizer
below the conference scale.

For the desired landing statement, (11) would require

```math
\Gamma(A_n)-\operatorname{cap}(A_n)=o(n^{3/2})       \tag{19}
```

on exact minimizers, with a power-saving version for a summable recurrence.
But (16) shows that (19) would itself imply

```math
M_n\ge {n\sqrt{n-1}\over2}-o(n^{3/2}),               \tag{20}
```

and hence the sharp asymptotic lower constant `1/2`.  Thus this additional
statistic is compressed but its required landing theorem is not a simpler
upper-bound obligation; it would settle the problem through the stronger
lower-bound route.

## 4. Exact finite separation and numerical audit

The obstruction is already strict on an exact minimizer.  The saved order-six
matrix is symmetric conference and has

```math
A_6^2=5I,
\qquad \operatorname{cap}(A_6)=5,
\qquad \Gamma(A_6)=3\sqrt5=6.7082039\ldots .         \tag{21}
```

Thus the state is demonstrably not full parent maximization: it is a
polynomial SDP and differs from the exact Boolean answer at order six.

After the theorem was fixed, the SDP was evaluated on the saved certified
minimizers through order 14.  Selected values are:

| `n` | exact `M_n` | numerical `Gamma(A)` | `Gamma-M_n` |
|---:|---:|---:|---:|
| 5 | 4 | 5.590170 | 1.590170 |
| 6 | 5 | 6.708204 | 1.708204 |
| 8, class 0 | 10 | 14.422205 | 4.422205 |
| 8, class 1 | 10 | 12.000000 | 2.000000 |
| 10 | 13 | 20.615528 | 7.615528 |
| 12 | 18 | 23.600782 | 5.600782 |
| 14 | 21 | 25.238859 | 4.238859 |

These SDP values are numerical diagnostics, not exact certificates.  The
caps and conference identities are independently checked exactly.  The
universal theorem (16) is analytic and does not depend on these solves.

## 5. Skew-Clifford evasion also stops at the conference floor

A skew microkernel evades the diagonal-range obstruction from the scalar
atom audit.  It does not evade the global quadratic-form floor.

Let `C_0` be a skew conference matrix of even order `k`:

```math
C_0^{\mathsf T}=-C_0,
\qquad C_0C_0^{\mathsf T}=(k-1)I.                   \tag{22}
```

Let `T` be any skew sign matrix of order `n`.  The joint two-channel matrix

```math
B_0=A\otimes I_k+T\otimes C_0                       \tag{23}
```

is symmetric.  Every cross-fibre block is full sign: its diagonal entries
come from `A` and its off-diagonal entries from `T tensor C_0`.  Adding
`I_n tensor D_0` fills the zero diagonal-fibre blocks.

Since `C_0^2=-(k-1)I`, real orthogonal block diagonalization of `C_0` gives
the exact spectral identity

```math
\boxed{
\|B_0\|_{\rm op}
=\left\|A+i\sqrt{k-1}\,T\right\|_{\rm op}.}          \tag{24}
```

Hence the full joint spectral certificate is

```math
\operatorname{cap}(B_0+I\otimes D_0)
\le {nk\over2}
\left\|A+i\sqrt{k-1}\,T\right\|_{\rm op}
+n\operatorname{cap}(D_0).                           \tag{25}
```

This keeps all Clifford channels inside one Hermitian norm.  Nevertheless,
every off-diagonal entry of the Hermitian matrix in (24) has squared modulus
`k`, so

```math
\left\|A+i\sqrt{k-1}\,T\right\|_F^2=kn(n-1).
```

The operator norm is at least the Frobenius norm divided by `sqrt(n)`, and
therefore

```math
{n\over2\sqrt k}
\left\|A+i\sqrt{k-1}\,T\right\|_{\rm op}
\ge {n\sqrt{n-1}\over2}.                            \tag{26}
```

Thus the effective `k^(3/2)`-scale state in (25) has exactly the same
conference floor as (16), even when `T` is chosen nonlocally from `A`.
Equality requires the complex conference identity

```math
\left(A+i\sqrt{k-1}T\right)^2=k(n-1)I.              \tag{27}
```

There is an independent quotient warning.  Because `C_0` is skew,
`1^T C_0 1=0`, so a cross block in (23) has total sum only `ka_ij`, not
`k^(3/2)a_ij`.  The construction replaces the missing child signal with a
new skew baseline at the conference scale; it does not amplify a possibly
better child optimizer.

## 6. Scope and surviving target

The proved no-go covers:

1. every symmetric full-sign microkernel controlled through one
   self-adjoint-contraction/signed-elliptope norm;
2. arbitrary symmetric Hadamard kernels, which are optimal for that bound;
3. skew-conference/Clifford encodings controlled through their complete
   Hermitian operator norm; and
4. arbitrary nonlocal choice of the skew macro channel `T` in (23).

It does not cover a genuinely Boolean higher-order state that cuts off the
universal degree-two pseudo-solutions while retaining polynomial or otherwise
subexponential complexity.  However, the existing exact duplication theorem
in `nonlocal_lift_joint_sos_inheritance.md` shows that every standard fixed
SOS level is inherited by an exact compressed lift with the full
`k^(3/2)` factor.  Such a relaxation must already have a project-scale-small
integrality gap on the child; the lift cannot repair it.

The next falsifiable theorem must therefore supply a new signing-specific
global inequality, not another universal operator norm.  A useful statement
would be a compressed Boolean state `Theta(A)` satisfying both

```math
\operatorname{cap}(A\otimes R+I\otimes D_0)
\le k^{3/2}\Theta(A)+O_k(n^{3/2-\delta}),
\qquad
\Theta(A)-\operatorname{cap}(A)=O(n^{3/2-\delta}),   \tag{28}
```

for exact minimizers and some fixed `delta>0`.  The first clause is joint
composition; the second is landing.  The conference floor proves that
`Theta` cannot be the signed elliptope, a spectral/Clifford norm, or any
standard fixed-level duplication-closed SDP without an additional
minimizer-specific theorem.

## Reproduction

```text
.venv/bin/python computations/audit_joint_finite_fibre_sdp.py \
  --output computations/results/joint_finite_fibre_sdp_audit.json
```

The result file records exact cap re-evaluation, exact conference-identity
checks, numerical SDP feasibility residuals, source hashes, and a canonical
payload hash.
