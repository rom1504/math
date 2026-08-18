# Independent audit: bounded-degree carriers for the actual product shadow

Status: **proved, with a strict scope qualification**.  This note was
derived independently of the bounded-degree separator draft and only then
compared with it.  It proves a different statement: every globally optimal
row factor itself is uniformly approximable by the positive part of a
bounded-degree Walsh polynomial, and the resulting restricted product
variational problem approximates the unrestricted product value at
`o(N)` normalized error.  The carrier has polynomial description size at
every fixed accuracy.  It is not yet a child-only algorithm: evaluating or
optimizing its objective may still query the full pressure oracle.

Throughout, logarithms are natural, `U_k` is the uniform law on
`{+-1}^k`, and

```math
D_2(P\Vert U_k)=\log E_{U_k}\left({dP\over dU_k}\right)^2.
```

## 1. Weak-coordinate densities have a bounded Walsh Dirichlet energy

Let `f=dP/dU_k>0`, `E_U f=1`.  Suppose

```math
D_2(P\Vert U_k)\le C,
\qquad
\left|\log f(x)-\log f(x^{(j)})\right|
\le {A\over\sqrt k}                                  \tag{LD.1}
```

for every one-bit flip.  Write

```math
f=\sum_{S\subseteq[k]}\widehat f(S)\chi_S,
\qquad g_d=\Pi_{\le d}f.
```

**Lemma LD.1 (uniform Walsh truncation).**  Put

```math
\tau_k=\tanh {A\over2\sqrt k}.
```

Then

```math
\boxed{
 \sum_S|S|\widehat f(S)^2\le k\tau_k^2e^C
 \le {A^2e^C\over4},
 \qquad
 \|f-g_d\|_2
 \le {Ae^{C/2}\over2\sqrt{d+1}}.}                  \tag{LD.2}
```

*Proof.*  For two adjacent values `a,b>0`, (LD.1) gives

```math
{|a-b|\over a+b}
=\tanh {|\log a-\log b|\over2}\le\tau_k.
```

Thus `(a-b)^2<=2tau_k^2(a^2+b^2)`.  Flip invariance of
`U_k` yields

```math
E_U(f-f^{(j)})^2\le4\tau_k^2E_Uf^2.
```

The exact cube Dirichlet identity is

```math
\sum_{j=1}^kE_U(f-f^{(j)})^2
=4\sum_S|S|\widehat f(S)^2.
```

Use `E_Uf^2<=e^C`, `tanh z<=z`, and then divide the
Dirichlet bound by `d+1` on levels strictly above `d`. `square`

This argument is about the density, not its logarithm.  In particular, it
really does approximate the probability factor rather than merely exposing
one direction that separates two factors.

## 2. Positivity and normalization do not destroy the approximation

The polynomial `g_d` need not be nonnegative.  Define

```math
z_d=E_U(g_d)_+,
\qquad
h_d={(g_d)_+\over z_d}.                             \tag{LD.3}
```

Because `E_Ug_d=E_Uf=1`,

```math
z_d=1+E_U(g_d)_-\ge1.
```

Let

```math
\eta_d={Ae^{C/2}\over2\sqrt{d+1}},
\qquad K=e^{C/2}.
```

**Lemma LD.2 (positive presented carrier).**  The function `h_d` is a
probability density and

```math
\boxed{
 \|h_d-f\|_2
 \le \rho_d:=\eta_d(1+K),
 \qquad
 \|h_d\|_2\le K.}                                 \tag{LD.4}
```

*Proof.*  Since `f>=0`, pointwise projection onto the positive half-line
gives

```math
\|(g_d)_+-f\|_2\le\|g_d-f\|_2\le\eta_d.
```

Also `z_d-1=E(g_d)_-<=eta_d`.  Orthogonal projection is an
`L^2` contraction, so `||g_d||_2<=||f||_2<=K`, and hence
`||(g_d)_+||_2<=K`.  Therefore

```math
\left\|{(g_d)_+\over z_d}-(g_d)_+\right\|_2
={z_d-1\over z_d}\|(g_d)_+\|_2
\le\eta_dK,
```

which proves (LD.4). `square`

Although taking a positive part generally restores high Walsh degrees,
`h_d` still has a finite presentation: it is the normalized positive part
of one degree-`d` polynomial.  At fixed `d` this uses

```math
\sum_{a=0}^d {k\choose a}=k^{O(d)}                 \tag{LD.5}
```

real coefficients and one normalization scalar.  Exact computation of the
normalizer can be hard; (LD.5) is a description-size statement, not an
algorithmic running-time statement.

## 3. A dimension-free entropy modulus on bounded `L^2` balls

The entropy term cannot be controlled by total variation with an alphabet-
size Fannes bound: that would reintroduce a factor `k`.  The `L^2` bound
instead supplies uniform integrability.

**Lemma LD.3 (entropy continuity).**  For every `K_0<\infty` there is a
modulus `\omega_{K_0}(s)\downarrow0` as `s\downarrow0` such that any two
densities `f,h` on any probability space satisfying

```math
\|f\|_2,\|h\|_2\le K_0,
\qquad \|f-h\|_2\le s
```

obey

```math
\boxed{|E f\log f-E h\log h|\le\omega_{K_0}(s).}   \tag{LD.6}
```

The modulus is independent of the size of the underlying space.

*Proof.*  Let `phi(x)=x log x`, with `phi(0)=0`, and fix `R>=e`.
For `x>=R`, `log x/x` is decreasing, so

```math
E[\phi(f)-\phi(f\wedge R)]
\le {\log R\over R}E f^2
\le K_0^2{\log R\over R},                         \tag{LD.7}
```

and the same holds for `h`.  On `[0,R]`, `phi` is uniformly
one-half-Holder.  Explicitly, the finite number

```math
c_R=\sup_{0\le x<y\le R}
 { |\phi(y)-\phi(x)|\over\sqrt{y-x}}
```

exists because `x|log x|=o(sqrt x)` at zero.  Since truncation is
one-Lipschitz,

```math
E|\phi(f\wedge R)-\phi(h\wedge R)|
\le c_RE|f-h|^{1/2}
\le c_Rs^{1/2}.                                    \tag{LD.8}
```

First choose `R` large and then `s` small.  For example, the infimum over
`R>=e` of

```math
2K_0^2{\log R\over R}+c_Rs^{1/2}
```

is a valid modulus. `square`

## 4. Stability of the full row-product variational objective

Let `L:\{\pm1\}^{m\times n}\to\mathbb R` be any potential whose
oscillation under one bridge-bit flip is at most `2u`.  For `lambda>0`
define

```math
\mathcal F(P_1,\ldots,P_m)
=E_{\otimes_iP_i}L
+{1\over\lambda}\sum_{i=1}^mD(P_i\Vert U_n).        \tag{LD.9}
```

Suppose every factor `P_i` has density `f_i` satisfying (LD.1), with common
`A,C`.  Let `H_i` have the positive degree-`d` carrier density `h_i`
from (LD.3).

**Lemma LD.4 (sequential product replacement).**  With `rho_d` from
(LD.4), `K=e^(C/2)`, and the entropy modulus from LD.3,

```math
\boxed{
 |\mathcal F(H_1,\ldots,H_m)-\mathcal F(P_1,\ldots,P_m)|
 \le m\left[u\sqrt n\,\rho_d
       +{1\over\lambda}\omega_K(\rho_d)\right].}  \tag{LD.10}
```

*Proof.*  Replace the factors in order.  At the step replacing row `i`,
average `L` over all other, already or not-yet replaced, factors and call
the resulting row function `G_i(b)`.  A one-bit flip still changes `G_i`
by at most `2u`.  The exact cube Poincare inequality gives

```math
\|G_i-E_UG_i\|_2^2
\le {1\over4}\sum_{j=1}^nE_U(G_i-G_i^{(j)})^2
\le nu^2.                                          \tag{LD.11}
```

The two row densities have equal integral, so the constant part cancels,
and Cauchy--Schwarz gives an energy change at most
`rho_d u sqrt(n)`.  Lemma LD.3 bounds the row entropy change by
`omega_K(rho_d)`.  Sum over all rows. `square`

The centering in (LD.11) is essential.  Bounding by the full row
oscillation would give `2un rho_d` per row and the wrong
`N^(3/2)` total scale.

## 5. Application to the actual optimizing-child law

Fix contracted-temperature optimizing children of orders `m,n`, let
`N=m+n`, and put `u=beta/sqrt(N)`.  In either orientation let
`p^*=tensor_i p_i^*` be a global minimizer of the exact row-product
variational objective.  Its coordinate score equation is

```math
{dp_i^*\over dU_n}(b)
\propto
\exp\{-\lambda E_{p_{-i}^*}L(b,B_{-i})\}.          \tag{LD.12}
```

The actual pressure has one-bit oscillation `2u`.  Consequently

```math
\operatorname {osc}_{b_j}\log {dp_i^*\over dU_n}
\le2\lambda u
\le {2\lambda\beta\over\sqrt n},                  \tag{LD.13}
```

where the final inequality uses `n<=N`.  Theorem 37.19 gives

```math
D_2(p_i^*\Vert U_n)
\le\lambda^2u^2n\le\lambda^2\beta^2.              \tag{LD.14}
```

Thus Lemmas LD.1--LD.4 apply with

```math
A=2\lambda\beta,
\qquad C=\lambda^2\beta^2                         \tag{LD.15}
```

uniformly in the child orders, their exact minimizing signings, the
orientation, and the row.

It is essential to impose the inherited norm bound on the carrier.  Let
`mathcal P_d(K)` be the class of row products whose row densities have the
form

```math
h={g_+\over E_Ug_+},\qquad
\deg g\le d,\qquad E_Ug=1,\qquad\|g\|_2\le K,       \tag{LD.16}
```

and let

```math
V_\lambda^{(d,K)}
=\inf_{P\in\mathcal P_d(K)}
 \left\{E_PL+{1\over\lambda}D(P\Vert U_B)\right\}. \tag{LD.17}
```

**Theorem LD.5 (uniform finite-degree recovery of the actual product
shadow).**  For every fixed `beta,lambda>0`,

```math
\boxed{
 \lim_{d\to\infty}\ \sup_{m,n,A,D,\text{orientation}}
 {V_\lambda^{(d,K)}-V_\lambda^{\rm row}\over m+n}=0,
 \qquad K=e^{\lambda^2\beta^2/2}.}                 \tag{LD.18}
```

where the supremum is over all positive child orders, actual
contracted-temperature minimizing children, and both orientations.  More
quantitatively, the left side at degree `d` is at most

```math
\beta\rho_d+{1\over\lambda}\omega_K(\rho_d),    \tag{LD.19}
```

with the constants (LD.15).  Indeed, `m<=N` and
`u sqrt(n)<=beta` in (LD.10).

If

```math
\mathcal I_\lambda^{\leftarrow,(d,K)}
=\lambda(V_\lambda^{(d,K)}-V_\lambda),            \tag{LD.20}
```

then

```math
\boxed{
0\le
\mathcal I_\lambda^{\leftarrow,(d,K)}
-\mathcal I_\lambda^{\leftarrow}
\le\lambda N
 \left[\beta\rho_d
 +{1\over\lambda}\omega_K(\rho_d)\right].}      \tag{LD.21}
```

Therefore the asymptotic rates of both remaining product branches are
recovered by a hierarchy of fixed-degree carriers:

```math
{\mathcal I_\lambda^{\leftarrow}\over N}
=\lim_{d\to\infty}
 {\mathcal I_\lambda^{\leftarrow,(d,K)}\over N},
\qquad
{\mathcal J-\mathcal I_\lambda^{\leftarrow}\over N}
=\lim_{d\to\infty}
 {\mathcal J-\mathcal I_\lambda^{\leftarrow,(d,K)}\over N}.    \tag{LD.22}
```

uniformly up to the explicit error in (LD.21).  In particular, any fixed
positive extensive rate in either branch is visible at some fixed degree:
if `J-I^leftarrow>=cN`, choose `d` so that the right side of (LD.21) is
less than `cN/2`; then
`J-I^(leftarrow,(d,K))>=cN/2`.  The analogous assertion for an extensive
`I^leftarrow` follows from the same two-sided approximation.

Equation (LD.18) is stronger than a mere bounded-degree separator: it
approximates the value of the globally optimal product shadow.  It also
does not decompose the objective into independently paid Walsh channels;
all coefficients of every positive polynomial density are optimized
jointly before the single expectation of `L` is taken.

## 6. What this proves, and what it does not

This is a genuine representation reduction for the actual children.
At fixed thermodynamic parameters and fixed response accuracy, the row
product search is reduced from `m2^n` density values to

```math
m\sum_{a=0}^d{n\choose a}=N^{O_{\beta,\lambda,\varepsilon}(1)}   \tag{LD.23}
```

real parameters.  It cannot be subsumed by the low-transport ceiling:
degree-`d` Walsh polynomials may use all row bits and have macroscopic
fibre transport.  Nor does it separately pay scalar channels: the
positive part couples all retained coefficients nonlinearly.

The inherited mean and norm conditions in (LD.16), especially the norm
bound, are indispensable.  If a carrier is defined merely as the
normalized positive part of an unrestricted polynomial, degree one already
contains every point mass: for a prescribed `x in {+-1}^n`, the positive
part of

```math
g_x(b)=\sum_jx_jb_j-(n-1)
```

is nonzero only at `b=x`.  That unrestricted class can silently recover
full bridge maximization.  In contrast, every density in (LD.16) obeys

```math
\|h\|_2\le K,
\qquad
P(E)\le K\sqrt{U_n(E)}.                            \tag{LD.24}
```

Thus an event carrying `P`-mass at least `a` has uniform mass at least
`a^2/K^2`; a hard support contains at least `2^n/K^2` row words.  The
normed carrier is quantitatively diffuse and cannot encode a deterministic
bridge row.

There is also a finite-bit version.  Keep the constant coefficient of
`g_d` equal to one and quantize its other `D=\sum_{a=1}^d{n\choose a}`
coefficients within the Euclidean ball of radius `sqrt(K^2-1)`, so that
their total `ell_2` error is at most `delta`.  Parseval
gives `||g_d-\widetilde g_d||_2<=delta`.  Positive projection is an
`L^2` contraction, both positive-part normalizers are at least one, and

```math
\left\|{(g_d)_+\over E(g_d)_+}
-{(\widetilde g_d)_+\over E(\widetilde g_d)_+}\right\|_2
\le\delta(1+K).                                    \tag{LD.25}
```

Applying LD.3--LD.4 again shows that a fixed sufficiently fine
quantization changes the normalized objective by arbitrarily little.
The Euclidean ball has a `delta`-net of cardinality at most
`(1+2K/delta)^D`, so this costs `O(D\log(K/delta))` bits per row at fixed
`d`.
Thus the reduction is finite-information as well as finite-dimensional,
although it is still not an efficient evaluation algorithm.

Three limitations remain.

1. `mathcal P_d(K)` is a polynomial-description class, not automatically a
   polynomial-time class.  Its normalizers and objective expectations can
   still be hard to compute exactly.
2. The theorem decides the branch only after optimizing the degree-`d`
   carrier against the actual pressure oracle.  It does not produce its
   coefficients from a small child statistic.  In particular, the pressure
   response surface can still contain exponential cross-row information;
   the present reduction compresses the factor/carrier representation, not
   that response oracle.
3. The convergence in (LD.18) is an iterated limit: first choose a fixed
   response accuracy (hence a fixed degree), then take large orders.  An
   explicit `o(N)` approximation along one sequence uses `d=d_N\to
   \infty`, and its parameter count need not remain polynomial with a
   uniform exponent.

Thus the rigorous new smallest missing lemma is narrower but still
nontrivial:

> **Finite-degree child closure.**  For one sufficiently large fixed
> degree `d`, evaluate or bound the jointly optimized carrier value in
> (LD.17) from an optimizer-specific child state that is strictly smaller
> than the full pressure landscape, with error small enough to decide one
> of the extensive product branches.

Theorem LD.5 says that no genuinely high Walsh degree is needed merely to
represent an extensive product retuning or reverse-product rate.  What is
missing is child-side closure of the finitely many retained modes, not
control of an exponential row table.  This is a factor-representation
reset only, not a proved child-closure theorem.
