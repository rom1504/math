# A norm-bounded low-degree carrier for the actual-child product shadow

Status: **rigorous uniform recovery theorem for the actual optimizing-child
row-product problem**.  The full optimal row factors are not needed even as
a representation: at fixed additive accuracy per row, the optimum is
recovered by normalized positive parts of bounded-degree Walsh polynomials
whose `L^2(U)` norms are uniformly bounded.  The degree depends on the
accuracy and fixed thermodynamic parameters, not on the child orders.

This is stronger than a post hoc separator.  The reduced carrier is defined
without knowing the optimal factors, and optimizing over it recovers the
full row-product variational value to `epsilon N`.  It is not claimed to be
a standard log-linear exponential family.  Projecting a log density and
then exponentiating can destroy tail control; the stable operation is to
project the density and apply a positive-part link.

## 1. The child-independent carrier

Let `U_n` be the fair law on the row cube and let `Pi_(<=d)` denote Walsh
projection onto degrees at most `d`, including the constant.  For `K>=1`,
define

```math
\mathcal C_{n,d,K}
=\left\{
 q_g:
 {dq_g\over dU_n}={g_+\over\mathbb E_Ug_+},\quad
 \deg g\le d,\quad \mathbb E_Ug=1,\quad
 \|g\|_{L^2(U)}\le K
 \right\}.                                         \tag{LDC.1}
```

Here `g_+=max(g,0)`.  Since `E_Ug=1`, its normalizer is at least one.  Thus
every carrier density obeys

```math
\left\|{dq_g\over dU_n}\right\|_2\le K,
\qquad D_2(q_g\Vert U_n)\le2\log K.                \tag{LDC.2}
```

The norm constraint is essential.  Without it, a degree-one affine
threshold with a huge coefficient can isolate one cube atom, so an
unrestricted positive-polynomial family would secretly reconstruct the
forbidden row response.  In (LDC.1), point masses have `L^2(U)` norm
`2^(n/2)` and are excluded at fixed `K`.

The carrier has

```math
D(n,d)=\sum_{a=0}^d{n\choose a}                    \tag{LDC.3}
```

real parameters before its two scalar constraints.  At fixed `d,K`, a
coefficient `ell_2` net of mesh `tau` has at most
`(1+2K/tau)^(D(n,d))` points.  Hence one row carrier has polynomial real
dimension and `N^{O(d)}log(1/tau)`-scale finite description, rather than a
table of length `2^n`.

## 2. Weak log scores have low-degree density carriers

**Lemma LDC.1 (uniform density projection).**  Fix `A,C>0`.  Let `p` be a
strictly positive law on `{+-1}^n`, put `f=dp/dU_n`, and suppose

```math
D_2(p\Vert U_n)\le C,                              \tag{LDC.4}
```

while every bit flip changes `log f` by at most `A/sqrt(n)`.  Put

```math
s_d=\Pi_{(\le d)}f,
\qquad
q_d=q_{s_d}.
```

Then `E_Us_d=1`, `||s_d||_2<=e^(C/2)`, so
`q_d in C_(n,d,e^(C/2))`, and

```math
\boxed{
\left\|{dq_d\over dU_n}-f\right\|_2
\le(1+e^{C/2})
 {A e^{A+C/2}\over\sqrt{2(d+1)}}
=:\varepsilon_d(A,C).}                            \tag{LDC.5}
```

In particular, `epsilon_d(A,C)->0` independently of `n`.

*Proof.*  If `B^(j)` is obtained by flipping bit `j`, then

```math
|f(B)-f(B^{(j)})|
\le e^{A/\sqrt n}{A\over\sqrt n}
     \max\{f(B),f(B^{(j)})\}.
```

Using `max(x,y)^2<=x^2+y^2`, invariance of `U_n` under a bit flip, and
(LDC.4),

```math
\sum_{j=1}^n\mathbb E_U
 |f(B)-f(B^{(j)})|^2
\le2A^2e^{2A+C}.                                   \tag{LDC.6}
```

The Walsh Dirichlet identity therefore gives

```math
\|f-s_d\|_2^2
\le{A^2e^{2A+C}\over2(d+1)}.                       \tag{LDC.7}
```

Orthogonal projection preserves the constant and decreases `L^2` norm,
so `E_Us_d=1` and `||s_d||_2<=||f||_2<=e^(C/2)`.  Put `h=(s_d)_+` and
`z=E_Uh`.  Because `f>=0`, pointwise projection onto the nonnegative ray
gives

```math
\|h-f\|_2\le\|s_d-f\|_2.
```

Moreover `z=1+E_U(s_d)_-`, so
`1<=z<=1+||s_d-f||_2`.  It follows that

```math
\left\|{h\over z}-f\right\|_2
\le\|h-f\|_2+(z-1)\|f\|_2,
```

which together with (LDC.7) proves (LDC.5). `square`

The factor `e^A` in the displayed constant is deliberately conservative.
At physical scale `A` is fixed, so only the `d^(-1/2)` decay matters.

## 3. A dimension-free entropy modulus

**Lemma LDC.2 (`L^2` entropy continuity).**  For every `K<infty` there is
a nondecreasing modulus `omega_K:[0,infty)->[0,infty)` with
`omega_K(t)->0` as `t->0` such that, on every finite probability space,
probability densities `f,g` satisfying

```math
\|f\|_2,\|g\|_2\le K,
\qquad \|f-g\|_2\le t
```

obey

```math
\left|\mathbb E f\log f-\mathbb E g\log g\right|
\le\omega_K(t).                                    \tag{LDC.8}
```

The modulus is independent of the number of atoms.

*Proof.*  Let `phi(x)=x log x`, with `phi(0)=0`.  For `R>=e`, replace
`phi` above `R` by the constant `phi(R)` to obtain a bounded continuous
function `phi_R`.  Since `log x/x` decreases for `x>=e`,

```math
\mathbb E|\phi(f)-\phi_R(f)|
\le 2K^2{\log R\over R},                            \tag{LDC.9}
```

and the same estimate holds for `g`.  (One term bounds `f log f` on
`{f>R}`; the other bounds the constant truncation there using
`Pr(f>R)<=K^2/R^2`.)

For fixed `R`, let `varpi_R` be the ordinary uniform-continuity modulus of
`phi_R`.  Split the probability space according to `|f-g|<=s`.  Since
`||f-g||_1<=t`,

```math
\mathbb E|\phi_R(f)-\phi_R(g)|
\le\varpi_R(s)+2\|\phi_R\|_\infty{t\over s}.       \tag{LDC.10}
```

First choose `R` large, then take `s=sqrt(t)` and let `t` tend to zero.
Equations (LDC.9)--(LDC.10) define the required dimension-free modulus.
`square`

One may take an explicit small-`t` rate

```math
\omega_K(t)\le C(1+K^2)t^{1/4}\log(1/t).            \tag{LDC.10a}
```

Indeed, in the preceding proof choose `R=t^(-1/4)` and `s=sqrt(t)`.
On `[0,R]`, direct integration of `1+log x` gives

```math
\varpi_R(s)\le Cs\{1+\log R+|\log s|\},
```

while `||phi_R||_infty<=C R log R`.  Substitution in
(LDC.9)--(LDC.10) proves (LDC.10a), after changing the universal constant.
Combined with (LDC.5), this gives the conservative explicit rate

```math
\eta_d=O_{\beta,\lambda}(d^{-1/8}\log d).            \tag{LDC.10b}
```

No optimality of this exponent is asserted.

## 4. Uniform recovery of the actual row-product optimum

Fix actual contracted-temperature minimizing children, either orientation
and row direction, and a split `m+n=N`.  Put `u=beta/sqrt(N)` and let
`L(B)` be the actual parent bridge pressure.  For a row product
`P=tensor_iP_i`, write

```math
\mathcal F(P)=\mathbb E_PL+{1\over\lambda}
 \sum_iD(P_i\Vert U_n).                            \tag{LDC.11}
```

Let

```math
V^{\rm row}=\min_{P\ {m row\ product}}\mathcal F(P)             \tag{LDC.12}
```

and, for a child-independent carrier,

```math
V^{(d,K)}=
\inf_{P_i\in\mathcal C_{n,d,K}}\mathcal F(\otimes_iP_i).          \tag{LDC.13}
```

**Theorem LDC.3 (fixed-degree recovery of the actual product shadow).**
Fix `beta,lambda>0` and put

```math
C_0=\lambda^2\beta^2,
\qquad K_0=e^{C_0/2},
\qquad A_0=2\lambda\beta.                          \tag{LDC.14}
```

For every `d`, define `epsilon_d=epsilon_d(A_0,C_0)` by (LDC.5), and set

```math
\eta_d=\beta\epsilon_d
 +{1\over\lambda}\omega_{K_0}(\epsilon_d).             \tag{LDC.15}
```

```math
\boxed{
0\le V^{(d,K_0)}-V^{\rm row}\le m\eta_d,
\qquad \eta_d\longrightarrow0.}                  \tag{LDC.16}
```

This holds uniformly over the actual optimizing children, orientations,
and splits.

*Proof.*  Let `p^*=tensor_i p_i^*` be a global optimal product shadow and
put `f_i=dp_i^*/dU_n`.  AC.17 and the physical bit-flip bound give

```math
\operatorname {osc}_{b_j}\log f_i\le2\lambda u
\le {A_0\over\sqrt n}.                             \tag{LDC.17}
```

The last inequality uses only `n<=N`.  AC.18 gives

```math
D_2(p_i^*\Vert U_n)
\le n\log(1+\tanh^2(\lambda u))
\le\lambda^2u^2n\le C_0.                           \tag{LDC.18}
```

Apply LDC.1 to each factor and call the resulting carrier law `q_i`.
It belongs to `C_(n,d,K_0)` and its density differs from `f_i` in `L^2(U)`
by at most `epsilon_d`.

Replace the factors of `p^*` by the `q_i` one at a time.  With all other
rows held at their current laws, the effective row pressure

```math
F_i(b)=\mathbb E[L(B)\mid B_i=b]
```

changes by at most `2u` under one bit flip.  The cube Poincare inequality
therefore gives

```math
\|F_i-\mathbb E_UF_i\|_{L^2(U)}\le u\sqrt n\le\beta.              \tag{LDC.19}
```

Since the density difference has uniform mean zero, Cauchy--Schwarz yields

```math
|\mathbb E_{p_i^*}F_i-\mathbb E_{q_i}F_i|
\le\beta\epsilon_d.                                \tag{LDC.20}
```

Summing the sequential replacements controls the energy term in
(LDC.11) by `m beta epsilon_d`.  Lemma LDC.2 controls each row entropy by
`omega_(K_0)(epsilon_d)`.  Hence

```math
\mathcal F(\otimes_iq_i)-\mathcal F(p^*)\le m\eta_d.
```

The left inequality in (LDC.16) follows because the carrier products form
a subclass of all row products. `square`

The use of pressure Poincare in (LDC.19) is the key physical-scale fact.
A crude row-range/total-variation estimate would lose `sqrt(N)` per row;
the weak bridge coordinates instead make the replacement cost
dimension-free.

## 5. Restricted information identities and branch consequence

Let `q_lambda` be the actual negative-disorder escort and `V_lambda` its
unrestricted Gibbs value.  Define the restricted reverse projection

```math
\mathcal I_{d,K}^{\leftarrow}
=\inf_{P_i\in\mathcal C_{n,d,K}}D(\otimes_iP_i\Vert q_\lambda).
                                                               \tag{LDC.21}
```

The Gibbs identity gives exactly

```math
\boxed{
\mathcal I_{d,K}^{\leftarrow}
=\lambda\{V^{(d,K)}-V_\lambda\}.}                 \tag{LDC.22}
```

Since `I^leftarrow=lambda(V^row-V_lambda)`, Theorem LDC.3 yields

```math
\boxed{
0\le\mathcal I_{d,K_0}^{\leftarrow}
      -\mathcal I^{\leftarrow}
\le\lambda m\eta_d=o_d(1)N.}                     \tag{LDC.23}
```

If `r` is the canonical product and

```math
\mathcal J-\mathcal I^{\leftarrow}\ge\alpha N,
```

choose one fixed `d=d(beta,lambda,alpha)` for which
`lambda m eta_d<=alpha N/2`.  The carrier product constructed above then
satisfies

```math
\boxed{
\lambda\{\mathcal F(r)-\mathcal F(\otimes_iq_i)\}
\ge{\alpha N\over2}.}                              \tag{LDC.24}
```

Thus branch (iii) has a feasible, child-independent, norm-bounded
fixed-degree carrier witness.  The family is operationally smaller than
the full product shadow: it has `mD(n,d)=N^{O(d)}` bounded coefficients,
all its factors retain uniformly bounded `D_2`, and a finite coefficient
net suffices at any fixed objective accuracy by the same sequential
replacement estimate.

Likewise, whenever the unrestricted product shadow has target excess
`o(N)`, a fixed-degree carrier has target excess at most
`o(N)+epsilon N` for arbitrarily prescribed fixed `epsilon>0`.  If one
chooses a degree schedule `d_N->infinity` slowly enough that
`eta_(d_N)=o(1)`, the carrier excess is `o(N)` and its factors still have
the same uniform `D_2` bound, so Theorem 37.19's product-basin mechanism
applies.  This is not a Level-6 recurrence: neither a summable error rate
nor an operationally controlled growing-degree schedule has been proved.

## 6. A sound-and-complete finite-description certificate for branch (iii)

The carrier reduction is not merely an existential change of variables.  A
declared carrier can be checked without storing either an exponential row
table or the complete bridge-pressure table.  Put

```math
D=D(n,d)=\sum_{a=0}^d{n\choose a},
\qquad M=K\sqrt D.                                  \tag{LDC.25}
```

**Proposition LDC.4 (oracle-verifiable carrier certificate).**  Fix
`d,K,beta,lambda` and suppose the coefficients of
`g_1,...,g_m` in (LDC.1) are given.  Then:

1. the constraints `E_Ug_i=1` and `||g_i||_2<=K` are checked directly
   from the Walsh coefficients;
2. exact independent samples from each `q_(g_i)` use rejection sampling
   from `U_n` with expected at most `M` row proposals;
3. given point-evaluation access to the actual pressure `B -> L(B)`,
   `E_(tensor_iq_(g_i))L` is estimated to additive `epsilon N` with
   failure probability at most `zeta` using

   ```math
   O\!\left(\beta^2\epsilon^{-2}\log {2\over\zeta}\right)             \tag{LDC.26}
   ```

   independent bridge samples; and
4. the row entropy terms are estimated to total additive `epsilon N`
   from polynomially many uniform row samples at every fixed `d,K,
   epsilon`.

Consequently, at any fixed extensive margin `alpha>0`, branch (iii) has a
sound-and-complete family of randomized finite-description certificates:

```math
\boxed{
\mathcal J-\mathcal I^{\leftarrow}\ge\alpha N
\quad\Longrightarrow\quad
\text{some fixed-degree carrier certifies gain at least }\alpha N/2,}
                                                               \tag{LDC.27}
```

while every carrier whose checked gain is `cN` proves
`J-I^leftarrow>=cN` by the exact variational identity.  One certificate
uses `mD(n,d)` bounded real coefficients (or their fixed-accuracy
quantizations), one common endpoint comparison, and scalar pressure
evaluations.  It cannot encode bridge point masses: every row factor obeys

```math
q_{g_i}(E)\le K\sqrt{U_n(E)},
\qquad
\max_b q_{g_i}(b)\le {K\sqrt D\over2^n}.            \tag{LDC.28}
```

The restricted projection is also a robust selector between the two linear
product alternatives.  If `J>=eta N`, choose `d` so that the right side of
(LDC.23) is at most `eta N/4`.  Then

```math
\boxed{
\begin{array}{ll}
\mathcal I_{d,K_0}^{\leftarrow}\ge3\eta N/4
 &\Longrightarrow\quad
   \mathcal I^{\leftarrow}\ge\eta N/2,\\[2mm]
\mathcal I_{d,K_0}^{\leftarrow}<3\eta N/4
 &\Longrightarrow\quad
   \mathcal J-\mathcal I^{\leftarrow}>\eta N/4.
\end{array}}                                             \tag{LDC.28a}
```

Thus no unbounded Walsh-degree information in the row factors is needed to
choose between alternatives (ii) and (iii), once a linear canonical phase
has been established.  The remaining oracle issue is evaluation from the
children, not representation of the product competitors.

*Proof.*  Parseval proves item 1.  The Walsh evaluation bound gives

```math
0\le(g_i(b))_+\le ||g_i||_2\sqrt D\le M.
```

Accept a fair row `b` with probability `(g_i(b))_+/M`.  The acceptance
probability is `E(g_i)_+/M>=1/M`, and the accepted law has density
`(g_i)_+/E(g_i)_+`; this proves item 2.

Changing a complete bridge row changes `L` by at most `2un`.  Hoeffding's
lemma over the `m` independent carrier rows therefore gives

```math
\log E_Pe^{s(L-E_PL)}
\le {s^2m u^2n^2\over2}
\le {s^2\beta^2N^2\over2}.                         \tag{LDC.29}
```

The standard subgaussian sample-mean bound proves item 3.  For item 4,
write `z_i=E_U(g_i)_+>=1` and

```math
D(q_{g_i}\Vert U_n)
={E_U[(g_i)_+\log(g_i)_+]\over z_i}-\log z_i.       \tag{LDC.30}
```

Both functions sampled under `U_n` are bounded in terms of `M`:
`0<=(g_i)_+<=M` and
`-e^(-1)<=(g_i)_+log(g_i)_+<=Mlog M` (with the evident
harmless adjustment when `M=1`).  Since `z_i>=1`, ordinary bounded
sampling and a union bound estimate every row contribution to fixed
additive accuracy using a number of samples polynomial in `M`,
`epsilon^(-1)`, and `log(m/zeta)`.  Summing `m<=N` such errors proves
item 4.  Completeness is (LDC.24), and soundness is (SM.1)--(SM.2).
Finally, Cauchy--Schwarz and the same evaluation bound prove (LDC.28).
`square`

This is an information and verification reduction, not a polynomial-time
construction from the children.  A point evaluation of `L` is itself a
finite parent partition-function computation, and finding the good
coefficients is not supplied by the proposition.  What (LDC.27) rules out
is a different obstruction: an extensive coherent retuning cannot require
an exponential row table merely to *represent or verify its witness*.

## 7. What changed in the SML

The former retuning obligation asked for a coherent direction without
solving `m` exponential row tables.  The exact replacement is now:

> **Low-degree carrier decision.**  For one fixed degree selected from the
> desired macroscopic accuracy, determine whether the norm-bounded carrier
> value `V^(d,K_0)` improves the canonical product by order `N`, or prove
> that it cannot.  Equivalently, control the polynomially parameterized
> restricted projection (LDC.21).

This is a strict factor-side and certificate-side reduction: the carrier is
declared before seeing `p^*`, has polynomial rather than exponential row
description, excludes atomic encodings by its norm bound, recovers the full
product value uniformly, and any extensive witness can be verified with a
bounded number of actual-pressure samples once its coefficients are given.
No sector--Gram or other child-only identity presently *generates* or
uniformly bounds those coefficients.  The remaining question is therefore
narrower but still optimizer-specific and nonconvex:

> **Finite-degree child closure.**  Produce or exclude a successful
> norm-bounded carrier from an actual-child state strictly smaller than the
> complete pressure response surface, with a power-saving or fixed-density
> margin.  Equivalently, prove that one fixed-degree carrier response
> functional closes under the optimizer identities, or exhibit a scalable
> actual-child obstruction to such closure.
