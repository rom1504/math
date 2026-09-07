# Independent four-sign Taylor and absolute-energy concentration audit

Date: 2026-09-07. Status: proved. This verifies the director's proposed
scalar-Gaussian concentration step; it is not a convergence theorem.

Let A be a hollow symmetric full-sign matrix of order n. Suppose
G=I+sA is positive semidefinite and |s|<=1/6. Let Z be a centered Gaussian
with covariance G, X_i=sign(Z_i), and H=H_A(X).

## 1. Uniform four-sign expansion, including singular global covariances

For any four distinct indices their covariance has the form R(s)=I_4+sB,
where B is hollow symmetric with off-diagonal signs. Its eigenvalues lie
between 1/2 and 3/2, because ||B||op<=3. This holds even if the full G is
singular. Define f(s)=E product_{i=1}^4 sign(Z_i) for this marginal.

Its density is differentiable three times under the integral, uniformly
for |s|<=1/6: all derivatives are Gaussian density times a bounded-degree
polynomial, dominated by a fixed integrable polynomial times exp(-||z||^2/3).
At zero the log-density expansion has linear term

```math
L_1(z)={1\over2}z^T Bz=\sum_{i<j}B_{ij}z_i z_j,
```

and quadratic term (Tr B^2)/4-z^T B^2z/2. Against the product of four
coordinate signs, constant, linear, and all degree-two polynomial terms
integrate to zero by parity. Only L_1^2/2 survives at quadratic order,
and only its disjoint-edge terms survive. Since E|N(0,1)|=sqrt(2/pi),

```math
f(s)=\left({2\over\pi}\right)^2s^2
 (B_{12}B_{34}+B_{13}B_{24}+B_{14}B_{23})+O(|s|^3),
```

with a universal remainder constant.

An explicit sufficient constant is 55296. To check this, write u=||z||^2
and L for log density. Direct differentiation and ||R^-1||<=2, ||B||<=3 give

```math
|L'|\le12+6u,\quad |L''|\le72+72u,\quad
|L'''|\le864+1296u.
```

The third density derivative is density times
(L')^3+3L'L''+L'''. Under R, the first three moments of u are at most
6,54,648. Substitution bounds the integrated absolute third derivative by
331776. Taylor's theorem divides this by 6, proving the stated remainder.

The exact two-sign formula is E X_iX_j=(2/pi)arcsin(sA_ij). Its product
for the first pairing equals (2/pi)^2s^2 A_ij A_kl+O(s^4). Therefore,
for four distinct i,j,k,l,

```math
\operatorname{Cov}(X_iX_j,X_kX_l)
=\left({2\over\pi}\right)^2s^2
 (A_{ik}A_{jl}+A_{il}A_{jk})+O(|s|^3).
```

The cancellation of the A_ij A_kl pairing is exact at quadratic order.
For example 55297 is a valid uniform covariance remainder constant:
|arcsin(s)-s|<=|s|^3 and |arcsin(s)|<=2|s| on this interval suffice.

## 2. Exact four-cycle bookkeeping

Sum over ORDERED pairs of disjoint undirected edges. The identity is

```math
\sum_{\{i,j\}\cap\{k,l\}=\varnothing}
 A_{ij}A_{kl}(A_{ik}A_{jl}+A_{il}A_{jk})
={1\over2}\operatorname{Tr}_{\rm distinct}(A^4).
```

Each unoriented four-cycle occurs four times in the edge-pair sum and
eight times in the trace. Repeated vertices in Tr A^4 have exactly the
total contribution n(n-1)(2n-3): the possible equalities are i=k or j=l,
with their intersection subtracted once. Thus

```math
\operatorname{Tr}_{\rm distinct}(A^4)
=\operatorname{Tr}(A^4)-n(n-1)(2n-3)\le\operatorname{Tr}(A^4).
```

Shared-edge covariances are bounded in absolute value by 3|s|, from
Cov(X_iX_j,X_iX_k)=E X_jX_k-E X_iX_j E X_iX_k. There are
n(n-1)(n-2) ordered pairs of this kind. Identical-edge variances are at
most one. Combining all three cases proves, for a universal constant C,

```math
\operatorname{Var}H\le {n(n-1)\over2}
 +{2\over\pi^2}s^2\operatorname{Tr}(A^4)
 +C\bigl(n^3|s|+n^4|s|^3\bigr).
```

For instance C=60000 follows from the explicit constants above. No
independence among the four-spin terms has been assumed.

## 3. Cap-only spectral control gives subleading fluctuations

Let lambda be an eigenvalue of maximal absolute value and v a unit
eigenvector. Row Cauchy--Schwarz gives |lambda| ||v||infinity<=sqrt(n-1).
The vector y=|lambda|v/sqrt(n-1) lies in the continuous cube. Independent
rounding of its coordinates shows |H_A(y)|<=Q(A). Consequently

```math
\|A\|_{\rm op}^3\le2(n-1)Q(A).
```

If Q(A)<=C_0 n^(3/2), then ||A||op=O(n^(5/6)). Since Tr A^2=n(n-1),
Tr A^4<=||A||op^2 Tr A^2=O(n^(11/3)). For |s|<=K/sqrt(n), the variance
bound becomes O_{C_0,K}(n^(8/3)). Therefore

```math
\left|\mathbb E|H|-|\mathbb EH|\right|
\le\sqrt{\operatorname{Var}H}=O(n^{4/3})=o(n^{3/2}).
```

The exact mean is n(n-1)arcsin(s)/pi. The bound is uniform over all
signings and admissible s under the two displayed hypotheses.

### 3.1 A stronger elementary interpolation bound

The exponent can be improved without changing the four-sign argument.
For a real matrix with maximum entry magnitude one, let beta be its real
infinity-to-one norm. Its complex infinity-to-one norm is at most 2beta:
write the complex input as a+ib and use
||Aa+iAb||_1<=||Aa||_1+||Ab||_1. Its complex one-to-infinity norm is one.
Interpolation therefore gives ||A||op^2<=2beta.

Here is a direct three-lines verification of that interpolation step.
For unit real vectors u,v, on the strip 0<=Re z<=1 consider

```math
F(z)=\sum_{i,j}A_{ij}\operatorname{sign}(u_i)\operatorname{sign}(v_j)
                  |u_i|^{2z}|v_j|^{2z},
```

omitting terms with a zero coordinate. On Re z=0 its absolute value is at
most 2beta; on Re z=1 it is at most
sum_i |u_i|^2 sum_j |v_j|^2=1. The three-lines inequality at z=1/2 gives
|u^T A v|<=sqrt(2beta). For hollow symmetric A, put a=(x+y)/2 and
b=(x-y)/2. Both lie in the continuous cube and
x^T A y=2H_A(a)-2H_A(b). Independent rounding places each partial energy
in [-R,P], so beta<=2(P+R)=4w(A)<=4Q(A).
Consequently

```math
\|A\|_{\rm op}^2\le8w(A)\le8Q(A),\qquad
\operatorname{Tr}(A^4)\le8w(A)n(n-1).
```

Under the hypotheses of Section 3 the same variance formula thus actually
gives Var H=O(n^(5/2)) and
E|H|=|E H|+O(n^(5/4)). The earlier n^(8/3) variance estimate is valid
and does not require this interpolation strengthening.

## 4. Application to the actual near-minimal interval

Suppose w(A)>=c n^(3/2)-o(n^(3/2)) and Q(A)<=u n^(3/2)+o(n^(3/2)),
where b=2c-u>0. Then each of P,R is at least b n^(3/2)-o(n^(3/2)).
Positivity of I+sA, tested on the appropriate extremal Boolean spin,
gives |s|<=n/(2 min(P,R))=O(n^(-1/2)). Hence the preceding result yields

```math
\limsup {\mathbb E|H_A(\operatorname{sign}Z)|\over n^{3/2}}
\le {1\over2\pi(2c-u)}.
```

For the campaign's certified c and u this is approximately .427687,
strictly below c=.4333221116640807. Thus scalar covariances I+sA cannot
provide a Boolean mean-absolute-energy witness above the width of these
near-minimizers. This excludes this particular witness method, not other
Gaussian covariances, ramp certificates, or convergence mechanisms.
