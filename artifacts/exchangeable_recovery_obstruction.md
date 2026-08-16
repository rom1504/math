# Projective exchangeability cannot realize an extremal signed action object

Date: 2026-08-16.

Status: proved obstruction to the projectively consistent exchangeability
architecture. This does not obstruct nonprojective, order-dependent recovery
laws or deterministic realizers.

## 1. Statement

Let ​\((X_{ij})_{1\le i<j<\infty}\) be a jointly exchangeable symmetric
​\(\{\pm1\}\)-array, and let (X^{(n)}) be its hollow order-(n)
restriction. Suppose

```math
\frac{\lVert X^{(n)}\rVert_{\rm op}}{\sqrt n}
```

is tight. Then the off-diagonal entries are i.i.d. uniform signs. Moreover,
almost surely,

```math
\boxed{
\liminf_{n\to\infty}
\frac{Q(X^{(n)})}{n^{3/2}}
\ge
\frac23\sqrt{\frac2\pi}
=0.5319230405\ldots .}                                    \tag{EX.1}
```

In particular, a single infinite exchangeable signing cannot have
asymptotically extremal finite restrictions, since the rigorous all-order
upper limit for the minima is (1/2).

## 2. Exchangeability plus spectral tightness forces i.i.d. signs

The Aldous--Hoover representation for an exchangeable binary graph array can
be expressed as a random graphon. Conditional on a random symmetric
​\(p:[0,1]^2\to[0,1]\) and i.i.d. uniform vertex labels (U_i), the edges
are independent and

```math
\mathbb P(X_{ij}=1\mid p,(U_k)_k)=p(U_i,U_j).
```

See Diaconis--Janson, *Graph limits and exchangeable random graphs*,
[Sections 5--6](https://arxiv.org/abs/0712.2749).

Put ​\(w=2p-1\). We show that spectral tightness forces (w=0) almost
everywhere, almost surely. If the bounded self-adjoint kernel operator
​\(K_w\) is nonzero, there is a bounded real function (f) such that

```math
c:=\iint f(u)w(u,v)f(v)\,du\,dv\ne0.                       \tag{EX.2}
```

Indeed, a nonzero self-adjoint operator has a nonzero quadratic form, and
bounded functions are dense in (L^2). Set (v_i=f(U_i)). Conditional on
​\(p,(U_i)_i\), the centered edge variables are independent, so the strong
law for bounded (U)-statistics and a conditional variance estimate give

```math
\frac1{n^2}v^{\mathsf T}X^{(n)}v\longrightarrow c,
\qquad
\frac1n\lVert v\rVert_2^2\longrightarrow\lVert f\rVert_2^2. \tag{EX.3}
```

The Rayleigh quotient is therefore of order (n), contradicting tightness
of ​\(\lVert X^{(n)}\rVert_{\rm op}/\sqrt n\).  To make the random-mixture
quantifier explicit, take a countable dense family of bounded simple
functions.  If \(K_w\ne0\) on an event of positive probability, one member
of that family and one rational \(c_0>0\) satisfy
\(|\langle f,K_wf\rangle|>c_0\) on an event of positive probability.
Equation (EX.3) then violates tightness on that fixed event.  Hence (w=0)
almost surely.
For a sign variable, conditional mean zero means conditional probabilities
​\(1/2,1/2\). Conditional edge independence then shows that every finite
set of edges consists of independent uniform signs. Thus the entire array is
i.i.d. Rademacher.

This proof uses only the representation theorem and does not assume that the
random graphon is deterministic or ergodic.

## 3. An elementary (0.5319\ldots) lower bound for i.i.d. signs

Expose an i.i.d. signing one new vertex at a time. Set (x_1=1), and after
​\(x_1,\ldots,x_{j-1}\) have been chosen put

```math
S_j=\sum_{i<j}X_{ij}x_i,
\qquad
x_j=\operatorname{sign}(S_j),                              \tag{EX.4}
```

with either sign at zero. The resulting energy is exactly

```math
H_{X^{(n)}}(x)=\sum_{j=2}^n|S_j|.                          \tag{EX.5}
```

Conditional on all previously exposed rows, (S_j) is a simple random walk
of length (j-1). Its conditional law does not depend on the past; inductively,
the variables ​\(|S_j|\) are independent. The central limit theorem and uniform
integrability give

```math
\mathbb E|S_j|
=\sqrt{\frac{2(j-1)}\pi}+o(\sqrt j).                       \tag{EX.6}
```

Consequently,

```math
\frac1{n^{3/2}}\sum_{j=2}^n\mathbb E|S_j|
\longrightarrow\frac23\sqrt{\frac2\pi}.                  \tag{EX.7}
```

Also ​\(\operatorname{Var}|S_j|\le j-1\) and ​\(|S_j|\le j-1\).
Bernstein's inequality makes the probability of a fixed
​\(\varepsilon n^{3/2}\) downward deviation at most
​\(\exp[-c_\varepsilon\sqrt n]\), which is summable. Equations
(EX.5)--(EX.7), Borel--Cantelli, and (Q\ge H(x)) prove (EX.1).

## 4. Consequence for action recovery

A purified liminf cluster can be chosen with

```math
\Phi(T)\le1+\eta
```

for every fixed ​\(\eta>0\), because ​\(\liminf M_n/n^{3/2}\le1/2\).
The i.i.d. exchangeable restrictions instead obey

```math
\liminf_n\Phi(T_{X^{(n)}})
\ge\frac43\sqrt{\frac2\pi}
=1.063846081\ldots .                                      \tag{EX.8}
```

Thus, for sufficiently small ​\(\eta\), they cannot satisfy directed
one-profile recovery with the quantitative error required by action
continuity.

The obstruction is exact at the architectural level:

- **ruled out:** choose one projectively consistent exchangeable infinite
  signing and use all of its finite restrictions;
- **not ruled out:** an order-dependent microcanonical ensemble, a globally
  conditioned exchangeable law at each separate order, or deterministic
  nonprojective realizers.

The latter mechanisms lose Kolmogorov consistency, which is precisely the
feature that would otherwise have solved the all-order requirement for free.
