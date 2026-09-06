# Inhomogeneous Gaussian responses retain a nonlinear return

2026-09-06. Director proof of the trace-rank mechanism proposed by the seed
researcher; independently audited by the adversarial researcher. This is a
finite Gaussian-kernel theorem. Applying it to actual
Boolean feedback still requires the separate mixed-covariance comparison;
that comparison is not assumed proved in this file.

## 1. Statement

Let T be an n-by-n positive semidefinite matrix with `q_i=T_ii`, and
let `g_i:R -> [-1,1]` be bounded odd nondecreasing functions. Coordinates
with `q_i=0` are ignored. Expand in orthonormal probabilists' Hermite
polynomials, `g_i(t)=sum_(k odd) alpha_ik h_k(t)` on the retained coordinates;
set every `alpha_ik=0` on ignored zero-q rows. Define

```math
\bar g(t)=\frac1n\sum_i\sqrt{q_i}\,g_i(t),\qquad
a=\langle\bar g,h_1\rangle,
M=\frac1n\sum_i\sqrt{q_i}.
```

Assume `a>=c>0` and `M<=M0<infinity`. For odd k set
`D_k=diag(alpha_ik/q_i^(k/2))`, with zero entries when `q_i=0`.
Then

```math
\frac1n\sum_{k\ge3\ {\rm odd}}
 \operatorname{tr}(T D_k T^{\circ k}D_k)
\ge \|\bar g-a h_1\|_{L^2(\gamma)}^2
\ge \delta(c,M0)>0,                                 (1)
\delta(c,M0)=\mathbb E(c|G|-M0)_+^2,
\quad G\sim N(0,1).
```

There is a finite `K0=K0(c,M0)`, independent of n, and an odd
`3<=k<=K0` such that

```math
\frac1n\operatorname{tr}(T D_kT^{\circ k}D_k)
\ge\delta(c,M0)/(2K0).                              (2)
```

If another PSD matrix Q satisfies `Q>=T/C`, then T can be replaced by Q
on the left of (2), at the additional factor `1/C`. Also

```math
\operatorname{tr}(Q D_kT^{\circ k}D_k)
\le\|Q\|_{op}\sum_i\alpha_{ik}^2
\le\|Q\|_{op} n.                                   (3)
```

No common sign of the rowwise higher coefficients is needed. This is the
point: a selected-row argument would change the response being measured,
whereas (1)--(3) retain every response row.

## 2. Rank, rather than entrywise positivity

Represent `T_ij=<v_i,v_j>`. For `k=2r-1` let

```math
S_k=\sum_i (D_k)_{ii}
           (v_i^{\otimes r})(v_i^{\otimes r})^T.
```

This self-adjoint operator may be indefinite, but has rank at most n.
Consequently

```math
\operatorname{tr}(T D_kT^{\circ k}D_k)
=\sum_{i,j}(D_k)_{ii}(D_k)_{jj}T_{ij}^{2r}
=\|S_k\|_F^2
\ge\frac{(\operatorname{tr}S_k)^2}{n}
=n\left(\frac1n\sum_i\sqrt{q_i}\alpha_{ik}\right)^2.
```

The last quantity is n times the squared kth Hermite coefficient of
`bar g`. Sum over odd k at least three. Boundedness and oddness give
`|bar g(t)|<=M0`; pointwise
`|bar g(t)-at|>=(c|t|-M0)_+`. This proves (1). Monotonicity is not
needed for this infinite-tail assertion.

All sums are legitimate in finite dimension. Each Schur matrix
`D_k T^{circ k}D_k` is PSD and has trace `sum_i alpha_ik^2`.
The total trace is finite by Parseval; multiplication by fixed T or Q
therefore permits summation in trace norm. Singular q causes no problem
because the corresponding Gram row is zero.

## 3. A uniform finite Hermite degree

A bounded nondecreasing function with range inside `[-M0,M0]` is a
positive mixture of threshold functions with total variation at most `2M0`.
For standard Gaussian X,Y with correlation `0<=rho<1`, Jensen gives

```math
\mathbb E[\bar g(X)-\bar g(Y)]^2
\le\frac{4M0^2}{\pi}\arccos\rho.                     (4)
```

Here threshold-disagreement probability is maximal at threshold zero,
where it equals `arccos(rho)/pi`. To check maximality, differentiate
`2[Phi(u)-Phi_2(u,u;rho)]`: its derivative is
`2 phi(u)[1-2Phi(u sqrt((1-rho)/(1+rho)))]`, positive below zero and
negative above zero. The threshold integral follows by monotone limits
when the response has jumps or does not attain its endpoint bounds.

If `a_k` are the Hermite coefficients of `bar g`, (4) and Mehler's
identity show, for integer `K>=2` and `rho=1-1/K`,

```math
\sum_{k>K}a_k^2
\le\frac{2M0^2\arccos(1-1/K)}
          {\pi[1-(1-1/K)^{K+1}]}
\le\frac{3M0^2}{\sqrt K}.                           (5)
```

For the last bound use `arccos(1-1/K)<=pi/sqrt(2K)` and
`1-(1-1/K)^(K+1)>1-e^(-1)>1/2`. Taking
`K0>=max(2,36 M0^4/delta(c,M0)^2)` makes (5) at most delta/2.
Equation (1), followed by pigeonhole over at most K0 coefficients,
proves (2). The selected k may depend on the finite input, but belongs
to this fixed finite set; all-order comparison errors must be uniform
only over that set, not over growing degree.

For explicit evaluation,
`delta(c,M0)=2[(c^2+M0^2)barPhi(M0/c)-c M0 phi(M0/c)]`.
Its positivity also follows directly from its defining integral.

## 4. Exact prospective implication for actual feedback

Suppose B is symmetric and an actual feedback vector C and Gaussian-style probes U_k have
the already-normalized mixed identity

```math
\mathbb E[(BC)^T B D_k U_k]
=\operatorname{tr}(Q D_k T^{\circ k}D_k)+o(n),
\quad Q=B^2\succeq T/C_0,
```

and the matching second-moment bound (3), uniformly for the finite degrees
up to K0. Then one retained actual return correlation is at least
`n delta/(2 C0 K0)+o(n)`, with probe second moment at most
`||Q||op n+o(n)`. This is the explicit quantitative payoff of the lemma.
Turning that correlation into a spin-flip gain additionally requires the
probe's joint local Gaussian law relative to the LITERAL old coherent
queries. None of those Boolean comparison statements follows merely from
Gaussian marginal convergence or from the theorem above.

Even if these comparisons are proved for every fixed operator cap L,
their positive constants must still be compared with the O(1/L) spectral
deletion cost before asserting a new unrestricted signing bound.
