# Odd transported variance: a sharp reciprocal-moment audit

Date: 2026-09-06. The director supplied the positive-tensor subset argument;
the audit agent independently verified it and strengthened it using the
range projector. The director independently reconstructed the strengthening.
No operator cap or entry-flatness is needed for the variance theorem itself.

## 1. General positive-matrix partition inequality

Let L be a finite-dimensional positive semidefinite matrix of rank r, and
let P_1,...,P_n be mutually orthogonal projections summing to the identity.
Write w_i=Tr(P_i L) and v_i=Tr(P_i L²)=||P_i L||_F². Then

```math
\boxed{\sum_{i:v_i>0}\frac{w_i^2}{v_i}\le r.} \tag{1}
```

If v_i=0, then w_i=0, and that summand is interpreted as zero.

Indeed let Pi be the orthogonal projector onto range(L). Since LPi=L,

```math
\langle P_i\Pi,P_iL\rangle_F
=\operatorname{Tr}(\Pi P_iL)
=\operatorname{Tr}(LP_i)=w_i.
```

Moreover ||P_i Pi||_F²=Tr(P_i Pi)=:d_i. Cauchy--Schwarz gives
w_i²<=d_i v_i. Summing w_i²/v_i<=d_i proves (1), because sum_i d_i=r.
All traces are finite, and no inverse of L is used.

## 2. Odd tensor specialization

Let B be a real symmetric n-by-n matrix with every row Euclidean norm one.
Consequently every column also has norm one. Put Q=B² and, for odd p>=1,

```math
T_p=BQ^{\circ p}B.
```

Let v_j be row j of B, set r_p=(p+1)/2, and form the positive matrix

```math
L_p=\sum_{j=1}^n
(v_j^{\otimes r_p})(v_j^{\otimes r_p})^T
```

on the r_p-fold tensor space. Its rank is at most n. Let P_i select first
tensor coordinate i and leave all other tensor coordinates unrestricted.
The unit row and column norms give

```math
\operatorname{Tr}(P_iL_p)
=\sum_j B_{ji}^2\|v_j\|^{2r_p-2}=1.
```

The squared row-block norm is exactly

```math
\begin{aligned}
\|P_iL_p\|_F^2
&=\sum_{jk}B_{ji}B_{ki}
 \langle v_j,v_k\rangle^{r_p-1+r_p}\\
&=\sum_{jk}B_{ij}B_{ik}Q_{jk}^{p}
=(T_p)_{ii}.
\end{aligned}
```

Thus every diagonal is strictly positive, and (1) gives the stronger
reciprocal bound

```math
\boxed{\sum_i\frac1{(T_p)_{ii}}
\le\operatorname{rank}(L_p)\le n.} \tag{2}
```

The director's subset inequality follows either by Cauchy--Schwarz from
(2), or directly by setting P=sum_{i in S}P_i and using
Tr(PL_pP)=|S| and rank(PL_pP)<=n:

```math
\sum_{i\in S}(T_p)_{ii}=\|PL_p\|_F^2
\ge\|PL_pP\|_F^2\ge |S|^2/n. \tag{3}
```

The oddness of p is used precisely in forming the positive integer tensor
power with 2r_p-1=p. It is not dispensable.

## 3. Mixtures and consequences

For any nonnegative finite or countable weights c_p on odd p, with
tau²=sum_p c_p in (0,infinity), put T=sum_p c_p T_p and sigma_i²=T_ii.
At each fixed matrix order the sum converges: |Q_ij|<=1, and its entries
are bounded by the summable weights times finite matrix-dependent bounds.
Scalar convexity of x->1/x gives

```math
\frac{\tau^2}{T_{ii}}
\le\sum_p\frac{c_p}{\tau^2}\frac1{(T_p)_{ii}}.
```

Nonnegative summation and (2) therefore imply

```math
\boxed{\frac1n\sum_i\frac{\tau^2}{\sigma_i^2}\le1.} \tag{4}
```

In particular,

```math
\#\{i:\sigma_i^2\le\epsilon\tau^2\}\le\epsilon n,
\qquad \frac1n\sum_i\sigma_i^{-1}\le\tau^{-1},
\qquad \frac1n\sum_i\sigma_i\ge\tau. \tag{5}
```

The second assertion is Cauchy--Schwarz; the last follows from Jensen for
x->x^{-2}. Each sigma_i² is also at least tau²/n. This order-dependent
floor is compatible with the actual apex counterexample tending to zero;
there is still no positive order-independent pointwise floor.

For independent standard normal N, (5) directly yields the averaged
small-ball estimate

```math
\frac1n\sum_i\mathbb P(|\sigma_iN|\le\eta)
\le\sqrt{2/\pi}\,\eta/\tau. \tag{6}
```

Equations (4)--(6) apply to every zero-first odd residual covariance in
the audited feedback theorem, including countable Hermite expansions.

## 4. Zero-first hard-sign passage under the feedback operator cap

Now additionally assume the actual-signing bounded-operator setting of
the audited zero-first masked covariance theorem. Let f be the fixed odd
zero-first input, tau²=E f(N)²>0, Y=Bf(BS), and H bounded even with the
stated marginal continuity. Set mu=E H(N). The hard response
C_i=H((BS)_i) sign(Y_i) is obtained by smooth odd approximations to sign.
The value assigned at zero does not affect the limit: the averaged
actual mass near zero tends to zero by the one-root comparison and (6).

For example use tanh(y/eta). Its Gaussian comparison L² error from sign
is O(eta/sigma_i); averaging gives O(eta/tau). Transfer the actual error
using fixed-width bounded continuous small-ball majorants, taking matrix
order first and eta second. Matrix Cauchy--Schwarz then transfers both
normalized covariance in nuclear norm and normalized energy. For the
linear-noise projection, write a_i^eta=E[N tanh(sigma_iN/eta)]/sigma_i.
Although 1/sigma_i is not uniformly bounded, the relevant error is

```math
\sigma_i^2\left(a_i^\eta-\sqrt{2/\pi}/\sigma_i\right)^2
\le\mathbb E\left[\tanh(\sigma_iN/\eta)
                  -\operatorname{sign}(N)\right]^2.
```

Thus the same averaged estimate controls the Gaussian weighted-vector
energy, with no inverse-variance operator bound. The final formula is

```math
\boxed{\frac{\mathbb E C^TBC}{2n}
=\frac{\mu^2}{\pi n}\operatorname{Tr}
 \left(BD_{\sigma^{-1}}TD_{\sigma^{-1}}\right)+o(1).} \tag{7}
```

The actual cross gain is

```math
\frac{\mathbb E f(BS)^TBC}{n}
=\mu\sqrt{2/\pi}\,\frac1n\sum_i\sigma_i+o(1). \tag{8}
```

For a nonnegative feasible mask, (5) bounds (8) below by
mu sqrt(2/pi) tau. The old zero-first self-energy is o(1), so the exact
endpoint identity adds the absolute value of (7). No sign of that
self-energy is assumed.

This passage is deliberately stated for zero-first feedback. For a
nonzero-first response, a raw-Z cross with unbounded inverse-variance
coefficients must keep its own cutoff argument; averaged unweighted
raw-covariance convergence alone does not justify removing those weights.

## 5. Even powers fail, including on actual bounded-op signings

Already B=[[1,-1],[-1,1]]/sqrt(2) has unit row and column norms,
Q=[[1,-1],[-1,1]], Q^{circ2}=J, and BQ^{circ2}B=0. Thus neither (2)
nor a positive variance conclusion holds for arbitrary even p.

There is also a scalable ACTUAL hollow-signing counterexample. Let H_m
be a symmetric Hadamard matrix, let R=[[1,-1],[-1,1]], and set

```math
F=R\otimes H_m,\quad D=\operatorname{diag}(F),\quad
A=F-D,\quad d=2m-1,\quad B=A/\sqrt d.
```

Then A is symmetric hollow with every off-diagonal entry a sign, and
||B||<=(2sqrt(m)+1)/sqrt(2m-1)->sqrt(2). Writing Q=B² and
c=(d-1)/d, direct entrywise squaring gives

```math
Q^{\circ2}=J_2\otimes W+(1-c^2)I,
\qquad W_{aa}=c^2.
```

Since R J_2=0, F annihilates J_2 tensor W on both sides. Consequently
A(J_2 tensor W)A=D(J_2 tensor W)D, whose diagonal is c². Every root
therefore has the exact even-transport variance

```math
\boxed{[BQ^{\circ2}B]_{ii}
=\frac{c^2}{d}+1-c^2
=1-\left(\frac{d-1}{d}\right)^3
=\frac{3d^2-3d+1}{d^3}\longrightarrow0.} \tag{9}
```

Oddness is essential even under all the original signing and fixed
operator-cap hypotheses. The script
`computations/continued_audit_tensor_variance_leverage_2026_09_06.py`
replays (9) with integer arithmetic and checks the odd reciprocal bounds
on random, Steiner, and apex/opposite-twin signings. Numerical odd checks
are diagnostics; the positive-tensor proof is the certificate.

## 6. Independent check of the nonzero-first cutoff theorem

The feedback agent's
`continued_feedback_threshold_without_variance_floor_2026_09_06.md`
also passes independently. Its zero-first formulas agree with (7)--(8).
For nonzero first coefficient b, the coherent hard response is literally
`c0_i=H(G_i)[2 Phi(b(QS)_i/sigma_i)-1]`, retaining the Boolean QS law.
The derivative coefficient is

```math
a_i=\sqrt{2/\pi}\,\sigma_i^{-1}
\mathbb E_S[H(G_i)e^{-b^2(QS)_i^2/(2\sigma_i^2)}].
```

For each fixed epsilon, set a_i^epsilon=0 on sigma_i²<=epsilon tau².
The source's energy formula with this cutoff and error O(sqrt(epsilon))
is correctly ordered. For fixed smoothing width eta, a_i^eta is uniformly
bounded, so the raw-noise normalized-nuclear covariance approximation
transfers its weighted squared norm after the matrix limit. Gaussian
integration by parts gives |a_i^eta|sigma_i<=sqrt(2/pi)||H||infinity.
Consequently the discarded rows cost O(epsilon) normalized squared norm,
and O(sqrt(epsilon)) normalized cross/energy by the fixed operator cap.

On retained rows the coefficient convergence is uniform. Indeed the
difference is bounded by C eta/sigma_i², by integrating the Gaussian
IBP expression on the interval where the smoothed sign differs. Hence
it is at most C eta/(epsilon tau²). This validates passage to a^epsilon
in the raw-Z cross after fixing epsilon. The Gaussian trace can then
remove its cutoff, because its diagonal a_i²sigma_i² stays bounded.
The raw cross is correctly left in the ordered cutoff form: no unjustified
unbounded multiplication of the covariance error occurs.

The local absolute-value gain is continuous and uniformly integrable,
so it needs no such inverse-variance cutoff. This checks the source's
normalizations and every limit order in its formulas (7)--(9), but does
not enlarge the source class beyond scalar or fixed independent colors.

## 7. Anchored mixed-correlation extension and its parity boundary

The director's mixed-correlation generalization also passes independently.
Let P be ANY real correlation matrix and retain Q=B². Then

```math
T=B[Q\circ P^{\circ2}]B
\quad\Longrightarrow\quad
\sum_iT_{ii}^{-1}\le n. \tag{10}
```

Choose unit Gram vectors p_j for P and use w_j=b_j tensor p_j in the
positive matrix L=sum_j w_j w_j^T. Its rank is at most n. Projection on
first B-coordinate i gives Tr(P_i L)=sum_j B_ji²=1. Its squared block
norm is exactly

```math
\sum_{jk}B_{ji}B_{ki}
 \langle p_j,p_k\rangle
 \langle b_j,b_k\rangle\langle p_j,p_k\rangle
=[B(Q\circ P^{\circ2})B]_{ii}.
```

The general inequality (1) proves (10). Nonnegative mixtures inherit the
same reciprocal bound after rescaling by their total coefficient mass.
In particular, odd p and even q give a correlation matrix
`P_a=Q^{circ((p-1)/2)} circ P^{circ(q/2)}` by the Schur product theorem,
where exponent zero means J. Thus mixtures of
`Q^{circ p} circ P^{circ q}` have the bound. This covers an odd anchored
Gaussian field and even dependence on an independent field of a different
correlation. The anchored Q factor and the squared correlation are real
requirements, not merely total odd parity.

Here is a bounded-response falsifier for arbitrary globally odd dependence
on three different covariance channels. Use the actual opposite-twin
signings of Section 5, and let three independent Gaussian fields have
covariances Q,Q,P respectively, where P=J_2 tensor I_m. All three
correlation operators are uniformly bounded. Set
`F_i=sin(G_i^(1)) sin(G_i^(2)) sin(U_i)`; this is bounded and globally odd.
Its exact covariance is

```math
R=e^{-3}(\sinh^{\circ}Q)^{\circ2}\circ\sinh^{\circ}P
=\beta(J_2\otimes I_m)+(\tau^2-\beta)I,
```

where `tau²=(e^-1 sinh(1))³` and
`beta=e^-3 sinh(c)² sinh(1)`, c=(d-1)/d. The same exact annihilation
calculation as Section 5 gives, for every i,

```math
[BRB]_{ii}=\beta/d+\tau^2-\beta\longrightarrow0.
```

Thus total oddness with three independent distinct correlation channels
does not imply (10), even for a bounded response and an actual hollow
bounded-op signing B. The third field is a literal Gaussian linear field
with covariance P; no claim that it is another copy of BS is made. This
does not contradict the equal-correlation independent-color theorem.
