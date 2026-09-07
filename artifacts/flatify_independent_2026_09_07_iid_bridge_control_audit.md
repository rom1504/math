# Independent audit: a dense iid bridge already exceeds the child target

Date: 2026-09-07. Verdict: PASS for the adversary's stopped-Brownian Parisi
control and iid-bridge obstruction. This is a construction-class obstruction,
not a lower bound on arbitrary sign bridges or on the original minimum.

## 1. Primary theorem and normalization checked

[Auffinger--Chen, Theorem 1 and Corollary 2](https://arxiv.org/html/1606.05335v2)
identify the SK ground-state constant P for covariance N q^2/2 with

```math
P=inf_gamma [Psi_gamma(0,0)-(1/2)int_0^1 s gamma(s) ds],
Psi_gamma(0,0)=sup_{|u|<=1, adapted}
 [E|B_1+int_0^1 gamma(s)u_s ds|
  -(1/2)int_0^1 gamma(s) E u_s^2 ds].
```

The source permits nonnegative nondecreasing integrable gamma, extending its
step-function control formula by L1 continuity. Its SK convention is exactly
xi(s)=s^2/2. This includes an irrelevant common Gaussian relative to the
usual hollow i<j Hamiltonian.

## 2. A genuinely admissible common control cancels every gamma

Let W be Brownian motion, tau its first exit from (-1,1), and
p(t)=P(tau>t). Set

```math
a(t)=E min(t,tau)=int_0^t p(r)dr.
```

Optional stopping gives E tau=1, so a maps [0,infinity) onto [0,1).
It is strictly increasing, with derivative p(t)>0. Put t(s)=a^(-1)(s),

```math
B_s=int_0^{t(s)} sqrt(p(r))dW_r,
u_s=W_{t(s) wedge tau},  0<=s<1,
u_1=sigma=W_tau in {+-1}.
```

B has quadratic variation s and is a standard Brownian motion. The
deterministic positive integrand and time change preserve the filtration
before time 1, so u is adapted to B, not to an illicit larger filtration.
It is a bounded continuous martingale, with

```math
E u_s^2=E min(t(s),tau)=s,
du_s=1_{tau>t(s)}p(t(s))^(-1/2)dB_s.
```

The coefficient need not be bounded near 1, but its expected squared
integral is 1. Thus the stochastic integrals and endpoint L2 limits are
legitimate. Martingale conditioning gives E[sigma u_s]=E u_s^2=s.
For every admissible gamma, use |z|>=sigma z in its control objective.
The two gamma contributions and the external Parisi penalty cancel exactly,
yielding

```math
P>=E[sigma B_1].
```

The latter covariance is the expected quadratic covariation of u and B:

```math
E[sigma B_1]
 =int_0^1 sqrt(p(t(s)))ds
 =int_0^infinity p(t)^(3/2)dt.                       (1)
```

This is a lower bound on EVERY Parisi trial functional, not an upper bound
obtained by choosing one gamma. That direction is essential.

## 3. Explicit survival estimate and exact constant

The killed heat equation on (-1,1), with initial survival value 1, gives

```math
p(t)=(4/pi)sum_{j>=0} (-1)^j/(2j+1)
 exp[-(2j+1)^2 pi^2 t/8].
```

The alternating terms decrease in magnitude for t>=0; therefore, with
c=pi^2/8,

```math
p(t)>=(4/pi)e^(-ct)[1-e^(-8ct)/3].
```

Convexity gives (1-z)^(3/2)>=1-(3/2)z for 0<=z<=1/3. Substituting into
(1) and integrating the two exponentials yields

```math
P >= (4/pi)^(3/2) c^(-1)[2/3-1/19]
   =2240/(57 pi^(7/2))
   =.7150700826077228... .                           (2)
```

In particular P>1/sqrt(2). This inequality can be certified without decimal
arithmetic using pi<22/7 and the integer inequality
2*2240^2*7^7 > 57^2*22^7.

## 4. Gaussian comparison to a square bridge

For an m by m iid standard Gaussian matrix J, define

```math
X_(x,y)=m^(-1/2)sum_ij J_ij x_i y_j.
```

Its variance is m and its covariance is m q_x q_y. On the same pair-index
set take a 2m-spin SK process Y with covariance

```math
Cov(Y_(x,y),Y_(x',y'))=m(q_x+q_y)^2/4.
```

Both have variance m, and Cov(Y)-Cov(X)=m(q_x-q_y)^2/4>=0. Gaussian
comparison (or its direct Gaussian integration-by-parts softmax proof)
therefore gives E max X>=E max Y. Consequently

```math
liminf_m E ||J||_(infinity->1)/m^(3/2) >= 2P.        (3)
```

The absolute value is automatic for a bilinear form by reversing x.
For the conventional hollow SK process, add one common N(0,1/2) variable
to obtain the displayed covariance. This does not change its expected
maximum, and avoids a missing-diagonal normalization mistake.

## 5. Bernoulli transfer and probability strength

For independent symmetric Bernoulli entries replace J one coordinate at a
time in

```math
f_beta(J)=(beta m)^(-1)log sum_(x,y)
 exp[(beta/sqrt(m))sum_ij J_ij x_i y_j].
```

Each third derivative has magnitude <=8 beta^2/m^(5/2), by its bounded
third centered spin cumulant. First and second moments match, and the sum
of third-moment Taylor remainders is O(beta^2/sqrt(m)). Softmax contributes
at most 2log2/beta. With beta=m^(1/6), the normalized expectation transfer
error is O(m^(-1/6)); the separate SK limiting error is merely o(1).

Changing one Bernoulli entry changes the raw bilinear cap by at most 2.
Bounded differences therefore imply

```math
P{||J||_(infinity->1) <= E||J||_(infinity->1)-epsilon m^(3/2)}
 <= exp(-epsilon^2 m/2).
```

Together with (2)--(3), with probability tending to 1 (indeed exponentially
for any fixed slack after a sufficiently large order),

```math
||J||_(infinity->1)/(2m)^(3/2)
 >= 2240/(57 sqrt(2) pi^(7/2))-o(1)
 = .5056309044355455...-o(1).                        (4)
```

## 6. Arbitrary diagonal children cannot cancel the obstruction

Complete J by ANY hollow symmetric diagonal blocks D_1,D_2, even chosen
after seeing J. For fixed x,y the full energy is a+b, with
a=H_D1(x)+H_D2(y), b=x^T J y. Reversing the entire first fibre changes
this to a-b while leaving a unchanged. Hence the full cap is at least
max(|a+b|,|a-b|)>=|b|. Optimizing x,y proves pointwise

```math
Q([[D_1,J],[J^T,D_2]])>=||J||_(infinity->1).
```

Thus the dense iid balanced bridge exceeds .50563 after parent normalization,
regardless of its child completion. This is strictly above the proved
selected-child target .493608094. Independent unbiased cross edges therefore
cannot provide the desired favorable balanced flatification, even if their
within-child signs are chosen adaptively.
