# Actual full signings: high energy alone does not imply cheap physical response

2026-09-17. Bernoulli-track discriminator for the adaptive response
theorem's operator-norm hypothesis. This is an ACTUAL full-sign
construction with bounded normalized cap, not a weighted surrogate.
Its hard code is explicitly a macroscopic distance BELOW the absolute
ground cap. It therefore does NOT answer the exact-ground dual question.

## 1. The statement

There are hollow symmetric full signings A_N, at every sufficiently
large EVEN order N, such that

```
Q(A_N)=O(N^(3/2)),
for every fixed 0<c<1/2,
 inf_(ALL physical sign laws nu)
 max_(|H_A(x)|>=c N^(3/2)) E_nu |h dot x|/sqrt(N)
       ->kappa=sqrt(2/pi).                         (1)
```

The same limit holds if the infimum is restricted to centered isotropic
sign laws, or further to uniformly subGaussian ones: the lower bound
applies to ALL laws, and independent fair signs give the matching upper.

In fact a subset C_N witnessing the lower bound has all its energies
within N^(11/8) of two values +-q_0, where q_0~N^(3/2)/2, while

```
Q(A_N)>=q_0+(kappa/4-o(1))N^(3/2).                 (2)
```

Thus every word of this witnessing code has absolute deficit at least
(kappa/8)N^(3/2), eventually. Exact or vanishing-deficit ground support
is NOT being smuggled into (1).

## 2. Genuine full-sign rounding with a controlled cap

Take even m~N^(3/4), l=N-m even, and

```
a=m(m-2)/[l(l-2)],
q_0=[m(m-1)+a l]/2=[m+a l(l-1)]/2.
```

For large N, 0<a<1. Use a positive clique on the m core vertices.
Choose independently every large-block edge B_ij in {+-1} with mean
-a, and independently every cross edge C_ij as a fair sign. All edges
are physical unit signs; there are no zero or weighted edges and no
constant term in the actual energy.

The EXPECTED energy matrix is precisely the two-block weighted model
from [the preceding scope test](paper_bernoulli_centered_energy_scope_2026_09_17.md),
whose cap is q_0. For any fixed physical word z, the centered energy
Z_z=H_A(z)-E H_A(z) is a sum of independent variables of interval
length two. Its subGaussian variance proxy is at most

```
V=binom(l,2)+ml<=N(N-1)/2.
```

Thus, for all u>0,

```
P(|Z_z|>u)<=2exp[-u^2/(2V)]<=2exp[-u^2/N^2].        (3)
```

A union over the ENTIRE 2^N cube yields, for every fixed D>sqrt(log2),
with probability tending to one,

```
Q(A)<=q_0+D N^(3/2).                               (4)
```

For example D=1 gives the explicit cap constant 3/2+o(1). The latent
mean matrix is used only for probabilistic analysis, never substituted
for the resulting full signing.

## 3. Almost all of the balanced slice remains at the prescribed energy

Let B_l be the balanced l-coordinate slice. Define G subset B_l to
contain y exactly when both words (+1_m,y) and (-1_m,y) satisfy

```
|H_A(+-1_m,y)-q_0|<=u_N,   u_N=N^(11/8).
```

The set G is antipodal. For each fixed y, (3) and a two-event union
bound give failure probability at most 4exp(-N^(3/4)). Therefore
Markov's inequality gives, with probability tending to one,

```
|B_l\G|/|B_l|<=d_N:=exp[-N^(3/4)/2].               (5)
```

Fix ANY balanced word b on the m core coordinates. The four words
(+-b,+-1_l) each have expected energy -q_0. Another finite union of
(3) shows that all four satisfy

```
|H_A(+-b,+-1_l)+q_0|<=u_N                         (6)
```

with probability tending to one. On the common event, put

```
C_N=({+-1_m} times G) union ({+-b} times {+-1_l}).
```

Every word of C_N has |H_A|>=q_0-u_N, hence belongs to the high-energy
set in (1) for every fixed c<1/2 eventually.

## 4. The response game is robust under the exponentially small deletion

Let f_l be the absolute overlap mean of two independent uniform
balanced l-words. The exact preceding scope proof gives

```
f_l=kappa sqrt(l)+O(l^(-1/2)),
E_(Y uniform B_l)|h_l dot Y|
       >=f_l[1-|sum h_l|/l].                       (7)
```

Since an overlap is at most l, (5) implies, uniformly in EVERY physical
h_l,

```
E_(Y uniform G)|h_l dot Y|
 >=f_l[1-|sum h_l|/l]-d_N l/(1-d_N).               (8)
```

This is a deterministic consequence of missing fraction, so no union
over h is hidden. For an arbitrary law nu let b_nu=E_nu|sum h_l|.
Averaging the four queries in the second sector of C_N gives a lower
bound b_nu: averaging the two core signs removes any possible favorable
cancellation with their overlap. Averaging the first sector similarly
gives at least the mean in (8), since |u+v|+|-u+v|>=2|v|. Consequently

```
max_(z in C_N) E_nu|h dot z|
 >=max{b_nu, f_l(1-b_nu/l)-d_N l/(1-d_N)}
 >=l f_l/(l+f_l)-d_N l/(1-d_N)
 =(kappa-o(1))sqrt(N).                             (9)
```

Here l/N->1. Independent fair signs, on the other hand, have response
E|S_N|=(kappa+o(1))sqrt(N) at EVERY Boolean word, proving (1).

An explicit near-optimal dual law for (9) mixes the first query sector
with the second with weights 1-alpha and alpha, where
alpha=f_l/(l+f_l). Its covariance norm is at least (1-alpha)m, tested
on the normalized all-ones direction of the core. This diverges while
the coherent core occupies only o(N) coordinates. Most of its marked
energy is carried by that core. The bounded-query-covariance lemma is
therefore not contradicted.

## 5. Actual full-sign bulk energy makes the hard code nonground

This part is included to quantify the exact-extremality gap, not merely
leave it unverified. Pair the l large-block vertices arbitrarily and
assign opposite spins within each pair. There are l/2 pair spins. For
two distinct pairs i,j the effective coefficient is

```
D_ij=B_(i1,j1)-B_(i1,j2)-B_(i2,j1)+B_(i2,j2).
```

These coefficients are independent over unordered pair pairs, centered,
bounded by four, and have variance 4(1-a^2). Split the pairs into I,J,
each of size l/4+O(1). Give all I spins +1, and set each J spin to
the sign of its field F_j=sum_(i in I)D_ij. The fields are independent
over j, and

```
E|F_j|=(kappa+o(1))sqrt(l),
sum_(j in J)|F_j|=(kappa/4+o_P(1))l^(3/2).          (10)
```

For the first assertion, the elementary bounded-summand Lindeberg
theorem applies since a->0 and the normalized largest summand tends
to zero; second moments give uniform integrability. Alternatively its
characteristic-function proof uses only the fourth-order expansion
of the finite product. For the second assertion, the variance of the
sum is at most O(l^2), so Chebyshev suffices.

Internal I-pair and J-pair interaction coefficients are independent of
all the cross fields used for selection. Their respective energies
have conditional mean zero and variance O(l^2), including after the
selected J signs are fixed. The original within-pair energies have
absolute total at most l/2. Hence the balanced word y constructed in
this way satisfies

```
H_B(y)>=(kappa/4-o_P(1))l^(3/2).                   (11)
```

The original core/large-block bridge is independent of B and y. With
the core aligned, its conditional variance is ml=o(N^2), and therefore
its value is o_P(N^(3/2)). Since the clique energy is q_0-a l/2 and
a l/2=o(N^(3/2)), (11) proves (2) with probability tending to one.
All events (4)--(6) and (11) thus occur together with probability
tending to one, proving deterministic existence of A_N and C_N.

## 6. Precise implication

Removing the actual-A operator-norm hypothesis from the uniform
MACROSCOPIC-HIGH-ENERGY response theorem would make it false: bounded
cap, full unit coefficients, and such high-energy membership are not
enough. Nor is a large signed average energy by itself enough to
control an arbitrary dual query covariance.

This construction does NOT refute a theorem for exact minimizing
signings, exact absolute ground codes, or vanishing absolute deficits.
It explicitly demonstrates that residual full-sign bulk variation can
produce a higher cap outside an apparently useful code. That is the
escape information absent from the marked covariance identity alone.

## 7. Replay and verification scope

`computations/paper_bernoulli_2026_09_17_actual_high_energy_obstruction.py`
passes 4,800 exact rational deleted-slice response inequalities and
exhaustively computes the actual full parent cap in 12 finite cases.
It verifies the paired-spin energy decomposition exactly and separately
records large-order greedy diagnostics, explicitly labeled as Monte
Carlo rather than asymptotic proof. The finite tolerance used for those
small deleted-slice tests is not promoted to the asymptotic probability
estimate; that estimate is established by (3)--(5).

Output: `tmp/paper_portfolio_2026_09_17/bernoulli/actual_high_energy_obstruction.json`.
Python compilation also passes.

The localization track independently reconstructed the full proof,
including uniform deletion control before averaging any physical law,
the paired-spin conditional independence, and the macroscopic gap
between the witness code and the actual cap: PASS. In particular this
audit confirms the explicitly nonground scope.
