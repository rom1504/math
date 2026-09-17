# Complete-sign noise stability from amplified preparation and mixed-tail chaining

2026-09-17. Director proof combining the actual cloned-block regularizer
with Dirksen's mixed-tail chaining mechanism. **Verified: all three
researchers independently reconstructed the proof and the primary
mixed-tail theorem.** Unlike a Gaussian neighborhood alone, the perturbed matrices
here are EXACT hollow full signings. No convergence or new cap constant
is claimed; external novelty is unestablished.

## 1. The theorem and a concrete parameter instance

Fix 0<gamma<2/3. Given ANY full signing A_N, prepare W_N by the
[cloned-block theorem](paper_director_cloned_block_regularization_2026_09_17.md).
There are universal positive constants c,C such that, for sufficiently
large N (depending on gamma), put

```math
 T=C N^{3/2-\gamma},\qquad
 \rho=c\,N^{-3\gamma/2}/\log N,
 \mathcal C=\{x:Q(W)-|H_W(x)|\le T\}.
```

Flip EVERY edge of W independently with probability rho, obtaining W'.
Then

```math
 Q(W)\le Q(A)+C N^{3/2-\gamma/2},\qquad
 \log|\mathcal C|\le C_\gamma N^{1-\gamma/4}\log N,
 \Pr\{\text{some absolute maximizer of }W'
                              \text{ lies outside }\mathcal C\}
       \le C\exp[-c_\gamma N^{1-\gamma/2}\log N].       \tag{1}
```

The set C is determined BEFORE the fresh edge flips. Every tie and both
energy polarities are covered. For gamma=1/2, preparation costs O(N^(5/4)),
the witness window is O(N), its log size is O(N^(7/8)log N), and an
independent flip rate c N^(-3/4)/log N is allowed. The expected number
of changed edges is Theta(N^(5/4)/log N), while witness failure is at
most C exp(-c N^(3/4)log N).

## 2. Primary-source reconstruction: mixed tails retain TWO metrics

The primary source is S. Dirksen,
[*Tail bounds via generic chaining*](https://arxiv.org/abs/1309.3522),
Theorem3.5, equations(12)--(15), and Remark3.3(i), arXiv v2 pages4,6--8.
The director read the statement and full proof of this theorem. The PDF
is preserved in the dated research archive. The published paper is
Electronic Journal of Probability20(2015), paper53.

For a finite centered process X_t with increments satisfying

```math
 \Pr\{|X_t-X_s|>\sqrt u\,d_2(s,t)+u d_1(s,t)\}
                                  \le2e^{-u},
```

the theorem gives, for a fixed t0 and u>=1,

```math
 \sup_t|X_t-X_{t0}|
 \le C[\gamma_2(T,d_2)+\gamma_1(T,d_1)
                +\sqrt u\,\operatorname{diam}_{d_2}T
                +u\,\operatorname{diam}_{d_1}T]       \tag{2}
```

except on an event of probability e^(-u). Constants are universal.
The relevant proof mechanism is as follows. Choose admissible refining
partitions approximating each gamma functional separately; intersect
their level j-1 cells to obtain at most2^(2^j) cells. Representatives
telescope along this common tree. Above level floor(log2 u), a union
bound assigns increment budgets proportional to2^(j/2)d2+2^j d1;
their sums are bounded by the two gamma functionals. The lower levels
are controlled by the diameter budgets sqrt(u)d2+u d1. Integrating
the same tail, or truncating the tree at floor(log2 p), gives the stated
moment estimate. This argument requires NO independence among increments.
The arXiv proof's displayed good-event inequality on page7 has its sign
reversed; the subsequent argument uses the correct upper-bound event.
The elementary union bound and telescope verify the intended direction.

Two further inputs are explicit. First, for a finite set with d1 diameter
D, a partition that stays trivial until2^(2^j)>=|T| and then becomes
discrete proves gamma1(T,d1)<=C D log(2|T|). Second, Talagrand's Gaussian
majorizing-measure theorem, recalled in the primary paper Remark3.3(i),
identifies gamma2 of the canonical Gaussian metric up to universal
constants with its expected anchored absolute supremum. That classical
theorem is IMPORTED, not claimed to be re-proved here. No exact constant1
Gaussian comparison is assumed in(2).

## 3. The actual edge-flip process has the required metrics

Let v_(sigma,x)=(sigma x_i x_j)_(i<j), and choose an old signed ground
state i0=(sigma0,x0), so sigma0 H_W(x0)=Q(W). Define

```math
 d_i=Q(W)-\sigma H_W(x)\ge0,\qquad
 Z_e=-2W_e(B_e-\rho),\quad B_e\sim\operatorname{Bernoulli}(\rho),
 X_i=\langle Z,v_i-v_{i0}\rangle.
```

Then W'=(1-2rho)W+Z EXACTLY entrywise; W' has entries plus or minus1.
The independent centered Z_e have variance4rho(1-rho) and magnitude<=2.
Scalar Bernstein, proved by the usual bounded-variable exponential
series, implies mixed increments with

```math
 d_2(i,j)=C\sqrt\rho\,\|v_i-v_j\|_2,\qquad
 d_1(i,j)=C\|v_i-v_j\|_\infty.                       \tag{3}
```

The diameters are at most C sqrt(rho)N and C, respectively. Include the
zero anchor i0 in every set under consideration. Applying(2), the simple
gamma1 bound, and the Gaussian majorizing-measure comparison gives

```math
 \sup_{d_i\le D}|X_i|
 \le C[\sqrt\rho\,\mathcal W_{edge}(D)
       +\log(2|E_W(D)|)+\sqrt\rho N\sqrt u+u],         \tag{4}
```

with probability at least1-e^(-u). Here W_edge is the expected ANCHORED
absolute supremum of independent Gaussian edge disorder. The signed
level has at most2|E_W(D)| elements. Reversal duplicates do not matter.

## 4. Whole-energy profiles of ONE prepared actual signing

Write p=N^(gamma/2), r=N^(1/2-3gamma/4), with the integer floors of
the construction, B0=N/sqrt(p), and ell=log(ep). Its SAME realization
satisfies at every D>=0

```math
 b(E_W(D))\le C(D/r+B_0),\qquad
 \log(2|E_W(D)|)\le C\ell(D/r+B_0).                 \tag{5}
```

The first bound adds the pr new physical coordinates to the core width;
pr<=C B0 for these parameters. For the second, the core VC bound is
B=min(m,D/r+C B0); since B>=c B0 before saturation,
log(em/B)<=C ell. Above saturation the bound only weakens. Physical
magnetization multiplicities add at most pr log2 and are already absorbed.

Clipping independent Gaussian coordinates at sqrt(2 log(ep)), followed
by the Bernoulli contraction inequality, proves

```math
 w_G(E_W(D))\le C\sqrt\ell(D/r+B_0).
```

Indeed the clipped part costs at most its clipping threshold times
Bernoulli width, and the expected total tail is at most C N/p<=C B0.
Within one polarity the quadratic edge metric obeys
||v_x-v_y||2<=sqrt(N)||x-y||2. Gaussian increment comparison, plus the
two polarity union and the anchor, therefore gives

```math
 \mathcal W_{edge}(D)
       \le C\sqrt{N\ell}(D/r+B_0)+CN.                \tag{6}
```

All costs in(5)--(6) are for the FULL nearcode, not only a preferred
branch, an average query, or a compressed magnetization list.

## 5. Absorb both slopes, then exclude EVERY outside energy shell

Combining(4)--(6), its deterministic complexity term is at most

```math
 C[\sqrt{\rho N\ell}+\ell](D/r+B_0)+C\sqrt\rho N.    \tag{7}
```

For fixed gamma<2/3, r/ell tends to infinity. Choose the absolute
constant in rho small enough that sqrt(rho N ell)<=c r, and take
N large enough that C ell/r is sufficiently small. Thus the coefficient
of D in(7) is, say, at most1/32. Choose T=C0 rB0 with sufficiently
large fixed C0 to absorb the intercept. At the upper endpoint
D=2^(j+1)T of shell[2^jT,2^(j+1)T], (7) is at most2^(j-3)T.

An outside state can match or beat the old ground-state anchor under
W' only if

```math
 X_i\ge(1-2\rho)d_i\ge d_i/2.                       \tag{8}
```

Use(4) on each larger level, with
u_j=c min((2^jT)^2/(rho N^2),2^jT). Its remaining diameter budgets
are at most2^(j-3)T for sufficiently small c. The dyadic union then
has failure at most

```math
 C\exp\left[-c\min\left\{T^2/(\rho N^2),T\right\}\right].       \tag{9}
```

This includes all ties by making the displayed tail threshold strictly
smaller than d_i/2. With the prescribed parameters, T~N^(3/2-gamma)
and T^2/(rho N^2)~N^(1-gamma/2)log N. The latter is the smaller term
for every fixed gamma<2/3. Equation(9) proves(1).

## 6. General tolerance and exact scope

The proof also works with eta=o(1), p~eta^(-1/2),
r~sqrt(N)eta^(3/4), ell=log(ep), provided r>=C ell. Then

```math
 Q(W)\le Q(A)+C\sqrt\eta N^{3/2},\quad
 T=C\eta N^{3/2},\quad \rho=c\eta^{3/2}/\ell,
 \log|E_W(T)|\le C N\eta^{1/4}\ell,
 \Pr(\text{escape})\le C\exp[-cN\sqrt\eta\ell].       \tag{10}
```

A convenient sufficient condition is
eta >> N^(-2/3)(log N)^(4/3). This substantially weakens the
eta >> N^(-1/3) restriction of the independent third-moment
Lindeberg proof; the latter remains valid and preserved. Mixed tails
pay their jump geometry through gamma1 instead of a global replacement
error. This is the decisive combination, not a renamed Gaussian bound.

The resulting random neighborhood is genuinely made of signings, not
weighted Gaussian matrices. But the code is near-ground for W, not
necessarily for the starting A, and the preparation cost can dominate
the perturbative value gain. No propagation to comparable larger orders,
favorable parent value, or improved original interval follows. The theorem
removes a realization gap in local stability; it does not remove the
global cross-order optimization gap.

## 7. A quantitative cap envelope, not a favorable original-value theorem

On the same exponential-probability scale one also has

```math
 |Q(W')-(1-2\rho)Q(W)|\le C T.                       \tag{11}
```

For the upper bound apply(4) on the inner level d<=T, then the same
all-shell argument outside, with half-deficit budgets. The anchored
response there is at most C T. The unanchored ground response itself
has variance at most C rho N^2 and jumps at most2, so scalar Bernstein
bounds its absolute value by C T with failure at most the right side
of(9). The ground state alone proves the corresponding lower bound.
Integrating the shell and anchor tails also gives
E|Q(W')-(1-2rho)Q(W)|<=C T.

Here T~eta N^(3/2), whereas rho~eta^(3/2)/log(1/eta). In particular
the error is LARGER than the apparent shrinkage2rho Q(W) at bounded
normalized cap. Equation(11) cannot be advertised as an improvement to
M_N or a sublinear cross-order recurrence. It is genuine local sign
stability, with its value error explicitly paid.
