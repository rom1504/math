# Recursive Hadamard bases: a uniform orbital bound and a type recursion

Date: 2026-09-06. These are one-row upper-bound tools for the restricted
weave, not an original cap or convergence theorem.

## 1. Gaussian orbital norm

Let `G_m` be the full signed-permutation group on `R^m`. For `t>0`, put

\[
P_t(v)=\mathbb E_{g\in G_m}e^{-t\|v-gv\|^2},
\qquad L_t(v)=\sqrt{P_t(v)}.
\tag{1}
\]

Then `P_t(v)=per K_t[|v|]/m!`, with
`K_t(a,b)=(e^{-t(a-b)^2}+e^{-t(a+b)^2})/2`.

**Uniform orbit theorem.** For fixed `t,C>0`, every orthogonal `U`, and
every `v` with `||v||^2<=Cm`,

\[
\mathbb E_g L_t(Ugv)
\le \sqrt2\exp\!\left\{\frac\pi2\sqrt{D/3}\right\}L_t(v),
\qquad D=\left\lceil e^2(2tC+1)m\right\rceil.
\tag{2}
\]

In particular its loss is `exp(O_{t,C}(sqrt(m)))`, uniformly in `U,v`.
The orthogonal matrix may be random provided it is independent of the
fresh uniform `g`, conditionally on `v`.

A useful energy-sensitive version takes `D=ceil(8t||v||^2)+1`. The
Poisson tail is again at most `exp(-4t||v||^2)<=P_t(v)` (the zero-energy
case is immediate). This gives logarithmic loss `O(1+sqrt(t||v||^2))`.
For `b` terminal blocks sharing total energy at most `Cm`, their total
loss is `O_t(b+sqrt(bCm))`. This observation does not remove the separate
finite-alphabet type-error obligation when the depth grows.

### Proof

Use the Gaussian symmetric-Fock feature

\[
\phi(v)=e^{-t\|v\|^2}
 \bigoplus_{d\ge0}\frac{(2t)^{d/2}}{\sqrt{d!}}v^{\otimes d},
\qquad \langle\phi(v),\phi(w)\rangle=e^{-t\|v-w\|^2}.
\]

Orthogonal transformations act unitarily on each degree. Write `Q_G`
for projection to the `G_m`-invariant subspace. Then
`P_t(v)=||Q_G phi(v)||^2`. The orbit covariance operator

\[
C_v=\mathbb E_g |\phi(gv)\rangle\langle\phi(gv)|
\]

has operator norm exactly `P_t(v)`. Indeed its nonzero eigenvalues agree
with those of the normalized finite orbit Gram matrix. That matrix has
nonnegative entries and constant row sum `P_t(v)`, so this row sum is
its largest eigenvalue. Repetitions in the orbit do not affect this argument.

Let `Pi_D` project to degrees at most `D`. The squared norm of the omitted
feature is the upper tail of a Poisson variable with mean `2t||v||^2`.
Since degrees are orthogonal, commute with both group actions, and the
invariant projection is contractive,

\[
\mathbb E_gL_t(Ugv)^2
\le \operatorname{rank}(Q_G\Pi_D)P_t(v)
 +\Pr\{\operatorname{Poisson}(2t\|v\|^2)>D\}.
\tag{3}
\]

An invariant homogeneous polynomial has even degree `2j` and is a
symmetric polynomial in `v_1^2,...,v_m^2`. Its dimension is the number
of partitions of `j` into at most `m` parts. Hence

\[
\operatorname{rank}(Q_G\Pi_D)
\le\sum_{j\le D/2}p(j)\le e^{\pi\sqrt{D/3}}.
\]

The last estimate follows directly from the partition generating function:
for `s>0`, its logarithm is at most `pi^2/(6s)`, and optimizing
`s(D/2)+pi^2/(6s)` gives the displayed bound.

The Poisson Chernoff bound and the choice of `D` make the tail at most
`exp(-(4tC+1)m)`. On the other hand each term in (1) is at least
`exp(-4tCm)`. Thus the tail is at most `P_t(v)`, and taking the square
root in (3) proves (2).

## 2. Two exact permanent projection inequalities

For PSD Gram matrices `A,B` of common order `n`,

\[
\frac{\operatorname{per}(A\circ B)}{n!}
\ge \frac{\operatorname{per}A}{n!}
       \frac{\operatorname{per}B}{n!}.
\tag{4}
\]

If `A` and `B` are the Gram matrices of vectors `v_i` and `w_i`, set
`v=otimes_i v_i`, `w=otimes_i w_i`. The diagonal symmetric-group
projection `(1/n!)sum_sigma P_sigma otimes P_sigma` contains the range
of the separate projection `P_sym otimes P_sym`. Their difference is
PSD; its quadratic form on `v otimes w` is exactly (4).

For concatenated coordinate lists `v=(v_+,v_-)`,

\[
L_t(v)\le L_t(v_+)L_t(v_-).
\tag{5}
\]

One proof is symmetrization of the tensor product of their Gaussian-kernel
feature tensors. Equivalently the `G_m`-invariant subspace is contained
in the `G_{m_+} times G_{m_-}`-invariant subspace; compare the corresponding
projection norms of the Gaussian Fock feature.

Finally diagonal-coordinate deletion in the weave incurs only a polynomial
factor:

\[
\max_i L_t(v\setminus v_i)\le\sqrt{2m}\,L_t(v).
\tag{6}
\]

Indeed the permanent terms fixing coordinate `i` give
`P_t(v)>=K_t(v_i,v_i)P_t(v\setminus v_i)/m`, and `K_t(a,a)>=1/2`.
No delocalization of the removed coordinate is needed.

## 3. A finite-depth randomized Hadamard ensemble

Write `U_s=H_s/sqrt(s)`. At a node of order `s=2q`, independently sample
the two child Hadamard bases and a fresh uniform `g in G_s`, and set

\[
U_s=\operatorname{diag}(U_q^{(1)},U_q^{(2)})
 \frac1{\sqrt2}\begin{pmatrix}I&I\\I&-I\end{pmatrix}g.
\tag{7}
\]

Every entry has magnitude `1/sqrt(s)` and the matrix is orthogonal.
At each terminal node use any fixed Hadamard of the required order,
followed on its input side by an independent uniform signed permutation.
All these random choices are independent. Thus, conditionally on a node's
input vector, its signed permutation is fully uniform and independent of
both child bases. No adaptive choice of a child matrix is being assumed.

For a vector of signed empirical type, let `nu` denote its symmetrized
law: its atom at zero is unchanged and its mass at each positive magnitude
is split equally between the two signs. For a symmetric finite law `nu`,
define

\[
F_t(\nu)=\inf_{\eta\in\Pi(\nu,\nu)}
 \{D(\eta\Vert\nu\otimes\nu)+t\mathbb E_\eta(X-Y)^2\},
\quad \Phi_t(\nu)=-F_t(\nu)/2.
\tag{8}
\]

This is the folded-kernel permanent exponent, because optimizing the
relative signs conditionally on the two magnitudes produces `K_t`.

Let `pi` range over signed pair laws whose average **absolute** marginal
is the absolute law of `nu`. Let `nu_+` and `nu_-` be the symmetrized
laws of `(A+B)/sqrt(2)` and `(A-B)/sqrt(2)`, respectively. Define

\[
(\mathcal B f)(\nu)=\sup_\pi
 \left\{\tfrac12 f(\nu_+)+\tfrac12 f(\nu_-)
       -\tfrac12D(\pi\Vert\nu\otimes\nu)\right\}.
\tag{9}
\]

The input signed-word count is
`m! 2^(#nonzero)/prod_magnitude count!`. An ordered pair table has
`(m/2)!/prod_pair count!` realizations. Thus the probability of its
empirical law `pi` has normalized log probability
`-D(pi||nu otimes nu)/2+o(1)`. The identity holds even if its two signed
marginals are unequal: their average absolute marginal is prescribed,
and `nu` has equal weights on opposite signs.

For each **fixed** depth `r` and finite initial alphabet, the exact
recursion (7), (5), the terminal estimate (2), and the method of types give

\[
\limsup_{m\to\infty}\frac1m
 \log\mathbb E L_t(U_m v_m)
\le (\mathcal B^r\Phi_t)(\nu),
\tag{10}
\]

whenever the symmetrized type of `v_m` tends to `nu` and its squared norm
is `O(m)`. The limit uses orders divisible by `2^r` with terminal Hadamard
order `m/2^r`. After fixed depth `r`, every alphabet is still finite.
Every terminal norm is bounded by the root norm; relative to its terminal
dimension the constant in (2) is at most `2^r` times the root constant.
There are only `2^r` terminal estimates, whose total logarithmic loss is
`O_{r,t,C}(sqrt(m))=o(m)`. All intermediate type counts and factorial
errors are polynomial in `m` for this fixed-depth calculation. No limit
uniform in a growing depth is asserted in (10).

## 4. The soft type potential is a genuine supersolution

For a finite law `mu`, put
`C_t(mu)=2H(mu)-F_t(mu)=max_eta[H(eta)-t E(X-Y)^2]`.
The following three properties hold.

1. `F_t` is subadditive under dependent coordinate marginals. Apply (4)
   to repeated empirical samples of a joint law with Gaussian coordinate
   kernels, and take the fixed-type permanent limit.
2. `C_t` is concave in its marginal law: mix two feasible self-couplings
   and use concavity of their entropy. It is also subadditive under
   coordinate marginals, by subadditivity of the coupling entropy.
3. `F_t` is concave in its marginal law. This follows from (5), by
   concatenating repeated empirical words and taking the fixed-type limit.
   In particular symmetrizing a law cannot decrease `F_t`.

For Property 3 on unsymmetrized signed laws, use the unsigned-permutation
version of the same symmetrization inequality, with raw Gaussian kernel
`exp(-t(a-b)^2)`. This proves concavity before the folding operation.

For an admissible pair law `pi` and its orthogonal output law `rho`,
orthogonal invariance gives `F_t(pi)=F_t(rho)` and `C_t(pi)=C_t(rho)`.
Properties 1 and 3 give
`F_t(nu_+)+F_t(nu_-)>=F_t(pi)`. Property 2, reflection invariance, and
the average absolute-marginal constraint give `C_t(pi)<=2C_t(nu)`.
Since `D(pi||nu otimes nu)=2H(nu)-H(pi)`, these imply exactly

\[
\tfrac12\Phi_t(\nu_+)+\tfrac12\Phi_t(\nu_-)
 -\tfrac12D(\pi\Vert\nu\otimes\nu)
\le \Phi_t(\nu).
\tag{11}
\]

The gap is one quarter of the sum of the two nonnegative quantities
`F_t(nu_+)+F_t(nu_-)-F_t(pi)` and `2C_t(nu)-C_t(pi)`.
Thus `B Phi_t<=Phi_t`, and the iterates in (10) decrease pointwise.

At `t=infinity`, `F_infinity=H`, and the gap is

\[
\tfrac14\{D(\pi\Vert\nu\otimes\nu)
 + I(C;D)
 +H(\nu_+)-H(C)+H(\nu_-)-H(D)\}.
\tag{12}
\]

For symmetric `nu=(delta_{-1}+delta_1)/2`, the first Bellman iterate
has gain exactly `(log 2)/4`. This follows by parameterizing the fraction
of pairs with equal signs; the input coupling cost and output dependence
sum to at least `log 2`, with equality for balanced signs.

## 5. Precise remaining cap obligation

For retention `p`, the input type is

\[
\nu_p=(1-p)\delta_0+\tfrac p2\delta_{1/\sqrt p}
                         +\tfrac p2\delta_{-1/\sqrt p}.
\]

In the weave use `H_i=sqrt(m) U_m^{(i)T}`, so that its row spectrum
`H_i[T,:]^T x/sqrt(k)` is precisely `U_m^{(i)}(1_Tx)/sqrt(k/m)`.
Sample an **independent** recursive base in each weave fibre, in addition
to the fresh output-column permutations and edge signs from the soft
Finner theorem. That theorem applies to different fixed fibre bases;
independence then produces `(E Z)^m`, not `E[Z^m]`.

By input sign/permutation invariance, (6), and (10), the expected one-row
weave sum for the ensemble (7) satisfies

\[
\limsup\frac1m\log\mathbb E_H Z_T(t)
\le p\log2+(\mathcal B^r\Phi_t)(\nu_p).
\tag{13}
\]

This holds for any fixed selector `T` of size `k/m -> p`; the new input
permutations randomize its position inside each base matrix.
A sufficient finite variational certificate is therefore

\[
p\log2+(\mathcal B^r\Phi_t)(\nu_p)+t(1-\sqrt p)<0.
\tag{14}
\]

No such inequality is claimed yet. A strict margin would survive a small
positive `eta` in the weave deficit target. It would yield an original
sub-half cap along the constructed orders. The director's separate artifact
`continued_director_recursive_weave_all_order_implication_2026_09_06.md`
supplies the all-order extension using Paley terminal orders and their
asymptotic density. Even an all-order upper improvement would not by itself
settle original convergence.

## 6. A rigorous shallow-depth falsifier

One must retain energy-condensation policies in (9). Choose `A=B` with
law `nu`; then `D(pi||nu otimes nu)=H(nu)`, the plus child has law
`sqrt(2) nu`, and the minus child is identically zero. Repeating only in
the nonzero branch gives

\[
(\mathcal B^r\Phi_t)(\nu)
\ge -(1-2^{-r})H(\nu)+2^{-r}\Phi_{2^rt}(\nu)
\ge -(1-2^{-r-1})H(\nu).
\tag{15}
\]

The last inequality follows by choosing the identity self-coupling in
(8), which gives `F_t(nu)<=H(nu)`. It is valid at every positive `t`,
not just at zero temperature.

For `p=15/16,t=4,r=2`, substituting
`H(nu_p)=h(p)+p log2` in (15) gives a strictly positive lower bound
`0.003676...` on the candidate exponent (14). Thus depth at most two
cannot certify this parameter choice. At depth three the same lower
bound is negative and is inconclusive. A selected policy always gives
a **lower** Bellman bound; it cannot prove the upper inequality (14).

The one-step value at `p=15/16,t=4` is numerically about
`-0.60606065765`, with candidate exponent `+0.17078147792`.
The diagnostic is reproducible by
`computations/continued_convergence_recursive_pairing_2026_09_06.py`.
No numerical optimizer value is being used as a rigorous upper certificate.

## 7. Exactly which input symmetries may be averaged

For arbitrary depth, global reversal `(A,B)->(-A,-B)` and input interchange
`(A,B)->(B,A)` preserve both output **absolute-value** laws individually.
Averaging over these two operations increases pair entropy and hence never
decreases the Bellman objective. Independent reversal of just one input is
different: it interchanges the two output laws, and averaging it can destroy
a beneficial separation of the two children.

At the first level only, the latter averaging is valid. Indeed

`R_t(mu)=H(mu)-F_t(mu)`

is concave: it is the supremum of `H(Y|X)-t E(X-Y)^2` over self-couplings,
and conditional entropy is jointly concave. Consequently, for any mixture
of laws, the entropy Jensen gain dominates the `F_t` Jensen gain. Since
each symmetrized output law is a fixed stochastic pushforward of the pair
law, its entropy gain is at most the pair entropy gain. Thus the first-level
objective has Jensen gain

`(1/2) Delta H(pair) - (1/4)(Delta F(plus)+Delta F(minus)) >= 0`.

It is concave in the pair law. The sign and swap averages therefore reduce
the first-level problem to a symmetric magnitude coupling and independent
relative signs. For ternary `nu_p`, write its off-diagonal magnitude mass
as `u`, with `0<=u<=min(p,1-p)`. The magnitude coupling is

`theta = [[1-p-u,u],[u,p-u]]`.

Both child magnitude laws have weights

`(1-p/2-3u/2, 2u, (p-u)/2)`

on `0, 1/sqrt(2p), sqrt(2/p)`, respectively. This proves the one-dimensional
first-step reduction used by the numerical script.

### A finite-temperature counterexample at the next level

Let `f_1=B Phi_t`, and initially take `nu` uniform on `{-1,1}`. At
`t=infinity`, the exact first-step values are

`f_1(delta_0)=0`,
`f_1(uniform{-sqrt(2),sqrt(2)})=-(3/4)log 2`,
`f_1((1/2)delta_0+(1/4)delta_{-sqrt(2)}+(1/4)delta_{sqrt(2)})=-log 2`.

For the last formula, the one-dimensional reduction above with `p=1/2`
gives, on writing `z=2u`,

`-log 2 + (1-z)[(3/4)log 2-h(1/4)]/2`,

whose maximum is `-log 2` because the bracket is negative.
The pair law `A=-B` therefore gives second-step value `-(7/8)log 2`,
whereas its independent-sign average gives `-log 2`. The average loses
`(log 2)/8`. Here the fully sign-symmetric pair law is unique, so this
actually falsifies the restricted optimization, not merely its proposed proof.

This persists at `t=8`. All child alphabets needed for these first-level
calculations contain at most five points and have minimum spacing at least
one. The conditional-entropy representation gives

`0 <= H(mu)-F_t(mu) <= log(1+4 exp(-t))`.

Thus the averaged policy can improve over its infinite-temperature-limit
value by at most `(1/2)log(1+4 exp(-8))`, still smaller than `(log 2)/8`.
The fully sign-symmetric restriction is therefore also a false upper-bound
reduction at this finite temperature. It remains valid for selected lower
policies. The audit agent independently reconstructed these calculations.

## 8. Exact finite-alphabet asymmetric tree objective

At level `d` the input magnitudes are
`i/sqrt(p 2^d)`, `i=0,...,2^d`. After the always-valid global-reversal and
swap averaging, a pair class is `(i,j,s)` with `i<=j` and relative sign
`s in {+1,-1}`; only one sign class is needed when `i=0`. Its multiplicity
is one for `(0,0)`, two for `i=j>0`, and four for `i<j`.
If the total class masses are `q_c`, the signed pair entropy is

`H_pair(q) = -sum_c q_c log(q_c/multiplicity_c)`.

The input magnitude law is obtained by giving half of each class mass to
each of its endpoints. The two output magnitude maps send that mass to
`|i+s j|` and `|i-s j|`, respectively. Hence every consistency constraint
between adjacent levels is linear. With `H_signed(p)=H(p)+(1-p_0)log 2`,
the exact depth-`r` finite objective is

`sum_{nodes v at depth d<r} 2^{-d} [H_pair(q_v)/2-H_signed(p_v)]`
`+ 2^{-r} sum_{leaves w} Phi_t(p_w)`.

Maximizing this expression subject to the linear marginal/child constraints
is precisely the Bellman iterate; independent relative signs have not been
imposed. At depth three there are 122 class-mass variables and 28 independent
linear equations. This is the frozen finite variational problem, not an
asymptotic unknown row partition function.

`computations/continued_convergence_asymmetric_bellman_2026_09_06.py`
implements the full objective, an analytic gradient, and SLSQP lower-policy
diagnostics. Its permanent evaluation additionally repairs the numerical
Sinkhorn coupling to an exactly self-marginal primal coupling by uniform
scaling and adding the nonnegative residual to the diagonal. This gives a
one-sided feasible-coupling evaluation at each leaf, although floating-point
policy feasibility is still not an interval proof. No numerical local optimum
is asserted to upper-bound the Bellman maximum.

### An interval-certified depth-three falsifier

At `p=15/16,t=4`, the full asymmetric depth-three objective has an explicitly
admissible rational policy with

`p log 2 + B^3 Phi_4(nu_p) + 4(1-sqrt(p))`
`>= 0.00152671438833974810585956660916941383395981761617838`.

The class masses have denominator `2^30`. Each terminal magnitude law is
equipped with an explicit rational self-coupling of denominator `2^50`.
Substituting these feasible couplings in the infimum defining `F` gives an
upper bound on `F`, hence a **lower** bound on `Phi` and the Bellman value.
All marginal and parent/child constraints are verified by integer arithmetic.
The remaining logarithms, exponentials and square root are enclosed using
60-decimal-digit interval arithmetic. Thus the positive sign does not rely
on floating-point optimization or a small Sinkhorn residual.

The certificate is
`computations/continued_convergence_bellman_depth3_certificate_2026_09_06.json`.
Verification, without running an optimizer, is

```
.venv/bin/python computations/continued_convergence_bellman_certificate_2026_09_06.py computations/continued_convergence_bellman_depth3_certificate_2026_09_06.json
```

This excludes depth three at this particular `p,t`; by monotonicity it also
excludes all shallower depths. It does not exclude other temperatures,
retentions, deeper recursions, or the actual recursive ensemble. In particular,
the more favorable greedy depth-three value was not a valid upper estimate:
optimizing earlier pairings for their eventual descendants changes the answer
enough to reverse the sign of the candidate exponent. The best policies found
at this particular parameter happen to have nearly independent relative
signs at every node; the need to permit asymmetric policies is established
by the separate exact counterexample in Section 7, not by this numerical
optimizer's observed shape.

## 9. Uniform exclusion of every depth at most two

The strict finite criterion (14) is impossible for **every** fixed
`0<p<1`, `t>0`, and `r<=2`. This statement concerns the criterion, not the
actual ensemble. It follows from two explicit lower policies and a small
interval certificate.

The first policy uses independent signed inputs at both levels, paying
zero KL at every node. If `X_1,...,X_4` are iid with law `nu_p`, it gives

`B^2 Phi_t(nu_p) >= Phi_t(law((X_1+...+X_4)/2))`.

Writing `delta=1-p`, the absolute law on `i/(2sqrt(p))`, `i=0,...,4`, is

```
w0 = delta^4 + 3p^2 delta^2 + 3p^4/8,
w1 = 4p delta^3 + 3p^3 delta,
w2 = 3p^2 delta^2 + p^4/2,
w3 = p^3 delta,
w4 = p^4/8.
```

The second policy is (15), giving `B^2 Phi_t(nu_p)>=-7H(nu_p)/8`.
Consequently the criterion's left side is at least the maximum of

`I(p,t)=p log2-F_t(w_p)/2+t(1-sqrt(p))`

and

`C(p,t)=p log2/8-7h(p)/8+t(1-sqrt(p))`.

Here `F_t(w_p)` uses the folded magnitude kernel, so it equals the signed
cost for the displayed fourfold-sum law.

### Analytic reduction to a compact rectangle

The earlier Gaussian-profile obstruction applies to every Hadamard base,
uniformly in its choice. The recursive ensemble has fresh input signed
permutations, so a fixed selector can be averaged to a uniform selector.
Combining its actual one-row pressure lower bound with (13) shows that
the Bellman criterion also obeys that obstruction. In particular,

`p log2 + (1/4)log(1-p) > 0` for `0<p<=0.92`:

this function is concave, vanishes at zero, and its value at `0.92` is
interval-enclosed above `0.006263245`. Hence those retentions are excluded.

For `p>=0.981`, `C(p,0)>0`: it is positive at `0.981` (above
`0.0026409869`) and its derivative is positive for `p>1/2`.
For the remaining `p in [0.92,0.981]`, the function `C(p,5)` is strictly
convex, because

`d^2 C(p,5)/dp^2 = 7/[8p(1-p)] + 5/[4p^(3/2)] > 0`.

Its tangent at `p=0.945` is interval-enclosed above `0.03488078`
throughout this entire retention interval. Since `C` increases in `t`,
all `t>=5` are excluded. Only

`[0.92,0.981] x [0,5]`

remains.

### The finite interval certificate

The remaining rectangle is partitioned into 237 rectangles with exact
rational endpoints. On each, either interval evaluation gives `C>0`, or
an explicit self-coupling gives `I>0`. For the latter, symmetric nonnegative
coefficients `c_ij` with denominator `2^40` define

`Q_ij(p)=c_ij w_i(p)w_j(p)` for `i!=j`,
`Q_ii(p)=w_i(p)[1-sum_{j!=i} c_ij w_j(p)]`.

Interval arithmetic verifies the bracket is nonnegative throughout its
rectangle. Thus `Q(p)` has exactly the required self-marginals for every
retention in the box. Evaluating its entropy-plus-kernel cost upper-bounds
`F_t(w_p)`, giving the required lower bound on `I`. The kernel is evaluated
without any square-root interval dependency through

`log K_t(i/(2sqrt(p)),j/(2sqrt(p)))`
`= -t(i-j)^2/(4p)+log(1+exp(-tij/p))-log2`.

The smallest certified rectangle lower bound is greater than
`0.0000185621595342`. A standalone verifier checks every saved coefficient,
all interval inequalities, exact rational total area, and pairwise interior
disjointness of the boxes. The area and containment checks ensure the closed
rectangles cover the target rectangle, not merely a finite grid of points.
The verifier does not call a numerical optimizer or import its libraries.

```
.venv/bin/python computations/continued_convergence_depth2_uniform_exclusion_2026_09_06.py --verify --precision 60 --output computations/continued_convergence_depth2_uniform_exclusion_2026_09_06.json
```

Both generation at 35 digits and independent replay at 60 digits pass.
Finally, monotonicity `B^r Phi_t >= B^2 Phi_t` for `r<=2` proves the claimed
uniform shallow-depth exclusion. Deeper Bellman iterates and the underlying
construction remain open; no conclusion about original convergence follows.
