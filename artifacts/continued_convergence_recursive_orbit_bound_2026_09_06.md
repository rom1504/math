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
sub-half cap along the constructed orders. An all-order upper realization
requires a separately stated dense family of terminal Hadamard orders;
even that upper improvement would not by itself settle original convergence.

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
