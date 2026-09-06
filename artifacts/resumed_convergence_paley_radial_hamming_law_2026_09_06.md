# Exact Paley ground-state radial law and explicit Hamming stability

Date: 2026-09-06. A Paley-specific sharpening of the general Cayley
Hamming theorem. Its finite symmetry identity is elementary; the
all-order asymptotic consequence uses the newly audited Paley saturation
theorem. No radial ground-state law is asserted for arbitrary minimizers.

Let `q=1 mod 4` be an odd prime power and let the hollow Paley core be

\[
(A_q)_{i,j}=\chi_q(i-j),\qquad \chi_q(0)=0.
\]

Write `H_A(f)=f^T A f/2` and `Q(A)=max_f |H_A(f)|`.

## 1. Exact radial law from one spin

Fix any Boolean `f`. Choose `b` uniformly in `F_q^*` and `t` uniformly
in `F_q`, independently. Set

\[
x_i=f(bi+t),\qquad \sigma=\chi_q(b).
\]

Then

\[
\mathbb E[\sigma x_i x_j]
=\frac{f^T A_qf}{q(q-1)}(A_q)_{ij}\quad(i\ne j),
\qquad \mathbb E\sigma=0.                           \tag{1}
\]

For fixed `i!=j`, the map `(b,t)->(u,v)=(bi+t,bj+t)` is a bijection onto
the ordered distinct pairs of field elements. Since
`b=(u-v)/(i-j)`, the expectation is

\[
\frac{\chi_q(i-j)}{q(q-1)}
\sum_{u\ne v}\chi_q(u-v)f(u)f(v),
\]

which is (1). Also a change of variables gives the exact energy identity

\[
H_{A_q}(x)=\chi_q(b)H_{A_q}(f).                     \tag{2}
\]

If `f` is an absolute ground state, multiply the displayed orientation
by `sign(H_A(f))`. The resulting law is supported entirely on exact
oriented ground states and has covariance

\[
\boxed{\quad
\mathbb E[\sigma x x^T]
=\rho_q A_q,
\qquad \rho_q=\frac{2Q(A_q)}{q(q-1)}.\quad}         \tag{3}
\]

Its diagonal is zero. This is an exact finite statement for every
admissible `q`; it does not require a Boolean eigenvector or the
asymptotic saturation theorem.

A nonsquare dilation in (2) also shows that the positive and negative
Paley-core extrema are exactly equal in magnitude at every finite order.

## 2. Exact comparison with an arbitrary target matrix

For any real symmetric hollow `B` of the same order, (3) gives

\[
\boxed{\quad
Q(B)\ge \rho_q
\left|\sum_{i<j}(A_q)_{ij}B_{ij}\right|.\quad}      \tag{4}
\]

Indeed the expression without the absolute value is the expected
oriented half-energy of `B`; it lies in `[-Q(B),Q(B)]`.

If `B` is a signing differing from `A_q` on exactly `d` unordered edges,
then

\[
\boxed{\quad
Q(B)\ge Q(A_q)
\left|1-\frac{4d}{q(q-1)}\right|.\quad}            \tag{5}
\]

This also treats a target near `-A_q`, by the absolute value. No operator
bound, regularity, or near-optimality assumption on `B` is used.

There is also an exact discrete local-optimality consequence. Every
oriented energy of an order-`q` signing has parity `binom(q,2)`, hence
distinct possible cap values differ by at least two. If

\[
d<\frac{q(q-1)}{2Q(A_q)},                           \tag{5a}
\]

then (5) gives `Q(B)>Q(A_q)-2`, so necessarily `Q(B)>=Q(A_q)`.
Thus fewer than `(1-o(1))sqrt(q)` coefficient flips cannot lower a
Paley core's cap. This is a radius of discrete local optimality; it
asserts neither global optimality nor suboptimality of the Paley core.

The law and inequality are invariant under simultaneously switching or
relabeling the Paley signing and its target. Thus the comparison applies
to every switched and relabeled copy, not just the displayed field
indexing.

## 3. Explicit asymptotic exclusion radius for actual minimizers

The all-order Paley saturation theorem proves

\[
Q(A_q)/q^{3/2}\longrightarrow\frac12
\]

along all odd prime powers `q=1 mod 4`. It is independently proved in
`resumed_convergence_paley_squarefree_audit_2026_09_06.md`.
Combining it with (5), uniformly for `0<=d<=q(q-1)/4`, gives

\[
\frac{Q(B)}{q^{3/2}}
\ge\frac12-\frac{2d}{q(q-1)}-o(1).                 \tag{6}
\]

In particular a sequence of actual minimizers, or any signings at all,
with normalized cap at most `1/2-delta` must satisfy

\[
\boxed{\quad d\ge(\delta/2-o(1))q^2.\quad}         \tag{7}
\]

The same lower bound applies to distance from `-A_q`, using the other
endpoint in (5). More generally, if the limiting target coefficient is
at most `c<1/2`, the asymptotic edge-edit fraction lies between
`1/4-c/2` and `1/4+c/2` relative to a fixed Paley core.

For the bordered Paley conference of order `q+1`, take the target's
finite core and apply (6). Its cap is at most the full target's cap by
independent averaging over the border spin. The core edit count is no
larger than the total edit count, and changing the normalization from
`q` to `q+1` is `o(1)`. Therefore the same asymptotic exclusion radius
`(delta/2-o(1))(q+1)^2` holds around every bordered Paley conference.

## 4. Relation to the general theorem and remaining boundary

The broader result in
`resumed_bound_audit_cayley_hamming_stability_2026_09_06.md` works for
every spectrally near-Frobenius-optimal finite-field additive-Cayley
family, using a spread polynomial covariance and Grothendieck transfer.
The explicit coefficient in (7) is special to Paley affine symmetry:
it eliminates the Ramsey-host and rounding losses entirely.

Neither result proves that arbitrary minimizers can be recovered into
one of these families. Instead, they show that a hypothetical persistent
sub-`1/2` minimizing sequence must remain a positive edge-edit fraction
away. They therefore obstruct one substantial class of proposed
upper-preserving landing maps, while leaving original convergence open.

The integer verifier
`computations/resumed_convergence_paley_radial_verify_2026_09_06.py`
enumerates all projective spins at prime orders `5,13,17`, then every
affine sample of a positive ground. It checks both the constant oriented
ground energy and the exact unnormalized covariance identity
`sum_(b,t) chi(b) xx^T=(f^T A f) A`. All checks pass; the core caps are
respectively `4,20,32`.
