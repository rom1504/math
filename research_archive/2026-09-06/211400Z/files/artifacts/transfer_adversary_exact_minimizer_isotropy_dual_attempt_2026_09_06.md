# Frozen attempt: exact global minimality and isotropic mean slack

2026-09-06. Status: **unresolved, bounded attempt closed without a new
global-minimizer theorem.** This note preserves the algebra tested during
the campaign rather than leaving it only in agent messages.

For a hollow signing A let

```
I(A)=max_{mu: E_mu xx^T=I} E_mu |H_A(x)|,
Q(A)=max_x |H_A(x)|.
```

The two distinct possible targets were:

1. Every exact minimizer A of order n obeys `I(A)>=Q(A)-O(n)`, with a
   universal constant.
2. At every order there exists at least one exact minimizer with that
   property.

Neither was proved or refuted. The finite rational LP screen in
`transfer_adversary_exact_minimizer_isotropy_finite_2026_09_06.md` is
compatible with both. The scalable planted-clique and balanced near-
minimizer counterexamples are NOT exact minimizers and do not settle
either quantifier.

## Exact finite LP dual and the sign obstruction to a naive descent

Finite-dimensional LP duality gives

```
I(A)=min t
subject to t + sum_(i<j) h_ij x_i x_j >= |H_A(x)| for every Boolean x,
```

where the h_ij are unrestricted real numbers. If `I(A)=Q(A)-g`, an
optimal dual therefore has `H_h(x)>=g` at BOTH positive and negative
absolute ground states of A.

This is a common unsigned-slack certificate, not a signed energy descent
certificate. The continuous move `A -> A-epsilon*h` lowers H_A at a
positive ground but makes H_A more negative at a negative ground. It
does not lower the absolute cap on both orientations. Reversing the move
exchanges the difficulty. In addition, such a move generally does not
preserve the required edge set {-1,+1}.

The order-four exact example from the finite screen already displays
the dual without approximation. Take all edges positive except edge
{2,3}. Then

```
|H_A(x)|=2+x_0 x_1+x_2 x_3,
Q(A)=4, I(A)=2.
```

Thus t=2 and `H_h=x_0x_1+x_2x_3` form an exact dual optimum. It is a
mean-slack witness, not a simultaneous oriented descent theorem.

## What a discrete perturbation would actually have to prove

Flipping an edge set E changes the Hamiltonian exactly by

```
H_{A^E}(x)=H_A(x)-2 sum_({i,j} in E) A_ij x_i x_j.
```

A global-minimizer contradiction requires a SINGLE admissible edge set
whose resulting absolute Hamiltonian is below Q(A) at EVERY spin. A
dual inequality against an isotropic average does not supply that
uniform signed conclusion. Randomizing a matrix perturbation also does
not cure this automatically: expectation and maximum have the wrong
inequality direction for an upper-cap certificate.

The continuous box relaxation has the zero matrix available and hence
does not retain the discrete minimization problem. John-type continuous
stationarity or a local edge condition cannot be substituted for exact
global optimality. The scalable Hadamard local traps in
`transfer_adversary_exact_minimizer_orientation_gap_2026_09_06.md`
explicitly warn against that substitution.

No workable actual-sign perturbation was obtained from the LP dual.
This closes this particular route, not the two original global-
minimizer isotropy questions.
