# Exact scope of the XOR-game and Fourier-type imports

Date: 2026-09-06. Bounded primary-source pass for the stabilized Hadamard
seed norm. No imported theorem in this pass establishes `R=T` or a strict
`R<T` seed. This file records exact failed mappings so that recent results
are not confused with the required catalyst theorem.

## 1. The exact catalyst quantity

For a real symmetric seed `B`, write

```math
\beta(B)=\max_{x,y\in\{-1,1\}^n}|x^TBy|,\qquad
R(B)=\frac12\sup_{r\ge0}\frac{\beta(H_4^{\otimes r}\otimes B)}{(4^r)^{3/2}},
```

where the order-four generator is a regular symmetric Hadamard. In the
archived stabilization theorem the same-spin and bilinear limits agree.
Define separately

```math
2T(B)=\min_{P\succeq B,\ P\succeq-B}\max_{x\in\{-1,1\}^n}x^TPx.
```

The inequality `R<=T` has a direct proof independent of the name of any
operator ideal. For a feasible `P`, write `B=P^(1/2) J P^(1/2)` with a
self-adjoint contraction `J` on the support of `P`. For an outer Hadamard
of order `s`, its normalized version is orthogonal. Blockwise Cauchy--Schwarz
therefore bounds every bilinear Boolean objective by

```math
\sqrt{s}\left(\sum_a x_a^TPx_a\right)^{1/2}
\left(\sum_b y_b^TPy_b\right)^{1/2}
\le s^{3/2}\max_x x^TPx.
```

This works for every orthogonal outer, flat or not. Equality with the
supremum over PRESCRIBED flat Walsh outers is the missing theorem; an
arbitrary orthogonal polar optimizer is not a Walsh realization.

## 2. The recent repeated-CHSH theorem is about a different product

Primary source checked directly: Andris Ambainis,
[Optimal bounds on the classical value of the repeated CHSH game,
arXiv:2608.16439v1](https://arxiv.org/html/2608.16439v1), 17 August 2026,
Sections 2--4.

Its parallel product requires the players to output answer TUPLES and win
EVERY component. The theorem concerns its exponential winning rate and a
nonnegative answer-vector relaxation. The matrix product in `R`, by
contrast, is an XOR/parity product with one output bit per player.

This distinction is already exact at two copies. Let
`C=[[1,1],[1,-1]]`, the unnormalized CHSH sign matrix. For the parity game,
the two-copy classical bias is `beta(C tensor C)/16=8/16=1/2`, and its
winning probability is `3/4`. The all-win two-copy probability in the
source is `10/16=5/8`. These are different optimization problems with
different answer sets, despite both being denoted by a tensor-product
symbol in their respective settings.

More decisively, at every even number `2r` of parity copies, a bent Boolean
vector saturates the Hadamard spectral bound. Its classical bias is exactly
`2^(-r)`, so its winning probability tends to `1/2`; the all-win probability
instead decays exponentially to zero. The new all-win asymptotic separation
therefore supplies neither `R=T` nor a strict inequality for `R`.

The Dinur--Steurer relaxation in that source uses nonnegative answer vectors
with coordinatewise disjoint supports for different answers to a question.
That is not the uniform Boolean-row cut covariance and orthogonal polar
coupling defining the majorant `T`. No equality of these norms is inferred.

## 3. Ordinary quantum XOR bias is already separated from the target

The exact standard quantum-bias formula is the supremum of
`sum_ij B_ij <u_i,v_j>` over arbitrary families of unit vectors. It appears,
with the probability-weighted sign-matrix convention made explicit, in
Section 4 of the same primary source. Let its unnormalized value be `beta_*`.

An elementary five-cycle signing prevents substituting `beta_*/2` for `R`.
Let `B` have diagonal zero, entry `-1` on the five-cycle and `+1` on its
complement. Then

```math
B\mathbf1=0,\quad B^2=5I-J,\quad
|B|=\sqrt5\left(I-\frac J5\right).
```

Taking `P=|B|` gives

```math
R(B)\le T(B)\le\frac12\max_x x^TPx
=\frac{12\sqrt5}{5}<\frac{5\sqrt5}{2}.
```

The ordinary quantum value is `beta_*(B)=5sqrt5`. The spectral inequality
gives its upper bound. For equality, put `d=4sqrt5/5` and use the unit
vectors
`u_i=d^(-1/2)|B|^(1/2)e_i`,
`v_j=d^(-1/2)sign(B)|B|^(1/2)e_j`.
Their objective is `tr B^2/d=5sqrt5`. Thus the strict separation above is
proved, not merely a convention warning.

This does NOT establish `R<T`; the displayed upper bound for `T` may be
attained by `R`. It establishes only that the usual quantum bias is not
the desired norm, even after taking arbitrarily large prescribed catalysts.

## 4. Fourier-type factorization does not currently close the endpoint gap

Primary source freshly located and its exact theorem statements checked:
Aicke Hinrichs,
[Hilbert space factorization and Fourier type of operators,
Studia Mathematica 145 (2001), 199--212](https://www.impan.pl/shop/en/publication/transaction/download/product/89984).
The Walsh gradations in Theorem 6.1 are vector-valued `L2`-based ideal norms,
and their nonuniform comparison is at a prescribed gradation. The present
quantity uses endpoint `L-infinity`-to-`L1` amplification and a supremum
over unbounded catalyst order for each fixed seed.

Neither replacing the endpoints nor replacing a gradation by the fixed-seed
infinite supremum is justified by that theorem. The earlier detailed audit
in `fresh_fixed_walsh_factorization_final_audit_2026_09_05.md` remains correct
about this scope. No newer exact endpoint catalyst theorem was found in this
bounded pass; this negative search finding is not a claim that none exists.

## Remaining discriminating theorem

Either prove the prescribed Walsh polar realization needed for `R(B)=T(B)`,
or exhibit a symmetric full sign seed with a uniform strict catalyst cap
`R(B)<n^(3/2)/2`. Because the original all-order upper constant is now
strictly below `1/2`, the first outcome would rigorously rule out asymptotically
lossless Hadamard tensor transfer on actual optimizing full seeds.

A possible exact test of the second outcome uses a finite random seed
ensemble: the finite-depth catalyst caps are nondecreasing, so their good
events are nested. A positive lower bound on good-event probability UNIFORM
over every catalyst depth yields one seed good at all depths. Existing
scalar restricted-weave proofs do not provide this uniform bound: tensoring
a fixed seed reuses its random edge signs and fibre bases, and the previous
independent-fibre product moment no longer applies.
