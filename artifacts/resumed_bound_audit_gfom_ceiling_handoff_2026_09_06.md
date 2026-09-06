# Independent bounded check of the involution ceiling and finite-tree handoff

Date: 2026-09-06. Scope: the original-class fixed-algorithm ceiling, not
the Boolean optimum or a growing-depth algorithm.

## 1. Primary hypotheses checked

I read Wang--Zhong--Fan, *Universality of approximate message passing
algorithms and tensor networks*, arXiv:2206.13037v5, Definition 2.6,
Proposition 2.7(b)(2), Theorem 2.8, and Lemma 2.14/Remark 2.15 directly
in the [primary PDF](https://arxiv.org/pdf/2206.13037).

The explicit power test in Proposition 2.7(b)(2) applies to a flat
deterministic symmetric involution after independent signed-permutation
conjugation: every even power is I, every odd power is B, and its limiting
spectral law is balanced Bernoulli because the diagonal is uniformly o(1).
Its norm is exactly one. Theorem 2.8 requires nondegenerate AMP covariance
and Lipschitz state maps; it is not a theorem for arbitrary bounded
functions of polynomial iterates. The convergence audit addresses these
requirements by exact GFOM-to-AMP serialization and fresh Gaussian
regularization, with the perturbation removed after the matrix limit.

I reconstructed those reductions and the Haar comparator calculation in
`resumed_convergence_involution_local_ceiling_2026_09_06.md`. In particular,
the positive-definite covariance regularization uses a fresh independent
Gaussian input before each matrix multiplication, not noise recycled from
an earlier register. The zero spectral mean makes the diagonal Onsager
coefficient zero; a fresh orthogonal image contributes independent
variance delta^2 in the new residual. Fixed Lipschitz recursions are stable
in normalized mean square as delta tends to zero.

For Haar O diag(I,-I) O^T, a fixed Boolean energy divided by n is
2 Beta(n/4,n/4)-1. Its exponential tail rate is
one quarter times log(1-t^2). The union bound over 2^n Boolean vectors
therefore gives the normalized half-energy upper bound sqrt(15)/8.
Removing the uniformly vanishing diagonal correctly extends this to
cube-valued outputs. Thus the numerical factor and its applicability to
fixed equivariant Lipschitz GFOM outputs survive this independent check.

## 2. Finite injective trees really have the required GFOM approximation

The warning about exact injective formulas in the convergence artifact
is appropriate: it needs a proof, supplied here for this particular tree
family. Work on a bounded-operator hollow flat-sign sequence, including
the hollow symmetric Hadamard subsequence.

For an old tree T let h_T be its normalized even child Hermite product.
The one-level target is

    Y_T=B[S h_T(X_children)].                         (1)

Under sign inputs, the local difference between h_T(X_children) and the
fully injective child Wick forest has L2 norm o(1), uniformly at the
local root. This follows directly by expanding the squared difference:
each of the two norms and their cross moment has the same leading
whole-tree Hermite matching patterns; extra identifications or nonempty
free parity edges vanish. This is the sign-input local Wick lemma used
in the full nonlinear channel proof, not an inference from Gaussian L2
error for an arbitrary nonmultilinear polynomial.

Both local terms exclude their own root spin. Reduce them on the cube
to their fixed-degree Walsh expansions. For own-coordinate-free e_j,

    E|sum_j c_j S_j e_j|^2
       <=(D+1) sum_j c_j^2 E|e_j|^2,                (2)

where D bounds their degree. A resulting Walsh coefficient can have at
most D+1 choices for the newly inserted root label; Cauchy--Schwarz
proves (2). Because each B row has squared norm one, (2) transports the
local Wick error in (1). Deleting collisions with the output root has
the usual vanishing fixed-root squarefree kernel norm. Consequently

    avg E|Y_T-X_T|^2=o(1).                           (3)

Now approximate the polynomial h_T by a bounded globally Lipschitz even
map h_T,R, agreeing on a large box. The exact child fields have uniformly
bounded moments of every fixed order, so their averaged L2 truncation
error tends to zero as R tends to infinity. At fixed R, replace each
exact child field by an already constructed GFOM approximation. The
Lipschitz constant and ||B||op bound control the resulting normalized
L2 error. Multiplication by the root sign is made globally Lipschitz
on all real seed inputs by clipping that seed to [-1,1]; it does not
change its actual Rademacher value. All local maps remain jointly odd
in the spin registers. Thus (3), induction over the finite ancestor
tree bank, and successively chosen finite truncation levels yield any
prescribed normalized L2 accuracy.

The choices are made with the target accuracy fixed before n. One does
NOT recursively substitute unbounded polynomials and assume that a
small averaged L2 error survives arbitrary polynomial multiplication.
Instead each parent is truncated first and its finite Lipschitz constant
determines the child accuracy required. There are only finitely many
such choices.

## 3. Feasibility and the discontinuous final sign

Approximate a finite feasible odd/even pair F,H by smooth bounded
coordinate maps preserving |F|+H<=1 and 0<=H<=1; conditional expectation
and joint Ornstein--Uhlenbeck smoothing give this, followed if necessary
by a symmetric globally Lipschitz approximation and projection onto the
pointwise feasible diamond. Realize their finite old fields by Section 2.

Do not assume that sign(BF) has a normalized L2 approximation when BF
has mass near zero. Use the globally Lipschitz odd softsign

    tau_eta(x)=clip(x/eta,-1,1),
    x tau_eta(x)>=|x|-eta.

The two feasible outputs are

    m_plus= F+H tau_eta(BF),
    m_minus=-F+H tau_eta(BF).                        (4)

They are fixed finite equivariant GFOM computations after all finite
truncation choices. Their half energy difference is exactly
F^T B[H tau_eta(BF)], whose expectation divided by n differs from
avg E H|BF| by at most eta. The normalized energy is stable under the
preceding L2 approximations by bounded operator norm and feasibility.

Apply the fixed-GFOM ceiling to each output in (4). Its conclusion needs
no approximation of an exact injective formula on the Haar comparator:
the realized Lipschitz outputs themselves are feasible on the comparator,
so its cube cap applies directly. Remove the fixed approximation errors
and then eta. Finally take the supremum over finite Gaussian pairs.

Therefore the full-response variational supremum proved by this campaign
obeys

    sup_(feasible finite F,H) E H E_N|K_F+tau_F N|
       <=sqrt(15)/8=0.4841229182759271... .           (5)

The infinite Gaussian closure cannot evade (5), because its functional
is the L2 limit of finite feasible constructions. This blocks an attempt
to reach 1/2 merely by increasing the finite response bank and taking
its variational supremum. It does not address algorithms whose size
grows with n, matrix-adapted seeds, nonequivariant rules, nonlocal
constructions, or the actual Boolean ground state of Hadamard signings.
