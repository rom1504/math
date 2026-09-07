# Full-sign obstruction to optimizer-independent integrated temperature payment

This is an explicit asymptotic counterexample to paying the contracted-child
temperature loss from the complete Hubbard identities plus a uniform cap or
operator-norm bound alone. It is NOT a sequence of globally minimizing
signings, and does not disprove a comparison using that stronger hypothesis.

## 1. Actual sign matrices with mesoscopic coherent blocks

Let k=2^r, ell=64k and N=k ell. Thus ell=8 sqrt(N). Write H_d for the
symmetric Sylvester Hadamard matrix of order d. Define

    S=H_k-diag(H_k),
    B=[[H_(ell/2),-H_(ell/2)],[-H_(ell/2),H_(ell/2)]],
    A=I_k tensor (J_ell-I_ell) + S tensor B.

Here J_ell is the all-one matrix. A is hollow symmetric and EVERY
off-diagonal entry is +/-1: its diagonal blocks are complete positive
cliques, and its off-diagonal blocks are signed copies of the full-sign B.
No zero or weighted edges have been substituted into the actual signing.

The matrix B annihilates the constant vector, and

    ||B||op=sqrt(2 ell),
    ||S||op <= sqrt(k)+1,
    Loff:=||S tensor B||op <= sqrt(2N)(1+1/sqrt(k)).       (1)

In particular ||A||op/sqrt(N) is bounded, as is Q(A)/N^(3/2).

## 2. Exact coherent maximum and a bound for the opposite polarity

Let m_i be the sum of the spins in clique i, and decompose each spin
block into its constant part and orthogonal part. The off-diagonal
interaction annihilates all constant parts. Set

    Delta=(1/2) sum_i (ell^2-m_i^2).

The squared norm of the total orthogonal part is 2 Delta/ell. Therefore

    H_A(x) <= (N ell-N)/2 - Delta +(Loff/ell) Delta.

For all k>=2, Loff<ell. Hence the positive maximum is EXACTLY

    max_x H_A(x)=(N ell-N)/2,                         (2)

attained by making each clique constant, with arbitrary independent
clique signs. For the other polarity the positive-semidefinite clique
part and (1) give

    -min_x H_A(x) <= N/2 + N Loff/2.                  (3)

The maximum absolute cap is therefore O(N^(3/2)); no spectral truncation
or probabilistic exceptional event is being used.

## 3. A leading-order failure of the integrated width payment

For an order-d signing C define the width log pressure

    L_C(beta)=(1/2)[log sum_x exp(beta H_C(x)/sqrt(d))
                       +log sum_x exp(-beta H_C(x)/sqrt(d))].

Split A into two principal children by assigning k/2 entire cliques to
each child. Each child has d=N/2 vertices. Equations (2)--(3), and the
elementary log-sum upper bound, imply

    L_A(beta) <= N log 2 + beta N (8+Loff/sqrt(N))/4.   (4)

For either child, the coherent spin configuration gives the positive
maximum at least (d ell-d)/2. The negative-branch partition is at least
2^d by Jensen, since the uniform-spin quadratic mean is zero. Summing
these lower bounds gives

    sum_children L_child(beta)
      >= (N/2)log 2
           + (beta N/4)(8 sqrt(2)-sqrt(2/N)).          (5)

Consequently at beta=1,

    L_A(1)-sum_children L_child(1)
      <= N[ (log 2)/2
          +(sqrt(2)(1+1/sqrt(k))-8(sqrt(2)-1)
                                     +sqrt(2/N))/4 ]
       = (-0.1283001438...+o(1)) N.                   (6)

Now decompose the left side exactly as G_A-T_A: G_A is the full
cross-bridge insertion gain at the parent normalization, and T_A is
the reheating loss of the fixed actual principal children. G_A has the
exact integrated Gibbs-variance expression, including ALL Hubbard
degree-one, degree-two and latent-mean variance terms. Thus (6) is a
counterexample to the optimizer-independent integrated inequality

    G_A >= T_A-o(N).

Every matrix along the fixed bridge-deletion path is uniformly
cap-bounded: deletion is an average of two switching conjugates, and
the path is their convex interpolation. The root's delocalized quadratic
variance theorem therefore applies uniformly along the entire path.
Its positive gain is genuine, but need not pay reheating.

## 4. Scope

This construction contains mesoscopic ferromagnetic cliques and is not
claimed to minimize the homogeneous or anisotropic pressure. Its purpose
is to locate exactly the indispensable extra hypothesis. An actual
globally minimizing-signing theorem could exclude this coherent geometry;
the full Gaussian product-mixture variance identity and a uniform cap
bound, by themselves, cannot do so. This is an integrated asymptotic
obstruction in actual full signs, not a finite-size conditional-PSD test
or a counterexample obtained by allowing large edge magnitudes.
