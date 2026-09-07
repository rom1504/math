# Scalable exceptional clouds inside the actual rank-two sign ensemble

Date: 2026-09-07. Exact witness, exact architecture-preserving inflation,
and an exact lift to actual optimizing children. This disproves dropping
the profile restriction from a universal all-cloud .74 bound. It does NOT
disprove the existence of a favorable bridge for fixed children.

## 1. Independently checkable finite witness

The solver-free script
`computations/flatify_independent_2026_09_07_exceptional_cloud_exact.py`
and matching results JSON contain four left and four right order-eight
Hadamard bases, their literal rank-two tiles, and four Boolean vectors.
All arithmetic in the witness verification is integer arithmetic. They give

    C in {+1,-1}^{32 by 32},  C C^T=32 I,
    x0^T C y0=0,             x^T C y=152,
    d_H(x,x0)=d_H(y,y0)=3.                               (1)

The normalized target is 152/(128sqrt(2))=.839689302659... . The checker
rebuilds C from every tile (1/2)F_i,pair H_2 G_j,pair^T, rather than merely
checking orthogonality of an arbitrary sign matrix.

Two additional exact identities use zero-based coordinates a=9,b=2:

    x_a=x0_a,  y_b=y0_b,
    x_a(Cy)_a=0,  y_b(C^T x)_b=2,  x_a C_ab y_b=1.       (2)

The discovery pilot and its unsuccessful/random-center cases are retained
separately in the local-cloud and adversarial-cloud pilot scripts/results.
The proof uses only the exact extracted witness and checker.

## 2. Inflation stays inside the SAME rank-two architecture

For any Hadamard H of order q, define an order-q^2 matrix

    R[(s,t),(u,v)]=H[t,u] H[v,s],
    w[(s,t)]=H[t,s].

Direct orthogonality sums give

    R R^T=q^2 I,     Rw=q w,     R^T=R.                  (3)

This is not merely an arbitrary tensor enlargement. In the new rank-two
ensemble use 4q fibres of size 8q on each side. The basis in new left
fibre (i,s) is F_i tensor H, and in right fibre (j,u) is G_j tensor H.
Assign the pair of columns (2j,u),(2j+1,u) to opposite fibre (j,u), and
reciprocally the pair (2i,s),(2i+1,s). The resulting tile entry is exactly

    C_ij[a,b] H[t,u] H[v,s].

After reindexing ((i,s),(a,t)) as ((i,a),(s,t)), the full sign bridge is
C tensor R. Thus the inflated matrix is a literal output in the same
rank-two construction with k'=8q, m'=4q, n=32q^2. The checker constructs
all these new tiles explicitly for q=2 and q=4 and verifies their equality
to the tensor matrix under the stated permutation. The displayed formula
proves it for every q.

Tensoring all four spins in (1) with w gives center energy zero, target
152q^3, and distances 3q^2 on both sides.

## 3. Increase the distances to floor(n/10) without losing the counter

Put r=floor(q^2/5). In the tensor target x tensor w, flip r coordinates
within base coordinate a from (2); in y tensor w flip the SAME subset S
of r microcoordinates within base coordinate b. These coordinates were
agreements with their respective centers. Therefore the new distances are

    3q^2+r=floor(n/10).                                 (4)

Choose S with w_S^T R[S,S] w_S>=0. Such a subset exists for every r:
the switched matrix diag(w)Rdiag(w) has diagonal entries one and every
row sum q, so a uniform r-subset has expected quadratic sum

    r+[r(r-1)/(q^2(q^2-1))](q^3-q^2)>=0                 (5)

for q>1. The cases r=0 and q=1 are immediate. Conditional expectation,
or finite exhaustive selection, gives a deterministic choice if desired.

The first target's linear energy change is zero by (2). The second costs
exactly 4qr. Their bilinear correction is 4w_S^T R[S,S]w_S>=0, again by
(2). The resulting target energy is consequently at least

    152q^3-4qr >=(756/5)q^3.

The center was unchanged. Hence, along q=2^j for example, actual rank-two
sign bridges and actual Boolean center/target pairs satisfy

    center energy =0,
    both target distances =floor(n/10),
    target energy/n^(3/2) >= b_*:=189/(160sqrt(2))>.83.   (6)

The directed interval checker certifies b_*>83/100. This is a scalable
counterexample to an unrestricted .74 bound, including asymptotically
exact .1 noise level. It does not rely on extrapolating finite numerics.

## 4. The exceptional pair can be aligned with ACTUAL optimal children

For each inflated order n choose any actual minimizer A with Q(A)=M_n.
Orient its global sign so some ground spin attains +M_n, then switch
that ground spin to the prescribed bridge center x0. The analogous
operation is allowed independently for the right child and center y0.
These operations preserve actual optimality.

In the ground-switched coordinates, randomly permute the child matrix.
The prescribed target relative to its center has r_n=floor(n/10) negative
coordinates. Under a uniform permutation this is a uniform Hamming sphere.
For any distinct coordinate pair the exact sign-product expectation is

    kappa_n=((n-2r_n)^2-n)/(n(n-1)).                      (7)

Consequently the expected target child energy is kappa_n M_n. Some
permutation attains at least this mean. Choose such permutations for the
two children independently, and undo their respective ground switches.
The prescribed centers remain their actual maximizing spins and the
prescribed targets have child energies at least kappa_n M_n each.

Together with the fixed positive bridge energy in (6), this gives actual
optimal children A_n,D_n and actual bridges C_n for which

    Q([[A_n,C_n],[C_n^T,D_n]])
       >= b_* n^(3/2)+2 kappa_n M_n.                    (8)

Since kappa_n->16/25 and M_n/n^(3/2)<=U+o(1), U=.493608094,

    liminf [Q(parent)-2sqrt(2) M_n]/n^(3/2)
       >=b_*-(2sqrt(2)-32/25)U > .07.                   (9)

The last strict margin is also checked by directed arithmetic. No
unproved convergence of M_n was used. The child permutations are chosen
by finite averaging; no approximate optimizer or finite optimum table
is substituted for the actual minimizers in this argument.

## 5. Exact logical scope

This proves that even child optimality does not make all pairs belong to
the controlled noisy-profile sector, or make a uniformly centered bridge
automatically satisfy a constant-preserving parent bound. There are bad
alignments at every inflated scale, with an explicit positive excess.

The quantifiers are existential in BOTH the bad bridge and its child
alignment. We have NOT shown that every bridge is bad for a fixed pair
of children, that a selected bridge cannot avoid this pair, or that the
exceptional set has enough entropy to defeat a more refined construction.
Thus (9) is not a nonconvergence theorem and not an obstruction to the
desired favorable all-seed composition operation.
