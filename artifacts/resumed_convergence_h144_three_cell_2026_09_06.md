# An exact three-point reflection inside the prescribed H144 catalyst

Date: 2026-09-06. Status: explicit integer certificate and exact norm
consequence. This settles the finite invariant-partition question left
open in `fresh_limit_variational_2026_09_05.md`; it does not establish
`R=T` for general seeds or convergence of the original minima.

## 1. Construction and independently executable certificate

Let `F` be the Paley order-12 Hadamard indexed by infinity and F_11:

    F_(infinity,j)=1,
    F_(i,infinity)=-1 for finite i,
    F_(i,i)=1,
    F_(i,j)=quadratic_character(j-i) for distinct finite i,j.

Thus `F^T F=12I`. The prescribed regular symmetric square Hadamard has

    H_((i,j),(l,k))=F_ij F_ik F_lk F_lj.

Its normalized action, on 12-by-12 arrays P, is

    U(P)=F circ [F (F circ P)^T F/12],       U=H/12.              (1)

It has `U^2=I` and `U(1)=1`. Row and column sign changes of F cancel
out of its four-entry product, and row/column permutations merely
reindex H, so this description is invariant under the usual equivalences
of the underlying order-12 Hadamard.

The baseline generator uses the opposite finite-field orientation
`chi(i-j)`. If that baseline matrix is `F_old` and the displayed one is
`F_new`, then `F_new=D F_old^T D` with
`D=diag(1,-1,...,-1)`. The diagonal gauges cancel in H, while replacing
F by its transpose conjugates H by the coordinate swap `(i,j)->(j,i)`.
Thus this certificate is for exactly the prescribed H144 up to a
coordinate permutation; it is not a change of Hadamard generator.

Here is a three-color partition of its 144 coordinates (row order is
infinity,0,...,10 and likewise for columns):

    0 1 2 1 1 2 2 2 1 0 0 0
    0 2 0 2 0 1 1 2 2 1 0 1
    2 2 0 1 0 1 0 1 2 2 1 0
    2 2 1 0 0 1 0 0 2 2 1 1
    2 0 0 1 1 2 2 0 0 1 1 2
    0 1 2 1 0 2 2 2 0 1 0 1
    2 0 0 0 2 2 1 1 1 2 1 0
    1 0 1 2 2 1 0 0 1 0 2 2
    1 0 2 0 2 0 1 1 0 2 2 1
    1 1 1 2 2 0 1 0 0 0 2 2
    1 1 2 0 1 0 0 2 2 1 2 0
    0 2 1 2 1 0 2 1 1 0 0 2

For its indicator masks `P_c` and `Z_c=F circ P_c`, direct integer
multiplication gives

    Z_c^T F+F^T Z_c=8I,       |P_c|=48,       c=0,1,2.          (2)

An exact standalone verifier, requiring no solver, is
`computations/resumed_convergence_h144_partition_verify_2026_09_06.py`.
It constructs F and H, verifies (2), symmetry, `H^2=144I`, regularity,
the quotient identity below, and the exact energy in Section 2.
The search script in `tmp/resumed_convergence_h144_partition.py` used
CP-SAT only to discover the coloring. Its solver status is not used as
the mathematical certificate.

Multiply (2) on the left by F. Equation (1) immediately yields

    U P_c=(2/3)1-P_c.                                          (3)

Thus the space of functions constant on the three cells is an invariant
probability algebra, and its induced orthogonal action is exactly

    Q3=2E3-I3,

where E3 is averaging for the uniform measure on three points. This is
a genuine non-Walsh alphabet inside the existing allowed catalyst; no
new outer matrix has been added to the original regularization.

## 2. Exact landing of a non-Walsh full sign seed

Take

    B=J3-2I3,

whose eigenvalues are 1,-2,-2. Its absolute value is `2I3-J3/3`.
This matrix is a feasible absolute PSD majorant, and for Boolean x

    x^T |B| x=6-(sum_i x_i)^2/3 <=17/3.

Therefore `T(B)<=17/6`. Conversely the cut covariance `C=B^2/3`,
realized by the three columns of B with equal probability, gives

    2T(B)>=Tr|B|^3/3=(1+8+8)/3=17/3.

So `T(B)=17/6` exactly.

Let `X_(a,i)=B_(color(a),i)` be the Boolean array on the H144 coordinates
and the three seed coordinates. Its columns are functions in the
invariant algebra (3). Consequently

    X^T U X/144=B Q3 B/3,

and the normalized same-spin quadratic energy is

    X_vec^T (H tensor B) X_vec / 144^(3/2)
        =Tr(B^3 Q3)/3=17/3.                                   (4)

In exact integers the unhalved energy is 9792 and the halved energy is
4896. Hence already at this single outer order

    q(H144 tensor B)/144^(3/2)=17/6,

where the upper inequality follows from `R<=T`. In particular

    R(J3-2I3)=T(J3-2I3)=17/6.                                  (5)

This is a finite exact landing, not merely a limiting witness. It is
not a useful small original-minimum upper construction: dividing (5)
by `3^(3/2)` gives about 0.5453, above 1/2. Its positive contribution
is a coherent non-Walsh realization step, extending the previously
known scalar/H2 landing and exposing a new three-point reflection
operation available for further catalyst constructions.

## 3. Tensor consequences, without an overclaim

Tensoring the partition realizes `Q3^(tensor t)` inside `U144^(tensor
t)` for every t. Combining it with the four-point outer realizes
products of the uniform three-point and four-point averaging
reflections. Boolean vector-valued functions on these finite product
probability spaces can therefore be used directly as exact witnesses
in the prescribed R norm.

The raw tensor-product partition does not itself give the single
reflection `2E_(3^t)-I`: it assigns signs by the parity of nonconstant
factors. The following additional coding construction does overcome
that particular obstacle.

## 4. Simplex compression gives every ternary-power alphabet

For each integer r>=1, let `t=(3^r-1)/2` and form an r-by-t matrix L
over F_3 with one nonzero column representative from each projective
line in F_3^r. The linear map `L:F_3^t -> F_3^r` is onto, so every
fiber has the same size. For every nonzero `a in F_3^r`, the word
`a^T L` has precisely `3^(r-1)` nonzero entries: the number of all
projective lines minus the number in the kernel hyperplane is

    (3^r-1)/2 - (3^(r-1)-1)/2 = 3^(r-1).

This number is ODD. On the product three-point space, `Q3^(tensor t)`
multiplies a character by (-1) to its number of nonconstant coordinate
characters. Therefore on EVERY nonconstant character pulled back from
F_3^r along L its eigenvalue is -1; its eigenvalue on constants is +1.
Characters span all functions on this quotient. The algebra of
functions of L is consequently invariant, and its exact action is

    Q_(3^r)=2E_(3^r)-I.                                        (6)

Combining (6) with (3), the prescribed outer `H144^(tensor t)` contains
an exact invariant uniform probability algebra of size `3^r` with the
single averaging-reflection action (6). This is a quotient by a linear
map, not an assertion that raw product parity equals global parity.
The proof is elementary character counting and needs no coding-theory
existence theorem.

## 5. Arbitrary finite probability reflections are limiting tests

The averaging reflection preserves every subalgebra defined by a
partition of its atoms. Partition the `3^r` quotient atoms into k
groups of sizes `m_1,...,m_k`. The induced probability vector is
`p_j=m_j/3^r`, and the induced action on functions on k atoms is

    (Q_p f)(j)=2 sum_l p_l f(l)-f(j).                            (7)

For any fixed probability vector p, choose integer multiplicities
with `m_j/3^r -> p_j`. For fixed finite Boolean vector labels and a
fixed seed B, every bilinear or quadratic objective of this test is
continuous in p. Thus R includes all finite probability-space
averaging reflections as limiting tests, and by tensoring includes
every finite product of such reflections. No prescribed generator was
changed or enlarged.

For example a single such reflection gives the concrete lower test

    R(B) >= (1/2) sup_(law of X in {+/-1}^k)
                   |2(EX)^T B(EX)-E[X^T B X]|.                 (8)

More generally choose Boolean vector functions on any finite product
probability space and apply the product of its averaging reflections.
The resulting objective is a rigorous lower bound on R(B). Formula
(8) alone is NOT the whole R norm; for the H2 seed its optimum is
5/4, below the known limiting value sqrt(2), so tensor products can
be essential.

This realizes arbitrary finite atom weights but not arbitrary
orthogonal operators on those atoms. In particular it still does not
prove the optimizing polar coupling in the T dual can be realized.

## 6. A second full-sign landing using a permissible one-sided permutation

The same three-cell quotient settles the distinct seed

    B2 = [-1  1  1; 1 -1  1; 1  1  1].

Its eigenvalues are 2,-1,-2. Put `v=(1,1,-1)^T`. Since
`|B2|=2I-vv^T/3`, the cube-majorant argument and the column covariance
`B2^2/3` again prove `T(B2)=17/6`.

Define the Boolean row-label matrix and a permutation matrix by

    F0 = [1  1  1; 1 -1 -1; -1  1 -1],
    P e_1=e_1,   P e_2=e_3,   P e_3=e_2.

Direct 3-by-3 multiplication gives

    (1/3) Tr[B2 F0^T (-Q3 P) F0] = 17/3.                       (9)

Choose the input Boolean labels `F0` on the three cells and the output
labels `-P F0`. Their normalized bilinear objective with outer H144 is
the left side of (9). A permutation of the labels on ONE side of a
bilinear witness is always allowed. This is not a composition lemma
for overlapping invariant quotient operators.

To convert the bilinear witness x,y for `C=H144 tensor B2` into an
allowed same-spin witness, use `R4=J4-2I4` and

    z=(x,y,y,x),
    z^T (R4 tensor C) z=8 x^T C y.

After the outer normalization this gives the required half-bilinear
value. Hence

    R(B2)=T(B2)=17/6,                                          (10)

with a finite same-spin certificate at outer order `4*144=576`.
Together with switching, permutations, and global negation, Sections
2 and 6 cover two different indefinite order-three full-sign classes.
The class `J3-2e_1e_1^T` is NOT claimed here: its T value is
`5/sqrt(2)` and its two repeated coordinates reduce it to the weighted
two-by-two seed `[-1,2;2,4]`, an additional realization problem.
For the stated T value, the majorant is the block matrix with first
entry sqrt(2), zero first-to-last cross entries, and last 2-by-2 block
sqrt(2) times J2. It has cube maximum `5sqrt(2)` and dominates both
signs of the seed. The matching dual uses spins `(X,Y,Y)` with
`E[XY]=3/4`; its compressed nuclear norm is `5sqrt(2)`.
