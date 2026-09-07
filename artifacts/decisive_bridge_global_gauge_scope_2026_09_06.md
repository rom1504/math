# Why signed-transpose conjugacy does not alone erase joint sources

The exact algebraic observation is due to the transfer-seeds agent. For a
symmetric sign seed A, choose a coordinate diagonal gauge D with
D_ij D_ji=A_ij for i!=j. Then D R_A D is ordinary transpose off the
diagonal. Fixing the seed diagonal fixes the remaining action. This DOES
erase the seed for any source class invariant under all such coordinate
gauges. Product symmetric coordinate sources are an important such class.

It need not preserve joint Hadamard-spectrum sources. Here is an exact small
test, not an all-order conclusion. In the Sylvester H4 basis the sign vector
(1,1,1,-1) has spectrum h=(2,2,2,-2). The coordinate multiplier
(-1,1,1,1) sends this to h'=(-2,2,2,-2), whose inverse transform is
(0,0,0,-2), not a sign vector. This multiplier occurs in one row of the
triangular coordinate gauge for a seed with a single negative offdiagonal
edge (choose the lower-triangular entry negative).

The example survives allowing a COMMON column-sign multiplier for the two
rows h,h': each is flat and their spectral inner product is 8. If both were
Hadamard transforms of sign vectors for a COMMON H4, those sign vectors
would be flat-spectrum sign vectors with inner product 2. Direct enumeration
gives possible inner products only -4,0,4. Every order-four Hadamard is
equivalent to H4, so common signed row/column permutations do not repair this.
The short exact check is in the companion computations file.

This is deliberately a narrow scope result. Tensoring a fixed basis extends
the elementary failure of its cube automorphism property. It does NOT prove
that arbitrary larger randomized Hadamard ensembles have the same support
obstruction, nor that the final coherent-state supremum retains a quantitative
seed advantage. Those require a separate proof. The global truncation theorem
keeps precisely the joint source polynomials where this question lives.
