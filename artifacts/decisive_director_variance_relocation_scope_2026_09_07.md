# Variance relocation: two distinct obstructions and the remaining question

Status: elementary counterexamples and a clearly labelled open comparison.
This is not a convergence theorem or a universal impossibility theorem.

The normalized direct-sum endpoint for two children is

    B=diag(sqrt(N/m) A_m, sqrt(N/n) A_n), N=m+n.

Its normalized cap is at most the size-weighted average of the child
normalized caps (with equality if both children are globally oriented to
attain their absolute cap positively). Its squared coefficient mass is
(N/2)(N-2), within O(N) of a flat
parent, and each row squared norm is N-O(1) at comparable splits. Its
coefficients exceed one internally and vanish across the split.

Thus the needed operation is NOT unbiased rounding inside the coefficient
cube. The independent soft-flatness proof establishes that no signing can
approximate this endpoint in the Boolean difference norm to o(N^(3/2)):
the missing flat cross block alone forces a leading error. This does not
exclude a favorable comparison of the two TOTAL caps.

## Total squared mass alone is inadequate, even at bounded amplitude

Fix L=2, take k=floor(N/L), and put B=L A_k on that principal block and
zero elsewhere, with A_k chosen from the proved all-order upper construction.
Then sum_(i<j) B_ij^2=(1+o(1))N^2/2, while

    Q(B)/N^(3/2) <= (0.494515125+o(1))/sqrt(2) <0.350.

Every flat signing has asymptotic normalized cap at least 0.4333221116640807.
Hence no general variance-mass-only relocation theorem can preserve total
cap at vanishing error, even when all coefficients have magnitude <=2.
Rows outside the block have variance zero: this does NOT falsify a theorem
requiring balanced row variance.

## Balanced row variance alone is inadequate without amplitude control

Use an order-14 signing with Q=21, whose upper witness is part of the exact
finite archive. For N=14t, take t disjoint copies, each multiplied by
sqrt((N-1)/13). Every row has squared norm exactly N-1. The triangle
inequality (no assumption about matching positive/negative extrema) gives

    Q(B)/N^(3/2)
       <= 21 sqrt((N-1)/N)/(14 sqrt(13)) <0.417.

This is again strictly below the verified universal flat-sign lower bound.
The amplitudes now grow as sqrt(N), so this does NOT falsify balanced-row
relocation with a fixed amplitude bound. Only the existence of the finite
upper witness, not its global optimality, is needed for this example.

## Open, not selected as a proved reduction

One could ask for a cap-preserving flat replacement for bounded-amplitude,
asymptotically row-regular matrices. The two-child endpoint lies in that
class for comparable splits. No such theorem is proved here. It is much
stronger than simply smoothing signs, and its target-order conclusion must
not be assumed in a variational-limit argument. Even a comparison with an
unspecified o(N^(3/2)) error would require checking the accumulated error
before asserting a convergence theorem.
