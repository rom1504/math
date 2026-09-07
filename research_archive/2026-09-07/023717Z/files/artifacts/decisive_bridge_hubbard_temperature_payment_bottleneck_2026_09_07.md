# What full Hubbard information must supply to pay the child-temperature loss

This records an exact sufficient comparison, not an assumption that the new
quadratic fluctuation theorem proves it. The positive bridge-pressure theorem
is proved separately. Here the desired normalization payment remains open.

## 1. The full variance decomposition

Use the Gaussian augmentation K_s=D+sJ, and condition on (g,s,Y). Put
m_i=tanhY_i and v_i=1−m_i². For a hollow quadratic F=Σ_{i<j}b_ij xi xj,
the product expansion gives the exact identity

    E_g Var_(μ_g)(sF)
      =E[Σ_{i<j} b_ij² v_i v_j
           +Σ_i v_i(Σ_{j≠i}b_ij m_j)²]
         +E_g Var_(s,Y|g)(sΣ_{i<j}b_ij m_i m_j).       (1)

The first term is the degree-two product chaos retained in the root's
delocalized variance bound. The second is its degree-one chaos; the third
is fluctuation of the conditional mean. These two additional terms are
nonnegative, but depend on how the actual bridge acts on the latent mean
vectors, not only on coordinatewise conditional variances.

An exact energy identity is

    E[Y_i tanhY_i]=d_i+E[s xi(Jx)_i],
    2E[sH_J]=Σ_i(E[Y_i tanhY_i]−d_i).                (2)

Indeed Y_i has conditional mean (D+sJ)x_i in the Gaussian channel, while
E[x_i|s,Y,g]=tanhY_i. There is no independence assumption on the Y_i.

For a fixed branch or prescribed mixture of branches, the corresponding
entropy identity is H(X|s)=Σ_i E h₂((1+m_i)/2)+I(X;Y|s). Thus conditional
product entropy does not equal the full spin entropy; dropping its mutual
information loses precisely another potentially extensive quantity.

Equations (1) and (2) do not currently imply that the bridge variance pays
the temperature loss. Such a conclusion needs a lower bound on the
bridge-alignment or latent-variance terms in (1), or a sharper comparison
using their joint law. The cap/diagonal-majorant estimates give a positive
constant factor but not the required matching coefficient.

## 2. Exact integral payment condition for the width objective

Let N=m+n, r=m/N, and let A be an actual parent signing. In this section
Z_A^±(t)=Σ_x exp(±tH_A(x)) uses the RAW temperature t. Define

    L_A(β)=[logZ_A^+(β/√N)+logZ_A^−(β/√N)]/2,

where the displayed arguments here denote raw coupling temperatures.
Let A_1,A_2 be its principal children. With B its cross bridge scaled by
β/√N, define the fixed-parent bridge gain

    G_A = L_A(β)
          −L_(A_1)(β√r)−L_(A_2)(β√(1−r)).

In this line the child L uses its own order normalization, so both children
at the subtraction point have raw temperature β/√N. The exact formula is

    G_A=1/2 Σ_{s=±1} ∫_0^1(1−t)
                       Var_(s,J_0+tB)(H_B) dt.       (3)

Define the fixed-child reheating loss

    T_A=Σ_{i=1,2}[L_(A_i)(β)−L_(A_i)(β√r_i)].      (4)

It is an integral of the child energy expectation (2), with the corresponding
temperature derivative. A sufficient fixed-parent theorem for a same-β width
comparison would be: for every large N, some Ψ_N(β)-minimizer A and a
comparable split satisfy

    G_A ≥ T_A−C_β N^α,       α<1.                    (5)

Then Ψ_N(β)≥Ψ_m(β)+Ψ_n(β)−C_βN^α. This is a one-sided comparison involving
actual globally minimizing parent signings. It is stronger than the version
using only optimized child reheating losses, but is noncircular and explicit.

The root's pressure theorem proves G_A≥κ_β,r β²mn/(2N)>0. It does not prove
(5), since κ_β,r is not shown to match (4). Formula (1) states exactly which
discarded positive quantities could in principle strengthen that estimate.

## 3. Original absolute objective and the rate requirement

Even a proof of (5) for width pressure is not automatically an original
absolute-cap result. At a block-diagonal endpoint the absolute partition is
Z_1^+Z_2^+ + Z_1^-Z_2^-, rather than the product of the child absolute
partitions. The width average avoids this polarity obstruction by changing
the objective; it must not subsequently be identified with the original
absolute minimum without an additional theorem.

Likewise, an unspecified o(N) comparable-block error does not by itself
imply convergence. For example f(N)=N sin(log log(N+2)) has o(N) two-block
additive defects but f(N)/N does not converge. A power-saving error as in
(5), or another summable almost-additivity modulus, is the appropriate
target. This is why the exact scale of any latent-field comparison matters.

## 4. A separate squared-degree profile route would need a new theorem

Let w_ij be bounded nonnegative magnitudes with every squared row sum n−1.
A useful sufficient profile statement would compare the globally optimized
homogeneous pressure from above with the optimized block-normalized profile,
to error C_β n^α, α<1, uniformly for comparable block sizes. The block
endpoint then provides an almost-subadditive upper comparison for the
original absolute partition, since its common-polarity sum is at most the
product of the child absolute sums.

This statement is not a consequence of convexity of logZ in its couplings.
Within one signing branch, differentiation with respect to a squared
magnitude direction d gives

    d²/du² logZ(λ_e(u)=β√(v_e+u d_e)A_e/√n)
      =Var(Σ_e λ'_e χ_e)+Σ_e λ''_e Eχ_e,

and the second term can be negative. More importantly, taking the minimum
over signings adds downward derivative jumps at branch crossings. The new
positive quadratic variance theorem controls neither those jumps nor the
negative normalization term at its required coefficient. A proof of this
profile comparison must address both explicitly.
