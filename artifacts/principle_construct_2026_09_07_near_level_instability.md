# A deterministic near-level instability inequality

2026-09-07. **Proved.** Let A be any nonzero hollow symmetric matrix, not
necessarily a signing, and put

    H(x)=sum_{i<j} A_ij x_i x_j,   Q=max_x |H(x)|,
    I(x)=sum_i [-x_i(Ax)_i]_+.

Then every Boolean word satisfies

    I(x)^2 <= 4Q [Q-H(x)].

This is a deterministic link between the actual energy deficit and total
one-coordinate instability. It does not assert that small instability
forces near-optimal energy.

## Proof

Let U be the negative-field coordinates of x. Flip each coordinate in U
independently with probability p, leaving every other coordinate fixed.
Writing x_U for the partial word supported on U, exact expansion gives

    E[H(new)-H(x)] = 2p I(x)+4p^2 H(x_U)
                  >=2p I(x)-4p^2 Q.

The last step is principal cap monotonicity, obtained by independent
unbiased completion of a partial Boolean word. Every resulting full word
has energy at most Q, so the displayed lower bound cannot exceed Q-H(x).

The standard polarization inequality beta(A)<=4Q, together with
`2I(x)=||Ax||_1-2H(x)`, gives I(x)<=2Q-H(x)<=3Q. Thus
`p=I(x)/(4Q)` belongs to [0,1]. Substituting this value gives the claim.
The zero matrix case is immediate.

## Actual low-cap near-level consequence

If Q(A)<=c n^(3/2), every positive-polarity word of deficit at most
eta n^(3/2) has

    I_A(x) <= 2 sqrt(c eta) n^(3/2).

The corresponding statement for a negative-polarity word uses -A. Hence,
for eta<1/(8 pi c), the sharp universal lower-tail theorem for actual
bounded-cap signings applies to both parts of the absolute near-level:
its cardinality is at most `2^(n+1) exp(-delta(c,eta)n)` for some positive
delta. One should use strict cap envelopes and strict threshold margins
when applying the asymptotic lower-tail result.

This gives an energy-specific consequence of universal instability
concentration, without assuming a spectral bound on the full matrix.
The entropy deficit is only a fixed positive fraction; it does not
control the factorization norm of the near-level code, so it does not
verify the hypothesis of hierarchy dilution or solve convergence.

## Quantitative relative-gap count

I reread the diagonal-majorant and quantitative zero-field parts of
`principle_director_universal_stability_entropy_2026_09_07.md`.
Their exponential estimate gives the following explicit dependence.
Let c0=1/sqrt(2 pi), C0=max(C,1), and

    E_rel(rho)={x: |H(x)| >= (1-rho)Q},
    Delta=min(c0/2, c0-2C sqrt(rho)).

For rho<1/(8 pi C^2), Delta is positive. Universal constants K,c>0 give

    |E_rel(rho)| <= 2^(n+1) K exp(-c Delta^3 n/C0^2).

Indeed the positive part is contained in
`{I_A <=2Q sqrt(rho)}`, and the negative part uses -A. Each falls below
the universal instability mean by at least Delta in n^(3/2) units.
The harmless endpoint Delta=c0/2 is obtained by taking the limit in the
tail theorem or shrinking Delta by an arbitrary fixed factor. If one
prefers avoiding endpoints altogether, replace c0/2 here by c0/3.

For an absolute normalized deficit eta, the same formula holds with
`Delta=min(c0/2,c0-2sqrt(C eta))`. These are bounds on actual energy
superlevels, not on a Gibbs law after an unproved entropy-preserving
cleanup. The constants are uniform over full-sign A with the given cap
envelope; the relative-gap range is not asserted to be sharp.

## A sharper bound in thin energy windows

The same exact flip calculation also gives, with d=Q-H(x),

    I(x) <= min(2sqrt(Qd), d+sqrt(2Qd)).

For the second bound the case I<=d is immediate. Otherwise write
H(x_U)=-alpha. The full flip cannot exceed Q, so
`2I-4alpha<=d`; since I>d this forces alpha>0 and the optimizing
probability I/(4alpha) to be at most one. Hence the random-subset
inequality gives `I^2<=4alpha d`. On the other hand the FULL flipped
word has energy `H+2I-4alpha>=-Q`, giving
`4alpha<=2Q-d+2I`. Combining yields

    (I-d)^2 <=2Qd,

as claimed. No field condition on the complement is needed.
This improves the leading small-deficit coefficient from 2 to sqrt(2).
In the relative-gap count, `2C sqrt(rho)` can therefore be replaced by
`C min(2sqrt(rho),rho+sqrt(2rho))`. For an absolute normalized deficit,
replace `2sqrt(C eta)` by
`min(2sqrt(C eta),eta+sqrt(2C eta))`.
