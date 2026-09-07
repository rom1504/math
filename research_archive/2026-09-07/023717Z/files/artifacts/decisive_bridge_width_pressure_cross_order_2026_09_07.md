# A genuine cross-order lower bound for optimized width pressure

This applies the cap-controlled fluctuation theorem in
decisive_director_delocalized_quadratic_fluctuations_2026_09_07.md.
I independently checked its Gaussian mixture, product-chaos covariance,
coefficient-delocalization normalization, and quenched-variance distinction.
The present recurrence has contracted child temperatures and does not prove
convergence of the absolute cap.

## Definitions and a uniform actual-cap bound

Let

    Z_A^±(β)=Σ_{x∈{±1}^N} exp[±β H_A(x)/√N],
    Ψ_N(β)=min_A [log Z_A^+(β)+log Z_A^−(β)]/2.

Ψ is an **unnormalized** log pressure. At zero temperature its leading
coefficient is the optimized half-width (max H−min H)/2, not necessarily
the optimized absolute cap.

For any signing A and J=βA/√N,

    [log Z_A^++log Z_A^−]/2
       ≥[max H_J−min H_J]/2 ≥Q(J)/2,

because H_J has uniform-spin mean zero. On the other hand averaging all
independent edge signs and applying Jensen yields

    Ψ_N(β)≤N log2+binom(N,2)log cosh(β/√N)
           ≤N log2+β²(N−1)/4.

Thus every Ψ-minimizer satisfies

    Q(J)≤C_β N,       C_β=2log2+β²/2.                 (1)

Optionally, the all-order random-sign cap candidate Q(A)≤√log2 N^(3/2)
also gives Ψ_N≤Nlog2+β√log2 N. Thus one may replace C_β everywhere by

    C_β=min{2log2+β²/2, 2log2+2β√log2}.

The minimizer remains deterministic; no expectation-adaptive signing was used.

## Fixed-parent bridge comparison

Split N=m+n, put r=m/N, and let J_0 be the block-diagonal restriction of
J, B=J−J_0. Switching one whole block gives

    J_0=(J+D_S J D_S)/2,

so convexity and switching invariance of Q imply Q(J_0)≤Q(J). The whole
deletion path J_t=(1−t)J_0+tJ therefore obeys Q(J_t)≤C_βN.

For this flat bridge,

    V_B=β²mn/N,       η_B=N/min(m,n).

Define the explicit positive constant

    c(C,η)=exp{−2η[min(4C,2C+1)+4K_G C+√(8K_G C/π)]}.

The root fluctuation theorem applies separately to either fixed energy
orientation. Along each branch the bridge derivative vanishes at t=0 by
block spin switching, and the second derivative is at least c(C_β,η_B)V_B.
Hence

    log Z_J^± ≥ log Z_(J_0)^± + c(C_β,η_B)β²mn/(2N). (2)

This inequality holds for each fixed bounded-cap parent before any
minimization, so changing optimizers does not enter the integration.

## Exact optimized recurrence

At the block-diagonal endpoint, each fixed branch factors. Averaging its
two logarithms therefore yields the sum of the two child width pressures,
with temperatures β√r and β√(1−r). Applying (2) to a Ψ_N-minimizer and then
minimizing its child restrictions proves

    Ψ_N(β) ≥ Ψ_m(β√r)+Ψ_n(β√(1−r))
                 +[c(C_β,η_B)/2]β²mn/N.             (3)

Equivalently, for ψ_N=Ψ_N/N,

    ψ_N(β) ≥ r ψ_m(β√r)+(1−r)ψ_n(β√(1−r))
                 +[c(C_β,η_B)/2]β²r(1−r).

There is no common-polarity factorization error here: averaging the two
branch logarithms commutes exactly with summing the child logarithms.
For the original absolute partition, replacing its log-sum by this average
is a relaxation and must not be reversed.

The positive coefficient in (3) has not been shown to compensate the
contracted temperatures. Thus this is a genuine quantitative cross-order
inequality for an actual-signing optimization, but not ordinary or
almost-superadditivity at fixed β, and not an original-limit theorem.
