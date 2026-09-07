# Absolute-pressure reveal generator and the unpaid selection loss

Date: 2026-09-07. This is an original-cap-preserving version of the
partition-reveal diagnostic. No asymptotic drift sign is claimed.

## Joint-orientation generator

Let F_A(v)=log sum_x cosh(beta sum_e A_e sqrt(v_e)x_i x_j/sqrt(N-1))
and F_*=min_A F_A. In the joint Gibbs law on (sigma,x), set
Y_e=sigma A_e x_i x_j, d_e=E Y_e, K_ef=Cov(Y_e,Y_f). Then

    partial_e F_A = beta d_e/[2sqrt(N-1)sqrt(v_e)],
    partial_ef F_A = beta^2 K_ef/[4(N-1)sqrt(v_e v_f)]
                    -1_(e=f) beta d_e/[4sqrt(N-1)v_e^(3/2)].

For a continuous variance martingale dv_e=sum_j Gamma_ej dW_j, define
L_j=sum_e Gamma_ej Y_e/sqrt(v_e). The smooth active-branch Ito drift is

    beta^2/[8(N-1)] sum_j Var(L_j)
       - beta/[8sqrt(N-1)] sum_e d_e ||Gamma_e||^2/v_e^(3/2).

This is the joint-law covariance, not the average of the two separate
orientation covariances used for paired width. The minimum envelope adds
nonpositive switching local time. Neither term may be silently omitted.

For a discrete posterior step v -> V', choose any current optimal A.
There is an exact decomposition

    E F_*(V')-F_*(v)
      = [E F_A(V')-F_A(v)] - E[F_A(V')-F_*(V')].

The last quantity is nonnegative. Any lower-Jensen/submartingale argument
must pay this optimizer-selection loss as well as control frozen drift.
Successor profiles must keep the same vertex labels for this identity.

## Full-sign numerical audit

`computations/flatify_adversary_2026_09_07_absolute_reveal_drift.py` exhausts
all signings modulo switching for N=4,6,7, at beta=.7,1.5,3. It constructs
the true hidden-partition posterior profiles with common labels, verifies
their martingale identity, computes optimized and frozen pressures, and
checks telescoping. These are floating-point diagnostics, not certificates.
All output is saved in the corresponding computations/results JSON.

At beta=3 the cumulative results are:

| N | optimized endpoint drift | frozen drift sum | selection-loss sum |
|---|---:|---:|---:|
| 4 | -0.2927264031 | -0.2927264031 | 0 |
| 6 | 1.6698663924 | 2.6290430442 | 0.9591766518 |
| 7 | 1.6431163246 | 4.9759594915 | 3.3328431668 |

Both conditional drift signs occur for every listed N at beta=3. The N=4
counterexample survives the switch from paired width to original absolute
pressure, without needing optimizer switching. At N=7 the selection loss
is a substantial fraction of frozen positive drift. Neither observation
is an asymptotic obstruction on actual optimal children.

## An exact normalization that removes the independent-edge baseline

At an active optimum, the edge-flip ratio implies d_e<=tanh(lambda_e),
lambda_e=beta sqrt(v_e)/sqrt(N-1). Thus

    R(v)=F_*(v)-sum_e log cosh(lambda_e)

is coordinatewise nonincreasing, in the integrated one-sided sense,
including zero-coordinate limits. This is distinct from claiming F_*
itself is coordinatewise monotone. It is adjacent to the archived
logcosh-interpolation method, not a replacement for its missing sign.

For any row-regular profile with max_e v_e<=C_delta,

    beta^2 N/4 - O_delta(beta^4)
       <= sum_e log cosh(lambda_e) <= beta^2 N/4.

Indeed x^2/2-x^4/12<=log cosh x<=x^2/2 for all real x, and
sum_e v_e^2<=C_delta N(N-1)/2. Therefore the independent-edge baseline
changes by only O_delta(beta^4) along the reveal profiles. At beta=N^(1/4)
this fits within the newly proposed beta^2 sqrt(N) pressure error budget.
This does NOT control the residual R or its optimizer-selection loss.

The high-temperature expansion makes that residual concrete:

    exp(R(v)-Nlog2)
      = min_A sum_(G Eulerian, |G| even) product_(e in G) A_e tanh(lambda_e),

including the empty graph. The orientation average removes graphs with
odd edge count. Thus the unaccounted part is a globally optimized sum of
even Eulerian correlations, beginning with four-cycles; it is not merely
independent-edge variance entropy. No sign inequality for its martingale
drift has been established.

Finally, equal-child terminal absolute optimization equals twice optimized
paired-width pressure, as proved in the companion reveal-count audit.
Keeping absolute pressure at intermediate states is necessary to preserve
the parent cap, but does not remove this stronger endpoint target.
