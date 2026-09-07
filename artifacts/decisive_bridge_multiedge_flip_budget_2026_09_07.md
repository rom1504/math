# Actual multi-edge optimality and the log-cosh slack budget

This is a bounded follow-up audit of the optimized WIDTH interpolation, not
a convergence proof and not a comparison for the original absolute cap.

## 1. Exact constraints at one actual anisotropic minimizer

Use the notation of decisive_bridge_width_logcosh_interpolation_2026_09_07.md.
Fix one globally minimizing actual signing A at its fixed positive magnitudes
lambda. For an arbitrary edge set E write V_E(x)=sum_(e in E)lambda_e A_e x_e.
The simultaneous-flip cost is exactly

    Delta_E = (log <exp(-2V_E)>_+ + log <exp(2V_E)>_-)/2 >= 0. (1)

No reoptimization or change of partition is used here. For a vertex set T,
let cut(T) be all crossing edges. Gauge switching gives Delta_cut(T)=0.
If F and G partition cut(T), then A^F and A^G are gauge-equivalent and

    Delta_F=Delta_G.                                  (2)

In particular F may be the within-interpolation-block part of a cut and G
the cross-interpolation-block part. Equality (2) does not say their sums of
singleton costs agree: interactions among the flipped edges are essential.

## 2. A uniform, summably small singleton logarithm conversion

The exact nonnegative singleton slack is

    a_e=(exp(2Delta_e)-1)/sinh(2lambda_e)^2.

Along a comparable two-block log-cosh path at fixed beta, one can replace
exp(2Delta_e)-1 by 2Delta_e with integrated error O_beta(N^(3/2)) in the
unscaled signed slack sum, hence O_beta(sqrt(N)) in total pressure.

Here is a proof that does not assume typical spin correlations. Convexity
of the even radial width pressure gives sum_e lambda_e d_e>=0. Since
2d_e=tanh(2lambda_e)(1-P_e-a_e), positivity of a_e and |P_e|<=1 imply

    sum_e lambda_e tanh(2lambda_e) a_e
       <= 2 sum_e lambda_e tanh(2lambda_e)=O_beta(N).   (3)

Also 0<=Delta_e<=2lambda_e, and 2Delta_e<=sinh(2lambda_e)^2 a_e.
For max lambda_e=O_beta(N^-1/2), the ratio
sinh(2lambda)^2/[lambda tanh(2lambda)] is uniformly bounded. Thus

    sum_e Delta_e=O_beta(N).                          (4)

In either constant-magnitude edge class with magnitude lambda_c,

    sum_class |a_e-2Delta_e/sinh(2lambda_c)^2|
      <= C_beta sum_class Delta_e^2/lambda_c^2
      <= C_beta N/lambda_c.

Within magnitudes are bounded below by c_beta/sqrt(N); the cross magnitude
is comparable to sqrt(1-u)/sqrt(N). Multiplying by the bounded block
weights D_e and integrating (1-u)^(-1/2) proves the asserted O(N^(3/2))
bound, including the zero-cross endpoint by integration. No stationarity
claim is made at an exactly zero edge.

Consequently the unresolved integrated slack can equivalently, up to an
already summable O(sqrt(N)) pressure error, be expressed in terms of the
ACTUAL logarithmic costs Delta_e, not just exponential flip ratios.

## 3. Biased random cuts give the right weights but no new inequality

For equal blocks let S_w=sinh(2lambda_within)^2 and
S_c=sinh(2lambda_cross)^2. Give independent vertex signs mean +sqrt(v)
on the first block and -sqrt(v) on the second, where
v=(S_w-S_c)/(S_w+S_c). Their random cut has within inclusion probability
S_c/(S_w+S_c), and cross inclusion probability S_w/(S_w+S_c).
These are exactly the inverse-variance proportions needed in the slack,
not an approximate homogeneous partition average.

Let F,G be its within/cross portions, and define the nonadditivity remainder

    R_E=Delta_E-sum_(e in E)Delta_e.

Then (2) yields the exact identity

    sum_within Delta_e/S_w - sum_cross Delta_e/S_c
      = (S_w+S_c)/(S_w S_c) E[R_G-R_F].              (5)

This identifies a precise candidate payment: control the difference of
the two multi-edge nonadditivity remainders. It does NOT establish that
payment. Equality (5) is algebraically a restatement of the singleton
budget once gauge equivalence is used. Global optimality (1) gives
R_E>=-sum_E Delta_e separately, but does not order R_F and R_G. Applying
(1) again to their union gives equality zero and still no ordering.

The analogous multiway decomposition only creates more gauge-equivalent
partial-cut costs; averaging them cannot, by itself, remove these
remainders. Replacing the log moment generating functions in (1) by their
first derivatives would discard row/cut-sized cumulants of leading order
at fixed beta, and is not justified by lambda_e=O(N^-1/2).

## 4. Outcome

The singleton exponential-to-log conversion is a genuine quantified
lemma at actual minimizers. The cut-based multi-edge extension does not
close the interpolation: after its exact cancellations it asks for the
same leading-order nonadditivity payment in a different form. A usable
next theorem would have to exploit non-cut simultaneous moves or an
additional Ising-specific comparison of these remainders. None is proved
here. In particular no generic positive-semidefinite or free-relabeling
argument is substituted for actual anisotropic global optimality.
