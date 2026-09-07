# Wave 51 main audit: spectral-excess coarea below lambda_1

Status: independently derived from (10.1220), pending cross-audit against the
coarea agent and exact common-core normalizations.

Let `D_C>0`, let `P` be the aggregate triple retention in (10.1235), and put

    delta = 1-lambda_1,
    gap   = lambda_1-lambda_2,
    kappa = gap/delta,
    M     = max_{z in C} a_z.

The definitions give the exact identity

    R_C = E[1_C a_z B_z]/(delta D_C) = (1-P)/delta.    (A)

Indeed the triple-event numerator is
`E[1_C a_z <f_z,K f_z>]`, while `B_z=a_z-<f_z,K f_z>`.
Average (10.1220) under the double-incidence law proportional to
`1_C a_z^2`.  Since its weighted mean of `a_z` is at most `M`,

    R_C >= 1+kappa-(1+kappa n)M.

Substitution of (A) and exact cancellation yield

    M >= [P-lambda_2]_+ /
         [(1-lambda_1)+n(lambda_1-lambda_2)].          (B)

Thus the previously targeted non-strict condition `P>=lambda_1` with a
constant spectral ratio is much stronger than convergence requires.  At a
project-row cap it suffices to prove, for some kernel,

    P-lambda_2 >= exp{-O(n^(3/4-c))},                  (C)

provided the denominator in (B) is at most polynomial (it is at most `n+1`
for a Markov kernel with spectrum in `[0,1]`).  Then (B) supplies a center of
degree `exp{-O(n^(3/4-c))}`, directly meeting the bare-tail mass scale; no
slice FKN or constant degree is needed.

Scope and remaining obstruction:

- At positive constant spectral gap, merely `P>lambda_2` by inverse
  polynomial already suffices, even if `P<lambda_1` and hard coarea fails.
- For the common-core `ell=1` kernel, `lambda_2=0`.  Positivity of `P_1` is
  automatic but can be exponentially too small.  Exact Cauchy on element
  loads gives `P_1>=r_C=sum a_z^3/sum a_z^2`; feeding this into (B) may be
  tautological because `r_C<=M`.
- Hence box positivity alone still does not prove (C).  The repaired target is
  quantitative spectral excess above `lambda_2`, not necessarily retention
  above `lambda_1`.

For comparison, if one insists on mixing a failing positive-core kernel with
`K_0`, let `Delta=lambda_1-P>0` and `r=r_C`.  The largest admissible positive-
core weight is `theta*=r/(r+Delta)`, and its best spectral ratio is

    kappa_max = r(lambda_1-lambda_2) /
                [Delta+r(1-lambda_1)].                (D)

Consequently constant-kappa endpoint repair with constant `Delta` already
forces constant `r<=M`; this is exactly the degree-circularity noted in Wave
50.  Formula (B) shows why direct spectral excess can be the cleaner target.
