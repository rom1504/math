# Independent audit: quantitative hierarchy compression and dilution

2026-09-07. **PASS.** Complete independent read of
`principle_invent_2026_09_07_hierarchy_complexity_dilution.md`.
This is an actual sign extension preserving every old edge, conditional
on the stated quantitative factorization complexity of absolute near-levels.

## Finite resource and exact bridge payment

The concatenated feature columns have squared norm at most
`1/2+sum_j a_j^2 <=1`. The primary Gram--Schmidt Walk normalization was
previously read and audited: variance proxy 40 times squared test norm.
Thus a level-j word has proxy `40 g_j n/a_j^2`, and the identity block
provides proxy 80n for every word. Independent new columns multiply
these proxies by q. The thresholds in (7) pay the two-sided union over
every old word, new word, and all the finitely many levels, with success
probability at least one-half.

The factorization-to-VC inequality is valid: a shattered coordinate set
of size d gives all 2^d labelings; factorization and averaging the squared
signed column sum force gamma squared at least d. Sauer's count therefore
gives log code size at most n H(g)+o(n), using the monotone entropy cap.
All the hierarchy's chosen g_j are eventually below one-half.

Substitution of the weights gives limiting level thresholds equal to
the specified deficit allowance divided by sqrt(2). The strict slack
absorbs each fixed hierarchy's finite entropy and union terms. For an
annular word, its old absolute deficit exceeds eta_j/2; for the terminal
level, no deficit is used. The outside region is paid by the universal
identity-feature bound once epsilon is small compared with fixed eta_0.
Finally the new principal block contributes at most q^(3/2). No cross
term is canceled by an uncharged choice of new-spin reversal.

## Summing the hierarchy

Choose K<kappa<tau/4096. The definition of K allows strictly larger
factorization bounds g(eta)>e(eta) with gH(g)<=kappa eta on all sufficiently
small fixed levels. Since H(g)>=2(log 2)g on [0,1/2],
`g <= sqrt(kappa eta/(2 log 2))`.

For u_j=eta_j/(2 tau epsilon), the main resource sum is at most

    (320 kappa/tau) sum_j u_j/(1+u_j)^2
    <=1280 kappa/tau.

The bound four on this dyadic sum follows by separate geometric bounds
above and below one. The remaining part of the resource is
`O(sqrt(kappa epsilon)/tau^(3/2))`; the corresponding dyadic sum
`sum sqrt(u)/(1+u)^2` is bounded by five. The terminal eta_J in
[epsilon^2/2,epsilon^2] contributes
`O(kappa eta_J/(tau^2 epsilon)+g_J/tau^2)=o(1)`.
Consequently the total is below 5/16+o(1), hence below one-half.

Crucially epsilon is fixed when the hierarchy is chosen. Its depth is
finite before n grows. Each level's limsup factorization bound therefore
holds simultaneously for sufficiently large n. No convergence at
n-dependent vanishing levels is assumed. Epsilon is sent to zero only
after the order limit.

## Consequences and limits of the result

The parent cap bound is exactly
`Q(parent)-Q(old) <= [tau epsilon+epsilon^(3/2)+o(1)] n^(3/2)`.
For a sequence realizing c_*, choosing
`4096K<tau<3c_*/2` would give a parent subsequence below the global liminf.
Thus `K>=3c_*/8192`. Inverting gH(g) near zero gives the stated
square-root-over-log complexity scale along arbitrarily small windows;
it does not give that bound at every window.

The cover estimate follows from gamma2 subadditivity: the repeated-center
matrix has the centers' factorization norm, while the difference matrix
has identity factorization with row norm at most 2 sqrt(r n). Optimizing
the two column weights gives the sum of these norms. The near-half
hierarchical examples have a linear transverse-energy cost after the
coarse-cut part absorbs the R correction. Their near-levels are therefore
within O_c(eta)n of the fibre-constant code, whose normalized gamma square
vanishes. This gives e(eta)=O_c(eta) and K=0 as claimed.

The theorem rules out a quantitative compressibility pattern for actual
liminf realizers. It does not prove that arbitrary selectable minimizers
satisfy that pattern, and it does not establish the original comparable-
child flatification inequality or convergence.

## Golden-ratio resource refinement

The subsequent Section 7 refinement also **PASS**es. For r>1,

    u/(1+u)^2 <= (r+1)/(4(r-1))
                  [1/(1+u/r)-1/(1+ru)].

After cancellation this is exactly the bound
`(1+u/r)(1+ru)/(1+u)^2 <= (r+1)^2/(4r)`, maximized at u=1.
Taking hierarchy ratio R=r^2 makes the right side telescope along the
levels. With coefficient 80(1+zeta) in the weights, the principal resource
is therefore at most

    [20(1+zeta) kappa/tau] r^2(r+1)/(r-1).

Its unique minimum occurs at r=phi=(1+sqrt(5))/2 and equals
`20(1+zeta) kappa phi^5/tau`. An arbitrarily small but fixed identity
weight theta leaves feature budget 1-theta; the outside universal bound
remains finite and is paid by shrinking epsilon after theta is chosen.
The factor 1+zeta supplies strict finite-level threshold slack. The other
resource terms and terminal weight still vanish with epsilon. Thus the
extension holds for every tau>20 phi^5 K, and the liminf necessity sharpens
to `K >= 3c_*/(40 phi^5)`. No additional limiting uniformity is used.
