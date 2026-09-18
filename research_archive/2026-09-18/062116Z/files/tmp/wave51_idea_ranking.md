# Wave 51 independent idea ranking

Evidence basis: README, the Wave 50 Updated frontier, and STEERING at the
Wave 51 boundary. This ranking was made before assigning agents.

1. **Positive-margin abundance from exact-minimizer averaging.** Average the
   condition `g=(1-p^2)e-(Y-t)>0` in (10.1244) over oriented child grounds and
   subsets, and seek an exact identity or one-sided inequality forcing enough
   positive-margin states for the positive Cantelli branch of (10.1245) to
   meet (10.795). State the precise saved-mass exponent and audit circularity.
2. **Near-Parseval abundance on negative-margin states.** On `g<=0`, study
   `D=q^2-e^2-V` jointly with the negative margin. Seek a minimizer-specific
   inequality forcing the reverse two-moment numerator in (10.1245) to be
   positive on `exp(-O(n^(3/4-c)))` mass, or exhibit a scalable obstruction.
3. **Direct lower tail from the linear-plus-quadratic completion structure.**
   Write `Z(w)=L(w)+Q(w)` and pair `w` with `-w`; derive a lower-tail theorem
   that uses more than support and variance. A useful statement must operate
   at deviation `O(n^(3/2-c))` and lose at most `exp(O(n^(3/4-c)))`.
4. **Higher-moment reverse-chaos bound specialized to sign quadratics.**
   Compute conditional third/fourth moments or use a verified primary-source
   reverse-small-ball theorem to lower-bound the completion CDF. Check all
   hypotheses and reject generic upper-tail concentration as irrelevant.
5. **Spectral-profile dichotomy for completions.** Split the conditional
   completion matrix according to operator/Frobenius norm. In the spiked case
   seek an explicit sign-chaos lower tail; in the flat case derive a
   minimizer-specific abundance or restriction consequence. The two branches
   must cover all exact minimizers without assuming (10.795).
6. **Cap-inflated coarea repair.** Starting from a box witness at cap `R0`,
   prove that adding rows up to some `R <= C(R0+n^2)` necessarily raises an
   exact triple-retention profile above `lambda_1` while keeping a uniformly
   positive harmonic gap. Quantify the required increments rather than just
   their monotonicity.
7. **Independent-resampling mixture without degree circularity.** Analyze the
   exact `K0`/positive-core tradeoff (10.1249)--(10.1250) and determine whether
   any choice of mixing weight and scale has uniform `kappa>0` without already
   assuming a constant point mass or bounded degree.
8. **Joint optimization of cap and kernel.** Couple the monotone cap profile
   with the mixture parameter, looking for a telescoping or pigeonhole theorem
   that forces either pure retention after bounded cap inflation or a
   noncircular mixed-kernel certificate.
9. **Finite exact-minimizer falsification atlas.** Extend the exact searches
   for the distributions in ideas 1--5 and the cap/mixing profiles in ideas
   6--8, emphasizing invariant quantities and scalable patterns rather than
   treating small cases as an asymptotic proof.
10. **Non-switching tight-decomposition certificate.** Search for a
    higher-order incompatibility outside the switching polytope exposed in
    (10.1251)--(10.1252). This is deliberately low-ranked and should be
    reopened only if a concrete non-switching invariant appears.

Selected independent attacks: 1--2 as a completion-abundance program; 3--5
as a direct/higher-moment negative-tail program; and 6--8 as a controlled
coarea-repair program. Ideas 9 and 10 are supporting falsification and a
dormant alternative, respectively.
