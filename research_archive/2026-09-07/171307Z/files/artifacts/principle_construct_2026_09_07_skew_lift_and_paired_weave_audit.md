# Skew-cross Boolean objective and paired-weave audit

2026-09-07, second campaign. Independent construct-agent work. No Boolean half-floor proof or scalable subhalf result is claimed yet. The synthesis agent owns the candidate paired-weave construction and its uniform scalar certificate.

## 1. Exact objective

For symmetric hollow full signs A and skew hollow signs C, let

    L=[ A   C ]
      [-C  -A ].

For Boolean x,y put u=(x+y)/2 and v=(x-y)/2. Their supports are disjoint and partition the n coordinates. Direct expansion gives

    H_L(x,y)=2u^T(A-C)v.

Swapping u and v preserves u^TAv and negates u^TCv. Thus

    Q(L)=2 max_disjoint (|u^TAv|+|u^TCv|).

Partial disjoint supports give the same maximum: missing coordinates can be filled by multilinearity without decreasing a chosen signed objective. Equivalently T=(A-C)/2 is a signed tournament matrix, with exactly one nonzero signed entry on each unordered pair, and

    Q(L)=4 max_disjoint |u^T T v|.

The proposed Boolean half-floor is therefore a universal asymptotic lower bound n^(3/2)/(2sqrt(2)) for this directed signed cut norm. The Frobenius/operator argument in the prior skew-lift audit proves only a spectral-certificate floor. It does not prove this Boolean statement.

I did not find a direct lower argument reaching that constant. Greedy random-partition/one-shore sign selection gives a lower constant below the already available original lower bound, and tensoring a small signed tournament with Hadamards does not automatically preserve its Boolean directed-cut optimum. Those are not scalable falsifiers or proofs.

## 2. A conditioning shortcut that is invalid

Independent balanced incidence rows induce a nearly uniform law on their edge products, up to exp(O(m log m)) density factors. This does not make prescribing all products eta_ij eta_ji=-1 free. Under independent balanced rows that particular product pattern has probability

    2^(-binom(m,2)) exp(O(m log m)).

At the weave scale N comparable to m^2 this is a leading entropy cost. A density comparison to uniform EDGE PRODUCTS must not be identified with a comparison of a conditioned compatible incidence law to the original independent ROW law.

The fixed-tournament construction below avoids this shortcut completely.

## 3. Fixed-tournament product-orbit construction: audit PASS

Fix an almost-regular tournament eta. At even Hadamard order m, each off-diagonal row has either m/2 or m/2-1 positive entries; choose its diagonal sign so the full row has exactly m/2 positive and m/2 negative entries. In each fibre use a full Hadamard root

    H_i=[ H_i,+   H_i,- ]
        [ H_i,+  -H_i,- ],

with independent recursive child bases. Permute the two column groups independently onto eta_i=+1 and eta_i=-1 positions. Select a fixed collection of k/2 paired rows, where k/m tends to p. The involution exchanging the two rows in every selected pair multiplies a column by eta_ij. Since eta_ij eta_ji=-1, every cross-fibre tile changes sign under the simultaneous involution. The outer symmetric sign on the tile does not affect this relation.

The tournament is FIXED before the independent fibre bases and within-group permutations are sampled. Consequently the row randomizations remain independent; no tournament-conditioning probability is paid.

For the folded PSD Gaussian kernel, the graph Cauchy-Schwarz/Finner estimate only needs the Hilbert tensor norm of each row-orbit average. For the orbit group G_plus times G_minus, the squared tensor norm is exactly

    P_t(v_plus) P_t(v_minus).

Thus its contribution is L_t(v_plus)L_t(v_minus). Complementary group labels at the two ends of an edge do not invalidate this: each edge still contracts the same PSD feature inner product, while the row norm factors according to that row's own product group. There is no requirement that the endpoint group names agree.

Deleting one diagonal coordinate affects only one group and costs the usual polynomial factor. Repeated coordinate values do not change the finite orbit argument.

## 4. Exact paired-spin counting and source scale: audit PASS

Among the k/2 selected row pairs let a have equal spins, and set theta=2a/k. There are exactly

    2^(k/2) binom(k/2,a)

such row configurations. After the normalized H2 root transform, the plus and minus child sources have nonzero densities p theta and p(1-theta), common nonzero magnitude sqrt(2/p), and second moments 2theta and 2(1-theta).

The scale relation E_t(cX)=E_(t c^2)(X) therefore gives the candidate row exponent

    R_(t,p)(theta)
      =p/2 [log 2+h(theta)]
        +1/2 E_(2t theta)(nu_(p theta))
        +1/2 E_(2t(1-theta))(nu_(p(1-theta))).

Here nu_q is the symmetric ternary unit-variance law with nonzero probability q. The zero-density endpoint is interpreted through the zero source, with zero contribution. The common physical selector is legitimate because the earlier recursive row theorem is uniform over fixed selectors, and the two child bases are independent even for the same input support.

The full row-spin sum has only O(k) possible theta values. Uniform fixed-depth estimates and a uniform strict bound over theta, if supplied, would pay this polynomial overhead. Checking isolated theta values is insufficient.

## 5. Off-diagonal energy and full-sign completion: audit PASS

Let E_off be the ordered cross-fibre quadratic sum and h_i the transformed fibre spin. The exact identity is

    D_off=2[m^2 k-sum_i h_i(i)^2-sigma E_off].

Therefore the event sigma E_off>=(1-gamma)m^2 k implies D_off<=2gamma m^2 k. The diagonal square term has the favorable sign and may be dropped in this implication. This uses the off-edge kernel directly; it does NOT bound removal of macrodiagonal energy by a triangle inequality.

Each fibre can be filled with a small anti-invariant signing of its own: choose its symmetric/skew components randomly and use a union bound to obtain cap O(k^(3/2)). The m fibre caps total O(m k^(3/2))=O(N^(5/4)) in the proportional regime. Each paired involution leaves one matching edge forced to zero; all such missing edges together number N/2 and may be filled at O(N) cap cost. Thus a proved uniform strict scalar certificate for the preceding row exponent would yield a genuine scalable strict-subhalf member of the prescribed anti-invariant family before matching fill.

## Current conclusion

The Boolean half-floor remains unresolved by the lower arguments here. The paired-weave realization and product-orbit reduction have passed independent algebraic/probabilistic audit. The remaining obligation is a uniform certified strict upper bound on R_(t,p)(theta), with the same fixed p,t and a legitimate finite recursive depth. It is not a free-conditioning assumption or an operator-norm claim.

## 6. Uniform density interpolation and common depth

The proposed test p=24/25, t=7/2 has fixed child amplitude squared 25/12. Let mu_r=(1-r)delta_0+(r/2)(delta_a+delta_-a), with a=sqrt(25/12). Then r->E_t(mu_r) is convex. Here is a channel-level proof, avoiding an assumption about the recursive Bellman iterates.

For a fixed channel X->L and a mixture input with hidden label J, the mean conditional variance satisfies

    E Var(X|L) >= E Var(X|L,J),

so it is concave in the input law. The function g_t is decreasing and convex: its Gaussian correlation representation is a supremum of affine functions of the variance, with nonpositive slopes. Also

    I(X;L)=I(J;L)+I(X;L|J)

for the Markov chain J->X->L, proving that mutual information of a fixed channel is concave in the input law. Therefore g_t(E Var(X|L))-I(X;L) is convex in the input law. Taking the supremum over channels preserves convexity.

It follows that certified upper values at r_j=(24/25)j/20 give upper chords for every intermediate r. On the corresponding theta intervals the sum of the two E chords is affine, because the grid is symmetric under theta->1-theta. Thus the row upper envelope has the form

    (p/2)h(theta)+b theta+c.

Its maximum on each closed interval is at an endpoint or at the clamped critical point theta=1/(1+exp(-2b/p)). This reduces continuum verification to finitely many explicit exponential/logarithmic interval checks.

Finally the actual row construction requires ONE fixed recursive depth for all theta. Pointwise depths at the gridpoints alone do not automatically suffice, unless convexity of the finite-depth iterates is separately proved. Instead use the archived explicit stopping estimate for B^r Phi-E: all root moments here are at most 2, and all initial budgets Phi-G are bounded by -g_t(2). Choose a common moment cutoff and gap tolerance, use the common positive compact-set drift, then choose one finite depth. This provides the needed uniform convergence on the entire source family. The fixed-depth type alphabets and their Stirling errors are uniform, including zero-count endpoint types.

## 7. Scalar verifier source audit and independent chord replay

Read `computations/principle_synthesis_2026_09_07_anti_weave_interval.py` and its imported base verifier completely. The fixed-amplitude generalization passes: A=exp(a^2 lambda z^2)-1, B=cosh(2a^2 lambda z)-1, variance=r a^2, and both envelope partial derivatives have the correct factors of r and a^2. The clipped weight interval contains the maximizing weight throughout each box, including active-set changes. The monotone bound c_t(upper lambda)-(lower lambda)*variance is safe, and all acceptance bounds are outward rational interval endpoints. The low-precision cutoff 1/(2a^2)=6/25 follows from the bounded centered source's a^2-subgaussian estimate; that branch is bounded above by g_t(variance).

The existing three-atom reproduction reduction is valid for every 0<r<1 and scales exactly to this source, not just the old r=24/25 test. Endpoint r=0 is treated as a zero source separately.

An independent chord verifier was implemented in `computations/principle_construct_2026_09_07_anti_chord_audit.py`. It does not use the synthesis code's optimized entropy-affine critical-point formula: it covers each theta cell by 200 subcells and applies the concave entropy tangent at each midpoint. At 55 interval digits it certified, CONDITIONAL ON ALL 20 PROPOSED RATIONAL E POINT BOUNDS, a maximum cap coefficient below

    24460867079760851396003115536506909043664574468748570519
    /49039857307708443467467104868809893875799651909875269632
    <0.498795642<0.499.

The synthesis critical-point chord implementation was also rerun at 55 digits and gives the consistent slightly sharper upper 0.498795640216. No point certificate was duplicated in this audit. Until every point rectangle has closed, these chord conclusions remain explicitly conditional and do not establish the asymptotic family.

## 8. Full canonical theorem audit

Read `principle_synthesis_2026_09_07_anti_invariant_weave.md` completely. Its construction, normalization, exact-family landing, filler, and all-order restriction arguments PASS independently.

The explicit filler variance bound was checked: for order k=2r and a fixed test, an unordered original pair contributes two coefficients of magnitude 2 exactly when its agreement labels differ, and otherwise both vanish. Hence the total squared coefficient sum is 8|S||T|<=2r^2=k^2/2. The stated threshold k sqrt((k+2)log 2) therefore gives total union failure at most one half.

Paired principal restriction preserves the anti-invariant family, with one zero matching edge per retained pair. Reordering one representative from each pair yields exactly a symmetric full-sign A and a skew full-sign C off their diagonals. Ordinary arbitrary-vertex restriction would not preserve this form; the canonical proof correctly restricts pairs.

The p_m->p amplitude passage is also uniform: for any source of second moment v, |E_t(mu)-E_s(mu)|<=v|t-s|, since the Gaussian correlation representation gives t-slopes between -v and zero. Scaling the fixed-amplitude source converts its amplitude change into a time change of order |p_m-p|, with moments uniformly bounded. This covers all theta, including zero-density endpoints.

At the time of this full-file audit, synthesis reports completed directed covers for indices 1,10,20 (4465,3041,26097 checked boxes respectively), while its all-20 run is in progress. This audit does not relabel that partial computational completion as a full certificate.

## 9. Completion update, 16:42 UTC

The canonical result file `computations/results/principle_synthesis_2026_09_07_anti_weave_interval.json` now reports ALL TWENTY completed directed covers, with 100466 checked boxes and 50243 accepted leaves. Independently checked that the indices are exactly 1 through 20, every status is a completed certificate, every worst accepted bound and Gaussian branch lies strictly below its rational target, and every binary-tree count obeys checked=2*leaves-1. The verifier source and the exact mathematical bounds had already been independently audited above; the continuum chord bound was independently implemented and replayed.

Final audit conclusion: **PASS**. The completed scalar proof closes the previously conditional construction. There are all-order skew anti-invariant lifts L(A,C) with

    limsup Q(L(A,C))/(2n)^(3/2) < 499/1000 < 1/2.

Therefore the proposed universal Boolean half-floor for this exact family is false. This is a scalable family counterexample, not a finite cap observation. The A children produced here are not asserted to be exact original minimizers or to preserve an input seed, so no favorable recurrence or convergence theorem follows.
