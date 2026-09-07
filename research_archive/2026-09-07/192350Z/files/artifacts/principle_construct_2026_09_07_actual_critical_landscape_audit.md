# Independent audit: actual critical-scale landscape

2026-09-07. **PASS for the deterministic graph conversion, actual-sign compiler, and all-spin cap/landscape proof.** Audited source: `principle_invent_2026_09_07_actual_critical_landscape.md`. I read its complete text; I did not independently reopen the cited Friedman/Bordenave primary in this audit, so the existence of the fixed-degree expander pieces remains the author's verified external input.

The prescribed matching switches preserve every degree, use each switched vertex at most once, and introduce no multiple edges. Deleting matching edges loses at most one crossing edge per minority vertex, so the stated retained internal expansion follows from the input expansion. The dyadic ceiling errors remain negligible uniformly through the last scale because M psi(r_J) tends to infinity.

For a coarse set S, majority rounding gives internal cut at least 38t and changes the external cut by at most t. Thus 37t+cut_external(U) is valid. When t<s/4, both U and its complement have size at least 3s/4. The one excluding the core is a satellite union, whose largest dyadic component has size at least half its total. This gives r_*>=3r/8 and exactly the factor 9a/32 in the quadratic-log profile.

The clique blowup cut identity is exact, including its internal-clique term. A rounded-crossing coarse edge contributes at least one half to f_u+f_v-2f_uf_v. Thus the same two-case rounding argument gives c_G=min(1/8,9c_H/32). There is no assumption that a physical set is a union of whole cliques.

The paired physical matrix has energy

    H=z^T(2G+I)z+w^T(2R-I)w,

where z,w have disjoint supports and z_i^2+w_i^2=1. The within-pair +1 edge is responsible for the displayed +I and -I terms; these signs check.

The complement signing has operator norm at most 8sqrt(m) by the stated real-sphere net argument. For each fixed unit vector, the sign-sum coefficient-square total is at most two, giving the tail 2exp(-t^2/4). A 1/4-net of size at most 9^m with threshold 4sqrt(m) proves the asserted operator bound. This estimate is uniform over the nonedge support.

For D=2q+1 the positive deficit is EXACTLY

    Dm-H=2z^T L_G z+w^T[(2q+2)I-2R]w.

It is at least 2sum_edges(z_i-z_j)^2+q||w||^2 for the stated large-order parameters. Equality at zero deficit forces w=0 and, by connectedness, z constant, so the two constant physical words are the only positive grounds.

The absolute-polarity issue is correctly paid, not inferred from positive maximality. The lexicographic clique spectrum gives lambda_min(G)>=-(d-1)M-1. Therefore

    Dm+H >=(4M-2)||z||^2+(2q-2||R||)||w||^2>=2Mm.

Generic regularity alone would not give this leading negative gap; the clique term is essential. Thus the only absolute grounds are indeed the two positive-polarity constants.

For arbitrary physical x, after a global reversal let s count negative aligned pairs and t count misaligned pairs. Then r=(2s+t)/(2m)<=1/2. If t>=rm, the misaligned sector pays a linear-in-r cost. Otherwise rm/2<=s<=rm. Every edge leaving {z=-1} contributes at least one to (z_i-z_j)^2, including edges ending at z=0. Applying the full cut bound and psi(s/m)>=psi(r)/4 yields the claimed c_G M^3 psi(r)/2 lower bound. This avoids any unpaid subtraction for edges ending in the misaligned sector.

Whole-community reversals have w=0 and exact positive drop 8cut_G(C_j times [M]). This gives the matching upper critical scale, stays strictly positive by connectedness, and remains far below the total cap, so the absolute deficit agrees with the positive deficit there.

Finally n=2M^2 and Q=(202M-1)M^2 give cap coefficient 101/sqrt2 as stated. The family proves actual full-sign nonvacuity of the critical landscape phenomenon. Its large coefficient makes it neither an asymptotic minimizer nor a sharpness proof for the near-minimizer augmentation constant.
