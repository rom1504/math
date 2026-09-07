# Actual port routing for sparse coherent spectra

2026-09-07. **Proved derivation; audit requested.** This closes a non-uniformly-integrable sector not covered by the bounded-profile pressure theorem. It is an actual low-deficit port assignment and partition lower bound, not an original parent-cap upper bound.

Let S be any hollow full-sign macro matrix of order m with Q(S)<=C_seed m^(3/2). At each row prescribe a real multiset of m-1 ports, with at most R=R_m nonzero values and squared norm at most C_row m. Independently uniformly permute the row multisets. Fix t<infinity and set

    Z_S(t)=E product_(i<j) exp[-t(u_ij-S_ij u_ji)^2].

If R=o(m/log m), then uniformly over all such row profiles and all these seeds,

    log Z_S(t)=o(m^2).                                    (1)

The profiles may be asymmetric, have arbitrary real amplitudes and arbitrarily small atom masses, and carry all their energy in diverging spikes. In particular no uniform-integrability premise is used.

## 1. Greedy compatible matching of nonzero stubs

Fix a small delta>0. Initially set aside all ports with magnitude at most one. There are at most mR such ports, and their total energy is at most mR.

Partition the remaining magnitudes into geometric bins of ratio 1+delta. Since every magnitude is at most sqrt(C_row m), the number L of nonempty possible bins is O_(delta,C_row)(log m). Regard every nonzero port as a stub at its row vertex, retaining its sign and actual magnitude.

Process each bin. Greedily match two stubs at distinct row vertices i,j whenever their signs sigma_i,sigma_j satisfy S_ij=sigma_i sigma_j and the edge ij has not already been used by an earlier stub pair. Each matched pair is placed on the two ports of that edge. Stop when no such pair remains. Every row uses at most R edges throughout all bins.

At the end of a bin, let U be the vertices with an unmatched stub in that bin; choose one such stub sign sigma_i per vertex. By maximality, every UNUSED edge inside U has S_ij=-sigma_i sigma_j. Used edges have maximum degree at most R. Switching by sigma on U, therefore,

    H_(S_U switched)(1)
       <=-binom(|U|,2)+R|U|.

Restriction monotonicity of the absolute cube cap gives

    |U|<=2R+1+sqrt(2Q(S)).                               (2)

This is the only seed input. Restriction monotonicity follows directly by independently randomizing all spins outside U, so it does not assume an induced submatrix is itself optimal.

The unmatched energy in each bin is at most C_row m times its number of remaining vertices. Summing (2) over all bins bounds total unmatched large-port energy by

    C_row m L [2R+1+sqrt(2Q(S))].                        (3)

## 2. Complete the actual port array and pay its deficit

Place all unmatched stubs, including the small ones, arbitrarily into the remaining distinct ports of their row. This is always possible because every row has at most R nonzeros. No already matched edge can be damaged: both of its ports were filled together and are unavailable to later assignments.

On a matched edge, the sign constraint makes the two signed endpoints agree after multiplying by S, while their magnitudes differ by at most a relative delta. Hence its squared mismatch is at most delta^2 times the sum of its endpoint energies. On all other edges use (a-Sb)^2<=2(a^2+b^2). The completed ACTUAL array has total edge deficit at most

    delta^2 C_row m^2
      +2C_row m L[2R+1+sqrt(2Q(S))]+2mR.                 (4)

No real-valued relaxation or fractional degree realization occurs in this construction.

## 3. Its actual probability and pressure

A row with at most R nonzeros has at most m^R distinct multiset arrangements, even if all nonzero magnitudes differ. The constructed array therefore has probability at least exp(-mR log m) under the independent uniform row-permutation law. Since all kernel products are at most one,

    0>=log Z_S(t)
      >=-mR log m
       -t[delta^2 C_row m^2
          +2C_row m L(2R+1+sqrt(2Q(S)))+2mR].            (5)

At fixed delta, R=o(m/log m) and Q(S)=O(m^(3/2)) make every term except delta^2 C_row m^2 negligible after division by m^2. Let delta decrease afterward. This proves (1), uniformly over the stated class.

For example, pure coherent row spectra with support at most m^(1-eta), for any fixed eta>0, have zero leading pressure for every low-cap seed, despite potentially having nonuniformly integrable squared tails. A one-spike row is the simplest case: nearly all spikes can be put on compatible reciprocal matches, with only O(m^(3/4)) unmatched row vertices.

This theorem does not cover a row containing BOTH many bounded diffuse coordinates and sparse large coordinates: retaining and counting the diffuse part during the stub routing requires an additional conditioning argument. Nor does it cover support as large as m/log m with the displayed error, or turn a sector partition statement into a bound on all parent spins. Those obligations are not silently included in (1).

## 4. Stronger descending-bin matching lemma for dense heavy tails

The following refinement is the matching input to the joint full pressure compiler being developed with synthesis. It does not by itself pay the remaining light-coordinate conditioning.

Assume only row energy at most C_row m, and consider all stubs with magnitude at least V. Bin them geometrically with lower endpoints a_b=V(1+delta)^b and process bins in DECREASING magnitude. At any stage of bin b, every previously used edge at a vertex has consumed a stub of magnitude at least a_b. The used degree at every vertex is therefore at most

    d_b=C_row m/a_b^2.

Apply exactly the compatible maximal matching above. Choosing one remaining sign at each residual vertex and using the almost-negative-clique argument gives

    |U_b|<=2C_row m/a_b^2+1+sqrt(2Q(S)).                 (6)

Discard the residual stubs from this bin before proceeding, so they do not block any later edge. The energy discarded in bin b is at most C_row m |U_b|. With L=O_(delta,C_row,V)(log m) bins, the total discarded energy is at most

    E_drop <= 2C_row^2 m^2/[V^2(1-(1+delta)^(-2))]
                +C_row m L[1+sqrt(2Q(S))].             (7)

For bounded-cap seeds this is O_(C_row,delta)(m^2/V^2)+O_(C_seed,C_row,delta,V)(m^(7/4)log m). Every remaining heavy stub has been assigned to an actual reciprocal edge with a sign-compatible endpoint in its relative-amplitude bin. Its matched squared mismatch is at most delta^2 times its endpoint energy.

This avoids random edge-color hosts, degree-repair theorems, or a global support bound. The residual induced submatrix cap is controlled by Q(S) for every switching, so the argument remains valid for arbitrary signed mixtures in each bin. The finite number of bin categories can later be prescribed at much lower entropy than the one-arrangement bound in Section 3, while retaining free permutations of actual values within each bin; that separate counting and light-coupling step is not included in (6)--(7).
