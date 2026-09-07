# Actual low-cap signings with many stable states

2026-09-07. Two explicit constructions delimit the universal stability entropy theorem. The second works arbitrarily close to the original liminf optimum, but does not produce exact minimizers or favorable dilution.

## 1. Balanced fibres and the dependence on the cap bound

Fix even t and any hollow full-sign seed A of order m. Make an actual signing F of order N=mt with off-fibre blocks A_ij J_t and within-fibre blocks -J_t+I_t. If every fibre spin word sums to zero, all off-fibre fields cancel and the within-fibre field at a vertex is exactly its spin. Every such word is therefore strictly one-spin stable, and

    #stable(F) >=binom(t,t/2)^m.                              (1)

For arbitrary fibre sums s_i in [-t,t], the cross energy is H_A(s). Hollow multilinearity gives |H_A(s)|<=t^2 Q(A). The total within-fibre energy has absolute value at most mt(t-1)/2. Conversely take fibre-constant spins attaining an absolute seed extremum. Hence

    |Q(F)-t^2Q(A)|<=mt(t-1)/2.                               (2)

If Q(A)/m^(3/2) tends to c, then Q(F)/N^(3/2) tends to c sqrt(t), while its stable entropy rate is at least t^-1 log binom(t,t/2). Thus a universal stable-count deficit delta(C) for bounded-cap matrices cannot exceed

    log2 - t^-1 log binom(t,t/2)
       = [log(pi t/2)]/(2t)+O(t^-2)                         (3)

along cap coefficients C just above c sqrt(t). This is O(log C/C^2). The seed can be an actual liminf-minimizing sequence, but the amplified family with t>=2 is not near-minimal at its new orders.

## 2. A stable paired extension of an arbitrary actual seed

Let A be of order m, and orient it globally so it has a positive ground x0 with H_A(x0)=Q(A). Its one-spin local fields at x0 are nonnegative.

Fix a pairing of all but at most one old coordinate. For each of q new pairs independently, choose an old cross vector v_j as follows. In coordinates switched by x0, put opposite independent fair signs on every old pair; on an unpaired old coordinate, use an independent fair sign. Thus v_j is a genuine full-sign vector and

    |v_j dot x0|<=1.                                        (4)

Connect BOTH new vertices in pair j to the old vertices with the identical column v_j. Put a negative edge within pair j. Between distinct new pairs j,l put the full 2 by 2 constant tile (D_q)_jl J_2, where D_q is any low-cap hollow full signing with Q(D_q)<=q^(3/2).

For each of the 2^q choices of a Boolean word z, extend x0 by (z_j,-z_j) on new pair j. The new-pair sums are zero. Therefore:

- every old local field remains its original nonnegative ground field;
- the field at a new vertex, multiplied by its spin, equals 1 plus or minus v_j dot x0, and is nonnegative by (4);
- all cross energies cancel and every negative pair edge contributes +1.

These are 2^q distinct exactly one-spin-stable words, all of energy

    Q(A)+q.                                                 (5)

This includes odd m; some new local fields can then be zero. For even m all new local fields are strictly positive.

## 3. Uniform paid cap bound

For every old Boolean x, a random balanced column obeys

    E exp(lambda v_j dot x)<=exp(lambda^2 m).

Indeed its independent pair coefficients are in {0,+2,-2}, giving variance proxy at most 2m; the possible singleton does not increase this upper bound. For arbitrary new spins y, the bridge energy is

    sum_j (y_(j,1)+y_(j,2)) v_j dot x.

The q columns are independent and the squared new coefficients sum to at most 4q. The bridge energy therefore has variance proxy at most 8mq. Union over both tails and all 2^(m+2q) parent words gives a realization with

    beta(C)<=4 sqrt[mq(m+2q+2)log2].                          (6)

For example, the right side makes the union failure probability at most 1/2. Here C is the m by 2q bridge with duplicated columns, so its bilinear cap is exactly the relevant cross energy.

The new principal block has cap at most 4Q(D_q)+q. Consequently the ACTUAL parent A^+ of order N=m+2q satisfies

    Q(A^+)<=Q(A)+4 sqrt[mq(m+2q+2)log2]+4q^(3/2)+q.           (7)

Every old edge is preserved. There is no leading old-edge deletion or unpaid cancellation in this cap estimate. The cancellation in (5) is used only to exhibit the stable words, not to upper-bound other parent words.

## 4. Arbitrarily near-minimal stable multiplicity

Let c_inf be the positive finite liminf of M_m/m^(3/2), and take an actual liminf-realizing seed sequence. For q=floor(rho m), (7) has normalized incremental cost at most

    e(rho)=4 sqrt[rho(1+2rho)log2]+4rho^(3/2),                (8)

before dividing by the parent normalization (1+2rho)^(3/2). This is O(sqrt(rho)).

For every fixed sufficiently small eta>0, choose a fixed rho=c eta^2 with a sufficiently small absolute constant c. Along sufficiently large selected orders, the parents satisfy

    Q(A^+)/N^(3/2)<=c_inf+eta,

and all the words in (5) have deficit at most eta N^(3/2). Their number is

    2^q=exp[(rho log2/(1+2rho)+o(1))N]
        >=exp[c' eta^2 N].                                (9)

The constants may be fixed uniformly for eta in a bounded small interval. The seed convergence error is taken small after eta and rho are fixed. Thus exponential stable near-ground multiplicity persists arbitrarily close to the optimum coefficient. This does not assert the parents are exact minimizers, and it does not contradict a fixed positive entropy deficit from the entire cube.

The extension cost in (8) is O(sqrt(rho)), whereas profitable dilution would need O(rho) with a sufficiently small coefficient or o(rho). This construction therefore supplies a rigorous boundary and actual near-minimal examples, not the missing convergence mechanism.
