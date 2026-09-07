# Independent reconstruction: selectors close the exceptional stability sector

2026-09-07. **Analytic audit PASS**, conditional only on the already proved full-row certificate used below. This is an actual restricted-weave upper-construction improvement, not a convergence theorem or a transfer from exact minimizing seeds. Root proposed the selector argument; this note reconstructs its net count, tensor placement, reciprocal filtering, and diagonal payment independently.

**Stronger same-session update:** Section 6 conditions on selectors having NO concentrated word at all. This preserves every old row exponent and removes the phase-gap premise from Sections 2-3. The original sector argument is retained as an independently valid first proof.

## 1. Uniform sparse-spectrum count for every orthogonal matrix

Fix k, m, p_m=k/m, and ANY real orthogonal O of order m. Uniformly choose a selector T of size k and sum over its Boolean words. This is exactly the uniform-selector average of all constant-type vectors

    f in {0,+1/sqrt(p_m),-1/sqrt(p_m)}^m,
    |supp f|=k,

each counted with weight 1/binom(m,k). Each vector has norm sqrt(m). Its normalized Hadamard spectrum is u=Of.

Call a row bad if sum_(|u_j|<=V) u_j^2<=epsilon m. Then its heavy support S has size at most m/V^2, and its distance to R^S is at most sqrt(epsilon m). The radius-sqrt(m) Euclidean ball in R^S admits a sqrt(epsilon m)-net of cardinality at most (1+2/sqrt(epsilon))^|S|. Thus every bad u lies within 2sqrt(epsilon m) of a net point.

Pull the net point back by O^T and round each coordinate to the nearest of 0,+1/sqrt(p_m),-1/sqrt(p_m). If the rounded coordinate differs from f, the coordinate error before rounding is at least 1/(2sqrt(p_m)). Consequently at most 16 p_m epsilon m coordinates differ. This argument does not require that the rounded vector have the prescribed type.

For 16epsilon<1/2 and V>sqrt(2), union over supports and altered coordinates therefore bounds the number of bad f by exp[(xi+o(1))m], where

    xi=h(V^-2)+V^-2 log(1+2/sqrt(epsilon))
          +h(16epsilon)+16epsilon log2.              (1)

Here h is binary entropy in natural logarithms. We used p_m<=1 to give a simple uniform formula. The number of alternatives at a changed ternary coordinate is two, so log2 suffices; using log3 would merely weaken the bound. Taking epsilon small and then V large makes xi arbitrarily small. All constants are fixed before m grows. The count is uniform in O, hence survives any independent distribution of exact Hadamard frames.

The old orbital row factor L is at most one because its kernel entries are at most one. Thus the selector-averaged bad-row weighted count is at most

    exp[(-h(p)+xi+o(1))m]                            (2)

when p_m tends to p. No probability-conditioning fee other than the EXACT division by binom(m,k) is hidden here.

## 2. The strict row-exponent gap is concrete

At p=24/25 and t=97/20, the existing full-row certificate gives exponent

    a=p log2-5151/6250.

Its difference from the limiting sparse-row exponent -h(p) is positive:

    d=h(p)+a=.00920544107172045138...>0.

For example, epsilon=10^-6 and V=100 in (1) give

    xi=.00198494624674918495...,
    d-xi=.00722049482497126642...>0.

These displays were independently checked with 40-digit outward interval arithmetic. No optimized value is needed: equivalently choose epsilon and V by continuity so xi<d/2. Let b=-h(p)+xi<a.

## 3. The row indicator passes the actual graph contraction

Badness depends only on the absolute spectrum multiset, hence is unchanged by the independent output column permutations. At a tensor vertex its indicator is therefore a fixed scalar multiplying the row-orbit tensor. Graph Cauchy-Schwarz gives that scalar times L, because its square inside the squared norm equals itself. This is not a claim that arbitrary correlated stability indicators factorize across vertices.

For a prescribed set of bad rows, independence of the selectors and of the row bases now gives the product of the bad-row bound (2) and the full-row bound exp[(a+o(1))m]. If fewer than delta m rows are light (not bad), at least (1-delta)m are bad. Union over their index set costs at most 2^m=exp(o(m^2)); the total exponent is at most

    [a-(1-delta)(a-b)+o(1)]m^2.                       (3)

This is the missing strict exceptional-profile gap, not merely a small unweighted fraction of codewords.

## 4. Positive light-row mass forces covered physical columns

Suppose at least delta m rows satisfy sum_(|u_i(j)|<=V)u_i(j)^2>epsilon m. Put ell=delta epsilon, so the total directed light energy is at least ell m^2. Every row has total energy m.

Choose U and C so V^2/U^2<=ell/8 and V^2/C<=ell/8. Removing coordinates for which the reciprocal outgoing magnitude exceeds U loses at most V^2 m^2/U^2: there are at most m^2/U^2 such outgoing coordinates. Removing the diagonal loses at most V^2 m, since only LIGHT incoming coordinates are counted.

The total incoming energy over all columns is m^2. Hence at most m/C columns have incoming energy greater than Cm. Deleting them loses at most V^2 m^2/C of the light energy. For large m, at least ell m^2/2 light energy remains on off-diagonal pairs with incoming magnitude at most V, outgoing magnitude at most U, and acceptable incoming column energy.

Under the energy tilt, each such pair has thermal variance multiplier at least

    a0=sech^2(2tUV)>0.

The sum of bounded incoming thermal variances over acceptable columns is therefore at least a0 ell m^2/2. Each column contributes at most V^2 m. With

    kappa=a0 ell/4,     eta=a0 ell/(4V^2),

at least eta m columns have bounded incoming thermal variance at least kappa m AND total incoming energy at most Cm. These conclusions hold for EVERY output column permutation, not merely on average.

The audited Talagrand row theorem, with k/m bounded below and physical deterministic term D=2, gives row probability at most exp[-c k], where c=g^2/(256V^2)>0. Finner over the independent tilted edge signs supplies the uniform factor exp[-c eta mk/2] on this entire light-row sector. Because it is uniform in the column permutations and bases, it can be extracted BEFORE applying the old graph contraction.

## 5. Full diagonal accounting and actual cap consequence

Keep the full woven sign matrix W, including its within-fibre rank-one blocks, and put C=W-diag(W). The exact defect identity is

    D_sigma=2(m^2 k-sigma x^T W x).

The tilted stability law has the same off-diagonal edge means whether the tilt is written using C or W, since their difference is a deterministic trace for fixed S_ii. One first applies the uniform stability penalty under that law. Only then may the nonnegative diagonal terms of D_sigma be discarded for the old PSD graph-kernel bound. No macrodiagonal block is deleted by a triangle inequality.

For a stable witness Q(C)>=q, one has sigma x^T W x>=2q-N, since |Tr W|<=N. This costs only exp(O(m)) in the moment at tilt t/k. It is negligible compared with the strict m^2 savings. Thus the physical diagonal remains paid exactly at its lower scale.

Combining (3) with the light-sector stability penalty improves the full moment exponent t+a by the fixed positive quantity

    Delta=min{(1-delta)(a-b), c eta p/2}>0.

The actual all-order cap coefficient supplied by the same construction is consequently

    (t+a-Delta)/(2t sqrt(p)),                         (4)

strictly below the preceding certificate (t+a)/(2t sqrt(p)). One may take delta=1/2 throughout. Finite-depth approximation losses and the m-to-all-orders passage are handled as in the existing row certificate, at fixed epsilon,V,U,C,kappa,eta,Delta; no uniform growing-depth claim is needed.

The proof gives a potentially extremely small but fixed improvement. A useful numerical decimal improvement would require quantitative optimization not performed here. It does not identify the optimal signing constant or prove convergence of M_n/n^(3/2).

## 6. Stronger selector conditioning: no phase-gap premise

Fix ANY p in (0,1). Choose epsilon>0 and V<infinity with the exponent xi in (1) strictly smaller than h(p). This is always possible. For a fixed orthogonal O and a uniform size-k selector T, let B_O(T) be the INTEGER number of its Boolean words having light spectral energy at most epsilon m. The same counting argument gives, uniformly in O,

    E_T B_O(T)<=exp[-(h(p)-xi+o(1))m].

Consequently

    Pr_T{B_O(T)>=1}<=E_T B_O(T)<=exp(-c0 m)           (5)

for some c0>0 and all sufficiently large m. No factor 2^k is missing: B_O(T) already counts ALL words, and the original ternary net count was divided only by the number of selectors.

For each independently sampled row basis O_i, choose its selector uniformly conditioned on B_(O_i)(T_i)=0. Every such actual selected row submatrix has the required positive light energy for EVERY one of its 2^k Boolean words. This is a property of the matrix, not a restriction artificially imposed on the parent's spin optimization.

The conditioning leaves the old row exponent unchanged. If Z(O,T)>=0 is the orbital row weight and q(O)=Pr_T(B_O(T)=0), then

    E_O E_T[Z(O,T) | B_O(T)=0]
       =E_O [E_T(Z 1_(B_O=0))/q(O)]
       <=(1-exp(-c0 m))^-1 E_O E_T Z.                (6)

The denominator bound is uniform in O. This formula specifies the conditioning order: sample the original basis law, then condition the selector given that basis. There is no unaccounted reweighting of the basis law. The logarithmic cost of (6) is o(1), much smaller than the o(m) allowed in a row exponent.

The no-bad-word event is invariant under final output signed column permutations, since it depends only on absolute spectrum values. Thus those permutations remain uniform and independent conditional on the accepted row data. Sampling the accepted basis-selector pairs independently in distinct fibres preserves all independence used in the graph contraction.

Now EVERY parent spin has m light rows. Apply Section 4 with delta=1. The resulting eta m covered columns, and hence the factor exp(-c eta mk/2), hold uniformly over all spins and all output permutations. The bad-sector split is unnecessary. For any finite-temperature full-row certificate exp[(a+o(1))m] at this fixed p<1 and fixed t>0, the modified ACTUAL ensemble therefore improves the cap coefficient by the strictly positive amount

    [c eta p/2]/[2t sqrt(p)] = c eta sqrt(p)/(4t),    (7)

up to arbitrarily small fixed-depth approximation losses. Neither a>-h(p) nor the numerical phase-gap check of Section 2 is required.

The boundary p=1 is excluded for a substantive reason: h(p)=0, so (5) cannot follow. For a full Hadamard row set, its column-sign words give exactly concentrated spectra, making the no-bad-word property false. Thus the stronger statement does not contradict the full-retention matching eigenvectors.
