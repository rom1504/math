# Dephased Hadamard bridge and the limitation of scalar overlap envelopes

## 1. The bridge estimate is exact and elementary

Let H be an n-by-n dephased sign Hadamard:
H H^T=nI, H1=n e1, H^T1=n e1. For arbitrary chosen centers x0,y0 use
C=D_x0 H D_y0. Write u=D_x0 x, v=D_y0 y, alpha=1^Tu/n,
beta=1^Tv/n, u_perp=u-alpha1, v_perp=v-beta1. Exactly,

 u^T H v=u_perp^T H v_perp+n(alpha v1+beta u1-alpha beta).

Consequently, uniformly over ALL spins,

 |x^T C y|<=n^(3/2)sqrt((1-alpha²)(1-beta²))+3n.       (1)

No probabilistic Gaussian assumption, typical-cloud restriction, or
special property of the selected centers is used.

## 2. Low-cap full signs retain energy oscillation at zero overlap

Define L0=2/(3sqrt(3pi))=0.217156671956853... . Let A_n be ANY hollow
symmetric full-sign sequence with Q(A_n)=O(n^(3/2)). For every chosen
Boolean center v_n, on the exactly zero-overlap slice (n even),

 max H_A(x)-min H_A(x)>=(2L0-o(1)) n^(3/2).           (2)

In particular the absolute cap on that slice is at least
(L0-o(1))n^(3/2). The error is uniform over matrices with a fixed cap
constant and over centers. For odd n use the nearest balanced slice,
of overlap magnitude1/n, with an additional O(n) error.

This theorem requires no convergence assumption and no actual-optimality
assumption. Actual minimizers satisfy its hypotheses by the banked upper
bound. The proof is a paired-spin greedy construction; optimizing its
split improves the naive equal-split constant sqrt(2/pi)/4.

### A. A spectral consequence used only to control exceptional pairs

We need only the elementary bound ||A||op=O(n^(5/6)), avoiding any
interpolation theorem. Let lambda=||A||op and choose a real unit eigenvector
z for an eigenvalue of magnitude lambda. The entry bound gives
lambda||z||infinity<=||z||1<=sqrt n. On the other hand, the bilinear cube
norm beta(A) gives lambda=|z^T A z|<=beta(A)||z||infinity². Combining yields
lambda³<=n beta(A). Finally beta(A)<=4Q(A) by polarization and multilinear
cube contraction: for Boolean x,y, set u=(x+y)/2,w=(x-y)/2 and use
x^TAy=2H_A(u)-2H_A(w). Thus lambda=O(n^(5/6)).

Thus

 Tr A^4<=||A||op² ||A||F²=O(n^(11/3)),
 ||A1||²<=n||A||op²=O(n^(8/3)).                      (3)

Switch the arbitrary center v to1; these estimates still hold. Let
d=A1. For a uniformly selected distinct vertex pair u,v, put
r_i=A_iu-A_iv for i outside that pair. With M=n-2,

 R2=sum r_i²=2(n-2)-2(A²)_uv,
 R1=sum r_i=d_u-d_v.

Equation(3) implies

 E[(A²)_uv²]=O(n^(5/3)), E[R1²]=O(n^(5/3)).

Hence, outside an o(1) fraction of pairs, R2=2n+o(n) and R1²=o(n²).
For example thresholds |(A²)_uv|<=n^(11/12), R1²<=n^(11/6) leave an
exceptional fraction O(n^(-1/6)).

### B. Random pairing produces a Gaussian-sized field

Fix p in(0,1). Randomly perfectly match the n vertices, and designate
l=pn/2+O(1) pairs as I; all remaining pairs are J. Assign opposite spins
inside EVERY pair, so all resulting full spins are balanced. On I use
independent signs epsilon_i; later optimize the J signs.

For a J pair(u,v), and an I pair(i,i'), its cross coefficient is

 b=A_iu-A_i'u-A_iv+A_i'v=r_i-r_i', with |b|<=4.

Conditional on(u,v), the2l vertices of I are a uniform subset of the
remaining M vertices, uniformly paired. Therefore their field variance
S=sum_I b² has the exact conditional mean

 E S=[2l/(M-1)] [R2-R1²/M].                          (4)

For the good pairs in A this is2pn+o(n). Also Var(S)=O(n), uniformly in
the bounded vector r: realize the subset and matching by a random
permutation. Swapping two entries changes at most two paired terms,
by a bounded constant. Exposing the permutation successively and coupling
completions by a transposition bounds the martingale increments by a
constant, giving the asserted variance bound. Thus S/n converges in
probability to2p for a random J pair and its random I matching.

Conditional on the pairing, the field is sum_I b epsilon_i. Its maximal
summand is bounded4 and its variance is S. On the good event S is of
order n, the elementary Lindeberg central limit theorem applies uniformly.
One may verify it directly by expanding product cos(t b/sqrt S);
sum b^4/S²=O(1/n). Its normalized second moment is1, which gives uniform
integrability for the absolute first moment. Consequently

 E_epsilon |field|=(sqrt(2/pi)+o(1))sqrt S
                 =(2sqrt(p/pi)+o(1))sqrt n.

Exceptional pairings contribute nonnegatively and can be discarded for
the lower bound. Summing over the (1-p)n/2+O(1) J pairs proves existence
of a pairing and I signs such that the optimized cross energy is at least

 [(1-p)sqrt(p)/sqrt(pi)-o(1)]n^(3/2).                 (5)

### C. The two balanced completions give oscillation, not just one cap

Choose each J-pair sign to align its field. Let L denote the resulting
cross energy. Reversing ALL J spins preserves all internal I and internal
J quadratic contributions, so the two balanced full energies are T+L
and T-L. Their difference is2L regardless of T. Taking p=1/3 maximizes
(1-p)sqrt p, giving L>= (L0-o(1))n^(3/2) and proving(2).

For odd n, discard one vertex, apply the even construction, and give the
discarded vertex a fixed sign. Its incident edges alter either energy
by at most n-1; the leading oscillation bound is unchanged.

## 3. Even exact signed energy profiles cannot close (1) by overlap alone

Consider two equal-order low-cap children A,D, arbitrary centers, and
arbitrary global child polarities. On their balanced slices write their
energy intervals as [l_A,u_A] and [l_D,u_D]. Equation(2) gives each width
at least(2L0-o(1))n^(3/2). Hence

 max_{balanced x,y}|H_A(x)+H_D(y)|
 >=[(u_A-l_A)+(u_D-l_D)]/2
 >=(2L0-o(1))n^(3/2).                               (6)

At alpha=beta=0 the bridge term in(1) is n^(3/2). Thus an upper certificate
which adds that bridge bound to even the EXACT separately optimized signed
child energy profiles has value at least

 [1+2L0-o(1)]n^(3/2)= [1.4343133439...-o(1)]n^(3/2).

The favorable equal-child leading budget is2sqrt2 c n^(3/2), where
Q(A),Q(D) approximately c n^(3/2). The scalar-overlap certificate cannot
reach this budget for c<(.5+L0)/sqrt2=0.5071063459...; in particular it
fails throughout the known regime c<1/2. Opposite child polarities do not
repair it, because(6) uses interval widths rather than assumed one-sided
ground energies.

CRITICAL SCOPE: this does NOT prove that a dephased Hadamard bridge has
such a large actual parent cap. The bridge spectral envelope need not be
attained at the balanced child-energy witnesses. It proves that optimizing
child energies separately at a fixed overlap and then adding (1) loses
too much information. A successful construction would have to retain
joint bridge/child structure beyond these two scalar overlaps. It also
does not constrain global changes to the old child edges.

## 4. Extension to any fixed finite partition

The same oscillation constant survives balancing inside every atom of
any partition into a fixed number R of vertex classes. Initially suppose
all atom sizes are at least epsilon n and even. Pair only within atoms,
and designate a fraction p of each atom's pairs as I. There are still
(1-p)n/2+O_R(1) greedy J pairs.

Let P be the orthogonal projection onto the atom-constant vectors. For a
J pair(u,v), the loss in its expected field variance, in addition to the
column correlation term, is asymptotically

 ||P A(e_u-e_v)||²=sum_atoms [sum_{i in atom}(A_iu-A_iv)]²/|atom|.

Excluding u,v only changes bounded terms. Because
||P A||F²<=R||A||op²=O_R(n^(5/3)), the average of this loss over J pairs
is O_{R,epsilon}(n^(2/3))=o(n). Here each vertex has uniform marginal,
and pair distributions within macroscopic atoms have density at most a
constant relative to uniformly selected pairs. Likewise their averaged
squared column correlation remains O_{R,epsilon}(n^(5/3))=o(n²).

Applying the exact formula(4) within each atom and summing gives
E S=2pn+o(n) for almost all J pairs. Independent within-atom permutation
exposure still has total variance O(n). The bounded-summand CLT and greedy
completion are unchanged, so(2) holds with the same constant while both
witnesses are balanced on EVERY atom.

Small atoms do not require an uncontrolled edge-deletion estimate. For
fixed epsilon, omit atoms smaller than epsilon n from the greedy step;
their total size is at most R epsilon n. Extend both resulting witnesses
by the SAME randomly chosen balanced spin on each omitted even atom.
Each omitted coordinate has mean zero. The expected difference of the
two full energies is exactly the difference already obtained: omitted
internal terms cancel, and all added cross terms have zero expectation.
Some common extension therefore retains at least that difference. Let
n tend to infinity first, then epsilon tend to zero. This proves the
same leading constant for arbitrary atom sizes when R is fixed.

If an atom is odd, leave one vertex aside and fix its spin. There are at
most R such vertices, whose incident terms cost O_R(n); each atom then
has imbalance at most1. This is the appropriate near-balanced statement
when exact atom balance is parity-impossible.

In particular, finitely many chosen Boolean centers define a fixed finite
partition by their joint coordinate labels. The constructed witnesses
have overlap O_R(1/n) with every center, yet retain oscillation
(2L0-o(1))n^(3/2). Thus merely adding a bounded number of overlap labels
does not remove the residual energy obstruction. A residual-operator
certificate that continues to charge a full unit bridge term on that
joint balanced slice still has the same failure as Section3. No theorem
for growing R is asserted; its constants and atom-size dependence must
be tracked separately.
