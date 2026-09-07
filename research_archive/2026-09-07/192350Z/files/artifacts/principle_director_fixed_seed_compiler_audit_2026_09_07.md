# Director reconstruction: fixed-seed mixed-orbit compiler

2026-09-07. **PASS, with a strict value-versus-information distinction.**
I reconstructed `principle_synthesis_2026_09_07_fixed_seed_mixed_orbit_weave.md`
after the separate construction-role audit, rather than treating that verdict
as a proof step.

## Imported dependencies actually needed

Only the old finite signed-input Bellman type recursion, its Gaussian-boundary
closure, terminal-order realization, and the replayed ternary E certificate
are imported. The new terminal orbital inequality is reconstructed directly
and does not import an unproved unsigned-input supersolution.

## Terminal orbital comparison

For the Gaussian Fock feature phi_t, the signed-input orbit covariance has
operator norm exactly P_signed(v). This follows from its finite nonnegative
Gram matrix's constant row sum, including multiplicities if the orbit has
a stabilizer. Ordinary-permutation output invariants of homogeneous degree d
have dimension at most the partition number p(d), not an exponential-in-d
dimension. The partition generating function gives

    sum_(d<=D)p(d)<=exp[pi sqrt(2D/3)].

For D=ceil(e^2(2tC+1)m), the feature tail is a Poisson tail with mean at
most2tCm, hence at most exp[-(4tC+1)m]. Since P_signed(v)>=exp(-4tCm),
it can be absorbed multiplicatively. Cauchy--Schwarz/Jensen therefore give
the mixed ordinary-output/signed-input inequality at exp(O(sqrt m)) cost,
uniformly in terminal orthogonal U. Degree preservation by U is essential.

At each recursive node the ordinary-output invariant norm is submultiplicative
under concatenation because its invariant subspace is contained in the
within-child invariant subspace. Its deleted-coordinate factor is at most
sqrt(m), obtained by retaining the permanent terms fixing that coordinate.
Every signed random permutation in the induction is an INPUT permutation.
Thus the same signed-pair empirical-type Bellman operator applies, with
ordinary output norms until the terminal comparison. Fixed depth precedes m.

## Actual sign construction and twisted contraction

For fixed S, the block weave remains an actual full signing, with exact
defect identity D_sigma=2(m^2k-sigma x^T W x). Since each off-diagonal
pair occurs twice in D, exponentiating -tD/(2k) produces exactly the raw
kernel exp[-t(a-sigma S_ij b)^2]. It is a Gaussian feature contraction
with a unitary reflection on that edge. Absorbing all incident reflections
in the relevant tensor coordinates preserves its Hilbert norm. The
edge-variable Finner/Cauchy--Schwarz bound gives the product of these norms.
There is no unjustified assertion that the reflected kernel itself is PSD.

The squared vertex norm is the ordinary-permutation orbit kernel of its
SIGNED spectrum. After independent basis averaging, this yields the old
exponent uniformly for EVERY fixed macro signing S. Both polarities are
counted. Only physical diagonal entries are removed at O(N) cap cost;
macrodiagonal blocks are not removed by a leading triangle estimate.

## The seed is not erased by hidden final signs

For U_s=diag(U_1,U_2)H2 g and physical H_s=sqrt(s)U_s^T, the top g is
a physical ROW gauge/permutation. There is no independent final COLUMN
sign. Such a column sign would indeed be fatal to the information claim:
S_ij sigma_i(j)sigma_j(i) would become independent fair macro edge signs,
erasing every S exactly. The proposed construction does not do that.

A physical row decomposes into independent signed-uniform terminal rows.
For a dephased order-q physical terminal Hadamard, the averaged distinct
second moment is zero and the averaged distinct fourth moment is1/(q-3).
The latter is the sum of the all-positive row contribution1/q and the
balanced-row contribution3/[q(q-3)]. For b terminal blocks and m=bq,

    mu4=b q(q-1)(q-2)/(m)_4.

Independence across five fibres gives exactly mu4^5 times the macro K5
sign product. At fixed depth this is of order m^-5, whereas conditioning
on cap at most(U0+epsilon)N^(3/2) changes any bounded expectation by only
O(exp(-c_epsilon m^2)). Thus the seed remains distinguishable in the
successful low-cap law. The epsilon margin is fixed before depth/order.

For m>=7, choose a triple T disjoint from a,b,c,d. The mod-two sum of
K5(T+a+c), K5(T+a+d), K5(T+b+c), K5(T+b+d) is the four-cycle
ac+ad+bc+bd. Four-cycle products determine off-diagonal edge signs modulo
switching and global negation. No diagonal information is claimed.

## What is and is not removed

The theorem removes the claim that independent macro-edge averaging, or
complete seed-law erasure, is necessary for the old strict upper compiler.
It gives a positive cycle-sensitive actual construction for EVERY fixed seed.

It does NOT give Q(parent)<=Q(seed) times the required normalization. The
present graph norm bound discards every reflection before estimating the
cap, so its quantitative upper bound is independent of Q(S). A nonzero
polynomially small observable is not a leading n^(3/2) seed-value transfer.
All-order principal restriction may also lose the full seed interface;
the exact fixed-interface theorem is at the constructed orders mk.

The remaining useful question is a genuinely phase-sensitive estimate of
this twisted contraction that yields a smaller parent cap for favorable
seeds. Merely renaming that estimate, or retaining additional cycle bits
without an objective inequality, does not complete the convergence argument.
