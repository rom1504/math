# Favorable high-rank directions need not themselves be regular

Date: 2026-09-07. Status: exact asymptotic construction, pending independent
audit. This is a diagnostic for one proposed route, not an obstruction to
all favorable high-rank sign recovery.

There are actual bounded-cap full signings A_n and signed contractions G_n
of rank o(n) such that Q(A_n-5sqrt(n)offdiag(G_n)) is strictly below Q(A_n)
by a fixed multiple of n^(3/2), although the correction has entries tending
to infinity on coherent blocks. Thus favorable descent alone cannot force
the given direction to have small entries or low maximum leverage.

Take n=2^(4k), and let H_n be the symmetric Sylvester Hadamard matrix. Put
F_n=H_n-diag(H_n), so

    Q(F_n)<=n^(3/2)/2+n/2.

Let u=1/sqrt(n) be the unit constant vector. Independently round the
entrywise-feasible mean

    Y=(1-4/sqrt(n))F_n+(4/sqrt(n))offdiag(J_n)

to actual hollow signs A_n. For n>=16 these are convex combinations of
signs. Total variance is O(n^(3/2)), so Bernstein over all Boolean spins
gives a realization with

    Q(A_n-Y)=O(n^(5/4)).

Consequently Q(A_n)<= (5/2+o(1))n^(3/2). On the constant spin vector,

    H_Y(1)>=2n^(3/2)-Q(F_n)-O(sqrt(n)),

so Q(A_n)>=(3/2-o(1))n^(3/2).

Partition the vertices into r=2^(3k)=n^(3/4) blocks of size s=2^k=n^(1/4).
Let P be the orthogonal projector onto vectors constant on each block.
Then Pu=u, rank(P)=r, and P_ij=1/s inside a block, zero across blocks.
Define

    G=(4/5)uu^T+(1/5)P.

This is PSD, has operator norm1 and rank r=o(n). The proposed perturbation
is

    Delta=5sqrt(n)G=(4/sqrt(n))J_n+sqrt(n)P.

After subtracting its off-diagonal part,

    A_n-offdiag(Delta)
      =(1-4/sqrt(n))F_n-sqrt(n)offdiag(P)+E_n,
    Q(E_n)=O(n^(5/4)).

Since 0<=x^T P x<=n and Tr(P)=r,

    Q(sqrt(n)offdiag(P))<=n^(3/2)/2.

Therefore

    Q(A_n-offdiag(Delta))<=(1+o(1))n^(3/2),

a strict favorable descent of at least (1/2-o(1))n^(3/2). Yet inside
each coherent block Delta_ij=4/sqrt(n)+n^(1/4), which diverges; no
entrywise feasibility or small maximum-leverage conclusion is possible.

The important limitation is that this example ALSO has the obvious
implementable rank-one descent Delta_0=(4/sqrt(n))J_n, with nuclear norm
4sqrt(n)=o(n), and target cap at most (1/2+o(1))n^(3/2). Hence it does not
refute a theorem that REPLACES an arbitrary favorable high-rank direction
by another implementable one. It only rules out deducing regularity of
the given direction solely from the favorable cap inequality. The A_n
here are deliberately not near-minimizers.

## Follow-up: what a replacement argument still has to prove

For any op-bounded rank-r G with r=o(n), Delta=sqrt(n)G satisfies
||Delta||_F²<=nr=o(n²). Thus, if epsilon_n tends to zero sufficiently
slowly, the exceptional edge set {|Delta_ij|>epsilon_n} has o(n²) edges.
The small-entry part has a uniform rank-one magnitude envelope of trace
n epsilon_n=o(n), hence is implementable by target contraction and
independent rounding. This does NOT license deleting the exceptional
large-entry correction: deleting an arbitrary o(n²)-edge graph can change
cap by order n^(3/2), even though random refilling of that graph has
o(n^(3/2)) noise. Low-rank geometry must enter that missing deletion step.

For literal coherent cluster frames, the exceptional correction often looks
like a diagonal quadratic term in the coarse block magnetizations. It is
tempting to discard it because actual within-block edge energy is only
O(n²/r). This argument is invalid. Fix the internal spin pattern in each
block and reverse entire blocks. The cross-block energy becomes a coarse
hollow quadratic K(y), whereas the diagonal magnetization correction is a
constant shift on this reversal orbit. Such a shift can center the interval
[min K,max K], reducing its absolute cap from W(K)+|I(K)| to W(K).
Therefore removing the coarse diagonal can lose a genuine midpoint benefit.
This connects the proposed replacement operation to the unresolved
width-versus-absolute-cap problem rather than resolving it by clipping.

No general higher-rank favorable replacement theorem, cap-safe exceptional
edge deletion, or midpoint repair is proved by these observations.
