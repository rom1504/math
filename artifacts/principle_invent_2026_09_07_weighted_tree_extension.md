# Universal bounded-profile extension: the two sensitive tree modules

Date: 2026-09-07. Status: **independently audited PASS by the director**;
see `principle_director_weighted_extension_audit_2026_09_07.md`, including
the approximate-diagonal scalar step and fixed numerical dependency.
This note supplies the profile-sensitive changes in
the already audited marked-response lower proof. It does not prove
favorable flatification, equality of optimized profile and sign values,
or convergence of the original extremal sequence.

Let W be real, symmetric and hollow of order n, and set B=W/sqrt(n).
Assume, for fixed finite K,L,

    max_ij |W_ij| <= K,       ||B||op <= L,
    max_i |sum_j B_ij^2 - 1| = eta_n -> 0.

Write Q=B^2. The bilinear Boolean norm satisfies
beta(W)<=n||W||op<=L n^(3/2). All diagrams below are fixed before n tends
to infinity; their constants may depend on their sizes and on K,L.

## 1. A weighted single-edge cancellation rule

Suppose a normalized finite diagram contains an edge occurring exactly
once between two labels that are both summed. Other edges may be squared
or have higher fixed multiplicities. If no other formal edge uses that
same pair of vertices, fixing every other label makes all remaining
factors incident to the two chosen vertices into bounded unary functions.
Consequently the normalized average is bounded by

    C_diagram,K beta(W)/n^2 = O(n^(-1/2)).

This uses the actual W entry on the chosen edge, not its sign. Bounded
real unary tests are allowed in beta by scaling and multilinear cube
rounding. Squared W entries on OTHER edges must be retained, but they
do not prevent this cancellation. Collision exclusions cost O(1/n) in
each bounded normalized fixed diagram.

## 2. Exact old-input Gram without constant-amplitude cancellation

Use the fully injective old input trees of
`fresh_limit_injective_input_gram_2026_09_05.md`: a tree R with odd vertex
count d, even degree at its marked root, and odd degree at every other
marked vertex. Its coefficient is the product of its d-1 B entries,
divided by the square root of its rooted automorphism number.

Different marked degrees have exactly zero covariance. Equal degrees
force a bijection of the two marked vertex sets. Every resulting edge
has multiplicity at most two. Retain each common edge with weight B^2;
the remaining edges form the parity graph.

At equal roots the parity graph is Eulerian. If nonempty, it contains a
free-free edge and Section 1 applies. If empty, the two trees are rooted
isomorphic, and the common graph is a doubled tree. Unrestricted leaf
summation gives 1+O_R(eta_n), since every row square sum is 1+O(eta_n).
Removing label collisions costs O_R,K(1/n). The automorphism factors
cancel exactly as before. Thus diagonal Gram errors are
O(eta_n+n^(-1/2)).

At distinct roots i,j there are d-2 free labels and normalization
n^(-(d-1)), giving a prefactor 1/n times a bounded weighted diagram.
The parity graph has odd degree at i,j and even degree elsewhere.
In addition its NUMBER OF EDGES IS EVEN: total edge occurrences are
2(d-1), and deleting common pairs preserves edge-count parity. In
particular the apparent exception consisting solely of the edge ij
is impossible.

There is therefore a parity edge other than ij. In the squared
Frobenius norm, take two copies with only i,j glued. That edge remains
single, because it meets a free vertex that is not glued between copies.
All doubled weights remain bounded and Section 1 applies. With the 1/n
matrix prefactor the Frobenius square is O(n^(-1/2)); consequently the
off-diagonal operator error is O(n^(-1/4)). Finite-pattern collision
errors have entries O(n^(-2)), hence operator norm O(1/n). Altogether,

    ||Cov(V_R,V_U)-1_(R iso U) I||op
        <= C(eta_n+n^(-1/4)).                         (1)

The output-root collision covariance proof in Section 6 of that source
uses only coefficient magnitudes and label counts. It still gives
O(1/n) in operator norm for bounded K. Hence

    ||Cov(X_T,X_U)-1_(T iso U) Q||op = o(1).          (2)

The same parity-even observation removes the sole-root-edge case in
the raw even response input Gram argument: every constituent old output
tree has an odd number of edges, while an even response has an even
number of such factors. The total edge-occurrence count is even.
One must still use that source's leading-label classification before
applying the single-edge rule; this observation does not license deleting
arbitrary squared factors in unrestricted collision diagrams.

## 3. Common-tree lemma for the exceptional old-input/new-forest pairing

Let V be an input tree on q marked vertices, rooted at the marked vertex
a. Let P be a tree on those same q marked vertices plus an external
unmarked root j, whose degree is the odd number k>=3. V has q-1 edges;
P has q edges. All marked P vertices have odd degree; V has even degree
at a and odd degree elsewhere. The identification is a marked-set
bijection, so j is absent from V, and a!=j.

Consider a covariance pattern having NO free-free parity edge, where
free excludes a,j. Its parity graph can only consist of p paths
a-u_r-j and an optional edge a-j. At j every P edge survives, because
j is absent from V. Thus k=p+e, with e in {0,1}. The total number of
edge occurrences, 2q-1, is odd, so its parity graph has an odd edge
count. Since that count is 2p+e, necessarily e=1 and p=k-1>=2.

Let F be the common-edge graph on the q marked vertices. It is a forest,
as a subgraph of V. It has

    |E(F)| = q-1-p,

and hence p+1 components. The edges a-j and j-u_r are P-only. The edge
a-u_r cannot also be in P, since that would create the P triangle
a-j-u_r-a. Thus every a-u_r is V-only.

The p+1 anchors a,u_1,...,u_p lie in distinct F components:

* An F path from a to u_r would form a cycle in V with a-u_r.
* An F path from u_r to u_s would form a cycle in P with u_r-j-u_s.

Since the component count equals the anchor count, EVERY F component is
a tree rooted at exactly one anchor. There are no unrooted components,
no component meeting two anchors, and no residual common cycle.

## 4. Positive unary weights control the exceptional covariance

First allow unrestricted labels. Sum the nonanchor vertices in every
F component. A component rooted at v produces a nonnegative function
f(v), obtained by summing a product of B-edge squares over a rooted
tree. If row square sums are bounded by r, a tree with h edges has

    0 <= f(v) <= r^h.

Thus no exact row normalization is required for the following smallness
bound. Let f_0 be the component weight rooted at a and f_r the one
rooted at u_r. The remaining pattern matrix is precisely

    diag(f_0) [ B circ Q_1 circ ... circ Q_p ],
    Q_r = B diag(f_r) B,                 p>=2.       (3)

Each Q_r is PSD, has bounded diagonal, and has bounded operator norm
under the fixed K,L hypotheses. In particular its entries and the
Euclidean norms of all rows are uniformly bounded. Apply Cauchy--Schwarz
to two of the p factors and bound the others entrywise. Since
max|B_aj|<=K/sqrt(n), every absolute row sum in (3) is O(n^(-1/2)).
Alternatively direct entrywise squaring and the same two-factor argument
gives Frobenius norm O(1). Multiplication by diag(f_0) is harmless.

For exact row square sums one, all f_r=1 and (3) reduces to the original
B circ Q^(circ p). For approximate row normalization, it is unnecessary
to approximate (3) by that expression; the positive unary factorization
already proves the needed O(1) Frobenius estimate.

Finite injectivity exclusions contribute O(n^(-3/2)) per matrix entry:
the original fixed-root covariance has nominal prefactor n^(-1/2), and
each extra collision loses one free label. Their Frobenius norm is
O(n^(-1/2)). Patterns with a free-free parity edge also have entries
O(1/n), by Section 1 and the nominal n^(-1/2) prefactor. Therefore

    ||Cov(V_T,P)||F = O(1).                         (4)

This preserves the all-root mixed covariance estimate used in the
root-not-hit full-contraction argument, not merely its diagonal average.

## 5. Why the other weighted modules keep their hypotheses

The finite local old-tree moment proof and its marked bridge identity
use the single-edge rule above or doubled-tree leaf summations. The
latter acquire only O_T(eta_n) errors. The two-output-root classification
similarly leaves actual cross factors Q_ij, not an averaged variance
profile. Doubled rooted attachments contribute 1+O_T(eta_n).

Global and fixed-root tree flattenings use only fixed operator norm,
bounded row square sums, maxentry O(n^(-1/2)), and their combinatorial
edge counts. Replacing the sign entry bound by K changes their fixed
constants only. Proper forest cuts and partial nonlinear matchings
therefore retain their smallness. Whole pairings use the genuine (2),
which preserves the Q-Schur covariance main term up to o(n) nuclear norm.

The root-hit doubled diagram retains a single edge at the new root:
its two P copies have disjoint leading free labels, while at most their
common root-to-a edge can double. The other edges at the P root occur
only once even with weighted amplitudes. Section 1 gives the same
E||J||F^2=O(sqrt(n))+O(1)=o(n) bound. Squared weights on the remaining
formal edges are retained as bounded unary tests, not erased.

The derivative, collision, sign-input replacement, and fixed-polynomial
approximation steps use the same bounds. Their constants are fixed
before n tends to infinity. Uniform eta_n=o(1) contributes additional
o(1) finite-tree moment and normalization errors; it need not have a
specific rate. The positive unary factorization in Section 4 avoids
amplifying eta_n in the sensitive all-root mixed covariance bound.

## 6. Consequence and boundary

Subject to the director's audit of this dependency transfer, the existing
finite marked-response certificate gives for every such W_n

    liminf_n width(W_n)/n^(3/2) >=
      0.433322111664080753415812928897579346558634648033693413106996.

Here width(W)=(max_x H_W(x)-min_x H_W(x))/2 and H_W(x)=x^T W x/2.
The numerical certificate is unchanged. The only theorem extension is
from full signs to bounded-entry, uniformly approximately row-regular
real matrices under a fixed normalized operator bound.

Removing that operator assumption is a SEPARATE step. The director's
variance-completion construction after spectral deletion appears to
provide it; its details are recorded separately rather than being
smuggled into this fixed-operator statement.

## 7. Explicit approximate-diagonal repair for the scalar variance step

The final variance inequality is robust without symmetrically rescaling
B itself. Put D=diag(Q), Qhat=D^(-1/2) Q D^(-1/2). At fixed L,

    ||Q-Qhat||op=O_L(eta_n).

For any fixed finite odd-Schur mixture R=sum w_p Q^(circ p), of total
weight one, put Rhat=sum w_p Qhat^(circ p). Telescoping each Schur power
and using the bounded diagonal Schur-multiplier norms gives

    ||R-Rhat||op=O_L,degree(eta_n).

The correlation theorem applies to Qhat and yields
Tr[Qhat sqrt(Rhat)]>=n. Since Tr sqrt(Rhat)<=n, positivity gives

    Tr[Q sqrt(Rhat)] >= n-O_L(eta_n)n.

Every column b_i of B has norm sqrt(D_ii)<=sqrt(1+eta_n). Hence

    sum_i sqrt((B Rhat B)_ii)
      >= (1+eta_n)^(-1/2) Tr[Q sqrt(Rhat)]
      >= n-O_L(eta_n)n.

Scalar square-root continuity and the operator estimate for R-Rhat
then give

    n^(-1) sum_i sqrt((B R B)_ii) >= 1-o(1).

This makes the approximate-row hypothesis in the fixed-operator theorem
explicit at the only matrix-correlation step that originally had an exact
unit diagonal. General total Schur weight tau^2 follows by scaling.

## 8. Audited variance completion and universal consequence

I independently checked
`principle_director_variance_completion_2026_09_07.md`, including its exact
deficit formula, positivity, simultaneous concentration, and objective
control. Let W_n now have maxentry bounded by fixed K and row square
sums n-1+o(n), uniformly, with NO operator hypothesis. If its normalized
width is potentially below the displayed constant, then Q(W)<=2 width(W)
supplies the bounded-cap hypothesis for spectral deletion.

At fixed epsilon, the director's construction produces a principal order
k>=(1-epsilon)n and a bounded-entry weighted V_k with fixed normalized
operator bound, uniformly approximate row square sums k-1+o(k), and

    width(V_k) <= (1+O_K(epsilon)+o(1)) width(W_n)
                    +O_K(sqrt(epsilon)n^(3/2)).

Apply the fixed-operator result to V_k, then take n to infinity, and only
then epsilon to zero. The resulting universal lower bound is

    liminf_n width(W_n)/n^(3/2) >=
      0.433322111664080753415812928897579346558634648033693413106996.

This is a substantive extension of the marked lower theorem to EVERY
bounded-amplitude, uniformly row-regular real variance profile. It rules
out any scalable sub-c_* counterexample in that class. It still does not
show optimized profile universality above c_* and does not manufacture
full signings from a weighted profile.
