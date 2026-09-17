# Full-sign regularization from exchangeable vertex information

2026-09-17. **Verified:** director derivation and independent reconstruction
by all three researchers, including the read-two inequality and tie rule. This
strengthens the fixed-pattern Stein regularizer by including the entire
random new child in the information argument. The random variables are
ordinary exact signs; Gaussian fields are unnecessary for the entropy
conclusion. No external novelty claim is established.

## 1. An overlapping-coordinate information lemma

Let V_1,...,V_d be independent random variables and U any finite-valued
random variable coupled to them (additional randomness is allowed).
Let S_1,...,S_k be coordinate subsets, each coordinate appearing in at
most b of them. Then

```math
 \sum_j I(U;V_{S_j})\le b I(U;V_1,\ldots,V_d)\le b H(U). \tag{1}
```

Proof: expand I(U;V_S) in increasing coordinate order. Independence
makes I(U;V_i|V_{S intersect [i-1]}) no larger than
I(U;V_i|V_[i-1]): both equal H(V_i) minus the appropriate conditional
entropy. Sum and use the multiplicity bound. This argument also works
for continuous V by the chain rule and nonnegative conditional mutual
information, rather than subtracting possibly infinite differential
entropies. Only discrete fair signs are needed below.

## 2. Exact increment of the randomly extended signing

Fix ANY old full signing A of order m. For k>=0, form W_k by adjoining
k vertices and choosing independently every edge with at least one new
endpoint. Write P_k=E Q(W_k), where Q is the absolute Boolean cap.
Put N=m+k. Then

```math
 0\le P_{k+1}-P_k
 \le \sqrt N+2\sqrt{\frac{N(N+1)\log2}{k+1}}.        \tag{2}
```

To prove the upper bound, select Z uniformly among all absolute maximizing
spin words of W_(k+1). This tie rule is equivariant under permutations
of the k+1 new vertices. Let R_j be the incident random edge row of
new vertex j; it has N independent sign entries. Every random edge
belongs to at most two of these rows, so (1) and exchangeability imply

```math
 I(Z;R_j)\le 2H(Z)/(k+1)\le 2(N+1)\log2/(k+1).       \tag{3}
```

For each fixed word z, |R_j dot z_-j| has mean
mu_N=E|sum_(i=1)^N epsilon_i|<=sqrt N, independently of z.
Its centered MGF is bounded by exp(Nt^2/2), by bounded differences
with coordinate changes at most2. The entropy variational inequality
therefore gives

```math
 E|R_j\cdot Z_{-j}|\le\mu_N+\sqrt{2N I(Z;R_j)}.     \tag{4}
```

For clarity, compare the joint law of (Z,R_j) with its product
marginals; the latter has the same uniform MGF bound for EVERY z.
The KL divergence is I(Z;R_j). Optimize the elementary bound
(I+Nt^2/2)/t. No conditional product law for the actual optimizer is
assumed. Pointwise deletion gives

```math
 Q(W_{k+1})-Q(W_{k+1}\setminus j)
    \le |R_j\cdot Z_{-j}|.
```

The deleted graph has law W_k, proving (2). Nonnegativity follows
from principal-restriction cap monotonicity. Internal edges of the
new child are INCLUDED in this calculation; no Q(D) error is paid.

## 3. A single conditional increment controls every near-level set

For a fixed full signing W of order N define

```math
 \Delta(W)=E_h\max_z\{|H_W(z)|+|h\cdot z|\}-Q(W),
```

where h is a fresh uniform sign row. This is exactly the conditional
mean increase in the absolute cap on adding ONE vertex: global reversal
of all old spins leaves H_W unchanged and reverses the row energy.
Thus E Delta(W_k)=P_(k+1)-P_k, with no cancellation assumption.

Let E_W(T)={z:Q(W)-|H_W(z)|<=T}. It is antipodal. Its Bernoulli
width b(C)=E_h max_(z in C) h dot z satisfies, for every T>=0,

```math
 b(E_W(T))\le T+\Delta(W).                           \tag{5}
```

Indeed restrict the fresh-row maximum to this code, then average.
Crucially the SAME Delta controls ALL T simultaneously.

If C is a nonempty Boolean code and d its VC dimension, then d<=b(C).
Choose a shattered d-set; realize the d corresponding random signs by
a codeword depending only on them. The remaining fair signs contribute
zero in expectation, proving this inequality. Sauer's induction gives
|C|<=sum_(j=0)^d binom(N,j)<=(eN/d)^d for d>=1.
Since x log(eN/x) is increasing on[0,N],

```math
 \log|C|\le B\log(eN/B),\qquad B=\min\{N,b(C)\}.     \tag{6}
```

The convention at B=0 is zero. A code of VC dimension zero is a
singleton; the formula also covers that case. All logarithms are natural.

## 4. Same-order actual-sign theorem

There is a universal constant C such that for EVERY full signing A_N
and every 1<=q<=N/2, there is a full signing W_N obtained by keeping
an arbitrary principal (N-q)-block and randomly refilling all remaining
edges, with

```math
 Q(W_N)\le Q(A_N)+C N\sqrt q,\qquad
 \Delta(W_N)\le C N/\sqrt q.                         \tag{7}
```

It follows SIMULTANEOUSLY for every T>=0 that

```math
 b(E_W(T))\le T+C N/\sqrt q,
 \log|E_W(T)|\le B_T\log(eN/B_T),\quad
 B_T=\min\{N,T+C N/\sqrt q\}.                       \tag{8}
```

For the proof, (2) and Markov put Delta below eight times its mean
with probability at least7/8. The rewritten edge count is
d=(N-q)q+binom(q,2)<=Nq. For each spin, its random rewritten energy
has subGaussian proxy d. Union over 2^N spins and both signs puts its
absolute maximum below C N sqrt(q) with probability at least7/8 for
a universal C. These two events intersect. The retained principal cap
is no larger than Q(A_N), by averaging omitted old spins. This proves
(7), and (5)--(6) give (8). At most Nq edges are changed.

For q->infinity, q=o(N), the cap loss is o(N^(3/2)) and EVERY window
T_N=o(N) has subexponential cardinality for this same W_N. In particular
q=floor(N/log N) gives cap loss O(N^(3/2)/sqrt(log N)), and at the
explicit window T=sqrt(N log N), entropy O(sqrt(N)(log N)^(3/2)).
The corresponding constant factors are fixed by (7), not silently
sent to zero. This is a statement at every large order.

Define the class R_(N,q) by Delta(W)<=C N/sqrt(q), independently of
the unknown optimum. Then

```math
 M_N\le\min_{W\in R_{N,q}}Q(W)\le M_N+C N\sqrt q.    \tag{9}
```

This controls the mean cost of ONE random new vertex, not a favorable
linear-density extension or an accumulated recurrence. Iterating (7)
does not give a summable normalized defect for free.

## 5. Combination with favorable frame selection

Choose q0=floor(N^(4/5)(log N)^(2/5)) in (7), and let
T=N/sqrt(q0). Then a regularized actual minimizer has

```math
 Q(W)\le M_N+O(N^{7/5}(\log N)^{1/5}),
 \log|E_W(T)|=O(N^{3/5}(\log N)^{4/5}).              \tag{10}
```

Apply [entropy-frame selection](paper_director_entropy_frame_selection_2026_09_17.md)
to this actual code, for a prospective new child of comparable order and
a mode size k=Theta(sqrt N). Its all-offset restricted-parent
sign/Gaussian comparison costs

```math
 O(N^{7/5}(\log N)^{1/5}).                           \tag{11}
```

Indeed W4=O(N^2 log|E|+N^(5/2)log N); take its fourth root and
multiply by O(N^(3/4)). This balances the regularization cost in (10).
All vertices and signs are genuine, all new-spin words and both
polarities remain in the comparison, and the child offsets are exact.

The stronger codewise-field comparison now improves the balanced combined
cost to O(N^(4/3)(log N)^(1/3)); see
[the later theorem](paper_director_codewise_field_universality_2026_09_17.md).
The fourth-mass result above remains valid and records a distinct route.

Two obligations remain: control of the Gaussian parent VALUE and of
all old words outside E_W(T). Neither is removed by the information
bound or by fourth-order universality. The reported original interval
and convergence status are therefore unchanged.

## 6. Proof provenance

The initial localization/Stein construction regularized fixed new-spin
patterns. The Bernoulli researcher then combined exchangeable-sample
mutual information with exact new-spin maximization to remove that
pattern union, but initially paid the entire new child separately.
The director's overlapping-row argument (1)--(4) incorporates its
internal edges and yields the unrestricted q=o(N) result here.
Classical ingredients are entropy selection, Sauer's lemma, and
elementary concentration; their combination and exact signing
consequences, not the ingredients themselves, are the campaign result.
