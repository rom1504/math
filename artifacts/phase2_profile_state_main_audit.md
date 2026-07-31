# Main-agent audit of the six-vertex switching profile

Status: proved reductions and limitations, followed by a scalable universal
falsifier. This note is part of the second sustained campaign.

## 1. The state and its exact landing implication

For an order-`n` signing `A`, let `phi_6(A)` be the 16 counts of induced
switching/permutation/global-negation classes on vertex sets of sizes four,
five, and six. There are respectively 2, 4, and 10 such classes. Every
coordinate is an integer between zero and `binom(n,k)`, so the number of
possible states is polynomial in `n` (with a large fixed exponent).

Choose, without reference to `M_n`, one canonical representative of every
realizable state and call the resulting family `F_n`. If

~~~math
\phi_6(A)=\phi_6(B)
\quad\Longrightarrow\quad
|\operatorname{cap}(A)-\operatorname{cap}(B)|=o(n^{3/2})               \tag{Q1}
~~~

uniformly over all order-`n` pairs, then the representative sharing the state
of an exact minimizer has cap `M_n+o(n^(3/2))`. On the known
`Theta(n^(3/2))` scale this is equivalent to

~~~math
\min_{S\in F_n}\operatorname{cap}(S)^{2/3}-M_n^{2/3}=o(n).             \tag{Q2}
~~~

Thus (Q1) really would prove structured landing for a polynomial-cardinality,
noncircular family. It would not supply the separate composition/update
clause, and a canonical representative defined by lexicographic search is not
yet an efficient generator. Section 4 records that (Q1) is false.

## 2. What the profile determines algebraically

Any switching-invariant signed monomial supported on at most six vertices,
and invariant under global negation, has a sum over embeddings determined by
`phi_6`. In particular, equal profiles imply equal even spectral traces

~~~math
\operatorname{tr}A^2,\quad\operatorname{tr}A^4,
\quad\operatorname{tr}A^6.                            \tag{Q3}
~~~

For example, if `N_4^+(A)` is the number of four-sets on which the three
Hamilton-cycle products are all positive, direct closed-walk enumeration
gives

~~~math
\operatorname{tr}A^4
=n(n-1)(2n-3)
 +8\left(4N_4^+(A)-{n\choose4}\right).               \tag{Q4}
~~~

The same profile also determines the first three even moments of the Boolean
energy

~~~math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j:
\qquad
\mathbb E H_A^2,\quad\mathbb E H_A^4,
\quad\mathbb E H_A^6.                                \tag{Q5}
~~~

To see (Q5), expand `H_A^(2r)`. A nonzero Walsh expectation is an even-degree
multigraph with `2r` edge occurrences; after isolated vertices are removed,
it uses at most `2r<=6` vertices. Its signed embedding sum is one of the local
invariants recorded by the profile.

These facts explain why the state is a strong finite predictor. They also
identify its limitation: it fixes finitely many polynomial moments, while cap
is an exponentially rare maximum over `2^(n-1)` projective states.

## 3. Fixed-moment norm bounds miss the required scale

The strongest immediate spectral consequence of (Q3) is

~~~math
\|A\|_{op}\le(\operatorname{tr}A^6)^{1/6}.
~~~

Even at the minimal trace scale `tr(A^6)=Theta(n^4)`, this and the ordinary
spectral cap estimate give only

~~~math
\operatorname{cap}(A)\le {n\over2}\|A\|_{op}
=O(n^{5/3}),                                         \tag{Q6}
~~~

not `O(n^(3/2))`. Finite energy moments give no upper bound on a maximum at
the desired scale without a tail/entropy theorem. The existing full-spectrum
cap collisions inside the balanced `PC(26)` restrictions prove that adding
all spectral moments does not repair this at finite order.

Therefore a proof of (Q1) cannot consist only of converting the 16 counts to
spectral traces or low energy moments. It needs one of:

1. an exact profile-preserving transformation coupling the full Boolean
   energy landscapes;
2. a uniform high-moment/tail theorem derived from the local counts; or
3. a rigidity theorem showing that every realizable low-cap profile has only
   `o(n^(3/2))` cap variation.

The second option must explicitly avoid the earlier common-active-face
entropy loss; controlling six moments alone is the same obstruction in new
language.

## 4. Exact falsifier and current classification

The falsification track found exact order-ten signings `A` and `B` with the
same oriented switching/permutation restriction counts at every order at
most six, but

~~~math
\operatorname{cap}(A)=19,\qquad \operatorname{cap}(B)=21.              \tag{Q7}
~~~

Replace every vertex by `L` positive twins. The induced profile of a set of
at most six vertices is a fixed transform of the oriented base profiles,
indexed by its occupancy composition, so the two blowups still have equal
profiles. Exact positive-definiteness certificates for `A+5I` and `B+5I`
and coordinatewise convexity give, for every `L>=3`,

~~~math
\operatorname{cap}(T_L(A))=24L^2-5L,\qquad
\operatorname{cap}(T_L(B))=26L^2-5L.                \tag{Q8}
~~~

The full proof and reproducible matrices are in
`phase2b_phi6_collision_report.md`. I independently reconstructed the pair,
checked every spin, checked the common profile, and enumerated the reduced
blowup energies for `L=2,3`. Thus the universal stability statement (Q1) is
false by a `Theta(N^2)` gap.

There is also a stronger correct-scale falsifier. With the symmetric
Sylvester Hadamard matrix `H_k` of order `k=4^r`, replace every base entry by
that common microblock and fill each diagonal macroblock with
`H_k-diag(H_k)`. The oriented profiles still agree. An exact rational
positive-semidefinite dual certificate and a Boolean Hadamard eigenvector
give

~~~math
\operatorname{cap}(S_A(k))\le {207\over8}k^{3/2},\qquad
\operatorname{cap}(S_B(k))\ge 26k^{3/2}.             \tag{Q9}
~~~

Thus both families have the correct `O(N^(3/2))` cap scale and their caps
differ by at least `N^(3/2)/(8*10^(3/2))`. I independently checked the SDP
normalization and reproduced its exact positive leading principal minors.
This is a scalable falsifier of fixed-profile control on the entire
correct-scale class, not merely a quadratic spike.

The certified constants in (Q9) remain substantially above the best known
upper constant for `M_N`, so the construction does not by itself produce two
near-minimizers. It does show that choosing an arbitrary canonical
representative of a profile cannot be justified by the profile alone.
Defining the representative to minimize cap inside its fiber merely hides
the original optimization and is not a reduction.

The exhaustive order-eight and Paley-half successes are now explained as
finite predictor behavior, not evidence for a universal mechanism. Fixed
`phi_6` is inactive as a standalone landing state. A low-cap-only repair
would need an independently checkable *near-optimal* certificate and a
composition rule; the phrase "among near-minimizers" is circular because
membership already invokes the objective being approximated.
