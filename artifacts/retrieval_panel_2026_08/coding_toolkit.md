# Coding theory / covering-radius retrieval toolkit

Status: retrieval packet checked through 2026-08-15; not a proposed solution.

Prepared 2026-08-15. This is a theorem-level literature packet, not a proposed attack. It uses the problem's native binary-code and signed-complete-graph objects, and it separates results that quantify over **every ambient word/coset** from packing results that constrain only **pairs of codewords**.

## 0. Executive classification

The decisive distinction is not LP versus SDP. It is what is being quantified.

| Result family | Native data | Arbitrary ambient root/coset? | What it can actually certify |
|---|---|---:|---|
| Sphere-covering bound | code size and Hamming-ball volume | yes | a lower bound on covering radius |
| Signed-graph switching/frustration | one signed graph modulo switching | yes, exactly | the minimum weight in that coset; maximizing gives the covering radius |
| Delsarte external distance | outer distribution / dual support | yes | `rho(C) <= s(C)` for every coset |
| Tietavainen dual-distance bound | Krawtchouk polynomial and dual distance | yes | an all-coset upper bound on radius |
| Bazzi small-code estimates | distributions of coset weights | almost all initially; all after an existential extension | almost-covering, then covering after adding a small auxiliary space |
| Completely regular codes | equitable distance partition / coset graph | yes, conditional on complete regularity | exact reduction of all cosets to `rho+1` layers |
| Graham-Sloane amalgamation | syndromes and a normal coordinate | yes | deterministic length/radius transfer |
| Schrijver Terwilliger SDP | triples rooted at a codeword after averaging | no | packing/minimum-distance upper bounds |
| Gijswijt-Polak covering SDP | rooted covering inequalities plus triple moments | yes | lower bounds for the size of covering codes |
| Moving-projection binary/spherical bounds (OpenAI 2026) | code-point-dependent stabilizer representations | no | stronger packing exponents for every fixed relative distance/inner product |
| Exact spherical SDP and AI-found kissing configurations | pair/triple separation constraints | no | packing optimality or improved packing constructions |

The sources below use the following tags:

- **ALL-COSET**: the theorem itself has a universal quantifier over ambient words, cosets, or syndromes.
- **ALMOST-ALL**: it controls the measure of bad ambient words but not necessarily the deepest word.
- **CONDITIONAL ALL-COSET**: it becomes exact after a structural hypothesis such as complete regularity or normality is proved.
- **PACKING ONLY**: it concerns pairwise distance or inner product and has no direct deep-hole conclusion.

## 1. Native dictionary for cut codes and the augmented problem

Let `Gamma=(V,E)` be a graph, `N=|E|`, and identify an edge signing with `x in F_2^E`; write

`C*(Gamma) = { delta(S) : S subseteq V }`

for its binary cut code. Switching the signed graph at `S` sends `x` to `x+delta(S)`. Hence:

- switching classes of signings are exactly the cosets of `C*(Gamma)`;
- a coset leader is a minimum-frustration representative;
- `min_S wt(x+delta(S))` is the frustration index of the signing `x`;
- `rho(C*(Gamma))` is the maximum frustration index over all signings.

For a complete graph, put `a_ij=(-1)^{x_ij}` and encode a switch by `s_i in {+1,-1}`. Then

`wt(x+delta(S)) = N/2 - (1/2) sum_{i<j} a_ij s_i s_j`,

so the unaugmented distance is

`d(x,C*) = N/2 - (1/2) max_s sum_{i<j} a_ij s_i s_j`.

For the augmented code

`C^pm = C*(K_n) + <1>`,

the complement of every switched signing is also allowed. Therefore the exact native objective is

`d(x,C^pm) = N/2 - (1/2) max_s |sum_{i<j} a_ij s_i s_j|`.                         (A)

Thus augmentation identifies a switching class with its global sign reversal. This absolute value is not a cosmetic modification: the all-negative class, which is the unique deepest switching class for the ordinary cut code of `K_n`, belongs to `C^pm` and has augmented distance zero.

For `n>=4`, with `N=binom(n,2)`:

- `C*(K_n)` has dimension `n-1`; its dual is the cycle space and has minimum distance `3`.
- `C^pm` has dimension `n`; its dual is the even-weight part of the cycle space and has minimum distance `4`, witnessed by a 4-cycle.
- The sphere-covering inequality gives

  `2^n sum_{j=0}^rho binom(N,j) >= 2^N`,

  and consequently

  `rho(C^pm) >= N H_2^{-1}(1-n/N)
             = N/2 - sqrt(N n ln(2)/2) + lower-order terms`.                     (B)

- Tietavainen's dual-distance theorem, applied only through `d((C^pm)^perp)=4`, gives the universal but structure-blind upper bound

  `rho(C^pm) <= xi_2^N = (N-sqrt(N))/2` up to the necessary integer rounding,    (C)

  where `xi_2^N` is the smaller zero of the binary degree-2 Krawtchouk polynomial.

Equations (B)-(C) place the deficit `N/2-rho` generically between order `sqrt(N)=Theta(n)` and order `sqrt(Nn)=Theta(n^(3/2))`. Neither endpoint uses the special higher-order cut/cocycle structure. Petersdorf's exact ordinary-cut result does **not** determine (A).

## 2. Primary sources and theorem-level extraction

### 2.1 Signed graphs, switching classes, and highly symmetric cut codes

#### 1. Patrick Sole and Thomas Zaslavsky, “A Coding Approach to Signed Graphs” (1994)

- Primary source: [author-hosted paper](https://people.math.binghamton.edu/zaslav/Tpapers/cas.sidma1994.pdf); [DOI](https://doi.org/10.1137/S0895480189174374).
- Status: **ALL-COSET**, for the ordinary cut code.
- Normalization and theorem. For a graph with `m` edges, `n` vertices, and `c` connected components, the cut code has length `m` and dimension `n-c`; its dual is the cycle code. Lemmas 1-2 identify switching classes with cut-code cosets and identify frustration with minimum coset weight. In particular the maximum frustration `D(Gamma)` equals `rho(C*(Gamma))`.
- Quantitative result. Their Theorem 1 applies the covering-volume inequality to obtain

  `D(Gamma) >= m/2 - sqrt(m(n-c) ln(2)/2)`

  in its displayed elementary form. The underlying entropy statement is `D >= m H_2^{-1}(1-(n-c)/m)` when the radius is at most `m/2`. They also derive `D(K_{t,t})=t^2/2-Theta(t^(3/2))` and record Petersdorf's exact formula

  `D(K_n)=floor((n-1)^2/4)`.
- Proof mechanism. Translate switching to addition of cut vectors, then apply the binary sphere-covering bound; special graph families use discrepancy/probabilistic estimates.
- Relevance boundary. This is the exact dictionary needed for coset leaders, but all of its signed-graph radius statements concern `C*`, not `C*+<1>`. Petersdorf's extremizer is the all-negative signing and is annihilated by augmentation.

#### 2. Maximilien Gadouleau and Huiying Zeng, “On the Maximum and Negative Frustration Indices of Graphs” (2026)

- Primary source: [arXiv:2606.11108](https://arxiv.org/abs/2606.11108).
- Status: **ALL-COSET** comparison for ordinary signed graphs; not an augmented-code theorem.
- Exact structural statements used. If `l(Sigma)` is the minimum number of negative edges after switching, then a signature is switching-minimal iff every vertex cut has at most as many negative as positive crossing edges:

  `d^-_Sigma(S) <= d^+_Sigma(S)` for every `S subseteq V`.

  For the all-negative signing, `l(-G)=|E|-MaxCut(G)`. The paper restates/proves the complete-graph extremal result

  `l_max(K_n)=l(-K_n)=floor((n-1)^2/4)`,

  with the maximizing signatures precisely the switching class of `-K_n`.
- Proof mechanism. Switching inequalities turn minimality into an all-cuts local-optimality condition. For the classified graph families the authors combine switching, Menger-type arguments, and packings of negative cycles/triangles.
- Closest no-go theorem. The paper constructs chordal and non-chordal families in which all-negative is nonmaximal, maximal but nonunique, or uniquely maximal, refuting three proposed generalizations. Thus “antibalanced is deepest” is a complete-graph phenomenon, not a robust signed-graph principle.
- Relevance boundary. Global sign reversal is not quotiented out. Its complete-graph theorem therefore gives the deepest coset that disappears in the augmented code, not the next-deepest augmented class.

#### 3. D. E. Taylor, “Regular 2-Graphs” (1977)

- Primary source: [DOI](https://doi.org/10.1112/plms/s3-35.2.257).
- Status: structural switching theory; **not by itself ALL-COSET**.
- Native object. Two-graphs are switching classes of signed complete graphs (equivalently Seidel matrices modulo conjugation by diagonal `+-1` matrices). A two-graph is regular exactly when every edge lies in the same number of negative triangles.
- Spectral theorem. A regular two-graph has a Seidel matrix `S` with two eigenvalues `rho_1>rho_2` satisfying `rho_1 rho_2=1-n`. After switching so a selected root row is all `+1` and writing the remaining block as `B`,

  `(B-rho_1 I)(B-rho_2 I)=-J`, and `B 1=(rho_1+rho_2)1`.

  The descendant at the root is therefore regular and, in the standard integral cases, strongly regular. This is a precise root/stabilizer reduction.
- Proof mechanism. Diagonal switching conjugacy preserves the Seidel spectrum; the two-eigenvalue minimal polynomial, separated into a root row and its orthogonal block, gives the descendant identities.
- Relevance boundary. The spectrum is switching invariant and the theorem controls negative-triangle incidence, not the maximum absolute switching quadratic form in (A). It supplies the natural representation language for moving roots, but not a covering-radius bound.

### 2.2 Krawtchouk transforms, external distance, and true all-coset bounds

#### 4. F. J. MacWilliams, “A Theorem on the Distribution of Weights in a Systematic Code” (1963)

- Primary source: [author-hosted scan](https://www.terpconnect.umd.edu/~abarg/ECC/macwilliams1963.pdf); [DOI](https://doi.org/10.1002/j.1538-7305.1963.tb04003.x).
- Status: pair/weight-distribution foundation; **not ALL-COSET on its own**.
- Theorem. For a `q`-ary linear code `C`,

  `W_{C^perp}(X,Y)=|C|^{-1} W_C(X+(q-1)Y, X-Y)`.

  Equivalently, if `A_i` and `B_j` are the weight distributions of `C` and `C^perp`, then

  `B_j=|C|^{-1} sum_i A_i K_j(i)`.
- Proof mechanism. Additive-character orthogonality/Fourier transform over the finite field.
- Relevance boundary. The ordinary weight enumerator is based at zero. It neither supplies the weight enumerator of every translate nor identifies a deepest coset. A use of “MacWilliams” that never introduces outer distributions or translates remains pair-sensitive.

#### 5. Philippe Delsarte, “Four Fundamental Parameters of a Code and Their Combinatorial Significance” (1973)

- Primary source: [DOI](https://doi.org/10.1016/S0019-9958(73)80007-5).
- Status: **ALL-COSET**.
- Theorem and normalization. The external distance `s(C)` is the number of nonzero entries in the MacWilliams transform of the distance distribution; for a linear code it is the number of distinct nonzero weights occurring in `C^perp`. Delsarte proves

  `rho(C) <= s(C)`.

  Equivalently, every ambient word has distance at most `s(C)` from `C`. In outer-distribution language, the number of nonzero dual eigenspaces supporting the characteristic vector controls the number of distance layers.
- Proof mechanism. Build an annihilator polynomial in the Krawtchouk basis / use the outer distribution module. Evaluating it at a translate forces one of the first `s` distance coordinates to be nonzero.
- Relevance boundary. This is genuinely universal over cosets, but can be very coarse when the dual has many distinct weights. It bounds the worst coset without describing or classifying it.

#### 6. Aimo Tietavainen, “An Upper Bound on the Covering Radius as a Function of the Dual Distance” (1990)

- Primary source: [DOI](https://doi.org/10.1109/18.59949).
- Status: **ALL-COSET**.
- Binary exact root form. Express the orthogonal-array strength as

  `tau=d^perp-1=2k-1+epsilon`, with `epsilon in {0,1}`.

  Then the covering radius is at most the smallest real zero `xi_k^{n-epsilon}` of the degree-`k` binary Krawtchouk polynomial for length `n-epsilon` (with integer rounding implicit for the radius). In the binary length-`N`, `d^perp=4` case, `tau=3`, so `k=2`, `epsilon=0`, and

  `K_2^N(x)=((N-2x)^2-N)/2`, hence `rho <= (N-sqrt(N))/2`.
- Proof mechanism. A Delsarte LP dual polynomial, constructed from Krawtchouk roots, is evaluated against every translated distance distribution. Orthogonality and the vanishing of low dual moments exclude an uncovered root past the stated zero.
- Relevance boundary. It sees only dual distance. Since the augmented complete-graph cut code has fixed `d^perp=4`, this method alone cannot exploit its expanding family of 4-cycles or higher cycle/cocycle incidence.

#### 7. Louay Bazzi, “On the Tightness of Tietavainen's Bound for Distributions with Limited Independence” (2017)

- Primary source: [arXiv:1707.00552](https://arxiv.org/abs/1707.00552).
- Status: theorem-level **NO-GO** for low-degree moment methods; its witness need not be a linear code.
- Theorem. There are absolute constants `k_0,n_0` such that, for `n>=n_0` and

  `k_0 <= k <= n^(1/3)/log^2 n`,

  there exists a `k`-wise independent distribution on `{0,1}^n` whose support has covering radius at least

  `n/2-sqrt(k n)`.

  The dual obstruction is also explicit: a real polynomial of degree at most `k` that is nonpositive at all integral weights in the central interval `|w-n/2|<=sqrt(kn)` must have nonpositive expectation under the binomial distribution.
- Proof mechanism. Exact LP duality between limited-independence distributions and low-degree polynomial witnesses, followed by discrete approximation/quantized Chebyshev analysis.
- Relevance boundary. This shows that a bounded or slowly growing number of Krawtchouk moments cannot generally push an all-coset radius much below the Tietavainen scale. The constructed distribution is not asserted to be uniform on a linear subspace; the paper explicitly leaves that stronger linear-code realization open. It therefore is a methodological barrier, not a lower bound for cut codes.

#### 8. Louay Bazzi, “On the Covering Radius of Small Codes Versus Dual Distance” (2019)

- Primary sources: [arXiv:1707.06628](https://arxiv.org/abs/1707.06628); [DOI](https://doi.org/10.1109/TIT.2018.2857495); [author/institutional PDF](https://scholarworks.aub.edu.lb/server/api/core/bitstreams/31db7970-04a9-4f73-82ee-e40cda3ef6f1/content).
- Status: **ALMOST-ALL**, upgraded to **ALL-COSET after an existential extension**.
- Precise representative consequence. For odd dual distance `d>=7`, in the regime `d=o(n)` covered by the paper, set

  `R = n/2 - sqrt(((d-5)n/13) log(n/(d-1)))`.

  The fraction of words farther than `R` from the code is at most

  `((d-1)/n)^((d-5)/13)`.

  Thus the first theorem controls almost every coset, not the deepest one. A later corollary shows that there exists a binary linear space `D` with `dim D <= ceil(log_2 n)` such that the enlarged code `C+D` has (full, all-word) covering radius at most the same displayed `R`.
- Proof mechanism. Compare coset-weight distributions with the binomial law in `L_1`, use Markov concentration, create bilateral dual distance by appending independent coordinates, and then use Cohen's iterative completion argument to square the density of the uncovered set.
- Relevance boundary. The all-coset conclusion changes the code by an existential, generally unstructured extension `D`. For a prescribed one-dimensional augmentation `<1>`, the almost-all statement cannot simply be promoted to a deepest-coset statement.

### 2.3 Deterministic length transfer, amalgamation, and normal codes

#### 9. R. L. Graham and N. J. A. Sloane, “On the Covering Radius of Codes” (1985)

- Primary sources: [author-hosted paper](https://mathweb.ucsd.edu/~ronspubs/85_01_covering_radius.pdf); [DOI](https://doi.org/10.1109/TIT.1985.1057039).
- Status: **CONDITIONAL ALL-COSET** length transfer.
- Native definitions. Split a binary code at coordinate `i` into the shortened halves `C_0^(i)` and `C_1^(i)`. Its coordinate norm is the maximum, over ambient words, of the sum of the distances to these two halves; an acceptable coordinate is one realizing a norm at most `2R(C)+1`, and a code with such a coordinate is normal.
- Amalgamated-direct-sum theorem (Theorem 19). If codes `B,C` have acceptable coordinates, glue those coordinates to obtain `B dot+ C`, of length `n_B+n_C-1` and (in the linear independent case) dimension `k_B+k_C-1`. The norm satisfies

  `N(B dot+ C) <= N(B)+N(C)-1`.

  In particular, if both codes are normal,

  `R(B dot+ C) <= R(B)+R(C)`.

  If equality holds in the radius inequality, the amalgam is itself normal. Their Theorem 20 glues copies of the repetition code `{000,111}` and produces normal `[n+2i,k]` codes of radius `R+i` from a normal `[n,k]` code.
- Proof mechanism. In parity-check/syndrome language, split a target syndrome between the two factors. The shared check column permits the two half-distances to be matched with a one-coordinate saving.
- Relevance boundary. This is a genuine every-syndrome transfer, not a volume heuristic. It requires a compatible acceptable coordinate and the exact amalgamated structure; ordinary direct sums or an arbitrary augmentation do not inherit the conclusion automatically.

#### 10. Karen Kilby and N. J. A. Sloane, “On the Covering Radius Problem for Codes I: Bounds on Normalized Covering Radius” (1987)

- Primary source: [DOI](https://doi.org/10.1137/0608049).
- Status: asymptotic/fixed-dimension **CONDITIONAL ALL-COSET** transfer.
- Main stabilization theorem. For fixed dimension `k` and all sufficiently large lengths `n`, an optimal normalized covering-radius code can be chosen normal and with all but one projective column type of multiplicity one. Consequently their extremal normalized parameter satisfies

  `t[n+2,k]=t[n,k]+1`

  for all sufficiently large `n` at fixed `k`.
- Proof mechanism. Record multiplicities of projective parity-check columns, subtract the forced half-multiplicity contribution to define normalized radius, contract column types, and prove stabilization/normality.
- Relevance boundary. The theorem's quantifiers are `k` fixed and `n -> infinity`. For complete-graph cut codes, both the dimension and the set of projective column types grow with graph order, so this recurrence does not transfer without new uniform input.

### 2.4 Completely regular codes and coset graphs

#### 11. J. Borges, J. Rifa, and V. A. Zinoviev, “On Completely Regular Codes” (2017/2019)

- Primary sources: [arXiv:1703.08684](https://arxiv.org/abs/1703.08684); [DOI](https://doi.org/10.1134/S0032946019010010).
- Status: **CONDITIONAL ALL-COSET**.
- Definition. A code `C` with covering radius `rho` is completely regular when its distance partition

  `C(0),C(1),...,C(rho)`

  is equitable: every vertex in layer `ell` has fixed numbers `c_ell` and `b_ell` of neighbors in layers `ell-1` and `ell+1`. Equivalently, the weight/distance distribution of a translate depends only on its distance from `C`.
- Structural theorems. The rank of the outer distribution matrix is `s+1`, where `s` is external distance; one always has `rho<=s`. For a completely regular code, `rho=s`. The converse is false: the review gives the extended binary quadratic-residue `[48,24,12]` code with `rho=s=8` that is not completely regular. For a linear completely regular code, the coset graph is distance-regular of diameter `rho` with the same intersection array. Complete transitivity—`Aut(C)` having exactly `rho+1` orbits on cosets—implies complete regularity.
- Proof mechanism. The Bose-Mesner outer distribution module turns equitable layers into a tridiagonal quotient action; for a linear code, quotienting the Hamming graph by cosets gives the distance-regular coset graph.
- Relevance boundary. Complete regularity would collapse the extreme-coset problem to finitely described layers, but it must be proved. A large automorphism group, a known weight enumerator, or equality `rho=s` alone is insufficient. Existing classification/nonexistence results for completely transitive linear codes at large minimum distance do not automatically apply to a growing family of dual cut/cycle codes.

### 2.5 Rooted LP/SDP, Terwilliger algebras, and growing degree

#### 12. Alexander Schrijver, “New Code Upper Bounds from the Terwilliger Algebra and Semidefinite Programming” (2005)

- Primary sources: [author-hosted paper](https://homepages.cwi.nl/~lex/files/codes.pdf); [DOI](https://doi.org/10.1109/TIT.2005.851748).
- Status: **PACKING ONLY**.
- Finite theorem. For a binary code of length `n` and minimum distance at least `d`, introduce normalized triple-orbit variables `x^t_{i,j}`, where a rooted triple has distances `i,j` from the root and intersection parameter `t`. Two symmetrized moment matrices—one averaging over roots in the code and one over roots outside it—are positive semidefinite. Together with orbit symmetry, nonnegativity, normalization, and forbidden-distance equations, maximizing

  `sum_i binom(n,i) x^0_{i,0}`

  gives an upper bound on `A_2(n,d)` that strengthens Delsarte's LP. The number of orbit variables is `O(n^3)`, and the Terwilliger algebra block diagonalization makes the program polynomial-size.
- Proof mechanism. Stabilize the zero word, average rank-one matrices over the automorphism group of the Hamming cube, and decompose the resulting noncommutative Terwilliger algebra into irreducible blocks.
- Rooting boundary. The root is made a codeword through isometric averaging of a packing configuration. No constraint says that every arbitrary ambient root lies within a specified radius of the code. Therefore these triple variables are richer than a pair distribution but remain unsuitable, without new localizing constraints, for deep holes.

#### 13. Dion Gijswijt and Sven Polak, “Semidefinite Lower Bounds for Covering Codes” (2025; revised 2026)

- Primary source: [arXiv:2504.01932v2](https://arxiv.org/abs/2504.01932).
- Status: **ALL-ROOT COVERING**, but it optimizes code size at a prescribed radius rather than the radius of a prescribed linear code.
- Problem normalization. `K_q(n,r)` is the minimum cardinality of a subset of `{0,...,q-1}^n` whose radius-`r` Hamming balls cover the full space.
- Binary cubic SDP bound (Theorem 4.9). Fix any valid family of all-root inequalities `sum_i lambda_i A_i(u)>=beta` for every radius-`r` covering code. For binary triple-orbit variables `x^t_{i,j}` satisfying the paper's basic symmetry constraints, Terwilliger PSD constraints, Lasserre localizing constraints induced by that inequality, and matrix cuts (Propositions 4.2, 4.3, 4.5, and 4.8),

  `K_2(n,r)^3 >= min_x 2^n sum_{i,j,t}
       multinomial(n; i-t, j-t, t, n-i-j+t) x^t_{i,j}`.                       (D)

  The `q`-ary theorem uses variables `x^{t,p}_{i,j}` and the corresponding orbit multiplicity

  `(q-1)^(i+j-t) (q-2)^(t-p)
   multinomial(n; p,t-p,i-t,j-t,n-i-j+t)`.

  The relevant invariant spaces have `binom(n+3,3)` binary and `binom(n+4,4)` nonbinary orbit dimensions before block decomposition.
- The genuinely covering ingredient. For every ambient root `u`, the local distance counts `A_i(u)` obey inequalities of the form

  `sum_i lambda_i A_i(u) >= beta`.

  These are inserted as scalar and matrix localizers; unlike packing Terwilliger constraints, they retain the universal root.
- Proof mechanism. Average code/complement three-point moment matrices, add Lasserre localizing matrices encoding ball coverage, strengthen with matrix cuts, then block-diagonalize under the stabilizer of zero.
- Relevance boundary. This is the closest current SDP framework to a coset-sensitive program. Published computations give record finite-parameter lower bounds on the *size* of unrestricted covering codes. The theorem does not impose linearity, a prescribed cut-code generator, or an asymptotic theorem for the extreme coset of a growing graph family.

#### 14. Cordian Riener, Jan Rolfes, and Frank Vallentin, “A Semidefinite Programming Hierarchy for Covering Problems in Discrete Geometry” (2023; revised 2025)

- Primary source: [arXiv:2312.11267](https://arxiv.org/abs/2312.11267).
- Status: genuine covering hierarchy; supplies a symmetry **NO-GO** at its first level.
- Theorem. For a compact metric space `X`, their moment/localizing hierarchy lower-bounds the minimum number `N(X,r)` of radius-`r` balls required to cover `X`, and has finite convergence in the finite-space case. If the isometry group acts transitively and every radius-`r` ball has normalized measure `omega_r`, Theorem 5.1 shows that the first primal and dual levels both equal exactly

  `1/omega_r`.
- Proof mechanism. Form measures on finite subsets, impose PSD moment and coverage-localizing conditions, and average feasible measures/functions over the transitive isometry group. At level one, averaging reduces the dual certificate to the constant volume witness.
- Relevance boundary. On a homogeneous Hamming space, first-level symmetry reduction cannot beat the sphere-covering/volume bound. Higher levels or code-specific constraints are necessary. The hierarchy concerns minimum cover size, not directly the radius of a fixed code.

#### 15. OpenAI, “Improved Bounds for Binary and Spherical Codes,” Chapter 2 of *Ten Advances in Mathematics and Theoretical Computer Science* (2026)

- Primary sources: [official research page](https://openai.com/index/ten-advances-in-mathematics/); [full paper, Chapter 2](https://cdn.openai.com/pdf/ten-proofs-oai.pdf); [Lean 4 certificates (`MetricCodes.lean`)](https://github.com/openai/ten-proofs).
- Status: **PACKING ONLY**, despite its especially relevant moving-root representation machinery. Institutional preprint updated 2026-08-06; the official page links a Lean formalization. It is not presented here as a covering theorem.
- Binary normalization and main theorem. Let `A_2(n,d)` be the maximum size of a binary length-`n` code with minimum distance `d`, and

  `R_2(delta)=limsup_{n->infinity} n^{-1} log_2 A_2(n,ceil(delta n))`.

  If `M_2(delta)` is the fully optimized second McEliece-Rodemich-Rumsey-Welch exponent, Theorem 1.1 constructs an explicit variational exponent `kappa_bin` and proves, for every fixed `0<delta<1/2`,

  `R_2(delta) <= kappa_bin(delta) < M_2(delta)`.                               (E)

  Thus the improvement is exponential and strict at every relative distance in the open range.
- Finite moving-projection theorem. Choose integers `0<=k<L<=n-k`, let `J_H(n,k,L)` be the symmetric tridiagonal transition matrix on Fourier levels `k,...,L`, with off-diagonal entry

  `c_{i,H}^{(k)} = ((i-k+1)(n-i-k)/n) / sqrt((i+1)(n-i))`,

  and let `lambda` be its largest eigenvalue. If `s=1-2d/n` and `lambda>s`, Theorem 2.1 states

  `A_2(n,d) <= ((1-s)/(d_k^square (lambda-s)))
                sum_{i=k}^L binom(n,i)`,                                      (F)

  where `d_k^square=binom(n,k)-binom(n,k-1)` is the Boolean harmonic dimension. At `k=0` this is the classical rank-one Krawtchouk path; positive `k` divides the ambient rank by an exponentially large stabilizer representation.
- Spherical theorem. For `A(n,s)` the largest code on `S^{n-1}` with pairwise inner products at most `s`, define `R_sph(s)` by exponential growth rate. Theorem 1.2 builds a hierarchy indexed by interlacing normalized row lengths

  `a_1>b_1>a_2>...>b_r>a_{r+1}>=0`.

  With the explicitly defined transition limit `Gamma_r(a,b)` and dimension exponent `Phi_r(a,b)`,

  `2 Gamma_r(a,b)>s  =>  A(n,s)<=2^((Phi_r(a,b)+o(1))n)`.

  Corollary 1.3 proves strict improvement at every hierarchy step and, for every `0<s<1`, a strict improvement over the spherical-cap-optimized Kabatianskii-Levenshtein exponent.
- Proof mechanism / why “moving projection” is literal. At each binary word `x`, the construction embeds a copy of the nontrivial Boolean harmonic module `E_k` across several Fourier levels and projects onto this `x`-dependent copy. In the spherical setting it uses `E_{k,x}=H_k(x^perp)`, harmonics on the tangent sphere at `x`, embedded across ambient harmonic degrees. Strong-Gelfand multiplicity-free branching makes overlaps scalar two-point kernels even though the subspace moves with the code point. A Perron eigenvector of the block-Jacobi/transition path produces the positive-definite certificate.
- Growing-degree lesson. The harmonic degree `k` and cutoff `L` scale linearly with `n`; fixed-degree association-scheme information would miss the entropy subtraction `H_2(k/n)` (and its spherical analogue). This is a concrete, verified example in which nontrivial stabilizer modules at growing degree beat the classical fixed-line LP.
- Relevance boundary. Every inequality in (E)-(F) assumes pairwise separation among codewords. The moving roots are **code points**, not arbitrary words that must be covered. The formal certificate validates a new packing LP construction; it does not imply an augmented-cut covering deficit.

### 2.6 Recent exact covering and machine-assisted packing results

#### 16. Minjia Shi, Shitao Li, Tor Helleseth, and Ferruh Ozbudak, “Determining the Covering Radius of All Generalized Zetterberg Codes in Odd Characteristic” (2024)

- Primary source: [arXiv:2411.14087](https://arxiv.org/abs/2411.14087).
- Status: **ALL-COSET**, for a special algebraic family.
- Family and theorem. For odd prime-power `q_0` and integer `s>=1`, `C_s(q_0)` is the generalized Zetterberg code of length `q_0^s+1` over `F_{q_0}`. Earlier work had settled all cases except `q_0^s == 7 (mod 8)`. In this remaining congruence class, Theorem 4.2 gives an exact finite-field criterion (`NP_i` for an explicitly indexed family in the paper): the radius is `3` iff one of the `NP_i` conditions holds, and otherwise it is `2`. The authors prove radius `2` for `s=1` and radius `3` beyond an explicit threshold in the relevant extension degrees, thereby resolving all generalized Zetterberg codes. They also define twisted half codes of length `(q_0^s+1)/2` and prove the same radius dichotomy.
- Proof mechanism. Represent every syndrome by at most three parity-check columns. The two-column obstruction and three-column existence become equations over finite fields; the latter is converted to rational points on algebraic curves and controlled with character/Weil estimates. A missing two-column representation supplies the matching lower bound.
- Relevance boundary. This is a model example of an exact proof that really checks every syndrome, rather than an average coset. Its constant radius and cyclic finite-field geometry are highly special; the theorem supplies no generic statement for binary cut/cocycle codes.

#### 17. Henry Cohn, David de Laat, and Nando Leijenhorst, “Optimality of Spherical Codes via Exact Semidefinite Programming Bounds” (2024)

- Primary source: [arXiv:2403.16874](https://arxiv.org/abs/2403.16874).
- Status: **PACKING ONLY**, with an especially strong exact-verification workflow.
- Theorem. The paper proves optimality for a list of spherical codes, uniqueness for those in dimensions at most `64`, and universal optimality of the `288`-point configuration in `R^16`. New optimal cases include `56` points in `R^20`, `50` and `77` points in `R^21`, and Kerdock-derived spherical configurations of sizes `4224`, `66048`, and `1050624` in dimensions `64`, `256`, and `1024`. It obtains corresponding optimality consequences for binary Kerdock codes, including uniqueness in the length-64 case.
- Proof mechanism. Start from the Bachoc-Vallentin three-point SDP, solve numerically, recognize the output in exact rational/algebraic number fields, and verify every PSD and sign condition exactly. Complementary-slackness identities restrict possible inner products and prove uniqueness. Exact matrices and verifier code accompany the paper.
- Relevance boundary. “Exact SDP” here means a rigorous upper bound for pairwise-separated spherical configurations. It is excellent evidence that large symbolic certificates can be extracted from numerics, but it has no arbitrary-root coverage constraint and no deep-hole conclusion.

#### 18. Alexander Novikov et al., “AlphaEvolve: A Coding Agent for Scientific and Algorithmic Discovery” (2025)

- Primary source: [arXiv:2506.13131](https://arxiv.org/abs/2506.13131).
- Authors: Alexander Novikov, Ngan Vu, Marvin Eisenberger, Emilien Dupont, Po-Sen Huang, Adam Zsolt Wagner, Sergey Shirobokov, Borislav Kozlovskii, Francisco J. R. Ruiz, Abbas Mehrabian, M. Pawan Kumar, Abigail See, Swarat Chaudhuri, George Holland, Alex Davies, Sebastian Nowozin, Pushmeet Kohli, and Matej Balog.
- Status: AI-assisted **PACKING CONSTRUCTION ONLY**.
- Checkable coding-theory result. The paper gives an explicit `593`-point kissing configuration in `R^11`. Its conversion lemma is elementary and exact: if a finite set `C subset R^d\{0}` satisfies

  `min_{x != y in C} ||x-y|| >= max_{x in C} ||x||`,

  then the normalized centers `2x/||x||` form a kissing configuration. The supplied `593` integral-coordinate vectors satisfy the required norm and pair-distance inequalities, yielding `K(11)>=593`.
- Proof/verification mechanism. An evolutionary coding agent searches over construction programs and scores exact finite configurations. Once the integral witness is emitted, norms and pairwise distances are independently checkable by integer arithmetic; the geometric lemma proves validity.
- Relevance boundary. This is discovery plus exact witness validation, not an optimality theorem. The objective is minimum pair distance after normalization, not covering radius. It demonstrates a reproducible certificate pattern but supplies no all-coset inequality.

## 3. Cross-source synthesis in the problem's native objects

### 3.1 What data are pair-sensitive, and what data retain an extreme coset?

For a binary linear code `C` and ambient root `u`, the useful hierarchy of data is:

1. **Weight enumerator of `C`**: counts `wt(c)` based at zero. MacWilliams transforms this into the weight enumerator of `C^perp`. It does not retain `u`.
2. **Coset/outer distribution at `u`**: counts `A_i(u)=|{c in C:d(u,c)=i}|`. The first nonzero coordinate is `d(u,C)`. A theorem uniform in `u` is therefore genuinely all-coset.
3. **Rooted pair/triple moments**: record intersections among displacement supports around `u`. These retain more than the distance distribution only when `u` remains an arbitrary ambient root in the constraints.
4. **Stabilizer modules at a moving point**: decompose rooted functions under `Stab(u)`. Their quantifier still matters: moving the module with code points proves a packing bound; moving it with every ambient point and imposing coverage localizers is covering-sensitive.

This yields the following exact audit.

| Source/method | Root that survives symmetrization | Extreme-coset information? |
|---|---|---:|
| MacWilliams | zero only | no |
| Delsarte external distance | arbitrary translate `u+C` | yes, coarse |
| Tietavainen | arbitrary translate, through low dual moments | yes |
| Signed switching | the signing/coset itself | yes, exact native objective |
| Complete regularity | every coset layer | yes, if the hypothesis holds |
| Schrijver SDP | a codeword root in a packing | no |
| Gijswijt-Polak SDP | every ambient root through localizers | yes |
| OpenAI moving projections | each code point in a packing | no |

For the augmented cut code there is an additional bilateral requirement. If `w_u` is the minimum ordinary cut-coset weight, augmentation replaces it by `min(w_u,N-w_u)` after optimizing the same switch. At the signed-matrix level this is exactly the absolute value in (A). A one-sided coset-weight estimate can therefore lose the relevant antipodal identification.

### 3.2 Closest structural theorems to highly symmetric binary cut/cocycle codes

- **Exact equivalence, no loss:** Sole-Zaslavsky. Cosets are switching classes and coset-leader weight is frustration.
- **Complete graph, ordinary code:** Petersdorf as recorded by Sole-Zaslavsky and Gadouleau-Zeng. The ordinary deepest class is uniquely antibalanced, with radius `floor((n-1)^2/4)`.
- **Antipodal quotient:** formula (A). It explains exactly why the preceding extremizer is irrelevant after adding the all-one word.
- **Root descendants:** Taylor. Regular two-graphs have two-eigenvalue Seidel matrices and strongly regular descendants after choosing a root. This is the cleanest theorem connecting switching symmetry to a stabilizer representation.
- **Distance-layer collapse:** Borges-Rifa-Zinoviev. If the relevant code were completely regular, its coset graph would be distance-regular and `rho=s`; symmetry alone does not establish this.
- **Universal generic bounds:** volume and Tietavainen give (B)-(C). Their gap is large because dimension and dual minimum distance discard almost all cycle/cut incidence.

No retrieved structural theorem identifies the deepest class of `C*(K_n)+<1>` or turns regular-two-graph spectra alone into the maximum in (A).

### 3.3 What length-transfer theorems actually transfer

| Operation | Hypotheses | Guaranteed radius statement | Failure mode outside hypotheses |
|---|---|---|---|
| Direct product/direct sum | independent coordinate blocks | radii add in the elementary Cartesian case | does not model a shared cut vertex/edge constraint |
| Graham-Sloane amalgamated sum | acceptable glued coordinates; normal factor codes | `R(B dot+ C)<=R(B)+R(C)` | a non-normal coordinate can have incompatible nearest halves |
| Repetition-code gluing | initial normal code | length `+2i`, radius `+i` | special parity-check column, not arbitrary augmentation |
| Kilby-Sloane stabilization | dimension fixed, length sufficiently large | `t[n+2,k]=t[n,k]+1` | no uniformity when `k` and projective support grow |
| Bazzi completion | code with required dual-distance regime; auxiliary `D` chosen existentially | `dim D<=ceil(log_2 n)` and full covering at the almost-covering radius | `D` is not prescribed and need not preserve graph structure |

These are all-syndrome results when their hypotheses hold. None states that adding the single vector `1`, deleting a graph vertex, or passing from `K_n` to `K_{n+1}` changes covering radius by a fixed amount.

### 3.4 Association schemes when degree grows with length

Three retrieved facts delimit the issue.

1. Delsarte/MacWilliams identities are finite and exact at every degree, but an asymptotic conclusion needs uniform control of Krawtchouk roots or coefficients as the degree changes with length.
2. Bazzi's 2017 duality shows a real barrier for degree `k<=n^(1/3)/log^2 n`: low-degree moment witnesses cannot generally rule out support radius `n/2-sqrt(kn)`. This is not a cut-code lower bound, but it rules out a universal moment-only shortcut in that range.
3. The 2026 moving-projection theorem uses `k=Theta(n)` and `L=Theta(n)`. The nontrivial stabilizer dimension contributes a full entropy term, and this is exactly what makes its packing exponent strictly better than MRRW. Fixed-degree limits would erase the improvement.

Schrijver and Gijswijt-Polak show the other axis: increasing *point order* from pairs to triples. Growing polynomial degree and increasing point order are distinct enrichments. A three-point SDP can still be packing-only, while a two-point construction can use sophisticated growing stabilizer modules. All-coset sensitivity enters through the universal root/localizing constraint, not automatically through either enrichment.

### 3.5 Explicit no-go and non-transfer statements

1. **Petersdorf is not the augmented answer.** Its unique ordinary extremizer `-K_n` becomes a word of `C^pm`.
2. **A MacWilliams transform is not a coset theorem.** It controls the zero-based distribution; external/outer distributions are the extra object needed for arbitrary translates.
3. **Fixed dual distance throws away the expanding graph structure.** For `C^pm`, `d^perp=4` for all `n>=4`; Tietavainen therefore sees only degree-2 Krawtchouk information.
4. **Almost every coset does not bound the deepest coset.** Bazzi 2019 needs an additional existential space to make that transition.
5. **Low-degree moments have a quantified limitation.** Bazzi 2017 constructs limited-independent supports meeting the `n/2-sqrt(kn)` obstruction.
6. **First-level transitive covering SDP is only volume.** Riener-Rolfes-Vallentin prove it equals `1/omega_r` exactly.
7. **Large automorphism groups do not imply complete regularity.** Even `rho=s` is not a converse; the `[48,24,12]` example witnesses the gap.
8. **A Terwilliger root need not be a covering root.** Schrijver's root is a member of the packed code; Gijswijt-Polak's localizers quantify over all ambient roots.
9. **Moving projections do not change the quantifier.** OpenAI's stabilizer module moves with every code point, hence remains a pairwise-separation certificate.
10. **Normal-code recurrences are conditional.** Amalgamation requires acceptable coordinates, and fixed-dimension stabilization is not uniform in growing dimension.
11. **Exact computational certification is objective-specific.** Exact PSD matrices or integer witnesses establish precisely their packing inequalities; they do not transfer to covering without coverage constraints in the certificate.

### 3.6 Recent capabilities and their verification level (2024-2026)

| Work | Mathematical capability | Verification artifact / status | Scope boundary |
|---|---|---|---|
| Shi-Li-Helleseth-Ozbudak 2024 | exact radius 2/3 for every generalized Zetterberg code in odd characteristic | theorem proof via finite-field equations and algebraic-curve point bounds | all syndromes, special family |
| Cohn-de Laat-Leijenhorst 2024 | new exact optimal spherical codes and uniqueness | exact rational/algebraic SDP matrices plus verifier code | packing |
| Gijswijt-Polak 2025/2026 | stronger finite covering-code lower bounds | reproducible symmetry-reduced SDP formulation and reported computations | covering size, unrestricted code |
| AlphaEvolve 2025 | explicit improved 11-dimensional kissing construction | integral coordinate list; norm and distance checks are exact | packing construction, not optimality |
| OpenAI 2026 Chapter 2 | strict asymptotic improvement over optimized MRRW and KL for every interior parameter | full proof plus official Lean formalization link; institutional preprint, not a substitute for independent peer review | packing exponents |

The verification ladder matters. An explicit integer configuration is easy to validate but proves only feasibility. An exact SDP dual proves an upper bound for the encoded optimization problem. A formalized asymptotic theorem checks the stated derivation relative to its formal definitions. None of these artifacts broadens a packing objective into an all-coset covering objective.

## 4. Compact answer sheet

- **Native covering object:** a switching class/coset, not a codeword pair. For the augmented complete-graph cut code it is the absolute switching quadratic form (A).
- **Best foundational all-coset tools retrieved:** signed frustration equivalence; Delsarte external distance; Tietavainen Krawtchouk-root bound; Graham-Sloane normal amalgamation; complete-regular distance partitions.
- **Best current all-root optimization framework:** Gijswijt-Polak's 2025/2026 covering SDP with explicit localizing inequalities for every ambient root.
- **Closest current moving-stabilizer theorem:** OpenAI 2026, with nontrivial Boolean/tangent harmonic modules of degree proportional to length. It is a rigorous packing theorem only.
- **Exact recent all-syndrome exemplar:** generalized Zetterberg codes (Shi et al. 2024), solved through parity-check-column representations and algebraic curves.
- **Extreme-coset warning:** ordinary complete-graph frustration has an exact unique extremizer, but augmentation removes it; almost-covering estimates and ordinary weight enumerators do not identify the replacement.
- **Generic augmented-code window from retrieved theorems:**

  `N/2-O(n^(3/2)) <= rho(C*(K_n)+<1>) <= N/2-Omega(n)`, with `N=binom(n,2)`.

  This window is a normalization check, not a claimed sharp asymptotic.
- **Main literature fault line:** pair/triple richness, semidefinite positivity, and moving representations are not sufficient indicators of covering relevance. The certificate must retain an arbitrary ambient root or an equivalent all-syndrome statement.
