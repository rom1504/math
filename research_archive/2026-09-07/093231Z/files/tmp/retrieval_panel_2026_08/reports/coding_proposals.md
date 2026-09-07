# Coding theory / association-scheme domain report

**Freeze and independence statement.**  The three architectures below were
frozen before any project-archive exposure.  I did not inspect the project
archive, repository history, repository computations, or any steering/ledger
file.  I used the naked packet, the coding-theory toolkit, and the narrowly
relevant primary sources linked below.

## 1. Native translation, with all normalizations checked

Put `N_n=binom(n,2)`.  Map a sign `a_e` to the bit
`u_e=(1-a_e)/2`.  The binary cut code and its prescribed antipodal
augmentation are

```math
C_n^*=\{(s_i+s_j)_{i<j}:s\in\mathbb F_2^n\},\qquad
C_n^+=C_n^*+\langle\mathbf 1\rangle.
```

For `n>=3`, `dim C_n^*=n-1`, `dim C_n^+=n`, and hence
`|C_n^+|=2^n`.  If a switch is represented by `x_i in {+-1}` and
the prescribed antipodal bit by `sigma in {+-1}`, then the number of
disagreements between `A` and the corresponding codeword is

```math
{1\over2}\left(N_n-\sigma\sum_{i<j}a_{ij}x_ix_j\right).
```

Minimizing over `sigma` and `x` gives, with no asymptotic error,

```math
d(u,C_n^+)={N_n-Q(A)\over2},\qquad
Q(A)=N_n-2d(u,C_n^+).
```

Consequently

```math
rho_n:=rho(C_n^+)=\max_u d(u,C_n^+),\qquad
M_n=N_n-2rho_n=\min_u\{N_n-2d(u,C_n^+)\}.
```

Thus a signing `A`, or equivalently its coset `u+C_n^+`, is an
**arbitrary ambient/covering root**.  A deepest coset maximizes `d` and
minimizes the defect

```math
Delta_n(u):=N_n-2d(u,C_n^+)=Q(A).
```

This is not pairwise code data.  The distance/weight distribution among
codewords, MacWilliams data based at zero, and a Terwilliger root chosen to be
a codeword are **packing roots**.  They do not retain the maximization over
arbitrary `u`.  The outer distribution

```math
A_i(u)=|\{c\in C_n^+:d(u,c)=i\}|
```

does retain `u`; its first nonzero coordinate is `d(u,C_n^+)`.

The antipodal augmentation is also essential.  For the ordinary cut code the
objective is the one-sided maximum of `H_A`; for `C_n^+` it is
`max_x|H_A(x)|`.  In particular, the all-negative ordinary deepest class is a
codeword after adding `1`.  Finally, the augmentation here is the **prescribed**
space `\langle 1\rangle`.  It is not the existential auxiliary space in
Bazzi's completion theorem.  In Architecture 1 below we choose the new
coordinates of an ambient root; we never enlarge the prescribed code.

The signed-graph/coset identification is the theorem-level content of
[Sole--Zaslavsky, *A Coding Approach to Signed Graphs*](https://doi.org/10.1137/S0895480189174374).
The extra absolute value and the displayed factors of two are the direct
calculation above.

## 2. Frozen architectures

| # | Mechanism | Intended conclusion | Assessment |
|---|---|---|---|
| 1 | Uniform terminal-band extension of every ambient root by one new vertex | convergence, without identifying the constant | **Best; viable missing lemma, presently unproved** |
| 2 | Degree `Theta(n)=Theta(sqrt(N_n))` all-root moment/Terwilliger certificate | convergence to `1/2` | arbitrary roots destroy the known compression; very low viability |
| 3 | Asymptotic terminal intersection number in the coset graph | convergence to the unique zero of a drift function | exact implication, but uniform lumpability is probably too rigid |

No credible nonconvergence architecture emerged.  Conference orders or other
special subsequences supply constructions, not the strict
`liminf<limsup` and matching universal inequalities required for genuine
nonconvergence.

## 3. Architecture 1 (best): all-root star extension

### 3.1 Exact extension identity

For a hollow signing `B` of order `n` and a star signing
`b in {+-1}^n`, let

```math
E_b(B)=\begin{pmatrix}B&b\\b^T&0\end{pmatrix}.
```

Writing the last spin as `t`,

```math
H_{E_b(B)}(y,t)=H_B(y)+t\,b\mathbin\cdot y.
```

For real `r,s`, `max_{t=+-1}|r+ts|=|r|+|s|`.  Therefore

```math
Q(E_b(B))=\max_{y\in\{+-1\}^n}
             \bigl(|H_B(y)|+|b\mathbin\cdot y|\bigr).                 \tag{1}
```

This is a simultaneous two-channel identity, not a triangle-inequality loss.
In code language, `b` specifies the `n` new star coordinates of one ambient
word on `K_{n+1}`; the target code remains exactly `C_{n+1}^+`.

### 3.2 The exact missing lemma

The threshold `n^{3/2}` below merely isolates the entire terminal band; the
supplied upper frontier ensures that every minimizing root lies in it for all
sufficiently large `n`.

> **Boxed lemma `L_star` (uniform terminal-band star discrepancy).**  There
> exist absolute constants `0<eta<1/2`, `C<infinity`, and `n_0` such that for
> every integer `n>=n_0` and every hollow symmetric
> `B in {+-1}^{n times n}`, there is a vector
> `b=b(B) in {+-1}^n` with the following properties:
>
> ```math
> Q(B)\le n^{3/2}\quad\Longrightarrow\quad
> \max_y\bigl(|H_B(y)|+|b\mathbin\cdot y|\bigr)
> \le \left(1+{3\over2n}\right)Q(B)+C n^{1/2-eta},              \tag{L*}
> ```
>
> and, when `Q(B)>n^{3/2}`, the same selected `b` need only obey the always
> valid bound `Q(E_b(B))<=Q(B)+n`.

Equivalently, for **every ambient coset/root** `u+C_n^+` in the terminal band,
there is an extension `u\oplus z` to the prescribed next code such that

```math
Delta_{n+1}(u\oplus z)
\le\left(1+{3\over2n}\right)Delta_n(u)+C n^{1/2-eta}.            \tag{2}
```

The quantifier is not “there exists a favorable near-minimizer”: it is uniform
over every root in the stated band.  Roots outside the band are also quantified
over, but only the trivial extension assertion is demanded because they can
never be the minimizing root used in the implication.

### 3.3 Imported results plus `L_star` imply convergence

The imported switching/coset theorem above, the exact identity (1), and
`L_star` give the following implication:

```math
\boxed{\text{Sole--Zaslavsky switching/coset equivalence}
+(1)+L_{star}}
\quad\Longrightarrow\quad
M_n/n^{3/2}\text{ converges}.                                  \tag{3}
```

Here is the complete proof of the last arrow.  Choose a deepest root `u_n`, so
`Delta_n(u_n)=M_n`.  From the supplied upper frontier,
`M_n<=n^{3/2}` for all sufficiently large `n`.  Apply (2) and then minimize
over all order-`n+1` roots:

```math
M_{n+1}\le\left(1+{3\over2n}\right)M_n+C n^{1/2-eta}.            \tag{4}
```

Set `a_n=M_n/n^{3/2}`.  Convexity gives

```math
(n+1)^{3/2}\ge n^{3/2}\left(1+{3\over2n}\right),
```

so (4) yields

```math
a_{n+1}\le a_n+C n^{-1-eta}.                                   \tag{5}
```

Let `p_n=max(a_{n+1}-a_n,0)`.  Then `sum_n p_n<infinity`.  Since
`a_n>=0`, finite total upward variation forces convergence: for `m>n`,

```math
a_m-a_n\le\sum_{j=n}^{m-1}p_j,
```

and the right side tends uniformly to zero as `n` tends to infinity.  Taking
`n` along a liminf subsequence proves `limsup a_n<=liminf a_n`.

The coefficient `3/2` is forced by the derivative of `n^{3/2}`.  Merely proving
an `O(sqrt(n))` remainder is not enough: after normalization it gives an
`O(1/n)` error, whose series diverges.  A power saving (or any explicitly
summable normalized error) is the substantive point of `L_star`.

### 3.4 Why `L_star` is strictly weaker than the full histogram or optimization

This is an information projection, not a description-length claim.

1. For a root `u`, `L_star` uses only the endpoint
   `d(u,C_n^+)` (equivalently `Delta_n(u)`) and the endpoint of one chosen
   extension.  It never uses any multiplicity `A_i(u)` above the first nonzero
   coordinate.  This projection is noninjective even on actual cosets of this
   code.  For `n>=6`, let `u_P` be the two-edge path and `u_M` two disjoint
   edges.  Both have coset-leader weight `2`: every nonzero cut has weight at
   least `n-1`, and every complemented cut remains farther after changing two
   coordinates.  But `A_{n-3}(u_P)=1`, contributed by the star at the path's
   degree-two vertex, whereas `A_{n-3}(u_M)=0`.  Indeed a cut producing
   distance `n-3` must have weight `n-1` and meet both chosen edges, hence must
   be exactly that star; all other cuts have weight at least `2(n-2)`, and
   complemented cuts have weight greater than `n-1` for `n>=6`.  Thus equal
   endpoint data can have genuinely different full outer histograms.
2. The only extremal consequence extracted from `L_star` is the one-sided
   recurrence (4).  That recurrence cannot recover `M_n` or even its limiting
   constant.  For every `c>=0`, the rounded sequence
   `m_n=round(c n^{3/2})` satisfies the same form of (4), with a remainder much
   smaller than `n^{1/2-eta}` for any fixed `eta<1/2`, yet these sequences have
   different limits.  Hence the lemma supplies smooth extendability, not the
   original minimizers, their values, or the answer's constant.
3. A witness `b(B)` is one star among `2^n`; it does not select a Boolean
   maximizer `y` in (1).  In particular, a proof may be a discrepancy or
   separation certificate for all `y` simultaneously without returning the
   maximizing `y` or the outer histogram.

### 3.5 Decisive falsifier and finite normalization check

Define the exact extension excess

```math
D_n(B):=\min_{b\in\{+-1\}^n}\max_y
          (|H_B(y)|+|b\mathbin\cdot y|)
        -\left(1+{3\over2n}\right)Q(B).                         \tag{6}
```

A structural family `B_n` with `Q(B_n)<=n^{3/2}` and

```math
D_n(B_n)\ge c\sqrt n
```

for one fixed `c>0` and infinitely many `n` decisively falsifies `L_star` for
every `eta>0`.  This is the cleanest exact target for SAT/ILP or a constructive
counterexample.  For proposed numerical values of `(C,eta)`, a single finite
`B` violating (L*) is already decisive.

The supplied exact anchors do not falsify the necessary extremal recurrence.
They do show that a proof cannot silently demand a nonpositive remainder: at
`n=6`,

```math
M_7-\left(1+{3\over12}\right)M_6=9-{5\over4}\,5={11\over4},
```

and at `n=10` the analogous excess is `41/20`.  These are normalization checks
only; `L_star` is uniform over all roots, which the anchor list does not test.

### 3.6 Nearest imported theorem and every unmatched hypothesis

The nearest deterministic covering transfer is Graham--Sloane's
[amalgamated-direct-sum theorem](https://neilsloane.com/doc/Me114.pdf): normal
codes with acceptable coordinates can be glued with
`R(B dot+ C)<=R(B)+R(C)`.  It does not prove `L_star`.  The unmatched points
are:

- its operation changes/glues two codes, whereas (2) extends one arbitrary
  ambient root inside the fixed family `C_n^+ -> C_{n+1}^+`;
- it assumes acceptable coordinates/normality, neither established here;
- `K_n -> K_{n+1}` adds an entire `n`-edge star, not one amalgamated coordinate;
- it supplies an additive radius bound, not the sharp defect coefficient
  `3/(2n)` with a power-saving `o(sqrt n)` error;
- it has no choice of a star uniformly for every old coset.

Standard vector-discrepancy theorems also do not immediately apply: the
constraint set in (1) has `2^{n-1}` antipodal spin vectors, and the allowed
correlation depends on the individual slack `Q(B)-|H_B(y)|`.  Treating the two
terms separately at their worst values introduces a fixed leading loss.

### 3.7 Circularity audit

A valid proof of `L_star` may use `Q(B)` as the endpoint defining the slack, but
it must not do any of the following:

- choose `B` by assuming a classification of minimizers or the convergence of
  `M_n/n^{3/2}`;
- choose `b` by comparing with the already minimized value `M_{n+1}`;
- enumerate the full outer distribution of every coset and call that a
  compression;
- replace `max_y(|H_B(y)|+|b.y|)` by two separate maxima and absorb a fixed
  `Theta(sqrt n)` coefficient loss;
- prove the claim only for almost every root, a random root, or a packing root;
- invoke an existential code enlargement.  The only existential object allowed
  in (L*) is the new **ambient star** `b`; the code is prescribed.

### 3.8 Confidence

- Truth of `L_star`: **0.30**.  The terminal-band restriction makes it
  plausible, but uniform control of all near-ground-state spin vectors with a
  power saving below `sqrt n` is strong.
- Tractability of `L_star`: **0.18**.  A weighted discrepancy/minimax proof is
  conceivable; ordinary union bounds and off-the-shelf discrepancy constants
  do not reach it.
- Correctness of implication (3), conditional on `L_star`: **0.99**.

## 4. Architecture 2: growing-degree all-root moment/Terwilliger certificate

Let `X` be uniform on `{+-1}^n` and, for an integer `k>=1`, set

```math
S_{n,k}(A)=2^{-n}\sum_x |H_A(x)|^{2k}.
```

The elementary moment-to-supremum inequality gives
`Q(A)>=S_{n,k}(A)^{1/(2k)}` for every arbitrary signing `A`.  Hence the
following lemma would combine with the supplied `limsup<=1/2` to prove
convergence to `1/2`.

> **Candidate `L_mom`.**  For every `epsilon>0` there are
> `alpha=alpha(epsilon)>0` and `n_0` such that, with
> `k_n=ceil(alpha n)`, for every `n>=n_0` and every hollow signing `A`,
>
> ```math
> S_{n,k_n}(A)^{1/(2k_n)}
>       \ge (1/2-epsilon)n^{3/2}.                               \tag{7}
> ```
>
> Moreover, (7) has a sign-uniform positive-semidefinite proof in the rooted
> Boolean harmonic algebra at degrees at most `2k_n` with: (i) exact closure
> before taking limits; (ii) at most `exp(o(n))` total block state; and
> (iii) error `o(n^{3/2})` after taking the `2k_n`-th root, uniformly in `A`.

The closure clause is part of the lemma, not an implementation preference.
Expanding the moment gives the signed Eulerian-multigraph sum

```math
S_{n,k}(A)=
\sum_{(e_1,...,e_{2k}):\,\deg_v(e_1+\cdots+e_{2k})\text{ even }\forall v}
     \prod_{j=1}^{2k}a_{e_j}.                                  \tag{8}
```

Thus an alleged compression must be closed under the contractions in (8) and
must control their cancellations; merely listing a small set of representations
does not suffice.  The desired radius displacement is `Theta(n^{3/2})` out of
length `N_n=Theta(n^2)`, i.e. relative scale `n^{-1/2}`.  The natural
Krawtchouk/Terwilliger degree is therefore `k=Theta(n)=Theta(sqrt N_n)`, since
the generic root scale is `sqrt(kN_n)`.  Fixed degree cannot see this scale.
The explicit `exp(o(n))` and `o(n^{3/2})` clauses prevent an exponential state
or a fixed leading loss from being hidden in “growing degree.”

The implication is exact:

```math
L_{mom}+\|f\|_infinity\ge\|f\|_{2k}
+\limsup M_n/n^{3/2}\le1/2
\Longrightarrow M_n/n^{3/2}\to1/2.
```

It is strictly weaker than a histogram: it uses only a lower bound on one
power sum.  The moment map is explicitly noninjective.  For any `p=2k`, the
two equal-length nonnegative lists `(1,0,0,...)` and
`(2^{-1/p},2^{-1/p},0,...)` have the same `p`-th power sum and different
maxima, much less different full histograms.  Nor does (7) identify any
minimizing signing.

**Nearest imported theorems and gap.**  Schrijver's
[Terwilliger SDP](https://homepages.cwi.nl/~lex/files/codes.pdf) has polynomial
blocks but roots at codewords and proves packing bounds.  Gijswijt--Polak's
[covering SDP](https://arxiv.org/abs/2504.01932) retains arbitrary ambient
roots through localizers, but optimizes the size of an unrestricted covering
code at prescribed radius, not the radius of this fixed cut code.  Bazzi's
[limited-independence obstruction](https://arxiv.org/abs/1707.00552) shows why
low-degree univariate moments cannot be treated as a universal shortcut.  The
unmatched hypotheses are: a prescribed linear cut code; arbitrary rather than
codeword roots; degree `Theta(n)`; uniform signed-cycle cancellation; exact
algebraic closure with subexponential state; and `o(n^{3/2})` error.

The classical all-root theorems also stop before the needed scale.
[Tietavainen's dual-distance bound](https://doi.org/10.1109/18.59949), using
`d((C_n^+)^perp)=4`, gives only
`rho_n<=(N_n-sqrt(N_n))/2`, hence `M_n>=sqrt(N_n)=Theta(n)`.
Delsarte external distance is universal over cosets but discards the cycle
incidence once reduced to the set of dual weights.  Bazzi's obstruction is not
a cut-code counterexample (its witness need not be linear), but in ambient
length `N_n` it rigorously blocks a generic univariate-moment shortcut through
degrees up to `N_n^{1/3}/log^2 N_n`; the desired `Theta(n)=Theta(sqrt N_n)`
degree lies beyond that range and still needs the new closure asserted above.

There is a further structural warning.  The coordinate automorphism group of
the code is essentially `S_n`, but a generic signing has trivial stabilizer in
`S_n`.  If one instead uses the full Hamming-coordinate group `S_{N_n}`, the
cut/cycle incidence is lost.  Thus the familiar polynomial Terwilliger orbit
count does not survive merely by declaring the root to move.

**Decisive falsifier.**  Let `P_j` be Boolean Fourier-level projections and
`D_A=diag(H_A(x))` on the antipodal Boolean space.  An infinite family for
which the commutant of `{P_j,D_A}` is scalar on exponentially large invariant
sectors, together with an `Omega(n)` log-moment separation between signings
identified by any proposed `exp(o(n))` state, falsifies the closure clause.
At finite `n`, the commutant calculation is exact; Burnside's theorem then
certifies full matrix blocks.  Independently, one signing violating (7) for a
proposed `(epsilon,alpha,n_0)` falsifies its numerical certificate.

**Circularity audit.**  A block containing `Q(A)`, all `2^{n-1}` values
`H_A(x)`, or all graph types in (8) is the original optimization/histogram in
disguise.  Averaging over codewords changes a covering root into a packing
root.  Dropping negative terms in (8) is invalid.  An `exp(Theta(n))` block
count or an `exp(Theta(n))` multiplicative uncertainty in the moment produces
a fixed loss after the `Theta(n)`-th root and fails the requested scale.

Confidence: truth of the numerical moment statement **0.20**; truth of the
required closure/compression **0.03**; tractability **0.03**; conditional
implication **0.99**.

## 5. Architecture 3: terminal coset-graph drift

Let `Gamma_n` be the Cayley coset graph on
`F_2^{E_n}/C_n^+`, generated by adding one edge coordinate.  For a coset `U`,
write

```math
r_n(U)=d(U,C_n^+),\qquad
z_n(U)={N_n-2r_n(U)\over n^{3/2}},
```

and define its outward intersection number

```math
b_n(U)=|\{e\in E_n:r_n(U+e)=r_n(U)+1\}|.                         \tag{9}
```

An exact completely regular code would make `b_n(U)` depend only on `r_n(U)`.
The following much weaker terminal statement alone would settle convergence.

> **Candidate `L_drift`.**  There are a compact interval `I` containing
> `[0.33,0.51]`, a continuous function `beta:I->[0,1]`, a point `c in I`, and
> a sequence `epsilon_n->0` such that:
>
> 1. `beta(z)=0` if and only if `z=c`;
> 2. for every sufficiently large `n` and **every coset** `U` with
>    `z_n(U) in I`,
>
> ```math
> \left|{b_n(U)\over N_n}-\beta(z_n(U))\right|\le\epsilon_n.      \tag{10}
> ```

Choose a deepest coset `U_n`.  It has `b_n(U_n)=0`, since no neighbor can lie
farther than the covering radius, and
`z_n(U_n)=M_n/n^{3/2}` lies in `I` by the supplied frontier.  Equation (10)
gives `beta(z_n(U_n))->0`.  Compactness and the unique-zero condition give
`z_n(U_n)->c`.  This proves convergence with no reconstruction of a full
intersection array.

The lemma retains only the first support point of a root and a one-bit
classification of its `N_n` neighboring minima.  It is invariant under every
change to all nonminimal coefficients of the `N_n+1` corresponding outer
distributions.  Hence it cannot recover any of those histograms.  It also does
not identify the deepest root or exact finite radii.  This is a genuine
terminal quotient, strictly less than complete regularity, which requires all
three intersection numbers at every layer and makes every translate's full
distance distribution depend only on its layer.

The nearest theorem is the completely-regular-code result summarized in
[Borges--Rifa--Zinoviev](https://arxiv.org/abs/1703.08684): a linear completely
regular code has a distance-regular coset graph, and complete transitivity is a
sufficient condition.  None of its structural hypotheses is known here.
Specifically missing are exact or approximate equitability, uniformity over
arbitrary terminal roots, convergence to one function `beta`, and uniqueness
of its zero.

Large symmetry does not fill this gap.  For `n>=6`, the minimum nonzero
codewords of `C_n^+` are exactly the `n` vertex-star cuts of weight `n-1`.
A coordinate automorphism must permute those stars; each edge is the unique
intersection of its two endpoint stars.  Thus the effective coordinate
automorphism group on cosets is `S_n`.  There are `2^{N_n-n}` cosets, so the
number of its coset orbits is at least `2^{N_n-n}/n!`, which exceeds
`N_n+1>=rho_n+1` for all sufficiently large `n`.  Consequently the family is
not completely transitive in large order.  Complete regularity is weaker, so
this is a no-go for the symmetry shortcut, not a disproof of (10).

**Decisive falsifier.**  Either of the following kills `L_drift`:

- two infinite root families `U_n,V_n` with
  `z_n(U_n)-z_n(V_n)->0` but
  `|b_n(U_n)-b_n(V_n)|/N_n` bounded away from zero;
- an infinite family of nondeep dead ends, `b_n(U_n)=0`, whose `z_n(U_n)`
  stays separated from the deepest-coset value.

For a proposed error sequence, the exact finite scatter plot of pairs
`(z_n(U),b_n(U)/N_n)` is a decisive test.  Packing distributions cannot perform
this test because (9) asks for the distances of arbitrary neighboring cosets.

**Circularity audit.**  The proof may not define `beta` by first inserting the
unknown deepest values, assume all local maxima of distance are global, infer
complete regularity from `S_n` symmetry, or calculate (10) by tabulating every
outer distribution.  A pairwise or codeword-root intersection number is not
`b_n(U)`.

Confidence: truth **0.08**; tractability **0.06**; conditional implication
**0.99**.

## 6. Specialist verdict and candidate card

The only proposal I recommend retaining is `L_star`.  It has the exact
arbitrary-root quantifier needed for a covering problem, treats the absolute
two-channel objective jointly, extends the prescribed family rather than the
code, and asks for a power-saving local error rather than the limit itself.
The growing-degree route currently has no credible algebraic compression once
the root is an arbitrary signing, and the terminal-drift route likely asks for
too much uniformity.

| Field | Entry for the preferred proposal |
|---|---|
| Domain | Binary covering codes, signed complete graphs, coset graphs |
| Imported theorem(s) | Sole--Zaslavsky switching/coset equivalence; nearest transfer comparison is Graham--Sloane normal amalgamation |
| Problem translation | `Q(A)=N_n-2d(u,C_n^+)`, `M_n=N_n-2rho(C_n^+)` |
| Proposed mechanism | Extend every terminal-band ambient root by a chosen new star and obtain summably small normalized upward variation |
| Exact missing lemma | `L_star`, equation (L*) / (2) |
| Why strictly weaker | Uses two coset-leader endpoints; endpoint projection is noninjective on histograms, and recurrence admits every limiting constant |
| Archive collisions | Not assessed: archive exposure was forbidden for this domain report |
| Falsification test | A family with `D_n(B_n)>=c sqrt(n)` in (6) |
| Specialist confidence | Truth `0.30`; tractability `0.18`; implication `0.99` |
| Recommendation | Retain as a sharply isolated lemma; reject importing packing SDP as if it proved it |
