# Contrarian nonconvergence report — initial, archive-blind

I read only `NAKED_PROBLEM.md` and `REPORT_SCHEMA.md`.  I did not inspect the
archive, ledger, steering files, saved computations, Git history, or any other
repository material.  The three architectures below were frozen before any
archive exposure.

## Executive judgment

I did **not** find a currently plausible proof of nonconvergence.  I found one
arithmetically explicit separation mechanism worth testing, one coding-theory
certificate version of it, and one genuinely different tensor mechanism.  In
all three, the decisive unsolved step is not the production of an unusual
matrix.  It is a lower bound holding for **every** signing on the proposed high
subsequence.  Existing conference, finite-field, excess, bent-function, and
spin-glass theorems do not perform that quantifier reversal.

The strongest frozen mechanism is Candidate 1: the towers

\[
 n_r^- = 5^{2r+1}+1,\qquad n_r^+=5^{2r}+1.
\]

On the even-degree field tower there is an exact trace-fibre Ising state in the
Paley conference matrix whose normalized energy tends to

\[
 \frac12\left(1-\frac1{5^2}\right)=\frac{12}{25}=0.48.
\]

Odd extension degrees do not have this state because the quadratic character
restricts nontrivially to \(\mathbb F_5^*\).  A fixed gap between a glassy odd
tower and a crystalline even tower is therefore a real, exact arithmetic
phenomenon for the **selected Paley models**.  What is wholly unproved, and in
my view unlikely, is that the even-tower obstruction applies to every signing.

My confidence that the full Candidate 1 lemma is true is **0.2%**; confidence
that it is tractable with presently identifiable methods is **0.5%**.  The
conditional implication to strict nonconvergence is elementary and has
confidence **99.9%**.

## 1. Native translation and normalization

### 1.1 Seidel switching

Let \(D_x=\operatorname{diag}(x_1,\ldots,x_n)\).  The switched Seidel matrix is
\(D_xAD_x\), and its (ordered) excess is

\[
 \mathbf 1^TD_xAD_x\mathbf 1=x^TAx=2H_A(x).
\]

Consequently

\[
 Q(A)=\frac12\max_{x\in\{\pm1\}^n}|x^TAx|.
\]

The absolute value is essential.  Standard signed-graph frustration is usually
one-sided, and standard matrix excess often permits independent row and column
signs, i.e. \(y^TAz\).  Neither is automatically the diagonal, antipodal
quantity above.

### 1.2 Exact code object

Put \(E=\binom n2\), identify an edge signing with a word
\(a\in\mathbb F_2^E\), and let \(B_n\) be the binary cut code of \(K_n\).  Then

\[
 C_n^+=B_n+\langle\mathbf1_E\rangle
\]

is a binary linear \([E,n]\) code for \(n\ge3\).  For the sign word
\(c=(\sigma x_ix_j)_{i<j}\),

\[
 \langle a,c\rangle_{\{\pm1\}}=E-2d(a,c)=\sigma H_A(x).
\]

Thus, exactly,

\[
 Q(A)=E-2d(a,C_n^+),\qquad
 M_n=E-2\rho(C_n^+).
\]

The dual has a useful concrete description:

\[
 (C_n^+)^\perp
 =\{F\subseteq E(K_n):\deg_F(v)\equiv0\pmod2\ \forall v,
 \ |F|\equiv0\pmod2\}.
\]

So it is the even-weight part of the Eulerian/cycle space.  Solé and Zaslavsky
give the signed-graph/cocycle-code correspondence, but their usual frustration
parameter is one-sided; the added all-one word above is what pays for the
absolute value here: [Solé–Zaslavsky, *A Coding Approach to Signed Graphs*](https://doi.org/10.1137/S0895480189174374).

### 1.3 A necessary density filter

There is a simple monotonicity that rules out most residue-class stories.
If \(S\subseteq[n]\), fix signs on \(S\), extend them independently and
uniformly to \([n]\setminus S\), and average.  All terms meeting the complement
average to zero, so

\[
 H_{A[S]}(x_S)=\mathbb E[H_A(x_S,X_{[n]\setminus S})].
\]

It follows that \(Q(A[S])\le Q(A)\), hence

\[
 \boxed{m\le n\implies M_m\le M_n.}
\]

In particular, if \(N_k/N\to1\) from above and
\(M_{N_k}\le(c+o(1))N_k^{3/2}\), then the same upper constant holds at order
\(N\).  Paley conference orders coming from primes \(q\equiv1\pmod4\) are
asymptotically dense by the prime number theorem in progressions.  Therefore:

* a fixed congruence class cannot by itself support nonconvergence;
* existence versus nonexistence of ordinary conference or Hadamard matrices at
  nearby orders cannot create a fixed leading gap;
* any viable low subsequence must be multiplicatively sparse, or its special
  estimate must fail to persist under principal restriction.

This filter is why all three frozen proposals use genuine towers rather than
orders modulo \(4\) or \(8\).

## 2. Candidate 1 — fixed-characteristic Paley crystal/glass alternation

### Candidate card

| Field | Entry |
|---|---|
| Domain | Finite fields, quadratic characters, conference matrices, orthogonal mean-field Ising models |
| Imported theorems | Paley conference construction; quadratic Gauss-sum diagonalization; Hasse–Davenport lifting; random-orthogonal-model ground-state heuristics |
| Exact subsequences | low: \(n_r^-=5^{2r+1}+1\); high: \(n_r^+=5^{2r}+1\) |
| Target constants | \(c_-=239/500=0.478\), \(c_+=479/1000=0.479\), gap \(10^{-3}\) |
| Proposed mechanism | Odd and even extension degrees have different restrictions of the quadratic character to the prime field.  Even degrees contain a macroscopic Boolean trace-fibre state in one spectral half; odd degrees do not. |
| Main falsifier | Any infinite even-tower family of arbitrary sign matrices with \(Q/n^{3/2}<0.479\), or odd Paley energies with limsup above \(0.478\) |
| Confidence | truth 0.2%; tractability 0.5%; implication 99.9% |

### 2.1 Imported leverage and exact arithmetic calculation

Paley's theorem constructs, for each odd prime power \(q\equiv1\pmod4\), a
symmetric conference matrix of order \(q+1\):
[Paley, *On Orthogonal Matrices*](https://doi.org/10.1002/sapm1933121311).
Write its core, indexed by \(\mathbb F_q\), as

\[
 B_q(a,b)=\chi_q(a-b),\qquad B_q\mathbf1=0,
 \qquad B_q^2=qI-J.
\]

Additive characters diagonalize \(B_q\); every nonconstant eigenvalue is a
quadratic Gauss sum of magnitude \(\sqrt q\).  The Hasse–Davenport lifting
formula controls the sign under field extension; a modern theorem statement is
given in [Zheng, *Davenport–Hasse's Theorem for Polynomial Gauss Sums over
Finite Fields*](https://doi.org/10.1016/j.jnt.2017.04.005), Theorem 1.1 (whose
classical specialization is all that is used here).

Take \(q=5^m\), a nonzero \(t\in\mathbb F_q\), and a function
\(\phi:\mathbb F_5\to\{\pm1\}\) with \(\sum_{u\in\mathbb F_5}\phi(u)=1\).
Set

\[
 x_a=\phi(\operatorname{Tr}_{\mathbb F_q/\mathbb F_5}(ta)).
\]

If \(m\) is even, then for every \(k\in\mathbb F_5^*\),

\[
 \chi_q(k)=\chi_5(k)^m=1.
\]

All four nonconstant Fourier modes of \(\phi\circ\operatorname{Tr}\) therefore
lie in the same \(B_q\)-eigenspace.  Since \(\sum_a x_a=q/5\), Parseval gives
the exact identity

\[
 |x^TB_qx|
 =\sqrt q\left(q-\frac1q\left(\frac q5\right)^2\right)
 =\frac{24}{25}q^{3/2}.
\]

Adding the conference border changes the quadratic value by only
\(2q/5\), with a freely chosen sign.  Hence the selected even-degree Paley
matrix satisfies

\[
 Q(C_q)\ge \frac{12}{25}q^{3/2}+O(q).
\]

For odd \(m\), in contrast,
\(\chi_q|_{\mathbb F_5^*}=\chi_5\); the four trace-line modes split between the
two spectral halves.  The exact Boolean eigenstate disappears.

This is consistent with the distinction between crystalline and glassy states
in orthogonal Ising models.  Parisi and Potters study random and deterministic
orthogonal interaction matrices in
[Parisi–Potters, *Mean-Field Equations for Spin Models with Orthogonal
Interaction Matrices*](https://arxiv.org/abs/cond-mat/9503009).  The physics is
motivation only, not a rigorous extreme-value theorem.

### 2.2 Exact missing lemma \(L_{\mathrm{FF}}\)

> **\(L_{\mathrm{FF}}\) (fixed-characteristic crystal/glass separation).**
> There is \(r_0\) such that for every integer \(r\ge r_0\):
>
> 1. for the symmetric Paley conference matrix \(C_{5^{2r+1}}\) of order
>    \(n_r^-=5^{2r+1}+1\),
>    \[
>    Q(C_{5^{2r+1}})\le \frac{239}{500}(n_r^-)^{3/2};
>    \]
> 2. for **every** symmetric hollow off-diagonal sign matrix \(A\) of order
>    \(n_r^+=5^{2r}+1\),
>    \[
>    Q(A)\ge \frac{479}{1000}(n_r^+)^{3/2}.
>    \]
>
> Both estimates are uniform in \(r\), and no exceptional signing is allowed
> in part 2.

Named imported results plus this lemma give

\[
 \liminf_n\frac{M_n}{n^{3/2}}\le0.478
 <0.479\le
 \limsup_n\frac{M_n}{n^{3/2}}.
\]

The imported finite-field theorems explain why \(12/25\) is a natural barrier,
but they prove neither numbered assertion of \(L_{\mathrm{FF}}\).

### 2.3 Why this would affect the optimum

Part 1 is an explicit admissible signing, so it is a legitimate upper bound on
\(M_{n_r^-}\).  Part 2 quantifies over the entire feasible set, so it is a
legitimate lower bound on \(M_{n_r^+}\).  No claim about a selected even Paley
matrix is substituted for part 2.

### 2.4 Why the lemma is strictly weaker than full minimization

The lemma supplies two scalar inequalities on two sparse sequences.  It does
not identify any minimizer, determine \(M_n\) exactly even on those sequences,
determine any other order, or give the weight enumerator of a single extremal
coset.  For example, every assignment of values
\(M_{n_r^-}\in[0.3365,0.478](n_r^-)^{3/2}\) and
\(M_{n_r^+}\in[0.479,0.5](n_r^+)^{3/2}\), together with arbitrary admissible
values elsewhere, is logically consistent with the lemma.  Those assignments
have mutually different minimizers and coset histograms.  Thus the lemma cannot
recover full minimization or a full histogram.

### 2.5 Decisive falsifiers

1. **Upper-side falsifier:** show
   \[
   \limsup_{r\to\infty}
   Q(C_{5^{2r+1}})/(5^{2r+1}+1)^{3/2}>0.478.
   \]
2. **Universal-side falsifier:** construct sign matrices \(A_r\) of orders
   \(5^{2r}+1\) with
   \(Q(A_r)/(5^{2r}+1)^{3/2}<0.479\) for infinitely many \(r\).
3. **Structural falsifier:** prove a restriction/rounding theorem that transfers
   the odd-tower Paley bound to the even tower with \(o(n^{3/2})\) loss.

The second is the fastest meaningful test.  A large value for the selected even
Paley matrix is not evidence for the universal statement.

### 2.6 Hypotheses not verified and circularity audit

Unverified hypotheses:

* a strict, uniform odd-tower Boolean ground-state gap below \(0.48\);
* a universality or rigidity theorem saying arbitrary even-tower signings retain
  a finite-field obstruction despite having no field structure;
* stability under the conference border and the absolute two-channel maximum at
  the advertised \(10^{-3}\) accuracy;
* exclusion of nonconference signings with much smaller \(Q\).

Quantifier reversals and channel risks:

* “the even Paley matrix has a high-energy state” is \(\exists A\,\exists x\),
  whereas the required lower theorem is \(\forall A\,\exists x\);
* the spectral norm bound concerns real unit vectors, not Boolean vectors;
* maximum conference excess often permits independent row and column operations
  (a bilinear \(y^TCz\)), whereas this problem requires the same switching on
  both sides (\(x^TCx\)); see
  [Momihara–Suda, *Conference Matrices with Maximum Excess and
  Two-Intersection Sets*](https://arxiv.org/abs/1611.01305);
* one-sided maximum cut or frustration does not pay for the antipodal absolute
  channel;
* a classification saying every *near-conference* matrix is Paley-like would
  already contain a substantial inverse theorem and must not silently assume
  near-optimality forces spectral flatness;
* relabelling vertices by \(\mathbb F_q\) gives no algebraic structure to an
  arbitrary signing.

The likely fatal point is the second clause of \(L_{\mathrm{FF}}\): no imported
theorem supplies even a weak route from extension-field parity to all
\(2^{\binom n2}\) signings.

## 3. Candidate 2 — Arf/Krawtchouk sector cancellation in the antipodal cut code

This is mathematically different from Candidate 1 on the lower-bound side.  It
seeks a universal covering certificate in the Hamming association scheme rather
than a rigidity classification of real matrices.

### Candidate card

| Field | Entry |
|---|---|
| Domain | Binary codes, Krawtchouk linear programming, quadratic refinements over \(\mathbb F_2\) |
| Imported theorems | MacWilliams transform; Delsarte covering-radius method; Arf classification/Gauss sums |
| Exact subsequences | low: \(5^{2r+1}+1\equiv6\pmod8\); high: \(5^{2r}+1\equiv2\pmod8\) |
| Target constants | again \(0.478<0.479\) |
| Proposed mechanism | A root-of-unity sector of the coset enumerator changes Arf sign between orders \(6\) and \(2\bmod8\); a uniform Krawtchouk witness would turn the high sector into an all-cosets central-hole theorem. |
| Main falsifier | The best one-dimensional Krawtchouk certificate has the same asymptotic threshold in the two sectors, or arbitrary coset phases erase the Arf sign. |
| Confidence | truth 0.05%; tractability 1%; implication 99.9% |

### 3.1 Imported leverage

MacWilliams' theorem converts a linear code weight enumerator to its dual
enumerator: [MacWilliams, *A Theorem on the Distribution of Weights in a
Systematic Code*](https://doi.org/10.1002/j.1538-7305.1963.tb04003.x) and the
nonlinear extension
[MacWilliams–Sloane–Goethals](https://doi.org/10.1002/j.1538-7305.1972.tb01947.x).
Delsarte's external-distance theorem and related covering arguments turn
Krawtchouk positivity into covering-radius bounds; a convenient primary
covering formulation is
[Navon–Samorodnitsky, *Linear Programming Bounds for Codes via a Covering
Argument*](https://arxiv.org/abs/math/0702425).

There is a real mod-eight sector to exploit.  For even \(n\),

\[
 q_n(z)=\sum_{i<j}z_iz_j\quad(z\in\mathbb F_2^n)
\]

is nondegenerate, with polar form \(I+J\).  Its Gauss sum is

\[
 \sum_z(-1)^{q_n(z)}
 =2^{n/2}\left(\cos\frac{n\pi}{4}+\sin\frac{n\pi}{4}\right).
\]

It is positive at \(n\equiv2\pmod8\) and negative at
\(n\equiv6\pmod8\).  These are precisely the two residues of
\(5^m+1\) as \(m\) alternates parity.  Bent/Kerdock theory shows that Arf and
rank sectors can produce exact square-root gaps in related Reed–Muller
covering problems; the relevant construction language appears in
[Calderbank et al., *The \(\mathbb Z_4\)-linearity of Kerdock, Preparata,
Goethals, and related codes*](https://doi.org/10.1109/18.312154).

For the present code, however, an edge signing adds a linear term to the
quadratic refinement at a root of unity and can flip its Arf sign.  The magnitude
is stable; the useful sign is not automatically uniform over cosets.

### 3.2 Exact certificate form

Let \(D_n=(C_n^+)^\perp\), and let \(K_j^{(E)}(w)\) be the binary
Krawtchouk polynomial.  For a radial function

\[
 P_n(w)=\sum_{j=0}^E\widehat P_{n,j}K_j^{(E)}(w),
\]

Fourier inversion on \(\mathbb F_2^E\) gives, for every coset representative
\(a\),

\[
 \sum_{c\in C_n^+}P_n(\operatorname{wt}(a+c))
 =|C_n^+|\sum_{u\in D_n}
   \widehat P_{n,\operatorname{wt}(u)}(-1)^{a\cdot u}.
\]

Thus the following domination condition is genuinely uniform in the coset:

\[
 \widehat P_{n,0}>
 \sum_{u\in D_n\setminus\{0\}}
 \left|\widehat P_{n,\operatorname{wt}(u)}\right|.
\tag{KP}
\]

If additionally \(P_n(w)\le0\) throughout a proposed central interval, the
positive sum forces every coset to have a word outside that interval.

### 3.3 Exact missing lemma \(L_{\mathrm{KP}}\)

> **\(L_{\mathrm{KP}}\) (sector-separated covering certificate).**  There is
> \(r_0\) such that for every \(r\ge r_0\):
>
> 1. at \(n=5^{2r+1}+1\), there is an explicitly specified edge word
>    \(a_r^-\) (for example the Paley word) such that every
>    \(c\in C_n^+\) satisfies
>    \[
>    |E-2\operatorname{wt}(a_r^-+c)|
>    \le0.478\,n^{3/2};
>    \]
> 2. at \(n=5^{2r}+1\), there is a radial function \(P_n\), with rational
>    Krawtchouk coefficients, satisfying (KP) and
>    \[
>    P_n(w)\le0\quad\text{whenever}\quad
>    |E-2w|<0.479\,n^{3/2}.
>    \]
>
> The coefficient domination and pointwise sign condition hold exactly, not
> merely after averaging over cosets.

Part 1 yields \(M_n\le0.478n^{3/2}\).  For part 2, if a coset were entirely in
the central interval, the left side of the Fourier identity would be
nonpositive, while (KP) makes it positive for every phase vector
\(((-1)^{a\cdot u})_u\).  Therefore every signing has an absolute correlation
at least \(0.479n^{3/2}\).

### 3.4 Strictness

The high-side certificate uses only the dual weight counts
\(|\{u\in D_n:\operatorname{wt}(u)=j\}|\), through an absolute coefficient
sum.  It never requests any coset weight distribution.  The low side requests
only the support interval of one coset, not its multiplicities.  Hence
\(L_{\mathrm{KP}}\) cannot reconstruct the full coset histogram.  It also gives
no exact covering radius: any radius within the certified interval remains
possible.

### 3.5 Decisive falsifier

Compute or asymptotically solve the one-dimensional feasibility problem (KP)
using the known weight enumerator of \(D_n\).  If its best certified constant
has the same limit for \(n\equiv2\) and \(6\pmod8\), this architecture is dead.
More structurally, if adding arbitrary linear terms realizes both Arf signs in
every coset sector, then a uniform witness cannot retain the sign alternation.

The magnitude of the Arf Gauss sum is only \(2^{n/2}\), exponentially smaller
than the \(2^n\)-state partition sum.  A fixed \(n^{3/2}\) endpoint gap would
therefore require an exceptionally sharp positivity amplification.  This is
the strongest mathematical reason I assign very low confidence.

### 3.6 Circularity and quantifier audit

* Using the complete coset enumerator in the witness would simply restate the
  optimization.  The certificate must remain radial and use only the dual
  weight enumerator.
* Cancellation depending on \((-1)^{a\cdot u}\) is not uniform in \(a\).
  Absolute coefficient domination pays for every phase but may destroy the
  proposed mod-eight advantage.
* Bent or semibent results for \(RM(1,m)\) concern a different code of length
  \(2^m\); importing their parity gap without an explicit map is invalid.
* The ordinary cut code omits the antipodal word.  A witness for it can prove
  only a one-sided statement.
* A polynomial that is negative only at one endpoint, rather than on the whole
  central interval, does not force an absolute extreme.

## 4. Candidate 3 — odd/even tensor depth of a nonsquare symmetric Hadamard seed

This is independent of the finite-field extension parity above.  Its arithmetic
is the parity of tensor depth and the appearance of an exact Boolean eigenvector
after squaring a symmetric Hadamard transform.

### Candidate card

| Field | Entry |
|---|---|
| Domain | Hadamard matrices, tensor norms, symmetric XOR/Ising games |
| Imported theorems | Existence of a symmetric Hadamard matrix of order 12; Kronecker closure of Hadamards |
| Exact subsequences | low: \(12^{2r+1}\); high: \(12^{2r}\) |
| Target constants | \(c_-=(1-\delta_H/2)/2\), \(c_+=(1-\delta_H/4)/2\), gap \(\delta_H/8\) |
| Proposed mechanism | A nonsquare seed has a finite Boolean/spectral gap.  Even tensor powers acquire an exact sign eigenvector via vectorization; odd powers might retain a uniform tensor gap. |
| Main falsifier | Sign vectors in odd tensor powers with Rayleigh quotient tending to one, or any low-\(Q\) family on the even tower |
| Confidence | truth 0.02%; tractability 0.5%; implication 99.9% |

### 4.1 Exact tensor phenomenon

Fix a symmetric Hadamard matrix \(H\) of order \(12\); Paley's second
construction supplies such a matrix from \(q=5\).  Symmetric Hadamard
constructions and further references are collected in
[Balonin et al., *Construction of Symmetric Hadamard
Matrices*](https://arxiv.org/abs/1708.05098).

Define the finite seed gap

\[
 \beta_H=max_{x\in\{\pm1\}^{12}}
 \frac{|x^THx|}{12\sqrt{12}},\qquad
 \delta_H=1-\beta_H>0.
\]

The strict inequality needs no conjecture: equality in the spectral bound would
give \(Hx=\pm\sqrt{12}\,x\), impossible because the left side is integral and
the right side has nonzero irrational coordinates.

Let

\[
 K_k=H^{\otimes k},\qquad
 S_k=K_k-\operatorname{diag}(K_k),\qquad n_k=12^k.
\]

Then \(S_k\) is an admissible hollow sign matrix.  For even depth \(k=2r\), put
\(G=H^{\otimes r}\), of order \(d=12^r\).  With a consistent vectorization,

\[
 (G\otimes G)\operatorname{vec}(G)
 =d\operatorname{vec}(G).
\]

Since \(\operatorname{vec}(G)\) is a sign vector and \(d=\sqrt{n_{2r}}\),
the even tensor has an exact Boolean spectral eigenvector.  Removing the
diagonal costs at most \(n\), so

\[
 Q(S_{2r})\ge\frac12n_{2r}^{3/2}-O(n_{2r}).
\]

At odd depth no analogous vectorization is forced.  This is a genuine parity
effect in selected tensors.

### 4.2 Exact missing lemma \(L_\otimes\)

> **\(L_\otimes\) (tensor-depth gap plus universal even-depth obstruction).**
> For all sufficiently large \(r\):
>
> 1. every sign vector \(x\in\{\pm1\}^{12^{2r+1}}\) obeys
>    \[
>    |x^TK_{2r+1}x|
>    \le(1-\delta_H/2)\,12^{3(2r+1)/2};
>    \]
> 2. every admissible hollow sign matrix \(A\) of order \(12^{2r}\) obeys
>    \[
>    Q(A)\ge
>    \frac12(1-\delta_H/4)(12^{2r})^{3/2}-o((12^{2r})^{3/2}),
>    \]
>    with a uniform little-oh term.

Part 1, after deleting the diagonal, gives the low-tower construction with
constant \((1-\delta_H/2)/2\).  Part 2 gives the universal high-tower constant
\((1-\delta_H/4)/2\), separated by \(\delta_H/8\).

### 4.3 Imported leverage and missing hypotheses

The Kronecker product of Hadamard matrices is Hadamard, so the spectral and
even-depth eigenvector statements are exact.  What is not imported is any
parallel-repetition theorem for the **same-strategy quadratic Boolean norm**.
Classical values of tensor games can exhibit nonmultiplicative behavior, and
multiplayer XOR games can fail parallel repetition; see
[Pérez-García et al., *Unbounded Violation of Tripartite Bell Inequalities*](https://doi.org/10.1007/s00220-008-0529-7)
for the tensor-norm setting and
[Briët–Vidick, *Explicit Lower and Upper Bounds on the Entangled Value of
Multiplayer XOR Games*](https://doi.org/10.1007/s00220-012-1461-z).
Neither paper proves the desired odd-depth gap.

More seriously, the exact sign eigenvector of the selected even tensor says
nothing about the minimum over all even-order sign matrices.  An inverse theorem
would have to say that any matrix beating the proposed high threshold is
impossible specifically at \(12^{2r}\).  No Hadamard theorem has this force.

### 4.4 Strictness

The odd clause evaluates one explicit family; the even clause gives one scalar
inequality for all matrices at one sparse family of orders.  It identifies no
optimizer and no coset distribution.  Infinitely many distinct minima and
histograms are consistent with both clauses, so the lemma is strictly weaker
than full minimization.

### 4.5 Decisive falsifiers and quantifier audit

* Find \(x_r\in\{\pm1\}^{12^{2r+1}}\) with
  \(|x_r^TK_{2r+1}x_r|/n_r^{3/2}\to1\).  This kills the retained seed gap.
* Construct any family \(A_r\) at orders \(12^{2r}\) below the stated universal
  threshold.  This kills the high side.
* A bilinear bound for \(y^TK_kz\) is not the same-strategy quadratic bound.
* Product sign vectors are not all sign vectors on a product set; proving the
  estimate only for \(x=x_1\otimes\cdots\otimes x_k\) reverses the maximization
  quantifier.
* The even tensor's high value is an \(\exists A\) statement, while the needed
  lower theorem is \(\forall A\).

This candidate has the cleanest exact parity algebra and the weakest connection
to the global optimum.

## 5. Cross-candidate no-go observations

1. **Conference existence alone is asymptotically too dense.**  Paley orders
   from primes in \(1\pmod4\) approach every large order multiplicatively, and
   principal restriction preserves \(Q\).  A conference/nonconference residue
   dichotomy cannot yield a fixed gap.
2. **Witt or parity integrality by itself is too small.**  Off-diagonal Gram
   products being odd instead of even changes fourth spectral moments by lower
   order terms.  Mod-four energy congruences change extrema by \(O(1)\), not
   \(\Theta(n^{3/2})\).
3. **Maximum-excess constructions usually optimize the wrong channel.**
   Independent row/column switches solve a bilinear problem.  Symmetric Seidel
   switching requires the diagonal \(x=y\), and the absolute value requires both
   antipodal channels.
4. **Selected high-energy matrices never lower-bound \(M_n\).**  The most
   tempting arithmetic facts above produce a crystalline state in one chosen
   matrix.  Since \(M_n\) is a minimum over matrices, that is evidence in the
   wrong direction unless accompanied by a universal rigidity or covering
   theorem.
5. **Approximate nonexistence needs a macroscopic gap.**  Bruck–Ryser–Chowla,
   determinant congruences, and exact ETF/conference nonexistence generally
   exclude equality.  Deleting \(o(n)\) vertices from a nearby construction
   shows that exclusion of equality alone cannot supply a fixed normalized
   separation.

## 6. Ranked conclusion

1. **Fixed-characteristic Paley crystal/glass alternation** — strongest explicit
   mechanism.  It has an exact \(12/25\) even-extension Boolean state and a
   physically meaningful odd/even distinction.  Confidence in the full lemma:
   **0.2%**.
2. **Arf/Krawtchouk covering certificate** — the only proposal here whose high
   side is natively an all-cosets theorem.  The Arf signal is probably erased by
   arbitrary coset phases or is too small for a fixed endpoint gap.  Confidence:
   **0.05%**.
3. **Nonsquare Hadamard tensor-depth parity** — exact selected-family parity but
   essentially no leverage toward a universal high-subsequence bound.
   Confidence: **0.02%**.

Accordingly, my ledger-blind recommendation is **do not treat any candidate as a
probable nonconvergence route**.  If one test is funded, test Candidate 1 in this
order: (i) estimate the odd-degree Paley ground-state constant along
\(q=5^{2r+1}\); (ii) independently search unrestricted matrices at the first
even tower orders for values below \(0.479\).  The second test addresses the
actual optimum and is more decisive than further analysis of selected Paley
examples.
