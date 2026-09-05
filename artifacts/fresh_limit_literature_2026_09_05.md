# Fresh primary-literature audit, 5 September 2026

## Outcome and normalization

No imported theorem found here proves convergence or nonconvergence of

\[
q(A)=\max_{x\in\{\pm1\}^n}\left|\sum_{i<j}a_{ij}x_ix_j\right|,
\qquad M_n=\min_A q(A),\qquad c_n=M_n/n^{3/2}.
\]

The external architecture selected **before inspecting the local archive** was
near-order Hadamard dimension filling. Its missing ingredient is preservation
of an arbitrary low scalar optimum, not density of available orders. The
closest recent theorem actually proving an optimal-value limit is Hilbertian
and trilinear, with different quantifiers and spin domain. An explicit
five-cycle calculation below independently rejects the simplest proposed
scalar tensor transfer. The archive subsequently showed that the broad
obstruction and older KSZ/Sidon scope issues were already known locally.

The exact code identity supplied for this search is
\(M_n=\binom n2-2\rho(C_n^+)\), where
\(C_n^+=\{(c+b_i+b_j)_{i<j}\}\). Neither an ordinary cut-code radius nor
the radius of a cycle code of a triangular graph can be substituted.

## Primary-source map: actual theorem scope

1. **Kalai--Schulman, quasi-random multilinear polynomials.** Theorem 1,
   pp. 1--2, gives, for fixed maximum degree \(d\) and support \(U\),
   \(\inf_{h\in L_{U,\mathbb R}}\|h\|_\infty
   =\Theta_d(\sum_i\sqrt{u_i})\), where \(u_i\) counts incident monomials.
   Taking all pairs is exactly our coefficient-sign minimization and gives
   the right scale. It neither identifies a leading constant nor relates
   optimal values at different orders. [Primary preprint](https://arxiv.org/pdf/1804.04828).

2. **Pellegrino--Raposo, JFA 282 (2022), 109293.** Theorem 1.1 fixes
   multilinear order \(d\) and \(\epsilon>0\); for all sufficiently large
   side dimensions \(n_1,\ldots,n_d\), it constructs a unimodular form on
   independent \(\ell_\infty^{n_i}\) inputs with norm at most
   \((1+\epsilon)\max_i\sqrt{n_i}\prod_i\sqrt{n_i}\).
   The proof first chooses nearby Hadamard orders, builds an orthogonal-chain
   form, and restricts coordinates. It does not accept a prescribed
   low-norm seed. Even its bilinear switching application retains a gap
   between lower and upper constants. [Theorem and proof](https://arxiv.org/pdf/2006.12892).

3. **Pellegrino--Raposo, August 2026 Hilbertian KSZ preprint.** Theorem 5.8
   says that if \(1\le m_n\le n\) and
   \(\limsup r_n/\log_2 n<2\), the minimum norm \(\Lambda\) of real
   unimodular trilinear forms on
   \(\ell_2^{r_n}\times\ell_2^{m_n}\times\ell_2^n\) satisfies
   \(\Lambda(r_n,m_n,n)=(1+o(1))\sqrt n\).
   Theorem 3.1 characterizes exact equality in the square case by a
   Hadamard matrix of order \(n\) and \(r\le\rho(n)\), with
   \(\rho\) the Hurwitz--Radon function. Theorem 3.3's arithmetic gap is
   \((1+2/n^3)^{1/4}\), which vanishes. Consequently it supplies neither
   a scalar diagonal limit nor separated subsequential constants here.
   [Actual Theorem 5.8 and proof](https://arxiv.org/html/2608.08246v1#S5.SS3).

4. **Defant--Galicer--Mansilla--Mastylo--Muro, IMRN 2024.** Theorems 4.1
   and 4.3 compute projection-constant limits, in particular
   \(\lambda(\mathcal B^n_{=2})/n\to\sqrt{2/(\pi e)}\).
   This is an average absolute value of the all-positive Walsh kernel,
   not a minimum signed supremum. Theorems 5.1 and 5.4 compare ordinary
   Sidon/Gordon--Lewis/unconditional constants only within fixed factors.
   Our equal-modulus quantity is
   \(U_n=\binom n2/M_n\le\operatorname{Sid}(\mathcal B^n_{=2})\).
   Fixed-factor comparisons cannot transfer a limit or flatten arbitrary
   coefficient magnitudes. [Primary paper](https://arxiv.org/pdf/2302.00233).

5. **Friedland--Lim, symmetric Grothendieck inequality.** Theorem 4.3
   applies to every real symmetric matrix and vectors in dimension at least
   \(n\); the absolute vector quadratic maximum is at most a universal
   constant times the same-spin Boolean maximum. The stated real upper
   constant is \(\sinh(\pi/2)\). This formalizes the two different norms,
   but is not asymptotic equality on minimizing sign matrices.
   [Primary full paper](https://www.stat.uchicago.edu/~lekheng/work/sgi.pdf).

6. **Christoph--Gishboliner--Krivelevich, February 2026.** Theorem 1.2
   fixes \(\epsilon>0\), takes an \(n\)-vertex \(d\)-regular graph
   \(F\) with \(d\le(1-\epsilon)n\), and forces a copy with discrepancy
   at least \(c\sqrt{\epsilon d}\,n\) in every two-coloring of
   \(K_n\), for absolute \(c>0\). Balanced complete bipartite \(F\)
   gives the relevant order of magnitude, not a sharp constant. Their
   exact factor constants concern fixed clique factors, whose discrepancy
   is linear rather than our dense quadratic-cut scale.
   [Primary theorem](https://arxiv.org/html/2602.04069v1).

7. **Best, 1976 report / 1977 paper.** The excess \(\sigma(n)\) is a
   **maximum over all Hadamard matrices**, not a minimum over switching
   classes. Theorems 1--2 give \(\sigma(n)\le n\sqrt n\), with equality
   precisely when a regular Hadamard matrix exists. The random-row proof
   in Section 5 also gives a per-class independent-switch lower bound
   asymptotic to \(\sqrt{2/\pi}\,n^{3/2}\). It explicitly does not prove
   an approximately flat Boolean image for every class. Its tensor
   inequality is supermultiplicativity of maximum excess, the wrong
   direction and extremization for a minimax transfer.
   [Original CWI report](https://ir.cwi.nl/pub/6858/6858D.pdf).

8. **Momihara--Suda, Integers 17 (2017), A30.** Theorem 2 assumes a
   particular two-intersection set in the Paley translation configuration
   for \(q=4m^2+1\), and constructs an order-\(q+1\) conference matrix
   with row sums \(2m-1\) and \(2m+1\). Theorem 4 supplies these sets
   when the underlying prime is \(1\pmod4\). In the actual Theorem 2
   proof, columns in \(D\) are negated, then rows in the dual class
   \(D_\alpha^\perp\) are negated. These sets need not agree. This is
   not a theorem about \(DCD\), and cannot prove universal same-switch
   Boolean saturation of symmetric conferences.
   [Primary proof, printed pp. 5--6](https://math.colgate.edu/~integers/r30/r30.pdf).

9. **Haemers, MATCH 68 (2012), 653--659.** Seidel energy is exactly
   \(\|A\|_*\). Theorem 3.1 proves
   \(\|A\|_*\le n\sqrt{n-1}\), with equality iff \(A\) is symmetric
   conference; Corollary 3.2 fills orders to obtain the asymptotic maximum
   nuclear norm. Its proof is trace-square plus Cauchy--Schwarz. There is
   no comparison there with the maximum signed quadratic excess.
   [Original published paper](https://match.pmf.kg.ac.rs/electronic_versions/Match68/n3/match68n3_653-659.pdf).

Two additional leads were not imported: the actual Esmaeili--Zaghian
2009 full text could not be obtained from its official volume entry;
Craigen--Kharaghani's 2004 weaving paper was located at its publisher but
the full theorem text was unavailable. Abstracts and secondary sequence
tables are not theorem authority. The former is a known direct predecessor,
not a newly validated formula for \(\rho(C_n^+)\).

## Frozen architecture and the exact missing transfer

The near-order idea is mathematically adequate if scalar values can be
retained. For example, let
\(\mathcal S=\{4^a144^b:a,b\ge0\}\). Its consecutive ratios tend to
one: \(\log144/\log4\) is irrational, and sufficiently long finite sets
of its multiples are dense modulo one. Symmetric regular Hadamard matrices
exist at these orders. To see the nontrivial base directly, start with a
Hadamard \(H\) of order \(h=12\), put
\(K_{(i,j),(k,l)}=H_{il}H_{kj}\), and set \(z_{ij}=H_{ij}\).
Then \(K=K^T\), \(K^2=h^2I\), and \(Kz=hz\). Thus
\(\operatorname{diag}(z)K\operatorname{diag}(z)\) is symmetric regular
Hadamard of order 144. Tensor products give \(\mathcal S\).

A concrete sufficient missing theorem is the following **near-minimizer
block theorem**, with an absolute constant \(C\): for every sufficiently
large \(n\) and an exact minimizing seed \(A_n\), one can choose, for all
large \(s\in\mathcal S\), exact sign cross-blocks \(B_{ij}\) with
\(B_{ji}=B_{ij}^T\) such that the hollow-block cross Hamiltonian obeys

\[
q(B_{\rm cross})\le s^{3/2}M_n+Cns^{3/2}
                 +o_s((ns)^{3/2}).                       \tag{T}
\]

Filling the \(n\) diagonal blocks with independent minimizing order-\(s\)
signings costs at most \(nM_s=O(ns^{3/2})\). Hence (T) yields
\(\limsup_{s\in\mathcal S}c_{ns}\le c_n+O(n^{-1/2})\).
For arbitrary large \(N\), take the next \(ns\ge N\); the relative gap
tends to zero. Principal restriction cannot increase \(q\): fill removed
spins uniformly and use conditional expectation. Therefore
\(\limsup_Nc_N\le c_n+O(n^{-1/2})\). Taking \(n\) along a liminf
sequence proves convergence.

This is a conditional deduction, not a proof of (T). Published KSZ
dimension filling does not supply (T), even when the order arithmetic is
already solved. If regular blocks are required to preserve the seed on
constant fiber spins, their row and column sums must be
\(a_{ij}\sqrt s\); that constraint alone does not bound other fiber spins.

## Direct tests of the missing scalar step

### A five-cycle defeats the zero-error common tensor claim

Index \(A\) by \(\mathbb Z/5\mathbb Z\), with \(a_{ij}=-1\) at cyclic
distance one and \(+1\) at distance two. All row sums are zero. After
global spin negation at most two spins are negative. The energy is zero
when at most one is negative and is \(4a_{ij}\) when exactly \(i,j\)
are negative. Thus \(q(A)=4\).

Let \(H_4(i,j)=(-1)^{\operatorname{popcount}(i\mathbin\&j)}\), indexed
from zero, and use the five fiber spins

\[
(-,+,-,-),\quad (+,-,-,+),\quad (+,-,+,-),\quad
(-,-,+,+),\quad (+,+,+,+).
\]

For **every** macro-edge \(i<j\), direct integer multiplication gives
\(a_{ij}x_i^TH_4x_j=4\). Consequently

\[
q(A\otimes H_4)\ge40>4^{3/2}q(A)=32.                \tag{1}
\]

This is a partial-sign cross-block matrix, not a claimed order-20
admissible signing. The failure already occurs before filling diagonal
blocks. Since \(H_4\) is congruence-switchable to a symmetric regular
Hadamard, regularity does not repair (1). This finite counterexample does
not refute (T), which permits \(O(n)\) seed overhead and only needs
large minimizing seeds.

### Every vector relaxation already pays the half-constant benchmark

Define
\(V(A)=\max_{\|u_i\|=1}|\sum_{i<j}a_{ij}\langle u_i,u_j\rangle|\).
Set \(r=\sqrt{n-1}\) and
\(X_\pm=\tfrac12(I\pm A/r)^2\). Both are positive semidefinite with
diagonal one, so they are feasible Gram matrices. Since
\(\operatorname{Tr}A=0\) and \(\operatorname{Tr}A^2=n(n-1)\),

\[
\frac12\operatorname{Tr}(AX_\pm)
=\pm\frac n2\sqrt{n-1}
 +\frac{\operatorname{Tr}A^3}{4(n-1)}.
\]

The larger absolute value is at least \(n\sqrt{n-1}/2\). Thus

\[
V(A)\ge\frac n2\sqrt{n-1}\qquad\hbox{for every signing }A. \tag{2}
\]

Conference matrices attain equality by the spectral upper bound.
Therefore an amplification argument controlled only through this relaxation
cannot preserve an eventual constant below \(1/2\). For the five-cycle,
\(V(A)=5\sqrt5/2>4\): its nonzero eigenvalues are
\(\pm\sqrt5\), each twice, and the spectral projectors have constant
diagonal. This does not rule out a genuinely scalar near-minimizer theorem.

## Follow-up: nuclear norm and universal conference saturation

Root and the variational agent independently obtained the stronger
polar-Gram/nuclear inequalities during this search. Targeted primary
searches for Seidel energy versus cut excess, polar decomposition versus
quadratic rounding, and nuclear norm versus arcsine bounds did not locate
their displayed inequality. This is **not** a novelty certification.
Haemers' exact nuclear-norm extremality theorem above is relevant but does
not contain the Boolean tradeoff.

Nor was a primary theorem found asserting
\(\max_x|x^TCx|=(1-o(1))n\sqrt{n-1}\) for **every** symmetric
conference matrix \(C\). Regular/square-field constructions establish
special families; independent row/column excess results do not answer that
same-switch question. Even the Hadamard-class analogue must not be inferred
from a theorem maximizing excess over all matrices.

## Archive check after freezing the mechanism

The selected mechanism is not a new route. The local
`artifacts/unimodular_walsh_sidon_literature.md` already records the
projection/Sidon mismatch and the 2022 KSZ scope. Ledger Sections 3.4,
3.6, and 3.8 already reject universal Hadamard/vector amplification and
tensor submultiplicativity. `artifacts/independent_regular_hadamard_lift.md`
has a stronger growing-order obstruction to making the orthogonal fiber
channel negligible. The five-cycle witness above is an independently
checked compact instance, not grounds to reopen the route.

No existing citation to `2608.08246` was found in the searched ledger,
artifacts, drafts, and audits. Its exact Hilbertian all-order theorem is a
useful current scope check, but supplies no missing scalar transfer. The
rigorous output remains: audited source boundaries, the conditional
dimension-filling deduction, and explicit route tests; no changed bound
and no conclusion about convergence or nonconvergence.

## Targeted second phase: same-spin tensor regularization

The algebra agent fixed the regular generators
\(H_4=J_4-2I_4\) and a symmetric regular \(H_{144}\), and defined
\[
R(A)=\sup_{a,b\ge0}\frac{q(H_4^{\otimes a}\otimes
 H_{144}^{\otimes b}\otimes A)}{(4^a144^b)^{3/2}}.
\]
The following source checks were made specifically against this definition.

### What the actual XOR and completely bounded theorems say

* Helton--Mousavi--Nezhadi--Paulsen--Russell,
  [*Synchronous Values of Games*](https://arxiv.org/pdf/2109.14741),
  Theorems 2.2 and 6.7, identify the synchronous classical optimization
  with one Boolean assignment and the synchronous quantum optimization
  with one correlation matrix. A symmetric hollow cost matrix therefore
  maps exactly to our one-sided quadratic objective; taking both cost
  signs gives the absolute objective. Definition 6.10 requires the
  optimal diagonal dual majorant to dominate both signs of the cost.
  Theorem 6.12 proves tensor multiplicativity for the synchronous quantum
  bias under this **balancedness** hypothesis. Example 6.9 shows that
  multiplicativity fails without it. These are vector-SDP statements,
  not Boolean-Hadamard realization theorems. Also, the negative diagonal
  of our \(H_4\) is not itself a synchronous-game cost matrix under their
  definition; any application to that factor must use a separately
  justified algebraic extension.
* Briët--Escudero Gutiérrez--Gribling,
  [*Grothendieck Inequalities Characterize Converses to the Polynomial
  Method*](https://homepages.cwi.nl/~jop/qpoly_add_Version5.pdf),
  Proposition 3.5, proves that a real quadratic form has the completely
  bounded norm of its unique symmetric coefficient matrix. The matrix
  norm still optimizes noncommuting contractions, equivalently the
  ordinary independent-family vector bias. Proposition 3.2's
  same-contraction formulation does not turn those operators into
  commuting Boolean spins. For \(p_A(x)=x^TAx/2\), this norm is
  \(\Gamma(A)/2\), where
  \(\Gamma(A)=\max_{u_i,v_j}|\sum_{ij}a_{ij}\langle u_i,v_j\rangle|\).
* Cleve--Slofstra--Unger--Upadhyay,
  [*Perfect Parallel Repetition Theorem for Quantum XOR Proof
  Systems*](https://arxiv.org/pdf/quant-ph/0608146), Theorems 1 and 4,
  give multiplicativity and the vector description for the ordinary
  **two-family** quantum bias. Combining this with the real
  Grothendieck inequality gives
  \[
  \Gamma(A)^k/K_G\le\|A^{\otimes k}\|_{\infty\to1}
       \le\Gamma(A)^k,
  \]
  hence a self-tensor, \(k\)-th-root limit equal to \(\Gamma(A)\).
  Neither the fixed ancillary Hadamard factors nor the unrooted
  normalization defining \(R\) occur in this conclusion.

The first source's relevant hypotheses and proof were read on pages
29--32; the second source's matrix/polynomial identification was read on
pages 12--14. No theorem identifying \(R\) was imported.

### A strict obstruction to either standard vector identification

Let \(A\) be the five-cycle Seidel matrix used above. Then
\[
q(A)=4,\quad \|A\|_{\infty\to1}=8,\quad
V(A)=\Gamma(A)/2=5\sqrt5/2.
\]
The vector equalities follow from \(A^2=5I-J\), constant-diagonal
spectral projectors, and the spectral upper bound. In contrast, the
absolute PSD-majorant functional
\[
T(A)=\frac12\min_{B\succeq\pm A}\max_xx^TBx
\]
satisfies, by the feasible majorant \(|A|=\sqrt5(I-J/5)\),
\[
R(A)\le T(A)\le\frac{12\sqrt5}{5}
 =\frac{24}{25}\,V(A)<V(A).
\]
In fact equality holds in the displayed bound for \(T\), either by
dihedral averaging plus the self-complementing permutation, or by the
algebra agent's sign-column lower bound. For the **actual** generator
\(J_4-2I_4\), exhaustive integer evaluation gives
\(q(H_4\otimes A)=36\), so
\[
4.5\le R(A)\le12\sqrt5/5<5\sqrt5/2.
\]
Thus \(R\) is neither the standard same-spin vector norm nor the
quadratic polynomial completely bounded norm. This is an exact
same-spin separation, not an inference from bipartite factorization.
The earlier Walsh-four witness in this note uses a different congruence
class of symmetric regular Hadamard matrix and has cap 40, not 36.

### An exact dual problem for the remaining realization question

Finite-dimensional SDP duality gives the following useful formulation:
\[
T(A)=\frac12\max\{\operatorname{Tr}(AW):
 Z\in\operatorname{conv}\{xx^T:x\in\{\pm1\}^n\},\quad
 -Z\preceq W\preceq Z\}.                 \tag{3}
\]
To check constants, minimize \(t/2\) subject to
\(B\succeq\pm A\) and \(\langle xx^T,B\rangle\le t\).
The latter multipliers sum to \(1/2\); the two PSD dual variables
\(X,Y\) satisfy \(X+Y=Z/2\), and \(W=2(X-Y)\).
Strict primal feasibility is immediate from sufficiently large \(B=t_0I\)
and a larger epigraph variable, so there is no duality gap.

Writing a tensor spin assignment as Boolean functions \(f_i\) on the
uniform outer index set, its empirical Gram matrix belongs to the cut
polytope in (3), and
\(W_{ij}=\langle f_i,(H_s/\sqrt s)f_j\rangle\) obeys
\(-Z\preceq W\preceq Z\). Consequently, an identity \(R=T\)
would require an asymptotic Boolean-function realization of these
cut-Gram/contraction pairs by this particular family of Hadamard
operators. Ordinary elliptope realization is insufficient. For the
four-letter generator, \(H_4/2=2\mathbb E-I\); its product action
changes the sign of odd Efron--Stein levels. This is not a positive
Markov noise operator. No universality theorem for this exact
realization problem was found or assumed.

### The square Gale--Berlekamp calibration

Pellegrino--Raposo,
[*Upper Bounds for the Constants of Bennett's Inequality and the
Gale--Berlekamp Switching Game*](https://arxiv.org/pdf/2111.00445),
Section 5, explicitly records for the exact square bilinear minimum
\(G_n=\min_{A\in\{\pm1\}^{n\times n}}\|A\|_{\infty\to1}\)
the asymptotic bounds
\[
\sqrt{2/\pi}+o(1)\le G_n/n^{3/2}\le1+o(1).
\]
The paper proves improved uniform upper estimates; it does not assert
convergence of the normalized minima. No primary theorem proving that
limit was located. This is a search result, not a claim that every
possible publication has been excluded. Quantum self-tensor
multiplicativity does not fill this all-order minimum gap.

An exact archive search after freezing these source mappings found the
Briët--Escudero Gutiérrez--Gribling paper already listed in
`artifacts/retrieval_panel_2026_08/banach_toolkit.md`; it also found the
general warning about full self-tensor regularizations there. No local
hit for the synchronous-values paper's identifier or balanced-XOR
mechanism was found. Neither source changes the original bounds.

## Adversarial audit of the regularized-limit and nuclear proofs

The full proofs in `artifacts/fresh_limit_algebra_2026_09_05.md` and
`artifacts/fresh_limit_variational_2026_09_05.md` were independently
checked after the literature phase. The following qualifications and
improvements matter.

1. **The hollow regularized-limit theorem is valid.** Constant outer
   spins make the normalized tensor values coordinatewise increasing;
   the spectral bound makes them finite. Shifting a cofinal pair of
   indices proves exact tensor stability, with no interchange of a
   minimum and supremum. Principal restriction is legitimate on the
   hollow domain by conditional expectation. Irrationality of
   \(\log144/\log4\), proved by prime factorization, gives a finite
   residue net and hence asymptotically vanishing multiplicative gaps.
   One fixes the seed order before using this mesh and only afterwards
   takes a liminf sequence of seed orders. The resulting limsup/liminf
   deduction is sound.
   At any **fixed** seed order, the finite minimum does commute with
   the directed supremum: if finitely many functions \(f_s(A)\) are
   coordinatewise increasing, choose for each \(A\) an index where it
   is within \(\varepsilon\) of its supremum and then a common larger
   index. Thus \(\min_A\sup_s f_s(A)=\sup_s\min_A f_s(A)\).
   The remaining exchange obstruction concerns unbounded seed orders
   or the restricted tensor family, not this finite-order identity.
2. **Diagonal extension requires care, but is valid.** On all symmetric
   real matrices, \(q\) is merely a seminorm: trace-zero diagonal
   matrices are invisible. Nevertheless \(R\ge\|A\|_{\infty\to1}/2\)
   still holds by the same two Boolean eigenvectors of \(H_4\), so
   \(R\) is a genuine norm. In particular
   \(R(\operatorname{diag}d)=\|d\|_1/2\), with the upper bound supplied
   by the majorant \(\operatorname{diag}|d|\). Ordinary \(q\)-principal
   monotonicity fails with diagonals. For \(R\), its conditional
   expectation error is
   \(|\operatorname{Tr}H_s|\,|\operatorname{Tr}A_{I^c}|/(2s^{3/2})\),
   which tends to zero cofinally: every \(H_s\) has constant diagonal
   \((-1)^a\), so \(|\operatorname{Tr}H_s|=s\).
3. **The old completion estimate was correct but avoidable.** For the
   particular fill \((H_s-\operatorname{diag}H_s)\otimes I_n\), its
   leading \(ns^{3/2}/2\) norm is already forced by the all-one witness.
   Sharpening its operator bound cannot remove the leading
   \(1/(2\sqrt n)\) normalized seed payment. Root's full-sign-seed
   upgrade does remove it: if \(B\) is a symmetric full-sign seed,
   \(H_s\otimes B\) is already full-sign, and erasing only its diagonal
   costs exactly \(sn/2\) in \(R\). This is negligible after dividing
   by \((sn)^{3/2}\). Hollow and full-sign normalized minimum values
   differ by at most \(1/(2\sqrt n)\), so the strengthened identity is
   \[
   c_R=\inf_{n\ge1}\min_{B\text{ symmetric full-sign}}
          R(B)/n^{3/2}.
   \]
   The lower direction is the definition of infimum plus convergence;
   the upper direction fixes a full-sign seed, tensors, erases the
   diagonal, and deletes to the nearest target order. No min--sup
   exchange is involved.
4. **The full-sign majorant floor is also valid.** If
   \(D\succeq\pm B\), every column of full-sign \(B\) is Boolean.
   With \(C=\max_xx^TDx\),
   \(nC\ge\operatorname{Tr}(DB^2)\ge\sum_j|\lambda_j(B)|^3\ge
   n^{5/2}\), since \(\sum\lambda_j^2=n^2\).
   Thus \(T(B)\ge n^{3/2}/2\). An identification \(R=T\) would force
   \(c_R=1/2\), but would still not identify the original liminf.
5. **The polar/nuclear theorem checks out, including the arcsine
   strengthening.** All \(h_i=|A|_{ii}\) are positive, and
   \(D^{-1/2}(|A|\pm A)D^{-1/2}\) really have diagonal one. Their
   paired vector objectives give the required factor \(1/2\).
   For Gaussian rounding, \(p\pm q\in[-1,1]\) implies
   \(|p|+|q|\le1\), and
   \(\arcsin(p+t)-\arcsin(p-t)\ge2\arcsin t\) on that domain.
   Applying this before taking the difference of the two expected
   Boolean energies gives
   \[
   q(A)\ge\frac{n(n-1)}\pi
       \arcsin\!\left(\frac n{\|A\|_*}\right),\qquad
   V(A)\|A\|_*\ge\frac{n^2(n-1)}2.
   \]
   Both Jensen steps have the correct direction: \(\arcsin\) is
   convex on \([0,1]\), and \((u,v)\mapsto(uv)^{-1/2}\) is jointly
   convex. The argument includes the endpoint \(n=2\). The improved
   spectral stability constant \(2\varepsilon/(1+\varepsilon)\)
   follows immediately. These checks do not certify novelty.

### A finite, noncircular landing statement weaker than identifying R

Let \(m_n\) minimize \(q\) over symmetric **full-sign** matrices and
let \(\mathcal M_n\) be its finite set of minimizers. For
\(j\in\{4,144\}\), define the finite one-step excess
\[
e_j(n)=\min_{B\in\mathcal M_n}
 \frac{q(H_j\otimes B)}{(jn)^{3/2}}-
 \frac{m_n}{n^{3/2}}\ge0.
\]
It would suffice to prove
\[
\max\{e_4(n),e_{144}(n)\}\le\omega(n),\qquad
\sum_{k\ge0}\omega(4^k n)\longrightarrow0,             \tag{4}
\]
for some nonincreasing nonnegative \(\omega\). Indeed, each finite
minimum supplies \(m_{jn}/(jn)^{3/2}\le m_n/n^{3/2}+\omega(n)\).
Reselect a minimizer at each new order; telescope along any word in the
two generators; bound the total loss by the tail in (4); then use the
dense mesh and diagonal erasure/principal deletion. Taking a liminf
sequence of starting orders proves convergence of the original sequence.
Rates \(O(n^{-\delta})\) or
\(O((\log n)^{-1-\delta})\), with \(\delta>0\), satisfy (4).

This demands only two finite lifts of an appropriately selected
minimizer at each order. The selected minimizers may differ between
the two generators and between successive orders. No one seed must
remain good through every tensor power, no Boolean-Hadamard universality
is assumed, and no regularized norm need be computed. Thus (4) is a
concrete weaker landing target than requiring \(R-q=o(n^{3/2})\) along
fixed seeds. It is a sufficient statement, **not** a theorem proved in
this campaign and not a claim of logical necessity.

## Final targeted source checks: factorization, AMP, odd Walsh

### The correct Hilbert-factorization convention

Root and the variational agent subsequently derived the exact identities
\[
R(B)=\frac12\sup_s
  \frac{\|H_s\otimes B\|_{\infty\to1}}{s^{3/2}},\qquad
2T(B)=\Gamma_2(B:\ell_\infty^n\longrightarrow\ell_1^n),
\]
where \(\Gamma_2\) denotes factorization **as an operator through a
Hilbert space**, not the entrywise matrix factorization norm commonly
denoted \(\gamma_2\) in discrepancy or communication complexity.
The first follows by applying \(R\ge\|\cdot\|_{\infty\to1}/2\)
to each tensor descendant, then using tensor stability; its opposite
inequality is immediate from \(q\le\|\cdot\|_{\infty\to1}/2\).
For the second, a factorization \(B=P^TQ\) with balanced factor norms
gives the majorant \((P^TP+Q^TQ)/2\succeq\pm B\). Conversely,
\(D\succeq\pm B\) gives \(B=D^{1/2}CD^{1/2}\) with a
self-adjoint contraction \(C\). These are legitimate exact same-spin
reductions; they do not make \(T\) the standard vector SDP.

A directly relevant primary warning is Aicke Hinrichs,
[*Hilbert Space Factorization and Fourier Type of Operators*](https://www.impan.pl/shop/en/publication/transaction/download/product/89984),
Studia Mathematica 145 (2001), 199--212.
His \(\kappa(T\mid U_N)\) is the outer
\(\ell_2^N(X)\to\ell_2^N(Y)\) norm of \(U_N\otimes T\).
Theorem 1.2 characterizes Hilbert-factorable operators using all
contraction matrices. Theorem 1.3 proves that for **any** prescribed
sequence of \(N\)-by-\(N\) contractions, there are operators \(T_N\)
with \(\kappa_N(T_N)=1\) but
\(N^\alpha\kappa(T_N\mid U_N)\to0\) for every \(\alpha<1/10\).
Theorem 6.1 includes Walsh tests; Theorem 5.2 requires exponentially
many fixed-size tests for uniform equivalence. The counteroperators
act \(\ell_1^N\to\ell_\infty^N\), and the test dimension is tied
to \(N\). Thus these results do **not** refute our fixed-seed,
unbounded-ancilla, outer \(L_\infty\to L_1\) identity conjecture.
They also supply no such identity. The paper's Fourier-type/operator
factorization open question should not be reported as a verified
current open problem without a later-literature status check.

No local archive hit for this exact Hinrichs paper was found after
freezing its mechanism. This is an audited boundary, not a new
convergence route.

### Averaged fourth-power AMP replacement not supplied

The following actual hypotheses were checked in response to root's
question whether \(n^{-1}\sum_{ij}U_{ij}^4\) could control an AMP
error linearly in spectral/nuclear defect.

* Dudeja--Lu--Sen,
  [*Universality of Approximate Message Passing with Semi-Random
  Matrices*](https://arxiv.org/pdf/2204.04281), Definition 1 and
  Theorem 1, assume entrywise \(O(n^{-1/2+\varepsilon})\)
  delocalization for every fixed \(\varepsilon>0\), bounded operator
  norm, and uniform approximate row orthogonality and row norms.
  Lemma 1 includes sign-conjugated symmetric orthogonal matrices with
  that same maximum-entry bound. The theorem does not replace it by
  an averaged fourth-power condition.
* Wang--Zhong--Fan,
  [*Universality of Approximate Message Passing Algorithms and
  Tensor Networks*](https://arxiv.org/pdf/2206.13037), version 5,
  Definition 2.6 and Theorem 2.8, use maximum off-diagonal bounds for
  all fixed diagonal monomials, sign/permutation conjugation, and a
  matching limiting diagonal distribution. Proposition 2.7 gives
  simpler sufficient conditions, still of maximum-entry type.
* Gorini--Jones--Kunisky--Pesenti,
  [*Universality of First-Order Methods on Random and Deterministic
  Matrices*](https://arxiv.org/html/2604.11729v1), Theorem 5.3,
  bounds each fixed noncactus diagram by
  \(O_\alpha(\varepsilon+n^{-1/2})\), but \(\varepsilon\)
  controls **uniform open-cactus entries** under Assumption 5.2.
  The full traffic bound also contains
  \(n^{-1/2}(1+\varepsilon\sqrt n)^{O_\alpha(1)}\).
  Theorem 5.1 for orthogonal matrices requires
  \(\max_{ij}|U_{ij}|\le n^{-1/2+o(1)}\) and puncturing.
  No averaged-fourth-power replacement is stated.

Small average squared distance to a flat sign matrix does not directly
verify these maximum-entry conditions. Deriving a diagram estimate in
an averaged fourth-power parameter remains an additional proof task.
The first two matrix-universality families and the 2026 paper were
already discussed in
`artifacts/finite_temperature_universality_relaxation_audit.md`.
A downloaded primary PDF used for this check is
`tmp/fresh_gjkp_2026.pdf`; no local PDF-to-text utility was available,
so the actual theorem was read in the arXiv HTML instead.

### Exact primary match for the odd-Walsh same-spin problem

Carlet--Danielsen--Parker--Solé,
[*Self-Dual Bent Functions*](https://doi.org/10.1504/IJICOT.2010.032864),
International Journal of Information and Coding Theory 1 (2010),
384--399, Section 3, studies exactly
\(S(F)=F^TW_rF\), where \(W_r=H_2^{\otimes r}\).
The [author-uploaded full manuscript](https://www.researchgate.net/publication/210230164_Self-dual_bent_functions)
was read because the indexed university copies failed to load.
Theorem 3.1 identifies spectral equality at even \(r\) with
self-dual or anti-self-dual bent functions. For odd \(r\),
Theorem 3.3 gives normalized value \(1/\sqrt2\).
Theorem 3.4 proves monotonicity under
\[
F'=(F,G,G,-F),\qquad G=\operatorname{sign}(W_rF),
\]
or the swapped variant, using the exact identity
\(S(F')=8\|W_rF\|_1\ge8|S(F)|\).
The paper does not prove that the resulting limit is one or determine
the maximum odd-arity asymptote. Its Table 1 gives non-exhaustive
same-spin witnesses, reaching normalized value 0.933386 at \(r=25\);
only its \(r=3\) search is exhaustive. These cannot be promoted to
optimal caps or a limit theorem. This paper was already cited in
`artifacts/walsh_bent_stability_literature_audit.md`.

The exact odd-Walsh problem is therefore a genuine published neighbor
of the full seed \(H_2\) regularization, with a directly relevant
monotone construction. No primary theorem establishing asymptotic
spectral saturation was found. A gap for this one Walsh trajectory
would still need control of the additional order-144 factors before
it could bound the actual two-generator \(R(H_2)\).
