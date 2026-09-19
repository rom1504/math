# Contrarian-short initial report (archive-blind)

## Verdict first

The best retrospectively simple definition is the **adversarial finite-temperature
pressure of the antipodal cut code**:

\[
 p_n(\beta):=\min_{a\in\{\pm1\}^{\binom n2}}
 {1\over \beta n}\log\left(
 2^{-n}\sum_{c\in\mathcal C_n^+}
 \exp\!\left\{{\beta\langle a,c\rangle\over\sqrt n}\right\}
 \right),\qquad \beta>0.
\]

The one missing statement is merely that (p_n(\beta)) has a thermodynamic
limit at each fixed temperature, locally uniformly in (\beta).  No uniformity
as (\beta\to\infty) is needed.  The exact soft-max inequality then gives

\[
 {M_n\over n^{3/2}}-{\log2\over\beta}
 \le p_n(\beta)\le {M_n\over n^{3/2}},
\]

so fixed-temperature convergence followed by (\beta\to\infty) forces the
desired convergence.  This is the only candidate below for which the final
step is both coefficient-exact and one line long.

My confidence that the pressure-limit lemma is true is **0.36**; confidence
that it admits a genuinely short proof is **0.28**.  The conditional
implication itself has confidence **0.995**.  The main danger is that the
minimum over deterministic sign environments destroys the covariance
interpolation that proves the analogous quenched spin-glass theorem.

I freeze exactly three architectures.  They are ordered by preference and
are not intended as a menu of nearby variants.

## 1. Native translation and normalization audit

Write (E=\binom n2), regard (a=(a_{ij})_{i<j}\) as an element of
\(\{\pm1\}^E\), and put

\[
 q(a):={Q(a)\over n^{3/2}},\qquad m_n:={M_n\over n^{3/2}}.
\]

For (c=(\sigma x_ix_j)_{i<j}\in\mathcal C_n^+\),

\[
 \langle a,c\rangle
 =\sigma\sum_{i<j}a_{ij}x_ix_j.
\]

Consequently

\[
 Q(a)=\max_{c\in\mathcal C_n^+}\langle a,c\rangle
     =E-2d(a,\mathcal C_n^+),
\]

where (d(a,\mathcal C_n^+)\) is the minimum Hamming distance to the code.
Minimizing over (a) therefore gives

\[
 M_n=E-2\max_a d(a,\mathcal C_n^+)=E-2\rho(\mathcal C_n^+).
\]

The antipodal coordinate (\sigma) is essential: it turns the absolute
quadratic maximum into a single one-sided support function.  The code has
(2^n) distinct words at the orders under discussion, which is why the
entropy penalty in Candidate 1 is exactly (\log2/\beta), not
(2\log2/\beta).

For the matrix convention,

\[
 x^TAx=2H_A(x).
\]

Thus a conference matrix has

\[
 Q(A)={1\over2}\max_x|x^TAx|
 \le {1\over2}\|A\|_{\mathrm{op}}\|x\|_2^2
 ={n\over2}\sqrt{n-1}.
\]

This checks both factors of two and the (n^{3/2}) scale.

---

## 2. Architecture I: adversarial finite-temperature pressure

### Correct definition

For (a\in\{\pm1\}^E\), define its normalized antipodal pressure

\[
 P_{n,a}(\beta):={1\over\beta n}\log\left(
 2^{-n}\sum_{c\in\mathcal C_n^+}
 e^{\beta\langle a,c\rangle/\sqrt n}
 \right),
\qquad p_n(\beta):=\min_aP_{n,a}(\beta).
\tag{2.1}
\]

This is not the usual SK quenched pressure.  The disorder is chosen
adversarially *outside* the logarithm, and it is constrained to be a hollow
symmetric sign matrix.

### Exact missing lemma

\[
\boxed{
\begin{minipage}{0.89\linewidth}
**(FT) Fixed-temperature adversarial thermodynamic limit.**  There is a
function (p:(0,\infty)\to\mathbb R\) such that, for every (B\ge1),
\[
 \sup_{\beta\in[B^{-1},B]}|p_n(\beta)-p(\beta)|=o_B(1)
 \qquad(n\to\infty).
\]
Equivalently, for every (B\ge1\) there is a deterministic modulus
\(\omega_B(N)\downarrow0\) such that for every (N\), every (n\ge N\), and
every (\beta\in[B^{-1},B]),
\[
 |p_n(\beta)-p(\beta)|\le\omega_B(N).
\]
No estimate uniform in unbounded (\beta) is asserted.
\end{minipage}}
\tag{FT}
\]

Local uniformity is a natural interpolation target, but pointwise convergence
for every fixed (\beta>0\) would already suffice for the implication below.

### Exact short implication

For every (a), the largest summand in (2.1) is
(\exp(\beta Q(a)/\sqrt n)\).  An average of (2^n) nonnegative terms is at
most its largest term and at least (2^{-n}) times its largest term.  Hence

\[
 q(a)-{\log2\over\beta}\le P_{n,a}(\beta)\le q(a).
\]

Taking the minimum over (a) preserves both inequalities:

\[
 m_n-{\log2\over\beta}\le p_n(\beta)\le m_n.
\tag{2.2}
\]

Let (u=\limsup_nm_n\) and (\ell=\liminf_nm_n\).  Under (FT), for each
fixed (\beta\), (2.2) gives

\[
 u\le p(\beta)+{\log2\over\beta}
 \le \ell+{\log2\over\beta}.
\]

Letting (\beta\to\infty) yields (u\le\ell\), hence convergence.
Therefore

\[
 \boxed{\text{fixed-temperature interpolation proving (FT)}
 \Longrightarrow m_n\text{ converges}.}
\]

The order of limits is the point: first (n\to\infty) with (\beta) fixed,
then (\beta\to\infty).  There is no zero-temperature uniformity bill.

### Why (FT) contains strictly less information

The full coset-weight data for (a) are the (2^n) distances
\(d(a,c)\), with multiplicity.  At a fixed (\beta), (P_{n,a}(\beta)) is
only one normalized exponential moment of those distances, since

\[
 \langle a,c\rangle=E-2d(a,c).
\]

More importantly, (FT) does not provide even the finite-(n) values of those
moments: it asserts only convergence of their minimum, up to an (o(1))
error after division by (n).  Altering (e^{o(n)}) entries of any coset
histogram, changing subexponential multiplicities, changing every
non-minimizing histogram, or changing optimizer identities can leave all
statements in (FT) unchanged.  Such alterations plainly change the full
histogram and can change exact finite-(n) maxima.  Thus neither a histogram
nor an optimizer is recoverable from (FT).

The limiting pressure curve can encode the eventual leading support edge
after (\beta\to\infty); that is exactly the one scalar the implication is
supposed to recover.  It does not encode exact (M_n), the maximizing spins,
or the other coset weights.  This is an information statement, not a claim
that directly evaluating (p_n(\beta)) is computationally easy.

### Decisive falsifier

One asymptotic obstruction kills (FT): exhibit a fixed rational
\(\beta_0>0\), a (\delta>0), and two increasing sequences (r_j,s_j\) such
that

\[
 |p_{r_j}(\beta_0)-p_{s_j}(\beta_0)|\ge\delta
 \quad\text{for every }j.
\tag{2.3}
\]

Exact enumeration or certified bounds for the two scalar minima at a single
temperature are enough to establish (2.3); no ground-state histogram is
needed.  A structural version would be two order classes with disjoint
finite-temperature variational lower and upper bounds.

### Imported leverage and unverified hypotheses

The closest theorem is Guerra--Toninelli, *The Thermodynamic Limit in Mean
Field Spin Glass Models* (2002),
[primary preprint](https://arxiv.org/abs/cond-mat/0204280).  Their smooth
interpolation makes the quenched free energy subadditive by comparing one
Gaussian system with two independent Gaussian subsystems.

Every hypothesis that is missing here matters:

1. their disorder is Gaussian and supports integration by parts;
2. expectation over disorder is outside the logarithm, whereas here there is
   a minimum over deterministic sign environments;
3. Gaussian block interpolation remains inside the admissible disorder
   class, while convex interpolation of two sign matrices does not remain
   sign-valued;
4. the covariance is an explicit convex function of replica overlap;
5. the minimizing environment may change discontinuously with (n) and
   (\beta);
6. the antipodal (\sigma)-channel must be interpolated jointly, not bounded
   as two separate one-sided problems.

Thus Guerra--Toninelli is a proof template, not an imported solution.

### Circularity audit

- Proving (FT) by first taking (\beta\to\infty) uniformly would simply
  reintroduce the desired ground-state convergence.
- Choosing a (Q)-minimizing matrix and showing only that its pressure
  converges does not control the matrix minimizing pressure.
- Interchanging (\min_a) with an interpolation, expectation, or block
  decomposition is unjustified and can reverse the needed inequality.
- A Gaussian universality estimate for each fixed (a) is useless unless its
  error is uniform over (2^E) environments; a union bound at that scale can
  consume the leading constant.
- Rounding a Gaussian interpolant to signs changes (\Theta(n^2)) couplings.
  A per-edge Lipschitz estimate pays (\Theta(n^{1/2})) after normalization,
  not (o(1)).
- Replacing the antipodal partition function by two separately optimized
  one-sided pressures can pay a fixed leading loss and changes the problem.

### Confidence

- Truth of (FT): **0.36**.
- Tractability by a short interpolation/cavity proof: **0.28**.
- Correctness of (FT) (\Rightarrow) convergence: **0.995**.

---

## 3. Architecture II: the squared row-normalized completion norm

### Correct definition

For a sign matrix (A) of order (n), define

\[
 \lambda(A):=\left({Q(A)\over n}\right)^2,
 \qquad \lambda_n:=\min_A\lambda(A)={M_n^2\over n^2}.
\tag{3.1}
\]

The target sequence is now

\[
 m_n^2={\lambda_n\over n}.
\]

Thus the natural composition law is additive for (\lambda), not for (Q).
This is the unexpected normalization: divide energy by the number of
vertices, then square.

Given (A) of order (n), (B) of order (m), and
(D\in\{\pm1\}^{n\times m}), write

\[
 A\oplus_D B:=
 \begin{pmatrix}A&D\\D^T&B\end{pmatrix}.
\]

### Exact missing lemma

\[
\boxed{
\begin{minipage}{0.89\linewidth}
**(PG) Pythagorean sign completion.**  There exist universal constants
(C<\infty), (\eta>0), and (N_0\) such that for all (n,m\ge N_0\) and
all hollow symmetric sign matrices (A,B) satisfying
\[
 Q(A)\le2n^{3/2},\qquad Q(B)\le2m^{3/2},
\]
there is a sign matrix (D\in\{\pm1\}^{n\times m}) for which
\[
 \lambda(A\oplus_DB)
 \le \lambda(A)+\lambda(B)+C(n+m)^{1-\eta}.
\tag{3.2}
\]
The constants and error are uniform in the size ratio (n/m), in (A,B),
and in both signs of the quadratic energy.
\end{minipage}}
\tag{PG}
\]

The harmless constant 2 only restricts the lemma to the (n^{3/2}) regime;
the supplied upper frontier puts all sufficiently large minimizers there.

### Exact short implication

Apply (PG) to minimizers at orders (n,m).  For all sufficiently large
orders,

\[
 \lambda_{n+m}\le\lambda_n+\lambda_m+C(n+m)^\alpha,
 \qquad \alpha:=1-\eta<1.
\tag{3.3}
\]

The generalized subadditive lemma then gives existence of
\(\lim_n\lambda_n/n\).  For completeness, the special case needed here is
short.  Choose a large block size (k).  Merge (q) copies of (k) along a
balanced binary tree.  At level (j), there are (O(q/2^j)) merges of size
(O(2^jk)), so the total defect over all levels is

\[
 O\!\left(qk^\alpha\sum_{j\ge0}2^{j(\alpha-1)}\right)
 =O_\alpha(qk^\alpha).
\]

Hence

\[
 {\lambda_{qk}\over qk}
 \le {\lambda_k\over k}+O_\alpha(k^{\alpha-1}).
\]

For an arbitrary order (N=qk+r), (0\le r<k), use (q-1) blocks of
size (k) and one terminal block of size (k+r\in[k,2k)).  All leaves are
then above (N_0), and the exceptional terminal cost is fixed once (k) is
fixed.  The same balanced-tree estimate, divided by (N), differs by
\(o_{q\to\infty}(1)).  Choose (k\to\infty) along a subsequence realizing
\(\liminf\lambda_k/k\).  Because
(k^{\alpha-1}\to0), limsup equals liminf.  Nonnegativity then gives

\[
 m_n=\sqrt{\lambda_n/n}\longrightarrow
 \sqrt{\lim_n\lambda_n/n}.
\]

Therefore

\[
 \boxed{\text{Hammersley-type almost subadditivity}+(PG)
 \Longrightarrow m_n\text{ converges}.}
\]

### Why (PG) contains strictly less information

For each pair (A,B), (PG) asks only for one cross-block witness and records
only the three scalar norms (Q(A),Q(B),Q(A\oplus_DB)).  It never asks for a
maximizing spin, any non-extreme energy, a coset distance multiplicity, or the
best (D).  The induced recurrence (3.3) is one-sided.  Abstract sequences
\(\lambda_n=cn\) for different (c\ge0) all satisfy its zero-error form, so
the recurrence itself cannot recover the limit constant, much less the exact
values (M_n).  Full histograms can be changed while retaining the three
suprema in (3.2).

The lemma does use (Q(A)) and (Q(B)) as scalar inputs.  Its strictness is
relative to full minimization/histograms, not to the task of certifying one
matrix norm.  A successful proof must upper-bound the output norm without
enumerating its maximizing spins; otherwise the architecture loses its
claimed practical advantage, even though it remains logically sufficient.

### Decisive falsifier

Define the exact completion defect

\[
 \Delta(A,B):=\min_{D\in\{\pm1\}^{n\times m}}
 \left[
 \left({Q(A\oplus_DB)\over n+m}\right)^2
 -\left({Q(A)\over n}\right)^2
 -\left({Q(B)\over m}\right)^2
 \right].
\tag{3.4}
\]

A decisive obstruction is a (delta>0) and admissible matrices
(A_j,B_j\), with both orders tending to infinity, such that

\[
 \Delta(A_j,B_j)\ge\delta(n_j+m_j)
\quad\text{for all }j.
\tag{3.5}
\]

The linear defect in (3.5) rules out every (o(n+m)) version of (PG), not
just a proposed exponent.  At any fixed pair of orders, (3.4) is a finite
Boolean optimization and can be certified exactly; a scalable dual or
algebraic family of certificates would settle the route.

### Imported leverage and unverified hypotheses

The limit step is the discrete special case of Hammersley's
*Generalization of the Fundamental Theorem on Subadditive Functions* (1962),
[primary DOI](https://doi.org/10.1017/S030500410003646X).  Its role begins
only after (3.3); it supplies no cross block.

The hypotheses still to be proved are precisely the difficult ones:

1. existence of a (\{\pm1\}) cross block for arbitrary near-minimal inputs;
2. a root-sum-square, rather than triangle, combination of the within and
   cross energies;
3. an (o(n+m)) defect after squaring;
4. uniformity for every size ratio and every large order;
5. joint control of positive and negative energy channels;
6. preservation of symmetry and the hollow diagonal.

### Circularity audit

- Choosing (D) after identifying all maximizing spins of every candidate
  completion reconstructs the parent Boolean maximum.
- The triangle bound (Q(A\oplus_DB)\le Q(A)+Q(B)+\|D\|_{\infty\to1})
  does not imply (3.2); its cross term is itself of order
  ((n+m)^{3/2}) and produces a fixed loss.
- Bounding the two diagonal blocks and cross block separately forfeits the
  joint cancellation encoded by the square in (3.2).
- An (O(n+m)) defect is not enough.  It can shift
  \(\lambda_n/n=m_n^2\) by a fixed amount.
- A random cross block has the correct *scale* but not a vanishing leading
  coefficient.  Calling it lower order is a normalization error.
- A lemma only for a preferred infinite subsequence or only for balanced
  (n=m), without an independently proved all-order interpolation, is
  insufficient.

### Confidence

- Truth of (PG): **0.14**.
- Tractability by a short completion argument: **0.13**.
- Correctness of (PG) (\Rightarrow) convergence: **0.99**.

---

## 4. Architecture III: critically biased lifts and all-order realization

### Correct definition

Ordinary blow-up is normalized incorrectly: a constant sign on a
(t\times t) block has sum (t^2), whereas an (n^{3/2})-scale coarse edge
must contribute (t^{3/2}).  The correct dilation therefore has block bias
(t^{-1/2}).

Let (r_t) be the integer nearest (t^{3/2}) that is congruent to (t^2)
modulo 2 (break ties upward).  Then (r_t=t^{3/2}+O(1)), and (r_t) is the
sum of some (t^2) signs.

For a hollow sign matrix (A=(a_{ij})) of order (k), define
(\operatorname{Lift}_t(A)) to be the set of hollow symmetric sign matrices
(B) on a partition

\[
 V=V_1\sqcup\cdots\sqcup V_k,\qquad |V_i|=t,
\]

such that

\[
 \sum_{u\in V_i,v\in V_j}b_{uv}=a_{ij}r_t
 \quad(i<j),
\tag{4.1}
\]

and

\[
 \left|\sum_{\{u,v\}\subset V_i}b_{uv}\right|\le1
 \quad(1\le i\le k).
\tag{4.2}
\]

The parity choices make these fibers nonempty at the level of their linear
constraints.  If (y_u=x_i) for (u\in V_i), then (4.1)--(4.2) give

\[
 H_B(y)=r_tH_A(x)+O(k).
\]

Thus every such lift obeys the normalization check

\[
 Q(B)\ge r_tQ(A)-k.
\tag{4.3}
\]

Lift-tightness would match this forced coarse lower bound from above.

### Exact missing lemma

\[
\boxed{
\begin{minipage}{0.89\linewidth}
**(CBL) Critical-bias lift-tightness.**  For every (\varepsilon>0\) there
is (K(\varepsilon)\) such that, for every (k\ge K(\varepsilon)) and every
hollow symmetric sign matrix (A) of order (k) with
\(Q(A)\le2k^{3/2}\), there is (T(A,\varepsilon)<\infty) such that for every
integer (t\ge T(A,\varepsilon)) there exists
\(B\in\operatorname{Lift}_t(A)\) satisfying
\[
 Q(B)\le r_tQ(A)+\varepsilon(kt)^{3/2}.
\tag{4.4}
\]
The same (B) controls both energy signs.  The error is relative (\varepsilon)
at order (kt); (T) may depend on the fixed seed (A), but not on (t).
\end{minipage}}
\tag{CBL}
\]

The lower cutoff on (k) is essential: a finite anomalously good seed should
not be assumed to propagate below the genuine asymptotic infimum.

### Exact short implication

Let (\ell=\liminf_nm_n\), fix (\varepsilon>0\), and choose
(k\ge K(\varepsilon)) from a liminf subsequence with
(m_k\le\ell+\varepsilon).  Let (A) minimize (Q) at order (k).
The supplied upper frontier ensures the harmless (Q(A)\le2k^{3/2})
condition for large (k).

For (N=kt), (4.4) and (r_t/t^{3/2}\to1) imply

\[
 \limsup_{t\to\infty}m_{kt}\le m_k+\varepsilon.
\tag{4.5}
\]

For an arbitrary large (N), write (N=kt+s) with (0\le s<k).  Extend a
lift on (kt) vertices by (s) vertices with arbitrary incident signs.  The
number of newly added edges is at most (sN\le kN), hence

\[
 M_N\le Q(B)+kN.
\]

Because (k) is fixed before (N\to\infty), (kN/N^{3/2}=o(1)).  Thus
(4.5) holds for all orders, and

\[
 \limsup_Nm_N\le m_k+\varepsilon\le\ell+2\varepsilon.
\]

Let (\varepsilon\downarrow0\).  Therefore

\[
 \boxed{(CBL)+\text{fixed-seed padding}\Longrightarrow m_n\text{ converges}.}
\]

### Why (CBL) contains strictly less information

For each fixed seed, (CBL) is a feasibility assertion in one explicitly
defined affine fiber: it asks for one of the exponentially many matrices
having only the (O(k^2)) prescribed block sums (4.1)--(4.2).  It neither
minimizes over all order-(kt) matrices nor identifies the best lift.  Its
input from the seed is only (A) and the scalar (Q(A)), and its output is
one upper-bound witness.  It contains no multiplicities of energies and no
coset histogram.

Indeed, (CBL) supplies only upper envelopes for (M_{kt}).  Arbitrarily many
different sequences of smaller minima are compatible with the same lift
witnesses, so exact (M_{kt}) cannot be recovered.  Conversely, mere
convergence of (m_n) would not imply (CBL), because an unrestricted good
matrix need not realize the prescribed coarse block biases.  The lemma is
structural rather than a disguised recurrence for the optimum.

### Decisive falsifier

For a fixed seed define the exact excess

\[
 \Lambda_t(A):={1\over(kt)^{3/2}}
 \left(\min_{B\in\operatorname{Lift}_t(A)}Q(B)-r_tQ(A)\right).
\tag{4.6}
\]

The route is dead if there are an (\varepsilon_0>0\), arbitrarily large
orders (k), and matrices (A_k) in the stated (n^{3/2}) regime such that

\[
 \liminf_{t\to\infty}\Lambda_t(A_k)\ge\varepsilon_0.
\tag{4.7}
\]

Equation (4.7) says that microscopic spin patterns impose a fixed leading
energy beyond the forced coarse energy.  For fixed (k,t), the inner minimum
in (4.6) is a finite constrained discrepancy problem.  A family of dual
certificates assigning nonnegative weights to spin constraints and proving
the lower bound in (4.7) would be a decisive structural falsifier.

### Imported leverage and unverified hypotheses

The closest general result is Spencer, *Six Standard Deviations Suffice*,
Transactions of the AMS 289 (1985), 679--706,
[primary DOI](https://doi.org/10.2307/2000258).  Applied without the affine
block constraints, general discrepancy methods naturally see
\(\Theta(N^2)) sign variables and (2^N) spin constraints and hence the
(N^{3/2}) scale.

What they do not provide here is exactly what is needed:

1. exact prescribed block sums of size (t^{3/2});
2. a leading coefficient inherited from (Q(A)), rather than an unspecified
   universal discrepancy constant;
3. cancellation of the forced coarse component with microscopic degrees of
   freedom;
4. a relative error tending to zero for every fixed large seed;
5. simultaneous positive/negative control with no separate-channel loss;
6. every sufficiently large lift factor (t), not a favorable subsequence.

### Circularity audit

- Constant-sign cloning has block sums (t^2) and is at the wrong scale.
  Calling it a lift of this problem is a normalization error.
- Unbiased blocks erase the seed.  The (t^{-1/2}) bias is not cosmetic; it
  is the only bias that makes the coarse energy scale as (t^{3/2}Q(A)).
- Prescribed block sums alone do not control microscopic spin patterns.
  Assuming that they do is precisely assuming (CBL).
- Starting from an unrestricted minimizer at order (kt) and then asserting
  that a partition with all the signed biases exists smuggles in a strong
  large-deviation statement.
- A fixed additional microscopic discrepancy constant cannot be absorbed
  into (\varepsilon); it would destroy (4.4).
- Padding is lower order only because (k) is fixed first and (t\to\infty)
  second.  Reversing these limits can make the (kN) term leading.
- The lift must control the antipodal maximum jointly.  Constructing one lift
  for the positive channel and another for the negative channel is invalid.

### Confidence

- Truth of (CBL): **0.11**.
- Tractability by a short partial-coloring/lift argument: **0.10**.
- Correctness of (CBL) (\Rightarrow) convergence: **0.99**.

---

## 5. Candidate cards

| Field | Architecture I: pressure | Architecture II: squared completion | Architecture III: critical lift |
|---|---|---|---|
| Domain | Mean-field statistical mechanics / antipodal cut code | Discrepancy norm composition / almost subadditivity | Critical-scale graph lifts / constrained discrepancy |
| Imported theorem(s) | Guerra--Toninelli thermodynamic-limit interpolation, [primary](https://arxiv.org/abs/cond-mat/0204280) | Hammersley generalized subadditive theorem, [primary](https://doi.org/10.1017/S030500410003646X) | Spencer's six-deviations theorem, [primary](https://doi.org/10.2307/2000258) |
| Problem translation | (p_n(\beta)) is the minimum normalized log-Laplace transform of (\langle a,c\rangle=E-2d(a,c)) | (m_n^2=\lambda_n/n), where (\lambda(A)=(Q(A)/n)^2) | A coarse edge must have block sum (t^{3/2}), hence bias (t^{-1/2}) |
| Proposed mechanism | Prove a thermodynamic limit at every fixed temperature; take zero temperature only after the size limit | Choose a cross sign block that makes squared row-normalized energies add with sublinear defect | Realize every large near-minimizer as the coarse response of tight lifts at all dilation factors |
| Exact missing lemma | (FT) | (PG) | (CBL) |
| Why strictly weaker in information | Only limiting scalar pressures; no finite values, optimizers, or histograms | One witness and one-sided scalar inequality; it does not select the slope or exact minima | One witness in a prescribed affine fiber; only upper envelopes, no unrestricted minimization |
| Ledger collisions | Not assessed: archive access was prohibited until this report was frozen | Not assessed: archive access was prohibited until this report was frozen | Not assessed: archive access was prohibited until this report was frozen |
| Falsification test | Separated subsequences of (p_n(\beta_0)) at one fixed rational (\beta_0) | A family with linear completion defect (3.5) | A persistent positive lift excess (4.7) |
| Specialist confidence | Truth .36 / tractability .28 / implication .995 | Truth .14 / tractability .13 / implication .99 | Truth .11 / tractability .10 / implication .99 |
| Verifier confidence | Deferred to archive verifier | Deferred to archive verifier | Deferred to archive verifier |
| Director judgment | **Execute first**: exact entropy squeeze; test whether deterministic minimization admits interpolation | Hold unless a root-sum-square completion identity appears | Hold; first test for a vector-spin or microscopic excess obstruction |

## 6. Final ranking and fastest tests

1. **Adversarial pressure** is the strongest candidate.  Its decisive virtue
   is that a fixed-temperature theorem suffices; no rate as
   (\beta\to\infty), optimizer stability, or all-order ground-state
   construction is needed.  The first technical test is whether a two-block
   interpolation has a definite sign *after* minimizing over sign
   environments.  If the min destroys the sign, stop.
2. **Squared completion** is algebraically the cleanest normalization but
   asks for very strong joint cancellation.  Compute or certify the defect
   (3.4) on selected near-minimal pairs; a stable linear positive defect ends
   it.
3. **Critical lift** is the cleanest all-order realization mechanism.  Its
   first test is whether the constrained fibers exhibit a fixed microscopic
   excess (4.6); ordinary random blocks are not evidence because they carry
   their own leading ground-state constant.

No fourth architecture is frozen.  In particular, (L_p) soft maxima are
not listed separately from pressure, and generic graphon compactness is not
listed because ordinary dense convergence loses the entire (n^{3/2})
fluctuation scale.
