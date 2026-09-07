# Dense graph limits, action convergence, exchangeability, and \(\Gamma\)-convergence

Independent domain report, 2026-08-15.  The proposals below were frozen before
any archive exposure.  I used only the naked packet, the supplied graph-limit
toolkit, and the cited primary literature.

## Executive verdict

Ordinary graphon compactness cannot see this problem at its decisive scale.
The scale-correct object is the action of \(A/\sqrt n\), not the bounded graphon
\(W_A\).  The best conditional route is a full-index action-recovery theorem for
spectrally regular finite sign matrices.  It leads to convergence without
identifying the limiting operator or the limiting constant.

There is, however, a real missing theorem rather than a routine application of
compactness: one must both remove concentrated spectral spikes from
near-minimizers at \(o(n^{3/2})\) cost and realize every resulting action limit by
exact symmetric hollow sign matrices at **every sufficiently large order**.
Backhausz--Szegedy give neither assertion.  My truth/tractability estimates for
the combined missing lemma are \(0.35/0.12\).  Conditional on it, the implication
to convergence has confidence \(0.99\).

## 1. Native translation and normalization

Let

\[
\mathfrak A_n=\{A=A^{\mathsf T}:a_{ii}=0,
                 \ a_{ij}\in\{-1,1\}\ (i\ne j)\}.
\]

Put the uniform probability measure on \([n]\) and let the \(P\)-operator

\[
T_A:=\frac{A}{\sqrt n},\qquad
(T_Af)(i)=\frac1{\sqrt n}\sum_{j=1}^n a_{ij}f(j).
\]

For a \(P\)-operator \(T\), define

\[
\Phi(T):=\sup_{\|f\|_\infty\le1}
          \left|\mathbb E\big[f(Tf)\big]\right|.
\]

The normalization is exact:

\[
\mathbb E[f(T_Af)]
  =\frac{f^{\mathsf T}Af}{n^{3/2}},
\qquad
\Phi(T_A)=\frac{2Q(A)}{n^{3/2}},
\qquad
\frac{M_n}{n^{3/2}}
  =\frac12\min_{A\in\mathfrak A_n}\Phi(T_A).
\tag{1}
\]

The second equality is not an assumption about continuity.  Since (A) is
hollow, (f^{\mathsf T}Af) is affine in each coordinate separately.  The
absolute value of an affine function on an interval is maximized at an
endpoint.  Rounding coordinates one at a time therefore proves

\[
\sup_{f\in[-1,1]^n}|f^{\mathsf T}Af|
=\max_{x\in\{\pm1\}^n}|x^{\mathsf T}Ax|=2Q(A).
\]

This also checks the code normalization.  If (a=(a_{ij})_{i<j}) and
(c=(\sigma x_ix_j)_{i<j}\in\mathcal C_n^+), then

\[
\langle a,c\rangle=\binom n2-2d(a,c),
\quad
Q(A)=\binom n2-2d(a,\mathcal C_n^+).
\]

Near-minimizers are bounded in the bare action compactness norm.  Indeed, if

\[
R(A)=\max_x|x^{\mathsf T}Ax|=2Q(A),\qquad
B(A)=\max_{s,t\in\{\pm1\}^n}|s^{\mathsf T}At|,
\]

then polarization gives

\[
s^{\mathsf T}At
=\left(\frac{s+t}{2}\right)^{\!\mathsf T}
 A\left(\frac{s+t}{2}\right)
-\left(\frac{s-t}{2}\right)^{\!\mathsf T}
 A\left(\frac{s-t}{2}\right),
\]

so (B(A)\le2R(A)=4Q(A)).  Consequently

\[
\|T_A\|_{\infty\to1}=\frac{B(A)}{n^{3/2}}
\le\frac{4Q(A)}{n^{3/2}}.
\tag{2}
\]

This gives subsequential action compactness for any near-minimizing sequence,
but not continuity of \(\Phi\).

### Why ordinary graphons collapse

For the bounded signed step graphon (W_A),

\[
\|W_A\|_\square
\le \frac{B(A)}{n^2}
\le \frac{4Q(A)}{n^2}.
\tag{3}
\]

Thus every sequence with (Q(A)=O(n^{3/2})) converges in cut distance to the
zero signed graphon.  Bounded-graphon ground-state theorems correctly return
zero after (n^2) normalization and cannot distinguish the coefficient in
(1).  On the two-point weight space, every decoration is affine in the sign;
hence the raw probability-graphons likewise converge to the constant edge law
\(\tfrac12(\delta_{-1}+\delta_1)\).  Scaling to the action kernel
\(\sqrt n W_A\) places all edge mass at \(\pm\sqrt n\), so the weight laws are
not tight.  Finally, the Braides--Cermelli--Dovetta Gamma-limit at the zero
graphon is the zero functional.  Its full-index chattering recovery does not
recover a subleading (n^{3/2}) coefficient.

## 2. Frozen architectures

### Architecture A (best): spectral action compactness plus full-index recovery

Choose near-minimizers whose rescaled matrices have a common (2\to2) bound.
Backhausz--Szegedy compactness then gives an action limit (T).  Prove that
every such finite-sign action limit has exact symmetric hollow sign realizers
of every sufficiently large order, still with a uniform (2\to2) bound.
Continuity of \(\Phi\), proved in Section 4 below, transfers the liminf
coefficient to all large orders.

In Gamma-convergence language, compactness supplies the liminf cluster point
and the all-order realization is the recovery sequence.  A complete Hausdorff
limit of all feasible sets would be stronger than needed; the boxed lemma below
asks only for the spectral regularization and recovery facts actually used.

### Architecture B: exchangeable induced-restriction transfer

Randomly relabel a large low-energy sign matrix and restrict to \(n\) sampled
vertices.  The exact missing statement would be the following hereditary
fluctuation transfer.  For an order-\(N\) matrix put
\(q_N(A)=Q(A)/N^{3/2}\).  For every \(R<\infty\) and
\(\varepsilon>0\),

\[
\begin{split}
&\exists n_0\ \forall n\ge n_0\ \exists N_0(n,R,\varepsilon)\
 \forall N\ge N_0\ \forall A\in\mathfrak A_N:\\
&\hspace{25mm}q_N(A)\le R\Longrightarrow
 \exists S\subset[N],\ |S|=n,\quad
 \frac{Q(A[S])}{n^{3/2}}\le q_N(A)+\varepsilon.
\end{split}
\tag{IR}
\]

If (IR) holds, take arbitrarily large optimal orders along a subsequence
attaining \(\ell=\liminf M_N/N^{3/2}\).  For each sufficiently large fixed
order (n), choose a still larger host from that subsequence and apply (IR).
Then (M_n/n^{3/2}\le\ell+2\varepsilon), so limsup equals liminf.

This is genuinely projective/exchangeable rather than an action-realization
argument.  Diaconis--Janson only control fixed dense restrictions and their
graphon limits.  They do not give (IR), where (n\to\infty), the normalization
changes from (N^{-3/2}) to (n^{-3/2}), and a supremum over (2^n) spins must
be transferred with no fixed loss.  A structural falsifier is a double sequence
\(n_j\ll N_j\), \(A_j\in\mathfrak A_{N_j}\), and \(\delta>0\) such that every
\(n_j\)-vertex principal submatrix satisfies
\[
\frac{Q(A_j[S])}{n_j^{3/2}}
\ge \frac{Q(A_j)}{N_j^{3/2}}+\delta.
\]
My truth/tractability estimates for (IR) are
\(0.25/0.08\).

There is no third viable architecture from ordinary graphons or bounded
Young-measure Gamma-convergence: equation (3) makes their common limit exactly
zero.

## 3. Best implication

Let \(\mathscr L_C\) be the collection of weak-equivalence classes of all
action limits of sequences

\[
T_{A_j}=A_j/\sqrt{n_j},\qquad n_j\to\infty,qquad
A_j\in\mathfrak A_{n_j},\qquad \|A_j\|_{\rm op}\le C\sqrt{n_j}.
\]

The exact implication is

\[
\boxed{\begin{gathered}
\text{Backhausz--Szegedy Theorems 2.14 and 2.16}\cr
+\ \text{their self-adjoint closure result}\cr
+\ \text{the continuity proposition proved below}\cr
+\ \boxed{L_{\rm FA}:\text{ full-order spectral action recovery}}
\end{gathered}}
\quad\Longrightarrow\quad
\frac{M_n}{n^{3/2}}\ \text{converges}.
\tag{4}
\]

The primary source is Backhausz--Szegedy,
[Action convergence of operators and graphs](https://arxiv.org/abs/1811.00626),
published in the Canadian Journal of Mathematics
([DOI](https://doi.org/10.4153/S0008414X2000070X)).

## 4. Exact missing lemma and proof of the implication

### Boxed lemma \(L_{\rm FA}\): full-order spectral action recovery

There exist \(C<\infty\), an integer \(n_0\), and numbers
\(\epsilon_n\downarrow0\) such that:

1. **Spectral regularization at fluctuation accuracy.**  For every
   (n\ge n_0), there is (A_n\in\mathfrak A_n) satisfying

   \[
   Q(A_n)\le M_n+\epsilon_n n^{3/2},
   \qquad
   \|A_n\|_{\rm op}\le C\sqrt n.
   \tag{SR}
   \]

2. **Full-index exact realization.**  For every (T\in\mathscr L_C), there
   is a finite constant (C_T) such that, for every \(\varepsilon>0\), there
   is (N(T,\varepsilon)) for which

   \[
   \forall m\ge N(T,\varepsilon)\ \exists B_m\in\mathfrak A_m:
   \quad
   \|B_m\|_{\rm op}\le C_T\sqrt m,
   \qquad d_M(T_{B_m},T)\le\varepsilon.
   \tag{AR}
   \]

The quantifier \(\forall m\ge N\) is essential: a realizing subsequence is not
enough.  Each (B_m) has exact order (m), exact symmetry, exact zero
diagonal, and exact off-diagonal signs.  Equivalently, (AR) provides a recovery
modulus (r_T(m)\downarrow0) and realizers with
\(d_M(T_{B_m},T)\le r_T(m)\) for all sufficiently large (m).  Uniformity in
(T) is not needed for (4); uniformity over all large orders for each fixed
(T) is needed.

### Continuity of the Boolean supremum (proved, not assumed)

**Proposition.**  Suppose \(S,T\) are \(P\)-operators with
\(\|S\|_{2\to2},\|T\|_{2\to2}\le C\).  Let

\[
\delta=d_H^{\rm LP}(\overline{\mathcal S_1(S)},
                     \overline{\mathcal S_1(T)}).
\]

For (0<\delta\le1),

\[
|\Phi(S)-\Phi(T)|\le 5C\sqrt\delta+\delta.
\tag{5}
\]

Since \(\delta\le2d_M(S,T)\), this gives the explicit action modulus

\[
|\Phi(S)-\Phi(T)|
\le5C\sqrt{2d_M(S,T)}+2d_M(S,T).
\tag{6}
\]

**Proof.**  Every law \(\mu\) in either closed one-profile is the law of, or a
weak limit of laws of, a pair \((X,Y)=(f,Tf)\) with \(|X|\le1\) and
\(\mathbb E Y^2\le C^2\).  The latter bound passes to weak limits by lower
semicontinuity.  If \(d_{\rm LP}(\mu,\nu)\le\delta\), Strassen coupling gives
coupled \((X,Y)\sim\mu\), \((X',Y')\sim\nu\) for which their Euclidean distance
exceeds \(\delta\) with probability at most \(\delta\) (use
\(\delta+o(1)\) and pass to the limit if necessary).

Let \(\theta_R(y)=\max(-R,\min(y,R))\) and
\(g_R(x,y)=x\theta_R(y)\).  On the good coupling event,
\(|g_R(X,Y)-g_R(X',Y')|\le(R+1)\delta\); on the exceptional event the
difference is at most (2R).  Moreover,

\[
\mathbb E|XY-g_R(X,Y)|
\le\mathbb E[|Y|1_{\{|Y|>R\}}]\le C^2/R,
\]

and likewise for \(\nu\).  Hence

\[
\left|\int xy\,d\mu-\int xy\,d\nu\right|
\le \frac{2C^2}{R}+(3R+1)\delta.
\]

Taking (R=C/\sqrt\delta) gives (5) (the case (C=0) is immediate).
Hausdorff matching in both directions, followed by taking the two suprema of
the absolute integral, proves (5).  Finally the (k=1) term has weight (1/2)
in (d_M), giving (6).  This proof is also why a mere
\(\infty\to1\) bound is insufficient.  \(\square\)

### Deduction of convergence

Let

\[
\ell=\liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\]

and choose (n_j\to\infty) along which the normalized minima tend to
\(\ell\).  Choose (A_{n_j}) from (SR).  Theorem 2.16 of
Backhausz--Szegedy (with (p=q=2)) gives an action-convergent subsequence, and
Theorem 2.14 represents its limit by a \(P\)-operator \(T\) with
\(\|T\|_{2\to2}\le C\).  Their closure proposition preserves
self-adjointness.  By (1), (SR), and (6),

\[
\Phi(T)=\lim_j\Phi(T_{A_{n_j}})=2\ell.
\tag{7}
\]

Apply (AR) to this fixed (T).  It gives exact (B_m\in\mathfrak A_m) for
every sufficiently large (m), with a common (2\to2) bound and
\(T_{B_m}\to T\).  Another application of (6) and (1) yields

\[
\lim_{m\to\infty}\frac{Q(B_m)}{m^{3/2}}=\ell.
\]

Since (M_m\le Q(B_m)), \(\limsup M_m/m^{3/2}\le\ell\).  The reverse
inequality is the definition of \(\ell\), proving convergence.  Notice that no
uniqueness or identification of (T) is used.

## 5. Why \(L_{\rm FA}\) is strictly weaker than full Boolean/coset optimization

This is a stability statement, not a description-length assertion.  Suppose
(A,A'\in\mathfrak A_n) differ on (r) unordered edges.  Then

\[
|H_A(x)-H_{A'}(x)|\le2r\quad\text{for every }x,
\qquad
|Q(A)-Q(A')|\le2r,
\tag{8}
\]

while

\[
\|T_A-T_{A'}\|_{\infty\to1}\le\frac{4r}{n^{3/2}},
\qquad
d_M(T_A,T_{A'})
\le6\sqrt{\frac r{n^{3/2}}}
\tag{9}
\]

by Backhausz--Szegedy Lemma 2.19.  Thus changes on
\(r=o(n^{3/2})\) edges are invisible both to the action-limit assertion and to
the normalized (o(n^{3/2})) comparison in (SR).

Full Boolean/coset data are not invariant under even one such change.  Flipping
one edge changes the labeled energy (H_A(x)) by \(\pm2\) for every Boolean
vector.  In code language it changes every labeled distance
\(d(a,c)\), (c\in\mathcal C_n^+\), by \(\pm1\).  It can change the optimum
itself: at \(n=4\), the all-\(+1\) matrix has \(Q=6\), whereas flipping one edge
gives \(Q=4\).  Hence the information retained in (SR)--(AR) cannot reconstruct
the labeled Boolean table, the labeled coset-distance table, or even exact
finite (Q).  It retains only a stable asymptotic action class and one
one-sided recovery property.  That is a genuine quotient of the full
optimization information.

## 6. Decisive falsifiers and a mandatory obstruction test

### Exact falsifiers for \(L_{\rm FA}\)

The spectral clause is false if one proves that for every (C<\infty) there
are a \(\delta_C>0\) and infinitely many (n) such that

\[
\min_{\substack{A\in\mathfrak A_n\\
                 \|A\|_{\rm op}\le C\sqrt n}}Q(A)
\ge M_n+\delta_C n^{3/2}.
\tag{10}
\]

The all-order clause is false if, for some (C), an action limit
(T\in\mathscr L_C), \(\varepsilon_0>0\), and infinite set of orders (S),
every (m\in S) and every spectrally bounded exact sign matrix of order (m)
remain at action distance at least \(\varepsilon_0\) from (T).  A congruence
class obstruction for a conference/design-like limit would therefore kill the
route immediately.

### Bare action compactness fails decisively

Let (k_n=\lfloor n^{3/4}\rfloor), let (D_n) be the hollow adjacency matrix
of a (k_n)-clique padded by zeros, and put (S_n=D_n/\sqrt n).  For every
fixed profile size, (S_n) and the zero operator have profile laws that can be
coupled to agree off a set of mass (k_n/n\).  In both Hausdorff directions one
sets all test functions to zero on the exceptional clique when necessary.
Consequently

\[
d_M(S_n,0)\le k_n/n\longrightarrow0.
\]

Nevertheless, with (f=1) on the clique and zero elsewhere,

\[
\Phi(S_n)\ge\mathbb E[f(S_nf)]
=\frac{k_n(k_n-1)}{n^{3/2}}\longrightarrow1,
\]

and \(\|S_n\|_{2\to2}\sim n^{1/4}\).  Thus any proposal that invokes only
(2), Lévy--Prokhorov profile convergence, or ordinary action subsequential
compactness is structurally falsified.  Uniform integrability (the spectral
bound in \(L_{\rm FA}\) is one sufficient form) must be paid for explicitly.

## 7. Imported leverage and every missing hypothesis

### Closest theorem for Architecture A

Backhausz--Szegedy prove:

- Lemma 2.11: uniform \(\infty\to1\) boundedness gives an
  action-convergent subsequence.
- Theorem 2.14: an action-convergent sequence uniformly bounded in
  (p\to q), with (p<\infty), has a representing (P)-operator limit with
  the same bound.
- Theorem 2.16: for (p<\infty), (q>1), the weak-equivalence classes with a
  fixed (p\to q) bound form a compact metric space.
- Proposition 3.3: under the relevant moment bound, self-adjointness is closed.

Their random-matrix Proposition 11.1 is only subsequential, uses nonsymmetric
iid matrices as stated, and explicitly does not prove full-sequence
convergence.

The hypotheses/conclusions still missing here are all of the following:

1. a common (2\to2) bound for (o(n^{3/2}))-near-minimizers;
2. a proof that concentrated high-eigenvalue pieces can be removed without a
   fixed loss in the absolute two-channel objective;
3. characterization of which abstract action limits retain finite
   symmetric-hollow-sign realizations;
4. exact realizers, not real-weight or randomized approximants;
5. realizers at every sufficiently large integer order, not merely a
   subsequence or a set of admissible congruence classes;
6. a uniform moment bound along the realizing sequence;
7. compatibility of the (1/\sqrt m) normalization with whatever sampling,
   blow-up, or rounding constructs those realizers;
8. joint preservation of the positive and negative quadratic channels through
   the single absolute functional \(\Phi\), with no factor loss;
9. a proof independent of the desired convergence that the relevant finite
   closure is not order-dependent.

The Boolean-supremum continuity is not left as a missing hypothesis; it is
proved in (5)--(6).

### Closest Gamma theorem

Braides--Cermelli--Dovetta,
[Gamma-limit of the cut functional on dense graph sequences](https://arxiv.org/abs/1806.03436)
([DOI](https://doi.org/10.1051/cocv/2019029)), construct recovery sequences at
every index for a fixed bounded cut-convergent kernel sequence and a fixed
finite label set.  Missing here are boundedness of \(\sqrt nW_A\), nonzero
fluctuation-scale compactness, variation of the disorder (A) itself, and an
all-index exact sign realization.  At the only available ordinary graphon
limit, (W=0), their limit energy is identically zero.

### Closest exchangeability theorem for Architecture B

Diaconis--Janson,
[Graph limits and exchangeable random graphs](https://arxiv.org/abs/0712.2749),
identify dense graph limits with exchangeable infinite graph laws and prove
almost-sure convergence of finite restrictions to the graph limit.  To obtain
(IR), one still needs projective consistency for symmetrized near-minimizers,
a tangent/fluctuation-scale tightness theory, uniform control of a supremum over
(2^n) labels as (n\to\infty), no fixed concentration or entropy loss,
preservation of both signs of the Hamiltonian jointly, and deterministic
existence for every target order.

## 8. Circularity audit

1. **Do not define the limit using a chosen minimizing subsequence and then
   call subsequential approximation “recovery.”**  (AR) is universal for every
   (T\in\mathscr L_C) and quantifies over all sufficiently large orders.
2. **Do not assume a spectral bound from (2).**  An \(\infty\to1\) bound does
   not imply uniform integrability; the padded-clique example proves this.
3. **A proof of (SR) must be structural.**  Acceptable forms are a
   transformation taking any near-minimizer to a spectrally regular competitor,
   or a necessary optimality argument.  “Choose a bounded-norm minimizer” is
   exactly the missing assertion.
4. **A proof of (AR) may not use convergence of (M_n/n^{3/2}).**  It must be
   an action/sampling/rounding theorem for finite sign operators.
5. **Fixed profile size is not growing-spin control.**  Action convergence of
   all fixed profiles plus the proved moment bound happens to control \(\Phi\)
   through the one-profile; ordinary fixed motif convergence does not.
6. **The cube-to-Boolean reduction uses hollowness.**  A diagonal contribution
   must be removed before invoking coordinatewise affinity.
7. **Symmetry and signs must survive recovery exactly.**  An arbitrary
   (P)-operator discretization, a nonsymmetric iid matrix, or a real-weight
   approximation is not a valid (B_m\).
8. **No separate-channel payment is allowed.**  Proving bounds for the maximum
   and minimum with separate approximants can introduce a fixed loss.  The
   functional \(\Phi\) and each (B_m) handle the absolute maximum jointly.
9. **A subsequence of admissible orders is insufficient.**  Conference-matrix
   existence or a design tower on special orders does not establish (AR).
10. **The ordinary zero graphon cannot be promoted to a tangent limit without
    new data.**  Any such promotion is precisely the fluctuation compactness
    theorem being requested, not a consequence of dense graphon convergence.

## 9. Confidence

| Claim | Truth | Tractability |
|---|---:|---:|
| Normalization (1), graphon collapse (3), and continuity (5)--(6) | 0.99 | 0.95 |
| Spectral regularization (SR) alone | 0.55 | 0.25 |
| Universal full-index action realization (AR) alone | 0.45 | 0.18 |
| Combined \(L_{\rm FA}\) | 0.35 | 0.12 |
| \(L_{\rm FA}\) plus imported theorems implies convergence | 0.99 | 0.95 |
| Exchangeable restriction lemma (IR) | 0.25 | 0.08 |

The best exact lemma is therefore **full-order spectral action recovery
\(L_{\rm FA}\), clauses (SR) and (AR)**.  It is mathematically sufficient and
strictly coarser than full Boolean/coset optimization, but it is not presently
available from the cited graph-limit literature.
