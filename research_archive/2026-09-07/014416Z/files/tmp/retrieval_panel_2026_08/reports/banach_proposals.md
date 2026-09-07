# Frozen Banach-space / tensor-norm proposals

## Verdict

The native object is a unimodular, tetrahedral element of the second
symmetric injective tensor power of \(\ell _1^n\), or equivalently a column
discrepancy problem for the degree-two Walsh evaluation matrix.  The cleanest
strict reduction I can isolate is a **uniform principal-restriction
self-similarity lemma**.  It would force convergence without determining any
finite optimum or any coset-weight histogram.

No theorem I found proves that lemma.  The closest random-submatrix theorem
leaves an additive term of order \(m^{3/2}\), with a fixed constant, precisely
where the lemma needs \(o(m^{3/2})\).  Grothendieck, vector-valued
factorization, and full tensor-power regularization all change the problem or
retain a fixed leading loss.  Thus the first architecture is a genuine
conditional reduction, but I do **not** regard any of the three architectures
below as presently viable without a new theorem of comparable depth to the
original problem.

The proposals below were frozen before any archive exposure.

## 1. Native translation and normalization audit

Let \(\mathscr A_n\) be the set of symmetric, hollow matrices whose
off-diagonal entries are in \(\{\pm1\}\).  For \(A\in\mathscr A_n\), write

\[
 P_A(z)=\sum_{i<j}a_{ij}z_i z_j=\frac12z^{\mathsf T}Az,
 \qquad z\in[-1,1]^n.
\]

This polynomial is tetrahedral (multiaffine), so its norm on the cube is
attained at a vertex:

\[
 \|P_A\|_{C([-1,1]^n)}
 =\max_{x\in\{\pm1\}^n}|P_A(x)|=Q(A).
\]

Put \(e_i\vee e_j=(e_i\otimes e_j+e_j\otimes e_i)/2\) and

\[
 u_A=\sum_{i<j}a_{ij}e_i\vee e_j
       \in \otimes^{2,s}\ell_1^n.
\]

With the symmetric injective norm

\[
 \varepsilon_s(u)=\sup_{\varphi\in B_{\ell_\infty^n}}
 |\langle u,\varphi^{\otimes2}\rangle|,
\]

one has exactly

\[
 \boxed{Q(A)=\varepsilon_s(u_A)},\qquad
 \boxed{M_n=\min_{A\in\mathscr A_n}\varepsilon_s(u_A)}.
\]

Thus this is the minimum symmetric injective norm among the tetrahedral
sign tensors with *every* off-diagonal coefficient present.  It is not the
norm of a generic tensor, a Sidon constant, or a coefficient inequality.

### Walsh-discrepancy form

Let \(E_n=\binom{[n]}2\), let
\(\Omega_n=\{\pm1\}^n/\{x\sim-x\}\), and define the
\(2^{n-1}\times\binom n2\) Walsh evaluation matrix

\[
 W_n([x],\{i,j\})=x_ix_j.
\]

Then

\[
 M_n=\min_{a\in\{\pm1\}^{E_n}}\|W_na\|_\infty
     =\operatorname{disc}(W_n).
\]

This is ordinary column discrepancy of this one matrix.  It is not
\(\operatorname{herdisc}(W_n)\): arbitrary subsets of edge-columns are not
the principal vertex restrictions used below.

For the antipodal cut code \(\mathcal C_n^+\),

\[
 \langle a,\sigma(x_ix_j)_{i<j}\rangle
   =E-2d_H\bigl(a,\sigma(x_ix_j)_{i<j}\bigr),
 \qquad E=\binom n2.
\]

Since \(\mathcal C_n^+\) is antipodal,

\[
 Q(a)=E-2d_H(a,\mathcal C_n^+),\qquad
 M_n=E-2\rho(\mathcal C_n^+).
\]

### Every factor from the diagonal polynomial to bilinear norms

The symmetric bilinear form associated with \(P_A\) is

\[
 \check P_A(s,t)
 =\frac12\sum_{i<j}a_{ij}(s_it_j+s_jt_i)
 =\frac12s^{\mathsf T}At.
\]

Consequently

\[
 \|\check P_A\|=\frac12\|A\|_{\infty\to1}.
\]

Diagonal evaluation and real degree-two polarization give, respectively,

\[
 \|P_A\|\le \|\check P_A\|,
 \qquad
 \|\check P_A\|\le2\|P_A\|.
\]

In the normalization of the problem this is

\[
 \boxed{Q(A)\le\frac12\|A\|_{\infty\to1}\le2Q(A)},
 \qquad
 \boxed{2Q(A)\le\|A\|_{\infty\to1}\le4Q(A)}. \tag{1}
\]

The second inequality can also be seen directly.  For signs \(s,t\), put
\(u=(s+t)/2\) and \(v=(s-t)/2\); their supports are disjoint and
\(s^{\mathsf T}At=u^{\mathsf T}Au-v^{\mathsf T}Av\).  Each partial-sign
quadratic is the conditional expectation of a full-sign quadratic, so each
has absolute value at most \(2Q(A)\).

If

\[
 \operatorname{SDP}(A)=
 \sup_{\|u_i\|=\|v_j\|=1}
 \left|\sum_{i,j}a_{ij}\langle u_i,v_j\rangle\right|,
\]

then real Grothendieck gives

\[
 2Q(A)\le\|A\|_{\infty\to1}
 \le\operatorname{SDP}(A)
 \le K_G\|A\|_{\infty\to1}
 \le4K_GQ(A). \tag{2}
\]

Equivalently, the vector relaxation of \(\check P_A\) lies between \(Q(A)\)
and \(2K_GQ(A)\).  Neither the factor two from real polarization nor the two
independent sign families in Grothendieck can be silently removed.

For comparison with the sharp real bilinear Littlewood inequality, the
ordered coefficients of \(\check P_A\) are \(a_{ij}/2\) for \(i\ne j\).
Hence

\[
 \left(\sum_{i,j}|\check P_A(e_i,e_j)|^{4/3}\right)^{3/4}
   =\frac12[n(n-1)]^{3/4}.
\]

The sharp inequality with constant \(\sqrt2\), followed by polarization,
only yields

\[
 \frac12[n(n-1)]^{3/4}
 \le\sqrt2\,\|\check P_A\|
 \le2\sqrt2\,Q(A),
\]

or \(Q(A)\ge[n(n-1)]^{3/4}/(4\sqrt2)\).  The sharp bilinear constant is
therefore not a sharp constant for the same-spin problem.  In the direct
Boolean polynomial inequality the \(\binom n2\) coefficients have
\(\ell_{4/3}\)-norm \(\binom n2^{3/4}\), with no ordered-pair factor, but its
Boolean polynomial constant is a different constant.

## 2. Three frozen architectures

### Architecture I (best): hereditary principal-restriction self-similarity

For \(S\subset[n]\), the principal restriction satisfies

\[
 u_{A[S]}=(\pi_S\otimes_s\pi_S)u_A,
\]

where \(\pi_S:\ell_1^n\to\ell_1^S\) is a norm-one coordinate projection.
Injectivity gives only \(Q(A[S])\le Q(A)\).  The proposed mechanism is that
among all restrictions of a low-norm complete sign tensor, one restriction
at every smaller order has the sharper, scale-correct contraction

\[
 Q(A[S])\lesssim (|S|/n)^{3/2}Q(A).
\]

The exact uniform version is Lemma PRS in Section 4.  It is the best proposal
because it operates in the same symmetric injective norm, preserves
complete support and one shared spin vector, and its conclusion is enough by
itself to compare every large order with an arbitrary low subsequence.

### Architecture II: block gluing with genuinely joint cancellation

For \(A\in\mathscr A_n\), \(B\in\mathscr A_m\), and
\(D\in\{\pm1\}^{n\times m}\), define

\[
 C_D=\begin{pmatrix}A&D\\D^{\mathsf T}&B\end{pmatrix}.
\]

Its same-spin polynomial is exactly

\[
 P_{C_D}(x,y)=P_A(x)+P_B(y)+x^{\mathsf T}Dy. \tag{3}
\]

A sufficient gluing statement would be: for some \(c,\delta>0\), for every
\(n,m\) and optimal \(A,B\), some \(D\) obeys

\[
 Q(C_D)^{2/3}
 \le Q(A)^{2/3}+Q(B)^{2/3}
       +c(n+m)^{1-\delta}. \tag{4}
\]

Then \(a_n=M_n^{2/3}=O(n)\) is almost subadditive with a summable dyadic
error, and the standard approximate-Fekete argument gives convergence of
\(a_n/n\), hence of \(M_n/n^{3/2}\).

This is mathematically distinct from PRS: it composes two prescribed
tensors and chooses only the bipartite coefficient block.  However, the
ordinary rectangular KSZ estimate controls the cross term in (3) only by
\(O(\sqrt{nm(n+m)})\), which is a leading \((n+m)^{3/2}\) channel when
\(n\asymp m\).  Triangle inequality therefore gives an order-\(n+m\), not
an \(o(n+m)\), error after taking the \(2/3\) power.  Vector-valued
Grothendieck factorization does not choose unimodular \(D\) so as to cancel
the two offsets simultaneously for every \((x,y)\).  This route is frozen
but rejected as substantially more circular than PRS.

### Architecture III: diagonal completion and full tensor-power regularization

Choose a diagonal sign matrix \(D_0\) and put \(T=A+D_0\).  Then \(T\) is a
full symmetric sign matrix, so \(T^{\otimes k}\) is again full and symmetric
of order \(N=n^k\).  Hollowing it changes the same-spin polynomial by the
constant

\[
 \frac12\sum_{\alpha=1}^N(T^{\otimes k})_{\alpha\alpha},
\]

whose absolute value is at most \(N/2=o(N^{3/2})\).  Thus diagonal zeros are
not the obstruction to this particular tensorization.

Let \(\phi_T:\ell_\infty^n\to\ell_1^n\) have matrix \(T\).  The
Aubrun--Müller-Hermes tensor-power theorem gives

\[
 \lim_{k\to\infty}
 \|\phi_T^{\otimes k}:
   (\ell_\infty^n)^{\otimes_\varepsilon k}
   \to(\ell_1^n)^{\otimes_\pi k}\|^{1/k}
 =\gamma_2^*(\phi_T). \tag{5}
\]

Here the tensor spaces identify with \(\ell_\infty^{n^k}\) and
\(\ell_1^{n^k}\), so (5) is genuinely an independent-block norm of the
Kronecker powers.  It still does not yield the desired result: it gives only
a \(k\)-th-root exponential rate for one constructed subsequence; it gives
neither a sharp multiplicative prefactor nor a lower bound for the minimum;
and (1) leaves a fixed same-spin/bilinear interval.  An
\(\exp(o(k))=N^{o(1)}\) uncertainty is far larger than the
\(1+o(1)\) precision needed for a normalized constant.  Transferring from
orders \(n^k\) to every order is another missing theorem.  Any repair that
asserts the sharp prefactor uniformly over all sign tensors simply restates
the required optimization in tensor language.  I therefore reject this
architecture.

## 3. Best exact implication

The complete-support Boolean KSZ theorem of Defant--Galicer--Mansilla--
Mastyło--Muro gives a finite constant \(C_{\rm KSZ}\) with
\(M_n\le C_{\rm KSZ}n^{3/2}\).  In fact, the supplied frontier gives the
stronger \(\limsup M_n/n^{3/2}\le1/2\).  With Lemma PRS below,

\[
 \boxed{
 \begin{gathered}
 \text{complete-support Boolean KSZ (boundedness at scale }n^{3/2}\text{)}\\
 {}+\ \text{Lemma PRS}
 \end{gathered}}
 \quad\Longrightarrow\quad
 \boxed{\displaystyle\lim_{n\to\infty}\frac{M_n}{n^{3/2}}
 \text{ exists}.}
\]

Proof of the implication is short but uses the full uniformity in PRS.  Put
\(b_n=M_n/n^{3/2}\) and \(\ell=\liminf b_n\).  Choose a constant
\(C>C_{\rm KSZ}\), and a subsequence \(n_k\to\infty\) with
\(b_{n_k}\to\ell\).  Let \(A_k\) be an optimizer of order \(n_k\).  For any
fixed sufficiently large \(m\le n_k\), PRS gives an \(m\)-set \(S_k\) with

\[
 \frac{M_m}{m^{3/2}}
 \le\frac{Q(A_k[S_k])}{m^{3/2}}
 \le b_{n_k}+\omega_C(m).
\]

Letting \(k\to\infty\) gives \(b_m\le\ell+\omega_C(m)\).  Now let
\(m\to\infty\).  Since \(\omega_C(m)\to0\),
\(\limsup b_m\le\ell=\liminf b_m\).

## 4. Exact missing lemma PRS

> **Lemma PRS (uniform principal-restriction self-similarity).** For every
> finite \(C>0\) there is a function
> \(\omega_C:\mathbb N\to[0,\infty)\), with
> \(\lim_{m\to\infty}\omega_C(m)=0\), such that for all integers
> \(n\ge m\ge1\) and every \(A\in\mathscr A_n\) satisfying
> \(Q(A)\le Cn^{3/2}\), there is an \(S\subset[n]\), \(|S|=m\), for which
> \[
> \boxed{
> Q(A[S])
> \le m^{3/2}\left(\frac{Q(A)}{n^{3/2}}+\omega_C(m)\right).}
> \tag{PRS}
> \]

Equivalently, define the exact defect

\[
 \Delta_C(m)=
 \sup_{n\ge m}\ 
 \sup_{\substack{A\in\mathscr A_n\\Q(A)\le Cn^{3/2}}}
 \left[
 \min_{\substack{S\subset[n]\\|S|=m}}
 \frac{Q(A[S])}{m^{3/2}}
 -\frac{Q(A)}{n^{3/2}}
 \right]_+ . \tag{6}
\]

PRS is exactly \(\Delta_C(m)\to0\) for each fixed \(C\).  Thus the error is
\(o(m^{3/2})\), uniformly in the ambient order \(n\), the ratio \(m/n\),
and every low-norm complete signing \(A\).  No polynomial rate is needed;
the absence of a fixed positive leading error is essential.

The exponent \(3/2\) is forced.  If \(Q(A)=b\,n^{3/2}\), preserving its
normalized norm under restriction means
\(Q(A[S])\le(b+o(1))m^{3/2}\), i.e. a contraction by
\((m/n)^{3/2}\), not by the edge-density factor \((m/n)^2\) and not merely
by norm-one coordinate contraction.

## 5. Why PRS is strictly less information than full optimization

This is a structural, not a description-length, distinction.

1. **PRS is conditional.** It never constructs the input signing and gives
   no upper bound at any order unless a low-norm signing has already been
   supplied.  It also gives no lower bound and no value of \(M_n\).

2. **The coefficient choices are not reoptimized.** The conclusion keeps
   the coefficients inherited from one supplied \(A\) and selects only a
   vertex set \(S\).  Full minimization chooses all \(\binom m2\) signs
   freely.  PRS does not assert that \(A[S]\) is optimal, only that one
   scalar upper inequality holds.

3. **An existential principal restriction cannot recover the energy or
   coset histogram.** Even the stronger datum \(Q(A)\) does not determine
   that histogram.  Here is an exact six-vertex witness.  Let \(A_1\) be
   negative on every edge except \(\{4,6\},\{5,6\}\), and let \(A_2\) be
   negative on every edge except \(\{3,6\},\{4,5\}\).  Direct enumeration
   of the 64 spins gives \(Q(A_1)=Q(A_2)=11\), but the absolute-energy
   histograms are

   | \(|P_A(x)|\) | 1 | 3 | 5 | 7 | 9 | 11 |
   |---:|---:|---:|---:|---:|---:|---:|
   | \(A_1\) count | 28 | 20 | 10 | 2 | 2 | 2 |
   | \(A_2\) count | 28 | 16 | 14 | 4 | 0 | 2 |

   Moreover, the complete list of principal-minimum scalars appearing in
   PRS is the same for the two matrices:
   \[
   \left(\min_{|S|=m}Q(A_r[S])\right)_{m=2}^6=(1,3,4,6,11),
   \qquad r=1,2.
   \]
   The displayed energy tables are also the corresponding distinct
   antipodal coset-correlation histograms.  Thus even all finite
   principal-minimum norm data used by PRS are non-injective with respect to
   the histogram.  The conclusion of PRS supplies no equations for the
   remaining correlations, so its asserted inequalities do not determine
   the omitted histogram.

4. **Many possible limiting constants remain compatible with PRS.** The
   proof in Section 3 determines only equality of liminf and limsup.  It
   does not select a value inside the supplied interval, classify an
   optimizer, or recover any finite \(M_n\).

Thus PRS is not full Boolean/coset optimization in compressed notation.
It is a hereditary compactness assertion about tensors already known to
have the correct order of norm.

## 6. Decisive falsifier

The exact structural falsifier is a constant \(C<\infty\), a number
\(\delta>0\), integers \(n_k\ge m_k\to\infty\), and
\(A_k\in\mathscr A_{n_k}\) such that \(Q(A_k)\le Cn_k^{3/2}\) and

\[
 \min_{|S|=m_k}\frac{Q(A_k[S])}{m_k^{3/2}}
 \ge \frac{Q(A_k)}{n_k^{3/2}}+\delta
 \qquad\text{for every }k. \tag{7}
\]

Equation (7) is precisely \(\limsup_m\Delta_C(m)>0\), so it decisively
kills PRS rather than merely one proof method.  A finite exact search should
compute (6), not merely sample restrictions or compare average cut norms.

The supplied finite anchors already impose a useful normalization check.  If
\(A\) is optimal at order 14, then \(Q(A)=21\), while every 11-vertex
restriction has norm at least \(M_{11}=17\).  Hence

\[
 \Delta_1(11)\ge
 \frac{17}{11^{3/2}}-\frac{21}{14^{3/2}}
 =0.0650802151\ldots. \tag{8}
\]

Thus any proposed quantitative proof with
\(\omega_1(11)<0.0650802151\ldots\) is already false.  The qualitative
lemma allows such finite defects; what must be ruled out is a positive
defect along unbounded \(m\).

For Architecture II, the corresponding decisive obstruction is a sequence
of optimal (or \(o(n^{3/2})\)-near-optimal) blocks \(A_k,B_k\) for which every
cross signing \(D\) violates (4) by \(\delta(n_k+m_k)\).  For Architecture
III, failure of a sharp \(1+o(1)\) prefactor, even when the tensor-power root
limit exists, is already fatal.

## 7. Closest imported theorem and all unmatched hypotheses

The closest primary theorem to PRS is Rudelson--Vershynin,
[*Sampling from large matrices: an approach through geometric functional
analysis*](https://arxiv.org/abs/math/0503442), specifically their
random-principal-submatrix cut-norm estimate.  If each coordinate is retained
with probability \(p=q/n\), their hollow-matrix specialization is

\[
 \mathbb E\|A[Q]\|_{\mathrm C}
 \le K\left(
 p^2\|A\|_{\mathrm C}
 +p^{3/2}(\|A\|_{\mathrm{Col}}+\|A^{\mathsf T}\|_{\mathrm{Col}})
 \right). \tag{9}
\]

For a hollow sign matrix,
\(\|A\|_{\mathrm{Col}}=n\sqrt{n-1}\), so the second term in (9) is
\(\Theta(q^{3/2})\).  Using

\[
 \tfrac14\|A\|_{\infty\to1}
 \le\|A\|_{\mathrm C}\le\|A\|_{\infty\to1}
\]

and (1), (9) yields at best a same-spin statement of the form

\[
 \mathbb E Q(A[Q])
 \le K'\bigl(p^2Q(A)+q^{3/2}\bigr), \tag{10}
\]

with fixed constants.  Their proof uses Banach-valued decoupling,
symmetrization, and a Rademacher Slepian inequality, and they show their
general estimate is optimal.

Every unmatched requirement is material:

- PRS needs an **exact-size** principal set for every \(m\), whereas (9)
  uses Bernoulli size and expectation.  Conditioning can address size but
  not the following losses.
- PRS needs the coefficient of the inherited normalized norm to be one;
  (9)--(10) contain fixed comparison constants.
- The \(q^{3/2}\) term in (10) must improve from a fixed leading term to
  \(o(q^{3/2})\).  This is the central unmatched hypothesis.
- PRS is uniform for every ambient \(n\ge m\), including \(m/n\to0\).
- (9) controls a bilinear cut norm.  Passing to the diagonal same-spin norm
  incurs the factors audited in (1); those factors cannot be absorbed in an
  \(o(1)\) error.
- The proposed improvement may use symmetry, complete unimodular support,
  hollowness, and the extra hypothesis \(Q(A)=O(n^{3/2})\).  The imported
  theorem does not extract any advantage from this joint special class.
- The ordinary injectivity of \(\varepsilon_s\) proves only
  \(Q(A[S])\le Q(A)\).  The Boolean noise operator contracts \(L_\infty\)
  and multiplies a homogeneous level-two polynomial by \(\rho^2\), but it
  averages fractional restrictions; it does not produce an integral
  principal submatrix with retained coefficients in \(\{\pm1\}\).

The direct scale inputs are Defant--Mastyło--Pérez,
[*On the Fourier Spectrum of Functions on Boolean Cubes*](https://arxiv.org/abs/1706.03670),
for Boolean Bohnenblust--Hille lower bounds, and
Defant--Galicer--Mansilla--Mastyło--Muro,
[*Asymptotic Insights for Projection, Gordon--Lewis, and Sidon Constants in
Boolean Cube Function Spaces*](https://arxiv.org/abs/2302.00233), for the
complete-support Boolean KSZ signing.  Neither is a restriction theorem.

For Architecture II, the closest vector-valued result is Defant--Junge,
[*A Vector-Valued Grothendieck Inequality with an Application to
\((p,q)\)-Completely Bounded Operators*](https://doi.org/10.1512/iumj.1999.48.1692).
It assumes independent bilinear factors and produces a norm factorization;
it neither selects a full unimodular cross block nor cancels the two scalar
offset landscapes in (3), and its universal constant would be a leading
loss.

For Architecture III, the exact imported theorem is Aubrun--Müller-Hermes,
[*Limit Formulas for Norms of Tensor Power Operators*](https://arxiv.org/abs/2410.23063).
Its unmatched hypotheses/conclusions are: full rather than symmetric
diagonal tensor powers, a fixed operator rather than minimization over
complete sign tensors, independent input/output blocks, only a root limit,
one multiplicative order subsequence, and no sharp prefactor or all-order
transfer.

## 8. Circularity audit

- A proof of PRS may not choose \(S\) by first inserting an independently
  known optimal order-\(m\) pattern.  That would use the optimization it is
  meant to bypass.
- The error modulus may not be defined or bounded through
  \(\sup_{r\ge m}(b_m-b_r)_+\).  Formula (6) is stated solely in terms of a
  supplied tensor and its actual principal restrictions; replacing it by a
  recurrence for \(M_n\) would erase the structural reduction.
- Convergence alone does not imply PRS.  Assuming that optimal signings can
  be chosen as a nested or projectively consistent family is essentially
  assuming the missing lemma.
- Random-submatrix expectation cannot be turned into PRS by dropping the
  column-norm term in (9).  That term is exactly on the target scale and is
  known to be necessary for general matrices.
- Cut norm, \(\ell_\infty\to\ell_1\), Hilbert-vector value, and same-spin
  norm cannot be identified.  Equations (1)--(2) account for every fixed
  factor.
- A proof that conditions simultaneously on all spin energies risks using
  the full coset-weight histogram.  It must explain why a smaller statistic
  controls the selection of \(S\).
- Architecture II cannot invoke Fekete before proving the sublinear error in
  (4).  Separate bounds on the two internal blocks and the cross block pay a
  fixed leading channel and do not constitute joint cancellation.
- Architecture III cannot replace an \(N^{o(1)}\) tensor-power estimate by
  \(1+o(1)\), or pass from bilinear to diagonal evaluation, without a new
  theorem.  Taking \(k\)-th roots hides exactly the prefactor needed here.
- A generic Boolean BH inequality supplies a lower bound for every
  coefficient array; a generic KSZ argument supplies one upper signing.
  Neither theorem compares different orders and neither can be cited as the
  missing compactness law.

## 9. Confidence and specialist recommendation

| Item | Truth | Tractability with current Banach tools | Conditional implication |
|---|---:|---:|---:|
| PRS (best) | 0.30 | 0.08 | 0.995 |
| Joint-canceling block gluing (4) | 0.15 | 0.03 | 0.98 |
| Tensor-power route as a route to this limit | 0.05 | 0.02 | 0.25 |

The high confidence in the PRS implication is confidence in the deduction,
not in the lemma.  Its main risk is the uniformity over all ambient ratios
and all low-norm signings: multiscale low-norm tensors could have every
intermediate principal restriction separated by a fixed normalized gap.
The finite anchors show that nontrivial gaps already occur.

**Recommendation:** retain PRS as a clean, falsifiable `B`-type conditional
architecture, but do not launch a broad Grothendieck/tensorization program.
The first serious test should be exact or certified lower bounds on
\(\Delta_C(m)\) for structured low-norm families, aimed specifically at the
structural falsifier (7).  Unless those defects decrease, this domain has no
viable architecture.

## Candidate card: PRS

| Field | Entry |
|---|---|
| Domain | Symmetric injective tensor norms; Boolean polynomial discrepancy |
| Imported theorem(s) | Boolean complete-support KSZ ([Defant et al.](https://arxiv.org/abs/2302.00233)); random cut-norm sampling ([Rudelson--Vershynin](https://arxiv.org/abs/math/0503442)) as the closest but insufficient theorem |
| Problem translation | \(M_n=\min\{\varepsilon_s(\sum_{i<j}a_{ij}e_i\vee e_j):a_{ij}=\pm1\}=\operatorname{disc}(W_n)\), with \(P_A(x)=\frac12x^{\mathsf T}Ax\) |
| Proposed mechanism | A low-norm complete sign tensor contains, at every smaller order, one principal restriction with asymptotically no increase in normalized symmetric injective norm |
| Exact missing lemma | Lemma PRS, equivalently \(\Delta_C(m)\to0\) in (6), uniformly in \(n\ge m\) and all \(Q(A)\le Cn^{3/2}\) |
| Why strictly weaker | Conditional on an input signing; selects only one inherited principal submatrix; gives no finite optimum, optimizer, or energy/coset histogram; the explicit order-six pair proves the max norm does not determine the histogram |
| Archive collisions | Not assessed by design; proposals were frozen without archive exposure |
| Falsification test | Produce (7), equivalently \(\limsup_m\Delta_C(m)>0\); (8) is the finite normalization anchor |
| Specialist confidence | Truth 0.30 / tractability 0.08 / implication 0.995 |
| Verifier confidence | To be assigned by archive verifier |
| Director judgment | Hold for targeted falsification; reject broad execution absent decreasing restriction defects |
