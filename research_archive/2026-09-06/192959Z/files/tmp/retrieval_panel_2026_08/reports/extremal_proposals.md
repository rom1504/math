# Extremal/probabilistic-combinatorics proposal report

The three architectures below were frozen before inspecting any primary paper beyond the supplied retrieval toolkit.  No project archive, repository computation, or prior route history was consulted.

## 1. Native translation and normalization

Let \(\mathcal A_n\) be the set of symmetric hollow \(n\times n\) sign matrices.  Put \(E_n=\binom n2\).  For \(A\in\mathcal A_n\),

\[
 H_A(x)=\sum_{i<j}a_{ij}x_ix_j=\frac12x^TAx,
 \qquad Q(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|.
\]

Thus \(Q(A)\), not \(2Q(A)\), is the maximum correlation of the edge-sign vector \(a=(a_{ij})_{i<j}\) with an antipodal cut character.  If

\[
 c_x=(x_ix_j)_{i<j},\qquad \mathcal C_n^+=\{\sigma c_x:\sigma\in\{\pm1\},\ x\in\{\pm1\}^n\},
\]

then

\[
 \langle a,c_x\rangle=H_A(x),
 \qquad \langle a,c\rangle=E_n-2d_H(a,c),
\]

and antipodality gives the exact identities

\[
 Q(A)=E_n-2d_H(a,\mathcal C_n^+),
 \qquad M_n=E_n-2\rho(\mathcal C_n^+).
\]

In signed-graph language, if \(T_A=\sum_{i<j}a_{ij}\) and \(C_A(S)=\sum_{i\in S,j\notin S}a_{ij}\), choose \(x=1\) on \(S\) and \(-1\) off \(S\).  Then

\[
 H_A(x)=T_A-2C_A(S).
\]

Equivalently, switch the edge signing across a vertex cut and measure the absolute total imbalance; \(M_n\) is the least possible largest imbalance in a switching class.  This is not induced-subgraph discrepancy.  In particular, principal restriction does **not** automatically decrease \(Q\): spins outside a principal submatrix cannot be set to zero.  Constant-factor cut/induced-discrepancy bridges therefore cannot be used to preserve a putative leading constant.

The exact block formula is central.  For \(A\in\mathcal A_n\), \(B\in\mathcal A_m\), and \(R\in\{\pm1\}^{n\times m}\), let

\[
 C(A,B;R)=\begin{pmatrix}A&R\\R^T&B\end{pmatrix}.
\]

Then

\[
 H_C(x,y)=H_A(x)+H_B(y)+x^TRy.
\]

Since replacing \(y\) by \(-y\) fixes \(H_B(y)\) and negates \(x^TRy\),

\[
 \boxed{\quad
 Q(C(A,B;R))=
 \max_{x,y}\bigl(|H_A(x)+H_B(y)|+|x^TRy|\bigr).
 \quad} \tag{1}
\]

Thus the internal and interface channels cannot be optimized separately without a fixed leading loss.

## 2. Frozen architectures

### Architecture A (best): compatible interface completion and nearly subadditive linearization

The \(3/2\)-homogeneous parameter becomes linear after the transform

\[
 a(n):=M_n^{2/3}.
\]

Given two low-discrepancy diagonal blocks, choose their rectangular sign interface so that the completed block matrix has, up to a summable sublinear error, the \(\ell_{2/3}\)-sum suggested by homogeneity.  Formula (1) makes this a concrete two-block extremal problem rather than an informal gluing principle.  The exact missing lemma is stated in Section 4.

The interface itself is not a small error.  For every \(R\in\{\pm1\}^{n\times m}\), if \(S_r\) denotes a sum of \(r\) independent Rademacher variables, then

\[
 \begin{aligned}
 \|R\|_{\infty\to1}
 &=\max_{x,y}|x^TRy|\\
 &\ge \mathbb E_x\max_y x^TRy
 =m\,\mathbb E|S_n|,
 \end{aligned} \tag{2}
\]

and symmetrically it is at least \(n\mathbb E|S_m|\).  At \(n=m\), this is \((\sqrt{2/\pi}+o(1))n^{3/2}\).  Hence the proposed gain must come from joint energy/interface cancellation inside (1), not from claiming that the interface is \(o(n^{3/2})\).

### Architecture B: hereditary proportional thinning of every low-discrepancy signing

A direct alternative would be a uniform function \(r(t)\to0\) such that for all \(N\ge n\ge N_0\) and every \(A\in\mathcal A_N\) with \(Q(A)\le0.51N^{3/2}\), some \(n\)-set \(S\) satisfies

\[
 Q(A[S])\le \left(\frac nN\right)^{3/2}Q(A)+r(n)n^{3/2}. \tag{3}
\]

This would force convergence directly.  If \(N_j\) realizes the liminf and \(A_j\) minimizes at order \(N_j\), then for every fixed large \(n\), apply (3) and let \(j\to\infty\):

\[
 \frac{M_n}{n^{3/2}}
 \le \liminf_j\frac{M_{N_j}}{N_j^{3/2}}+r(n).
\]

Then let \(n\to\infty\).  This route is mathematically distinct from block completion, but (3) is very strong: it requires every scale ratio, arbitrary low-discrepancy inputs, and an error independent of the ambient order.  Neither principal restriction nor spectral interlacing proves it.  A sequence of low-discrepancy signings for which every selected scale has a fixed positive normalized hereditary gap would decisively falsify it.

### Architecture C: finite-template random switching plus exact absorption

Decompose \(K_N\) into copies of \(K_k\), place a fixed order-\(k\) signing on each copy, independently switch its vertices, and absorb the final leave.  Exact design theory makes the edge bookkeeping plausible for fixed \(k\), but it does not transmit the local minimax constant.

Indeed, the triangle inequality gives, for \(b=E_N/E_k\) blocks,

\[
 Q(\text{global})\le bM_k\asymp \frac{N^2}{\sqrt k}, \tag{4}
\]

which is larger than \(N^{3/2}\) by \(\sqrt{N/k}\) unless \(k=\Omega(N)\).  Random switching changes the calculation but not in the needed way.  For each fixed global spin vector, a switched block has mean zero and exact variance \(E_k\), so the total variance is \(bE_k=E_N=\Theta(N^2)\).  A Bernstein estimate based only on \(M_k=O(k^{3/2})\), at deviation \(t\asymp N^{3/2}\), has exponent only

\[
 \Theta\!\left(
 \frac{N^3}{N^2+k^{3/2}N^{3/2}}
 \right)
 =\Theta\!\left(
 \frac{N}{1+k^{3/2}/N^{1/2}}
 \right). \tag{5}
\]

It can pay a union bound over \(2^N\) spins using this information only when \(k=O(N^{1/3})\).  In that regime the leading bound is variance-driven and forgets the local optimum \(M_k\).  For larger \(k\), one needs the full local log-moment profile

\[
 \log\left(2^{-k}\sum_x e^{\lambda H_{A_k}(x)}\right),
\]

not merely \(M_k\); at zero temperature this approaches the full energy/coset histogram obligation.

For fixed \(k\), exact absorption can remove a leave, and moving to a nearby admissible order costs only \(O_k(N)=o(N^{3/2})\) edges.  But convergence would require \(k\to\infty\), while the fixed-template design theorem has no uniform growing-\(k\) absorber or all-order \(o(N^{3/2})\) retouch bound.  More importantly, an upper construction from selected templates cannot prove the universal lower comparison needed for convergence.  I therefore do not regard this as a viable standalone architecture.

## 3. Best exact implication

Let \(a(n)=M_n^{2/3}\).  The exact implication is

\[
 \boxed{
 \begin{gathered}
 \text{Füredi--Ruzsa, Theorem 5 on comparable-pair nearly subadditive sequences}\\
 {}+\ \text{Interface Completion Lemma }L_{\mathrm{IC}}
 \end{gathered}
 }
 \quad\Longrightarrow\quad
 \boxed{\displaystyle \lim_{n\to\infty}\frac{M_n}{n^{3/2}}\text{ exists}.} \tag{6}
\]

The primary source for the sequence theorem is Z. Füredi and I. Z. Ruzsa, *Nearly subadditive sequences*, Theorem 5, [arXiv:1810.11723](https://arxiv.org/abs/1810.11723), revision dated 2026-03-22.  Its hypotheses are: a real sequence on every positive integer; a fixed \(\mu>1\) and threshold \(N\); a nonnegative nondecreasing error \(f\) with \(\sum_t f(t)/t^2<\infty\); and

\[
 a(n+m)\le a(n)+a(m)+f(n+m)
\]

for all \(N\le n\le m\le\mu n\).  It then gives convergence of \(a(n)/n\).

## 4. Exact missing lemma \(L_{\mathrm{IC}}\)

Fix \(\kappa=51/100\), and use natural logarithms.  The boxed missing assertion is

\[
\boxed{L_{\mathrm{IC}}}.
\]

> There exist constants \(K<\infty\) and an integer \(N_0\) such that the following holds for every pair of integers
>
> \[
> N_0\le n\le m\le2n.
> \]
>
> For every \(A\in\mathcal A_n\) and \(B\in\mathcal A_m\) satisfying
>
> \[
> Q(A)\le\kappa n^{3/2},\qquad Q(B)\le\kappa m^{3/2},
> \]
>
> there is a single matrix \(R\in\{\pm1\}^{n\times m}\), chosen after \(A,B\) but before \(x,y\), for which
>
> \[
> \left[
> \max_{x\in\{\pm1\}^n,\,y\in\{\pm1\}^m}
> \left(|H_A(x)+H_B(y)|+|x^TRy|\right)
> \right]^{2/3}
> \le
> Q(A)^{2/3}+Q(B)^{2/3}
> +K\frac{n+m}{\log^2(e+n+m)}.
> \tag{L_IC}
> \]
>
> The assertion is deterministic and holds at every sufficiently large comparable pair of orders, with no divisibility or arithmetic restriction.

The error is genuinely below the target scale.  Write \(N=n+m\), \(s=Q(A)^{2/3}+Q(B)^{2/3}=O(N)\), and \(f(N)=KN/\log^2(e+N)=o(N)\).  Then

\[
 (s+f)^{3/2}-s^{3/2}
 \le \frac32(s+f)^{1/2}f
 =O\!\left(\frac{N^{3/2}}{\log^2 N}\right)
 =o(N^{3/2}). \tag{7}
\]

This is stronger than an unspecified \(o(N^{3/2})\).  After the \(2/3\) transform, the error is summable in the exact sense required by the imported theorem:

\[
 \sum_{N\ge1}\frac{f(N)}{N^2}
 \ll\sum_{N\ge2}\frac1{N\log^2 N}<\infty. \tag{8}
\]

### Proof of the implication

The supplied upper frontier implies that every exact minimizer satisfies \(Q(A_n)=M_n\le0.51n^{3/2}\) for all sufficiently large \(n\).  Apply \(L_{\mathrm{IC}}\) to minimizers \(A_n,A_m\).  By (1), its interface constructs \(C\in\mathcal A_{n+m}\), so

\[
 \begin{aligned}
 a(n+m)=M_{n+m}^{2/3}
 &\le Q(C)^{2/3}\\
 &\le a(n)+a(m)+K\frac{n+m}{\log^2(e+n+m)}
 \end{aligned} \tag{9}
\]

for every sufficiently large \(n\le m\le2n\).  The error in (9) is eventually nondecreasing, nonnegative, and satisfies (8).  Füredi--Ruzsa Theorem 5 therefore gives \(a(n)/n\to\gamma\).  The supplied two-sided frontier makes \(0<\gamma<\infty\).  Finally,

\[
 \frac{M_n}{n^{3/2}}=\left(\frac{a(n)}n\right)^{3/2}\longrightarrow\gamma^{3/2}. \tag{10}
\]

No identification of \(\gamma\) is needed.

## 5. Why \(L_{\mathrm{IC}}\) is strictly weaker than full signing/coset optimization

The information used by \(L_{\mathrm{IC}}\) is precisely:

1. two scalar upper envelopes \(Q(A),Q(B)\) for inputs already known to lie in the low-discrepancy region;
2. one interface \(R\) for each such pair; and
3. one joint upper certificate for the maximum in (1).

It does not request an optimal interface, an optimal order-\(n+m\) signing, the maximizing spins, the values of \(Q(D)\) for arbitrary \(D\), or the coset-weight/energy histogram.

There is a proper information separation, not merely a description-length comparison.  For sufficiently large \(r\), the all-positive signing \(J_r-I_r\) has \(Q=E_r>0.51r^{3/2}\).  It can never be an input block of \(L_{\mathrm{IC}}\).  A full all-positive signing at order \(n+m\) can never be one of its completed candidates either, since both of its diagonal blocks are all-positive and fail the premises.  Consequently all distance and histogram data on this entire high-discrepancy family are outside every query and every witness required by the lemma.  Altering those entries in a formal optimization oracle leaves all \(L_{\mathrm{IC}}\) answers unchanged while changing the full coset table.  Thus the projection from full optimization data to the data needed by the lemma is non-injective.

The scalar consequence is also strictly non-identifying.  For every \(c>0\), the abstract sequence \(q_n=cn^{3/2}\) has \(q_n^{2/3}=c^{2/3}n\) and satisfies (9) with zero error.  Hence the same recurrence and the supplied frontier are compatible with every limiting constant in the frontier interval.  They neither determine the constant nor recover any optimizer.  This many-to-one compatibility is a mathematical separation, independent of how many bits are used to state the lemma.

The lemma is nevertheless not computationally cheap: a certificate must control all \(2^{n+m}\) spin pairs.  “Strictly weaker” here means less information than the parent minimization/histogram, not polynomial-time verifiability.

## 6. Decisive falsifier

Define the exact finite interface value

\[
 \Gamma(A,B):=\min_{R\in\{\pm1\}^{n\times m}}
 \max_{x,y}\bigl(|H_A(x)+H_B(y)|+|x^TRy|\bigr), \tag{11}
\]

and the transformed excess

\[
 \Delta(A,B):=\Gamma(A,B)^{2/3}-Q(A)^{2/3}-Q(B)^{2/3}. \tag{12}
\]

For any fixed \(A,B\), (11) is an exact finite min--max problem and can be encoded as SAT/MILP or exhaustively evaluated.  For proposed numerical \(K,N_0\), one low-discrepancy pair with

\[
 \Delta(A,B)>K\frac{n+m}{\log^2(e+n+m)} \tag{13}
\]

falsifies that instance of the lemma.  A decisive structural falsifier, independent of the unknown constants, is a sequence of comparable low-discrepancy pairs \((A_j,B_j)\) such that

\[
 \frac{\log^2(e+n_j+m_j)}{n_j+m_j}\,\Delta(A_j,B_j)\longrightarrow\infty. \tag{14}
\]

A positive linear excess \(\Delta(A_j,B_j)\ge c(n_j+m_j)\) is an especially strong obstruction.  The elementary floor (2) is a fast necessary-condition screen:

\[
 \Gamma(A,B)\ge
 \max\{m\mathbb E|S_n|,\ n\mathbb E|S_m|\}. \tag{15}
\]

If the right side of (15) already exceeds the \(3/2\)-power of the proposed right side of \(L_{\mathrm{IC}}\), no interface can work.

The finite anchors also rule out an all-order zero-error version before any interface search.  With \(a(n)=M_n^{2/3}\),

\[
 a(10)-2a(5)=13^{2/3}-2\cdot4^{2/3}=0.489090614\ldots>0,
\]

and

\[
 a(12)-2a(6)=18^{2/3}-2\cdot5^{2/3}=1.020249979\ldots>0.
\]

Thus exact subadditivity is false as an all-order statement; these small cases do not by themselves exclude eventual exact subadditivity.  A quantitative, summable error is therefore not cosmetic in the proposed all-order lemma.

## 7. Imported leverage and every missing hypothesis

### Exact conclusion engine

The closest theorem that actually completes the implication is Füredi--Ruzsa Theorem 5.  Relative to that theorem:

- **All-integer sequence:** verified by \(a(n)=M_n^{2/3}\).
- **Fixed comparable range:** supplied by \(n\le m\le2n\) in \(L_{\mathrm{IC}}\).
- **Threshold only:** allowed by the theorem.
- **Nonnegative nondecreasing error:** \(KN/\log^2(e+N)\) is eventually increasing; enlarge \(N_0\).
- **Summability:** verified in (8).
- **Finite limiting slope:** follows from the supplied positive finite frontier.
- **Missing hypothesis:** only recurrence (9), exactly the content furnished by \(L_{\mathrm{IC}}\), is unproved.

Mere \(o(N^{3/2})\) error on the original scale is insufficient.  If the original excess is \(e(N)\), its linearized size is \(\Theta(e(N)/\sqrt N)\), and the relevant condition is essentially

\[
 \sum_N\frac{e(N)}{N^{5/2}}<\infty. \tag{16}
\]

For example, \(e(N)=N^{3/2}/\log N\) is \(o(N^{3/2})\) but fails (16).  This is precisely why the rate in \(L_{\mathrm{IC}}\) is explicit.

Astashkin--Lykov, *Random unconditional convergence of Rademacher chaos in \(L_\infty\) and sharp estimates for discrepancy of weighted graphs and hypergraphs*, [arXiv:2412.20107](https://arxiv.org/abs/2412.20107), is the closest direct order-of-magnitude theorem in the packet.  It gives \(\Theta(n^{3/2})\) for a minimum signing followed by an induced-vertex-set maximum.  It is not an exact theorem for the present spin maximum.  Writing \(x_i=2z_i-1\) introduces, besides the induced quadratic term, the total sum and a degree-weighted linear term:

\[
 H_A(2z-\mathbf1)
 =4\sum_{i<j}a_{ij}z_iz_j
 -2\sum_i z_i\sum_{j\ne i}a_{ij}
 +\sum_{i<j}a_{ij}. \tag{17}
\]

Its universal implicit constants, different state space, and absence of a composition recurrence leave every hypothesis of \(L_{\mathrm{IC}}\) open.  It supplies scale validation, not the limiting mechanism.

### Closest mechanism, but not a transferable theorem

Guerra--Toninelli, *The thermodynamic limit in mean field spin glass models*, [arXiv:cond-mat/0204280](https://arxiv.org/abs/cond-mat/0204280), is the closest native mechanism: interpolate between one system and two blocks, use the covariance/overlap identity to obtain sub/superadditivity, and pass to a limit.  Every hypothesis missing here is substantial:

- their disorder is independent Gaussian randomness, not an adversarial deterministic signing;
- their conclusion concerns quenched expectation/free energy (and then a random ground state), not an outer minimum over sign matrices;
- Gaussian interpolation supplies a sign-definite covariance remainder; no analogue is known for choosing one Boolean interface \(R\);
- \(L_{\mathrm{IC}}\) must be uniform over every low-discrepancy \(A,B\), not typical disorder;
- the maximum here is absolute, and (1) shows the extra joint channel exactly;
- the interface entries must be exactly \(\pm1\), the diagonal must remain hollow, and the construction must work at every comparable integer pair;
- an explicit summable zero-temperature error is required.

Carmona--Hu universality does not fill these gaps: it compares selected independent disorder laws and does not commute a random maximum with deterministic minimization over the disorder.

### Imported theorems relevant to rejected architectures

Glock--Kühn--Lo--Osthus, *The existence of designs via iterative absorption*, [arXiv:1611.06827](https://arxiv.org/abs/1611.06827), supplies exact \(F\)-decompositions for fixed \(F\) at sufficiently large admissible orders.  Missing for Architecture C are: growing-template uniformity, an all-order \(o(N^{3/2})\) retouch estimate, simultaneous cancellation for every spin vector, and a universal lower comparison rather than a selected upper construction.

Paley's conference construction ([DOI](https://doi.org/10.1002/sapm1933121311)) supplies selected matrices with spectral constant \(1/2\).  If an order-\(N\) conference matrix is restricted to \(n\) coordinates, then

\[
 Q(B)\le\frac n2\sqrt{N-1}. \tag{18}
\]

This preserves \(1/2\) only when \(N/n\to1\); a fixed ratio produces the fixed factor \(\sqrt{N/n}\).  It supplies no universal statement about minimizers and no lower bound.  It therefore cannot substitute for Architecture B or \(L_{\mathrm{IC}}\).

## 8. Scale and quantifier stress test

### Absorption

An uncovered edge contributes at most one if it is being filled from coefficient zero; flipping an already signed edge changes every Hamiltonian value by at most two.  Thus \(L\) uncovered edges or retouched signs cost at most \(L\) or \(2L\), respectively.  The required condition is \(L=o(N^{3/2})\).  A standard \(o(N^2)\) leave is not enough.  Fixed-template exact absorption passes this bookkeeping test on admissible orders but fails the growing-template and joint-cancellation tests described above.

### Random and random-algebraic construction

Random sign matrices, Paley matrices, and random switching of design blocks are selected families.  They can improve an upper bound, even uniformly over all spins after sampling, but cannot furnish a universal lower bound for arbitrary signings.  Random SK universality identifies a typical-disorder maximum, not \(\min_A Q(A)\).  No exchange of these quantifiers is valid.

### Restriction and extension

For this spin functional, arbitrary principal restriction has no monotonicity theorem.  The conference bound (18) works because operator norm is inherited, and it has a fixed loss unless the ambient order is \((1+o(1))n\).  Extension from an order \(m<n\) signing creates \(m(n-m)+\binom{n-m}{2}\) new edges.  Without cancellation the cost is \(O(n(n-m))\), so preserving a leading constant requires \(n-m=o(\sqrt n)\), not merely \(o(n)\).

### Stability

If two signings differ on \(L\) edges, then \(|Q(A)-Q(B)|\le2L\).  Therefore a stability theorem stated only as \(o(N^2)\)-edge closeness is too coarse; one needs \(o(N^{3/2})\)-edge closeness or directly \(o(N^{3/2})\) switching-norm closeness.  Moreover, stability to Paley/conference examples would still be a selected-structure assertion until proved for every near-minimizer at every large order.

### Near-subadditivity

The correct linearized object is \(M_n^{2/3}\), and the error must satisfy the summability threshold (16).  Bounds, restriction along an arithmetic subsequence, or an \(o(N^{3/2})\) recurrence with no rate do not force a limit.  This is an actual obstruction, not a technical preference: the converse part of Füredi--Ruzsa constructs maximally nonconvergent nearly subadditive slopes when the error series diverges.

## 9. Circularity and quantifier audit

1. **No assumed limit.**  The threshold \(0.51\) uses only the supplied \(\limsup\le1/2\), not convergence or a guessed constant.
2. **Uniform inputs.**  \(L_{\mathrm{IC}}\) applies to every signing below the explicit threshold, hence in particular to arbitrary exact or near minimizers.  A proof for Paley, conference, random, or any other selected family would not prove the lemma.
3. **One interface, all witnesses.**  \(R\) may depend on \(A,B\) but not on \(x,y\).  Choosing a different interface for each spin pair reverses the minimax order and is invalid.
4. **Joint payment.**  Formula (1) handles the absolute value and the internal/interface channels exactly.  Separately bounding \(|H_A|\), \(|H_B|\), and \(|x^TRy|\) pays a fixed main-scale loss.
5. **Every comparable order.**  Arithmetic subsequences or admissible design orders do not meet the hypothesis of the sequence theorem.
6. **Summable error.**  Replacing the displayed rate by an unnamed \(o(N^{3/2})\) is not sufficient.
7. **No parent optimizer hidden in the construction.**  Defining \(R\) as the off-diagonal block of an order-\(n+m\) optimizer would not work: its diagonal blocks need not equal the prescribed \(A,B\).  Requiring an optimizer with all prescribed principal blocks would be at least as strong as the parent problem and would be circular.
8. **No histogram smuggling.**  A proof that tabulates the full energy distribution of every low signing, or all distances in every coset, may establish the statement but loses the advertised information reduction.  The intended proof must use a coarser certificate, such as an entropy bound only for near-extremal spin pairs plus a discrepancy construction for their interface constraints.
9. **Selected constructions are only upper bounds.**  No conference/random/design example is used as a universal lower bound.
10. **Nonconvergence standard.**  None of the architectures currently implies genuine nonconvergence.  Failure of \(L_{\mathrm{IC}}\), unusual finite values, or failure of subadditivity would not show \(\liminf<\limsup\).

## 10. Candidate card and confidence

| Field | Entry |
|---|---|
| Domain | Extremal discrepancy/coding theory with a nearly-subadditive sequence theorem |
| Imported theorem(s) | Füredi--Ruzsa Theorem 5 ([primary](https://arxiv.org/abs/1810.11723)); Guerra--Toninelli interpolation is mechanism-only ([primary](https://arxiv.org/abs/cond-mat/0204280)) |
| Problem translation | Covering-radius defect of the antipodal cut code; equivalently least maximum switching-class total imbalance; \(H=x^TAx/2\) |
| Proposed mechanism | Choose one Boolean rectangular interface jointly adapted to two arbitrary low-discrepancy diagonal blocks; linearize by the \(2/3\) power and invoke summable near-subadditivity |
| Exact missing lemma | \(L_{\mathrm{IC}}\) in Section 4, all large \(n\le m\le2n\), excess \(KN/\log^2(e+N)\) after the \(2/3\) transform |
| Why strictly weaker | It only certifies one extension for low blocks, ignores all high-block cosets, and the resulting recurrence is compatible with every possible limiting constant |
| Falsification test | Compute \(\Gamma(A,B)\) in (11); a sequence satisfying (14), or a positive linear transformed excess, kills the route |
| Specialist confidence | Truth of \(L_{\mathrm{IC}}\): **20%**; tractability with current discrepancy methods: **8%**; conditional implication: **99%**; overall route: **about 7%** |

The low confidence is not caused by the sequence theorem, whose use is exact.  It reflects the thin margin left by the unavoidable bipartite interface lower bound (2) and the need to control exponentially many joint near-ground-state pairs with one deterministic Boolean interface.  Among the three frozen proposals, however, \(L_{\mathrm{IC}}\) is the only one that simultaneously has the correct quantifiers, an all-order realization, a genuinely \(o(N^{3/2})\) and summable error, and a conclusion mechanism that does not require identifying the limiting constant.
