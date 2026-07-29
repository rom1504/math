# Projective compression: exact width transfer and the integrality wall

## Status

This note tests the following proposed convergence mechanism.  Start
with a competitive signing \(B\) of order \(N\), partition its vertices
into \(n\) equal fibres, project each fibre onto one signed Boolean
mode, and round the resulting \(n\times n\) real matrix to a signing.
One would like to obtain
\[
 \frac{M(A)}{n^{3/2}}
 \le
 \frac{M(B)}{N^{3/2}}+o_n(1)
 \tag{1}
\]
uniformly for \(N\gg n\).

There is an exact centered-width compression inequality, but the
usual random projection destroys all of the information at the
relevant scale.  More precisely:

* the centered width contracts with exactly the desired
  \(s^{3/2}\) normalization;
* a random one-mode-per-fibre projection produces a fractional
  matrix converging to zero, even when \(B\) is an optimal signing;
* every sign rounding which stays close in centered width then pays
  a compulsory \(\Omega(n^{3/2})\) integrality cost;
* avoiding that cost requires an exceptionally coherent projection:
  the compressed coefficients must be \(L^2\)-close to signs;
* for a conference input, such coherence is equivalent to finding an
  almost invariant subspace spanned by disjointly supported Boolean
  vectors, with a flat compressed involution.

Thus random equipartition plus ordinary sign rounding is a rigorous
no-go.  A successful projective-compression proof would need a new
structural theorem producing the exceptional almost-invariant
Boolean block subspace.  Low Boolean quadratic norm alone is not
currently known to provide one.

The midpoint is independent of this issue: it is a fixed internal
energy offset of the block code.  Centered-width transfer does not
control it.

## 1. Notation

For a symmetric zero-diagonal real matrix \(D\), write
\[
 H_D(x)=\sum_{i<j}d_{ij}x_ix_j,
 \qquad
 M(D)=\max_{x\in\{\pm1\}^m}|H_D(x)|,
\]
and
\[
 W(D)=\frac{\max_xH_D(x)-\min_xH_D(x)}2,\qquad
 d(D)=\frac{\max_xH_D(x)+\min_xH_D(x)}2.
\]
Thus
\[
 M(D)=W(D)+|d(D)|.
 \tag{2}
\]

Let \(N=ns\), and partition \([N]\) into fibres
\[
 V_1\sqcup\cdots\sqcup V_n,\qquad |V_a|=s.
\]
Choose a sign vector
\(\sigma^{(a)}\in\{\pm1\}^{V_a}\) in each fibre.  Define
\[
 I_\sigma
 =
 \sum_{a=1}^n\ \sum_{\{i,j\}\subset V_a}
 b_{ij}\sigma_i\sigma_j
 \tag{3}
\]
and, for \(a<b\),
\[
 S_{ab}
 =
 \sum_{i\in V_a,\ j\in V_b}
 b_{ij}\sigma_i\sigma_j,
 \qquad
 C_{ab}=\frac{S_{ab}}{s^{3/2}}.
 \tag{4}
\]
For a macrospin \(x\in\{\pm1\}^n\), let
\[
 z_i(x)=\sigma_i x_a\qquad(i\in V_a).
 \tag{5}
\]

## 2. Exact compression theorem

For every macrospin \(x\),
\[
 \boxed{
 H_B(z(x))=I_\sigma+s^{3/2}H_C(x).
 }
 \tag{6}
\]
This gives the exact centered-width contraction
\[
 \boxed{
 W(C)\le \frac{W(B)}{s^{3/2}}.
 }
 \tag{7}
\]
Indeed, the energy range on the \(2^n\)-point block code
\(\{z(x)\}\) is a subinterval of the full energy range of \(B\), and
the translation \(I_\sigma\) cancels from the range.

For the absolute objective one only gets
\[
 \boxed{
 M(C)\le \frac{M(B)+|I_\sigma|}{s^{3/2}}.
 }
 \tag{8}
\]
If \(d_{\rm code}(B)\) denotes the midpoint of the energy interval of
\(B\) restricted to the block code, then the exact midpoint identity
is
\[
 \boxed{
 d(C)=\frac{d_{\rm code}(B)-I_\sigma}{s^{3/2}}.
 }
 \tag{9}
\]
Although averaging (6) over macrospins shows
\[
 |I_\sigma|\le M(B),
 \tag{10}
\]
there is no general \(o(N^{3/2})\) bound.  Consequently (7) transfers
centered width at the correct scale, while (9) leaves a separate
midpoint problem.

## 3. Random projection collapses to zero

Now take the fibre signs independently and uniformly at random.  For
every cross-fibre pair,
\[
 \mathbb E S_{ab}=0,\qquad
 \mathbb E S_{ab}^2=s^2.
 \tag{11}
\]
The second identity is exact: in the expansion of \(S_{ab}^2\), the
only surviving terms have the same endpoint in each of the two
fibres.  Therefore
\[
 \boxed{
 \mathbb E\sum_{a<b}C_{ab}^2
 =\frac1s\binom n2.
 }
 \tag{12}
\]
Likewise, distinct degree-two characters within a fibre are
orthogonal, so
\[
 \boxed{
 \mathbb E I_\sigma=0,\qquad
 \mathbb E I_\sigma^2=n\binom s2.
 }
 \tag{13}
\]

In particular, there is a deterministic choice of the fibre signs
for which, up to absolute constants,
\[
 \sum_{a<b}C_{ab}^2\le \frac{n^2}{s},
 \qquad
 \frac{|I_\sigma|}{s^{3/2}}\le \sqrt{\frac{2n}{s}}.
 \tag{14}
\]
For that choice, Cauchy--Schwarz gives
\[
 M(C)
 \le
 \sum_{a<b}|C_{ab}|
 \le
 O\!\left(\frac{n^2}{\sqrt s}\right).
 \tag{15}
\]
Hence whenever \(s/n\to\infty\),
\[
 \boxed{
 \frac{M(C)}{n^{3/2}}\longrightarrow0.
 }
 \tag{16}
\]
If \(s\gg n^2\), the second-moment bound and a union bound also give a
choice with
\[
 \max_{a<b}|C_{ab}|\le1,
 \tag{17}
\]
so the collapse occurs inside the fractional signing cube
\([-1,1]^{\binom n2}\).

This conclusion is independent of the quality of \(B\).  It applies
equally to an exact minimizer, a conference matrix, or an arbitrary
signing.

The scale mismatch has a simple interpretation.  A random cross-block
bilinear form has size \(s\), whereas scale-preserving compression
requires size \(s^{3/2}\).  A one-dimensional random coarse mode loses
a factor \(\sqrt s\).

## 4. A quantitative centered-width integrality obstruction

The collapse in Section 3 cannot be repaired by a sign rounding with
subleading centered-width error.

### Lemma 4.1

For every symmetric zero-diagonal real \(n\times n\) matrix \(D\),
\[
 \boxed{
 W(D)\ge
 \frac1{4\sqrt2}
 \sum_{i=1}^n
 \left(\sum_{j\ne i}d_{ij}^2\right)^{1/2}.
 }
 \tag{18}
\]

#### Proof

Use the exact cut identity
\[
 W(D)=\max_{S\subset[n]}
 \|D_{S,S^c}\|_{\infty\to1}.
 \tag{19}
\]
For fixed \(S\), choose independent random signs \(y_j\) on \(S^c\),
and then choose each row sign on \(S\) to agree with
\(\sum_{j\in S^c}d_{ij}y_j\).  The sharp \(L^1\) Khintchine inequality
gives
\[
 \|D_{S,S^c}\|_{\infty\to1}
 \ge
 \frac1{\sqrt2}
 \sum_{i\in S}
 \left(\sum_{j\in S^c}d_{ij}^2\right)^{1/2}.
 \tag{20}
\]
Choose \(S\) by putting every vertex in \(S\) independently with
probability \(1/2\).  If
\(R_i^2=\sum_{j\ne i}d_{ij}^2\) and
\(X_i=\sum_{j\notin S}d_{ij}^2\), then
\[
 \sqrt{X_i}\ge X_i/R_i,\qquad
 \mathbb E X_i=R_i^2/2.
 \]
Conditioning on \(i\in S\) therefore gives
\(\mathbb E\sqrt{X_i}\ge R_i/2\), while
\(\Pr(i\in S)=1/2\).  Averaging (20) proves (18). \(\square\)

If \(|d_{ij}|\le2\), then every row norm is at most \(2\sqrt n\), and
(18) yields
\[
 \boxed{
 W(D)\ge
 \frac{\|D\|_F^2}{8\sqrt{2n}}.
 }
 \tag{21}
\]

Now let \(A\) be any signing and let \(C\in[-1,1]^{\binom n2}\).
Putting \(D=A-C\), (21) shows
\[
 W(A-C)=o(n^{3/2})
 \quad\Longrightarrow\quad
 \|A-C\|_F^2=o(n^2).
 \tag{22}
\]
In particular,
\[
 \sum_{a\ne b}(1-|C_{ab}|)^2=o(n^2).
 \tag{23}
\]
Thus a fractional matrix can be rounded with subleading
centered-width error only if almost all of its entries are already
close to \(\pm1\) in mean square.

For the random compression in Section 3,
\(\|C\|_F=o(n)\) when \(s\to\infty\).  If, as we may in a proposed
uniform \(N\gg n\) compression theorem, we take \(s\gg n^2\), the same
choice can also be made with \(\|C\|_{\max}\le1\) by (17).
Consequently, uniformly over all signings \(A\),
\[
 \|A-C\|_F^2=(1-o(1))n^2
 \tag{24}
\]
and hence, along this regime,
\[
 \boxed{
 W(A-C)\ge
 \left(\frac1{8\sqrt2}-o(1)\right)n^{3/2}.
 }
 \tag{25}
\]
This is a deterministic leading-order loss.  No independence
assumption on the rounding is used.

There is an even simpler objective-level formulation.  Since every
signing satisfies the established universal lower bound
\[
 M(A)\ge(c_*+o(1))n^{3/2},
 \qquad
 c_*=0.336493364431\ldots,
 \tag{26}
\]
while (16) gives \(M(C)=o(n^{3/2})\), every sign purification of the
collapsed projection opens an integrality gap of at least
\((c_*-o(1))n^{3/2}\).

## 5. Why using the inherited budget is circular

Equation (7) only says
\[
 W(C)\le
 \frac{W(B)}{s^{3/2}}
 =
 \frac{W(B)}{N^{3/2}}\,n^{3/2}.
 \tag{27}
\]
For a random projection, the left side is much smaller than the right
side: it tends to zero after normalization.  One might try to round
the near-zero \(C\) using the unused budget on the right side.

But after the collapse, the projected object carries no
leading-order information about \(B\).  If \(B\) is taken along a
liminf sequence with normalized value \(c\), a theorem asserting that
one can round this collapsed object to a signing of order \(n\) with
value at most
\[
 (c+o_n(1))n^{3/2}
 \tag{28}
\]
is already the desired all-order realization theorem.  The
compression identity supplies no additional constraint from which
(28) follows.

Thus the random projection has not transferred the good large-order
geometry; it has erased it and left the original integrality problem.

## 6. The coherent alternative and its spectral meaning

To avoid (25), one must choose the partition and fibre signs
nonrandomly so that the normalized block sums
\[
 C_{ab}=s^{-3/2}\sigma_a^\top B_{ab}\sigma_b
 \tag{29}
\]
are already close to signs.  This means cross-block sums of order
\(s^{3/2}\), rather than their generic size \(s\).

For a conference input this condition has an exact invariant-subspace
interpretation.  Suppose
\[
 B^2=(N-1)I,\qquad U=\frac{B}{\sqrt{N-1}}.
 \tag{30}
\]
Let \(v_a=s^{-1/2}\sigma^{(a)}\), let \(V\) be the \(N\times n\)
matrix with columns \(v_a\), and put \(P=VV^\top\).  The columns of
\(V\) are orthonormal and supported on disjoint fibres.  If the
diagonal block energies are included in a full compressed matrix
\(\widetilde C\), then
\[
 \widetilde C
 =
 \sqrt{\frac{N-1}{s}}\,V^\top U V.
 \tag{31}
\]
Since \(U\) is orthogonal,
\[
 \boxed{
 \|\widetilde C\|_F^2
 =
 \frac{N-1}{s}
 \left(
 n-\|(I-P)UV\|_F^2
 \right).
 }
 \tag{32}
\]

The maximum possible right side is
\((1+o(1))n^2\).  If the off-diagonal entries of \(C\) are
\(L^2\)-close to signs, then
\(\|C\|_F^2=(1-o(1))n^2\), and (32) forces
\[
 \boxed{
 \|(I-P)UV\|_F^2=o(n).
 }
 \tag{33}
\]
Thus the block-Boolean subspace \({\rm ran}(V)\) must be almost
invariant under the conference involution.  In addition, the
compression \(V^\top U V\) must be asymptotically flat, with
off-diagonal entries of magnitude \(1/\sqrt n\).

For independent random fibre signs,
\[
 \mathbb E\|V^\top U V\|_F^2=O(n^2/N)=O(n/s),
 \tag{34}
\]
so
\[
 \mathbb E\|(I-P)UV\|_F^2=n-O(n/s).
 \tag{35}
\]
Random block modes therefore lose essentially all invariant energy.
Successful compression requires the opposite extreme.

In the exact case, (33) says that for every block pair the matrices
\(B_{ab}\) share compatible Boolean singular directions; schematically,
\[
 B_{ab}\sigma_b\approx \sqrt{s}\,a_{ab}\sigma_a.
 \tag{36}
\]
This is precisely the rigid compressed-lift geometry encountered in
the centered-width amplification program.

## 7. Exact sufficient target

A projective-compression proof of convergence would follow from the
following structural theorem.

> **Coherent compression target.**  For every competitive order-\(N\)
> signing \(B\), and every sufficiently large target order
> \(n=o(N)\), there are an equipartition, fibre signs, and an
> order-\(n\) signing \(A\) such that
> \[
> W(A-C)=o(n^{3/2}),
> \qquad
> |d_{\rm code}(B)-I_\sigma|=o(N^{3/2}).
> \tag{37}
> \]

The first condition, together with (7), transfers centered width.
The second condition, together with (9), removes the midpoint.

Lemma 4.1 shows that the first condition necessarily implies
\[
 \|A-C\|_F^2=o(n^2),
 \tag{38}
\]
so it cannot be obtained from a generic projection or a generic
rounding theorem.  In the flat-spectrum case it further implies the
almost-invariant Boolean block-subspace condition (33).

No theorem currently derives (33), (36), or (37) from the sole
hypothesis \(M(B)=O(N^{3/2})\).  Conversely, no argument confined to
random equipartitions and one random Boolean mode per fibre can prove
(37), because Sections 3 and 4 give an explicit leading-order
obstruction.

## 8. Verdict

The projective route clarifies the scale-transfer problem but does not
yet solve it.

* **Proved:** exact \(s^{3/2}\) transfer for centered width.
* **Proved:** exact, separate midpoint-offset identity.
* **Proved:** random one-mode block projections collapse to zero.
* **Proved:** any subleading-width sign purification requires
  \(L^2\)-near-saturation of the projected coefficients.
* **Proved for conference inputs:** near-saturation forces an almost
  invariant disjoint-support Boolean subspace with a flat compressed
  involution.
* **Decisive no-go:** random equipartition/projection plus ordinary
  sign rounding cannot yield convergence.
* **Remaining possible route:** prove a new coherent-compression
  theorem of the form (37); its content is the existence of the
  exceptional subspace in (33), not a routine rounding step.
