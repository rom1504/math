# Opposite-face rigidity: exact uncrossing and localization

## Status

This note attacks the remaining same-signing cap-face compatibility
problem from `common_gibbs_cap_law_rank.md`.

The main output is an exact four-point/uncrossing package for the two
opposite faces, followed by a finite-dimensional localization theorem.
It does **not** prove linear effective rank.  It proves that bounded
affine dimension of the two cap faces forces every positively aligned
edge into a bounded collection of exact face-difference cuts, and that
each such cut gives an exact two-block principal closure.  Thus the
low-dimensional obstruction really is a quotient/descent obstruction,
not merely a generic low-rank Boolean cloud.

The remaining gap is quantitative: sublinear matrix effective rank
does not presently imply bounded (or even suitably small) cap-face
affine dimension, and a balanced two-block closure does not by itself
preserve the \(n^{3/2}\) normalization.

Throughout,
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
P=\max H_A,\qquad Q=-\min H_A,\qquad
W=\frac{P+Q}{2}.
\]

## 1. One signed cut function contains both faces

Fix a positive ground state \(x^+\), switch it to the all-one vector,
and write
\[
B_{ij}=a_{ij}x_i^+x_j^+.
\]
For \(S\subseteq[n]\), let
\[
c(S)=\sum_{ij\in\delta(S)}B_{ij}.
\tag{1.1}
\]
If \(z^S\) is the spin vector which is \(-1\) on \(S\) and \(+1\)
outside \(S\), then
\[
c(S)=\frac{P-H_B(z^S)}2.
\tag{1.2}
\]
Consequently
\[
\boxed{0\le c(S)\le W\quad\text{for every }S.}
\tag{1.3}
\]
The positive face is represented by the zero-cut family
\[
\mathcal Z=\{S:c(S)=0\},
\tag{1.4}
\]
and the negative face is represented by the maximum-cut family
\[
\mathcal M=\{S:c(S)=W\}.
\tag{1.5}
\]

Thus opposite-face compatibility is exactly the assertion that the
two marginal laws live on the zero and maximum level sets of the
*same* signed cut function \(c\).

## 2. Exact four-point intersection inequalities

For two vertex sets \(S,T\), put
\[
I_B(S,T)
=
\sum_e B_e\,\mathbf1_{\delta(S)}(e)\mathbf1_{\delta(T)}(e).
\tag{2.1}
\]
The pointwise cut-vector identity
\[
\mathbf1_{\delta(S)}+\mathbf1_{\delta(T)}
=
\mathbf1_{\delta(S\triangle T)}
+2\mathbf1_{\delta(S)}\mathbf1_{\delta(T)}
\]
gives
\[
\boxed{
c(S)+c(T)=c(S\triangle T)+2I_B(S,T).
}
\tag{2.2}
\]
Combining (2.2) with (1.3) gives three exact opposite-face laws:
\[
\boxed{
\begin{array}{ll}
S,T\in\mathcal Z
&\Longrightarrow I_B(S,T)\le0,\\[1mm]
S,T\in\mathcal M
&\Longrightarrow I_B(S,T)\ge W/2,\\[1mm]
S\in\mathcal Z,\ T\in\mathcal M
&\Longrightarrow 0\le I_B(S,T)\le W/2.
\end{array}}
\tag{2.3}
\]
No PSD relaxation or first-moment cut-cone argument contains (2.3);
it uses simultaneous membership in the two opposite faces.

Let \(S,S'\) be independent samples from a law on \(\mathcal Z\), and
let \(T,T'\) be independent samples from a law on \(\mathcal M\).
Define the crossing probabilities
\[
p_e=\Pr(e\in\delta(S)),\qquad
q_e=\Pr(e\in\delta(T)).
\tag{2.4}
\]
Averaging (2.3) yields
\[
\boxed{
\sum_eB_ep_e^2\le0,\qquad
\sum_eB_eq_e^2\ge\frac W2,\qquad
0\le\sum_eB_ep_eq_e\le\frac W2.
}
\tag{2.5}
\]
The face equations themselves are
\[
\boxed{
\sum_eB_ep_e=0,\qquad
\sum_eB_eq_e=W.
}
\tag{2.6}
\]

In the fixed top gauge, the two correlation matrices obey
\[
R^+_e=1-2p_e,\qquad R^-_e=1-2q_e,
\tag{2.7}
\]
so their half-difference is \(q-p\).

One further useful consequence is
\[
\sum_eB_e(q_e^2-p_e^2)
=
\sum_e g_e(p_e+q_e)
\ge\frac W2,
\qquad
g_e:=B_e(q_e-p_e).
\tag{2.8}
\]
This says that a positive proportion of the width is carried on
coordinates which are not only aligned but also active in at least one
of the two marginals.

## 3. Common balance forces a dense stochastic edge set

Assume an exact common cap law has the coordinate balance
\[
\boxed{g_e=B_e(q_e-p_e)\le\lambda<1\quad\text{for every }e.}
\tag{3.1}
\]
Since
\[
\sum_eg_e=W,
\tag{3.2}
\]
the positive aligned set
\[
E_+=\{e:g_e>0\}
\]
satisfies
\[
\boxed{|E_+|\ge W/\lambda.}
\tag{3.3}
\]

Call an edge deterministic in one marginal when its cut indicator is
constant almost surely, equivalently when its crossing probability is
in \(\{0,1\}\).  If an edge is deterministic in both marginals, then
\[
g_e\in\{-1,0,1\}.
\]
Under \(\lambda<1\), such an edge cannot belong to \(E_+\).
Therefore
\[
\boxed{
E_+\subseteq E_{\rm var}^+\cup E_{\rm var}^-,
\qquad
|E_{\rm var}^+\cup E_{\rm var}^-|\ge W/\lambda.
}
\tag{3.4}
\]

At the natural common-Gibbs scale
\(\lambda=b/\sqrt n\) and \(W\ge c n^{3/2}\), (3.4) becomes
\[
|E_{\rm var}^+\cup E_{\rm var}^-|
\ge(c/b)n^2.
\tag{3.5}
\]
Thus exact opposite-face balance cannot be supported on a sparse set
of cap-varying coordinates.

## 4. Low cap-face dimension localizes all variation to cuts

Let
\[
\mathcal V_+
=
\{(x_ix_j)_{i<j}:x\text{ is a positive ground state}\}
\]
and define \(\mathcal V_-\) analogously.  Write
\[
d_\pm=\dim\operatorname{aff}\mathcal V_\pm.
\tag{4.1}
\]

### Proposition 4.1 (face-difference cut cover)

There are \(d_+\) differences of positive ground states and \(d_-\)
differences of negative ground states whose cut supports cover every
edge which is variable in the corresponding face.  In particular,
\[
\boxed{
E_{\rm var}^+\cup E_{\rm var}^-
\subseteq
\bigcup_{r=1}^{d_++d_-}\delta(D_r).
}
\tag{4.2}
\]

#### Proof

Choose a reference vertex \(v_0\) of \(\mathcal V_+\), and choose
\(d_+\) further vertices \(v_1,\ldots,v_{d_+}\) such that
\(v_r-v_0\) form a basis of the affine tangent space.  If coordinate
\(e\) is constant on these \(d_++1\) vertices, its coordinate
functional vanishes on the tangent basis and hence is constant on the
whole affine hull.  Therefore every variable coordinate belongs to
the support of some \(v_r-v_0\).

For two cut-correlation vertices, the support of their difference is
exactly the edge cut associated with their coordinatewise spin
difference.  This proves the positive-face assertion.  Apply the same
argument to the negative face. \(\square\)

Each \(D_r\) in the positive family is a zero cut after switching its
positive reference state to \(1\).  Each \(D_r\) in the negative
family is a zero cut after negating the signing and switching its
negative reference state to \(1\).  Thus every cut in (4.2) is an
exact same-face zero cut in its natural endpoint gauge.

Combining (3.3) and (4.2), one of the face-difference cuts obeys
\[
|\delta(D_r)|
\ge
\frac{W}{\lambda(d_++d_-)}.
\tag{4.3}
\]
If \(s_r=\min\{|D_r|,n-|D_r|\}\), then
\[
\boxed{
s_r
\ge
\frac{W}{\lambda(d_++d_-)n}.
}
\tag{4.4}
\]
For fixed \(b\), \(\lambda=b/\sqrt n\), \(W\ge cn^{3/2}\), and
\(d_++d_-=O(1)\), this is a genuinely macroscopic cut.

There is also a bridge from exact marginal matrix rank.  If a cap
law has correlation matrix of rank \(r\), all spin vectors in its
support lie in an \(r\)-dimensional linear subspace.  Their quadratic
edge features therefore lie in a space of dimension at most
\[
\binom{r+1}{2}.
\tag{4.5}
\]
Hence
\[
d_\pm\le\binom{r_\pm+1}{2}-1.
\tag{4.6}
\]
This only converts \(r_\pm=o(\sqrt n)\) into sublinear face dimension;
it does not settle the square-root-rank obstruction.

## 5. Every localized cut is an exact principal two-block closure

The following elementary lemma converts (4.2) into genuine principal
structure.

### Proposition 5.1 (zero-cut principal closure)

Let \(C\) be any signing whose all-one vector is a positive ground
state, and suppose
\[
\sum_{ij\in\delta(S)}C_{ij}=0.
\tag{5.1}
\]
Then
\[
\boxed{
p(C)=p(C[S])+p(C[S^c]),
}
\tag{5.2}
\]
the all-one vector is a positive ground state of both principal
blocks, and
\[
\boxed{
\nu(C)\ge\nu(C[S])+\nu(C[S^c]),
\qquad
W(C)\ge W(C[S])+W(C[S^c]).
}
\tag{5.3}
\]

#### Proof

The all-one energy is
\[
p(C)
=
H_{C[S]}(\mathbf1)+H_{C[S^c]}(\mathbf1)
\]
because the cross-block total is zero.  On the other hand, take
positive ground states independently in the two principal blocks.
Flipping every spin in one block reverses the cross term and leaves
both internal energies fixed, so one orientation has nonnegative
cross energy.  Hence
\[
p(C)\ge p(C[S])+p(C[S^c])
\ge
H_{C[S]}(\mathbf1)+H_{C[S^c]}(\mathbf1)
=p(C).
\]
All inequalities are equalities, proving (5.2) and principal
all-one maximality.

Apply the same block-flip argument to negative ground states of the
two principal blocks, now choosing the orientation with nonpositive
cross energy.  This proves the first inequality in (5.3); adding
(5.2) proves the width inequality. \(\square\)

Switching and global negation preserve principal widths.  Therefore
every cut in the face-difference cover (4.2) gives an exact
two-block principal closure of the original signing, in its natural
endpoint gauge.

The equality case contains more structure than (5.2)--(5.3).

### Proposition 5.2 (Cartesian ground closure and cross annihilation)

In the setting of Proposition 5.1, let
\(\mathcal G_S,\mathcal G_{S^c}\) be the positive ground-state sets of
the two principal blocks, and let
\[
L_S=\operatorname{span}\mathcal G_S,\qquad
L_{S^c}=\operatorname{span}\mathcal G_{S^c}.
\tag{5.4}
\]
If \(C_{S,S^c}\) denotes the cross block, then
\[
\boxed{
u^\mathsf TC_{S,S^c}v=0
\quad
(u\in\mathcal G_S,\ v\in\mathcal G_{S^c}),
}
\tag{5.5}
\]
and hence
\[
\boxed{
P_{L_S}C_{S,S^c}P_{L_{S^c}}=0.
}
\tag{5.6}
\]
Every concatenation \((u,v)\), with \(u,v\) chosen independently from
the two principal positive-ground sets, is a positive ground state of
the full signing.

More generally, for arbitrary Boolean \(u,v\), the same two block
orientations give the exact capped-bilinear inequality
\[
\boxed{
\left|u^\mathsf TC_{S,S^c}v\right|
\le
\bigl[p(C[S])-H_{C[S]}(u)\bigr]
+
\bigl[p(C[S^c])-H_{C[S^c]}(v)\bigr].
}
\tag{5.7}
\]
Every further zero cut which crosses \(S\) is an equality profile in
(5.7): its positive full-ground energy forces its cross-block gain to
equal the sum of the two principal energy deficits.  Thus the
nonlaminar zero-cut problem is equivalently an inverse problem for
the equality cases of (5.7).

#### Proof

By (5.2), the two internal energies already sum to \(p(C)\).
The two full states \((u,v)\) and \((u,-v)\) have energies
\[
p(C)\pm u^\mathsf TC_{S,S^c}v.
\]
Neither may exceed \(p(C)\), so the cross term is zero and both states
are positive grounds.  Bilinearity then gives (5.6). \(\square\)

Writing \(n_1=|S|\), \(n_2=|S^c|\), and
\(r_i=\dim L_i\), (5.6) implies the exact rank obstruction
\[
\boxed{
\operatorname{rank}C_{S,S^c}
\le(n_1-r_1)+(n_2-r_2).
}
\tag{5.8}
\]
Indeed, \(C_{S,S^c}\) maps \(L_{S^c}\) into \(L_S^\perp\);
the complementary domain has dimension \(n_2-r_2\).

Thus a macroscopic zero-cut split has a sharp dichotomy:

* a full-rank cross block forces a linear total deficit in the two
  principal ground-state spans;
* small ground-span deficit forces the cross sign block itself to
  have low matrix rank, which is a concrete quotient/lift signature.

Repeated zero-cut splitting also multiplies positive-ground sets.
After a laminar recursion into \(\ell\) leaves, independently changing
the signs of leaf ground states gives a type-constant
\(2^{\ell-1}\)-state positive-ground cube.  The affine type-closure
theorem then applies to the leaf partition.  In particular, a deep
balanced zero-cut tree cannot remain hidden: it creates a large
explicit affine ground quotient.  What is not yet proved is that an
arbitrary *crossing* zero-cut family admits a laminar refinement with
comparable macroscopic atoms.

This proves the following finite-dimensional form of the desired
opposite-face alternative:

> If one common exact cap law has bounded total cap-face affine
> dimension, then its aligned cancellation is supported on a bounded
> family of same-face zero cuts; at the natural balance scale (with
> fixed \(b\)) one of those cuts is macroscopic, and it yields an exact
> principal two-block closure.

## 6. Precise obstruction to completion

Two gaps remain.

1. **Effective rank versus face dimension.**  A Boolean correlation
   matrix of rank \(r\) can support edge features of affine dimension
   \(\Theta(r^2)\).  At the sharp obstruction \(r\asymp\sqrt n\),
   Proposition 4.1 can require \(\Theta(n)\) cuts.  PSD, the Boolean
   triangle identities, and (2.5) have not yet reduced this cover to
   bounded depth.

2. **Balanced closure versus scale-preserving descent.**  Proposition
   5.1 gives exact superadditivity
   \(W(C)\ge W(C[S])+W(C[S^c])\), but for a balanced split
   \[
   |S|^{3/2}+|S^c|^{3/2}<n^{3/2}.
   \]
   Thus this scalar budget alone loses a fixed normalization factor.
   A successful continuation must either show that the zero-cut
   decomposition has a nearly macroscopic atom, or extract an
   additional cross-block payment.

The four-point inequalities (2.3) are the strongest exact
opposite-face information found in this audit.  They should be kept
in any subsequent PSD, laminarity, or principal-descent attack.

### 6.1 Uniform near-cap versions

The uncrossing statements are stable at exactly the error scale needed
for the conditioned common Gibbs law.  Suppose, in the fixed exact-top
gauge,
\[
\mathcal Z_\alpha=\{S:c(S)\le\alpha\},\qquad
\mathcal M_\alpha=\{S:c(S)\ge W-\alpha\}.
\tag{6.1}
\]
Then (2.2) gives
\[
\boxed{
\begin{array}{ll}
S,T\in\mathcal Z_\alpha
&\Longrightarrow I_B(S,T)\le\alpha,\\[1mm]
S,T\in\mathcal M_\alpha
&\Longrightarrow I_B(S,T)\ge W/2-\alpha,\\[1mm]
S\in\mathcal Z_\alpha,\ T\in\mathcal M_\alpha
&\Longrightarrow
-\alpha/2\le I_B(S,T)\le W/2+\alpha/2.
\end{array}}
\tag{6.2}
\]

There is also a stable principal-closure lemma.  If
\[
0\le c(S)\le\alpha,
\tag{6.3}
\]
then, writing \(C_1=C[S]\), \(C_2=C[S^c]\),
\[
\boxed{
P-\alpha\le p(C_1)+p(C_2)\le P.
}
\tag{6.4}
\]
Moreover, for every pair of principal positive grounds \(u,v\),
\[
\boxed{
\left|u^\mathsf TC_{S,S^c}v\right|\le\alpha.
}
\tag{6.5}
\]
Indeed, the all-one internal energy is \(P-c(S)\), while the universal
block-flip argument gives \(p(C_1)+p(C_2)\le P\).  Once (6.4) is
known, the two extensions \((u,v)\) and \((u,-v)\) prove (6.5).

For the common shell in `common_gibbs_cap_law_rank.md`,
\(\alpha=O(n^{3/2}/b_n)=o(n^{3/2})\).  Thus every localized
near-face difference cut gives an \(o(n^{3/2})\)-approximate Cartesian
ground product.  The remaining stability issue is linear-algebraic:
small bilinear values on Boolean ground generators do not control the
whole span unless those generators have a quantitative frame lower
bound.

## 7. Hadamard-moment face inequalities

There is an infinite extension of the four-point law which is useful
for auditing proposed PSD arguments.  Let \(R^+\) and \(R^-\) be the
two cap correlation matrices in the fixed top gauge.  For every
positive integer \(m\), the entrywise powers
\[
(R^\pm)^{\circ m}
\]
are again Boolean correlation matrices: sample \(m\) independent cap
states and multiply them coordinatewise.  The product state need not
remain in the cap, but its energy remains in the global interval.
Therefore
\[
\boxed{
-Q\le
\left\langle B,(R^\pm)^{\circ m}\right\rangle
\le P.
}
\tag{7.1}
\]
For odd \(m\), define the nonnegative divided difference
\[
K_m(r,s)=
\begin{cases}
\dfrac{r^m-s^m}{r-s},&r\ne s,\\[2mm]
mr^{m-1},&r=s.
\end{cases}
\tag{7.2}
\]
With
\[
g_e=\frac12B_e(R^+_e-R^-_e),
\]
subtracting the two appropriate bounds in (7.1) gives
\[
\boxed{
\sum_e g_eK_m(R^+_e,R^-_e)\le W
\qquad(m\ \text{odd}).
}
\tag{7.3}
\]
Since \(\sum_eg_e=W\),
\[
\boxed{
\sum_e g_e\bigl(K_m(R^+_e,R^-_e)-1\bigr)\le0.
}
\tag{7.4}
\]

For \(m=3\), (7.3) is the sum of the two nonnegative closure defects
\[
\begin{aligned}
\Delta_+
&=
P-\langle B,(R^+)^{\circ3}\rangle,\\
\Delta_-
&=
\langle B,(R^-)^{\circ3}\rangle+Q,
\end{aligned}
\qquad
\Delta_++\Delta_-
=2W-2\sum_eg_eK_3(R^+_e,R^-_e).
\tag{7.5}
\]
If \(\Delta_+=0\) and the positive marginal has full support on its
face, the product of every three positive support states is again a
positive ground state.  The support is then a torsor under
coordinatewise multiplication, hence an affine \(\mathbb F_2\)
ground family; the exact affine type-closure theorem applies.
The analogous statement holds on the negative face.

Thus the third-product defect is a precise, falsifiable distance from
the already-understood affine quotient branch.  The unresolved case
has both a substantial product-closure defect and low marginal
effective rank.  Bare PSD cannot exclude that case: entrywise powers
remain in the correlation polytope and (7.1) is exactly the full
information supplied by independent product probes.

## 8. Two affine faces and the signed-incidence insertion model

There is a clean insertion reduction when both cap families are
affine.  Write
\[
x_i(w)=\alpha_i(-1)^{\tau_i\cdot w},\qquad
y_i(u)=\beta_i(-1)^{\sigma_i\cdot u}.
\tag{8.1}
\]
Make a bipartite multigraph whose left vertices are the occupied
\(\tau\)-types, whose right vertices are the occupied \(\sigma\)-types,
and whose edge \(i\) joins \(\tau_i\) to \(\sigma_i\).  Give edge \(i\)
the sign
\[
\gamma_i=\alpha_i\beta_i.
\tag{8.2}
\]
For a proposed new-vertex signing \(b_i\), put
\(c_i=b_i\alpha_i\).  Its Fourier coefficients on the two cap
families are exactly
\[
\widehat f_+(\tau)=\sum_{i:\tau_i=\tau}c_i,\qquad
\widehat f_-(\sigma)=\sum_{i:\sigma_i=\sigma}\gamma_ic_i.
\tag{8.3}
\]
Thus simultaneous insertion is the discrepancy problem for a signed
bipartite incidence matrix.

If every signed cycle is unfrustrated,
\[
\prod_{i\in C}\gamma_i=+1,
\tag{8.4}
\]
then \(\gamma\) is a vertex coboundary and row switching reduces
(8.3) to ordinary bipartite incidence.  Alternating signs on Euler
cycles make every even type sum zero.  The affine parity theorem says
that each type partition has at most one odd class, so only the
possible two trail endpoints have residual magnitude one.  In this
unfrustrated case one obtains
\[
\boxed{
\sup_w\left|\sum_i b_ix_i(w)\right|=O(1),\qquad
\sup_u\left|\sum_i b_iy_i(u)\right|=O(1).
}
\tag{8.5}
\]

For a single signed cycle, (8.4) is also necessary for zero
discrepancy: propagating the alternating equations around the cycle
returns consistently exactly when the cycle sign is \(+1\).  A
frustrated cycle forces a residual of magnitude \(2\) somewhere.
The unresolved same-\(A\) question in this language is whether many
independent frustrated components can occur without producing either
a width witness or one of the zero-cut quotients in Sections 4--5.

There is, however, an important obstruction to applying this exact
affine insertion picture directly to the common Gibbs law.  If the
two marginals are **uniform** on the full affine families (8.1), then
\[
R^+_{ij}
=
\alpha_i\alpha_j\mathbf1_{\{\tau_i=\tau_j\}},
\qquad
R^-_{ij}
=
\beta_i\beta_j\mathbf1_{\{\sigma_i=\sigma_j\}}.
\tag{8.6}
\]
Consequently
\[
g_{ij}
=
\frac12a_{ij}(R^+_{ij}-R^-_{ij})
\in\{0,\pm\tfrac12,\pm1\}.
\tag{8.7}
\]
If the common edge balance has \(\lambda<1/2\), every positive value
in (8.7) is forbidden.  This contradicts
\(\sum_eg_e=W>0\).  Therefore:
\[
\boxed{
\text{the two exact uniform affine faces cannot themselves be the
common }\lambda=o(1)\text{ cap marginals.}
}
\tag{8.8}
\]
At least one marginal must be a genuinely thick or nonuniform
near-cap law.  The signed-cycle insertion idea is therefore a model
for the limiting quotient branch, not a direct completion of the
finite-temperature proof.

## 9. Audited theorem and remaining obstruction

The proved same-signing chain is:
\[
\begin{gathered}
\text{opposite cap faces}
\Longrightarrow
\text{zero/max levels of one cut function}
\Longrightarrow
\text{intersection laws (2.3)},\\
\text{common balance}
\Longrightarrow
\Omega(W/\lambda)\text{ stochastic aligned edges},\\
\text{bounded cap-face affine dimension}
\Longrightarrow
\text{a bounded face-difference cut cover},\\
\text{a face-difference zero cut}
\Longrightarrow
\text{exact principal closure, Cartesian grounds, and
cross annihilation (5.6)}.
\end{gathered}
\tag{9.1}
\]

This is a genuine opposite-face rigidity theorem, but not yet the
desired linear-rank theorem.  The sharp remaining obstruction can be
stated without ambiguity:

> Can a competitive same-signing common cap law have
> square-root-scale marginal rank, substantial cubic product-closure
> defect on both faces, and a face-difference cut cover whose
> crossing equality profiles in (5.7) admit no nearly macroscopic
> principal atom?

Neither PSD nor the full independent Hadamard-moment hierarchy rules
this out.  Resolving it requires an inverse theorem for the equality
profiles of the capped bilinear inequality (5.7), or a mechanism
turning their crossing structure into an additional
\(n^{3/2}\)-scale cross-block payment.
