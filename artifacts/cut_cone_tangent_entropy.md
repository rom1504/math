# Cut-cone tangent entropy and feature-span degeneracy

## Status

This note audits the proposed entropy-versus-tangent route under the
*full endpoint structure* of a signing.  It proves:

1. the exact cut-cone tangent and capped-cut identities;
2. a pointwise conversion from low \(c\)-slack to high \(h\)-energy;
3. a covariance/entropy inequality for every low-slack cut family;
4. an exact feature-covariance kernel theorem at zero slack, and its
   quantitative low-eigenvalue version at positive slack;
5. a Schur-complement form of the degeneracy, suitable for block
   replacement;
6. an exact global-minimality obstruction phrased as a discrepancy
   problem on the low-slack feature cloud.

The desired traffic entropy bound is **not** proved.  The constants show
that ordinary vertex covariance is intrinsically too weak.  What survives
is a precise dichotomy target: either the shell is small enough for traffic,
or its internal edge features have a flat, conditionally predictable
low-variance direction aligned with the midpoint signing \(h\).  Turning
that predictable direction into a sign-valued block replacement remains
open.

---

## 1. Endpoint gauge and the tangent pair

Let \(A\) be a signing, let \(x^+\) and \(x^-\) be top and bottom states,
and write

\[
P=\max_xH_A(x),\qquad m=\min_xH_A(x),
\]

\[
W=\frac{P-m}{2},\qquad d=\frac{P+m}{2}.
\]

Switch \(x^+\) to the all-plus state and put

\[
u_i=x_i^+x_i^-.
\]

Define the two all-plus-maximized signings

\[
b^+_{ij}=a_{ij}x_i^+x_j^+,
\qquad
b^-_{ij}=-a_{ij}x_i^-x_j^-,
\]

and

\[
c=\frac{b^++b^-}{2},\qquad
h=\frac{b^+-b^-}{2}.
\tag{1.1}
\]

Then

\[
c\pm h=b^\pm\in \operatorname{CUT}_n^*,
\tag{1.2}
\]

\[
c\cdot\mathbf1=W,\qquad h\cdot\mathbf1=d.
\tag{1.3}
\]

Moreover,

\[
c_{ij},h_{ij}\in\{-1,0,1\},\qquad c_{ij}h_{ij}=0,
\qquad c_{ij}^2+h_{ij}^2=1.
\tag{1.4}
\]

If

\[
U=\{i:u_i=1\},\qquad V=\{i:u_i=-1\},
\]

then \(c\) is supported exactly on \(U\times V\), while \(h\) is
supported exactly on

\[
E_h=\binom U2\sqcup\binom V2.
\tag{1.5}
\]

Thus \(c\) is a rectangular signing and \(h\) is the union of the two
internal signings.

For a cut \(S\), write \(\delta(S)\) for its edge-incidence vector.  From
(1.2),

\[
\boxed{
|h\cdot\delta(S)|\le c\cdot\delta(S)
\quad\text{for every }S\subseteq[n].
}
\tag{1.6}
\]

This is the exact two-sided tangent inequality.

---

## 2. Capped-cut form

Represent a cut by \(z\in\{\pm1\}^n\), modulo global negation, and put

\[
H_c(z)=\sum_ec_ez_iz_j,\qquad
H_h(z)=\sum_eh_ez_iz_j.
\]

Since

\[
c\cdot\delta(z)=\frac{W-H_c(z)}2,
\qquad
h\cdot\delta(z)=\frac{d-H_h(z)}2,
\]

(1.6) is equivalent to

\[
\boxed{
|d-H_h(z)|\le W-H_c(z).
}
\tag{2.1}
\]

Flipping every spin in \(U\) changes \(H_c\) to \(-H_c\) and leaves
\(H_h\) unchanged.  Applying (2.1) to whichever of these two
orientations has nonnegative cross energy gives the sharper symmetric
form

\[
\boxed{
|d-H_h(z)|\le W-|H_c(z)|.
}
\tag{2.2}
\]

Writing \(z=(s,t)\in\{\pm1\}^U\times\{\pm1\}^V\), with rectangular
matrix \(C=(c_{ij})_{U\times V}\), this is

\[
\left|
d-H_{h_U}(s)-H_{h_V}(t)
\right|
\le
W-|s^\top Ct|.
\tag{2.3}
\]

In particular, \(W=\|C\|_{\infty\to1}\), with the all-plus pair chosen
as one maximizer.

Equation (2.3) is the useful dual cap identity: on every near-maximizer
of the rectangular form, the internal quadratic form is forced into a
thin slab about \(d\).

---

## 3. Low slack forces high tangent energy

For \(r\ge0\), define

\[
\mathcal L_r
=
\left\{
z:\ c\cdot\delta(z)\le r n^{3/2}
\right\}.
\tag{3.1}
\]

For every \(z\in\mathcal L_r\),

\[
|d-H_h(z)|\le2r n^{3/2}.
\tag{3.2}
\]

If \(d\ge\varepsilon n^{3/2}\), then

\[
\boxed{
H_h(z)\ge(\varepsilon-2r)n^{3/2}
\quad(z\in\mathcal L_r).
}
\tag{3.3}
\]

This is stronger than an average statement: the entire low-slack cloud
lies in one quadratic halfspace.

Let \(\mu_r\) be any probability measure supported on
\(\mathcal L_r\), symmetrized under \(z\mapsto-z\), and let

\[
R=\mathbb E_{\mu_r}zz^\top.
\tag{3.4}
\]

Then \(R\succeq0\), \(R_{ii}=1\), and, with
\(\Delta=(\varepsilon-2r)_+\),

\[
\sum_{i<j}h_{ij}R_{ij}\ge\Delta n^{3/2}.
\tag{3.5}
\]

Since \(|E_h|\le\binom n2\), Cauchy--Schwarz gives

\[
\sum_{i<j}R_{ij}^2
\ge
\frac{\Delta^2n^3}{|E_h|}
\ge
(2\Delta^2-o(1))n.
\tag{3.6}
\]

Equivalently,

\[
\boxed{
\operatorname{tr}R^2-n
\ge
(4\Delta^2-o(1))n.
}
\tag{3.7}
\]

Thus a macroscopic midpoint forces at least linear Frobenius excess in
the vertex covariance of every sufficiently low-slack family.

---

## 4. A vertex covariance/entropy lemma

### Lemma 4.1

Let \(Z\in\{\pm1\}^n\) have a globally sign-symmetric law, covariance
\(R=\mathbb EZZ^\top\), and Shannon entropy in natural units
\(\mathsf H(Z)\).  Then

\[
\boxed{
\mathsf H(Z)
\le
n\log2
-
\frac{\operatorname{tr}R^2-n}
{8\|R\|_{\mathrm{op}}}.
}
\tag{4.1}
\]

### Proof

Choose a vertex bipartition \(U\sqcup V=[n]\).  Conditional
subadditivity gives

\[
\mathsf H(Z)
\le
|U|\log2+\sum_{i\in V}\mathsf H(Z_i\mid Z_U).
\tag{4.2}
\]

Put \(m_i(Z_U)=\mathbb E(Z_i\mid Z_U)\).  The binary entropy inequality

\[
h_{\mathrm b}\!\left(\frac{1+t}{2}\right)
\le\log2-\frac{t^2}{2}
\]

implies

\[
\mathsf H(Z_i\mid Z_U)
\le\log2-\frac12\mathbb E m_i(Z_U)^2.
\tag{4.3}
\]

The conditional expectation is the best \(L^2\) predictor.  Comparing
it with the best linear predictor from \(Z_U\),

\[
\mathbb E m_i(Z_U)^2
\ge
R_{iU}R_{UU}^{\dagger}R_{Ui}
\ge
\frac{\|R_{iU}\|_2^2}{\|R\|_{\mathrm{op}}}.
\tag{4.4}
\]

The covariance vector \(R_{Ui}\) lies in the range of \(R_{UU}\), so
the pseudoinverse formula is legitimate.

For a uniformly random bipartition,

\[
\mathbb E_{U,V}
\sum_{i\in V,j\in U}R_{ij}^2
=
\frac14(\operatorname{tr}R^2-n).
\]

Some bipartition therefore has at least this cross mass.  Combining
(4.2)--(4.4) proves (4.1). \(\square\)

For the uniform cut law on a family \(\mathcal L\), the symmetric spin
lift has entropy \(\log|\mathcal L|+\log2\).  Combining (3.7) and
(4.1) gives

\[
\frac1n\log|\mathcal L_r|
\le
\log2-\frac{\Delta^2}{2\|R\|_{\mathrm{op}}}+o(1).
\tag{4.5}
\]

### Quantitative audit

The \(s=4\) traffic criterion needs shell entropy below

\[
\frac5{64}r^2.
\tag{4.6}
\]

Even in the best possible covariance regime
\(\|R\|_{\mathrm{op}}=1\), (4.5) only subtracts
\(\Delta^2/2\) from \(\log2\).  It cannot approach the small right side
of (4.6) in the relevant parameter cage.  Therefore:

\[
\boxed{
\text{ordinary vertex covariance cannot prove the traffic bound.}
}
\tag{4.7}
\]

This is a genuine scale obstruction, not merely a missing constant.

---

## 5. Internal edge-feature covariance

The sharper object is the internal edge feature

\[
\phi(z)=(z_iz_j)_{ij\in E_h}\in\{\pm1\}^{E_h}.
\tag{5.1}
\]

Let

\[
\bar\phi=\mathbb E_{\mu_r}\phi,
\qquad
\Sigma_r
=
\mathbb E_{\mu_r}
(\phi-\bar\phi)(\phi-\bar\phi)^\top.
\tag{5.2}
\]

Equation (3.2) says that the scalar \(h\cdot\phi(z)\) lies in an
interval of length \(4rn^{3/2}\).  Popoviciu's variance inequality
therefore gives

\[
\boxed{
h^\top\Sigma_r h
\le4r^2n^3.
}
\tag{5.3}
\]

At zero slack this becomes an exact kernel identity:

\[
\boxed{
\Sigma_0h=0.
}
\tag{5.4}
\]

Indeed, \(h\cdot\phi(z)=d\) on \(\mathcal L_0\).

Let \(P_{>\lambda}\) be the spectral projector of \(\Sigma_r\) on
eigenvalues greater than \(\lambda\).  From (5.3),

\[
\boxed{
\|P_{>\lambda}h\|_2^2
\le
\frac{4r^2n^3}{\lambda}.
}
\tag{5.5}
\]

Since

\[
|E_h|
=
\binom{|U|}{2}+\binom{|V|}{2}
\ge\frac{n^2}{4}-O(n)
\tag{5.6}
\]

and \(\|h\|_2^2=|E_h|\), choosing
\(\lambda=Kr^2n\) yields

\[
\|P_{\le Kr^2n}h\|_2^2
\ge
\left(1-\frac{16}{K}-o(1)\right)\|h\|_2^2.
\tag{5.7}
\]

Thus a low-slack cloud with a macroscopic midpoint always carries an
almost-flat, low-variance internal edge direction.  This is the exact
covariance degeneracy missed by vertex covariance.

---

## 6. Feature span and entropy

Let \(r_\phi=\operatorname{rank}\Sigma_r\), the affine dimension of the
feature cloud.

An \(r_\phi\)-dimensional affine subspace of \(\mathbb R^{E_h}\)
intersects \(\{\pm1\}^{E_h}\) in at most \(2^{r_\phi}\) points: choose
\(r_\phi\) coordinates on which the affine projection is injective.
The internal features determine the spins within \(U\) and within
\(V\), each up to a sign.  Modulo global negation, each internal
feature vector has at most two cut preimages.  Hence

\[
\boxed{
|\mathcal L_r|\le2^{r_\phi+1}.
}
\tag{6.1}
\]

Consequently, if

\[
\log|\mathcal L_r|>\frac5{64}r^2n,
\]

then necessarily

\[
r_\phi
\ge
\frac{5}{64\log2}r^2n-1.
\tag{6.2}
\]

At zero slack, (5.4) says simultaneously that this potentially
high-dimensional affine cloud lies in the hyperplane normal to the
flat sign vector \(h\).  Equations (5.4) and (6.2) are the precise
feature-span formulation of the surviving obstruction.

The rank lower bound is only order \(n\), while
\(|E_h|=\Theta(n^2)\).  Thus dimension counting alone leaves an
enormous polar space and cannot force a contradiction.

---

## 7. Conditional block degeneracy

For an internal edge block \(T\subseteq E_h\), let

\[
\Sigma_{T\mid T^c}
=
\Sigma_{TT}
-
\Sigma_{T,T^c}\Sigma_{T^c,T^c}^{\dagger}
\Sigma_{T^c,T}
\tag{7.1}
\]

be the Schur complement, interpreted on the supported range.  It is
the covariance left in \(\phi_T\) after optimal linear prediction from
\(\phi_{T^c}\).

The variational characterization of the Schur complement gives

\[
h_T^\top\Sigma_{T\mid T^c}h_T
\le h^\top\Sigma_rh.
\]

Therefore

\[
\boxed{
h_T^\top\Sigma_{T\mid T^c}h_T
\le4r^2n^3
\quad\text{for every }T\subseteq E_h.
}
\tag{7.2}
\]

At \(r=0\), the signed block statistic \(h_T\cdot\phi_T\) is determined
affinely by the features outside \(T\), on the entire exact cap
family.  Equivalently,

\[
\Sigma_{T\mid T^c}h_T=0.
\tag{7.3}
\]

This is the rigorous block-localization statement currently
available:

> every candidate replacement block carrying \(h\)-mass is either
> conditionally degenerate in the \(h_T\) direction, or its low-slack
> shell has positive quadratic variation large enough to contradict
> (7.2).

A restricted-invertibility argument can select large well-conditioned
coordinate sets when the Schur complements have large stable rank.
Equation (7.2) shows that no such set can be well-conditioned *in the
flat direction \(h_T\)* at zero slack.  What is still missing is an
inverse theorem converting this conditional edge-feature
predictability into an induced vertex block on which a
\(\{\pm1\}\)-valued replacement can be made.

---

## 8. Exact replacement obstruction

Let \(\beta\in\{\pm1\}^{E_h}\) be a new internal signing.  Flipping all
spins in \(U\) sends \(H_c\) to \(-H_c\) while preserving
\(H_\beta\).  Consequently,

\[
\boxed{
\max_z|H_c(z)+H_\beta(z)|
=
\max_z\left\{
|H_c(z)|+|H_\beta(z)|
\right\}.
}
\tag{8.1}
\]

Thus the continuous midpoint \(c\) has absolute norm exactly \(W\),
and a sign rounding \(c+\beta\) stays below \(W+e\) provided

\[
|H_\beta(z)|
\le
W-|H_c(z)|+e
\quad\text{for every }z.
\tag{8.2}
\]

If the original signing is globally width-minimizing, no internal
sign replacement can lower the width.  In particular, every
\(\beta\) must have an obstructing cap profile at the corresponding
accuracy:

\[
\boxed{
\max_z
\left\{
|H_\beta(z)|-[W-|H_c(z)|]
\right\}
\ge0.
}
\tag{8.3}
\]

This is the exact polar/discrepancy form of global minimality.

For random independent \(\beta_e\), \(H_\beta(z)\) is a Rademacher
process indexed by the internal feature cloud \(\phi(z)\).  A union
bound recovers the traffic entropy threshold.  If that bound fails,
the correct replacement statistic is not cardinality but the
multiscale metric entropy of

\[
d_\phi(z,z')^2
=
\|\phi(z)-\phi(z')\|_2^2.
\tag{8.4}
\]

Highly correlated shell states may have exponentially large
cardinality but small Rademacher width; in that case generic chaining,
rather than scalar entropy, could still produce a valid replacement.
Conversely, failure of both union-bound traffic and chaining forces a
large separated feature packing.  Combining such a packing with the
flat low-variance direction (5.5) is the remaining restricted
invertibility/polar-volume problem.

---

## 9. What has and has not been proved

The following implication is rigorous:

\[
\begin{aligned}
d\ge\varepsilon n^{3/2},\quad
\mathcal L_r\ne\varnothing
\quad\Longrightarrow\quad
&
\operatorname{tr}R^2-n
\ge(4(\varepsilon-2r)_+^2-o(1))n,\\
&
h^\top\Sigma_rh\le4r^2n^3.
\end{aligned}
\tag{9.1}
\]

At \(r=0\), \(h\) is an exact flat kernel vector of the internal
feature covariance and of every block Schur complement.

What is **not** proved is

\[
\log|\mathcal L_r|
<
\frac5{64}r^2n.
\tag{9.2}
\]

Nor has (7.3) yet been converted into a concrete sign replacement
contradicting (8.3).  Standard covariance entropy, Hanson--Wright,
one-replica log determinants, and bare dimension counting all lose a
fixed amount at the exponent scale and cannot close (9.2).

The sharpened next lemma is:

> **Flat conditional-kernel inverse theorem.**  
> For a cut-feature cloud satisfying the triangle identities, if a
> dense \(\{\pm1\}\) internal vector \(h\) lies in the approximate
> kernel of every relevant conditional covariance as in (7.2), then
> either the cloud has traffic entropy below (9.2), or there is an
> induced vertex block and a sign replacement \(\beta\) satisfying
> (8.2) with \(e=o(n^{3/2})\).

This formulation uses all information currently verified: the cut
triangle structure, the flatness and disjoint support of \(c,h\), the
two-sided cap identity, and global block replacement.  It is strictly
stronger than an ordinary covariance or spectral statement.
