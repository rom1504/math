# Biased margin channels: principal escape and paired-star structure

Checkpoint: 2026-07-26.

## 1. Scope

Let \(A\) be a symmetric zero-diagonal sign matrix of order \(n\)
satisfying

\[
A^2=(n-1)I.
\tag{1.1}
\]

Suppose its vertices are partitioned into exact affine positive-ground
types \(V_\alpha\), of sizes \(k_\alpha\), every internal type block
is all positive, and every cross-type block has total zero.  Write

\[
d_\alpha=k_\alpha-1,\qquad
u_\alpha=k_\alpha^{-1/2}\mathbf 1_{V_\alpha},
\qquad
r_{j,\alpha}=\sum_{i\in V_\alpha}a_{ji}.
\tag{1.2}
\]

The exact affine closure theorem gives, for \(j\in V_\beta\),

\[
\sum_{\alpha\ne\beta}|r_{j,\alpha}|\le d_\beta.
\tag{1.3}
\]

This note attacks the inverse problem left by the quantitatively
biased margin mass in Section 10 of
`margin_preserving_affine_refill.md`.

The main new outputs are:

1. the complete margin-channel system is an orthogonal family and
   each channel spans an exact invariant two-plane with its type
   indicator;
2. high normalized bias has a rigorous sequence-level dichotomy:
   it either escapes through \(o(\sqrt n)\)-sized source types,
   whose \(o(n)\) vertices can be removed by an asymptotically
   scale-preserving principal restriction, or a positive part of
   the bias lies on mesoscopic source types and has bounded incidence
   depth;
3. at the extremal compensated constant
   \(c=1/(2\sqrt2)\), the latter structure sharpens to an
   asymptotic depth-one paired-star quotient on all but \(o(n)\)
   vertices;
4. the 2026 bounded-\(\gamma_2\) inverse theorem applies to the
   threshold incidence layer, but its guaranteed homogeneous
   rectangle is necessarily an all-zero rectangle.  It does not by
   itself yield the desired positive Ferrers block.

The depth-one endpoint theorem is a genuine quotient statement.  At
a fixed moderate sub-half constant, bounded depth is obtained, but
global Ferrers nesting remains unproved.

## 2. Exact orthogonal margin channels

Let \(U\) be the matrix with columns \(u_\alpha\), and define

\[
g_\alpha=(I-UU^\mathsf T)Au_\alpha.
\tag{2.1}
\]

### Lemma 2.1 (channel coordinates and orthogonality)

For every occupied type \(\alpha\),

\[
Au_\alpha=d_\alpha u_\alpha+g_\alpha,
\tag{2.2}
\]

where

\[
g_\alpha(j)=
\begin{cases}
0,&j\in V_\alpha,\\[2mm]
r_{j,\alpha}/\sqrt{k_\alpha},&j\notin V_\alpha.
\end{cases}
\tag{2.3}
\]

Moreover,

\[
\boxed{
\langle g_\alpha,g_\beta\rangle=0\quad(\alpha\ne\beta),
\qquad
\|g_\alpha\|_2^2=n-1-d_\alpha^2,
}
\tag{2.4}
\]

and

\[
\boxed{
Ag_\alpha=
\bigl(n-1-d_\alpha^2\bigr)u_\alpha-d_\alpha g_\alpha.
}
\tag{2.5}
\]

Thus the mutually orthogonal planes
\(\operatorname{span}\{u_\alpha,g_\alpha\}\) are \(A\)-invariant
(when \(g_\alpha\ne0\)), and in the orthonormal basis
\((u_\alpha,g_\alpha/\|g_\alpha\|)\), the restriction of \(A\) is

\[
\begin{pmatrix}
d_\alpha&\sqrt{n-1-d_\alpha^2}\\
\sqrt{n-1-d_\alpha^2}&-d_\alpha
\end{pmatrix}.
\tag{2.6}
\]

#### Proof

The compression \(U^\mathsf TAU\) is diagonal: its diagonal entries
are \(d_\alpha\), since the internal block is all positive, and its
off-diagonal entries are zero, since every cross block has total
zero.  This gives (2.2)--(2.3).

Taking inner products of \(A^2u_\alpha=(n-1)u_\alpha\) with
\(u_\beta\) gives

\[
\langle Au_\alpha,Au_\beta\rangle
=(n-1)\mathbf1_{\{\alpha=\beta\}}.
\]

Substitute (2.2) and use \(g_\alpha\perp\operatorname{span}U\).
This proves (2.4).  Finally,

\[
A(d_\alpha u_\alpha+g_\alpha)
=(n-1)u_\alpha,
\]

and (2.2) gives (2.5).  Equation (2.6) follows. \(\square\)

In unnormalized margin coordinates, (2.4) is

\[
\boxed{
\sum_jr_{j,\alpha}r_{j,\beta}=0\quad(\alpha\ne\beta),
\qquad
\sum_jr_{j,\alpha}^2
=k_\alpha\bigl[n-1-(k_\alpha-1)^2\bigr].
}
\tag{2.7}
\]

Also \(d_\alpha\le\sqrt{n-1}\), so

\[
\boxed{k_\alpha\le1+\sqrt{n-1}.}
\tag{2.8}
\]

This is the exact signed-incidence information that a descent or
quotient theorem must use.  The scalar \(L^1\)-budget alone misses
both (2.7) and the invariant-plane relation (2.5).

The uncompressed block equations are also exact.  If
\(B_{\alpha\beta}=A[V_\alpha,V_\beta]\), then

\[
\boxed{
\sum_{\beta\ne\alpha}
B_{\alpha\beta}B_{\alpha\beta}^{\mathsf T}
=(n-2)I_{k_\alpha}-(k_\alpha-2)J_{k_\alpha},
}
\tag{2.9}
\]

and, for \(\alpha\ne\beta\),

\[
\boxed{
(J-I)B_{\alpha\beta}
+B_{\alpha\beta}(J-I)
+\sum_{\gamma\ne\alpha,\beta}
B_{\alpha\gamma}B_{\gamma\beta}=0.
}
\tag{2.10}
\]

These follow by taking the corresponding diagonal and off-diagonal
blocks of \(A^2=(n-1)I\).  In particular, arbitrary bounded-degree
threshold supports are not enough: their signs and column defects
must simultaneously satisfy the tight-frame identity (2.9) and the
three-type compatibility equation (2.10).

## 3. Principal monotonicity

### Lemma 3.1

For every principal set \(S\subseteq[n]\),

\[
p(A[S])\le p(A),\qquad
\nu(A[S])\le\nu(A),\qquad
W(A[S])\le W(A).
\tag{3.1}
\]

#### Proof

Fix a spin vector on \(S\), and extend it by independent uniform
spins on \(S^c\).  The expected full energy is exactly the energy
on \(S\).  Hence the full maximum is at least every energy on \(S\),
and the full minimum is at most every energy on \(S\).  Maximize and
minimize on \(S\), respectively. \(\square\)

Consequently, if \(|S|=n-o(n)\) and \(W(A)=O(n^{3/2})\), then

\[
\boxed{
\frac{W(A[S])}{|S|^{3/2}}
\le
\frac{W(A)}{n^{3/2}}+o(1).
}
\tag{3.2}
\]

When \(S\) is a union of complete affine types, the exact
positive-ground closure is inherited as well.

## 4. Threshold incidence and its exact paired-star form

Fix \(0<\theta<1\).  Define the directed threshold incidence set

\[
\mathcal E_\theta
=\{(j,\alpha):j\notin V_\alpha,\ 
|r_{j,\alpha}|\ge\theta k_\alpha\}.
\tag{4.1}
\]

For \((j,\alpha)\in\mathcal E_\theta\), put
\(\sigma_{j,\alpha}=\operatorname{sgn}r_{j,\alpha}\).  The column
\((a_{ij})_{i\in V_\alpha}\) then has the exact representation

\[
a_{ij}=
\begin{cases}
\sigma_{j,\alpha},&i\in V_\alpha\setminus D_{j,\alpha},\\
-\sigma_{j,\alpha},&i\in D_{j,\alpha},
\end{cases}
\qquad
|D_{j,\alpha}|
=\frac{k_\alpha-|r_{j,\alpha}|}{2}
\le\frac{1-\theta}{2}k_\alpha.
\tag{4.2}
\]

Equivalently, pair every minority edge with a majority edge.  The
pairs cancel in the column sum, leaving exactly
\(|r_{j,\alpha}|\) unpaired edges, all of sign
\(\sigma_{j,\alpha}\).  This is an exact local paired-star
decomposition; no regularity or probability is used.

### Lemma 4.1 (bounded depth on mesoscopic sources)

For \(\kappa>0\), restrict (4.1) to source types satisfying

\[
k_\alpha\ge\kappa\sqrt n.
\tag{4.3}
\]

Then every target vertex belongs to at most

\[
\boxed{
D_{\theta,\kappa}
=\left\lceil\frac{2}{\theta\kappa}\right\rceil
}
\tag{4.4}
\]

threshold incidences, for all sufficiently large \(n\).

#### Proof

Each such incidence consumes at least
\(\theta\kappa\sqrt n\) from the left side of (1.3), while
\(d_\beta\le\sqrt{n-1}\) by (2.8). \(\square\)

Thus the mesoscopic threshold layer is a bounded-depth signed
incidence quotient: every target vertex is attached to at most
\(D_{\theta,\kappa}\) source types, and every attachment is the
paired-star column (4.2).  This is weaker than common Ferrers
nesting because the defect sets \(D_{j,\alpha}\) may depend
arbitrarily on \(j\).

## 5. Escape-to-zero versus bounded-depth mass

Assume

\[
p(A)=c\,n^{3/2}+o(n^{3/2}),
\qquad
\frac1{2\sqrt2}\le c<\frac12.
\tag{5.1}
\]

Put \(q=y\sqrt n+o(\sqrt n)\).  Exact flatness and the affine
compression give

\[
\frac1{2c}\le y\le4c,
\tag{5.2}
\]

because \(q\sum k_\alpha^2\ge n^2\) and

\[
2\sum k_\alpha^2\ge q(n-2)+3n.
\]

In particular \(q=\Theta(\sqrt n)\).

Let

\[
\gamma(c)=\frac1{4c^2}-1,\qquad
\theta=\frac{\gamma(c)}2.
\tag{5.3}
\]

The compensated-margin calculation gives

\[
\sum_\alpha\eta_\alpha
=(y-2c)n^{3/2}+o(n^{3/2})
\ge
\left(\frac1{2c}-2c-o(1)\right)n^{3/2},
\tag{5.4}
\]

and the \(L^1\)-weighted mean normalized bias is at least
\(\gamma(c)-o(1)\).  Therefore

\[
\boxed{
\sum_{(j,\alpha)\in\mathcal E_\theta}|r_{j,\alpha}|
\ge b(c)n^{3/2}
}
\tag{5.5}
\]

for some \(b(c)>0\) and all sufficiently large \(n\).  One may take

\[
b(c)
=\frac{\gamma(c)}{2-\gamma(c)}
\left(\frac1{2c}-2c\right)/2
\tag{5.6}
\]

after absorbing the asymptotic error.

### Theorem 5.1 (principal escape or bounded-depth quotient)

Along every sequence satisfying (5.1), one of the following occurs
after passage to a subsequence.

1. **Escape through small source types.**  There is
   \(\kappa_n\downarrow0\) such that all but \(o(n^{3/2})\) of the
   mass in (5.5) comes from source types
   \(k_\alpha<\kappa_n\sqrt n\).  Their union \(T_n\) has

   \[
   |T_n|\le q\kappa_n\sqrt n=o(n).
   \tag{5.7}
   \]

   The principal signing on \(S_n=[n]\setminus T_n\) inherits exact
   affine positive-ground closure and satisfies

   \[
   \boxed{
   \frac{W(A[S_n])}{|S_n|^{3/2}}
   \le\frac{W(A)}{n^{3/2}}+o(1).
   }
   \tag{5.8}
   \]

2. **Mesoscopic bias.**  There are fixed
   \(\kappa,\varepsilon>0\) such that the part of (5.5) with
   \(k_\alpha\ge\kappa\sqrt n\) has mass at least
   \(\varepsilon n^{3/2}\).  It contains at least
   \((\varepsilon/2)n\) incidences, every target vertex has depth at
   most \(D_{\theta,\kappa}\), and each incidence has the exact
   paired-star form (4.2).

#### Proof

Apply the standard tightness dichotomy at zero to the finite measures
on normalized source sizes \(k_\alpha/\sqrt n\), weighted by
\(|r_{j,\alpha}|/n^{3/2}\) on \(\mathcal E_\theta\).

In the escaping branch choose \(\kappa_n\downarrow0\) by a diagonal
argument.  Equation (5.2) gives (5.7), and Lemma 3.1 gives (5.8).

In the tight branch, (2.8) implies
\(|r_{j,\alpha}|\le(1+o(1))\sqrt n\), so mass
\(\varepsilon n^{3/2}\) requires at least
\((\varepsilon/2)n\) incidences.  Lemma 4.1 and (4.2) finish the
proof. \(\square\)

The first branch is a genuine asymptotically scale-preserving
principal restriction, not merely deletion of a covariance model.
It also removes the margin channels responsible for the biased mass.
The second branch is the promised bounded-depth paired quotient,
although at moderate \(\theta\) it need not be globally Ferrers.

## 6. Endpoint rigidity: depth one at \(1/(2\sqrt2)\)

The compensated lower endpoint has much stronger equality
information.

### Theorem 6.1 (asymptotic depth-one paired quotient)

Assume

\[
p(A)=\left(\frac1{2\sqrt2}+o(1)\right)n^{3/2}.
\tag{6.1}
\]

Then, after deleting \(o(n)\) exceptional vertices (including all
vertices in size-exceptional types) and discarding
\(o(n^{3/2})\) absolute margin mass, the following hold.

1. Every remaining type has size

   \[
   k_\alpha=(1+o(1))\sqrt{n/2}.
   \tag{6.2}
   \]

2. Every retained nonzero margin satisfies

   \[
   |r_{j,\alpha}|=(1-o(1))k_\alpha.
   \tag{6.3}
   \]

3. Every remaining target vertex \(j\) has at most one retained
   source type \(\alpha\), and all but \(o(n)\) remaining target
   vertices have one.

4. The column from \(j\) to that source type is
   \(o(k_\alpha)\)-Hamming-close to a constant signed column.
   For each ordered type pair \((\beta,\alpha)\), the positive and
   negative such columns can be paired except for total
   \(o(n)\) exceptional columns (in the margin-weighted sense).

Hence the extremal compensated branch is an asymptotic depth-one
paired-star/Ferrers quotient.

#### Proof

Let \(S_2=\sum k_\alpha^2\).  From (6.1),

\[
S_2=\left(\frac1{\sqrt2}+o(1)\right)n^{3/2}.
\]

Cauchy and the exact flat inequality squeeze

\[
q=(\sqrt2+o(1))\sqrt n,\qquad
\mu:=n/q=(1+o(1))\sqrt{n/2},
\tag{6.4}
\]

and

\[
\sum_\alpha(k_\alpha-\mu)^2
=S_2-\frac{n^2}{q}=o(n^{3/2}).
\tag{6.5}
\]

Thus, for every fixed \(\varepsilon>0\), all but \(o(n)\) vertices
belong to types with
\((1-\varepsilon)\mu\le k_\alpha\le(1+\varepsilon)\mu\).

Next, the exact compensation identity gives

\[
\sum_\alpha\eta_\alpha
=\left(\frac1{\sqrt2}+o(1)\right)n^{3/2}
=2p(A)+o(n^{3/2}).
\]

Since

\[
\sum_\alpha\eta_\alpha\le L\le2p(A),
\]

we have

\[
\boxed{
L-\sum_\alpha\eta_\alpha
=\sum_{j,\alpha}|r_{j,\alpha}|
\left(1-\frac{|r_{j,\alpha}|}{k_\alpha}\right)
=o(n^{3/2}).
}
\tag{6.6}
\]

For any fixed \(\varepsilon>0\), (6.6) discards only
\(o(n^{3/2})\) mass when we require
\(|r_{j,\alpha}|\ge(1-\varepsilon)k_\alpha\).
Channels belonging to size-exceptional source types also have
negligible mass: their total is at most
\(\sum_{\alpha\ {\rm exceptional}}\eta_\alpha+o(n^{3/2})\);
each \(\eta_\alpha\le n\), while (6.5) leaves only
\(o(\sqrt n)\) exceptional source types.
Exceptional target vertices have capacity at most \(\sqrt n\), so
their mass is negligible as well.

Choose \(\varepsilon< (5-\sqrt{17})/4\).  Two retained incidences at
one target vertex would consume more than its domination budget,
because

\[
2(1-\varepsilon)^2\mu>(1+\varepsilon)\mu.
\]

This proves depth one.  Formula (4.2) gives the
\(o(k_\alpha)\)-defect claim after letting
\(\varepsilon\downarrow0\) diagonally.

The sum of all target domination budgets is \(2p(A)\), while the
retained margin mass is
\(L-o(n^{3/2})=2p(A)-o(n^{3/2})\).  Every good target has budget
\((1+o(1))\mu\).  Hence only \(o(n)\) good target vertices can have
no retained incidence.

Finally, every cross block has total zero:
\(\sum_{j\in V_\beta}r_{j,\alpha}=0\).  After the discarded margin
mass is removed, retained margins have a common magnitude
\((1+o(1))k_\alpha\).  Thus their positive and negative counts in
each ordered block balance up to the discarded weighted error;
pair them greedily. \(\square\)

## 7. The invariant plane produces half-hard principal faces

The depth-one quotient does not automatically give a favorable
principal descent.  In fact, the exact invariant-plane equation
forces the most natural half-size restrictions to sit at the
spectral half scale.

### Theorem 7.1 (endpoint crossover faces)

Under the hypotheses of Theorem 6.1, there exist Boolean vectors
\(x,y\in\{\pm1\}^n\) such that

\[
\langle x,y\rangle=o(n),
\tag{7.1}
\]

\[
x^\mathsf TAx
=\left(\frac1{\sqrt2}+o(1)\right)n^{3/2},
\qquad
y^\mathsf TAy
=-\left(\frac1{\sqrt2}+o(1)\right)n^{3/2},
\tag{7.2}
\]

and

\[
x^\mathsf TAy
=\left(\frac1{\sqrt2}+o(1)\right)n^{3/2}.
\tag{7.3}
\]

If

\[
I=\{j:x_j=y_j\},\qquad
J=\{j:x_j=-y_j\},
\tag{7.4}
\]

then

\[
|I|=|J|=n/2+o(n)
\tag{7.5}
\]

and the two principal signings obey

\[
\boxed{
\frac{M(A[I])}{|I|^{3/2}}\ge\frac12-o(1),
\qquad
\frac{M(A[J])}{|J|^{3/2}}\ge\frac12-o(1).
}
\tag{7.6}
\]

#### Proof

Choose arbitrary type signs \(t_\alpha\), and put

\[
x=\sum_\alpha\sqrt{k_\alpha}\,t_\alpha u_\alpha.
\tag{7.7}
\]

This is a type-constant Boolean vector.  Normalize the margin
channels by

\[
h_\alpha
=\frac{g_\alpha}{s_\alpha},
\qquad
s_\alpha=\sqrt{n-1-d_\alpha^2},
\tag{7.8}
\]

and define the real vector

\[
y_0=\sum_\alpha\sqrt{k_\alpha}\,t_\alpha h_\alpha.
\tag{7.9}
\]

Theorem 6.1 implies that, outside \(o(n)\) coordinates,
\(h_\alpha\) has disjoint depth-one support and its nonzero
coordinates are

\[
(1+o(1))\sigma_{j,\alpha}/\sqrt{k_\alpha}.
\]

Round \(y_0\) coordinatewise on the good support, using precisely
the corresponding signs, and fill the exceptional coordinates
arbitrarily.  Equation (6.6), the size concentration (6.5), and
the fact that all retained margins are
\((1-o(1))k_\alpha\) give

\[
\|y-y_0\|_2=o(\sqrt n).
\tag{7.10}
\]

For completeness, the normalization behind (7.10) is

\[
\sqrt{k_\alpha}\,h_\alpha(j)
=\frac{r_{j,\alpha}}{s_\alpha},
\qquad
s_\alpha=(1+o(1))\sqrt{n/2}.
\]

On a good target coordinate there is one retained term, equal to
\((1+o(1))\sigma_{j,\alpha}\), while the sum of all discarded
terms has total \(\ell_1\)-mass \(o(n)\) after division by
\(\sqrt{n/2}\).  Row domination bounds the coordinatewise error by
\(O(1)\), so its squared \(\ell_2\)-mass is \(o(n)\).  The
size-exceptional source channels contribute
\(\sum_{\alpha\ {\rm exceptional}}k_\alpha=o(n)\) to
\(\|y_0\|_2^2\); exceptional targets are controlled in the same
way by their \(o(n^{3/2})\) margin capacity.  This proves (7.10)
without a hidden factor of \(n^{1/4}\).

The columns \(h_\alpha\) are orthonormal and perpendicular to every
\(u_\beta\), so

\[
\langle x,y_0\rangle=0.
\]

Equation (7.10) proves (7.1).

By (2.6),

\[
x^\mathsf TAx=\sum_\alpha k_\alpha d_\alpha,
\qquad
y_0^\mathsf TAy_0=-\sum_\alpha k_\alpha d_\alpha,
\tag{7.11}
\]

and

\[
x^\mathsf TAy_0=\sum_\alpha k_\alpha s_\alpha.
\tag{7.12}
\]

At the endpoint, (6.2) gives

\[
d_\alpha=s_\alpha
=(1+o(1))\sqrt{n/2}
\]

off a negligible set of types.  Equations (7.10)--(7.12), together
with \(\|A\|_{\rm op}=\sqrt{n-1}\), prove (7.2)--(7.3).

Put

\[
u=\frac{x+y}{2},\qquad v=\frac{x-y}{2}.
\]

These are signed indicator vectors supported on \(I,J\),
respectively.  Equations (7.1)--(7.3) give

\[
u^\mathsf TAu
=\left(\frac1{2\sqrt2}+o(1)\right)n^{3/2},
\qquad
v^\mathsf TAv
=-\left(\frac1{2\sqrt2}+o(1)\right)n^{3/2}.
\tag{7.13}
\]

Since one-copy energy is half the doubled quadratic form and
\(|I|,|J|=n/2+o(n)\), (7.6) follows. \(\square\)

Thus the endpoint paired quotient is not merely a proof artifact:
its crossover faces are themselves asymptotically half-hard.
Any favorable descent must use a different, more global selection
than the agreement/disagreement faces of one invariant channel
cube.

## 8. Audit of the bounded-\(\gamma_2\) inverse theorem

Balla--Hambardzumyan--Tomon (Math. Ann. 394 (2026), article 52,
DOI `10.1007/s00208-026-03355-2`) prove that a Boolean matrix of
bounded \(\gamma_2\)-norm contains a linear-by-linear homogeneous
submatrix.

Apply this to the \(0/1\) incidence matrix \(M_{\theta,\kappa}\)
whose rows are target vertices and whose columns are mesoscopic
source types, with a one on \(\mathcal E_\theta\).  Lemma 4.1 says
every row has at most \(D_{\theta,\kappa}\) ones.  The associated
bipartite graph is \(D_{\theta,\kappa}\)-degenerate, and the same
paper's degeneracy estimate gives

\[
\gamma_2(M_{\theta,\kappa})
\le2\sqrt{D_{\theta,\kappa}}.
\tag{8.1}
\]

The inverse theorem therefore yields a
\(\delta n\times\delta q\) homogeneous rectangle, with
\(\delta>0\) depending only on \(\theta,\kappa\).  For large \(n\)
it cannot be an all-one rectangle, because each row contains at most
\(D_{\theta,\kappa}\) ones while \(\delta q\to\infty\).
Consequently the theorem yields

\[
\boxed{
\text{a linear set of target vertices with no high-bias incidence
to a linear set of mesoscopic source types.}
}
\tag{8.2}
\]

This is a useful macroscopic low-bias rectangle, but it is the
opposite polarity from the positive Ferrers block needed for direct
quotient descent.  The normalization therefore closes, but the
conclusion does not finish the target.

## 9. Sharp remaining obstruction

At a moderate fixed \(c\in(1/(2\sqrt2),1/2)\), Theorem 5.1 gives a
linear-sized bounded-depth paired-star layer, and (2.7) supplies
exact orthogonality of the complete signed margin columns.  What is
still missing is a theorem upgrading these **local** paired stars
to common Ferrers nesting or to a principal restriction of a fixed
smaller proportion with no normalized-width loss.

The gap is real at the level of threshold incidence: bounded row
degree permits arbitrary sparse expanders, and the bounded
\(\gamma_2\) inverse theorem can force only a large zero rectangle.
Therefore any further upgrade must use the signs and exact
orthogonality in (2.7), or the invariant-plane equation (2.5), not
only the support graph and row-domination.
