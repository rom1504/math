# Deficit-kernel inverse campaign

Checkpoint date: 2026-07-26.

## 1. Setup

Let \(\mathcal F\subseteq\binom{[n]}m\) be a family of positive-good
sets for a symmetric zero-diagonal signing \(A\).  For \(S,T\in
\mathcal F\), write

\[
X=S\cap T,\qquad Y=S\setminus T,\qquad Z=T\setminus S
\]

and

\[
t(S,T)=e(X,Y)=e(X,Z)=-e(Y,Z)\geq 0.
\]

The purpose of this note is to isolate what the low-rank/Fourier
structure of the kernel \(t\) proves, and exactly where zero-deficit
pairs remain an obstruction.

## 2. The deficit kernel has polynomial rank

### Proposition 2.1

The matrix

\[
K=(t(S,T))_{S,T\in\mathcal F}
\]

has

\[
\boxed{\operatorname{rank}K\leq n+\binom n2.}
\tag{2.1}
\]

#### Proof

For arbitrary \(S,T\), put \(s_i=\mathbf1_{\{i\in S\}}\) and
\(u_i=\mathbf1_{\{i\in T\}}\), and define the symmetrized expression

\[
\widetilde t(S,T)
=\frac12\bigl(e(S\cap T,S\setminus T)
e(S\cap T,T\setminus S)\bigr).
\]

The crossing law says \(\widetilde t=t\) on
\(\mathcal F\times\mathcal F\).  Expanding with an ordered sum gives

\[
\widetilde t(S,T)
=\frac12\sum_{i,j}a_{ij}
\left(s_is_j u_i+s_i u_i u_j-2s_is_j u_i u_j\right).
\tag{2.2}
\]

Every term in (2.2) is bilinear in the feature vector

\[
\Phi(S)=\left((s_i)_{i\in[n]},(s_is_j)_{i<j}\right)
\in\mathbb R^{\,n+\binom n2}
\]

and the corresponding feature vector \(\Phi(T)\).  Hence
\(K=\Phi^\mathsf TB\Phi\) for one fixed symmetric matrix \(B\), and
(2.1) follows. \(\square\)

This is an association-scheme-type reduction which does not require
\(\mathcal F\) itself to be distance regular.

## 3. A finite-value polynomial bound

### Proposition 3.1

Let \(\mathcal G\subseteq\mathcal F\) have the property that

\[
1\leq t(S,T)\leq L\qquad(S\neq T,\ S,T\in\mathcal G).
\tag{3.1}
\]

Put \(r=n+\binom n2\).  Then

\[
\boxed{|\mathcal G|\leq\binom{r+L}{L}.}
\tag{3.2}
\]

#### Proof

Factor the symmetric rank-\(r\) kernel as

\[
t(S,T)=v_S^\mathsf TJv_T
\]

in \(\mathbb R^r\), with \(J\) a fixed symmetric bilinear form.
For each \(S\in\mathcal G\), consider the polynomial

\[
P_S(x)=\prod_{\ell=1}^L(v_S^\mathsf TJx-\ell).
\]

Since \(t(S,S)=0\),

\[
P_S(v_S)=(-1)^L L!\neq0,
\]

whereas (3.1) gives \(P_S(v_T)=0\) for \(T\neq S\).
The polynomials \(P_S\), restricted to the finite set
\(\{v_T:T\in\mathcal G\}\), are therefore linearly independent.
They lie in the space of polynomials of total degree at most \(L\)
in \(r\) variables, whose dimension is \(\binom{r+L}{L}\).
\(\square\)

At the balanced incidence threshold
\(|\mathcal G|=2^{n/2-o(n)}\), (3.2) implies

\[
\boxed{\max_{S\neq T}t(S,T)
\geq \left(\frac{\log 2}{2}+o(1)\right)\frac n{\log n}}
\tag{3.3}
\]

provided every off-diagonal deficit is positive.  Indeed, for
\(L=o(r)\),

\[
\log\binom{r+L}{L}
\leq L\log\frac{e(r+L)}L,
\]

and \(r=(1/2+o(1))n^2\).

This is rigorous but far below the tangent scale \(n^{3/2}\).
More importantly, the proof fails exactly when \(t(S,T)=0\):
the isolating polynomial takes the same nonzero value at a
zero-deficit neighbour as on the diagonal.  Thus the low rank of
\(K\) by itself does not turn a dense zero-deficit graph into an
affine family.

## 4. The pairwise-pseudorandom branch already gives a tangent witness

Let \(S\) be uniform on \(\mathcal F\), and put

\[
p_i=\Pr(i\in S),\qquad p_{ij}=\Pr(i,j\in S),\qquad
d_{ij}=p_i+p_j-2p_{ij}.
\]

The exact mean identity from the main inverse note is

\[
\mathbb E_{S,T}t(S,T)
=\sum_{i<j}a_{ij}p_{ij}d_{ij}.
\tag{4.1}
\]

Every member is an absolute positive child ground, so

\[
\overline E:=\sum_{i<j}a_{ij}p_{ij}
=\mathbb E H_{A[S]}(\mathbf1)
\geq M_m.
\tag{4.2}
\]

### Proposition 4.1

For every real \(\delta\),

\[
\boxed{
\mathbb E t(S,T)
\geq \delta\,\overline E
-\sum_{i<j}p_{ij}|d_{ij}-\delta|.
}
\tag{4.3}
\]

Consequently, if \(\delta\geq\delta_0>0\) and

\[
\sum_{i<j}p_{ij}|d_{ij}-\delta|=o(n^{3/2}),
\tag{4.4}
\]

then some pair \(S,T\) has

\[
t(S,T)\geq \delta_0M_m-o(n^{3/2})
=\Omega(n^{3/2}).
\tag{4.5}
\]

#### Proof

Subtract \(\delta\overline E\) from (4.1), and use
\(|a_{ij}|=1\).  Since the mean of the nonnegative kernel is at
most its maximum, (4.5) follows. \(\square\)

Thus a pairwise-pseudorandom incidence family is not the difficult
case: it immediately creates a tangent-scale principal witness.
Any counterexample to the desired inverse theorem must correlate the
signed child energy very strongly with the separation probabilities
\(d_{ij}\), or must put substantial mass on exact zero-deficit pairs.

## 5. Precise wall and next target

The polynomial-rank argument proves that a large clique of uniformly
positive, low integer deficits is impossible.  It does **not** prove
that a large family contains a large zero-deficit clique: a dense
zero graph need not have a large clique, and the rank factorization is
indefinite (so Euclidean-distance or PSD arguments are unavailable).

The next useful statement would be one of the following.

1. A same-\(A\) theorem showing that the zero-deficit graph of
positive absolute child grounds has a large clique or many additive
rectangles.
2. A weighted decomposition proving either (4.4) on a large
subfamily or a coordinate/junta descent from its failure.
3. A polynomial which also separates diagonal pairs from
zero-deficit off-diagonal pairs by using one of the strict
complementarity fields, rather than the scalar kernel \(t\) alone.

No Balog--Szemerédi--Gowers conclusion follows from the scalar
condition \(t=0\) without such an additional bridge.

## 6. Tangent deficit or a large fixed-intersection fibre

There is a second exact reduction which uses the **absolute** child
ground condition and therefore survives the boundary-only
counterexamples.

### Proposition 6.1

Assume \(m\asymp n\) and \(M(A)\leq Cn^{3/2}\).  Fix
\(S\in\mathcal F\), and put

\[
p_S=H_{A[S]}(\mathbf1)=M(A[S]).
\]

There is a constant \(c=c(C,m/n)>0\) such that either

\[
\boxed{\max_{T\in\mathcal F}t(S,T)\geq p_S/4}
\tag{6.1}
\]

or some \(X\subseteq S\) has

\[
\boxed{
\left|\{T\in\mathcal F:T\cap S=X\}\right|
\geq
\frac{|\mathcal F|}{2^m\exp(-c n^{3/4})}.
}
\tag{6.2}
\]

In particular, at \(m=n/2\) and
\(|\mathcal F|\geq2^{-m}\binom nm=2^{m-o(n)}\), failure of (6.1)
produces a fixed-intersection fibre of size

\[
\boxed{\exp(\Omega(n^{3/4})).}
\tag{6.3}
\]

#### Proof

Let \(B=A[S]\).  For \(Y=S\setminus T\), flipping precisely the
vertices of \(Y\) inside the all-one child gives

\[
H_B(\mathbf1_{S\setminus Y},-\mathbf1_Y)
=p_S-2e(Y,S\setminus Y)
=p_S-2t(S,T).
\tag{6.4}
\]

If (6.1) fails, every such spin has energy at least \(p_S/2\).

We use the following elementary spectral bootstrap.

### Lemma 6.2 (quadratic norm versus spectral norm)

For every symmetric zero-diagonal signing \(B\),

\[
\boxed{
Q(B):=\max_{x\in\{\pm1\}^m}|x^\mathsf TBx|
\geq \frac12\|B\|_{\rm op}(\|B\|_{\rm op}+1).
}
\tag{6.5}
\]

#### Proof

Let \(Bv=\lambda v\), where
\(|\lambda|=\|B\|_{\rm op}\), \(\|v\|_2=1\), and put
\(\alpha=\|v\|_\infty\).  At a coordinate where
\(|v_i|=\alpha\), the zero diagonal and \(|b_{ij}|=1\) give

\[
|\lambda|\alpha
=\left|\sum_{j\ne i}b_{ij}v_j\right|
\leq\|v\|_1-\alpha.
\]

Thus \(\|v\|_1/\alpha\geq|\lambda|+1\).  For
\(z=v/\alpha\in[-1,1]^m\),

\[
\|B\|_{\infty\to1}
\geq\|Bz\|_1
=|\lambda|\frac{\|v\|_1}{\alpha}
\geq|\lambda|(|\lambda|+1).
\tag{6.6}
\]

For sign vectors \(x,y\), put
\(u=(x+y)/2\) and \(w=(x-y)/2\).  Polarization gives

\[
x^\mathsf TBy=u^\mathsf TBu-w^\mathsf TBw.
\]

The absolute value of a zero-diagonal quadratic form on
\([-1,1]^m\) is maximized at a cube vertex, so both terms have
absolute value at most \(Q(B)\).  Hence
\(\|B\|_{\infty\to1}\leq2Q(B)\), proving (6.5). \(\square\)

The principal-submatrix inequality \(M(B)\leq M(A)\) and Lemma 6.2
give

\[
\|B\|_{\rm op}=O(n^{3/4}).
\tag{6.7}
\]

This improves, in the present spectral range, the earlier bootstrap

\[
\|B\|_{\rm op}^3\leq m\,Q(B),
\qquad Q(B)=2M(B),
\tag{6.8}
\].

Also \(\|B\|_F^2=m(m-1)\).  The Rademacher Hanson--Wright inequality,
applied to \(x^\mathsf TBx=2H_B(x)\), therefore yields

\[
\Pr_x\!\left(H_B(x)\geq p_S/2\right)
\leq
2\exp\left[
-c_0\min\left\{
\frac{p_S^2}{\|B\|_F^2},
\frac{p_S}{\|B\|_{\rm op}}
\right\}\right].
\tag{6.9}
\]

The universal child lower bound \(p_S=M(B)\geq M_m=\Omega(m^{3/2})\)
makes the two exponents in (6.9) respectively
\(\Omega(n)\) and \(\Omega(n^{3/4})\).  Hence at most

\[
2^m\exp(-c n^{3/4})
\tag{6.10}
\]

different subsets \(Y\subseteq S\) can occur.  Since
\(Y=S\setminus(T\cap S)\), pigeonholing the members of
\(\mathcal F\) over these traces proves (6.2). \(\square\)

The fibre in (6.2) has additional exact structure.  Write its members
as

\[
T=X\cup Z,\qquad Z\subseteq V:=S^c,
\]

with \(X\) fixed.  The crossing law gives one fixed number
\[
t_0=e(X,S\setminus X)=e(X,Z)=-e(S\setminus X,Z)
\]
for every \(Z\) in the fibre.  If
\[
b=A_{V,X}\mathbf1_X,\qquad D=A[V],
\]
then strict complementarity becomes the affine system

\[
\boxed{
(D\mathbf1_Z+b)_v=0\quad(v\notin Z),\qquad
(D\mathbf1_Z+b)_v\geq1\quad(v\in Z).
}
\tag{6.11}
\]

Thus failure of a tangent deficit does not leave an arbitrary
code: it leaves a superpolynomial family of constant-weight solutions
of one common **affine** symmetric LCP.  The remaining gap is an
inverse theorem for (6.11) at entropy \(\Omega(n^{3/4})\), below the
linear-entropy affine-subspace theorem already proved.

## 7. Zero deficit is exact top-face factorization

The scalar equality \(t=0\) has a stronger consequence once the
global child-ground inequalities are retained.

### Proposition 7.1 (quantitative face splitting)

Let \(S,T\in\mathcal F\), write

\[
X=S\cap T,\qquad Y=S\setminus T,
\]

and put \(B=A[X]\), \(C=A[Y]\), \(R=A[X,Y]\).  Write

\[
f_X=H_B(\mathbf1),\quad f_Y=H_C(\mathbf1),\quad
p_X=\max_xH_B(x),\quad p_Y=\max_yH_C(y).
\]

Then

\[
\boxed{
0\leq p_X-f_X\leq t(S,T),\qquad
0\leq p_Y-f_Y\leq t(S,T).
}
\tag{7.1}
\]

Moreover, for every positive maximizer \(x\) of \(B\) and every
positive maximizer \(y\) of \(C\),

\[
\boxed{
|x^\mathsf TRy|
\leq t(S,T)-(p_X-f_X)-(p_Y-f_Y)
\leq t(S,T).
}
\tag{7.2}
\]

In particular, if \(t(S,T)=0\), then

\[
\boxed{
p_X=f_X,\qquad p_Y=f_Y,\qquad
x^\mathsf TRy=0
}
\tag{7.3}
\]

for every pair of positive component grounds.  Every concatenation
\((x,y)\) of component grounds is then a positive ground of
\(A[S]\).

#### Proof

The crossing law gives

\[
H_{A[S]}(\mathbf1)=f_X+f_Y+t.
\tag{7.4}
\]

For arbitrary \(x\in\{\pm1\}^X\), the two child configurations
\((x,\mathbf1_Y)\) and \((x,-\mathbf1_Y)\) have energies

\[
H_B(x)+f_Y\pm x^\mathsf TR\mathbf1.
\]

Their average is \(H_B(x)+f_Y\).  Since the all-one child is a
positive maximum of value (7.4), taking \(x\) to be a positive
ground of \(B\) gives

\[
p_X+f_Y\leq f_X+f_Y+t.
\]

This is the first inequality in (7.1); the second is symmetric.

Now take positive component grounds \(x,y\).  The two configurations
\((x,y)\) and \((x,-y)\) have energies

\[
p_X+p_Y\pm x^\mathsf TRy.
\]

Their larger value cannot exceed (7.4), which proves (7.2).
When \(t=0\), all inequalities are equalities and the cross term
vanishes, so every \((x,y)\) has energy \(p_X+p_Y\), equal to the
child maximum. \(\square\)

Thus a zero-deficit edge is not merely a scalar cut cancellation: it
is an exact Cartesian factorization of the positive top face, with
the cross block annihilating the two component ground clouds.
Low deficit is the corresponding quantitative near-factorization.
The remaining combinatorial problem is to turn many such compatible
pairwise factorizations into a common recursive partition.

## 8. Geometry of the affine strict-LCP fibre

The affine system (6.11) has an exact secant identity that is absent
from a generic sign-pattern arrangement.

### Proposition 8.1 (strict quadratic separation)

Let \(\eta,\theta\in\{0,1\}^V\) be two solutions of

\[
h^\eta=D\eta+b,\qquad
\operatorname{supp}h^\eta=\operatorname{supp}\eta,\qquad
h^\eta_i\geq1\quad(\eta_i=1),
\tag{8.1}
\]

and the analogous equations for \(\theta\).  Put
\[
P=\operatorname{supp}\eta\setminus\operatorname{supp}\theta,
\qquad
Q=\operatorname{supp}\theta\setminus\operatorname{supp}\eta.
\]
Then

\[
\boxed{
(\eta-\theta)^\mathsf TD(\eta-\theta)
=\sum_{i\in P}h^\eta_i+\sum_{i\in Q}h^\theta_i
\geq |P|+|Q|
=\|\eta-\theta\|_2^2.
}
\tag{8.2}
\]

#### Proof

The affine field cancels on taking differences:
\[
D(\eta-\theta)=h^\eta-h^\theta.
\]
On \(P\), the two fields are respectively at least one and zero; on
\(Q\), they are respectively zero and at least one.  Multiplying by
\(\eta-\theta\) and summing gives (8.2). \(\square\)

Thus every secant of the solution code lies in the positive cone of
\(D-I\).  In particular, a large fibre cannot be treated as an
arbitrary collection of cells of a hyperplane arrangement.

There is also an exact averaged form.  If \(\eta\) is uniform on the
fibre, let

\[
p_i=\Pr(\eta_i=1),\qquad
\mu_i=\mathbb E(h^\eta_i\mid\eta_i=1),\qquad
\Sigma=\operatorname{Cov}(\eta).
\]

For two independent solutions, averaging (8.2) before its final
inequality gives

\[
\boxed{
\operatorname{tr}(D\Sigma)
=\sum_i p_i(1-p_i)\mu_i
\geq\sum_i p_i(1-p_i)
=\operatorname{tr}\Sigma.
}
\tag{8.3}
\]

This identity exposes the remaining dichotomy.  Entropy carried by
coordinates with nontrivial variance consumes positive spectral mass
of \(D-I\); energy carried by coordinates with \(p_i\) extremely close
to one instead forms an almost-fixed principal core.  A quantitative
proof that one of these two resources is macroscopic would yield,
respectively, a spectral Boolean witness through Lemma 6.2 or a
principal-core descent.

The exponent \(3/4\) in Proposition 6.1 is natural for this
dichotomy: a competitive signing has
\(\|D\|_{\rm op}=O(n^{3/4})\), while a fibre of entropy
\(\Omega(n^{3/4})\) can in principle concentrate all of its
variability on that many coordinates.  Therefore an arrangement- or
rank-only argument cannot close the gap without also using the
component ground energies in Proposition 7.1.

## 9. Exact endpoint-face cover on \(n\equiv1\pmod4\)

The favorable-edge replacement theorem from the insertion campaign
interfaces especially cleanly with Proposition 7.1.  Switch a positive
absolute ground to \(\mathbf1\), and let

\[
c(R)=e(R,R^c),\qquad H_A(\mathbf1^R)=M-2c(R).
\]

For \(n\equiv1\pmod4\), every cut has even signed size and \(M\) is
even.  Therefore a deficit-two replacement witness is automatically
an **exact** endpoint witness:

\[
c(R)\leq1\Longrightarrow c(R)=0,\qquad
c(R)\geq M-1\Longrightarrow c(R)=M.
\tag{9.1}
\]

After switching a negative endpoint ground and replacing \(A\) by
\(-A\), the second alternative is again a zero cut of a positive
endpoint face.  Hence every favorable edge in the saturated opposite-
ground cross block is separated by a zero cut in one of the two exact
endpoint gauges.

### Proposition 9.1 (zero-cut Cartesian factorization)

If \(c(R)=0\), put \(B=A[R]\), \(D=A[R^c]\), and
\(C=A[R,R^c]\).  Then

\[
\boxed{
p(B)=H_B(\mathbf1),\qquad
p(D)=H_D(\mathbf1),
}
\tag{9.2}
\]

and for every positive ground \(x\) of \(B\) and \(y\) of \(D\),

\[
\boxed{x^\mathsf TCy=0.}
\tag{9.3}
\]

Every concatenation \((x,y)\) is a positive ground of \(A\).

This is Proposition 7.1 with \(S=[n]\) and \(t=c(R)=0\).  Thus the
edge cover is a cover by exact Cartesian factorizations of one of the
two endpoint faces, not merely by low-energy cuts.

There is also a cut-level uncrossing law which does not require the
rowwise traffic equations.

### Proposition 9.2 (endpoint zero-cut uncrossing)

If \(c(R)=c(Q)=0\), put

\[
P=R\setminus Q,\qquad W=Q\setminus R.
\]

Then

\[
\boxed{
t_0(R,Q):=-e(P,W)
=\frac{c(R\cap Q)+c(R\cup Q)}2\geq0.
}
\tag{9.4}
\]

If \(t_0(R,Q)=0\), both \(R\cap Q\) and \(R\cup Q\) are again exact
zero cuts.  If \(t_0>0\), then

\[
\boxed{M(A[P\cup W])\geq t_0(R,Q).}
\tag{9.5}
\]

Pairwise-disjoint symmetric-difference supports pack additively,
exactly as in Theorem 4.3 of the main inverse note.

#### Proof

The standard signed-cut identity is

\[
c(R)+c(Q)
=c(R\cap Q)+c(R\cup Q)+2e(P,W).
\]

All cut values are nonnegative because \(\mathbf1\) is a positive
ground, proving (9.4).  Equation (9.5) follows by comparing the two
relative sign choices on \(P,W\).  The disjoint packing proof uses
independent global signs on the disjoint supports and is unchanged.
\(\square\)

This gives an exact finite closure procedure for the insertion
branch:

1. repeatedly uncross pairs with \(t_0=0\), adjoining their unions
   and intersections, and adjoining complements (which are
   automatically zero cuts);
2. a positive \(t_0\) creates a principal deficit witness, and
   disjoint such witnesses add;
3. if closure never creates a positive deficit, it ends in a Boolean
   algebra of zero cuts.

In the last branch, let \(P_1,\ldots,P_k\) be the atoms of that Boolean
algebra.  Since every union of atoms has zero cut,

\[
\boxed{e(P_i,P_j)=0\quad(i\ne j).}
\tag{9.6}
\]

Indeed,
\[
c(P_i\cup P_j)=c(P_i)+c(P_j)-2e(P_i,P_j).
\]
Proposition 9.1 then iterates: the positive endpoint face contains the
Cartesian product of the positive component ground clouds, and all
cross blocks annihilate those clouds.

In fact this is exactly the affine-type closure already analyzed in
`affine_type_closure_recursion.md`.  Every union of atoms is a zero
cut, so flipping any union preserves the endpoint value \(M\).
Consequently **every atom-constant spin is a positive absolute ground**.
The atom partition therefore obeys all of the affine closure
conclusions: intertype block totals vanish, the positive maxima add
hereditarily, and the vertexwise row-domination inequalities hold.
Thus the zero-deficit closure branch is not a new exception; it lands
directly in the established macro-type versus mesoscopic-residual
dichotomy.  In particular, the square-root type scale is the genuine
spectral \(1/2\) wall of that branch.

The unresolved quantitative step is now precise.  Starting from an
\(\Omega(n^2)\)-edge incidence cover, prove that the closure procedure
either packs total \(t_0=\Omega(n^{3/2})\), or produces a Boolean
algebra and hence triggers the existing affine-type recursion.
Coverage alone does not yet force the former, because a first positive
uncrossing deficit can have weight \(O(1)\) on a support of linear
size.

## 10. Ground-span constraints on an exact factorization

Cut totals alone are insufficient for quantitative laminarization.
The additional same-face information can be recorded as a rank
constraint.

For a signing \(B\), let
\[
\mathcal G_+(B)=\{x:H_B(x)=p(B)\},\qquad
d_+(B)=\dim\operatorname{span}_{\mathbb R}\mathcal G_+(B).
\]

### Proposition 10.1

For an endpoint zero cut \(R\), use the block notation of Proposition
9.1 and put \(r=|R|\), \(s=n-r\).  Then

\[
\boxed{
C\bigl(\operatorname{span}\mathcal G_+(D)\bigr)
\subseteq
\bigl(\operatorname{span}\mathcal G_+(B)\bigr)^\perp
}
\tag{10.1}
\]

and consequently

\[
\boxed{
\operatorname{rank}C
\leq (r-d_+(B))+(s-d_+(D)).
}
\tag{10.2}
\]

Moreover,

\[
\boxed{d_+(A)\geq d_+(B)+d_+(D).}
\tag{10.3}
\]

#### Proof

Equation (10.1) is exactly \(x^\mathsf TCy=0\) from Proposition 9.1
for all component grounds, extended bilinearly to their spans.
The image of the \(d_+(D)\)-dimensional subspace has dimension at most
\(r-d_+(B)\), while the complementary domain has dimension
\(s-d_+(D)\), proving (10.2).

Proposition 9.1 supplies both \((x,y)\) and \((x,-y)\) as top grounds
for every component-ground pair.  Their sum and difference show that
the global ground span contains the direct sum of the two component
ground spans, proving (10.3). \(\square\)

Lemma 6.2 gives \(\|A\|_{\rm op}=O(n^{3/4})\) for a competitive
signing, hence \(\|C\|_{\rm op}=O(n^{3/4})\).  Since \(C\) is a flat
sign matrix,

\[
\operatorname{rank}C
\geq\frac{\|C\|_F^2}{\|C\|_{\rm op}^2}
=\Omega\!\left(\frac{rs}{n^{3/2}}\right).
\tag{10.4}
\]

For a balanced zero cut this forces a combined component ground-span
codimension \(\Omega(\sqrt n)\).  This is the first quantitative
constraint unavailable to the signed cut totals alone, and it again
identifies the square-root mesoscopic scale as the natural residual
wall.

The laminarization obstruction is sharp at the level of totals:
between the two opposite crossing cells \(P,W\), one may have
\(N_-=N_++1=\Theta(|P||W|)\).  Then \(t_0=1\), although uncrossing
loses \(\Theta(|P||W|)\) favorable edges; row-balanced realizations
also have low first-moment traffic.  Any successful weighted packing
theorem must therefore charge overlap to the ground-span
codimensions in (10.2), rather than to \(t_0\) alone.

## 11. Sharp wall for affine-LCP counting

The entropy \(\exp(\Omega(n^{3/4}))\) in Proposition 6.1 cannot by
itself force a spectral contradiction, even if the common affine
matrix has optimal \(O(\sqrt n)\) operator norm.

Let \(V\) be partitioned into \(q\) pairs
\(P_r=\{(r,+),(r,-)\}\).  Put \(+1\) on the edge inside each pair and,
for \(r\ne s\), set

\[
D[P_r,P_s]
=c_{rs}
\begin{pmatrix}1&-1\\-1&1\end{pmatrix},
\qquad c_{rs}=c_{sr}\in\{\pm1\}.
\tag{11.1}
\]

For every union \(Z\) of whole pairs,

\[
\boxed{
(D\mathbf1_Z)_i=
\begin{cases}
1,&i\in Z,\\
0,&i\notin Z.
\end{cases}}
\tag{11.2}
\]

Thus, with \(b=0\), every union of pairs is a strict affine-LCP
solution, and the fixed-weight layer selecting \(q/2\) pairs has

\[
\binom q{q/2}=2^{q-o(q)}
\tag{11.3}
\]

solutions.  On pair-sum vectors \(D\) acts as \(+I\), while on the
pair-difference subspace it acts as \(-I+2C\).  Choosing \(C\) to be
conference-like gives

\[
\|D\|_{\rm op}=O(\sqrt q)=O(\sqrt{|V|}).
\tag{11.4}
\]

This construction is precisely the paired quotient, so it is
structured rather than generic.  It shows, however, that
sign-pattern, determinant, P-matrix, operator-norm, or solution-count
bounds applied only to (6.11) cannot close the argument.

The absolute child-ground hypothesis eliminates the complete model:
if \(Z\) contains at least two pairs, assigning opposite spins within
each selected pair exposes the quotient quadratic form \(C[Z]\), and
the all-one state is not an absolute ground (the energy becomes
\(-|R|+4H_{C[R]}(y)\)).  A successful inverse theorem must therefore
prove that every large strict-LCP family is close to such a paired
quotient and then invoke the child-ground inequality to rule it out.
