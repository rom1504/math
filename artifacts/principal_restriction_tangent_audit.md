# Tangent-scale principal restriction: exact identities and the Walsh audit

Checkpoint date: 2026-07-26.

## 1. Status

Write

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_x|H_A(x)|.
\]

For a vertex partition \(V=S\sqcup T\), write

\[
A=\begin{pmatrix}B&C\\ C^\mathsf T&D\end{pmatrix}.
\]

The hoped-for scale-transfer statement is that, for a global or
competitive signing of order \(N\), one can choose \(|S|=\alpha N+o(N)\)
such that

\[
M(B)\le \alpha^{3/2}M(A)+o(N^{3/2}). \tag{1.1}
\]

This audit does not prove (1.1).  It does isolate the obstruction more
sharply.

* There is an exact principal-extension variance inequality, valid for
  every signing and every partition.
* A nonadaptive random partition always creates a leading
  \(\Theta(N^{3/2})\) cross field.
* The selector expansion shows exactly why taking a supremum before
  averaging loses that gain: its first ANOVA term is itself of leading
  order for every balanced partition.
* The scalable Bush/Walsh obstruction annihilates an
  \(\exp(\Theta(\sqrt N))\) child ground family on prescribed
  macroscopic splits of every limiting density.
* Uniform random restriction rigorously destroys this planted ground
  family, but it does not control new restriction-adaptive ground
  states.
* The same Walsh family has carefully chosen symplectic principal
  subspaces on which (1.1) holds with only an \(O(N)\) error.  Thus it
  is not a counterexample to optimized principal restriction; it is a
  counterexample to any proof using a prescribed split or only the
  inherited ground layers.

The missing assertion is therefore an **anti-emulation theorem**:
restriction ground states must not be able to choose a different
cross-field cancellation for almost every candidate subset.

## 2. Exact extension-variance inequality

Fix \(x\in\{\pm1\}^S\), put

\[
b=H_B(x),\qquad f=C^\mathsf T x,
\]

and let \(Y\) be uniform on \(\{\pm1\}^T\).  Then

\[
Z=H_D(Y)+f\cdot Y
\]

has mean zero and

\[
\operatorname{Var}Z
=\binom{|T|}{2}+\|f\|_2^2. \tag{2.1}
\]

The linear and quadratic Walsh levels are orthogonal.  Since

\[
|b+Z|\le M(A)
\]

pointwise, \(Z\) is supported on

\[
[-M(A)-b,M(A)-b].
\]

The elementary Bhatia--Davis inequality for a centered variable
supported on \([-u,v]\),

\[
\operatorname{Var}Z\le uv,
\]

therefore gives

\[
\boxed{
H_B(x)^2+\binom{|T|}{2}+\|C^\mathsf T x\|_2^2
\le M(A)^2.
} \tag{2.2}
\]

In particular, for every absolute ground state \(x\) of \(B\),

\[
\boxed{
M(B)^2+\binom{|T|}{2}+\|C^\mathsf T x\|_2^2
\le M(A)^2.
} \tag{2.3}
\]

There is also the elementary \(\ell_1\) form

\[
\boxed{
M(A)\ge
H_B(x)-M(D)+\|C^\mathsf T x\|_1
} \tag{2.4}
\]

when \(H_B(x)\ge0\), obtained by choosing \(y\) maximizing
\(|f\cdot y|\) in

\[
M(A)\ge |H_B(x)+H_D(y)|+|x^\mathsf TCy|.
\]

Equation (2.3) is sharp information for nearly flat deletion, but its
variance terms are normally only \(O(N^2)\), whereas
\(M(A)^2=\Theta(N^3)\).  It cannot by itself pay a fixed-ratio
\(3/2\)-homogeneity deficit.  Equation (2.4) is at the correct scale
provided one can prevent the \(B\)-ground state from making its cross
field atypically small.

## 3. Random partitions: the nonadaptive cross field is large

Fix a full spin \(x\in\{\pm1\}^N\), and choose \(S\) uniformly among
the \(m\)-subsets, with \(T=S^c\), \(|T|=k=N-m\).  For \(j\in T\),

\[
(C^\mathsf T x_S)_j
=\sum_{i\in S}a_{ij}x_i.
\]

Conditional on \(j\notin S\), this is a sample-without-replacement sum
of \(m\) signs from a population of \(N-1\) signs.  Convex ordering of
hypergeometric laws shows that its expected absolute value is
minimized by the most balanced population.  The local central limit
theorem then gives, uniformly in the population,

\[
\mathbb E\left[
\left|\sum_{i\in S}a_{ij}x_i\right|
\;\middle|\;j\notin S
\right]
\ge
\left(\sqrt{\frac2\pi}+o(1)\right)
\sqrt{\frac{m(k-1)}{N-2}}. \tag{3.1}
\]

Summing over \(j\) yields

\[
\boxed{
\mathbb E_S\|C^\mathsf T x_S\|_1
\ge
\left(\sqrt{\frac2\pi}+o(1)\right)
k\sqrt{\frac{m(k-1)}{N-2}}.
} \tag{3.2}
\]

For \(m/N\to\alpha\in(0,1)\), the right side is

\[
\left[
(1-\alpha)\sqrt{\frac{2\alpha(1-\alpha)}{\pi}}
+o(1)
\right]N^{3/2}. \tag{3.3}
\]

The exact second-moment companion is also useful.  With switched row
sums

\[
r_j=x_j\sum_{i\ne j}a_{ij}x_i,
\]

direct hypergeometric expansion gives

\[
\boxed{
\begin{aligned}
\mathbb E_S\|C^\mathsf T x_S\|_2^2
={}&
\frac{km(k-1)}{N-2}\\
&+
\frac{km(m-1)}
{N(N-1)(N-2)}
\sum_{j=1}^N r_j^2 .
\end{aligned}
} \tag{3.4}
\]

Thus randomization really does create the missing tangent-scale
boundary traffic for every **fixed** spin.  The unresolved issue is
the order of quantifiers:

\[
\mathbb E_S\|C^\mathsf T x_S\|_1
\quad\hbox{versus}\quad
\mathbb E_S\min_{x\in\mathcal G(B_S)}
\|C^\mathsf T x\|_1. \tag{3.5}
\]

The first is controlled by (3.2); the second is the quantity needed
in (2.4), and the restriction ground state is allowed to depend on
\(S\).

## 4. Exact selector expansion and its leading obstruction

For independent selectors \(\eta_i\in\{0,1\}\), with
\(\mathbb E\eta_i=\alpha\), put \(\xi_i=\eta_i-\alpha\).  For every
full spin \(x\),

\[
\boxed{
\begin{aligned}
H_{A[S]}(x_S)
={}&\alpha^2H_A(x)
+\alpha\sum_i\xi_i x_i(Ax)_i\\
&+\sum_{i<j}a_{ij}x_ix_j\xi_i\xi_j .
\end{aligned}
} \tag{4.1}
\]

At \(\alpha=1/2\), write \(\sigma_i=2\eta_i-1\).  Then the three terms
are exactly

\[
\frac14H_A(x),\qquad
\frac12\bigl(H_B(x_S)-H_D(x_T)\bigr),\qquad
\frac14H_{A^\sigma}(x). \tag{4.2}
\]

The middle term cannot be discarded or bounded as lower order.  Its
supremum is

\[
\frac12
\max\{P(B)+N(D),N(B)+P(D)\}
\ge \frac12\bigl(W(B)+W(D)\bigr), \tag{4.3}
\]

where \(P=\max H\), \(N=-\min H\), and
\(W=(P+N)/2\).  Universal range discrepancy makes (4.3)
\(\Theta(N^{3/2})\) for every balanced macroscopic split.

Consequently a triangle-inequality treatment of (4.1) necessarily
loses at leading scale.  Any successful selector proof must preserve
the cancellation among all three terms at the optimizing,
selector-dependent spin.

## 5. The Bush/Walsh prescribed-split obstruction

Let \(b=2^r\), \(N=b^2\), and index vertices by

\[
V=\mathbb F_2^r\times\mathbb F_2^r.
\]

Define the symplectic Walsh matrix

\[
K_{(u,v),(s,t)}
=(-1)^{v\cdot s+u\cdot t},
\qquad A=K-I. \tag{5.1}
\]

Then

\[
K=K^\mathsf T,\qquad K^2=NI,
\]

and \(A\) is a symmetric zero-diagonal signing.  For every Boolean
function \(g:\mathbb F_2^r\to\mathbb F_2\), put

\[
X_g(u,v)=(-1)^{u\cdot v+g(u)}.
\]

Direct summation gives

\[
KX_g=bX_g. \tag{5.2}
\]

Now choose an arbitrary set \(R\subseteq\mathbb F_2^r\) and the
macroscopic fibre union

\[
S_R=R\times\mathbb F_2^r.
\]

For \(u\in R\), the same calculation restricted to \(S_R\) gives

\[
K[S_R]\,X_g|_{S_R}=b\,X_g|_{S_R}, \tag{5.3}
\]

while for \(u\notin R\),

\[
K[S_R^c,S_R]\,X_g|_{S_R}=0. \tag{5.4}
\]

Because \(\|K[S_R]\|_{\rm op}\le b\), the restricted vectors in
(5.3) are positive ground states, with

\[
P(A[S_R])=\frac12|S_R|(b-1). \tag{5.5}
\]

The complement has the analogous formula, and hence

\[
P(A[S_R])+P(A[S_R^c])
=\frac12N(b-1)=P(A). \tag{5.6}
\]

Thus every density \(|R|/b\) supports a prescribed macroscopic split
on which an \(\exp(\Theta(\sqrt N))\) inherited ground family has
identically zero cross field.  Taking more states from that same
ground family cannot recover a boundary excess: the whole inherited
ground span is annihilated.

Moreover,

\[
P(A[S_R])
\sim\frac{\alpha}{2}N^{3/2},
\qquad \alpha=\frac{|R|}{b},
\]

which exceeds the desired inherited scale
\(\frac12\alpha^{3/2}N^{3/2}\) by a leading amount whenever
\(0<\alpha<1\).  This is the exact bent/Walsh obstruction to a
prescribed-split or ground-layer-only proof.

## 6. Uniform random restriction removes the planted Walsh layer

The obstruction in Section 5 is highly adaptive to the fibre
partition.  Fix one \(X_g\) and switch \(A\) by it.  The switched
matrix has constant row sum \(b-1\), Frobenius norm squared
\(N(N-1)\), and operator norm at most \(b+1\).

Let \(S\) be a uniform \(m\)-subset with \(m/N\to\alpha\in(0,1)\).
The fixed-cardinality selector expansion has no linear fluctuation,
because all row sums are equal and
\(\sum_i(\mathbf1_S(i)-m/N)=0\).  Hanson--Wright on the slice
(or independent Bernoulli Hanson--Wright followed by conditioning on
the cardinality) gives, for every fixed \(t>0\),

\[
\Pr\left(
\left|
H_{A[S]}(X_g|_S)
-\frac{(m)_2}{(N)_2}P(A)
\right|\ge tN^{3/2}
\right)
\le \exp(-c_{\alpha,t}N). \tag{6.1}
\]

Since the displayed Walsh family has only \(2^b=\exp(O(\sqrt N))\)
members, a union bound is effective.  In particular, because

\[
\alpha^{3/2}-\alpha^2>0,
\]

with probability \(1-\exp(-\Omega_\alpha(N))\), none of the planted
positive full-ground restrictions reaches

\[
\alpha^{3/2}P(A)-o(N^{3/2}). \tag{6.2}
\]

Thus uniform randomization rigorously bypasses the **known inherited
ground layer**.  What it does not do is bound \(M(A[S])\): the
principal submatrix may create a new maximizing spin chosen after
seeing \(S\).  This is precisely the adaptivity gap in (3.5).

## 7. Carefully chosen symplectic subsets satisfy scale transfer

The Walsh family is not a counterexample to optimized restriction.
Let \(U\le V\) be a nondegenerate symplectic subspace of dimension
\(2s\), and put

\[
m=|U|=2^{2s},\qquad \alpha=m/N.
\]

After choosing a symplectic basis of \(U\), the principal matrix
\(K[U]\) is the same symplectic Walsh matrix of order \(m\).  Quadratic
forms of both Arf signs supply Boolean eigenvectors with eigenvalues
\(\pm\sqrt m\).  Therefore

\[
\boxed{
M(A[U])=\frac12(m^{3/2}+m).
} \tag{7.1}
\]

The full matrix similarly satisfies

\[
M(A)=\frac12(N^{3/2}+N).
\]

Consequently

\[
\boxed{
M(A[U])
=\alpha^{3/2}M(A)
+\frac N2(\alpha-\alpha^{3/2})
=\alpha^{3/2}M(A)+O(N).
} \tag{7.2}
\]

For every fixed codimension \(2j\), this proves the desired
\(o(N^{3/2})\) transfer at density \(\alpha=4^{-j}\).

This is not enough for convergence: a single geometric ratio (or its
powers) still permits log-periodic scalar oscillation.  It does show
that the fibre-union obstruction is a failure of subset choice, not a
failure of optimized principal restriction.

## 8. Surviving theorem target

Equations (2.4) and (3.2) suggest the concrete missing lemma.  For
each fixed \(\alpha\), prove that every global minimizer has some
\(m\)-subset \(S\), \(m=\alpha N+o(N)\), and an absolute \(B_S\)-ground
state \(x_S\) for which

\[
\|A_{S,S^c}^\mathsf T x_S\|_1
\]

retains the nonadaptive tangent scale, unless \(M(B_S)\) already
satisfies (1.1).  An averaged version would be

\[
\mathbb E_S
\min_{x\in\mathcal G(B_S)}
\|A_{S,S^c}^\mathsf T x\|_1
\ge c_\alpha N^{3/2}-o(N^{3/2}), \tag{8.1}
\]

after excluding the small structured family of tight partitions.
The Bush/Walsh calculation shows that the exclusion or an inverse
classification is essential: (8.1) is false on its fibre-union
partitions, even though it is true for every fixed nonadaptive spin
on average.

The alternatives are now exact:

1. boundary traffic survives the adaptive ground-state choice and
   (2.4) yields a tangent decrement; or
2. many subsets admit restriction-specific cross-field cancellation,
   in which case one needs an inverse theorem showing that these
   cancellations organize into a replaceable affine/symplectic
   structure.

## 9. An exact incidence reduction for the adaptivity gap

There is a useful way to remove the selector-dependent spin at the
price of a sharp subset-counting problem.  For a full spin
\(z\in\{\pm1\}^N\), define

\[
\mathcal E_z^{\rm gs}(L)=
\left\{
S\in\binom{[N]}m:
z_S\in\mathcal G(A[S]),\
\left\|A_{S,S^c}^\mathsf Tz_S\right\|_1\le L
\right\}. \tag{9.1}
\]

Suppose that every \(m\)-subset \(S\) possesses at least one
restriction ground state \(x_S\) whose cross field is at most \(L\).
Choose one oriented representative \(x_S\), and count all full
extensions \(z\) satisfying \(z_S=x_S\).  There are exactly

\[
\binom Nm2^{N-m}
\]

pairs \((S,z)\).  Since there are \(2^N\) full spins, some \(z\)
occurs for at least

\[
\boxed{
2^{-m}\binom Nm
} \tag{9.2}
\]

different subsets.  Every one of those subsets belongs to
\(\mathcal E_z^{\rm gs}(L)\).  Consequently,

\[
\boxed{
\left[
\forall S\ \exists x_S\in\mathcal G(A[S]):
\|A_{S,S^c}^\mathsf T x_S\|_1\le L
\right]
\Longrightarrow
\max_z|\mathcal E_z^{\rm gs}(L)|
\ge2^{-m}\binom Nm .
} \tag{9.3}
\]

At \(m/N\to\alpha\), the required exceptional-family entropy is

\[
\frac1N\log\!\left(2^{-m}\binom Nm\right)
=h(\alpha)-\alpha\log2+o(1). \tag{9.4}
\]

For a balanced restriction this is \((\log2)/2+o(1)\), so the
remaining statement is genuinely a speed-\(N\) lower-tail problem,
not a polynomial counting problem.

The contrapositive supplies a particularly concrete inverse target:
prove, for the relevant tangent threshold \(L=c_\alpha N^{3/2}\),

\[
\boxed{
\max_z|\mathcal E_z^{\rm gs}(L)|
<2^{-m}\binom Nm
} \tag{9.5}
\]

unless \(A\) has a structured affine/symplectic quotient that can
itself be punctured at the correct scale.  The fibre-union Walsh
exceptions have only

\[
2^{\sqrt N}=\exp(o(N))
\]

members for each planted \(z\), far below (9.2).  Thus the known
Walsh obstruction does not falsify (9.5); it identifies the type of
inverse exception that a proof must classify.

The ground-state condition in (9.1) is indispensable.  There is a
spectrally tame competitive counterexample to the same count with
that condition omitted.  Let \(N=2d\), pair the vertices, choose any
order-\(d\) signing \(C=(c_{pq})\), put \(+1\) on each internal pair
edge, and put

\[
A_{\{p\},\{q\}}
=c_{pq}
\begin{pmatrix}1&-1\\-1&1\end{pmatrix}
\tag{9.6}
\]

between distinct pairs.  If a spin is antiuniform on the pair set
\(R\), with induced signs \(y_R\), then exactly

\[
\boxed{
H_A(x)=d-2|R|+4H_{C[R]}(y_R).
} \tag{9.7}
\]

Hence \(M(A)\le d+4M(C)=O(N^{3/2})\) whenever \(C\) is competitive,
and

\[
\|A\|_{\rm op}\le1+2\|C\|_{\rm op}.
\]

For the all-uniform full spin \(z\), every subset made of complete
pairs has zero cross traffic.  At \(m=N/2\) there are

\[
\binom{N/2}{N/4}
=\left(\sqrt2+o(1)\right)
2^{-N/2}\binom N{N/2} \tag{9.8}
\]

such subsets, just beyond the threshold (9.2).  Thus no universal
large-deviation estimate for low traffic alone can prove (9.5), even
under \(O(\sqrt N)\) operator norm.

The reason this example does not immediately kill (9.5) is exactly
the ground requirement.  On a union of \(s\) complete pairs, the
uniform restriction has energy \(s\), whereas an all-antiuniform
competitor has absolute energy at least

\[
4M(C[R])-s=\Theta(s^{3/2}).
\]

It is therefore not a ground state for macroscopic \(s\).  More
generally, (9.7) shows that equality candidates recursively encode the
same quadratic problem on the antiuniform pair coordinates.  A valid
inverse theorem must use this ground-state condition to either expose
a large Boolean witness or descend to a smaller competitive signing;
row-sum large deviations without the energy condition are provably
insufficient.

## 10. Exact quotient descent for the paired equality model

The paired construction is not merely an obstruction; it has an exact
descent.  Denote it by \(A=\mathcal L(C)\), globally negate \(C\) if
necessary so that

\[
-\min_yH_C(y)=M(C).
\]

From (9.7),

\[
|H_A(x)|\le d+4M(C)
\]

for every spin.  Taking every pair antiuniform and taking a negative
ground state of \(C\) gives equality.  Therefore

\[
\boxed{
M(\mathcal L(C))=d+4M(C).
} \tag{10.1}
\]

Choosing one fixed vertex from each pair gives a principal submatrix
equal to \(C\).  Hence the paired low-traffic model automatically has
a half-order principal restriction satisfying

\[
\boxed{
M(C)=\frac{M(A)-d}{4}
<2^{-3/2}M(A)
} \tag{10.2}
\]

for all sufficiently large competitive instances.  In normalized
form,

\[
\boxed{
\frac{M(C)}{d^{3/2}}
=\frac1{\sqrt2}
\frac{M(A)}{(2d)^{3/2}}
+O(d^{-1/2}).
} \tag{10.3}
\]

Thus the apparent equality case in the traffic-only incidence count
is strictly easier after taking one representative per type.

The joint ground-and-traffic data also descend exactly.  Let
\(Q\subseteq[d]\), let \(S_Q\) be the union of its full vertex pairs,
and let \(z\) be antiuniform on every pair, with pair signs \(y\).
Then

\[
A[S_Q]=\mathcal L(C[Q]) \tag{10.4}
\]

and

\[
\boxed{
\left\|A_{S_Q,S_Q^c}^\mathsf Tz_{S_Q}\right\|_1
=4\left\|C_{Q,Q^c}^\mathsf Ty_Q\right\|_1.
} \tag{10.5}
\]

If \(y_Q\) is a negative absolute ground state of \(C[Q]\), then
\(z_{S_Q}\) is an absolute ground state of \(A[S_Q]\), and

\[
M(A[S_Q])=|Q|+4M(C[Q]). \tag{10.6}
\]

Consequently, an incidence family concentrated on the paired
equality model either:

* uses uniform pair modes, in which case the restrictions are not
  ground states at macroscopic order; or
* uses antiuniform pair modes, in which case both the ground
  certificate and the low-traffic certificate descend to the
  half-order signing \(C\), whose normalized objective is strictly
  smaller by (10.3).

This proves the desired inverse/descent assertion for the exact
two-point affine model.  The general affine-ground theorem supplies
the natural extension: a competitive exact affine ground family has
a type class of size \(\Omega(\sqrt N)\), and all nonzero intertype
block sums vanish.  What is still missing is a purification theorem
turning those weighted type sums into an unweighted principal
quotient with the analogue of (10.2).  The paired calculation shows
the correct conclusion and normalization when every type has size
two.

## 11. Verdict on this route

The raw random-selector route should not continue as a standalone
large-deviation argument:

* its first ANOVA term is leading scale;
* low traffic alone has a spectrally tame speed-\(N\) counterfamily;
* and the needed fixed-spin tail must be intersected with the full
  child-ground condition.

The route remains viable only as an **inverse/descent program**:

1. use the exact incidence threshold
   \(2^{-m}\binom Nm\);
2. prove that a joint ground-plus-low-traffic family above that
   threshold has an affine/type quotient;
3. purify the quotient into a smaller principal signing with no worse
   normalized objective.

Step 3 is exact for paired types by (10.1)--(10.6).  Extending it from
two-point fibres to the \(\Omega(\sqrt N)\) type classes forced by
affine ground rigidity is the next nonredundant lemma.  Without that
new structural input, further generic concentration estimates will
repeat a falsified traffic-only argument.

## 12. The exact paired exception has depth at most one

Let

\[
c(A)=\frac{M(A)}{|A|^{3/2}}.
\]

For a paired lift of an order-\(d\) signing, (10.1) gives

\[
\boxed{
c(\mathcal L(C))
=\sqrt2\,c(C)+O(d^{-1/2}).
} \tag{12.1}
\]

Every signing sequence has the universal lower bound

\[
c(C)\ge c_* -o(1),
\qquad
c_*=0.336493364431\ldots ,
\]

while an asymptotic global minimizer has

\[
c(A)\le\frac12+o(1).
\]

It follows that a single exact paired exception inside a minimizing
sequence is confined to the narrow interval

\[
\boxed{
0.475873479627\ldots
=\sqrt2\,c_*
\le c(\mathcal L(C))
\le\frac12+o(1),
} \tag{12.2}
\]

and its quotient satisfies

\[
c_*-o(1)\le c(C)
\le\frac1{2\sqrt2}+o(1)
=0.353553390593\ldots+o(1). \tag{12.3}
\]

Two consecutive exact pair lifts are impossible in a minimizing
sequence, since

\[
\boxed{
c(\mathcal L^2(D))
=2c(D)+o(1)
\ge 2c_* -o(1)
=0.672986728863\ldots-o(1)
>\frac12+o(1).
} \tag{12.4}
\]

Thus the exact paired inverse class is closed at bounded depth:
one quotient descent always produces a strictly better normalized
principal problem, and the descended problem cannot itself carry the
same exact paired obstruction.  Conditional on an inverse theorem
saying that every incidence-saturating joint ground/traffic family is
of paired type, the next scale must fall on the nonexceptional
traffic branch.  The unresolved extension is to prove the analogous
bounded-depth growth factor for general affine type classes after
purification.
