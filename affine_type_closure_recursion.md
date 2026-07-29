# Affine positive-ground closure and exact type recursion

## 1. Scope and normalization

Let \(A=(a_{ij})\) be a symmetric signing of \(K_n\), with zero
diagonal, and use the one-copy energy

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j.
\]

Write

\[
p(A)=\max_x H_A(x),\qquad
\nu(A)=-\min_x H_A(x),\qquad
M(A)=\max\{p(A),\nu(A)\},
\]

and let

\[
R(A)=p(A)+\nu(A),\qquad W(A)=\frac{R(A)}2.
\]

Assign to every vertex \(i\) a type
\(\tau_i\in\mathbb F_2^d\), and put

\[
z_i(w)=(-1)^{\tau_i\cdot w}
\qquad(w\in\mathbb F_2^d).
\]

For an occupied type \(\phi\), write

\[
V_\phi=\{i:\tau_i=\phi\},\qquad k_\phi=|V_\phi|.
\]

The hypothesis in the main theorem below is exact and one-sided:

\[
\boxed{H_A(z(w))=p(A)\quad\text{for every }w.} \tag{1.1}
\]

Thus every point of the affine cloud is a positive global ground
state.  The conclusion is substantially stronger than the
character-fibre parity constraints obtained from constancy of the
affine energy alone.

## 2. The type-closure theorem

For distinct occupied types \(\phi,\psi\), define the total of the
corresponding cross block by

\[
b_{\phi\psi}
=\sum_{\substack{i\in V_\phi\\j\in V_\psi}}a_{ij}.
\tag{2.1}
\]

For a type \(\phi\), define its internal all-one energy

\[
e_\phi
=\sum_{\{i,j\}\subset V_\phi}a_{ij}.
\tag{2.2}
\]

### Theorem 2.1 (affine grounds close to the full type cube)

Under (1.1), the following statements hold.

1. Every cross-type block has zero total:

   \[
   \boxed{b_{\phi\psi}=0\quad(\phi\ne\psi).} \tag{2.3}
   \]

2. Every spin vector which is constant on each type class is a
   positive global ground state.  In particular,

   \[
   \boxed{p(A)=\sum_\phi e_\phi.} \tag{2.4}
   \]

3. Every union of type classes is positive-ground closed.  More
   precisely, for every set \(\mathcal I\) of occupied types and
   \(U=\bigcup_{\phi\in\mathcal I}V_\phi\),

   \[
   \boxed{
   p(A[U])=\sum_{\phi\in\mathcal I}p(A[V_\phi])
   =\sum_{\phi\in\mathcal I}e_\phi.
   } \tag{2.5}
   \]

4. Negative extrema and ranges obey the exact budgets

   \[
   \boxed{
   \nu(A[U])\ge
   \sum_{\phi\in\mathcal I}\nu(A[V_\phi]),
   } \tag{2.6}
   \]

   \[
   \boxed{
   R(A[U])\ge
   \sum_{\phi\in\mathcal I}R(A[V_\phi]),
   \qquad
   W(A[U])\ge
   \sum_{\phi\in\mathcal I}W(A[V_\phi]).
   } \tag{2.7}
   \]

5. If

   \[
   r_{i,\psi}=\sum_{j\in V_\psi}a_{ij}
   \quad(\psi\ne\phi),\qquad
   r_{i,\phi}=\sum_{\substack{j\in V_\phi\\j\ne i}}a_{ij}
   \quad(i\in V_\phi),
   \]

   then the exact row-domination inequality

   \[
   \boxed{
   r_{i,\phi}\ge
   \sum_{\psi\ne\phi}|r_{i,\psi}|
   } \tag{2.8}
   \]

   holds for every vertex.

#### Proof

Fix a union \(U\) of complete type classes.  For every \(w\), define
the switched cut

\[
C_U(w)=
\sum_{\substack{i\in U\\j\notin U}}
a_{ij}z_i(w)z_j(w).
\tag{2.9}
\]

Flipping all spins in \(U\) changes the energy by \(-2C_U(w)\).
Since \(z(w)\) is a positive global ground state,

\[
C_U(w)\ge0. \tag{2.10}
\]

Every edge crossing \(U\) joins two distinct types.  Its character
\(\tau_i+\tau_j\) is therefore nonzero, so

\[
\mathbb E_w C_U(w)=0. \tag{2.11}
\]

Equations (2.10)--(2.11) imply \(C_U(w)=0\) for every \(w\).
At \(w=0\),

\[
\sum_{\substack{\phi\in\mathcal I\\\psi\notin\mathcal I}}
b_{\phi\psi}=0
\quad\text{for every }\mathcal I. \tag{2.12}
\]

Taking singleton cuts first gives zero row sums in the weighted type
graph.  Taking the two-type cut \(\{\phi,\psi\}\) next gives

\[
0=C(\{\phi\})+C(\{\psi\})-C(\{\phi,\psi\})
=2b_{\phi\psi}.
\]

This proves (2.3).

If a spin vector has constant value \(s_\phi\) on \(V_\phi\), then
its energy is

\[
\sum_\phi e_\phi
+\sum_{\phi<\psi}b_{\phi\psi}s_\phi s_\psi
=\sum_\phi e_\phi.
\]

The all-one vector is \(z(0)\), so this common value is \(p(A)\).
This proves (2.4) and the full-type-cube assertion.

Now fix one type \(\phi\), a vector
\(y\in\{\pm1\}^{V_\phi}\), and choose independent uniform signs
\(s_\psi\) on every other type.  Extend \(y\) by setting all spins
on \(V_\psi\) equal to \(s_\psi\).  Cross terms have mean zero, and
therefore

\[
\mathbb E H_A(x)
=H_{A[V_\phi]}(y)+\sum_{\psi\ne\phi}e_\psi
\le p(A).
\]

Using (2.4) gives
\(H_{A[V_\phi]}(y)\le e_\phi\).  Equality holds at the all-one
vector, so

\[
p(A[V_\phi])=e_\phi. \tag{2.13}
\]

The same averaging argument, now with an arbitrary vector on a union
\(U\) and independent type signs outside \(U\), proves

\[
p(A[U])\le\sum_{\phi\in\mathcal I}e_\phi.
\]

The reverse inequality is attained by a type-constant vector, proving
(2.5).

For (2.6), choose a negative ground vector \(y_\phi\) of every
\(A[V_\phi]\), and independently multiply each whole vector
\(y_\phi\) by a uniform sign \(\sigma_\phi\).  Internal energies
remain \(-\nu(A[V_\phi])\), while every cross term has mean zero.
Thus

\[
\mathbb E_\sigma H_{A[U]}((\sigma_\phi y_\phi)_\phi)
=-\sum_{\phi\in\mathcal I}\nu(A[V_\phi]).
\]

Some realization is no larger than its mean, which proves (2.6).
Adding (2.5) and (2.6) proves (2.7).

Finally, every type-constant state is a positive ground state.
Singleton-flip stability at \(i\in V_\phi\) gives

\[
r_{i,\phi}
+\sum_{\psi\ne\phi}r_{i,\psi}s_\phi s_\psi
\ge0
\]

for all independent choices of the type signs.  Minimizing the left
side over those choices proves (2.8). \(\square\)

### Corollary 2.2 (parity rigidity)

For distinct occupied types,

\[
k_\phi k_\psi\equiv b_{\phi\psi}=0\pmod2.
\]

Hence at most one occupied type has odd multiplicity.

This eliminates the odd-parallelogram leakage from the earlier
constant-energy argument: additive collisions are relevant without
ground maximality, but exact positive-ground maximality forces every
cross-type total to vanish separately.

## 3. Collision mass and exclusion of bounded multiplicities

Define the type collision mass

\[
S_0=\sum_\phi\binom{k_\phi}{2}.
\tag{3.1}
\]

Equation (2.4) gives the elementary but decisive bound

\[
0<p(A)\le S_0. \tag{3.2}
\]

The strict positivity follows because a nonzero quadratic polynomial
has mean zero on the Boolean cube and therefore cannot have maximum
zero.

The one-sided discrepancy product, translated to the present
one-copy normalization, is

\[
p(A)\bigl(p(A)+\nu(A)\bigr)
\ge
\frac{(1-\rho^2)n^3}{6400},
\qquad
\rho=\frac{2p(A)}{n(n-1)},
\tag{3.3}
\]

whenever \((1-\rho^2)/4\ge1/n\).

### Theorem 3.1 (competitive affine grounds require
\(\sqrt n\)-scale type collisions)

Put

\[
\alpha_0=
1-\left(\frac{2S_0}{n(n-1)}\right)^2.
\]

Whenever \(\alpha_0/4\ge1/n\),

\[
\boxed{
M(A)\ge
\frac{\alpha_0 n^3}{6400S_0}-S_0.
} \tag{3.4}
\]

In particular, if \(M(A)\le Cn^{3/2}\), then

\[
\boxed{
S_0\ge
\frac{(1-o(1))n^{3/2}}{12800C}.
} \tag{3.5}
\]

If \(K=\max_\phi k_\phi\), then

\[
S_0\le\frac{(K-1)n}{2},
\]

and consequently

\[
\boxed{
K\ge
1+\frac{(1-o(1))\sqrt n}{6400C}.
} \tag{3.6}
\]

Thus no affine positive-ground family with
\(K=o(\sqrt n)\)—in particular, no bounded-multiplicity family—can
occur in a signing competitive at the \(n^{3/2}\) scale.

#### Proof

Equation (3.3) gives

\[
\nu(A)\ge
\frac{(1-\rho^2)n^3}{6400p(A)}-p(A).
\]

The right side is decreasing in \(p(A)>0\).  Use
\(p(A)\le S_0\) and \(1-\rho^2\ge\alpha_0\) to obtain (3.4).

Alternatively, if \(M(A)\le Cn^{3/2}\), then

\[
p(A)(p(A)+\nu(A))
\le 2S_0M(A)
\le2CS_0n^{3/2}.
\]

Combining this with (3.3) proves (3.5), and the bound on \(S_0\)
then proves (3.6). \(\square\)

The factor \(6400\) in (3.3) is the correct one-copy factor.  The
doubled energy theorem is
\(P(P+N)\ge(1-\rho^2)n^3/1600\); substituting
\(P=2p\), \(N=2\nu\) divides its left side by four.

The closure theorem is hereditary, so the collision estimate can be
applied to every union of type classes.

### Corollary 3.2 (hereditary small-type mass bound)

Assume

\[
M(A)\le Cn^{3/2}.
\tag{3.7}
\]

Let \(U\) be a union of type classes, put \(m=|U|\), and suppose all
classes in \(U\) have size at most \(K\).  If \(m/K\to\infty\), then

\[
\boxed{
m\le
(80\sqrt C+o(1))\,n^{3/4}\sqrt K.
} \tag{3.8}
\]

Equivalently, if \(K=t\sqrt n\) and \(m=\beta n\), then

\[
\boxed{
\beta\le\sqrt{6400Ct}+o(1).
} \tag{3.9}
\]

#### Proof

By Theorem 2.1, the principal signing \(A[U]\) is itself
positive-ground closed under its inherited type partition.  Put

\[
S_U=\sum_{V_\phi\subset U}\binom{k_\phi}{2}.
\]

Then \(p(A[U])\le S_U\le(K-1)m/2\).  A principal energy is the
expectation of the full energy over independent outside spins, so

\[
M(A[U])\le M(A)\le Cn^{3/2}. \tag{3.10}
\]

Apply the one-sided product theorem to \(A[U]\).  Since
\[
\frac{2p(A[U])}{m(m-1)}
\le\frac{K}{m-1}=o(1),
\]
it gives

\[
\frac{(1-o(1))m^3}{6400}
\le
p(A[U])R(A[U])
\le
2S_U M(A)
\le
(K-1)mCn^{3/2}.
\]

Canceling \(m\) and taking square roots proves (3.8), and (3.9) is
the same inequality after substitution. \(\square\)

### Corollary 3.3 (exact affine ground entropy is sublinear)

Under (3.7), the number \(q\) of occupied types satisfies

\[
\boxed{q=O_C(n^{5/6}).} \tag{3.11}
\]

Consequently the number of distinct states in the affine cloud is at
most

\[
\boxed{2^{O_C(n^{5/6})}=\exp(o(n)).} \tag{3.12}
\]

#### Proof

Take \(K=n^{1/6}\).  Corollary 3.2 says that the total number of
vertices in types of size at most \(K\) is
\(O_C(n^{5/6})\), and hence there are at most that many such types.
There are at most \(n/K=n^{5/6}\) larger types.  This proves (3.11).

The vectors \(z(w)\) are determined by the evaluations of the
occupied type characters.  Their effective affine dimension is at
most the number of occupied types, proving (3.12). \(\square\)

Thus exact affine positive-ground clouds of linear entropy cannot
occur at the competitive scale.  Any entropy-based replacement
argument must work with a thick near-ground layer rather than an
exact affine subspace of grounds.

## 4. The exact recursion branch

Let

\[
\Lambda=\sum_\phi k_\phi^{3/2}.
\tag{4.1}
\]

By (2.7), some occupied type satisfies

\[
\boxed{
\frac{W(A[V_\phi])}{k_\phi^{3/2}}
\le
\frac{W(A)}{\Lambda}.
} \tag{4.2}
\]

Thus, if

\[
\Lambda\ge(1-\varepsilon)n^{3/2},
\]

there is a smaller principal signing with

\[
\boxed{
\frac{W(A[V_\phi])}{k_\phi^{3/2}}
\le
\frac1{1-\varepsilon}
\frac{W(A)}{n^{3/2}}.
} \tag{4.3}
\]

This is an actual scale-preserving descent for centered width.  It is
not merely a covariance statement.

The condition describes concentration into a nearly macroscopic type:
if \(K=\max k_\phi\), then

\[
\Lambda
\le \sqrt K\sum_\phi k_\phi
=n\sqrt K.
\tag{4.4}
\]

Hence \(\Lambda\ge(1-\varepsilon)n^{3/2}\) implies
\(K\ge(1-\varepsilon)^2n\).

When (4.1) is appreciably smaller than \(n^{3/2}\), the recursion
loses a constant.  Section 5 gives the exact remaining alternatives.

## 5. Replacement entropy versus a pure mesoscopic core

For a type \(\phi\), put

\[
m_\phi=\binom{k_\phi}{2},\qquad
d_\phi=\frac{m_\phi-p(A[V_\phi])}{2}.
\tag{5.1}
\]

Here \(d_\phi\) is exactly the number of negative internal edges,
because the all-one vector is a positive ground state of
\(A[V_\phi]\).

Replacing the signing on the induced block \(E(V_\phi)\) preserves
the complete affine cap profile exactly when its total sign sum is
preserved.  Therefore its exact fibre replacement entropy is

\[
\boxed{
\mathscr R_\phi
=\log\binom{m_\phi}{d_\phi}.
} \tag{5.2}
\]

Fix \(0<\delta<1/2\).  Call a type internally mixed if
\(d_\phi\ge\delta m_\phi\), and internally pure otherwise.
For mixed classes,

\[
\mathscr R_\phi
\ge m_\phi h(\delta)-O(\log(m_\phi+1)),
\tag{5.3}
\]

where
\[
h(\delta)=-\delta\log\delta-(1-\delta)\log(1-\delta).
\]

Consequently, if mixed classes carry a fixed fraction of a collision
mass \(S_0=\Omega(n^{3/2})\), their induced blocks have total exact
replacement entropy \(\Omega(n^{3/2})\).

Otherwise the pure classes carry nearly all of \(S_0\), and

\[
\sum_{\phi\ {\rm pure}}p(A[V_\phi])
\ge
(1-2\delta)
\sum_{\phi\ {\rm pure}}m_\phi.
\tag{5.4}
\]

This is a genuine positive-energy core, not just a second-moment
surrogate.  For every threshold \(t\),

\[
\sum_{k_\phi<t}m_\phi
\le\frac{(t-1)n}{2}.
\tag{5.5}
\]

Thus if the pure collision mass is \(s n^{3/2}\), then classes of
size at least \(s\sqrt n\) carry at least
\(\frac{s}{2}n^{3/2}\) of that mass, and there are at most
\(\sqrt n/s\) such classes.  Their union \(U\) is an induced
positive-ground-closed core satisfying

\[
p(A[U])
=\sum_{V_\phi\subset U}p(A[V_\phi])
\ge
(1-2\delta)\frac{s}{2}n^{3/2}.
\tag{5.6}
\]

The fragmented alternative is therefore not an amorphous affine
cloud.  It reduces to at most \(O(\sqrt n)\) mesoscopic,
almost-all-positive internal signings, coupled only through
mean-zero ANOVA residuals.

Indeed, for a cross block \(B_{\phi\psi}\), let
\(Q_\phi=k_\phi^{-1}\mathbf1\mathbf1^\mathsf T\) and
\(P_\phi=I-Q_\phi\).  Equation (2.3) is exactly

\[
Q_\phi B_{\phi\psi}Q_\psi=0,
\]

so

\[
\boxed{
B_{\phi\psi}
=P_\phi B_{\phi\psi}P_\psi
+P_\phi B_{\phi\psi}Q_\psi
+Q_\phi B_{\phi\psi}P_\psi.
} \tag{5.7}
\]

There is no macro--macro channel.  Moreover the one-sided channels
obey the vertexwise \(\ell_1\) domination (2.8).  Equation (5.7) is
the exact ANOVA residual left for a further recursion or correlated
rounding argument.

Combining Sections 3--5 gives the proved trichotomy:

1. a nearly macroscopic type gives the scale-preserving descent
   (4.3);
2. internally mixed collision mass gives large exact induced-block
   replacement entropy through (5.2)--(5.3);
3. otherwise an \(O(\sqrt n)\)-type mesoscopic pure core carries
   \(\Theta(n^{3/2})\) actual positive energy, with all cross
   interaction confined to the residual decomposition (5.7).

The third branch is now the only unresolved affine positive-ground
model.  In particular, the formerly open bounded-multiplicity branch
has been eliminated.

## 6. A sharp square-order residual construction

The mesoscopic pure-core alternative in Section 5 is genuine.  It
cannot be removed using only natural-scale competitiveness.

Let \(k\) be an even Hadamard order, put \(n=k^2\), and partition the
vertices into \(k\) type classes \(V_1,\ldots,V_k\), each of size
\(k\).  Fix, in each \(V_a\), an orthogonal basis of balanced sign
vectors

\[
\{v_{a,b}:b\ne a\}\subset\{\pm1\}^k,
\qquad
v_{a,b}^{\mathsf T}\mathbf1=0,
\qquad
v_{a,b}^{\mathsf T}v_{a,c}
=k\,\mathbf1_{\{b=c\}}.
\tag{6.1}
\]

These are the nonconstant rows of a normalized Hadamard matrix.
Put \(+1\) on every edge internal to a type, and, for \(a\ne b\),
define the cross blocks by

\[
\boxed{
A_{ab}=v_{a,b}v_{b,a}^{\mathsf T},
\qquad A_{ba}=A_{ab}^{\mathsf T}.
} \tag{6.2}
\]

Every cross block has zero row and column sums.  Thus every
type-constant state has energy

\[
p(A)=k\binom{k}{2}=\frac{n(k-1)}2.
\tag{6.3}
\]

### Proposition 6.1 (exact spectrum and Boolean extrema)

The signing (6.2) satisfies

\[
\boxed{
p(A)=\frac{n(k-1)}2,\qquad
\nu(A)=\frac{n(k+1)}2,\qquad
M(A)=\frac{n(k+1)}2.
} \tag{6.4}
\]

Consequently,

\[
\boxed{
\frac{M(A)}{n^{3/2}}
=\frac12+\frac1{2k}
\longrightarrow\frac12.
} \tag{6.5}
\]

Every type-constant vector is a positive global ground state, every
type has zero internal replacement entropy, and
\(\max_\phi k_\phi=\sqrt n\).

#### Proof

The class-constant subspace is annihilated by every cross block, and
the internal matrix \(J_k-I_k\) acts on it by \(k-1\).

On the direct sum of the balanced subspaces, use the orthonormal
basis

\[
e_{a,b}=\frac{v_{a,b}}{\sqrt k}
\quad(a\ne b).
\]

The cross-block operator \(T\) satisfies

\[
T e_{b,a}=k e_{a,b}. \tag{6.6}
\]

It is therefore \(k\) times the involution which swaps the two
directed-edge basis vectors associated with each unordered pair
\(\{a,b\}\).  The internal matrix acts as \(-I\) on this residual
space.  Hence the full spectrum is contained in

\[
\{k-1,-k-1\}. \tag{6.7}
\]

For every Boolean vector \(x\),

\[
-n(k+1)\le x^{\mathsf T}Ax\le n(k-1).
\]

The upper endpoint is attained by every type-constant vector, proving
positive global maximality.

To attain the lower endpoint, choose a perfect matching of the
\(k\) type classes.  For every matched pair \(\{a,b\}\), set

\[
x|_{V_a}=v_{a,b},\qquad
x|_{V_b}=-v_{b,a}.
\]

Orthogonality kills every unmatched cross channel, while (6.6) gives
\(Tx=-kx\).  Since the internal action on a balanced vector is
\(-x\), this Boolean vector obeys

\[
Ax=-(k+1)x.
\]

It attains the lower spectral endpoint.  Dividing the doubled
quadratic energies by two proves (6.4)--(6.5). \(\square\)

For powers of two \(k\), Sylvester Hadamard matrices give an infinite
sequence.  This proves that the \(\Omega(\sqrt n)\) multiplicity
conclusion (3.6) has the correct order and that the third branch of
the trichotomy is indispensable.

It also identifies its quotient exactly: the ANOVA residual becomes
a signed involution on directed type pairs.  Any theorem intended to
eliminate this branch must use near-minimality below the \(1/2\)
spectral scale, rather than merely \(M(A)=O(n^{3/2})\).

### Proposition 6.2 (variable-type spectral envelope)

Let \(V_1,\ldots,V_q\) be a partition with sizes \(k_i\).  Put
\(+1\) on every internal edge, and let \(B\) be the symmetric matrix
of cross-block signs.  Assume every cross block has zero row and
column sums, and

\[
\|B\|_{\rm op}\le L\sqrt n.
\tag{6.8}
\]

Then the resulting signing \(A\) satisfies

\[
\boxed{
M(A)\le\frac12\max\left\{
\sum_i k_i\max\{k_i,L\sqrt n\}-n,\,
Ln^{3/2}+n
\right\}.
} \tag{6.9}
\]

In particular, if

\[
L=1+o(1),\qquad
\max_i k_i\le(1+o(1))\sqrt n,
\]

then

\[
\boxed{M(A)\le(1/2+o(1))n^{3/2}.} \tag{6.10}
\]

#### Proof

For a Boolean vector \(x\), put

\[
s_i=\mathbf1^{\mathsf T}x|_{V_i},\qquad
u_i=P_{k_i}x|_{V_i},\qquad d_i=\|u_i\|_2^2
=k_i-\frac{s_i^2}{k_i}.
\]

Because every cross block has zero row and column sums,

\[
x^{\mathsf T}Bx=u^{\mathsf T}Bu.
\]

The internal energy is

\[
\frac12\sum_i(s_i^2-k_i)
=p_0-\frac12\sum_i k_id_i,
\qquad
p_0=\frac12\left(\sum_i k_i^2-n\right).
\]

Therefore

\[
H_A(x)
\le
p_0+\frac12\sum_i(L\sqrt n-k_i)d_i.
\]

Using \(0\le d_i\le k_i\) only where the coefficient is positive
gives

\[
p(A)
\le
\frac12\left(
\sum_i k_i\max\{k_i,L\sqrt n\}-n
\right).
\tag{6.11}
\]

In the other orientation,

\[
-H_A(x)
\le
-p_0+\frac12\sum_i(k_i+L\sqrt n)d_i
\le
\frac12(Ln^{3/2}+n).
\tag{6.12}
\]

Equations (6.11)--(6.12) prove (6.9). \(\square\)

The spectral target in (6.8) is sharp.  Since

\[
\|B\|_F^2=n^2-\sum_i k_i^2
\]

and \(B\) acts on a space of dimension at most \(n-q\),

\[
\|B\|_{\rm op}
\ge
\sqrt{\frac{n^2-\sum_i k_i^2}{n-q}}
=(1-o(1))\sqrt n
\tag{6.13}
\]

whenever \(\max k_i=O(\sqrt n)\) and \(q=o(n)\), as holds in the
competitive affine setting by Corollary 3.3.  The construction in
Proposition 6.1 attains this floor exactly for equal square Hadamard
fibres.
Producing (6.8) for more general type profiles is therefore a sharp
flat fusion-frame or weighing-matrix design problem, not a matter of
improving the spectral estimate.

## 7. Stability with a nonzero cap slack

Assume now only

\[
p(A)-r\le H_A(z(w))\le p(A)
\quad\text{for every }w. \tag{7.1}
\]

For a union \(U\) of type classes, (2.9) satisfies

\[
C_U(w)\ge-\frac r2,\qquad
\mathbb E_w C_U(w)=0.
\tag{7.2}
\]

Hence

\[
\mathbb E_w|C_U(w)|\le r. \tag{7.3}
\]

For a fixed nonzero character \(\lambda\), the unordered type pairs
\(\{\phi,\psi\}\) satisfying \(\phi+\psi=\lambda\) form a matching.
The \(\lambda\)-Fourier coefficient of \(C_U\) is

\[
\widehat C_U(\lambda)
=
\sum_{\substack{\phi\in\mathcal I,\ \psi\notin\mathcal I\\
\phi+\psi=\lambda}}
b_{\phi\psi}. \tag{7.4}
\]

Since \(C_U+r/2\) is nonnegative and has mean \(r/2\),

\[
|\widehat C_U(\lambda)|\le\frac r2. \tag{7.5}
\]

Every subset of a matching can be realized as the set of its edges
crossing a suitable type cut.  Applying (7.5) to the positive and
negative subsets separately yields the sharp fibrewise quotient
bound

\[
\boxed{
\sum_{\substack{\{\phi,\psi\}\\\phi+\psi=\lambda}}
|b_{\phi\psi}|
\le r.
} \tag{7.6}
\]

In particular,

\[
\boxed{|b_{\phi\psi}|\le r/2.} \tag{7.7}
\]

Finally, averaging (6.1) over \(w\) kills every cross-type edge and
gives

\[
p(A)-r
\le\sum_\phi e_\phi
\le p(A). \tag{7.8}
\]

For every type,

\[
\boxed{
e_\phi\le p(A[V_\phi])\le e_\phi+r.
} \tag{7.9}
\]

The upper bound follows by extending an arbitrary internal vector by
independent uniform type signs outside the class and comparing the
mean full energy with \(p(A)\).

Equations (7.6), (7.8), and (7.9) are a stable version of type
closure.  They show exactly how an \(r\)-thick affine cap can leak:
at most \(r\) total quotient weight per character fibre, and at most
\(r\) positive-ground error per individual type.  They do not by
themselves control the sum over all character fibres, so no global
Frobenius claim is made.

## 8. Exact recursion-or-quadratic-entropy dichotomy

The cross-block conclusion (2.3) supplies the replacement branch
without any internal-purity assumption.

For two distinct types \(\phi,\psi\), the block
\(V_\phi\times V_\psi\) contains \(k_\phi k_\psi\) signs and has
total zero.  Hence \(k_\phi k_\psi\) is even, exactly half of its
signs are positive, and the block has

\[
\boxed{
\binom{k_\phi k_\psi}{k_\phi k_\psi/2}
} \tag{8.1}
\]

alternative signings with the same total.  Holding every other edge
fixed, all these alternatives preserve the complete affine cap
profile exactly.

More generally, every cross-type block can be rebalanced
independently.  Define

\[
\mathscr R_{\rm cross}
=
\sum_{\phi<\psi}
\log\binom{k_\phi k_\psi}{k_\phi k_\psi/2}.
\tag{8.2}
\]

This is a lower bound for the logarithm of the exact
profile-preserving replacement set.  Its integer affine dimension is
at least

\[
\sum_{\phi<\psi}(k_\phi k_\psi-1).
\tag{8.3}
\]

### Theorem 8.1 (recursion or quadratic replacement entropy)

Assume \(M(A)\le Cn^{3/2}\), and let
\(K=\max_\phi k_\phi\).  Then, after passing to a subsequence, one
of the following asymptotic branches occurs.

1. If \(K=n-o(n)\), the largest type \(V_*\) gives a
   scale-preserving centered-width descent:

   \[
   \boxed{
   \frac{W(A[V_*])}{K^{3/2}}
   \le
   \frac{W(A)}{n^{3/2}}+o(1)
   } \tag{8.4}
   \]

   whenever \(W(A)=O(n^{3/2})\).

2. If \(K\le(1-\delta)n\) along a subsequence, for some fixed
   \(\delta>0\), then

   \[
   \boxed{
   \mathscr R_{\rm cross}
   =
   (\log2)E_{\rm cross}+o_C(n^2)
   \ge
   \left(\frac{\delta\log2}{2}-o_C(1)\right)n^2,
   } \tag{8.5}
   \]

   and the integer replacement dimension is
   \(\Omega_\delta(n^2)\).

#### Proof

Theorem 2.1 gives \(W(A[V_*])\le W(A)\).  If \(K/n\to1\), division
by \(K^{3/2}\) proves (8.4).

The number of cross-type edges is

\[
E_{\rm cross}
=\sum_{\phi<\psi}k_\phi k_\psi
=\frac12\left(n^2-\sum_\phi k_\phi^2\right).
\tag{8.6}
\]

Since \(\sum k_\phi^2\le K\sum k_\phi=Kn\),

\[
E_{\rm cross}\ge\frac{\delta n^2}{2}. \tag{8.7}
\]

For every even \(m\ge2\),

\[
\log\binom m{m/2}
=m\log2+O(\log(m+1)),
\]

and also

\[
\binom m{m/2}\ge2^{m/2};
\]

the latter follows from
\(\binom{2r}{r}=\prod_{j=1}^r(r+j)/j\ge2^r\).
Corollary 3.3 gives \(q=O_C(n^{5/6})\) occupied types, so the total
error in the first estimate is

\[
O(q^2\log n)=O_C(n^{5/3}\log n)=o_C(n^2).
\]

Together with (8.7), this proves (8.5).  Similarly,
\(m-1\ge m/2\) in every nonempty block, so (8.3) is at least
\(E_{\rm cross}/2=\Omega_\delta(n^2)\). \(\square\)

Thus the exact affine model now has the requested
rigidity-or-recursion dichotomy:

\[
\boxed{
\text{near-total type core}
\quad\text{or}\quad
\exp(\Omega(n^2))\text{ exact cap-preserving replacements}.
} \tag{8.8}
\]

The Hadamard construction of Section 6 lies in the second branch:
its internal replacement entropy is zero, but its balanced cross
blocks supply quadratic replacement entropy.  The remaining analytic
question is whether this replacement pool contains a signing whose
traffic on the rest of the Boolean cube stays below the cap.

### 8.2 Exact law of the simultaneous balanced refill

Choose every cross block independently and uniformly from all
balanced sign matrices of its dimensions, while keeping all internal
blocks fixed.  For a fixed spin vector \(x\), put

\[
s_\phi=\sum_{i\in V_\phi}x_i,
\qquad
m_{\phi\psi}=k_\phi k_\psi.
\]

The refilled cross energy has mean zero.  Its block contributions are
independent, and the exact variance is

\[
\boxed{
V(x)
=
\sum_{\phi<\psi}
\frac{m_{\phi\psi}^2-s_\phi^2s_\psi^2}
{m_{\phi\psi}-1}.
} \tag{8.9}
\]

#### Proof

Flatten one block to a balanced random vector
\(\beta\in\{\pm1\}^m\), and put
\(u=(x_ix_j)_{(i,j)\in V_\phi\times V_\psi}\).  Then

\[
\sum_eu_e=s_\phi s_\psi.
\]

For uniform balanced \(\beta\),

\[
\mathbb E\beta_e=0,\qquad
\mathbb E\beta_e\beta_f=-\frac1{m-1}\quad(e\ne f).
\]

Consequently

\[
\operatorname{Var}\langle\beta,u\rangle
=m-\frac{(\sum_eu_e)^2-m}{m-1}
=\frac{m^2-(\sum_eu_e)^2}{m-1}.
\]

Independence of the blocks proves (8.9). \(\square\)

In particular, every type-constant state has \(V(x)=0\), so the
affine cap is preserved deterministically.  A generic state has
\(V(x)=\Theta(E_{\rm cross})=\Theta(n^2)\), and its
\(n^{3/2}\)-scale tail therefore has speed \(n\).

For internally all-positive types, if
\[
D(x)=\sum_\phi
\frac{k_\phi^2-s_\phi^2}{2}
\]
is the loss of internal energy from the type-constant cap, then

\[
V(x)
\le
\sum_{\phi<\psi}
\frac{k_\phi k_\psi}
{1-1/(k_\phi k_\psi)}
\left[
\left(1-\frac{s_\phi^2}{k_\phi^2}\right)
+
\left(1-\frac{s_\psi^2}{k_\psi^2}\right)
\right].
\tag{8.10}
\]

If every represented type has size at least \(\kappa\sqrt n\), then
(8.10) simplifies to

\[
\boxed{
V(x)\le\frac{4}{\kappa}\sqrt n\,D(x).
} \tag{8.11}
\]

Thus the exact stationarity question becomes a joint
entropy--variance problem in the type magnetizations.  The
\(\exp(\Theta(n^2))\) size of the refill cloud does not alone settle
it: each bad event indexed by a Boolean state has probability at
speed \(n\), and there are \(2^n\) possible witnesses.  A proof still
needs a sharp stratified union bound or a chained/common-law argument
using (8.9), not merely replacement counting.

### 8.3 Doubly balanced refill removes all one-sided channels

There is a stronger profile-preserving refill.  Suppose
\(k,\ell\) are even, and let \(\boldsymbol B\) be uniform over all
\(k\times\ell\) sign matrices satisfying

\[
\boldsymbol B\mathbf1=0,\qquad
\mathbf1^{\mathsf T}\boldsymbol B=0.
\tag{8.12}
\]

Its block total is zero, so it can replace any cross-type block
without changing a single affine cap energy.

Let

\[
P_k=I_k-\frac1k\mathbf1\mathbf1^{\mathsf T},
\qquad
P_\ell=I_\ell-\frac1\ell\mathbf1\mathbf1^{\mathsf T}.
\]

### Proposition 8.2 (exact doubly balanced covariance)

\[
\boxed{
\operatorname{Cov}(\operatorname{vec}\boldsymbol B)
=
\frac{k\ell}{(k-1)(\ell-1)}
(P_k\otimes P_\ell).
} \tag{8.13}
\]

Consequently, for \(x\in\{\pm1\}^k\) and
\(y\in\{\pm1\}^\ell\),

\[
\boxed{
\operatorname{Var}(x^{\mathsf T}\boldsymbol By)
=
\frac{
\left(k^2-(\mathbf1^{\mathsf T}x)^2\right)
\left(\ell^2-(\mathbf1^{\mathsf T}y)^2\right)}
{(k-1)(\ell-1)}.
} \tag{8.14}
\]

#### Proof

The law is invariant under row permutations, column permutations,
and global sign reversal.  Hence the covariance of two entries has
four possible values, according as they are equal, share only a row,
share only a column, or share neither.

The diagonal value is \(1\).  Taking covariance with a deterministic
zero row sum gives the same-row value \(-1/(\ell-1)\); the zero
column sum gives the same-column value \(-1/(k-1)\).  Taking covariance
with a different zero row then gives
\(1/((k-1)(\ell-1))\) for entries sharing neither coordinate.
These are exactly the entries of (8.13).

For the rank-one feature \(x\otimes y\), (8.13) gives

\[
\frac{k\ell}{(k-1)(\ell-1)}
(x^{\mathsf T}P_kx)(y^{\mathsf T}P_\ell y),
\]

which simplifies to (8.14). \(\square\)

For internally all-positive types, retain the deficit notation

\[
D_\phi(x)
=
\frac{k_\phi^2-s_\phi^2}{2},
\qquad s_\phi=\sum_{i\in V_\phi}x_i.
\]

Independent doubly balanced refills then have exact total variance

\[
\boxed{
V_{\rm db}(x)
=
4\sum_{\phi<\psi}
\frac{D_\phi(x)D_\psi(x)}
{(k_\phi-1)(k_\psi-1)}.
} \tag{8.15}
\]

In particular,

\[
\boxed{
V_{\rm db}(x)
\le
2\left(
\sum_\phi\frac{D_\phi(x)}{k_\phi-1}
\right)^2.
} \tag{8.16}
\]

If all types have size at least \(\kappa\sqrt n\), then, with
\(D=\sum_\phi D_\phi\),

\[
\boxed{
V_{\rm db}(x)
\le
\left(\frac{2+o(1)}{\kappa^2}\right)\frac{D^2}{n}.
} \tag{8.17}
\]

This is qualitatively stronger than (8.11): variance is now
quadratic, rather than linear, in the cap deficit.  It vanishes when
either endpoint of a block is type-constant.

The covariance constraints removed by one block have dimension
\(k+\ell-1\).  Across \(q\) even type classes their total count is

\[
\boxed{
\sum_{\phi<\psi}(k_\phi+k_\psi-1)
=(q-1)n-\binom q2.
} \tag{8.18}
\]

Thus even two macroscopic types supply \(\Theta(n)\) independent
low-variance directions, precisely the rank scale missing from the
single conditioned half-flip.

By Corollary 2.2, at most one occupied type is odd.  If a
\(k\times\ell\) block has \(\ell\) even but \(k\) odd, one may still
choose its rows independently and uniformly among balanced sign
vectors.  This preserves the block total and has exact covariance

\[
I_k\otimes\frac{\ell}{\ell-1}P_\ell,
\]

so

\[
\operatorname{Var}(x^{\mathsf T}\boldsymbol By)
=
\frac{k(\ell^2-(\mathbf1^{\mathsf T}y)^2)}{\ell-1}.
\tag{8.19}
\]

The outstanding analytic input is now narrower: establish a
uniform exponential tail for the doubly balanced ensemble strong
enough to combine (8.15) with the exact entropy of type
magnetization shells.  Covariance alone is not asserted to provide
that tail.

Tikhomirov--Youssef, *The spectral gap of dense random regular
graphs* (arXiv:1610.01765), prove a Bennett/Freedman concentration
inequality for arbitrary linear forms of uniform matrices with
prescribed row and column degrees.  Their theorem is directly
relevant to the square equal-size version of (8.12), but its constants
and exceptional conditioning event still have to be audited against
the speed-\(n\) Boolean traffic exponent.

### 8.4 The balanced shell kills every direct first-moment bound

The improved covariance does **not** make the direct traffic union
bound viable.  This can be proved without importing a concentration
theorem, using an explicit subensemble of doubly balanced matrices.

In an even-by-even block, choose independent uniform perfect matchings
of its rows and columns.  On every matched \(2\times2\) cell put

\[
\varepsilon
\begin{pmatrix}
1&-1\\
-1&1
\end{pmatrix},
\qquad \varepsilon\in\{\pm1\},
\tag{8.20}
\]

with independent signs.  Every resulting matrix has all row and
column sums zero.

Assume the internal type blocks are all positive, and consider the
shell in which every type spin vector is balanced (or differs from
balanced by one coordinate for the unique possible odd type).  If a
row matching has \(a\) discordant pairs and a column matching has
\(b\), the block energy is exactly

\[
4\sum_{j=1}^{ab}\varepsilon_j.
\tag{8.21}
\]

For a uniform matching of a balanced \(k\)-set,

\[
\mathbb E a=\frac{k^2}{4(k-1)}
=\frac k4+O(1),
\tag{8.22}
\]

and \(a=k/4+o(k)\) with exponentially high probability whenever
\(k\to\infty\).  Small types contribute only \(o(n^2)\) cross edges
by Corollary 3.2.  Since \(p(A)\le Cn^{3/2}\) in the pure case forces

\[
\max_\phi k_\phi=O_C(n^{3/4})=o(n),
\tag{8.23}
\]

the total checkerboard energy on a fixed balanced-shell state has

\[
\sigma_n^2
=\left(\frac12+o(1)\right)n^2.
\tag{8.24}
\]

Conditional on the matchings it is a sum of independent variables
of magnitude at most \(4\).  Standard binomial moderate deviations
therefore give, for every fixed \(B>0\),

\[
\boxed{
\Pr\!\left[
Y\ge Bn^{3/2}
\right]
=
\exp\left(-(B^2+o(1))n\right).
} \tag{8.25}
\]

The internal energy of a balanced state is only \(-n/2+O(q)\), which
does not change the exponent.

On the other hand, the number of balanced-shell states is

\[
\prod_\phi
\binom{k_\phi}{\lfloor k_\phi/2\rfloor}
=
\exp\left((\log2-o(1))n\right).
\tag{8.26}
\]

Indeed Corollary 3.3 gives \(q=O_C(n^{5/6})\), so the total Stirling
correction \(O(\sum_\phi\log(k_\phi+1))\) is \(o(n)\).

Combining (8.25)--(8.26), the contribution of this one shell to a
direct union-bound sum is

\[
\boxed{
\exp\left(
(\log2-B^2-o(1))n
\right).
} \tag{8.27}
\]

It diverges exponentially whenever

\[
\boxed{
B<\sqrt{\log2}=0.832554611\ldots.
} \tag{8.28}
\]

Thus row--column balancing does not bring a first-moment refill
certificate to \(1/2\); its direct threshold is already at least
\(\sqrt{\log2}\).  The route is stopped in that form.  Any use of
(8.13)--(8.17) must exploit overlap, chaining, or a quenched
second-moment argument rather than summing one-state tail
probabilities.

## 9. Exact microprofile-to-quotient variational formula

The residual branch also admits an exact smaller type-level
description.  For a microprofile

\[
r=(r_\phi)_\phi,\qquad
r_\phi\in\{\pm1\}^{V_\phi},
\]

define

\[
\Delta_\phi(r_\phi)
=p(A[V_\phi])-H_{A[V_\phi]}(r_\phi)
\ge0,
\tag{9.1}
\]

and let \(C(r)\) be the weighted signing on the occupied types with
coefficients

\[
c_{\phi\psi}(r)
=r_\phi^{\mathsf T}A_{\phi\psi}r_\psi.
\tag{9.2}
\]

Write

\[
\nu(C(r))
=-\min_{s_\phi=\pm1}
\sum_{\phi<\psi}c_{\phi\psi}(r)s_\phi s_\psi.
\tag{9.3}
\]

Then

\[
\boxed{
R(A)
=
\max_r\left\{
\sum_\phi\Delta_\phi(r_\phi)+\nu(C(r))
\right\},
} \tag{9.4}
\]

or equivalently

\[
\boxed{
W(A)
=\frac12
\max_r\left\{
\sum_\phi\Delta_\phi(r_\phi)+\nu(C(r))
\right\}.
} \tag{9.5}
\]

#### Proof

Every full spin vector can be written

\[
x|_{V_\phi}=s_\phi r_\phi
\]

for a microprofile \(r\) and type signs \(s\).  Using
\(p(A)=\sum_\phi p(A[V_\phi])\),

\[
p(A)-H_A(x)
=
\sum_\phi\Delta_\phi(r_\phi)
-
\sum_{\phi<\psi}
c_{\phi\psi}(r)s_\phi s_\psi.
\]

Maximizing first over \(s\) and then over \(r\) gives
\[
p(A)-\min_xH_A(x)
=R(A),
\]
which proves (9.4)--(9.5). \(\square\)

For internally all-positive types,

\[
\Delta_\phi(r_\phi)
=\frac{k_\phi^2-(\mathbf1^{\mathsf T}r_\phi)^2}{2}.
\tag{9.6}
\]

Thus the formerly vague ANOVA branch is an exact finite-level
optimization: local type deficits plus the negative one-sided energy
of a weighted quotient.  The Hadamard construction realizes the
quotient term by a perfect matching of directed residual channels.

## 10. Status

The exact affine positive-ground frontier is now structurally
settled:

\[
\boxed{
\begin{array}{c}
\text{a near-total recursive type core}\\
\text{or}\\
\exp(\Omega(n^2))\text{ exact cross-block replacements}.
\end{array}
}
\]

In addition, competitive exact affine clouds have only
\(O_C(n^{5/6})\) effective type dimension, bounded or
o\((\sqrt n)\) multiplicities are impossible, and the
\(\sqrt n\)-scale is sharp by the Hadamard residual construction.

What remains is analytic rather than algebraic: use minimizer
stationarity to show that the quadratic balanced-refill cloud contains
a low-traffic signing, or extract a principal descent from failure.
Equations (8.9)--(8.11) and (9.4) are the exact two formulations of
that last affine-model problem.
