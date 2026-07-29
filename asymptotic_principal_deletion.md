# Principal deletion, puncturing, and sparse edge repair

## 1. Normalization

For a symmetric zero-diagonal signing \(A\) of order \(n\), write

\[
E_A(x)=x^\top A x,\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|E_A(x)|,
\qquad q_n=\min_A Q(A).
\]

It is often cleaner to use the one-copy edge energy

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j=\frac12E_A(x),
\qquad M(A)=\frac12Q(A).
\]

For the augmented cuts

\[
\mathcal V_n=\{\sigma(x_ix_j)_{i<j}:
  \sigma=\pm1,\ x\in\{\pm1\}^n/\{\pm\mathbf1\}\},
\]

put

\[
h_v=a\cdot v,\qquad
M=M(A)=\max_{v\in\mathcal V_n}h_v,\qquad
g_v=\frac{M-h_v}{2}.
\]

If \(L=\binom n2\) and

\[
N_v=\{e:a_ev_e=-1\},
\]

then

\[
|N_v|=\frac{L-M}{2}+g_v.
\]

All gaps are integral.

## 2. Exact one-vertex deletion identities

Fix a positive ground state \(x\):

\[
x^\top Ax=Q(A),
\qquad
\ell_i=x_i(Ax)_i.
\]

One-spin optimality gives

\[
0\le \ell_i\le \frac{Q(A)}2,
\qquad
\sum_i\ell_i=Q(A).
\]

Let

\[
Q_i=Q(A_{-i}),\qquad
\Delta_i=Q(A)-Q_i
\]

and define the replenishment of the old ground state by

\[
g_i=Q_i-\bigl(Q(A)-2\ell_i\bigr).
\]

Then

\[
\boxed{2\ell_i=\Delta_i+g_i}
\tag{2.1}
\]

and hence

\[
\boxed{\sum_i\Delta_i=2Q(A)-\sum_i g_i.}
\tag{2.2}
\]

Thus the proposed averaged contraction

\[
\frac1n\sum_iQ(A_{-i})
\le
\left(1-\frac{3}{2n}+o(n^{-1})\right)Q(A)
\]

is exactly the assertion

\[
\sum_i g_i\le\left(\frac12+o(1)\right)Q(A).
\]

There is a second useful interpretation.  Choose an oriented absolute
ground state of \(A_{-i}\).  Its two extensions to vertex \(i\) have
oriented energies

\[
Q_i+2r_i,\qquad Q_i-2r_i.
\]

Both are at most \(Q(A)\), so

\[
|r_i|\le\frac{\Delta_i}{2}.
\]

Their two gaps below \(Q(A)\) add to \(2\Delta_i\).  Therefore a flat
principal deletion is equivalent to a pair of near-ground augmented
cuts, differing by the \(i\)-th vertex switch, with total gap
\(2\Delta_i\).

Finally, if \(A\) is an exact order-\(n\) minimizer, then its \(i\)-th
star is an optimal insertion row for its core.  Indeed replacing that
star by any other sign row gives another order-\(n\) signing, whose
quadratic norm is at least \(q_n=Q(A)\).

## 3. A conference obstruction to raw principal contraction

The coefficient \(3/2\) cannot follow from spectral regularity,
pseudorthogonality, or conference structure.

### Parity lemma

If \(n\) is even, then for every signing \(A\) and every vertex \(i\),

\[
\boxed{Q(A_{-i})\le Q(A)-2.}
\tag{3.1}
\]

For a spin \(y\) on the \(n-1\) surviving vertices, write

\[
h=y^\top A_{-i}y,\qquad
r=\sum_{j\ne i}a_{ij}y_j.
\]

The two extensions have energies \(h\pm2r\), whence

\[
Q(A)\ge \max_{\varepsilon=\pm1}|h+2\varepsilon r|
=|h|+2|r|.
\]

The integer \(r\) is a sum of the odd number \(n-1\) of signs, so
\(|r|\ge1\).  Maximizing over \(y\) proves (3.1).

Equality in (3.1) holds if a ground state has oriented local field
\(1\) at \(i\).  Consequently, if the switching-automorphism group is
vertex transitive and one ground-state orbit contains a unit local
field, then

\[
Q(A_{-i})=Q(A)-2
\quad\text{for every }i.
\tag{3.2}
\]

This flatness is fully compatible with exact optimality of every
one-star replacement.  In the one-copy normalization, if a core \(D\)
has odd order \(m\), every inserted sign row \(b\) obeys

\[
\max_x\bigl(|H_D(x)|+|b\cdot x|\bigr)\ge M(D)+1.
\tag{3.3}
\]

Indeed, evaluate at a core ground state; \(b\cdot x\) is an odd
integer.  Therefore any even-order parent satisfying

\[
M(A)=M(A_{-i})+1
\]

already uses an optimal insertion row at vertex \(i\).  If this holds
for every \(i\), no argument based only on separately reoptimizing one
star can distinguish it from a global minimizer.

This is the exact mechanism behind the flat principal minors in the
small Paley examples.  Direct exact enumeration gives:

| parent order | \(Q(A)\) | every \(Q(A_{-i})\) | \(\sum_i\Delta_i/Q(A)\) |
|---:|---:|---:|---:|
| \(6\) | \(10\) | \(8\) | \(6/5\) |
| \(14\) | \(42\) | \(40\) | \(2/3\) |
| \(18\) | \(66\) | \(64\) | \(6/11\) |

The order-six conference signing is an exact global minimizer.  At the
larger two orders global optimality is not known, so these rows are an
optimizer-compatible obstruction, not a counterexample involving
proved large-order optimizers.  They show, however, that the desired
principal theorem would rule out the leading nonresonant-conference
candidate family for a reason not visible to any spectral argument.

The order-\(13\) and order-\(17\) principal Paley cores above are also
strict one-edge local minima: every single surviving-edge flip raises
their one-copy norms from \(20\) to \(22\), and from \(32\) to \(34\),
respectively.  A repair theorem must therefore permit genuinely
nonlocal batches.

## 4. Exact fixed-cardinality repair formula

Let \(D\) be any signing of order \(m\), let

\[
L=\binom m2,\qquad M=M(D),
\]

and retain the notation \(g_v,N_v\) from Section 1.  Choose an edge set
\(S\) uniformly among the \(k\)-subsets of the \(L\) edges and flip all
edges in \(S\).  For a fixed augmented cut \(v\), put

\[
X_v=|S\cap N_v|.
\]

Then

\[
X_v\sim
\operatorname{Hyp}\left(
L,\frac{L-M}{2}+g_v,k
\right)
\]

and the exact edge-flip formula is

\[
\boxed{
M(D^S)=M-2\min_v\bigl(g_v+k-2X_v\bigr).
}
\tag{4.1}
\]

Let \(0<R<2(k/L)M\).  The event

\[
M(D^S)>M-R
\]

is contained in the union over \(v\) of

\[
X_v>
\frac{g_v+k-R/2}{2}.
\tag{4.2}
\]

Since \(X_v\le k\), only layers

\[
\boxed{g_v<k+\frac R2}
\tag{4.2a}
\]

can occur in this union.

Put

\[
p=\frac{k}{L},\qquad
\alpha=
\frac{M}{2L}-\frac{R}{4k},\qquad
\beta=\frac{1-2p}{2k}.
\]

The threshold in (4.2) exceeds the population mean by the fraction

\[
\alpha+\beta g_v.
\]

More sharply, put

\[
a_v=\frac12+\frac{g_v}{2k}-\frac{R}{4k},
\qquad
\theta_v=\frac12-\frac{M}{2L}+\frac{g_v}{L}.
\]

The hypergeometric Chernoff bound gives

\[
\Pr[(4.2)]
\le
\exp[-kD(a_v\Vert\theta_v)],
\tag{4.3}
\]

where \(D(\,\cdot\,\Vert\,\cdot\,)\) is binary relative entropy.
Pinsker's inequality gives

\[
D(a_v\Vert\theta_v)\ge2(a_v-\theta_v)^2.
\]

Consequently

\[
\Pr[(4.2)]
\le
\exp\left[-2k(\alpha+\beta g_v)^2\right].
\]

We obtain the following deterministic repair theorem.

### Theorem 4.1 (weighted-layer repair)

If

\[
\boxed{
\mathcal W_D(k,R):=
\sum_{\substack{v\in\mathcal V_m\\g_v<k+R/2}}
\exp\left[
-2k\left(
\frac{M}{2L}-\frac{R}{4k}
+\frac{(1-2k/L)g_v}{2k}
\right)^2
\right]
<1,
}
\tag{4.4}
\]

then there is a \(k\)-edge set \(S\) such that

\[
\boxed{M(D^S)\le M(D)-R.}
\tag{4.5}
\]

This is a direct use of the full augmented-cut/deep-hole profile.  It
does not assume independence between the energy values of different
spin configurations.  Replacing each quadratic exponent in (4.4) by
\(kD(a_v\Vert\theta_v)\) gives a strictly sharper criterion.

## 5. The mesoscopic scale-transfer consequence

Let \(A\) be an exact order-\(N\) minimizer in the one-copy
normalization:

\[
M(A)=M_N.
\]

Delete any \(h\) vertices, put \(m=N-h\), and call the resulting core
\(D\).  Since \(M(D)\le M_N\), Theorem 4.1 proves

\[
M_m\le M(D)-R
\]

whenever its criterion holds.

The target is

\[
T_N=
\left(\frac{N-h}{N}\right)^{3/2}M_N.
\]

If \(M(D)>T_N\), set

\[
R=M(D)-T_N.
\]

Take

\[
h\to\infty,\qquad h=o(N),\qquad
p=\frac{k}{L}=\lambda\frac hN+o(h/N),
\qquad \lambda>\frac34.
\]

Write

\[
\frac{M_N}{N^{3/2}}=c_N
\]

and parameterize the loss already supplied by deletion as

\[
M_N-M(D)=\rho\,h\sqrt N+o(h\sqrt N).
\]

Then

\[
R=
\left(\frac32c_N-\rho\right)h\sqrt N
+o(h\sqrt N),
\]

and, for gaps on the mesoscopic scale

\[
g_v=u\,h\sqrt N+o(h\sqrt N),
\]

the exponent in (4.4) is

\[
\boxed{
2k(\alpha+\beta g_v)^2
=
\frac{
\left[
u+c_N(\lambda-\tfrac34)+\rho/2
\right]^2
}{\lambda}\,h
+o(h).
}
\tag{5.1}
\]

In the worst flat-deletion case \(\rho=0\), the ground-layer rate is

\[
\boxed{
\Gamma_\lambda(c_N)
=
c_N^2
\left(
\lambda-\frac32+\frac{9}{16\lambda}
\right).
}
\tag{5.2}
\]

Thus the desired \(3/2\)-coefficient transfer is rigorously reduced to
a weighted entropy estimate for the energy layers of one principal
core.  The edge batch has size

\[
k=(1+o(1))\frac{\lambda}{2}hN.
\]

The size \(O(h\sqrt N)\) suggested by the raw energy decrement is not
enough for a diffuse active layer: a pseudorandom near-ground mixture
has edge bias only \(O(N^{-1/2})\), and therefore needs
\(\Theta(hN)\) edge trials to obtain an \(h\sqrt N\) common drift.

## 6. Exact entropy alternative when transfer fails

Theorem 4.1 has a converse in the only sense needed for a no-go
certificate.  If

\[
M_m>T_N,
\]

then every order-\(m\) signing, including every \(D^S\), has norm
larger than \(T_N\).  Hence the union in (4.2) has probability one and

\[
\boxed{\mathcal W_D(k,M(D)-T_N)\ge1.}
\tag{6.1}
\]

This forces a quantitative energy-layer population.  Fix a band width
\(\Delta>0\), put

\[
\mathcal L_j=
\{v:j\Delta\le g_v<(j+1)\Delta\},
\qquad
0\le j\le J:=\left\lceil\frac{M}{\Delta}\right\rceil,
\]

and let

\[
I(t)=2k(\alpha+\beta t)^2.
\]

Since \(I\) is increasing, (6.1) implies that some \(j\) satisfies

\[
\boxed{
\log|\mathcal L_j|
\ge
I(j\Delta)-\log(J+1).
}
\tag{6.2}
\]

For

\[
\Delta=\varepsilon h\sqrt N,\qquad h\gg\log N,
\]

equations (5.1) and (6.2) say that failure of scale transfer forces,
at some explicit normalized gap \(u=j\varepsilon\),

\[
\log|\mathcal L_j|
\ge
\frac{
\left[
u+c_N(\lambda-\tfrac34)+\rho/2
\right]^2
}{\lambda}\,h
-o(h).
\tag{6.3}
\]

Because every \(g_v\) is an integer, one may instead take
\(\Delta=1\).  If \(h\gg\log N\), the term
\(\log(J+1)=O(\log N)\) is still \(o(h)\).  Thus the obstruction can be
placed on one **exact equal-energy layer**:

\[
\boxed{
\log|\{v:g_v=g\}|
\ge
\frac{
\left[
g/(h\sqrt N)+c_N(\lambda-\tfrac34)+\rho/2
\right]^2
}{\lambda}\,h
-o(h)
}
\tag{6.4}
\]

for some integer \(g<k+R/2\).  This exact-layer form is what makes the
two-replica identities in Section 8 lossless.

This is the clean frozen-versus-entropic dichotomy:

1. if every band violates (6.3), a batch repair gives the
   \(3/2\)-coefficient scale transfer;
2. if the transfer fails, a mesoscopic near-ground layer with the
   explicit entropy rate (6.3) must exist.

## 7. What exact global optimality itself forces

Apply (4.1) to an exact order-\(n\) minimizer \(A\), without deleting
vertices.  No edge set can lower \(M(A)\).  Therefore, for every
\(k<L/2\), the hypergeometric rescue events

\[
X_v\ge\frac{g_v+k}{2}
\]

cover the entire set of \(k\)-edge flips.  The same calculation gives
the exact necessary weighted certificate

\[
\boxed{
\sum_{\substack{v\in\mathcal V_n\\g_v\le k}}
\exp\left[
-2k\left(
\frac{M}{2L}
+\frac{(1-2k/L)g_v}{2k}
\right)^2
\right]
\ge1.
}
\tag{7.1}
\]

For

\[
\frac{M}{n^{3/2}}\to c,\qquad
\frac{k}{L}=\lambda\frac hn,\qquad
h\to\infty,\ h=o(n),
\]

the ground term in (7.1) has exponent

\[
\lambda c^2h+o(h).
\]

So a low-entropy premise cannot be inserted for free: resistance to
all coefficient flips already forces an entropic rescue certificate
at every mesoscopic scale.  This is the large-order form of the
high-degeneracy obstruction seen at orders five and six.  Deletion can
still break the certificate because it changes every gap by its
boundary energy, but proving that it does so is a genuinely
overlap-sensitive statement.

## 8. Two-replica geometry of an obstructing layer

For two augmented cuts \(v=\sigma(xx^\top)\) and
\(w=\tau(yy^\top)\), their edge overlap is determined by the spin
overlap:

\[
v\cdot w
=
\sigma\tau\,
\frac{\left(\sum_i x_i y_i\right)^2-m}{2}.
\tag{8.1}
\]

If two same-orientation configurations become \(x=\mathbf1\) and
\(y=\mathbf1_{S^c}-\mathbf1_S\) after switching, then

\[
H_A(x)=I+C,\qquad H_A(y)=I-C,
\]

where \(I\) is the total signed energy internal to \(S\) and \(S^c\),
and \(C\) is the signed cross cut.  Thus two configurations in gaps
\(g_x,g_y\) from the same positive maximum satisfy

\[
I=M-(g_x+g_y),
\qquad
|C|=|g_x-g_y|\le g_x+g_y.
\tag{8.2}
\]

The complete block statement is as follows.  Switch by \(x\) and
partition the matrix according to the difference set \(S\):

\[
A^x=
\begin{pmatrix}
D&B\\
B^\top&E
\end{pmatrix}.
\]

Let \(u\) and \(z\) be the all-one restrictions to the two blocks, and
write

\[
I_D=H_D(u),\qquad I_E=H_E(z),\qquad C=u^\top Bz.
\]

Then

\[
\boxed{
I_D+I_E=M-(g_x+g_y),\qquad
C=g_y-g_x.
}
\tag{8.3}
\]

Put

\[
P_D=\max_s H_D(s),\qquad P_E=\max_t H_E(t).
\]

For arbitrary positive block maximizers \(s,t\), changing the relative
global sign of the two blocks gives

\[
P_D+P_E+|s^\top Bt|\le M.
\]

Since \(P_D+P_E\ge I_D+I_E\), it follows that

\[
\boxed{
M-(g_x+g_y)
\le P_D+P_E\le M,
\qquad
|s^\top Bt|\le g_x+g_y.
}
\tag{8.4}
\]

In particular, two exact positive ground states with difference set
\(S\) force

\[
P_D+P_E=M
\]

and the cross block is identically zero on the Cartesian product of
the two positive block-ground-state layers.

This statement also gives an exact factorization of every additive
difference multiplicity in a positive ground layer.  Represent
projective spins by fixing one vertex, and write

\[
\mathcal Z=\{z\in\mathbb F_2^{m-1}:H_A(x^z)=M\},\qquad
r_{\mathcal Z}(d)=|\mathcal Z\cap(\mathcal Z+d)|.
\]

For \(d\ne0\), choose one representation \(d=z+(z+d)\), switch by
\(z\), and partition the vertices into the two shores \(P,Q\) of
\(d\), with the fixed vertex in \(P\).  Let
\(\gamma_P,\gamma_Q\) be the numbers of projective positive ground
states of the two principal blocks.  Then

\[
\boxed{r_{\mathcal Z}(d)=2\gamma_P\gamma_Q.}
\tag{8.4a}
\]

Indeed, a spin with block restrictions \(s,t\) and its product with
\(d\) have energies

\[
H_D(s)+H_E(t)\pm s^\top Bt.
\]

Both equal \(M\) precisely when \(s,t\) are positive block ground
states and \(s^\top Bt=0\); the preceding Cartesian-annihilation
statement says that the last condition is then automatic.  The factor
two in (8.4a) comes from the two representatives of the projective
ground state on the block not containing the fixed vertex.

Consequently the additive energy obeys the sharp universal bound

\[
\boxed{
E(\mathcal Z):=\sum_d r_{\mathcal Z}(d)^2
\ge 3|\mathcal Z|^2-2|\mathcal Z|.
}
\tag{8.4b}
\]

This is just the trivial lower bound in an elementary two-group:
\(r(0)=|\mathcal Z|\), every realized nonzero difference has
\(r(d)\ge2\), and
\(\sum_{d\ne0}r(d)=|\mathcal Z|(|\mathcal Z|-1)\).
Equality holds exactly when \(\mathcal Z\) is Sidon, or, equivalently
in the present setting, when every realized nonzero difference gives
a tight principal decomposition whose two blocks each have a unique
projective positive ground state.  Thus the additive-combinatorial
problem is exactly a counting problem for tight principal
decompositions.

There is also an exact local-field decomposition.  For \(i\) in either
block let \(d_i\) be its same-block row sum and \(r_i\) its cross-block
row sum in the switching by \(x\).  If

\[
\ell_i^{(x)}=x_i(Ax)_i,\qquad
\ell_i^{(y)}=y_i(Ay)_i,
\]

then

\[
\boxed{
d_i=\frac{\ell_i^{(x)}+\ell_i^{(y)}}2,\qquad
r_i=\frac{\ell_i^{(x)}-\ell_i^{(y)}}2.
}
\tag{8.5}
\]

Moreover,

\[
\sum_i d_i=2[M-(g_x+g_y)],
\qquad
\sum_i r_i=2(g_y-g_x),
\tag{8.6}
\]

and one-spin comparison with the global maximum gives

\[
-g_x\le\ell_i^{(x)}\le M-g_x,\qquad
-g_y\le\ell_i^{(y)}\le M-g_y.
\tag{8.7}
\]

Thus an equal-gap pair has zero total cross field, while its
coordinatewise cross fields may still have a large cancelling
\(\ell_2\) mass.

For later uncrossing attempts, the exact signed-cut identities are
also worth recording.  For two vertex sets \(S,T\), put

\[
P=S\cap T,\quad Q=S\setminus T,\quad
R=T\setminus S,\quad U=(S\cup T)^c,
\]

and write \(w(X,Y)\) for the signed edge sum between \(X\) and \(Y\).
Then

\[
\boxed{
C(S)+C(T)-C(S\triangle T)
=2\bigl(w(P,U)+w(Q,R)\bigr),
}
\tag{8.8}
\]

\[
\boxed{
C(S)+C(T)-C(S\cap T)-C(S\cup T)
=2w(Q,R).
}
\tag{8.9}
\]

The error terms have no fixed sign.  This is the exact obstruction to
ordinary cut uncrossing.

If the reference switching is a true ground state, every signed cut is
nonnegative.  For two disjoint zero cuts \(S,T\), (8.9) reduces to

\[
w(S,T)\le0,\qquad
C(S\cup T)=-2w(S,T).
\tag{8.10}
\]

For pairwise disjoint zero cuts \(S_1,\ldots,S_r\),

\[
C\left(\bigcup_{j\in J}S_j\right)
=-2\sum_{\{i,j\}\subseteq J}w(S_i,S_j).
\tag{8.11}
\]

Hence all pair weights are nonpositive and

\[
\sum_{i<j}-w(S_i,S_j)\le\frac M2.
\tag{8.12}
\]

Since every nonzero aggregate is an integer of magnitude at least one,
at most \(M/2\) pairs can have strictly negative aggregate interaction.
The Caro--Wei bound therefore gives a subfamily of size at least

\[
\frac{r^2}{r+M}
\]

with all pair aggregates zero; every union of that subfamily is again
a zero cut.  This is a genuine structure conclusion, but it only starts
to bite beyond the \(M^{1/2}=\Theta(n^{3/4})\) scale and does not by
itself control the mesoscopic regime needed in Section 5.

There is a continuous-box form of the same rigidity.  In a positive
ground switching, let \(\ell_i=\sum_{j\ne i}a_{ij}\) and define the
multilinear cut polynomial

\[
F(z)=\sum_i\ell_i z_i-2\sum_{i<j}a_{ij}z_iz_j.
\tag{8.13}
\]

Independent Bernoulli sampling with means \(z_i\) shows that

\[
\boxed{F(z)\ge0\quad\text{for every }z\in[0,1]^m.}
\tag{8.14}
\]

If \(s,t\in\{0,1\}^m\) are two zero cuts and \(d=t-s\), the restriction
to their segment is exactly

\[
\boxed{
F((1-u)s+ut)
=
2u(1-u)\sum_{i<j}a_{ij}d_id_j
\quad(0\le u\le1).
}
\tag{8.15}
\]

Consequently every pair of zero cuts obeys

\[
\sum_{i<j}a_{ij}d_id_j\ge0.
\tag{8.16}
\]

For a probability distribution on zero cuts, with coordinate means
\(p_i\), averaging (8.13) gives the covariance identity

\[
\boxed{
F(p)=2\sum_{i<j}a_{ij}
\operatorname{Cov}(\mathbf1_{\{i\in S\}},
\mathbf1_{\{j\in S\}})
\ge0.
}
\tag{8.17}
\]

These are exact pairwise constraints, but they are not a PSD statement
for the full Hessian: (8.16) is only known on differences of actual
zero-cut vertices.  Small exact minimizers already show that zero cuts
need not be closed under symmetric difference.

In fact the obstruction is stronger than failure of closure.  The
following exact order-five minimizer has

\[
A=\begin{pmatrix}
0&1&1&1&1\\
1&0&1&-1&1\\
1&1&0&1&-1\\
1&-1&1&0&-1\\
1&1&-1&-1&0
\end{pmatrix},
\qquad M(A)=4.
\tag{8.17a}
\]

After fixing one vertex, its positive ground/zero-cut masks are

\[
\mathcal Z=\{0,4,6,8,9\}\subset\mathbb F_2^4.
\tag{8.17b}
\]

They form a Sidon set:

\[
E(\mathcal Z)=65=3\cdot5^2-2\cdot5,
\]

so there is no nontrivial additive quadruple at all.  The family also
fails symmetric-difference closure (\(4+6=2\notin\mathcal Z\)) and
fails the delta-matroid symmetric-exchange axiom (take the pair
\(4,9\) and the first coordinate of their difference).  The signed
Laplacian is indefinite; numerically its eigenvalues are

\[
-1.828\ldots,\ 0,\ 1,\ 3.828\ldots,\ 5.
\]

Exhaustive switching-gauge enumeration gives an even more systematic
warning:

| order | exact one-copy optimum | positive-ground sizes | all positive-ground families |
|---:|---:|---:|---|
| \(5\) | \(4\) | \(5\) | Sidon; fail symmetric exchange |
| \(6\) | \(5\) | \(6\) | Sidon; fail symmetric exchange |
| \(7\) | \(9\) | \(3,4,\) or \(7\) | Sidon |

For order seven, 2280 of the 3240 oriented switching-gauge minimizer
records also fail symmetric exchange.  These computations use all
\(2^{\binom{m-1}{2}}\) gauge-fixed signings and all projective spins,
so they are exact rather than heuristic.  Therefore neither
cardinality nor global nonnegativity of the cut polynomial supplies
the high additive energy needed to invoke Balog--Szemerédi--Gowers.
Any successful additive argument must prove that a *large
mesoscopic* obstruction cannot remain Sidon, or else bound the number
of tight principal decompositions in (8.4a).

For an equal-energy layer at gap \(g_0\), switch by one member of the
layer.  Its difference sets are again zeros of (8.13), but now the
global maximum gives only

\[
F(z)\ge-g_0\quad(z\in[0,1]^m).
\]

The same segment calculation yields the robust pair inequality

\[
\boxed{
\sum_{i<j}a_{ij}d_id_j\ge-2g_0.
}
\tag{8.18}
\]

This is the exact orientation-even quadratic constraint carried by a
large equal-gap obstruction.

There is also an elementary rank consequence that is useful for
structure-versus-randomness arguments.  From an exact layer
\(\mathcal L\), retain the larger common-orientation half and switch by
one of its members.  The resulting family \(\mathcal Z\) consists of
zero cuts and contains the empty set.  Its binary linear span has
dimension

\[
\boxed{
\dim_{\mathbb F_2}\langle\mathcal Z\rangle
\ge\log_2|\mathcal Z|
\ge\log_2|\mathcal L|-1.
}
\tag{8.19}
\]

Therefore the failure certificate (6.4) supplies
\(\Omega(h)\) linearly independent zero-cut directions in one
switching whose energy is only \(2g\) below the global maximum.
Independence alone is not closure: the other vectors in their span
need not be zero cuts.

The edge-flip rescue events also have an exact two-replica reduction.
For augmented cuts \(v,w\), let

\[
r_{ab}
=
|\{e:\mathbf1_{\{e\in N_v\}}=a,\
\mathbf1_{\{e\in N_w\}}=b\}|
\qquad(a,b\in\{0,1\}).
\]

Then

\[
\boxed{
r_{11}
=
\frac{L-h_v-h_w+v\cdot w}{4}
=
\frac{L-2M+2g_v+2g_w+v\cdot w}{4},
}
\tag{8.20}
\]

with

\[
r_{10}=|N_v|-r_{11},\qquad
r_{01}=|N_w|-r_{11},\qquad
r_{00}=L-r_{11}-r_{10}-r_{01}.
\]

If

\[
v=\sigma(xx^\top),\qquad w=\tau(yy^\top),
\]

then

\[
\boxed{
v\cdot w
=
\sigma\tau\,
\frac{\langle x,y\rangle^2-m}{2}.
}
\tag{8.21}
\]

For a uniform \(k\)-edge repair set, the four selected category counts
have the multivariate hypergeometric law with cell sizes
\((r_{00},r_{01},r_{10},r_{11})\).  Consequently the joint probability
that \(v\) and \(w\) both rescue a failed repair is an explicit
function only of

\[
(g_v,g_w,\sigma\tau,\langle x,y\rangle).
\]

This is the exact Franz--Parisi/two-replica object for the repair
problem; no further matrix statistics enter its pair law.  A
second-moment or covering argument must therefore control the overlap
distribution inside the exact layer supplied by (6.4).

If additive closure could be obtained, its consequence is completely
explicit.  Suppose a binary subspace

\[
W\le\mathbb F_2^m
\]

has the property that every \(z\in W\) is a zero cut in a fixed
switching of energy \(R\).  For each vertex \(i\), let
\(\phi_i\in W^*\) be coordinate evaluation on \(W\), and let
\(n_\phi=|\{i:\phi_i=\phi\}|\).  Fourier expansion on \(W\) gives

\[
\boxed{
\sum_{\substack{i<j\\\phi_i+\phi_j=\psi}}a_{ij}=0
\quad(\psi\ne0),
\qquad
\sum_{\substack{i<j\\\phi_i=\phi_j}}a_{ij}=R.
}
\tag{8.22}
\]

In particular,

\[
\sum_\phi\binom{n_\phi}{2}\ge R,
\qquad
\boxed{\max_\phi n_\phi\ge1+\frac{2R}{m}.}
\tag{8.23}
\]

At \(R=\Theta(m^{3/2})\), full subspace closure therefore forces a
frozen vertex type of size \(\Omega(\sqrt m)\), while every nonzero
quotient-difference class has exact signed cancellation.  This is the
precise finite-type structure promised by the zero-cut heuristic.
The unresolved additive-combinatorial step is to promote the large
set in (8.19) to a nontrivial subspace on which the cut energy remains
zero; a set of size \(\exp(\Theta(h))\) in an ambient
\(\mathbb F_2^m\) can be Sidon-like when \(h=o(m)\), so cardinality
alone cannot do this.

Consequently a large obstructing layer either clusters in Hamming
space or produces many partitions with anomalously small signed cross
energy and almost all of \(M\) stored inside the two principal blocks.
Equation (8.2) is the precise two-replica target left by (6.3).

There is a quantitative, but insufficient, consequence.  Suppose a
single band

\[
\mathcal L=\{v:g\le g_v<g+\Delta\}
\]

contains \(L_0\) augmented cuts.  At least \(L_0/2\) of them have the
same orientation sign.  If

\[
\sum_{j=0}^r\binom mj<\frac{L_0}{2},
\tag{8.24}
\]

two of those configurations have projective Hamming distance greater
than \(r\).  Applying (8.2) to that pair produces a set \(S\), with

\[
r<|S|\le\frac m2,
\]

such that

\[
\boxed{
|c_A(S,S^c)|<\Delta,\qquad
I_A(S)+I_A(S^c)>M-2(g+\Delta).
}
\tag{8.25}
\]

For example, if

\[
\log L_0\ge\kappa h,\qquad h=o(m),
\]

one may take

\[
r=\left\lfloor
\frac{\kappa h}{4\log(em/h)}
\right\rfloor
\]

for all sufficiently large \(m/h\).  Thus the entropy alternative
forces a growing low-cut mode, of size
\(\Omega(h/\log(m/h))\), whenever the obstructing band itself is
narrow.  The proof is just the Hamming-ball packing bound followed by
the exact two-replica identity.

A mere count does not finish the argument.  A family of
\(\exp(\Theta(h))\) spins can be contained in a Hamming ball of radius
\(\Theta(h/\log(N/h))\), and deleting \(h\) coordinates has
\(2^h\) fibers.  Therefore neither a first-moment count nor a
pair-distance bound alone converts (6.3) into a contradiction.
One needs an overlap-sensitive assertion: for example, a bound on the
Franz--Parisi profile of a near-minimizing signing, or a theorem saying
that a layer meeting (6.3) forces a Boolean witness strictly above its
declared maximum.

## 9. Why currently verified regularity does not close the gap

Suppose a regularized near-minimizer satisfies

\[
\|D\|_{\mathrm{op}}\le K\sqrt m.
\]

Hanson--Wright gives, for a uniform Boolean spin,

\[
\Pr(|x^\top Dx|\ge t)
\le
2\exp\left[
-c\min\left(
\frac{t^2}{m^2},
\frac{t}{K\sqrt m}
\right)
\right].
\]

At \(t=\Theta(m^{3/2})\) this only removes
\(\exp(\Omega(m/K))\) from the ambient \(2^m\) states.  It does not
give the \(\exp(O(h))\) control needed when \(h=o(m)\).  The
orientation-even \(A^2\) defect theorem likewise controls a moment of
the witness law, not the overlap profile of the deterministic
near-ground layer.

Thus the exact regularization theorem is compatible with Sections
4--7 but does not prove their entropy hypothesis.  The remaining
statement is genuinely a low-temperature, two-replica estimate.

## 10. Status

The raw averaged principal-deletion route is stopped: conference
signings give a sharp algebraic obstruction, and one of them is already
an exact finite optimizer.

The delete-then-re-sign route survives in the exact form (4.4).
It gives either:

- the required mesoscopic \(3/2\)-scale contraction, or
- the explicit weighted layer entropy certificate (6.3).

Global signing optimality itself supplies the companion necessary
certificate (7.1), explaining why unqualified low-ground-state-entropy
arguments fail.  The unresolved bridge is to use the two-replica
geometry (8.1)--(8.2) to show that the entropy rate required by (6.3)
is incompatible with asymptotic global minimality.
