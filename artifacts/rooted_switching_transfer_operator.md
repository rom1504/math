# Rooted switching-class transfer: exact state and projective obstruction

## 1. The rooted switching-class tree

Let

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
M(A)=\max_x|H_A(x)|.
\]

Given a row \(b\in\{\pm1\}^n\), write

\[
A_b=
\begin{pmatrix}
A&b\\
b^\top&0
\end{pmatrix}.
\]

Flipping the new spin gives the exact affine norm identity

\[
\boxed{
M(A_b)=
\max_x\bigl(|H_A(x)|+|b\cdot x|\bigr).
}
\tag{1.1}
\]

Switch the old vertices by \(b\).  The new row becomes all positive
and the old core becomes \(A^b\).  If \(T(A)\) denotes extension by an
all-positive row, then

\[
[A_b]=[T(A^b)]
\]

as labeled vertex-switching classes.

For a fixed labeled switching class \([A]\), the vectors
\(b\in\{\pm1\}^n/\{\pm\mathbf1\}\) give exactly \(2^{n-1}\) distinct
children of order \(n+1\).  Thus the rooted class tree has the exact
transition

\[
\boxed{
[A]\longmapsto
\bigl\{[T(A^b)]:b\in\{\pm1\}^n/\{\pm\mathbf1\}\bigr\}.
}
\tag{1.2}
\]

In particular,

\[
\boxed{
F_{n+1}
=
\min_A
\max_x\left(
|H_A(x)|+\left|\sum_i x_i\right|
\right).
}
\tag{1.3}
\]

The issue is whether (1.2) admits a compact state and a contractive
max-plus transfer.

## 2. The signed magnetization profile closes for one rooted child

For a fixed gauge, define

\[
U_A(m)=\max_{\sum_i x_i=m}H_A(x),
\qquad
L_A(m)=\min_{\sum_i x_i=m}H_A(x),
\tag{2.1}
\]

for \(m=-n,-n+2,\ldots,n\).  This signed profile determines the norm
of the distinguished all-positive child:

\[
\boxed{
M(T(A))
=
\max_m
\left(
\max\{U_A(m),-L_A(m)\}+|m|
\right).
}
\tag{2.2}
\]

It also has an exact max-plus/min-plus transition.  If \(q\) is the
total magnetization after adding the new spin, then

\[
\boxed{
\begin{aligned}
U_{T(A)}(q)
&=
\max_{\substack{s=\pm1\\m=q-s}}
\bigl(U_A(m)+sm\bigr),\\
L_{T(A)}(q)
&=
\min_{\substack{s=\pm1\\m=q-s}}
\bigl(L_A(m)+sm\bigr).
\end{aligned}
}
\tag{2.3}
\]

Indeed, the extended energy is \(H_A(x)+s\sum_i x_i\).

The map in (2.3) is order preserving and additively homogeneous.
On the \(U\)-coordinate, for example,

\[
(\mathcal TU)(q)
=
\max\{U(q-1)+(q-1),\,U(q+1)-(q+1)\}.
\tag{2.4}
\]

It is nonexpansive in the sup norm and in the additive Hilbert
seminorm

\[
d_{\rm H}(U,V)
=
\max_m(U(m)-V(m))
-\min_m(U(m)-V(m)).
\tag{2.5}
\]

It is not strictly contractive.  At the two boundary magnetizations
there is only one predecessor:

\[
(\mathcal TU)(n+1)=U(n)+n,
\qquad
(\mathcal TU)(-n-1)=U(-n)+n.
\tag{2.6}
\]

Therefore boundary differences are copied exactly.  The projective
diameter is not finite and the Birkhoff contraction coefficient is
one.

More importantly, (2.1) is Markovian only for the one distinguished
child \(T(A)\).  The true class transition (1.2) requires the profiles
of every switching \(A^b\).

## 3. Small exact counterexample to rooted-profile closure

The failure occurs even when the entire signed profile (2.1), not just
the norm or the absolute energy histogram, is retained.

Consider

\[
A_1=
\begin{pmatrix}
0&-1&1&-1&1&-1&-1\\
-1&0&1&-1&-1&1&1\\
1&1&0&1&-1&-1&-1\\
-1&-1&1&0&1&1&1\\
1&-1&-1&1&0&-1&1\\
-1&1&-1&1&-1&0&-1\\
-1&1&-1&1&1&-1&0
\end{pmatrix}
\]

and

\[
A_2=
\begin{pmatrix}
0&-1&-1&-1&1&1&-1\\
-1&0&-1&-1&-1&1&1\\
-1&-1&0&1&-1&1&1\\
-1&-1&1&0&1&-1&-1\\
1&-1&-1&1&0&1&1\\
1&1&1&-1&1&0&-1\\
-1&1&1&-1&1&-1&0
\end{pmatrix}.
\]

Both have \(M=13\), and both have exactly the same signed
energy-versus-magnetization profile:

\[
\begin{array}{c|rrrrrrrr}
m&-7&-5&-3&-1&1&3&5&7\\ \hline
U(m)&-1&3&7&7&7&7&3&-1\\
L(m)&-1&-5&-9&-13&-13&-9&-5&-1.
\end{array}
\tag{3.1}
\]

For \(b\in\{\pm1\}^7\), put

\[
\Delta_A(b)=M(A_b)-M(A).
\]

Exact enumeration of the \(128\) rows gives

\[
\begin{array}{c|rrrr}
\Delta&1&3&5&7\\ \hline
A_1&70&42&14&2\\
A_2&40&60&24&4.
\end{array}
\tag{3.2}
\]

Hence the two cores have the same state (2.1) and the same
distinguished-child norm \(M(T(A_i))=14\), but different multisets of
child norms in the rooted switching-class tree.

The entries in (3.1)--(3.2) are direct finite checks.  Exhaustive
enumeration of every signing through order six found no earlier pair
with the same signed profile and different child-norm multiset; thus
order seven is computationally minimal for this exact obstruction.

## 4. The one-step repair and why it fails at the next step

To repair (2.1) for one branching step, retain the profile in every
gauge:

\[
\boxed{
\mathcal P_A(b;m)
=
\left(
\max_{b\cdot x=m}H_A(x),
\min_{b\cdot x=m}H_A(x)
\right),
}
\tag{4.1}
\]

where \(b\in\{\pm1\}^n/\{\pm\mathbf1\}\).  This determines every child
norm in (1.2), because

\[
M(T(A^b))
=
\max_m
\left(
\max\mathcal P_A(b;m)_{\rm abs}+|m|
\right).
\tag{4.2}
\]

But \(\mathcal P_A\) is not closed under another branching step.  Let
\(C=T(A^b)\), switch the child by \((s,\tau)\), and write the old spin
as \(z\).  Its energy and total magnetization reduce to expressions of
the form

\[
H_A(z)+\tau u\,b\cdot z,
\qquad
(bs)\cdot z+u.
\tag{4.3}
\]

Computing the next signed profile therefore requires the joint
extrema of \(H_A(z)\) under simultaneous constraints on
\(b\cdot z\) and \((bs)\cdot z\).  The one-center profile (4.1)
contains only either constraint separately.

Adding the two-center profile repairs two steps; the next extension
introduces a third center.  Exact closure generates the full
multi-overlap hierarchy

\[
\max\left\{
H_A(x):
b_1\cdot x=m_1,\ldots,b_r\cdot x=m_r
\right\},
\tag{4.4}
\]

and its minimum analogue, for arbitrary \(r\).  At \(r=n\), choosing
independent centers recovers the complete labeled energy word.

## 5. The natural exact closed state: external-field support functions

The entire hierarchy (4.4) has a concise exact encoding.  Define the
two orientation-specific support functions

\[
\Phi_A^\alpha(h)
=
\max_{x\in\{\pm1\}^n}
\left(\alpha H_A(x)+h\cdot x\right),
\qquad
\alpha\in\{\pm1\},\ h\in\mathbb R^n.
\tag{5.1}
\]

They are exactly Markovian under vertex extension.  For a new row
\(b\), old external field \(h\), and new field \(t\),

\[
\boxed{
\Phi_{A_b}^\alpha(h,t)
=
\max_{u=\pm1}
\left\{
tu+\Phi_A^\alpha(h+\alpha u b)
\right\}.
}
\tag{5.2}
\]

Switching acts by coordinate reflection:

\[
\Phi_{A^s}^\alpha(h)=\Phi_A^\alpha(sh).
\tag{5.3}
\]

Within the cap/support-function formulation, this is the natural exact
repair, and it admits no information reduction: it is injective.  For
every \(x_0\),

\[
\boxed{
H_A(x_0)
=
\lim_{t\to\infty}
\left(\Phi_A^+(t x_0)-tn\right).
}
\tag{5.4}
\]

For sufficiently large \(t\), \(x_0\) is the unique maximizer, so the
limit is eventually constant.  Thus an exact state which supports all
future fields retains the original \(2^{n-1}\)-entry energy word; the
apparent profile compression is only a transform of the original
problem.

## 6. The exact transfer is a projective isometry

Let \(\mathcal T_b\) denote the max-plus operator in (5.2).  It is
order preserving and additively homogeneous:

\[
f\le g\Longrightarrow \mathcal T_bf\le\mathcal T_bg,
\qquad
\mathcal T_b(f+c)=\mathcal T_bf+c.
\tag{6.1}
\]

Consequently it is nonexpansive in the additive Hilbert metric

\[
d_{\rm H}(f,g)
=
\sup_h(f(h)-g(h))
-\inf_h(f(h)-g(h)).
\tag{6.2}
\]

In fact it is an isometry.  If

\[
(\mathcal T_bf)(h,t)
=
\max\{t+f(h+b),-t+f(h-b)\},
\tag{6.3}
\]

then, for every fixed \(h\), the first branch is selected for all
sufficiently large positive \(t\).  Hence

\[
(\mathcal T_bf-\mathcal T_bg)(h,t)
=
(f-g)(h+b)
\]

there.  As \(h+b\) ranges over the whole old field space, both the
supremum and infimum in (6.2) are inherited by the image.  Together
with nonexpansiveness this proves

\[
\boxed{
d_{\rm H}(\mathcal T_bf,\mathcal T_bg)
=d_{\rm H}(f,g).
}
\tag{6.4}
\]

The same argument applies separately to both orientations in (5.2).
Thus the exact rooted transfer has Birkhoff contraction coefficient
one, not a coefficient below one.

The obstruction is not an artifact of unbounded fields.  In the
magnetization reduction, the boundary states in (2.6) copy old data
exactly.  Removing those states destroys exactness because an extreme
magnetization can itself be a ground state.

## 7. Additive eigenvalue verdict

Finite-state max-plus Perron--Frobenius theorems obtain a unique
additive eigenvalue from irreducibility and projective contraction (or
finite projective diameter).  None of those hypotheses holds here:

1. the state dimension grows with \(n\);
2. the exact external-field state is injective and exponentially
   informative;
3. the transfer embeds the previous state isometrically on its
   boundary field faces;
4. the gauge branching forces an unbounded multi-overlap hierarchy.

An additive eigenfunction may exist on a specially selected invariant
family, but it cannot attract all exact states by a Hilbert-metric
argument.  Iterating (5.2) preserves projective differences along
successive boundary faces.

Therefore the rooted switching-class transfer does **not** yield a
convergence proof through ordinary max-plus contraction or a universal
additive eigenvalue.  The minimal exact repair is the full
external-field support function, which is equivalent to the original
energy word.  Any viable transfer theorem must deliberately coarse
grain the cap/overlap hierarchy and prove that the discarded boundary
information has only \(o(n^{3/2})\) effect; exact projective dynamics
cannot supply that loss automatically.
