# Affine cap clouds: exact conditional kernels and the missing replacement condition

## 1. Status and purpose

This note audits the proposed implication

\[
\text{flat conditional kernel}
\quad\Longrightarrow\quad
\text{localized sign replacement}.
\]

It gives three rigorous conclusions.

1. On an affine cut-feature cloud, every conditional Schur complement
   can be computed exactly by grouping edges according to their
   character type.
2. Exact conditional degeneracy can come entirely from remote copies
   of the same character.  It need not provide any local freedom to
   alter signs.
3. There is a cut-triangle-compatible example with feature entropy
   \(\Theta(n)\), a dense flat kernel vector, and proportional induced
   blocks on which the exact cap-preserving replacement is unique.

Thus (7.2), the triangle identities, flatness, and large shell entropy
are not by themselves sufficient for the inverse theorem proposed in
`cut_cone_tangent_entropy.md`.  The full theorem may still be true,
but its proof must use the cap envelope/global extremality, or an
additional integer "replacement entropy" condition defined below.

Throughout this note,

\[
H_a(z)=\sum_{e=\{i,j\}}a_ez_iz_j
\]

is in the one-copy normalization.

## 2. Affine cut clouds and character fibres

Let \(V\) be a vertex set, let \(D=\mathbb F_2^d\), and assign a type

\[
\tau_i\in D
\qquad(i\in V).
\]

For \(w\in D\), define the spin vector

\[
z_i(w)=(-1)^{\tau_i\cdot w}.
\tag{2.1}
\]

The feature on an edge \(e=\{i,j\}\) is

\[
\phi_e(w)=z_i(w)z_j(w)
=(-1)^{\lambda(e)\cdot w},
\qquad
\lambda(e)=\tau_i+\tau_j.
\tag{2.2}
\]

For \(\lambda\in D\), write

\[
\mathcal E_\lambda
=\{e:\lambda(e)=\lambda\}.
\tag{2.3}
\]

These are the character fibres of the edge set.  They automatically
satisfy every cut-triangle identity, since

\[
\phi_{\{i,j\}}\phi_{\{j,k\}}\phi_{\{k,i\}}=1.
\tag{2.4}
\]

Let \(w\) be uniform on \(D\), and let
\(\Sigma=\operatorname{Cov}(\phi(w))\).  Nonzero characters are
orthonormal, while the zero character is constant.  Consequently,

\[
\Sigma_{ef}
=
\begin{cases}
1,&\lambda(e)=\lambda(f)\ne0,\\
0,&\text{otherwise}.
\end{cases}
\tag{2.5}
\]

Thus \(\Sigma\) is a direct sum of all-ones matrices, one for every
nonzero character fibre, together with a zero block on
\(\mathcal E_0\).

## 3. Exact Schur-complement formula

Let \(T\) be any edge set and let \(u\in\mathbb R^T\).  Put

\[
b_\lambda(u;T)
=\sum_{e\in T\cap\mathcal E_\lambda}u_e.
\tag{3.1}
\]

### Theorem 3.1 (conditional variance is unsupported fibre mass)

For the affine cloud (2.1), the conditional covariance of the
features in \(T\), after optimal linear prediction from \(T^c\),
satisfies

\[
\boxed{
u^\top\Sigma_{T\mid T^c}u
=
\sum_{\substack{\lambda\ne0\\
\mathcal E_\lambda\cap T^c=\varnothing}}
b_\lambda(u;T)^2.
}
\tag{3.2}
\]

#### Proof

Fix a nonzero \(\lambda\).  If its fibre has \(q\) edges in \(T\) and
\(r\) edges in \(T^c\), the corresponding covariance block is

\[
\begin{pmatrix}
J_q&J_{q,r}\\
J_{r,q}&J_r
\end{pmatrix}.
\]

If \(r>0\), then \(J_r^\dagger=J_r/r^2\), and hence

\[
J_q-J_{q,r}J_r^\dagger J_{r,q}=0.
\]

One outside copy determines the common character exactly.  If
\(r=0\), the residual block is \(J_q\), whose quadratic form at
\(u\) is \(b_\lambda(u;T)^2\).  Different nonzero characters are
orthogonal, and the zero-character block has zero covariance.  Adding
the blocks proves (3.2). \(\square\)

The formula shows that Schur-complement degeneracy is a statement
about *availability of duplicate characters*, not necessarily about
local geometry in the vertex set.

### Corollary 3.2 (balanced fibres annihilate every Schur complement)

Suppose a flat signing \(h\in\{\pm1\}^{E}\) obeys

\[
\sum_{e\in\mathcal E_\lambda}h_e=0
\qquad(\lambda\ne0).
\tag{3.3}
\]

Then, for every edge block \(T\),

\[
\boxed{
h_T^\top\Sigma_{T\mid T^c}h_T=0.
}
\tag{3.4}
\]

Indeed, a fibre represented outside \(T\) contributes zero to (3.2),
while a fibre wholly contained in \(T\) contributes its balanced
total, also zero.

Moreover,

\[
H_h(z(w))
=
\sum_{\lambda\in D}
\left(\sum_{e\in\mathcal E_\lambda}h_e\right)
(-1)^{\lambda\cdot w}.
\tag{3.5}
\]

Thus (3.3) makes the energy constant on the affine cloud; if the
zero-fibre sum also vanishes, that constant is zero.

## 4. The exact integer replacement law

Let \(\beta\in\{\pm1\}^T\) replace \(h_T\), leaving \(h_{T^c}\)
unchanged.  Define the fibre-sum changes

\[
\Delta_\lambda
=
\sum_{e\in T\cap\mathcal E_\lambda}(\beta_e-h_e).
\tag{4.1}
\]

The change of energy on the affine cloud is the Fourier polynomial

\[
D_\beta(w)
=H_{h^{T\to\beta}}(z(w))-H_h(z(w))
=\sum_{\lambda\in D}\Delta_\lambda
(-1)^{\lambda\cdot w}.
\tag{4.2}
\]

Character orthogonality gives the exact Parseval identity

\[
\boxed{
\mathbb E_wD_\beta(w)^2
=\sum_{\lambda\in D}\Delta_\lambda^2.
}
\tag{4.3}
\]

Consequently, a replacement preserves the full affine cap profile
exactly if and only if

\[
\boxed{
\sum_{e\in T\cap\mathcal E_\lambda}\beta_e
=
\sum_{e\in T\cap\mathcal E_\lambda}h_e
\quad\text{for every }\lambda.
}
\tag{4.4}
\]

This is an integer constraint which is invisible to the covariance
kernel.  If

\[
m_\lambda(T)=|T\cap\mathcal E_\lambda|,
\qquad
p_\lambda(T)=|\{e\in T\cap\mathcal E_\lambda:h_e=1\}|,
\]

then the number of exact profile-preserving replacements is

\[
\boxed{
\#\mathcal R_h(T)
=
\prod_{\lambda\in D}
\binom{m_\lambda(T)}{p_\lambda(T)}.
}
\tag{4.5}
\]

This motivates the **fibre replacement entropy**

\[
\boxed{
\mathscr R_h(T)
=
\sum_{\lambda\in D}
\log\binom{m_\lambda(T)}{p_\lambda(T)}.
}
\tag{4.6}
\]

The conditional variance in (3.2) and the replacement entropy in
(4.6) are logically independent.  An outside copy of a character
makes its conditional variance zero, but it creates no alternative
signing inside \(T\).

There is also a quantitative approximate form.  If

\[
\sup_w|D_\beta(w)|\le e,
\]

then

\[
\boxed{
\sum_\lambda\Delta_\lambda^2\le e^2.
}
\tag{4.7}
\]

Since every \(\Delta_\lambda\) is an even integer, at most \(e^2/4\)
character sums can be changed nontrivially.

## 5. A linear-entropy proportional-block obstruction

The following example shows that the distinction above persists at
linear entropy and on proportional induced vertex blocks.

Let \(d\ge4\) be even, put \(n=2d\), and index the vertices by

\[
V=[d]\times\{0,1\}.
\]

Assign the repeated types

\[
\tau_{(i,0)}=\tau_{(i,1)}=e_i\in\mathbb F_2^d.
\tag{5.1}
\]

For \(i<j\), sign the \(2\times2\) block between the two vertex pairs
by

\[
h_{(i,a),(j,b)}=(-1)^a.
\tag{5.2}
\]

Its four signs have sum zero.  Sign the within-pair edge
\(\{(i,0),(i,1)\}\) by \(\sigma_i\), where
\(\sigma_i\in\{\pm1\}\) and

\[
\sum_i\sigma_i=0.
\tag{5.3}
\]

For every \(w\in\mathbb F_2^d\), use the pair-constant state

\[
z_{(i,0)}(w)=z_{(i,1)}(w)=(-1)^{w_i}.
\tag{5.4}
\]

Every nonzero character fibre \(e_i+e_j\) consists of the four edges
between pairs \(i\) and \(j\), and (5.2) balances it.  The
zero-character sum is zero by (5.3).  Therefore

\[
\boxed{
H_h(z(w))=0
\quad\text{for all }w,
}
\tag{5.5}
\]

and Corollary 3.2 gives

\[
\boxed{
h_T^\top\Sigma_{T\mid T^c}h_T=0
\quad\text{for every edge block }T.
}
\tag{5.6}
\]

The edge-feature cloud has \(2^{d-1}=2^{n/2-1}\) distinct points:
the characters \(e_i+e_j\) span the even-weight subspace.  Hence this
is not a low-cardinality or \(o(n)\)-entropy artefact.

Now choose the proportional transversal

\[
S=\{(i,0):i\in[d]\},
\qquad |S|=n/2,
\qquad T=E(S).
\tag{5.7}
\]

Every edge of \(T\) has a different character \(e_i+e_j\), while
each such character has three additional copies outside \(T\).
Thus \(T\) is perfectly conditionally predictable, but every fibre
inside \(T\) has size one.  Equations (4.4)--(4.6) give

\[
\boxed{
\mathcal R_h(T)=\{h_T\},
\qquad
\mathscr R_h(T)=0.
}
\tag{5.8}
\]

More quantitatively, if a replacement changes \(k\) edges of \(T\),
then all affected Fourier characters are distinct and

\[
\boxed{
\mathbb E_wD_\beta(w)^2=4k,
\qquad
\sup_w|D_\beta(w)|\ge2\sqrt{k}.
}
\tag{5.9}
\]

The same conclusion holds for every transversal choosing one vertex
from each repeated pair.  There are \(2^d\) such proportional induced
blocks.

This is the promised obstruction:

\[
\boxed{
\begin{array}{c}
\text{cut-triangle identities}
+\text{ dense flat }h
+\text{ linear feature entropy}\\
+\text{ exact kernel for every conditional covariance}
\end{array}
\quad\not\Longrightarrow\quad
\begin{array}{c}
\text{nontrivial exact cap-preserving}\\
\text{replacement on a proportional induced block.}
\end{array}
}
\tag{5.10}
\]

## 6. Why the example does not refute the full global theorem

The signing in Section 5 is deliberately not competitive.  It has a
quadratic Boolean witness.  Put the first \(d/2\) vertex pairs in
antiuniform mode \((1,-1)\), and the last \(d/2\) pairs in uniform
mode \((1,1)\).  Every cross block from the first group to the second
contributes \(4\), while same-group cross blocks contribute zero.
The within-pair edges cost at most \(d\).  Hence

\[
\boxed{
Q(h)\ge d^2-d=\frac{n^2}{4}-\frac n2.
}
\tag{6.1}
\]

Likewise, Proposition 5.1 of
`tight_principal_decomposition_count.md` proves the stronger
structural statement: if all pair-constant states were *positive
ground states* of a signing, then

\[
Q(A)\ge2d^2-d=\frac{n^2-n}{2}.
\tag{6.2}
\]

Thus competitive \(O(n^{3/2})\) scale and genuine endpoint
extremality can rule out this particular repeated-pair obstruction.
The example does **not** disprove the full proposed inverse theorem,
which also assumes the two-sided cap envelope and global minimality.
It does prove that those hypotheses must be used essentially; they
cannot be discarded after deriving (7.2).

### Proposition 6.1 (even-type pair rigidity)

There is a partial extension of the paired-coordinate theorem which
uses only typewise local-field domination.

Suppose that all states (2.1) are positive ground states of a signing
\(A\).  Assume every occupied type class has even cardinality, and
let

\[
\mathcal P=\{\tau:|V_\tau|=2\},
\qquad m=|\mathcal P|.
\]

Then

\[
\boxed{
Q(A)\ge 2m^2-m.
}
\tag{6.3}
\]

In particular, \(Q(A)=O(n^{3/2})\) implies

\[
m=O(n^{3/4}).
\tag{6.4}
\]

#### Proof

For \(i\in V_\phi\), the signed local field on the affine ground
family is

\[
\ell_i(w)
=z_i(w)\sum_{j\ne i}a_{ij}z_j(w)
=\sum_\tau r_{i,\tau}
(-1)^{(\phi+\tau)\cdot w},
\tag{6.5}
\]

where

\[
r_{i,\tau}
=\sum_{\substack{j\in V_\tau\\j\ne i}}a_{ij}.
\]

Every ground state is stable under a singleton flip, so
\(\ell_i(w)\ge0\) for every \(w\).  Every Fourier coefficient of a
nonnegative function has absolute value at most its mean.  Therefore

\[
|r_{i,\tau}|
\le r_{i,\phi}
\le |V_\phi|-1.
\tag{6.6}
\]

If \(|V_\phi|=2\), then \(r_{i,\phi}\) is one sign.  Its
nonnegativity forces the within-pair edge to be \(+1\).  For every
other occupied \(\tau\), the row sum \(r_{i,\tau}\) is a sum of an
even number of signs.  By (6.6) its absolute value is at most one, so

\[
r_{i,\tau}=0.
\tag{6.7}
\]

Let \(U\) be the union of the \(m\) two-vertex type classes.  Fix any
affine ground state outside \(U\).  Since that state is constant
inside each type class, (6.7) says that its external field on every
vertex of \(U\) is zero.  Hence changing the spins arbitrarily inside
\(U\) changes the full energy by exactly the change of the principal
energy \(H_{A[U]}\).

Every pair-constant state on \(U\) has principal energy \(m\):
within-pair edges contribute \(m\), while every \(2\times2\)
cross-block has zero total sum by (6.7).  Fixing the outside of one
global ground state shows that no state on \(U\) can have principal
energy greater than \(m\).  Thus all \(2^m\) pair-constant states are
positive ground states of \(A[U]\).

Proposition 5.1 of `tight_principal_decomposition_count.md` gives

\[
Q(A[U])\ge2m^2-m.
\]

Finally, for a spin state on \(U\), the expectation of the full
energy over independent uniform outside spins equals its principal
energy.  Some extension therefore has absolute full energy at least
the absolute principal energy, so \(Q(A)\ge Q(A[U])\).  This proves
(6.3). \(\square\)

The parity hypothesis exposes two apparent leakage channels.  An odd
outside type allows a row sum \(\pm1\) from a two-vertex class, while
classes of size at least four permit still larger row sums.  The next
proposition removes the first channel; larger type multiplicities
remain open.

The odd-type leakage can in fact be removed by combining ground-state
stability with the one-sided discrepancy product.

### Proposition 6.2 (pair rigidity with arbitrary outside types)

Suppose again that all states (2.1) are positive ground states of
\(A\), but make no parity assumption on the occupied type classes.
If exactly \(m\) type classes have size two, then, as \(m\to\infty\),

\[
\boxed{
Q(A)\ge(1-o(1))\frac{m^2}{4000}.
}
\tag{6.8}
\]

The numerical constant only records the current
Bollobás--Scott discrepancy constant; the structural content is
\(Q(A)=\Omega(m^2)\).

#### Proof

Let \(U\) be the union of the two-vertex type classes.  As in
Proposition 6.1, (6.6) forces every within-pair edge to be \(+1\).
Between two size-two classes, every row sum is an even integer of
absolute value at most one, and is therefore zero.  The same holds
for column sums.  Thus every \(2\times2\) cross block has the form

\[
B_{ij}=c_{ij}vv^\top,
\qquad
v=(1,-1),
\qquad
c_{ij}\in\{\pm1\}.
\tag{6.9}
\]

The coefficients \(c_{ij}\) form an auxiliary signing \(C\) of
\(K_m\).

For a ground parameter \(w\), let \(s_i(w)\) be the common ground
spin on pair \(i\), and let

\[
f_{i,a}(w)
=\sum_{j\notin U}a_{(i,a),j}z_j(w)
\qquad(a=0,1)
\tag{6.10}
\]

be its field from outside \(U\).  Define its signed version

\[
b_{i,a}(w)=s_i(w)f_{i,a}(w).
\tag{6.11}
\]

The within-pair edge contributes \(1\) to the signed local field,
and all other pair-constant classes in \(U\) contribute zero by
(6.9).  Hence singleton-flip stability of every affine ground gives

\[
1+b_{i,a}(w)\ge0.
\tag{6.12}
\]

There are no vertices of type \(\tau_i\) outside its complete
two-vertex class.  Every term in (6.10)--(6.11) is therefore a
nonzero character of \(w\), and

\[
\mathbb E_wb_{i,a}(w)=0.
\tag{6.13}
\]

For any mean-zero random variable \(b\ge-1\),

\[
\mathbb E|b|
=2\mathbb E(-b)_+
\le2.
\tag{6.14}
\]

Consequently there is a ground parameter \(w_0\) for which

\[
\sum_{i=1}^m\sum_{a=0}^1|b_{i,a}(w_0)|
\le4m.
\tag{6.15}
\]

Hold the spins outside \(U\) fixed at \(z(w_0)\).  Replace the
uniform ground spin \(s_i\mathbf1\) on every pair by an antiuniform
spin \(t_iv\), and put \(u_i=t_is_i\).  Relative to the ground
energy, the within-pair contribution changes by \(-2m\), the pair
cross-block contribution becomes

\[
4\sum_{i<j}(c_{ij}s_is_j)u_iu_j,
\]

and the outside-field change on pair \(i\) is

\[
u_i(b_{i,0}-b_{i,1})-(b_{i,0}+b_{i,1}).
\]

Global maximality of the ground state therefore implies, for every
\(u\in\{\pm1\}^m\),

\[
\begin{aligned}
4\sum_{i<j}(c_{ij}s_is_j)u_iu_j
&\le
2m+
\sum_i\left[
(b_{i,0}+b_{i,1})
-u_i(b_{i,0}-b_{i,1})
\right]\\
&\le
2m+2\sum_{i,a}|b_{i,a}|
\le10m.
\end{aligned}
\tag{6.16}
\]

Switching by \(s\) does not change either one-sided extremum of
\(C\).  In the one-copy normalization, (6.16) says

\[
P(C)\le\frac52m.
\tag{6.17}
\]

The one-sided discrepancy theorem recorded in
`one_sided_energy_product.md` uses the doubled normalization
\(x^\top Cx\).  Translating its equation (6) to the present one-copy
normalization gives

\[
P(C)\bigl(P(C)+N(C)\bigr)
\ge(1-o(1))\frac{m^3}{6400},
\tag{6.18}
\]

because \(P(C)=O(m)\) makes the density correction \(o(1)\).
Combining (6.17)--(6.18) yields

\[
N(C)
\ge(1-o(1))\frac{m^2}{16000}.
\tag{6.19}
\]

On the all-antiuniform pair modes, the principal energy is exactly

\[
H_{A[U]}((t_iv)_i)
=-m+4H_C(t).
\tag{6.20}
\]

Taking a negative ground state of \(C\) gives

\[
Q(A[U])
\ge m+4N(C)
\ge(1-o(1))\frac{m^2}{4000}.
\tag{6.21}
\]

Finally \(Q(A)\ge Q(A[U])\) by the random-extension argument used in
Proposition 6.1.  This proves (6.8). \(\square\)

Proposition 6.2 proves the conjectural \(Q\gtrsim n\,d\) scale when a
linear fraction of the affine directions is represented by repeated
two-vertex coordinate types.  The surviving case is not odd
leakage, but larger and highly unequal type multiplicities.

### Proposition 6.3 (odd type classes force additive fibre collisions)

Let

\[
\mathcal O=\{\tau:|V_\tau|\ \text{is odd}\}.
\]

If the affine energy is constant (ground-state maximality is not
needed), then for every nonzero \(\lambda\in D\),

\[
\boxed{
\#\bigl\{\{\phi,\tau\}\subset\mathcal O:
\phi+\tau=\lambda\bigr\}
\equiv0\pmod2.
}
\tag{6.22}
\]

#### Proof

The nonzero Fourier coefficient of the affine energy at
\(\lambda\) is

\[
0=
\sum_{\substack{\{\phi,\tau\}\\\phi+\tau=\lambda}}
I_{\phi,\tau},
\qquad
I_{\phi,\tau}
=
\sum_{\substack{i\in V_\phi\\j\in V_\tau}}a_{ij}.
\tag{6.23}
\]

Modulo two,

\[
I_{\phi,\tau}
\equiv |V_\phi||V_\tau|.
\]

Thus precisely the pairs for which both multiplicities are odd
contribute one modulo two to (6.23), proving (6.22). \(\square\)

In particular, the occupied odd types cannot form a Sidon set once
there are at least two of them: every represented pair-sum must occur
at least twice.  Four types forming an affine parallelogram realize
the smallest nontrivial pattern.  Hence the obstruction to extending
Proposition 6.2 to bounded odd multiplicities is necessarily
additive: odd blocks cannot be isolated character by character, and
their cancellations are routed through repeated pair-sum fibres.

This gives a clean bounded-multiplicity dichotomy but not yet a norm
theorem:

\[
\boxed{
\begin{array}{ll}
\text{many size-two types}
&\Longrightarrow Q(A)=\Omega(m^2)
\quad\text{by Proposition 6.2},\\[2mm]
\text{many odd types of size }\ge3
&\Longrightarrow\text{pair-sum collisions by (6.22)}.
\end{array}}
\tag{6.24}
\]

Turning the second line into a replacement requires showing that the
mod-two collisions contain *mixed signs* inside an induced block.
Equation (6.23) balances the total signed block sums, but by itself it
does not guarantee positive fibre replacement entropy (4.6).
Likewise, even type classes of size at least four allow nonzero row
and column modes of size up to \(|V_\tau|-1\), so the auxiliary
signing extracted in Proposition 6.2 becomes a weighted matrix with
linear fields.  These are the two precise obstructions to the
bounded-\(K\) generalization.

## 7. The sharpened missing condition

For affine cap clouds, the exact missing statistic is now explicit:
the Schur kernel controls duplicate-character predictability, whereas
a sign replacement requires integer freedom inside the character
fibres.

A viable inverse theorem therefore needs at least one of the
following additional conclusions.

1. **Fibre unlocking.**  Find an induced block \(T=E(S)\) for which
   \(\mathscr R_h(T)\) is large.  Equivalently, many character fibres
   must contain both signs of \(h\) inside \(T\), so signs can be
   exchanged while preserving the active affine cap profile.
2. **Envelope-aware slack.**  If \(\mathscr R_h(T)\) is small, use
   the actual slack \(W-|H_c(z)|\), not only covariance, to pay for
   the Fourier changes in (4.7).
3. **Local predictor control.**  Strengthen Schur predictability to a
   predictor supported near the vertex boundary of \(S\), with a
   controlled \(\ell_1\) or discrepancy norm.  Formula (3.2) permits
   prediction from a single arbitrarily remote copy and therefore
   has no relation to the boundary term in
   \[
   W(A[S])\le W_{|S|}
   +2\|A_{S,S^c}\|_{\infty\to1}.
   \]
4. **Competitive affine-ground rigidity.**  Prove that a
   linear-entropy affine cap cloud with low fibre replacement entropy
   necessarily produces an \(\Omega(n^2)\) Boolean witness, extending
   the paired-coordinate theorem (6.2) to general type
   multiplicities.  This would eliminate the obstruction at the
   \(O(n^{3/2})\) scale.

The fourth route is the cleanest model-case inverse theorem.  In type
language it asks for a quantitative dichotomy:

\[
\boxed{
\begin{array}{c}
\text{many independent affine cap directions}\\
+\text{mostly rigid induced character fibres}
\end{array}
\Longrightarrow
\begin{array}{c}
\text{large Boolean quadratic witness, or}\\
\text{a block with large fibre replacement entropy.}
\end{array}
}
\tag{7.1}
\]

This is strictly stronger than a covariance statement and is aligned
with the exact global replacement law (10.128): the relevant object
is the integer null lattice of the active cut-character evaluation
matrix, not merely its real kernel.

## 8. Bottom line

The flat conditional-kernel program remains viable only after the
following correction:

\[
\boxed{
\text{conditional }L^2\text{ predictability}
\ \text{must be supplemented by}\
\text{integer replacement entropy or cap-envelope control}.
}
\]

The affine calculation identifies precisely where a proof based only
on Schur complements would fail, supplies a proportional
high-entropy obstruction, and isolates the next rigorous model lemma:
extend paired-coordinate quadratic rigidity from repeated pairs to
arbitrary low-replacement-entropy type systems.
