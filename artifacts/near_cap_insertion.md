# Near-cap insertion: the full deficit hierarchy

## Status

This note keeps the entire near-cap hierarchy in the one-vertex
insertion problem.  It gives:

1. an exact reformulation in the cut gauge of one absolute ground
   state;
2. a deletion/replacement optimality theorem forced by global
   minimality of the signing;
3. exact entropy and matching/chaining criteria which control all
   energy layers simultaneously; and
4. a quantitative audit of the example in which one row balances all
   exact grounds but has insertion excess at least \(n-4\).

For globally minimizing signings it also derives a sharp final
dichotomy.  Either a deficit-two witness has enough signed boundary
traffic to recover the scale-correct \(3/2\) deletion coefficient, or
the favorable edges of a saturated rectangular block are covered by
low-traffic cap cuts.  At prime orders \(p\equiv1\pmod4\), parity
upgrades every cut in this cover to an exact top- or bottom-face zero
cut.

The audit is instructive: after using the near-cap hierarchy to choose
the row, the same example has insertion excess at most \(2\) for all
even block sizes \(k\ge 12\).  Thus that example obstructs endpoint-only
selection, but it is not an obstruction to hierarchy-aware insertion.

Throughout,
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_x|H_A(x)|,
\]
\[
g_A(x)=M(A)-|H_A(x)|,\qquad
\Delta_A(b)=\max_x\bigl(|b\cdot x|-g_A(x)\bigr).
\]

## 1. Exact cut-gauge form of every near-cap constraint

After replacing \(A\) by \(-A\), if necessary, and switching by an
absolute ground state, assume
\[
H_A(\mathbf 1)=M(A)=M.
\]
For \(S\subseteq[n]\), put
\[
c(S)=\sum_{\substack{i\in S\\j\notin S}}a_{ij}.
\]
If \(\mathbf1^S\) denotes the spin obtained from \(\mathbf1\) by
flipping precisely \(S\), then
\[
H_A(\mathbf1^S)=M-2c(S).
\]
The absolute-ground inequality gives
\[
\boxed{0\le c(S)\le M\quad(S\subseteq[n]).}
\tag{1.1}
\]
Consequently, with
\[
q(S)=\min\{c(S),M-c(S)\},
\]
the complete deficit hierarchy is
\[
\boxed{g_A(\mathbf1^S)=2q(S).}
\tag{1.2}
\]

Switch the proposed insertion row by the same ground state and still
call it \(b\).  Write \(b(S)=\sum_{i\in S}b_i\) and
\(B=b([n])\).  Then
\[
b\cdot\mathbf1^S=B-2b(S).
\]
It follows that
\[
\boxed{
\Delta_A(b)\le d
\iff
|B-2b(S)|\le d+2q(S)
\quad\hbox{for every }S\subseteq[n].
}
\tag{1.3}
\]

Thus thick-cap insertion is exactly a vertex-sign discrepancy problem
with the symmetric signed-cut penalty
\(\min\{c(S),M-c(S)\}\).  No passage from the exact endpoints to nearby
states is needed: (1.3) is the full problem.

## 2. What global optimality forces

Let \(A\) be a globally optimal signing of order \(n\), so
\(M(A)=M_n\).  Delete vertex \(i\), writing the principal signing as
\(B_i\) and the deleted row as \(r_i\).

### Proposition 2.1 (optimal deletion rows)

For every \(i\),
\[
\boxed{
\min_{b\in\{\pm1\}^{n-1}}\Delta_{B_i}(b)
=M_n-M(B_i),
}
\tag{2.1}
\]
and the original row \(r_i\) realizes the minimum.

#### Proof

Every row \(b\) extending \(B_i\) produces an order-\(n\) signing, so
its norm is at least \(M_n\).  The original row reconstructs \(A\) and
has norm \(M_n\).  The exact insertion identity
\[
\min_b\max_{x,y}|H_{B_i}(x)+y\,b\cdot x|
=M(B_i)+\min_b\Delta_{B_i}(b)
\]
therefore gives (2.1), including attainment by \(r_i\). \(\square\)

This has an exact covering consequence.  Put
\[
d_i=M_n-M(B_i).
\]
For every \(d<d_i\), the variable-radius caps
\[
\left\{
b:\ |b\cdot x|>d+g_{B_i}(x)
\right\},
\qquad x\in\{\pm1\}^{n-1},
\tag{2.2}
\]
cover the entire row cube.  Hence any proposed entropy or chaining
upper bound for the cap union must fail below \(d_i\).  This is a
necessary near-cap richness statement forced by global optimality,
not an assumption about ground-state entropy.

There is also a deletion budget.  Let \(x\) be a positive absolute
ground of \(A\), and set
\[
\ell_i=x_i\sum_{j\ne i}a_{ij}x_j.
\]
Single-spin optimality and the opposite absolute bound give
\[
0\le\ell_i\le M_n,\qquad \sum_i\ell_i=2M_n.
\tag{2.3}
\]
The energy of the restricted spin in \(B_i\) is \(M_n-\ell_i\), so
\[
d_i\le\ell_i,\qquad
\boxed{\sum_i d_i\le2M_n.}
\tag{2.4}
\]
The leading derivative needed for normalized monotonicity is only
\((3/2)M_n/n\).  Thus (2.4), whose average is \(2M_n/n\), displays the
precise factor \(4/3\) still missing from a deletion-only argument.

It is essential here that \(B_i\) need not be optimal at order
\(n-1\).  Define its nonoptimality defect
\[
e_i=M(B_i)-M_{n-1}\ge0.
\]
Then Proposition 2.1 gives the exact compensation identity
\[
\boxed{
e_i+d_i=M_n-M_{n-1}
\quad\text{for every }i.
}
\tag{2.4a}
\]
Consequently, averaging the \(d_i\)'s, even if it reached
\((3/2)M_n/n\), would not by itself bound the global increment:
a small hierarchy-aware reinsertion cost is compensated by a large
child defect.  A genuine two-step exchange theorem must control
\(e_i\) and \(d_i\) jointly after replacing \(B_i\) by a near-optimal
order-\(n-1\) signing.  This is a sharp obstruction to using (2.4) as
a scalar recurrence.

There is an exact metric formulation of the required exchange.  For
two signings \(B,C\) of the same order, define their deficit-profile
distance
\[
\eta(B,C)=\max_x|g_B(x)-g_C(x)|.
\tag{2.4b}
\]
Then, for every row \(r\),
\[
|\Delta_B(r)-\Delta_C(r)|\le\eta(B,C).
\tag{2.4c}
\]
Applying this with the deleted row \(r_i\) and an optimal signing
\(C\) of order \(n-1\) gives the two-step exchange bound
\[
\boxed{
M_n-M_{n-1}
\le d_i+\eta(B_i,C).
}
\tag{2.4d}
\]
Indeed, insert \(r_i\) into \(C\) and use Proposition 2.1.  A coarser
but sometimes easier metric is
\[
R(B,C)=\max_x|H_B(x)-H_C(x)|,
\]
for which
\[
\eta(B,C)\le2R(B,C).
\tag{2.4e}
\]
Thus the missing exchange theorem can be stated cleanly: after finding
a deletion with small \(d_i\), find an optimal child close to \(B_i\)
in the full deficit-profile metric.  Edge-Hamming closeness is far
stronger than needed, while equality of norms alone is far weaker.

Combining (2.4a) and (2.4d) also gives the forced separation
\[
\boxed{
\eta(B_i,C)\ge e_i
\quad\text{for every optimal order-\((n-1)\) signing }C.
}
\tag{2.4f}
\]
This simply records, in the correct metric, why reusing the deleted
row cannot secretly improve the global minimizer.  Any successful
exchange argument must therefore prove an *upper* stability estimate
on this distance (and hence on \(e_i\)); scalar closeness of
\(M(B_i)\) to \(M_{n-1}\) is not enough.

If \(A\) has both a positive and a negative absolute ground, with
local loss profiles \(\ell_i^+,\ell_i^-\), then
\[
d_i\le\min\{\ell_i^+,\ell_i^-\}
\tag{2.5}
\]
and therefore
\[
\boxed{
\sum_i d_i
\le
2M_n-\frac12\sum_i|\ell_i^+-\ell_i^-|.
}
\tag{2.6}
\]
In particular, an \(L_1\)-separation of at least \(M_n\) between the
two local loss profiles improves the average deletion loss to the
scale-correct leading derivative.  The unresolved complementary case
is that all opposite endpoint profiles are strongly aligned.

That aligned case has an exact block description.  Switch the positive
ground to \(\mathbf1\), let \(z\) be a negative ground, and put
\[
S=\{i:z_i=-1\},\qquad
A=
\begin{pmatrix}
B&C\\ C^\mathsf T&D
\end{pmatrix}
\tag{2.7}
\]
on \(S\sqcup S^c\).  Comparing the two endpoint energies gives
\[
\mathbf1^\mathsf B\mathbf1+
\mathbf1^\mathsf D\mathbf1=0,
\qquad
\mathbf1^\mathsf C\mathbf1=M_n.
\tag{2.8}
\]
For \(i\in S\), the two local losses are
\[
\ell_i^+=(B\mathbf1)_i+(C\mathbf1)_i,\qquad
\ell_i^-=-(B\mathbf1)_i+(C\mathbf1)_i,
\tag{2.9}
\]
and the analogous formula holds on \(S^c\), with \(D\) in place of
\(B\).  In particular,
\[
C\mathbf1\ge|B\mathbf1|,
\qquad
C^\mathsf T\mathbf1\ge|D\mathbf1|
\tag{2.10}
\]
coordinatewise, and
\[
\boxed{
\frac12\|\ell^+-\ell^-\|_1
=
\|B\mathbf1\|_1+\|D\mathbf1\|_1.
}
\tag{2.11}
\]
Thus (2.6) is equivalently
\[
\boxed{
\sum_i d_i
\le
2M_n-\|B\mathbf1\|_1-\|D\mathbf1\|_1.
}
\tag{2.12}
\]

Exact equality of the opposite local profiles is therefore equivalent
to
\[
B\mathbf1=D\mathbf1=0.
\tag{2.13}
\]
All endpoint energy then lies in a bipartite cross block with
nonnegative row and column sums, while both diagonal blocks are
row-Eulerian.  More quantitatively, if every deletion loss exceeds
the leading scale-correct value \(3M_n/(2n)\), then (2.12) forces
\[
\boxed{
\|B\mathbf1\|_1+\|D\mathbf1\|_1<\frac{M_n}{2}.
}
\tag{2.14}
\]
So failure to recover the missing \(3/2\) from opposite endpoint
profiles is confined to a concrete near-bipartite/Eulerian regime.

The cross block in (2.7) is in fact exactly saturated in the
\(\ell_\infty\!\to\!\ell_1\) norm.  For arbitrary signs \(u\) on
\(S\) and \(v\) on \(S^c\), compare the two full spins
\((u,v)\) and \((-u,v)\).  Their internal energy is unchanged and
their cross energy changes sign, so
\[
|H_B(u)+H_D(v)|+|u^\mathsf TCv|\le M_n.
\tag{2.15}
\]
Together with (2.8), this yields
\[
\boxed{\|C\|_{\infty\to1}
=\max_{u,v}|u^\mathsf TCv|
=\mathbf1^\mathsf TC\mathbf1
=M_n.}
\tag{2.16}
\]
Thus the aligned-profile wall is simultaneously a saturated
rectangular switching problem.  Any proof through (2.14) must use
more than second moments or the rectangular norm alone.

Global minimality also forces the first nontrivial near-cap layer in
this regime.

### Proposition 2.2 (favorable-edge replacement witnesses)

In the gauge (2.7), let \(e=ij\) be a cross edge with \(a_{ij}=+1\).
Then there is a spin \(w\) such that \(g_A(w)\le2\), and either
\[
H_A(w)\ge M_n-2,\quad w_iw_j=-1,
\tag{2.17}
\]
or
\[
H_A(w)\le-M_n+2,\quad w_iw_j=+1.
\tag{2.18}
\]

#### Proof

Flip only the coefficient \(a_{ij}\).  Both distinguished endpoint
energies move from \(\pm M_n\) to \(\pm(M_n-2)\).  The modified
signing still has norm at least \(M_n\), by global minimality.  A spin
witnessing that norm has modified energy
\[
H_A(w)-2w_iw_j.
\]
If this is at least \(M_n\), then \(w_iw_j=-1\) and
\(H_A(w)\ge M_n-2\).  If it is at most \(-M_n\), then
\(w_iw_j=+1\) and \(H_A(w)\le-M_n+2\). \(\square\)

Equivalently, every favorable cross edge is covered either by a
top-side cut \(T\) with \(c(T)\le1\) which separates its endpoints,
or by a bottom-side cut with \(c(T)\ge M_n-1\) which does not separate
them.  Since
\[
\#\{e\in S\times S^c:a_e=+1\}
=\frac{|S||S^c|+M_n}{2},
\tag{2.19}
\]
this is a large, globally forced deficit-two incidence system.  It is
the concrete near-cap structure that a hierarchy-aware insertion
argument may exploit in the aligned-profile branch.

At orders \(n\equiv1\pmod4\) this statement sharpens from near-cap to
exact-cap.  Indeed \(|T|(n-|T|)\) is even when \(n\) is odd, so the
signed cut sum \(c(T)\) is even.  Also
\(\binom n2\), and hence \(M_n\), is even when \(n\equiv1\pmod4\).
Thus \(0\le c(T)\le1\) forces \(c(T)=0\), while
\(M_n-1\le c(T)\le M_n\) forces \(c(T)=M_n\).  Therefore
\[
\boxed{
n\equiv1\pmod4
\Longrightarrow
\text{every witness in Proposition 2.2 can be chosen on an exact
top or bottom face.}
}
\tag{2.19b}
\]
This is relevant because convergence may be tested on the odd prime
orders \(p\equiv1\pmod4\).  On those orders, the remaining
low-traffic incidence cover is governed by exact zero-cut face
factorization, with no deficit-two approximation error.

The same parity observation also covers the case in which only one
orientation attains the absolute cap.  Switch a positive ground to
\(\mathbf1\).  For every edge with \(a_{ij}=+1\), flipping that
coefficient lowers the displayed top energy by two.  Minimality and
parity force an exact ground \(w\) of the original signing with either
\[
H_A(w)=M_n,\quad w_iw_j=-1,
\]
or
\[
H_A(w)=-M_n,\quad w_iw_j=+1.
\tag{2.19c}
\]
If the negative face does not attain \(-M_n\), only the first
alternative is possible.  Since a top gauge has
\[
\#\{ij:a_{ij}=+1\}
=\frac{\binom n2+M_n}{2}=\Theta(n^2),
\]
one-sided minimizers at prime orders also force a dense exact
zero-cut cover, now entirely inside the positive face.

This one-sided cover already gives a canonical block quotient.  Define
\[
i\sim_+j
\iff
w_i=w_j\quad\hbox{for every positive exact ground }w
\tag{2.19d}
\]
in the top gauge, and let \(V_1,\ldots,V_q\) be the equivalence
classes.  Every positive ground is constant on each \(V_a\).  More
strongly,
\[
\boxed{
\text{if the negative face is inactive, then }
a_{ij}=-1
\quad(i\ne j,\ i,j\in V_a).
}
\tag{2.19e}
\]
Indeed, a \(+1\) edge inside a class would have to be separated by a
positive exact ground, contradicting (2.19d).  Thus the diagonal
blocks of the canonical quotient are negative cliques.

This quotient is directly visible to insertion.  Choose the row signs
to sum to zero in every even class and to a residual \(\pm1\) in every
odd class.  Its field on every positive exact ground is at most the
number of odd classes.  Hence
\[
\#\{a:|V_a|\text{ odd}\}=o(\sqrt n)
\tag{2.19f}
\]
would already annihilate the one active exact face at the required
scale.  Failure leaves a quotient with \(\Omega(\sqrt n)\) odd
negative-clique types.  This is an explicit structured alternative,
although controlling its thick cap still requires the hierarchy
criteria of Sections 3--4.

The size of the endpoint cut gives a further clean split.  If
\(\min\{|S|,|S^c|\}=o(n)\), the larger diagonal block \(D\) is a
principal signing on \(n-o(n)\) vertices and
\[
M(D)\le M_n,\qquad
\frac{M_{|D|}}{|D|^{3/2}}
\le
\frac{M_n}{n^{3/2}}(1+o(1)).
\tag{2.19a}
\]
Thus this branch already supplies a scale-preserving principal
descent.  If both sides have linear size, (2.19) contains
\(\Theta(n^2)\) favorable edges.  Proposition 2.2 then forces a
quadratically large incidence cover by the first deficit-two layers.
Accordingly, the genuinely new case is a balanced saturated
rectangular block together with a dense low-traffic cap cover.

There is a quantitative high-traffic/low-traffic dichotomy behind
these witnesses.  For any spin \(x\), orient it so that
\[
sH_A(x)=M_n-g(x),\qquad s\in\{\pm1\},
\]
and define its oriented deletion-loss profile
\[
u_i(x)=M_n-sH_{A[-i]}(x_{-i}).
\tag{2.20}
\]
Single-spin comparison with the absolute cap gives \(u_i(x)\ge0\),
while
\[
\sum_i u_i(x)=2M_n+(n-2)g(x).
\tag{2.21}
\]
Moreover \(M(B_i)\ge sH_{A[-i]}(x_{-i})\), so
\[
d_i\le u_i(x).
\]
Consequently, for any two spins \(x,y\),
\[
\boxed{
\sum_i d_i
\le
\frac{
\sum_i u_i(x)+\sum_i u_i(y)
-\|u(x)-u(y)\|_1
}{2}.
}
\tag{2.22}
\]

Apply this to the top ground \(\mathbf1\) and a top-near spin
\(\mathbf1^T\) with \(c(T)\le1\).  Let
\[
\partial_T(i)
=
\sum_{\substack{j:\,|\{i,j\}\cap T|=1}}a_{ij}
\tag{2.23}
\]
be the signed boundary degree at vertex \(i\).  Direct calculation
gives
\[
u_i(\mathbf1^T)-u_i(\mathbf1)
=2c(T)-2\partial_T(i).
\tag{2.24}
\]
Using (2.21) in (2.22) yields the exact estimate
\[
\boxed{
\sum_i d_i
\le
2M_n+(n-2)c(T)
-\sum_i|\partial_T(i)-c(T)|.
}
\tag{2.25}
\]
Thus a deficit-two witness with
\[
\sum_i|\partial_T(i)-c(T)|
\ge\frac{M_n}{2}+O(n)
\tag{2.26}
\]
recovers a deletion/reinsertion cost
\[
\min_i d_i\le\frac{3M_n}{2n}+O(1),
\tag{2.27}
\]
which is accurate enough at the summable-error scale.  If (2.26)
fails for every witness supplied by Proposition 2.2 (in the
appropriate top or bottom gauge), all favorable cross edges are
covered by deficit-two cuts whose signed boundary-degree traffic is
less than \(M_n/2+O(n)\).  This is precisely the non-affine
low-traffic inverse branch, now derived from global optimality.

As emphasized by (2.4a), (2.27) is a scale-correct step relative to
the particular child \(B_i\); a global recurrence additionally needs
to control its defect \(e_i\).  The new content is the sharp
dichotomy: high boundary traffic gives the \(3/2\) coefficient,
whereas failure forces a large low-traffic incidence cover.

## 3. An exact balanced-row entropy criterion

Assume \(n\) is even and retain the cut gauge of Section 1.  Choose
\(b\) uniformly from the balanced slice
\[
\mathcal B_n=\{b\in\{\pm1\}^n:b([n])=0\}.
\]
For \(|S|=s\), the number \(K\) of \(+1\) coordinates of \(b\) in
\(S\) has the hypergeometric law
\[
\Pr(K=k)
=
\frac{\binom sk\binom{n-s}{n/2-k}}{\binom n{n/2}},
\]
and
\[
b\cdot\mathbf1^S=-2b(S)=-2(2K-s).
\]
Define the exact tail
\[
\Pi_{n,s}(t)
=
\frac1{\binom n{n/2}}
\sum_k
\binom sk\binom{n-s}{n/2-k}
\mathbf1_{\{\,2|2k-s|>t\,\}}.
\tag{3.1}
\]

### Proposition 3.1 (balanced cap-union criterion)

If
\[
\boxed{
\sum_{S\subseteq[n]}
\Pi_{n,|S|}\bigl(d+2q(S)\bigr)<1,
}
\tag{3.2}
\]
then some balanced row \(b\) satisfies \(\Delta_A(b)\le d\).

#### Proof

For each \(S\), (3.1) is exactly the probability that its inequality
in (1.3) fails.  The union bound leaves positive probability that no
inequality fails. \(\square\)

This criterion is finite and directly checkable from the joint
cut-size/deficit census.  It does not require few exact or near-exact
grounds: an exponentially large layer is allowed whenever its
hypergeometric large-deviation cost is larger than its entropy.
Conversely, Proposition 2.1 says that for every deletion of a global
minimizer, (3.2) must fail for all \(d<d_i\).

A dyadic version only needs the counts
\[
N_{s,j}
=
\#\{S:|S|=s,\ 2^j\le q(S)<2^{j+1}\},
\tag{3.3}
\]
together with the separate \(q=0\) layer.  Monotonicity of
\(\Pi_{n,s}\) turns (3.2) into an explicit entropy-versus-tail
inequality involving the \(N_{s,j}\).

## 4. Matching compression and a chaining target

The balanced-slice criterion still counts highly correlated cuts
separately.  A perfect matching compresses those correlations.
Assume \(n=2m\), fix a perfect matching \(P\), and give the two
endpoints of every matching edge opposite signs.  Its orientation is
encoded by \(\varepsilon\in\{\pm1\}^m\).

For a subset \(S\), define
\[
\sigma_P(S)\in\{0,\pm1\}^m
\]
by putting zero on a matching edge if \(S\) contains both or neither
endpoint, and putting \(+1\) or \(-1\) according to which endpoint is
in \(S\).  Then
\[
b([n])=0,\qquad b(S)=\langle\varepsilon,\sigma_P(S)\rangle.
\tag{4.1}
\]
For a realized signature \(\sigma\), retain only the most stringent
deficit
\[
q_P(\sigma)
=
\min\{q(S):\sigma_P(S)=\sigma\}.
\tag{4.2}
\]

### Proposition 4.1 (oriented-matching criterion)

If
\[
\boxed{
\sum_{\substack{\sigma\in\Sigma_P\\
                 \|\sigma\|_0>0}}
2\exp\left(
-\frac{(d+2q_P(\sigma))^2}{8\|\sigma\|_0}
\right)<1,
}
\tag{4.3}
\]
where \(\Sigma_P=\{\sigma_P(S):S\subseteq[n]\}\), then some orientation
of \(P\) gives an insertion row with \(\Delta_A(b)\le d\).

#### Proof

For fixed \(\sigma\), Hoeffding's inequality gives
\[
\Pr_\varepsilon\left(
2|\langle\varepsilon,\sigma\rangle|>d+2q_P(\sigma)
\right)
\le
2\exp\left(
-\frac{(d+2q_P(\sigma))^2}{8\|\sigma\|_0}
\right).
\]
Signatures with zero support never violate the constraint.  Apply the
union bound over distinct signatures, rather than over all cuts.
\(\square\)

There is a stronger deterministic but more restrictive certificate:
write
\[
m_P(S)=\#\{e\in P:|e\cap S|=1\}.
\]
Then
\[
\boxed{
m_P(S)\le \frac d2+q(S)\quad(S\subseteq[n])
\Longrightarrow
\Delta_A(b)\le d
}
\tag{4.4}
\]
for every orientation of \(P\).  Indeed
\(|b(S)|\le m_P(S)\).  Thus the search for a hierarchy-aware row
contains a concrete thin-matching problem relative to the symmetric
cut penalty \(q\).

One can replace the union bound in (4.3) by generic chaining.  For
layers
\[
\Sigma_j
=
\{\sigma:q_P(\sigma)\in[s_j,s_{j+1})\},
\quad
w_j
=
\mathbb E_\varepsilon
\sup_{\sigma\in\Sigma_j}
|\langle\varepsilon,\sigma\rangle|,
\tag{4.5}
\]
bounded differences gives
\[
\Pr\left\{
\sup_{\sigma\in\Sigma_j}
|\langle\varepsilon,\sigma\rangle|
>
\frac{d+2s_j}{2}
\right\}
\le
\exp\left(
-\frac{\bigl((d+2s_j)/2-w_j\bigr)^2}{2m}
\right)
\tag{4.6}
\]
whenever the displayed gap is positive.  Hence summing the
right-hand side of (4.6) over the layers gives another sufficient
criterion.  Each \(w_j\) is controlled by the metric entropy of
\(\Sigma_j\) in the canonical metric
\(\|\sigma-\tau\|_2\), so (4.6) is a precise, checkable chaining
property whose proof would settle the thick-cap transfer.  It uses
overlap geometry and does not assume low cardinality of a ground
layer.

For odd \(n\), a near-perfect matching leaves one coordinate
unmatched.  The same proof has the deterministic bound
\[
|b([n])-2b(S)|\le1+2m_P(S),
\]
so only an additive \(1\) is lost.

## 5. The \(n-4\) endpoint counterexample is repaired by its hierarchy

Use the construction in `cap_discrepancy_insertion.md`.  Thus
\(n=2k\), \(k\) is even, \(B\) is the balanced circulant cross block,
\[
\rho=\|B\|_{\rm op}=\frac2{\sin(\pi/k)},
\]
and one cross coefficient of
\[
A_0=
\begin{pmatrix}
J_k-I_k&B\\
B^\mathsf T&J_k-I_k
\end{pmatrix}
\]
is flipped to obtain \(A\).  Its absolute norm is
\[
M(A)=k(k-1)+2.
\]
The endpoint-only row exhibited there has zero field on the unique
absolute ground but satisfies
\[
\Delta_A(b)\ge2k-4=n-4.
\]

Now choose instead
\[
b=(u,v),
\qquad
u,v\in\{\pm1\}^k,\quad
u\cdot\mathbf1=v\cdot\mathbf1=0.
\tag{5.1}
\]
This row balances both block-constant orientations, including the
slack-four state which defeated the endpoint-only row.

For arbitrary spins \(x,y\), put
\[
s=\mathbf1^\mathsf Tx,\quad
t=\mathbf1^\mathsf Ty,\quad
X=k^2-s^2,\quad Y=k^2-t^2.
\]
Orthogonality to the constants gives
\[
|u\cdot x+v\cdot y|\le\sqrt X+\sqrt Y.
\tag{5.2}
\]
If \(H_A(x,y)\ge0\), the spectral estimate for \(A_0\), together with
the fact that the coefficient flip changes energy by at most two,
gives
\[
g_A(x,y)
\ge
\frac{X+Y}{2}-\frac{\rho}{k}\sqrt{XY}
\ge
\frac{1-\rho/k}{2}(X+Y).
\tag{5.3}
\]
Therefore
\[
|b\cdot(x,y)|-g_A(x,y)
\le
\frac1{1-\rho/k}.
\tag{5.4}
\]
If \(H_A(x,y)<0\), the lower spectral bound gives
\[
g_A(x,y)\ge k(k-2-\rho),
\]
and hence
\[
|b\cdot(x,y)|-g_A(x,y)
\le
2k-k(k-2-\rho).
\tag{5.5}
\]
Combining (5.4) and (5.5),
\[
\boxed{
\Delta_A(b)
\le
\max\left\{
\frac1{1-\rho/k},
\ 2k-k(k-2-\rho)
\right\}.
}
\tag{5.6}
\]

For every even \(k\ge12\), one has \(\rho\le k-4\), so the second term
is nonpositive, while the first is less than \(3\).  Every field and
deficit in this even-order example is even.  Thus
\[
\boxed{
\min_b\Delta_A(b)\le2
\qquad(k\ge12\text{ even}).
}
\tag{5.7}
\]

The same signing therefore has a ground-balancing row with excess at
least \(n-4\) and hierarchy-aware balanced rows with excess at most
\(2\).  This sharply isolates the lesson of the example: the
near-cap hierarchy must participate in choosing the row, but it can
also completely remove the apparent obstruction.

## 6. Remaining theorem

The direct route would close if, for every globally minimizing
signing \(A\) of order \(n\), one could find a matching \(P\), or more
generally a balanced row ensemble, for which the chained criterion of
Section 4 holds with
\[
d
\le
M_n\left[\left(1+\frac1n\right)^{3/2}-1\right]+r_n,
\qquad
\sum_n\frac{r_n}{n^{3/2}}<\infty.
\tag{6.1}
\]
Proposition 2.1 supplies the exact optimality input: every deleted row
of a global minimizer is itself optimal for the full weighted
near-cap problem of its child, and the corresponding variable-radius
caps cover the row cube below its optimum.  What is still missing is
an inverse theorem converting such repeated covering, across all
deletions, into either:

1. the entropy/chaining bound (4.6) for the parent hierarchy; or
2. a principal restriction retaining the normalized optimum at the
   scale-correct \(3/2\) rate.

The factor \(4/3\) gap in (2.4) and the aligned-profile alternative in
(2.6) are the precise numerical and structural barriers.

The strongest form of the remaining aligned branch is now more
specific.  On orders \(p\equiv1\pmod4\):

1. the opposite endpoint partition has a saturated cross block
   \(C\) with \(\|C\|_{\infty\to1}=M_p\);
2. unless (2.27) holds, every favorable edge of \(C\) is separated by
   a low-traffic exact zero cut in one of the two endpoint gauges; and
3. every such zero cut \(R\) has exact Cartesian principal-ground
   factorization and
   \[
   \operatorname{rank}A_{R,R^c}
   \le
   \operatorname{codim}L_R+\operatorname{codim}L_{R^c},
   \tag{6.2}
   \]
   where the \(L\)'s are the two principal positive-ground spans.

Thus near-full component ground spans produce a low-rank sign block,
and hence a finite row/column-pattern quotient; large cross rank
forces a large ground-span deficit.  What is not yet proved is that
the many potentially crossing zero cuts in the favorable-edge cover
admit one compatible laminar refinement or bounded-depth quotient.
That is the exact block-quotient target usable by an amplification
argument.

Nor does zero-cut factorization alone control the child defect
\(e_i\).  Its inequalities compare principal pieces of the *same*
signing and give lower/additive width information, while \(e_i\)
compares \(B_i\) with a separately optimized order-\((n-1)\) signing.
The profile-distance separation (2.4f) shows that an additional
stability theorem across optimal signings is logically necessary.
This is the clean stopping frontier for the direct insertion route.
