# Ground-state entropy under products, blow-ups, and twin substitutions

## 1. Exact seed audit

Use

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
M(A)=\max_x|H_A(x)|.
\]

After fixing the first row by switching, an order-five minimizer is

\[
A_5=
\begin{pmatrix}
0&1&1&1&1\\
1&0&-1&1&-1\\
1&-1&0&-1&1\\
1&1&-1&0&1\\
1&-1&1&1&0
\end{pmatrix}.
\]

It has

\[
M(A_5)=4,
\qquad
\#\{x:|H_{A_5}(x)|=4\}=20.
\]

Its complete energy histogram is

\[
\#\{H=0\}=12,\qquad
\#\{H=4\}=10,\qquad
\#\{H=-4\}=10.
\]

Every positive ground state has the same sorted oriented-local-field
profile

\[
(0,0,2,2,4).
\tag{1.1}
\]

Thus the large finite degeneracy is accompanied by exactly two
zero-field coordinates.

An order-six minimizer is

\[
A_6=
\begin{pmatrix}
0&1&1&1&1&1\\
1&0&1&1&-1&-1\\
1&1&0&-1&1&-1\\
1&1&-1&0&-1&1\\
1&-1&1&-1&0&1\\
1&-1&-1&1&1&0
\end{pmatrix}.
\]

It has

\[
M(A_6)=5,
\qquad
\#\{x:|H_{A_6}(x)|=5\}=24,
\]

and energy histogram

\[
\#\{H=3\}=\#\{H=-3\}=20,\qquad
\#\{H=5\}=\#\{H=-5\}=12.
\]

Every positive ground state has oriented-local-field profile

\[
(1,1,1,1,1,5).
\tag{1.2}
\]

In particular it is strict under every single spin flip.  There are
ten balanced ground states, all at energy \(-5\).

These counts were re-enumerated directly over all spins and all
gauge-fixed minimizers: every order-five minimizer has ground count
20, and every order-six minimizer has ground count 24.

## 2. Constant-block blow-ups preserve the scale but not the grounds

Let \(S\) be an order-\(r\) signing and \(B\) an order-\(m\) signing.
The standard rank-one lexicographic substitution puts \(B\) on each
diagonal block and \(s_{ab}J_m\) on off-diagonal block \(ab\).
For fiber spins \(u^{(a)}\), let

\[
z_a=\sum_{i=1}^m u_i^{(a)}.
\]

Then its Hamiltonian is exactly

\[
H_G(u^{(1)},\ldots,u^{(r)})
=
\sum_{a=1}^rH_B(u^{(a)})
+\sum_{a<b}s_{ab}z_az_b.
\tag{2.1}
\]

Since the last expression is multi-affine in
\((z_1,\ldots,z_r)\in[-m,m]^r\), its cross part has maximum

\[
m^2M(S)
\tag{2.2}
\]

at corner magnetizations \(z_a=\pm m\).  Hence

\[
m^2M(S)-rM(B)
\le M(G)\le
m^2M(S)+rM(B).
\tag{2.3}
\]

If \(r\) is fixed and \(M(B)=O(m^{3/2})\), the construction has
\(M(G)=O((rm)^{3/2})\) only because its leading term is the rank-one
cross energy.  Balanced ground states of the order-six seed have
zero cross energy in (2.1), only \(O(r)\) internal energy, and therefore
are not grounds of a large blow-up: corner fibers create the unavoidable
\(\Theta(m^2)\) leading energy.

This rules out the simplest attempt to multiply the ten balanced
order-six grounds across blocks.

## 3. A sign block has no hidden saturated microstate

The natural repair would be to replace \(J_m\) by a sign block \(R\)
having several Boolean microstates which are externally
indistinguishable and all saturate the cross norm.  This is impossible.

### Theorem 3.1

Let \(R\in\{\pm1\}^{m\times m}\) and

\[
L=\max_{x,y\in\{\pm1\}^m}|x^\top Ry|.
\]

Suppose \(U\subseteq\{\pm1\}^m\) satisfies

\[
u^\top Rv=L
\qquad\text{for every }u,v\in U.
\tag{3.1}
\]

Then

\[
\boxed{|U|=1.}
\]

#### Proof

Assume that \(u,v\in U\) are distinct.  Gauge both sides by \(u\), so
that \(u=\mathbf1\).  Let

\[
F=\{i:v_i=-1\},\qquad C=[m]\setminus F.
\]

Write the sums of the four rectangular blocks of \(R\) as

\[
a=\mathbf1_F^\top R_{FF}\mathbf1_F,\quad
b=\mathbf1_F^\top R_{FC}\mathbf1_C,\quad
c=\mathbf1_C^\top R_{CF}\mathbf1_F,\quad
d=\mathbf1_C^\top R_{CC}\mathbf1_C.
\]

The four equalities

\[
\mathbf1^\top R\mathbf1
=
\mathbf1^\top Rv
=
v^\top R\mathbf1
=
v^\top Rv
=L
\]

give

\[
a=b=c=0,\qquad d=L.
\tag{3.2}
\]

Now keep all coordinates in \(C\) equal to \(1\), but choose arbitrary
Boolean vectors \(p,q\) on \(F\).  The resulting bilinear value is

\[
L+
p^\top R_{FF}q+
p^\top R_{FC}\mathbf1_C+
\mathbf1_C^\top R_{CF}q.
\tag{3.3}
\]

Because \(L\) is the global positive maximum, the last three terms in
(3.3) are at most zero for every \(p,q\).  Their uniform average over
\((p,q)\) is zero, so they must vanish identically.  But the Fourier
coefficient of \(p_iq_j\) is \((R_{FF})_{ij}\in\{\pm1\}\), a
contradiction.  \(\square\)

The theorem is valid without symmetry.  Thus an off-diagonal sign
gadget cannot carry even one hidden degeneracy bit if every pair of
microstates is required to realize its extremal logical interaction.
The brute-force search through all blocks of sizes \(m\le4\) found no
exception, as the theorem requires.

## 4. At most one fiber can remain macroscopically free

There is a complementary continuous obstruction which does not assume
pairwise saturation.

### Lemma 4.1

Let

\[
q(z)=\sum_{i<j}s_{ij}z_iz_j,
\qquad z\in[-1,1]^r,
\]

with every \(s_{ij}\ne0\).  At a global maximum or minimum of \(q\), at
most one coordinate can lie strictly inside \((-1,1)\).

#### Proof

If two coordinates \(i,j\) are interior, first-order optimality gives
\(\partial_iq=\partial_jq=0\).  Perturb only these coordinates.  The
linear terms vanish and the change is

\[
s_{ij}\delta_i\delta_j.
\]

Choosing the sign of \(\delta_i\delta_j\) increases \(q\), or decreases
it in the minimum case, a contradiction.  \(\square\)

There is also a quantitative probabilistic form.  Let \(M\) be the
maximum of \(q\) on the cube, and suppose

\[
q(z)\ge M-\varepsilon.
\]

Round coordinates independently to signs \(Z_i\) with
\(\mathbb E Z_i=z_i\).  Since \(q\) is multi-affine,

\[
\mathbb E q(Z)=q(z).
\]

The cube energies are integral and a nonground corner has gap at least
two.  Hence

\[
\Pr[Z\text{ is not a ground corner}]\le\varepsilon/2.
\tag{4.1}
\]

Put \(p_i=(1-|z_i|)/2\), the minority probability in coordinate \(i\).
For any pair \(i\ne j\), after all other coordinates are fixed, not all
four assignments of \((Z_i,Z_j)\) can be ground corners: their mixed
second difference is \(4s_{ij}\ne0\).  Therefore

\[
\boxed{p_ip_j\le\varepsilon/2\qquad(i\ne j).}
\tag{4.2}
\]

All but at most one coordinate consequently have

\[
p_i\le\sqrt{\varepsilon/2}.
\tag{4.3}
\]

Apply this to (2.1).  If \(M(B)=O(m^{3/2})\), any ground state of the
substitution must have cross energy within \(O(m^{3/2})\) of the
leading \(m^2M(S)\).  Thus \(\varepsilon=O(m^{-1/2})\) after dividing
magnetizations by \(m\).  In all but at most one fiber, the minority
spin fraction is \(O(m^{-1/4})\).  The number of possibilities in those
fibers is

\[
\exp\!\bigl(O(m^{3/4}\log m)\bigr)=e^{o(m)}.
\tag{4.4}
\]

Only one fiber can carry macroscopic entropy.  A fixed
lexicographic substitution therefore cannot multiply an old
ground-state entropy independently across a positive fraction of its
fibers.

For recursive substitution by a fixed order-\(r\) seed, this means
that any entropy inherited from the previous level appears in at most
one of the \(r\) new fibers.  Its entropy density is divided by \(r\),
up to the \(o(1)\) contribution in (4.4).  The high finite degeneracies
of \(A_5\) and \(A_6\) consequently do not bootstrap to a positive
asymptotic ground-state entropy through this product class.

## 5. Why linear twin towers also miss the scale

Duplicating one zero-field vertex rather than every vertex avoids the
\(\Theta(m^2)\) interaction between large fibers, and it can preserve a
constant degeneracy factor at a single step.  But the new twin has a
full sign row to the existing \(m\)-vertex core.  Against unrestricted
spins that row has a Boolean field as large as \(m\).  Repeating a
constant-size twin shell for a linear number of steps therefore has a
potential cumulative cost

\[
\sum_{j\le n}\Theta(j)=\Theta(n^2),
\]

not \(O(n^{3/2})\).  Restricting to \(O(\sqrt n)\) such steps keeps the
quadratic norm at the right scale but creates only
\(\exp(O(\sqrt n))\) states.

Turning this heuristic into a lower bound for every adaptive twin
tower would require a joint discrepancy theorem for the inherited
ground family.  The existing max-plus insertion-tower obstruction is
exactly this missing point.  Thus the twin route is not fully ruled out,
but it no longer offers a free entropy amplification.

## 6. Conference benchmark

For comparison, the symmetric Paley conference matrix of order ten
(constructed over \(\mathbb F_9\)) satisfies

\[
C^2=9I.
\]

Direct enumeration gives exactly

\[
12\ \text{Boolean }(+3)\text{-eigenvectors}
\quad\text{and}\quad
12\ \text{Boolean }(-3)\text{-eigenvectors}.
\]

Hence it has 24 spectral-cap ground states out of \(2^{10}\), matching
the order-six seed count but giving no evidence of exponential
eigenspace intersection.  The hidden-saturated-microstate theorem
explains why a fixed sign-block product cannot simply multiply these
Boolean eigenvectors.

## 7. Verdict

No infinite \(e^{cn}\)-ground-state family with
\(M=O(n^{3/2})\) emerged.

The following broad mechanisms are now ruled out:

1. constant rank-one blow-ups of the balanced \(A_6\) grounds;
2. any fixed sign-block product whose microstates pairwise saturate the
   block's rectangular norm;
3. recursive fixed lexicographic substitutions which attempt to make
   many fibers independently free.

The only survivor is an adaptive linear twin/insertion tower in which
each new row has exceptionally small discrepancy simultaneously on the
entire inherited ground family.  Producing such a tower would itself
be the sought Sidon-like construction; proving it impossible is
equivalent to a strong ground-family discrepancy theorem, not to an
ordinary graph-product calculation.
