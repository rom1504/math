# Capped-bilinear inverse: cyclic monotonicity and cubic-defect localization

## Status

This note studies the exact equality profiles left by the
zero-cut principal closure in `opposite_face_rigidity.md`.

The main proved outputs are:

1. equality profiles form a cyclically monotone transport support for
   the cross-block bilinear form;
2. every profile pays not only its two internal deficits, but also
   its interaction with the complete Cartesian family of principal
   ground states;
3. the global cubic product-closure defect has an exact
   disintegration into these capped-bilinear payments;
4. if the cap correlation has diverging effective rank, a
   macroscopic balanced zero cut carries any substantial cubic
   defect.

This gives the requested macroscopic cut and an
\(n^{3/2}\)-scale cross-block certificate.  It does **not** yet close
the convergence proof, because the exact centered-width bound takes
the maximum of the child-width budget and the cross payment, rather
than their sum.

Throughout, write a zero-cut principal split as
\[
A=\begin{pmatrix}B&C\\ C^{\mathsf T}&D\end{pmatrix},
\qquad
p(A)=p(B)+p(D),
\]
where the all-one vector is a positive ground state in both
principal blocks.  Put
\[
f(u)=p(B)-H_B(u),\qquad
g(v)=p(D)-H_D(v).
\tag{1}
\]
Then
\[
|u^{\mathsf T}Cv|\le f(u)+g(v).
\tag{2}
\]

## 1. Equality profiles are cyclically monotone

After changing \(v\) to \(-v\) when necessary, define the positive
equality set
\[
\mathcal E
=
\{(u,v):u^{\mathsf T}Cv=f(u)+g(v)\}.
\tag{3}
\]
These are exactly the positive full ground states in the chosen
block orientation.

For arbitrary \(u,v\), define the nonnegative slack
\[
s(u,v)=f(u)+g(v)-u^{\mathsf T}Cv.
\tag{4}
\]

### Proposition 1.1 (two-cycle monotonicity)

For any \((u,v),(u',v')\in\mathcal E\),
\[
\boxed{
(u-u')^{\mathsf T}C(v-v')
=s(u,v')+s(u',v)\ge0.
}
\tag{5}
\]
Moreover, equality in (5) holds if and only if both crossed pairs
\((u,v')\) and \((u',v)\) also belong to \(\mathcal E\).

More generally, for every cyclic permutation of
\((u_1,v_1),\ldots,(u_k,v_k)\in\mathcal E\),
\[
\boxed{
\sum_{i=1}^k u_i^{\mathsf T}Cv_i
\ge
\sum_{i=1}^k u_i^{\mathsf T}Cv_{i+1}.
}
\tag{6}
\]

#### Proof

For (5), expand its left side and replace the two diagonal terms by
\(f(u)+g(v)\) and \(f(u')+g(v')\).  The result is exactly the two
slacks in (4).  Their nonnegativity proves the assertion and also
the equality characterization.  Summing (4) around a cycle proves
(6). \(\square\)

Thus a zero-gap clique of profiles automatically completes to a
Cartesian rectangle of full ground states.  On such a rectangle,
all mixed differences are annihilated by \(C\).  The two-cycle gap
is the precise nonnegative obstruction to the affine/quotient
branch.

## 2. Interaction with the principal Cartesian ground family

Let \(\mathcal G_B,\mathcal G_D\) be the positive principal
ground-state sets.  Proposition 5.2 of
`opposite_face_rigidity.md` says that
\[
\mathcal G_B\times\mathcal G_D\subseteq\mathcal E,
\qquad
u_0^{\mathsf T}Cv_0=0
\quad(u_0\in\mathcal G_B,\ v_0\in\mathcal G_D).
\tag{7}
\]
Both ground families are invariant under global sign.

### Proposition 2.1 (projection payment)

Every \((u,v)\in\mathcal E\) and every
\(u_0\in\mathcal G_B,\ v_0\in\mathcal G_D\) obey
\[
\boxed{
f(u)+g(v)=u^{\mathsf T}Cv
\ge
|u^{\mathsf T}Cv_0|+|u_0^{\mathsf T}Cv|.
}
\tag{8}
\]

#### Proof

Apply (5) to \((u,v)\) and
\((\varepsilon u_0,\delta v_0)\), where
\(\varepsilon,\delta\in\{\pm1\}\) are independent choices.  Using
(7), the resulting inequality is
\[
u^{\mathsf T}Cv
\ge
\delta\,u^{\mathsf T}Cv_0
+\varepsilon\,u_0^{\mathsf T}Cv.
\]
Maximizing over the two signs proves (8). \(\square\)

Let \(R_B,R_D\) be covariance matrices of arbitrary symmetric laws
on the two principal ground families, and let \(P_B,P_D\) project
onto their spans.  If
\[
R_B\succeq \kappa_B P_B,\qquad
R_D\succeq \kappa_D P_D
\tag{9}
\]
on those spans, then (8) and root-mean-square domination give
\[
\boxed{
f(u)+g(v)
\ge
\sqrt{\kappa_D}\,\|P_D C^{\mathsf T}u\|_2
+
\sqrt{\kappa_B}\,\|P_B Cv\|_2.
}
\tag{10}
\]
This is a quantitative version of cross annihilation: an equality
profile outside the two projected kernels must spend internal
deficit.  The remaining limitation is explicit: without a frame
lower bound, span dimension alone gives no useful \(\kappa\).

## 3. Cubic defect is exactly capped-bilinear payment

Let \(\mu\) be a symmetric law supported on the positive ground face
of \(A\), with correlation matrix
\[
R=\mathbb E_\mu XX^{\mathsf T}.
\tag{11}
\]
Define its cubic closure defect
\[
\Delta_3(\mu)
=
p(A)-\mathbb E H_A(XYZ)
=
p(A)-\langle A,R^{\circ3}\rangle,
\tag{12}
\]
where \(X,Y,Z\) are independent with law \(\mu\), and the matrix
inner product uses the same off-diagonal normalization as \(H_A\).

For fixed \(X,Y\), let
\[
S=S(X,Y)=\{i:X_iY_i=-1\}.
\tag{13}
\]
After switching by \(X\), \(S\) is a zero cut because both \(X\)
and \(Y\) are positive grounds.  Write the switched restriction of
\(Z\) as \((u,v)\) on \(S,S^c\), and let \(C_S\) be the switched
cross block.

### Proposition 3.1 (exact defect disintegration)

For every fixed \(X,Y\) and every positive ground \(Z\),
\[
\boxed{
p(A)-H_A(XYZ)
=2u^{\mathsf T}C_Sv
=2\bigl[f_S(u)+g_S(v)\bigr].
}
\tag{14}
\]
Consequently, if
\[
\pi(S)
=
\mathbb E_Z\,u^{\mathsf T}C_Sv
=
\sum_{i\in S,\ j\notin S}a_{ij}R_{ij},
\tag{15}
\]
then
\[
\boxed{
\pi(S)\ge0,\qquad
\Delta_3(\mu)=2\mathbb E_{X,Y}\pi(S(X,Y)).
}
\tag{16}
\]

#### Proof

In the switched gauge, \(X=(\mathbf1,\mathbf1)\) and
\(Y=(-\mathbf1,\mathbf1)\).  The positive ground \(Z=(u,v)\)
saturates (2), so
\[
u^{\mathsf T}C_Sv=f_S(u)+g_S(v)\ge0.
\]
The product \(XYZ=(-u,v)\) reverses only the cross term and leaves
both internal energies unchanged.  Its deficit is therefore twice
the displayed cross term, proving (14).  Averaging over \(Z\), and
then over \(X,Y\), proves (15)--(16). \(\square\)

There is an equivalent edgewise check:
\[
\mathbb P\bigl(ij\in\delta(S(X,Y))\bigr)
=\frac{1-R_{ij}^2}{2},
\]
so the right side of (16) is
\(\sum_{i<j}a_{ij}R_{ij}(1-R_{ij}^2)\), exactly (12).

## 4. Substantial cubic defect localizes on a balanced zero cut

Define the correlation effective rank
\[
r_{\rm eff}(R)=\frac{(\operatorname{tr}R)^2}{\|R\|_F^2}
=\frac{n^2}{\|R\|_F^2}.
\tag{17}
\]
Because \(\mu\) is symmetric,
\[
|S(X,Y)|=\frac{n-X^{\mathsf T}Y}{2},
\qquad
\mathbb E(X^{\mathsf T}Y)^2=\|R\|_F^2.
\tag{18}
\]
Thus
\[
\mathbb P\left(\left||S|-\frac n2\right|>\varepsilon n\right)
\le
\frac{1}{4\varepsilon^2r_{\rm eff}(R)}.
\tag{19}
\]

Every payment in (15) is capped by the full centered width:
\[
0\le\pi(S)\le W(A).
\tag{20}
\]
Indeed, each integrand is an equality value
\(f_S(u)+g_S(v)\), and the exact range formula below gives
\(f_S(u)+g_S(v)\le W(A)\).

### Corollary 4.1 (macroscopic payment cut)

Suppose
\[
W(A)\le Cn^{3/2},\qquad
r_{\rm eff}(R)\longrightarrow\infty,\qquad
\Delta_3(\mu)\ge\delta n^{3/2}.
\tag{21}
\]
Then, for every fixed \(0<\varepsilon<1/2\) and all sufficiently
large \(n\), there is a zero cut \(S\) generated by two support states
such that
\[
\varepsilon n\le |S|\le(1-\varepsilon)n
\tag{22}
\]
and
\[
\boxed{
\pi(S)\ge\frac{\delta}{2}n^{3/2}-o(n^{3/2}).
}
\tag{23}
\]
For a square-root-scale effective rank, the exceptional contribution
in (19)--(20) is already \(O(n^{5/4})\).

## 5. The exact normalization wall

For the zero-cut split, put
\[
R_B=p(B)+\nu(B)=2W(B),\qquad
R_D=p(D)+\nu(D)=2W(D).
\]
A direct maximization over the two block orientations gives
\[
\boxed{
2W(A)
=
\max_{u,v}
\left\{
f(u)+g(v)+|u^{\mathsf T}Cv|
\right\}.
}
\tag{24}
\]
It follows that
\[
\boxed{
W(A)\ge W(B)+W(D),
\qquad
W(A)\ge\|C\|_{\infty\to1},
\qquad
W(A)\ge\pi(S).
}
\tag{25}
\]
The second inequality uses (2): at a bilinear maximizer, the
quantity in (24) is at least twice the bilinear value.

Equations (23)--(25) are a genuine macroscopic principal
localization with an \(n^{3/2}\)-scale cross certificate.  The sharp
remaining wall is that (24) yields
\[
W(A)\ge
\max\{W(B)+W(D),\,\pi(S)\},
\tag{26}
\]
not their sum.  A convergence proof still needs either:

* a lower bound on the cross term near the **negative principal
  ground layer**, where \(f+g=2W(B)+2W(D)\); or
* a theorem converting the cyclic-monotonicity gaps in (5) and the
  projection payment (10) into a common profile carrying both the
  child-width budget and a fixed part of \(\pi(S)\).

This max-versus-sum obstruction is the precise normalization loss in
the present equality-case method.

## 6. The max-versus-sum wall is already sharp for an optimum

The missing additive inequality cannot hold for arbitrary optimal
signings, even at the first nontrivial conference order.  Consider
\[
A=
\begin{pmatrix}
0&-1&-1& 1& 1&-1\\
-1&0& 1&-1& 1&-1\\
-1&1&0&1&-1&-1\\
1&-1&1&0&-1&-1\\
1&1&-1&-1&0&-1\\
-1&-1&-1&-1&-1&0
\end{pmatrix}.
\tag{27}
\]
Direct multiplication gives \(A^2=5I\), and exhaustive evaluation
of its \(64\) Boolean states gives
\[
p(A)=\nu(A)=5,\qquad W(A)=5=M_6.
\tag{28}
\]
There are twelve positive ground states.  Switch by
\[
x_0=(1,-1,1,1,-1,-1)
\]
and use the zero cut \(S=\{1,2\}\) (indices numbered from \(1\)).
For this \(2+4\) principal split,
\[
\boxed{W(A[S])+W(A[S^c])=5=W(A).}
\tag{29}
\]
Nevertheless, under the uniform law on the twelve positive grounds,
the equality-profile cross payments are
\[
0,0,0,0,4,4,4,4,4,4,4,4
\]
up to ordering, and hence
\[
\boxed{\pi(S)=\frac83>0.}
\tag{30}
\]
Moreover, the cross block annihilates every pair of negative
principal grounds, so the negative-layer payment is exactly zero.

Thus neither a positive cubic/equality payment nor even an optimal
conference signing forces a term to be *added* to the exact child
width budget.  Any successful use of (23) must be asymptotic and
structural: it must show that repeated sharp max saturation creates
a quotient/recursion, or that the payment is large enough by itself
to reach the target normalization.  A universal inequality of the
form
\[
W(A)\ge W(A[S])+W(A[S^c])+c\,\pi(S)
\]
is false for every \(c>0\).

## 7. Exact max saturation forces two-sided Cartesian closure

Although the payment cannot be added, equality in the child-width
budget has a rigid second endpoint which was absent from the
one-sided principal-closure statement.

Put
\[
R_0=2W(B)+2W(D)
\]
and retain the positive deficits \(f,g\) from (1).  Define the
negative-end deficits
\[
\bar f(u)=p(B)+\nu(B)-f(u)=\nu(B)+H_B(u),
\qquad
\bar g(v)=p(D)+\nu(D)-g(v)=\nu(D)+H_D(v).
\tag{31}
\]

### Proposition 7.1 (two-sided capped bilinear law)

If the zero-cut split saturates centered-width superadditivity,
\[
W(A)=W(B)+W(D),
\tag{32}
\]
then, for every Boolean \(u,v\),
\[
\boxed{
|u^{\mathsf T}Cv|
\le
\min\{f(u)+g(v),\,\bar f(u)+\bar g(v)\}.
}
\tag{33}
\]
Consequently the cross block annihilates both endpoint Cartesian
families:
\[
\boxed{
P_{L_B^+}CP_{L_D^+}=0,
\qquad
P_{L_B^-}CP_{L_D^-}=0,
}
\tag{34}
\]
where \(L_i^\pm\) is the span of the positive/negative principal
ground states.  Every independent concatenation of two negative
principal grounds is a negative full ground, just as every
concatenation of positive grounds is a positive full ground.

#### Proof

The first term in the minimum is (2).  By (24) and (32),
\[
f(u)+g(v)+|u^{\mathsf T}Cv|\le R_0.
\]
Rearranging gives the second term in (33).  On a pair of negative
principal grounds, \(\bar f=\bar g=0\), so the cross term vanishes.
Their full energy is then
\(-\nu(B)-\nu(D)=-\nu(A)\), proving negative Cartesian closure and
(34). \(\square\)

There is a useful combined rank form.  Write
\[
r_i^\pm=\dim L_i^\pm,\qquad
r_i^\cup=\dim(L_i^++L_i^-).
\]
From (34),
\[
\boxed{
\begin{aligned}
\operatorname{rank}C
&\le
(n_2-r_2^\cup)+(n_1-r_1^+)+(n_1-r_1^-),\\
\operatorname{rank}C
&\le
(n_1-r_1^\cup)+(n_2-r_2^+)+(n_2-r_2^-).
\end{aligned}
}
\tag{35}
\]
Indeed, on \(L_2^++L_2^-\), the image is contained in the sum of
\((L_1^+)^\perp\) and \((L_1^-)^\perp\); the complementary domain
costs \(n_2-r_2^\cup\) further dimensions.  Transposition gives the
second inequality.

Thus repeated **sharp** max saturation is not structureless.  Along
a laminar zero-cut recursion, both endpoint ground faces factor as
Cartesian products at every saturated node.  After \(\ell\) leaves,
both endpoint faces contain independent type-constant
\(2^{\ell-1}\)-cubes, so the affine type-closure machinery applies
twice.  The remaining asymptotic issue is stability: an
\(o(n^{3/2})\) slack only gives \(o(n^{3/2})\) cross values on
negative ground pairs, and turning that scalar estimate into
approximate projection annihilation again requires quantitative
frame lower bounds.

## 8. Repeated sharp saturation: exact laminar classification

The exact induction can be completed.  Let \(\mathcal T\) be a
binary laminar tree of principal vertex sets.  At every internal node
\(V=V_1\sqcup V_2\), assume that \(V_1\) is a zero cut in a positive
ground gauge for \(A[V]\), and that
\[
W(A[V])=W(A[V_1])+W(A[V_2]).
\tag{36}
\]
Let the leaves be \(L_1,\ldots,L_\ell\).

### Theorem 8.1 (double endpoint torsor)

Under (36),
\[
\boxed{
W(A)=\sum_{a=1}^{\ell}W(A[L_a]),
\quad
p(A)=\sum_{a=1}^{\ell}p(A[L_a]),
\quad
\nu(A)=\sum_{a=1}^{\ell}\nu(A[L_a]).
}
\tag{37}
\]
Choose independently a positive reference ground
\(\alpha_a^+\) and a negative reference ground \(\alpha_a^-\) in
each leaf.  Then the full positive and negative ground faces contain
the two Cartesian torsors
\[
\boxed{
\left\{(\varepsilon_a\alpha_a^+)_{a\le\ell}:
\varepsilon\in\{\pm1\}^{\ell}\right\},
\qquad
\left\{(\varepsilon_a\alpha_a^-)_{a\le\ell}:
\varepsilon\in\{\pm1\}^{\ell}\right\}.
}
\tag{38}
\]
After switching by the positive references,
\[
\boxed{
\sum_{i\in L_a,\ j\in L_b}
a_{ij}\alpha_{a,i}^+\alpha_{b,j}^+=0
\qquad(a\ne b),
}
\tag{39}
\]
and after switching \(-A\) by the negative references,
\[
\boxed{
\sum_{i\in L_a,\ j\in L_b}
a_{ij}\alpha_{a,i}^-\alpha_{b,j}^-=0
\qquad(a\ne b).
}
\tag{40}
\]
Thus repeated max saturation is exactly a common leaf partition
carrying affine torsors on both opposite endpoint faces.  The
uniform laws on the torsors have correlation rank precisely
\(\ell\).

#### Proof

Proposition 7.1 factors both endpoint ground sets across the two
children and makes both endpoint values additive.  Induction down
the tree proves (37)--(38).  The type characters
\(\varepsilon_a\) in (38) are a full affine family.  Applying the
exact affine type-closure theorem to the positive family gives
(39); applying it to \(-A\) and the negative family gives (40).
The uniform torsor correlation matrix is block diagonal with one
rank-one block on every nonempty leaf, proving the rank assertion.
\(\square\)

At each node, if \(r_V^\pm\) denotes the span dimension of the
Cartesian endpoint family inherited from its leaves, then
\[
r_V^\pm=\sum_{L_a\subseteq V}r_{L_a}^\pm
\]
and the node cross block obeys
\[
\boxed{
\operatorname{rank}C_V\le |V|-r_V^+,
\qquad
\operatorname{rank}C_V\le |V|-r_V^-.
}
\tag{41}
\]
Conversely, a sign matrix of real rank \(k\) has at most \(2^k\)
distinct rows and at most \(2^k\) distinct columns: choose \(k\)
coordinates on which projection is injective on its row (or column)
space.  Hence a genuinely low-rank node cross block is an explicit
bounded-pattern sign quotient.

The exact repeated-saturation branch therefore has the requested
classification:

* \(\ell=\Omega(n)\) supplies a linear-rank endpoint cap law;
* \(\ell=o(n)\) supplies a common low-type affine torsor on both
  faces;
* low-rank cross blocks further collapse to bounded-pattern sign
  quotients.

What does **not** follow is bounded depth or a normalized blow-up.
For \(\ell\) leaves, (37) and the universal lower bound only give
\[
W(A)\ge c_*\sum_{a=1}^{\ell}|L_a|^{3/2},
\tag{42}
\]
which can be as small as
\(c_*n^{3/2}/\sqrt{\ell}\) for equal leaves.  Cross-block
discrepancy restores an order-\(n^{3/2}\) lower bound but currently
only with the known rectangular constant, not an additive term.
The square-root leaf count is therefore a genuine remaining
mesoscopic quotient, rather than something excluded by the exact
induction.

For near saturation, the same proof yields only
\[
|u_-^{\mathsf T}C_Vv_-|
\le 2\bigl[
W(A[V])-W(A[V_1])-W(A[V_2])
\bigr]
\tag{43}
\]
on negative child-ground pairs.  Summing (43) over a tree controls
the scalar cross values by the telescoping total slack, but no
uniform frame lower bound converts this into operator or rank
control.  This is the precise stability wall.
