# Finite-fibre renormalization: exact audit

## Status

Let
\[
Q(B):=\max_{x\in\{\pm1\}^{N}}|x^\top Bx|
\]
for a symmetric zero-diagonal sign matrix \(B\) of order \(N\).  This
note audits the fixed three-fibre lift
\[
\mathcal T(C):=C\otimes R+I_k\otimes D,\qquad
R=J_3-2I_3,\qquad D=J_3-I_3. \tag{1}
\]
The main conclusions are negative but exact:

1. the lift has an exact signed four-state variational formula;
2. for the order-six Paley conference matrix,
   \[
   Q(\mathcal T(C_6))=78;
   \]
   its original half-energy normalization is
   \(78/(2\cdot18^{3/2})=0.5106882308\ldots>1/2\);
3. fibre pairs compose associatively, but \(Q\) is not multiplicative
   under that composition;
4. \(Q(C)\) is not enough state information for even one lift: two
   order-three signings with the same \(Q(C)=6\) have lifted values
   \(24\) and \(36\).

Thus this particular gadget cannot yield scale interpolation, and no
fixed-fibre scalar recursion in \(Q\) can describe the general lift.

## 1. Exact four-state formula

Every fibre spin belongs to one of four antipodal pairs.  Choose
\[
u_0=(1,1,1),\quad
u_1=(-1,1,1),\quad
u_2=(1,-1,1),\quad
u_3=(1,1,-1).
\]
Write the spin in fibre \(i\) uniquely as
\[
x_i=z_i u_{t_i},\qquad z_i\in\{\pm1\},\quad
t_i\in\{0,1,2,3\}.
\]
Direct multiplication gives
\[
K=(u_a^\top Ru_b)_{a,b=0}^3
=
\begin{pmatrix}
3&1&1&1\\
1&-5&3&3\\
1&3&-5&3\\
1&3&3&-5
\end{pmatrix},
\qquad
\ell=(u_a^\top Du_a)_{a=0}^3=(6,-2,-2,-2).
\tag{2}
\]
Consequently,
\[
\boxed{
Q(\mathcal T(C))
=
\max_{\substack{t_i\in\{0,1,2,3\}\\z_i\in\{\pm1\}}}
\left|
\sum_{i,j=1}^k c_{ij}z_i z_jK_{t_i,t_j}
+\sum_{i=1}^k\ell_{t_i}
\right|.
}
\tag{3}
\]
This is exact; no asymptotic or probabilistic approximation enters.

There is also a useful three-replica form.  If the three columns of the
fibre array are \(a,b,c\in\{\pm1\}^k\), then
\[
\begin{aligned}
x^\top\mathcal T(C)x
={}&-a^\top Ca-b^\top Cb-c^\top Cc\\
&+2a^\top Cb+2a^\top Cc+2b^\top Cc\\
&+2a\cdot b+2a\cdot c+2b\cdot c .
\end{aligned}
\tag{4}
\]
Equation (4) already shows why the lift needs a joint multi-replica
profile rather than the single scalar \(Q(C)\).

## 2. Conference spectral decomposition

Suppose \(C^2=(k-1)I_k\).  The uniform fibre direction is an
eigenvector of \(R,D\) with eigenvalues \(1,2\), respectively.  Each
of the two transverse directions has eigenvalues \(-2,-1\).
Therefore (1) is orthogonally equivalent to
\[
C+2I_k
\quad\oplus\quad
(-2C-I_k)
\quad\oplus\quad
(-2C-I_k). \tag{5}
\]
In particular,
\[
\|\mathcal T(C)\|_{\mathrm{op}}=2\sqrt{k-1}+1
\]
and hence
\[
Q(\mathcal T(C))
\le 3k(2\sqrt{k-1}+1). \tag{6}
\]
Thus
\[
\limsup_{k\to\infty}
\frac{Q(\mathcal T(C))}{(3k)^{3/2}}
\le \frac2{\sqrt3}. \tag{7}
\]
For the original half-energy objective the corresponding ceiling is
\(1/\sqrt3\), which is already worse than \(1/2\).  The spectral
bound does not rescue this gadget.

## 3. Exact order-six value

Take the Paley conference matrix
\[
C_6=
\begin{pmatrix}
0&1&1&1&1&1\\
1&0&1&-1&-1&1\\
1&1&0&1&-1&-1\\
1&-1&1&0&1&-1\\
1&-1&-1&1&0&1\\
1&1&-1&-1&1&0
\end{pmatrix},
\qquad C_6^2=5I_6. \tag{8}
\]
Exhaustive enumeration gives \(Q(C_6)=10\) and
\[
\boxed{Q(\mathcal T(C_6))=78.} \tag{9}
\]
The following gives a short certificate for (9), independent of the
full \(8^6\)-state enumeration.

Let \(m\) be the number of uniform fibres.  The squared norm of the
uniform fibre projection and its transverse complement are
\[
a_m=2+\frac83m,\qquad
b_m=16-\frac83m. \tag{10}
\]
From (5),
\[
-x^\top\mathcal T(C_6)x
\le(\sqrt5-2)a_m+(2\sqrt5+1)b_m. \tag{11}
\]
For every \(m\ge1\), the right side is less than \(78\).  Therefore a
state with energy at most \(-78\) must have all six fibres
nonuniform.

For nonuniform fibres \(t_i\in\{1,2,3\}\), put
\[
d_{ij}=c_{ij}z_i z_j,\qquad
S=\sum_{i<j}d_{ij},\qquad
T=\sum_{\substack{i<j\\t_i=t_j}}d_{ij}. \tag{12}
\]
Since \(K_{aa}=-5\), \(K_{ab}=3\) for \(a\ne b\), and every internal
fibre contributes \(-2\), formula (3) becomes
\[
x^\top\mathcal T(C_6)x=6S-16T-12. \tag{13}
\]
Gauge-fix \(z_1=1\).  The 32 switchings and \(3^6\) colourings give
the following exact finite table:

| \(S\) | maximum \(T\) |
|---:|---:|
| \(-5\) | \(2\) |
| \(-3\) | \(3\) |
| \(3\) | \(4\) |
| \(5\) | \(5\) |

Thus the minimum in (13) is \(-78\), attained at \(S=-3,T=3\).
For the positive side, (5) and (10) give
\[
x^\top\mathcal T(C_6)x
\le(2+\sqrt5)a_m+(2\sqrt5-1)b_m
\le18(2+\sqrt5)<78. \tag{14}
\]
This proves (9).  One lexicographic fibre witness, with
\(\{\pm1\}^3\) ordered as in `itertools.product((-1,1),repeat=3)`,
is
\[
(1,1,3,2,4,5),
\]
and has energy \(-78\).

The two normalizations are
\[
\frac{78}{18^{3/2}}=1.0213764617\ldots,\qquad
\frac{78}{2\cdot18^{3/2}}=0.5106882308\ldots . \tag{15}
\]
The second is already above the conference upper constant \(1/2\).

## 4. The algebraic semigroup and its failure as a scalar recursion

For any fixed full sign matrix \(R\) and zero-diagonal sign matrix
\(D\), define
\[
\mathcal T_{R,D}(C)=C\otimes R+I\otimes D.
\]
Composition is exact:
\[
\mathcal T_{R_2,D_2}(\mathcal T_{R_1,D_1}(C))
=\mathcal T_{R_1\otimes R_2,\,
D_1\otimes R_2+I\otimes D_2}(C). \tag{16}
\]
Thus fibre pairs form an associative semigroup under
\[
(R_1,D_1)\star(R_2,D_2)
=
(R_1\otimes R_2,\,
D_1\otimes R_2+I\otimes D_2). \tag{17}
\]
The second component in (17) remains a valid zero-diagonal sign
filling: off-diagonal entries receive exactly one nonzero summand.

After \(r\) repetitions of the same \(s\)-fibre gadget,
\[
R_r=R^{\otimes r},\qquad
D_r=\sum_{j=0}^{r-1}
I_{s^j}\otimes D\otimes R^{\otimes(r-1-j)}. \tag{18}
\]
However, the Boolean norm does not respect this semigroup.  For the
present \(R\),
\[
Q(R)=5,\qquad Q(R\otimes R)=33>25=Q(R)^2. \tag{19}
\]
Its normalized value actually jumps:
\[
\frac{Q(R)}{3^{3/2}}=0.9622504\ldots,\qquad
\frac{Q(R^{\otimes2})}{9^{3/2}}=\frac{33}{27}
=1.2222222\ldots . \tag{20}
\]
So iteration amplifies transverse Boolean choices instead of
contracting them.

There is an even sharper obstruction to a scalar recursion.  Let
\(C^-\) and \(C^\triangle\) be the order-three signings whose upper
triangular edge lists are
\[
(-1,-1,-1),\qquad(-1,-1,1),
\]
in the order \((12),(13),(23)\).  Exhaustive enumeration gives
\[
Q(C^-)=Q(C^\triangle)=6,
\]
but
\[
Q(\mathcal T(C^-))=24,\qquad
Q(\mathcal T(C^\triangle))=36. \tag{21}
\]
Therefore \(Q(\mathcal T(C))\) is not a function of \(Q(C)\), even at
one fixed finite scale.

One can see the missing statistic directly.  In (4), set the three
replicas to \((a,b,a)\).  Then
\[
x^\top\mathcal T(C)x
=4a^\top Cb-b^\top Cb+4a\cdot b+2k. \tag{22}
\]
Choosing \(a=-\operatorname{sgn}(b^\top Cb)\operatorname{sgn}(Cb)\)
shows
\[
Q(\mathcal T(C))
\ge
4\|Cb\|_1+|b^\top Cb|-6k. \tag{23}
\]
The lift therefore sees a joint local-field/quadratic profile, not
only the diagonal Boolean norm \(Q(C)\).

## 5. Consequence for convergence

The pair semigroup (17) is useful algebraically, but it does not by
itself interpolate the normalized minima:

- the explicit three-fibre gadget misses the desired \(1/2\) upper
  constant already at order \(18\);
- the transverse spectral multiplier is \(2/\sqrt3>1\);
- the kernel norm is nonmultiplicative at depth two;
- and equal scalar input norms can produce different lifted norms.

A viable renormalization theory would therefore need a profile-valued
state recording all finite-replica overlaps and mixed energies, plus
a separate compactness or contraction theorem for the minimizing
profiles.  Ordinary scalar submultiplicativity is decisively ruled
out by (19) and (21).

## Reproducibility

`verify_finite_fibre.py` checks (2), (8)--(15), (19), and (21) in
exact integer arithmetic.
