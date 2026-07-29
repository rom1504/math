# Multiscale deterministic conversion of heavy local fields

## Goal and verdict

Let \(A\) be a symmetric zero-diagonal sign matrix of order \(n\),
\[
Q(A)=\max_{z\in\{\pm1\}^n}|z^\top A z|,
\]
and fix a Boolean witness \(x\).  Put
\[
u=Ax,\qquad h_i=\frac{|u_i|}{\sqrt n}.
\]

The aim was to replace a fixed cap on the fields by a dyadic,
tail-free conversion of \(\sum_i h_i^2\) into a lower bound on \(Q(A)\).

The campaign produced exact block inequalities which give the desired
cross-versus-internal recursion.  They are useful and seem not to have
been recorded in the existing notes in this form.  The same calculation
also locates the unresolved step exactly:

* a sparse heavy level whose field is genuinely external gives an
  additive, uncapped gain;
* a level whose field is generated internally must be passed to the
  induced principal block;
* replacing the old restricted witness by a ground state of that block
  loses precisely the replenishment/adaptivity gap from the earlier
  peeling route.

An exact universal-vertex construction proves that no increasing
functional of the squared fields of one witness can be a tail-free
leading-order lower bound.  Thus the deterministic multiscale route
does not by itself close the defect argument.  Its surviving target is
a **grouped-scale replenishment theorem**, not another scalar moment
inequality.

All statements below use the doubled normalization \(Q(A)\).

## 1. Global polarization and the first-moment bound

For a symmetric zero-diagonal matrix,
\[
\max_{y,z\in\{\pm1\}^n}|y^\top A z|\le 2Q(A). \tag{1.1}
\]

Indeed, put \(v=(y+z)/2\), \(w=(y-z)/2\).  Then
\[
y^\top A z=v^\top Av-w^\top Aw.
\]
Both \(v,w\) lie in \([-1,1]^n\), and coordinatewise randomized
completion (the diagonal of \(A\) is zero) gives
\[
|v^\top Av|,|w^\top Aw|\le Q(A).
\]

Taking \(y=\operatorname{sign}(Ax)\) yields
\[
\boxed{\quad
Q(A)\ge\frac12\|Ax\|_1
=\frac{\sqrt n}{2}\sum_i h_i.
\quad} \tag{1.2}
\]

This is sharp enough to show that, on a sequence with
\(Q(A)=O(n^{3/2})\), every one-witness level above \(H\sqrt n\) has at
most \(O(n/H)\) vertices.  It cannot control its squared-field mass:
one field of order \(n\) contributes order \(n^2\) to that mass but
only order \(n\) to (1.2).

## 2. Exact two-block inequalities

Partition the coordinates into \(S,T\), \(|S|=s\), and write
\[
A=
\begin{pmatrix}
D&B\\
B^\top&E
\end{pmatrix},
\qquad x=(x_S,x_T).
\]
Set
\[
e_T=x_T^\top E x_T,\qquad
C_T=\|Bx_T\|_1,\qquad
L_T=\|Bx_T\|_2,
\qquad Q_S=Q(D).
\]

### 2.1 Additive \(L_2\) extraction

For a uniform \(z\in\{\pm1\}^S\), compare the two full vectors
\((z,x_T)\) and \((-z,x_T)\).  Their energies are
\[
z^\top Dz+e_T\pm2z^\top Bx_T.
\]
For real \(a,b\),
\[
\max\{|a+b|,|a-b|\}=|a|+|b|.
\]
Average in \(z\), use
\(\mathbb E z^\top Dz=0\), Jensen, and the sharp lower Khintchine
constant \(1/\sqrt2\).  This gives
\[
\boxed{\quad
Q(A)\ge |e_T|+2\mathbb E_z|z^\top Bx_T|
\ge |e_T|+\sqrt2\,L_T.
\quad} \tag{2.1}
\]

This is an exact additive inequality; no ground-state assumption is
made on \(x_T\).

### 2.2 Biased \(L_1\) extraction

Let \(z_0=\operatorname{sign}(Bx_T)\), with arbitrary signs on zero
coordinates.  Orient the whole matrix by
\(\operatorname{sign}(e_T)\), and independently round a vector \(Z\)
with
\[
\mathbb E Z_i=t(z_0)_i,\qquad 0\le t\le1.
\]
The expected oriented energy is at least
\[
|e_T|+2tC_T-t^2Q_S.
\]
Consequently
\[
\boxed{\quad
Q(A)\ge |e_T|+\Phi(C_T,Q_S),
\quad} \tag{2.2}
\]
where
\[
\Phi(C,d)=\max_{0\le t\le1}(2tC-t^2d)
=
\begin{cases}
2C,&d=0,\\[2mm]
C^2/d,&0<C\le d,\\[2mm]
2C-d,&C\ge d.
\end{cases} \tag{2.3}
\]

The \(L_2\) estimate is strongest for diffuse cross fields.  The biased
\(L_1\) estimate is strongest for sparse coherent fields.

## 3. Conversion from the full fields on a selected set

On \(S\),
\[
u_S=Dx_S+Bx_T.
\]
Since every row of \(D\) has \(s-1\) signs,
\[
\|Dx_S\|_1\le s(s-1),\qquad
\|Dx_S\|_2\le\sqrt{s}(s-1).
\]
Therefore (2.1)--(2.2) imply the completely explicit bounds
\[
\boxed{\quad
Q(A)\ge |e_T|
+\sqrt2\left(
\|u_S\|_2-\sqrt{s}(s-1)
\right)_+,
\quad} \tag{3.1}
\]
and
\[
\boxed{\quad
Q(A)\ge |e_T|
+\Phi\left(
\bigl(\|u_S\|_1-s(s-1)\bigr)_+,\,
Q(A[S])
\right).
\quad} \tag{3.2}
\]

Replacing \(Q(A[S])\) by the elementary upper bound \(s(s-1)\) gives
a version depending only on the field level and its cardinality.

There is also a cross-only inequality.  Put
\[
y_S=\operatorname{sign}(u_S).
\]
Then
\[
\|u_S\|_1
=y_S^\top Dx_S+y_S^\top Bx_T.
\]
The first term has magnitude at most \(s(s-1)\).  Flipping the relative
global sign between \(S\) and \(T\) shows that
\[
Q(A)\ge2|y_S^\top Bx_T|.
\]
Hence
\[
\boxed{\quad
Q(A)\ge
2\left(\|u_S\|_1-s(s-1)\right)_+.
\quad} \tag{3.3}
\]

## 4. Dyadic level-set recursion

Fix \(H>0\), let \(L_\ell=2^\ell H\), and define
\[
S_\ell
=\left\{i:
L_\ell\sqrt n<|u_i|
\le2L_\ell\sqrt n
\right\},
\qquad s_\ell=|S_\ell|.
\]
For \(T_\ell=S_\ell^c\), (3.1)--(3.3) give
\[
Q(A)\ge
|x_{T_\ell}^\top A[T_\ell]x_{T_\ell}|
+\sqrt2\sqrt{s_\ell}
\left(L_\ell\sqrt n-(s_\ell-1)\right)_+, \tag{4.1}
\]
\[
Q(A)\ge
2s_\ell
\left(L_\ell\sqrt n-(s_\ell-1)\right)_+, \tag{4.2}
\]
and the stronger biased form
\[
Q(A)\ge
|x_{T_\ell}^\top A[T_\ell]x_{T_\ell}|
+\Phi\left(
s_\ell\left(L_\ell\sqrt n-(s_\ell-1)\right)_+,
Q(A[S_\ell])
\right). \tag{4.3}
\]

In particular, if
\[
s_\ell-1\le\frac13L_\ell\sqrt n, \tag{4.4}
\]
then the argument of \(\Phi\) is at least twice the elementary bound
on \(Q(A[S_\ell])\).  Equations (2.3)--(4.3) yield the clean sparse
branch
\[
\boxed{\quad
Q(A)\ge
|x_{T_\ell}^\top A[T_\ell]x_{T_\ell}|
+s_\ell L_\ell\sqrt n.
\quad} \tag{4.5}
\]

Keeping the constants before the simplification gives
\[
Q(A)\ge
|x_{T_\ell}^\top A[T_\ell]x_{T_\ell}|
+2s_\ell L_\ell\sqrt n-3s_\ell(s_\ell-1).
\tag{4.6}
\]

Thus a sparse extreme level is not lost by capping: its first-moment
field contribution can be added to the energy of the old restricted
state.  If (4.4) fails, then \(s_\ell\gtrsim L_\ell\sqrt n\); the level
is large enough that its fields may be generated internally by a
dense principal block, and the recursion must pass to \(A[S_\ell]\).

This is the exact cross-versus-internal dichotomy requested by the
multiscale program.

## 5. Positive fields at a global maximizer

Suppose now that \(x\) is switched to \(\mathbf1\), the energy is
positive, and
\[
q=\mathbf1^\top A\mathbf1=Q(A).
\]
Write \(r=A\mathbf1\).  One-flip optimality gives \(r_i\ge0\), and
global maximality gives, for every \(S\),
\[
0\le b_S:=\mathbf1_S^\top A_{S,T}\mathbf1_T\le q/2.
\tag{5.1}
\]
Put
\[
R_S=\sum_{i\in S}r_i,\qquad
a_S=\mathbf1_S^\top A[S]\mathbf1_S,\qquad
e_T=\mathbf1_T^\top A[T]\mathbf1_T.
\]
Then exactly
\[
R_S=a_S+b_S,\qquad
e_T=q-2R_S+a_S. \tag{5.2}
\]
In particular,
\[
\boxed{\quad
Q(A[S])\ge(R_S-q/2)_+.
\quad} \tag{5.3}
\]

The continuous cut test is also exact.  Select every vertex of \(S\)
independently with probability \(t\).  The expected signed cut is
\[
tR_S-t^2a_S.
\]
Since every cut has weight at most \(q/2\),
\[
\boxed{\quad
tR_S-t^2a_S\le q/2
\qquad(0\le t\le1).
\quad} \tag{5.4}
\]
For \(R_S\le q\), optimizing (5.4) recovers
\(a_S\ge R_S-q/2\).

These inequalities show why scalar multiscale updates stop on the
positive branch.  For a nonnegative field profile with
\(\frac1n\sum_i r_i/\sqrt n=q/n^{3/2}=c\), every capped local
perturbation sees at most this same first moment.  A two-point profile
\[
h=
\begin{cases}
L,&\text{with probability }c/L,\\
0,&\text{otherwise}
\end{cases}
\]
has mean \(c\), arbitrarily large second moment \(cL\), and saturates
every cap:
\[
\frac1H\mathbb E[h\min(h,H)]\le c
\quad\text{for all }H>0. \tag{5.5}
\]
Additional graph structure, rather than a better choice of scalar cap,
is necessary.

## 6. Exact no-go: one universal positive vertex

Let \(D\) have order \(m\), and switch a positive maximizer to
\(\mathbf1\):
\[
R=\mathbf1^\top D\mathbf1=Q(D),\qquad r=D\mathbf1.
\]
Adjoin \(k\) universally positive vertices:
\[
\widetilde D_k=
\begin{pmatrix}
J_k-I_k&J_{k,m}\\
J_{m,k}&D
\end{pmatrix}.
\]
Then exactly
\[
Q(\widetilde D_k)=R+2km+k(k-1). \tag{6.1}
\]
At the all-one maximizer, the new local fields are
\[
\widetilde r_i=m+k-1\quad(i\le k),\qquad
\widetilde r_{k+j}=r_j+k,
\]
so
\[
\boxed{\quad
\sum_i\widetilde r_i^2
=\sum_jr_j^2+2kR+mk^2+k(m+k-1)^2.
\quad} \tag{6.2}
\]

Already for \(k=1\), (6.2) adds \(m^2+O(m^{3/2})\) to the squared
fields, while (6.1) adds only \(2m\) to \(Q\).  Hence
\[
\frac{Q(\widetilde D_1)}{(m+1)^{3/2}}
-\frac{Q(D)}{m^{3/2}}\longrightarrow0,
\]
whereas the normalized one-witness squared-field mass changes by
\(1+o(1)\).

Therefore no universally calibrated stability inequality can force a
strict leading improvement merely because
\(\|Ax\|_2^2/n^2\) increases by a fixed amount.  More precisely, if a
putative bound
\[
\frac{Q(A)}{n^{3/2}}
\ge G\left(\frac{\|Ax\|_2^2}{n^2}\right)-o(1)
\]
is asymptotically tight on the base sequence \(D_m\), then it cannot
also have a fixed positive increase
\(G(s+1)-G(s)>0\) at the corresponding squared-field values.  This
remains true when \(x\) is a positive global maximizer and all its
oriented fields are nonnegative.

The sparse formula (4.6) behaves correctly on this example: it charges
only order \(m\), exactly the lower-order cost of the universal
vertex.  Thus (4.6) cannot be iterated into a leading gain without
using what happens to the successive principal cores.

## 7. Fixed-scale internal obstructions

The positive-clique tower gives the corresponding obstruction at many
finite scales.  Partition \(n=N^2\) vertices into blocks of size
\[
k_j=K_j\sqrt n
\]
occupying vertex fractions \(p_j\), put \(+1\) on every block interior,
and choose the between-block edges randomly.  There are realizations
with
\[
Q(A)=O(n^{3/2})
\]
provided
\[
\sum_jp_jK_j<\infty.
\]
At the all-one state, the internal contribution to the fields on a
type-\(j\) block is \(K_j\sqrt n+O(1)\).  The random exterior part is
only of the natural \(\sqrt n\) scale for a typical vertex.  Thus the
construction supports arbitrarily many dyadic field levels while its
leading deterministic cost is precisely
\[
\left(\sum_jp_jK_j\right)n^{3/2}.
\]

This stress test shows that any multiscale theorem based only on level
cardinalities and first-moment field charges can have an integrable but
non-uniformly bounded tower.  To distinguish a near-minimizer, the
tower must be charged to *excess above the asymptotic optimum*, or its
high-scale vertices must be peeled.

## 8. Why the recursion meets the replenishment gap

Formula (4.5) adds the heavy-level gain to
\[
|x_T^\top A[T]x_T|,
\]
the energy of the **old restricted state**.  A recursive theorem wants
instead the new core norm \(Q(A[T])\).  The difference is the
replenishment gap.

For a positive ground state and a peeled block \(S\), write
\[
d=Q(A)-Q(A[T]),\qquad
g=Q(A[T])-x_T^\top A[T]x_T.
\]
Together with (5.2), one obtains the exact block identity
\[
\boxed{\quad
2R_S=d+a_S+g.
\quad} \tag{8.1}
\]

Across disjoint peeling steps, the decrements telescope and the
induced-block norms have the known bound
\[
\sum_t d_t\le Q(A_0),\qquad
\sum_tQ(A[S_t])\le2Q(A_0).
\]
Thus every dyadic heavy-field recursion reduces to controlling
\(\sum_tg_t\).  This is not an artifact of the estimates above:
finite regular examples have a large \(g_t\) even when the deleted
cross row is orthogonal to every successor ground state.  Hence no
pointwise replacement of the old state by a successor ground state is
valid, even under \(O(\sqrt n)\) spectral regularity.

The exact surviving conjecture is consequently:

> For a peeling order formed by grouping current positive local fields
> into dyadic levels, is the cumulative replenishment
> \(\sum_tg_t\) bounded by \(O(Q(A_0))\), or by the normalized
> optimality excess plus \(o(n^{3/2})\)?

An affirmative excess-sensitive version would let (4.5)--(4.6) peel
the genuinely sparse tail, pass internally generated mass to principal
blocks, and preserve the leading normalized lower bound.

## 9. Bottom line

The multiscale attack does give a tail-free deterministic conversion,
but only in recursive form:

\[
\text{heavy level}
\quad\Longrightarrow\quad
\begin{cases}
\text{additive cross gain via (4.5)--(4.6),}\\
\text{or induced-block recursion.}
\end{cases}
\]

It cannot be collapsed to a scalar function of
\(\sum_i|(Ax)_i|^2\).  Universal positive vertices falsify such a
collapse at the exact leading scale.  The only unresolved loss in the
recursive form is the same replenishment/adaptivity gap already
isolated by principal peeling.  Thus further work should target a
**grouped dyadic replenishment theorem**, not another cap or
single-state moment estimate.
