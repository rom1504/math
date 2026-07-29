# Correlated internal fill: kernel rounding and Gaussian width

## Status

Independent internal resampling is ruled out by
`internal_resampling_traffic_no_go.md`.  This note formulates the
remaining conditioned problem in the exact cap kernel, gives a
Gaussian-width criterion for any sufficiently rich kernel rounding
law, derives the quantitative projection/rank threshold that such a
law must meet, and records an exact small-order audit.

No unconditional centered fill is proved.  The main new conclusions
are:

1. exact \(+W\) and \(-W\) cross caps impose the same linear
   constraints;
2. a conditional subgaussian law on the kernel sign points converts
   the raw union count into a Gaussian-width criterion;
3. even under ideal subgaussian rounding, a competitive
   \(W\le\frac12n^{3/2}\) needs a substantial projection gain:
   at balanced endpoints, the projected typical norm must fall by a
   factor below \(1/\sqrt{4\log2}=0.60056\ldots\);
4. in the most optimistic rank-only model, the exact-cap span must
   capture more than \(27.865\%\) of the internal edge dimension at
   the spectral ceiling, and more at smaller \(W\);
5. exact finite examples show that the complete thick-band feasible
   sign set can collapse to the original pair \(\{\pm h\}\).

## 1. Exact cap kernel

Use the top/bottom decomposition
\[
A=c+h
\]
with cross matrix \(C\), internal edge set
\[
E_{\rm int}=E(U)\cup E(V),\qquad m=|E_{\rm int}|,
\]
and
\[
\|C\|_{\infty\to1}=W.
\]
For an orientation \(\alpha=(p,q)\), write
\[
t_\alpha=p^\mathsf TCq,\qquad
v_\alpha=\delta(p,q)|_{E_{\rm int}},
\qquad
b_\alpha=\frac{W-|t_\alpha|}{2}.
\tag{1}
\]
The exact internal-fill problem is
\[
g\in\{\pm1\}^m,\qquad
|g\cdot v_\alpha|\le b_\alpha
\quad\text{for every }\alpha.
\tag{2}
\]

Let
\[
\mathcal Z_0=\{\alpha:|t_\alpha|=W\},\qquad
L_0=\operatorname{span}\{v_\alpha:\alpha\in\mathcal Z_0\},
\qquad
K=L_0^\perp.
\tag{3}
\]
Every feasible fill belongs to
\[
\mathcal S_0=K\cap\{\pm1\}^m.
\tag{4}
\]
The original \(h\) and \(-h\) lie in \(\mathcal S_0\).

The doubled cap supplies only one system of equations.  The map
\[
(p,q)\longmapsto(p,-q)
\]
changes \(t_\alpha\) to \(-t_\alpha\) but leaves \(v_\alpha\)
unchanged.  Hence the \(+W\) and \(-W\) constraints in (3) coincide.

The stronger dual-cone interval is
\[
\boxed{
|d-H_B(p)-H_D(q)|\le W-|p^\mathsf TCq|.
}
\tag{5}
\]
At an exact cap this gives
\[
h\cdot\phi(p,q)=d.
\tag{6}
\]
Thus \(h\) is an unavoidable null direction of the centered internal
feature Gram.  This null direction alone does not lower the width:
the two sign points \(\pm h\) simply interchange the two oriented
endpoint signings.

## 2. Why ordinary partial coloring is not automatic

The convex body
\[
\mathcal P=
\left\{
x\in K\cap[-1,1]^m:
|x\cdot v_\alpha|\le b_\alpha\ \forall\alpha
\right\}
\tag{7}
\]
contains \(0,h,-h\).  The fractional point \(0\) has perfect midpoint,
but the objective requires a cube vertex.  A partial-coloring walk in
\(K\) may freeze many coordinates and still stop with an
unroundable residual system.  In the extreme case,
\[
\mathcal P\cap\{\pm1\}^m=\{\pm h\}.
\tag{8}
\]
Therefore dimension of \(K\), by itself, is not a rounding theorem.

Equivalently, writing \(g=h\odot\varepsilon\), exact-cap preservation
becomes
\[
\sum_e h_ev_{\alpha,e}\varepsilon_e=0
\qquad(\alpha\in\mathcal Z_0).
\tag{9}
\]
Although \(\varepsilon=\pm\mathbf1\) solve (9), a large real nullspace
need not contain any other Boolean point.

## 3. Conditional Gaussian-width criterion

The useful invariant is not the number of cap constraints but the
quality of the best probability law supported on \(\mathcal S_0\).
Suppose there is a symmetric probability measure \(\nu\) on
\(\mathcal S_0\) such that the process
\[
X_z(g)=z\cdot g,\qquad g\sim\nu,
\]
has \(L\)-subgaussian increments on \(K\):
\[
\mathbb E_\nu
\exp\{\lambda(X_z-X_{z'})\}
\le
\exp\left\{
\frac{L^2\lambda^2}{2}\|z-z'\|_2^2
\right\}
\quad(z,z'\in K).
\tag{10}
\]

For every non-cap orientation, project its feature into \(K\) and
normalize by its available margin:
\[
\mathcal V=
\left\{
\pm\frac{P_Kv_\alpha}{b_\alpha}:
b_\alpha>0
\right\}.
\tag{11}
\]
To demand a midpoint bound \(|g\cdot\mathbf1|\le\tau\), add
\[
\pm\frac{P_K\mathbf1}{\tau}
\]
to \(\mathcal V\).

Generic chaining for the subgaussian process gives
\[
\mathbb E_\nu\sup_{z\in\mathcal V}X_z
\le C L\,\gamma_2(\mathcal V,\|\cdot\|_2),
\tag{12}
\]
where \(C\) is universal.  Consequently the explicit sufficient
condition
\[
\boxed{
C L\,\gamma_2(\mathcal V,\|\cdot\|_2)<1
}
\tag{13}
\]
produces a kernel sign point obeying all bands (2) and the requested
midpoint bound.  Up to universal constants, \(\gamma_2\) is the
Gaussian width
\[
\mathbb E\sup_{z\in\mathcal V}\langle Z,z\rangle,
\qquad Z\sim N(0,P_K).
\tag{14}
\]

This is the promised replacement of raw cardinality by projected
Gaussian width.  Its unresolved hypothesis is the existence of a
kernel-supported law with \(L=O(1)\).  The degenerate law on
\(\{\pm h\}\) has \(L\) of order \(\sqrt m\) in its active direction
and gives no improvement.

## 4. Explicit shell threshold

For a shell \(\mathcal A\) of \(e^{(\mathfrak h+o(1))n}\)
orientations, suppose
\[
b_\alpha=(\beta+o(1))n^{3/2},\qquad
\|P_Kv_\alpha\|_2^2
\le(\theta+o(1))n^2
\quad(\alpha\in\mathcal A).
\tag{15}
\]
Even the raw subgaussian estimate already gives the exact benchmark
\[
\boxed{
\mathfrak h<
\frac{\beta^2}{2L^2\theta}.
}
\tag{16}
\]
The Gaussian-width condition can improve on (16) when the shell has
strong metric correlation, but cannot be worse than the corresponding
entropy estimate.

For typical orientations of a balanced top/bottom partition,
\[
k(p,q)=\left(\frac18+o(1)\right)n^2,
\qquad
b_\alpha=\left(\frac c2+o(1)\right)n^{3/2},
\quad W=cn^{3/2}.
\tag{17}
\]
If projection onto \(K\) reduces the typical norm by a factor
\(\eta\), then
\[
\theta=\frac{\eta^2}{8},\qquad
\beta=\frac c2,
\]
and (16) becomes
\[
\boxed{
\mathfrak h<\frac{c^2}{L^2\eta^2}.
}
\tag{18}
\]
At the full orientation entropy \(\mathfrak h=\log2\) and the largest
competitive constant \(c=1/2\), an ideal \(L=1\) law still requires
\[
\boxed{
\eta^2<\frac1{4\log2},
\qquad
\eta<0.600561\ldots.
}
\tag{19}
\]
Thus conditioning must remove more than \(63.93\%\) of the typical
projected squared norm.  An arbitrary small gain over independent
rounding is not enough.

## 5. Optimistic rank threshold

Under uniform orientations, the internal cut-feature second moment is
exactly
\[
G_{\rm unif}
=\mathbb E[v_\alpha v_\alpha^\mathsf T]
=\frac14(I_m+J_m).
\tag{20}
\]
For a rank-\(r_0\) cap span \(L_0\),
\[
\mathbb E\|P_Kv_\alpha\|_2^2
=\operatorname{tr}(P_KG_{\rm unif}).
\tag{21}
\]
The most optimistic rank-\(r_0\) projection captures the all-one
eigenvector first, leaving
\[
\operatorname{tr}(P_KG_{\rm unif})
\ge\frac{m-r_0}{4}.
\tag{22}
\]
Since the unprojected typical squared norm is \(m/2\), this gives
\[
\eta^2\ge\frac{1-r_0/m}{2}.
\tag{23}
\]

Combining (18), (23), and \(\mathfrak h=\log2\), even this optimistic
rank-only model needs
\[
\boxed{
\frac{r_0}{m}
>
1-\frac{2c^2}{L^2\log2}.
}
\tag{24}
\]
At \(c=1/2,L=1\),
\[
\boxed{
\frac{r_0}{m}>1-\frac1{2\log2}
=0.278652\ldots.
}
\tag{25}
\]
At \(c=0.336493\ldots\), the threshold rises to approximately
\[
\frac{r_0}{m}>0.6733.
\tag{26}
\]
If the cap span misses the all-one direction, the true requirement is
stronger.  These are necessary benchmarks for a rank-based proof, not
sufficient conditions for Boolean kernel rounding.

## 6. Exact small-order audit

The verifier `audit_correlated_internal_fill.py` exhaustively
enumerates all switching classes through \(n=7\), selects
centered-width minimizers, and audits one exact representative of
every endpoint profile.  It counts both the cap-kernel sign points
and the sign points satisfying every thick band (2).

| \(n\) | \(W\) | \(|d|\) | split | \(m\) | cap rank | cap signs | all-band signs | best all-band \(|g\cdot1|\) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 1 | \(1+2\) | 1 | 0 | 2 | 2 | 1 |
| 4 | 4 | 2 | \(2+2\) | 2 | 0 | 4 | 4 | 0 |
| 4 | 4 | 0 | \(2+2\) | 2 | 0 | 4 | 4 | 0 |
| 5 | 4 | 0 | \(2+3\) | 4 | 1 | 8 | 4 | 0 |
| 5 | 4 | 0 | \(1+4\) | 6 | 0 | 64 | 12 | 0 |
| 6 | 5 | 0 | \(3+3\) | 6 | 2 | 10 | 2 | 0 |
| 6 | 5 | 0 | \(1+5\) | 10 | 0 | 1024 | 12 | 0 |
| 7 | 8 | 1 | \(3+4\) | 9 | 0 | 512 | 64 | 1 |
| 7 | 8 | 1 | \(3+4\) | 9 | 1 | 256 | 40 | 1 |
| 7 | 8 | 1 | \(2+5\) | 11 | 1 | 768 | 120 | 1 |

All arithmetic is integral, and the script independently verifies that
the original \(h\) satisfies every band.

Two observations matter:

- the \(n=4,|d|=2\) representative has a same-width centered fill;
- the balanced \(n=6\) representative has exactly two all-band sign
  points.  Since \(\pm h\) are always feasible, the full thick-band
  set is exactly \(\{\pm h\}\).

Thus collapse of the conditioned feasible set is a real finite
phenomenon, not merely a logical possibility.  Any asymptotic theorem
must use additional global-minimality or large-order structure.

## 7. Surviving dichotomy

The correlated-fill route is now reduced to two quantitative objects:

1. **kernel richness:** construct a symmetric \(O(1)\)-subgaussian
   law on \(\mathcal S_0\), or prove that failure supplies a discrete
   block replacement lowering \(W\);
2. **projected thick-cap width:** prove that the normalized projected
   feature set (11) has Gaussian width below the threshold (13).

Small eigenvalues of the cap Gram are useful only if they create
additional Boolean kernel points.  Large real nullity without Boolean
richness is insufficient, and full-rank thick features with margins
of order \(\sqrt n\) are insufficient.  The exact obstruction is the
integrality gap between the fractional center \(0\) and the sign
vertices of the band polytope (7).
