# Sparse perturbations of conference signings

## Purpose

This note studies the following proposed construction.  Start from a symmetric
conference signing \(C=(c_{ij})\), independently flip every edge with
probability \(\delta\), and hope that:

1. the sparse Paley square-wave resonances are damped by the factor
   \(1-2\delta\); while
2. the exponentially many ordinary configurations pay only their
   energy-layer entropy, rather than the full \(2^n\) union-bound penalty.

The calculation below is exact at exponential scale.  It gives:

- a **conditional all-orders upper constant**
  \[
  0.498605960816\ldots<\frac12
  \]
  from the conjectured flat energy-layer envelope with a subexponential
  exceptional resonant set; and
- a **sharp obstruction** to reaching the ROM value
  \(\sqrt{15}/8\) by independent flips.  Near-threshold ordinary layers gain
  more from the random perturbation than the resonant cap loses.

The first calculation reduced the desired strict upper bound to one concrete
energy-counting statement.  A subsequent exact Hamming-sphere calculation
shows that this statement is **false** for every sequence with a cap-near
resonance: one resonant vector creates an exponentially large cloud of
near-resonant vectors.  The corrected surviving route must therefore
stratify simultaneously by energy and distance from the resonant centers,
and exploit the strong correlation of the perturbation noise inside each
cloud.

## 1. Exact one-configuration law

Write
\[
H_C(x)=\sum_{i<j}c_{ij}x_ix_j,\qquad
r_C(x)=\frac{x^\top Cx}{n\sqrt{n-1}}
      =\frac{2H_C(x)}{n\sqrt{n-1}}.
\]
Let \(\eta_{ij}\) be independent signs with
\[
\Pr(\eta_{ij}=-1)=\delta,\qquad
\mu:=\mathbb E\eta_{ij}=1-2\delta,
\]
and put \(a_{ij}=c_{ij}\eta_{ij}\).  For fixed \(x\), set
\(b_{ij}=c_{ij}x_ix_j\).  If \(N=\binom n2\), then
\[
N_+(x)=\frac{N+H_C(x)}2,\qquad
N_-(x)=\frac{N-H_C(x)}2.
\]
Consequently the complete moment generating function is
\[
\mathbb E e^{\lambda H_A(x)}
=\big((1-\delta)e^\lambda+\delta e^{-\lambda}\big)^{N_+(x)}
 \big((1-\delta)e^{-\lambda}+\delta e^\lambda\big)^{N_-(x)}. \tag{1}
\]
In particular, the perturbed energy law depends on \(x\) only through its
original energy \(H_C(x)\).  No further structural assumption on \(C\) is
needed for this fact.

Let
\[
Z_x=H_A(x)-\mu H_C(x)
    =\sum_{i<j}b_{ij}(\eta_{ij}-\mu),\qquad
D=\delta(1-\delta).
\]
For fixed \(0<\delta<1\) and \(\lambda=O(n^{-1/2})\), Taylor expansion of
(1), uniformly in \(x\), gives
\[
\log\mathbb E e^{\lambda Z_x}
=2DN\lambda^2+O_\delta(N|\lambda|^3)
=D n^2\lambda^2+O_\delta(\sqrt n). \tag{2}
\]
Chernoff optimization therefore yields, uniformly for bounded \(z>0\),
\[
\Pr\{Z_x\ge z n^{3/2}\}
\le
\exp\left\{-\,\frac{z^2}{4D}\,n+o_\delta(n)\right\}. \tag{3}
\]
The same estimate holds for the lower tail.  Formula (3) is the exact
speed-\(n\) moderate-deviation rate relevant to this problem.

## 2. Conditional stratified perturbation theorem

Define the flat-spectrum rate and entropy functions
\[
R_f(r)=\frac14\log\frac1{1-r^2},\qquad
s_f(r)=\log2-R_f(r)
      =\log2+\frac14\log(1-r^2), \tag{4}
\]
and
\[
r_*=\frac{\sqrt{15}}4,\qquad
c_*=\frac{r_*}{2}=\frac{\sqrt{15}}8.
\]
Notice that \(s_f(r_*)=0\).

Consider a sequence of conference signings \(C_n\).  The precise missing
energy-layer hypothesis is:

> **Flat-envelope with sparse exceptions (FESE).** There are sets
> \(\mathcal E_n\subset\{\pm1\}^n/\{\pm1\}\), with
> \(|\mathcal E_n|=\exp(o(n))\), such that, uniformly over polynomially many
> bins \(I\subset[0,r_*+o(1)]\),
> \[
> \#\{x\notin\mathcal E_n: |r_{C_n}(x)|\in I\}
> \le
> \exp\left\{n\left(\sup_{r\in I}s_f(r)+o(1)\right)\right\}, \tag{5}
> \]
> and every nonexceptional configuration has
> \(|r_{C_n}(x)|\le r_*+o(1)\).

The exceptional vectors may have energy all the way up to the spectral cap
\(|r|=1\).

Under FESE, independently flipping edges at fixed rate \(\delta\) and applying
(3) to each layer proves that some perturbed signing satisfies
\[
\frac{\max_x|H_A(x)|}{n^{3/2}}
\le G(\delta)+o(1), \tag{6}
\]
where
\[
\begin{split}
B(\delta)
&=\max_{0\le r\le r_*}
\left[
\frac{1-2\delta}{2}\,r
+2\sqrt{\delta(1-\delta)s_f(r)}
\right],\\
G(\delta)&=\max\left\{\frac12-\delta,\ B(\delta)\right\}. \tag{7}
\end{split}
\]

Indeed, (5) and (3) give the term \(B(\delta)\).  For
\(\exp(o(n))\) exceptional configurations, a union bound in (3) makes their
centered noise \(o(n^{3/2})\), while their mean is at most
\((1-2\delta)n^{3/2}/2\), giving \(1/2-\delta\).

This is not a heuristic use of Gaussian noise: (1)--(3) prove the rate
directly for Bernoulli edge flips.

## 3. Optimizing the conditional bound

The minimum of (7) occurs where the exceptional and ordinary-layer terms
meet.  At the maximizing ordinary energy \(r\), stationarity and equality
\(B(\delta)=1/2-\delta\) reduce to the scalar equation
\[
s_f(r)=\frac{r}{4(1+r)}. \tag{8}
\]
It has the relevant solution
\[
r=0.947800854926\ldots .
\]
If
\[
q=\frac{2\sqrt{s_f(r)}}{1-r},\qquad
1-2\delta=\frac{q}{\sqrt{1+q^2}}, \tag{9}
\]
then
\[
\delta_*=0.001394039184177\ldots
\]
and
\[
\boxed{
\min_{0\le\delta\le1/2}G(\delta)
=G(\delta_*)
=0.498605960815823\ldots <\frac12.} \tag{10}
\]

Therefore FESE for one ratio-dense sequence of conference orders, together
with the already available padding/prime-gap transfer, would make the
strict upper bound unconditional.

## 4. Why independent perturbation cannot reach the ROM value

The obstruction occurs immediately below the zero-entropy endpoint.  From
(4),
\[
s_f(r_*-u)=8r_*u+O(u^2). \tag{11}
\]
For \(\delta\downarrow0\), substitute \(u=32r_*\delta+o(\delta)\) into the
ordinary-layer objective in (7).  This gives
\[
B(\delta)
=c_*+15r_*\delta+o(\delta)
=c_*+14.523687548\ldots\,\delta+o(\delta). \tag{12}
\]
Thus every sufficiently small positive flip rate makes the bulk bound
strictly *worse* than \(c_*\), even though it damps each isolated resonance.

More starkly, damping a spectral-cap resonance from \(1/2\) down to \(c_*\)
requires
\[
\delta\ge\frac12-c_*
=0.0158770817241\ldots . \tag{13}
\]
At equality in (13), the ordinary-layer variational value is
\[
B(\delta)=0.567909095566\ldots, \tag{14}
\]
far above both \(c_*\) and \(1/2\).

Hence independent homogeneous edge flips cannot prove the conjectural
\(\sqrt{15}/8\) upper constant under the flat entropy law.  Their best
possible role is the much smaller strict improvement (10).

## 5. The exact missing deterministic theorem

The whole route now hinges on proving FESE, or a quantitatively weaker
variant, for a concrete algebraic sequence.  For Paley matrices the Fourier
identity has the form
\[
x^\top C_px
=\frac1p\sum_{\xi\in\mathbb F_p}\lambda_\xi
  |\widehat x(\xi)|^2,\qquad
\lambda_\xi\in\{\pm\sqrt p\}
\]
up to the standard zero-frequency convention.  Equivalently, the energy is
the imbalance of Fourier mass between the quadratic-residue and
nonresidue eigenspaces.

The known square waves put anomalously large Fourier mass on one arithmetic
eigenspace.  FESE asks for a counting theorem saying that all such
delocalized failures are confined to only \(\exp(o(p))\) Boolean vectors,
while the rest obey the flat rate (4).

A sufficient version is:
\[
\#\left\{x\notin\mathcal E_p:
\frac{|x^\top C_px|}{p^{3/2}}\ge r\right\}
\le
\exp\{p(s_f(r)+o(1))\},\qquad 0\le r\le r_*, \tag{15}
\]
with \(|\mathcal E_p|=\exp(o(p))\), and no nonexceptional vector above
\(r_*+o(1)\).

The resonance theorem alone does **not** imply the required exceptional-set
bound: it constructs subexponentially many explicit square waves but does
not classify every high-energy Boolean vector.  This classification/counting
step is exactly where character-sum delocalization, an inverse theorem for
Boolean Fourier mass, or a moment-generating estimate with resonance
excision must enter.

## 6. A resonant Hamming-cloud theorem: FESE is false

The sparse-exception premise in FESE cannot hold in the presence of even one
cap-near vector.  This follows from an exact averaging identity that does not
use the conference property.

Fix any signing \(C\), a vector \(x\), and choose \(y\) uniformly from the
Hamming sphere of radius \(k\) around \(x\).  Put \(z_i=x_iy_i\), so exactly
\(k\) coordinates of \(z\) equal \(-1\).  For every \(i\ne j\),
\[
\mathbb E z_i z_j
=\rho_{n,k}
:=\frac{(n-2k)^2-n}{n(n-1)}. \tag{16}
\]
It follows immediately that
\[
\mathbb E_y H_C(y)=\rho_{n,k}H_C(x),\qquad
\mathbb E_y r_C(y)=\rho_{n,k}r_C(x). \tag{17}
\]

Because \(r_C(y)\le1\), for every \(t<\rho_{n,k}r_C(x)\), the proportion of
the Hamming sphere satisfying \(r_C(y)\ge t\) is at least
\[
\frac{\rho_{n,k}r_C(x)-t}{1-t}. \tag{18}
\]
Indeed, if that proportion is \(p\), then the sphere average is at most
\(p+(1-p)t\).

Suppose now that \(r_C(x)=1-o(1)\), take \(k=\alpha n+o(n)\), and let
\(t<(1-2\alpha)^2\) be fixed.  Equations (18) and Stirling's formula give
\[
\#\{y:r_C(y)\ge t\}
\ge
\exp\{n(h(\alpha)+o(1))\}, \tag{19}
\]
where
\[
h(\alpha)=-\alpha\log\alpha-(1-\alpha)\log(1-\alpha).
\]

At the flat-rate endpoint \(r_*=\sqrt{15}/4\), set
\[
\alpha_*=\frac{1-\sqrt{r_*}}2
=0.008002582184\ldots .
\]
Taking a threshold just below \(r_*\), (19) gives layer entropy arbitrarily
close to
\[
h(\alpha_*)=0.046606870654\ldots, \tag{20}
\]
whereas \(s_f(r_*)=0\).  Removing only \(\exp(o(n))\) configurations cannot
repair this discrepancy.

This applies directly to exact Boolean conference eigenvectors and, with
the evident \(r_C(x)=1-o(1)\) modification, to cap-approaching Paley
square-wave sequences.  Thus the proposed FESE theorem is not merely
unproved; it is false on precisely the resonant sequences it was intended
to handle.

## 7. Why the Hamming cloud does not yet kill perturbation

Applying the independent per-vector union bound (3) to the cloud (19) is
too pessimistic.  The centered perturbation process has the exact covariance
\[
\operatorname{Cov}(Z_x,Z_y)
=4D\sum_{i<j}x_ix_jy_iy_j
=2D\big((x^\top y)^2-n\big). \tag{21}
\]
At overlap \(q=x^\top y/n\), its correlation is \(q^2+o(1)\).  A small
Hamming cloud therefore shares almost all of its perturbation noise.

The corrected target is a two-parameter estimate for
\[
\max\{Z_y-Z_x:
  d_H(x,y)=\alpha n,\ r_C(y)\ \text{in a specified layer}\}, \tag{22}
\]
not a union bound on the absolute variables \(Z_y\).  For a flip set
\(S=\{i:y_i\ne x_i\}\),
\[
Z_y-Z_x
=-2\sum_{i\in S,\ j\notin S}
  c_{ij}x_ix_j(\eta_{ij}-\mu), \tag{23}
\]
so the local process is a signed random-cut process.  Its canonical metric
contains an extra factor determined by symmetric differences of vertex
sets; this is the correlation that the scalar entropy formula (7) discards.

Consequently:

- the numerical value (10) remains a correct consequence of FESE, but FESE
  cannot describe resonant Paley orders;
- the unconditional perturbation problem is now a **resonant-center plus
  random-cut chaining** problem;
- a valid proof must cover each cap-near center by Hamming shells, control
  (23) relative to the center, and separately apply the flat layer estimate
  away from all such clouds.

## 8. Exact canonical metric and Gaussian-width bounds

The correlated-noise penalty in (22) can be quantified sharply, independently
of \(C\).  Replace the centered edge variables temporarily by independent
Gaussians \(W_{ij}\sim N(0,4D)\), and define
\[
Y_S=-2\sum_{i\in S,\ j\notin S}W_{ij},
\qquad |S|=k.
\]
For two \(k\)-subsets \(S,T\), the symmetric difference of their cut-edge
sets is the cut generated by \(S\triangle T\).  If
\(d=|S\triangle T|\), then
\[
\boxed{
\mathbb E(Y_S-Y_T)^2=16D\,d(n-d).} \tag{24}
\]
This is the exact canonical metric of the resonant-cloud perturbation.

Assume \(k/n\to\alpha\le1/2\), let
\[
z_\alpha=\Phi^{-1}(1-\alpha),
\]
and let \(\phi,\Phi\) denote the standard normal density and distribution.
Compare \(Y\) with the two linear Gaussian processes
\[
L_S^+=4\sqrt{Dn}\sum_{i\in S}g_i,\qquad
L_S^-=4\sqrt{D(n-2k)}\sum_{i\in S}g_i.
\]
Their squared increment metrics are respectively
\[
16Dn\,d,\qquad 16D(n-2k)d.
\]
Since \(d\le2k\),
\[
\mathbb E(L_S^--L_T^-)^2
\le\mathbb E(Y_S-Y_T)^2
\le\mathbb E(L_S^+-L_T^+)^2.
\]
Sudakov--Fernique and the fact that the largest sum of \(k\) coordinates is
the sum of the \(k\) largest order statistics give
\[
\boxed{
4\sqrt{D(1-2\alpha)}\,\phi(z_\alpha)
\le
\liminf\frac{\mathbb E\max_{|S|=k}Y_S}{n^{3/2}}
\le
\limsup\frac{\mathbb E\max_{|S|=k}Y_S}{n^{3/2}}
\le
4\sqrt D\,\phi(z_\alpha).} \tag{25}
\]
For \(\alpha\downarrow0\), the ratio of the lower and upper bounds tends to
one.  Thus the small-cloud width is asymptotically determined:
\[
\frac{\mathbb E\max Y_S}{n^{3/2}}
=4\sqrt D\,\phi(z_\alpha)(1+O(\alpha)). \tag{26}
\]

This Gaussian answer is universal for the centered Bernoulli flip variables
at the \(n^{3/2}\) scale.  One elementary proof replaces the \(O(n^2)\) edge
variables one at a time in the soft maximum
\[
\beta^{-1}\log\sum_{|S|=k}e^{\beta Y_S}.
\]
The soft-max error is \(O(n/\beta)\), while the third-order Lindeberg
remainder is \(O(n^2\beta^2)\).  Taking
\(\beta=n^{-3/8}\) makes both \(o(n^{3/2})\).  The edge variables have
matching mean and variance and uniformly bounded third moments.  Therefore
(25) also holds, up to \(o(1)\), for the actual centered edge-flip process.
Bounded differences further shows concentration on the \(o(n^{3/2})\)
scale.

### Comparison with scalar entropy penalties

A per-vector union bound on the absolute noise, which led to (7), charges
\[
2\sqrt{D\,h(\alpha)}.
\]
Even a union bound on increments inside the shell charges
\[
4\sqrt{2D\alpha(1-\alpha)h(\alpha)}.
\]
The actual small-shell Gaussian width is instead
\[
4\sqrt D\,\phi(z_\alpha)(1+O(\alpha)).
\]
At the shell that reaches the ROM endpoint from a cap vector,
\(\alpha_* =0.008002582184\ldots\), the three coefficients after dividing
by \(\sqrt D\) are respectively
\[
0.4317724894,\qquad 0.1088103741,\qquad
0.0870\text{--}0.0877031470. \tag{27}
\]
Thus geometry removes about eighty percent of the original scalar
energy-layer penalty.  The Hamming cloud invalidates FESE, but it does not
invalidate sparse perturbation.

The exact limit in (25) is the zero-temperature constrained SK/random-cut
energy at magnetization \(1-2\alpha\).  A Parisi formula can characterize
it, but (25) is already nearly sharp at the resonant radii relevant here and
requires no spin-glass machinery.

## 9. The corrected variational frontier

For a cap-near center \(x\), write
\[
u_C(\alpha)=
\frac1{n^{3/2}}\max_{d_H(x,y)=\alpha n}H_C(y).
\]
Equations (23)--(26) imply the rigorous shell upper bound
\[
\frac1{n^{3/2}}\max_{d_H(x,y)=\alpha n}H_A(y)
\le
(1-2\delta)u_C(\alpha)
+4\sqrt{\delta(1-\delta)}\,\phi(z_\alpha)+o(1). \tag{28}
\]
This formulation isolates the remaining deterministic input.  The spectral
bound gives only \(u_C(\alpha)\le1/2\), which is too weak.  A strict
perturbative construction follows from any uniform shell profile
\(\bar u(\alpha)\) for which
\[
\inf_{\delta>0}\sup_{0\le\alpha\le1/2}
\left[
(1-2\delta)\bar u(\alpha)
+4\sqrt{\delta(1-\delta)}\,\phi(z_\alpha)
\right]<\frac12. \tag{29}
\]
The average identity (17) says that the shell mean is
\(\tfrac12(1-2\alpha)^2+o(1)\) around an exact cap vector, but an upper
profile for its maximum is not presently known.  This maximum-profile
theorem, rather than FESE, is the exact missing deterministic statement.

## 10. Paley square waves make every one-center shell flat

The deterministic profile \(u_C(\alpha)\) in (28) cannot decay when it is
defined relative to only one Paley square wave.  This is an exact symmetry
obstruction.

For \(p=4m+1\), write
\[
s(t)=\operatorname{sgn}\cos(2\pi t/p),\qquad
x_{k,a}(j)=s(kj+a).
\]
The positive set of \(s\) is the cyclic interval
\(\{-m,-m+1,\ldots,m\}\), of size \((p+1)/2\).  If
\(\chi_p(k)=1\), all \(x_{k,a}\) have the same Paley energy: translations
and multiplication by a quadratic residue are automorphisms of the Paley
signing.

For fixed \(k\), put
\[
\ell=\min(|a-b|_p,p-|a-b|_p).
\]
For \(0\le\ell\le(p-1)/2\), the two cyclic half-intervals have symmetric
difference
\[
d_H(x_{k,a},x_{k,b})=2\ell. \tag{30}
\]
After identifying a vector with its global negative, these equal-energy
translates occur at every normalized Hamming radius in \([0,1/2]\), with
mesh \(O(1/p)\).

Consequently, if a square-wave sequence has
\[
\frac{H_C(x_{k,a})}{p^{3/2}}\longrightarrow\frac12,
\]
then the one-center shell profile satisfies
\[
\liminf_p u_C(\alpha)\ge\frac12
\qquad\text{for every fixed }0\le\alpha\le\frac12. \tag{31}
\]
Since the spectral bound gives the reverse inequality, the profile is flat.
Thus no theorem of the form
\(u_C(\alpha)\le1/2-c\alpha\) can hold around a single resonant template.

This does **not** defeat random perturbation.  The correct resonant set is
the complete affine orbit
\[
\mathcal R_p=
\{\pm x_{k,a}:a\in\mathbb F_p,\ \chi_p(k)=1\},
\qquad |\mathcal R_p|=O(p^2)=\exp(o(p)). \tag{32}
\]
The maximum centered perturbation over \(\mathcal R_p\) is
\(o(p^{3/2})\), so all members of the orbit are damped together by
\(1-2\delta\).

The remaining deterministic target must therefore use distance to the
whole orbit:
\[
u_C^{\rm orb}(\alpha)=
\frac1{p^{3/2}}
\max_{\substack{x:\,
 d_H(x,\mathcal R_p)=\alpha p}}
H_C(x). \tag{33}
\]
A useful decay estimate for (33) is an inverse theorem: a Boolean vector
whose Fourier mass is almost entirely on the quadratic-residue eigenspace
must be close to an affine square wave.  In projection language, it is a
stability/restricted-projection theorem for the Paley Fourier subspace.
Ordinary spectral information cannot prove it, and exact additive
uncertainty only rules out a genuinely sparse vector lying *exactly* in the
eigenspace; it supplies no constant quantitative gap at linear sparsity.

## 11. Two resonant templates create a second exponential obstruction

Distance to the full affine orbit removes the trivial flatness (31), but it
still does not give the mean Hamming-shell decay.  Two cap-near vectors can
be patched coordinatewise to create exponentially many new high-energy
vectors.

Let \(U=C/\sqrt{n-1}\), \(P_-=(I-U)/2\), and suppose
\[
\frac{x^\top Ux}{n},\frac{z^\top Uz}{n}\ge1-\varepsilon.
\]
Let
\[
D_0=\{i:x_i\ne z_i\},\quad |D_0|=\gamma n,\quad
a=\frac{x+z}{2},\quad b=\frac{x-z}{2}.
\]
The vectors \(a,b\) have disjoint supports, and
\[
\|P_-a\|^2
\le\frac{\|P_-x\|^2+\|P_-z\|^2}{2}
\le\frac{\varepsilon n}{2}. \tag{34}
\]
Choose an arbitrary sign vector \(w\) supported on \(D_0\) and set
\(y=a+w\).  This gives all \(2^{\gamma n}\) Boolean coordinatewise
patchworks of \(x,z\).  Averaging uniformly over \(w\), the mixed term
vanishes and, since \(U\) has zero diagonal,
\[
\mathbb E_w w^\top Uw=\operatorname{tr}U[D_0,D_0]=0.
\]
Therefore
\[
\mathbb E_w\frac{y^\top Uy}{n}
=\frac{a^\top Ua}{n}
=\frac{\|a\|^2-2\|P_-a\|^2}{n}
\ge1-\gamma-\varepsilon. \tag{35}
\]
As every normalized energy is at most one, a constant fraction of the
patchworks have energy at least
\(1-\gamma-\varepsilon-o(1)\).

Take \(x,z\) to be translated square waves at distance \(\gamma n\).  A
uniform patchwork differs from each endpoint in
\(\gamma n/2+O(\sqrt n)\) coordinates.  More generally, its agreement with
any fixed member of the \(O(n^2)\)-sized affine orbit is a sum of independent
signs on \(D_0\).  Hoeffding plus a union bound shows that, except for an
exponentially small fraction, its distance from the complete orbit is
\(\gamma n/2+o(n)\), provided the endpoint pair is chosen at the relevant
orbit separation.  Thus, with \(\alpha=\gamma/2\),
\[
\boxed{
u_C^{\rm orb}(\alpha)
\ge\frac12-\alpha-o(1),} \tag{36}
\]
and the corresponding orbit-relative layer has entropy at least
\[
2\alpha\log2-o(1). \tag{37}
\]

This halves the slope suggested by the one-center shell average
\(\tfrac12(1-2\alpha)^2=\tfrac12-2\alpha+O(\alpha^2)\).
Accordingly, an orbit-stability theorem stronger than (36) is false.

The perturbation noise on this patchwork family is again much more correlated
than its entropy suggests.  Conditional on the coordinates outside \(D_0\),
the cross edges generate independent Gaussian local fields of variance
\(4D(n-|D_0|)\).  Their leading width over all \(2^{|D_0|}\) completions is
\[
2\gamma\sqrt{\frac{2D(1-\gamma)}{\pi}}\,n^{3/2}, \tag{38}
\]
with an additional internal constrained-SK term of order
\(O(\sqrt D\,\gamma^{3/2}n^{3/2})\).  Hence the patchwork penalty is linear
in \(\alpha\) (up to the smaller \(\alpha^{3/2}\) term), rather than the
\(\sqrt\alpha\) penalty from a scalar entropy union bound.

For very sparse \(\delta\), the deterministic loss in (36) still dominates
(38).  Thus (36) is a genuine obstruction to overly strong stability, but
does not by itself rule out a strict perturbed upper bound below \(1/2\).

## 12. Exact patchwork width, optimized model, and \(k\)-template hull

Let \(\gamma=2\alpha\).  The full two-template patchwork family is indexed by
all subsets of a disagreement set of size \(\gamma n\).  Applying the same
metric comparison as in Section 8, now without a cardinality constraint,
gives the sharper all-subcube bounds
\[
\boxed{
2\gamma\sqrt{\frac{2D(1-\gamma)}{\pi}}
\le W_\gamma(D)
\le
2\gamma\sqrt{\frac{2D}{\pi}},} \tag{39}
\]
where \(W_\gamma(D)\) is the limiting \(n^{-3/2}\) expected supremum of the
centered flip noise on that patchwork subcube.

Indeed, compare with
\[
4\sqrt{D(n-\gamma n)}\sum_{i\in S}g_i
\quad\text{and}\quad
4\sqrt{Dn}\sum_{i\in S}g_i,
\qquad S\subseteq D_0,
\]
and use
\[
\mathbb E\max_{S\subseteq D_0}\sum_{i\in S}g_i
=|D_0|\,\mathbb E(g_+)
=\frac{\gamma n}{\sqrt{2\pi}}.
\]

There is also an exact, though implicit, description.  If
\(\mathcal P_{\rm SK}(h)\) denotes the zero-temperature Parisi ground energy
per spin of the SK model with iid Gaussian external field of strength \(h\),
then
\[
W_\gamma(D)
=2\sqrt D\,\gamma^{3/2}
\mathcal P_{\rm SK}\left(\sqrt{\frac{1-\gamma}{\gamma}}\right). \tag{40}
\]
The external-field asymptotic of (40) agrees with (38) as
\(\gamma\downarrow0\).

### Exact optimization under the extremal mean profile

It is useful to optimize the model in which the deterministic shell profile
equals the patchwork obstruction
\[
u(\alpha)=\frac12-\alpha=\frac{1-\gamma}{2}. \tag{41}
\]
This is **not yet a proved upper profile for every vector in the shell**;
it is the sharp guaranteed mean/convex-hull profile.  The calculation below
therefore diagnoses whether patchwork geometry itself obstructs
perturbation.

Put
\[
a=\frac12-\delta,\qquad
b=2\sqrt{\frac{2\delta(1-\delta)}{\pi}}.
\]
By (39), the combined deterministic-plus-noise model is bracketed by
\[
a(1-\gamma)+b\gamma\sqrt{1-\gamma}
\quad\text{and}\quad
a(1-\gamma)+b\gamma. \tag{42}
\]
Choose \(a=b\).  This gives
\[
\delta_*=
\frac12-\sqrt{\frac{2}{\pi+8}}
=0.0763166866997\ldots . \tag{43}
\]
At this value, the upper expression in (42) is identically \(a\) for every
\(\gamma\), and \(\gamma=0\) attains it.  Therefore the exact patchwork model
optimization is
\[
\boxed{
\inf_\delta\sup_\gamma
\left[(1-2\delta)\frac{1-\gamma}{2}+W_\gamma(D)\right]
\le
\sqrt{\frac{2}{\pi+8}}
=0.423683313300\ldots<\frac12.} \tag{44}
\]
In particular, two-template clouds do not constitute a variational
obstruction to a strict improvement.

### Why more templates do not flatten the guaranteed mean profile

Let \(x^1,\ldots,x^k\) be cap-near templates, let
\(\lambda_j\ge0\), \(\sum_j\lambda_j=1\), and put
\[
m=\sum_j\lambda_jx^j.
\]
Generate a Boolean patchwork \(y\) independently coordinate by coordinate
with \(\mathbb Ey_i=m_i\).  Zero diagonal gives
\[
\mathbb E\frac{y^\top Uy}{n}
=\frac{m^\top Um}{n}
=\frac{\|m\|^2}{n}-o(1). \tag{45}
\]
Let
\[
q_{\max}=\max_j\frac{\langle m,x^j\rangle}{n}.
\]
Since \(m\) is a convex combination of the templates,
\[
\frac{\|m\|^2}{n}
=\sum_j\lambda_j\frac{\langle m,x^j\rangle}{n}
\le q_{\max}. \tag{46}
\]
A typical patchwork has distance
\(\alpha n+o(n)\) from the nearest template, where
\[
\alpha=\frac{1-q_{\max}}2.
\]
Equations (45)--(46) imply
\[
\mathbb E\frac{H_C(y)}{n^{3/2}}
\le\frac12-\alpha+o(1). \tag{47}
\]
Equal mixtures of two appropriately separated square-wave translates attain
equality.  Thus the limiting convex-hull/mean profile generated by arbitrary
\(k\)-template patching is exactly \(1/2-\alpha\); adding more templates
cannot create a flatter guaranteed mean-energy obstruction.

What remains possible is a sparse exceptional subset of a patchwork cloud
whose energy lies above its mean profile.  Controlling that upper tail,
rather than adding further templates, is the unresolved deterministic task.

## 13. Verdict

- **Rigorous and unconditional:** exact conditional law (1), uniform
  moderate-deviation rate (3), perturbation variational formula (7), the
  obstruction (12)--(14), the resonant-cloud theorem (16)--(20), the exact
  metric (24), and the Gaussian-width bounds (25).
- **Rigorous conditional consequence:** FESE implies the new upper constant
  \(0.498605960816\ldots\).
- **Rigorous negative result:** FESE is false for every sequence containing
  a cap-near Boolean resonance.
- **Still missing:** a deterministic shell-maximum bound strong enough to
  make (29) less than \(1/2\), formulated relative to the full affine orbit
  as in (33).  The correlated-noise side is now controlled within a
  \(1+O(\alpha)\) factor on small resonant shells.

Independent sparse perturbation is not settled.  Pure energy-layer counting
cannot handle arithmetic resonances, but the exact covariance (21) leaves a
more structured route open.  Any future strict improvement below \(1/2\)
must exploit this geometry rather than treating the exponentially many
nearby vectors as independent.
