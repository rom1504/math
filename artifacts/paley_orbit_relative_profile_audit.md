# Orbit-relative Paley profile: patchwork theorem and counterexamples

## Status

This note audits the proposed orbit-relative shell inequality

\[
u_p^{\rm orb}(\alpha)\stackrel{?}{\le}\frac12-\alpha+o(1)
\tag{1}
\]

around the affine square-wave orbit.

The two-template patchwork calculation is correct and its optimized
correlated-noise constant is

\[
\sqrt{\frac2{\pi+8}}=0.423683313300\ldots .
\]

However, (1) is false when the orbit contains only affine
single-frequency square waves.  It already fails strongly for the Paley
matrix at \(p=17\), and a finite-cyclic \(3\)-adic gadget gives an
asymptotic second structured orbit far from every ordinary square wave.

The correct target must therefore use a growing arithmetic template library,
not one affine orbit.

## 1. Translation flatness and the need to quotient the full orbit

For

\[
x_{k,a}(j)=\operatorname{sgn}
\cos\left(\frac{2\pi k(j-a)}p\right),
\]

translation invariance of the Paley matrix makes the energy independent of
\(a\).  A direct interval count gives

\[
d_H(x_{k,a},x_{k,a+\ell})=2\ell
\qquad
0\le\ell\le\frac{p-1}{4}. \tag{2}
\]

Thus a single resonant square wave has equal-energy translates at every
normalized Hamming radius in \([0,1/2]\), up to \(O(1/p)\) mesh.
A one-center shell-decay theorem is impossible.

The natural first quotient is

\[
\mathcal R_p^{\rm sq}
=\{\pm x_{k,a}:a\in\mathbb F_p,\ \chi_p(k)=1\},
\qquad
|\mathcal R_p^{\rm sq}|=O(p^2). \tag{3}
\]

A polynomial number of centers does not change the leading perturbation
width: unioning the shell bounds over (3) costs only
\(O(n\sqrt{\log n})=o(n^{3/2})\).

## 2. Independent verification of the two-template patchwork theorem

Let \(U=C/\sqrt{n-1}\), so \(U^2=I\), and put
\(P_-=(I-U)/2\).  Suppose two Boolean templates satisfy

\[
\frac{x^\top Ux}{n},\frac{z^\top Uz}{n}\ge1-\varepsilon.
\]

If

\[
D_0=\{i:x_i\ne z_i\},\quad |D_0|=\gamma n,\quad
a=\frac{x+z}{2},
\]

then

\[
\|P_-a\|^2
\le\frac{\|P_-x\|^2+\|P_-z\|^2}{2}
\le\frac{\varepsilon n}{2}. \tag{4}
\]

Choose independent signs \(w_i\) on \(D_0\), put \(w=0\) off \(D_0\),
and set \(y=a+w\).  Zero diagonal gives

\[
\mathbb E_w y^\top Uy=a^\top Ua
=\|a\|^2-2\|P_-a\|^2
\ge(1-\gamma-\varepsilon)n. \tag{5}
\]

This energy actually concentrates around its mean.  The mixed term has
variance at most \(4\|Ua\|^2=O(n)\), and

\[
\operatorname{Var}(w^\top Uw)
\le 2\sum_{i,j\in D_0}U_{ij}^2=O(n).
\]

Thus all but \(o(1)\) of the patchworks have doubled normalized energy at
least \(1-\gamma-\varepsilon-o(1)\).

They are at distance \(\gamma n/2+o(n)\) from the nearer endpoint and, after
a polynomial union bound, from the full affine orbit when the two endpoints
are chosen at the corresponding separation.  With \(\alpha=\gamma/2\),

\[
\boxed{\displaystyle
u_p^{\rm orb}(\alpha)\ge\frac12-\alpha-o(1).} \tag{6}
\]

So no orbit-relative upper profile strictly below this line can hold.

## 3. Exact patchwork-noise width and the \(0.4236833\) diagnostic

For a disagreement set of size \(\gamma n\), the complete patchwork family
is the subcube \(S\subseteq D_0\).  If \(d=|S\triangle T|\), the centered
edge-perturbation process has metric

\[
\mathbb E(Y_S-Y_T)^2=16D\,d(n-d),
\qquad D=\delta(1-\delta). \tag{7}
\]

Because \(d\le\gamma n\), comparison with linear Gaussian processes gives

\[
2\gamma\sqrt{\frac{2D(1-\gamma)}{\pi}}
\le W_\gamma(D)
\le2\gamma\sqrt{\frac{2D}{\pi}}. \tag{8}
\]

The internal SK term is already included in this metric comparison; it
must not be added a second time.

Under the diagnostic deterministic profile

\[
u(\alpha)=\frac12-\alpha=\frac{1-\gamma}{2},
\]

write

\[
a=\frac12-\delta,\qquad
b=2\sqrt{\frac{2\delta(1-\delta)}{\pi}}.
\]

The upper comparison in (8) gives the affine envelope

\[
a(1-\gamma)+b\gamma.
\]

Choosing \(a=b\) yields

\[
\delta_*=\frac12-\sqrt{\frac2{\pi+8}}
=0.0763166866997\ldots
\]

and

\[
\boxed{\displaystyle
\sup_\gamma\bigl[a(1-\gamma)+W_\gamma(D)\bigr]
\le
\sqrt{\frac2{\pi+8}}
=0.423683313300\ldots .} \tag{9}
\]

This verifies the constant, but only for the patchwork model or for a
future theorem classifying the whole dangerous shell by such subcubes.
The deterministic line (1) alone does not imply (9), because arbitrary
shell vectors have the larger fixed-cardinality random-cut width.

## 4. More templates do not worsen the convex-hull mean line

Let \(x^1,\ldots,x^k\) be cap-near templates, choose weights
\(\lambda_j\), and put

\[
m=\sum_j\lambda_jx^j.
\]

Generate \(y\) independently coordinatewise with \(\mathbb Ey_i=m_i\).
Convexity gives \(\|P_-m\|^2=o(n)\), so

\[
\mathbb E\frac{y^\top Uy}{2n}
=\frac{\|m\|^2}{2n}+o(1). \tag{10}
\]

Put

\[
q_{\max}=\max_j\frac{\langle m,x^j\rangle}{n}.
\]

Since \(m\) is their convex combination,

\[
\frac{\|m\|^2}{n}
=\sum_j\lambda_j\frac{\langle m,x^j\rangle}{n}
\le q_{\max}. \tag{11}
\]

The distance of a typical patchwork to its nearest template is

\[
\alpha=\frac{1-q_{\max}}2+o(1).
\]

Equations (10)--(11) imply

\[
\mathbb E\frac{y^\top Uy}{2n}
\le\frac12-\alpha+o(1). \tag{12}
\]

Two equally weighted, appropriately separated translates attain equality.
Thus \(k\)-template independent convexification cannot create a guaranteed
mean profile flatter than (6).  The unresolved danger is a sparse
above-mean subset of a patchwork family, not the number of templates in the
mixture.

## 5. Finite counterexample to the pure square-wave orbit

For the \(17\times17\) Paley core, exhaustive enumeration gives the Boolean
vector

```text
++ - + - + - +++ ------
```

with spaces removed, equivalently

```text
++-+-+-+++------
```

and

\[
H_{C_{17}}(x)=32,\qquad
\frac{H_{C_{17}}(x)}{17^{3/2}}
=0.456537647127\ldots . \tag{13}
\]

Its distance from the complete \(\pm\), translate, QR-frequency square-wave
orbit of 136 distinct vectors is \(4\).  Hence

\[
\frac{H_{C_{17}}(x)}{17^{3/2}}
>
\frac12-\frac4{17}
=0.264705882353\ldots . \tag{14}
\]

There are 2516 Boolean vectors violating the finite-\(p\) version of (1).
This does not by itself disprove an asymptotic theorem, but it shows that the
square-wave orbit misses a major structured energy family already at the
first informative order.

## 6. Asymptotic second orbit: the \(3\)-adic length-27 gadget

The depth-one \(3\)-adic gadget from `paley_resonance_gadget.md` is the
length-27 word

```text
+-++-++-++-++-++---+--+--+
```

viewed as a step function on 27 equal intervals.  On suitable fixed
arithmetic refinements of the primes, its doubled normalized Paley energy
tends

\[
\frac{80}{81}=0.987654320988\ldots>r_* . \tag{15}
\]

It is not close to an affine single-frequency square wave.  For the
continuous word \(f\),

\[
\sup_{k\ge1,\ \varphi\in\mathbb R}
\int_0^1
f(t)\operatorname{sgn}\cos(2\pi kt-\varphi)\,dt
=\frac49. \tag{16}
\]

Here is a finite exact verification:

- for fixed \(k\), the correlation as a function of phase is piecewise
  linear, so its maximum occurs when a cosine zero meets one of the 27 word
  boundaries;
- enumerating those rational events for \(1\le k\le30\) gives maximum
  \(4/9\), attained at \(k=8\);
- on each constant word interval, the integral of a frequency-\(k\) square
  wave has magnitude at most \(1/(2k)\), so for \(k\ge31\) the total
  correlation is at most
  \[
  \frac{27}{2k}\le\frac{27}{62}<\frac49.
  \]

Consequently its limiting Hamming distance from every continuous
single-frequency square wave is at least

\[
\boxed{\displaystyle
\frac{1-4/9}{2}=\frac5{18}=0.277777\ldots .} \tag{17}
\]

Passing (16) uniformly to every discrete affine frequency requires a
rank-one lattice discrepancy estimate when the frequency varies with
\(p\).  The fixed-frequency and high-frequency regimes above already show
the structural issue; a fully uniform discrete statement should be included
in any published version.

The finite-cyclic gadget therefore supplies a second natural arithmetic
template family, asymptotically cap-near but geometrically separated from
ordinary square waves.

## 7. Exact weakest profile needed by the general shell-width method

For a template library \(\mathcal T_p\), write

\[
\bar u(\alpha)=\frac12-\Delta(\alpha).
\]

The general correlated-shell upper width from
`paley_correlated_shell_optimization_audit.md` yields a strict improvement
whenever some fixed \(\delta>0\) satisfies

\[
\boxed{\displaystyle
(1-2\delta)\Delta(\alpha)+\delta
>
4\sqrt{\delta(1-\delta)}\,
\phi(\Phi^{-1}(1-\alpha))
\quad(0<\alpha\le1/2).} \tag{18}
\]

Equivalently, the weakest pointwise gap for that \(\delta\) is

\[
\Delta_\delta(\alpha)
=
\frac{\left[
4\sqrt{\delta(1-\delta)}\phi(z_\alpha)-\delta
\right]_+}{1-2\delta}. \tag{19}
\]

Any uniform linear rigidity

\[
\Delta(\alpha)\ge\kappa\alpha
\]

with fixed \(\kappa>0\) is sufficient after choosing \(\delta\) small
enough.  For practical rates, the required slopes are approximately

\[
\begin{array}{c|c|c}
\delta&\kappa_{\rm required}&\text{worst }\alpha\\ \hline
10^{-3}&0.2567005&0.0213645\\
10^{-4}&0.0971313&0.0075928\\
10^{-6}&0.0124736&0.0009093
\end{array} \tag{20}
\]

The ideal line \(\Delta(\alpha)=\alpha\), if proved relative to a complete
template library, gives through the **general** shell-width bound

\[
\delta_{\rm opt}=0.01587069\ldots,\qquad
\alpha_{\rm worst}=0.02637874\ldots,
\qquad
c_{\rm upper}=0.489147479915\ldots . \tag{21}
\]

The stronger \(0.4236833\) value in (9) additionally uses the exact subcube
classification and its smaller width.

## 8. Fourier formulation of the missing inverse theorem

For the Paley core, let \(M_-(y)\) be the normalized Fourier mass of \(y\)
on the quadratic-nonresidue eigenspace and let

\[
m_0(y)=|\widehat y(0)|^2.
\]

The exact identity is

\[
\frac12-\frac{H_C(y)}{p^{3/2}}
=M_-(y)+\frac12m_0(y). \tag{22}
\]

Thus a linear orbit-relative profile is precisely the Boolean inverse
inequality

\[
M_-(y)+\frac12m_0(y)
\ge
\kappa\,\frac{d_H(y,\mathcal T_p)}p-o(1). \tag{23}
\]

The pure square-wave choice
\(\mathcal T_p=\mathcal R_p^{\rm sq}\) fails by Sections 5--6.
The live target is to construct a subexponential arithmetic template
library containing square waves and finite-cyclic gadgets, then prove
(23) for some \(\kappa>0\).

No positive-entropy family of mutually orbit-separated cap-near templates
has been found.  The known gadget constructions use only
\(\exp(o(p))\) templates at any fixed target accuracy, so the enlarged
template-library strategy remains plausible.

