# Paley orbit stability: first counterexample and corrected two-phase model

## 1. The proposed one-phase inverse theorem is false

Let \(C_p\) be the \(p\times p\) Paley signing
\[
(C_p)_{ij}=\chi_p(i-j),\qquad (C_p)_{ii}=0,
\]
and let \(\mathcal R_p\) be the affine orbit of the QR square waves
\[
x_{k,a}(j)=\operatorname{sgn}\cos(2\pi(kj+a)/p),
\qquad \chi_p(k)=1,
\]
including global negatives.  Define
\[
u_p(d)=
\max_{\substack{x\in\{\pm1\}^p\\
 d_H(x,\mathcal R_p)=d}}
\frac{H_{C_p}(x)}{p^{3/2}}.
\]

The tempting inverse inequality
\[
u_p(d)\le\frac12-\frac d p+o(1) \tag{1}
\]
is already false at \(p=17\).  Exhaustive enumeration gives:

| orbit distance \(d\) | exact maximum \(H\) | \(H/17^{3/2}\) | \(1/2-d/17\) |
|---:|---:|---:|---:|
| 0 | 24 | 0.342403 | 0.500000 |
| 1 | 28 | 0.399470 | 0.441176 |
| 2 | 28 | 0.399470 | 0.382353 |
| 3 | 28 | 0.399470 | 0.323529 |
| 4 | 32 | 0.456538 | 0.264706 |

Thus the first strict violation is at \(d=2\), and the true Paley optimizer
is at distance \(4\) from the square-wave orbit.  By contrast, the complete
\(p=13\) cube satisfies (1), and its square wave is an optimizer.

The exhaustive implementation is `paley_orbit_profile.cpp`.

## 2. Classification of the \(p=17\) extremizer

One optimizer has minus set
\[
\{0,1,2,4,6,8,9,10\}\subset\mathbb F_{17}
\]
(bit mask \(1879\)).  All \(34\) positive-energy maximizers modulo global
sign form one affine-QR orbit; the orbit has \(68\) members after including
global negatives.

Its normalized QR Fourier mass is
\[
0.954807543321,
\]
its NQR mass is
\[
0.041732249067,
\]
and its zero-frequency mass is \(0.003460207612\).  The largest Fourier
coefficients occur at
\[
\{\pm2,\pm8\}\quad\text{and}\quad\{\pm1,\pm4\},
\]
all quadratic-residue frequencies.  It is therefore a multi-harmonic QR
threshold, rather than a one-frequency square wave.

If \(U=C_{17}/\sqrt{17}\) and \(P_+=(I+U)/2\), the optimizer satisfies
\[
x=\operatorname{sgn}(P_+x). \tag{2}
\]
Its switched local fields are multiples of \(2/\sqrt{17}\):
\[
x_i(Ux)_i\in
\left\{0,\frac2{\sqrt{17}},\frac4{\sqrt{17}},\frac6{\sqrt{17}}\right\}.
\]
This identifies the next finite-order template family: affine orbits of
multi-harmonic fixed points of the projection-sign map (2).

The full high-energy histogram modulo global sign is
\[
\#\{H=20,24,28,32\}=\{2516,1224,272,34\}.
\]
Thus the \(p=17\) optimizer is a small algebraic orbit sitting above a much
larger bulk.

## 3. Why the global linear profile had to fail

The bound \(1/2-\alpha\) was derived as a convex-hull **mean** profile for
coordinatewise patchworks of cap-near templates.  It was never an upper
bound for every Boolean vector at that orbit distance.

Far from a polynomial-sized arithmetic orbit, a ROM-like bulk can still
have energy of order \(c_{\rm bulk}p^{3/2}\).  Consequently the minimal
plausible corrected envelope is two-phase:
\[
\boxed{
u^{\rm corr}(\alpha)
=\max\left\{\frac12-\alpha,\ c_{\rm bulk}\right\}.} \tag{3}
\]
The arithmetic ridge controls only
\[
0\le\alpha\le\alpha_0:=\frac12-c_{\rm bulk};
\]
the bulk plateau controls larger distances.

For the ROM candidate
\[
c_{\rm bulk}=c_*=\frac{\sqrt{15}}8,
\qquad
\alpha_0=0.0158770817241\ldots . \tag{4}
\]
The \(p=17\) data have exactly this qualitative form: the multi-harmonic
optimizer is a bulk extremizer outside the square-wave ridge.

## 4. Correlated perturbation of the arithmetic ridge

Let
\[
\mu=1-2\delta,\qquad D=\delta(1-\delta),\qquad
z_\alpha=\Phi^{-1}(1-\alpha).
\]
The Gaussian-width calculation for a Hamming shell gives the ridge bound
\[
R(\delta)
=\sup_{0\le\alpha\le\alpha_0}
\left[
\mu\left(\frac12-\alpha\right)
+4\sqrt D\,\phi(z_\alpha)
\right]. \tag{5}
\]
Since
\[
\frac{d}{d\alpha}\phi(z_\alpha)=z_\alpha,
\]
the interior maximizer, when it lies below \(\alpha_0\), is
\[
\alpha_\delta
=1-\Phi\left(\frac{\mu}{4\sqrt D}\right). \tag{6}
\]
For the sparse rates relevant below this is exponentially tiny in
\(1/\delta\).  Therefore
\[
R(\delta)=\frac12-\delta+O(e^{-c/\delta}). \tag{7}
\]
The exponentially many Hamming descendants of a resonant vector do not
produce the scalar entropy penalty; their common noise makes the ridge cost
essentially just its center.

## 5. Flat-rate perturbation of the bulk

Assume the resonance-excised Paley bulk satisfies the flat layer entropy
\[
s_f(r)=\log2+\frac14\log(1-r^2),
\qquad
0\le r\le r_*:=\frac{\sqrt{15}}4. \tag{8}
\]
The exact Bernoulli moderate-deviation calculation then gives
\[
B(\delta)
=\max_{0\le r\le r_*}
\left[
\frac{\mu r}{2}
+2\sqrt{D\,s_f(r)}
\right]. \tag{9}
\]
The combined corrected two-phase upper bound is
\[
G_{\rm corr}(\delta)=\max\{R(\delta),B(\delta)\}. \tag{10}
\]

Numerical optimization gives
\[
\delta=0.001394039182\ldots,
\qquad
r=0.947800854944\ldots,
\]
and
\[
\boxed{
\inf_\delta G_{\rm corr}(\delta)
=0.498605960818\ldots<\frac12.} \tag{11}
\]
The local correlated-ridge correction changes the earlier scalar answer by
only about \(2\times10^{-12}\).  At the optimizer, the ridge maximum is at
\(\alpha_\delta\approx1.2\times10^{-11}\); the bulk and the damped center
equioscillate.

Ignoring the exponentially tiny term, the algebraic defining equations are
\[
s_f(r)=\frac{r}{4(1+r)},\qquad
\frac{1-2\delta}{2\sqrt{\delta(1-\delta)}}
=\frac{2\sqrt{s_f(r)}}{1-r}. \tag{12}
\]

## 6. Verdict

- The pure orbit inverse theorem (1) is false.
- The first counterexample is the unique affine orbit of multi-harmonic QR
  projection fixed points at \(p=17\).
- The corrected arithmetic-ridge plus ROM-bulk variational model still gives
  the strict constant \(0.498605960818\ldots\).
- This strict bound remains **conditional**, because the required
  resonance-excised flat bulk rate (8), or any sufficiently strong
  substitute, has not been proved for a deterministic Paley sequence.

The next useful theorem is no longer a one-orbit stability result.  It is a
two-phase counting theorem: classify only the cap-near arithmetic ridge,
and prove a flat-rate entropy bound for all remaining multi-harmonic/bulk
fixed points collectively.
