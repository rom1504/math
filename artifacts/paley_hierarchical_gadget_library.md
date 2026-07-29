# Hierarchical Paley gadget library: exact construction, counts, and entropy

## Status

The finite-cyclic resonance gadget admits a much larger exact closure than
the original product construction.  At depth \(r\), **any** Boolean function
of the even base-\(\ell\) digits has Fourier transform zero at every
frequency of odd \(\ell\)-adic valuation.

For nearly balanced active-digit functions this gives exact multiplier
energy

\[
1-\ell^{-2r-2}. \tag{1}
\]

The resulting family is large:

\[
\log |\mathcal G_{\ell,r}|=\Theta(\ell^{r+1}),
\]

but its period is \(q=\ell^{2r+1}\), so

\[
\log|\mathcal G_{\ell,r}|=\Theta(\sqrt{\ell q}).
\]

Consequently the entire proven macroscopic step-gadget closure remains
\(\exp(o(p))\) when sampled on \(\mathbb F_p\) with period \(q=o(p)\).
No positive-entropy family of orbit-separated cap-near templates arises
from this mechanism.

## 1. The full active-digit construction

Fix an odd prime \(\ell\), an integer \(r\ge0\), and put

\[
q=\ell^{2r+1}.
\]

Write

\[
j=\sum_{d=0}^{2r}j_d\ell^d,\qquad 0\le j_d<\ell.
\]

Let

\[
F:\{0,\ldots,\ell-1\}^{r+1}\longrightarrow\{\pm1\}
\]

be any Boolean function, and define

\[
\boxed{\displaystyle
v_F(j)=F(j_0,j_2,\ldots,j_{2r}).} \tag{2}
\]

Thus \(v_F\) is independent of all odd-position digits.  The earlier
product construction

\[
v(j)=\prod_{e=0}^r u_e(j_{2e})
\]

is only a special case.

Use the normalized DFT

\[
\widehat v_F(k)
=\frac1q\sum_{j\bmod q}v_F(j)e^{-2\pi i kj/q}.
\]

## 2. Exact support-annihilation theorem

Suppose \(k\ne0\bmod q\) and \(s=v_\ell(k)\) is odd.  The digit

\[
d=2r-s
\]

is odd.  In the exponential, its one-digit factor is

\[
\sum_{a=0}^{\ell-1}
\exp\left(
-\,\frac{2\pi i(k/\ell^s)a}{\ell}
\right)=0, \tag{3}
\]

because \(k/\ell^s\) is a unit modulo \(\ell\).  The function \(v_F\) is
independent of this odd digit, so the full DFT factorizes through (3).
Therefore

\[
\boxed{\displaystyle
\widehat v_F(k)=0
\quad\text{whenever }v_\ell(k)\text{ is odd}.} \tag{4}
\]

This proof uses no product structure and no condition on \(F\).

## 3. Exact energy

There are \(\ell^r\) choices of the inactive odd digits, so

\[
\sum_{j\bmod q}v_F(j)
=\ell^r\sum_zF(z). \tag{5}
\]

The number

\[
M=\ell^{r+1}
\]

of active-digit states is odd.  Choose \(F\) with

\[
\sum_zF(z)=\pm1. \tag{6}
\]

Then

\[
|\widehat v_F(0)|^2
=\left(\frac{\ell^r}{\ell^{2r+1}}\right)^2
=\ell^{-2r-2}. \tag{7}
\]

For the signed Fourier multiplier

\[
\eta_\ell(k)=(-1)^{v_\ell(k)},
\]

every nonzero frequency in the support of \(\widehat v_F\) has multiplier
\(+1\), by (4).  Parseval therefore gives the exact identity

\[
\boxed{\displaystyle
\sum_{k\ne0}\eta_\ell(k)|\widehat v_F(k)|^2
=1-\ell^{-2r-2}.} \tag{8}
\]

This is the exact finite-cyclic energy before the standard
Dirichlet-refinement and circle-sampling transfer to Paley primes.

## 4. The length-27 certificate

Take \(\ell=3\), \(r=1\), and

\[
u=(1,1,-1),\qquad
F(j_0,j_2)=u_{j_0}u_{j_2}.
\]

In the ordering \(j=j_0+3j_1+9j_2\), (2) is

```text
+-++-++-++-++-++---+--+--+
```

Its sum is \(3\), so its normalized DC mass is \(1/81\).
Its DFT vanishes at precisely the nonzero frequencies with
\(v_3(k)=1\).  Hence

\[
\boxed{\displaystyle
E_{\eta_3}(v)=1-\frac1{81}
=\frac{80}{81}
=0.987654320987\ldots .} \tag{9}
\]

This lies above

\[
r_*=\frac{\sqrt{15}}4=0.968245836552\ldots .
\]

At depth \(r=2\), the same mechanism has period \(3^5=243\), DC mass
\(3^{-6}=1/729\), and energy

\[
\frac{728}{729}=0.998628257888\ldots .
\]

## 5. Exact template count at one prime

Modulo global negation, impose \(\sum F=1\).  The number of active-digit
functions is

\[
\boxed{\displaystyle
N_{\ell,r}
=\binom{\ell^{r+1}}{(\ell^{r+1}+1)/2}.} \tag{10}
\]

Stirling's formula gives

\[
\log N_{\ell,r}
=(\log2)\ell^{r+1}
-\frac12\log(\ell^{r+1})+O(1). \tag{11}
\]

Since \(q=\ell^{2r+1}\),

\[
\ell^{r+1}=\sqrt{\ell q},
\]

and therefore

\[
\boxed{\displaystyle
\log N_{\ell,r}=\Theta(\sqrt q)
\quad(\ell\text{ fixed}).} \tag{12}
\]

For the length-27 family, \(M=9\), so there are

\[
N_{3,1}=\binom95=126 \tag{13}
\]

nearly balanced active-digit shapes modulo global negation.  The original
separable product family contains only \(3^2=9\) of them.

## 6. Several bad primes

Let \(B\) be a finite set of odd primes, use depths \(r_\ell\), and put

\[
Q=\prod_{\ell\in B}\ell^{2r_\ell+1},
\qquad
L_B=\prod_{\ell\in B}\ell.
\]

The active-state space has size

\[
M_B=\prod_{\ell\in B}\ell^{r_\ell+1}
=\sqrt{Q L_B}. \tag{14}
\]

An arbitrary nearly balanced Boolean function on this product active space
again annihilates every local odd-valuation frequency, so the number of
shapes is

\[
N_B=\binom{M_B}{(M_B+1)/2},
\qquad
\log N_B=(\log2)M_B+O(\log M_B). \tag{15}
\]

The retained Fourier mass and energy depend on the treatment of local zero
coordinates exactly as in `paley_resonance_gadget.md`; the counting statement
(15) is independent of that boundary loss.

Because \(L_B\le Q\),

\[
M_B\le Q. \tag{16}
\]

Thus every macroscopic step-gadget regime with

\[
Q=o(p) \tag{17}
\]

has

\[
\log N_B=o(p). \tag{18}
\]

Condition (17) is also the natural regime in which every step interval
contains \(p/Q\to\infty\) sample points and the Riemann-sum Paley transfer is
uniformly meaningful.

For a fixed target energy gap below the cap, the Fourier cutoff, bad-prime
set, and depths may all be fixed before \(p\to\infty\).  Then \(N_B\) is an
enormous but \(p\)-independent constant.

## 7. Affine orbit size

Let \(f_p\) be a sampled nonconstant gadget vector.  The signed affine group

\[
\{\pm1\}\times\operatorname{AGL}(1,p)
\]

has size

\[
2p(p-1).
\]

Therefore every gadget shape has affine signed orbit

\[
\boxed{\displaystyle
|\mathcal O_p(f)|
=\frac{2p(p-1)}
{|\operatorname{Stab}^{\pm}_{\rm aff}(f_p)|}
\le2p(p-1).} \tag{19}
\]

The translation part of the signed stabilizer is always trivial for a
nonconstant vector:

- a nonzero ordinary translation generates the additive group and would
  force the vector to be constant;
- a signed anti-translation is impossible for odd \(p\), since iterating it
  \(p\) times would give \(x=-x\).

A nontrivial multiplicative affine stabilizer is possible in principle for a
specially symmetric active-digit function.  For the explicit length-27 word,
direct exact checks at \(p=109,163,271,433\) find trivial full signed
stabilizer, and hence the maximal orbit size \(2p(p-1)\).

The upper bound in (19), not stabilizer triviality, is what matters for the
entropy and perturbation arguments.

## 8. Does hierarchical closure acquire positive entropy?

For all gadget shapes with periods at most \(Q_p=o(p)\), (15)--(18) imply

\[
\log|\mathcal T_p|
\le o(p)+O(\log p), \tag{20}
\]

even after multiplying by every affine signed orbit.  Hence

\[
\boxed{\displaystyle
|\mathcal T_p|=\exp(o(p)).} \tag{21}
\]

To obtain positive entropy from active-digit freedom alone, one would need

\[
M_B=\Theta(p),
\]

which by (14) requires \(Q L_B=\Theta(p^2)\).  In particular the period can
no longer be \(o(p)\) unless the bad-prime product has an extreme,
currently unsupported growth pattern.  This lies outside the proven
step-function transfer mechanism.

So the current verdict is:

- **proved subexponential:** every fixed-accuracy gadget library and every
  growing library with period \(Q_p=o(p)\);
- **not constructed:** any positive-entropy family of mutually
  orbit-separated cap-near arithmetic templates;
- **unresolved:** genuinely discrete period-\(\Theta(p)\) gadgets that do
  not arise as sampled macroscopic step functions.

## 9. Perturbation after quotienting the full gadget library

Let \(\mathcal T_p\) be any template library satisfying (21).  For a fixed
edge-flip rate \(\delta\), the centered perturbation at one template has a
speed-\(p\) tail at scale \(p^{3/2}\).  Therefore

\[
\max_{x\in\mathcal T_p}|Z_x|=o_{\Pr}(p^{3/2}). \tag{22}
\]

All template centers are simultaneously damped by \(1-2\delta\).
Likewise:

- unioning fixed-radius shell estimates over \(\mathcal T_p\) adds only
  \(o(p)\) to the logarithmic complexity;
- unioning two-template patchwork bounds over
  \(\mathcal T_p^2\) still adds only \(o(p)\);
- the exact patchwork width and its
  \(\sqrt{2/(\pi+8)}\) diagnostic constant are unchanged at leading order.

Thus arithmetic template multiplicity is not itself a perturbative
obstruction.

The remaining missing theorem is a **coverage** statement:

> Every Paley vector above a chosen threshold is either within a controlled
> Hamming tube of one gadget orbit in \(\mathcal T_p\), or lies in a
> coordinatewise patchwork hull whose centered random-cut width is governed
> by the subcube comparison.

No such inverse theorem is proved.  But the counting audit shows that the
known hierarchical arithmetic closure is small enough for this strategy to
be logically viable.

## 10. Bottom line

The length-27 resonance is not an isolated anomaly.  It is the first member
of an exact hierarchical family with energy \(1-\ell^{-2r-2}\) and
\(\exp(\Theta(\sqrt q))\) shapes at period \(q\).  This forces the Paley
inverse theorem to use an arithmetic template library richer than square
waves.

At the same time, the hierarchy does **not** yet create positive entropy at
the ambient order \(p\).  Quotienting all currently proved gadget orbits
remains compatible with sparse perturbation and with a strict upper bound
below \(1/2\).

