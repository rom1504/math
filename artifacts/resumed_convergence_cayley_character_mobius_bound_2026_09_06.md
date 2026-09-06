# Exact single-character Cayley tests and Möbius spectral control

Date: 2026-09-06. This gives an exact spectral-norm lower test for actual
Cayley signings, with full higher-harmonic and zero-mode accounting. It
does not prove saturation for every Cayley signing.

## 1. The exact character-square-wave formula

Let `G` be a finite abelian group of order `n`, and let `A` be any real
symmetric convolution matrix. Its real even Fourier eigenvalues are
`lambda(a)`, where `a` ranges over the dual group. Write
`Q(A)=max_(x Boolean)|x^T A x|/2`.

For a character `a`, choose a uniform phase `theta in [0,2pi)` and the
Boolean spin

\[
f_{a,\theta}(t)=\operatorname{sign}
\operatorname{Re}(e^{i\theta}\chi_a(t)).
\]

The finitely many tie phases have measure zero and can be assigned
arbitrarily. Define

\[
g(a)=\frac1n\mathbb E_\theta[f_{a,\theta}^T A f_{a,\theta}].
\]

Then, exactly,

\[
\boxed{\quad
g(a)=\frac8{\pi^2}\sum_{m\ge1,\ m\ \mathrm{odd}}
\frac{\lambda(ma)}{m^2}.\quad}                       \tag{1}
\]

In particular `g(0)=lambda(0)` and `|g(a)|<=2Q(A)/n`.

To prove (1), the square wave `sign(cos theta)` has Fourier coefficients
of squared magnitude `4/(pi^2 m^2)` at every nonzero odd integer `m`,
and zero at all other integers. Sampling it along the character `a`
may alias several harmonics into the same group frequency, but averaging
over `theta` eliminates every cross term between distinct integer
harmonics. Parseval in the phase variable yields (1); its series is
absolutely convergent. Real symmetry combines the positive and negative
harmonics using `lambda(-a)=lambda(a)`.

## 2. Möbius inversion and the universal factor three-halves

Let `mu(m)` be the integer Möbius function. Absolute convergence permits
rearrangement of the two series in

\[
\sum_{r\ge1,\ r\ \mathrm{odd}}\frac{\mu(r)}{r^2}g(ra).
\]

The coefficient of `lambda(ka)` is proportional to
`sum_(r|k)mu(r)`, which is zero except at `k=1`. Therefore

\[
\boxed{\quad
\lambda(a)=\frac{\pi^2}{8}
\sum_{m\ge1,\ m\ \mathrm{odd}}
\frac{\mu(m)}{m^2}g(ma).\quad}                       \tag{2}
\]

This remains true when multiplication by `m` annihilates `a`; those
terms are the zero-mode contribution, not terms to discard.

The Euler product gives

\[
\sum_{m\ge1,\ m\ \mathrm{odd}}\frac{|\mu(m)|}{m^2}
=\frac{\zeta(2)}{\zeta(4)(1+2^{-2})}
=\frac{12}{\pi^2}.
\]

Thus

\[
\|\lambda\|_\infty\le\frac32\|g\|_\infty,
\qquad
\boxed{\quad Q(A)\ge\frac n3\|A\|_{\mathrm{op}}.\quad}   \tag{3}
\]

There is no positivity or hollow assumption in (1)--(3), and no
assumption that the kernel entries are signs. The conclusion applies,
in particular, to every real symmetric additive-Cayley hollow signing.

An actual Cayley signing with normalized spectral norm greater than
`3/2+epsilon` consequently has cap coefficient greater than
`1/2+epsilon/3`. A small spectral excess above one is not enough for
this test alone.

## 3. Exact finite-characteristic weights

For `G=F_q` of odd characteristic `p` and `a!=0`, every nontrivial
character has order `p`. One can choose a Boolean function on `F_p`
that is positive on `(p+1)/2` consecutive residues and negative on the
rest. Its mean is `1/p`. A geometric sum computes the squared Fourier
weight on the pair of frequencies `+ell,-ell`, using the unique odd
representatives `ell=1,3,...,p-2`, as

\[
w_{p,\ell}=\frac{2}{p^2\sin^2(\pi\ell/(2p))}.
\]

These weights sum to `1-1/p^2`. They give the exact finite form

\[
\boxed{\quad
g(a)=\frac{\lambda(0)}{p^2}
+\sum_{\substack{1\le\ell\le p-2\\\ell\ \mathrm{odd}}}
w_{p,\ell}\lambda(\ell a).\quad}                   \tag{4}
\]

Up to translations and global negation, all generic phase square waves
on `F_p` give this same profile, so the phase average does not conceal a
loss in the finite-field case.

For example, in characteristic three,

\[
g(a)=\frac{8\lambda(a)+\lambda(0)}9.                \tag{5}
\]

Writing `M=max_(a!=0)|lambda(a)|` and `L=|lambda(0)|`, this gives

\[
\frac{2Q(A)}n\ge
\max\left\{L,\frac{8M-L}9\right\}\ge\frac45M.
\]

Hence `Q(A)>=(2/5)n M`; if `lambda(0)=o(sqrt(n))`, the stronger bound is

\[
Q(A)\ge\frac49 n M-o(n^{3/2}).                     \tag{6}
\]

For general `p`, set `w_p=w_(p,1)` and
`kappa_p=2w_p-1+1/p^2`. First-harmonic dominance alone gives

\[
\frac{2Q(A)}n\ge
\max\{L,\kappa_pM-L/p^2\}
\ge\frac{\kappa_p}{1+1/p^2}M.                     \tag{7}
\]

As `p` grows, `kappa_p->16/pi^2-1`, whose reciprocal is about `1.610`.
The exact Möbius inversion improves this limiting amplitude threshold
to `3/2` by using the tests at all character frequencies jointly.

## 4. Sharpness for this family of tests, not for Boolean caps

The factor `3/2` in `||lambda||infty<=(3/2)||g||infty` is asymptotically
sharp for arbitrary real even Fourier symbols on prime cyclic groups.
It can even be sharp within hollow **weighted** convolution matrices.
This does not claim sharpness of the bound on the full Boolean cap in
(3), and does not construct a sign kernel.

Here is a direct construction. Fix a large integer `R`, and then a prime
`p>4R`. On the quotient of `F_p^*` by overall sign, assign
`g([m])=sign(mu(m))` for odd squarefree `m<=R`. These classes are
distinct. Set `g` to one common value on all remaining classes so that
its average over nonzero frequencies is zero; this common value has
absolute value at most one once `p` is sufficiently large. Set `g(0)=0`.

Define `lambda` by (2). The forward transform (1) returns the prescribed
`g` by the same absolutely convergent divisor cancellation. Moreover

\[
\lambda(1)\ge\frac{\pi^2}{8}
\left(
\sum_{\substack{m\le R\\m\ \mathrm{odd}}}
\frac{\mu(m)^2}{m^2}
-\sum_{\substack{m>R\\m\ \mathrm{odd}}}
\frac{\mu(m)^2}{m^2}
\right).                                           \tag{8}
\]

The right side tends to `3/2` as `R` grows. Also `lambda(0)=0` and
`sum_a lambda(a)=0`: on nonzero frequencies, the forward transform
multiplies the average by `1-1/p^2`. Thus inverse Fourier transformation
gives a real symmetric matrix with zero diagonal. Its kernel values are
general real numbers, not necessarily `+1` and `-1`.

This is a precise reason that a spectral-amplitude argument based only
on the single-character profiles cannot force a `1/2` cap for every
normalized spectral amplitude just above one. Actual sign-kernel
structure must enter an additional lemma.

## 5. A scalable actual-signing obstruction to a naive lower-bulk lemma

Flat real-space entries and a bounded `n^(3/2)` cap do not by themselves
force the lower spectral bulk required by the Ramsey saturation theorem.

Let `q=3^(2e)` and identify `F_q x F_3` additively with the field vector
space of order `n=3q`. Set `b(0)=-1,b(+1)=b(-1)=+1`. Define an actual
hollow even Cayley kernel by

\[
a(t,s)=\chi_q(t)b(s)\quad(t\ne0),\qquad
a(0,\pm1)=1,\qquad a(0,0)=0.                       \tag{9}
\]

Its convolution matrix is exactly
`A_q^Paley tensor B_3 + I_q tensor A_3^complete`.
At frequencies with nonzero `F_q` component, its eigenvalues are

\[
\tau_q\sqrt q\,\chi_q(\alpha)\widehat b(\beta)
+\widehat c(\beta),
\quad
\widehat b=(1,-2,-2),\quad\widehat c=(2,-1,-1).
\]

Here `tau_q in {+1,-1}` is the real normalized quadratic Gauss-sum
sign for the chosen additive-character convention.

Therefore normalized spectral magnitudes approach `1/sqrt(3)` on one
third of all modes and `2/sqrt(3)` on two thirds. Its spectral norm is
`(2/sqrt(3)+o(1))sqrt(n)`, so its cap is `O(n^(3/2))`, but a positive
fraction of eigenvalues stay strictly below `sqrt(n)`.

This is **not** a sub-`1/2` construction: its zero eigenvalue is only two,
and the characteristic-three test (6) already gives cap coefficient at
least `8/(9sqrt(3))`, approximately `0.5132`.

Thus any next lower-bulk theorem for actual Cayley minimizers must use a
genuine deficit below `1/2`, not merely bounded cap or sign-valued
entries. The quantitative statement below makes the necessary spectral
defect uniform, rather than merely a positive limsup on a subsequence.

## 6. A uniform necessary spectral-defect density below one-half

Fix `0<delta<1/2`. There are finite constants `L_delta,Q_delta` such
that for every odd prime power `q>=Q_delta`, an actual hollow even
additive-Cayley signing satisfying

    Q(A)/q^(3/2) <= 1/2-delta

must obey BOTH

\[
\frac{\|A\|_{\rm op}}{\sqrt q}\le\frac32-3\delta,
\qquad
\#\{a\ne0:|\lambda(a)|<(1-\delta/2)\sqrt q\}
\ge\frac{q-1}{L_\delta}.                           \tag{10}
\]

The first statement is (3). For the second use the finite Ramsey-host
module in `resumed_convergence_fourier_involution_ramsey_saturation_2026_09_06.md`,
with `epsilon=delta/16`, and let `H` be its fixed nonzero frequency host
of size at most `L_delta`. Its guarantee, for every even two-coloring
of `bH`, is a Boolean `f` and a mean-zero polynomial `P` such that

    ||f-P||_2 <= epsilon,   ||P||_2 >= 1-epsilon,
    supp Fourier(P) lies in one color of bH.

The norms here are probability-normalized. The host depends on the
characteristic regime as in that proof, but its size and the threshold
can be bounded uniformly using the finite large-characteristic host
and the finitely many bounded-characteristic hosts.

Put `eta=delta/2` and call the modes below `(1-eta)sqrt(q)` bad.
If any nonzero dilation `b` avoids all bad modes on `bH`, color these
frequencies by the sign of their eigenvalue. For its polynomial choose
the monochromatic orientation `sigma in {+1,-1}`. With
`U=A/sqrt(q)` and `C=||U||op<=3/2`,

\[
\begin{aligned}
\sigma\langle f,Uf\rangle
&\ge(1-\eta)(1-\epsilon)^2
       -C\epsilon(2+\epsilon)\\
&\ge1-\frac{13}{16}\delta
       +\frac{31\delta^2-\delta^3}{512}
>1-\delta>1-2\delta.                              \tag{11}
\end{aligned}
\]

The transfer error follows directly from Cauchy--Schwarz and
`||P||_2<=1+epsilon`. Equation (11) contradicts the assumed cap.
Hence EVERY dilation `b` hits the bad set. Each fixed host element
ranges bijectively over the nonzero field under dilation, so the union
bound gives `q-1<=|H| |Bad|`, proving (10).

This is an actual-signing necessary condition with a uniform positive
defect density depending on `delta`. It does not contradict Parseval:
the identity `q^-1 sum_a lambda(a)^2/q=1-1/q` can coexist with a fixed
fraction of smaller modes and spectral norm below `3/2`. No result in
this artifact excludes that remaining regime or proves all-Cayley
saturation.

One geometric consequence of (10) is worth making explicit. With
`eta=delta/2`, a deficient Cayley signing is separated from EVERY scaled
real symmetric involution, not just from Fourier-diagonal ones:

\[
\begin{aligned}
\inf_{J=J^T,\ J^2=I}\|A-\sqrt q J\|_F^2
&=\sum_a(|\lambda(a)|-\sqrt q)^2
\ge\frac{\eta^2 q(q-1)}{L_\delta},\\
\frac1q\operatorname{Tr}\left[(A^2/q-I)^2\right]
&\ge\frac{q-1}{q L_\delta}(2\eta-\eta^2)^2.
\tag{12}
\end{aligned}
\]

For the first equality, expanding the square reduces the optimization
to maximizing `Tr(AJ)`, whose value is `sum_a|lambda(a)|`; the spectral
sign of `A`, with arbitrary signs on its kernel, attains it. Each bad
mode contributes the stated amount to both sums. Consequently the
sub-half Cayley regime requires a genuinely nonvanishing spectral
variance and distance from all conference/Hadamard-type involutions.
This statement remains restricted to a Cayley TARGET; it does not show
that every non-Cayley signing close to an arbitrary involution has a
one-half cap.

## 7. An unconditional four-ninths theorem in characteristic three

For EVERY hollow even additive-Cayley signing on `F_(3^r)`, with
`q=3^r>=81`,

\[
\boxed{\qquad Q(A)\ge\frac49 q^{3/2}.\qquad}       \tag{13}
\]

There is no spectral-flatness or small-row-sum hypothesis. In fact the
complete square-wave accounting in (5), together with the actual
sign-kernel identities

    sum_a lambda(a)=0,   sum_a lambda(a)^2=q(q-1),

gives the exact formula

\[
\sum_{a\ne0}g(a)^2
=\frac{64q(q-1)+(q-81)\lambda(0)^2}{81}.          \tag{14}
\]

For `q>=81`, the right side is at least `64q(q-1)/81`.
Thus some nonzero character has `|g(a)|>=8sqrt(q)/9`, proving (13).
This strengthens the general universal lower bound within a genuine
nonflat actual-signing class. It does not extend unchanged to other
odd characteristics, whose single-character profiles have additional
independent harmonic eigenvalues.

The four-ninths constant is quantitatively improved, without adding a
flatness assumption, in
`resumed_convergence_characteristic_three_incidence_gain_2026_09_06.md`:
the asymptotic lower bound there is
`2080/(9 sqrt(269441))`, approximately `0.4452346798594428`. It uses the actual
compatibility of projective-line color patterns inside three-dimensional
subspaces, not merely the spectral second moments.
