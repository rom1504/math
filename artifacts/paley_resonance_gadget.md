# A finite-cyclic resonance obstruction for Paley signings

## Statement

For a prime \(p\equiv 1\pmod 4\), let

\[
(A_p)_{ij}=\chi_p(i-j),\qquad i,j\in\mathbb F_p,
\]

where \(\chi_p\) is the quadratic character and \(\chi_p(0)=0\).  Put

\[
\mathcal E_p(x)=\frac{x^\top A_px}{p^{3/2}},
\qquad x\in\{\pm1\}^{p}.
\]

The following is the useful obstruction.

**Theorem (resonance in every finite character cylinder).**  
Fix the values \(\epsilon_\ell\in\{\pm1\}\) of \(\chi_p(\ell)\) for
finitely many rational primes \(\ell\), compatibly with
\(p\equiv1\pmod4\).  For every \(\delta>0\), there is a refinement by
finitely many further character conditions, hence a fixed arithmetic
progression containing infinitely many primes, such that for every
sufficiently large prime \(p\) in that progression there is
\(x^{(p)}\in\{\pm1\}^{p}\) satisfying

\[
\mathcal E_p(x^{(p)})
\ge 1-\delta.
\]

Equivalently, for the half-energy
\(H_p(x)=\frac12x^\top A_px\),

\[
\frac{H_p(x^{(p)})}{p^{3/2}}
\ge\frac12-\frac{\delta}{2}.
\]

Consequently, no condition involving only finitely many Legendre
symbols—and, more generally, no fixed admissible congruence class—can
have any eventual uniform gap below the Paley spectral cap \(1/2\).
In particular, it cannot be an eventually nonresonant Paley
subsequence at the proposed random-orthogonal value
\(\sqrt{15}/8=0.484122918276\ldots\).  Any such subsequence must use an
increasing, \(p\)-dependent amount of arithmetic information.

Equivalently, for every fixed admissible progression \(\mathcal A\) of
primes \(p\equiv1\pmod4\),

\[
\boxed{\displaystyle
\limsup_{\substack{p\to\infty\\p\in\mathcal A}}
\max_{x\in\{\pm1\}^p}
\frac{|x^\top A_px|}{2p^{3/2}}
=\frac12.}
\]

The theorem does **not** show that every prime is resonant, nor that the
Paley maxima fail to have a limit.  It rules out finite-congruence
criteria for proving a complementary upper bound.

## 1. Paley Fourier normalization

For \(x:\mathbb F_p\to\mathbb R\), use

\[
\widehat x_p(m)=\frac1p\sum_{j\in\mathbb F_p}
x_j e^{-2\pi i mj/p}.
\]

Then Parseval says

\[
\sum_{m\in\mathbb F_p}|\widehat x_p(m)|^2
=\frac1p\sum_j|x_j|^2.
\]

For Boolean \(x\), the right side is \(1\).  The quadratic Gauss-sum
identity, with a harmless global sign fixed by the Fourier convention,
gives

\[
\boxed{\displaystyle
\mathcal E_p(x)
=\sum_{m\in\mathbb F_p^\times}
\chi_p(m)|\widehat x_p(m)|^2.}
\tag{1}
\]

Thus \(\mathcal E_p(x)\) is the difference between the normalized
Fourier mass on quadratic-residue and nonresidue frequencies.

## 2. One odd-prime gadget at arbitrary depth

Fix an odd prime \(\ell\), an integer \(r\ge1\), put
\(q=\ell^{2r+1}\), and choose
\(u_0,\ldots,u_{\ell-1}\in\{\pm1\}\) with

\[
\sum_{a=0}^{\ell-1}u_a=1.
\]

Write the base-\(\ell\) expansion

\[
j=\sum_{d=0}^{2r}j_d\ell^d,
\qquad 0\le j_d<\ell,
\]

and define

\[
v^{(\ell,r)}_j=\prod_{\substack{0\le d\le2r\\d\ {\rm even}}}u_{j_d}.
\tag{2}
\]

Its sum is

\[
\sum_jv^{(\ell,r)}_j
=\ell^r.
\tag{3}
\]

Use the normalized finite DFT

\[
\widehat v^{(\ell,r)}(k)
=\frac1{\ell^{2r+1}}\sum_{j\bmod\ell^{2r+1}}
v^{(\ell,r)}_j e^{-2\pi i kj/\ell^{2r+1}}.
\]

The digit expansion makes this DFT a product of \(2r+1\) one-digit
sums.  Suppose \(k\ne0\) modulo \(\ell^{2r+1}\) and
\(s=v_\ell(k)\) is odd.  The digit

\[
d=2r-s
\]

is odd, so the corresponding digit in (2) is unweighted.  Its DFT
factor is a complete nontrivial \(\ell\)-th-root sum and vanishes.
Therefore

\[
\boxed{\widehat v^{(\ell,r)}(k)=0
\quad\text{whenever }v_\ell(k)\text{ is odd}.}
\tag{4}
\]

Hence the entire non-DC Fourier support has even \(\ell\)-adic
valuation.  By (3), its normalized DC mass is

\[
|\widehat v^{(\ell,r)}(0)|^2
=\left(\frac{\ell^r}{\ell^{2r+1}}\right)^2
=\ell^{-2r-2}.
\tag{5}
\]

Thus exactly \(1-\ell^{-2r-2}\) of its Fourier mass is non-DC and has
even \(\ell\)-adic valuation.  This loss tends to zero as the depth
\(r\) grows.

## 3. The prime \(2\)

If the finite prescription includes \(\chi_p(2)=-1\), the length-two
word

\[
v^{(2)}=(1,-1)
\tag{6}
\]

has zero DC and full Fourier mass at the odd frequency.  Thus this
gadget also has only even \(2\)-adic valuation on its support and
incurs no loss.  Higher-depth binary versions follow from the same
digit construction with \(u=(1,-1)\), but are unnecessary.

## 4. CRT tensor and its Fourier mass

Let \(B\) be the finite set of primes prescribed to have character
\(-1\).  For each odd \(\ell\in B\), choose a depth \(r_\ell\), use the
gadget above, and put \(q_\ell=\ell^{2r_\ell+1}\).  If \(2\in B\), use
the length-two gadget and put \(q_2=2\).  Let

\[
Q=\prod_{\ell\in B}q_\ell.
\]

Through the Chinese remainder isomorphism

\[
\mathbb Z/Q\mathbb Z
\cong\prod_{\ell\in B}\mathbb Z/q_\ell\mathbb Z,
\]

define the Boolean tensor word

\[
V(j)=
\prod_{\substack{\ell\in B\\\ell\ {\rm odd}}}
v^{(\ell,r_\ell)}(j_\ell)
\times
\begin{cases}
v^{(2)}(j_2),&2\in B,\\
1,&2\notin B.
\end{cases}
\tag{7}
\]

If \(B=\varnothing\), use instead any balanced square wave on the
circle (for example \(\operatorname{sgn}\cos(2\pi t)\)).  Its mean is
zero and, after the later finite character refinement makes all
relevant prime symbols positive, its entire controlled Fourier mass
has multiplier \(+1\).  Thus the empty-\(B\) case has no loss and is
strictly easier.  In the rest of the tensor argument assume
\(B\ne\varnothing\).

The finite DFT factors under the dual CRT isomorphism.  Multiplication
of a local frequency by the CRT unit does not alter its
\(\ell\)-adic valuation.  Let \(\Gamma\) be the part of the Fourier
support for which no odd local coordinate is zero.  Equations (4) and
(5), Parseval, and tensorization give

\[
\sum_{k\in\Gamma}|\widehat V(k)|^2
=G_B,
\qquad
G_B:=\prod_{\substack{\ell\in B\\ \ell\ {\rm odd}}}
(1-\ell^{-2r_\ell-2}).
\tag{8}
\]

At every \(k\in\Gamma\), each bad-prime valuation is even.
The \(2\)-coordinate, when present, is odd by (6), hence has valuation
zero.  Consequently

\[
\eta_B(k):=\prod_{\ell\in B}(-1)^{v_\ell(k)}=+1
\qquad(k\in\Gamma).
\tag{9}
\]

## 5. From the finite word to a circle step function

Define \(f_B:[0,1)\to\{\pm1\}\) by

\[
f_B(t)=V(j)
\quad\text{for}\quad
\frac jQ\le t<\frac{j+1}{Q}.
\tag{10}
\]

Let

\[
c_m=\int_0^1f_B(t)e^{-2\pi i mt}\,dt.
\]

For \(m\ne0\),

\[
c_m
=\frac{1-e^{-2\pi im/Q}}{2\pi i m}
\sum_{j=0}^{Q-1}V(j)e^{-2\pi imj/Q}.
\tag{11}
\]

For every residue \(k\bmod Q\), the classical partial-fraction identity

\[
\sum_{r\in\mathbb Z}
\frac{\sin^2(\pi k/Q)}
{\pi^2(k/Q+r)^2}=1
\]

implies the exact aliasing identity

\[
\boxed{\displaystyle
\sum_{r\in\mathbb Z}|c_{k+rQ}|^2
=|\widehat V(k)|^2.}
\tag{12}
\]

When \(k\in\Gamma\), every integer \(m\equiv k\pmod Q\) has the same
even valuation (strictly below \(2r_\ell+1\)) at each odd
\(\ell\in B\), and is odd if \(2\in B\).  Thus
\(\eta_B(m)=+1\) on all aliases of \(\Gamma\).
Equations (8) and (12) show that at least \(G_B\) of the Fourier mass
of \(f_B\) lies on \(\eta_B=+1\).  Treating every remaining nonzero
frequency adversarially gives

\[
\sum_{m\ne0}\eta_B(m)|c_m|^2
\ge 2G_B-1.
\tag{13}
\]

Because \(B\) is finite, the depths \(r_\ell\) can be chosen so large
that \(G_B\) is arbitrarily close to \(1\).  In particular, for any
\(\varepsilon>0\), one can arrange

\[
G_B>1-\varepsilon
\tag{14}
\]

and hence

\[
\sum_{m\ne0}\eta_B(m)|c_m|^2
\ge1-2\varepsilon.
\tag{15}
\]

As a depth-one uniform numerical benchmark, taking \(r_\ell=1\) for
every odd bad prime gives

\[
G_B\ge
\prod_{\ell\ {\rm odd\ prime}}(1-\ell^{-4})
=\frac{96}{\pi^4}
\]

and therefore the earlier explicit bound
\(192/\pi^4-1=0.971068592899\ldots\).  Arbitrary depth is what upgrades
this to the full cap.

## 6. Dirichlet refinement and sampling on \(\mathbb F_p\)

Start from any compatible finite prescription
\(\chi_p(\ell)=\epsilon_\ell\).  Its negative primes are the set \(B\)
used above.  Fix \(\delta>0\).  Since \(f_B\in L^2([0,1])\), choose a
finite Fourier cutoff \(L\) so that

\[
\sum_{|m|>L}|c_m|^2<\delta.
\tag{16}
\]

For every rational prime \(r\le L\) not already prescribed, impose

\[
\chi_p(r)=+1.
\tag{17}
\]

These conditions are simultaneously realizable.  For odd \(r\),
quadratic reciprocity and \(p\equiv1\pmod4\) give

\[
\left(\frac rp\right)=\left(\frac pr\right),
\]

so choose a residue or nonresidue class for \(p\bmod r\) as required.
For \(r=2\), choose \(p\bmod8\) according to the supplementary law.
The Chinese remainder theorem combines the choices into an admissible
class \(p\equiv a\pmod M\), and Dirichlet's theorem supplies infinitely
many primes in that class.  The prime number theorem in arithmetic
progressions also shows that this fixed refined class has positive
density among all primes and that its successive primes have ratio
tending to \(1\).

For these primes and every \(0<|m|\le L\), complete multiplicativity
gives

\[
\chi_p(m)=\eta_B(m).
\tag{18}
\]

Now sample the step function:

\[
x^{(p)}_j=f_B(j/p),\qquad 0\le j<p.
\tag{19}
\]

For each fixed integer \(m\), Riemann-sum convergence gives

\[
\widehat x_p(m)\longrightarrow c_m.
\tag{20}
\]

Using (1), split the Paley energy into \(|m|\le L\) and its complement.
The complementary signed mass is bounded below by minus its total
Fourier mass.  Parseval, (16), and (20) therefore give

\[
\liminf_{\substack{p\to\infty\\p\equiv a\pmod M}}
\mathcal E_p(x^{(p)})
\ge
\sum_{m\ne0}\eta_B(m)|c_m|^2-2\delta.
\tag{21}
\]

Because the cutoff in (16) may be chosen with \(\delta/2\) in place of
\(\delta\), equations (15) and (21) prove the theorem with any desired
final error.

If one starts with a fixed admissible progression rather than only a
finite list of Legendre symbols, retain its existing residue conditions
and add (17) only at new coprime prime moduli.  This produces a fixed
refined progression inside the original one.  The positive-density
claim is for that refined progression; its density may be extremely
small and depends on the cutoff.

## 7. Independent check for \(B=\{3\}\)

Take \(\ell=3\), \(u=(1,1,-1)\), and

\[
v_{b+3c+9a}=u_au_b.
\]

The resulting length-\(27\) word is

```text
++-++-++-++-++-++---+--+--+
```

It has sum \(3\), hence DC mass \(1/81\).  A direct FFT gives zero
(up to floating-point roundoff) at exactly the six frequencies

\[
3,6,12,15,21,24,
\]

which are those with \(3\)-adic valuation \(1\) modulo \(27\).  All
non-DC Fourier mass therefore has even \(3\)-adic valuation, and the
signed multiplier energy is exactly

\[
1-\frac1{81}=\frac{80}{81}
=0.987654320987\ldots .
\]

This numerical check independently matches the support proof and the
DC normalization.

At depth \(r=2\), the same construction has length \(3^5=243\).  A
separate FFT check gives zero at every frequency of \(3\)-adic
valuation \(1\) or \(3\), DC mass \(1/729\), and signed energy
\(728/729=0.998628257888\ldots\), confirming the arbitrary-depth
formula.

## 8. Relation to the simpler square-wave resonance

For \(p=4r+1\), the balanced interval square wave
\(s_j=\operatorname{sgn}\cos(2\pi j/p)\) has

\[
\widehat s(m)
=(-1)^{(m-1)/2}\csc\!\frac{\pi m}{2p}
\quad(m\ {\rm odd}),
\]

with the paired even formula giving only \(O(p^{-1})\) normalized mass.
Consequently

\[
\mathcal E_p(s)
=\frac8{\pi^2}
\sum_{\substack{h\ge1\\h\ {\rm odd}}}
\frac{\chi_p(h)}{h^2}
+O(p^{-1}).
\tag{22}
\]

If all odd \(h\le25\) are residues, then even an adversarial tail gives

\[
\mathcal E_p(s)
\ge
\frac{16}{\pi^2}
\sum_{\substack{h\le25\\h\ {\rm odd}}}\frac1{h^2}
-1-O(p^{-1})
=0.968839592155\ldots-O(p^{-1})
>\frac{\sqrt{15}}4.
\]

It is enough to require
\(3,5,7,11,13,17,19,23\) to be residues.  These conditions have
relative density \(2^{-8}\) among primes \(p\equiv1\pmod4\).
Thus even the elementary square wave already disproves a
“density-one nonresonant primes” claim.  The finite-cyclic gadget is
stronger: it adapts to arbitrary finitely many prescribed nonresidues
and shows that every finite character cylinder contains its own
resonant refinement.

In fact, (22) gives a clean limiting-distribution statement.  Let
\((X_\ell)_{\ell\ {\rm odd\ prime}}\) be independent Rademacher signs
and extend them completely multiplicatively to odd integers.  As
\(p\to\infty\) through primes \(p\equiv1\pmod4\), the square-wave
statistics converge in distribution (with primes counted by natural
density) to

\[
T=\frac8{\pi^2}\sum_{\substack{h\ge1\\h\ {\rm odd}}}
\frac{X(h)}{h^2}
=\frac8{\pi^2}\prod_{\ell\ {\rm odd\ prime}}
\left(1-\frac{X_\ell}{\ell^2}\right)^{-1}.
\tag{23}
\]

To prove this, first truncate the absolutely convergent series.
Quadratic reciprocity, CRT, and the prime number theorem in arithmetic
progressions make the finitely many symbols at the primes in the
truncation independent uniform signs in the limit.  The omitted tail
is bounded uniformly by \(O(1/L)\), so the truncation may then tend to
infinity.

The law in (23) is nondegenerate.  If \(X_3=-1\), all remaining Euler
factors are maximized by taking the other signs positive, and the
result is at most \(0.8\).  Thus

\[
\mathbb P(T\le0.8)\ge\frac12.
\]

On the other hand, if
\(X_3=X_5=X_7=X_{11}=X_{13}=X_{17}=X_{19}=X_{23}=+1\),
the same worst-tail estimate as above gives

\[
T\ge0.968839592155\ldots>\frac{\sqrt{15}}4,
\]

so

\[
\mathbb P\!\left(T>\frac{\sqrt{15}}4\right)\ge2^{-8}.
\]

Therefore the square-wave Paley statistic itself does not converge in
density to any constant.  This still does not, by itself, prove that
the **maximum over all Boolean vectors** has no limit.
