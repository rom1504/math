# Quadratic Signing Limit — Research Ledger

Last updated: 2026-07-25 (America/Los_Angeles)

## Purpose

This is the durable checkpoint for the problem

\[
M_n=\min_{a_{ij}\in\{\pm1\}}\max_{x_i\in\{\pm1\}}
\left|\sum_{1\le i<j\le n}a_{ij}x_ix_j\right|,
\qquad
\text{determine whether }\lim_{n\to\infty}\frac{M_n}{n^{3/2}}
\text{ exists.}
\]

It separates:

- **Verified:** proof has been reconstructed and checked in the current work.
- **Pending audit:** plausible result from an earlier research wave, but its full proof has not yet been reconstructed.
- **Numerical:** computational evidence only.
- **Falsified:** a precise gap or counterexample has been found.
- **Open target:** a lemma that would materially advance or settle the problem.

## Notation

Let \(A\) be the symmetric zero-diagonal matrix with off-diagonal entries
\(a_{ij}\), and define

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top A x|.
\]

Then

\[
x^\top A x=2H_A(x),\qquad Q(A)=2\max_x|H_A(x)|,
\qquad M_n=\frac12\min_A Q(A).
\]

The conjecturally simplest outcome is

\[
\frac{M_n}{n^{3/2}}\longrightarrow \frac12,
\]

equivalently \(Q(A)\ge(1-o(1))n^{3/2}\) for every signing \(A\).
This is not proved.

---

## 1. Verified results

### 1.1 Current rigorous asymptotic interval

\[
\boxed{
0.336493364431\ldots
\le
\liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\le
\limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le
\frac12
}
\]

The former verified lower constant \(1/\pi=0.318309886\ldots\) is
superseded by the field-plus-spin theorem in §1.7.

### 1.2 Finite-\(n\) Gaussian-rounding lower bound

For every \(n\ge2\),

\[
\boxed{
M_n\ge
\frac{n(n-1)}{\pi}
\arcsin\!\frac1{\sqrt{n-1}}
}
\]

and hence

\[
M_n\ge
\frac1\pi n^{3/2}
-\frac1{3\pi}\sqrt n
+O(n^{-1/2}).
\]

Proof checkpoint:

1. Fix \(A\), put \(m=n-1\), and let
   \[
   P=\max_x H_A(x),\qquad Q=-\min_x H_A(x).
   \]
2. For \(s=\pm1\), take \(g\sim N(0,I)\) and set
   \[
   X^{(s)}=\operatorname{sgn}\bigl((sA+\sqrt m\,I)g\bigr).
   \]
3. The pre-sign coordinates have variance \(2m\). For \(i\ne j\), after
   multiplication by \(a_{ij}\), their correlation is
   \[
   s\,a+b_{ij},
   \qquad
   a=\frac1{\sqrt m},
   \qquad
   b_{ij}=\frac{a_{ij}(A^2)_{ij}}{2m}.
   \]
4. The Gaussian arcsine identity gives
   \[
   P+Q\ge
   \frac2\pi\sum_{i<j}
   \left[
   \arcsin(a+b_{ij})+\arcsin(a-b_{ij})
   \right].
   \]
5. On the admissible domain,
   \[
   \frac{\arcsin(a+b)+\arcsin(a-b)}2\ge\arcsin a.
   \]
6. Since \(\max(P,Q)\ge(P+Q)/2\), the displayed finite bound follows.

### 1.3 Conference-matrix upper bound

If \(C\) is a symmetric conference matrix of order \(N\), then

\[
C^2=(N-1)I
\]

and therefore

\[
\max_x\left|\sum_{i<j}c_{ij}x_ix_j\right|
=\frac12\max_x|x^\top Cx|
\le\frac12N\sqrt{N-1}.
\]

Symmetric Paley conference matrices exist at orders \(N=q+1\) for prime
powers \(q\equiv1\pmod4\). Primes in this progression can be chosen with
\(q=n+o(n)\), and principal submatrices handle intermediate orders. Thus

\[
\limsup_{n\to\infty}\frac{M_n}{n^{3/2}}\le\frac12.
\]

### 1.4 Exact augmented-cut-code formulation

Let \(N=\binom n2\) and

\[
\mathcal C_n=
\left\{
(t+z_i+z_j)_{i<j}:t,z_i\in\mathbb F_2
\right\}.
\]

This is the cut code of \(K_n\) together with its complement. If
\(\rho(\mathcal C_n)\) is its covering radius, then

\[
\boxed{M_n=N-2\rho(\mathcal C_n).}
\]

So the problem asks whether the \(n^{3/2}\)-scale deficit of this covering
radius has a limit.

### 1.5 Elementary finite-size facts

\[
M_n\le M_{n+1}\le M_n+n.
\]

- The first inequality follows by restricting an \((n+1)\)-vertex signing
  and maximizing over the last spin.
- The second follows by extending an optimal \(n\)-vertex signing with any
  new signed row; the new linear term has absolute value at most \(n\).

These bounds are too weak for convergence: \(O(n)\) increments permit
\(O(1)\) changes after normalization across windows of only \(O(\sqrt n)\)
vertices.

### 1.6 Exact block-gluing obstruction

For

\[
G=\begin{pmatrix}A&B\\B^\top&D\end{pmatrix},
\]

maximizing over the relative global sign of the two blocks yields the exact
identity

\[
\max_{x,y}|H_A(x)+H_D(y)+x^\top By|
=
\max_{x,y}
\left(
|H_A(x)+H_D(y)|+|x^\top By|
\right).
\]

Thus cross edges cannot cancel internal energy. A successful gluing theorem
must anti-align the large bilinear values of \(B\) with the high-energy
layers of both blocks; a scalar inequality involving only \(M_n\) cannot
express this.

### 1.7 Universal field-plus-spin lower bound

For every sequence of symmetric zero-diagonal sign matrices,

\[
\boxed{
\liminf_{n\to\infty}\frac{Q(A_n)}{n^{3/2}}
\ge c_*
=0.672986728863\ldots
}
\]

and therefore

\[
\boxed{
\liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\ge\frac{c_*}{2}
=0.336493364431\ldots .
}
\]

#### Spectral bootstrap

Let \(\lambda=\|A\|_{\mathrm{op}}\), and choose a unit eigenvector \(v\)
whose eigenvalue has absolute value \(\lambda\). The eigenvector equation
gives

\[
\lambda\|v\|_\infty\le\|v\|_1\le\sqrt n.
\]

Set \(c=\|v\|_\infty^{-1}\), and choose independent signs \(X_i\) with
\(\mathbb E X_i=cv_i\). Since \(A\) has zero diagonal,

\[
\mathbb E X^\top AX=c^2v^\top Av,
\]

so

\[
\boxed{Q(A)\ge\frac{\|A\|_{\mathrm{op}}^3}{n}.}
\]

Thus, on a sequence with \(Q(A)=O(n^{3/2})\),
\[
\|A\|_{\mathrm{op}}=O(n^{5/6}).
\]
Writing \(m=n-1\) and
\[
q_{ij}=\frac{(A^2)_{ij}}m,
\]
we obtain
\[
\frac1{n(n-1)}\sum_{i\ne j}q_{ij}^2
\le
\frac{\operatorname{tr}A^4}{m^2n(n-1)}
\le
\frac{\|A\|_{\mathrm{op}}^2}{m^2}
=O(n^{-1/3}).
\]
Hence distinct rows are asymptotically orthogonal in mean square.

#### Smoothed opposite-orientation rounding

Fix \(t,\tau>0\), take independent Rademachers \(\xi_i\) and Gaussians
\(Z_i\), put
\[
h=\frac{A\xi}{\sqrt m},
\qquad
Y_i^\sigma
=\operatorname{sgn}(\sigma h_i+t\xi_i+\tau Z_i),
\quad \sigma\in\{\pm1\}.
\]

For a pair \(i\ne j\), isolate its direct coupling
\(\varepsilon=m^{-1/2}\), integrate out \(Z\), and Taylor-expand the
exact smoothed pair response before making any CLT approximation.
The two orientations have identical baselines, which cancel exactly.
The first-order response is

\[
a_{ij}\bigl(F_+(\varepsilon)-F_-(\varepsilon)\bigr)
=
4\varepsilon\,\mathbb E[\alpha_\tau(u)\beta_\tau(v)]
+O_{t,\tau}(\varepsilon^2),
\]

where

\[
\psi_\tau(z)=2\Phi(z/\tau)-1,
\]
\[
\alpha_\tau(u)
=\frac{\psi_\tau'(u+t)+\psi_\tau'(u-t)}2,
\qquad
\beta_\tau(v)
=\frac{\psi_\tau(v+t)-\psi_\tau(v-t)}2.
\]

A smooth two-dimensional Lindeberg replacement has error
\(O_{t,\tau}(m^{-1/2})\) after this Taylor extraction. Gaussian
covariance interpolation and the mean-square estimate for \(q_{ij}\)
then give, on average over pairs,

\[
\mathbb E[\alpha_\tau(u)\beta_\tau(v)]
=
2\phi_{1+\tau^2}(t)
\left[
2\Phi\!\left(\frac{t}{\sqrt{1+\tau^2}}\right)-1
\right]
+o(1).
\]

Summing over ordered pairs and using
\(2Q(A)\ge\max x^\top Ax-\min x^\top Ax\) yields

\[
\frac{Q(A)}{n\sqrt{n-1}}
\ge
4\phi_{1+\tau^2}(t)
\left[
2\Phi\!\left(\frac{t}{\sqrt{1+\tau^2}}\right)-1
\right]
-o(1).
\]

Take \(n\to\infty\) first, then \(\tau\downarrow0\). The unique positive
optimizer \(t_*\) solves

\[
2\phi(t_*)=t_*\bigl(2\Phi(t_*)-1\bigr),
\qquad
t_*=0.876902\ldots,
\]

and

\[
c_*=4\phi(t_*)\bigl(2\Phi(t_*)-1\bigr)
=0.672986728863\ldots .
\]

The order of operations is essential: a direct Berry–Esseen estimate on
the unexpanded pair response has leading-order error.

### 1.8 Local continuity in the order

For \(N=n+h\), random cross edges give

\[
\boxed{
M_{n+h}\le M_n+M_h+
\sqrt{2nh(n+h+2)\log2}.
}
\]

Together with monotonicity, this implies, uniformly for \(h=o(n)\),

\[
\boxed{
\frac{M_{n+h}}{(n+h)^{3/2}}
-
\frac{M_n}{n^{3/2}}
=o(1).
}
\]

Consequently, if \(d_{k+1}/d_k\to1\) and
\(M_{d_k}/d_k^{3/2}\to c\), then the full sequence converges to \(c\).
The missing fact is convergence of the minima on any such ratio-dense
subsequence.

### 1.9 Conference heavy-row dichotomy

Let \(D\) be a switching of a symmetric conference matrix for which
\[
R=\mathbf1^\top D\mathbf1=\max_x x^\top Cx,
\qquad r=D\mathbf1,
\qquad m=n-1.
\]
Then \(r_i>0\), \(\|r\|_2^2=nm\), and \(Dr=m\mathbf1\).

For any \(p\in[0,1]^n\), put each vertex independently in a random set
\(S\) with probabilities \(p_i\). All-cut positivity gives
\[
\mathbb E\,c_D(S,S^c)
=r\cdot p-p^\top Dp
\le \frac R2.
\]
Taking \(p=tr\) gives the exact inequality
\[
mt(n-tR)\le\frac R2
\qquad
\left(0\le t\le\frac1{r_{\max}}\right).
\]

Consequently:

- If \(r_{\max}\le 2R/n\), then \(t=n/(2R)\) is feasible and
  \[
  \boxed{R\ge n\sqrt{\frac{n-1}{2}}.}
  \]
- In general, writing
  \[
  h=\frac{r_{\max}}{\sqrt{n-1}},
  \qquad
  \rho=\frac{R}{n\sqrt{n-1}},
  \]
  the choice \(t=1/r_{\max}\) gives
  \[
  \boxed{\rho\ge\frac{2h}{h^2+2}.}
  \]

Thus any conference sequence with \(\rho<1/\sqrt2-o(1)\) must have a
genuinely heavy row. Converting that heavy row into a stronger Boolean
witness remains open.

For a conference matrix of order \(n=4\ell+2\), all row sums in any
switching are congruent modulo \(4\). At a maximizing switching they are
positive odd integers in one common residue class. This supplies useful
arithmetic rigidity but has not yet closed the heavy-row case.

### 1.10 Audited small-order insertion obstruction

Exhaustive enumeration after fixing the first row by switching gives:

\[
M_5=4,\qquad M_6=5.
\]

For every gauge-class optimizer found:

- at \(n=5\), \(20\) of the \(32\) Boolean vectors are exact extremizers;
- at \(n=6\), \(24\) of the \(64\) Boolean vectors are exact extremizers.

Define the best one-vertex extension profile
\[
E(A)=
\min_{b\in\{\pm1\}^n}
\max_{x\in\{\pm1\}^n}
\left(|H_A(x)|+|b\cdot x|\right).
\]
All twelve gauge-fixed order-\(6\) minimizers satisfy
\[
\boxed{E(A)=9.}
\]
Thus their best one-vertex extension jumps by \(4\), whereas the
derivative-scale target \(3M_6/(2\cdot6)\) is only \(1.25\). This does not
disprove an asymptotic insertion theorem, but it rules out the hoped-for
uniform finite theorem and shows that low entropy of the extremal layer is
not automatic.

### 1.11 Reduction to primes \(1\bmod 4\)

Let
\[
L_n=\frac{M_n}{n^{3/2}}.
\]
The full limit exists if and only if \(L_p\) converges as \(p\to\infty\)
through primes \(p\equiv1\pmod4\), and the two limits then agree.

Indeed, the prime number theorem in arithmetic progressions implies that
the consecutive primes in this progression have ratio tending to \(1\).
For every \(n\), choose consecutive such primes
\[
p_-(n)\le n\le p_+(n).
\]
Then \(p_\pm(n)/n\to1\), while monotonicity gives
\[
M_{p_-(n)}\le M_n\le M_{p_+(n)}.
\]
After division by \(n^{3/2}\), convergence on the prime subsequence
squeezes the full sequence. The converse is immediate.

This is a useful localization of the problem, but it does not itself
compare different prime orders.

### 1.12 Paley square-wave resonance gives the spectral limsup

Let \(p\equiv1\pmod4\) be prime and let
\[
S_{jk}=\chi_p(j-k),\qquad j,k\in\mathbb F_p,
\]
with zero diagonal. In the unitary Fourier normalization,
\[
x^\top Sx
=
p^{-1/2}\sum_{m\ne0}\chi_p(m)|\widehat x(m)|^2.
\]

For the square wave
\[
x_j=\operatorname{sgn}\cos(2\pi j/p),
\]
its unnormalized Fourier coefficients satisfy, on symmetric
frequencies,
\[
\widehat x(m)
=
\frac{(-1)^{(m-1)/2}}
{\sin(\pi m/(2p))}
\quad(m\ {\rm odd}),
\]
while the total even-frequency contribution is asymptotically
negligible. Consequently, for each fixed odd \(m\), the pair of
frequencies \(\pm m\) carries asymptotic Fourier mass
\[
\frac{8}{\pi^2m^2}.
\]

Fix \(L\). Choose a residue class
\[
p\equiv1
\pmod{8\prod_{\substack{\ell\le2L+1\\\ell\ {\rm odd\ prime}}}\, .
\]
Quadratic reciprocity makes every positive odd
\(m\le2L+1\) a quadratic residue modulo \(p\), and Dirichlet's theorem
supplies infinitely many primes in this class. Since
\[
\frac8{\pi^2}\sum_{m\ {\rm odd}>0}\frac1{m^2}=1,
\]
first taking \(p\to\infty\) in the class and then \(L\to\infty\) gives
\[
\frac{x^\top Sx}{p^{3/2}}\to1.
\]
Adding the extra row and column to form the Paley conference matrix
changes the quadratic energy by only \(O(p)\). Hence, in the original
half-quadratic normalization,
\[
\boxed{
\limsup_{\substack{p\to\infty\\p\equiv1(4)}}
\frac{\max_x|H_{\rm Paley,p}(x)|}{p^{3/2}}
=\frac12.
}
\]

This theorem shows that exact spectral flatness does not prevent sparse
arithmetic resonance. It does **not** show that the Paley values fail to
converge; a rigorously controlled nonresonant subsequence is still
missing.

The resonance is not a zero-density curiosity. The exact square-wave
formula is
\[
\frac{x^\top Sx}{p^{3/2}}
=
\frac8{\pi^2}
\sum_{\substack{h\ge1\\h\ {\rm odd}}}
\frac{\chi_p(h)}{h^2}
+O(p^{-1}).
\]
If every odd \(h\le25\) is a quadratic residue, then even an adversarial
tail gives
\[
\frac{x^\top Sx}{p^{3/2}}
\ge
\frac{16}{\pi^2}
\sum_{\substack{h\le25\\h\ {\rm odd}}}\frac1{h^2}
-1-o(1)
=0.9688395921\ldots-o(1),
\]
strictly above the Haar doubled benchmark
\(\sqrt{15}/4=0.9682458366\ldots\).
It suffices to force
\[
3,5,7,11,13,17,19,23
\]
to be quadratic residues. Among primes \(p\equiv1\pmod4\), these
conditions have relative Dirichlet density \(2^{-8}=1/256\).
Therefore the external claim that Paley values converge to the Haar
constant on a density-one set of primes is false. A useful upper-bound
route needs only a ratio-dense *good* prime subsequence, but that still
requires an all-Boolean nonresonance theorem.

There is now a stronger obstruction. Fix any admissible arithmetic
progression of primes \(p\equiv1\pmod4\), or equivalently any finite
compatible prescription of Legendre symbols together with congruence
conditions. For every \(\delta>0\), a fixed refined progression inside
it has Boolean Paley witnesses with doubled normalized energy at least
\(1-\delta\). Consequently,
\[
\boxed{
\limsup_{\substack{p\to\infty\\p\ \mathrm{in\ any\ fixed\
admissible\ progression}}}
\frac{\max_x|x^TA_px|}{2p^{3/2}}
=\frac12.
}
\]

The construction is explicit. For each prescribed nonresidue prime
\(\ell\), choose \(u\in\{\pm1\}^{\ell}\) with \(\sum u=1\), write
\[
j=\sum_{d=0}^{2r}j_d\ell^d,
\qquad
v_j=\prod_{d\ {\rm even}}u_{j_d}
\quad
\left(j\bmod\ell^{2r+1}\right).
\]
Its normalized DFT vanishes whenever the frequency has odd
\(\ell\)-adic valuation, while its DC mass is only
\(\ell^{-2r-2}\). CRT-tensor these gadgets over the finitely many
prescribed nonresidues. The resulting circle step function has
arbitrarily close to all of its Fourier mass on frequencies where
every bad-prime valuation is even. After prescribing the finitely many
remaining small prime symbols to be positive, Dirichlet's theorem and
Riemann sampling transfer the witness to \(\mathbb F_p\).

This proves that no fixed congruence class, and no finite-character
inverse theorem, can yield a Paley upper bound below \(1/2\). A
hypothetical low-valued Paley subsequence must use an increasing,
\(p\)-dependent amount of arithmetic information. The full proof and
normalization audit are in `paley_resonance_gadget.md`.

### 1.13 Exact opposite-orientation \(A^2\)-energy theorem

For every signing \(A\), every \(t,\tau\), and the smoothed witnesses
\[
X_i^\sigma
=\operatorname{sgn}\!\left(
\sigma\frac{(A\xi)_i}{\sqrt{n-1}}+t\xi_i+\tau Z_i
\right),
\qquad \sigma\in\{\pm1\},
\]
one has the exact finite-\(n\) inequality
\[
\boxed{
\frac12\sum_{\sigma=\pm1}
\mathbb E\|AX^\sigma\|_2^2
\ge n(n-1).
}
\]

To prove it, expand
\[
\frac12\sum_\sigma\mathbb E\|AX^\sigma\|_2^2
=\operatorname{tr}A^2+
\sum_{i\ne j}(A^2)_{ij}C_{ij},
\qquad
C_{ij}=\frac12\sum_\sigma
\mathbb E[X_i^\sigma X_j^\sigma].
\]
For a fixed pair \(i,j\), switch coordinates so \(a_{ij}=1\) and the
common coefficients in row \(i\) are all \(+1\). Write the corresponding
row-\(j\) coefficients as \(c_k\in\{\pm1\}\). Then
\[
d=\sum_{k\ne i,j}c_k=(A^2)_{ij}
\]
in this gauge. Conditioning on every variable except one common
Rademacher shows that \(C_{ij}\) is nondecreasing in each \(c_k\): the
change from an opposite to an equal coupling is a product of two
nonnegative smoothed-threshold increments. A measure-preserving
transformation that flips all common Rademachers, \(\sigma\), the
\(j\)-spin, and the \(j\)-dither shows
\[
C(-c)=-C(c).
\]
Hence \(C(c)\) has the sign of \(\sum c_k=d\), so every off-diagonal
summand satisfies
\[
(A^2)_{ij}C_{ij}\ge0.
\]
Switching preserves this product, completing the proof.

This supplies the full conference-scale squared local field for the
two opposite-orientation witnesses **on average**, without a CLT,
pseudorthogonality, or a near-minimizer hypothesis. The remaining
selection problem is to combine high oriented energy and high
\(A^2\)-energy in one witness and then exploit the only surviving case:
cut-stable, positive, heavy local fields.

### 1.14 Improved spectral bootstrap by asymmetric Boolean rounding

For every symmetric zero-diagonal matrix \(A\), not only sign matrices,
\[
\boxed{Q(A)\ge\frac12\|A\|_{\rm op}^2.}
\]

First let
\[
B(A)=\max_{x,y\in\{\pm1\}^n}|x^TAy|.
\]
Polarization writes \(x^TAy\) as the difference of the quadratic forms
of the two disjointly supported vectors
\((x+y)/2,(x-y)/2\in\{0,\pm1\}^n\). Randomly completing every zero
coordinate to a sign shows that the absolute quadratic form of either
partial sign vector is at most \(Q(A)\). Hence
\[
B(A)\le2Q(A).
\]

Choose a unit eigenvector \(v\) with eigenvalue
\(|\lambda|=\|A\|_{\rm op}\), take \(x=\operatorname{sign}v\), and choose
an independent Boolean vector \(Y\) with
\[
\mathbb EY_i=\frac{v_i}{\|v\|_\infty}.
\]
Then
\[
B(A)\ge
\left|\mathbb E\,x^TAY\right|
=|\lambda|\,\frac{\|v\|_1}{\|v\|_\infty}.
\]
At a coordinate where \(|v_i|=\|v\|_\infty\), the eigenvector equation
and \(|a_{ij}|\le1\) give
\[
|\lambda|\|v\|_\infty\le\|v\|_1.
\]
Combining the last three displays proves the claim.

Thus any competing sequence with \(Q(A)=O(n^{3/2})\) satisfies
\[
\boxed{\|A\|_{\rm op}=O(n^{3/4}),}
\]
improving the \(O(n^{5/6})\) bootstrap in §1.7. Consequently, for
\(q_{ij}=(A^2)_{ij}/(n-1)\),
\[
\frac1{n(n-1)}\sum_{i\ne j}q_{ij}^2
\le\frac{\|A\|_{\rm op}^2}{(n-1)^2}
=O(n^{-1/2}).
\]
The remaining gap is qualitative: \(O(n^{3/4})\) is still much larger
than the \(O(\sqrt n)\) operator-norm regime where squared local fields
are uniformly \(O(n^2)\) for every Boolean witness.

### 1.15 Joint selection, localized spectral anomalies, and capped
profiles

For the opposite-orientation field-plus-spin law, orient the sampled
energy as
\[
R=\sigma X^TAX,\qquad S=X^TA^2X.
\]
If a competing sequence has
\[
Q=(c_*+o(1))n^{3/2},
\qquad
\mathbb ER\ge(c_*-o(1))n^{3/2},
\]
then \(R\le Q\) pointwise implies
\[
\mathbb E(Q-R)=o(n^{3/2}).
\]
A near-\(Q\) sample automatically has
\[
S\ge R^2/n=(c_*^2-o(1))n^2
\]
and negligible negative local-field mass by the correction in §3.22.
The exact average \(\mathbb ES\ge n(n-1)\) upgrades this to
\[
S\ge(1-o(1))n^2
\]
for the same near-maximal sample whenever \(S/n^2\) is uniformly
integrable. In particular this holds if
\[
\|A\|_{\rm op}=O(\sqrt n).
\]

Every failure of this regular condition is localized. If
\(Av=\lambda v\), \(\|v\|_2=1\), and \(Q=Q(A)\), then for every
\(\theta>1\), with
\[
s^2=\frac{\lambda}{\theta Q},
\qquad
T=\{i:|v_i|>s\},
\]
one has
\[
\boxed{
|T|\le\frac{\theta Q}{\lambda},
\qquad
\|v_T\|_2^2\ge\frac{1-\theta^{-1}}3.
}
\]
Indeed, writing \(v=u+w\) on \(T\cup T^c\),
\[
w^TAw=\lambda(1-2\|u\|_2^2)+u^TAu
\ge\lambda(1-3\|u\|_2^2),
\]
whereas \(w/s\in[-1,1]^n\) gives
\[
|w^TAw|\le Qs^2=\lambda/\theta.
\]
Thus a mode \(\lambda=L_n\sqrt n\), \(L_n\to\infty\), has a fixed
fraction of its mass on \(O(n/L_n)=o(n)\) vertices.

There is a quantitative theorem after a successful peeling/capping
step. Switch a near-maximal witness to \(\mathbf1\), write
\[
r=D\mathbf1,\quad
c=\frac{\mathbf1^TD\mathbf1}{n^{3/2}},\quad
s=\frac{\|r\|_2^2}{n^2},
\]
and assume
\[
\|D\|_{\rm op}\le K\sqrt n,\qquad
\max_i|r_i|\le H\sqrt n,
\]
with negligible negative-field mass. The spectral measure of
\(\mathbf1\) has first three normalized moments \(c,s,t\) on
\([-K,K]\). From
\[
(K-z)(z-a)^2\ge0
\]
and optimization in \(a\),
\[
t\le
M_K(c,s)
:=
Ks-\frac{(s-Kc)^2}{K-c}.
\]
Applying all-cut positivity to independent probabilities
\(p_i=\alpha(r_i)_+/\sqrt n\) gives
\[
\boxed{
\frac c2\ge
\max_{0\le\alpha\le1/H}
\left[\alpha s-\alpha^2M_K(c,s)\right],
\qquad
c\ge\frac{s}{H}.
}
\]
For the coefficient-one, flat-spectrum case \(s=K=1\), a sequence
saturating \(c=c_*\) must therefore retain a local field at least
\[
\boxed{(1.941916296\ldots-o(1))\sqrt n.}
\]
The improved cap remains nontrivial only up to
\[
K<1.0220798875\ldots.
\]
This reduces the lower-bound problem to peeling localized heavy
coordinates without allowing a succession of unrelated maximizing
states. The detailed proofs are in
`joint_selection_and_spectral_localization.md`.

### 1.16 Exact orientation-even \(A^2\)-energy gain

The unsmoothed opposite-orientation Gaussian witnesses from §1.2 obey
a stronger finite identity than §1.13 detects. Put
\[
m=n-1,\qquad
X^\sigma=\operatorname{sign}\bigl((\sigma A+\sqrt m\,I)g\bigr).
\]
Then
\[
\boxed{
\frac12\sum_{\sigma=\pm1}
\mathbb E\,(X^\sigma)^TA^2X^\sigma
\ge
nm+
\frac{\|A^2-mI\|_F^2}
{\pi m\sqrt{1-1/m}}
}
\qquad(m>1).
\]

For \(i\ne j\), write
\[
c=(A^2)_{ij},\qquad
u=\frac{c}{2m},\qquad
v=\frac{a_{ij}}{\sqrt m}.
\]
The pair correlation averaged over the two orientations is
\[
C_{ij}
=
\frac1\pi
\left[
\arcsin(u+v)+\arcsin(u-v)
\right].
\]
If
\[
G_v(u)=\arcsin(u+v)+\arcsin(u-v),
\]
then \(G_v(0)=0\), and convexity of
\[
g(z)=(1-z^2)^{-1/2}
\]
gives
\[
G_v'(u)=g(u+v)+g(u-v)
\ge2g(v)
=\frac2{\sqrt{1-1/m}}.
\]
The inequality reverses in the correct way when integrating to
negative \(u\), so
\[
cC_{ij}
\ge
\frac{c^2}{\pi m\sqrt{1-1/m}}.
\]
Summing the off-diagonal terms and adding
\(\operatorname{tr}A^2=nm\) proves the theorem.

This is a leading-scale orientation-even statistic:
a Wigner-scale defect
\[
\|A^2-mI\|_F^2=\Theta(n^3)
\]
raises the averaged squared local field by \(\Theta(n^2)\). The
remaining task is joint selection with high oriented energy for these
same witnesses, followed by control of heavy positive coordinates.

### 1.17 Exact capped-field conversion

The joint-selection loss in §1.16 can be isolated in one explicit tail
term.  For a Boolean witness \(x\), orient \(A\) so that
\[
q=x^TAx,\qquad r_i=x_i(Ax)_i,
\]
and put
\[
S_K=\sum_i r_i\,\operatorname{clip}
\bigl(r_i,[-K\sqrt m,K\sqrt m]\bigr),\qquad m=n-1.
\]
Switching \(x\) to \(\mathbf1\), define
\[
u_i=\operatorname{clip}
\bigl(r_i,[-K\sqrt m,K\sqrt m]\bigr),\qquad
\mu=(1-\alpha)\mathbf1+\frac{\alpha u}{K\sqrt m}.
\]
For \(0\le\alpha\le1\), \(\mu\in[-1,1]^n\).  Multilinear Boolean
rounding and \(|u^TAu|\le K^2mQ(A)\) give the exact inequality
\[
\boxed{
(1+\alpha^2)Q(A)
\ge
(1-\alpha)^2q+
\frac{2\alpha(1-\alpha)}{K\sqrt m}S_K .
}
\]
This inequality may be averaged over a witness distribution without
selecting one sample having both large \(q\) and large \(S_K\).

Write
\[
C=\frac{Q(A)}{n\sqrt m},\qquad
c=\frac{\mathbb E q}{n\sqrt m},\qquad
z=\frac{\mathbb E S_K}{Knm}.
\]
Optimization in \(\alpha\) yields
\[
\boxed{
C\ge{\cal F}(c,z),\qquad
{\cal F}(c,z)=
\begin{cases}
c,&z\le c,\\[1mm]
c-z+\sqrt{z^2+(z-c)^2},&z>c.
\end{cases}
}
\]
For the opposite-orientation Gaussian witnesses of §§1.2 and 1.16,
\[
c\ge\frac2\pi+o(1).
\]
If
\[
\delta_n=\frac{\|A^2-mI\|_F^2}{nm^2}
\]
and
\[
\Psi_n(K)=\frac1{nm}\,
\mathbb E\sum_i r_i^2
\mathbf1_{\{|r_i|>K\sqrt m\}},
\]
then §1.16 gives
\[
\boxed{
C\ge
{\cal F}\left(
\frac2\pi+o(1),
\frac{1+\delta_n/(\pi\sqrt{1-1/m})-\Psi_n(K)}K
\right).
}
\]
Thus the orientation-even defect already converts into a stronger
quadratic witness whenever a fixed cap retains enough squared local
field.  The remaining loss is exactly the positive heavy-field tail;
universal-positive-vertex examples show that it cannot be bounded from
the first two moments alone.

### 1.18 Regularized asymptotic near-minimizers

Raw near-minimizers need not be spectrally regular, but regular
near-minimizers can always be constructed with an explicit two-limit
tradeoff.  Let
\[
q_n=\min_AQ(A).
\]
For every \(K\ge1\), Grothendieck--Pietsch deletion followed by a
conference-type refill and random cross signing gives a full
order-\(n\) signing \(A'_{n,K}\) satisfying
\[
\boxed{
Q(A'_{n,K})
\le q_n+O(K^{-1/2}n^{3/2}),\qquad
\|A'_{n,K}\|_{\rm op}=O(K\sqrt n).
}
\]
Indeed, delete at most \(n/K\) vertices so that the retained principal
matrix has norm \(O(K\sqrt n)\).  The refill has internal norm
\(O((n/K)^{3/2})\), and a random cross block has Boolean bilinear norm
\(O(n^{3/2}/\sqrt K)\) and operator norm \(O(\sqrt n)\).

Letting \(K\to\infty\) arbitrarily slowly produces spectrally controlled
asymptotic near-minimizers.  This does not by itself control local-field
tails uniformly in \(K\), but it legitimizes proving lower bounds first
for fixed \(K\) and then tracking the dependence on \(K\).

There is also an exact block visibility inequality.  For
\[
A=\begin{pmatrix}D&B\\B^T&E\end{pmatrix}
\]
and any Boolean \(y\),
\[
\boxed{
Q(A)\ge |y^TEy|+\sqrt2\,\|By\|_2.
}
\]
If \(y\) is an oriented \(Q(E)\)-ground state and
\(x=\operatorname{sign}(By)\), then the sharper deterministic form is
\[
\boxed{
Q(A)-Q(E)
\ge2\max_{y\in{\rm GS}(E)}\|By\|_1-Q(D).
}
\]
Consequently, a replenishment tower can persist only if successive
cross-block singular directions avoid every ground-state frame of the
retained cores.

---

## 2. Important exact reformulations

### 2.1 Switching and signed cuts

For a signing \(A\), replacing \(A\) by
\(\operatorname{diag}(x)A\operatorname{diag}(x)\) is Seidel switching.
The value \(x^\top A x\) is the total excess after switching.

For a conference matrix \(C\), choose a switching \(D\) for which
\[
R=\mathbf1^\top D\mathbf1=\max_x x^\top Cx.
\]
If \(c_D(S,S^c)\) denotes the signed sum across the cut, global maximality
implies
\[
\boxed{
0\le c_D(S,S^c)\le R/2
\quad\text{for every }S\subseteq[n].
}
\]

Also, with \(r=D\mathbf1\),
\[
D^2=(n-1)I,\qquad
\|r\|_2^2=n(n-1),\qquad
Dr=(n-1)\mathbf1.
\]

A strong lower bound on \(R=\sum_i r_i\) from these cut constraints would
be a direct route to the conference case.

### 2.2 Projection form for conference matrices

For
\[
U=\frac{C}{\sqrt{n-1}},\qquad
P=\frac{I+U}{2},
\]
we have \(U^2=I\) and \(P\) is a rank-\(n/2\) projection. Moreover,
\[
\frac{x^\top Cx}{n\sqrt{n-1}}
=
\frac{x^\top Ux}{n}
=
\frac{2\|Px\|_2^2}{n}-1.
\]

Thus the spectral ceiling \(1\) is equivalent to finding a Boolean vector
at \(o(\sqrt n)\) Euclidean distance from one eigenspace.

---

## 3. Falsified or stopped routes

### 3.1 Falsified: naive repeated conference AMP reaches \(1\)

The earlier claim that repeated projection has state evolution
\[
0.67299\to0.68512\to0.69539\to\cdots\to1
\]
is false for a fixed conference matrix.

Let \(U=C/\sqrt{n-1}\), so \(U^2=I\). Start with
\[
G\sim N(0,I+U),\qquad X_0=\operatorname{sgn}G,
\]
and put
\[
a=\sqrt{2/\pi},\qquad \sigma^2=1-a^2,
\qquad
R=\frac{X_0-aG}{\sigma},\qquad W=UR.
\]

The first update is
\[
X_1=\operatorname{sgn}(X_0+UX_0).
\]

In the scalar limit, write
\[
Y=\operatorname{sgn}\bigl(\operatorname{sgn}Z+aZ+\sigma W\bigr),
\quad Z,W\stackrel{\mathrm{iid}}{\sim}N(0,1).
\]
Define
\[
b=\mathbb E[ZY],\qquad
c=\mathbb E[WY],\qquad
d=\mathbb E[RY],
\quad
R=\frac{\operatorname{sgn}Z-aZ}{\sigma}.
\]
Numerically,
\[
b=0.7920592175,\quad
c=0.05415383275,\quad
d=0.5726603594.
\]

The decisive algebraic obstruction is
\[
UW=R.
\]
If
\[
\eta=X_1-bG-cW,
\]
then \(\langle \eta,R\rangle\to d\ne0\), and hence
\[
\langle W,U\eta\rangle
=\langle UW,\eta\rangle
=\langle R,\eta\rangle
\to d.
\]
Therefore \(U\eta\) is not a fresh Gaussian residual. The omitted
backtracking term is order one.

The corrected paired prediction after one step is
\[
b^2+2cd=0.68938131\ldots,
\]
not the naive
\[
b^2+cd=0.65836956\ldots.
\]

A legitimate recursion must retain Gram–Schmidt pairs
\((R_s,W_s)\) satisfying
\[
UR_s=W_s,\qquad UW_s=R_s.
\]
Population recursion then gives approximately
\[
0.63662,\ 0.68938,\ 0.69937,\ 0.70254,\ 0.70377,\ldots
\]
and appears to converge near \(0.7054\), not \(1\).

### 3.2 Numerical confirmation of the AMP obstruction

For direct iteration
\[
x\mapsto\operatorname{sgn}((I+U)x)
\]
on Paley conference matrices, terminal ratios
\[
\frac{|x^\top Cx|}{n\sqrt{n-1}}
\]
cluster around \(0.70\), and fixed points usually appear after \(3\)–\(5\)
steps.

Independent current-work computations:

- Orders from \(102\) through \(402\): typical terminal means approximately
  \(0.69\)–\(0.73\).
- Single-spin greedy local search does substantially better, around
  \(0.91\)–\(0.94\) at orders \(102,194,402\), but this is still numerical
  and does not prove approach to \(1\).

### 3.3 Falsified/stopped: scalar Fekete argument

Splitting into comparable blocks introduces a rectangular minimax term of
order \((n+m)^{3/2}\), exactly the leading scale. Ordinary subadditivity
therefore does not imply convergence.

At free-energy scale \(t=\beta/\sqrt n\), splitting \(2n\) vertices changes
the internal-block inverse temperature from \(\beta\) to
\(\beta/\sqrt2\). Fixed-temperature Fekete subadditivity sees the wrong
diagonal scaling.

### 3.4 Stopped: universal local tensor/Hadamard gadget

Common Hadamard lifts immediately relax scalar Boolean choices into vector
choices. On the order-\(6\) conference seed, the relaxed profile already
jumps away from the scalar optimum. Fixed local block templates cannot
provide the required scale-preserving amplification.

### 3.5 Stopped: ordinary spin-glass interpolation

After random vertex switching, the pairwise-overlap covariance is the same
for every signing. Yet the Boolean optimum can range from order
\(n^{3/2}\) to order \(n^2\). A method using only that covariance is
information-theoretically unable to distinguish the relevant structures.

### 3.6 Stopped: spectral certificate alone

The operator-norm bound gives the upper construction constant \(1/2\), but
standard SDP/vector relaxations lose the scalar Boolean geometry and have a
hard spectral floor. They do not prove the matching universal lower bound.

### 3.7 Stopped: \(n=5\to10\) two-copy optimizer lift

The \(n=10\) optimizer can be represented as a signed two-copy lift of an
\(n=5\) optimizer, explaining its repeated spectrum, but iterating the lift
causes the normalized objective to grow:
\[
0.716\to0.822\to1.073
\]
for the corresponding \(Q/(n\sqrt n)\) sequence. It is not an asymptotic
minimizing construction.

### 3.8 Falsified: tensor submultiplicativity

There are explicit symmetric zero-diagonal \(5\times5\) signings

\[
A=\begin{pmatrix}
0&1&1&-1&-1\\
1&0&-1&-1&1\\
1&-1&0&-1&1\\
-1&-1&-1&0&-1\\
-1&1&1&-1&0
\end{pmatrix},
\quad
B=\begin{pmatrix}
0&-1&1&-1&-1\\
-1&0&-1&-1&-1\\
1&-1&0&1&-1\\
-1&-1&1&0&-1\\
-1&-1&-1&-1&0
\end{pmatrix}
\]

with
\[
Q(A)=Q(B)=12,
\]
but a Boolean witness gives
\[
\boxed{Q(A\otimes B)\ge176>144=Q(A)Q(B).}
\]

Thus the equal-modulus quadratic norm is not tensor-submultiplicative.
Ordinary Sidon-constant tensor machinery leaves the full-support,
equal-modulus subclass and cannot supply the missing scale transfer.

### 3.9 Falsified: generic centered-free-energy axioms suffice

Fixed-raw-temperature centered subadditivity, evenness, real analyticity,
convexity, the correct variance, and a correct zero-temperature slope do
not force convergence on \(t=\beta/\sqrt n\).

A concrete countermodel is obtained by choosing an oscillating
\[
a_n=c_n n^{3/2},
\qquad
c_n=1+\varepsilon\sin(\log\log(n+n_0)),
\]
with \(a_n/\binom n2\) decreasing, and defining
\[
F_n(t)=\log2+
\frac{a_n^2}{\binom n2}
\log\cosh\!\left(\frac{\binom n2}{a_n}t\right).
\]

Then \(F_n\) is even, analytic and convex,
\[
F_n''(0)=\binom n2,
\qquad
\lim_{t\to\infty}\frac{F_n(t)}t=a_n,
\]
and the centered functions
\[
G_n(t)=F_n(t)-\binom n2\log\cosh t
\]
are subadditive. Nevertheless,
\[
\frac{F_n(\beta/\sqrt n)}n
=
2c_n^2\log\cosh\!\left(\frac{\beta}{2c_n}\right)+o(1)
\]
oscillates. A successful free-energy proof must use discrete
higher-cumulant identities specific to quadratic Ising energies.

### 3.10 Stopped: ground-state entropy alone gives sharp insertion

For fixed \(A\), a random new signed row \(b\) can achieve extension
overhead at most \(t\) if a weighted energy-layer sum of the form
\[
\sum_x
\exp\!\left(
-\frac{(M-|H_A(x)|+t)^2}{2n}
\right)
<\frac12.
\]
An \(O(\sqrt n)\) conclusion from this criterion requires only \(O(1)\)
effective entropy near the ground-state layer. The audited \(n=5,6\)
optimizers in §1.10 have respectively \(20/32\) and \(24/64\) exact
extremizers, so this premise fails badly even at minimizing signings.

After switching a positive maximizer, the relevant slack layer is a family
of low signed cuts. Ordinary Karger cut-counting does not apply because a
signed cut function is not submodular. No asymptotic replacement is known.

### 3.11 Stopped: orientation-odd higher-cycle stability

A higher-Hermite refinement of field-plus-spin rounding produces signed
odd-cycle corrections. In normalized notation \(B=A/\sqrt{n-1}\) and
\(q_{ij}=(B^2)_{ij}\), its correction has the schematic form
\[
S_t(B)=
\sum_{i\ne j}B_{ij}K_t(q_{ij}),
\qquad
K_t(q)=\sum_{\ell\ {\rm odd}}k_\ell(t)^2q^\ell.
\]
The leading term is proportional to \(\operatorname{tr}(B^3)\); the
remaining terms are longer signed theta-cycle profiles.

This entire family has a decisive null class. If
\[
PBP^\top=-B
\]
for a permutation matrix \(P\), then \(PB^2P^\top=B^2\) and hence
\[
S_t(B)=-S_t(B)=0
\]
for every \(t\). Self-complementary signings of this kind can be
arbitrarily far from conference structure (random examples have
Wigner-like fourth moment). Therefore no orientation-odd, one-channel
higher-cycle correction can force \(B^2\approx I\). A viable successor
must be orientation-even and detect
\(\|B^2-I\|_F^2\) or equivalent information.

### 3.12 Exact-computation checkpoint at order \(10\)

A complete mixed-integer optimization gives
\[
M_{10}=13,\qquad Q_{10}=26.
\]
One optimizer has spectrum
\[
\pm3.933464\ldots\ (\text{twice}),\quad
\pm2.554969\ldots\ (\text{twice}),\quad
\pm1,
\]
and \(A^2-9I\) is supported on two signed \(5\)-cycles with nonzero
entries \(\pm4\). This explains the earlier observed repeated spectral
polynomial, but its natural two-copy continuation was already falsified
as an asymptotic minimizing lift.

### 3.13 Optimality ceiling for one-probe coordinatewise rounding

The field-plus-spin constant \(c_*\) in §1.7 is optimal within the full
class of one-probe coordinatewise Boolean rules.

Let \(G\sim N(0,1)\), \(S\) be an independent sign, and
\(f(G,S)\in\{\pm1\}\). Put
\[
a=\mathbb E[Gf(G,S)],
\qquad
b=\mathbb E[Sf(G,S)].
\]
The direct-edge response of the associated one-probe rounding is
\(2ab\). Write
\[
g(z)=\frac{f(z,+1)+f(z,-1)}2,
\qquad
h(z)=\frac{f(z,+1)-f(z,-1)}2.
\]
Pointwise,
\[
|g(z)|+|h(z)|\le1.
\]
After choosing signs so \(a,b\ge0\), rearrangement shows that for fixed
\(\mathbb E h(G)\), the product is maximized by using the spin on the
smallest values of \(|G|\) and the field on the largest values:
\[
f(z,s)=\operatorname{sgn}(z+ts).
\]
For this threshold rule,
\[
a=2\phi(t),
\qquad
b=2\Phi(t)-1.
\]
Therefore
\[
\sup_f 2ab
=
\max_{t\ge0}
4\phi(t)(2\Phi(t)-1)
=c_*.
\]
Any improvement over \(0.3364933644\ldots\) for \(M_n/n^{3/2}\) must
use multiple dependent probes, a genuinely nonlocal rule, or additional
structure of \(A\).

### 3.14 Stopped: correlated flat-Fourier lifts

For an abelian fiber group \(G\) of order \(k\), a broad class of
correlated block lifts has flat Fourier kernels
\[
\widehat b_{ij}(\chi)=\sqrt{k}\,a_{ij}^{(\chi)}.
\]
If
\[
m_{i,\chi}=\frac1k\sum_{g\in G}x_i(g)\chi(g),
\qquad
\sum_\chi |m_{i,\chi}|^2=1,
\]
then the cross energy is exactly
\[
k^{3/2}
\sum_\chi\sum_{i<j}
a_{ij}^{(\chi)}
m_{i,\chi}m_{j,\chi}.
\]
Thus the lift exposes a row-sphere, multi-channel relaxation.

If the channel signings form a balanced full signed-permutation orbit of
a seed \(A\), this relaxation is exactly the spectral value
\[
\frac n2\|A\|_{\rm op}.
\]
Top eigenvectors of the conjugates attain the lower bound after their
squared coordinates are averaged; the spectral inequality gives the
reverse bound. Hence the natural orbit-symmetrized lift necessarily
jumps to the spectral ceiling and cannot preserve a scalar optimum below
\(1/2\).

A finite audit at fiber order \(4\) gives the same obstruction: all
\(768\) Fourier-compatible quadruples of order-\(5\) minimizers have
lifted cross value \(36\) or \(40\), while exact preservation would
require \(32\). This computational datum should be treated as an audited
finite no-go, not as a theorem covering every nonabelian lift.

### 3.15 Stopped: radial signed cut-code dual certificates

Let \(C=\mathcal C_n\), \(N=\binom n2\), and
\[
D=C^\perp
=
\{\text{even-cardinality Eulerian edge sets of }K_n\}.
\]
For a signing \(A\), define
\[
T_w(A)=
\sum_{\substack{F\in D\\|F|=w}}
\prod_{e\in F}a_e.
\]
The complete signed dual enumerator satisfies the exact identity
\[
\sum_wT_w(A)z^w
=
(\cosh\beta)^{-N}2^{-n}
\sum_x\cosh(\beta q_A(x)),
\qquad z=\tanh\beta,
\]
where \(q_A(x)=\sum_{i<j}a_{ij}x_ix_j\).

Equivalently, the coset weight enumerator is
\[
\sum_{c\in C}u^{d(y,c)}
=
\frac12\sum_x
\left[
u^{(N-q_A(x))/2}
+
u^{(N+q_A(x))/2}
\right].
\]
Thus the full signed Eulerian/Krawtchouk data is invertibly identical to
the original switching-energy histogram.

For a ball indicator, the exact signed certificate uses
\[
\sum_{j\le R}K_j^N(w)=K_R^{N-1}(w-1),
\]
but demanding its positivity for every translate is term-for-term
equivalent to the original covering-radius statement. Fixed dual degree
sees only fixed moments of the energy distribution; conference
sequences have the same fixed Gaussian-chaos limits even when sparse
resonant vectors change their maxima. A resonance-sensitive certificate
must therefore have genuinely growing degree, where it again becomes
the full energy-histogram problem.

### 3.16 Scalar Eulerian-pressure axioms still permit oscillation

The correct two-sided finite-temperature object is
\[
\Gamma_n(\rho)
=
\min_A\log
\left[
(\cosh t)^{-\binom n2}
2^{-n}\sum_x\cosh(tH_A(x))
\right],
\qquad
\rho=\tanh t.
\]
Exact random-edge arguments show that \(\Gamma_n\) is:

- subadditive in \(n\) at fixed \(\rho\);
- nonincreasing in \(n\);
- nonincreasing in \(\rho\).

The tempting scaling inequality
\[
\Gamma_n(\lambda\rho)\le\lambda^2\Gamma_n(\rho)
\]
is false already at \(n=4\), where
\[
\Gamma_4(\rho)=\log(1-\rho^4).
\]

Moreover, all the listed scalar properties, even together with even
analyticity and a first signing-dependent term of order \(\rho^4\),
allow diagonal oscillation. For a slowly oscillating positive \(c_n\)
chosen so \(\theta_n=c_n/\sqrt n\) decreases, the abstract family
\[
\Gamma_n(\rho)
=
-n\,\frac{(\rho/\theta_n)^4}{1+(\rho/\theta_n)^4}
\]
has all those properties, while
\[
\frac1n\Gamma_n(\beta/\sqrt n)
=
-\frac{(\beta/c_n)^4}{1+(\beta/c_n)^4}
\]
oscillates. A free-energy proof must use coefficient-level Eulerian
constraints, an overlap hierarchy, or another discrete feature—not only
the scalar pressure inequalities.

### 3.17 Fixed conference diagrams are universal, but resonance is
nonperturbative

For a conference matrix put
\[
U=\frac C{\sqrt{n-1}},
\qquad U^2=I,
\qquad |U_{ij}|=(n-1)^{-1/2}.
\]
In the linked-cluster expansion of
\[
\log\mathbb E_x
\exp\!\left(\frac\beta2x^\top Ux\right),
\]
every spin-index vertex has even degree. Repeatedly contracting
degree-\(2\) vertices by \(U^2=I\) leaves either a fully contractible
cycle/cactus diagram or a core of minimum degree at least \(4\). For each
fixed core with \(r\) vertices and \(e\) edges,
\[
e\ge2r,
\]
so its tensor sum is \(O(n^r n^{-e/2})=O(1)\). Only the contractible
diagrams contribute \(O(n)\), and their values are fixed by \(U^2=I\)
and flatness.

Therefore every fixed Taylor coefficient of the pressure divided by
\(n\) is switching-independent and agrees with the Haar
half-involution/Random-Orthogonal-Model coefficient.

The conclusion fails nonperturbatively. If a conference family has a
Boolean \(+1\) eigenvector, then its normalized pressure is at least
\[
\frac\beta2-\log2+o(1).
\]
For the Haar annealed pressure
\[
p_H(\beta)
=
\sup_{|u|<1}
\left\{
\frac{\beta u}{2}
+\frac14\log(1-u^2)
\right\},
\]
one has \(p_H(8)=3.22233\ldots\), whereas the Boolean eigenvector gives
\(4-\log2=3.30685\ldots\). Thus a zero-entropy resonant state changes
the pressure by \(O(n)\) while being invisible to every fixed diagram.
Fixed-order perturbation theory cannot justify the zero-temperature ROM
transfer.

### 3.18 Stopped: speed-\(n^2\) disorder LDP as a shortcut

Let
\[
K_n(T)=
\#\left\{
A:\max_x|H_A(x)|\le T
\right\},
\qquad
p_n(c)=2^{-\binom n2}K_n(cn^{3/2}).
\]
If the feasible set is empty, then
\[
-n^{-2}\log p_n(c)=+\infty.
\]
If it is nonempty, switching and global sign already produce at least
\(2^{n-1}\) feasible signings, so
\[
-\frac1{n^2}\log p_n(c)
\le
\frac{\binom n2-n+1}{n^2}\log2
\longrightarrow\frac{\log2}{2}.
\]
Consequently, if \(c\) lies strictly between the liminf and limsup of
\(M_n/n^{3/2}\), the proposed rate has an infinite subsequence and a
bounded subsequence. Proving its extended-real limit for every \(c\)
would already prove the original convergence by support; the LDP is not
an easier preliminary.

There is also a finite-profile planting obstruction. Start with a
low-norm signing \(A\), choose a Boolean vector \(s\) with
\(|H_A(s)|\le\sqrt{\binom n2}\), and flip
\[
h=\lfloor\delta n^{3/2}\rfloor
\]
edges on which \(a_{ij}s_is_j=-1\). The new signing \(B\) satisfies
\[
H_B(s)=2\delta n^{3/2}+o(n^{3/2}),
\]
so its maximum changes at leading scale. But for a uniform spin \(X\),
\[
\mathbb E_X\bigl(H_B(X)-H_A(X)\bigr)^2=4h.
\]
Thus for every fixed replica number \(k\), the normalized joint
\(k\)-replica energy laws of \(A\) and \(B\) have Wasserstein distance
\[
O_k(n^{-3/4}),
\]
and every fixed signed-subgraph density also agrees asymptotically.
Fixed-replica and bounded-degree profiles cannot detect a planted
zero-entropy extremizer.

### 3.19 Stopped: deterministic Haar transfer from flatness alone

The Haar half-subspace constant cannot be transferred from any of the
usual deterministic isotropy assumptions.

1. Block copies of \(C_6/\sqrt5\) give constant-diagonal projections
   whose maximum projection fraction is
   \[
   \frac12\left(1+\frac{10}{6\sqrt5}\right)
   =0.872678\ldots<\beta_*=0.984122\ldots .
   \]
   Degeneracy allows an eigenbasis to be mixed across blocks so that
   its coherence is \(O(n^{-1/2})\). Thus constant diagonal plus an
   incoherent eigenbasis is insufficient.
2. A balanced Hadamard eigenspace has constant diagonal and a perfectly
   flat eigenbasis but contains Boolean columns, so its maximum
   projection fraction is \(1\).
3. Most strongly, the exact Paley conference projections in §1.12 have
   constant diagonal \(1/2\) and every off-diagonal magnitude exactly
   \(1/(2\sqrt p)\), yet the resonant subsequence has Boolean projection
   fraction tending to \(1\).

Any Haar comparison must control anti-alignment with all \(2^n\) cube
points or explicitly classify arithmetic resonance; flatness,
equiangularity, and delocalization do not suffice.

### 3.20 Independent multi-probe rounding has the same exact ceiling

Let
\[
G\sim N(0,I_k),
\qquad
S\sim{\rm Unif}\{\pm1\}^k,
\]
independently, and let \(f(G,S)\in[-1,1]\) be arbitrary. Define
\[
a=\mathbb E[Gf(G,S)],
\qquad
b=\mathbb E[Sf(G,S)].
\]
The first-order direct-edge response of the associated independent
multi-probe rounding is \(2a\cdot b\).

Put
\[
m(s)=\mathbb E_Gf(G,s),
\qquad
r^2=\mathbb E_Sm(S)^2.
\]
The coordinate functions \(S_j\) are orthonormal, so Bessel's inequality
gives
\[
\|b\|_2\le r.
\]
For fixed \(s\), the sharp Gaussian centroid inequality gives
\[
\left\|\mathbb E_G[Gf(G,s)]\right\|_2\le J(m(s)),
\]
where
\[
J(u)=
2\phi\!\left(
\Phi^{-1}\!\left(\frac{1+u}{2}\right)
\right),
\qquad 0\le u\le1.
\]
The function \(K(v)=J(\sqrt v)\) is concave on \([0,1]\). With
\(u=2\Phi(z)-1\), this is equivalent to
\[
2\Phi(z)-1\ge2z\phi(z),
\]
whose two sides agree at \(z=0\) and whose difference has derivative
\(2z^2\phi(z)\ge0\). Hence
\[
\|a\|_2
\le\mathbb E J(m(S))
\le J(r).
\]
It follows that
\[
2a\cdot b
\le2rJ(r)
\le\max_{0\le r\le1}2rJ(r)
=c_*.
\]
Equality is attained by a one-spin threshold/dictatorship
\[
f(G,S)=\operatorname{sgn}(G_1+t_*S_1).
\]
Thus arbitrarily many fresh independent channels, arbitrary response
boundaries, and their mixtures cannot improve the verified
\(0.3364933644\ldots\) lower bound. An improvement must use dependent
matrix probes and retain their backtracking/Onsager structure.

### 3.21 Exact edge-flip/deep-hole certificate

For the augmented cut vectors
\[
\mathcal V=\{\pm(x_ix_j)_{i<j}:x\in\{\pm1\}^n\},
\]
write
\[
h_v=a\cdot v,\qquad
M=\max_{v\in\mathcal V}h_v,
\qquad
g_v=\frac{M-h_v}{2},
\]
and
\[
N_v=\{e:a_ev_e=-1\}.
\]
If \(a^S\) is obtained by flipping the edge set \(S\), then
\[
\boxed{
M(a^S)
=
M-2\min_{v\in\mathcal V}
\left[
g_v+|S|-2|S\cap N_v|
\right].
}
\]
Therefore \(a\) is stable under every flip set of size at most \(k\) if
and only if, for every such \(S\), some \(v\) satisfies
\[
\boxed{
2|S\cap N_v|\ge |S|+g_v.
}
\]
Only energy layers \(g_v\le|S|\) can certify a \(k\)-edge flip. For a
single edge, the negative supports from the top two energy layers
\(g=0,1\) must cover every coordinate.

For the exact order-\(10\) optimizer, the \(g=0\) layer has \(40\)
gauge-fixed states with \(|N_v|=16\), and the \(g=1\) layer has \(80\)
with \(|N_v|=17\). Every set of at most four edges is certified by the
active layer. Thirteen five-edge sets fail active-only certification;
they are perfect matchings, and the \(g=1\) layer rescues all of them at
equality.

Sparse-flip stability is not sufficient for global optimality. The
following order-\(11\) signing has \(M=19\) and is stable under every
set of at most four edge flips:
\[
\begin{pmatrix}
0&1&-1&-1&-1&-1&1&1&1&1&-1\\
1&0&1&1&-1&1&1&1&1&1&1\\
-1&1&0&-1&-1&-1&1&-1&-1&1&-1\\
-1&1&-1&0&-1&-1&1&1&1&-1&1\\
-1&-1&-1&-1&0&1&-1&-1&1&1&-1\\
-1&1&-1&-1&1&0&-1&1&-1&-1&1\\
1&1&1&1&-1&-1&0&1&-1&1&-1\\
1&1&-1&1&-1&1&1&0&-1&-1&1\\
1&1&-1&1&1&-1&-1&-1&0&-1&-1\\
1&1&1&-1&1&-1&1&-1&-1&0&-1\\
-1&1&-1&1&-1&1&-1&1&-1&-1&0
\end{pmatrix}.
\]
The flip set
\[
\{(0,3),(0,7),(1,6),(2,7),(5,9)\}
\]
lowers its value to \(17\), and an independent order-\(11\) signing
with value \(17\) is
\[
\begin{pmatrix}
0&1&-1&-1&1&1&-1&-1&-1&1&1\\
1&0&1&1&-1&-1&1&-1&-1&1&-1\\
-1&1&0&-1&-1&-1&-1&-1&-1&-1&1\\
-1&1&-1&0&-1&-1&1&1&-1&-1&1\\
1&-1&-1&-1&0&-1&1&-1&-1&-1&-1\\
1&-1&-1&-1&-1&0&1&1&-1&-1&1\\
-1&1&-1&1&1&1&0&-1&-1&-1&-1\\
-1&-1&-1&1&-1&1&-1&0&1&-1&-1\\
-1&-1&-1&-1&-1&-1&-1&1&0&1&-1\\
1&1&-1&-1&-1&-1&-1&-1&1&0&-1\\
1&-1&1&1&-1&1&-1&-1&-1&-1&0
\end{pmatrix}.
\]
Exhaustive verification over
\(\sum_{j=1}^4\binom{55}{j}\) flip sets gives minimum post-flip values
\(21,19,19,19\) for sizes \(1,2,3,4\), respectively. Thus even a
fourth-order local certificate can be trapped above the true minimum.

### 3.22 Dependent negative-field correction and its stable branch

Fix a signing \(A\), a Boolean vector \(x\), and write
\[
q=x^TAx,\qquad
r_i=x_i(Ax)_i,\qquad
L_-=\sum_i(-r_i)_+,\qquad Q=Q(A).
\]
Let \(B=\{i:r_i<0\}\), and independently flip each coordinate in \(B\)
with probability \(p\). For the resulting Boolean vector \(y\),
\[
\mathbb E\,y^TAy
=q+4pL_-+4p^2x_B^TA_{BB}x_B.
\]
Randomly completing \(x_B\) outside \(B\) shows
\[
|x_B^TA_{BB}x_B|\le Q,
\]
and hence
\[
Q\ge q+\max_{0\le p\le1}
\left(4pL_--4p^2Q\right).
\]
In the branch \(L_-\le2Q\), the optimizer is
\(p=L_-/(2Q)\), so
\[
\boxed{\quad
Q\ge q+\frac{L_-^2}{Q},
\qquad
L_-^2\le Q(Q-q).
\quad}
\]
Thus any witness with \(q=Q-o(n^{3/2})\) that cannot be improved by this
dependent second stage has \(L_-=o(n^{3/2})\).

There is also an exact quantitative description of this stable branch.
Put
\[
R=\sum_i r_i=q,\qquad
H=\max_i(r_i)_+,\qquad
S=\sum_i r_i^2=x^TA^2x.
\]
Since \(r_i^2\le Hr_i\) on the positive fields and
\(r_i^2\le(n-1)(-r_i)\) on the negative fields,
\[
\boxed{\quad
S\le H(R+L_-)+(n-1)L_-,
\qquad
H\ge\frac{S-(n-1)L_-}{R+L_-}.
\quad}
\]
Consequently, if \(q\sim c n^{3/2}\), \(Q-q=o(n^{3/2})\), and this
adaptive witness additionally satisfies \(S=(1-o(1))n^2\), then either
the negative-field correction improves the leading constant or
\[
H\ge\left(\frac1c-o(1)\right)\sqrt n.
\]
At \(c=c_*=0.672986728863\ldots\), this is
\(H>1.485\sqrt n\), more than \(2.20\) times the average field
\(R/n\).

Section 1.13 has now replaced the proposed four-point estimate with the
stronger exact opposite-orientation average
\[
\frac12\sum_\sigma\mathbb E[(X^\sigma)^TA^2X^\sigma]
\ge n(n-1).
\]
The unresolved issue is joint sample selection and the heavy positive
branch: a single heavy positive field is not harvested by the
negative-coordinate flip.

A broader exact dependent rule reaches the same stopping point. Clip
\[
u_i=\operatorname{clip}
\left(r_i,[-K\sqrt{n-1},K\sqrt{n-1}]\right),
\]
and, conditional on \(x\), generate independent \(Y_i\) with means
\[
\mathbb E[Y_i\mid x]
=x_i\left[(1-\gamma K)+\frac{\gamma u_i}{\sqrt{n-1}}\right],
\qquad 0\le\gamma K\le1.
\]
Writing
\[
S_K=\sum_i
\min\!\left(r_i^2,K\sqrt{n-1}|r_i|\right),
\]
random-sign completion of the clipped vector gives
\[
\boxed{
Q\ge
(1-\gamma K)^2q+
\frac{2\gamma(1-\gamma K)}{\sqrt{n-1}}S_K
-\gamma^2K^2Q.
}
\]
This harvests unstable or negative local-field mass. If all \(r_i\ge0\),
however, \(S_K\le Kq\) pointwise, so the first variation cannot improve
the witness. Coordinatewise local-field feedback therefore stops
exactly at a globally cut-stable positive-heavy configuration. A
successful successor must use a genuinely multi-vertex cut inequality.

### 3.23 Abstract flip certificates have no coercivity without cut
triangle identities

The edge-flip formula in §3.21 is exact, but its aggregate
majority-with-gap condition alone cannot force a useful lower bound on
\(M\). In an abstract antipodal state space containing only a vector
\(w\) and its negative \(-w\), the \(k\)-flip stability inequalities can
be satisfied even at the parity floor. Thus no argument using only
active-layer cardinality, antipodality, parity, and intersection sizes
can yield leading-order coercivity.

Any successful use of the certificate must exploit the defining
triangle relations of genuine cut vectors:
\[
v_{ij}v_{jk}v_{ki}=\text{the same global sign}
\quad\text{for every triangle }ijk.
\]
That test has now been completed in §3.24.

### 3.24 Local triangle rigidity reduces exactly to covering radius

Let \(H=(V,S)\) be a support graph, with \(c(H)\) components and cycle
rank
\[
\beta(H)=|S|-|V|+c(H).
\]
The restriction of the augmented cut code to \(S\) has dimension
\[
\dim(\mathcal C_n|_S)=
\begin{cases}
|V|-c(H),&H\text{ bipartite},\\
|V|-c(H)+1,&H\text{ nonbipartite},
\end{cases}
\]
and therefore codimension
\[
\boxed{
q(H)=
\begin{cases}
\beta(H),&H\text{ bipartite},\\
\beta(H)-1,&H\text{ nonbipartite}.
\end{cases}
}
\]
Indeed, ordinary cuts are the image of the binary vertex-edge incidence
map. The global complement bit adds the all-one edge vector, which
already lies in the cut space exactly when \(H\) is bipartite.

The dual local constraints are precisely even-cardinality Eulerian
subgraphs. Consequently every forest, matching, star, single triangle,
and connected odd unicyclic graph sees the full local cube. The first
genuine constraints occur on a \(4\)-cycle or between two independent
odd cycles. Moreover, every arbitrary local pattern is within
\(q(H)\) edge changes of a consistent augmented-cut pattern, by
pivoting a full-rank parity-check matrix. Bounded-cycle-excess tests
therefore have only bounded coercive power.

The global reduction is also exact. For an edge cochain
\(\alpha\in\mathbb F_2^E\), define
\[
(\delta\alpha)_{ijk}
=\alpha_{ij}+\alpha_{jk}+\alpha_{ki}.
\]
Then
\[
\mathcal C_n=\delta^{-1}(\langle\mathbf1\rangle),
\]
and the image of \(\delta\) is the space of two-graphs satisfying the
tetrahedron parity identity. Thus
\[
\rho(\mathcal C_n)
=
\max_{\tau\ {\rm two\text{-}graph}\bmod\mathbf1}
\min\{|\alpha|:\delta\alpha\in\{\tau,\tau+\mathbf1\}\}.
\]
The complete triangle system plus the all-flip majority certificate is
therefore exactly the original maximum cofilling/covering-radius
problem, not an extra regularity condition.

Rooting makes the obstruction transparent. After switching
\(a_{1i}=1\), the rooted triangle signs
\(\tau_{1ij}=a_{ij}\) are an arbitrary signing of \(K_{n-1}\), and
with \(x_1=1\) the energy is
\[
\sum_{i>1}x_i+
\sum_{2\le i<j}\tau_{1ij}x_ix_j.
\]
This is precisely the already unclosed affine recurrence.

The exact \(n=10\) audit reinforces the theorem: every triangle and
every \(4\)-cycle is certified by the top energy layer, whereas the
first top-layer failures occur on perfect matchings, which are forests
and carry no triangle constraint. The missing nonlocal statistic is the
conditional extension-gap profile
\[
\Gamma_H(y)=\min\{g_v:N_v|_{E(H)}=y\};
\]
triangle parity determines its domain but gives no control of its
values. The full proof is in `triangle_rigidity_reduction.md`.

### 3.25 Exact insertion profiles either fail closure or contain the
whole landscape

With the half-energy normalization
\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\]
the gauge-fixed insertion recurrence is
\[
\boxed{
M_{n+1}
=\min_A\max_x
\left(
|H_A(x)|+\left|\sum_i x_i\right|
\right).
}
\]
The complete signed magnetization-extrema profile
\[
U_A(m)=\max_{\sum x_i=m}H_A(x),
\qquad
L_A(m)=\min_{\sum x_i=m}H_A(x)
\]
does not close under choosing the next signed row. Two explicit
order-\(7\) minimizers have identical \((L_A(m),U_A(m))\) for every
magnetization \(m\), and indeed identical radial external-field
supports for every real field, but their best extension values are
\[
\boxed{E(A_1)=12,\qquad E(A_2)=10.}
\]
The two matrices and an exhaustive \(128\times128\) integer
reproducer are in `scale_transfer_profile_no_go.md`.

The exact closed state is the pair of external-field supports
\[
F_A^\sigma(h)=
\max_x\bigl(\sigma H_A(x)+h\cdot x\bigr),
\qquad \sigma=\pm1.
\]
For
\[
C=\begin{pmatrix}A&b\\b^T&0\end{pmatrix},
\]
it obeys
\[
F_C^\sigma(h,t)
=\max_{y=\pm1}
\left[
ty+F_A^\sigma(h+\sigma yb)
\right].
\]
But this closure has no compression: at every cube vertex \(x\),
\[
\boxed{
H_A(x)
=\lim_{\lambda\to\infty}
\left(F_A^+(\lambda x)-\lambda n\right),
\qquad
(F_A^+)^*(x)=-H_A(x).
}
\]
Thus the smallest evident exact insertion profile contains the entire
switching-energy word. A useful convergence proof needs a genuinely
lossy asymptotic regularity theorem for high-energy level sets, not
another exact finite recursion.

### 3.26 Nonuniform random-cut moments below \(n^{1/3}\) are
asymptotically tautological

Let \(D\) be switched so
\[
R=\mathbf1^TD\mathbf1=Q(D),
\]
choose independent cut indicators \(Z_i\sim{\rm Bernoulli}(p_i)\), put
\(v=1-2p\), and let \(C=C_D(Z)\). Then
\[
\boxed{
\mathbb EC=\frac{R-v^TDv}{4}.
}
\]
Hence \(0\le\mathbb EC\le R/2\) is exactly the existing cube inequality
\(|v^TDv|\le R\), even when each \(p_i\) is an arbitrary function of
the row sums.

Writing \(\eta=Z-p\), the exact Hoeffding decomposition is
\[
C-\mathbb EC
=(Dv)\cdot\eta
-2\sum_{i<j}D_{ij}\eta_i\eta_j.
\]
Orthogonality of the Bernoulli monomials and \(D_{ij}^2=1\) give
\[
\boxed{
\operatorname{Var}C
=
\sum_i p_i(1-p_i)(Dv)_i^2
+4\sum_{i<j}p_i(1-p_i)p_j(1-p_j).
}
\]
Therefore
\[
\operatorname{Var}C
\le\frac n4\|D\|_{\rm op}^2+\frac{n(n-1)}8.
\]
For a competing sequence \(R=c_n n^{3/2}\), \(c_n\ge c_0>0\), the
spectral bootstrap gives
\[
\frac{\operatorname{Var}C}{R^2}
\le
\frac14c_n^{-4/3}n^{-1/3}
+\frac18c_n^{-2}n^{-1}
=o(1).
\]
A Hanson--Wright/Khintchine bound upgrades this: every normalized
moment of order \(k=o(n^{1/3})\) collapses to the corresponding power
of the mean. Fixed-order cut moments and cycle statistics therefore
see only a point in \([0,1/2]\).

This loss is real. If the dense constraint is relaxed, take \(m\)
disjoint edges of weight
\[
w=\frac n{\sqrt{2m}},
\qquad m=\left\lfloor\frac{c^2n}{2}\right\rfloor.
\]
This weighted matching (not a sign matrix) satisfies
\[
Q=(c+o(1))n^{3/2},\qquad
\sum r_i^2=n^2,\qquad
\|W\|_F^2=n^2,
\]
all row sums are nonnegative, every cut lies in \([0,Q/2]\), and the
spectral bootstrap holds, for any fixed \(c>0\). Thus row moments,
cut positivity, and all \(o(n^{1/3})\)-order cut statistics cannot use
the essential condition \(|D_{ij}|=1\). The surviving target is a
growing-order endpoint statistic, such as the overlap-resolved number
of pairs of cuts within \(o(R)\) of \(R/2\).

### 3.27 Heavy positive rows can be manufactured at negligible leading
cost

Let \(D\) have order \(m\), switched so
\[
\mathbf1^TD\mathbf1=Q(D)=R.
\]
Adjoin \(k\) universally positive vertices:
\[
\widetilde D_k=
\begin{pmatrix}
J_k-I_k&J_{k,m}\\
J_{m,k}&D
\end{pmatrix}.
\]
Then the value is exact:
\[
\boxed{
Q(\widetilde D_k)=R+2km+k(k-1).
}
\]
For Boolean block sums \(s,t\), the energy is
\[
s^2-k+2st+y^TDy,
\]
whose absolute value is at most the displayed quantity, with equality
at the all-one vector.

Consequently \(k=o(\sqrt m)\) leaves the leading normalized value
unchanged. Already \(k=1\) creates a row of size \(m\), and changes the
row-square mass from \(S=\sum r_i^2\) to
\[
\boxed{
\widetilde S=m^2+S+2R+m,
}
\]
at energy cost only \(2m=o(m^{3/2})\). This is an actual dense sign
matrix, not a weighted relaxation. It proves that global cut
positivity, mean row sum, and a lower bound on the second row moment
cannot by themselves turn exceptional positive-heavy rows into a
leading-order gain.

For an exact conference switching, the full two-dimensional row-threshold
consequence can also be solved. Put
\[
U=\frac D{\sqrt{n-1}},\qquad
u=U\mathbf1,\qquad
\rho=\mathbb Eu.
\]
Then \(\mathbb Eu^2=1\), \(U\mathbf1=u\), \(Uu=\mathbf1\), and choosing
\(x_i=\operatorname{sign}(1-u_i)\) in the two-dimensional spectral
projection gives
\[
\boxed{
\mathbb E|1-u|\le\sqrt{1-\rho^2}.
}
\]
But among distributions with \(0\le u\le L\),
\(\mathbb Eu=\rho\), and \(\mathbb Eu^2=1\),
\[
\inf\mathbb E|1-u|
=(1-\rho)\left(1+\frac2L\right),
\]
attained on \(\{0,1,L\}\). Since \(L=\sqrt{n-1}\) is possible, the
threshold inequality is asymptotically vacuous.

The correct surviving statistic is the uniformly-integrated row-square
tail
\[
\Psi_n(K)=
\frac1n\sum_i
\left(\frac{r_i}{\sqrt n}\right)^2
\mathbf1_{\{r_i>K\sqrt n\}}.
\]
A useful continuation needs a regular-versus-peeling theorem: prove a
stronger result when \(\Psi_n\) is uniformly integrable, and show that
the exceptional vertices can otherwise be peeled or structurally
converted without losing the \(n^{3/2}\) objective. The complete proof
is in `multicut_heavy_field_note.md`.

### 3.28 Signed-Johnson endpoint hierarchy and the flat-Sidon no-go

Switch a signing \(D\) so
\[
R=\mathbf1^TD\mathbf1=Q(D),
\qquad
C_D(S)=\sum_{i\in S,j\notin S}d_{ij}
=\frac{R-q_D(\mathbf1_S)}4.
\]
For \(k\)-subsets define the signed Johnson operator
\[
(T_k)_{S,S-\{i\}+\{j\}}=d_{ij}.
\]
Restricting the Fourier convolution matrices of the nonnegative cube
polynomials \(R\pm q_D\) to level \(k\) gives the exact hierarchy
\[
\boxed{-\frac R2I\preceq T_k\preceq\frac R2I.}
\]
It also satisfies
\[
T_k\mathbf1(S)=C_D(S),
\qquad
\frac{\langle\mathbf1,T_k\mathbf1\rangle}{\binom nk}
=\frac{k(n-k)}{n(n-1)}R.
\]
For every Boolean \(x\), with
\(\psi_x(S)=\prod_{i\in S}x_i\),
\[
\frac{\langle\psi_x,T_k\psi_x\rangle}{\|\psi_x\|_2^2}
=\frac{k(n-k)}{n(n-1)}q_D(x).
\]

Entrywise flatness fixes only
\[
\binom nk^{-1}\operatorname{tr}T_k^2=k(n-k).
\]
Together with the spectral interval this yields merely
\[
R\ge2\sqrt{k(n-k)}=O(n),
\]
even at \(k\asymp n\). Thus Frobenius/trace information misses the
required \(n^{3/2}\) scale. A successful growing-order hierarchy must
use the correlated cycle structure inherited from the same base-edge
signing.

Likewise,
\[
\mathbb EC_D(S)=\frac{k(n-k)}{n(n-1)}R
\]
for a uniform \(k\)-set. Markov endpoint counts produce exponentially
many near-maximizers, but every resulting moment inequality is
scale-free because \(R\) cancels from both sides.

Finally, \(M_n\) is the minimum sup norm of a *flat* degree-two
tetrahedral Walsh polynomial. This is not the ordinary level-two Sidon
constant, which optimizes over arbitrary coefficient magnitudes.
Already at \(n=4\), the best flat ratio is
\[
\binom42/M_4=3/2,
\]
whereas the four-cycle polynomial with coefficients of magnitude
\(1/2\) has coefficient \(\ell_1\)-norm \(2\) and sup norm \(1\).
Ordinary Sidon asymptotics therefore do not determine the flat
minimum. Direct sums introduce a leading rectangular cross block, and
ordinary tensors raise degree two to degree four, so standard Sidon
stabilization does not give scale transfer.

### 3.29 Sharp insertion is endpoint-weighted discrepancy

With the half-energy normalization, define
\[
M=M(A)=\max_x|H_A(x)|,
\qquad
E(A)=\min_b\max_x\bigl(|H_A(x)|+|b\cdot x|\bigr).
\]
Then exactly
\[
E(A)-M
=
\min_b\max_x
\left[
|b\cdot x|-\bigl(M-|H_A(x)|\bigr)
\right].
\]
After switching a positive absolute maximizer to \(\mathbf1\), put
\[
C_D(S)=\sum_{i\in S,j\notin S}d_{ij}.
\]
Since
\[
H_D(y_S)=M-2C_D(S),
\qquad
0\le C_D(S)\le M,
\]
the recurrence becomes the exact weighted-discrepancy identity
\[
\boxed{
E(D)-M
=
\min_b\max_S
\left[
|b\cdot y_S|
-2\min\{C_D(S),M-C_D(S)\}
\right].
}
\]

There is a useful mixed, but not pure, insertion theorem. If
\[
d=M-|H_A(x)|
\]
and the local fields are oriented toward the sign of \(H_A(x)\), then
\[
-d/2\le\ell_j,
\qquad
\sum_j\ell_j=2(M-d),
\]
and hence
\[
\sum_j|\ell_j|\le2M+(n-2)d.
\]
Choosing a uniformly random duplicate row of \(A\), with its missing
diagonal entry filled by either sign, gives against every fixed state
\[
\boxed{
\mathbb E_j\bigl[
|H_A(x)|+|b^{(j)}\cdot x|
\bigr]
\le
M+\frac{2M}{n}+1-\frac{2d}{n}.
}
\]
The quantifier order is the entire gap: insertion requires one pure
row working for every endpoint cut. This isolates a purification
problem, not a scalar entropy estimate.

The exact ground-state obstruction is
\[
E(A)-M
\ge
\min_b\max_{\{|H_A(x)|=M\}}|b\cdot x|.
\]
A full orthogonal Boolean basis of extremizers forces overhead at least
\(\sqrt n\), although a nontrivial conference matrix cannot have such a
Boolean eigenbasis by a simple integrality argument. The complete note
is `sharp_insertion_weighted_discrepancy.md`.

### 3.30 Conditional regular peeling and replenishment towers

Let \(A_t\) be nested principal cores, \(q_t=Q(A_t)\), and at each
level switch a \(+q_t\) maximizer to \(\mathbf1\). Peel a set \(H_t\),
and define
\[
R_t=\sum_{i\in H_t}r_i,\qquad
h_t=\mathbf1_{H_t}^T(A_t)_{H_t,H_t}\mathbf1_{H_t},
\]
\[
d_t=q_t-q_{t+1}.
\]
If \(e_t\) is the energy in \(A_{t+1}\) of the old maximizer restricted
to the new core, put
\[
g_t=q_{t+1}-e_t\ge0.
\]
Direct deletion gives the exact identity
\[
\boxed{2R_t=d_t+h_t+g_t.}
\]

For disjoint peeled blocks,
\[
\sum_t|h_t|
\le\sum_tQ(A[H_t])
\le2q_0.
\]
The last inequality follows by choosing a ground state on every block,
retaining a common energy-sign class carrying at least half the total,
and randomizing the block-global signs to cancel cross terms. Also
\[
\sum_td_t=q_0-q_L\le q_0.
\]
Consequently,
\[
\boxed{
2\sum_tR_t\le3q_0+\sum_tg_t.
}
\]
If every peeled coordinate has \(r_i>K\sqrt n\), then
\[
\boxed{
\left|\bigcup_tH_t\right|
\le
\frac{3q_0+\sum_tg_t}{2K\sqrt n}.
}
\]
Thus \(\sum_tg_t=O(q_0)\) would make the heavy set \(O(n/K)\), and
\(K\to\infty\) would leave an \(n-o(n)\) principal regular core at no
leading cost.

The sole exact obstruction is now a **replenishment tower**:
successive deletions reveal substantially different maximizers and
\[
\sum_tg_t\gg Q(A_0).
\]
The spectral localization theorem in §1.15 has the same obstruction:
deleting one \(o(n)\)-vertex anchor may reveal a new high mode. A
convergence proof along this route must bound cumulative replenishment,
or assemble a large-gap tower into a Boolean witness exceeding
\(Q(A_0)\). The full checkpoint is in
`regular_peeling_replenishment.md`.

### 3.31 Falsified: a uniformly regular \(n-o(n)\) principal core

The following plausible strengthening of spectral localization is
false:

> \(Q(A_n)=O(n^{3/2})\) implies that deleting \(o(n)\) vertices leaves
> a principal matrix with operator norm \(O(\sqrt n)\), with one
> uniform implied constant.

For \(n=N^2\), partition the vertices into positive cliques. Let type
\(j\) occupy a fraction
\[
p_j\asymp j^{-3}
\]
of the vertices in blocks of size
\[
k_j=K_jN,\qquad K_j=aj,
\]
and put independent signs between distinct blocks. The total internal
quadratic norm is bounded by
\[
\sum_B|B|^2
=
\left(\sum_jp_jK_j+o(1)\right)n^{3/2}
=O(n^{3/2}),
\]
while a Hoeffding union bound gives a realization of the between-block
signing with quadratic norm \(O(n^{3/2})\).

For every fixed \(C\), choose a fixed type with \(K_j>2(C+1)\).
That type occupies \(p_jn+o(n)\) vertices. Any \(o(n)\) deletion leaves
one of its cliques with more than \((C+1)\sqrt n\) vertices, whose
principal Rayleigh quotient exceeds \(C\sqrt n\). Hence no
\(o(n)\)-deletion can produce a uniform \(O(\sqrt n)\) core.

The obstruction can be implanted into any \(O(n^{3/2})\) signing at an
arbitrarily small fixed normalized cost by overwriting only the clique
interiors. It cannot retain the obstruction at \(o(n^{3/2})\) cost;
this indicates that the correct optimizer-specific theorem must charge
spectral-tail mass to the excess above the asymptotic optimum.

There is an optimal qualitative replacement. Grothendieck--Pietsch
factorization and
\[
\|A\|_{\infty\to1}\le2Q(A)
\]
imply that for every fixed \(\varepsilon>0\), some principal set
\(U\) satisfies
\[
\boxed{
|U|\ge(1-\varepsilon)n,
\qquad
\|A[U]\|_{\rm op}
\le
\frac{4K_GQ(A)}{\varepsilon n}
=O_\varepsilon(\sqrt n).
}
\]
Thus regularization is always possible after a fixed fractional
deletion, but its constant must diverge as
\(\varepsilon\downarrow0\). The full construction and proof are in
`spectral_peeling_counterexample.md`.

### 3.32 Falsified: the sharp replenishment bound

The attractive conjecture
\[
\sum_tg_t\le2Q(A_0)
\]
for singleton suffix peeling is false.  An explicit \(15\times15\)
signing has
\[
Q(A_0)=62,\qquad
\max_{\text{ground-state ties}}\sum_tg_t=128>124=2Q(A_0).
\]
The suffix norms, from orders \(15\) down to \(1\), are
\[
62,58,52,44,42,30,28,28,18,14,12,8,6,2,0,
\]
and one maximizing tie choice gives replenishment gaps
\[
20,16,12,16,8,16,16,0,8,8,4,4,0,0.
\]
These values were independently re-enumerated over every Boolean state
of every suffix; their sum is \(128\).

There is an exact online-learning interpretation.  Encode the chosen
suffix ground state at time \(t\) as an augmented cut \(Z_t\), and form
the triangular mosaic \(W\) whose edge \((i,j)\), \(i<j\), is read from
the last leader before \(i\) is deleted.  Edgewise telescoping gives
\[
\boxed{\sum_tg_t=\langle A,W\rangle-Q(A)}
\]
with the corresponding ordered-matrix normalization.  Thus a constant
replenishment theorem is precisely a constant adaptivity-gap theorem
for suffix follow-the-leader over augmented cuts.  The example proves
that any such constant must exceed \(2\); it does not yet show that no
absolute constant exists.

### 3.33 Falsified: raw optimality excess controls regularity

Even an asymptotically minimizing sequence may have divergent
\(\|A\|_{\rm op}/\sqrt n\) and non-uniformly-integrable row squares.
Take an exact minimizer \(B_m\) along a liminf subsequence and adjoin
\(k=k_m\to\infty\), \(k=o(\sqrt m)\), universally positive vertices
after switching a positive maximizer of \(B_m\) to \(\mathbf1\).  Then
\[
\widetilde B_m=
\begin{pmatrix}
J_k-I_k&J_{k,m}\\
J_{m,k}&B_m
\end{pmatrix}
\]
satisfies the exact identity
\[
\boxed{
Q(\widetilde B_m)=Q(B_m)+2km+k(k-1).
}
\]
Hence its normalized value has the same liminf, while
\[
\frac{\|\widetilde B_m\|_{\rm op}}{\sqrt{m+k}}
\ge(1-o(1))\sqrt k\longrightarrow\infty.
\]
The exceptional set has \(k=o(n)\) vertices, so this kills only an
unpeeled excess-to-regularity theorem.  A post-peeling or
spectral-tail-charge theorem remains possible.

### 3.34 Proved: cumulative ground-state-frame visibility

For a block peeling step
\[
A_t=\begin{pmatrix}D_t&B_t\\B_t^\top&A_{t+1}\end{pmatrix},
\qquad
V_t=\max_{y\in\operatorname{GS}(A_{t+1})}\|B_ty\|_1,
\]
comparison with a ground state of the new core gives
\[
2V_t\le
Q(A_t)-Q(A_{t+1})+Q(D_t).
\]
For disjoint peeled blocks the decrements telescope, while the
same-sign-class/random-block-sign argument gives
\[
\sum_tQ(D_t)\le2Q(A_0).
\]
Therefore
\[
\boxed{\sum_tV_t\le\frac32Q(A_0).}
\]
Thus a long replenishment tower can exist only by repeatedly hiding its
cross modes from the exact ground-state frames of the successive cores.

### 3.35 Falsified: regular pointwise visibility inverse

Spectral regularity does not imply \(g_t\le C(K)V_t\) at a single step.
The order-nine suffix of the order-15 replenishment counterexample has
\[
Q(T)=Q(E)=28,\qquad g=16,
\]
but its order-eight core has exactly one absolute ground-state pair and
the deleted cross row is orthogonal to both, so \(V=0\).  Moreover the
integer matrix \(81I-4T^2\) is positive definite, hence
\[
\|T\|_{\rm op}<\frac92=\frac32\sqrt9.
\]
The exact matrix, Boolean witnesses, and Sylvester certificate are in
`optimality_excess_and_regularization.md`.  Any viable inverse must
group scales, use a near-ground-state layer, or charge a statistic other
than exact ground-state visibility.

---

## 4. Pending audit — do not cite as proved yet

### 4.1 Resolved: finite arcsine stability refinement

The earlier claimed quantitative remainder has now been audited and
proved with the stronger denominator
\[
8\pi(n-1)^{5/2}
\left(1-\frac1{n-1}\right)^{3/2}.
\]
Its proof and its factor-\(n\) scale limitation are recorded in
`orientation_even_stability_audit.md`; the stronger leading-scale
\(A^2\)-energy theorem is in §1.16.

### 4.2 External Claude campaign: useful leads, proof files unavailable

A separate campaign supplied a detailed summary, but its named technical
note and certificate bundles were not available in the shared file store.
The following are therefore research leads rather than accepted inputs:

1. **Haar half-subspace benchmark.** For a Haar-random rank-\(p/2\)
   projection \(P\), the claimed limit is
   \[
   \max_{x\in\{\pm1\}^p}\frac{\|Px\|^2}{p}
   \longrightarrow
   \beta_*=\frac12\left(1+\frac{\sqrt{15}}4\right).
   \]
   The union-bound side is immediate from the beta tail:
   \[
   \beta_*(1-\beta_*)=\frac1{64}.
   \]
   The asserted matching second-moment lower bound, especially the
   microscopic-overlap truncation, still requires reconstruction.
   In the quadratic-energy normalization this suggests
   \[
   c_{\rm Haar}=\frac{\sqrt{15}}8
   =0.4841229182\ldots .
   \]
2. **Paley resonance.** This claim has now been reconstructed and moved
   to the verified results in §1.12.
3. **Important logical caveat.** That limsup theorem alone does not prove
   that the Paley values fail to converge. The same summary labels
   convergence to \(\sqrt{15}/8\) on density-one nonresonant primes as a
   conjecture, and §1.12 now disproves the density-one part. A proved
   low-valued second subsequence is still missing.
4. **Entropic alternative.** The summary proposes that rare random
   signings might undercut \(c_{\rm Haar}\), governed by a lower-tail or
   Franz--Parisi rate. A speed-\(n^2\) rate can locate the onset of
   exponentially many good signings, but by itself cannot exclude a
   single algebraic switching orbit, whose probability already has the
   maximal cost \((\log2)n^2/2+o(n^2)\).

---

## 5. Current ten-route research cycle

The workspace permits three subagents concurrently, so routes run in waves.
A route survives only if it produces a scale-transfer inequality, a
proof-grade structural lemma, or a decisive obstruction.

| # | Route | Concrete target | Status |
|---:|---|---|---|
| 1 | Explicit conference AMP | Valid Onsager recursion reaching the spectral ceiling | Old claim falsified; paired recursion capped near \(0.705\) |
| 2 | Conference cut positivity | From all-cut positivity and \(C^2=(n-1)I\), prove \(R\ge(1-o(1))n^{3/2}\) | Heavy-row dichotomy proved; near-\(1\) conclusion still open |
| 3 | Spectral-moment dichotomy | Either direct Boolean witness \(\ge n^{3/2}\), or strong conference-like structure | Produced verified \(0.336493\) theorem; sharp dichotomy still open |
| 4 | Higher-order stability | Upgrade near-optimality to a two-eigenvalue/ETF approximation strong enough for rounding | Orientation-odd hierarchy obstructed by self-complementary signings |
| 5 | Near-ground-state entropy | Prove an \(O(\sqrt n)\) vertex-insertion inequality | Uniform theorem falsified at \(n=6\); asymptotic entropy lemma unknown |
| 6 | Correlated nonlocal lifts | Preserve scalar geometry under amplification without vector-relaxation loss | Flat-Fourier/orbit mechanism falsified |
| 7 | Scale-preserving free energy | Prove convergence on \(t=\beta/\sqrt n\), not merely fixed \(t\) | Fixed diagrams universal, but nonperturbative resonance kills all-\(\beta\) transfer |
| 8 | Signed coding dual | Build a signed high-degree certificate for the augmented cut-code covering radius | Exact dual is the original histogram; stopped |
| 9 | Optimizer mining | Discover or falsify scalable algebraic structure from exact/heuristic optima | Two-copy recurrence and Paley-minor hypothesis falsified |
| 10 | Energy–entropy compactness | Define a full profile closed under gluing and prove a unique asymptotic value | Fixed profiles/LDP shortcut obstructed by sparse planting |

Current successor routes:

- regular-versus-peeling control of the positive-heavy row-square tail;
- joint selection of high energy and high \(A^2\)-energy in the
  opposite-orientation rounding;
- growing-order endpoint/overlap statistics for the dense cut
  polynomial.

Recently completed successor routes:

- **Paley nonresonance:** stopped with a stronger negative theorem:
  every fixed admissible progression has Paley limsup \(1/2\).
- **Dependent local-field rounding:** exact \(A^2\)-energy and
  instability corrections proved; all coordinatewise variants stop at
  the cut-stable positive-heavy branch.
- **Triangle rigidity:** exact local-rank and two-graph reductions prove
  that bounded local triangle tests add no coercion beyond the original
  covering-radius problem.
- **Insertion profiles:** magnetization profiles fail closure, while
  the exact external-field closure is injective and contains the full
  energy word.
- **Nonuniform random cuts:** all \(o(n^{1/3})\)-order moments collapse
  to the cube-norm first moment and cannot see dense sign flatness.
- **Positive heavy fields:** universal-positive-vertex extensions show
  that exceptional heavy rows and an \(n^2\) second row moment can be
  created at only \(O(n)\) energy cost; first/second moments are
  insufficient without a peeling theorem.
- **Raw optimality-excess regularity:** falsified even for asymptotically
  minimizing sequences; the anomaly is always localized in the explicit
  construction.
- **Sharp replenishment coefficient \(2\):** falsified at order \(15\);
  the surviving question is an absolute or regularized adaptivity-gap
  bound.

---

## 6. Most useful sufficient lemmas

Any one of the following would settle convergence.

### 6.1 Uniform amplification

\[
\frac{M_{kn}}{(kn)^{3/2}}
\le
\frac{M_n}{n^{3/2}}+\varepsilon_n,
\qquad
\varepsilon_n\to0,
\]
uniformly in \(k\).

### 6.2 Sharp one-vertex insertion

\[
M_{n+1}
\le
M_n+\left(\frac32+o(1)\right)\frac{M_n}{n}.
\]

This would make \(M_n/n^{3/2}\) asymptotically nonincreasing.

### 6.3 Universal sharp lower bound

\[
Q(A)\ge(1-o(1))n^{3/2}
\quad\text{for every signing }A.
\]

Together with Paley conference upper bounds, this would prove
\[
\frac{M_n}{n^{3/2}}\to\frac12.
\]

### 6.4 Convergent finite-temperature minima

Define
\[
F_n(\beta)=
\frac1n\min_A
\log\mathbb E_x
\exp\!\left(\frac{\beta H_A(x)}{\sqrt n}\right).
\]

If \(F_n(\beta)\) converges for every fixed \(\beta>0\), uniformly enough to
send \(\beta\to\infty\), then
\[
\frac{M_n}{n^{3/2}}
\]
converges, because the zero-temperature approximation error is at most
\(\log2/\beta\).

The unresolved issue is convergence along the changing raw temperature
\(t=\beta/\sqrt n\).

---

## 7. Immediate next actions

1. Combine the smoothed field-plus-spin energy theorem with an
   orientation-even \(A^2\)-energy gain and the exact capped-field
   conversion, with constants uniform after fixed-\(K\) regularization.
2. Decide whether suffix augmented-cut follow-the-leader has an absolute
   adaptivity gap, or whether a scalable replenishment tower disproves it.
3. Prove a ground-state-frame visibility theorem for spectrally regular
   block decompositions; this would charge replenishment to \(Q\).
4. Build a closed scale-transfer profile, or prove that a proposed
   energy--magnetization profile is still insufficient.
5. Revisit Paley/nonabelian/random correlated upper constructions only
   if they include a moving \(p\)-dependent condition and a mechanism
   that controls all \(2^p\) Boolean vectors.
6. Revisit other correlated lifts only if they include a
   mechanism that controls zero-entropy resonant cube points.
7. Save a new checkpoint after every agent wave or any material
   proof/counterexample.

## 8. Proof acceptance policy

A result is moved to **Verified** only when:

1. all definitions and factors of \(2\) are explicit;
2. asymptotic error terms are uniform in the required matrix class;
3. any Gaussian/cavity step states the dependency structure and Onsager
   terms;
4. a second derivation, independent audit, or falsifying computation has
   been attempted;
5. the result advances the original convergence question, or is clearly
   labeled as a bound/obstruction only.
