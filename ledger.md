# Quadratic Signing Limit — Research Ledger

Last updated: 2026-07-29 (UTC)

## Purpose

This is the durable checkpoint for the problem

```math
M_n=\min_{a_{ij}\in\{\pm1\}}\max_{x_i\in\{\pm1\}}
\left|\sum_{1\le i<j\le n}a_{ij}x_ix_j\right|,
\qquad
\text{determine whether }\lim_{n\to\infty}\frac{M_n}{n^{3/2}}
\text{ exists.}
```

It separates:

- **Verified:** proof has been reconstructed and checked in the current work.
- **Pending audit:** plausible result from an earlier research wave, but its full proof has not yet been reconstructed.
- **Numerical:** computational evidence only.
- **Falsified:** a precise gap or counterexample has been found.
- **Open target:** a lemma that would materially advance or settle the problem.

## Notation

Let $`A`$ be the symmetric zero-diagonal matrix with off-diagonal entries
$`a_{ij}`$, and define

```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
\qquad
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top A x|.
```

Then

```math
x^\top A x=2H_A(x),\qquad Q(A)=2\max_x|H_A(x)|,
\qquad M_n=\frac12\min_A Q(A).
```

The conjecturally simplest outcome is

```math
\frac{M_n}{n^{3/2}}\longrightarrow \frac12,
```

equivalently $`Q(A)\ge(1-o(1))n^{3/2}`$ for every signing $`A`$.
This is not proved.

---

## 1. Verified results

### 1.1 Current rigorous asymptotic interval

```math
\boxed{
0.336493364431\ldots
\le
\liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\le
\limsup_{n\to\infty}\frac{M_n}{n^{3/2}}
\le
\frac12
}
```

The former verified lower constant $`1/\pi=0.318309886\ldots`$ is
superseded by the field-plus-spin theorem in §1.7.

### 1.2 Finite-$`n`$ Gaussian-rounding lower bound

For every $`n\ge2`$,

```math
\boxed{
M_n\ge
\frac{n(n-1)}{\pi}
\arcsin\!\frac1{\sqrt{n-1}}
}
```

and hence

```math
M_n\ge
\frac1\pi n^{3/2}
-\frac1{3\pi}\sqrt n
+O(n^{-1/2}).
```

Proof checkpoint:

1. Fix $`A`$, put $`m=n-1`$, and let
   $`\displaystyle P=\max_x H_A(x),\qquad Q=-\min_x H_A(x).`$
2. For $`s=\pm1`$, take $`g\sim N(0,I)`$ and set
   $`\displaystyle X^{(s)}=\operatorname{sgn}\bigl((sA+\sqrt m\,I)g\bigr).`$
3. The pre-sign coordinates have variance $`2m`$. For $`i\ne j`$, after
   multiplication by $`a_{ij}`$, their correlation is
   $`\displaystyle s\,a+b_{ij}, \qquad a=\frac1{\sqrt m}, \qquad b_{ij}=\frac{a_{ij}(A^2)_{ij}}{2m}.`$
4. The Gaussian arcsine identity gives
   $`\displaystyle P+Q\ge \frac2\pi\sum_{i<j} \left[ \arcsin(a+b_{ij})+\arcsin(a-b_{ij}) \right].`$
5. On the admissible domain,
   $`\displaystyle \frac{\arcsin(a+b)+\arcsin(a-b)}2\ge\arcsin a.`$
6. Since $`\max(P,Q)\ge(P+Q)/2`$, the displayed finite bound follows.

### 1.3 Conference-matrix upper bound

If $`C`$ is a symmetric conference matrix of order $`N`$, then

```math
C^2=(N-1)I
```

and therefore

```math
\max_x\left|\sum_{i<j}c_{ij}x_ix_j\right|
=\frac12\max_x|x^\top Cx|
\le\frac12N\sqrt{N-1}.
```

Symmetric Paley conference matrices exist at orders $`N=q+1`$ for prime
powers $`q\equiv1\pmod4`$. Primes in this progression can be chosen with
$`q=n+o(n)`$, and principal submatrices handle intermediate orders. Thus

```math
\limsup_{n\to\infty}\frac{M_n}{n^{3/2}}\le\frac12.
```

### 1.4 Exact augmented-cut-code formulation

Let $`N=\binom n2`$ and

```math
\mathcal C_n=
\left\{
(t+z_i+z_j)_{i<j}:t,z_i\in\mathbb F_2
\right\}.
```

This is the cut code of $`K_n`$ together with its complement. If
$`\rho(\mathcal C_n)`$ is its covering radius, then

```math
\boxed{M_n=N-2\rho(\mathcal C_n).}
```

So the problem asks whether the $`n^{3/2}`$-scale deficit of this covering
radius has a limit.

### 1.5 Elementary finite-size facts

```math
M_n\le M_{n+1}\le M_n+n.
```

- The first inequality follows by restricting an $`(n+1)`$-vertex signing
  and maximizing over the last spin.
- The second follows by extending an optimal $`n`$-vertex signing with any
  new signed row; the new linear term has absolute value at most $`n`$.

These bounds are too weak for convergence: $`O(n)`$ increments permit
$`O(1)`$ changes after normalization across windows of only $`O(\sqrt n)`$
vertices.

### 1.6 Exact block-gluing obstruction

For

```math
G=\begin{pmatrix}A&B\\B^\top&D\end{pmatrix},
```

maximizing over the relative global sign of the two blocks yields the exact
identity

```math
\max_{x,y}|H_A(x)+H_D(y)+x^\top By|
=
\max_{x,y}
\left(
|H_A(x)+H_D(y)|+|x^\top By|
\right).
```

Thus cross edges cannot cancel internal energy. A successful gluing theorem
must anti-align the large bilinear values of $`B`$ with the high-energy
layers of both blocks; a scalar inequality involving only $`M_n`$ cannot
express this.

### 1.7 Universal field-plus-spin lower bound

For every sequence of symmetric zero-diagonal sign matrices,

```math
\boxed{
\liminf_{n\to\infty}\frac{Q(A_n)}{n^{3/2}}
\ge c_*
=0.672986728863\ldots
}
```

and therefore

```math
\boxed{
\liminf_{n\to\infty}\frac{M_n}{n^{3/2}}
\ge\frac{c_*}{2}
=0.336493364431\ldots .
}
```

#### Spectral bootstrap

Let $`\lambda=\|A\|_{\mathrm{op}}`$, and choose a unit eigenvector $`v`$
whose eigenvalue has absolute value $`\lambda`$. The eigenvector equation
gives

```math
\lambda\|v\|_\infty\le\|v\|_1\le\sqrt n.
```

Set $`c=\|v\|_\infty^{-1}`$, and choose independent signs $`X_i`$ with
$`\mathbb E X_i=cv_i`$. Since $`A`$ has zero diagonal,

```math
\mathbb E X^\top AX=c^2v^\top Av,
```

so

```math
\boxed{Q(A)\ge\frac{\|A\|_{\mathrm{op}}^3}{n}.}
```

Thus, on a sequence with $`Q(A)=O(n^{3/2})`$,
```math
\|A\|_{\mathrm{op}}=O(n^{5/6}).
```
Writing $`m=n-1`$ and
```math
q_{ij}=\frac{(A^2)_{ij}}m,
```
we obtain
```math
\frac1{n(n-1)}\sum_{i\ne j}q_{ij}^2
\le
\frac{\operatorname{tr}A^4}{m^2n(n-1)}
\le
\frac{\|A\|_{\mathrm{op}}^2}{m^2}
=O(n^{-1/3}).
```
Hence distinct rows are asymptotically orthogonal in mean square.

#### Smoothed opposite-orientation rounding

Fix $`t,\tau>0`$, take independent Rademachers $`\xi_i`$ and Gaussians
$`Z_i`$, put
```math
h=\frac{A\xi}{\sqrt m},
\qquad
Y_i^\sigma
=\operatorname{sgn}(\sigma h_i+t\xi_i+\tau Z_i),
\quad \sigma\in\{\pm1\}.
```

For a pair $`i\ne j`$, isolate its direct coupling
$`\varepsilon=m^{-1/2}`$, integrate out $`Z`$, and Taylor-expand the
exact smoothed pair response before making any CLT approximation.
The two orientations have identical baselines, which cancel exactly.
The first-order response is

```math
a_{ij}\bigl(F_+(\varepsilon)-F_-(\varepsilon)\bigr)
=
4\varepsilon\,\mathbb E[\alpha_\tau(u)\beta_\tau(v)]
+O_{t,\tau}(\varepsilon^2),
```

where

```math
\psi_\tau(z)=2\Phi(z/\tau)-1,
```
```math
\alpha_\tau(u)
=\frac{\psi_\tau'(u+t)+\psi_\tau'(u-t)}2,
\qquad
\beta_\tau(v)
=\frac{\psi_\tau(v+t)-\psi_\tau(v-t)}2.
```

A smooth two-dimensional Lindeberg replacement has error
$`O_{t,\tau}(m^{-1/2})`$ after this Taylor extraction. Gaussian
covariance interpolation and the mean-square estimate for $`q_{ij}`$
then give, on average over pairs,

```math
\mathbb E[\alpha_\tau(u)\beta_\tau(v)]
=
2\phi_{1+\tau^2}(t)
\left[
2\Phi\!\left(\frac{t}{\sqrt{1+\tau^2}}\right)-1
\right]
+o(1).
```

Summing over ordered pairs and using
$`2Q(A)\ge\max x^\top Ax-\min x^\top Ax`$ yields

```math
\frac{Q(A)}{n\sqrt{n-1}}
\ge
4\phi_{1+\tau^2}(t)
\left[
2\Phi\!\left(\frac{t}{\sqrt{1+\tau^2}}\right)-1
\right]
-o(1).
```

Take $`n\to\infty`$ first, then $`\tau\downarrow0`$. The unique positive
optimizer $`t_*`$ solves

```math
2\phi(t_*)=t_*\bigl(2\Phi(t_*)-1\bigr),
\qquad
t_*=0.876902\ldots,
```

and

```math
c_*=4\phi(t_*)\bigl(2\Phi(t_*)-1\bigr)
=0.672986728863\ldots .
```

The order of operations is essential: a direct Berry–Esseen estimate on
the unexpanded pair response has leading-order error.

### 1.8 Local continuity in the order

For $`N=n+h`$, random cross edges give

```math
\boxed{
M_{n+h}\le M_n+M_h+
\sqrt{2nh(n+h+2)\log2}.
}
```

Together with monotonicity, this implies, uniformly for $`h=o(n)`$,

```math
\boxed{
\frac{M_{n+h}}{(n+h)^{3/2}}
-
\frac{M_n}{n^{3/2}}
=o(1).
}
```

Consequently, if $`d_{k+1}/d_k\to1`$ and
$`M_{d_k}/d_k^{3/2}\to c`$, then the full sequence converges to $`c`$.
The missing fact is convergence of the minima on any such ratio-dense
subsequence.

### 1.9 Conference heavy-row dichotomy

Let $`D`$ be a switching of a symmetric conference matrix for which
```math
R=\mathbf1^\top D\mathbf1=\max_x x^\top Cx,
\qquad r=D\mathbf1,
\qquad m=n-1.
```
Then $`r_i>0`$, $`\|r\|_2^2=nm`$, and $`Dr=m\mathbf1`$.

For any $`p\in[0,1]^n`$, put each vertex independently in a random set
$`S`$ with probabilities $`p_i`$. All-cut positivity gives
```math
\mathbb E\,c_D(S,S^c)
=r\cdot p-p^\top Dp
\le \frac R2.
```
Taking $`p=tr`$ gives the exact inequality
```math
mt(n-tR)\le\frac R2
\qquad
\left(0\le t\le\frac1{r_{\max}}\right).
```

Consequently:

- If $`r_{\max}\le 2R/n`$, then $`t=n/(2R)`$ is feasible and
  $`\displaystyle \boxed{R\ge n\sqrt{\frac{n-1}{2}}.}`$
- In general, writing
  $`\displaystyle h=\frac{r_{\max}}{\sqrt{n-1}}, \qquad \rho=\frac{R}{n\sqrt{n-1}},`$
  the choice $`t=1/r_{\max}`$ gives
  $`\displaystyle \boxed{\rho\ge\frac{2h}{h^2+2}.}`$

Thus any conference sequence with $`\rho<1/\sqrt2-o(1)`$ must have a
genuinely heavy row. Converting that heavy row into a stronger Boolean
witness remains open.

For a conference matrix of order $`n=4\ell+2`$, all row sums in any
switching are congruent modulo $`4`$. At a maximizing switching they are
positive odd integers in one common residue class. This supplies useful
arithmetic rigidity but has not yet closed the heavy-row case.

### 1.10 Audited small-order insertion obstruction

Exhaustive enumeration after fixing the first row by switching gives:

```math
M_5=4,\qquad M_6=5.
```

For every gauge-class optimizer found:

- at $`n=5`$, $`20`$ of the $`32`$ Boolean vectors are exact extremizers;
- at $`n=6`$, $`24`$ of the $`64`$ Boolean vectors are exact extremizers.

Define the best one-vertex extension profile
```math
E(A)=
\min_{b\in\{\pm1\}^n}
\max_{x\in\{\pm1\}^n}
\left(|H_A(x)|+|b\cdot x|\right).
```
All twelve gauge-fixed order-$`6`$ minimizers satisfy
```math
\boxed{E(A)=9.}
```
Thus their best one-vertex extension jumps by $`4`$, whereas the
derivative-scale target $`3M_6/(2\cdot6)`$ is only $`1.25`$. This does not
disprove an asymptotic insertion theorem, but it rules out the hoped-for
uniform finite theorem and shows that low entropy of the extremal layer is
not automatic.

### 1.11 Reduction to primes $`1\bmod 4`$

Let
```math
L_n=\frac{M_n}{n^{3/2}}.
```
The full limit exists if and only if $`L_p`$ converges as $`p\to\infty`$
through primes $`p\equiv1\pmod4`$, and the two limits then agree.

Indeed, the prime number theorem in arithmetic progressions implies that
the consecutive primes in this progression have ratio tending to $`1`$.
For every $`n`$, choose consecutive such primes
```math
p_-(n)\le n\le p_+(n).
```
Then $`p_\pm(n)/n\to1`$, while monotonicity gives
```math
M_{p_-(n)}\le M_n\le M_{p_+(n)}.
```
After division by $`n^{3/2}`$, convergence on the prime subsequence
squeezes the full sequence. The converse is immediate.

This is a useful localization of the problem, but it does not itself
compare different prime orders.

### 1.12 Paley square-wave resonance gives the spectral limsup

Let $`p\equiv1\pmod4`$ be prime and let
```math
S_{jk}=\chi_p(j-k),\qquad j,k\in\mathbb F_p,
```
with zero diagonal. In the unitary Fourier normalization,
```math
x^\top Sx
=
p^{-1/2}\sum_{m\ne0}\chi_p(m)|\widehat x(m)|^2.
```

For the square wave
```math
x_j=\operatorname{sgn}\cos(2\pi j/p),
```
its unnormalized Fourier coefficients satisfy, on symmetric
frequencies,
```math
\widehat x(m)
=
\frac{(-1)^{(m-1)/2}}
{\sin(\pi m/(2p))}
\quad(m\ {\rm odd}),
```
while the total even-frequency contribution is asymptotically
negligible. Consequently, for each fixed odd $`m`$, the pair of
frequencies $`\pm m`$ carries asymptotic Fourier mass
```math
\frac{8}{\pi^2m^2}.
```

Fix $`L`$. Choose a residue class
```math
p\equiv1
\pmod{8\prod_{\substack{\ell\le2L+1\\\ell\ {\rm odd\ prime}}}}\, .
```
Quadratic reciprocity makes every positive odd
$`m\le2L+1`$ a quadratic residue modulo $`p`$, and Dirichlet's theorem
supplies infinitely many primes in this class. Since
```math
\frac8{\pi^2}\sum_{m\ {\rm odd}>0}\frac1{m^2}=1,
```
first taking $`p\to\infty`$ in the class and then $`L\to\infty`$ gives
```math
\frac{x^\top Sx}{p^{3/2}}\to1.
```
Adding the extra row and column to form the Paley conference matrix
changes the quadratic energy by only $`O(p)`$. Hence, in the original
half-quadratic normalization,
```math
\boxed{
\limsup_{\substack{p\to\infty\\p\equiv1(4)}}
\frac{\max_x|H_{\rm Paley,p}(x)|}{p^{3/2}}
=\frac12.
}
```

This theorem shows that exact spectral flatness does not prevent sparse
arithmetic resonance. It does **not** show that the Paley values fail to
converge; a rigorously controlled nonresonant subsequence is still
missing.

The resonance is not a zero-density curiosity. The exact square-wave
formula is
```math
\frac{x^\top Sx}{p^{3/2}}
=
\frac8{\pi^2}
\sum_{\substack{h\ge1\\h\ {\rm odd}}}
\frac{\chi_p(h)}{h^2}
+O(p^{-1}).
```
If every odd $`h\le25`$ is a quadratic residue, then even an adversarial
tail gives
```math
\frac{x^\top Sx}{p^{3/2}}
\ge
\frac{16}{\pi^2}
\sum_{\substack{h\le25\\h\ {\rm odd}}}\frac1{h^2}
-1-o(1)
=0.9688395921\ldots-o(1),
```
strictly above the Haar doubled benchmark
$`\sqrt{15}/4=0.9682458366\ldots`$.
It suffices to force
```math
3,5,7,11,13,17,19,23
```
to be quadratic residues. Among primes $`p\equiv1\pmod4`$, these
conditions have relative Dirichlet density $`2^{-8}=1/256`$.
Therefore the external claim that Paley values converge to the Haar
constant on a density-one set of primes is false. A useful upper-bound
route needs only a ratio-dense *good* prime subsequence, but that still
requires an all-Boolean nonresonance theorem.

There is now a stronger obstruction. Fix any admissible arithmetic
progression of primes $`p\equiv1\pmod4`$, or equivalently any finite
compatible prescription of Legendre symbols together with congruence
conditions. For every $`\delta>0`$, a fixed refined progression inside
it has Boolean Paley witnesses with doubled normalized energy at least
$`1-\delta`$. Consequently,
```math
\boxed{
\limsup_{\substack{p\to\infty\\p\ \mathrm{in\ any\ fixed\
admissible\ progression}}}
\frac{\max_x|x^TA_px|}{2p^{3/2}}
=\frac12.
}
```

The construction is explicit. For each prescribed nonresidue prime
$`\ell`$, choose $`u\in\{\pm1\}^{\ell}`$ with $`\sum u=1`$, write
```math
j=\sum_{d=0}^{2r}j_d\ell^d,
\qquad
v_j=\prod_{d\ {\rm even}}u_{j_d}
\quad
\left(j\bmod\ell^{2r+1}\right).
```
Its normalized DFT vanishes whenever the frequency has odd
$`\ell`$-adic valuation, while its DC mass is only
$`\ell^{-2r-2}`$. CRT-tensor these gadgets over the finitely many
prescribed nonresidues. The resulting circle step function has
arbitrarily close to all of its Fourier mass on frequencies where
every bad-prime valuation is even. After prescribing the finitely many
remaining small prime symbols to be positive, Dirichlet's theorem and
Riemann sampling transfer the witness to $`\mathbb F_p`$.

This proves that no fixed congruence class, and no finite-character
inverse theorem, can yield a Paley upper bound below $`1/2`$. A
hypothetical low-valued Paley subsequence must use an increasing,
$`p`$-dependent amount of arithmetic information. The full proof and
normalization audit are in `paley_resonance_gadget.md`.

### 1.13 Exact opposite-orientation $`A^2`$-energy theorem

For every signing $`A`$, every $`t,\tau`$, and the smoothed witnesses
```math
X_i^\sigma
=\operatorname{sgn}\!\left(
\sigma\frac{(A\xi)_i}{\sqrt{n-1}}+t\xi_i+\tau Z_i
\right),
\qquad \sigma\in\{\pm1\},
```
one has the exact finite-$`n`$ inequality
```math
\boxed{
\frac12\sum_{\sigma=\pm1}
\mathbb E\|AX^\sigma\|_2^2
\ge n(n-1).
}
```

To prove it, expand
```math
\frac12\sum_\sigma\mathbb E\|AX^\sigma\|_2^2
=\operatorname{tr}A^2+
\sum_{i\ne j}(A^2)_{ij}C_{ij},
\qquad
C_{ij}=\frac12\sum_\sigma
\mathbb E[X_i^\sigma X_j^\sigma].
```
For a fixed pair $`i,j`$, switch coordinates so $`a_{ij}=1`$ and the
common coefficients in row $`i`$ are all $`+1`$. Write the corresponding
row-$`j`$ coefficients as $`c_k\in\{\pm1\}`$. Then
```math
d=\sum_{k\ne i,j}c_k=(A^2)_{ij}
```
in this gauge. Conditioning on every variable except one common
Rademacher shows that $`C_{ij}`$ is nondecreasing in each $`c_k`$: the
change from an opposite to an equal coupling is a product of two
nonnegative smoothed-threshold increments. A measure-preserving
transformation that flips all common Rademachers, $`\sigma`$, the
$`j`$-spin, and the $`j`$-dither shows
```math
C(-c)=-C(c).
```
Hence $`C(c)`$ has the sign of $`\sum c_k=d`$, so every off-diagonal
summand satisfies
```math
(A^2)_{ij}C_{ij}\ge0.
```
Switching preserves this product, completing the proof.

This supplies the full conference-scale squared local field for the
two opposite-orientation witnesses **on average**, without a CLT,
pseudorthogonality, or a near-minimizer hypothesis. The remaining
selection problem is to combine high oriented energy and high
$`A^2`$-energy in one witness and then exploit the only surviving case:
cut-stable, positive, heavy local fields.

### 1.14 Improved spectral bootstrap by asymmetric Boolean rounding

For every symmetric zero-diagonal matrix $`A`$, not only sign matrices,
```math
\boxed{Q(A)\ge\frac12\|A\|_{\rm op}^2.}
```

First let
```math
B(A)=\max_{x,y\in\{\pm1\}^n}|x^TAy|.
```
Polarization writes $`x^TAy`$ as the difference of the quadratic forms
of the two disjointly supported vectors
$`(x+y)/2,(x-y)/2\in\{0,\pm1\}^n`$. Randomly completing every zero
coordinate to a sign shows that the absolute quadratic form of either
partial sign vector is at most $`Q(A)`$. Hence
```math
B(A)\le2Q(A).
```

Choose a unit eigenvector $`v`$ with eigenvalue
$`|\lambda|=\|A\|_{\rm op}`$, take $`x=\operatorname{sign}v`$, and choose
an independent Boolean vector $`Y`$ with
```math
\mathbb EY_i=\frac{v_i}{\|v\|_\infty}.
```
Then
```math
B(A)\ge
\left|\mathbb E\,x^TAY\right|
=|\lambda|\,\frac{\|v\|_1}{\|v\|_\infty}.
```
At a coordinate where $`|v_i|=\|v\|_\infty`$, the eigenvector equation
and $`|a_{ij}|\le1`$ give
```math
|\lambda|\|v\|_\infty\le\|v\|_1.
```
Combining the last three displays proves the claim.

Thus any competing sequence with $`Q(A)=O(n^{3/2})`$ satisfies
```math
\boxed{\|A\|_{\rm op}=O(n^{3/4}),}
```
improving the $`O(n^{5/6})`$ bootstrap in §1.7. Consequently, for
$`q_{ij}=(A^2)_{ij}/(n-1)`$,
```math
\frac1{n(n-1)}\sum_{i\ne j}q_{ij}^2
\le\frac{\|A\|_{\rm op}^2}{(n-1)^2}
=O(n^{-1/2}).
```
The remaining gap is qualitative: $`O(n^{3/4})`$ is still much larger
than the $`O(\sqrt n)`$ operator-norm regime where squared local fields
are uniformly $`O(n^2)`$ for every Boolean witness.

### 1.15 Joint selection, localized spectral anomalies, and capped
profiles

For the opposite-orientation field-plus-spin law, orient the sampled
energy as
```math
R=\sigma X^TAX,\qquad S=X^TA^2X.
```
If a competing sequence has
```math
Q=(c_*+o(1))n^{3/2},
\qquad
\mathbb ER\ge(c_*-o(1))n^{3/2},
```
then $`R\le Q`$ pointwise implies
```math
\mathbb E(Q-R)=o(n^{3/2}).
```
A near-$`Q`$ sample automatically has
```math
S\ge R^2/n=(c_*^2-o(1))n^2
```
and negligible negative local-field mass by the correction in §3.22.
The exact average $`\mathbb ES\ge n(n-1)`$ upgrades this to
```math
S\ge(1-o(1))n^2
```
for the same near-maximal sample whenever $`S/n^2`$ is uniformly
integrable. In particular this holds if
```math
\|A\|_{\rm op}=O(\sqrt n).
```

Every failure of this regular condition is localized. If
$`Av=\lambda v`$, $`\|v\|_2=1`$, and $`Q=Q(A)`$, then for every
$`\theta>1`$, with
```math
s^2=\frac{\lambda}{\theta Q},
\qquad
T=\{i:|v_i|>s\},
```
one has
```math
\boxed{
|T|\le\frac{\theta Q}{\lambda},
\qquad
\|v_T\|_2^2\ge\frac{1-\theta^{-1}}3.
}
```
Indeed, writing $`v=u+w`$ on $`T\cup T^c`$,
```math
w^TAw=\lambda(1-2\|u\|_2^2)+u^TAu
\ge\lambda(1-3\|u\|_2^2),
```
whereas $`w/s\in[-1,1]^n`$ gives
```math
|w^TAw|\le Qs^2=\lambda/\theta.
```
Thus a mode $`\lambda=L_n\sqrt n`$, $`L_n\to\infty`$, has a fixed
fraction of its mass on $`O(n/L_n)=o(n)`$ vertices.

There is a quantitative theorem after a successful peeling/capping
step. Switch a near-maximal witness to $`\mathbf1`$, write
```math
r=D\mathbf1,\quad
c=\frac{\mathbf1^TD\mathbf1}{n^{3/2}},\quad
s=\frac{\|r\|_2^2}{n^2},
```
and assume
```math
\|D\|_{\rm op}\le K\sqrt n,\qquad
\max_i|r_i|\le H\sqrt n,
```
with negligible negative-field mass. The spectral measure of
$`\mathbf1`$ has first three normalized moments $`c,s,t`$ on
$`[-K,K]`$. From
```math
(K-z)(z-a)^2\ge0
```
and optimization in $`a`$,
```math
t\le
M_K(c,s)
:=
Ks-\frac{(s-Kc)^2}{K-c}.
```
Applying all-cut positivity to independent probabilities
$`p_i=\alpha(r_i)_+/\sqrt n`$ gives
```math
\boxed{
\frac c2\ge
\max_{0\le\alpha\le1/H}
\left[\alpha s-\alpha^2M_K(c,s)\right],
\qquad
c\ge\frac{s}{H}.
}
```
For the coefficient-one, flat-spectrum case $`s=K=1`$, a sequence
saturating $`c=c_*`$ must therefore retain a local field at least
```math
\boxed{(1.941916296\ldots-o(1))\sqrt n.}
```
The improved cap remains nontrivial only up to
```math
K<1.0220798875\ldots.
```
This reduces the lower-bound problem to peeling localized heavy
coordinates without allowing a succession of unrelated maximizing
states. The detailed proofs are in
`joint_selection_and_spectral_localization.md`.

### 1.16 Exact orientation-even $`A^2`$-energy gain

The unsmoothed opposite-orientation Gaussian witnesses from §1.2 obey
a stronger finite identity than §1.13 detects. Put
```math
m=n-1,\qquad
X^\sigma=\operatorname{sign}\bigl((\sigma A+\sqrt m\,I)g\bigr).
```
Then
```math
\boxed{
\frac12\sum_{\sigma=\pm1}
\mathbb E\,(X^\sigma)^TA^2X^\sigma
\ge
nm+
\frac{\|A^2-mI\|_F^2}
{\pi m\sqrt{1-1/m}}
}
\qquad(m>1).
```

For $`i\ne j`$, write
```math
c=(A^2)_{ij},\qquad
u=\frac{c}{2m},\qquad
v=\frac{a_{ij}}{\sqrt m}.
```
The pair correlation averaged over the two orientations is
```math
C_{ij}
=
\frac1\pi
\left[
\arcsin(u+v)+\arcsin(u-v)
\right].
```
If
```math
G_v(u)=\arcsin(u+v)+\arcsin(u-v),
```
then $`G_v(0)=0`$, and convexity of
```math
g(z)=(1-z^2)^{-1/2}
```
gives
```math
G_v'(u)=g(u+v)+g(u-v)
\ge2g(v)
=\frac2{\sqrt{1-1/m}}.
```
The inequality reverses in the correct way when integrating to
negative $`u`$, so
```math
cC_{ij}
\ge
\frac{c^2}{\pi m\sqrt{1-1/m}}.
```
Summing the off-diagonal terms and adding
$`\operatorname{tr}A^2=nm`$ proves the theorem.

This is a leading-scale orientation-even statistic:
a Wigner-scale defect
```math
\|A^2-mI\|_F^2=\Theta(n^3)
```
raises the averaged squared local field by $`\Theta(n^2)`$. The
remaining task is joint selection with high oriented energy for these
same witnesses, followed by control of heavy positive coordinates.

### 1.17 Exact capped-field conversion

The joint-selection loss in §1.16 can be isolated in one explicit tail
term.  For a Boolean witness $`x`$, orient $`A`$ so that
```math
q=x^TAx,\qquad r_i=x_i(Ax)_i,
```
and put
```math
S_K=\sum_i r_i\,\operatorname{clip}
\bigl(r_i,[-K\sqrt m,K\sqrt m]\bigr),\qquad m=n-1.
```
Switching $`x`$ to $`\mathbf1`$, define
```math
u_i=\operatorname{clip}
\bigl(r_i,[-K\sqrt m,K\sqrt m]\bigr),\qquad
\mu=(1-\alpha)\mathbf1+\frac{\alpha u}{K\sqrt m}.
```
For $`0\le\alpha\le1`$, $`\mu\in[-1,1]^n`$.  Multilinear Boolean
rounding and $`|u^TAu|\le K^2mQ(A)`$ give the exact inequality
```math
\boxed{
(1+\alpha^2)Q(A)
\ge
(1-\alpha)^2q+
\frac{2\alpha(1-\alpha)}{K\sqrt m}S_K .
}
```
This inequality may be averaged over a witness distribution without
selecting one sample having both large $`q`$ and large $`S_K`$.

Write
```math
C=\frac{Q(A)}{n\sqrt m},\qquad
c=\frac{\mathbb E q}{n\sqrt m},\qquad
z=\frac{\mathbb E S_K}{Knm}.
```
Optimization in $`\alpha`$ yields
```math
\boxed{
C\ge{\cal F}(c,z),\qquad
{\cal F}(c,z)=
\begin{cases}
c,&z\le c,\$$1mm]
c-z+\sqrt{z^2+(z-c)^2},&z>c.
\end{cases}
}
```
For the opposite-orientation Gaussian witnesses of §§1.2 and 1.16,
```math
c\ge\frac2\pi+o(1).
```
If
```math
\delta_n=\frac{\|A^2-mI\|_F^2}{nm^2}
```
and
```math
\Psi_n(K)=\frac1{nm}\,
\mathbb E\sum_i r_i^2
\mathbf1_{\{|r_i|>K\sqrt m\}},
```
then §1.16 gives
```math
\boxed{
C\ge
{\cal F}\left(
\frac2\pi+o(1),
\frac{1+\delta_n/(\pi\sqrt{1-1/m})-\Psi_n(K)}K
\right).
}
```
Thus the orientation-even defect already converts into a stronger
quadratic witness whenever a fixed cap retains enough squared local
field.  The remaining loss is exactly the positive heavy-field tail;
universal-positive-vertex examples show that it cannot be bounded from
the first two moments alone.

### 1.18 Regularized asymptotic near-minimizers

Raw near-minimizers need not be spectrally regular, but regular
near-minimizers can always be constructed with an explicit two-limit
tradeoff.  Let
```math
q_n=\min_AQ(A).
```
For every $`K\ge1`$, Grothendieck--Pietsch deletion followed by a
conference-type refill and random cross signing gives a full
order-$`n`$ signing $`A'_{n,K}`$ satisfying
```math
\boxed{
Q(A'_{n,K})
\le q_n+O(K^{-1/2}n^{3/2}),\qquad
\|A'_{n,K}\|_{\rm op}=O(K\sqrt n).
}
```
Indeed, delete at most $`n/K`$ vertices so that the retained principal
matrix has norm $`O(K\sqrt n)`$.  The refill has internal norm
$`O((n/K)^{3/2})`$, and a random cross block has Boolean bilinear norm
$`O(n^{3/2}/\sqrt K)`$ and operator norm $`O(\sqrt n)`$.

Letting $`K\to\infty`$ arbitrarily slowly produces spectrally controlled
asymptotic near-minimizers.  This does not by itself control local-field
tails uniformly in $`K`$, but it legitimizes proving lower bounds first
for fixed $`K`$ and then tracking the dependence on $`K`$.

There is also an exact block visibility inequality.  For
```math
A=\begin{pmatrix}D&B\\B^T&E\end{pmatrix}
```
and any Boolean $`y`$,
```math
\boxed{
Q(A)\ge |y^TEy|+\sqrt2\,\|By\|_2.
}
```
If $`y`$ is an oriented $`Q(E)`$-ground state and
$`x=\operatorname{sign}(By)`$, then the sharper deterministic form is
```math
\boxed{
Q(A)-Q(E)
\ge2\max_{y\in{\rm GS}(E)}\|By\|_1-Q(D).
}
```
Consequently, a replenishment tower can persist only if successive
cross-block singular directions avoid every ground-state frame of the
retained cores.

### 1.19 Relative-invariance defect gain for the same field-plus witnesses

**Independently audited.**  The derivative-product replacement,
orientation factors, parity anchors, and exact small-order pair laws
have all been checked.

The Gaussian first-Hermite defect coefficient transfers to the actual
Rademacher row-field witnesses.  Fix $`t\in\mathbb R,\tau>0`$, put
```math
X_i^\sigma=\operatorname{sign}\left(
\sigma\frac{(A\xi)_i}{\sqrt{n-1}}+t\xi_i+\tau Z_i
\right),
```
and let
```math
f_{t,\tau}(u)
=\frac{\psi_\tau(u+t)+\psi_\tau(u-t)}2,\qquad
K_{t,\tau}(q)=\mathbb E f_{t,\tau}(G)f_{t,\tau}(H).
```
For a pair, switch its endpoints so that $`a_{ij}=1`$.  In the original
gauge put
```math
d=a_{ij}(A^2)_{ij}.
```
Its opposite-orientation correlation $`C_m(d)`$, $`m=n-1`$, satisfies
uniformly
```math
\boxed{
\left|C_m(d)-K_{t,\tau}(d/m)\right|
\le C_{t,\tau}\frac{|d|}{m^{3/2}}.
}
```
The proof telescopes in the number of common row signs.  A
one-coordinate change factors exactly as a product of two threshold
increments.  Each increment and its first three derivatives is
$`O_{t,\tau}(m^{-1/2})`$, so a two-dimensional Lindeberg replacement
of the remaining $`m-2`$ common Rademachers costs
$`O_{t,\tau}(m^{-3/2})`$ per step.  The Gaussian step is
```math
\frac2mK_{t,\tau}'(d/m)+O_{t,\tau}(m^{-3/2}).
```
Oddness supplies the parity anchor at $`d=0`$, or across
$`-1\to1`$, and summing the errors remains relative even for
$`|d|=\Theta(m)`$.

Since
```math
qK_{t,\tau}(q)
\ge4\phi_{1+\tau^2}(t)^2q^2,
```
one obtains
```math
\boxed{
\frac12\sum_{\sigma=\pm1}
\mathbb E (X^\sigma)^TA^2X^\sigma
\ge nm+
\left(
4\phi_{1+\tau^2}(t)^2-O_{t,\tau}(m^{-1/2})
\right)
\frac{\|A^2-mI\|_F^2}{m}.
}
```
This uses the same witnesses whose oriented energy is asymptotically
```math
e(t,\tau)n\sqrt m,\qquad
e(t,\tau)=4\phi_{1+\tau^2}(t)
\left[2\Phi\!\left(t/\sqrt{1+\tau^2}\right)-1\right].
```
Consequently, with
```math
\delta=\frac{\|A^2-mI\|_F^2}{nm^2},
```
the capped-field conversion becomes
```math
\boxed{
\frac{Q(A)}{n\sqrt m}
\ge{\cal F}\left(
e(t,\tau)-o(1),
\frac{
1+\left(4\phi_{1+\tau^2}(t)^2-o(1)\right)\delta
-\Psi_{t,\tau,n}(H)}
H
\right).
}
```
At $`t=t_*=0.8769009856\ldots`$, then $`\tau\downarrow0`$, the defect
coefficient is $`0.2950713629\ldots`$.  A tail-free sequence saturating
$`c_*=0.672986728863\ldots`$ must therefore have
```math
H\ge
\frac{1+0.2950713629\ldots\,\delta}{c_*}.
```
For $`\delta=1`$, this is $`H\ge1.92436\ldots`$.  The only loss left in
this scalar theorem is the sparse heavy-positive-field tail; the
negative-field update controls the opposite tail but not this one.
Full details and the finite monotonicity fallback are in
`orientation_even_stability_audit.md`.

This theorem has now been independently reconstructed.  The audit
checked the switching gauge, exact $`1/4`$ increment, the uniform
two-dimensional Lindeberg estimate even at singular covariance,
direct-edge and variance corrections, parity telescoping, and the
ordered-pair Frobenius normalization.  Exact small-order enumeration
also verifies the pair increment and full matrix identity.  The
certificate is `relative_lindeberg_independent_audit.md`.

### 1.20 Conditional-independence $`2\to4`$ tail theorem

For the same smoothed field-plus witnesses, condition on
$`(\xi,\sigma)`$ and put
```math
\mu_j=\psi_\tau\left(
\sigma(A\xi)_j/\sqrt m+t\xi_j
\right),\qquad b=A\mu,\qquad m=n-1.
```
The dithers make the output coordinates conditionally independent, so
```math
(AX^\sigma)_i=b_i+\eta_i,
```
where $`\eta_i`$ is a centered sum of $`m`$ independent variables of
range length $`2`$.  If $`\|A\|_{\rm op}\le K\sqrt m`$, then for
$`H\ge1`$,
```math
\boxed{
\Psi_{t,\tau,n}(H)
\le
\frac{20}{H^2nm^2}\mathbb E_\xi\|A\mu(\xi)\|_4^4
+(4K^2+H^2+8)e^{-H^2/8}.
}
```

The remaining $`2\to4`$ quantity is uniformly bounded for fixed
$`(K,t,\tau)`$.  With
```math
L_1=\|\psi_\tau'\|_\infty,\qquad
L_2=\|\psi_\tau''\|_\infty,
```
a pointwise flip-influence calculation gives
```math
\sum_k(\Delta_kb_i)^2
\le
\left[8K^2L_1^2+2(2+2L_2)^2\right]n
=:\Gamma_{K,t,\tau}n.
```
Hypercube log-Sobolev and oddness of $`b_i`$ then yield
```math
\boxed{
\frac{\mathbb E\|A\mu\|_4^4}{nm^2}
\le C\Gamma_{K,t,\tau}^2.
}
```
Thus squared local fields are uniformly integrable on every fixed
spectrally regular, fixed-dither branch.  Quantitatively this does not
close the problem: $`L_1\asymp\tau^{-1}`$,
$`L_2\asymp\tau^{-2}`$, and the product-noise term is already too
large at the decisive $`H\approx1.9`$.  It is therefore a rigorous
tail localization theorem, but not the fixed-cap estimate needed in
the scalar criterion of §1.19.  Full constants and proof are in
`orientation_even_stability_audit.md`.

---

## 2. Important exact reformulations

### 2.1 Switching and signed cuts

For a signing $`A`$, replacing $`A`$ by
$`\operatorname{diag}(x)A\operatorname{diag}(x)`$ is Seidel switching.
The value $`x^\top A x`$ is the total excess after switching.

For a conference matrix $`C`$, choose a switching $`D`$ for which
```math
R=\mathbf1^\top D\mathbf1=\max_x x^\top Cx.
```
If $`c_D(S,S^c)`$ denotes the signed sum across the cut, global maximality
implies
```math
\boxed{
0\le c_D(S,S^c)\le R/2
\quad\text{for every }S\subseteq[n].
}
```

Also, with $`r=D\mathbf1`$,
```math
D^2=(n-1)I,\qquad
\|r\|_2^2=n(n-1),\qquad
Dr=(n-1)\mathbf1.
```

A strong lower bound on $`R=\sum_i r_i`$ from these cut constraints would
be a direct route to the conference case.

### 2.2 Projection form for conference matrices

For
```math
U=\frac{C}{\sqrt{n-1}},\qquad
P=\frac{I+U}{2},
```
we have $`U^2=I`$ and $`P`$ is a rank-$`n/2`$ projection. Moreover,
```math
\frac{x^\top Cx}{n\sqrt{n-1}}
=
\frac{x^\top Ux}{n}
=
\frac{2\|Px\|_2^2}{n}-1.
```

Thus the spectral ceiling $`1`$ is equivalent to finding a Boolean vector
at $`o(\sqrt n)`$ Euclidean distance from one eigenspace.

---

## 3. Falsified or stopped routes

### 3.1 Falsified: naive repeated conference AMP reaches $`1`$

The earlier claim that repeated projection has state evolution
```math
0.67299\to0.68512\to0.69539\to\cdots\to1
```
is false for a fixed conference matrix.

Let $`U=C/\sqrt{n-1}`$, so $`U^2=I`$. Start with
```math
G\sim N(0,I+U),\qquad X_0=\operatorname{sgn}G,
```
and put
```math
a=\sqrt{2/\pi},\qquad \sigma^2=1-a^2,
\qquad
R=\frac{X_0-aG}{\sigma},\qquad W=UR.
```

The first update is
```math
X_1=\operatorname{sgn}(X_0+UX_0).
```

In the scalar limit, write
```math
Y=\operatorname{sgn}\bigl(\operatorname{sgn}Z+aZ+\sigma W\bigr),
\quad Z,W\stackrel{\mathrm{iid}}{\sim}N(0,1).
```
Define
```math
b=\mathbb E[ZY],\qquad
c=\mathbb E[WY],\qquad
d=\mathbb E[RY],
\quad
R=\frac{\operatorname{sgn}Z-aZ}{\sigma}.
```
Numerically,
```math
b=0.7920592175,\quad
c=0.05415383275,\quad
d=0.5726603594.
```

The decisive algebraic obstruction is
```math
UW=R.
```
If
```math
\eta=X_1-bG-cW,
```
then $`\langle \eta,R\rangle\to d\ne0`$, and hence
```math
\langle W,U\eta\rangle
=\langle UW,\eta\rangle
=\langle R,\eta\rangle
\to d.
```
Therefore $`U\eta`$ is not a fresh Gaussian residual. The omitted
backtracking term is order one.

The corrected paired prediction after one step is
```math
b^2+2cd=0.68938131\ldots,
```
not the naive
```math
b^2+cd=0.65836956\ldots.
```

A legitimate recursion must retain Gram–Schmidt pairs
$`(R_s,W_s)`$ satisfying
```math
UR_s=W_s,\qquad UW_s=R_s.
```
Population recursion then gives approximately
```math
0.63662,\ 0.68938,\ 0.69937,\ 0.70254,\ 0.70377,\ldots
```
and appears to converge near $`0.7054`$, not $`1`$.

### 3.2 Numerical confirmation of the AMP obstruction

For direct iteration
```math
x\mapsto\operatorname{sgn}((I+U)x)
```
on Paley conference matrices, terminal ratios
```math
\frac{|x^\top Cx|}{n\sqrt{n-1}}
```
cluster around $`0.70`$, and fixed points usually appear after $`3`$–$`5`$
steps.

Independent current-work computations:

- Orders from $`102`$ through $`402`$: typical terminal means approximately
  $`0.69`$–$`0.73`$.
- Single-spin greedy local search does substantially better, around
  $`0.91`$–$`0.94`$ at orders $`102,194,402`$, but this is still numerical
  and does not prove approach to $`1`$.

### 3.3 Falsified/stopped: scalar Fekete argument

Splitting into comparable blocks introduces a rectangular minimax term of
order $`(n+m)^{3/2}`$, exactly the leading scale. Ordinary subadditivity
therefore does not imply convergence.

At free-energy scale $`t=\beta/\sqrt n`$, splitting $`2n`$ vertices changes
the internal-block inverse temperature from $`\beta`$ to
$`\beta/\sqrt2`$. Fixed-temperature Fekete subadditivity sees the wrong
diagonal scaling.

### 3.4 Stopped: universal local tensor/Hadamard gadget

Common Hadamard lifts immediately relax scalar Boolean choices into vector
choices. On the order-$`6`$ conference seed, the relaxed profile already
jumps away from the scalar optimum. Fixed local block templates cannot
provide the required scale-preserving amplification.

### 3.5 Stopped: ordinary spin-glass interpolation

After random vertex switching, the pairwise-overlap covariance is the same
for every signing. Yet the Boolean optimum can range from order
$`n^{3/2}`$ to order $`n^2`$. A method using only that covariance is
information-theoretically unable to distinguish the relevant structures.

### 3.6 Stopped: spectral certificate alone

The operator-norm bound gives the upper construction constant $`1/2`$, but
standard SDP/vector relaxations lose the scalar Boolean geometry and have a
hard spectral floor. They do not prove the matching universal lower bound.

### 3.7 Stopped: $`n=5\to10`$ two-copy optimizer lift

The $`n=10`$ optimizer can be represented as a signed two-copy lift of an
$`n=5`$ optimizer, explaining its repeated spectrum, but iterating the lift
causes the normalized objective to grow:
```math
0.716\to0.822\to1.073
```
for the corresponding $`Q/(n\sqrt n)`$ sequence. It is not an asymptotic
minimizing construction.

### 3.8 Falsified: tensor submultiplicativity

There are explicit symmetric zero-diagonal $`5\times5`$ signings

```math
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
```

with
```math
Q(A)=Q(B)=12,
```
but a Boolean witness gives
```math
\boxed{Q(A\otimes B)\ge176>144=Q(A)Q(B).}
```

Thus the equal-modulus quadratic norm is not tensor-submultiplicative.
Ordinary Sidon-constant tensor machinery leaves the full-support,
equal-modulus subclass and cannot supply the missing scale transfer.

### 3.9 Falsified: generic centered-free-energy axioms suffice

Fixed-raw-temperature centered subadditivity, evenness, real analyticity,
convexity, the correct variance, and a correct zero-temperature slope do
not force convergence on $`t=\beta/\sqrt n`$.

A concrete countermodel is obtained by choosing an oscillating
```math
a_n=c_n n^{3/2},
\qquad
c_n=1+\varepsilon\sin(\log\log(n+n_0)),
```
with $`a_n/\binom n2`$ decreasing, and defining
```math
F_n(t)=\log2+
\frac{a_n^2}{\binom n2}
\log\cosh\!\left(\frac{\binom n2}{a_n}t\right).
```

Then $`F_n`$ is even, analytic and convex,
```math
F_n''(0)=\binom n2,
\qquad
\lim_{t\to\infty}\frac{F_n(t)}t=a_n,
```
and the centered functions
```math
G_n(t)=F_n(t)-\binom n2\log\cosh t
```
are subadditive. Nevertheless,
```math
\frac{F_n(\beta/\sqrt n)}n
=
2c_n^2\log\cosh\!\left(\frac{\beta}{2c_n}\right)+o(1)
```
oscillates. A successful free-energy proof must use discrete
higher-cumulant identities specific to quadratic Ising energies.

### 3.10 Stopped: ground-state entropy alone gives sharp insertion

For fixed $`A`$, a random new signed row $`b`$ can achieve extension
overhead at most $`t`$ if a weighted energy-layer sum of the form
```math
\sum_x
\exp\!\left(
-\frac{(M-|H_A(x)|+t)^2}{2n}
\right)
<\frac12.
```
An $`O(\sqrt n)`$ conclusion from this criterion requires only $`O(1)`$
effective entropy near the ground-state layer. The audited $`n=5,6`$
optimizers in §1.10 have respectively $`20/32`$ and $`24/64`$ exact
extremizers, so this premise fails badly even at minimizing signings.

After switching a positive maximizer, the relevant slack layer is a family
of low signed cuts. Ordinary Karger cut-counting does not apply because a
signed cut function is not submodular. No asymptotic replacement is known.

### 3.11 Stopped: orientation-odd higher-cycle stability

A higher-Hermite refinement of field-plus-spin rounding produces signed
odd-cycle corrections. In normalized notation $`B=A/\sqrt{n-1}`$ and
$`q_{ij}=(B^2)_{ij}`$, its correction has the schematic form
```math
S_t(B)=
\sum_{i\ne j}B_{ij}K_t(q_{ij}),
\qquad
K_t(q)=\sum_{\ell\ {\rm odd}}k_\ell(t)^2q^\ell.
```
The leading term is proportional to $`\operatorname{tr}(B^3)`$; the
remaining terms are longer signed theta-cycle profiles.

This entire family has a decisive null class. If
```math
PBP^\top=-B
```
for a permutation matrix $`P`$, then $`PB^2P^\top=B^2`$ and hence
```math
S_t(B)=-S_t(B)=0
```
for every $`t`$. Self-complementary signings of this kind can be
arbitrarily far from conference structure (random examples have
Wigner-like fourth moment). Therefore no orientation-odd, one-channel
higher-cycle correction can force $`B^2\approx I`$. A viable successor
must be orientation-even and detect
$`\|B^2-I\|_F^2`$ or equivalent information.

### 3.12 Exact-computation checkpoint at order $`10`$

A complete mixed-integer optimization gives
```math
M_{10}=13,\qquad Q_{10}=26.
```
One optimizer has spectrum
```math
\pm3.933464\ldots\ (\text{twice}),\quad
\pm2.554969\ldots\ (\text{twice}),\quad
\pm1,
```
and $`A^2-9I`$ is supported on two signed $`5`$-cycles with nonzero
entries $`\pm4`$. This explains the earlier observed repeated spectral
polynomial, but its natural two-copy continuation was already falsified
as an asymptotic minimizing lift.

### 3.13 Optimality ceiling for one-probe coordinatewise rounding

The field-plus-spin constant $`c_*`$ in §1.7 is optimal within the full
class of one-probe coordinatewise Boolean rules.

Let $`G\sim N(0,1)`$, $`S`$ be an independent sign, and
$`f(G,S)\in\{\pm1\}`$. Put
```math
a=\mathbb E[Gf(G,S)],
\qquad
b=\mathbb E[Sf(G,S)].
```
The direct-edge response of the associated one-probe rounding is
$`2ab`$. Write
```math
g(z)=\frac{f(z,+1)+f(z,-1)}2,
\qquad
h(z)=\frac{f(z,+1)-f(z,-1)}2.
```
Pointwise,
```math
|g(z)|+|h(z)|\le1.
```
After choosing signs so $`a,b\ge0`$, rearrangement shows that for fixed
$`\mathbb E h(G)`$, the product is maximized by using the spin on the
smallest values of $`|G|`$ and the field on the largest values:
```math
f(z,s)=\operatorname{sgn}(z+ts).
```
For this threshold rule,
```math
a=2\phi(t),
\qquad
b=2\Phi(t)-1.
```
Therefore
```math
\sup_f 2ab
=
\max_{t\ge0}
4\phi(t)(2\Phi(t)-1)
=c_*.
```
Any improvement over $`0.3364933644\ldots`$ for $`M_n/n^{3/2}`$ must
use multiple dependent probes, a genuinely nonlocal rule, or additional
structure of $`A`$.

### 3.14 Stopped: correlated flat-Fourier lifts

For an abelian fiber group $`G`$ of order $`k`$, a broad class of
correlated block lifts has flat Fourier kernels
```math
\widehat b_{ij}(\chi)=\sqrt{k}\,a_{ij}^{(\chi)}.
```
If
```math
m_{i,\chi}=\frac1k\sum_{g\in G}x_i(g)\chi(g),
\qquad
\sum_\chi |m_{i,\chi}|^2=1,
```
then the cross energy is exactly
```math
k^{3/2}
\sum_\chi\sum_{i<j}
a_{ij}^{(\chi)}
m_{i,\chi}m_{j,\chi}.
```
Thus the lift exposes a row-sphere, multi-channel relaxation.

If the channel signings form a balanced full signed-permutation orbit of
a seed $`A`$, this relaxation is exactly the spectral value
```math
\frac n2\|A\|_{\rm op}.
```
Top eigenvectors of the conjugates attain the lower bound after their
squared coordinates are averaged; the spectral inequality gives the
reverse bound. Hence the natural orbit-symmetrized lift necessarily
jumps to the spectral ceiling and cannot preserve a scalar optimum below
$`1/2`$.

A finite audit at fiber order $`4`$ gives the same obstruction: all
$`768`$ Fourier-compatible quadruples of order-$`5`$ minimizers have
lifted cross value $`36`$ or $`40`$, while exact preservation would
require $`32`$. This computational datum should be treated as an audited
finite no-go, not as a theorem covering every nonabelian lift.

### 3.15 Stopped: radial signed cut-code dual certificates

Let $`C=\mathcal C_n`$, $`N=\binom n2`$, and
```math
D=C^\perp
=
\{\text{even-cardinality Eulerian edge sets of }K_n\}.
```
For a signing $`A`$, define
```math
T_w(A)=
\sum_{\substack{F\in D\\|F|=w}}
\prod_{e\in F}a_e.
```
The complete signed dual enumerator satisfies the exact identity
```math
\sum_wT_w(A)z^w
=
(\cosh\beta)^{-N}2^{-n}
\sum_x\cosh(\beta q_A(x)),
\qquad z=\tanh\beta,
```
where $`q_A(x)=\sum_{i<j}a_{ij}x_ix_j`$.

Equivalently, the coset weight enumerator is
```math
\sum_{c\in C}u^{d(y,c)}
=
\frac12\sum_x
\left[
u^{(N-q_A(x))/2}
+
u^{(N+q_A(x))/2}
\right].
```
Thus the full signed Eulerian/Krawtchouk data is invertibly identical to
the original switching-energy histogram.

For a ball indicator, the exact signed certificate uses
```math
\sum_{j\le R}K_j^N(w)=K_R^{N-1}(w-1),
```
but demanding its positivity for every translate is term-for-term
equivalent to the original covering-radius statement. Fixed dual degree
sees only fixed moments of the energy distribution; conference
sequences have the same fixed Gaussian-chaos limits even when sparse
resonant vectors change their maxima. A resonance-sensitive certificate
must therefore have genuinely growing degree, where it again becomes
the full energy-histogram problem.

### 3.16 Scalar Eulerian-pressure axioms still permit oscillation

The correct two-sided finite-temperature object is
```math
\Gamma_n(\rho)
=
\min_A\log
\left[
(\cosh t)^{-\binom n2}
2^{-n}\sum_x\cosh(tH_A(x))
\right],
\qquad
\rho=\tanh t.
```
Exact random-edge arguments show that $`\Gamma_n`$ is:

- subadditive in $`n`$ at fixed $`\rho`$;
- nonincreasing in $`n`$;
- nonincreasing in $`\rho`$.

The tempting scaling inequality
```math
\Gamma_n(\lambda\rho)\le\lambda^2\Gamma_n(\rho)
```
is false already at $`n=4`$, where
```math
\Gamma_4(\rho)=\log(1-\rho^4).
```

Moreover, all the listed scalar properties, even together with even
analyticity and a first signing-dependent term of order $`\rho^4`$,
allow diagonal oscillation. For a slowly oscillating positive $`c_n`$
chosen so $`\theta_n=c_n/\sqrt n`$ decreases, the abstract family
```math
\Gamma_n(\rho)
=
-n\,\frac{(\rho/\theta_n)^4}{1+(\rho/\theta_n)^4}
```
has all those properties, while
```math
\frac1n\Gamma_n(\beta/\sqrt n)
=
-\frac{(\beta/c_n)^4}{1+(\beta/c_n)^4}
```
oscillates. A free-energy proof must use coefficient-level Eulerian
constraints, an overlap hierarchy, or another discrete feature—not only
the scalar pressure inequalities.

### 3.17 Fixed conference diagrams are universal, but resonance is
nonperturbative

For a conference matrix put
```math
U=\frac C{\sqrt{n-1}},
\qquad U^2=I,
\qquad |U_{ij}|=(n-1)^{-1/2}.
```
In the linked-cluster expansion of
```math
\log\mathbb E_x
\exp\!\left(\frac\beta2x^\top Ux\right),
```
every spin-index vertex has even degree. Repeatedly contracting
degree-$`2`$ vertices by $`U^2=I`$ leaves either a fully contractible
cycle/cactus diagram or a core of minimum degree at least $`4`$. For each
fixed core with $`r`$ vertices and $`e`$ edges,
```math
e\ge2r,
```
so its tensor sum is $`O(n^r n^{-e/2})=O(1)`$. Only the contractible
diagrams contribute $`O(n)`$, and their values are fixed by $`U^2=I`$
and flatness.

Therefore every fixed Taylor coefficient of the pressure divided by
$`n`$ is switching-independent and agrees with the Haar
half-involution/Random-Orthogonal-Model coefficient.

The conclusion fails nonperturbatively. If a conference family has a
Boolean $`+1`$ eigenvector, then its normalized pressure is at least
```math
\frac\beta2-\log2+o(1).
```
For the Haar annealed pressure
```math
p_H(\beta)
=
\sup_{|u|<1}
\left\{
\frac{\beta u}{2}
+\frac14\log(1-u^2)
\right\},
```
one has $`p_H(8)=3.22233\ldots`$, whereas the Boolean eigenvector gives
$`4-\log2=3.30685\ldots`$. Thus a zero-entropy resonant state changes
the pressure by $`O(n)`$ while being invisible to every fixed diagram.
Fixed-order perturbation theory cannot justify the zero-temperature ROM
transfer.

### 3.18 Stopped: speed-$`n^2`$ disorder LDP as a shortcut

Let
```math
K_n(T)=
\#\left\{
A:\max_x|H_A(x)|\le T
\right\},
\qquad
p_n(c)=2^{-\binom n2}K_n(cn^{3/2}).
```
If the feasible set is empty, then
```math
-n^{-2}\log p_n(c)=+\infty.
```
If it is nonempty, switching and global sign already produce at least
$`2^{n-1}`$ feasible signings, so
```math
-\frac1{n^2}\log p_n(c)
\le
\frac{\binom n2-n+1}{n^2}\log2
\longrightarrow\frac{\log2}{2}.
```
Consequently, if $`c`$ lies strictly between the liminf and limsup of
$`M_n/n^{3/2}`$, the proposed rate has an infinite subsequence and a
bounded subsequence. Proving its extended-real limit for every $`c`$
would already prove the original convergence by support; the LDP is not
an easier preliminary.

There is also a finite-profile planting obstruction. Start with a
low-norm signing $`A`$, choose a Boolean vector $`s`$ with
$`|H_A(s)|\le\sqrt{\binom n2}`$, and flip
```math
h=\lfloor\delta n^{3/2}\rfloor
```
edges on which $`a_{ij}s_is_j=-1`$. The new signing $`B`$ satisfies
```math
H_B(s)=2\delta n^{3/2}+o(n^{3/2}),
```
so its maximum changes at leading scale. But for a uniform spin $`X`$,
```math
\mathbb E_X\bigl(H_B(X)-H_A(X)\bigr)^2=4h.
```
Thus for every fixed replica number $`k`$, the normalized joint
$`k`$-replica energy laws of $`A`$ and $`B`$ have Wasserstein distance
```math
O_k(n^{-3/4}),
```
and every fixed signed-subgraph density also agrees asymptotically.
Fixed-replica and bounded-degree profiles cannot detect a planted
zero-entropy extremizer.

### 3.19 Stopped: deterministic Haar transfer from flatness alone

The Haar half-subspace constant cannot be transferred from any of the
usual deterministic isotropy assumptions.

1. Block copies of $`C_6/\sqrt5`$ give constant-diagonal projections
   whose maximum projection fraction is
   $`\displaystyle \frac12\left(1+\frac{10}{6\sqrt5}\right) =0.872678\ldots<\beta_*=0.984122\ldots .`$
   Degeneracy allows an eigenbasis to be mixed across blocks so that
   its coherence is $`O(n^{-1/2})`$. Thus constant diagonal plus an
   incoherent eigenbasis is insufficient.
2. A balanced Hadamard eigenspace has constant diagonal and a perfectly
   flat eigenbasis but contains Boolean columns, so its maximum
   projection fraction is $`1`$.
3. Most strongly, the exact Paley conference projections in §1.12 have
   constant diagonal $`1/2`$ and every off-diagonal magnitude exactly
   $`1/(2\sqrt p)`$, yet the resonant subsequence has Boolean projection
   fraction tending to $`1`$.

Any Haar comparison must control anti-alignment with all $`2^n`$ cube
points or explicitly classify arithmetic resonance; flatness,
equiangularity, and delocalization do not suffice.

### 3.20 Independent multi-probe rounding has the same exact ceiling

Let
```math
G\sim N(0,I_k),
\qquad
S\sim{\rm Unif}\{\pm1\}^k,
```
independently, and let $`f(G,S)\in[-1,1]`$ be arbitrary. Define
```math
a=\mathbb E[Gf(G,S)],
\qquad
b=\mathbb E[Sf(G,S)].
```
The first-order direct-edge response of the associated independent
multi-probe rounding is $`2a\cdot b`$.

Put
```math
m(s)=\mathbb E_Gf(G,s),
\qquad
r^2=\mathbb E_Sm(S)^2.
```
The coordinate functions $`S_j`$ are orthonormal, so Bessel's inequality
gives
```math
\|b\|_2\le r.
```
For fixed $`s`$, the sharp Gaussian centroid inequality gives
```math
\left\|\mathbb E_G[Gf(G,s)]\right\|_2\le J(m(s)),
```
where
```math
J(u)=
2\phi\!\left(
\Phi^{-1}\!\left(\frac{1+u}{2}\right)
\right),
\qquad 0\le u\le1.
```
The function $`K(v)=J(\sqrt v)`$ is concave on $`[0,1]`$. With
$`u=2\Phi(z)-1`$, this is equivalent to
```math
2\Phi(z)-1\ge2z\phi(z),
```
whose two sides agree at $`z=0`$ and whose difference has derivative
$`2z^2\phi(z)\ge0`$. Hence
```math
\|a\|_2
\le\mathbb E J(m(S))
\le J(r).
```
It follows that
```math
2a\cdot b
\le2rJ(r)
\le\max_{0\le r\le1}2rJ(r)
=c_*.
```
Equality is attained by a one-spin threshold/dictatorship
```math
f(G,S)=\operatorname{sgn}(G_1+t_*S_1).
```
Thus arbitrarily many fresh independent channels, arbitrary response
boundaries, and their mixtures cannot improve the verified
$`0.3364933644\ldots`$ lower bound. An improvement must use dependent
matrix probes and retain their backtracking/Onsager structure.

### 3.21 Exact edge-flip/deep-hole certificate

For the augmented cut vectors
```math
\mathcal V=\{\pm(x_ix_j)_{i<j}:x\in\{\pm1\}^n\},
```
write
```math
h_v=a\cdot v,\qquad
M=\max_{v\in\mathcal V}h_v,
\qquad
g_v=\frac{M-h_v}{2},
```
and
```math
N_v=\{e:a_ev_e=-1\}.
```
If $`a^S`$ is obtained by flipping the edge set $`S`$, then
```math
\boxed{
M(a^S)
=
M-2\min_{v\in\mathcal V}
\left[
g_v+|S|-2|S\cap N_v|
\right].
}
```
Therefore $`a`$ is stable under every flip set of size at most $`k`$ if
and only if, for every such $`S`$, some $`v`$ satisfies
```math
\boxed{
2|S\cap N_v|\ge |S|+g_v.
}
```
Only energy layers $`g_v\le|S|`$ can certify a $`k`$-edge flip. For a
single edge, the negative supports from the top two energy layers
$`g=0,1`$ must cover every coordinate.

For the exact order-$`10`$ optimizer, the $`g=0`$ layer has $`40`$
gauge-fixed states with $`|N_v|=16`$, and the $`g=1`$ layer has $`80`$
with $`|N_v|=17`$. Every set of at most four edges is certified by the
active layer. Thirteen five-edge sets fail active-only certification;
they are perfect matchings, and the $`g=1`$ layer rescues all of them at
equality.

Sparse-flip stability is not sufficient for global optimality. The
following order-$`11`$ signing has $`M=19`$ and is stable under every
set of at most four edge flips:
```math
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
```
The flip set
```math
\{(0,3),(0,7),(1,6),(2,7),(5,9)\}
```
lowers its value to $`17`$, and an independent order-$`11`$ signing
with value $`17`$ is
```math
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
```
Exhaustive verification over
$`\sum_{j=1}^4\binom{55}{j}`$ flip sets gives minimum post-flip values
$`21,19,19,19`$ for sizes $`1,2,3,4`$, respectively. Thus even a
fourth-order local certificate can be trapped above the true minimum.

### 3.22 Dependent negative-field correction and its stable branch

Fix a signing $`A`$, a Boolean vector $`x`$, and write
```math
q=x^TAx,\qquad
r_i=x_i(Ax)_i,\qquad
L_-=\sum_i(-r_i)_+,\qquad Q=Q(A).
```
Let $`B=\{i:r_i<0\}`$, and independently flip each coordinate in $`B`$
with probability $`p`$. For the resulting Boolean vector $`y`$,
```math
\mathbb E\,y^TAy
=q+4pL_-+4p^2x_B^TA_{BB}x_B.
```
Randomly completing $`x_B`$ outside $`B`$ shows
```math
|x_B^TA_{BB}x_B|\le Q,
```
and hence
```math
Q\ge q+\max_{0\le p\le1}
\left(4pL_--4p^2Q\right).
```
In the branch $`L_-\le2Q`$, the optimizer is
$`p=L_-/(2Q)`$, so
```math
\boxed{\quad
Q\ge q+\frac{L_-^2}{Q},
\qquad
L_-^2\le Q(Q-q).
\quad}
```
Thus any witness with $`q=Q-o(n^{3/2})`$ that cannot be improved by this
dependent second stage has $`L_-=o(n^{3/2})`$.

There is also an exact quantitative description of this stable branch.
Put
```math
R=\sum_i r_i=q,\qquad
H=\max_i(r_i)_+,\qquad
S=\sum_i r_i^2=x^TA^2x.
```
Since $`r_i^2\le Hr_i`$ on the positive fields and
$`r_i^2\le(n-1)(-r_i)`$ on the negative fields,
```math
\boxed{\quad
S\le H(R+L_-)+(n-1)L_-,
\qquad
H\ge\frac{S-(n-1)L_-}{R+L_-}.
\quad}
```
Consequently, if $`q\sim c n^{3/2}`$, $`Q-q=o(n^{3/2})`$, and this
adaptive witness additionally satisfies $`S=(1-o(1))n^2`$, then either
the negative-field correction improves the leading constant or
```math
H\ge\left(\frac1c-o(1)\right)\sqrt n.
```
At $`c=c_*=0.672986728863\ldots`$, this is
$`H>1.485\sqrt n`$, more than $`2.20`$ times the average field
$`R/n`$.

Section 1.13 has now replaced the proposed four-point estimate with the
stronger exact opposite-orientation average
```math
\frac12\sum_\sigma\mathbb E[(X^\sigma)^TA^2X^\sigma]
\ge n(n-1).
```
The unresolved issue is joint sample selection and the heavy positive
branch: a single heavy positive field is not harvested by the
negative-coordinate flip.

A broader exact dependent rule reaches the same stopping point. Clip
```math
u_i=\operatorname{clip}
\left(r_i,[-K\sqrt{n-1},K\sqrt{n-1}]\right),
```
and, conditional on $`x`$, generate independent $`Y_i`$ with means
```math
\mathbb E[Y_i\mid x]
=x_i\left[(1-\gamma K)+\frac{\gamma u_i}{\sqrt{n-1}}\right],
\qquad 0\le\gamma K\le1.
```
Writing
```math
S_K=\sum_i
\min\!\left(r_i^2,K\sqrt{n-1}|r_i|\right),
```
random-sign completion of the clipped vector gives
```math
\boxed{
Q\ge
(1-\gamma K)^2q+
\frac{2\gamma(1-\gamma K)}{\sqrt{n-1}}S_K
-\gamma^2K^2Q.
}
```
This harvests unstable or negative local-field mass. If all $`r_i\ge0`$,
however, $`S_K\le Kq`$ pointwise, so the first variation cannot improve
the witness. Coordinatewise local-field feedback therefore stops
exactly at a globally cut-stable positive-heavy configuration. A
successful successor must use a genuinely multi-vertex cut inequality.

### 3.23 Abstract flip certificates have no coercivity without cut
triangle identities

The edge-flip formula in §3.21 is exact, but its aggregate
majority-with-gap condition alone cannot force a useful lower bound on
$`M`$. In an abstract antipodal state space containing only a vector
$`w`$ and its negative $`-w`$, the $`k`$-flip stability inequalities can
be satisfied even at the parity floor. Thus no argument using only
active-layer cardinality, antipodality, parity, and intersection sizes
can yield leading-order coercivity.

Any successful use of the certificate must exploit the defining
triangle relations of genuine cut vectors:
```math
v_{ij}v_{jk}v_{ki}=\text{the same global sign}
\quad\text{for every triangle }ijk.
```
That test has now been completed in §3.24.

### 3.24 Local triangle rigidity reduces exactly to covering radius

Let $`H=(V,S)`$ be a support graph, with $`c(H)`$ components and cycle
rank
```math
\beta(H)=|S|-|V|+c(H).
```
The restriction of the augmented cut code to $`S`$ has dimension
```math
\dim(\mathcal C_n|_S)=
\begin{cases}
|V|-c(H),&H\text{ bipartite},\\
|V|-c(H)+1,&H\text{ nonbipartite},
\end{cases}
```
and therefore codimension
```math
\boxed{
q(H)=
\begin{cases}
\beta(H),&H\text{ bipartite},\\
\beta(H)-1,&H\text{ nonbipartite}.
\end{cases}
}
```
Indeed, ordinary cuts are the image of the binary vertex-edge incidence
map. The global complement bit adds the all-one edge vector, which
already lies in the cut space exactly when $`H`$ is bipartite.

The dual local constraints are precisely even-cardinality Eulerian
subgraphs. Consequently every forest, matching, star, single triangle,
and connected odd unicyclic graph sees the full local cube. The first
genuine constraints occur on a $`4`$-cycle or between two independent
odd cycles. Moreover, every arbitrary local pattern is within
$`q(H)`$ edge changes of a consistent augmented-cut pattern, by
pivoting a full-rank parity-check matrix. Bounded-cycle-excess tests
therefore have only bounded coercive power.

The global reduction is also exact. For an edge cochain
$`\alpha\in\mathbb F_2^E`$, define
```math
(\delta\alpha)_{ijk}
=\alpha_{ij}+\alpha_{jk}+\alpha_{ki}.
```
Then
```math
\mathcal C_n=\delta^{-1}(\langle\mathbf1\rangle),
```
and the image of $`\delta`$ is the space of two-graphs satisfying the
tetrahedron parity identity. Thus
```math
\rho(\mathcal C_n)
=
\max_{\tau\ {\rm two\text{-}graph}\bmod\mathbf1}
\min\{|\alpha|:\delta\alpha\in\{\tau,\tau+\mathbf1\}\}.
```
The complete triangle system plus the all-flip majority certificate is
therefore exactly the original maximum cofilling/covering-radius
problem, not an extra regularity condition.

Rooting makes the obstruction transparent. After switching
$`a_{1i}=1`$, the rooted triangle signs
$`\tau_{1ij}=a_{ij}`$ are an arbitrary signing of $`K_{n-1}`$, and
with $`x_1=1`$ the energy is
```math
\sum_{i>1}x_i+
\sum_{2\le i<j}\tau_{1ij}x_ix_j.
```
This is precisely the already unclosed affine recurrence.

The exact $`n=10`$ audit reinforces the theorem: every triangle and
every $`4`$-cycle is certified by the top energy layer, whereas the
first top-layer failures occur on perfect matchings, which are forests
and carry no triangle constraint. The missing nonlocal statistic is the
conditional extension-gap profile
```math
\Gamma_H(y)=\min\{g_v:N_v|_{E(H)}=y\};
```
triangle parity determines its domain but gives no control of its
values. The full proof is in `triangle_rigidity_reduction.md`.

### 3.25 Exact insertion profiles either fail closure or contain the
whole landscape

With the half-energy normalization
```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,
```
the gauge-fixed insertion recurrence is
```math
\boxed{
M_{n+1}
=\min_A\max_x
\left(
|H_A(x)|+\left|\sum_i x_i\right|
\right).
}
```
The complete signed magnetization-extrema profile
```math
U_A(m)=\max_{\sum x_i=m}H_A(x),
\qquad
L_A(m)=\min_{\sum x_i=m}H_A(x)
```
does not close under choosing the next signed row. Two explicit
order-$`7`$ minimizers have identical $`(L_A(m),U_A(m))`$ for every
magnetization $`m`$, and indeed identical radial external-field
supports for every real field, but their best extension values are
```math
\boxed{E(A_1)=12,\qquad E(A_2)=10.}
```
The two matrices and an exhaustive $`128\times128`$ integer
reproducer are in `scale_transfer_profile_no_go.md`.

The exact closed state is the pair of external-field supports
```math
F_A^\sigma(h)=
\max_x\bigl(\sigma H_A(x)+h\cdot x\bigr),
\qquad \sigma=\pm1.
```
For
```math
C=\begin{pmatrix}A&b\\b^T&0\end{pmatrix},
```
it obeys
```math
F_C^\sigma(h,t)
=\max_{y=\pm1}
\left[
ty+F_A^\sigma(h+\sigma yb)
\right].
```
But this closure has no compression: at every cube vertex $`x`$,
```math
\boxed{
H_A(x)
=\lim_{\lambda\to\infty}
\left(F_A^+(\lambda x)-\lambda n\right),
\qquad
(F_A^+)^*(x)=-H_A(x).
}
```
Thus the smallest evident exact insertion profile contains the entire
switching-energy word. A useful convergence proof needs a genuinely
lossy asymptotic regularity theorem for high-energy level sets, not
another exact finite recursion.

### 3.26 Nonuniform random-cut moments below $`n^{1/3}`$ are
asymptotically tautological

Let $`D`$ be switched so
```math
R=\mathbf1^TD\mathbf1=Q(D),
```
choose independent cut indicators $`Z_i\sim{\rm Bernoulli}(p_i)`$, put
$`v=1-2p`$, and let $`C=C_D(Z)`$. Then
```math
\boxed{
\mathbb EC=\frac{R-v^TDv}{4}.
}
```
Hence $`0\le\mathbb EC\le R/2`$ is exactly the existing cube inequality
$`|v^TDv|\le R`$, even when each $`p_i`$ is an arbitrary function of
the row sums.

Writing $`\eta=Z-p`$, the exact Hoeffding decomposition is
```math
C-\mathbb EC
=(Dv)\cdot\eta
-2\sum_{i<j}D_{ij}\eta_i\eta_j.
```
Orthogonality of the Bernoulli monomials and $`D_{ij}^2=1`$ give
```math
\boxed{
\operatorname{Var}C
=
\sum_i p_i(1-p_i)(Dv)_i^2
+4\sum_{i<j}p_i(1-p_i)p_j(1-p_j).
}
```
Therefore
```math
\operatorname{Var}C
\le\frac n4\|D\|_{\rm op}^2+\frac{n(n-1)}8.
```
For a competing sequence $`R=c_n n^{3/2}`$, $`c_n\ge c_0>0`$, the
spectral bootstrap gives
```math
\frac{\operatorname{Var}C}{R^2}
\le
\frac14c_n^{-4/3}n^{-1/3}
+\frac18c_n^{-2}n^{-1}
=o(1).
```
A Hanson--Wright/Khintchine bound upgrades this: every normalized
moment of order $`k=o(n^{1/3})`$ collapses to the corresponding power
of the mean. Fixed-order cut moments and cycle statistics therefore
see only a point in $`[0,1/2]`$.

This loss is real. If the dense constraint is relaxed, take $`m`$
disjoint edges of weight
```math
w=\frac n{\sqrt{2m}},
\qquad m=\left\lfloor\frac{c^2n}{2}\right\rfloor.
```
This weighted matching (not a sign matrix) satisfies
```math
Q=(c+o(1))n^{3/2},\qquad
\sum r_i^2=n^2,\qquad
\|W\|_F^2=n^2,
```
all row sums are nonnegative, every cut lies in $`[0,Q/2]`$, and the
spectral bootstrap holds, for any fixed $`c>0`$. Thus row moments,
cut positivity, and all $`o(n^{1/3})`$-order cut statistics cannot use
the essential condition $`|D_{ij}|=1`$. The surviving target is a
growing-order endpoint statistic, such as the overlap-resolved number
of pairs of cuts within $`o(R)`$ of $`R/2`$.

### 3.27 Heavy positive rows can be manufactured at negligible leading
cost

Let $`D`$ have order $`m`$, switched so
```math
\mathbf1^TD\mathbf1=Q(D)=R.
```
Adjoin $`k`$ universally positive vertices:
```math
\widetilde D_k=
\begin{pmatrix}
J_k-I_k&J_{k,m}\\
J_{m,k}&D
\end{pmatrix}.
```
Then the value is exact:
```math
\boxed{
Q(\widetilde D_k)=R+2km+k(k-1).
}
```
For Boolean block sums $`s,t`$, the energy is
```math
s^2-k+2st+y^TDy,
```
whose absolute value is at most the displayed quantity, with equality
at the all-one vector.

Consequently $`k=o(\sqrt m)`$ leaves the leading normalized value
unchanged. Already $`k=1`$ creates a row of size $`m`$, and changes the
row-square mass from $`S=\sum r_i^2`$ to
```math
\boxed{
\widetilde S=m^2+S+2R+m,
}
```
at energy cost only $`2m=o(m^{3/2})`$. This is an actual dense sign
matrix, not a weighted relaxation. It proves that global cut
positivity, mean row sum, and a lower bound on the second row moment
cannot by themselves turn exceptional positive-heavy rows into a
leading-order gain.

For an exact conference switching, the full two-dimensional row-threshold
consequence can also be solved. Put
```math
U=\frac D{\sqrt{n-1}},\qquad
u=U\mathbf1,\qquad
\rho=\mathbb Eu.
```
Then $`\mathbb Eu^2=1`$, $`U\mathbf1=u`$, $`Uu=\mathbf1`$, and choosing
$`x_i=\operatorname{sign}(1-u_i)`$ in the two-dimensional spectral
projection gives
```math
\boxed{
\mathbb E|1-u|\le\sqrt{1-\rho^2}.
}
```
But among distributions with $`0\le u\le L`$,
$`\mathbb Eu=\rho`$, and $`\mathbb Eu^2=1`$,
```math
\inf\mathbb E|1-u|
=(1-\rho)\left(1+\frac2L\right),
```
attained on $`\{0,1,L\}`$. Since $`L=\sqrt{n-1}`$ is possible, the
threshold inequality is asymptotically vacuous.

The correct surviving statistic is the uniformly-integrated row-square
tail
```math
\Psi_n(K)=
\frac1n\sum_i
\left(\frac{r_i}{\sqrt n}\right)^2
\mathbf1_{\{r_i>K\sqrt n\}}.
```
A useful continuation needs a regular-versus-peeling theorem: prove a
stronger result when $`\Psi_n`$ is uniformly integrable, and show that
the exceptional vertices can otherwise be peeled or structurally
converted without losing the $`n^{3/2}`$ objective. The complete proof
is in `multicut_heavy_field_note.md`.

### 3.28 Signed-Johnson endpoint hierarchy and the flat-Sidon no-go

Switch a signing $`D`$ so
```math
R=\mathbf1^TD\mathbf1=Q(D),
\qquad
C_D(S)=\sum_{i\in S,j\notin S}d_{ij}
=\frac{R-q_D(\mathbf1_S)}4.
```
For $`k`$-subsets define the signed Johnson operator
```math
(T_k)_{S,S-\{i\}+\{j\}}=d_{ij}.
```
Restricting the Fourier convolution matrices of the nonnegative cube
polynomials $`R\pm q_D`$ to level $`k`$ gives the exact hierarchy
```math
\boxed{-\frac R2I\preceq T_k\preceq\frac R2I.}
```
It also satisfies
```math
T_k\mathbf1(S)=C_D(S),
\qquad
\frac{\langle\mathbf1,T_k\mathbf1\rangle}{\binom nk}
=\frac{k(n-k)}{n(n-1)}R.
```
For every Boolean $`x`$, with
$`\psi_x(S)=\prod_{i\in S}x_i`$,
```math
\frac{\langle\psi_x,T_k\psi_x\rangle}{\|\psi_x\|_2^2}
=\frac{k(n-k)}{n(n-1)}q_D(x).
```

Entrywise flatness fixes only
```math
\binom nk^{-1}\operatorname{tr}T_k^2=k(n-k).
```
Together with the spectral interval this yields merely
```math
R\ge2\sqrt{k(n-k)}=O(n),
```
even at $`k\asymp n`$. Thus Frobenius/trace information misses the
required $`n^{3/2}`$ scale. A successful growing-order hierarchy must
use the correlated cycle structure inherited from the same base-edge
signing.

Likewise,
```math
\mathbb EC_D(S)=\frac{k(n-k)}{n(n-1)}R
```
for a uniform $`k`$-set. Markov endpoint counts produce exponentially
many near-maximizers, but every resulting moment inequality is
scale-free because $`R`$ cancels from both sides.

Finally, $`M_n`$ is the minimum sup norm of a *flat* degree-two
tetrahedral Walsh polynomial. This is not the ordinary level-two Sidon
constant, which optimizes over arbitrary coefficient magnitudes.
Already at $`n=4`$, the best flat ratio is
```math
\binom42/M_4=3/2,
```
whereas the four-cycle polynomial with coefficients of magnitude
$`1/2`$ has coefficient $`\ell_1`$-norm $`2`$ and sup norm $`1`$.
Ordinary Sidon asymptotics therefore do not determine the flat
minimum. Direct sums introduce a leading rectangular cross block, and
ordinary tensors raise degree two to degree four, so standard Sidon
stabilization does not give scale transfer.

### 3.29 Sharp insertion is endpoint-weighted discrepancy

With the half-energy normalization, define
```math
M=M(A)=\max_x|H_A(x)|,
\qquad
E(A)=\min_b\max_x\bigl(|H_A(x)|+|b\cdot x|\bigr).
```
Then exactly
```math
E(A)-M
=
\min_b\max_x
\left[
|b\cdot x|-\bigl(M-|H_A(x)|\bigr)
\right].
```
After switching a positive absolute maximizer to $`\mathbf1`$, put
```math
C_D(S)=\sum_{i\in S,j\notin S}d_{ij}.
```
Since
```math
H_D(y_S)=M-2C_D(S),
\qquad
0\le C_D(S)\le M,
```
the recurrence becomes the exact weighted-discrepancy identity
```math
\boxed{
E(D)-M
=
\min_b\max_S
\left[
|b\cdot y_S|
-2\min\{C_D(S),M-C_D(S)\}
\right].
}
```

There is a useful mixed, but not pure, insertion theorem. If
```math
d=M-|H_A(x)|
```
and the local fields are oriented toward the sign of $`H_A(x)`$, then
```math
-d/2\le\ell_j,
\qquad
\sum_j\ell_j=2(M-d),
```
and hence
```math
\sum_j|\ell_j|\le2M+(n-2)d.
```
Choosing a uniformly random duplicate row of $`A`$, with its missing
diagonal entry filled by either sign, gives against every fixed state
```math
\boxed{
\mathbb E_j\bigl[
|H_A(x)|+|b^{(j)}\cdot x|
\bigr]
\le
M+\frac{2M}{n}+1-\frac{2d}{n}.
}
```
The quantifier order is the entire gap: insertion requires one pure
row working for every endpoint cut. This isolates a purification
problem, not a scalar entropy estimate.

The exact ground-state obstruction is
```math
E(A)-M
\ge
\min_b\max_{\{|H_A(x)|=M\}}|b\cdot x|.
```
A full orthogonal Boolean basis of extremizers forces overhead at least
$`\sqrt n`$, although a nontrivial conference matrix cannot have such a
Boolean eigenbasis by a simple integrality argument. The complete note
is `sharp_insertion_weighted_discrepancy.md`.

### 3.30 Conditional regular peeling and replenishment towers

Let $`A_t`$ be nested principal cores, $`q_t=Q(A_t)`$, and at each
level switch a $`+q_t`$ maximizer to $`\mathbf1`$. Peel a set $`H_t`$,
and define
```math
R_t=\sum_{i\in H_t}r_i,\qquad
h_t=\mathbf1_{H_t}^T(A_t)_{H_t,H_t}\mathbf1_{H_t},
```
```math
d_t=q_t-q_{t+1}.
```
If $`e_t`$ is the energy in $`A_{t+1}`$ of the old maximizer restricted
to the new core, put
```math
g_t=q_{t+1}-e_t\ge0.
```
Direct deletion gives the exact identity
```math
\boxed{2R_t=d_t+h_t+g_t.}
```

For disjoint peeled blocks,
```math
\sum_t|h_t|
\le\sum_tQ(A[H_t])
\le2q_0.
```
The last inequality follows by choosing a ground state on every block,
retaining a common energy-sign class carrying at least half the total,
and randomizing the block-global signs to cancel cross terms. Also
```math
\sum_td_t=q_0-q_L\le q_0.
```
Consequently,
```math
\boxed{
2\sum_tR_t\le3q_0+\sum_tg_t.
}
```
If every peeled coordinate has $`r_i>K\sqrt n`$, then
```math
\boxed{
\left|\bigcup_tH_t\right|
\le
\frac{3q_0+\sum_tg_t}{2K\sqrt n}.
}
```
Thus $`\sum_tg_t=O(q_0)`$ would make the heavy set $`O(n/K)`$, and
$`K\to\infty`$ would leave an $`n-o(n)`$ principal regular core at no
leading cost.

The sole exact obstruction is now a **replenishment tower**:
successive deletions reveal substantially different maximizers and
```math
\sum_tg_t\gg Q(A_0).
```
The spectral localization theorem in §1.15 has the same obstruction:
deleting one $`o(n)`$-vertex anchor may reveal a new high mode. A
convergence proof along this route must bound cumulative replenishment,
or assemble a large-gap tower into a Boolean witness exceeding
$`Q(A_0)`$. The full checkpoint is in
`regular_peeling_replenishment.md`.

### 3.31 Falsified: a uniformly regular $`n-o(n)`$ principal core

The following plausible strengthening of spectral localization is
false:

> $`Q(A_n)=O(n^{3/2})`$ implies that deleting $`o(n)`$ vertices leaves
> a principal matrix with operator norm $`O(\sqrt n)`$, with one
> uniform implied constant.

For $`n=N^2`$, partition the vertices into positive cliques. Let type
$`j`$ occupy a fraction
```math
p_j\asymp j^{-3}
```
of the vertices in blocks of size
```math
k_j=K_jN,\qquad K_j=aj,
```
and put independent signs between distinct blocks. The total internal
quadratic norm is bounded by
```math
\sum_B|B|^2
=
\left(\sum_jp_jK_j+o(1)\right)n^{3/2}
=O(n^{3/2}),
```
while a Hoeffding union bound gives a realization of the between-block
signing with quadratic norm $`O(n^{3/2})`$.

For every fixed $`C`$, choose a fixed type with $`K_j>2(C+1)`$.
That type occupies $`p_jn+o(n)`$ vertices. Any $`o(n)`$ deletion leaves
one of its cliques with more than $`(C+1)\sqrt n`$ vertices, whose
principal Rayleigh quotient exceeds $`C\sqrt n`$. Hence no
$`o(n)`$-deletion can produce a uniform $`O(\sqrt n)`$ core.

The obstruction can be implanted into any $`O(n^{3/2})`$ signing at an
arbitrarily small fixed normalized cost by overwriting only the clique
interiors. It cannot retain the obstruction at $`o(n^{3/2})`$ cost;
this indicates that the correct optimizer-specific theorem must charge
spectral-tail mass to the excess above the asymptotic optimum.

There is an optimal qualitative replacement. Grothendieck--Pietsch
factorization and
```math
\|A\|_{\infty\to1}\le2Q(A)
```
imply that for every fixed $`\varepsilon>0`$, some principal set
$`U`$ satisfies
```math
\boxed{
|U|\ge(1-\varepsilon)n,
\qquad
\|A[U]\|_{\rm op}
\le
\frac{4K_GQ(A)}{\varepsilon n}
=O_\varepsilon(\sqrt n).
}
```
Thus regularization is always possible after a fixed fractional
deletion, but its constant must diverge as
$`\varepsilon\downarrow0`$. The full construction and proof are in
`spectral_peeling_counterexample.md`.

### 3.32 Falsified: the sharp replenishment bound

The attractive conjecture
```math
\sum_tg_t\le2Q(A_0)
```
for singleton suffix peeling is false.  An explicit $`15\times15`$
signing has
```math
Q(A_0)=62,\qquad
\max_{\text{ground-state ties}}\sum_tg_t=128>124=2Q(A_0).
```
The suffix norms, from orders $`15`$ down to $`1`$, are
```math
62,58,52,44,42,30,28,28,18,14,12,8,6,2,0,
```
and one maximizing tie choice gives replenishment gaps
```math
20,16,12,16,8,16,16,0,8,8,4,4,0,0.
```
These values were independently re-enumerated over every Boolean state
of every suffix; their sum is $`128`$.

There is an exact online-learning interpretation.  Encode the chosen
suffix ground state at time $`t`$ as an augmented cut $`Z_t`$, and form
the triangular mosaic $`W`$ whose edge $`(i,j)`$, $`i<j`$, is read from
the last leader before $`i`$ is deleted.  Edgewise telescoping gives
```math
\boxed{\sum_tg_t=\langle A,W\rangle-Q(A)}
```
with the corresponding ordered-matrix normalization.  Thus a constant
replenishment theorem is precisely a constant adaptivity-gap theorem
for suffix follow-the-leader over augmented cuts.  The example proves
that any such constant must exceed $`2`$; it does not yet show that no
absolute constant exists.

### 3.33 Falsified: raw optimality excess controls regularity

Even an asymptotically minimizing sequence may have divergent
$`\|A\|_{\rm op}/\sqrt n`$ and non-uniformly-integrable row squares.
Take an exact minimizer $`B_m`$ along a liminf subsequence and adjoin
$`k=k_m\to\infty`$, $`k=o(\sqrt m)`$, universally positive vertices
after switching a positive maximizer of $`B_m`$ to $`\mathbf1`$.  Then
```math
\widetilde B_m=
\begin{pmatrix}
J_k-I_k&J_{k,m}\\
J_{m,k}&B_m
\end{pmatrix}
```
satisfies the exact identity
```math
\boxed{
Q(\widetilde B_m)=Q(B_m)+2km+k(k-1).
}
```
Hence its normalized value has the same liminf, while
```math
\frac{\|\widetilde B_m\|_{\rm op}}{\sqrt{m+k}}
\ge(1-o(1))\sqrt k\longrightarrow\infty.
```
The exceptional set has $`k=o(n)`$ vertices, so this kills only an
unpeeled excess-to-regularity theorem.  A post-peeling or
spectral-tail-charge theorem remains possible.

### 3.34 Proved: cumulative ground-state-frame visibility

For a block peeling step
```math
A_t=\begin{pmatrix}D_t&B_t\\B_t^\top&A_{t+1}\end{pmatrix},
\qquad
V_t=\max_{y\in\operatorname{GS}(A_{t+1})}\|B_ty\|_1,
```
comparison with a ground state of the new core gives
```math
2V_t\le
Q(A_t)-Q(A_{t+1})+Q(D_t).
```
For disjoint peeled blocks the decrements telescope, while the
same-sign-class/random-block-sign argument gives
```math
\sum_tQ(D_t)\le2Q(A_0).
```
Therefore
```math
\boxed{\sum_tV_t\le\frac32Q(A_0).}
```
Thus a long replenishment tower can exist only by repeatedly hiding its
cross modes from the exact ground-state frames of the successive cores.

### 3.35 Falsified: regular pointwise visibility inverse

Spectral regularity does not imply $`g_t\le C(K)V_t`$ at a single step.
The order-nine suffix of the order-15 replenishment counterexample has
```math
Q(T)=Q(E)=28,\qquad g=16,
```
but its order-eight core has exactly one absolute ground-state pair and
the deleted cross row is orthogonal to both, so $`V=0`$.  Moreover the
integer matrix $`81I-4T^2`$ is positive definite, hence
```math
\|T\|_{\rm op}<\frac92=\frac32\sqrt9.
```
The exact matrix, Boolean witnesses, and Sylvester certificate are in
`optimality_excess_and_regularization.md`.  Any viable inverse must
group scales, use a near-ground-state layer, or charge a statistic other
than exact ground-state visibility.

### 3.36 Proved: existential regularized adaptivity-gap bound

If
```math
\|A\|_{\rm op}\le K\sqrt n,
```
there is a singleton principal-deletion order for which
```math
\boxed{
\sum_tg_t\le2Kn^{3/2}-Q(A).
}
```
Indeed, on a current core $`S`$, choose a positive ground state $`x`$
and write its oriented local fields as
$`\ell_i=x_i(A[S]x)_i\ge0`$.  Uniformly deleting $`i\in S`$ gives
```math
\mathbb E_i\ell_i^2
=|S|^{-1}\|A[S]x\|_2^2
\le K^2n.
```
The exact singleton identity is $`2\ell_i=d_i+g_i`$; summing conditional
expectations and using $`\sum_td_t=Q(A)`$ proves the claim.  A
deterministic order follows by taking a coordinate no larger than the
root-mean-square field at every step.  With any universal
$`Q(A)\ge c_0n^{3/2}`$, this is a $`C(K)Q(A)`$ bound.

This settles the existential regularized suffix-FTL question, but not
the heavy-core version: the low-field order is not constrained to
delete coordinates above the peeling threshold.

### 3.37 Stopped upgrades: capped bias and average deletion

For capped field-biased deletion
```math
p_i\propto\min\left(1,\frac{\ell_i}{H\sqrt m}\right),
```
scale-local regularity gives only
```math
\sum_iw_i\le Km/H,\qquad
\mathbb E_p\ell_i\le K^2m^{3/2}/H
```
while a heavy coordinate exists.  Across a linear number of deletions
this costs $`O(K^2n^{5/2}/H)`$, a polynomial factor above the target,
and the hazard for a lone heavy coordinate is only $`H/(Km)`$.
Moving or late-created heavy fields remain uncontrolled.

For uniform one-vertex deletion, the exact first-step identity is
```math
\mathbb E_i d_i+\mathbb E_i g_i=\frac{2Q(A)}n.
```
Thus the desired average contraction coefficient $`3/2`$ is equivalent
to the new estimate $`\mathbb E_i g_i\le Q(A)/(2n)`$.  The regularized
cumulative theorem does not imply it.  Regularity alone is insufficient:
the order-nine matrix in §3.35 has
```math
Q(T_{-i})=Q(T)=28
\quad\text{for all }i,
```
despite $`\|T\|_{\rm op}<1.5\sqrt9`$.  Even an order-seven exact
minimizer has this all-deletions-flat property.  Any asymptotic theorem
must exploit large-order optimality beyond spectral regularity.

### 3.38 Max-plus insertion towers require new asymptotic types

For a core $`B`$ of order $`m`$, adjoining a row $`x`$ makes
$`(1,x)`$ a parent ground state exactly when
```math
E_B(x)
=
\max_z\left(|E_B(z)|-4\delta_H(x,z)\right).
```
After switching $`x`$ to $`\mathbf1`$, this is equivalent to the
all-cut condition
```math
\boxed{
-|S|\le c_B(S,S^c)\le |S|+\frac{E_B(x)}2
\quad (|S|\le m/2).
}
```
If $`E_B(x)=Q(B)-4r`$, an exact layer form is
```math
\boxed{
\delta_H(x,Z_s)\ge r-s
\quad\text{for every higher absolute layer}\quad
Z_s=\{z:|E_B(z)|=Q(B)-4s\}.
}
```

A certified recursive tower satisfying this max-plus condition exists
through order $`25`$.  It has
```math
Q(A_{25})=228,\qquad \sum g_t=332,
```
so the finite data do not indicate a diverging adaptivity ratio.

There is a rigorous obstruction to periodic extrapolation.  If $`u_m`$
is a preceding ground state and $`d_m=Q(A_m)-Q(A_{m-1})`$, then
```math
d_m\ge2|\langle x_m,u_m\rangle|.
```
Consequently
```math
\sum_m|\langle x_m,u_m\rangle|
\le\frac{Q(A_N)-Q(A_{\rm start})}{2}.
```
Any positive-density rule with linear overlap therefore forces
$`Q(A_N)=\Omega(N^2)`$.  More generally, every fixed finite
vertex-type kernel $`a_{ij}=K(c_i,c_j)`$ has a nonzero Boolean Fourier
mode and hence quadratic norm.  A subquadratic max-plus tower, if one
exists, must continually generate genuinely new near-orthogonal
types; bounded-period and bounded-type constructions are stopped.
The exact order-25 certificate is in `regular_peeling_tower25.md`.

---

## 4. Pending audit — do not cite as proved yet

### 4.1 Resolved: finite arcsine stability refinement

The earlier claimed quantitative remainder has now been audited and
proved with the stronger denominator
```math
8\pi(n-1)^{5/2}
\left(1-\frac1{n-1}\right)^{3/2}.
```
Its proof and its factor-$`n`$ scale limitation are recorded in
`orientation_even_stability_audit.md`; the stronger leading-scale
$`A^2`$-energy theorem is in §1.16.

### 4.2 External Claude campaign: useful leads, proof files unavailable

A separate campaign supplied a detailed summary, but its named technical
note and certificate bundles were not available in the shared file store.
The following are therefore research leads rather than accepted inputs:

1. **Haar half-subspace benchmark.** For a Haar-random rank-$`p/2`$
   projection $`P`$, the claimed limit is
   $`\displaystyle \max_{x\in\{\pm1\}^p}\frac{\|Px\|^2}{p} \longrightarrow \beta_*=\frac12\left(1+\frac{\sqrt{15}}4\right).`$
   The union-bound side is immediate from the beta tail:
   $`\displaystyle \beta_*(1-\beta_*)=\frac1{64}.`$
   The asserted matching second-moment lower bound, especially the
   microscopic-overlap truncation, still requires reconstruction.
   In the quadratic-energy normalization this suggests
   $`\displaystyle c_{\rm Haar}=\frac{\sqrt{15}}8 =0.4841229182\ldots .`$
2. **Paley resonance.** This claim has now been reconstructed and moved
   to the verified results in §1.12.
3. **Important logical caveat.** That limsup theorem alone does not prove
   that the Paley values fail to converge. The same summary labels
   convergence to $`\sqrt{15}/8`$ on density-one nonresonant primes as a
   conjecture, and §1.12 now disproves the density-one part. A proved
   low-valued second subsequence is still missing.
4. **Entropic alternative.** The summary proposes that rare random
   signings might undercut $`c_{\rm Haar}`$, governed by a lower-tail or
   Franz--Parisi rate. A speed-$`n^2`$ rate can locate the onset of
   exponentially many good signings, but by itself cannot exclude a
   single algebraic switching orbit, whose probability already has the
   maximal cost $`(\log2)n^2/2+o(n^2)`$.

---

## 5. Current ten-route research cycle

The workspace permits three subagents concurrently, so routes run in waves.
A route survives only if it produces a scale-transfer inequality, a
proof-grade structural lemma, or a decisive obstruction.

| # | Route | Concrete target | Status |
|---:|---|---|---|
| 1 | Explicit conference AMP | Valid Onsager recursion reaching the spectral ceiling | Old claim falsified; paired recursion capped near $`0.705`$ |
| 2 | Conference cut positivity | From all-cut positivity and $`C^2=(n-1)I`$, prove $`R\ge(1-o(1))n^{3/2}`$ | Heavy-row dichotomy proved; near-$`1`$ conclusion still open |
| 3 | Spectral-moment dichotomy | Either direct Boolean witness $`\ge n^{3/2}`$, or strong conference-like structure | Produced verified $`0.336493`$ theorem; sharp dichotomy still open |
| 4 | Higher-order stability | Upgrade near-optimality to a two-eigenvalue/ETF approximation strong enough for rounding | Orientation-odd hierarchy obstructed by self-complementary signings |
| 5 | Near-ground-state entropy | Prove an $`O(\sqrt n)`$ vertex-insertion inequality | Uniform theorem falsified at $`n=6`$; asymptotic entropy lemma unknown |
| 6 | Correlated nonlocal lifts | Preserve scalar geometry under amplification without vector-relaxation loss | Flat-Fourier/orbit mechanism falsified |
| 7 | Scale-preserving free energy | Prove convergence on $`t=\beta/\sqrt n`$, not merely fixed $`t`$ | Fixed diagrams universal, but nonperturbative resonance kills all-$`\beta`$ transfer |
| 8 | Signed coding dual | Build a signed high-degree certificate for the augmented cut-code covering radius | Exact dual is the original histogram; stopped |
| 9 | Optimizer mining | Discover or falsify scalable algebraic structure from exact/heuristic optima | Two-copy recurrence and Paley-minor hypothesis falsified |
| 10 | Energy–entropy compactness | Define a full profile closed under gluing and prove a unique asymptotic value | Fixed profiles/LDP shortcut obstructed by sparse planting |

Current successor routes:

- regular-versus-peeling control of the positive-heavy row-square tail;
- joint selection of high energy and high $`A^2`$-energy in the
  opposite-orientation rounding;
- growing-order endpoint/overlap statistics for the dense cut
  polynomial.

Recently completed successor routes:

- **Paley nonresonance:** stopped with a stronger negative theorem:
  every fixed admissible progression has Paley limsup $`1/2`$.
- **Dependent local-field rounding:** exact $`A^2`$-energy and
  instability corrections proved; all coordinatewise variants stop at
  the cut-stable positive-heavy branch.
- **Triangle rigidity:** exact local-rank and two-graph reductions prove
  that bounded local triangle tests add no coercion beyond the original
  covering-radius problem.
- **Insertion profiles:** magnetization profiles fail closure, while
  the exact external-field closure is injective and contains the full
  energy word.
- **Nonuniform random cuts:** all $`o(n^{1/3})`$-order moments collapse
  to the cube-norm first moment and cannot see dense sign flatness.
- **Positive heavy fields:** universal-positive-vertex extensions show
  that exceptional heavy rows and an $`n^2`$ second row moment can be
  created at only $`O(n)`$ energy cost; first/second moments are
  insufficient without a peeling theorem.
- **Raw optimality-excess regularity:** falsified even for asymptotically
  minimizing sequences; the anomaly is always localized in the explicit
  construction.
- **Sharp replenishment coefficient $`2`$:** falsified at order $`15`$;
  the surviving question is an absolute or regularized adaptivity-gap
  bound.

---

## 6. Most useful sufficient lemmas

Any one of the following would settle convergence.

### 6.1 Uniform amplification

```math
\frac{M_{kn}}{(kn)^{3/2}}
\le
\frac{M_n}{n^{3/2}}+\varepsilon_n,
\qquad
\varepsilon_n\to0,
```
uniformly in $`k`$.

### 6.2 Sharp one-vertex insertion

```math
M_{n+1}
\le
M_n+\left(\frac32+o(1)\right)\frac{M_n}{n}.
```

This would make $`M_n/n^{3/2}`$ asymptotically nonincreasing.

### 6.3 Universal sharp lower bound

```math
Q(A)\ge(1-o(1))n^{3/2}
\quad\text{for every signing }A.
```

Together with Paley conference upper bounds, this would prove
```math
\frac{M_n}{n^{3/2}}\to\frac12.
```

### 6.4 Convergent finite-temperature minima

Define
```math
F_n(\beta)=
\frac1n\min_A
\log\mathbb E_x
\exp\!\left(\frac{\beta H_A(x)}{\sqrt n}\right).
```

If $`F_n(\beta)`$ converges for every fixed $`\beta>0`$, uniformly enough to
send $`\beta\to\infty`$, then
```math
\frac{M_n}{n^{3/2}}
```
converges, because the zero-temperature approximation error is at most
$`\log2/\beta`$.

The unresolved issue is convergence along the changing raw temperature
$`t=\beta/\sqrt n`$.

---

## 7. Immediate next actions

1. Prove the positive-heavy-field tail or peeling dichotomy needed to
   turn the verified relative-invariance $`A^2`$-gain and capped-field
   conversion into a strictly stronger universal lower bound.
2. Use large-order global optimality—not spectral regularity alone—to
   prove an averaged principal-deletion contraction at coefficient
   $`3/2`$, or find its asymptotic obstruction.
3. Test degree-$`\Theta(n)`$ moment/Laplace methods for the flat
   Rademacher chaos; fixed-degree diagrams are already known to be
   insufficient.
4. Seek a grouped or near-ground-layer inverse for replenishment.
   Hard max-plus towers of bounded period/type are now ruled out.
5. Build a closed scale-transfer profile, or prove that a proposed
   energy--magnetization profile is still insufficient.
6. Revisit Paley/nonabelian/random correlated upper constructions only
   if they include a moving $`p`$-dependent condition and a mechanism
   that controls all $`2^p`$ Boolean vectors.
7. Revisit other correlated lifts only if they include a
   mechanism that controls zero-entropy resonant cube points.
8. Save a new checkpoint after every agent wave or any material
   proof/counterexample.

## 8. Proof acceptance policy

A result is moved to **Verified** only when:

1. all definitions and factors of $`2`$ are explicit;
2. asymptotic error terms are uniform in the required matrix class;
3. any Gaussian/cavity step states the dependency structure and Onsager
   terms;
4. a second derivation, independent audit, or falsifying computation has
   been attempted;
5. the result advances the original convergence question, or is clearly
   labeled as a bound/obstruction only.

---

## 9. Checkpoint: sparse repair and linear-order moments

### 9.1 Exact Bernoulli repair criterion (provisional, being audited)

Let $`D`$ be an $`m`$-vertex signing and put
```math
G_x=Q(D)-|x^\top D x|.
```
Independently flip each surviving edge with probability $`p`$, and set
$`\mu=1-2p`$.  A Bernstein bound followed by the union bound gives a
concrete sufficient condition for a realization $`B`$ satisfying
$`Q(B)\le \mu Q(D)+s`$:
```math
2\sum_{x/\{\pm1\}}
\exp\!\left(
-\frac{(s+\mu G_x)^2}
{16pm^2+\frac83(s+\mu G_x)}
\right)<1. \tag{9.1}
```
For a principal core of order $`m=N-h`$, taking $`p=\lambda h/N`$ and
```math
s<\left(2\lambda-\frac32\right)\frac hN\,q_N
```
would give the desired $`3/2`$-coefficient mesoscopic contraction.
Thus this route has been reduced to a weighted near-ground-layer entropy
bound.  A frozen layer with logarithmic size $`o(pN)`$ is sufficient;
an exponentially rich near-ground layer is the only remaining
obstruction.  Operator-norm regularity and a black-box Hanson--Wright
estimate do not by themselves reduce that entropy enough.

This criterion is not yet promoted to the verified theorem list: its
normalizations, the absolute-value union, and the full $`h,p,s`$ regime
are undergoing a second audit.

### 9.2 Exact convolution representation of high moments

Let $`f_A`$ be the signed indicator of the weight-two sphere in
$`\mathbb F_2^n`$.  For
```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j
```
and every $`k\ge1`$,
```math
\mathbb E H_A^{2k}
=\sum_{T\subset[n]}\bigl(f_A^{*k}(T)\bigr)^2. \tag{9.2}
```
The top boundary $`|T|=2k`$ equals $`k!\operatorname{haf}(A_T)`$.
At fourth order this gives
```math
\mathbb E H_A^4
=N^2+4\sum_{i<j}(A^2_{ij})^2
+4\sum_{|T|=4}\operatorname{haf}(A_T)^2,
\qquad N=\binom n2. \tag{9.3}
```

### 9.3 Positive-diagram strategy fails at $`k=\Theta(n)`$

The positive diagrams made of $`k`$ distinct edges, each repeated twice,
already contribute
```math
D_{N,k}=\binom Nk\frac{(2k)!}{2^k}.
```
For $`k=\alpha n`$,
```math
D_{N,k}^{1/(2k)}
\sim\sqrt{\frac{\alpha}{e}}\,n^{3/2}.
```
On conference signings the entire $`2k`$-th moment is at most
$`(\tfrac12+o(1))^{2k}n^{3k}`$.  Hence for
$`\alpha>e/4`$ this positive subfamily exceeds the complete moment:
signed Eulerian diagrams cancel it at exponential scale.  Consequently,
no proof at linear moment order may keep edge-even/pairing diagrams and
discard the rest.  Positivity survives only after the exact boundary
grouping in (9.2).

Superlinear moments are asymptotically equivalent to the original
$`L_\infty`$ problem, while linear moments retain a fixed entropy factor
and the same cancellation obstruction.  The top hafnian layer alone
also gives a constant far below the current $`0.336493`$ bound.

### 9.4 Sidon and Bohnenblust--Hille audit

The flat Sidon quantity is exactly
```math
\operatorname{Sid}_{\rm flat}(\mathcal B^n_{=2})
=\frac{\binom n2}{M_n}.
```
Existing unrestricted Sidon estimates for the degree-two Boolean
Walsh space are only sharp in order, not in the leading flat constant,
and do not prove that a Sidon extremizer has equal coefficient
magnitudes.

A tetrahedral degree-two Boolean Bohnenblust--Hille constant $`B`$ would
give
```math
\frac{M_n}{n^{3/2}}\ge \frac{2^{-3/4}}{B}+o(1).
```
It improves the current bound only if $`B<1.767059\ldots`$.
The often quoted $`1.83737\ldots`$ is a two-variable real-polynomial
constant with square terms, not a dimension-free constant for the
Boolean tetrahedral subclass; even if inserted formally it would give
only $`0.32362\ldots`$.

### 9.5 Current rotation

The three active proof routes are now:

1. turn failure of (9.1) into an entropic two-replica witness, thereby
   closing a repair-or-energy dichotomy;
2. derive a multiscale, tail-free conversion from the full local-field
   distribution, replacing the capped positive-heavy-field estimate;
3. derive a coefficient-level Eulerian/MacWilliams recursion whose
   algebra survives the scale change that defeats scalar free energy.

---

## 10. Checkpoint: the common replenishment/entropy frontier

### 10.1 Repair theorem promoted to verified

The fixed-cardinality version of the sparse-repair argument has now
been audited in the one-copy normalization.  Let $`D`$ have $`L`$ edges,
```math
M=M(D),\qquad h_v=a\cdot v,\qquad
g_v=\frac{M-h_v}{2},
```
where $`v`$ ranges over augmented cuts.  Flip a uniformly random
$`k`$-edge set, put $`p=k/L`$, and target a drop $`0<R<2pM`$.  Define
```math
\alpha=\frac{M}{2L}-\frac{R}{4k},
\qquad
\beta=\frac{1-2p}{2k}.
```
If
```math
\boxed{
\sum_v\exp\!\left[-2k(\alpha+\beta g_v)^2\right]<1,
}
\tag{10.1}
```
then some $`k`$-edge batch produces $`M(D^S)\le M-R`$.  Replacing the
Pinsker exponent by its exact hypergeometric binary-relative-entropy
form gives a sharper valid criterion.

For a core of order $`N-h`$, $`p=\lambda h/N`$, and the desired
$`3/2`$-scale drop, failure of (10.1) forces an exact equal-energy
layer with logarithmic multiplicity
```math
\log|\mathcal L_g|
\ge
\frac{\left[
g/(h\sqrt N)+c_N(\lambda-\tfrac34)+\rho/2
\right]^2}{\lambda}\,h-o(h), \tag{10.2}
```
where $`c_N=M_N/N^{3/2}`$ and
$`M_N-M(D)=\rho h\sqrt N+o(h\sqrt N)`$.
Integrality permits unit-width bands when $`h\gg\log N`$.

There is also an exact converse certificate for a global minimizer:
for every $`k<L/2`$,
```math
\boxed{
\sum_v
\exp\!\left[
-2k\left(
\frac{M}{2L}+\frac{(1-2k/L)g_v}{2k}
\right)^2
\right]\ge1.
}
\tag{10.3}
```
Thus coefficient optimality itself forces weighted entropic rescue at
every mesoscopic flip scale.  A frozen-layer premise is not available
without an additional structural theorem.

### 10.2 Exact two-replica geometry

Two same-orientation states in one exact energy layer become, after
switching by one of them, two exact ground states whose difference set
is an exact zero cut.  For a difference partition
```math
A^x=\begin{pmatrix}D&B\\B^\top&E\end{pmatrix}
```
and gaps $`g_x,g_y`$,
```math
I_D+I_E=M-(g_x+g_y),\qquad C=g_y-g_x. \tag{10.4}
```
If $`P_D,P_E`$ are the positive block maxima, then
```math
M-(g_x+g_y)\le P_D+P_E\le M.
```
For two exact parent ground states,
```math
P_D+P_E=M
```
and $`B`$ annihilates the Cartesian product of the two positive block
ground layers.

Failure of scale transfer therefore yields exponentially many exact
zero cuts in a single ground switching.  Their
$`\mathbb F_2`$-span has dimension $`\Omega(h)`$, but the span need not
itself consist of zero cuts.  Conditional on that missing additive
closure, an exact Fourier/type calculation gives
```math
\sum_{\phi_i+\phi_j=\psi}a_{ij}=0
\quad(\psi\ne0)
```
for vertex evaluation types $`\phi_i\in W^*`$, and forces a frozen type
of size at least $`1+2M/n=\Omega(\sqrt n)`$.  The live question is now
an additive-structure theorem for the zero set of a globally
nonnegative complete-sign cut polynomial; a Sidon-like zero-cut layer
is the precise obstruction.

### 10.3 Tail-free block inequalities, and why scalar recursion stops

For
```math
A=\begin{pmatrix}D&B\\B^\top&E\end{pmatrix},
\quad x=(x_S,x_T),
\quad e=x_T^\top Ex_T,
\quad C=\|Bx_T\|_1,\quad L=\|Bx_T\|_2,
```
one has the exact bounds
```math
\boxed{Q(A)\ge |e|+\sqrt2\,L}
\tag{10.5}
```
and
```math
\boxed{
Q(A)\ge |e|+\Phi(C,Q(D)),
\quad
\Phi(C,d)=
\begin{cases}
C^2/d,&C\le d,\\
2C-d,&C\ge d.
\end{cases}}
\tag{10.6}
```
The second inequality follows by biased independent rounding on the
peeled block.  They harvest sparse heavy local-field levels without a
cap and force dense heavy levels to recurse into their induced block.

The obstruction is exact: (10.5)--(10.6) retain the old restricted
energy $`x_T^\top Ex_T`$.  Replacing it by $`Q(E)`$ introduces precisely
the successor replenishment gap.  Therefore the multiscale
heavy-field route and the deletion route are now the same problem at
different resolutions.  Universal-positive-vertex extensions also
show that no scalar one-witness $`L_2`$-tail functional can close this
gap.

### 10.4 Eulerian free energy: exact channel form and finite-hierarchy
no-go

At $`\rho=\tanh t`$, the high-temperature factor
```math
W_A(\rho)
=\mathbb E_{\sigma,x}
\prod_{i<j}(1+\rho\sigma a_{ij}x_ix_j)
```
is exactly $`2^L`$ times the probability of output $`A`$ when a uniform
augmented cutword passes through independent binary symmetric noise
of mean $`\rho`$.  Thus the minimum centered free energy is a
$`D_\infty`$ (least-output-likelihood) problem for the noisy cut code.

The edge recursion is exact, as is the vertex boundary-sector
recursion, but deleting one vertex opens all $`2^{n-2}`$ even boundary
sectors.  Its $`k`$-th row moment requires the full $`k`$-replica
overlap array.  The first nontrivial coefficient is
```math
T_4(A)=
\frac{\|A^2\|_F^2-n(n-1)(2n-3)}8
\ge-\frac{n(n-1)(n-2)}8,
\tag{10.7}
```
with equality exactly at conference matrices.

This finite hierarchy does not order the relevant partition
functions.  Two explicit order-six polynomials cross at positive
$`\rho`$, and an exact local trap
```math
W(\rho)=(1-\rho^4)^3
```
has twelve of fifteen flat edge flips and only three active
correlations, while the conference polynomial is nevertheless lower.
Hence local anti-alignment, susceptibility, $`T_4`$, and every fixed
replica truncation are stopped.  The only surviving free-energy route
is a genuinely growing-replica $`D_\infty`$ large-deviation theorem.

### 10.5 Active routes after convergence of the obstructions

1. additive-energy/zero-cut rigidity for the exact layer forced by
   (10.2);
2. grouped-scale replenishment using cumulative ground-layer
   visibility, rather than a false pointwise inverse;
3. block purification of the exact mixed insertion theorem, with
   arbitrary new-block spins explicitly controlled.

### 10.6 Block purification audit

For
```math
G=\begin{pmatrix}A&B\\B^\top&D\end{pmatrix}
```
the exact half-energy identity is
```math
M(G)=\max_{x,y}
\left(
|H_A(x)+H_D(y)|+|x^\top By|
\right). \tag{10.8}
```
It yields a sufficient weighted cross-block criterion in terms of the
joint energy-gap/local-field exponential profile.  However, the
natural duplicate-row purification does not preserve the
$`3/2`$-scale:

* the mixed-insertion $`L^1`$ slack reverses sign once two or more
  duplicate rows are accumulated;
* taking all duplicate rows with the natural internal copy gives
  $`4M(A)\pm n`$, a normalized loss of $`\sqrt2`$;
* exhaustive optimization over all internal signings gives best
  order-doubling factors $`2.50,3.25,3.60`$ for optimal cores of
  orders $`4,5,6`$, respectively.

Thus simple block purification is stopped.  Its exact surviving
condition is a joint energy--local-field large-deviation estimate of
the form
```math
\sum_x e^{-a d_A(x)/\sqrt n}
\left[
\frac{\cosh(a/\sqrt n)}n
\sum_j\cosh\!\left(\frac{a(Ax)_j}{\sqrt n}\right)
\right]^{\theta n},
\tag{10.9}
```
coupled to the corresponding deficit partition function of the
second block.  Ground-state entropy or spectral regularity alone does
not control (10.9).

### 10.7 Exact zero-cut obstructions at the smallest orders

Enumeration through order seven shows that positive ground-state
families can already have minimum possible additive energy
```math
E(Z)=3|Z|^2-2|Z|,
```
so they can be maximally Sidon.  The order-five and order-six
families also fail delta-matroid symmetric exchange.  Therefore
cardinality alone supplies no Balog--Szemerédi--Gowers leverage.

Even the full absolute-ground condition
```math
0\le C(S)\le M\qquad\text{for every cut }S
\tag{10.10}
```
does not imply positive semidefiniteness of the signed Laplacian or
closure of its zero cuts.  The smallest counterexample is the
order-five signing
```math
A=\begin{pmatrix}
0&1&1&1&1\\
1&0&1&-1&1\\
1&1&0&1&-1\\
1&-1&1&0&-1\\
1&1&-1&-1&0
\end{pmatrix},
\qquad M=4.
```
Its signed Laplacian has eigenvalues
```math
-1.828427\ldots,\ 0,\ 1,\ 3.828427\ldots,\ 5,
```
and, in the gauge excluding vertex $`0`$, its zero-cut masks are
$`\{0,4,6,8,9\}`$, while $`4\mathbin\triangle6=2`$ is not a zero cut.

The surviving exact datum is pair factorization: every realized
ground-state difference factors into the two positive ground
degeneracies of its induced principal blocks.  A Sidon layer must
therefore encode exponentially many distinct tight principal
decompositions.  Whether that can coexist with $`M=O(n^{3/2})`$ is
the next structural question.

### 10.8 Exact factorization and entropy-product no-go

If $`Z`$ is a projective positive ground family and $`d\ne0`$ is a
realized difference, then
```math
\boxed{r_Z(d)=2\gamma_P(d)\gamma_Q(d),} \tag{10.11}
```
where $`\gamma_P,\gamma_Q`$ are the projective positive-ground
degeneracies of the two principal blocks in the tight decomposition
defined by $`d`$.  Consequently
```math
E(Z)\ge3|Z|^2-2|Z|,
```
with equality exactly when $`Z`$ is Sidon, equivalently when every
realized tight decomposition has unique projective block grounds.
Every oriented exact minimizer ground family enumerated at orders
five, six, and seven is Sidon.

Large finite ground-state entropy also fails to bootstrap through the
natural saturated block products.  For any sign block $`R`$, if a
projective set $`U`$ satisfies
```math
u^\top Rv=\|R\|_{\infty\to1}
\qquad\text{for every }u,v\in U,
```
then $`|U|=1`$.  More generally, a complete sign quadratic on
$`[-1,1]^r`$ has a maximizer and a minimizer with at most one
fractional coordinate; if its value is within $`\varepsilon`$ of an
extreme, all but at most one coordinate have minority probability at
most $`\sqrt{\varepsilon/2}`$.  Thus fixed rank-one or
edgewise-saturated recursive products cannot independently tensor the
large ground layers seen at orders five and six.

### 10.9 Adaptive ground-state closure

Fix one orientation and let
```math
P(C)=\max_xx^\top Cx.
```
Starting from a state $`y`$ on a current core $`C_j`$, choose a
positive ground state $`z_j`$ and delete the smaller disagreement set
$`D_j`$, so that $`y=z_j`$ on the surviving core.  If
```math
g_j=P(C_j)-y^\top C_jy,
```
then the exact decomposition gives
```math
g_{j+1}
=\frac12g_j+a_j-\bigl(P(C_j)-P(C_{j+1})\bigr),
\qquad
a_j\le P(D_j). \tag{10.12}
```
The blocks $`D_j`$ are disjoint and independent block flips give
```math
\sum_jP(D_j)\le P(A),\qquad
\sum_jN(D_j)\le N(A).
```
More strongly, positive maxima are block-superadditive:
```math
P\!\begin{pmatrix}D&B\\B^\top&C\end{pmatrix}
\ge P(D)+P(C). \tag{10.13}
```
Indeed, use positive ground states of the two diagonal blocks and
choose their relative global sign so that the cross term is
nonnegative.  Therefore
```math
P(C_j)-P(C_{j+1})\ge P(D_j)\ge a_j,
```
and (10.12) collapses to the error-free contraction
```math
\boxed{g_{j+1}\le\tfrac12g_j.} \tag{10.14}
```
In particular,
```math
\boxed{
\sum_jg_j\le2g_1\le4Q(A),
} \tag{10.15}
```
and the last bound improves to $`2Q(A)`$ when the inherited energy is
nonnegative.

There are compatible one-sided block inequalities
```math
P(A)\ge e+\sqrt2\|By\|_2,\qquad
P(A)\ge e+\Phi(\|By\|_1,N(D)).
```
They improve cumulative successor-ground visibility to
```math
\sum_jV_j^+\le\frac{P(A)+N(A)}2\le Q(A).
```
Moreover, the closure gap supplies
$`\|B_jy_{j+1}\|_1\ge g_j/4`$, and
```math
\Phi(b/4,d)+d\ge b/2.
```
Thus the former cumulative replenishment/adaptivity gap is solved
without an induced-block range charge.  The remaining issue for scale
transfer is geometric rather than energetic: the closure chain may
delete a macroscopic disagreement block.  In the certified
order-fifteen example, a singleton deletion with gap $`20`$ has
closest-ground disagreement size $`5`$ in a core of order $`14`$.
The next target is to exploit block superadditivity and
$`3/2`$-homogeneity when such a macroscopic split occurs.

The half-contraction was stress-tested in exact integer arithmetic on
907 adaptive closure steps from 800 random signings of orders
three through ten.  Every step satisfied block superadditivity, the
identity (10.12), and $`2g_{j+1}\le g_j`$.

### 10.10 Sharp algebraic floor for ground-family counting

Let $`N=2^{2m}`$, index rows and columns by
$`(u,v)\in\mathbb F_2^m\times\mathbb F_2^m`$, and define the symmetric
Hadamard matrix
```math
K_{(u,v),(x,y)}=(-1)^{v\cdot x+u\cdot y}.
```
Its diagonal is $`1`$, so $`A=K-I`$ is a complete sign quadratic.
For every Boolean function $`g:\mathbb F_2^m\to\mathbb F_2`$,
```math
X_g(x,y)=(-1)^{x\cdot y+g(x)}
```
is a Boolean $`+\sqrt N`$-eigenvector of $`K`$.  Hence $`A`$ has at
least
```math
2^{\sqrt N-1}
```
projective positive grounds and at least
$`2^{\sqrt N-1}-1`$ tight principal decompositions, while its absolute
norm remains at the spectral $`N^{3/2}+O(N)`$ scale.  Therefore no
polynomial or $`\exp(o(\sqrt N))`$ tight-decomposition bound is
possible.

There is a matching rigidity phenomenon for a linear-sized product
cube.  If $`2d`$ vertices are paired and all $`2^d`$ pair-constant
configurations are positive grounds, each cross $`2\times2`$ block is
$`c_{ij}vv^\top`$.  Three-pair tests force
```math
c_{ij}c_{ik}c_{jk}=-1,
```
so $`c_{ij}=-\sigma_i\sigma_j`$, and an antiuniform witness has
```math
|H|=2d^2-d=\frac{n^2-n}{2}.
```
Thus a bounded-block ground cube of dimension $`\Theta(n)`$ forces
quadratic norm, whereas the symplectic construction realizes
dimension $`\Theta(\sqrt n)`$ at the desired $`n^{3/2}`$ scale.  The
live generalization is an arbitrary-affine-subspace inequality of the
form $`Q(A)\gtrsim n\dim W`$.

### 10.11 Universal no-go for duplicate-row exponential purification

Use the half-energy normalization
```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_x|H_A(x)|,
```
and for an exact minimizer define its absolute gap partition function
```math
Z_A(\lambda)=
\sum_x\exp\!\left[-\lambda\bigl(M(A)-|H_A(x)|\bigr)\right].
```
If $`M(A)=c\,n^{3/2}+o(n^{3/2})`$, the exact sparse-flip converse
implies, for every fixed $`b\ge0`$,
```math
\boxed{
\liminf_{n\to\infty}\frac1n
\log Z_A(b/\sqrt n)
\ge
\psi_c(b):=\frac{(2c-b)_+^2}{8}.
} \tag{10.16}
```
There is also a simpler floor which was initially omitted.  Since
$`0\le M(A)-|H_A(x)|\le M(A)`$ for all $`2^n`$ states,
```math
\boxed{
\liminf_{n\to\infty}\frac1n
\log Z_A(b/\sqrt n)
\ge
\chi_c(b):=(\log2-bc)_+.
} \tag{10.17}
```
This already contradicts the formerly proposed ceiling
$`\log Z_A(c/\sqrt n)\le(c^2/8+o(1))n`$, since
$`\log Z_A(c/\sqrt n)/n\ge\log2-c^2-o(1)`$.

There is a stronger signing-independent local-field floor.  Put
```math
\ell(b)=\mathbb E_{G\sim N(0,1)}\log\cosh(bG).
```
Under a uniform Boolean spin, every row field $`(Ax)_j`$ is exactly a
sum of $`n-1`$ independent signs, regardless of $`A`$.  For the
balanced $`k`$-fold core factor $`S_{A,k}`$, Jensen and the CLT give
```math
\boxed{
\liminf\frac1n\log S_{A,k}(b/\sqrt n)
\ge \log2-bc+k\ell(b).
} \tag{10.18}
```
A ground state independently gives the floor
$`k\log\cosh(2bc)`$.

Let $`s=bc`$, $`L=\log2`$, and
$`\delta_k=(1+k)^{3/2}-1-k^{3/2}`$.  The logarithmic union-bound
prefactor minus the entire scale-preserving allowance is therefore at
least
```math
\boxed{
\begin{aligned}
R_k(s)={}&(kL-sk^{3/2})_+\\
&+\max\{0,L-s+k\ell(2s),k\log\cosh(2s)\}
-s\delta_k.
\end{aligned}
} \tag{10.19}
```
An exact cusp analysis proves
```math
\boxed{R_k(s)>0\quad\text{for every integer }k\ge1,\ s\ge0.}
\tag{10.20}
```
For $`k\le4`$, retain the ground-state branch and minimize at
$`s=L/\sqrt k`$.  For $`k\ge5`$, use
```math
u-u^3\le\ell'(u)\le u,\qquad
\ell(u)\ge u^2/2-u^4/4;
```
the same cusp is the unique minimum and has positive value.  In the
self-doubling case the gap is exactly
```math
\boxed{
\log\cosh(2\log2)
-(2\sqrt2-2)\log2
=\log(17/8)-(2\sqrt2-2)\log2
=0.1795498767\ldots .
} \tag{10.21}
```
For $`k\to\infty`$, the limiting gap is at least
```math
2(\log2)^2-\tfrac12\log2
=0.6143324376\ldots .
```

The independent duplicate-row criterion with $`h=o(n)`$ also fails:
new-spin entropy plus the core-ground local-field factor exceeds its
available allowance by at least
```math
\boxed{
\min_{s\ge0}
\{\log2+\log\cosh(2s)-\tfrac32s\}
=\log8-\tfrac78\log7
=0.3767701613\ldots .
} \tag{10.22}
```
Thus the entire duplicate-row exponential-union-bound family is
stopped, at every fixed or diverging block ratio and in the
sublinear-insertion regime.  Any replacement must exploit
correlations or clusters among bad spin pairs instead of counting
them separately.

### 10.12 Exact macroscopic-closure dichotomy

At a closure step
```math
C_j=\begin{pmatrix}D_j&B_j\\B_j^\top&C_{j+1}\end{pmatrix},
```
write
```math
r_j^+=P(C_j)-P(D_j)-P(C_{j+1}),\qquad
r_j^-=N(C_j)-N(D_j)-N(C_{j+1}).
```
Both are nonnegative.  If $`g_j`$ is the inherited-state positive
gap, $`a_j`$ its deleted-block energy, and
$`\delta_j=P(D_j)-a_j`$, then the adaptive-closure identity sharpens
to
```math
\boxed{\frac12g_j-g_{j+1}=\delta_j+r_j^+.} \tag{10.23}
```
For a complete closure partition with terminal core $`T`$, the range
decomposes exactly:
```math
\boxed{
R(C_1)=\sum_jR(D_j)+R(T)+\sum_j(r_j^++r_j^-).
} \tag{10.24}
```

Let $`h_j^\pm`$ be the largest absolute cross correlation of
same-orientation ground states of the two children, and let $`V_j^\pm`$
be the corresponding maximum successor-ground $`\ell_1`$ fields.
The exact interaction excess
```math
\omega_j=
\max\!\left\{
2(h_j^++h_j^-),\,
\left[
\Phi(V_j^+,N(D_j))+\Phi(V_j^-,P(D_j))-R(D_j)
\right]_+
\right\}
```
satisfies
```math
\boxed{
R(C_1)\ge\sum_jR(D_j)+R(T)+\sum_j\omega_j.
} \tag{10.25}
```
Using the universal range theorem
$`R(A_m)\ge(2c_*-o(1))m^{3/2}`$, with
$`c_*=0.672986728863\ldots`$ in doubled normalization, a closure with
block proportions $`\alpha_j`$ and terminal proportion $`\rho`$
forces a strict improvement exactly when child normalized-range
excess plus $`\sum_j\omega_j/(2n^{3/2})`$ pays
```math
c_*\left(1-\sum_j\alpha_j^{3/2}-\rho^{3/2}\right).
```
For one balanced split this unresolved concavity budget is
```math
c_*(1-1/\sqrt2)=0.197113\ldots .
```

A scalable bent/Walsh construction suppresses more than one chosen
pair.  For both complete positive and negative child ground layers,
the cross fields satisfy
```math
\Phi(V_+,N(D))+\Phi(V_-,P(D))
\le n^{3/2}<R(D)=2n^{3/2},
```
so the full field excess $`\Delta`$ is exactly zero, while the parent
still has $`Q=O(n^{3/2})`$ and operator norm $`O(\sqrt n)`$.
Consequently every proof based only on positive/negative ground-layer
assembly is stopped.  The surviving version must use a non-ground
energy--field profile or a dependent multi-replica quantity.

### 10.13 Correlated energy-layer reduction and the iid cross-block floor

For
```math
G=\begin{pmatrix}A&B\\B^\top&D\end{pmatrix}
```
and exact energy layers
$`\mathcal X_A(p)=\{x:H_A(x)=p\}`$, define
```math
S_{p,q}(B)=
\max_{x\in\mathcal X_A(p),\,y\in\mathcal X_D(q)}
|x^\top By|.
```
There are only $`O(n^2m^2)`$ nonempty layer pairs, and exactly
```math
M(G)=\max_{p,q}\{|p+q|+S_{p,q}(B)\}. \tag{10.26}
```
If $`B`$ is iid Rademacher, bounded differences over this polynomial
family yield
```math
\boxed{
\mathbb E M(G)
=
\max_{p,q}\left\{
|p+q|+\mathbb E S_{p,q}(B)
\right\}
+o((n+m)^{3/2}).
} \tag{10.27}
```
Thus one can retain all correlations inside each energy layer; an
exponential state-by-state union bound is unnecessary.

For Gaussian $`B`$, the canonical metric is
```math
d((x,y),(x',y'))^2
=2nm-2(x\!\cdot x')(y\!\cdot y'),
```
and Sudakov--Fernique gives the constant-one comparison
```math
\mathbb E\sup_{x\in X,y\in Y}x^\top By
\le\sqrt m\,w(X)+\sqrt n\,w(Y). \tag{10.28}
```

This correlated route nevertheless fails for an iid cross block.
Every complete quadratic Hamiltonian has an exact
$`O(n)`$-energy layer of inverse-polynomial cube density.  Such a
layer has Gaussian width
```math
\left(\sqrt{2/\pi}+o(1)\right)n.
```
By switching invariance, a maximizer of an iid $`B`$ has uniform cube
marginals.  Moving those marginals to the dense exact layers changes
the bilinear value by only $`o((n+m)^{3/2})`$.  Hence, uniformly in
the child signings,
```math
\boxed{
\mathbb E M(G)
\ge
\mathbb E\|B\|_{\infty\to1}
-o((n+m)^{3/2}).
} \tag{10.29}
```

For balanced $`n\times n`$ iid $`B`$, two exact greedy updates give
```math
\mathbb E\|B\|_{\infty\to1}
\ge(\kappa-o(1))n^{3/2},
\qquad
\kappa=\mathbb E|Z+\sqrt{2/\pi}|
=1.0391966601\ldots .
```
After normalization by $`(2n)^{3/2}`$, the cross floor is
```math
0.3674115027\ldots>0.3364933644\ldots .
```
Therefore iid-cross chaining cannot preserve the current extremal
constant.  The only surviving block variant is a structured flat
cross matrix with both a sufficiently low rectangular norm and
uniform restricted-norm control on every high child-energy layer.

### 10.14 Two-sided switching-class width, not switching-minimality

Switch an absolute positive ground state to $`\mathbf1`$, so that in
half normalization $`H_A(\mathbf1)=M`$.  Every signed cut then obeys
```math
\boxed{0\le C_A(S)\le M.} \tag{10.30}
```
If $`G`$ is the graph of negative edges, this is exactly
```math
\frac{|S|(n-|S|)-M}{2}
\le e_G(S,S^c)
\le\frac{|S|(n-|S|)}2. \tag{10.31}
```
The upper edge inequality alone is the usual definition of a
switching-minimal graph; it omits the essential upper cap on signed
cuts.

That omission is fatal.  For even $`n`$, the maximum number of edges
in an ordinary switching-minimal graph is exactly
```math
\frac{n(n-2)}4,
```
by the singleton degree bound, with equality for
$`K_{n/2}\dot\cup K_{n/2}`$.  Equivalently, the balanced rank-one
signing $`a_{ij}=-u_iu_j`$ has positive maximum only $`n/2`$, while
its negative energy is $`-\binom n2`$.  Thus one-sided
switching-minimality forces only an $`O(n)`$ deficit and cannot prove
the desired $`n^{3/2}`$ scale.

For odd $`n`$, the exact maximum is $`(n-1)^2/4`$, so the minimum
signed total is $`(n-1)/2`$.  Indeed, the nonnegative row sums are
even; the zero-row vertices form a clique by the two-cut inequality
and number at most $`(n+1)/2`$, forcing total row sum at least $`n-1`$.
The two-clique rank-one construction of sizes $`(n\pm1)/2`$ attains
equality.

The original problem is instead the minimum width of the edge-count
distribution within a Seidel switching class:
```math
\frac12\binom n2-\frac M2
\le |E(G^S)|
\le\frac12\binom n2+\frac M2
\quad\text{for every }S. \tag{10.32}
```

There are two useful exact recursions.  Delete a vertex $`v`$, let
```math
r_v=n-1-2d_G(v),
```
and let $`\Delta`$ be the edge reduction obtained by switching the
induced graph to a switching-minimal representative.  Comparing the
two complementary full switches gives
```math
\boxed{0\le\Delta\le r_v/2,} \tag{10.33}
```
so the repaired induced signed total lies in $`[M-r_v,M]`$.
For two cuts,
```math
\boxed{
D(S,T)=
\frac{C_A(S)+C_A(T)-C_A(S\triangle T)}2
\in[-M/2,M],
} \tag{10.34}
```
where $`D(S,T)`$ is the signed sum of edges crossing both cuts.  Every
weighted partition quotient inherits the same two-sided cap.

These identities do not yet improve the constant.  Fixed-degree
cut-cone/SOS localization yields only
$`-Q I/2\preceq A\preceq Q I/2`$; growing degree returns the original
covering-radius problem.  The surviving graph-theoretic target must
retain both sides of (10.30) across a growing correlated family of
cuts.

### 10.15 Nonnegative quadratic Fourier cone and signed-Johnson no-go

With doubled energy $`q_A(x)=x^\top A x`$ and
$`R=\max_x|q_A(x)|`$, switch an absolute positive ground state to
$`\mathbf1`$.  Then
```math
f_\pm(x)=\frac{R\pm q_A(x)}4\ge0
```
are nonnegative degree-two Walsh polynomials with constant coefficient
$`R/4`$ and every nonconstant coefficient of magnitude $`1/2`$.  If
$`H`$ is the graph of positive edges, this cone is exactly the class of
s-maximal graphs together with the essential opposite-endpoint cap
```math
\boxed{
0\le e_H(S,S^c)-\frac12|S||S^c|\le\frac R4
\quad\text{for all }S.
} \tag{10.35}
```
In particular, the original minimax problem is a two-sided width problem
inside a Seidel switching class, not merely one-sided switching
minimality.

Restricting the Fourier-convolution PSD matrix of $`f_\pm`$ to
$`k`$-subsets gives the signed Johnson, or hard-core-boson, operator
```math
(T_k)_{S,S-\{i\}+\{j\}}=a_{ij},\qquad
-\frac R2I\preceq T_k\preceq\frac R2I.
```
Its first two exact moments are
```math
\frac{\operatorname{tr}T_k^2}{\binom nk}=k(n-k)
```
and
```math
\boxed{
\frac{\operatorname{tr}T_k^3}{\binom nk}
=6\frac{k(n-k)}{n(n-1)}
\sum_{i<j<\ell}a_{ij}a_{j\ell}a_{\ell i}.
} \tag{10.36}
```
Even the endpoint inequality
$`\lvert\operatorname{tr}T_k^3\rvert\le(R/2)\operatorname{tr}T_k^2`$
only produces an $`O(n)`$ lower bound for $`R`$.  More generally, every
fixed-order trace or fixed subgraph statistic is blind at the
$`n^{3/2}`$ scale: (10.35) with $`R=o(n^2)`$ already makes the ordinary
graphon limit constant $`1/2`$.

The tempting comparison with the fermionic exterior-power lift also
fails under frustration.  Direct Paley-conference benchmarks make the
failure quantitative.  At half filling, the lower certificate
$`2\|T_k\|/n^{3/2}`$ is $`0.57646,0.60407,0.61172`$ for
$`n=6,14,18`$, respectively, all below the proved doubled constant
$`0.6729867\ldots`$.  Thus the bare many-particle norm cannot improve
the current bound even on the canonical flat-spectrum examples.

Some exact endpoint structure remains useful.  If $`\omega(H)`$ is the
clique number of the s-maximal representative, then
```math
R\ge\omega(H)(\omega(H)-1).
```
A triangle-free s-maximal graph is balanced complete bipartite and its
opposite orientation has quadratic norm $`n(n-1)`$.  Therefore a viable
SOS/cut-cone route must be a genuinely growing-order fluctuation
dichotomy: correlated cycles must either yield a Boolean witness
directly or force proximity to a switching extreme.  Degree-one PSD,
any fixed trace hierarchy, a plain boson--fermion comparison, and the
bare $`T_k`$ norm are now stopped as standalone routes.

### 10.16 Two-sided multicut hierarchies require linear entropy

Let $`\Gamma\le\{\pm1\}^n`$ be a subgroup of rank $`k`$, represented by
vertex labels $`\lambda_i\in\mathbb F_2^k`$.  For the edge-difference
fibers
```math
W_d=\sum_{\lambda_i+\lambda_j=d}a_{ij},
```
the complete energy profile on $`\Gamma`$ is exactly its Walsh transform:
```math
H_A(g_t)=\sum_dW_d(-1)^{t\cdot d},\qquad
2^{-k}\sum_tH_A(g_t)^2=\sum_dW_d^2. \tag{10.37}
```
Thus every symmetric-difference and two-cut identity among the selected
cuts is already present in this quotient.

This entire hierarchy has a sharp entropy obstruction.  For every
prescribed $`\Gamma`$, there is an actual flat signing $`B`$, switched so
that $`H_B(\mathbf1)=M_\Gamma`$, with
```math
\boxed{
|H_B(g)|\le M_\Gamma\quad(g\in\Gamma),\qquad
M_\Gamma\le\sqrt{2\binom n2\log(48|\Gamma|)}.
} \tag{10.38}
```
The same example can retain at least half of the quotient Parseval
energy, $`\sum_dW_d^2\ge\binom n2/2`$.  The proof is a direct random
signing union bound combined with Paley--Zygmund, followed by switching
at a maximizing group element.  Consequently, when $`k=o(n)`$, the
complete two-sided cut cap on $`\Gamma`$, all its triangle consequences,
and flat second moments are jointly compatible with
$`M_\Gamma=o(n^{3/2})`$.

The averaged identity is equally definitive:
```math
\mathbb E_\lambda\sum_dW_d^2
=(1-2^{-k})\binom n2+2^{-k}H_A(\mathbf1)^2. \tag{10.39}
```
For every $`k`$, combining it with the cap recovers only
$`M^2\ge\binom n2`$.  More generally, averaging any test function over
random labels merely reproduces the ordinary Rademacher-chaos law; it
does not create a new moment inequality.

There is also a deterministic local obstruction.  For balanced
$`u\in\{\pm1\}^n`$, the rank-one signing
$`a_{ij}=-u_i u_j`$ satisfies, after switching,
```math
M=n/2,\qquad
C_A(S)=\left(\sum_{i\in S}u_i\right)^2.
```
It therefore obeys the full two-sided cap for every
$`|S|\le\sqrt{n/2}`$, including all singleton and pair cuts, while the
remote cut $`x=u`$ has quadratic energy of order $`n^2`$.  Fixed-degree
and even $`o(\sqrt n)`$-radius local tests cannot see the obstruction.

Finally, if $`T`$ is a uniformly random $`t`$-vertex restriction, with
$`p_s=(t)_s/(n)_s`$, signed total $`M=H_A(\mathbf1)`$, and switched row
sums $`r_i`$, then
```math
\boxed{
\mathbb E H_{A[T]}(\mathbf1)^2
=p_4M^2+(p_2-p_4)\binom n2
 (p_3-p_4)\left(\sum_i r_i^2-n(n-1)\right).
} \tag{10.40}
```
The deletion repair can replenish the entire scalar loss, and the
rank-one obstruction has only $`\sum_i r_i^2=O(n)`$.  Hence scalar
restriction moments do not close the recursion.

The conclusion is sharp: a multicut, SOS, or deletion proof must retain
$`\exp(\Omega(n))`$ genuinely remote configurations, or information
equivalent to a linear-entropy energy profile.  Another finite quotient
or fixed moment cannot prove the $`n^{3/2}`$ scale, much less
convergence.

### 10.17 Flat anti-conjugate blocks: exact ellipse and failure of closure

Let $`B=\sqrt n\,U`$ be Hadamard, and suppose the two child signings
satisfy the exact anti-conjugacy relation
```math
D=-U^\top A U.
```
For Boolean $`x,y`$, put $`z=Uy`$ and
$`t=|x^\top z|/n`$.  Exact polarization gives the layer ellipse
```math
\boxed{
|H_A(x)+H_D(y)|
\le n\|A\|_{\rm op}\sqrt{1-t^2},\qquad
|x^\top By|=n^{3/2}t.
} \tag{10.41}
```
Equivalently,
```math
S_{p,q}(B)\le n^{3/2}
\sqrt{1-\left(\frac{|p+q|}{n\|A\|_{\rm op}}\right)^2},
```
and the full parent obeys
```math
M\!\begin{pmatrix}A&B\\B^\top&D\end{pmatrix}
\le n\sqrt{\|A\|_{\rm op}^2+n}. \tag{10.42}
```
An operator-norm perturbation $`E=D+U^\top A U`$ costs at most
$`(n/2)\|E\|_{\rm op}`$.  This completely prevents addition of the
internal and cross spectral costs, but the Frobenius floor
$`\|A\|_{\rm op}\ge\sqrt{n-1}`$ leaves the certificate at the
half-normalized constant $`1/2`$.

Hadamard cross blocks also satisfy the sharp conditional moment bound
```math
\mathbb E_x\exp\!\left(\frac t{\sqrt n}x^\top By\right)
\le(\cosh t)^n. \tag{10.43}
```
Independent row/column switching therefore yields an exact two-sided
free-energy recursion, but its optimized dyadic limit is
$`0.76668975\ldots`$, far above $`1/2`$.  Random relative orientation
without child--basis correlation is stopped.

Exact sign anti-conjugate pairs nevertheless form an infinite
phase-space family.  For $`V=\mathbb F_2^d`$, Walsh matrix
$`W_{u,v}=(-1)^{u\cdot v}`$, any permutation
$`\pi:V\setminus\{0\}\to V\setminus\{0\}`$ satisfying
$`r\cdot\pi(r)=0`$, and signs $`\varepsilon_r`$, define
```math
A_{u,v}=\varepsilon_{u+v}(-1)^{v\cdot\pi(u+v)}
\quad(u\ne v).
```
Then $`W^\top A W/n`$ is again symmetric, zero-diagonal, and
entrywise signed, so $`D=-W^\top A W/n`$ gives exact anti-conjugacy.
For every even $`d`$, $`\pi(r)=Jr`$ with invertible alternating $`J`$
provides such examples.

The small-order audit is decisive.  At $`d=2`$, every resulting
order-$`8`$ parent has $`M=10=F(8)`$, normalized
$`0.4419417\ldots`$.  At $`d=3`$, all $`3072`$ admissible parents have
$`M=32`$, normalized exactly $`1/2`$.  Moreover, the good order-$`8`$
parent has only twelve zero-energy Boolean rays and no Hadamard basis
among them satisfying the transformed off-diagonal sign condition.
Thus the isolated optimal lift cannot be iterated.

The sole surviving flat-block lemma is now a forbidden-band problem:
prove that Boolean pairs in an exact or approximate anti-conjugate sign
family cannot approach the support point of the continuous ellipse.
Without such a discrete improvement, anti-conjugacy is not a
scale-preserving recursion below $`1/2`$.

There is now a complete resonance theorem for the easiest infinite
subfamily.  If $`\pi(r)=Jr`$ with $`J`$ invertible alternating, the
associated Weyl summands commute.  A Boolean quadratic-phase basis
$`\xi_s(v)=(-1)^{q(v)+s\cdot v}`$, where $`q`$ has polar form $`J`$,
diagonalizes $`A`$.  Writing
$`\eta_0=0`$ and
$`\eta_r=\varepsilon_r(-1)^{q(r)}`$, one obtains exactly
```math
\boxed{
M(A)=\frac n2\max_s|\widehat\eta(s)|
\ge\frac12n^{3/2}\sqrt{1-\frac1n}.
} \tag{10.44}
```
Hence every linear exact anti-conjugate family is completely resonant
and asymptotically trapped at $`1/2`$.  A candidate below $`1/2`$ must
use nonlinear $`\pi`$, equivalently noncommuting Weyl operators.  For
two distinct displacements $`r,s`$, the two ordered contributions to
$`A^2`$ cancel exactly when
```math
r\cdot\pi(s)+s\cdot\pi(r)=1. \tag{10.45}
```
The remaining phase-space problem is therefore a precise
collision--resonance tradeoff: linear graphs have the cancellations
needed for flat spectra but Boolean eigenvectors, while generic
nonlinear graphs lose those eigenvectors but revert to Wigner-scale
fourth moments.

### 10.18 Finite arithmetic obstructions cannot cause nonconvergence

The normalized sequence is uniformly continuous under sublinear
padding.  For arbitrary $`n,h`$, random cross edges between optimal
children give
```math
\boxed{
F(n+h)\le F(n)+F(h)
+\sqrt{2nh(n+h+2)\log2}.
} \tag{10.46}
```
Together with monotonicity and the random-sign bound
```math
F(h)\le
\sqrt{2\binom h2(h+2)\log2},
```
this proves
```math
\boxed{
h=o(n)\quad\Longrightarrow\quad
\frac{F(n+h)}{(n+h)^{3/2}}-\frac{F(n)}{n^{3/2}}\to0.
} \tag{10.47}
```
More generally, every subsequence whose next-element gap is $`o(n)`$
has exactly the same cluster set as the full normalized sequence.
Therefore parity, every fixed residue class, every fixed $`2`$-adic
class, and any other $`o(n)`$-gap arithmetic family cannot carry a
different limiting constant.  Nonconvergence, if real, must occur on
multiplicatively separated scales.

The exact code arithmetic explains why the usual bent/semi-bent
analogy is misleading here.  For the augmented cut code
$`\mathcal C_n`$, the gcd of nonzero codeword weights is
```math
\boxed{
\Delta_n=
\begin{cases}
2,&n\equiv1\pmod4,\\
1,&\text{otherwise}.
\end{cases}
} \tag{10.48}
```
Its dual consists of even-cardinality Eulerian subgraphs and has
minimum distance $`4`$.  In contrast, first-order Reed--Muller codes
have divisibility growing proportionally to block length, which is what
supports their parity-dependent bent obstruction.  Here the universal
energy lattice has only bounded spacing: one residue modulo $`4`$ for
odd $`n`$, at most two for even $`n`$; triangle traces have spacing
$`12`$.  All disappear after $`n^{-3/2}`$ normalization.

Seidel characteristic-polynomial parity likewise depends only on the
ordinary parity of $`n`$:
```math
\det(\lambda I-A)\equiv
\begin{cases}
(\lambda+1)^n,&n\text{ even},\\
\lambda(\lambda+1)^{n-1},&n\text{ odd}
\end{cases}\pmod2.
```
Finally, principal submatrices of Paley conference matrices of order
$`n+o(n)`$ give the spectral $`1/2+o(1)`$ upper bound in every
congruence class.  Exact conference nonexistence or determinant
integrality cannot force a leading class-dependent gap.

This closes the finite-arithmetic branch.  The only surviving
nonconvergence mechanism would be a genuinely geometric oscillation
across multiplicative scales, not congruence or divisibility.

### 10.19 Signing-space entropy gives a genuine convergence criterion

Let
```math
Z_n(c)=
\#\{A:M(A)\le c n^{3/2}\}.
```
Switching and global negation imply that $`Z_n(c)`$ is either zero or
at least $`2^n`$.  More importantly, the edge-Hamming Lipschitz bound
```math
|M(A)-M(B)|\le2d_H(A,B)
```
gives a refined thickening scale.  If $`Z_n(c)>0`$, then for every
fixed $`\varepsilon>0`$,
```math
\boxed{
\log Z_n(c+\varepsilon)
\ge
\left(\frac{\varepsilon}{4}+o(1)\right)n^{3/2}\log n.
} \tag{10.49}
```
Consequently, if the refined microcanonical entropies
```math
\Sigma_n(c)=
\frac{\log(1+Z_n(c))}{n^{3/2}\log n}
```
converge pointwise for fixed $`c`$ throughout the bounded interval
containing the cluster set of $`F(n)/n^{3/2}`$, then the normalized
minimum converges.  A liminf subsequence thickens to a positive limiting
entropy at $`c+\varepsilon`$; pointwise convergence then forces
nonemptiness at every sufficiently large order.

There is an even cleaner canonical criterion.  Define the Gibbs
partition function on the space of edge signings
```math
\mathfrak Z_n(\beta)=
\sum_Ae^{-\beta\sqrt n\,M(A)},\qquad
\Phi_n(\beta)=\frac1{n^2}\log\mathfrak Z_n(\beta).
```
Writing $`f_n=F(n)/n^{3/2}`$, the exact entropy squeeze is
```math
\boxed{
-\frac{\Phi_n(\beta)}\beta
\le f_n\le
-\frac{\Phi_n(\beta)}\beta+
\frac{\binom n2}{n^2}\frac{\log2}{\beta}.
} \tag{10.50}
```
Therefore, if $`\Phi_n(\beta)`$ converges for any unbounded set of fixed
$`\beta`$'s, then
```math
\limsup f_n-\liminf f_n\le\frac{\log2}{2\beta}
```
for each such $`\beta`$, and the requested limit follows by
$`\beta\to\infty`$.  This avoids an exchange of the thermodynamic and
zero-temperature limits.

The canonical recursion has an exact one-row cavity form.  For an
order-$`n`$ core $`B`$ and new row $`b`$,
```math
M(B,b)=\max_x\bigl(|H_B(x)|+|b\cdot x|\bigr),
\qquad
\Delta_B(b)=M(B,b)-M(B).
```
Set
```math
R_{B,n}(\lambda)=\sum_b e^{-\lambda\Delta_B(b)},\quad
\beta'=\beta\sqrt{\frac{n+1}{n}},\quad
\lambda=\beta\sqrt{n+1}.
```
Then
```math
\boxed{
\mathfrak Z_{n+1}(\beta)=
\mathfrak Z_n(\beta')\,
\mathbb E_{\nu_{n,\beta'}}R_{B,n}(\lambda),
} \tag{10.51}
```
where $`\nu_{n,\beta'}`$ is the signing Gibbs measure.  Thus convergence
has been reduced to the row-cavity free energy
```math
\frac1n\log
\mathbb E_{\nu_{n,\beta'}}
R_{B,n}(\beta\sqrt{n+1}).
```
Closure in the scalar $`M(B)`$ is false: two explicit order-$`6`$
signings with $`M=9`$ have respectively
```math
R_1=36e^{-2\lambda}+24e^{-4\lambda}+4e^{-6\lambda}
```
and
```math
R_2=8+40e^{-2\lambda}+14e^{-4\lambda}+2e^{-6\lambda}.
```
The second has eight norm-preserving rows; the first has none.

The exact one-step cavity state can be identified geometrically.  Put
```math
g_B(x)=M(B)-|H_B(x)|,\qquad
\Delta_B(b)=\max_x\bigl(|b\cdot x|-g_B(x)\bigr),
```
and
```math
V_B(u)=|\{b:\Delta_B(b)\le u\}|.
```
Since $`|b\cdot x|=n-2d_\pm(b,x)`$, where
$`d_\pm(b,x)=\min\{d_H(b,x),d_H(b,-x)\}`$,
```math
\boxed{
\{\Delta_B\le u\}
=
\bigcap_x
\left\{
b:d_\pm(b,x)\ge\frac{n-g_B(x)-u}{2}
\right\}.
} \tag{10.51a}
```
Thus $`V_B`$ is an intersection-of-cap profile carrying the complete
overlap geometry of all near-ground configurations.  At
$`\lambda=\beta\sqrt n`$, its Legendre transform determines the row
factor up to $`o(1)`$:
```math
\max_d\left\{\frac1n\log V_B(d)-\beta\frac d{\sqrt n}\right\}
\le\frac1n\log R_{B,n}(\lambda)
\le
\max_d\left\{\frac1n\log V_B(d)-\beta\frac d{\sqrt n}\right\}
+\frac{\log(n+1)}n. \tag{10.51b}
```
Two explicit order-$`8`$ signings have the same norm $`16`$ and the
same complete $`|H_B|`$-histogram, but different row-increment
histograms:
```math
\begin{array}{c|rrrrr}
\Delta&0&2&4&6&8\\ \hline
B_1&60&110&66&18&2\\
B_2&68&112&58&16&2.
\end{array}
```
So even the full unlabeled energy histogram is not a closed cavity
state; the missing information begins with layer overlaps.

Switching invariance gives a useful gauge reduction:
```math
\Delta_{B^s}(b)=\Delta_B(bs),\qquad
\boxed{
\mathbb E_\nu R_{B,n}(\lambda)
=2^n\mathbb E_\nu e^{-\lambda\Delta_B(\mathbf1)}.
} \tag{10.51c}
```
In this gauge, the signed energy-versus-magnetization profile
```math
U_B(m)=\max_{\sum x_i=m}H_B(x),\qquad
L_B(m)=\min_{\sum x_i=m}H_B(x)
```
has an exact max-plus/min-plus transition under a universal-positive
vertex extension:
```math
\begin{aligned}
U_{T(B)}(q)&=\max_{\substack{s=\pm1\\m=q-s}}
\bigl(U_B(m)+sm\bigr),\\
L_{T(B)}(q)&=\min_{\substack{s=\pm1\\m=q-s}}
\bigl(L_B(m)+sm\bigr).
\end{aligned} \tag{10.51d}
```
The obstruction is that the next cavity step averages over every
gauge, and hence requires the whole switching orbit of this profile,
equivalently the higher cap intersections in (10.51a).

Ordinary concentration is quantitatively too weak.  A single core-edge
flip changes every $`\Delta_B(b)`$ by at most $`4`$, whence
```math
|\Delta\log R_{B,n}(\beta\sqrt n)|\le4\beta\sqrt n,
\qquad
\operatorname{Var}(\log R)=O(\beta^2n^3)
```
under the product signing law.  The target variable has scale $`n`$,
so the edgewise variance proxy misses self-averaging by a factor $`n`$.
Under the signing Gibbs law, conditional edge odds may vary by
$`\exp(\Theta(\beta\sqrt n))`$, so a dimension-free Dobrushin argument
is unavailable.

There is an exact hereditary inequality, but it follows the wrong
temperature characteristic.  If $`E_n=\binom n2`$, then
```math
\boxed{
\frac{\log Z_N(t)}{E_N}\le\frac{\log Z_m(t)}{E_m},
\qquad
\frac{\log\mathfrak Z_N(\lambda)}{E_N}
\le
\frac{\log\mathfrak Z_m(\lambda E_m/E_N)}{E_m}.
} \tag{10.51e}
```
At $`\lambda=\beta\sqrt N`$, the induced smaller-order inverse
temperature is
$`\beta_m=\beta(m/N)^{3/2}(1+o(1))`$, rather than fixed $`\beta`$.
This is a substantive barrier: the abstract thresholds
```math
T_n=n^{3/2}\bigl(c_0+\varepsilon\sin(\log\log n)\bigr)
```
can be chosen nondecreasing with $`0\le T_{n+1}-T_n\le n`$.  The
profiles $`\widehat Z_n(t)=0`$ below $`T_n`$ and $`2^{E_n}`$ above it
satisfy hereditary Shearer, switching-orbit quantization, every
Hamming-thickening lower bound, and the one-step count sandwich, while
$`T_n/n^{3/2}`$ oscillates.  They even match the universal scalar
cavity range exactly through
```math
\widehat R_n(\lambda)
=2^ne^{-\lambda(T_{n+1}-T_n)}
\in[2^ne^{-\lambda n},2^n].
```
Therefore scalar entropy, scalar increments, Shearer, and cavity-factor
magnitude cannot jointly prove convergence.  A proof must use a
speed-$`n`$ law for the cap-intersection profile, or an equivalent
non-scalar geometric state.

A speed-$`n^2`$ lower-tail LDP alone is insufficient.  A single
switching orbit already has the maximal finite probability rate
$`(\log2)/2`$, so the leading rate cannot distinguish emptiness from a
sparse algebraic feasible phase.  The refined entropy or the full
canonical pressure is essential.

Exact counting transfer under adding $`h`$ vertices is available, but
at proportional scale it pays
```math
\sqrt{\log2\,(1-\alpha^2)}
```
in the normalized threshold when $`n/(n+h)\to\alpha`$.  Hence local
entropy regularization is proved; the missing theorem is convergence of
the cavity pressure, or an equivalent cross-block anti-alignment
principle at multiplicative scales.

### 10.20 Nonlinear phase space: collision defect versus Boolean resonance

The exact Walsh anti-conjugate family admits a Weyl-operator
description.  For $`V=\mathbb F_2^d`$, $`n=2^d`$, and an orthogonal
permutation $`\pi`$ with $`r\cdot\pi(r)=0`$, write
```math
A=\sum_{r\ne0}\varepsilon_rP_r,\qquad
w_r=(r,\pi(r)).
```
Two summands commute exactly when their symplectic product vanishes.
If $`\mathcal P_t`$ is the set of commuting unordered pairs whose
phase-space sum is $`t`$, then
```math
\boxed{
\|A^2-(n-1)I\|_F^2
=4n\sum_t
\left(
\sum_{\{r,s\}\in\mathcal P_t}
\sigma_{rs}\varepsilon_r\varepsilon_s
\right)^2.
} \tag{10.52}
```
This is the exact signed collision equation that any nonlinear
near-conference construction must solve.

Finite geometry forces many collisions before signs are chosen.  The
points $`w_r`$ lie on the hyperbolic quadric $`Q^+(2d-1,2)`$.  The
least-eigenvalue bound in its collinearity graph yields
```math
\boxed{
C(\pi)\ge
\frac{(n-1)(n+4)(n-4)}{8(n+2)}
=\left(\frac18+o(1)\right)n^2,
} \tag{10.53}
```
where $`C(\pi)`$ is the number of commuting pairs.  Hence
```math
\mathbb E_\varepsilon
\|A^2-(n-1)I\|_F^2=4nC(\pi)
```
and the expected normalized fourth moment is at least
$`3/2-o(1)`$, rather than the conference/Haar value $`1`$.

There is an explicit fully nonlinear obstruction.  Identify $`V`$ with
$`\mathbb F_{2^d}`$, choose nonzero $`\alpha`$ with
$`\operatorname{Tr}\alpha=0`$, and set
```math
\pi(r)=\alpha/r.
```
The graph $`\{(r,\alpha/r)\}`$ is Sidon: equality of two distinct pair
sums fixes both the sum and product of the underlying field elements.
Thus every commuting pair occupies a singleton class and no choice of
$`\varepsilon`$ can cancel it.  A Kloosterman count gives
```math
C(\pi)=\left(\frac14+O(n^{-1/2})\right)n^2,\qquad
\boxed{
\frac{\operatorname{tr}A^4}{n(n-1)^2}
=2+O(n^{-1/2})
} \tag{10.54}
```
for every coefficient signing.  This natural nonlinear family is
irreducibly Wigner-like.

At the opposite extreme, an alternating linear map
$`\pi(r)=Lr`$ becomes a Cayley convolution after a quadratic vertex
switch.  Walsh characters are Boolean eigenvectors and
```math
\boxed{
Q(A)=n\max_k|\widehat a(k)|
\ge n\sqrt{n-1}.
} \tag{10.55}
```
Bent choices match this asymptotically, so the optimized original
one-copy constant in the entire linear family tends exactly to
$`1/2`$.  Linear structure supplies spectral flatness but forces
Boolean resonance.

Ordinary additive inverse theory cannot bridge the intermediate
regime.  Small conference defect forces only
$`\Theta(n^2)`$ unsigned additive energy for a size-$`n`$ set in a
group of size $`n^2`$.  In Balog--Szemerédi--Gowers notation this is
$`K=\Theta(n)`$, so the extracted structured subset is nonmacroscopic;
the loss is a full factor $`n`$.  A viable inverse theorem would have
to exploit the simultaneous *signed* cancellation equations in
(10.52), not unsigned parallelogram counts.  This closes the generic
nonlinear and affine endpoints but leaves a sharply defined signed
intermediate problem.

### 10.21 Exact transfer states and why scalar/projective dynamics fail

Gauge fixing the newest vertex gives an exact rooted tree of switching
classes.  An order-$`n`$ class $`[A]`$ has $`2^{n-1}`$ children
$`[T(A^b)]`$, indexed by $`b\in\{\pm1\}^n/\{\pm\mathbf1\}`$, and
```math
M(T(A^b))
=\max_x\bigl(|H_A(x)|+|b\cdot x|\bigr). \tag{10.56}
```
The signed magnetization envelope closes for the distinguished
all-positive child by the max-plus/min-plus recursion (10.51d), but it
is not Markov under the full gauge branching.  Two explicit order-$`7`$
signings have the same norm $`13`$ and the same complete signed
energy-versus-magnetization profile, yet their child-increment
multisets are respectively
```math
\begin{array}{c|rrrr}
\Delta&1&3&5&7\\ \hline
A_1&70&42&14&2\\
A_2&40&60&24&4.
\end{array} \tag{10.57}
```
The smallest exact closed support state is
```math
\Phi_A^\sigma(h)=\max_x\bigl(\sigma H_A(x)+h\cdot x\bigr),
```
with transfer
```math
\Phi_{A_b}^\sigma(h,t)
=\max_{u=\pm1}\{tu+\Phi_A^\sigma(h+\sigma ub)\}. \tag{10.58}
```
This state is injective: it recovers every $`H_A(x)`$ from a boundary
limit.  The transfer is order preserving and additively homogeneous,
but its Hilbert-projective contraction coefficient is exactly $`1`$:
on a boundary face one branch is selected and old differences are
copied isometrically.  Thus an ordinary nonlinear
Perron--Frobenius/unique-attractor proof cannot establish convergence.

The same obstruction appears in the canonical signing pressure.
For
```math
Y_n(\lambda)=2^{-\binom n2}\sum_Ae^{-\lambda M(A)}
```
one has the exact raw-parameter monotonicity
```math
e^{-\lambda n}Y_n(\lambda)\le Y_{n+1}(\lambda)\le Y_n(\lambda),
\tag{10.59}
```
but the diagonal $`\lambda=\beta\sqrt n`$ can still oscillate in the
abstract scalar countermodel following (10.51e).  Scalar pressure,
rooted scalar increments, and projective contraction are therefore
closed routes.

### 10.22 Correct fluctuation state and the finite-traffic no-go

For each signing, retain the support of every finite overlap-energy
tuple of Boolean rank-one projectors $`P_x=xx^\top/n`$.  The resulting
projective support state is compact and contains
```math
\frac{M(A)}{n^{3/2}}
=\frac12\max_x\left|\operatorname{Tr}
\left(\frac A{\sqrt n}P_x\right)\right|
```
continuously on bounded-energy sets.  A full speed-$`n^2`$ LDP on this
state would force convergence: a single signing has finite maximal
rate $`(\log2)/2`$, so an LDP lower bound propagates any liminf-state
neighborhood to every sufficiently large order.

This state cannot be replaced by graphons, spectral moments, or fixed
replicas.  For Walsh orders $`n=2^d`$, start from the off-diagonal
Walsh signing $`A_n`$, for which
```math
\frac{M(A_n)}{n^{3/2}}\le\frac12+o(1).
```
Flip a random controlled subset of
$`(1/3+o(1))n^{3/2}`$ negative edges.  A deterministic realization
$`B_n`$ satisfies
```math
\frac{M(B_n)}{n^{3/2}}\ge\frac23+o(1),\qquad
\left\|\frac{A_n-B_n}{\sqrt n}\right\|_{2,\tau}=o(1),
```
while both normalized operator norms remain bounded.  Every fixed
trace polynomial, fixed signed graph density, and fixed-replica
normalized energy law has the same limit for $`A_n`$ and $`B_n`$.
The planted Boolean direction is invisible to all such finite
traffic states.

### 10.23 Proportional restriction: a conditional scale-transfer theorem

The most direct surviving convergence mechanism is optimized
principal restriction.  If every competitive order-$`N`$ signing had,
for each fixed $`0<\alpha<1`$, a subset $`|S|=\alpha N+o(N)`$ with
```math
M(A[S])\le\alpha^{3/2}M(A)+o(N^{3/2}), \tag{10.60}
```
then a liminf sequence could be transferred to every multiplicative
scale, proving convergence.

The exact finite inequality is false, but conference data support the
$`o(N^{3/2})`$ version.  For symmetric conference orders $`6,14,18`$,
the best half-subset excess over the target is exactly the unavoidable
energy-lattice rounding; the normalized excess decreases rapidly.

There is a rigorous energy-layer criterion.  Let $`K=M(A)`$,
```math
\mathcal L_t^\pm=\{x:\pm H_A(x)\ge t\},\qquad
\delta=\frac{L-t}{K-t},
```
with $`0\le t<L<K`$.  If $`m<N`$ and
```math
|\mathcal L_t^+|+|\mathcal L_t^-|
<\delta\,2^{N-m}, \tag{10.61}
```
then some $`|S|=m`$ has $`M(A[S])<L`$.  More sharply, if $`D_\pm(L)`$
is the maximum, over a full spin $`x`$, of the fraction of $`m`$-sets
on which its restricted energy has the indicated sign and exceeds
$`L`$, it is enough that
```math
D_+(L)|\mathcal L_t^+|+
D_-(L)|\mathcal L_t^-|
<\delta\,2^{N-m}. \tag{10.62}
```
Both statements follow by averaging a bad restricted ground state over
all full extensions and double-counting $`(S,x)`$.

Raw layer entropy cannot prove (10.61).  For every homogeneous
quadratic $`f`$, $`K=\|f\|_\infty`$, and $`0<\theta<1`$, noise around a
maximizer gives
```math
\liminf_{n\to\infty}\frac1n\log_2
|\{x:|f(x)|\ge\theta K\}|
\ge
h_2\!\left(\frac{1-\sqrt\theta}{2}\right)
\ge1-\theta. \tag{10.63}
```
Thus exponentially large near-maximal Hamming clouds are universal.
After quotienting those clouds, a two-threshold packing theorem gives
a subcritical center rate, but a sharp entropy inequality shows that
center count times worst extension-cylinder intersection can still
fill the whole cylinder.  The multiplicity-weighted criterion
(10.62), not raw layer or cluster entropy, is the correct target.

### 10.24 Growing-degree triangle certificates and maximal hafnian cancellation

For every $`1\le k\le n/2`$, Parseval applied to $`H_A^k`$ gives a
genuine degree-$`\Theta(n)`$ certificate.  If $`|S|=2k`$, its top
Fourier coefficient is
```math
\widehat{H_A^k}(S)=k!\operatorname{haf}(A_S).
```
Because a signed hafnian is a sum of an odd number
$`(2k-1)!!`$ of signs, it never vanishes.  Hence
```math
\boxed{
M(A)\ge
\left[\binom n{2k}(k!)^2\right]^{1/(2k)}.
} \tag{10.64}
```
Optimizing (10.64) gives only order $`n`$, not $`n^{3/2}`$.

This loss is sharp for the entire top-Fourier/parity route.  On
vertices $`L_1,R_1,\ldots,L_k,R_k`$, take every edge positive except
```math
a_{R_rL_s}=-1\qquad(r<s).
```
Then
```math
\boxed{\operatorname{haf}(A)=1\quad\text{for every }k.} \tag{10.65}
```
An elementary recurrence expands at the final pair; the two cross
sums cancel under the reverse $`L\leftrightarrow R`$ involution, so
the hafnian remains its initial value $`1`$.  This was also verified
independently by a creation-annihilation calculation.  Thus induced
signed matchings can cancel all the way down to parity at linear
degree, although fixed $`k=2`$ anti-cancellation is strong by
Goodman's theorem.

### 10.25 Planted signing entropy and why fixed replicas still miss the hole

One good signing generates a rigorous speed-$`n^2`$ family of nearby
good signings.  If $`M(A_0)\le c_0n^{3/2}`$ and every edge is flipped
independently with probability $`\delta`$, then a moderate-deviation
union bound gives, with high probability,
```math
\frac{M(A_0\odot\xi)}{n^{3/2}}
\le
c_\delta+o(1),\qquad
c_\delta=(1-2\delta)c_0+
2\sqrt{\delta(1-\delta)\log2}. \tag{10.66}
```
Typical-shell counting consequently yields
```math
\liminf\frac1{n^2}\log Z_n(c_\delta+\varepsilon)
\ge\frac12h(\delta).
```
The function $`c_\delta`$ is concave on $`[0,1/2]`$, so it cannot
improve the seed constant.

For two oriented cut directions whose spin overlap is $`q`$, their
edge-direction correlation is
```math
\rho=\pm\frac{nq^2-1}{n-1}.
```
The joint upper-tail exponent is
```math
\Pr(A\cdot u,A\cdot v\ge cn^{3/2})
=
\exp\left\{-\frac{2c^2}{1+\rho}n+o(n)\right\},
```
and the fixed-replica second-moment exponent is
```math
\Delta_c(q)=
h\!\left(\frac{1-q}{2}\right)-\log2+
\frac{2c^2q^2}{1+q^2}\le0
\qquad(c\le1/2).
```
Thus fixed replicas are replica-symmetric in the relevant range but
control only the speed-$`n`$ violation count, not the speed-$`n^2`$
probability of a hole containing no violating direction.

### 10.26 Grothendieck aggregation and spectral regularization

A new exact norm inequality improves the spectral bootstrap.  Put
```math
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top Ax|=2M(A).
```
Zero diagonal makes $`z^\top Az`$ multiaffine, so its absolute maximum
on the box is $`Q(A)`$.  Polarization then gives
```math
\|A\|_{\infty\to1}\le2Q(A).
```
Since $`\|A\|_{1\to\infty}=\max_{ij}|a_{ij}|=1`$, interpolation yields
```math
\boxed{
\|A\|_{\mathrm{op}}\le\sqrt{2Q(A)}=2\sqrt{M(A)}.
} \tag{10.67}
```
Competitive signings therefore have operator norm $`O(n^{3/4})`$,
improving the earlier $`O(n^{5/6})`$ bootstrap.

There is also an aggregated mixed-sign spectral inequality.  For any
orthonormal eigenpairs $`Au_r=\lambda_ru_r`$ and weights $`w_r\ge0`$,
let
```math
D=\max_i\sum_rw_ru_r(i)^2.
```
Real Grothendieck, with the signs of $`\lambda_r`$ placed in one vector
family, gives
```math
\boxed{
Q(A)\ge
\frac{\sum_rw_r|\lambda_r|}{2K_GD}.
} \tag{10.68}
```
This rules out a uniformly diffuse high-eigenvalue subspace; all
spectral signs are handled simultaneously.

Most importantly, Grothendieck--Pietsch factorization gives a spectral
removal theorem.  For every $`\varepsilon>0`$, there is
$`R\subset[n]`$, $`|R|\ge(1-\varepsilon)n`$, such that
```math
\boxed{
\|A[R]\|_{\mathrm{op}}
\le\frac{4K_GQ(A)}{\varepsilon n}.
} \tag{10.69}
```
Indeed, factor the bilinear form using probability weights
$`\mu,\nu`$, delete coordinates where either weight exceeds
$`2/(\varepsilon n)`$, and restrict the factorization.  Equivalently,
for $`Q=O(n^{3/2})`$, every spectral scale
$`L/\sqrt n\to\infty`$ is supported on a common $`o(n)`$-vertex
exceptional set.

The theorem is not yet the desired scale transfer.  After deleting a
fixed $`\varepsilon n`$, Hanson--Wright has an exponential-in-$`n`$
selector tail, but its constant is too small to pay the full
$`\alpha n\log2`$ spin entropy.  Letting $`\varepsilon\to0`$ makes the
tail coefficient vanish.  A $`\sqrt n`$-hub construction shows this
loss is real for norm-only estimates and is repaired precisely by
deleting the hub.  Nonuniform Pietsch sampling neutralizes hub
variance but pays a leading bias; ordinary chaining over all full
extensions is also impossible because extension averaging forces a
deterministic fluctuation larger than the scale-transfer margin.

The current frontier is therefore exceptionally concrete: combine
(10.69) with either dependent rounding that cancels the first-order
selector field, a max-extension quotient, or a recursive Pietsch
hierarchy whose entropy losses are summable.  Any successful version
at each fixed proportional scale would prove the limit.

### 10.27 Exact obstructions after spectral regularization

The three natural ways of using (10.69) all have leading-scale
obstructions.

First, nonuniform Pietsch sampling can suppress the selector variance
at high-weight vertices, but moving inclusion probabilities away from
the uniform value $`\alpha`$ creates a deterministic bias of order
```math
\|A\|_{\infty\to1}
\bigl(2\alpha\|p-\alpha\mathbf1\|_\infty+
\|p-\alpha\mathbf1\|_\infty^2\bigr),
```
which consumes the entire scale-transfer margin.  At half density the
first ANOVA term cannot be canceled by any dependent rounding.  If
$`e=\mathbf1_S-\tfrac12\mathbf1`$, then
```math
D_eAD_{1/2}+D_{1/2}AD_e
=\frac12\bigl(A[S]\oplus(-A[S^c])\bigr), \tag{10.70}
```
whose Boolean norm is already $`\Theta(n^{3/2})`$ for every balanced
selector.

Second, recursive Pietsch peeling sees no hub at all on a flat
conference core: the uniform factorization measure is valid, so the
recursion terminates and gives only linear-in-$`\alpha`$, not
$`\alpha^{3/2}`$, control.

Third, the max-extension quotient has an orientation obstruction.  For
```math
A=\begin{pmatrix}B&C\\C^\top&D\end{pmatrix}
```
and $`y\in\{\pm1\}^{S}`$, the positive and negative extension gaps are
```math
\begin{aligned}
D_S^+(y)&=\max_z\{H_D(z)+|(C^\top y)\cdot z|\},\\
D_S^-(y)&=\max_z\{-H_D(z)+|(C^\top y)\cdot z|\}.
\end{aligned} \tag{10.71}
```
Their sum has a strong rectangular lower bound, but either oriented
gap can be small.  Cross columns can be chosen orthogonal to both a
positive and a negative extremizer of $`B`$, even while the cross block
remains spectrally and rectangularly regular.

### 10.28 Centered width and the midpoint conjecture

Let
```math
W(A)=\frac{\max H_A-\min H_A}{2}.
```
There is an exact cut-norm identity
```math
\boxed{
W(A)=\max_{S\subset[n]}
\|A_{S,S^c}\|_{\infty\to1}.
} \tag{10.72}
```
Consequently $`G_n=\min_AW(A)`$ is superadditive.  If
```math
R_{m,n}=\min_{B\in\{\pm1\}^{m\times n}}\|B\|_{\infty\to1},
```
then
```math
G_{n+m}\le G_n+G_m+R_{n,m},
```
and $`R`$ is symmetric and separately subadditive.  These exact
two-parameter axioms still do not force $`G_n/n^{3/2}`$ to converge:
smoothly rounded versions of
```math
g(x)=x^{3/2}\bigl(1+\varepsilon\sin(\log\log(x+e^e))\bigr)
```
together with $`r(m,n)=K\sqrt{mn(m+n)}`$ satisfy all of them and
oscillate.

Nor can internal width and cross discrepancy be combined by a
uniform Pythagorean inequality.  An exact order-$`5`$, $`1+4`$ block
example has
```math
W(A)=W(D)=4,\qquad \|C\|_{\infty\to1}=4.
```
In general, with endpoint slacks $`s_+,s_-`$,
```math
W(A)=W(B)+W(D)+\frac12\left[
\max(|\phi|-s_+)+\max(|\phi|-s_-)
\right], \tag{10.73}
```
so a large cross norm may lie entirely below both internal
energy-layer slack profiles.

There is nevertheless strong evidence for a separate midpoint theorem.
Write $`P=\max H_A`$, $`Q=-\min H_A`$.  Exact minimizers through order
$`10`$ have $`|P-Q|\le2`$; all order-$`8`$ minimizers are exactly
centered.  In coding language, if $`C`$ is the cut code and
$`D=C\cup(\mathbf1+C)`$, one asks whether a deepest $`D`$-hole $`a`$
can be chosen with
```math
\boxed{
|d(a,C)-d(a,\mathbf1+C)|\le1.
} \tag{10.74}
```
The statement holds for every binary linear code through ambient
length $`7`$, extensive random codes through length $`16`$, and more
than $`7.15\times10^8`$ tested cubelike center instances.  It is a
two-center bisector statement in a cubelike quotient graph and is
related to, but weaker than, the classical open normal-code
conjecture.  The deepest-hole plateau can be totally disconnected, so
geodesic intermediate-value proofs fail.

Even if (10.74) is proved, an additional scale-transfer theorem is
needed.  The most direct centered rectangular-separation target is
false on conference orders $`6,14,18`$; at order $`18`$, half-subsets
supply at most $`34`$ against the required $`42.666`$.

### 10.29 Purification and action concentration--compactness

Weak action convergence is fine enough to retain many operator
profiles but still loses non-uniformly-integrable Boolean spikes.  An
$`n^{3/4}`$-vertex dense clique can disappear in every finite action
profile while changing $`M/n^{3/2}`$ by a fixed constant.

For minimization, however, such bubbles can be removed.  Fix
$`\varepsilon>0`$.  Apply (10.69), replace the exceptional
$`\varepsilon n`$-vertex principal block by a flat signing, and replace
its cross edges by a rectangular signing of
$`\infty\to1`$ norm $`O(\sqrt\varepsilon\,n^{3/2})`$.  The repaired
signing $`A'`$ satisfies
```math
\boxed{
M(A')\le M(A)+O(\sqrt\varepsilon)\,n^{3/2},
\qquad
\|A'\|_{\mathrm{op}}\le C_\varepsilon\sqrt n.
} \tag{10.75}
```
Thus for every $`\eta>0`$, the infimum is uniformly
$`\eta n^{3/2}`$-approximable inside a $`2\to2`$-bounded action class,
where the Boolean objective is action-continuous by uniform
integrability.

This produces compact liminf objects but not all-order realization.
Any sign-block blow-up preserving a macro coefficient must leave a
residual Frobenius budget $`(1-o(1))k^2`$ in its orthogonal fiber
modes.  Random residuals add a Wigner component; Hadamard residuals
add a tensor component; constant blocks amplify the base operator.
Consequently no known rational blow-up realizes a bounded signed
action limit at every large order.  Compactness and continuity alone
still permit slowly drifting finite admissible sets.

### 10.30 Stratified conference perturbation

Sparse random edge perturbation gives a concrete conditional
improvement of the upper constant.  Start with a conference signing and
flip each edge independently at rate $`\delta`$.  Conditional on
```math
r=\frac{x^\top Cx}{n\sqrt{n-1}},
```
the normalized half-energy upper-tail rate at level $`b`$ is
```math
\frac{(b-(1-2\delta)r/2)^2}{4\delta(1-\delta)}. \tag{10.76}
```
Suppose all nonresonant conference energy layers obey the flat
entropy envelope
```math
s_f(r)=\log2+\frac14\log(1-r^2),
\qquad 0\le r\le r_*=\frac{\sqrt{15}}4, \tag{10.77}
```
while cap-near square-wave resonances form only $`\exp(o(n))`$ states.
Then the perturbed upper constant is
```math
G(\delta)=\max\left\{
\frac12-\delta,\ 
\max_{r\le r_*}\left[
\frac{(1-2\delta)r}{2}
+2\sqrt{\delta(1-\delta)s_f(r)}
\right]\right\}. \tag{10.78}
```
Numerical optimization, independently checkable from (10.78), gives
```math
\boxed{
\delta_*=0.001394039184\ldots,\qquad
G(\delta_*)=0.498605960816\ldots<\frac12.
} \tag{10.79}
```
This would be the first rigorous upper improvement if the stratified
entropy theorem is proved.  The same variational calculation shows
that independent perturbation cannot reach $`\sqrt{15}/8`$: enough
noise to suppress cap resonances raises the bulk term above the ROM
value.  The exact open lemma is therefore a Paley/conference
delocalization theorem with the known square-wave resonant families
explicitly excised.

### 10.31 Corrected two-step, finite-type, profile, and dependent-lift frontier

This section records the verified conclusions of the subsequent
multi-agent wave. The normalization is
```math
Q(A)=\max_{x\in\{\pm1\}^n}|x^\top Ax|=2M(A).
```
Population computations, finite computations, and open lemmas are
labeled explicitly.

#### 10.31.1 Corrected dependent two-step constant

For a normalized symmetric conference involution $`U`$, the dependent
two-step population law uses
```math
F=\operatorname{sign}(G+tS),\qquad
R=\frac{F-aG-bS}{\sqrt{1-a^2-b^2}},
```
where $`S`$ is a sign, $`G,W`$ are independent standard Gaussians, and
$`a=\mathbb E(GF), b=\mathbb E(SF)`$. For a second Boolean response
$`Y=g(S,G,R,W)`$, write
```math
(\alpha,\beta,\gamma,\delta)
=\bigl(\mathbb E(SY),\mathbb E(GY),
\mathbb E(RY),\mathbb E(WY)\bigr).
```
The paired population energy is exactly
```math
\mathcal E(g)=2(\alpha\beta+\gamma\delta). \tag{10.80}
```
Its nondegenerate stationary responses obey
```math
Y=\operatorname{sign}
(\beta S+\alpha G+\delta R+\gamma W).
```

**Computed population optimum, corrected.** Splitting the Gaussian
integral at $`G=-tS`$, then optimizing the exact finite-dimensional
functional, gives
```math
\boxed{
t_2=0.8414699114\ldots,\qquad
c_2=0.783387533648\ldots
}
\tag{10.81}
```
in doubled $`Q/n^{3/2}`$ normalization. The corresponding original
half-energy value would be $`c_2/2=0.391693766824\ldots`$.

**Exact-conference theorem.** Fixed dithers make the response smooth;
the Onsager residual has hybrid Hermite rank at least three. Truncation,
flat-involution contraction estimates, homogeneous-sum invariance, and
the multivariate fourth-moment theorem prove
```math
\boxed{
\liminf_{n\to\infty}
\frac{Q(C_n)}{n\sqrt{n-1}}\ge c_2
}
```
for every sequence of exact symmetric conference matrices $`C_n`$.
The infinite Hermite tail is removed before the dithers are sent to
zero. Thus (10.81) is a rigorous conference lower theorem; its decimal
is the numerical evaluation of the explicit population variational
integral. Extension from conference matrices to every competing
signing remains open. See `dependent_two_step_amp.md`.

#### 10.31.2 Corrected Schatten--3 inequality

Two exact inequalities replace the earlier weaker bootstrap:
```math
\boxed{
\|A\|_{\mathrm{op}}^2
\le\|A\|_{\infty\to1}
\le2Q(A)
}
\tag{10.82}
```
and
```math
\boxed{
Q(A)\ge
\frac{\operatorname{tr}|A|^3}{2K_G(n-1)}.
}
\tag{10.83}
```
For (10.83), the left Grothendieck vectors have spectral coordinates
$`|\lambda_r|u_r(i)`$, while the right vectors have coordinates
$`\operatorname{sgn}(\lambda_r)|\lambda_r|u_r(i)`$. Their norms are
one after division by $`\sqrt{n-1}`$, because
$`(A^2)_{ii}=n-1`$. The factor $`2K_G`$, including the factor $`2`$
needed to return from bilinear to same-sign Boolean norm, is the
corrected constant.

For
```math
B=A/\sqrt{n-1},\quad
C_A=Q(A)/(n\sqrt{n-1}),\quad
\mu_4=n^{-1}\operatorname{tr}B^4,
```
(10.82)--(10.83) imply
```math
\mu_4
\le2\sqrt2K_G C_A^{3/2}n^{1/4}(1+o(1)). \tag{10.84}
```
At $`C_A=c_2`$, this is
```math
\mu_4\le(3.49519+o(1))n^{1/4}. \tag{10.85}
```
It excludes the dense profile $`q_{ij}\asymp n^{-1/4}`$ on
$`\Theta(n^2)`$ row pairs, but not the sparse profile
$`|q_{ij}|\asymp n^{-1/8}`$ on $`\Theta(n^{3/2})`$ pairs. See
`multispike_grothendieck_spectral_removal.md`.

#### 10.31.3 Exact three-fibre obstruction

Let $`C_k`$ be a symmetric conference matrix and put
```math
R=J_3-2I_3,\qquad D=J_3-I_3,\qquad
A_k=C_k\otimes R+I_k\otimes D. \tag{10.86}
```
This is an order-$`3k`$ signing with $`Q(A_k)=O((3k)^{3/2})`$, and
```math
A_k^2
=I_k\otimes\bigl((k-1)R^2+D^2\bigr)
+4C_k\otimes I_3. \tag{10.87}
```
For $`B_k=A_k/\sqrt{3k-1}`$, distinct rows in the same fibre have
correlation tending to $`-1/3`$; correlations between different
fibres are zero or $`O(k^{-1})`$. Hence
```math
\boxed{
\Lambda_4(B_k)
=\frac1{3k}\sum_{i\ne j}(B_k^2)_{ij}^4
\longrightarrow\frac2{81}.
}
\tag{10.88}
```
This persists after every $`o(n)`$-vertex deletion. Thus bounded
$`Q/n^{3/2}`$, Schatten control, and small vertex deletion do not
force the conference cavity law.

There is an exact four-state formula. With antipodal fibre
representatives $`u_0=(1,1,1)`$, $`u_1=(-1,1,1)`$,
$`u_2=(1,-1,1)`$, $`u_3=(1,1,-1)`$, the kernels are
```math
K=
\begin{pmatrix}
3&1&1&1\\
1&-5&3&3\\
1&3&-5&3\\
1&3&3&-5
\end{pmatrix},
\qquad
\ell=(6,-2,-2,-2). \tag{10.89}
```
Thus
```math
Q(A_k)=
\max_{z_i,t_i}
\left|
\sum_{i,j}(C_k)_{ij}z_i z_jK_{t_i,t_j}
+\sum_i\ell_{t_i}
\right|. \tag{10.90}
```
For the order-six Paley conference seed, an exact finite certificate
gives
```math
\boxed{Q(A_6)=78.} \tag{10.91}
```
At order $`18`$, this is $`1.0213764617\ldots`$ in doubled
normalization, or $`0.5106882308\ldots`$ in the original
normalization. This lift is not an upper-bound amplification
mechanism. See `finite_fibre_renormalization.md`.

#### 10.31.4 Finite-type cavity law and mesoscopic gap

For a fixed sign fibre $`T=R/\sqrt s`$, the exact finite-channel
population functional is
```math
\mathcal E_T(Y)=\frac2s\left[
\operatorname{tr}(T\mathcal A T\mathcal H^\top)
+\operatorname{tr}(T\mathcal M T\mathcal J^\top)
\right]. \tag{10.92}
```
The first nonlinear residual has an exact variance-rigidity correction.
If $`C=T^2`$, then at the optimized scalar threshold
```math
\frac1s\operatorname{tr}(TK_0T)-s_t^2
\ge
0.00445420230\ldots\,
\frac1s\sum_{i\ne j}C_{ij}^4. \tag{10.93}
```
Equality requires $`T^2=I`$. Whitening the two receiving fields
produces orthogonal coupling matrices in both paired channels, so the
remaining finite-dimensional obstruction is cube-versus-rotation,
not loss of Gaussian information.

**Open finite-type lemma.** It remains to prove
```math
\sup_Y\mathcal E_T(Y)\ge c_2
\quad\text{for every fixed symmetric sign fibre }T. \tag{10.94}
```
For the three-fibre obstruction, the inherited scalar response gives
approximately $`0.78344`$, about $`5\cdot10^{-5}`$ above $`c_2`$;
this is a Sobol computation, not an interval certificate. All
switching/permutation classes through $`s=4`$, and sampled fibres
through $`s=7`$, were computationally above $`c_2`$.

The mesoscopic reduction needed after (10.94) is now proved. Fix
$`0<\varepsilon<1`$. Grothendieck--Pietsch removal gives a principal
set $`R`$, $`|R|=r\ge(1-\varepsilon)n`$, on which, for
```math
B_R=A_R/\sqrt{r-1},\qquad C_R=B_R^2,\qquad
q_{ij}=(C_R)_{ij},
```
```math
\|C_R\|_{\mathrm{op}}\le L_\varepsilon+o(1),\qquad
L_\varepsilon
=\frac{(4K_Gc_2)^2}{\varepsilon^2(1-\varepsilon)}
<\frac{31.1884}{\varepsilon^2(1-\varepsilon)}.
```
For any $`\eta>0`$, retain only edges satisfying
```math
|q_{ij}|\ge \sqrt{\eta/L_\varepsilon}.
```
Because $`C_R\succeq0`$, $`(C_R)_{ii}=1`$, and
$`(C_R^2)_{ii}\le\|C_R\|_{\mathrm{op}}`$, the retained graph has
```math
\boxed{
\Delta_{\varepsilon,\eta}
\le \frac{L_\varepsilon^2}{\eta}
<\frac{972.72}
{\varepsilon^4(1-\varepsilon)^2\eta},
}
```
while the discarded correlations obey
```math
\boxed{
\frac1r\sum_{\substack{i\ne j\\
|q_{ij}|<\sqrt{\eta/L_\varepsilon}}}q_{ij}^4
\le\eta+o(1).
}
```
Thus the fixed-$`\varepsilon,\eta`$ bounded-degree reduction is a
theorem, not an open target. There is also a one-subset
$`n-o(n)`$ version: one explicit choice deletes
$`O(n^{99/100})`$ vertices, leaves maximum threshold degree
$`O(n^{31/50})`$, and discards $`o(n)`$ unnormalized fourth mass.
The remaining global step is a bounded-degree
compatible-transport theorem, since the row-correlation kernel
$`C_R`$ alone does not determine its flat square root $`B_R`$. See
`mesoscopic_correlation_dichotomy.md`.

There is one newly audited local fact at the conference point. For
fixed dither $`\psi_\tau(u)=2\Phi(u/\tau)-1`$, take the stationary
scalar response
```math
y=\psi_\tau(U),\qquad
U=hS+aG+jR+mZ,
```
with the coefficients reversed from the response moments, and put
```math
\ell_\tau=\mathbb E\psi_\tau'(U)>0,\qquad
k_\tau=-\mathbb E[U\psi_\tau''(U)]>0.
```
Perturb only by neighboring initial spins,
```math
y_i^{(\lambda)}
=\psi_\tau\!\left(U_i+
\lambda h\sum_{j\ne i}q_{ij}S_j\right).
```
A complete second-order covariance calculation at $`q=0`$ gives
```math
\boxed{
\mathcal E(Q)
=c_{2,\tau}
+h^2(2\lambda\ell_\tau-\lambda^2k_\tau)
\frac1r\sum_{i\ne j}q_{ij}^2+O_3.
}
```
With
$`\lambda_\tau=\min\{1,\ell_\tau/k_\tau\}`$, its quadratic
coefficient is at least
```math
h^2\lambda_\tau\ell_\tau>0.
```
This Hessian is exact. The earlier proposed coefficient involving the
full neighbor field $`\beta S+\delta R`$, and the uniform
connected-cluster remainder inferred from it, are withdrawn: the
dense compatible transport $`B_RKB_R`$ can make the cubic remainder
nonlocal. Positivity is proved for each fixed finite instance, but no
radius uniform in $`r`$ is yet known. See
`mesoscopic_correlation_dichotomy.md` and
`finite_type_conference_fibres.md`.

#### 10.31.5 Paley stable-tail reductions and exact walls

For a conference involution $`U`$, let $`P=(I+U)/2`$. A Gaussian
Brascamp--Lieb argument using only that $`P`$ is a rank-$`n/2`$
projection with diagonal $`1/2`$ proves, for every $`t\ge0`$,
```math
\boxed{
\mathbb E_y e^{t y^\top Uy}
\le
e^{-tn}\left(\frac{1+e^{4t}}2\right)^{n/2}.
}
\tag{10.95}
```
Consequently,
```math
\begin{aligned}
\frac1n\log\#\{y:y^\top Uy\ge(1-\varepsilon)n\}
\le{}&
\frac{\log2}{2}
+\frac{\varepsilon}{4}\log\frac{2-\varepsilon}{\varepsilon}\\
&-\frac12\log(1-\varepsilon/2).
\end{aligned}
\tag{10.96}
```
The right side tends to $`(\log2)/2`$, not zero. Constant leverage
and geometric Brascamp--Lieb therefore have a $`2^{n/2}`$ wall.

After switching by one exact positive Boolean eigenvector, every
other such eigenvector is equivalent to an equitable switch set
$`S`$ in the conference strongly regular graph:
```math
\deg_S(v)=|S|/2\quad(v\notin S),\qquad
\deg_{G[S]}(v)
=\frac{|S|-1-\sqrt{n-1}}2\quad(v\in S). \tag{10.97}
```
Spectral mixing is equality on these sets. Exact cap counting is
therefore an equitable-partition enumeration problem, not a generic
projection problem. See `projection_sign_brascamp_lieb.md`.

For cyclic Paley matrices, two further reductions are exact:

1. leave-one-out prediction gives
   $`\displaystyle \Pr(\operatorname{sign}(UY)_i\ne Y_i)\le4\varepsilon,`$
   but small $`H(Y_i\mid Y_{-i})`$ does not bound joint entropy; the
   even-parity distribution is the decisive counterexample;
2. writing
   $`\displaystyle \alpha_y=\sum_jy_j\zeta^j\in\mathbb Z[\zeta],`$
   the negative-eigenspace condition is
   $`\displaystyle \sum_{k\in N}|\sigma_k(\alpha_y)|^2 \le\varepsilon p^2. \tag{10.98}`$
   The missing theorem is a Littlewood cyclotomic near-conjugate
   lattice-point bound.

Poisson summation isolates a favorable determinant term but produces
a dual theta series containing the same resonant near-eigenvectors.
Large-spectrum lemmas control $`O_\tau(1)`$ spiky frequencies but not
flat mass on $`\Theta(p)`$ frequencies. These are exact reductions
and failures, not a stable-tail proof. See
`paley_stable_tail_cyclotomic_entropy.md`.

The hierarchical resonance library also required a correction. For
$`q=\ell^{2r+1}`$, the recursive finite-cyclic family has exact
multiplier energy $`1-\ell^{-2r-2}`$ and cardinality
```math
T_{\ell,r}
=\binom{\ell}{(\ell+1)/2}^{
(\ell^{2r+2}-1)/(\ell^2-1)}
=\exp(\Theta(q)). \tag{10.99}
```
The proved transfer to an ambient Paley prime has error at most
```math
4\sqrt{q/p}. \tag{10.100}
```
Thus Paley defect $`\varepsilon`$ forces
$`q/p=O(\varepsilon^2)`$, and the transferred ambient entropy rate is
only $`O(\varepsilon^2)`$. The corrected family disproves
subexponential counting at fixed period scale, but does not disprove
vanishing stable-tail entropy. See
`paley_hierarchical_gadget_library.md`.

#### 10.31.6 Compact spectral and action profiles

On the purified class $`\|A\|_{\mathrm{op}}=O(\sqrt n)`$, define
```math
\mu_{A,x}(E)=n^{-1}\langle x,E_A(E)x\rangle,\qquad
\mathcal K_1(A)=\{\mu_{A,x}:x\in\{\pm1\}^n\}.
```
This compact set retains isolated maximizing vectors, and
```math
\frac{Q(A)}{n^{3/2}}
=\sup_{\mu\in\mathcal K_1(A)}
\left|\int t\,d\mu(t)\right|. \tag{10.101}
```
Hausdorff convergence in $`W_1`$ makes the objective continuous.
Matrix-valued spectral measures of all finite Boolean tuples give a
compact full profile retaining overlaps and bilinear actions.

For finite-fibre composition, a smaller exact quotient suffices. For
$`X\in\{\pm1\}^{n\times r}`$, set
```math
O=n^{-1}X^\top X,\qquad H=n^{-3/2}X^\top AX,
```
and let $`\mathcal P_r(A)`$ contain all such pairs. If
$`\widehat A=A\otimes R+I_n\otimes D`$, with fibre size $`s`$ and
$`\eta=n^{-1/2}`$, then
```math
\widetilde O_{pq}
=s^{-1}\sum_aO_{(p,a),(q,a)}, \tag{10.102}
```
```math
\widetilde H_{pq}
=s^{-3/2}\left[
\sum_{a,b}R_{ab}H_{(p,a),(q,b)}
+\eta\sum_{a,b}D_{ab}O_{(p,a),(q,b)}
\right]. \tag{10.103}
```
This map is exact, continuous, and associative on the all-level
profile $`(\mathcal P_r)_{r\ge1}`$. It also displays the obstruction:
the $`d`$-th iterate reads tuple level $`s^dr`$. The maps are not
equicontinuous in the product metric, so compactness, associativity,
or an idempotent does not force convergence.

All-order recovery remains missing. A block preserving one macro
coefficient has compulsory orthogonal-fibre variance
```math
\left\|L-\frac{a}{\sqrt s}J_s\right\|_F^2=s^2-s. \tag{10.104}
```
Standard graphon blow-ups therefore inject a leading microscopic
action component. A successful extremal profile must absorb that
component and be realizable at every sufficiently large order. See
`boolean_spectral_profile_compactification.md` and
`profile_renormalization_semigroup.md`.

#### 10.31.7 Exact dependent $`4`$-lift and the centering theorem

For square $`s`$, let $`\mathcal L_s(A)`$ contain the order-$`ns`$
signings partitioned into $`n`$ fibres whose cross-block sums satisfy
```math
\mathbf1^\top B_{ij}\mathbf1=a_{ij}s^{3/2},
```
and set $`G_s(A)=\min_{B\in\mathcal L_s(A)}Q(B)`$.

**Certified finite computation.** For the all-negative triangle
$`C^-`$, $`Q(C^-)=6`$, and an explicit order-$`12`$ signing has all
three cross-block sums $`-8`$ and
```math
\boxed{
Q(B)=40<48=4^{3/2}Q(C^-).
}
\tag{10.105}
```
All $`2^{11}`$ antipodal Boolean states were enumerated independently.
The normalized values are
```math
\frac{40}{12^{3/2}}
=0.962250448649\ldots,\qquad
\frac{40}{2\cdot12^{3/2}}
=0.481125224324\ldots . \tag{10.106}
```
This is a certified witness, not a proof that $`G_4(C^-)=40`$; the
current exact lower bound is $`32`$. An analogous $`Q=40`$ witness
exists for the other triangle switching class.

The interpretation is fixed by an exact range theorem. Let
```math
P(A)=\max_xx^\top Ax,\qquad m(A)=\min_xx^\top Ax.
```
Constant-fibre configurations in any $`B\in\mathcal L_s(A)`$ have
energies $`s^{3/2}x^\top Ax+d`$, where the diagonal-fibre contribution
$`d`$ is independent of $`x`$. Hence
```math
\boxed{
G_s(A)\ge
\frac{s^{3/2}}2\bigl(P(A)-m(A)\bigr).
}
\tag{10.107}
```
For $`C^-`$, $`(P,m)=(2,-6)`$, so (10.107) gives $`32`$. The
improvement $`48\to40`$ is finite-size recentering of an asymmetric
energy interval, not repeatable normalized contraction.

The order-$`12`$ witness is chiral: an exact signed permutation $`S`$
satisfies
```math
S^2=-I,\qquad SB=-BS. \tag{10.108}
```
Thus $`P(B)=40`$, $`m(B)=-40`$, and every exact compressed four-lift
obeys
```math
\boxed{G_4(B)\ge8Q(B)=320.} \tag{10.109}
```
Strict contraction cannot repeat; equality is the only useful target.

The constrained classes compose:
```math
B\in\mathcal L_s(A),\ C\in\mathcal L_t(B)
\Longrightarrow C\in\mathcal L_{st}(A),
```
so
```math
G_{st}(A)\le
\min_{B\in\mathcal L_s(A)}G_t(B). \tag{10.110}
```
This is profile-valued, not a closed scalar inequality. See
`dependent_profile_recovery.md`,
`dependent_profile_recovery_witness.json`, and
`verify_dependent_profile_recovery.py`.

#### 10.31.8 Chiral equality conditions and complete Clifford-mask no-go

After signed permutation of coordinates, (10.108) gives
```math
B=Z\otimes A_0+X\otimes C_0
=
\begin{pmatrix}A_0&C_0\\C_0&-A_0\end{pmatrix}, \tag{10.111}
```
or $`x^\top Bx=\operatorname{Re}(z^\top(A_0-iC_0)z)`$ on
$`z\in\{\pm1\pm i\}^6`$. The original $`4\times4`$ fibre blocks do
not generate a small algebra: two already generate
$`M_4(\mathbb R)`$. The witness has twelve distinct eigenvalues, so
it is not switching-equivalent to a conference matrix, a regular
two-graph, or a strongly regular Seidel construction.

The natural inherited-chiral equality ansatz is
```math
\widehat B
=B\otimes R+(SB)\otimes K+S\otimes E
+\operatorname{diag}_uD_u. \tag{10.112}
```
Here $`R,K`$ are complementary-support symmetric partial sign
matrices, $`R`$ contains the micro-diagonal,
```math
\sum R=8,\qquad\sum K=0,
```
$`E`$ is skew on the $`K`$-support, and a signed micro-involution
$`P`$ satisfies
```math
[P,R]=[P,K]=0,\qquad\{P,E\}=0,\qquad
D_{\pi(u)}=-PD_uP, \tag{10.113}
```
where $`\pi`$ is the underlying coordinate permutation of the signed
macro operator $`S`$.
These conditions make $`\widehat B`$ a valid exact compressed
four-lift with inherited chiral symmetry.

**Exact finite no-go.** There are $`1{,}008`$ compatible
$`(P,R,K,E)`$ masks. Repeated-fibre configurations and the exact joint
Boolean profile of $`(B,SB)`$ give
```math
\begin{array}{c|rrrr}
\text{lower bound}&320&400&416&448\\ \hline
\text{number of masks}&96&288&336&288.
\end{array}
\tag{10.114}
```
The $`96`$ survivors form four simultaneous-coordinate-permutation
orbits. For those four orbits, one explicit two-pattern fibre
configuration has minimum absolute energy
```math
392,\quad480,\quad392,\quad480 \tag{10.115}
```
over all $`2^{24}`$ compatible diagonal fillings. Exact integer
subset-sum dynamic programming certifies (10.115); no heuristic
optimization is used.

Thus all $`1{,}008`$ inherited-chiral Clifford masks are rigorously
ruled out at the necessary equality value $`320`$. This is a no-go
for (10.112), not for every nonlocal dependent microblock lift. A
scale-preserving continuation, if it exists, must leave this
complementary-support Clifford algebra. See
`chiral_scale_preserving_lift.md` and
`verify_chiral_clifford_no_go.py`; the verifier also certifies over
$`\mathbb Q`$ both the full $`M_4`$ fibre algebra and the square-free
characteristic polynomial used above.

#### 10.31.9 Static microblock filters: exact ANOVA no-go

There is a separate exact obstruction to repairing amplification with
a static $`A`$-dependent polynomial quotient. For an $`r\times r`$
cross block, let $`e_0=\mathbf1/\sqrt r`$ and decompose
```math
B_{ij}
=s_{ij}e_0e_0^\top+q_{ij}e_0^\top
+e_0c_{ij}^\top+R_{ij},
\qquad C_{ij}=s_{ij}/\sqrt r.
```
The four summands are Frobenius-orthogonal, so the sign constraint
gives the exact conservation law
```math
\boxed{
r^2
=rC_{ij}^2+\|q_{ij}\|_2^2+\|c_{ij}\|_2^2
+\|R_{ij}\|_F^2.
}
\tag{10.116}
```

Condition each fibre to magnetization $`m_i=t\xi_i`$, with the
$`\xi_i`$ independent signs, and set
```math
\alpha=\frac r{r-1}(1-t^2).
```
The exact Hamming-slice Hoeffding/ANOVA decomposition, averaged over
$`\xi`$, implies for $`t^2\ge1/r`$
```math
\boxed{
\mathbb E_\xi\operatorname{Var}(\mathcal H_B\mid t\xi)
\ge
\alpha^2\sum_{i<j}\bigl(r^2-rC_{ij}^2\bigr).
}
\tag{10.117}
```
Regular blocks attain the associated continuous allocation minimum.
Consequently every bounded natural-scale quotient leaves a
$`1-O(1/r)`$ fraction of its squared mass in microscopic modes.

In particular, an odd static filter
```math
C_p=\sqrt n\,p(A/\sqrt n),
\qquad
C_\eta=A-\frac{\eta}{n}A^3
```
can shrink the coarse spectral extremes, but it does not shrink the
leading residual variance. For regular blocks the global identity
```math
N(N-1)=r\|C_p\|_F^2+\|B_\perp\|_F^2
\tag{10.118}
```
also forces every purely continuous operator-norm certificate back to
the original half-energy floor $`1/2`$. This rules out static
polynomial filtering as a continuous or one-step second-moment
amplification mechanism. It does not rule out dependent Boolean
microblocks such as (10.105), which exploit the full action profile.
See `spectral_filter_microblock_recovery.md`.

### 10.32 Centered-width recovery, cut-cone tangents, and compatible graphings

This section is the post-v40 ledger.  It uses the one-copy
normalization
```math
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
W(A)=\frac{\max H_A-\min H_A}{2}.
```
Every statement below is labelled by its actual status.  In
particular, no convergence theorem is claimed.

#### 10.32.1 Centered-width monotonicity and the first equality no-go

For square fibre size $`s`$, let $`\mathcal L_s(A)`$ be the exact
compressed-lift class
```math
\mathbf1^\top B_{ij}\mathbf1=a_{ij}s^{3/2}.
```
Constant-fibre states retain the complete translated seed interval,
so translations cancel from its range and
```math
\boxed{
W(B)\ge s^{3/2}W(A)
\qquad(B\in\mathcal L_s(A)).
}
\tag{10.119}
```
Thus exact compressed lifts can only preserve or increase normalized
centered width.  An amplification proof in this class needs
asymptotic equality; finite strict contraction is impossible.

Equality already fails for the first nontrivial test.  For the
all-negative triangle $`C^-`$,
```math
[m(C^-),P(C^-)]=[-3,1],\qquad W(C^-)=2.
```
At $`s=4`$, equality in (10.119) would give $`W(B)=16`$.  Endpoint
local-field constraints reduce every equality candidate to $`24,000`$
integer block patterns.  Exact enumeration of all antipodal
microstates gives
```math
\boxed{
\min_{\text{equality candidates}}W(B)=20>16.
}
\tag{10.120}
```
The earlier order-$`12`$, doubled-$`Q=40`$ witness has interval
$`[-20,20]`$; its apparent absolute-norm gain is finite recentering,
not centered-width contraction.  See
`centered_width_amplification_reboot.md` and
`verify_centered_width_equality_no_go.py`.

#### 10.32.2 Exact traffic criterion and compulsory variance

For a microstate $`X=(x_i)_{i\le n}`$, let
$`\mu_i=s^{-1}\mathbf1^\top x_i`$, and separate an exact-sum
microblock into its coarse part and residual $`R_B(X)`$.  The upper
and lower residual excesses $`E_+,E_-`$ satisfy the exact identity
```math
\boxed{
W(B)-s^{3/2}W(A)=\frac{E_+(A,B)+E_-(A,B)}2.
}
\tag{10.121}
```

If an $`s\times s`$ sign block is sampled uniformly with prescribed
sum $`a_{ij}s^{3/2}`$, then for fibre spins of means $`\mu,\nu`$,
```math
\boxed{
\operatorname{Var}(x^\top B_{ij}y)
=\kappa_s(1-\mu^2\nu^2),
\qquad
\kappa_s=\frac{s^3}{s+1}.
}
\tag{10.122}
```
For a full microstate, with $`z_i=\mu_i^2`$,
```math
V_s(X)
=
\frac{\kappa_s}{2}
\left[
n(n-1)-\left(\sum_i z_i\right)^2+\sum_i z_i^2
\right].
\tag{10.123}
```
Bernstein plus a union bound therefore gives an exact recovery
criterion $`\mathcal T_{A,s}(e)<1`$.  At normalized slack $`r`$ and
squared-magnetization statistic $`\rho`$, its asymptotic shell
threshold is
```math
\boxed{
\text{entropy rate}
<
\frac{r^2}{\kappa_s(1-\rho^2)}.
}
\tag{10.124}
```
Among nontrivial square $`s\ge4`$,
```math
\max_s\kappa_s^{-1}=\kappa_4^{-1}=\frac5{64}.
\tag{10.125}
```

The variance in (10.122) is compulsory, not an artifact of the
chosen ensemble.  The exact block ANOVA identity
```math
r^2=rC_{ij}^2+\|q_{ij}\|_2^2+\|c_{ij}\|_2^2+\|R_{ij}\|_F^2
\tag{10.126}
```
and its Hamming-slice variance bound show that static polynomial
filters and independent block gadgets retain leading microscopic
mass.  Independent internal refill is even more decisively excluded:
at every competitive $`W\le(1/2+o(1))n^{3/2}`$, its cap union sum is
at least
```math
\exp\{(\log2-1/4-o(1))n\}.
\tag{10.127}
```
Correlated, cap-conditioned replacement remains open.  See
`asymptotic_centered_width_recovery.md`,
`spectral_filter_microblock_recovery.md`, and
`internal_resampling_traffic_no_go.md`.

#### 10.32.3 Global replacement and its minimax dual

For a pair $`p=(x,y)`$, put
```math
s_e(p)=\frac{x_ix_j-y_iy_j}{2},\qquad
\delta_A(p)=W(A)-\sum_ea_es_e(p).
```
Replacing an edge block $`T`$ by $`\beta\in\{\pm1\}^T`$ gives
```math
\boxed{
W(A^{T\to\beta})-W(A)
=
\max_p\{(\beta-a_T)\cdot s_T(p)-\delta_A(p)\}.
}
\tag{10.128}
```
This is the exact global replacement theorem.

For a global width minimizer and
```math
\eta_T=\sqrt{2|T|(2v(T)\log2+1)},
```
finite minimax plus sign rounding produces a pair-profile law
$`\mu_T`$ such that
```math
\mathbb E_{\mu_T}\operatorname{score}_A
\ge W(A)-\eta_T,
\qquad
2\sum_{e\in T}(a_e\mathbb E_{\mu_T}s_e)_+
\le\eta_T.
\tag{10.129}
```
Thus every block has a near-active dual measure whose mean gradient
is almost coordinatewise anti-aligned with the current signing.

For an induced vertex block $`S`$, with cross matrix $`C`$,
global minimality also gives the deterministic exchange inequality
```math
\boxed{
W(A[S])\le W_{|S|}
+2\|C\|_{\infty\to1}.
}
\tag{10.130}
```
It rules out a localized high-width block behind a weak boundary.
It does not rule out a diffuse $`n^{-1/2}`$-density midpoint spread
over $`\Theta(n^2)`$ edges.  See
`global_block_replacement_midpoint.md`.

#### 10.32.4 Local cut geometry: two exact counterexamples

Single-edge optimality and cut-triangle closure do not center the
energy interval.  The all-positive signing is an edgewise flat local
minimum of $`W`$, yet has midpoint $`\Theta(n^2)`$.  More sharply,
there are signings consisting of an all-positive block of order
$`k\asymp n^{3/4}`$, a chiral centered bulk, and a weak random
boundary, for which
```math
W(A_n)=O(n^{3/2}),\qquad
|d(A_n)|=\left(\frac14+o(1)\right)n^{3/2}.
\tag{10.131}
```
They satisfy all cut, symmetric-difference, triangle, and
maximum-cut domination identities.

The construction (10.131) is not a global minimizer: its planted
block violates (10.130).  Hence the audit cleanly separates the
claims:

* local endpoint coverage plus triangle geometry is insufficient;
* correct $`n^{3/2}`$ scale plus full cut domination is insufficient;
* global replacement excludes localized planted spikes;
* diffuse midpoint bias remains open.

See `cut_triangle_midpoint_balancing.md` and
`global_block_replacement_midpoint.md`.

#### 10.32.5 Dual capped cuts and the tangent-covariance obstruction

Switch a top state to $`\mathbf1`$, orient a bottom state as a top
state of the negated signing, and set
```math
c=\frac{A+\widetilde A}{2},\qquad
h=\frac{A-\widetilde A}{2}.
```
Then $`c`$ is the cross block of the relative endpoint partition,
$`h=B\oplus D`$ consists of its two internal blocks,
```math
c\cdot\mathbf1=W,\qquad h\cdot\mathbf1=d,
```
and every cut obeys
```math
\boxed{
|h\cdot\delta(S)|
\le
\min\{c\cdot\delta(S),W-c\cdot\delta(S)\}.
}
\tag{10.132}
```
Equivalently, if $`C`$ is the rectangular cross signing,
```math
\boxed{
|d-H_B(p)-H_D(q)|
\le W-|p^\top Cq|.
}
\tag{10.133}
```
Thus every near-cap orientation of $`C`$ lies in an equally thin
internal-energy tube.  Exact cap cuts annihilate $`h`$; in the capped
dual-cut polytope, $`c\pm h`$ lie in the minimal face of $`c`$.
Extreme $`c`$ would force $`h=0`$, but exact cap faces of actual small
minimizers are highly non-extreme.

There are two rigorous covariance consequences.  If
$`\mathcal L_r=\{c\cdot\delta\le rn^{3/2}\}`$,
$`\mu`$ is supported there, and $`R=\mathbb E_\mu zz^\top`$, then
for $`d\ge\varepsilon n^{3/2}`$,
```math
\operatorname{tr}R^2-n
\ge(4(\varepsilon-2r)_+^2-o(1))n.
\tag{10.134}
```
For any sign-symmetric cube law,
```math
\boxed{
\mathsf H(Z)
\le n\log2-
\frac{\operatorname{tr}R^2-n}{8\|R\|_{\rm op}}.
}
\tag{10.135}
```
The resulting exponent is far too weak for (10.124).

The correct covariance lives on internal degree-two features
$`\phi(z)=(z_iz_j)_{ij\in E_h}`$.  With
$`\Sigma_r=\operatorname{Cov}_\mu\phi`$,
```math
\boxed{
h^\top\Sigma_rh\le4r^2n^3.
}
\tag{10.136}
```
At $`r=0`$, $`h`$ is an exact flat kernel vector.  More generally,
```math
\|P_{>\lambda}h\|_2^2\le\frac{4r^2n^3}{\lambda},
\tag{10.137}
```
and the same bound holds for every conditional block Schur
complement.  This is a verified feature-span degeneracy, not yet a
replacement theorem.

Bare covariance cannot force a traffic entropy rate: a symmetrized
product cloud with coordinate bias $`n^{-1/4}`$ has precisely the
required $`\Theta(n)`$ degree-two barycenter mass at only
$`O(\sqrt n)`$ entropy cost.  A shifted Hubbard--Stratonovich
one-replica LDP is also quantitatively capped far below what traffic
needs.  See `dual_cut_cone_diffuse_midpoint.md`,
`cut_cone_tangent_entropy.md`, and `endpoint_slack_ldp.md`.

#### 10.32.6 Exact-cap rank audit and the thick-cap frontier

Exhaustive switching-class enumeration through order eight computes
the affine rank of internal cut features at
$`|p^\top Cq|=W`$.  The rank is usually zero or one; at $`n=8`$, all
$`67,200`$ endpoint pairs have only the two trivial cap
representatives and rank zero.  Exact-face extremality is therefore
falsified on the known minimizers.

Thickening the cap repairs rank but not scale.  If $`r_i`$ are the
cross local fields and $`L_T=\{i:r_i\le T\}`$, every singleton and
double flip in $`L_T`$ lies in the $`4T`$-near-cap layer.  Their
internal edge features span the full $`K_{|L_T|}`$ edge space, and
the associated feature matrix has
```math
\boxed{\sigma_{\min}=2\quad(|L_T|\ge5).}
\tag{10.138}
```
However, typical fields are $`T=\Theta(\sqrt n)`$, so inversion loses
exactly the factor $`\sqrt n`$ needed to control the flat
coefficients.  Thresholding gives
```math
|\{i:r_i>T\}|\le W/T.
\tag{10.139}
```
Making the exceptional set $`o(n)`$ forces
$`T\gg\sqrt n`$, which worsens the cap margin; taking
$`T=o(\sqrt n)`$ leaves a linear exceptional set.  This critical
tradeoff is the current thick-cap frontier.  See
`thick_near_cap_rank.md`.

#### 10.32.7 Bounded-degree compatible graphing: one local theorem

After deleting $`\varepsilon n`$ vertices by
Grothendieck--Pietsch and putting
```math
B=A_R/\sqrt{r-1},\qquad C=B^2=I+Q,
```
one has a dimension-free operator bound on $`C`$.  Thresholding the
row correlations $`q_{ij}`$ leaves a bounded-degree large-edge graph;
the discarded correlations have arbitrarily small normalized fourth
mass in the order
```math
r\to\infty,\qquad\eta\downarrow0,\qquad\varepsilon\downarrow0.
```

After smoothing both sign maps, the neighbor-spin perturbation has
an audited positive Hessian.  For
$`\widehat Q=(|q_{ij}|)`$,
$`\varrho=\|\widehat Q\|_{\rm op}\le\Delta q_{\max}`$, there is a
dimension-free constant such that
```math
\mathcal E(Q)
\ge c_{2,\sigma,\tau}
+\frac{a_\tau^{\rm safe}}2\frac{\operatorname{tr}Q^2}{r}
\tag{10.140}
```
whenever $`\Delta q_{\max}`$ is below an explicit fixed threshold.
This uniform small-edge theorem is proved.

The proposed global compatible homotopy
```math
C_u=I+uQ,\qquad
B_u=\operatorname{sgn}(B)C_u^{1/2}
\tag{10.141}
```
is not proved monotone.  Its exact stationary envelope derivative
contains an indefinite residual-regression term plus a skew polar
velocity,
```math
\mathscr I_u
=
\operatorname{tr}(T_u\widetilde M_u'O_{2,u}J_0^\top)
+
\operatorname{tr}(T_u\widetilde M_u
\Omega_uO_{2,u}J_0^\top),
\tag{10.142}
```
with no known positive-square representation.  A former size-four
numerical counterexample was caused by normalization and iteration
error and is withdrawn; the corrected census shows no coarse-grid
decrease but is numerical only.  The earlier neighbor-$`(S,R)`$
coefficient is also withdrawn because it omitted residual-channel
covariance terms.  See `bounded_degree_compatible_graphing.md`.

#### 10.32.8 Sparse flat squares: amenable and nonamenable models

For
```math
S=A^2-(n-1)I,
```
off-diagonal entries satisfy
```math
S_{ij}\equiv n-2\pmod2.
```
Hence exact sparse squares require even order.  Also $`A`$ and $`S`$
commute, and on every $`S`$-eigenspace,
```math
A=\sum_\theta\sqrt{n-1+\theta}\,J_\theta,
\tag{10.143}
```
with a self-adjoint involution $`J_\theta`$.  Nonsquare integral
eigenvalues force even multiplicity, and blocks belonging to
different connected components of $`S`$ intertwine only through
shared eigenvalues.

Finite-component reduction is false.  The explicit order-$`18`$
circulant signing
```math
(0,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,1,1,-1,-1,-1,-1)
```
satisfies
```math
A^2
=17I+4(P^2+P^5+P^{13}+P^{16})-8P^9,
\tag{10.144}
```
whose off-diagonal support is connected of degree five.

Fixed-threshold finite-type reduction is also false.  A
block-palindromic circulant family has limiting correlation kernel
```math
q_d=(1-|d|/L)_+,
\tag{10.145}
```
so its threshold graphs are connected, bounded degree, unbounded,
and hyperfinite.

Hyperfiniteness itself is false without an extra spectral hypothesis.
An inverse-orbit factor on bounded-degree Cayley expanders produces
symmetric group-circulant flat signings whose threshold correlation
graphs contain those expanders.  The present spectral certificate is
only
```math
\|B^2\|_{\rm op}=O_D(\log n),
\tag{10.146}
```
not $`O_D(1)`$.  Therefore this construction refutes purely
structural amenability, but not hyperfiniteness under the
spectral-tameness inherited from a Pietsch-regularized near-minimizer.
Removing the logarithm requires a noncommutative trace-moment or
discrepancy theorem.  See
`flat_square_root_sparse_classification.md` and
`compatible_graphing_hyperfiniteness.md`.

#### 10.32.9 Current frontier and withdrawn claims

The following steps are proved and reusable:

1. centered-width lift monotonicity (10.119);
2. the exact residual/traffic criterion (10.121)--(10.125);
3. global block replacement and induced-block exchange
   (10.128)--(10.130);
4. the capped-cut tangent identity (10.132)--(10.133);
5. the flat feature-covariance degeneracy
   (10.136)--(10.137);
6. bounded-degree reduction and the uniform small-edge theorem
   (10.140).

The remaining proof targets are:

1. an entropy-sensitive or chained sign replacement converting the
   conditional kernel (10.136) into a centered internal fill;
2. an inverse theorem localizing a flat conditional-kernel direction
   into a replaceable induced vertex block;
3. a positive compatible-homotopy identity, or a bounded-degree
   graphing proof that avoids homotopy;
4. a spectral-tame classification of compatible graphings.

Explicitly withdrawn or falsified:

* exact finite equality as a universal centered-width gadget;
* independent traffic/refill;
* edge-local and triangle-only midpoint balancing;
* exact-cap full-rank rigidity;
* covariance alone as an entropy-rate theorem;
* finite-component and unconditional hyperfinite graphing reduction;
* the old neighbor-$`(S,R)`$ Hessian and the apparent size-four
  homotopy counterexample.

None of the open bullets above is silently assumed in later work.

#### 10.32.10 Status table and sidecar index

| Object | Status | Exact output or obstruction |
|---|---|---|
| Compressed centered-width range | Proved | $`W(B)\ge s^{3/2}W(A)`$ |
| Triangle four-lift equality | Falsified exactly | best equality candidate $`W=20>16`$ |
| Exact-sum traffic variance | Proved | $`\kappa_s=s^3/(s+1)`$, best square coefficient $`5/64`$ |
| Independent internal refill | Falsified asymptotically | union sum $`\ge e^{(\log2-1/4-o(1))n}`$ |
| Global block replacement | Proved | identities (10.128)--(10.130) |
| Local/triangle midpoint route | Falsified | local and correct-scale planted counterexamples |
| Dual capped-cut tangent | Proved | identities (10.132)--(10.133) |
| Vertex covariance entropy | Proved but insufficient | bound (10.135) misses the traffic exponent |
| Internal feature covariance | Proved structural reduction | flat kernel/low-eigen direction (10.136)--(10.137) |
| Exact-cap rank rigidity | Falsified computationally | complete certified audit through $`n=8`$ |
| Thick-cap feature rank | Proved but insufficient | full rank and $`\sigma_{\min}=2`$, wrong slack scale |
| Bounded-degree small-edge rounding | Proved | positive uniform gain (10.140) |
| Compatible homotopy monotonicity | Open | exact derivative has indefinite term (10.142) |
| Finite-type/hyperfinite graphing | Falsified structurally | circulant and Cayley-expander constructions |
| Spectral-tame graphing classification | Open | nonamenable example still loses $`O(\log n)`$ |
| Limit of $`F_n/n^{3/2}`$ | Open | no step above proves convergence |

Primary sidecars for this block:

* `centered_width_amplification_reboot.md`
* `asymptotic_centered_width_recovery.md`
* `spectral_filter_microblock_recovery.md`
* `internal_resampling_traffic_no_go.md`
* `global_block_replacement_midpoint.md`
* `cut_triangle_midpoint_balancing.md`
* `dual_cut_cone_diffuse_midpoint.md`
* `cut_cone_tangent_entropy.md`
* `endpoint_slack_ldp.md`
* `thick_near_cap_rank.md`
* `bounded_degree_compatible_graphing.md`
* `flat_square_root_sparse_classification.md`
* `compatible_graphing_hyperfiniteness.md`

Verification programs and data explicitly cited above:

* `verify_centered_width_equality_no_go.py`
* `cut_cap_face_rank_audit.cpp`
* `_homotopy_finite_channel.py`
* `_sparse_flat_square_search.py`
* `_circulant_block_palindrome.py`

### 10.33 One-sided energy-product balance

Let
```math
P(A)=\max_xx^\top Ax,\qquad
N(A)=-\min_xx^\top Ax,\qquad R(A)=P(A)+N(A).
```
Switch a positive ground state to $`\mathbf1`$, and form the graph
$`G`$ whose edges are the $`+1`$ entries of $`A`$.  Its density is
```math
p=\frac12+\frac{P(A)}{2n(n-1)}.
```
For every $`U\subset[n]`$, randomizing all spins outside $`U`$ gives
```math
-\frac{N(A)}2
\le\sum_{\{i,j\}\subset U}a_{ij}
\le\frac{P(A)}2.
```
Consequently the positive and negative induced discrepancies of $`G`$
satisfy
```math
\operatorname{disc}^{+}(G)\le\frac{P(A)}4,\qquad
\operatorname{disc}^{-}(G)\le\frac{R(A)}4. \tag{10.147}
```

Bollobás--Scott's product theorem
```math
\operatorname{disc}^{+}(G)\operatorname{disc}^{-}(G)
\ge\frac{p(1-p)n^3}{6400}
```
therefore yields
```math
P(A)R(A)\ge
\frac{\left(1-\left(P(A)/n(n-1)\right)^2\right)n^3}{1600}.
\tag{10.148}
```
Applying the same argument to $`-A`$ gives the analogous inequality
with $`N(A)`$ in place of $`P(A)`$.

In particular, along every sequence with $`Q(A)=O(n^{3/2})`$,
```math
\boxed{
P(A)N(A)\ge(1-o(1))\frac{n^3}{3200}.
} \tag{10.149}
```
Thus every asymptotically competitive signing has both one-sided
extrema of order $`n^{3/2}`$; severe orientation imbalance is
impossible.  This is now an input to the macroscopic-closure,
repair--purification, affine-ground rigidity, and centered-width
routes.  It does not alone provide scale transfer, because it contains
no localization among principal blocks.

The complete derivation, including all factors of two, is in
`one_sided_energy_product.md`.

### 10.34 Post-product checkpoint: strict deletion, integer cap rigidity,
sparse stationarity, and a graphing counterexample

This section records the next audited wave.  None of its statements
proves convergence, but it both strengthens the scale inequalities and
removes several plausible false routes.

#### 10.34.1 Strict macroscopic deletion

Use the one-copy notation
```math
q(A)=\max_x|H_A(x)|,\qquad
\ell(A)=\min\{\max H_A,-\min H_A\}.
```
For a principal split
```math
A=\begin{pmatrix}B&C\\C^\mathsf T&D\end{pmatrix},
```
relative block switching gives the exact orientation-balanced
superadditivity
```math
q(A)\ge q(B)+\ell(D),\qquad
q(A)\ge q(D)+\ell(B).                 \tag{10.150}
```
The one-sided product theorem implies, uniformly for competitive
order-$`m`$ signings,
```math
\ell(A)\ge(1-o(1))\frac{m^3}{12800q(A)}. \tag{10.151}
```
Consequently every fixed-ratio split satisfies
```math
q(A)\ge q(B)+(1-o(1))
\frac{|D|^3}{12800q(D)},              \tag{10.152}
```
and symmetrically.  Applied to an order-$`n`$ minimizer, this gives
```math
\boxed{
M_{n-k}+(1-o(1))\frac{k^3}{12800M_n}\le M_n,
}
\qquad
\boxed{
M_n^2\ge M_{n-k}^2+(1-o(1))\frac{k^3}{12800}.
}                                                   \tag{10.153}
```
The $`o(1)`$ is uniform when both blocks have order at least
$`\delta n`$.

This is a genuine strict scale inequality, but it does not settle the
limit.  Its payment is cubic in the deleted fraction, whereas
$`n^{3/2}`$-homogeneity requires a linear tangent payment.  An explicit
slowly oscillating sequence satisfies (10.153), monotonicity, local
continuity, the rigorous cage, and the corresponding one-sided block
inequalities.  Thus a boundary-traffic or energy-layer theorem is still
essential.  The full proof and scalar obstruction are in
`macroscopic_closure_block_dichotomy.md`.

#### 10.34.2 Affine cap clouds and integer replacement entropy

For an affine cut cloud with vertex types
$`\tau_i\in\mathbb F_2^d`$, group an edge $`e=ij`$ by its character
$`\lambda(e)=\tau_i+\tau_j`$, and write $`\mathcal E_\lambda`$ for the
corresponding fibre.  For every edge block $`T`$ and vector $`u`$ on
$`T`$, the conditional feature covariance is exactly
```math
\boxed{
u^\mathsf T\Sigma_{T\mid T^c}u
=
\sum_{\substack{\lambda\ne0\\
\mathcal E_\lambda\cap T^c=\varnothing}}
\left(\sum_{e\in T\cap\mathcal E_\lambda}u_e\right)^2.
}                                                   \tag{10.154}
```
Thus a zero Schur complement can arise merely because a remote copy of
each character exists.  It need not create any local sign freedom.

If a replacement $`\beta`$ is required to preserve the complete affine
cap profile, its exact integer constraints are
```math
\sum_{e\in T\cap\mathcal E_\lambda}\beta_e
=
\sum_{e\in T\cap\mathcal E_\lambda}a_e
\quad\text{for every }\lambda,                      \tag{10.155}
```
and the number of replacements is
```math
\prod_\lambda
\binom{|T\cap\mathcal E_\lambda|}
{|\{e\in T\cap\mathcal E_\lambda:a_e=1\}|}.          \tag{10.156}
```
This defines the fibre replacement entropy missing from the real
covariance formulation.

There is a linear-entropy, cut-triangle-compatible proportional-block
example in which (10.154) vanishes for every $`T`$ but (10.155) has a
unique solution on every transversal.  The example has a quadratic
Boolean witness, so it does not refute a theorem using competitive
scale and global extremality.

A partial competitive rigidity theorem was proved.  If an affine
positive-ground family has $`m`$ occupied type classes of size exactly
two, then, in one-copy normalization,
```math
\boxed{
q(A)\ge(1-o(1))\frac{m^2}{4000}.
}                                                   \tag{10.157}
```
Hence a competitive signing has only $`O(n^{3/4})`$ such classes.
Odd-multiplicity occupied types also obey an exact parity constraint:
every represented nonzero pair sum occurs an even number of times, so
the odd types contain additive parallelogram collisions.  Converting
those collisions, or even classes of size at least four, into mixed-sign
replacement freedom remains open.  See
`affine_ground_subspace_rigidity.md`.

#### 10.34.3 Why one conditioned tangent direction is insufficient

The natural correlated refill chooses
$`\beta\in\{\pm1\}^m`$ uniformly subject to
$`\langle\beta,h\rangle=0`$.  For a cap feature $`\phi`$, put
$`u=\langle h,\phi\rangle`$.  Then
```math
K\sim\operatorname{Hyp}\left(m,\frac{m+u}{2},\frac m2\right),
\qquad
\langle\beta,\phi\rangle=m+u-4K,                    \tag{10.158}
```
and
```math
\operatorname{Var}\langle\beta,\phi\rangle
=\frac{m^2-u^2}{m-1}.                               \tag{10.159}
```
At $`m=\Theta(n^2)`$ and $`u,t=O(n^{3/2})`$, the exact
hypergeometric exponent satisfies
```math
mJ(u/m,t/m)=\frac{t^2}{2m}+O(1).                    \tag{10.160}
```
Thus rank-one half-flip conditioning has the same leading speed-$`n`$
tail exponent as independent refill.  It changes the known traffic
deficit only subexponentially.  A constant reduction of the
$`\Theta(n^2)`$ feature variance requires
```math
\boxed{\Omega(n)\ \text{essentially independent low-variance
directions}.}                                       \tag{10.161}
```
The exact variational bridge also shows that failure of sparse repair
lower-bounds the same one-replica gap partition function that
independent purification needs to upper-bound.  A successful synthesis
must therefore be genuinely multiconstraint or chained.  See
`repair_purification_duality.md`.

#### 10.34.4 Sparse cap shaving and Paley orbit stationarity

Let edges be flipped independently with
```math
p_e=\delta(1+\lambda a_e\phi_e).
```
Uniformly over the full Boolean cube, Bernstein's inequality gives
```math
H_B(x)
=(1-2\delta)H_A(x)
-2\delta\lambda\langle\phi,z(x)\rangle
+O\!\left(\sqrt{\delta n^3}+n\right).                \tag{10.162}
```
For $`\delta=\Theta(n^{-1/2})`$, only
$`\Theta(n^{3/2})`$ edges are changed and the stochastic error is
$`O(n^{5/4})=o(n^{3/2})`$.  Therefore the full-cube metric entropy is
not the obstruction to sparse shaving; the deterministic cap geometry
is.

For exact minimizers this yields a quantitative convex-balance
certificate for every oriented thick-cap cloud.  More decisively, let
$`A`$ be a resonant Paley signing of prime order $`p`$, let $`x`$ be a
square-wave cap state with energy $`E`$, and let $`B`$ differ from
$`A`$ on $`m`$ edges.  Averaging over the affine quadratic-residue
orbit gives the exact lower bound
```math
M(B)\ge
E-\frac{2m(E+|T|)}{\binom p2},                       \tag{10.163}
```
with $`T`$ the opposite orbit energy.  In particular,
```math
m=o(p^2)
\quad\Longrightarrow\quad
M(B)\ge(1/2-o(1))p^{3/2},                            \tag{10.164}
```
and at Hamming density $`\delta`$,
```math
\frac{M(B)}{p^{3/2}}\ge\frac12-\delta-o(1).           \tag{10.165}
```
Hence no sparse perturbation can remove Paley resonance; a successful
perturbative construction must change $`\Omega(p^2)`$ edges.  See
`sparse_cap_shaving.md`.

#### 10.34.5 Spectral tameness does not force hyperfiniteness

For inverse-orbit sign convolutions on a finite group $`G`$, Schur
orthogonality and a sphere-net argument remove the previous logarithmic
loss.  If the dependency degree is bounded and
```math
\frac{d_{\min}(G)}{\log|G|}\longrightarrow\infty,
```
then one realization satisfies
```math
\boxed{\|A\|_{\rm op}=O(\sqrt{|G|})}.                 \tag{10.166}
```
This applies to $`G=\mathrm{PSL}_2(q)`$.  The same realization retains
the prescribed positive correlations along the Cayley generators, so
the fixed-threshold graph of
```math
\left(A/\sqrt{|G|-1}\right)^2
```
contains a bounded-degree expander, while the normalized square has
uniformly bounded operator norm.

Therefore spectrally tame compatible graphings need not be hyperfinite,
even at the competitive $`O(n^{3/2})`$ scale.  Any positive graphing
classification must use near-optimality or another condition stronger
than flatness, compatibility, bounded threshold degree, and bounded
operator norm.  See
`quasirandom_convolution_spectral_tameness.md`.

#### 10.34.6 Active proof targets after this wave

The live routes are now:

1. a tangent-scale principal restriction theorem using cross-block
   energy on actual near-ground layers;
2. a common-law theorem forcing $`\Theta(n)`$ independent tangent
   directions from many overlapping block certificates;
3. a projective compression/rounding theorem transferring a large
   near-minimizer to every smaller macroscopic order.

The affine bounded-multiplicity generalization and any Paley
perturbation changing only $`o(n^2)`$ edges are stopped pending a new
idea.

### 10.35 Projective compression, adaptive restriction, and the
traffic--descent frontier

This checkpoint records the audited conclusions of the next wave.
The main outcome is not yet a convergence proof.  It is a reduction
of the surviving difficulty to a joint inverse theorem for child
ground states and boundary traffic.

#### 10.35.1 Exact projective width transfer and its integrality wall

Let $`B`$ have order $`N=ns`$, partition its vertices into $`n`$
fibres of size $`s`$, and choose a Boolean mode
$`\sigma^{(a)}`$ in each fibre.  If
```math
S_{ab}=\sum_{i\in V_a,\ j\in V_b}
b_{ij}\sigma_i^{(a)}\sigma_j^{(b)},\qquad
C_{ab}=S_{ab}/s^{3/2},
```
then the block code satisfies the exact identity
```math
H_B(z(x))=I_\sigma+s^{3/2}H_C(x).
\tag{10.167}
```
Consequently,
```math
\boxed{W(C)\le W(B)/s^{3/2}.}                       \tag{10.168}
```
This has exactly the normalization needed for scale transfer.

The ordinary random projection nevertheless collapses:
```math
\mathbb E\sum_{a<b}C_{ab}^2=\frac1s\binom n2.
\tag{10.169}
```
For $`s/n\to\infty`$, one can therefore choose the modes so that
```math
M(C)/n^{3/2}\longrightarrow0.
\tag{10.170}
```
Rounding this fractional matrix back to signs has a compulsory
leading-order cost.  For every symmetric zero-diagonal real matrix,
```math
\boxed{
W(D)\ge\frac1{4\sqrt2}
\sum_i\left(\sum_{j\ne i}d_{ij}^2\right)^{1/2},
}
\tag{10.171}
```
and, when $`|d_{ij}|\le2`$,
```math
\boxed{
W(D)\ge\frac{\|D\|_F^2}{8\sqrt{2n}}.
}
\tag{10.172}
```
Thus sign rounding after a random projection costs
$`\Omega(n^{3/2})`$.  For a conference input, avoiding this loss is
equivalent to finding an almost invariant subspace spanned by
disjoint-support Boolean vectors whose compressed involution is
flat.  No such structural theorem follows from low Boolean quadratic
norm alone.  See `projective_compression_integrality_gap.md`.

#### 10.35.2 Exact principal-extension inequalities

For a principal split
```math
A=\begin{pmatrix}B&C\\C^\mathsf T&D\end{pmatrix}
```
and $`x\in\{\pm1\}^{V(B)}`$, orthogonality of the linear and quadratic
Walsh levels plus the Bhatia--Davis inequality gives
```math
\boxed{
H_B(x)^2+\binom{|D|}{2}+\|C^\mathsf Tx\|_2^2
\le M(A)^2.
}
\tag{10.173}
```
For a positive $`B`$-ground state there is also the tangent-scale
$`\ell_1`$ form
```math
\boxed{
M(A)\ge H_B(x)-M(D)+\|C^\mathsf Tx\|_1.
}
\tag{10.174}
```
A nonadaptive random restriction makes the last term
$`\Theta(n^{3/2})`$ for every fixed spin.  The obstruction is exactly
the reversed quantifier: the ground state of the restriction can
depend on the selected vertex set.

For a full spin $`z`$, define
```math
\mathcal E_z^{\rm gs}(L)=
\left\{
S\in\binom{[N]}m:
z_S\in\mathcal G(A[S]),\
\|A_{S,S^c}^\mathsf Tz_S\|_1\le L
\right\}.
\tag{10.175}
```
If every $`m`$-subset has some child ground state with traffic at
most $`L`$, exact incidence counting yields
```math
\boxed{
\max_z|\mathcal E_z^{\rm gs}(L)|
\ge2^{-m}\binom Nm.
}
\tag{10.176}
```
At $`m=N/2`$, the exceptional family must therefore have entropy
$`(\log2)N/2+o(N)`$.

The ground condition in (10.175) is essential.  A paired signing
built from an order-$`d`$ signing $`C`$ has energy
```math
H_{\mathcal L(C)}(x)
=d-2|R|+4H_{C[R]}(y_R),
\tag{10.177}
```
and its unions of complete pairs saturate the traffic-only incidence
threshold.  Uniform pair states, however, are not child ground
states.  In the actual ground orientation the model descends exactly:
```math
\boxed{
M(\mathcal L(C))=d+4M(C).
}
\tag{10.178}
```
Taking one representative from each pair recovers $`C`$, and
```math
\boxed{
c(\mathcal L(C))
=\sqrt2\,c(C)+O(d^{-1/2}),
\qquad c(A)=M(A)/|A|^{3/2}.
}
\tag{10.179}
```
Because $`2c_*>1/2`$, where
$`c_*=0.336493364431\ldots`$, two consecutive exact pair exceptions
cannot occur along a minimizing sequence.  The paired inverse class
therefore has bounded depth.  See
`principal_restriction_tangent_audit.md`.

#### 10.35.3 One common near-cap Gibbs law

The previously separate blockwise dual measures can be replaced by
one law.  For
```math
s_e(x,y)=\frac{x_ix_j-y_iy_j}{2},\qquad
\mathcal Z_A(\lambda)
=\sum_{x,y}\exp\!\left(\lambda\sum_ea_es_e(x,y)\right),
```
choose $`A_\lambda`$ minimizing $`\mathcal Z_A(\lambda)`$.  Its Gibbs
law $`\mu_\lambda`$ obeys, simultaneously for every edge,
```math
\boxed{
\left(a_e\mathbb E_{\mu_\lambda}s_e\right)_+
\le\tanh\lambda\le\lambda.
}
\tag{10.180}
```
Hence the same law controls every edge block $`T`$:
```math
\sum_{e\in T}
\left(a_e\mathbb E_{\mu_\lambda}s_e\right)_+
\le\lambda|T|.
\tag{10.181}
```
Moreover,
```math
W(A_\lambda)\le W_n+\frac{2n\log2}{\lambda}.
\tag{10.182}
```
Taking $`\lambda=b_n/\sqrt n`$, with $`b_n\to\infty`$ arbitrarily
slowly, and conditioning on a common
$`O(n^{3/2}/b_n)`$-cap shell preserves (10.180) up to
$`O(2^{-n})`$.  Thus every $`O(n)`$-edge graphing has
$`o(n)`$ aligned gradient under one and the same near-cap law.

This does not yet imply tangent rank.  There are low-rank Boolean
clouds with the same scalar balance but no cut-triangle structure.
Even the triangle identities alone are insufficient: a partition
into $`\sqrt n`$ groups and the complementary group cuts gives a
genuine cut-compatible cloud of rank $`O(\sqrt n)`$, cap-scale
feature mass, and coordinate means $`O(n^{-1/2})`$.  What this model
lacks is simultaneous near-activity for one globally competitive
signing.  Global near-minimality, not merely cut compatibility, is
the indispensable remaining input.  See
`common_gibbs_cap_law_rank.md`.

#### 10.35.4 Exact affine positive-ground closure

Suppose
```math
z_i(w)=(-1)^{\tau_i\cdot w}
\quad(w\in\mathbb F_2^d)
```
is an affine family of positive global ground states.  Let $`V_\phi`$
be the occupied type classes and put
```math
b_{\phi\psi}
=\sum_{i\in V_\phi,\ j\in V_\psi}a_{ij}.
```
Ground maximality applied to every union of complete types, followed
by Fourier averaging, gives the exact quotient closure
```math
\boxed{b_{\phi\psi}=0\quad(\phi\ne\psi).}
\tag{10.183}
```
In fact every type-constant spin is then a positive ground state.
For every union $`U`$ of type classes,
```math
\boxed{
p(A[U])=\sum_{V_\phi\subset U}p(A[V_\phi]),
\qquad
\nu(A[U])\ge\sum_{V_\phi\subset U}\nu(A[V_\phi]),
}
\tag{10.184}
```
and hence centered widths obey the same superadditive budget.
Vertexwise singleton stability also gives
```math
\boxed{
r_{i,\phi}\ge\sum_{\psi\ne\phi}|r_{i,\psi}|.
}
\tag{10.185}
```

Let $`k_\phi=|V_\phi|`$ and
```math
S_0=\sum_\phi\binom{k_\phi}{2}.
```
Since $`p(A)\le S_0`$, the one-sided product theorem implies, for
$`M(A)\le Cn^{3/2}`$,
```math
\boxed{
S_0\ge\frac{(1-o(1))n^{3/2}}{12800C},
\qquad
\max_\phi k_\phi
\ge1+\frac{(1-o(1))\sqrt n}{6400C}.
}
\tag{10.186}
```
The conclusion is hereditary.  If $`U`$ is the union of all types
with $`k_\phi\le t\sqrt n`$, and $`|U|=\beta n`$, then applying
(10.186) to $`A[U]`$ yields
```math
\boxed{
\beta\le\sqrt{6400Ct}+o(1).
}
\tag{10.187}
```
Thus, as $`t\downarrow0`$, all but $`o(n)`$ vertices lie in
$`\sqrt n`$-scale or larger type classes.

There is an exact trichotomy:

1. a nearly macroscopic type gives a scale-preserving centered-width
   descent;
2. mixed internal classes carry $`\Omega(n^{3/2})`$ exact
   cap-preserving replacement entropy;
3. an $`O(\sqrt n)`$-type pure mesoscopic core carries
   $`\Theta(n^{3/2})`$ actual positive energy, with all cross
   interaction in the double-centered ANOVA residual.

The third branch is sharp at the spectral ceiling.  For Hadamard
order $`k`$, an explicit signing on $`n=k^2`$ vertices, split into
$`k`$ all-positive classes of size $`k`$, has
```math
p(A)=\frac{n(k-1)}2,\qquad
\nu(A)=\frac{n(k+1)}2,
\qquad
\boxed{\frac{M(A)}{n^{3/2}}=\frac12+\frac1{2k}.}
\tag{10.188}
```
Its residual operator is a signed involution on directed type pairs.
Therefore the $`\sqrt n`$ multiplicity threshold is sharp, and any
elimination of the mesoscopic branch must use a quantitative gap
below $`1/2`$.  See `affine_type_closure_recursion.md`.

#### 10.35.5 Current traffic--descent target

After these audits, the most concrete remaining statement is:

> For a competitive signing, either a macroscopic principal
> restriction has tangent-scale boundary traffic on all of its child
> ground states, or an incidence-saturating family from (10.176)
> forces an affine/type quotient which admits a normalized principal
> descent.  The pure mesoscopic quotient may persist only at the
> $`1/2`$ ceiling.

The exact-zero form is already rigid.  After switching the common
full spin to $`1`$, zero boundary traffic and the child-ground
condition say
```math
(A1_S)_{S^c}=0,\qquad (A1_S)_S\ge0.
\tag{10.189}
```
Thus a single sign matrix maps exponentially many half-support
indicators back into their own coordinate faces.  The paired model
does this and descends by (10.178); generic sign matrices do not.
The active inverse problem is to prove that every
incidence-saturating family with (10.189) has the same bounded-depth
quotient/descent behavior, and then stabilize the result for
traffic $`o(n^{3/2})`$.

### 10.36 The affine branch reaches a sharp $`1/2`$ threshold

The exact affine analysis has now been completed far enough to
separate its algebraic and analytic parts.

#### 10.36.1 Recursion or quadratic replacement freedom

For occupied affine types of sizes $`k_\phi`$, every intertype block
has total zero by (10.183).  It can therefore be replaced
independently by any balanced signing of the same dimensions without
changing a single type-constant cap energy.  The logarithm of the
complete replacement pool is
```math
\mathscr R_{\rm cross}
=
\sum_{\phi<\psi}
\log\binom{k_\phi k_\psi}{k_\phi k_\psi/2}.
\tag{10.190}
```
Competitive exact affine clouds have only $`O_C(n^{5/6})`$ occupied
types.  Hence, if $`K=\max_\phi k_\phi\le(1-\delta)n`$,
```math
\boxed{
\mathscr R_{\rm cross}
=
(\log2)E_{\rm cross}+o_C(n^2)
\ge
\left(\frac{\delta\log2}{2}-o_C(1)\right)n^2.
}
\tag{10.191}
```
If instead $`K=n-o(n)`$, the largest type is a centered-width
principal restriction with asymptotically no normalized loss.
Thus the exact algebraic dichotomy is
```math
\boxed{
\text{scale-preserving near-total type descent}
\quad\text{or}\quad
\exp(\Omega(n^2))\text{ exact cap-profile replacements}.
}
\tag{10.192}
```

The replacement entropy alone is insufficient.  A refill which
preserves only each block total has, for a fixed spin with type
magnetizations $`s_\phi`$, exact variance
```math
\boxed{
V(x)=
\sum_{\phi<\psi}
\frac{(k_\phi k_\psi)^2-s_\phi^2s_\psi^2}
{k_\phi k_\psi-1}.
}
\tag{10.193}
```
Its leading tails still have speed $`n`$, exactly the speed of the
$`2^n`$ Boolean witness set.

#### 10.36.2 Doubly balanced refills

For even $`k,\ell`$, let a block be uniform among all sign matrices
with every row and column sum zero.  With
```math
P_k=I_k-k^{-1}J_k,
```
its covariance is exactly
```math
\boxed{
\operatorname{Cov}(\operatorname{vec}B)
=
\frac{k\ell}{(k-1)(\ell-1)}
(P_k\otimes P_\ell).
}
\tag{10.194}
```
For internally all-positive types, put
```math
D_\phi(x)=\frac{k_\phi^2-s_\phi^2}{2}.
```
Independent doubly balanced blocks then have
```math
\boxed{
V_{\rm db}(x)
=4\sum_{\phi<\psi}
\frac{D_\phi(x)D_\psi(x)}
{(k_\phi-1)(k_\psi-1)}
\le
2\left(\sum_\phi\frac{D_\phi(x)}{k_\phi-1}\right)^2.
}
\tag{10.195}
```
This multiplicative deficit is the first refill law to remove every
one-sided local channel: the variance vanishes if either endpoint
type is constant.

Nevertheless, the direct first-moment route is rigorously false.  A
random checkerboard subensemble of doubly balanced blocks and the
shell in which every type spin is balanced give
```math
\Pr\!\left[Y\ge Bn^{3/2}\right]
=\exp(-(B^2+o(1))n),
\tag{10.196}
```
while that shell contains
```math
\exp((\log2-o(1))n)
\tag{10.197}
```
states.  Its direct traffic sum therefore diverges whenever
```math
\boxed{
B<\sqrt{\log2}=0.832554611\ldots.
}
\tag{10.198}
```
Thus neither total-balanced nor row--column-balanced refills can be
certified at $`1/2`$ by a union bound over individual Boolean states.
Any further probabilistic use must exploit overlaps or an inverse
description of the exceptional refill set.

#### 10.36.3 The conditional flat-residual envelope

There is an exact deterministic statement explaining why $`1/2`$
appears.  Split the vertices into all-positive types $`V_i`$ of
sizes $`k_i`$, put zero row and column sums on every cross block, and
write $`B`$ for the complete cross-block matrix.  If
```math
\|B\|_{\rm op}\le L\sqrt n,
```
then
```math
\boxed{
M(A)\le\frac12\max\left\{
\sum_i k_i\max\{k_i,L\sqrt n\}-n,\,
Ln^{3/2}+n
\right\}.
}
\tag{10.199}
```
In particular,
```math
L=1+o(1),\qquad
\max_i k_i\le(1+o(1))\sqrt n
\quad\Longrightarrow\quad
M(A)\le(1/2+o(1))n^{3/2}.
\tag{10.200}
```
The target $`L=1`$ is sharp: on the residual space,
```math
\|B\|_{\rm op}
\ge
\sqrt{\frac{n^2-\sum_i k_i^2}{n-q}}
=(1-o(1))\sqrt n
\tag{10.201}
```
in the competitive mesoscopic regime.  The square Hadamard
construction attains equality.  General existence at this floor is
a fusion-frame/weighing-matrix design problem.

The resulting exact microprofile formula is also useful.  If
```math
\Delta_\phi(r_\phi)
=p(A[V_\phi])-H_{A[V_\phi]}(r_\phi)
```
and
```math
c_{\phi\psi}(r)
=r_\phi^\mathsf TA_{\phi\psi}r_\psi,
```
then, for the weighted type quotient $`C(r)`$,
```math
\boxed{
p(A)+\nu(A)
=
\max_r\left\{
\sum_\phi\Delta_\phi(r_\phi)+\nu(C(r))
\right\}.
}
\tag{10.202}
```
This is the precise residual optimization left by the affine branch:
local type deficits plus a smaller weighted negative-ground problem.

The current target is an inverse $`1/2`$-threshold theorem.  A
mesoscopic affine signing avoiding all Boolean witnesses above
$`(1/2-\varepsilon)n^{3/2}`$ should be forced close to the flat
residual structure in (10.199)--(10.201); that structure should in
turn force either $`\sqrt n`$-scale ceiling behavior or a principal
descent.  The direct annealed refill argument is stopped.

### 10.37 Exact low-traffic child grounds

The exact-zero form of the adaptive principal-restriction
obstruction now has a proof-level complementarity and uncrossing
theory.

After fixing the common full spin and switching it to $`1`$, call an
even $`m`$-set $`S`$ positive good when $`1_S`$ is a positive
absolute ground of $`A[S]`$ and its outside field vanishes.  Parity
upgrades the weak local inequalities to
```math
\boxed{
A1_S\ge1_S,\qquad
\operatorname{supp}(A1_S)=S.
}
\tag{10.203}
```
Thus every good support is a strict coordinate face mapped into
itself by $`A`$.

If $`S,T`$ are good and
```math
X=S\cap T,\quad
Y=S\setminus T,\quad
Z=T\setminus S,\quad
W=(S\cup T)^c,
```
then exact boundary cancellation and child-ground stability give
```math
\boxed{
e(X,Y)=e(X,Z)=-e(Y,Z)=:t\ge0,
}
\tag{10.204}
```
and, row by row on $`W`$,
```math
\boxed{
e(\{w\},Y)=e(\{w\},Z)=-e(\{w\},X).
}
\tag{10.205}
```
Consequently,
```math
M(A)\ge M(A[Y\cup Z])\ge t.
\tag{10.206}
```
The dichotomy is exact: $`t>0`$ exposes a negative interaction
between the symmetric-difference cells; $`t=0`$ creates additional
exact child grounds by flipping $`X`$ inside both restrictions.

For adjacent good sets
```math
T=S\setminus\{i\}\cup\{j\},
```
this becomes literal twin rigidity:
```math
\boxed{
a_{ij}=-1,\qquad
a_{ki}=a_{kj}\quad(k\notin S\cup\{j\}),
}
\tag{10.207}
```
and both exchanged local fields equal $`1`$.  A speed-$`n`$
constant-weight code can avoid adjacent pairs, so a local argument
alone is insufficient.

There is also a quantitative purification statement.  If every
member has traffic at most $`L_n`$, then, for the conditional
covariance $`\Sigma_j`$ of $`1_S`$ given $`j\notin S`$,
```math
\boxed{
\sum_jq_j\,a_j^\mathsf T\Sigma_ja_j\le mL_n.
}
\tag{10.208}
```
At even $`m`$, only $`L_n/2`$ outside rows per member can have a
nonzero field; at odd $`m`$, parity forces traffic at least $`n-m`$.
Thus $`o(n)`$ traffic is an $`o(n)`$-row perturbation of the exact
system.

The incidence threshold produces an exact multivariate
Littlewood--Offord atom.  For some $`t`$-set $`J`$,
```math
\boxed{
\Pr_{\xi\in\{0,1\}^{J^c}}
\left(A_{J,J^c}\xi=0\right)
\ge
2^{-(n-t)}|\mathcal F|
\frac{\binom{n-m}{t}}{\binom nt}.
}
\tag{10.209}
```
At $`m=n/2`$, $`|\mathcal F|\ge2^{-m}\binom nm`$, and
$`t=o(n)`$, the right side is
```math
2^{-n/2}
\exp\!\left(-O(t^2/n+\log n)\right).
\tag{10.210}
```
Rank alone is too weak at this atom size; a genuine inverse
Littlewood--Offord theorem must use the compatible flat rows and the
positive inequalities in (10.203).

Finally, both complete entropy-saturating pair models are now
excluded by the ground condition.

* If every middle-layer union of complete pairs is good, all
  interpair blocks are checkerboards and internal pair edges are
  positive.  Antiuniform pair spins then force a nonzero Boolean
  quadratic polynomial to be nonnegative everywhere, a
  contradiction.
* If every pair transversal is good, every interpair block is
  constant and the induced order-$`d`$ quotient $`C`$ satisfies
  $`\displaystyle \boxed{M(C)=d/2,} \tag{10.211}`$
  contradicting the universal $`\Omega(d^{3/2})`$ lower bound for
  large $`d`$.

Thus the obvious pair-union and transversal saturators are not the
missing exception.  Any surviving family at the incidence entropy
must be genuinely code-separated and nonlinear.  The live exact
target is an extremal stability theorem using (10.203) **together
with the full absolute child-ground inequalities**
```math
\left|H_{A[S]}(x)\right|
\le H_{A[S]}(\mathbf1)
\qquad(x\in\{\pm1\}^{S}).
\tag{10.211a}
```
A near-saturating middle-layer family should force one of the
pair/affine quotient structures already known to be impossible or
to descend.

Neither boundary cancellation nor strict complementarity alone is
enough.  For a balanced sign vector $`s`$, the signing
```math
A=I-ss^\mathsf T
```
satisfies
```math
A1_S=1_S
```
for all $`\binom{n/2}{n/4}^2=2^{n-o(n)}`$ type-balanced middle
sets.  Yet the child spin $`s_S`$ has energy
$`-\binom{|S|}{2}`$, so (10.211a) fails maximally.  These supports,
over all even sizes, attain the sharp central-binomial bound for
strict local maxima.  Consequently the absolute-ground condition,
not just positivity/support, is indispensable.  See
`inverse_low_traffic_descent.md`.

### 10.38 Margin-preserving affine refills and the exact half threshold

The nonrecursive affine branch has now been analyzed inside the
entire fibre that preserves every cross-block row and column margin.
Row domination gives the global margin budget
```math
\sum_{\phi<\psi}\left(
\sum_{i\in V_\phi}|r_{i,\psi}|
+\sum_{j\in V_\psi}|r_{j,\phi}|
\right)
\le 2p(A)=O(n^{3/2}).
\tag{10.212}
```
An alternating-rectangle count then gives
$`\Omega_\delta(n^2)`$ independent $`2\times2`$ margin-preserving
trades whenever no type contains more than $`(1-\delta)n`$
vertices.  Thus this branch really does have quadratic refill
entropy; the obstruction is analytic rather than a shortage of
admissible signings.

For a uniformly random doubly balanced $`k\times\ell`$ sign block,
the exact covariance is
```math
\boxed{
\operatorname{Cov}(\operatorname{vec}B)
=\frac{k\ell}{(k-1)(\ell-1)}
 (P_k\otimes P_\ell).
}
\tag{10.213}
```
For two Boolean replicas of overlap $`\rho`$, the corresponding
checkerboard second-moment exponent at cap $`c n^{3/2}`$ is
```math
\Phi_c(\rho)
=h\!\left(\frac{1-\rho}{2}\right)-\log2
+\frac{2c^2\rho^2}{1+\rho^2}.
\tag{10.214}
```
Pinsker's inequality shows $`\Phi_c(\rho)\le0`$ for
$`c\le1/2`$, with the overlap-zero saddle losing local stability
exactly above $`1/2`$.  Consequently, a generic isotropic refill
already contains half-scale Boolean witnesses.  Random refill is
therefore a sharp no-go for constructing a sub-half cap, not a
route to one.

There is a complementary deterministic fourth-moment theorem.
For internally pure affine types of sizes $`k_\alpha`$, normalized
type means $`u_\alpha`$, $`q`$ occupied types, and
```math
d_\alpha
=u_\alpha^\mathsf T(A^2-(n-1)I)u_\alpha,\qquad
\Delta_{\rm type}=\sum_\alpha d_\alpha^2,
```
one has
```math
\boxed{
2\sum_\alpha k_\alpha^2
\ge q(n-2)+3n-\sqrt{q\Delta_{\rm type}}.
}
\tag{10.215}
```
Writing $`q=y\sqrt n`$ and
$`\Delta_{\rm type}\le\delta n^{5/2}`$ gives
```math
\frac{p(A)}{n^{3/2}}
\ge
\max\left\{
\frac1{2y},\
\frac y4-\frac{\sqrt{y\delta}}4
\right\}-o(1).
\tag{10.216}
```
In particular, projected defect $`o(n^{5/2})`$ gives the
$`1/(2\sqrt2)`$ floor; adding negligible margin compensation gives
the $`1/2`$ floor.

Under exact flatness, a sub-half exception must instead carry a
quantitative family of biased columns.  If
$`p(A)=c n^{3/2}+o(n^{3/2})`$, then the
$`L^1`$-margin-weighted mean normalized bias is at least
```math
\gamma(c)=\frac1{4c^2}-1-o(1).
\tag{10.217}
```
At least a fraction
```math
\frac{\gamma(c)}{2-\gamma(c)}
\tag{10.218}
```
of the absolute margin mass lies on columns with normalized bias at
least $`\gamma(c)/2`$.  At the ROM value
$`c=\sqrt{15}/8`$, this means at least $`1/29-o(1)`$ of the
margin mass has bias at least $`1/30`$.  At
$`c=1/(2\sqrt2)`$, essentially all margin mass is carried by
nearly constant columns.

The precise surviving inverse target is now: turn this biased
incidence mass, or an anisotropic kernel of the rectangle-trade
covariance, into a Ferrers/paired quotient or a scale-preserving
principal descent.  The known fourth-moment stability estimate is
too weak to infer (10.215) from a merely sub-half Boolean cap:
it permits full defect $`O(n^4)`$, whereas the needed projected
scale is $`o(n^{5/2})`$.  See
`margin_preserving_affine_refill.md`.

### 10.39 Certified order-eight test of absolute-ground low traffic

The first complete finite test of the corrected low-traffic inverse
target is now certified.  Every signing of order $`8`$ was gauged to
have positive first row, leaving $`2^{21}`$ switching classes.
Exactly $`4200`$ of these have the optimal value $`M(A)=10`$.
For every such class, every common gauge $`z`$, and every
four-subset $`S`$, we tested both

```math
z_S\text{ is a positive absolute ground of }A[S],
\qquad
A_{S^c,S}z_S=0.
\tag{10.219}
```

The largest family satisfying (10.219) for one fixed $`(A,z)`$
has size
```math
\boxed{4.}
\tag{10.220}
```
The middle-layer incidence threshold is
```math
2^{-4}\binom84=4.375,
```
so its integer form is $`5`$: no optimal order-eight signing reaches
the threshold.  Two independent exhaustive implementations agree
on all three certification statistics:

```math
\#\{A:M(A)=10\}=4200,\qquad
\max|\mathcal F|=4.
\tag{10.221}
```

This is evidence, not an asymptotic theorem, but it cleanly separates
global minimizers from arbitrary signings: a nonoptimal order-eight
example with $`M(A)=14`$ has eight good four-sets.  Hence the desired
entropy deficit must use global near-minimality together with strict
complementarity; neither condition alone is sufficient.

### 10.40 Opposite-face uncrossing and capped-bilinear equality

The common-law route now has an exact same-signing rigidity package.
Switch a positive ground to $`1`$, write the resulting signing as
$`B`$, and define
```math
c(S)=\sum_{e\in\delta(S)}B_e.
```
If $`P=\max H_B`$, $`Q=-\min H_B`$, and
$`W=(P+Q)/2`$, then
```math
0\le c(S)\le W.
\tag{10.222}
```
The top face is the zero-level family $`\mathcal Z`$, while the
bottom face is the maximum-level family $`\mathcal M`$.

For
```math
I_B(S,T)=
\sum_eB_e1_{\delta(S)}(e)1_{\delta(T)}(e),
```
the pointwise cut identity gives
```math
c(S)+c(T)=c(S\triangle T)+2I_B(S,T).
\tag{10.223}
```
Consequently,
```math
\boxed{
\begin{array}{ll}
S,T\in\mathcal Z&\Longrightarrow I_B(S,T)\le0,\\
S,T\in\mathcal M&\Longrightarrow I_B(S,T)\ge W/2,\\
S\in\mathcal Z,\ T\in\mathcal M
&\Longrightarrow0\le I_B(S,T)\le W/2.
\end{array}}
\tag{10.224}
```
These are genuine opposite-face constraints; they do not follow
from PSD or from a generic cut-compatible correlation model.

If $`p_e,q_e`$ are the two face crossing probabilities and the common
edge balance is
```math
g_e=B_e(q_e-p_e)\le\lambda<1,
\qquad
\sum_eg_e=W,
```
then at least $`W/\lambda`$ positively aligned edges are stochastic
in one of the two faces.  If the two edge-feature faces have affine
dimensions $`d_+,d_-`$, all such edges are covered by
$`d_++d_-`$ same-face difference cuts.  One covering cut therefore
has
```math
|\delta(S)|\ge
\frac{W}{\lambda(d_++d_-)}.
\tag{10.225}
```

Every zero difference cut gives more than scalar width
superadditivity.  For a principal split $`S,S^c`$, let
$`\mathcal G_S,\mathcal G_{S^c}`$ be the two positive-ground sets,
let $`L_S,L_{S^c}`$ be their linear spans, and let $`C`$ be the
cross sign block.  Then all principal grounds close
Cartesianly and
```math
\boxed{
P_{L_S}CP_{L_{S^c}}=0,\qquad
\operatorname{rank}C
\le\operatorname{codim}L_S+\operatorname{codim}L_{S^c}.
}
\tag{10.226}
```
More generally,
```math
\left|u^\mathsf TCv\right|
\le
\bigl[p(A[S])-H_{A[S]}(u)\bigr]
+
\bigl[p(A[S^c])-H_{A[S^c]}(v)\bigr].
\tag{10.227}
```
Thus a full-rank macroscopic cross block forces a linear ground-span
deficit; small deficit forces an actual low-rank sign quotient.

There is also an infinite Hadamard-moment audit.  For the two cap
correlation matrices $`R^\pm`$,
```math
-Q\le\langle B,(R^\pm)^{\circ m}\rangle\le P
\qquad(m\ge1).
\tag{10.228}
```
Vanishing cubic closure defect makes the face support a torsor under
coordinatewise ternary products, hence an affine
$`\mathbb F_2`$ ground family.  The exact unresolved case is now
square-root-scale effective rank together with substantial cubic
closure defect and many crossing equality profiles in (10.227).
See `opposite_face_rigidity.md`.

### 10.41 Incidence-scale affine child families are impossible

The affine branch of the exact low-traffic inverse problem has now
been eliminated completely.

First, let
```math
\mathcal V=v_0+U\subseteq\mathbb F_2^n
```
be an affine subspace of dimension $`d`$, all of whose vectors have
the same Hamming weight.  Parametrize by $`w\in\mathbb F_2^d`$.
Every nonconstant coordinate has the form
```math
v_i(w)=
\frac{1-\sigma_i(-1)^{\lambda_i\cdot w}}2.
```
Constant weight and Fourier uniqueness force, in every represented
nonzero character fibre,
```math
\sum_{i:\lambda_i=\lambda}\sigma_i=0.
\tag{10.229}
```
If $`q`$ is the number of represented nonzero characters and $`r`$
the number of nonconstant coordinates, then
```math
\boxed{d\le q\le r/2\le n/2.}
\tag{10.230}
```
Equality $`d=n/2`$ is exactly the full transversal cube of
$`n/2`$ hidden opposite-coordinate pairs.

At the incidence threshold
```math
|\mathcal V|
\ge2^{-n/2}\binom n{n/2},
```
one has
```math
d\ge\frac n2-\frac12\log_2n-O(1).
\tag{10.231}
```
Choose $`d`$ independent character fibres and one opposite-phase
coordinate pair $`P_r=\{p_r,q_r\}`$ from each.  All remaining
coordinates form an exceptional set $`E`$ of size
```math
h=n-2d=O(\log n).
\tag{10.232}
```

Now suppose every support in $`\mathcal V`$ is positive good: the
restricted all-one spin is a positive absolute ground and every
outside field is zero.  Mark a basis index dirty if it occurs in an
exceptional character of Hamming weight one or two.  There are at
most $`2h`$ dirty indices.  For two clean indices $`r\ne s`$, restrict
the outside-field equation of $`p_r`$ to the hyperplane where
$`p_r`$ is absent.  The contribution of pair $`s`$ has singleton
Fourier coefficient
```math
\frac{a_{p_rp_s}-a_{p_rq_s}}2\,\chi_s.
```
An exceptional character can cancel it only if its support is
$`\{s\}`$ or $`\{r,s\}`$, which would make an index dirty.
Repeating for $`q_r`$, then using symmetry, proves
```math
\boxed{
A[P_r,P_s]=c_{rs}J_2
\quad\text{for every two clean pairs.}
}
\tag{10.233}
```

The constant Fourier coefficient in the same outside-field equation
gives
```math
\left|\sum_{\substack{s\ {\rm clean}\\s\ne r}}c_{rs}\right|
\le1+3h.
\tag{10.234}
```
Every child contains one representative of each clean pair and only
$`O(h)`$ other vertices.  Equations (10.233)--(10.234) therefore
bound its all-one energy by
```math
H_{A[S]}(\mathbf1)=O(n\log n).
\tag{10.235}
```
But this spin is an absolute child ground, while the universal lower
bound gives
```math
M(A[S])\ge(c_*+o(1))(n/2)^{3/2},
\qquad
c_*=0.336493364431\ldots.
\tag{10.236}
```
This is a contradiction.

Hence, for all sufficiently large even $`n`$,
```math
\boxed{
\text{no incidence-threshold affine family of positive-good
middle sets exists.}
}
\tag{10.237}
```
The surviving exact inverse obstruction must be genuinely
non-affine and code-separated.  See
`inverse_low_traffic_descent.md`.

### 10.42 Orthogonal margin channels and the half-hard endpoint

The biased-column branch of the affine refill problem has an exact
invariant-plane structure under conference flatness
$`A^2=(n-1)I`$.  Let $`V_\alpha`$ be an internally positive affine
type, $`k_\alpha=|V_\alpha|`$,
```math
u_\alpha=k_\alpha^{-1/2}1_{V_\alpha},
\qquad
d_\alpha=k_\alpha-1,
```
and
```math
g_\alpha=(I-UU^\mathsf T)Au_\alpha.
```
Then
```math
\boxed{
\langle g_\alpha,g_\beta\rangle=0\ (\alpha\ne\beta),
\quad
\|g_\alpha\|^2=n-1-d_\alpha^2,
}
\tag{10.238}
```
and
```math
\boxed{
Ag_\alpha=
\bigl(n-1-d_\alpha^2\bigr)u_\alpha-d_\alpha g_\alpha.
}
\tag{10.239}
```
Thus every type mean and its signed margin channel span an exact
$`A`$-invariant two-plane.

For a fixed sub-half flat exception, the quantitatively biased margin
mass has a rigorous dichotomy.

1. If it escapes through $`o(\sqrt n)`$-sized source types, the
   union of those types has $`o(n)`$ vertices; deleting it is a
   principal restriction with asymptotically no increase in
   normalized centered width.
2. Otherwise, $`\Omega(n)`$ high-bias incidences lie on
   $`\sqrt n`$-scale source types.  Row domination gives bounded
   incidence depth, and every such column has an exact paired-star
   decomposition: a constant column plus canceled opposite-sign
   edge pairs.

At the compensated endpoint
```math
p(A)=\left(\frac1{2\sqrt2}+o(1)\right)n^{3/2},
\tag{10.240}
```
equality forces
```math
q=(\sqrt2+o(1))\sqrt n,\qquad
k_\alpha=(1+o(1))\sqrt{n/2},
\tag{10.241}
```
and all but $`o(n)`$ vertices form a depth-one paired-star quotient.
The nonzero margin columns are $`o(\sqrt n)`$-close to constant and
pair by opposite signs.

This endpoint does **not** yield a favorable half-size principal
descent.  Rounding the orthogonal channel cube and using
(10.239) produces Boolean vectors $`x,y`$ with
```math
\langle x,y\rangle=o(n),
```
```math
x^\mathsf TAx
=\left(\frac1{\sqrt2}+o(1)\right)n^{3/2},
\quad
y^\mathsf TAy
=-\left(\frac1{\sqrt2}+o(1)\right)n^{3/2},
\quad
x^\mathsf TAy
=\left(\frac1{\sqrt2}+o(1)\right)n^{3/2}.
\tag{10.242}
```
Their agreement and disagreement sets $`I,J`$ both have size
$`n/2+o(n)`$, but
```math
\boxed{
\frac{M(A[I])}{|I|^{3/2}}\ge\frac12-o(1),
\qquad
\frac{M(A[J])}{|J|^{3/2}}\ge\frac12-o(1).
}
\tag{10.243}
```
The natural crossover faces are therefore half-hard.  The remaining
moderate-bias route must use global nesting or a different principal
selection, not the obvious invariant-channel split.  See
`biased_margin_descent.md`.

### 10.43 Capped-bilinear payment and exact repeated saturation

For a positive-ground zero-cut split
```math
A=\begin{pmatrix}B&C\\C^\mathsf T&D\end{pmatrix},
\qquad p(A)=p(B)+p(D),
```
put
```math
f(u)=p(B)-H_B(u),\qquad g(v)=p(D)-H_D(v).
```
The full cap gives
```math
|u^\mathsf TCv|\le f(u)+g(v).
\tag{10.244}
```
On the positive equality set
$`\mathcal E=\{(u,v):u^\mathsf TCv=f(u)+g(v)\}`$, two profiles obey
the exact cyclic-monotonicity identity
```math
(u-u')^\mathsf TC(v-v')
=s(u,v')+s(u',v)\ge0,
\tag{10.245}
```
where $`s(u,v)=f(u)+g(v)-u^\mathsf TCv`$.  If
$`\mathcal G_B,\mathcal G_D`$ are the principal positive ground
families, every equality profile also pays
```math
f(u)+g(v)
\ge |u^\mathsf TCv_0|+|u_0^\mathsf TCv|
\quad
(u_0\in\mathcal G_B,\ v_0\in\mathcal G_D).
\tag{10.246}
```

There is an exact global disintegration.  Let $`\mu`$ be a symmetric
law on the positive ground face, let
```math
R=\mathbb E_\mu XX^\mathsf T,
\qquad
\Delta_3=p(A)-\langle A,R^{\circ3}\rangle,
```
and take independent $`X,Y,Z\sim\mu`$.  The disagreement set
$`S=X\triangle Y`$, after switching by $`X`$, is a zero cut.  If
$`Z=(u,v)`$ across that cut, then
```math
\boxed{
p(A)-H_A(XYZ)
=2u^\mathsf TC_Sv
=2\bigl(f_S(u)+g_S(v)\bigr).
}
\tag{10.247}
```
Consequently
```math
\Delta_3=2\mathbb E_{X,Y}\pi(S),
\qquad
\pi(S)=\sum_{i\in S,\ j\notin S}a_{ij}R_{ij}\ge0.
\tag{10.248}
```
If $`r_{\rm eff}(R)\to\infty`$ and
$`\Delta_3\ge\delta n^{3/2}`$, some balanced macroscopic zero cut
has
```math
\pi(S)\ge(\delta/2-o(1))n^{3/2}.
\tag{10.249}
```

The exact range formula at such a cut is
```math
2W(A)=\max_{u,v}
\left\{f(u)+g(v)+|u^\mathsf TCv|\right\}.
\tag{10.250}
```
It yields a maximum, not an additive payment.  A direct check of an
optimal order-$`6`$ conference signing gives a $`2+4`$ zero cut with
child widths summing to $`W(A)=5`$ and $`\pi(S)=8/3`$, while the
cross term vanishes on every pair of negative principal grounds.
Thus
```math
W(A)\ge W(A[S])+W(A[S^c])+c\,\pi(S)
```
is false for every universal $`c>0`$.

Exact saturation nevertheless has a complete two-sided
classification.  If
```math
W(A)=W(B)+W(D),
```
then
```math
|u^\mathsf TCv|
\le
\min\{f(u)+g(v),\bar f(u)+\bar g(v)\},
\tag{10.251}
```
where
$`\bar f(u)=\nu(B)+H_B(u)`$ and
$`\bar g(v)=\nu(D)+H_D(v)`$.
Hence $`C`$ annihilates both the positive-ground Cartesian family
and the negative-ground Cartesian family.  Along a laminar tree of
sharp zero-cut saturations with leaves $`L_1,\ldots,L_\ell`$,
```math
\boxed{
p(A)=\sum_a p(A[L_a]),\quad
\nu(A)=\sum_a\nu(A[L_a]),\quad
W(A)=\sum_aW(A[L_a]).
}
\tag{10.252}
```
Both endpoint faces contain independent leaf-sign torsors; all
interleaf block totals vanish in their corresponding endpoint
gauges.  Low-rank cross blocks therefore reduce to explicit
bounded-pattern sign quotients.

The remaining wall is stability.  Near saturation supplies only
scalar near-annihilation of negative child-ground pairs:
```math
|u_-^\mathsf TCv_-|
\le2\bigl[W(A)-W(B)-W(D)\bigr].
\tag{10.253}
```
Without quantitative frame lower bounds this does not imply
approximate rank collapse, and exact additivity over many leaves
loses a factor $`1/\sqrt\ell`$ in the scalar lower bound.  See
`capped_bilinear_inverse.md`.

### 10.44 Exact insertion scale and the thick-cap counterexample

Adjoining one vertex with incident row $`b`$ has the exact cost
```math
\boxed{
E(A)=M(A)+\min_b\Delta_A(b),\qquad
\Delta_A(b)=\max_x
\bigl\{|b\cdot x|-[M(A)-|H_A(x)|]\bigr\}.
}
\tag{10.254}
```
The derivative-scale target is
```math
D_n=M_n\left[\left(1+\frac1n\right)^{3/2}-1\right].
\tag{10.255}
```
An estimate
```math
\Delta_A(b)\le D_n+r_n,\qquad
\sum_n r_n/n^{3/2}<\infty,
\tag{10.256}
```
for optimal $`A`$ would force convergence.  Merely
$`r_n=o(\sqrt n)`$ does not: normalized oscillations of the form
$`\sin(\log\log n)`$ survive such a one-sided derivative error.
Moreover a total $`o(\sqrt n)`$ increment is impossible.  The
current rigorous cage implies that in every sufficiently large
dyadic interval some order satisfies
```math
M_{n+1}-M_n
\ge(0.319433\ldots-o(1))\sqrt n.
\tag{10.257}
```

When the two exact endpoint faces are affine, their simultaneous
insertion discrepancy is exactly a signed bipartite-incidence
problem
```math
\kappa(G,\gamma)
=\min_{c\in\{\pm1\}^{E(G)}}
\max\{\|L(c)\|_1,\|R(c)\|_1\}.
\tag{10.258}
```
If $`\tau(G,\gamma)`$ is the signed frustration index, gauged
alternating Euler tours prove
```math
\boxed{\kappa(G,\gamma)\le2\tau(G,\gamma)+1,}
\tag{10.259}
```
with the final $`1`$ absent in the all-even case.  Connected
incidence components are simultaneous positive/negative
type-unions and hence exact principal zero-cut closures.

Exact endpoint balancing alone cannot control the thick cap.  For
each even $`k\ge8`$, take two positive $`k`$-cliques joined by a
balanced circulant sign block $`B`$, then flip one cross edge.  The
resulting order-$`n=2k`$ signing has one absolute-ground pair.
There is a row $`b`$ orthogonal to both absolute grounds, but a
former affine ground lies only four below the cap and has field
$`n`$.  Therefore
```math
\boxed{\Delta_A(b)\ge n-4.}
\tag{10.260}
```
This construction was independently brute-force checked at
$`k=8,10`$: the flipped matrices have respective absolute maxima
$`58,92=k(k-1)+2`$, each with exactly one ground pair.
A viable insertion theorem must choose $`b`$ from the full near-cap
hierarchy rather than from exact endpoint rank or discrepancy alone.
See `cap_discrepancy_insertion.md`.

### 10.45 First bounds in the non-affine code-separated branch

For a positive-good middle-layer family $`\mathcal F`$, let
$`K_{S,T}=t(S,T)\ge0`$ be the pairwise uncrossing deficit.  Because
this kernel is represented by Boolean features of degrees one and
two,
```math
\boxed{\operatorname{rank}K\le n+\binom n2.}
\tag{10.261}
```
If a subfamily $`\mathcal G`$ has every off-diagonal deficit in
$`\{1,\ldots,L\}`$, polynomial interpolation in $`K`$ gives
```math
|\mathcal G|
\le
\binom{n+\binom n2+L}{L}.
\tag{10.262}
```
Thus an incidence-threshold clique already forces
```math
\max_{S\ne T}t(S,T)
\ge
\left(\frac{\log2}{2}+o(1)\right)\frac n{\log n}.
\tag{10.263}
```
This is still far below the required $`n^{3/2}`$ payment, and
off-diagonal zeros evade the isolating polynomial completely.

There is a sharper pseudorandom identity.  If
```math
p_i=\Pr(i\in S),\quad p_{ij}=\Pr(i,j\in S),\quad
d_{ij}=p_i+p_j-2p_{ij},
```
then for any $`\delta>0`$,
```math
\boxed{
\mathbb Et(S,T)
\ge
\delta\,\mathbb E M(A[S])
-
\sum_{i<j}p_{ij}|d_{ij}-\delta|.
}
\tag{10.264}
```
Hence weighted-$`L^1`$ pseudorandom separation probabilities force
one pair with an $`\Omega(n^{3/2})`$ deficit.  The genuinely hard
branch must instead have strong $`A`$-correlated two-coordinate
bias, a dense exact-zero deficit graph, or both.  The strict-field
polynomial/zero-graph bridge remains under active investigation in
`deficit_code_inverse.md`.

### 10.46 A sharper spectral bootstrap and the strict-LCP fibre

For every symmetric zero-diagonal signing $`B`$, with
```math
Q(B)=\max_{x\in\{\pm1\}^m}|x^\mathsf TBx|,
```
there is an exact spectral inequality
```math
\boxed{
Q(B)\ge
\frac12\|B\|_{\rm op}\bigl(\|B\|_{\rm op}+1\bigr).
}
\tag{10.265}
```
To see this, take a unit eigenvector $`v`$ of eigenvalue magnitude
$`\lambda=\|B\|_{\rm op}`$, and put
$`\alpha=\|v\|_\infty`$.  At a coordinate attaining $`\alpha`$,
the zero diagonal gives
```math
\|v\|_1/\alpha\ge\lambda+1.
```
For $`z=v/\alpha`$,
```math
\|B\|_{\infty\to1}
\ge\|Bz\|_1
=\lambda\|z\|_1
\ge\lambda(\lambda+1).
```
Polarization and cube extremality for zero-diagonal quadratic forms
give
```math
\|B\|_{\infty\to1}\le2Q(B),
```
proving (10.265).  This proof was independently audited.

Every principal submatrix of a competitive signing has
$`Q(B)=O(n^{3/2})`$, hence
```math
\boxed{\|B\|_{\rm op}=O(n^{3/4}),}
\tag{10.266}
```
improving the earlier $`O(n^{5/6})`$ bootstrap.

This strengthens the code-separated inverse branch.  Fix a
positive-good support $`S`$, of child maximum $`p_S`$.  If no family
member $`T`$ has
```math
t(S,T)\ge p_S/4,
```
then every trace $`S\setminus(T\cap S)`$ lies in the child upper tail
$`H_{A[S]}\ge p_S/2`$.  Hanson--Wright, (10.266), and
$`p_S=\Omega(n^{3/2})`$ bound that tail by
```math
2^{|S|}\exp[-\Omega(n^{3/4})].
\tag{10.267}
```
An incidence-size family therefore contains an
$`\exp[\Omega(n^{3/4})]`$-sized fibre with one fixed intersection
$`X=T\cap S`$.

Writing its members as $`T=X\cup Z`$, $`Z\subseteq V=S^c`$, the
fibre solves one common affine strict complementarity problem
```math
\boxed{
(D1_Z+b)_v=0\ (v\notin Z),\qquad
(D1_Z+b)_v\ge1\ (v\in Z),
}
\tag{10.268}
```
where $`D=A[V]`$ and $`b=A_{V,X}1_X`$.  For two solutions
$`\eta,\theta`$, with active fields $`h^\eta,h^\theta`$,
```math
\boxed{
(\eta-\theta)^\mathsf TD(\eta-\theta)
=
\sum_{\eta\setminus\theta}h^\eta_i
+
\sum_{\theta\setminus\eta}h^\theta_i
\ge\|\eta-\theta\|_2^2.
}
\tag{10.269}
```
Thus every code secant lies in the positive cone of $`D-I`$.  If
$`\eta`$ is uniform on the fibre, with covariance $`\Sigma`$,
coordinate densities $`p_i`$, and conditional active margins
$`\mu_i`$, then
```math
\boxed{
\operatorname{tr}(D\Sigma)
=\sum_i p_i(1-p_i)\mu_i
\ge\operatorname{tr}\Sigma.
}
\tag{10.270}
```
The remaining dichotomy is now spectral versus fixed core:
coordinate entropy either consumes positive spectral mass, or
concentrates on an almost-fixed principal set.  The exponent
$`n^{3/4}`$ is critical for the present bounds, so rank/arrangement
arguments alone do not yet force a macroscopic descent.

Finally, zero pair deficit has exact face content.  If
$`X=S\cap T`$, $`Y=S\setminus T`$, then deficit $`t=t(S,T)`$
implies
```math
0\le p_X-H_{A[X]}(1)\le t,\qquad
0\le p_Y-H_{A[Y]}(1)\le t,
\tag{10.271}
```
and the cross block has value at most $`t`$ on every pair of
positive component grounds.  At $`t=0`$, both component all-one
spins are exact tops, the cross block annihilates the complete
Cartesian ground family, and all component-ground concatenations
are child grounds.  Dense zero-deficit structure is therefore a
compatible face-factorization problem, not merely a scalar
zero graph.  See `deficit_code_inverse.md`.

For an endpoint zero cut with blocks $`B,C,D`$, let
$`d_+(B),d_+(D)`$ be the dimensions of the two principal positive
ground spans.  Exact Cartesian factorization gives
```math
\boxed{
\operatorname{rank}C
\le
(|B|-d_+(B))+(|D|-d_+(D)),
\qquad
d_+(A)\ge d_+(B)+d_+(D).
}
\tag{10.271a}
```
Since $`\|C\|_F^2=|B||D|`$ and every competitive principal block has
operator norm $`O(n^{3/4})`$, a balanced zero cut costs
$`\Omega(\sqrt n)`$ combined ground-span codimension.  This is the
first quantitative resource unavailable to scalar uncrossing, but
it lands exactly at the mesoscopic $`\sqrt n`$ wall.

That wall is genuine for affine-LCP methods.  A paired checkerboard
signing on hidden coordinate pairs has every union of pairs as a
strict solution of (10.268), hence an exponentially large
fixed-weight solution layer, while a conference-like quotient keeps
operator norm $`O(\sqrt n)`$.  Counting solutions, hyperplane
arrangements, determinants, and operator norm alone therefore cannot
prove the inverse theorem.  The complete countermodel is excluded
only by the absolute child-ground condition: opposite spins inside
the selected pairs expose an energy larger in absolute value than
the stipulated all-one child ground.  Any successful inverse theorem
must use precisely this ground-energy information.

### 10.47 Global optimality, exact endpoint covers, and insertion

In a positive absolute-ground gauge, write
```math
c(S)=\sum_{i\in S,j\notin S}a_{ij},\qquad
q(S)=\min\{c(S),M-c(S)\}.
```
Then
```math
0\le c(S)\le M,\qquad
M-|H_A(1^S)|=2q(S),
```
and the full insertion hierarchy is exactly
```math
\boxed{
|b([n])-2b(S)|\le d+2q(S)
\quad(S\subseteq[n]).
}
\tag{10.272}
```

If $`A`$ is globally optimal at order $`n`$, delete vertex $`i`$,
call the child $`B_i`$, and let $`r_i`$ be the deleted row.  Then
```math
\boxed{
\min_b\Delta_{B_i}(b)
=M_n-M(B_i)=:d_i,
}
\tag{10.273}
```
and $`r_i`$ is an exact minimizer.  Thus below $`d_i`$, the
variable-radius near-cap constraints cover the entire row cube.
For any positive absolute ground with local-loss profile $`\ell`$,
```math
\sum_i d_i\le\sum_i\ell_i=2M_n.
\tag{10.274}
```
With opposite absolute grounds,
```math
\sum_i d_i
\le
2M_n-\frac12\|\ell^+-\ell^-\|_1.
\tag{10.275}
```

Switch the positive ground to $`1`$, and let the negative ground
define $`S\sqcup S^c`$, with blocks $`B,C,D`$.  Then
```math
1^\mathsf TC1=M_n,\qquad
\|C\|_{\infty\to1}=M_n,
\tag{10.276}
```
and
```math
\frac12\|\ell^+-\ell^-\|_1
=\|B1\|_1+\|D1\|_1.
\tag{10.277}
```
Hence failure of the scale-correct $`3/2`$ average forces a
near-bipartite/Eulerian regime:
```math
\|B1\|_1+\|D1\|_1<M_n/2.
\tag{10.278}
```

Every favorable positive cross edge has, by global edge-flip
optimality, a replacement witness within two of one absolute
endpoint.  On
```math
n\equiv1\pmod4,
```
all cut sums and $`M_n`$ have the parity that sharpens this
completely:
```math
\boxed{
\text{every such witness is an exact top or bottom ground}
\quad(c=0\ \text{or}\ c=M_n).
}
\tag{10.279}
```
This matters because convergence of the full normalized sequence is
equivalent to convergence on primes $`1\bmod4`$.

The exact cover does not require both absolute orientations to be
active.  After switching any positive ground to $`1`$, flipping any
positive coefficient lowers that endpoint by two.  Global minimality
and parity force an exact original ground which either lies on the
top face and separates the edge, or lies on the bottom face and does
not separate it.  If the bottom face is inactive, only the first
alternative can occur.  Since there are
```math
\frac{\binom n2+M_n}{2}=\Theta(n^2)
```
positive coefficients in this gauge, every prime-order minimizer has
a dense exact endpoint-cut cover, even in the one-sided case.

In the one-sided case this already yields a canonical quotient.  Put
```math
i\sim_+j
\iff
w_i=w_j\quad\hbox{for every positive exact ground }w.
```
Every positive ground is constant on each equivalence class, and
every diagonal class block is an all-negative clique: a positive
edge inside a class would have to be separated by an exact top
ground.  An insertion row can be balanced to sum to zero on every
even class and to $`\pm1`$ on every odd class.  Its field on the
whole active exact face is therefore at most the number of odd
classes.  Consequently, fewer than $`o(\sqrt n)`$ odd classes gives
the required endpoint-scale balance; failure leaves
$`\Omega(\sqrt n)`$ explicit negative-clique types.  Thick-cap
control of that quotient remains necessary.

For a general near-cap witness $`x`$, define its oriented deletion
profile $`u_i(x)`$.  Two witnesses obey
```math
\sum_i d_i
\le
\frac{\sum_i u_i(x)+\sum_i u_i(y)
-\|u(x)-u(y)\|_1}{2}.
\tag{10.280}
```
For the top ground and a top-near cut $`T`$,
```math
\boxed{
\sum_i d_i
\le
2M_n+(n-2)c(T)
-
\sum_i|\partial_T(i)-c(T)|.
}
\tag{10.281}
```
Thus boundary traffic at least $`M_n/2+O(n)`$ recovers
```math
\min_i d_i
\le\frac{3M_n}{2n}+O(1),
\tag{10.282}
```
whose $`O(1)`$ error is summable after normalization.  Otherwise
the $`\Theta(n^2)`$ favorable cross edges are covered by low-traffic
endpoint cuts; at prime orders these are exact zero/max cuts, so the
Cartesian factorization theorem applies with no approximation loss.

There remains a precise exchange obstruction.  Put
```math
e_i=M(B_i)-M_{n-1}.
```
Then
```math
\boxed{
e_i+d_i=M_n-M_{n-1}.
}
\tag{10.283}
```
For an optimal child $`C`$, the deficit-profile metric
```math
\eta(B_i,C)=\max_x|g_{B_i}(x)-g_C(x)|
```
satisfies
```math
M_n-M_{n-1}\le d_i+\eta(B_i,C),
\qquad
\eta(B_i,C)\ge e_i.
\tag{10.284}
```
So the scale-correct $`d_i`$ estimate becomes a recurrence only
after controlling the nonoptimal-child defect in this full profile
metric.  The active prime-order route is to turn the exact
low-traffic endpoint incidence cover into such stability or into a
bounded laminar quotient.  See `near_cap_insertion.md` and
`prime_face_cover_quotient.md`.

### 10.48 Canonical pressure reduced to a cut-code down-set

Define
```math
\mathfrak Z_n(\beta)
=\sum_Ae^{-\beta\sqrt n\,M(A)},\qquad
\Phi_n(\beta)=n^{-2}\log\mathfrak Z_n(\beta).
```
If $`\Phi_n(\beta)`$ converges for arbitrarily large fixed
$`\beta`$, then the entropy squeeze gives
```math
\limsup\frac{M_n}{n^{3/2}}
-
\liminf\frac{M_n}{n^{3/2}}
\le\frac{\log2}{2\beta},
```
and hence proves the desired limit.

The exact Gibbs--Shearer restriction inequality is
```math
\log\mathfrak Z_n(\beta)
\le q^{-1}
\log\mathfrak Z_m\!\left(\beta q\sqrt{n/m}\right),
\qquad
q=\frac{(m)_2}{(n)_2}.
\tag{10.285}
```
Its characteristic is
$`\beta\mapsto\beta(m/n)^{3/2}`$, so it does not close at fixed
temperature.  The missing fixed-temperature strengthening is
equivalent to a $`3/2`$-coefficient average-deletion theorem, and
conference matrices already falsify that statement without
additional global-minimality structure.

There is nevertheless an exact ground-gauge reduction.  Let
```math
\mathcal G_n
=\{B:T(B)=M(B)\},\qquad T(B)=\sum_{i<j}b_{ij}.
```
Then
```math
B\in\mathcal G_n
\iff
T(B)\ge0,\quad
0\le c_B(S)\le T(B)\quad(S\subseteq[n]).
\tag{10.286}
```
For
```math
\mathcal Y_n(\beta)
=\sum_{B\in\mathcal G_n}e^{-\beta\sqrt n\,T(B)},
```
double-counting augmented-cut ground pairs gives
```math
\boxed{
2\mathcal Y_n(\beta)
\le\mathfrak Z_n(\beta)
\le2^n\mathcal Y_n(\beta).
}
\tag{10.287}
```
Thus the two pressures have identical subsequential limits.

If $`F`$ is the negative-edge set, this cone is the down-set
```math
\mathcal I_n
=
\{F:2|F\cap D|\le|D|
\text{ for every augmented cut-code support }D\},
\tag{10.288}
```
and $`\mathcal Y_n`$ is its high-fugacity independence polynomial
after a linear prefactor.  Exact checks show that the hoped-for
standard mechanisms are absent: signed cuts are not submodular,
$`\mathcal I_n`$ is not a matroid already at $`n=5`$, and it is not
hereditary under $`K_5\to K_4`$ restriction.  Its coefficients are
log-concave through $`n=7`$, but no cross-order theorem follows.
Block gluing reintroduces the rectangular Gale--Berlekamp pressure
and has no reverse inequality.  This scalar pressure route is
therefore stopped; any revival must retain a near-ground overlap or
cap-intersection profile.  See `signing_free_energy_limit.md`.

### 10.49 Prime-order face quotients and the exact crossing wall

At orders $`n\equiv1\pmod4`$, the exact replacement law can be
organized by the affine endpoint faces.  Choose affine bases of the
top and bottom ground-correlation faces, use their difference cuts as
binary vertex coordinates, and partition the vertices into equal
signature atoms.

In the one-sided case $`\nu(A)<M_n`$, every signature atom is an
all-negative clique.  Averaging products of a spanning set of
top-face difference spins kills every inter-atom correlation and
retains every intra-atom correlation.  If the nonempty atoms have
sizes $`m_\sigma`$ and their number is $`r`$, then
```math
\boxed{
\nu(A)\ge\sum_\sigma\binom{m_\sigma}{2}
\ge\frac12\left(\frac{n^2}{r}-n\right).
}
\tag{10.289}
```
Since $`r\le2^{d_+}`$, every one-sided
$`M(A)=O(n^{3/2})`$ competitor satisfies
```math
\boxed{
r=\Omega(\sqrt n),\qquad
d_+\ge\frac12\log_2n-O(1).
}
\tag{10.290}
```
Thus absence of the opposite absolute face forces a genuinely
growing exact top face.

In the balanced two-sided branch, let $`S`$ be the shore separating a
top ground from a bottom ground.  If a simultaneous signature atom
contains $`a_\sigma`$ vertices of $`S`$ and $`b_\sigma`$ vertices of
$`S^c`$, then its mixed block is all negative, whence
```math
a_\sigma b_\sigma\le M_n,\qquad
\sum_\sigma\min(a_\sigma,b_\sigma)
\le2^{d_++d_-}\sqrt{M_n}.
\tag{10.291}
```
Consequently, if
```math
d_++d_-\le(1/4-\varepsilon)\log_2n,
```
then $`S`$ is a union of signature atoms after deleting $`o(n)`$
vertices.  This is an exact finite quotient, but Boolean combinations
of its basis cuts need not themselves be ground differences.

Two crossing top zero cuts have a complete local disintegration.
Write their four regions as
```math
X=R\cap Q,\quad Y=R\setminus Q,\quad
Z=Q\setminus R,\quad T=[n]\setminus(R\cup Q),
```
and let $`w_{UV}`$ be the signed sum of the block $`U\times V`$.
For $`k=\frac12c(R\triangle Q)`$ one has
```math
\boxed{
w_{XY}+w_{ZT}=w_{XZ}+w_{YT}=k
}
\tag{10.292}
```
and
```math
\boxed{
\begin{aligned}
2w_{XY}&=c(X)+c(Y),&
2w_{ZT}&=c(Z)+c(T),\\
2w_{XZ}&=c(X)+c(Z),&
2w_{YT}&=c(Y)+c(T).
\end{aligned}}
\tag{10.293}
```
Thus $`k=0`$ forces all six inter-region block sums and all four
corner cuts to vanish: the two cuts refine to an exact four-way
Cartesian factorization.  For $`k>0`$, the two child-ground deficits
and the capped-bilinear payment are both exactly $`2k`$, while the
product ground has deficit $`4k`$.

This local identity is sharp but does not supply a size-dependent
payment.  The rank-one signing
```math
a_{ij}=-t_it_j,\qquad t_i\in\{\pm1\},\quad\sum_it_i=1,
```
has $`c(U)=T_U(T_U-1)`$, where $`T_U=\sum_{i\in U}t_i`$, and admits
four linear-sized regions supporting irreducibly crossing zero cuts
with constant $`k`$.  This family has norm $`\Theta(n^2)`$, so the
still-unused hypothesis is precisely the global $`O(n^{3/2})`$ cap.

Free laminarization is false even under global minimality.  An exact
order-nine minimizer with $`M_9=12`$ has nineteen nontrivial top- and
bottom-face difference cuts covering all twenty-four positive
coefficients, but exhaustive enumeration of all $`2^{19}`$
subfamilies shows that its 792 laminar subfamilies cover at most
twenty-three of them.  A viable continuation must therefore extract
mesoscopic energy from many crossing pairs using the global cap, not
replace the exact cover by a laminar one formally.  See
`prime_face_cover_quotient.md`.

### 10.50 Weight-two concatenation and association-scheme no-go

Every signing is the weight-two restriction of a quadratic Boolean
function.  Put
```math
b_{ij}=(1-a_{ij})/2\in\mathbb F_2,\qquad
q_B(z)=\sum_{i<j}b_{ij}z_iz_j.
```
Then
```math
(-1)^{q_B(e_i+e_j)}=a_{ij}.
\tag{10.294}
```
The ambient quadratic function is bent exactly when the alternating
matrix $`B`$ is nonsingular over $`\mathbb F_2`$.  This condition is
provably blind to the present norm: a perfect-matching negative set
is nonsingular and has $`M(A)=\Theta(n^2)`$, while random
nonsingular alternating matrices include signings with
$`M(A)=O(n^{3/2})`$.  The standard Schmidt direct-sum recursion
$`q\mapsto q+uv`$ restricts exactly to adding one negative matched
edge and making all old--new edges positive, and hence iterates into
the bad quadratic-scale family.

Tensor products do not remain in the problem: a product of two
weight-two characters has Fourier weight four and has zero
level-two projection.  Restricting the product to rank-one spin
matrices tests only the Segre subset, whereas $`K_{nm}`$ tests every
Boolean spin matrix; the missing configurations are exactly the
multichannel bilinear relaxation already known to lose the scalar
geometry.

The obstruction also holds at every degree of the ordinary Johnson
scheme.  Positive supports of any two $`d`$-regular graphs have the
same complete $`J(n,2)`$ inner distribution, since the numbers of
edge pairs sharing a vertex and disjoint edge pairs depend only on
$`(n,d)`$.  Yet two disjoint cliques give a rank-one signing with
$`M=\Theta(n^2)`$, while dense regular expanders give
$`M=O(n^{3/2})`$.  Hence the entire one-point
Johnson/Delsarte algebra, not merely bounded moments, is blind.

Finally, an exact triangle-aware pointwise algebra cannot remain
small.  The edge characters generate every even Boolean character,
so any unital pointwise algebra containing them has dimension
```math
2^{n-1}
\tag{10.295}
```
and is the full cut-character algebra.  Closure under products of at
most $`t`$ edge characters already has dimension at least
```math
\sum_{j=0}^t\binom n{2j}.
\tag{10.296}
```
Thus there is no polynomial-dimensional exact harmonic projection
which preserves the triangle identities and the relevant sup norm.
No scale-preserving bent/Reed--Muller/Johnson concatenation survives
this audit.  See `spherical_code_concatenation.md`.

### 10.51 Optimal punctures fail; mean puncture stability is the useful target

For an optimal parent $`A`$ and its principal children
$`B_i=A[-i]`$, recall
```math
d_i=M_n-M(B_i),\qquad
e_i=M(B_i)-M_{n-1},
\qquad
e_i+d_i=M_n-M_{n-1}.
\tag{10.297}
```
The strongest possible puncture statement is false at the first
relevant arithmetic order.  The explicit order-nine minimizer in
`optimal_child_stability.md` satisfies
```math
\boxed{
M(A)=12=M_9,\qquad
M(A[-i])=12\quad(i=1,\ldots,9),
}
\tag{10.298}
```
while $`M_8=10`$.  Hence
```math
\boxed{d_i=0,\qquad e_i=2\quad\text{for every }i.}
\tag{10.299}
```
This independently certified example shows that a globally optimal
order $`4k+1`$ signing need not have any optimal child, even when
every deletion has zero parent-cap cost.

In covering-radius notation,
```math
\rho_n=\frac{\binom n2-M_n}{2},\qquad
r_i=\frac{\binom{n-1}2-M(A[-i])}{2},
```
the defect is
```math
\boxed{e_i=2(\rho_{n-1}-r_i).}
\tag{10.300}
```
The exact surviving finite target is therefore
```math
\boxed{
\max_{A:\operatorname{dist}(A,\mathcal C_n)=\rho_n}
\min_i(\rho_{n-1}-r_i)\le1,
}
\tag{10.301}
```
possibly after choosing a compatible deepest hole.  Equivalently,
some child should have $`e_i\le2`$; the order-nine example shows the
constant would be sharp.  The elementary extension inequality
```math
\rho_n\le r_i+\left\lfloor\frac{n-1}{2}\right\rfloor
\tag{10.302}
```
is far too weak, so any proof must use compatibility among different
puncture approximants, exact endpoint optimality, or the triangle
relations of the augmented cut code.

There is a second correction: radius-one puncturing alone does not
close the convergence recurrence.  If
```math
\delta_n=M_n-M_{n-1},
```
then
```math
\boxed{e_i=\delta_n-d_i.}
\tag{10.303}
```
Thus the child with smallest defect is exactly the child with largest
deletion cost, while the near-cap argument selects a child with small
deletion cost.  Separate bounds on $`\min e_i`$ and $`\min d_i`$ need
not concern the same vertex.

The clean sufficient average statement is
```math
\boxed{
\overline e_n=\frac1n\sum_i e_i=O(1)
}
\tag{10.304}
```
when the traffic argument supplies
```math
\overline d_n\le\frac{3M_n}{2n}+O(1).
```
Indeed $`\delta_n=\overline d_n+\overline e_n`$, yielding the
scale-correct recurrence with a summable normalized error.  More
generally it is enough that
```math
\sum_n\frac{\overline e_n}{n^{3/2}}<\infty.
\tag{10.305}
```
For a specifically selected small-$`d_i`$ child, the corresponding
joint requirement is
```math
\sum_n\frac{e_i}{n^{3/2}}<\infty.
\tag{10.306}
```

The compatibility obstruction among closest child codewords can be
made exact.  For a deepest hole $`A`$, choose a closest
$`c_i\in\mathcal C_{n-1}`$ to each puncture and define
```math
\mathcal R
=
\min_{c\in\mathcal C_n}
\sum_i d(c[-i],c_i).
```
Counting every global edge in exactly $`n-2`$ punctures gives
```math
\boxed{
(n-2)M_n
\ge
nM_{n-1}+\sum_i e_i-2\mathcal R.
}
\tag{10.307}
```
If
```math
\Gamma
=
\sum_{i<j}
d\!\left(
c_i|_{E([n]\setminus\{i,j\})},
c_j|_{E([n]\setminus\{i,j\})}
\right),
```
a majority decoder followed by the augmented-cut triangle test gives,
for $`n\ge12`$,
```math
\boxed{\mathcal R\le100\,\Gamma/n.}
\tag{10.308}
```
Hence a scale-smooth minimizing sequence whose child defects do not
collapse must have
```math
\Gamma=\Omega(n^{5/2}).
\tag{10.309}
```
Ordinary local agreement is therefore intrinsically too coarse: the
difference between constant defects two and four is only an
$`O(n)`$ signal on top of a leading $`n^{3/2}`$ repair budget.  The
active continuation is to exploit the fact that these disagreements
are themselves cuts or complemented cuts and combine their massive
incompatibility with endpoint optimality.

Finally, parity alone cannot supply a universal completion theorem.
An explicit order-twelve signing $`B`$ has
```math
M(B)=18,\qquad
\min_b\max_x\bigl(|H_B(x)|+|b\cdot x|\bigr)=24.
\tag{10.310}
```
Thus even at orders divisible by four the fixed-core affine insertion
gap can be six.  The displayed matrix and an exhaustive deterministic
verifier are in `optimal_child_stability.md` and
`optimal_child_stability_verify.cpp`.

### 10.52 A linear zero-cut rank theorem and its sharp reuse obstruction

The previous $`\Omega(\sqrt n)`$ rank wall for balanced exact zero
cuts can be upgraded to a linear wall without any operator-norm
assumption.  For every flat $`r\times s`$ cross block $`C`$,
```math
\|C\|_{\infty\to1}\le W(A),
\tag{10.311}
```
because changing the relative global sign of the two shores produces
two full energies separated by $`2u^\mathsf TCv`$.  If
$`\operatorname{rank}C=k`$, $`\gamma_2`$-duality and real
Grothendieck give
```math
rs=\langle C,C\rangle
\le\gamma_2(C)\gamma_2^*(C)
\le\sqrt{k}\,K_G\|C\|_{\infty\to1}.
```
Consequently
```math
\boxed{
\operatorname{rank}A[U,U^c]
\ge
\left(
\frac{|U|(n-|U|)}{K_GW(A)}
\right)^2.
}
\tag{10.312}
```
For $`W(A)=O(n^{3/2})`$, every balanced vertex split therefore has
linear cross rank.

If $`U`$ is an exact endpoint zero cut, Cartesian ground
factorization converts (10.312) into
```math
\boxed{
\Phi(A[U])+\Phi(A[U^c])
\ge
\left(
\frac{|U|(n-|U|)}{K_GW(A)}
\right)^2,
}
\tag{10.313}
```
where
```math
\Phi(B)=|B|-\dim\operatorname{span}\mathcal G_+(B).
```
Thus every balanced exact zero cut consumes a linear combined
component-ground-span codimension.

Three top grounds give the three Klein-four bipartitions of their
four Venn cells.  If the cell sizes are $`n_\omega`$, the sum
$`\mathcal P`$ of their three component-codimension payments obeys
```math
\boxed{
\mathcal P
\ge
\frac{
\left(n^2-\sum_\omega n_\omega^2\right)^2
}{3K_G^2W(A)^2}.
}
\tag{10.314}
```
In particular, four equal cells at the spectral-width scale
$`W\le(1/2+o(1))n^{3/2}`$ cost at least
```math
\left(\frac{3}{4K_G^2}-o(1)\right)n.
```

If zero-payment uncrossing closes to a Boolean algebra with atoms
$`P_i`$, put
```math
\Psi=\sum_i\left(|P_i|-d_+(A[P_i])\right).
```
Every atom union $`U`$ then satisfies
```math
\boxed{
|U|(n-|U|)\le K_GW(A)\sqrt{\Psi}.
}
\tag{10.315}
```
Hence
```math
W(A)=O(n^{3/2}),\quad\Psi=o(n)
\quad\Longrightarrow\quad
\text{some atom has size }n-o(n).
\tag{10.316}
```

The linear budget cannot simply be summed.  For $`n=k^2`$, with
$`k`$ an even Hadamard order, partition the vertices into $`k`$
classes of size $`k`$, put $`+1`$ inside each class, and use
orthogonal rank-one Hadamard residual channels between each ordered
pair of classes.  This gives
```math
\boxed{
W(A)=\frac12n^{3/2},\qquad
M(A)=\frac12n^{3/2}+\frac n2,\qquad
\Psi=n-\sqrt n.
}
\tag{10.317}
```
Every union of the $`\sqrt n`$ classes is a zero cut, and a union of
$`a`$ classes has exact cross rank
```math
a(\sqrt n-a).
\tag{10.318}
```
Thus exponentially many zero cuts can reuse the same linear
codimension budget even at optimal spectral scale.  The potential
$`\Phi`$ is also neither submodular nor supermodular, so ordinary
Shearer summation cannot repair the loss.

For the puncture tangent family, a balanced pair $`i,j`$ has a
zero-triple closure set
```math
\mathcal K_{ij}
=
\{k:x_ix_jx_k\text{ is another top ground}\}
```
whose real span satisfies
```math
\operatorname{rank}\{x_k:k\in\mathcal K_{ij}\}
\le
n-\Phi(A[D_{ij}])-\Phi(A[D_{ij}^c]).
\tag{10.319}
```
A full-rank tangent family therefore forces linearly many nonclosed
cubic products per balanced pair.  Parity charges each such product
only a constant deficit, however, so averaging still loses the
leading scale.

The generic rank/uncrossing route is now stopped at a sharp infinite
model.  A continuation must use global-minimizer exchange, prove a
strict upper gap below centered constant $`1/2`$, or establish stable
near-ground frame bounds.  See `mesoscopic_zero_cut_rigidity.md`.

### 10.53 Endpoint slack, stable frames, and replenishment

This wave tested three independent continuations of the frontier in
Sections 10.47--10.52.  The endpoint-slack and frame calculations below
are verified analytically.  The displayed finite minima for
compatibility disagreement are numerical and are not used in any
proof.

#### 10.53.1 Extension slack controls signed cut weight, not agreement

Represent an augmented-cut word by
```math
q^{t,x}_{uv}=t x_ux_v,\qquad t\in\{\pm1\},
```
so that $`\langle A,q^{t,x}\rangle=tH_A(x)`$.  For a globally
optimal parent $`A`$, choose a closest word $`c_i=q^{t_i,x^{(i)}}`$
to each child $`B_i=A[-i]`$, and put
```math
h_i=t_i\sum_{j\ne i}a_{ij}x_j^{(i)},\qquad
d_i=M_n-M(B_i),\qquad
s_i=d_i-|h_i|.
```
Extending $`c_i`$ with the better deleted spin gives a full word
$`q_i`$.  Parent optimality implies $`|h_i|\le d_i`$, and the exact
extension identity is
```math
\boxed{
\langle A,q_i\rangle=M_n-s_i,\qquad s_i\ge0.
}
\tag{10.320}
```

Switch $`A`$ by $`q_j`$, write $`S_{ij}`$ for the difference shore
between $`q_i`$ and $`q_j`$, and let $`C_j(S)`$ be its signed cut
weight in that gauge.  Direct subtraction in (10.320) gives
```math
\boxed{
C_j(S_{ij})=
\begin{cases}
(s_i-s_j)/2,&t_i=t_j,\\[2mm]
M_n-(s_i+s_j)/2,&t_i=-t_j.
\end{cases}
}
\tag{10.321}
```
Relative to a minimum-slack reference, same-orientation differences
are small nonnegative signed cuts and opposite-orientation differences
are near-maximal signed cuts.

This identity does not control their cardinalities.  On the common
$`K_{n-2}`$, put $`k_{ij}=|S_{ij}\setminus\{i,j\}|`$ and $`r=n-2`$.
The pairwise contribution to the compatibility defect is exactly
```math
\boxed{
\Gamma_{ij}=
\begin{cases}
k_{ij}(r-k_{ij}),&t_i=t_j,\\[2mm]
\binom r2-k_{ij}(r-k_{ij}),&t_i=-t_j.
\end{cases}
}
\tag{10.322}
```
Thus $`s_i`$ controls a signed $`A`$-weighted cut sum, while
$`\Gamma_{ij}`$ counts unweighted overlap edges.  Cancellation can
make the former zero while the latter is large.

The order-nine minimizer from Section 10.51 makes this obstruction
analytic.  Every closest child word has
```math
d_i=h_i=s_i=0.
```
Nevertheless $`\Gamma`$ cannot vanish.  Otherwise the local words
would glue to a global augmented-cut word $`c`$.  Each child distance
is eight, so edge counting would give
```math
72=\sum_i d(A[-i],c[-i])=7d(A,c),
```
whereas $`d(A,c)\ge\rho_9=12`$ makes the right side at least $`84`$.
In fact a nonzero ordinary cut on the common $`K_7`$ has at least six
edges, and a complemented cut has at least nine, so
```math
\boxed{\Gamma\ge6.}
\tag{10.323}
```
Consequently every universal estimate
$`\Gamma\le F(s_1,\ldots,s_n)`$ with $`F(0,\ldots,0)=0`$ is
**falsified**.  Any surviving compatibility theorem must control the
cardinalities or crossing patterns of the exact zero/max cuts, not
only their extension slacks.

**Numerical.**  Two independent exact-integer exhaustive evaluators
found child-choice counts
```math
(2,4,2,4,2,2,4,4,4),
```
hence $`16{,}384`$ closest-word families, with
```math
\min\Gamma=228,\qquad \min\mathcal R=58.
```
The two minima need not occur for the same family, and neither value
is needed for (10.323).

#### 10.53.2 A stable endpoint-frame inequality and its scale wall

For
```math
A=\begin{pmatrix}B&C\\C^{\mathsf T}&D\end{pmatrix},
```
define
```math
e_+=p(A)-p(B)-p(D),\qquad
e_-=\nu(A)-\nu(B)-\nu(D).
```
Both endpoint excesses are nonnegative.  For $`\sigma\in\{+,-\}`$,
let $`\mu_{B,\sigma}`$ and $`\mu_{D,\sigma}`$ be independent laws on
arbitrary Boolean endpoint clouds.  Write
```math
R_{E,\sigma}=\mathbb E(xx^{\mathsf T})
\succeq\kappa_{E,\sigma}P_{E,\sigma},
```
where $`P_{E,\sigma}`$ projects onto the span of the law's support.
For the endpoint deficits
```math
d_E^+(x)=p(E)-H_E(x),\qquad
d_E^-(x)=\nu(E)+H_E(x),
```
put
```math
\eta_\sigma^2=
\mathbb E\bigl(e_\sigma+d_B^\sigma(u)+d_D^\sigma(v)\bigr)^2.
```

Changing the relative sign of $`u`$ and $`v`$ gives, with no hidden
factor of two,
```math
|u^{\mathsf T}Cv|
\le e_\sigma+d_B^\sigma(u)+d_D^\sigma(v).
```
Independence and the frame inequalities therefore imply
```math
\boxed{
\eta_\sigma^2
\ge
\kappa_{B,\sigma}\kappa_{D,\sigma}
\|P_{B,\sigma}CP_{D,\sigma}\|_F^2.
}
\tag{10.324}
```
If $`s_j(C)`$ are the singular values and
```math
a_\sigma=\operatorname{codim}L_{B,\sigma}
+\operatorname{codim}L_{D,\sigma},
```
then $`\operatorname{rank}(C-P_BCP_D)\le a_\sigma`$, so singular-value
interlacing sharpens (10.324) to
```math
\boxed{
\eta_\sigma^2
\ge
\kappa_{B,\sigma}\kappa_{D,\sigma}
\sum_{j>a_\sigma}s_j(C)^2
\ge
\kappa_{B,\sigma}\kappa_{D,\sigma}
[rs-\|C\|_{\mathrm{op}}^2a_\sigma]_+.
}
\tag{10.325}
```

For laws on exact endpoint grounds all deficits vanish.  Since
```math
W(A)-W(B)-W(D)=\frac{e_++e_-}{2},
```
one obtains the verified conditional superadditivity bound
```math
\boxed{
W(A)-W(B)-W(D)
\ge
\frac12\sum_{\sigma=\pm}
\sqrt{\kappa_{B,\sigma}\kappa_{D,\sigma}}
\sqrt{[rs-\|C\|_{\mathrm{op}}^2a_\sigma]_+}.
}
\tag{10.326}
```
For laws supported on $`\tau_B`$- and $`\tau_D`$-near clouds, replace
$`e_\sigma^2`$ by $`(e_\sigma+\tau_B+\tau_D)^2`$.

There is a useful near-zero-cut form.  Switch a positive full ground
to $`1`$ and let $`\alpha=1^{\mathsf T}C1`$.  The flipped full state
has endpoint deficit $`2\alpha`$, while $`0\le e_+\le\alpha`$.  Hence
```math
\boxed{
(\alpha+\tau_B+\tau_D)^2
\ge
\kappa_{B,+}\kappa_{D,+}
\sum_{j>a_+}s_j(C)^2.
}
\tag{10.327}
```

This theorem exposes a sharp scale limitation.  On a balanced split
with $`\|C\|_{\mathrm{op}}=O(\sqrt n)`$, frame product bounded below,
and $`a=o(n)`$, it forces only
```math
e_\sigma+\tau_B+\tau_D=\Omega(n),
```
not an $`n^{3/2}`$ payment.  This is intrinsic to second moments:
$`\operatorname{tr}R_E=|E|`$ forces $`\kappa_E=O(1)`$ on a
linear-dimensional span.

**Numerical finite audit.**  The globally optimal order-nine example
has a $`4+5`$ endpoint zero cut with cross rank four, but both child
positive-ground spans are only antipodal lines, of dimensions one and
one.  Thus global minimality does not itself force a useful exact-ground
frame.  The viable **open target** is a thick near-ground law with
```math
\tau_B+\tau_D=o(n),\qquad
\kappa_B\kappa_D\ge\kappa_0>0,\qquad
\sum_{j>a}s_j(C)^2\ge cn^2
```
for every balanced $`o(n)`$-traffic cut.  Equation (10.327) would then
give an immediate contradiction.

#### 10.53.3 Field-proportional peeling isolates replenishment

Use doubled normalization
```math
Q(A)=\max_x|x^{\mathsf T}Ax|=2M(A).
```
Switch a positive absolute ground to $`1`$, so
```math
q=1^{\mathsf T}A1=Q(A),\qquad r=A1\ge0.
```
For a deleted set $`H`$ with complement $`T`$, define
```math
\begin{aligned}
R_H&=\sum_{i\in H}r_i,&
h_H&=1_H^{\mathsf T}A[H]1_H,&
b_H&=1_H^{\mathsf T}A[H,T]1_T,\\
e_H&=1_T^{\mathsf T}A[T]1_T,&
d_H&=q-Q(A[T]),&
g_H&=Q(A[T])-e_H.
\end{aligned}
```
Random extension of a core ground proves $`d_H\ge0`$, and trivially
$`g_H\ge0`$.  Expanding $`q`$ and $`R_H`$ gives the exact identities
```math
\boxed{
2R_H=d_H+h_H+g_H,\qquad
d_H+g_H=R_H+b_H.
}
\tag{10.328}
```
Flipping $`H`$ in the parent ground gives $`0\le b_H\le q/2`$.
Thus positive heavy-field mass must appear as decrement or
replenishment, but (10.328) does not decide which.

Let $`S`$ delete vertices independently with probabilities $`p_i`$.
Averaging (10.328) yields
```math
\boxed{
2\sum_i p_ir_i
=\mathbb E d_S+p^{\mathsf T}Ap+\mathbb E g_S,\qquad
\mathbb E(d_S+g_S)
=\sum_i p_ir_i+\mathbb E b_S.
}
\tag{10.329}
```
Here $`0\le\mathbb E b_S\le q/2`$ and $`|p^{\mathsf T}Ap|\le q`$,
the latter by independent Boolean rounding with mean $`p`$.

For a fixed heavy set $`H`$, take
```math
p_i=\frac{r_i}{n-1}1_{\{i\in H\}}.
```
Then the row-square tail appears without loss:
```math
\boxed{
\mathbb E(d_S+g_S)
=\frac1{n-1}\sum_{i\in H}r_i^2+\mathbb E b_S
\ge\frac1{n-1}\sum_{i\in H}r_i^2.
}
\tag{10.330}
```
If $`A`$ is an exact order-$`n`$ minimizer and $`q_m=2M_m`$, then
```math
d_S=q_n-Q(A[S^c])\le q_n-q_{n-|S|}.
```
Consequently the exact minimality dichotomy is
```math
\boxed{
\frac1{n-1}\sum_{i\in H}r_i^2
\le
\mathbb E\bigl(q_n-q_{n-|S|}\bigr)+\mathbb E g_S.
}
\tag{10.331}
```
The sole uncontrolled positive term is successor replenishment.

Global signing minimality does not remove it pointwise.  For the
order-nine minimizer in Section 10.51, the absolute endpoint
```math
x=(1,-1,1,-1,-1,-1,-1,-1,1),\qquad x^{\mathsf T}Ax=-24,
```
has, after negative orientation and switching, row profile
```math
(6,4,8,2,0,2,2,0,0).
```
Every principal child has doubled norm $`24`$.  Deleting vertex three
therefore gives
```math
\boxed{
r_3=8,\qquad d_{\{3\}}=h_{\{3\}}=0,\qquad
e_{\{3\}}=8,\qquad g_{\{3\}}=16.
}
\tag{10.332}
```
Moreover the deleted star annihilates every absolute child ground:
otherwise the better extension would have parent norm greater than
$`24`$.  Hence exact successor-ground visibility is zero while the
replenishment gap is sixteen.

This **falsifies** every universal pointwise estimate such as
$`g_H\le C(d_H+|h_H|)`$, $`g_H\le CV_H`$, or a positive decrement
lower bound depending only on one heavy field, even for an exact
global minimizer.  It does not falsify a grouped asymptotic theorem.
The remaining **open target** is an excess-sensitive cumulative bound
on $`\mathbb E\sum_t g_t`$ for field-proportional deletions, using
successor near-ground layers up to deficit $`g_t`$ rather than only
exact successor grounds.

#### 10.53.4 Disposition of the three routes

- The slack-only compatibility route is stopped by (10.323).  Its
  successor is a prime-order crossing-energy theorem for the exact
  zero/max cut family forced by large $`\Gamma`$.

- Exact-ground frames and bare second moments are stopped at their
  scale wall.  Their successor is a common thick near-ground law with
  quantitative frame bounds at $`o(n)`$ deficit.

- Pointwise heavy-field peeling is stopped by (10.332).  Its successor
  is grouped multiscale replenishment charged to thick successor
  near-ground layers.

These are the three highest-priority independent routes for the next
wave.  None of the present results proves convergence by itself.

### 10.54 Crossing closure, SDP curvature, and grouped replenishment

The second wave tested the three successors in Section 10.53.4.  It
produced two general lemmas and three sharp walls.  All statements in
this section are verified analytically; none relies on numerical
optimization.

#### 10.54.1 Triple endpoint defect and an exact residue-class wall

Let $`q_1,\ldots,q_m`$ be oriented augmented-cut endpoint words for a
signing $`A`$, so
```math
\langle A,q_i\rangle=M(A).
```
For every edge $`e`$, put $`s_e=\sum_iq_i(e)`$.  Since augmented-cut
words form a group under entrywise multiplication, every
$`q_iq_jq_k`$ is again a valid augmented-cut word.  Direct expansion
gives the exact cubic identity
```math
\boxed{
\sum_{i,j,k}
\bigl(M(A)-\langle A,q_iq_jq_k\rangle\bigr)
=
\sum_e a_es_e(m^2-s_e^2).
}
\tag{10.333}
```
Every summand on the left is nonnegative.  The corresponding full
pairwise Hamming mass is
```math
\boxed{
\sum_{i<j}d(q_i,q_j)
=\frac14\sum_e(m^2-s_e^2).
}
\tag{10.334}
```
Thus (10.333) is the exact candidate for converting endpoint
disagreement into crossing payment.  It can vanish despite arbitrarily
large (10.334): this happens when the endpoint family is closed under
ternary products.

There is an exact flat signing showing that this wall survives at
orders $`1\bmod4`$ and at the competitive scale.  Let $`k`$ be an even
Hadamard order, $`N=k^2`$, and let $`A_0`$ be the square residual
signing from Section 10.52.  It has spectral values $`k-1`$ and
$`-k-1`$.  Add one universally positive vertex:
```math
A=
\begin{pmatrix}
A_0&1\\
1^{\mathsf T}&0
\end{pmatrix},
\qquad n=N+1\equiv1\pmod4.
```
If $`s_a`$ are the spin sums on the $`k`$ old vertex classes and
$`S=\sum_as_a`$, projection onto the class-constant positive
eigenspace gives
```math
x^{\mathsf T}A_0x
\ge-N(k+1)+2\sum_as_a^2.
```
Every $`s_a`$ is even, hence $`\sum_as_a^2\ge|S|`$.  Together with the
spectral upper bound this yields
```math
-N(k+1)
\le x^{\mathsf T}A_0x+2\varepsilon S
\le N(k+1).
```
A fixed perfect matching of the classes supplies Boolean vectors
attaining the lower endpoint, so
```math
\boxed{
Q(A)=N(k+1),\qquad
M(A)=\frac{N(k+1)}2
=\left(\frac12+o(1)\right)n^{3/2}.
}
\tag{10.335}
```

More explicitly, for every matched pair $`\{a,b\}`$, take
```math
x|_{V_a}=\eta_{ab}v_{a,b},\qquad
x|_{V_b}=-\eta_{ab}v_{b,a},
```
and fix the new spin.  The $`k/2`$ independent signs $`\eta_{ab}`$
give an affine family of $`2^{k/2}`$ exact negative absolute grounds.
It is ternary-closed.  For all large $`k`$, a standard random-code
selection gives $`n`$ members whose pairwise sign-code distances lie
between $`k/8`$ and $`3k/8`$.  Their difference shores therefore have
between $`N/4`$ and $`3N/4`$ old vertices.  Even after omitting the two
puncture stars, each pair contributes $`\Theta(n^2)`$, so an arbitrary
indexing of these endpoints has
```math
\Gamma=\Omega(n^4),
```
while every triple deficit in (10.333) is zero.

This **falsifies** any crossing-energy aggregation based only on
flatness, $`n\equiv1\bmod4`$, exact absolute endpoint cuts, and an
$`O(n^{3/2})`$ cap.  The example is not known to be an exact
order-$`n`$ minimizer, need not occur at prime $`n`$, and its indexed
endpoint restrictions are not proved to be closest child words.
Those global-minimizer, prime-order, and puncture-compatibility inputs
are exactly what a surviving theorem must use.

#### 10.54.2 The SDP dual slack controls every near-ground frame

First, generic isotropic thickening is too expensive.  Let $`G`$ have
any law on exact positive grounds of an $`m`$-vertex signing $`B`$,
and independently let exchangeable Boolean noise $`Y`$ satisfy
```math
\mathbb E(Y_iY_j)=\rho\qquad(i\ne j).
```
For $`X=G\circ Y`$,
```math
\boxed{
R_X=(1-\rho)I+\rho R_G,\qquad
\mathbb E[p(B)-H_B(X)]=(1-\rho)p(B).
}
\tag{10.336}
```
When $`0\le\rho\le1`$, an isotropic frame floor $`1-\rho`$ therefore
costs exactly $`(1-\rho)p(B)`$.  At
$`p(B)=\Theta(m^{3/2})`$, deficit $`o(m)`$ buys only
$`o(m^{-1/2})`$ isotropic floor.  Sampling many exact grounds before
adding the noise does not change this identity.

There is a general anisotropic replacement.  Define
```math
\operatorname{SDP}_+(B)
=
\frac12\max\{
\operatorname{tr}(BX):X\succeq0,\ \operatorname{diag}X=1
\},
\qquad
s_B=\operatorname{SDP}_+(B)-p(B).
```
Let
```math
L_B=\operatorname{Diag}(y)-B\succeq0,\qquad
\frac12\sum_i y_i=\operatorname{SDP}_+(B)
```
be an optimal dual slack.  Every Boolean vector obeys the exact
identity
```math
\boxed{
\frac12x^{\mathsf T}L_Bx
=s_B+p(B)-H_B(x).
}
\tag{10.337}
```
Hence any law of mean positive-endpoint deficit at most $`\tau`$, with
covariance $`R`$, satisfies
```math
\operatorname{tr}(L_BR)\le2(s_B+\tau).
```
If $`E_\gamma`$ is the spectral subspace of $`L_B`$ on
$`[\gamma,\infty)`$, and $`R\succeq\kappa P_L`$, then
```math
\boxed{
\operatorname{tr}(P_{E_\gamma}R)
\le\frac{2(s_B+\tau)}{\gamma},
\qquad
\kappa[\dim E_\gamma-\operatorname{codim}L]_+
\le\frac{2(s_B+\tau)}{\gamma}.
}
\tag{10.338}
```
Thus a nondegenerate near-ground frame must lie almost entirely in the
approximate kernel of the child SDP dual slack.

The square Hadamard residual model makes this obstruction sharp for
every law, not only isotropic ones.  Split its $`k`$ types into
$`a+b=k`$, with children $`B,D`$ and cross block $`C`$.  The cross
block has singular value $`k`$ on $`ab`$ orthogonal channels.  The
child SDPs are exact, and their dual slacks have eigenvalue $`k`$ on
the corresponding left and right channel spaces $`E_B,E_D`$.
Therefore laws of mean deficits $`\tau_B,\tau_D`$ obey
```math
\operatorname{tr}(P_{E_B}R_B)\le\frac{2\tau_B}{k},
\qquad
\operatorname{tr}(P_{E_D}R_D)\le\frac{2\tau_D}{k},
```
and the partial-isometry form $`C=kU:E_D\to E_B`$ gives
```math
\boxed{
\operatorname{tr}(R_BCR_DC^{\mathsf T})
\le4\tau_B\tau_D.
}
\tag{10.339}
```
Consequently $`\tau_B+\tau_D=o(n)`$ forces the cross second moment to
be $`o(n^2)`$ for every possible pair of laws.  Competitiveness,
spectral flatness, exact zero traffic, and thick sublinear caps do not
force the frame hypothesis from (10.327).

The construction is not known to be globally minimizing.  The new
**open target** is a block-exchange theorem: in an exact minimizer, a
balanced low-traffic cut cannot have most of $`C`$ coupling
$`\Omega(\sqrt n)`$-curvature subspaces of two nearly tight child SDP
slacks.  Equivalently, global minimality must force substantial cross
mass through both approximate kernels, after which a Boolean
near-ground realization theorem is still required.

#### 10.54.3 Exact grouped replenishment identities

Consider a nested random field-proportional peeling process.  At step
$`t`$, the current matrix has order $`m_t`$, a positive ground is
switched to $`1`$, and
```math
p_{t,i}
=\frac{r_{t,i}}{m_t-1}1_{\{i\in H_t\}},\qquad
a_t=p_t^{\mathsf T}r_t,\qquad
c_t=p_t^{\mathsf T}A_tp_t.
```
Delete vertices independently with probabilities $`p_{t,i}`$.  The
conditional identity (10.329) is
```math
2a_t=\mathbb E_t(d_t+g_t)+c_t.
```
Start from an exact order-$`n`$ minimizer and stop at random order
$`m_L`$.  With doubled minima $`q_m=2M_m`$ and terminal excess
```math
\varepsilon_L=Q(A_L)-q_{m_L}\ge0,
```
the decrements telescope pathwise.  Therefore
```math
\boxed{
\mathbb E\sum_{t<L}g_t
=
2\mathbb E\sum_{t<L}a_t
-\mathbb E\sum_{t<L}c_t
-q_n+\mathbb E q_{m_L}
+\mathbb E\varepsilon_L.
}
\tag{10.340}
```
All intermediate optimality excesses cancel.  Thus global minimality
contributes only the scalar endpoint term; it does not itself charge
leader switches.  The uncontrolled quantity is now the cumulative
signed internal term $`\sum_tc_t`$.

There is also a universal theorem using all successor energy layers.
For a disjoint block tower
```math
A_t=
\begin{pmatrix}
D_t&B_t\\
B_t^{\mathsf T}&A_{t+1}
\end{pmatrix},
```
let an oriented successor state $`z=(\sigma,y)`$ have deficit
```math
\delta_t(z)
=Q(A_{t+1})-\sigma y^{\mathsf T}A_{t+1}y,
```
and put
```math
\mathcal L_t
=
\sup_z[2\|B_ty\|_1-\delta_t(z)]_+.
```
Choosing the peeled-block signs to expose $`B_ty`$ proves
```math
\mathcal L_t\le d_t+Q(D_t).
```
The decrement and diagonal-block bounds telescope as before, giving
```math
\boxed{
\sum_t\mathcal L_t\le3Q(A_0).
}
\tag{10.341}
```
This already includes every near-ground layer, but it still does not
control replenishment.  For the predecessor restriction, the exact
signed relation is
```math
2b_t-g_t=d_t-h_t,
```
so the deficit penalty cancels $`g_t`$.  In the exact order-nine
minimizer from (10.332), a singleton step has
```math
d_t=h_t=0,\qquad g_t=16,\qquad \|B_ty\|_1=8,
```
and hence $`2\|B_ty\|_1-g_t=0`$.  The coefficient one on
$`\delta_t`$ in (10.341) is sharp even under global minimality.

Finally, Gibbs regularization gives an exact interpretation but no
automatic telescope.  Let $`\mu_A^\beta`$ be the Gibbs law on
oriented augmented states with score $`\sigma x^{\mathsf T}Ax`$.
For a split $`H\sqcup T`$, let $`\nu_H^\beta`$ be its core marginal,
$`\mu_E^\beta`$ the fresh-core Gibbs law, and
```math
K_H(\sigma,y)
=
\sum_u
\exp\!\left(
\beta\sigma(u^{\mathsf T}Du+2u^{\mathsf T}By)
\right).
```
Then
```math
\boxed{
\frac{d\nu_H^\beta}{d\mu_E^\beta}
=
\frac{K_H}{\mathbb E_{\mu_E^\beta}K_H},
\qquad
g_H(z)
=
\lim_{\beta\to\infty}
\frac1\beta
\log\frac{d\nu_H^\beta}{d\mu_E^\beta}(z)
}
\tag{10.342}
```
for the restriction $`z`$ of a parent absolute ground.  Replenishment
is exactly the zero-temperature information gain from the fresh-core
law to the parent restriction marginal.

Along a tower, the log partition functions and conditional entropies
in the associated KL divergences telescope.  The remaining term is
the expected soft block reward under the changing parent Gibbs laws;
at zero temperature it is precisely the FTL/mosaic payoff
$`\sum_t(d_t+g_t)`$.  Gibbs regularization therefore restates the
adaptivity gap rather than bounding it.

The surviving **open target** is a grouped layer-descent theorem:
produce successor states different from the predecessor such that
```math
2\|B_ty_t\|_1-\delta_t(y_t)
\ge\eta g_t-e_t,\qquad
\eta>0,\qquad
\sum_te_t=o(n^{3/2}).
```
Equation (10.341) would then bound cumulative replenishment at the
correct scale.  The order-nine equality shows that thickness,
inter-step coherence, or large-order exchange is indispensable.

#### 10.54.4 Updated frontier

This wave stops three generic mechanisms:

- endpoint disagreement does not aggregate without prime/global
  exchange and literal puncture compatibility;

- thick caps do not produce frames unless global exchange forces cross
  mass through child SDP-dual approximate kernels;

- all-layer visibility and scalar Gibbs entropy do not by themselves
  control replenishment.

The next wave should attack the common missing input directly:
global-minimizer block exchange.  Independent formulations are
cross-mass transfer through SDP kernels, grouped layer descent, and
prime puncture compatibility inside the endpoint affine coset.

### 10.55 Exact exchange and compatibility walls

The third block-exchange wave produced three exact statements.  Each
one narrows the missing input, but none yet yields the asymptotic
upper bound.  All claims in this section are **Verified** unless
explicitly labelled otherwise.

#### 10.55.1 Cross-only replacement is vacuous at an active zero cut

Write
```math
A(C)=
\begin{pmatrix}
B&C\\
C^{\mathsf T}&D
\end{pmatrix}.
```
Let `u` and `v` be positive grounds of `B` and `D`.  Changing the
relative sign of these two vectors proves, for every flat replacement
`C'`,
```math
\boxed{
p(A(C'))
\ge p(B)+p(D)+|u^{\mathsf T}C'v|
\ge p(B)+p(D).
}
\tag{10.343}
```
Consequently, if an active positive zero cut of a global minimizer has
`p(B)+p(D)=p(A)=M_n`, then every cross-only replacement already
satisfies `M(A(C'))\ge M_n`.  Global minimality imposes no further
condition on `C'`.  Away from exact activity, this comparison can
exploit at most the endpoint margin
`e_+=M_n-p(B)-p(D)`.

There is an exact SDP form of the same obstruction.  For
`\sigma\in\{+,-\}`, take optimal child SDP dual slacks
`L_{B,\sigma},L_{D,\sigma}` and put
```math
\begin{aligned}
S_\sigma
&=\operatorname{SDP}_+(\sigma B)
  +\operatorname{SDP}_+(\sigma D),\\
\Phi_\sigma(C)
&=\max_{u,v\in\{\pm1\}}
\left[
u^{\mathsf T}Cv
-\frac12u^{\mathsf T}L_{B,\sigma}u
-\frac12v^{\mathsf T}L_{D,\sigma}v
\right].
\end{aligned}
```
The dual identity (10.337), together with `v\mapsto-v`, gives
```math
\boxed{
p(A(C))=S_++\Phi_+(C),\qquad
\nu(A(C))=S_-+\Phi_-(C).
}
\tag{10.344}
```
Thus global minimality says only that every flat `C'` obeys
```math
M_n\le
\max\{S_++\Phi_+(C'),\,S_-+\Phi_-(C')\}.
```
A sufficient positive-endpoint completion would be
```math
\begin{pmatrix}
L_{B,+}&-C'\\
-C'^{\mathsf T}&L_{D,+}
\end{pmatrix}\succeq0,
```
because it forces `\Phi_+(C')\le0`.  But a positive semidefinite block
completion necessarily satisfies
```math
\boxed{
C'\ker L_{D,+}=0,\qquad
C'^{\mathsf T}\ker L_{B,+}=0.
}
\tag{10.345}
```
The analogous condition holds at the negative endpoint with the
appropriate sign.  Optimal SDP slacks are singular by complementary
slackness, so this kernel condition cannot be discarded.

The globally optimal order-nine signing from
`artifacts/prime_face_cover_quotient.md` gives an exact finite
obstruction.  After switching by positive-ground mask `247` and
using the zero cut
`U=\{0,1,2,4\}`, `V=\{3,5,6,7,8\}`, its blocks are
```math
B=J_4-I_4,\qquad
D=
\begin{pmatrix}
0&1&-1&1&1\\
1&0&1&1&-1\\
-1&1&0&1&1\\
1&1&1&0&1\\
1&-1&1&1&0
\end{pmatrix},
```
and
```math
C=
\begin{pmatrix}
1&1&-1&1&-1\\
-1&1&1&-1&-1\\
1&-1&1&-1&1\\
-1&-1&1&1&-1
\end{pmatrix}.
```
Here `p(B)=p(D)=6`, `1^{\mathsf T}C1=0`, and
`M(A)=M_9=12`.  Both child positive SDPs are exact.  Their canonical
slacks have spectra
```math
\operatorname{spec}(L_B)=\{0,4,4,4\},\qquad
\operatorname{spec}(L_D)=\{0,1,1,5,5\},
```
with all-one kernel lines.  Nevertheless
```math
C1_5=(1,-1,1,-1)^{\mathsf T},\qquad
1_4^{\mathsf T}C=(0,0,2,0,-2),
```
and, for `P_B=I-J_4/4` and `P_D=I-J_5/5`,
```math
\boxed{
\|P_BCP_D\|_F^2=\frac{86}{5}
=\frac{43}{50}\|C\|_F^2.
}
\tag{10.346}
```
Thus `86\%` of the cross mass couples the positive-curvature
subspaces despite exact child SDPs, global optimality, and zero
traffic.  More decisively, no flat `4\times5` matrix `C'` can satisfy
the positive completion: (10.345) would require `C'1_5=0`, whereas
every five-entry `\{\pm1\}` row has odd nonzero sum.

This falsifies a cross-only SDP completion or replacement proof.
A viable exchange must modify an internal block together with the
cross block, or exploit a quantitatively nonzero endpoint margin.

#### 10.55.2 Literal puncture compatibility has an exact criterion

Let `q` be an oriented endpoint word with
`\langle A,q\rangle=M(A)=M`.  Define its incident field and the
`i`th child defect by
```math
r_i(q)=\sum_{j\ne i}a_{ij}q_{ij},\qquad
d_i=M-M(A[-i]).
```
Since the punctured score is `M-r_i(q)`, one always has
```math
\boxed{
r_i(q)\ge d_i,
}
\tag{10.347}
```
with equality exactly when `q[-i]` is a closest oriented word for the
`i`th child.  Also every edge is counted twice, so
`\sum_i r_i(q)=2M`.  It follows that
```math
\boxed{
\text{a full endpoint word is closest on every puncture}
\quad\Longleftrightarrow\quad
\sum_i d_i=2M.
}
\tag{10.348}
```
Indeed, necessity follows by summing the equalities.  Conversely, if
the defect sum is `2M`, then the nonnegative quantities
`r_i(q)-d_i` sum to zero for every endpoint word `q`; hence every
endpoint word is simultaneously closest on every puncture.  The
equivalent compatibility identity is
```math
\sum_i M(A[-i])=(n-2)M.
```

For an affine endpoint coset `F=q_0L`, choosing a literal closest word
for puncture `i` is therefore exactly the minimization problem
```math
q_i\in\operatorname*{arg\,min}_{q\in F}r_i(q).
\tag{10.349}
```
After switching, the `r_i` are nonnegative Fourier sums on `L` and
their nonconstant Fourier coefficients cancel after summing over
`i`.  These facts do not force the coordinatewise minimizers to
intersect.

There is already a counterexample at the genuinely prime order
`n=5`:
```math
A=
\begin{pmatrix}
0&-1&1&-1&1\\
-1&0&-1&1&1\\
1&-1&0&1&1\\
-1&1&1&0&1\\
1&1&1&1&0
\end{pmatrix}.
\tag{10.350}
```
Its energies are `-4,0,4`, so `M(A)=4=M_5`.  The lower bound is
analytic: for a uniform Boolean spin `X`,
`\mathbb E H_A(X)^2=10`, while every energy is even, forcing
`M_5\ge4`.  Direct evaluation of the five displayed children gives
`M(A[-i])=4=M_4` for every `i`, hence every defect `d_i` is zero.

The endpoint face contains the affine plane
`q^{t,x}_{uv}=t x_ux_v` represented by
```math
\begin{array}{c|c}
t&x\\ \hline
-1&(1,1,1,1,-1)\\
-1&(1,1,-1,1,-1)\\
 1&(1,-1,1,-1,-1)\\
 1&(1,-1,-1,-1,-1).
\end{array}
```
The first three edge words multiply coordinatewise to the fourth.
Their incident-field profiles are
```math
\begin{pmatrix}
2&2&0&0&4\\
4&0&0&2&2\\
2&4&0&2&0\\
0&2&0&4&2
\end{pmatrix}.
\tag{10.351}
```
The zero sets cover all five punctures, so each child has a literal
closest representative in this same affine plane.  But
`\sum_i d_i=0<8=2M`, and (10.348) proves that no single representative
works for all punctures.

Thus prime order, global optimality, optimal children, literal
closeness, and a bounded affine endpoint face still do not give a
common compatible word.  The order-nine face used earlier has an
analogous puncture-covering plane, but `9` is composite and all of its
child defects equal `2`.  The surviving arithmetic target must use
large primes (beginning at `13`), growing affine dimension, bounded
endpoint reuse, or some additional large-order rigidity.

#### 10.55.3 A two-shore descent dichotomy

Use doubled normalization.  Let a positive ground of an order-`m`
signing `C` be switched to `1`, put `r=C1\ge0`, and let
`v=1^S` be the spin obtained by flipping the shore `S`.  Set
`T=S^c` and write
```math
C=
\begin{pmatrix}
D&B\\
B^{\mathsf T}&E
\end{pmatrix},\qquad
g=Q(C)-v^{\mathsf T}Cv.
```
Define
```math
H_D=1^{\mathsf T}D1,\quad
H_E=1^{\mathsf T}E1,\quad
\gamma=Q(E)-H_E,\quad
d=Q(C)-Q(E),
```
and let
`c_C(S)=1_S^{\mathsf T}B1_T` and
`R_S=\sum_{i\in S}r_i`.  Direct expansion gives the exact identities
```math
\boxed{
c_C(S)=\frac g4,\qquad
R_S=H_D+\frac g4,\qquad
d=H_D+\frac g2-\gamma.
}
\tag{10.352}
```
Use the predecessor restriction `y=1_T` as an all-layer successor.
Its layer payoff is
```math
\ell=2\|By\|_1-\gamma.
```
Since `\|By\|_1\ge1_S^{\mathsf T}By=g/4` and `d\ge0`,
```math
\boxed{
\ell\ge\frac g2-\gamma
=d-H_D
=d+\frac g4-R_S.
}
\tag{10.353}
```
Consequently, for every `\theta<1/4`,
```math
\boxed{
\ell<\theta g
\quad\Longrightarrow\quad
R_S>\left(\frac14-\theta\right)g.
}
\tag{10.354}
```
In particular, either `\ell\ge g/8` or `R_S>g/8`.  This is a genuine
descent dichotomy: a poorly visible predecessor gap must put positive
field mass on the flipped shore.

Field-proportional sampling exposes the remaining covariance problem.
With `p_i=r_i/(m-1)` on `S`,
```math
\sum_{i\in S}p_i=\frac{R_S}{m-1},\qquad
a_S=\frac1{m-1}\sum_{i\in S}r_i^2
\ge\frac{R_S^2}{(m-1)|S|}.
\tag{10.355}
```
Write
```math
\beta_i=\sum_{j\in T}C_{ij},\qquad
\alpha_i=\sum_{j\in S\setminus\{i\}}C_{ij},
```
so `r_i=\alpha_i+\beta_i\ge0` and
`\sum_{i\in S}\beta_i=g/4`.  The expected signed boundary removed
when sampling only `S` is
```math
\frac1{m-1}\sum_{i\in S}r_i\beta_i.
```
Neither `r_i\ge0` nor the positive total of the `\beta_i` controls
this covariance: vertices with positive `r_i` can have negative
`\beta_i`, while positive cross field can be cancelled internally so
that `r_i=0`.  The two-shore weight
`p_i+p_j-p_ip_j` has the same unresolved correlation.

There is also a recursive form of the wall.  If
`\ell<\eta g`, then (10.353) gives
`\gamma>(1/2-\eta)g`, while `d\ge0` gives
`\gamma\le g/2+Q(D)`.  Hence the gap can migrate to the next core at
roughly half size without ever paying a fixed visible fraction.  The
statement concerns the same augmented orientation; an
opposite-orientation transition is a complemented-cut problem.

The sharpened **open target** is a global-minimizer exchange estimate
on every bad shore, for some absolute `\kappa>0`:
```math
\boxed{
\sum_{i\in S}r_i\beta_i
\ge \kappa(m-1)c_C(S)-\operatorname{Err}_S,
\qquad
\sum_t\operatorname{Err}_{S_t}=o(n^{3/2}),
}
\tag{10.356}
```
or a two-shore analogue.  Such an estimate must exchange matched
positive cross edges together with internal edges; cross-only
minimality is vacuous by Section 10.55.1.

#### 10.55.4 Updated frontier

The common obstruction is now precise.  A global minimizer can have
an active zero cut whose cross block lies mostly in the curved child
SDP subspaces; a prime optimal signing can have compatible punctures
only one at a time; and a large cut gap can hide behind negative
field-boundary covariance.

The next proof wave should therefore prioritize genuinely coupled
internal/cross exchanges.  Independent secondary targets are a
large-prime puncture rigidity theorem (not a bounded-face theorem
already refuted at `n=5`), and an exchange using nonzero endpoint
margin or opposite orientation.  Any proposed lemma must survive the
exact order-nine zero-cut obstruction and the prime order-five affine
plane above.

### 10.56 Cut-space exchange, affine puncture families, and reset descent

The fourth wave found that two of the proposed exchanges were still
too close to gauge symmetries.  It also found the correct
orientation-safe way to distribute a leader gap.  All algebraic and
finite claims below are **Verified**; explicitly identified
enumeration counts are **Numerical**.

#### 10.56.1 Internal stars do not make a coupled exchange

Let `A^F` denote the signing obtained by flipping an edge set `F`,
and split the vertices as `S\sqcup T`.  For `R\subseteq S` and
`U\subseteq T`, let
```math
K_S=\delta_{A[S]}(R),\qquad
K_T=\delta_{A[T]}(U),
```
and define the cross rectangle
```math
X_{R,U}
=
\bigl[R\times(T\setminus U)\bigr]
\sqcup
\bigl[(S\setminus R)\times U\bigr].
```
The full cut `\delta_A(R\cup U)` is the disjoint union
`K_S\sqcup K_T\sqcup X_{R,U}`.  Hence, for every cross-edge set `G`,
```math
F=K_S\mathbin\triangle K_T\mathbin\triangle G
\quad\Longrightarrow\quad
F\mathbin\triangle\delta_A(R\cup U)
=G\mathbin\triangle X_{R,U}\subseteq E(S,T).
```
Flipping a full vertex cut is a switching equivalence, so
```math
\boxed{
M(A^F)
=
M\!\left(A^{\,G\mathbin\triangle X_{R,U}}\right).
}
\tag{10.357}
```
Thus arbitrary internal stars on both shores, even when coupled to
arbitrary cross rectangles, are exactly cross-only exchanges in
disguise.  At an active positive block cut
```math
p(A[S])+p(A[T])=M(A),
```
every signing on the right of (10.357) is already stable by
(10.343).  This whole exchange class is therefore vacuous.

Conversely, a coupled flip set has a cross-only representative modulo
switching if and only if both of its internal pieces are cuts.  In a
complete graph an edge set is a cut exactly when it meets every
triangle evenly.  Therefore every genuinely informative internal
residual has an odd-parity triangle witness.  On a shore of size `s`,
the residual quotient has dimension
```math
\boxed{
\binom{s}{2}-(s-1)=\binom{s-1}{2}.
}
\tag{10.358}
```
The next coupled exchange must operate in this triangle/cycle
quotient; internal stars cannot supply the missing covariance.

The prime order-five minimizer (10.350) also gives the exact wall for
an *unconditional* covariance inequality.  In doubled normalization,
use its all-one ground and `S=\{0,1\}`.  Then
```math
\begin{gathered}
q=8,\qquad
r=(0,0,2,2,4),\qquad
\beta=(1,1),\\
c=2,\quad g=8,\quad
H_D=-2,\quad H_E=6,\quad
\gamma=0,\quad d=2,\quad R_S=0,
\end{gathered}
```
while
```math
\boxed{
\sum_{i\in S}r_i\beta_i=0,\qquad
(m-1)c=8,\qquad
\ell=4=\frac g2.
}
\tag{10.359}
```
Consequently no positive `\kappa` can make
`\sum r_i\beta_i\ge\kappa(m-1)c` hold on every shore of every global
minimizer.  This does **not** falsify the narrowed bad-shore target:
the zero covariance in (10.359) occurs at maximal predecessor
visibility.

A useful opposite-endpoint-margin inequality survives.  In a
positive-ground gauge let
```math
c=c_A(S),\qquad
\Delta=\frac q2-c=\frac{2q-g}{4},\qquad
K=\sum_{i\in S}r_i\beta_i,\qquad
D=A[S].
```
Choose `R\subseteq S` independently with
`\Pr(i\in R)=t r_i`, where
`0<t\le1/\max_{i\in S}r_i`.  Direct edge expansion gives
```math
\boxed{
\mathbb E\,c_A(S\setminus R)
=
c+t\bigl(\|r_S\|_2^2-2K\bigr)
-t^2r_S^{\mathsf T}Dr_S.
}
\tag{10.360}
```
The absolute-ground cap says
`0\le c_A(H)\le q/2` for every `H`.  Applying its upper half to
(10.360) yields
```math
\boxed{
K\ge
\frac12\left(
\|r_S\|_2^2
-t\,r_S^{\mathsf T}Dr_S
-\frac{\Delta}{t}
\right).
}
\tag{10.361}
```
At an exact opposite endpoint, `\Delta=0`, letting `t\downarrow0`
proves
```math
K\ge\frac12\|r_S\|_2^2.
```
Equivalently, the pointwise cut comparison
`c_A(S\setminus\{i\})\le c_A(S)=q/2` gives
`\beta_i\ge r_i/2`.  Since
`r_S^{\mathsf T}Dr_S\le(|S|-1)\|r_S\|_2^2`, taking
`t=1/(m-1)` gives the stable form
```math
\boxed{
K\ge
\frac{|T|}{2(m-1)}\|r_S\|_2^2
-\frac{(m-1)(2q-g)}8.
}
\tag{10.362}
```
Its margin error is generally too large and is not known to
telescope, but it identifies the correct favorable regime: cuts
close to the opposite endpoint force positive field-boundary
covariance.

**Numerical exact audit.**  Across every oriented absolute-ground
gauge and rooted shore of the order-nine minimizer (5.1), no negative
`K` occurs.  There are `18` positive-traffic zero-covariance shores,
all with `\ell/g=1/2`.  Among the `636` shores with `\ell<g/8`,
```math
\min\frac{K}{(m-1)c}=\frac14;
```
among the four with `\ell<0`, the minimum is `5/8`.  These finite
counts support a clipped bad-shore theorem but do not prove an
asymptotic constant.

#### 10.56.2 Rank-two puncture obstructions persist at every prime residue order

There is a convenient normal form for an affine endpoint plane.
After switching one word to `q_0=1`, write
```math
q_h(i,j)
=
(-1)^{\langle h,\tau_i+\tau_j+\omega\rangle},
\qquad
h\in\mathbb F_2^2,
\tag{10.363}
```
where `\tau_i\in\mathbb F_2^2` is the vertex type and `\omega` is the
orientation character.  Put
```math
W_z
=
\sum_{\tau_i+\tau_j+\omega=z}a_{ij}.
```
Fourier inversion gives the exact endpoint criterion
```math
\boxed{
\langle A,q_h\rangle=M\quad(h\in\mathbb F_2^2)
\quad\Longleftrightarrow\quad
W_0=M,\quad W_z=0\ (z\ne0).
}
\tag{10.364}
```
For a vertex `i` of type `a`, let
```math
T_{i,z}
=
\sum_{\tau_j=a+\omega+z}a_{ij}.
```
Then its four fields are the Walsh transform
```math
r_i(h)=\sum_zT_{i,z}(-1)^{\langle h,z\rangle}.
```
The closest-puncture condition within this plane is exactly
`\min_h r_i(h)=0` for every `i`.  Once the four words are actual
endpoints, that condition automatically gives `d_i=0`: a zero field
retains score `M` after deletion, while a principal child cannot have
norm larger than its parent.

Prime arithmetic cannot rule out this configuration.  Let `m` be
even, put `n=2m+1`, and partition the vertices as
```math
X\sqcup\{z\}\sqcup Y,\qquad |X|=|Y|=m.
```
Choose balanced `s,t\in\{\pm1\}^m` and define
```math
A_m=
\begin{pmatrix}
J_m-I_m&s&J_m\\
s^{\mathsf T}&0&t^{\mathsf T}\\
J_m&t&I_m-J_m
\end{pmatrix}.
\tag{10.365}
```
Assign types
```math
\tau(X)=0,\qquad\tau(z)=1,\qquad\tau(Y)=2,\qquad
\omega=2.
```
The four oriented words are represented by
```math
\begin{array}{c|c}
\text{orientation}&\text{spin on }(X,z,Y)\\ \hline
+&(1_X, 1, 1_Y)\\
+&(1_X,-1, 1_Y)\\
-&(1_X, 1,-1_Y)\\
-&(1_X,-1,-1_Y).
\end{array}
```
They form an affine plane and all have oriented score `m^2`.

For a general spin `(x,\varepsilon,y)`, set
```math
S=1^{\mathsf T}x,\quad T=1^{\mathsf T}y,\quad
U=s^{\mathsf T}x,\quad V=t^{\mathsf T}y.
```
Its one-copy energy is
```math
H
=
\frac12(S^2-T^2)+ST+\varepsilon(U+V).
\tag{10.366}
```
Balance implies
`|S|+|U|\le m` and `|T|+|V|\le m`.  Put
`a=|S|`, `b=|T|` and write
```math
\max(a,b)=m-p,\qquad
\min(a,b)=m-q,\qquad q=p+2r.
```
Because `m` is even, `p,q` are even.  The exact gap computation is
```math
\begin{aligned}
m^2
&-
\left[
\frac12|a^2-b^2|+ab+2m-a-b
\right]\\
&=
p(2m-p-2)+2r(r-1)\ge0.
\end{aligned}
```
Type-constant states attain both signs, and hence
```math
\boxed{
M(A_m)=m^2=\frac{(n-1)^2}{4}.
}
\tag{10.367}
```

The complete field table is
```math
\begin{array}{c|c|c|c}
h&r_X(h)&r_z(h)&r_Y(h)\\ \hline
0&2m-1+s_i&0&1+t_j\\
1&2m-1-s_i&0&1-t_j\\
2&1-s_i&0&2m-1+t_j\\
3&1+s_i&0&2m-1-t_j.
\end{array}
\tag{10.368}
```
Every entry is nonnegative and every vertex column contains a zero.
Therefore
```math
\boxed{
M(A_m[-i])=m^2\qquad(i\in[n]).
}
```
This produces the puncture-covering rank-two obstruction at every
order `n\equiv1\pmod4`, including arbitrarily large primes.  Its norm
is quadratic, so it does not refute a theorem using global
minimality or an asymptotic `O(n^{3/2})` cap.

There is nevertheless an exact competitive-scale finite example at
the prime `n=13`:
```math
A=
\begin{pmatrix}
0&1&-1&1&1&-1&-1&1&1&-1&1&1&-1\\
1&0&1&1&1&1&-1&1&1&1&1&1&1\\
-1&1&0&1&1&1&1&1&1&1&1&-1&1\\
1&1&1&0&-1&1&-1&1&1&-1&1&1&1\\
1&1&1&-1&0&1&1&1&1&1&1&1&-1\\
-1&1&1&1&1&0&1&1&1&1&-1&1&1\\
-1&-1&1&-1&1&1&0&1&-1&1&1&-1&-1\\
1&1&1&1&1&1&1&0&-1&-1&-1&-1&-1\\
1&1&1&1&1&1&-1&-1&0&-1&-1&-1&-1\\
-1&1&1&-1&1&1&1&-1&-1&0&-1&1&1\\
1&1&1&1&1&-1&1&-1&-1&-1&0&-1&1\\
1&1&-1&1&1&1&-1&-1&-1&1&-1&0&-1\\
-1&1&1&1&-1&1&-1&-1&-1&1&1&-1&0
\end{pmatrix}.
\tag{10.369}
```
Use type list
```math
(0,0,0,0,0,0,1,2,2,2,2,2,2)
```
and `\omega=2`.  The four field profiles are
```math
\begin{pmatrix}
2&10&8&6&8&8&0&2&0&2&2&0&0\\
4&12&6&8&6&6&0&0&2&0&0&2&2\\
2&2&0&2&0&0&0&12&10&4&8&6&2\\
0&0&2&0&2&2&0&10&12&2&6&8&4
\end{pmatrix}.
```
Every column contains zero.  Exact enumeration of the `2^{12}`
projective spins and the thirteen `2^{11}` child spin spaces gives
```math
\boxed{
\min H_A=-24,\quad
\max H_A=24,\quad
M(A[-i])=24\ \text{for every }i.
}
\tag{10.370}
```
This example is not globally minimizing.  The quadratic-residue
circulant signing with `P_{ii}=0` and, for `i\ne j`,
```math
P_{ij}=1
\quad\Longleftrightarrow\quad
i-j\pmod {13}\in\{1,3,4,9,10,12\}
```
has exact enumerated norm `20`.  Thus `M_{13}\le20<24`.  The example
nevertheless falsifies every finite theorem based only on prime
order, flatness, four exact endpoints, puncture-zero coverage, and
`d_i=0`, even at ratio
`24/13^{3/2}=0.512\ldots`.  Global minimality or a genuinely
asymptotic competitive hypothesis is indispensable.

#### 10.56.3 Orientation reset and the clipped backward residual

Return to doubled normalization and the two-shore notation of
Section 10.55.3.  In addition to the forward core `E`, retain the
backward or peeled core `D`.  Put
```math
q_D=Q(D),\qquad q_E=Q(E),\qquad
\zeta=q_D+q_E-q,
```
and define the two same-orientation raw payoffs
```math
\ell_F
=2\|B1_T\|_1-(q_E-H_E),\qquad
\ell_B
=2\|B^{\mathsf T}1_S\|_1-(q_D-H_D).
```
Since each cross norm is at least `g/4` and
```math
(q_D-H_D)+(q_E-H_E)=\zeta+\frac g2,
```
one obtains
```math
\boxed{
\ell_F+\ell_B\ge\frac g2-\zeta.
}
\tag{10.371}
```
Consequently,
```math
\boxed{
\max\{(\ell_F)_+,(\ell_B)_+\}\ge\frac g8
\quad\text{or}\quad
\zeta>\frac g4.
}
\tag{10.372}
```
If `\zeta>0`, the two child absolute norms must occur in opposite
exclusive orientations.  If they shared a maximizing orientation,
the two child grounds and a choice of relative global sign would give
`Q(C)\ge q_D+q_E`.

There is also an orientation-safe small-shore estimate.  Let `X` be
either shore, of size `s`, choose an orientation `\tau_X` attaining
`Q(X)`, and put
```math
g_X=Q(X)-\tau_XH_X,\qquad
\lambda_X=2L_X-g_X,
```
where `L_X` is the corresponding cross `\ell_1` norm.  Since
`2L_X\ge g/2` and `g_X\le2Q(X)`,
```math
\boxed{
g\le2(\lambda_X)_++4Q(X)
\le2(\lambda_X)_++4s(s-1).
}
\tag{10.373}
```
Thus a gap invisible from the smaller shore forces that shore to
have size `\Omega(\sqrt g)`.

The augmented orientation cannot simply be ignored.  For a carried
`\sigma\in\{\pm1\}`, define
```math
g=q-\sigma v^{\mathsf T}Cv,\qquad
a=q-v^{\mathsf T}Cv.
```
Then
```math
a=
\begin{cases}
g,&\sigma=+,\\
2q-g,&\sigma=-.
\end{cases}
\tag{10.374}
```
With carried child gaps
`g_D=q_D-\sigma H_D` and `g_E=q_E-\sigma H_E`, the two directed
payoffs satisfy
```math
\boxed{
\begin{aligned}
\sigma=+:\quad&
\ell_D^++\ell_E^+
\ge\frac g2-\zeta,\\
\sigma=-:\quad&
\ell_D^-+\ell_E^-
\ge q-\frac{3g}{2}-\zeta.
\end{aligned}
}
\tag{10.375}
```
A pure orientation mismatch has `g=2q` but `a=0`.  It carries no
ordinary cut layer, so no scalar potential depending only on `g` can
survive arbitrary orientation changes.

There is a clean repair.  At every positive-oriented node, choose
independently in each child `X` an exact endpoint
`(\tau_X,w_X)` with
```math
\tau_Xw_X^{\mathsf T}Xw_X=Q(X).
```
Normalize the child by this endpoint, define its inherited gap `g_X`
and directed payoff `\lambda_X=2L_X-g_X`, and recurse on both
children.  At one node,
```math
g\le2L_D+2L_E
\le
(\lambda_D)_++(\lambda_E)_++g_D+g_E.
\tag{10.376}
```
All internal child gaps cancel on recursion, while singleton and
zero-gap leaves contribute nothing.  Therefore
```math
\boxed{
g_{\rm root}
\le
\sum_{\text{reset-tree edges}}(\lambda_e)_+.
}
\tag{10.377}
```
This is a fully orientation-safe binary distribution theorem.  It
does not yet bound the allocated tree capacity.

The smallest exact obstruction to a naive same-channel recursion is
```math
C=
\begin{pmatrix}
0&1&1&1&1\\
1&0&1&-1&1\\
1&1&0&1&-1\\
1&-1&1&0&-1\\
1&1&-1&-1&0
\end{pmatrix}.
\tag{10.378}
```
Its doubled energies are `-8,0,8`.  For
`v=(-1,-1,-1,1,1)` and the split
`D=\{0,1,2\}`, `E=\{3,4\}`,
```math
g=8,\quad
q_D=H_D=6,\quad
q_E=2,\quad H_E=-2,\quad\zeta=0.
```
The same-positive payoffs are `4` and `0`, with child gaps `0` and
`4`; the latter is a pure orientation mismatch and has no later cut
layer.  Thus the naive same-channel tree sees only `4<g`.
Orientation reset chooses the positive endpoint in `D` and the
negative endpoint in `E`; both reset gaps vanish and the payoffs are
`4+4=g`.

Finally, along a nested path on which the inherited orientation is
the positive channel, define
```math
f_t=(\ell_{F,t})_+,\qquad
h_t=
\left[
\frac{g_t}{2}-\zeta_t-f_t
\right]_+.
\tag{10.379}
```
Equation (10.371) gives
`h_t\le(\ell_{B,t})_+` and
```math
\frac{g_t}{2}\le f_t+h_t+\zeta_t.
```
The surplus telescopes *with its sign*:
```math
\sum_t\zeta_t
=
\sum_tQ(D_t)+Q(A_L)-Q(A_0).
\tag{10.380}
```
The peeled blocks together with `A_L` form one disjoint partition.
Put each block in a class according to an orientation attaining its
norm.  Multiblock endpoint superadditivity (choose the corresponding
block grounds and average over their independent relative signs)
shows that the positive-class norms sum to at most `P(A_0)` and the
negative-class norms to at most `N(A_0)`.  Hence
```math
\sum_tQ(D_t)+Q(A_L)
\le P(A_0)+N(A_0)\le2Q(A_0),
```
and so `\sum_t\zeta_t\le Q(A_0)`.  Combining this with the forward
path bound `\sum_t f_t\le3Q(A_0)` from (10.341) gives
```math
\boxed{
\sum_tg_t
\le
2\sum_th_t+8Q(A_0).
}
\tag{10.381}
```
This is the sharpest current reduction of gap migration.  Raw
backward capacity is the wrong object: in (10.378), the singleton
shore `S=\{3\}` has
```math
g=\zeta=\ell_F=0,\qquad\ell_B=8,
```
but (10.379) correctly assigns `h=0`.

The surviving **open target** is therefore the clipped triangular
estimate
```math
\boxed{
\sum_t h_t=O(Q(A_0)),
}
\tag{10.382}
```
or an asymptotically sharper version sufficient for the insertion
recurrence.  It must control only backward payoff not already paid by
forward visibility or signed child-norm surplus.

#### 10.56.4 Updated frontier

This wave makes the next three routes substantially narrower:

- a coupled exchange must contain a genuine triangle/cycle residual
  inside a shore; internal stars and rectangles are gauge-equivalent
  to the already-vacuous cross-only problem;

- prime order, a rank-two exact endpoint plane, literal puncture
  coverage, and zero child defects coexist at every relevant residue
  order.  Global minimality or an asymptotic competitive cap must do
  the work;

- leader gaps admit an exact orientation-reset binary distribution,
  and same-channel path migration reduces to the clipped residual
  (10.379).  Charging all backward or all tree capacity to the gaps
  is false at zero-gap nodes, while bounding it wholesale is
  unnecessarily strong.

The highest-priority next step is to combine the first and third
points: prove that a large clipped backward allocation forces a
non-cut internal triangle residual whose global flip certificate pays
for it.  Independent alternatives are to exclude the rank-two
puncture family under *global* minimality, or to prove a laminar
allocated-capacity theorem for the reset tree.

### 10.57 Fixed-block charging, residual no-go, and affine resonance

The fifth wave closes the clipped estimate on every compatible
fixed-orientation path.  It also shows that the proposed pointwise
triangle explanation is false, and that the matched-degree extension
of the affine puncture family cannot have competitive norm.  All
algebraic statements below are **Verified**.  Finite searches are
labelled **Numerical**.

#### 10.57.1 The clipped residual has an exact fixed-block charge

Let
`A_t=A_0[R_t]`.  At step `t`, let `y_t` be the inherited state and
let `x_t` be a fresh endpoint in one fixed positive orientation.  Put
```math
D_t=\{i:x_{t,i}=-y_{t,i}\},\qquad
R_{t+1}=\{i:x_{t,i}=y_{t,i}\},\qquad
y_{t+1}=x_t|_{R_{t+1}}.
```
In the `x_t` gauge, write
```math
c_t
=
x_t(D_t)^{\mathsf T}A[D_t,R_{t+1}]x_t(R_{t+1})
=\frac{g_t}{4},
\qquad
L_{F,t}=\|B_t1\|_1.
```
With
```math
\gamma_{D,t}=Q(D_t)-H_{D,t},\qquad
\gamma_{E,t}=Q(A_{t+1})-H_{E,t},
```
the definition (10.379) simplifies exactly to
```math
\boxed{
h_t=
\left[
g_t-\gamma_{D,t}
-\max\{\gamma_{E,t},2L_{F,t}\}
\right]_+.
}
\tag{10.383}
```
Indeed,
`\zeta_t=\gamma_{D,t}+\gamma_{E,t}-g_t/2`.
Since `L_{F,t}\ge c_t\ge0`, this immediately gives
```math
\boxed{
h_t
\le[4c_t-2L_{F,t}]_+
\le2c_t.
}
\tag{10.384}
```
The clipping is essential here: raw backward capacity need not vanish
when `c_t=0`.

There is a single global Boolean state which charges all the `c_t`.
The blocks
```math
D_0,D_1,\ldots,D_{L-1},D_L=R_L
```
partition `R_0`.  Define `s` by
```math
s|_{D_t}=x_t|_{D_t}\quad(t<L),
\qquad
s|_{D_L}=-y_L.
```
Compatibility gives, by induction, that `x_t=s` on `D_t` and
`x_t=-s` on every later block.  Consequently
```math
\boxed{
c_t
=-
\sum_{j>t}
s_{D_t}^{\mathsf T}A[D_t,D_j]s_{D_j}.
}
\tag{10.385}
```
Put
`H_j=s_{D_j}^{\mathsf T}A[D_j]s_{D_j}`.  Summing (10.385)
over all ordered block pairs yields the exact identity
```math
\boxed{
2\sum_t c_t
=
\sum_{j=0}^{L}H_j-s^{\mathsf T}A_0s.
}
\tag{10.386}
```
Independently randomizing the global sign on every fixed block state
has expected full energy `\sum_jH_j`.  Thus
`\sum_jH_j\le P(A_0)`, while
`-s^{\mathsf T}A_0s\le N(A_0)`.  Equations (10.384)--(10.386)
prove
```math
\boxed{
\sum_t h_t
\le2\sum_t c_t
\le P(A_0)+N(A_0)
\le2Q(A_0).
}
\tag{10.387}
```
This proves the open target (10.382) in its stated compatible
same-positive-orientation setting.  Combining (10.387) with (10.381)
gives the explicit bound
```math
\boxed{
\sum_tg_t
\le2\bigl(P(A_0)+N(A_0)\bigr)+8Q(A_0)
\le12Q(A_0).
}
\tag{10.388}
```
No global minimality, flatness, or asymptotics enter this argument.

There is a useful extension across orientation changes along one
path.  Write
```math
\mathcal R(A)=P(A)+N(A).
```
Suppose a fixed-orientation run `[a,b]` is followed at
`r=b+1` by an opposite-orientation reset.  If `q_r=Q(A_r)` and
`g_r` is the inherited augmented gap, let
```math
a_r=2q_r-g_r
```
be the reset's ordinary cut deficit.  Reserve an old-orientation
ground on the terminal core in the block-sign averaging argument.
Its improvement over the inherited terminal state is exactly
`\mathcal R(A_r)-a_r`.  Therefore
```math
\boxed{
2\sum_{t=a}^{b}c_t
\le
\mathcal R(A_a)-\mathcal R(A_r)+a_r.
}
\tag{10.389}
```
Both `P` and `N` are monotone under principal restriction, so the
range terms telescope between consecutive runs.  Along any compatible
path with resets,
```math
\boxed{
\sum_{\sigma_t=+}h_t
\le
\mathcal R(A_0)
+\sum_{\sigma_t=-}a_t.
}
\tag{10.390}
```
A pure orientation mismatch has `a_t=0`, as it should.  What remains
open is a bound for the cumulative nonzero reset layers and a
tensorization across the branches of the reset tree; (10.390) does
not by itself close that problem.

**Numerical audit.**  Exact brute force checked (10.383)--(10.387) on
2,100 random signings and compatible tower starts of orders `3`--`9`.
A separate audit deliberately encouraged orientation changes on 3,500
random paths, totaling 5,975 steps, and checked (10.389)--(10.390).
These computations only audit the normalization and implementation;
the proofs above are exact.

#### 10.57.2 Triangle parity cannot explain clipped mass pointwise

For one node, retain
```math
c=1_S^{\mathsf T}B1_T,\qquad
N_S=\sum_{i\in S}\bigl(-(B1_T)_i\bigr)_+.
```
Since `g=4c` and `\|B1_T\|_1=c+2N_S`, (10.383) becomes
```math
\boxed{
h=
\left[
4c-\gamma_D-
\max\{\gamma_E,2c+4N_S\}
\right]_+
\le[2c-\gamma_D-4N_S]_+.
}
\tag{10.391}
```
This does not force any non-cut internal exchange.  In fact, for
every singleton shore `S=\{i\}` in a positive-ground gauge,
```math
\boxed{h=2r_i.}
\tag{10.392}
```
To see this, put `d=Q(C)-Q(C[-i])\ge0`.  Then
```math
c=r_i,\qquad \gamma_D=0,\qquad
\gamma_E=2r_i-d,\qquad
\ell_F=d,\qquad \zeta=-d,
```
and (10.392) follows.  Yet a singleton has no internal edge, and the
internal quotient of a two-vertex shore is also zero.  Thus the
pointwise triangle-residual route proposed at the end of Section
10.56 is **Falsified**; any such exchange must first group many
nodes.

For completeness, the coupled-flip quotient has the following exact
normal form.  Choose a root in each shore.  After switching away all
root-incident internal flip bits, an arbitrary flip set is uniquely
represented by
```math
F^\circ=R_S\sqcup R_T\sqcup G,
```
where `G` is cross-shore and, for `i,j` away from the root `s_0`,
```math
1_{R_S}(ij)
=
1_F(ij)+1_F(s_0i)+1_F(s_0j)\pmod2.
```
Thus the internal coordinates are exactly odd parities on rooted
triangles.  Let `H` be their support graph and
`\beta(H)=|E(H)|-|V(H)|+\operatorname{comp}(H)`.  The codimension of
the augmented-cut restrictions on `H` is
```math
\boxed{
q(H)=
\begin{cases}
\beta(H),&H\text{ bipartite},\\
\beta(H)-1,&H\text{ nonbipartite}.
\end{cases}
}
\tag{10.393}
```
In particular, a forest and a single triangle both have codimension
zero.  The first locally coercive supports are a four-cycle or two
independent odd cycles.  Since an arbitrary edge pattern is within
`q(H)` changes of an allowed restriction, parity alone can change a
doubled residual reward by at most `8q(H)`.  One odd triangle is not
a quantitative payment certificate.

There is nevertheless an exact conditional certificate worth
retaining.  For an oriented augmented-cut word `z`, put
```math
\Delta_A(z)=Q(A)-\sigma x^{\mathsf T}Ax,
\qquad z_{ij}=\sigma x_ix_j.
```
For `R=R_S\sqcup R_T`, define
```math
Y_R(y)=-4\sum_{e\in R}a_ey_e,
```
and, with cross residual `G`,
```math
\Gamma_{G,H}(y)
=
\min_{z:z|_H=y}
\left[
\Delta_A(z)+4\sum_{e\in G}a_ez_e
\right].
```
Directly regrouping the exact edge-flip certificate gives
```math
\boxed{
Q(A^{G\sqcup R})-Q(A)
=
\max_{y\in\mathcal C_n|_H}
\bigl[Y_R(y)-\Gamma_{G,H}(y)\bigr].
}
\tag{10.394}
```
This isolates the genuine target: control the cheaply extensible
patterns of `\Gamma_{G,H}`.  Triangle parity determines its domain,
but not its values.

The small exact minimizers show that this is a real obstruction.
For (10.350), switch by `(-1,1,1,1,1)` and take `S=\{3\}`;
then
```math
q=8,\quad c=4,\quad\gamma_E=8,\quad
\ell_F=\zeta=0,\quad h=8,
```
with no internal residual.  In the all-one gauge of the order-nine
minimizer (5.1), `S=\{1,4\}` gives
```math
q=24,\quad c=10,\quad\gamma_E=20,\quad
\ell_F=\zeta=0,\quad h=20,
```
while its only internal edge is itself a cut.

**Numerical, exhaustive audit.**  At `n=5`, 66 of the 150 endpoint
gauge/rooted-shore nodes have positive `h`, and 54 of those have shore
size at most two.  At `n=9`, the corresponding counts are 3,636 of
5,355 and 624.  All 949 positive-`h` rooted three-shores at `n=9`
have every non-cut internal `K_3` flip certified by an original
zero-deficit endpoint.  For `S=\{0,1,4\}`, one has `h=20`, but the
best forced doubled endpoint payment is only `4`.

#### 10.57.3 The matched-degree affine extension is necessarily quadratic

The natural attempt to deform (10.365) toward competitive scale has
an exact Boolean resonance.  Work in the one-copy normalization of
Section 10.56.2.  Let `m` be even, let balanced
`s,t\in\{\pm1\}^m` index shores `X,Y`, and let `G,H` be simple
graphs on the two shores.  Let `K` be bipartite from `X` to `Y`, with
```math
\deg_K^{\rm row}(i)=\deg_G(i),\qquad
\deg_K^{\rm col}(j)=\deg_H(j).
```
Put
```math
B=J-I-2G,\qquad
D=I-J+2H,\qquad
C=J-2K,
```
and
```math
A=
\begin{pmatrix}
B&s&C\\
s^{\mathsf T}&0&t^{\mathsf T}\\
C^{\mathsf T}&t&D
\end{pmatrix}.
```
Degree matching gives
`|E(G)|=|E(H)|=:e` and `|E(K)|=2e`.  The four affine-plane
words from (10.365) have common oriented score
```math
\boxed{M_0=m^2-4e.}
\tag{10.395}
```
If `g_i=\deg_G(i)`, `h_j=\deg_H(j)`, and
```math
b_i=m-1-2g_i,\qquad a_j=m-1-2h_j,
```
their four field profiles are
```math
\begin{aligned}
X_i:&\quad
(2b_i+1+s_i,\ 2b_i+1-s_i,\ 1-s_i,\ 1+s_i),\\
z:&\quad(0,0,0,0),\\
Y_j:&\quad
(1+t_j,\ 1-t_j,\ 2a_j+1+t_j,\ 2a_j+1-t_j).
\end{aligned}
\tag{10.396}
```
Assume the proposed plane words are actual endpoints.  In an endpoint
gauge, if two vertices have zero field, flipping both changes the
score by four times their oriented edge.  Endpoint maximality forces
that edge to be `-1`.  Applying this to the two zero-field classes in
(10.396) shows that `G` is bipartite across `s^+\mid s^-` and `H`
is bipartite across `t^+\mid t^-`.  A forbidden within-class edge
would explicitly improve one plane word by `4`.

Every `G` edge therefore has `s_is_j=-1`, and balance gives
```math
H_B(s)
=
\sum_{i<j}B_{ij}s_is_j
=-\frac m2+2e.
```
Fix the singleton spin to `+1`, set `x=s`, and average the other
shore state `y` uniformly.  The `D` quadratic, `C` cross, and
`t`-linear terms all have mean zero, while `s^{\mathsf T}x=m`.
Hence
```math
\boxed{
\mathbb E_y H_A(s,1,y)
=
\frac m2+2e
=
\frac{m^2+m-M_0}{2}.
}
\tag{10.397}
```
Some Boolean `y` attains at least this average.  Therefore actual
endpoint status `M(A)=M_0` requires
```math
\boxed{
M_0\ge\frac{m^2+m-M_0}{2},
\qquad
M_0\ge\frac{m(m+1)}{3}.
}
\tag{10.398}
```
This rules out the entire matched-degree ansatz at competitive scale,
without regularity or pseudorandomness.  In the `d`-regular case,
`M_0=m(m-2d)`, so
```math
m-2d\ge\frac{m+1}{3},
\qquad
d\le\frac{2m-1}{6}.
```
In particular, `d=(m-C\sqrt m)/2` is impossible.  Equivalently, the
necessary bipartite closure creates the exact Boolean eigenvector
`Bs=(2d-1)s`; secondary spectral pseudorandomness cannot remove it.

**Numerical audit.**  Exhaustive enumeration at
`(m,d)=(4,1),(6,1),(8,2)` gives `(M_0,M(A))` equal to
`(8,16),(24,24),(32,44)`, respectively.  This is a normalization
check only; (10.398) is exact.

#### 10.57.4 Updated frontier

This wave changes the priority order:

- the clipped triangular target (10.382) is solved on every
  compatible fixed-orientation path by one fixed Boolean mosaic.  On
  paths with resets, the only remaining scalar cost is the sum of the
  ordinary reset layers in (10.390);

- clipped mass cannot pointwise force a triangle residual: singleton
  and two-vertex shores already falsify that implication.  A residual
  exchange would need grouped cycle excess and control of the
  conditional extension profile (10.394);

- the most natural matched-degree deformation of the rank-two affine
  puncture family is necessarily quadratic.  Any competitive affine
  obstruction must evade either degree matching or the zero-field
  pair-flip closure.

The highest-priority open problem is now a reset-layer or branch
tensorization theorem: charge
`\sum_{\sigma_t=-}(2q_t-g_t)` and the reset-tree allocation without
paying the same parent range on many branches.  Independent routes
are a grouped conditional-extension theorem for (10.394), and a
global-minimality argument outside the matched-degree affine ansatz.

### 10.58 Strict-reset walls and a conserved two-channel allocation

The sixth wave falsifies the natural constant-two bound for raw reset
layers, even when orientation changes are forced and ties preserve the
old channel.  It also gives a different orientation-safe construction:
one-sided endpoint capacities plus recursively conserved imbalance
obligations.  All algebraic claims and displayed finite certificates
below are **Verified**.  Broader finite searches are labelled
**Numerical**.

#### 10.58.1 Exact reset algebra and the bounded imbalance tax

At a forced reset, let `\rho` be the old carried orientation and put
`p=P(\rho A_t)` and `q=N(\rho A_t)>p`.  If `y_t` is the inherited
spin, write
```math
h=\rho\,y_t^{\mathsf T}A_ty_t,\qquad
b=q-p,\qquad u=p-h,\qquad \mathcal R=p+q.
```
There are two different gaps, and conflating them reverses the useful
sign.  The carried old-channel gap and the fresh-orientation ordinary
cut layer are respectively
```math
\boxed{
g=q-h=b+u,\qquad
a=q+h=2q-g=\mathcal R-u.
}
\tag{10.399}
```
Thus a small dominance imbalance `b` does not make the reset layer
`a` small: `a` is the *complement* of the old-channel replenishment
`u` inside the full range.

Let the fresh endpoint have orientation `\rho'=-\rho`, disagree with
`y_t` on `D_t`, and agree on `E_t`.  Direct cut expansion gives
```math
\boxed{
a_t
=
-4\rho'\,
y_t(D_t)^{\mathsf T}A[D_t,E_t]y_t(E_t).
}
\tag{10.400}
```
All internal terms cancel.  Raw reset accumulation is therefore an
alternating triangular cross-traffic problem.

The imbalance part is nevertheless completely controlled.  Along a
nested path, `P(A_t)` and `N(A_t)` separately decrease, while the
signs of strict dominance alternate at successive resets.  If a reset
entering the negative side at `r` is followed by one entering the
positive side at `s`, then
```math
\begin{aligned}
b_r+b_s
&=(N_r-P_r)+(P_s-N_s)\\
&=(N_r-N_s)-(P_r-P_s)
\le N_r-N_s.
\end{aligned}
```
Pairing consecutive resets, and charging a possible final unpaired
term to the terminal extremum, proves
```math
\boxed{
\sum_{\mathrm{forced\ resets}}b_t\le Q(A_0).
}
\tag{10.401}
```

The same fact tensorizes on a binary partition tree.  Orient every
node by its larger one-sided extremum, preserving its parent's
orientation on a tie.  Below a positive node `U`, the sum of
strict-reset imbalances is at most `N(U)`.  Indeed, a positive child
contributes at most `N(X)` by induction, while a negative child
contributes
```math
N(X)-P(X)+(\text{subtree contribution})\le N(X).
```
Now sum over the children and use `N`-superadditivity.  The negative
case is symmetric.  Hence the full tree also obeys (10.401), with the
sharper opposite-root one-sided extremum on its right-hand side.

There is a sharp one-step relation for the unbounded complement `u`.
Suppose the preceding step had fresh endpoint
`x=(x_D,x_E)` in the current old orientation, and let `z` attain the
best old-channel score on the retained core `E`.  Parent maximality,
applied with both `z` and `-z` while holding `x_D` fixed, gives
```math
\boxed{
u_t+
2\left|
\rho\,x_D^{\mathsf T}A[D,E]z
\right|
\le\frac{\alpha_{t-1}}2,
}
\tag{10.402}
```
where `\alpha_{t-1}` is the preceding step's ordinary cut deficit.
Thus `u_t\le\alpha_{t-1}/2`.  If that preceding step did not reset,
`\alpha_{t-1}` is its carried gap; if it did reset, it is its
potentially large reset layer.  This last distinction prevents
(10.402) from closing the raw reset recursion.

#### 10.58.2 A strict order-twelve reset tower beats two ranges

The following signing is an exact counterexample to the most natural
remaining scalar conjecture:
```math
A=
\begin{pmatrix}
0&-1&-1&-1&-1&-1&-1&1&-1&1&-1&-1\\
-1&0&-1&-1&-1&-1&1&-1&1&-1&-1&-1\\
-1&-1&0&1&-1&1&1&-1&1&1&1&-1\\
-1&-1&1&0&-1&1&1&1&1&1&1&1\\
-1&-1&-1&-1&0&1&1&1&1&-1&1&1\\
-1&-1&1&1&1&0&1&-1&-1&-1&-1&-1\\
-1&1&1&1&1&1&0&-1&1&-1&-1&-1\\
1&-1&-1&1&1&-1&-1&0&-1&1&1&1\\
-1&1&1&1&1&-1&1&-1&0&1&1&1\\
1&-1&1&1&-1&-1&-1&1&1&0&1&-1\\
-1&-1&1&1&1&-1&-1&1&1&1&0&-1\\
-1&-1&-1&1&1&-1&-1&1&1&-1&-1&0
\end{pmatrix}.
\tag{10.403}
```
Use the ordered blocks `(2,3,2,2,2,1)`, start in orientation
`\rho_0=-1`, and successively retain the suffix blocks.  Exact
enumeration gives

| step | `P` | `N` | inherited `H` | `\rho` | fresh `\tau` | reset `a` |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 44 | 40 | -4 | -1 | +1 | 48 |
| 1 | 34 | 38 | 22 | +1 | -1 | 60 |
| 2 | 26 | 22 | -6 | -1 | +1 | 32 |
| 3 | 8 | 16 | 8 | +1 | -1 | 24 |
| 4 | 6 | 2 | -2 | -1 | +1 | 8 |

At every row the fresh channel strictly dominates the carried channel,
and the prescribed block-cut spin is an absolute endpoint.  No reset
is a tie choice.  The root has
```math
P(A)=44,\qquad N(A)=40,\qquad Q(A)=44,\qquad\mathcal R(A)=84,
```
but
```math
\boxed{
\sum_ta_t=172
>168=2\mathcal R(A)
>132=3Q(A).
}
\tag{10.404}
```
The corresponding soft totals are
`\sum b_t=24`, `\sum u_t=64`, and `\sum g_t=88`.  Thus the bounded
imbalance theorem does not control the complementary reset layers.

This certificate also defeats the local coefficient-two estimate:
at the first step the range falls only from `84` to `72`, while
`a_0=48>2(84-72)`.  Smaller local examples exist, but are unnecessary
for the global falsification.

**Numerical search scope.**  Exact dynamic programming exhaustively
checked all signings, inherited orientations, and compatible shores
through order six.  The maximum tie-preserving ratios
`\sum a/\mathcal R` at orders `3,4,5,6` were `1,1,4/3,3/2`.
Exact-endpoint mixed-integer searches found `13/7` at order nine and
the verified `43/21` in (10.404) at order twelve.  These data do
**not** prove unbounded growth or disprove every constant-times-`Q`
theorem.

#### 10.58.3 A conserved two-channel successor allocation

Raw exact resets are avoidable at the algebraic allocation stage.  For
a nontrivial block `U`, choose spins `p,n` attaining `P(U)` and
`-N(U)`.  Split `U=S\sqcup T` according to the two values of
`p_i n_i`.  After gauging by `p`, the negative endpoint is exactly the
cut flip.  If
```math
c=p_S^{\mathsf T}A[S,T]p_T,
```
then
```math
\boxed{\mathcal R(U)=P(U)+N(U)=4c.}
\tag{10.405}
```
On either shore `X`, the restrictions of `p` and `n` differ only by a
global sign, so they have one common internal energy
```math
h_X=p_X^{\mathsf T}A[X]p_X=n_X^{\mathsf T}A[X]n_X.
```
Let `L_X` be the cross `\ell_1` norm obtained by fixing `p_X` and
exposing the opposite shore.  Then
`L_X\ge c=\mathcal R(U)/4`.

Write `C_+(X)=P(X)`, `C_-(X)=N(X)` and define the one-sided directed
payoffs
```math
\lambda_X^\sigma
=
2L_X-\bigl(C_\sigma(X)-\sigma h_X\bigr),
\qquad \sigma\in\{+,-\}.
```
For each shore,
```math
\lambda_X^++\lambda_X^-
=4L_X-\mathcal R(X)
\ge\mathcal R(U)-\mathcal R(X).
```
Since `P` and `N` are separately block-superadditive,
`\mathcal R(S)+\mathcal R(T)\le\mathcal R(U)`.  Hence
```math
\boxed{
\sum_{X\in\{S,T\}}\sum_{\sigma=\pm}
(\lambda_X^\sigma)_+
\ge\mathcal R(U).
}
\tag{10.406}
```

The payoffs required by the all-successor theorem (10.341) must use
the absolute child benchmark.  Put
```math
\mu_X^\sigma
=
2L_X-\bigl(Q(X)-\sigma h_X\bigr),
\qquad
b_X^\sigma=Q(X)-C_\sigma(X).
```
Then `\lambda_X^\sigma=\mu_X^\sigma+b_X^\sigma` and
```math
b_X^++b_X^-
=2Q(X)-\mathcal R(X)
=|P(X)-N(X)|.
```
Using `(v+b)_+\le v_++b` for `b\ge0` in (10.406) gives the exact
corrected coverage inequality
```math
\boxed{
\sum_{X,\sigma}(\mu_X^\sigma)_+
+|P(S)-N(S)|+|P(T)-N(T)|
\ge\mathcal R(U).
}
\tag{10.407}
```
Every `\mu_X^\sigma` is a valid oriented all-successor payoff:
```math
\mu_X^\sigma
=2L_X-\delta_X(\sigma,p_X),
\qquad
\delta_X(\sigma,p_X)=Q(X)-\sigma h_X.
```

Equation (10.407) has an exact conserved recursion.  Give `U` an
obligation `0\le w_U\le\mathcal R(U)`.  Distribute it among layer
allocations and two child obligations so that
```math
w_U
=
\sum_{X,\sigma}a_X^\sigma+w_S+w_T,
```
with
```math
0\le a_X^\sigma\le(\mu_X^\sigma)_+,\qquad
0\le w_X\le|P(X)-N(X)|.
```
The nonnegative buckets in (10.407) have enough total capacity for
this choice, and `w_X\le\mathcal R(X)` permits recursion on the child.
Use its own positive/negative endpoint cut.  Every nontrivial block
strictly splits; singleton imbalance is zero.  Summing the
conservation identities over the finite tree yields
```math
\boxed{
\sum_{\text{all chosen layers}}a_X^\sigma
=w_{\rm root}
\le\mathcal R(A)
\le2Q(A).
}
\tag{10.408}
```
If zero-obligation branches are completed by the same endpoint cuts,
every original edge lies in the cross block of exactly one lowest
common ancestor in the resulting full two-channel tree.

This is a genuine orientation-safe allocation theorem, but its scope
must remain explicit.  It bounds the **chosen fractional
allocations**, not the sum of all positive raw capacities and not the
exact-reset-tree sum falsified by (10.404).  It does not follow by
applying (10.341) independently to every branch.  The next step is a
sampling or purification theorem which harvests the conserved chosen
layers without reintroducing branch overcount.

#### 10.58.4 The imbalance obligations cannot be omitted

Immediate four-layer coverage fails if one replaces both one-sided
benchmarks by `Q`.  The smallest audited obstruction is
```math
A=
\begin{pmatrix}
0&-1&1&-1&-1&-1\\
-1&0&1&-1&1&1\\
1&1&0&-1&1&-1\\
-1&-1&-1&0&1&-1\\
-1&1&1&1&0&-1\\
-1&1&-1&-1&-1&0
\end{pmatrix}.
\tag{10.409}
```
Exact enumeration gives `P(A)=N(A)=10`.  Take endpoints
```math
p=(-1,-1,-1,1,-1,1),\qquad
n=(-1,-1,1,-1,1,1)
```
of energies `10` and `-10`.  Their endpoint cut is
```math
S=\{2,3,4\},\qquad T=\{0,1,5\}.
```
In the `p` gauge its cross block is
```math
B=
\begin{pmatrix}
1&1&1\\
1&1&-1\\
-1&1&1
\end{pmatrix},
```
and
```math
(P_S,N_S,P_T,N_T,h_S,h_T,L_S,L_T)
=(2,6,6,2,2,-2,5,5).
```
The four one-sided payoffs in (10.406) are `(10,2,2,10)`, but the
four `Q`-benchmarked payoffs are only `(6,2,2,6)`.  Therefore
```math
\boxed{
\sum_{X,\sigma}(\mu_X^\sigma)_+
=16<20=\mathcal R(A).
}
\tag{10.410}
```
Both child imbalances equal `4`, and they are precisely the recursive
obligations missing from the failed shortcut.

#### 10.58.5 External warning and updated frontier

The reset mosaic is triangular, so a generic matrix-norm shortcut
would be especially tempting.  It cannot be the whole proof:
[Mishura's triangular-truncation theorem
(arXiv:2110.06984)](https://arxiv.org/abs/2110.06984) proves that
triangular truncation on symmetric matrices has operator norm growing
without bound for the cut norm.  This theorem does **not** include the
endpoint-dominance constraints here, so it neither proves nor
disproves a reset bound.  It does verify that those extra constraints
must be used in any constant argument.

The frontier is now:

- dominance imbalance is harmless: its total is at most one opposite
  root extremum on a path or a full partition tree;

- raw ordinary reset layers are not harmless.  Even strict,
  tie-preserving resets beat `2(P+N)`, and no universal constant is
  currently proved or disproved;

- the conserved two-channel construction (10.408) replaces raw reset
  capacity by selected valid successor layers of total mass at most
  `2Q`.  Its remaining gap is no longer orientation algebra but
  harvesting/purification of those fractional layers in the original
  insertion or peeling recurrence.

The highest-priority next target is therefore a conserved-layer
harvesting theorem: turn the allocations in (10.408) into actual
successor choices, or a randomized choice with the correct expected
gain, while preserving the single global `3Q` capacity bound from
(10.341).  Independent routes are to use global minimality to exclude
long strict-reset mosaics such as (10.403), or to prove an
endpoint-constrained triangular theorem stronger than the generic
cut-norm setting.

### 10.59 Path congestion, branching rounding, and a half-factor no-go

The seventh wave identifies the exact difference between the
conserved allocation (10.408) and an actual successor chain.  The
difference is a normalized Carleson congestion, not integrality of
the allocation flow.  Lossless chain purification is false in a
genuine fixed-orientation recurrence, and even a universal factor
`1/2` for the most liberal chain is false.  On the positive side,
opposite directed layers admit an exact factor-two *branching*
rounding.  All algebraic statements and displayed finite witnesses
below are **Verified**.  Exhaustive and broader tree searches are
labelled **Numerical**.

#### 10.59.1 The exact path-cover norm

Fix one finite endpoint partition tree.  At node `v`, for child
`X` and orientation `\sigma`, write
```math
c_{vX}^{\sigma}=(\mu_{vX}^{\sigma})_+,
\qquad 0\le a_{vX}^{\sigma}\le c_{vX}^{\sigma}
```
for the available capacity and its chosen allocation.  An actual
successor chain chooses one pair `(X,\sigma)` at every visited node
and then continues only inside `X`.  Define the directed-edge load
```math
\ell_{vX}
=
\sum_{\sigma=\pm}\frac{a_{vX}^{\sigma}}{c_{vX}^{\sigma}},
\tag{10.411}
```
with a zero summand when `a=c=0` and value `+\infty` when
`a>0=c`.

The minimum total mass `K` of a finite measure on oriented
root-to-leaf chains whose integrated reward dominates every
allocation coordinate has two exact descriptions:
```math
\boxed{
\begin{aligned}
\theta_{\rm leaf}&=0,\\
\theta_v&=
\sum_{X\text{ child of }v}
\max\{\ell_{vX},\theta_X\},\\
K&=\theta_{\rm root}
=
\sup_{\mathcal F\text{ edge antichain}}
\sum_{e\in\mathcal F}\ell_e.
\end{aligned}
}
\tag{10.412}
```
Necessity of the antichain bound is immediate: one chain meets at
most one edge of an edge antichain.  For sufficiency, recursively
send `\max\{\ell_{vX},\theta_X\}` units of chain mass into `X`,
split enough of that mass between its two orientation labels to meet
the two lower bounds in (10.411), and pass the remaining mass through
the already constructed child flow.  This also proves the recursive
formula.

Thus `K\le1` is precisely the condition under which the allocations
can be averaged through (10.341) once.  In general, integrating
(10.341) against a chain measure of mass `K` gives only the multiplied
bound `3KQ(A)`.  Both the allocation polytope and the path polytope
have node--arc incidence matrices and are totally unimodular, but
they are different polytopes: for integral right-hand sides and
capacity bounds, integrality of (10.408) gives integral *energy units
which may branch*, not one
integral successor chain.

The scalar inequalities in (10.407)--(10.408) alone cannot control
`K`.  In an abstract depth-`d` binary system, internal nodes have
depths `0,\ldots,d-1`; give every depth-`j` internal node range and
obligation `2^{-j}`.  At depth `j<d-1`, give each child imbalance
`2^{-(j+1)}` and use zero local layer capacity.  At depth `d-1`,
give the two singleton children zero range and imbalance, and allocate
`2^{-d}` into one orientation layer on each of the two outgoing
edges.  Every coverage and conservation identity is exact.  There are
`2^d` bottom edges, so total allocated mass is one, but every path
receives only `2^{-d}`.  The bottom-edge antichain has unit normalized
load on every edge, and (10.412) gives `K=2^d`.  This is an abstract
capacity system, **not** a matrix realization.  It proves that any
constant theorem must use endpoint-matrix geometry.

#### 10.59.2 A genuine order-five recurrence defeats one path

Let `C=J_3-I_3`.  In doubled normalization,
```math
P(C)=6,\qquad N(C)=2,\qquad\mathcal R(C)=8.
```
Take positive endpoint `p=(1,1,1)` and negative endpoint
`n=(-1,-1,1)`.  Their cut has a two-vertex shore and a singleton.
The two oriented capacities toward the pair are `(4,0)`, and those
toward the singleton are `(4,4)`.  Both child imbalances are zero.
Consequently every conserved allocation of obligation eight is
entirely at the root, all its nonzero capacities equal four, and
```math
\boxed{K=2.}
\tag{10.413}
```

Even if unallocated descendant capacities are allowed for free, one
chain cannot recover the range.  Define the liberal Bellman value
```math
V(U)=
\max_{p,n,X,\sigma}
\left\{(\mu_X^\sigma)_++V(X)\right\},
\qquad V(\text{singleton})=0,
\tag{10.414}
```
where every child may optimize all fresh endpoint ties.  The pair
branch earns `4+2`, while the singleton branch earns `4`, so
```math
V(C)=6<8=\mathcal R(C).
```
Randomization cannot improve `V`, since expected path reward is a
linear average of deterministic path rewards.

This is not merely an abstract tree or an orientation tie.  Attach a
positive two-clique and use `n` as both cross rows:
```math
B=
\begin{pmatrix}
0&1&-1&-1&1\\
1&0&-1&-1&1\\
-1&-1&0&1&1\\
-1&-1&1&0&1\\
1&1&1&1&0
\end{pmatrix}.
\tag{10.415}
```
Exact enumeration gives `P(B)=N(B)=12`, and
`(1,1,-1,-1,1)` is a positive endpoint.  Equivalently, if `m` is
the two-clique spin sum, `s` the core spin sum, and `t=n\cdot y`,
every energy is
```math
m^2+s^2-5+2mt,
```
whose absolute value is at most twelve on the finite feasible set.
Deleting the two-clique leaves inherited core energy `-2`, while the
strictly dominant positive core endpoint has energy `6`.  The
same-orientation successor gap is exactly eight, but every liberal
chain has value only six.

#### 10.59.3 An exact family below the one-half path factor

The failure is asymptotic, not only finite.  Let
`m=k^2=4^r`, choose `s\in\{\pm1\}^m` with
`\mathbf1^{\mathsf T}s=k`, and put
```math
C_m=ss^{\mathsf T}-I_m,
\qquad
A_m=
\begin{pmatrix}
C_m&J_m\\
J_m&-C_m
\end{pmatrix}.
\tag{10.416}
```
This is a complete signing of order `2m`.  For a state `(x,y)`, set
```math
u=\mathbf1^{\mathsf T}x,\quad a=s^{\mathsf T}x,\quad
v=\mathbf1^{\mathsf T}y,\quad b=s^{\mathsf T}y.
```
Its energy is
```math
x^{\mathsf T}C_mx-y^{\mathsf T}C_my
+2(\mathbf1^{\mathsf T}x)(\mathbf1^{\mathsf T}y)
=a^2-b^2+2uv.
```
Simultaneous negation of `x,y` permits `u\ge0`.  The Boolean box
constraints give the exact fixed-`u` optimizations
```math
\max a^2=[m-|u-k|]^2,
\qquad
\max_y(-b^2+2uv)=
\begin{cases}
u^2+2u(m-k),&0\le u\le k,\\
-m+2mu,&k\le u\le m.
\end{cases}
\tag{10.417}
```
On the first interval their sum is increasing, and on the second it
is convex.  It is therefore enough to compare `u=k` and `u=m`.
The corresponding values are
```math
m^2+2mk-m
\quad\text{and}\quad
2m^2,
```
whose difference is `m(k-1)^2`.  For `k>1`, the positive endpoint is
therefore uniquely `(\mathbf1,\mathbf1)` modulo global sign.
Moreover the map `(x,y)\mapsto(y,-x)` negates every energy, so
```math
\boxed{P(A_m)=N(A_m)=2m^2,\qquad\mathcal R(A_m)=4m^2.}
\tag{10.418}
```

The unique endpoint cut is the displayed two-shore split.  On each
shore, `h=0`, `L=m^2`, and `Q(C_m)=m(m-1)`, so both root orientation
capacities equal `m^2+m`.  Switching and sign reversal reduce either
child Bellman value to that of the positive clique `K_m`.  For powers
of two its endpoint cut is balanced and
```math
V(K_m)=\frac{m^2}{2}+V(K_{m/2})
=\frac{2(m^2-1)}3.
```
Consequently
```math
\boxed{
V(A_m)=\frac{5m^2+3m-2}{3},
\qquad
\frac{V(A_m)}{\mathcal R(A_m)}
\longrightarrow\frac5{12}.
}
\tag{10.419}
```
At `m=4`, this is an order-eight signing with
`P=N=32` and `V=30<32=\mathcal R/2`.  Thus even the universal
one-half chain conjecture is **Falsified**.  This family does not tend
to zero; whether any positive universal chain factor exists remains
open.

#### 10.59.4 Orientation rounding is local, branching is global

There is a separate obstruction before any path choice.  Let
```math
D=
\begin{pmatrix}
0&1&-1&1&-1&-1\\
1&0&1&-1&-1&-1\\
-1&1&0&-1&-1&-1\\
1&-1&-1&0&-1&-1\\
-1&-1&-1&-1&0&-1\\
-1&-1&-1&-1&-1&0
\end{pmatrix}.
\tag{10.420}
```
Exact enumeration gives unique projective endpoints,
```math
P(D)=14,\qquad N(D)=18,\qquad\mathcal R(D)=32.
```
Their two shores both have zero imbalance.  The two oriented
`Q`-benchmarked capacities are `(12,16)` on one shore and `(8,8)`
on the other.  Hence fractional capacity is `44`, but selecting at
most one orientation on each directed shore gives only
```math
\boxed{16+8=24<32.}
\tag{10.421}
```
There is nothing to pass recursively.  This too occurs in an actual
strict-channel recurrence.  If
`p=(-1,-1,-1,-1,1,1)` is the positive endpoint, attach a negative
five-clique and use `-p` for every cross row.  The resulting order-
eleven parent has `(P,N)=(54,66)`, and `(\mathbf1_5,p)` has energy
`-66`.  After deleting the five-clique, the inherited negative-
oriented score is `-14`, whereas the fresh negative core endpoint has
score `18`; the exact carried gap is `32`.

There is nevertheless a sharp local branching repair.  Gauge a
positive endpoint of a node `U` to `\mathbf1`, so its negative
endpoint is the cut `U=S\sqcup T`.  Put
```math
c=\sum_{i\in S,j\in T}a_{ij}=\frac{\mathcal R(U)}4.
```
For `W\subseteq S`, let
```math
k_W=\sum_{i\in W,j\in S\setminus W}a_{ij},
\qquad
e_W=\sum_{i\in W,j\in T}a_{ij}.
```
Flipping `W` at the positive endpoint and at the negative endpoint
gives respectively `k_W+e_W\ge0` and `k_W-e_W\le0`.
Applying the same two inequalities to `S\setminus W` yields
```math
\boxed{
e_W\ge|k_W|,\qquad c-e_W\ge|k_W|,\qquad
c\ge2|k_W|.
}
\tag{10.422}
```
Also `L_S\ge c`.  Since every child spin is a cut flip and has
energy `h_S-4k_W`, (10.422) proves
```math
P(S)-h_S\le2L_S,qquad
N(S)+h_S\le2L_S.
```
Thus both one-sided payoffs are nonnegative:
```math
\boxed{
\lambda_S^+=2L_S-(P(S)-h_S)\ge0,qquad
\lambda_S^-=2L_S-(N(S)+h_S)\ge0,
}
\tag{10.423}
```
and likewise on `T`.  The `(\lambda_X^\sigma)_+` positive-part
operations in (10.406) were therefore unnecessary; this does not
remove the positive parts from the `\mu`-payoffs.

There is also a clipping-safe opposite-pair bound.  With
`C_+=P,C_-=N`, exposing the cross field and bounding the opposite
shore by its opposite one-sided extremum gives
```math
\lambda_X^\sigma
\le
C_\sigma(U)+C_{-\sigma}(Y)-C_\sigma(X),
\qquad \{X,Y\}=\{S,T\}.
```
Adding the two opposite directions, using
`(\mu_X^\sigma)_+\le\lambda_X^\sigma`, yields
```math
\boxed{
(\mu_S^\sigma)_+
+(\mu_T^{-\sigma})_+
\le\mathcal R(U).
}
\tag{10.424}
```

Now group the four allocations at every node as
`\{(S,+),(T,-)\}` and `\{(S,-),(T,+)\}`.  Keep the group carrying
more allocated mass and recurse into both children.  Conservation
and induction give the exact branching guarantee
```math
\boxed{
\text{retained allocated mass}
\ge\frac12\,w_{\rm root}.
}
\tag{10.425}
```
The selected capacity at node `U` is at most `\mathcal R(U)` by
(10.424).  This is not one actual successor tower: it selects two
directed witnesses at every lowest common ancestor, and their
exposing spins need not be nested-compatible.  Thus (10.341) cannot
be applied once.  The immediate estimate on the sum of selected
available capacities is only
```math
\sum_{U\text{ internal}}\mathcal R(U).
```
Separate `P/N` superadditivity bounds the range sum on every block
antichain by the root range, but it does not bound the sum over all
depths.  LCA uniqueness prevents repetition of original cross-block
edges; it does not supply this missing depth estimate.

#### 10.59.5 Exhaustive scope and updated frontier

**Numerical exhaustive audit.**  Every labelled signing, every
endpoint tie, and every recursive child choice were enumerated
through order six.  Among strict-dominance matrices, the minimum
liberal ratios `V/\mathcal R` at orders `3,4,5,6` are respectively
```math
\frac34,\qquad\frac58,\qquad\frac34,\qquad\frac59.
```
Every signing through order five has an endpoint choice for which
one orientation per shore plus the two child imbalances covers the
root range.  Matrix (10.420) is the first failure at order six, with
ratio `24/32`.  These finite statements were independently rerun in
exact integer arithmetic.

For a numerical upper-envelope audit, let `T(U)` optimize endpoint
ties independently at every node, select the better `\mu` orientation
independently on each shore, add both selected capacities, and recurse
into both children.  This is more liberal than the opposite-pair
selection in (10.425).  The exact computations found `T/Q=2` for
positive and negative cliques through order twelve, `4.4` for the
matrix in (10.409), and `188/44` for the strict-reset matrix (10.403).
These data neither prove nor disprove a universal branching constant.

The frontier is now exact:

- lossless path purification and even a one-half path factor are
  false.  No realizable matrix family with path factor tending to
  zero is known;

- for coordinatewise path domination of a fixed allocation, the
  missing hypothesis is precisely the antichain/Carleson norm
  (10.412).  Conservation controls raw energy mass, not normalized
  path congestion;

- opposite orientations admit the factor-two branching theorem
  (10.425).  LCA uniqueness prevents cross-edge repetition, but range
  superadditivity controls only one block antichain at a time.  A
  depth-summed endpoint-matrix Carleson estimate or nested exposing-
  spin theorem is still missing.

The highest-priority target is therefore to prove `K=O(1)` in
(10.412) for a suitably chosen conserved allocation arising from
actual endpoint matrices, or to construct a realizable family with
`K\to\infty`.  An independent route is to use global minimality to
exclude the recurrence embeddings above or to force compatibility
between the exposing spins on different levels.  Here `K=O(1)` would
give constant-loss path control; the literal single-probability,
single-`3Q` use of (10.341) requires `K\le1`.

### 10.60 Minimax congestion and an exact endpoint stopping target

The eighth wave jointly optimizes the conserved allocation and the
endpoint tree, rather than testing an arbitrary allocation.  It finds
two stronger realizable congestion obstructions, culminating in an
exact order-eight value `10/3`.  It also sharpens the endpoint geometry:
both cross-exposure norms are exactly one quarter of the parent range.
This gives an unconditional half-imbalance stopping theorem and reduces
a universal `K\le4` theorem to one explicit local inequality.  The local
inequality remains **Open**.  All LP values called exact below have
matching rational primal and dual certificates; finite searches without
such certificates are labelled **Numerical**.

#### 10.60.1 The joint allocation--path minimax LP

Fix an endpoint tree, let `c_{vX}^{\sigma}>0` denote its nonzero layer
capacities, and put
`d_X=|P(X)-N(X)|`.  For root obligation `W`, choose allocations and
obligations satisfying
```math
\boxed{
w_v=\sum_{X,\sigma}a_{vX}^{\sigma}+\sum_Xw_X,
\qquad
w_{\rm root}=W,
\qquad
0\le a_{vX}^{\sigma}\le c_{vX}^{\sigma},
\quad 0\le w_X\le d_X.
}
\tag{10.426}
```
After eliminating the obligations, this says that the total allocated
mass is `W` and that the mass allocated in every proper subtree `X` is
at most `d_X`.  Minimize the path-cover mass in (10.412) over this
polytope.  Denote the result by `K_*(\mathcal T,W)`, and also minimize
over every endpoint tie when the tree is not prescribed.

There is a useful exact dual certificate.  Let nonnegative edge weights
`q_e` lie in the fractional edge-antichain polytope
```math
q_e\ge0,
\qquad
\sum_{e\in\pi}q_e\le1
\quad
\text{for every root-to-leaf path }\pi.
\tag{10.427}
```
For unrestricted node potentials `y_v\in\mathbb R`, bucket slacks
`z_{vX\sigma}\ge0`, and obligation-cap slacks `s_X\ge0`, the dual is
```math
\boxed{
\begin{aligned}
\text{maximize}\quad&
W y_{\rm root}
-\sum_{v,X,\sigma}c_{vX}^{\sigma}z_{vX\sigma}
-\sum_{X\ne\mathrm{root}}d_Xs_X,\\
\text{subject to}\quad&
y_v-z_{vX\sigma}\le\frac{q_{vX}}{c_{vX}^{\sigma}},\\
&y_v-y_X-s_X\le0.
\end{aligned}
}
\tag{10.428}
```
This follows by exchanging the allocation minimum with the antichain
maximum in (10.412) and dualizing the resulting tree flow.  It is also a
convenient way to certify a lower bound without trusting a floating-point
LP value.

The following values are exact.

| matrix and tree convention | root obligation | `K_*` |
|---|---:|---:|
| positive triangle | `8=\mathcal R` | `2` |
| matrix (10.420), unique tree | `32=\mathcal R` | `5/2` |
| family (10.416), `m=4` | `64=\mathcal R` | `12/5` |
| family (10.416), `m\ge16` | `4m^2` | `\displaystyle\frac{12m^2}{5m^2+3m-8}` |
| reset matrix (10.403), lexicographic tree | `84=\mathcal R` | `22/7` |
| reset matrix, all endpoint ties | `84` | `28/13` |
| reset matrix, lexicographic tree | carried gap `40` | `58/51` |
| reset matrix, all endpoint ties | `40` | `19/23` |
| prescribed fresh reset endpoint, later ties optimized | `40` | `104/119` |
| same prescribed endpoint | `84` | `22/9` |

For the reset matrix, all twelve positive and fifteen negative
projective root endpoints and all recursive ties reduce to `250`
LP-distinct full-tree signatures.  Twelve attain `28/13`, and every
signature has an exact dual certificate at least that large.  Fixing the
fresh positive endpoint from the displayed reset chain leaves `43`
distinct signatures.  Thus the large raw reset total in (10.404) does
not itself force large optimized path congestion for its actual first
carried gap.

One optimized range allocation for (10.403) has root branch data
```math
((14,34),d=4),
\qquad
((26,42),d=0).
\tag{10.429}
```
Pass four units into the first child and allocate root masses `34,42,4`.
The two root edge maxima are `1` and `15/13`; the first child congestion
is smaller.  This gives the primal value `28/13`.

The family (10.416) also does not produce growing minimized congestion.
For a clique block of size `s\ge4`, the largest usable layer on either
child edge is `s^2/2`, its imbalance is `s(s-2)`, and the size-two
imbalance vanishes.  The total capacity usable along one canonical
chain is therefore
```math
\boxed{
B_m
=m^2+m+\frac12\sum_{s=m,m/2,\ldots,4}s^2
=\frac{5m^2+3m-8}{3}.
}
\tag{10.430}
```
Assigning antichain weight proportional to the largest capacity on each
edge gives the matching lower certificate.  A symmetric chain measure
attains it for `m\ge16`; the imbalance caps are respected because
`4m(s+2)\le5m^2+3m-8` at every encountered size.  Hence
```math
\boxed{
K_*(A_m)
=\frac{12m^2}{5m^2+3m-8}
\longrightarrow\frac{12}{5}.
}
\tag{10.431}
```
At `m=4`, an imbalance cap binds and the exact value is separately
`12/5`.

#### 10.60.2 Exact exposure and unconditional half-imbalance stopping

Return to one endpoint node `U`.  Gauge its positive endpoint to
`\mathbf1`, let the negative endpoint be the cut `U=S\sqcup T`, and
write
```math
c=\sum_{i\in S,j\in T}a_{ij}=\frac{\mathcal R(U)}4.
```
For `i\in S`, let `r_i` be its internal row sum in `S` and `e_i` its
cross row sum into `T`.  Flipping `i` at the positive endpoint gives
`r_i+e_i\ge0`; flipping it at the negative endpoint gives
`r_i-e_i\le0`.  Thus
```math
e_i\ge|r_i|\ge0.
```
Every cross row sum is nonnegative, and the symmetric argument in `T`
makes every cross column sum nonnegative.  Their totals are `c`, so the
inequality `L_X\ge c` used previously is in fact an equality:
```math
\boxed{
L_S=L_T=c=\frac{\mathcal R(U)}4.
}
\tag{10.432}
```

Put `H=h_S+h_T`, so
```math
P(U)=2c+H,
\qquad
N(U)=2c-H,
\qquad
I(U):=|P(U)-N(U)|=2|H|.
```
For a shore `X`, let its best single-orientation capacity be
```math
b_X:=\max_{\sigma=\pm}(\mu_X^\sigma)_+.
```
The orientation attaining `Q(X)` has `\mu=\lambda\ge0` by (10.423), so
the positive part does not alter the maximum.  Equation (10.432) gives
the exact formula
```math
\boxed{
b_X=2c-Q(X)+|h_X|.
}
\tag{10.433}
```

This already proves a nontrivial stopping theorem.  Since
`Q(X)\le\mathcal R(X)` and separate `P/N` superadditivity gives
`\mathcal R(S)+\mathcal R(T)\le4c`,
```math
\boxed{
b_S+b_T
=4c-Q(S)-Q(T)+|h_S|+|h_T|
\ge |h_S|+|h_T|
\ge |H|
=\frac{I(U)}2.
}
\tag{10.434}
```
Thus at least half of every inherited imbalance is always available in
the two best local orientation buckets.  This statement is
**Verified** and uses no global-minimality assumption.

#### 10.60.3 The precise missing local inequality and conditional `K\le4`

The finite evidence supports the exact factor-two upgrade
```math
\boxed{
\mathcal B(U):=b_S+b_T\ge I(U).
}
\tag{10.435}
```
By (10.433), this is equivalent to
```math
\boxed{
\sum_{X=S,T}\bigl(Q(X)-|h_X|\bigr)
\le 2\min\{P(U),N(U)\}
=4c-2|H|.
}
\tag{10.436}
```
Both (10.435) and (10.436) are **Open targets**, not proved claims.  A
stronger sufficient inequality, also supported by all current tests, is
```math
\boxed{
Q(S)+Q(T)\le4c-|H|.
}
\tag{10.437}
```
It is important that (10.437) is stronger, not merely a rewriting of
(10.436), when `h_S` and `h_T` have opposite signs.

There is a complete conditional consequence.  Assume (10.435) at every
endpoint node.  At the root put
```math
C=\sum_{X,\sigma}(\mu_X^\sigma)_+,
\qquad
B=C+I(S)+I(T),
\qquad
t=\frac{\mathcal R(A)}B\le1,
\tag{10.438}
```
where `B\ge\mathcal R(A)` is (10.407).  Allocate `t` times every root
capacity and pass obligations `w_X=tI(X)`.  At a child `X`, use only the
best orientation on each of its two shores.  If `I(X)=0`, set `s_X=0`
and allocate or pass nothing.  Otherwise allocate `s_Xb_Y` there, where
`s_X=w_X/\mathcal B(X)\le t`, and pass nothing farther.

Each root directed-edge load is at most `2t`, while every child subtree
has `\theta_X\le2s_X\le2t`.  The recursion (10.412) therefore gives
```math
\boxed{
K\le
\sum_{X=S,T}\max\{2t,2t\}
\le4.
}
\tag{10.439}
```
More directly, each of the two root summands is at most `2t`; no depth
sum or LCA estimate occurs.  Hence (10.435) would produce a conserved
allocation with universal path-cover mass four and would turn (10.341)
into a `12Q(A)` bound.

There is a sharp cut formulation of the unresolved hard case.  In the
positive-endpoint gauge, use undoubled edge weights.  Let `a_X` be the
total internal edge weight of `X`, and let `m_X,M_X` be its minimum and
maximum internal cut weights.  If `a=a_S+a_T\ge0` and both children
are positive-dominant, positive dominance says
`M_X\le a_X-m_X`, and (10.437) becomes exactly
```math
\boxed{
a_S+a_T-(m_S+m_T)\le c.
}
\tag{10.440}
```
For cuts `W\subseteq S`, `Z\subseteq T`, put
`k=k_S(W)+k_T(Z)` and let `q(W,Z)` be the separated cross weight.
The full parent endpoint conditions are
```math
0\le k+q\le c,
\qquad
0\le k+c-q\le c,
\qquad
|k|\le q\le c-|k|.
\tag{10.441}
```
One-shore use of (10.422) loses information; the simultaneous
constraints in (10.441) are essential.  For child maximum-cut spins
`x,y`, define the undoubled cross bilinear
`z=\sum_{i\in S,j\in T}a_{ij}x_i y_j`.  Their favorable relative sign
has parent energy
`2a-4(M_S+M_T)-2|z|`.  It is therefore sufficient to choose `x,y`
so that
```math
\boxed{
2(M_S+M_T)+|z|
\ge a_S+a_T-(m_S+m_T).
}
\tag{10.442}
```
Their favorable relative sign would then give a parent energy at most
`2(m_S+m_T)`, proving (10.440).  No valid averaging or cut-cone proof of
(10.440) or (10.442) is currently known.

**Numerical audit.**  The identities (10.432) and the two candidate
inequalities were checked for every switching class and every endpoint
pair through order seven.  At order seven this is `32,768` matrices and
`288,064` endpoint pairs; the minimum doubled-normalization slack in
both (10.435) and (10.437) is twelve.  An independent order-six audit
covered `1,024` switching classes and `2,242` endpoint pairs.  A further
`20,000` random order-eight draws contributed `41,278` endpoint
pairs without a failure.  Weighted cut-cone LPs also verify (10.437)
through shore size `4+4`.  These counts are evidence, not a proof.

#### 10.60.4 Exact optimized congestion `8/3` and `10/3`

The local obstruction can be stronger than (10.420), even after every
endpoint tie and allocation are jointly optimized.  At order seven,
let
```math
A_7=
\begin{pmatrix}
0&1&1&1&1&1&1\\
1&0&1&-1&1&1&-1\\
1&1&0&1&1&1&-1\\
1&-1&1&0&1&1&1\\
1&1&1&1&0&-1&1\\
1&1&1&1&-1&0&-1\\
1&-1&-1&1&1&-1&0
\end{pmatrix}.
\tag{10.443}
```
Exact enumeration gives
`(P,N,\mathcal R)=(22,18,40)`.  There are two projective endpoints of
each sign.  All four endpoint pairs have one shore with capacities
`(16,12)`, imbalance four, and child `(P,N)=(2,6)`; the other shore has
capacities `(12,12)`, zero imbalance, and child `P=N=8`.

At most four units can pass down the first branch, so at least `36`
units are root-local.  If the two edge maxima in (10.412) are `q_1,q_2`,
the first edge carries at most `12q_1+4` local units and the second at
most `12q_2`.  Hence
```math
36\le12(q_1+q_2)+4=12K+4,
\qquad
K\ge\frac83.
```
Equality passes four units to the first child, fills its capacity
sixteen, and allocates twenty units on the second edge.  The child
congestion is one.  Therefore
```math
\boxed{K_{\min}(A_7)=\frac83.}
\tag{10.444}
```

There is an even cleaner order-eight witness:
```math
A_8=
\begin{pmatrix}
0&1&1&1&1&1&1&1\\
1&0&1&-1&1&1&-1&-1\\
1&1&0&1&-1&1&-1&-1\\
1&-1&1&0&-1&-1&-1&1\\
1&1&-1&-1&0&-1&1&-1\\
1&1&1&-1&-1&0&1&1\\
1&-1&-1&-1&1&1&0&1\\
1&-1&-1&1&-1&1&1&0
\end{pmatrix}.
\tag{10.445}
```
Here `P=N=20` and `\mathcal R=40`.  There are four projective
endpoints of each sign.  Every one of the sixteen endpoint pairs splits
into two four-vertex shores, and every directed shore has
```math
(c^+,c^-)=(12,12),
\qquad
|P_X-N_X|=0,
\qquad
P_X=N_X=8,
\qquad
h_X=0,
\qquad
L_X=10.
```
No mass can pass recursively.  Every allocation coordinate therefore
has denominator twelve in (10.411), and every feasible allocation has
```math
\boxed{
K_{\min}(A_8)=\frac{40}{12}=\frac{10}{3}.
}
\tag{10.446}
```
These lower bounds are elementary integer certificates and do not rely
on solver tolerances.  In particular, every universal optimized-tree
theorem must allow `K\ge10/3`.

#### 10.60.5 Search scope and updated frontier

The exact optimized values on named finite examples are

| matrix | `K_{\min}` |
|---|---:|
| (10.409) | `2` |
| (10.420) | `5/2` |
| reset matrix (10.403) | `28/13` |
| (10.416), `m=4` | `12/5` |
| (10.443) | `8/3` |
| (10.445) | `10/3` |

The switching-gauge search is exhaustive through order six: its maximum
is two at order five and `5/2` at order six.  Among random draws, the
largest solver values found were `8/3` in `1,000` order-seven draws,
`10/3` in `11,000` order-eight draws, `8/3` in `500` order-nine
draws, and `26/9` in `500` order-ten draws.  No search found
`K>10/3`, `K>4`, or a realizable family with growing congestion.

A second structured attempt used
```math
C=ss^{\mathsf T}-I,
\qquad
B=J-I-\operatorname{diag}(s),
\qquad
\widetilde A_m=
\begin{pmatrix}C&B\\B^{\mathsf T}&-C\end{pmatrix},
\qquad
\mathbf1^{\mathsf T}s=\sqrt m.
\tag{10.447}
```
At `m=4` it has `P=N=20` and exact optimized `K=3`.  At `m=9`
(order eighteen), exact endpoint enumeration gives `P=N=138` and
unique projective endpoints.  A 600-second MILP stopped with the
following nonmatching **Numerical solver bounds**:
```math
\frac{74}{33}\lesssim K_{\min}\lesssim\frac{130}{57}.
```
No exact order-eighteen congestion claim is made.  This route currently
offers no evidence of growth.

The frontier is now narrower:

- optimized endpoint choice does not remove congestion: the exact
  realizable lower bound is `10/3`; bounds such as two or `5/2` are
  false;

- the rank-two family (10.416) tends to `12/5`, and the severe raw reset
  tower becomes mild after joint endpoint/allocation optimization.  No
  growing realizable family is known;

- endpoint geometry gives the exact identity (10.432) and the verified
  half-imbalance capacity (10.434).  The sole local upgrade needed for
  the depth-two construction is (10.435), equivalently (10.436);

- proving (10.435) would give `K\le4` and a `12Q` harvesting bound.
  Falsifying it should target the simultaneous cut system (10.441), not
  one shore in isolation.  A counterexample with `K>4`, or a family
  with `K\to\infty`, remains the independent alternative.

The highest-priority next target is therefore the endpoint stopping
inequality (10.435), with the hard same-positive reduction (10.440).
Independent routes are a conic/averaging proof of (10.442), an exact
MILP search for a violation of (10.435), and a structured construction
that drives the local four-bucket obstruction toward or past four.

### 10.61 The stopping inequality is false: a strict obstruction and a sign blow-up

The ninth wave falsifies the local stopping route from Section 10.60.
There is a strict rational weighted counterexample on nine vertices, and
a biased random clone construction lifts it to genuine complete
`\{\pm1\}` sign matrices.  Thus (10.435)--(10.437), the hard reduction
(10.440), and the universal proposed selection statement (10.442) are
all **Falsified**.  The unconditional half-imbalance bound (10.434)
remains **Verified**.  This result invalidates the particular depth-two
proof of `K\le4`; it does not disprove `K\le4` itself.

#### 10.61.1 An exact strict weighted certificate

Use the undoubled cut normalization.  Let
`S=\{0,1,2,3\}`, `T=\{4,5,6,7,8\}`, and `B=\{4,5,6\}`.  Put
`w_{ij}=M_{ij}/672`, where

```math
M=
\begin{pmatrix}
0&54&54&54&28&28&28&55&29\\
54&0&54&54&28&28&28&55&29\\
54&54&0&54&28&28&28&55&29\\
54&54&54&0&28&28&28&55&29\\
28&28&28&28&0&83&83&-55&-55\\
28&28&28&28&83&0&83&-55&-55\\
28&28&28&28&83&83&0&-55&-55\\
55&55&55&55&-55&-55&-55&0&113\\
29&29&29&29&-55&-55&-55&113&0
\end{pmatrix}.
\tag{10.448}
```

For `R\subseteq[9]`, write
`g(R)=\sum_{i\in R,j\notin R}w_{ij}`.  Exact enumeration of the 256
projective cuts gives

```math
g(\varnothing)=0,
\qquad
g(S)=g(T)=1,
\qquad
\frac1{112}\le g(R)\le\frac{111}{112}
\tag{10.449}
```

for every other projective cut.  Hence the constant spin and the
`S\mid T` spin are the unique projective parent endpoints, with `c=1`.
The strict gaps in (10.449) are important for the blow-up below.

Let `a_X` be the total internal edge weight of a shore, and let
`m_X,M_X` be its minimum and maximum internal cut weights.  A second
exact enumeration gives

| shore | `a_X` | `m_X` | `M_X` | `a_X-m_X-M_X` |
|---|---:|---:|---:|---:|
| `S` | `27/56` | `0` | `9/28` | `9/56` |
| `T` | `1/21` | `-55/112` | `19/112` | `31/84` |

Both last-column entries are positive, so both children are strictly
positive-dominant.  Nevertheless,

```math
\boxed{
a_S+a_T-(m_S+m_T)
=\frac{49}{48}
>1=c.
}
\tag{10.450}
```

This is already a strict weighted failure of (10.440).  To check the
quadratic-energy normalization directly, the child data are

```math
\begin{aligned}
h_S&=\frac{27}{28},&
(P(S),N(S),Q(S))&=
\left(\frac{27}{28},\frac9{28},\frac{27}{28}\right),\\
h_T&=\frac2{21},&
(P(T),N(T),Q(T))&=
\left(\frac{173}{84},\frac7{12},\frac{173}{84}\right).
\end{aligned}
\tag{10.451}
```

At the parent,

```math
P=\frac{257}{84},
\qquad
N=\frac{79}{84},
\qquad
I=|P-N|=\frac{89}{42}.
\tag{10.452}
```

The exact capacity identity (10.433) gives

```math
b_S=2,
\qquad
b_T=\frac1{28},
\qquad
\boxed{I-(b_S+b_T)=\frac1{12}.}
\tag{10.453}
```

Thus (10.435), and therefore its equivalent form (10.436), fail.
The stronger candidate (10.437) fails by the same amount:

```math
Q(S)+Q(T)-\bigl(4c-|h_S+h_T|\bigr)=\frac1{12}.
\tag{10.454}
```

The maximum-cut averaging target fails on the same certificate.  There
are three projective maximum-cut spins on `S` and six on `T`.  All rows
of the `S`--`T` block are identical, while every maximizing `S` spin is
balanced.  Hence `z=x^{\mathsf T}By=0` for all eighteen maximizing
pairs, and

```math
\boxed{
2(M_S+M_T)+|z|
=\frac{55}{56}
<\frac{49}{48}
=a_S+a_T-(m_S+m_T).
}
\tag{10.455}
```

The exact gap is `13/336`.  In the doubled favorable-energy check the
gap is `13/168`, confirming the factor of two in (10.442).

#### 10.61.2 Lifting the obstruction to complete sign matrices

The unequal magnitudes in (10.448) do not explain away the failure.
Fix `\kappa=5`; then

```math
\max_{i<j}|\kappa w_{ij}|=\frac{565}{672}<1.
\tag{10.456}
```

Replace every macro vertex `i` by a clone class `V_i` of size `L`.
Independently give every edge between `V_i` and `V_j` a sign with mean
`\kappa w_{ij}`, and give edges inside a clone class independent
unbiased signs.  Equivalently,

```math
\Pr(A_{uv}=1)=\frac{1+\kappa w_{ij}}2
\qquad(u\in V_i, v\in V_j, i\ne j).
\tag{10.457}
```

This assigns exactly one sign in `\{\pm1\}` to every unordered pair of
the `9L` vertices.

For a micro cut `R`, put `x_i=|R\cap V_i|/L`.  If independent Bernoulli
variables `Z_i` have parameters `x_i`, the expected undoubled cut value
is exactly

```math
\mathbb E C_L(R)
=\kappa L^2G(x),
\qquad
G(x)=\mathbb E\,g(\{i:Z_i=1\}).
\tag{10.458}
```

Thus `G` is the multilinear extension of the macro cut function.  The
following shell argument is what preserves exact endpoints; a bare
`o(L^2)` cut-norm estimate would not control cuts differing from an
endpoint in only a few clones.

Choose the complement of `R` if necessary, let
`d=\min\{|R|,9L-|R|\}`, and put `s=d/L\le9/2`.  The expected number of
disagreeing pairs among the nine Bernoulli variables satisfies

```math
\begin{aligned}
D(x)
&=8s-s^2+\sum_i x_i^2\\
&\ge8s-\frac89s^2
\ge4s.
\end{aligned}
\tag{10.459}
```

A nonconstant nine-bit vector has at most twenty disagreeing pairs.
With `\delta=1/112`, (10.449) therefore implies

```math
\mathbb E C_L(R)
\ge\frac{\kappa\delta}{5}Ld.
\tag{10.460}
```

Let `E=V_4\cup\cdots\cup V_8`.  Applying the same argument to
`1-g` after XOR with `E` gives, for projective Hamming distance `d_E`
from `E`,

```math
\mathbb E\bigl[C_L(E)-C_L(R)\bigr]
\ge\frac{\kappa\delta}{5}Ld_E.
\tag{10.461}
```

The difference between cuts at Hamming distance `d` uses at most
`9Ld` independent edge signs, with coefficients in `\{-1,0,1\}`.
Hoeffding's inequality, retaining half the mean margin, bounds either
failure probability by

```math
\exp\!\left(-\frac{\kappa^2\delta^2}{1800}Ld\right).
\tag{10.462}
```

There are at most two copies of each projective shell.  A union bound
for both endpoint inequalities is at most

```math
4\sum_{d=1}^{\lfloor9L/2\rfloor}
{9L\choose d}
\exp\!\left(-\frac{\kappa^2\delta^2}{1800}Ld\right)
=o(1).
\tag{10.463}
```

Indeed, `\binom{9L}{d}\le(9eL/d)^d`; the negative exponent is linear
in `Ld`, whereas the shell entropy is `O(d\log L)`.  With probability
tending to one, every other cut lies strictly between the two intended
endpoint cut values.  Hence the finite sign matrix has exactly the
required parent endpoints, not merely approximate ones.

It remains to transfer the child quantities.  For a fixed child cut,
Hoeffding gives an `\exp(-\Omega(L^2))` tail for an error of order
`L^2`.  A union bound over at most `2^{5L}` child cuts still tends to
zero.  Since a multilinear function attains its extrema at cube
vertices, simultaneously

```math
\frac{(a_X^{(L)},m_X^{(L)},M_X^{(L)})}{\kappa L^2}
\longrightarrow(a_X,m_X,M_X)
\qquad(X=S,T)
\tag{10.464}
```

in probability.  The strict positive-dominance gaps in the table above
therefore persist.  Also `c_L/(\kappa L^2)\to1`, and the exact formulas
`P_X=2a_X-4m_X` and `N_X=4M_X-2a_X` transfer all child benchmarks.
The intended parent endpoint energies retain opposite signs, since
`P_L/(\kappa L^2)\to257/84` and
`N_L/(\kappa L^2)\to79/84`.
Intersecting these events with the exact endpoint event from (10.463)
gives

```math
\boxed{
\frac{\mathcal B_L-I_L}{\kappa L^2}
\longrightarrow-\frac1{12},
\qquad
\frac{a_S^{(L)}+a_T^{(L)}-m_S^{(L)}-m_T^{(L)}-c_L}
     {\kappa L^2}
\longrightarrow\frac1{48}.
}
\tag{10.465}
```

Consequently, for every sufficiently large `L`, there exists a complete
zero-diagonal `\{\pm1\}` matrix of order `9L` that violates (10.435),
(10.437), and (10.440).  Any maximum-cut pair satisfying (10.442) would
imply (10.440), so its proposed universal existence form also fails on
these eventual blow-ups.  This is an existence proof; the crude shell
constants make no claim about the first order at which failure occurs.

#### 10.61.3 Audit scope and corrected frontier

Three independent exact scripts reconstructed the rational certificate,
the capacity values, and all eighteen maximum-cut pairs.  A separate
strictification of another weighted optimum gives the same `1/48` hard
cut gap and `1/12` stopping gap.  The probabilistic lifting argument was
also derived independently in both forms.

Small actual signings still conceal the obstruction.  Exhaustive
switching-gauge audits of (10.442) pass at order six (`1,024` classes and
`322` hard endpoint pairs) and order seven (`32,768` classes and `73,297`
hard endpoint pairs), with minimum undoubled slacks three and one.
Random audits of `2,000` matrices at each order eight through twelve also
found no failure; order nine reached equality.  These random counts are
**Numerical** only and do not conflict with the asymptotic existence
proof.

The corrected frontier is:

- (10.435)--(10.437), (10.440), and the universal selection proposal
  (10.442) are **Falsified**.  They must not be reused as proof targets;

- the exact exposure identity (10.432), the unconditional
  half-imbalance theorem (10.434), the coverage recursion (10.407), and
  all previously verified path-congestion results remain valid;

- the implication “(10.435) at every node implies `K\le4`” is still a
  correct conditional statement, but its hypothesis is false.  The
  blow-up does not supply `K_{\min}>4`, does not disprove a different
  `K\le4` argument, and does not settle the quadratic signing limit;

- the highest-priority next task is to optimize the full allocation and
  endpoint tree on the weighted obstruction and its blow-ups.  Either a
  different two-orientation or deeper stopping rule absorbs the missing
  `1/12`, or this family can be amplified into a genuine congestion lower
  bound.  Independent routes are a global tree dual, a graphon limit of
  (10.426), and new structured families whose optimized `K_{\min}` passes
  the current exact lower bound `10/3`.

### 10.62 Complemented half-stopping proves `K\le4`

The tenth wave finds the missing deeper stopping rule.  The false local
claim `\mathcal B(U)\ge I(U)` is unnecessary: the deficit is paid by
half-sized child-imbalance reservoirs.  This gives an unconditional
conserved allocation with path-cover mass at most four on every endpoint
tree.  The proof below is **Verified** and uses only (10.407),
(10.412), (10.432)--(10.434), and separate `P/N`
superadditivity.

#### 10.62.1 The exact complement inequality

At an endpoint node `U` with shores `S,T`, write

```math
I(X)=|P(X)-N(X)|,
\qquad
D(U)=I(S)+I(T),
\qquad
\mathcal B(U)=b_S+b_T.
```

For every block `X`, the definitions give

```math
Q(X)=\frac{\mathcal R(X)+I(X)}2.
```

The exact capacity formula (10.433) therefore becomes

```math
\begin{aligned}
\mathcal B(U)
&=\mathcal R(U)-Q(S)-Q(T)+|h_S|+|h_T|\\
&=\mathcal R(U)
-\frac{\mathcal R(S)+\mathcal R(T)+D(U)}2
+|h_S|+|h_T|.
\end{aligned}
\tag{10.466}
```

Separate `P/N` superadditivity gives
`\mathcal R(S)+\mathcal R(T)\le\mathcal R(U)`.  Also
`I(U)=2|h_S+h_T|`, so
`|h_S|+|h_T|\ge I(U)/2`.  Hence

```math
\boxed{
\mathcal B(U)
\ge\frac{\mathcal R(U)+I(U)-D(U)}2,
\qquad
\mathcal B(U)+\frac{D(U)}2
\ge\frac{\mathcal R(U)+I(U)}2
=Q(U)
\ge I(U).
}
\tag{10.467}
```

This is compatible with Section 10.61.  The strict weighted witness has
`\mathcal B<I`; its missing local capacity is supplied by the `D/2`
term.  Thus (10.467) complements rather than rehabilitates (10.435).

#### 10.62.2 The recursive reservoir lemma

Fix any endpoint tree below `U`.  Suppose `0\le t\le1` and `U`
receives an obligation

```math
0\le w_U\le tI(U).
\tag{10.468}
```

At `U`, make four continuous reservoirs of sizes

```math
tb_S,
\qquad
tb_T,
\qquad
\frac t2 I(S),
\qquad
\frac t2 I(T).
\tag{10.469}
```

For the first two, use one orientation attaining the best capacity
`b_X`; the last two are child-obligation reservoirs.  Their total is at
least `tI(U)` by (10.467), so amounts bounded by (10.469) can be chosen
to sum exactly to `w_U`.  In particular,

```math
w_X\le\frac t2I(X).
\tag{10.470}
```

Recurse into `X` with parameter `t/2`.  On the edge from `U` to `X`,
only one local orientation is used, so its normalized load is at most
`t`.  By induction, the child path-cover mass is at most
`2(t/2)=t`.  The exact recursion (10.412) now gives

```math
\boxed{
\theta_U
=\sum_{X=S,T}\max\{\ell_{UX},\theta_X\}
\le2t.
}
\tag{10.471}
```

This includes all boundary cases.  If `b_X=0`, its local reservoir and
allocation are zero and its load is zero.  If `I(X)=0`, its incoming
obligation is zero.  Singletons are therefore the induction base, and
every nonzero-range endpoint split strictly decreases block size.

#### 10.62.3 Root scaling and the universal bound

At the root `A`, put

```math
C_0=\sum_{X=S,T}\sum_{\sigma=\pm}(\mu_X^\sigma)_+,
\qquad
D_0=I(S)+I(T).
```

Coverage (10.407) says `C_0+D_0\ge\mathcal R(A)`.  For nonzero root
range, set

```math
t_0=\frac{\mathcal R(A)}{C_0+D_0}\le1.
\tag{10.472}
```

Allocate `t_0` times every root capacity and pass obligations
`w_X=t_0I(X)`.  The root allocations plus obligations total exactly
`\mathcal R(A)`, so this is a conserved root distribution.  Each root
edge has at most two nonzero orientation buckets and hence load at most
`2t_0`.  Applying (10.471) to each child with parameter `t_0` gives
`\theta_X\le2t_0`.  Therefore

```math
\boxed{
K=\theta_A
\le\sum_{X=S,T}\max\{2t_0,2t_0\}
=4t_0
\le4.
}
\tag{10.473}
```

If the root range is zero, the zero allocation handles the omitted
`0/0` case.  No endpoint-tie optimization is needed: (10.473) holds
for every fixed endpoint partition tree.

By (10.412), the chosen allocation is dominated coordinatewise by a
finite measure on oriented root-to-leaf chains of mass at most four.
Integrating the all-successor bound (10.341) therefore gives the
unconditional path-cover harvesting estimate

```math
\boxed{3KQ(A)\le12Q(A).}
\tag{10.474}
```

This closes the constant-congestion target from Sections 10.58--10.60.
It does not prove that four is sharp: the exact matrix (10.445) gives
the current lower bound `10/3`.  It also does not by itself settle the
original limit; the fourfold harvesting loss must still be connected
to a scale-sensitive globally minimizing deletion or insertion
recurrence.

#### 10.62.4 Exact audit and updated frontier

The construction was implemented independently with rational arithmetic,
checking every capacity, conservation identity, complement slack, and
path recursion.  It also passed all `1,024` switching classes at order
six; the largest constructed path mass there is `10/3`.  A separate
audit covered all `1,099` labelled signings through order five and all
`7,482` of their root endpoint-tie pairs, including zero capacities,
ties, and zero child imbalances.  Named fixed trees give:

| tree | `t_0` | constructed `\theta_A` |
|---|---:|---:|
| positive triangle | `2/3` | `2` |
| matrix (10.420) | `8/11` | `32/11` |
| `A_7` (10.443) | `5/7` | `20/7` |
| `A_8` (10.445) | `5/6` | `10/3` |
| reset matrix (10.403), lexicographic tree | `7/8` | `7/2` |
| family (10.416), `m=4` | `2/3` | `8/3` |

The endpoint-tree congestion interval is now

```math
\boxed{
\frac{10}{3}
\le \sup_A K_{\min}(A)
\le4.
}
\tag{10.475}
```

The lower bound is attained by the exact order-eight matrix (10.445);
the upper bound is (10.473).  The highest-priority next target is no
longer a local stopping inequality or a constant Carleson estimate.
It is to insert the verified fourfold chain harvest into the grouped
replenishment identities (10.340)--(10.342) and prove a
scale-preserving descent for globally minimizing matrices.  Improving
`4` toward `10/3` is secondary unless the constant loss becomes the
specific obstruction in that recurrence.

### 10.63 Critical-scale falsification and benign macro congestion

Two independent scope checks sharpen the interpretation of Sections
10.61--10.62.  First, regular pseudorandom blocks lift the strict
weighted obstruction to actual sign matrices with
`Q(A)=O(n^{3/2})`; the local stopping inequality remains false at the
critical scale, not only in dense biased blow-ups.  Second, the full
allocation/path minimax on the obstruction is mild: its exact macro
value is `175/99`, and even every fractional endpoint face stays below
`1.826`.  The regular-block proof and the rational macro certificates
are **Verified**.  Recursive amplification data are explicitly labelled
**Numerical**.

#### 10.63.1 A regular-block lift at `n^{3/2}` scale

Let `W=(w_{ij})` be (10.448), with strict projective endpoint gap
`\delta=1/112`.  Fix `\alpha=1/2` and `m=2` in Theorems A and B of
[Tikhomirov--Youssef, *The spectral gap of dense random regular
graphs*](https://arxiv.org/abs/1610.01765), and let `\Gamma` absorb
their two spectral constants and the fixed factors below.

Take `L\equiv1\pmod4`.  For every macro edge choose an odd integer

```math
r_{ij}=\kappa w_{ij}\sqrt L+O(1),
\qquad
d_{ij}=\frac{L-|r_{ij}|}{2}.
```

The cited directed theorem supplies an `L\times L` zero--one matrix
`H_{ij}` with every row and column sum `d_{ij}` and second singular
value `O(\sqrt L)`.  Put

```math
B_{ij}=\operatorname{sgn}(r_{ij})(J-2H_{ij}),
\qquad
E_{ij}=B_{ij}-\frac{r_{ij}}LJ.
\tag{10.476}
```

Then `B_{ij}` is a sign block with every row and column sum `r_{ij}`,
while `E_{ij}` annihilates constants and
`\|E_{ij}\|_{\rm op}\le\Gamma\sqrt L`.  Inside clone class `i`, use

```math
B_{ii}=J-I-2G_i,
```

where `G_i` is an `(L-1)/2`-regular simple graph supplied by the cited
undirected theorem.  The congruence on `L` makes this degree even;
`B_{ii}` has zero row sums and the same `O(\sqrt L)` norm.  These
theorems apply uniformly because every degree is
`L/2-O(\sqrt L)` and lies between `L^{1/2}` and `L/2` for large `L`.

Assemble the symmetric zero-diagonal sign matrix `A_L`.  If `P` is
orthogonal projection onto the nine block-constant vectors, write
`A_L=D_L+E_L`, where the off-diagonal block of `D_L` is
`(r_{ij}/L)J`.  Blockwise regularity gives

```math
E_LP=PE_L=0,
\qquad
\|E_L\|_{\rm op}\le9\Gamma\sqrt L.
\tag{10.477}
```

This exact annihilation is what removes the distance-one obstruction
that defeats independent critical biases.

Put `\beta_{ij}=r_{ij}/\sqrt L`; then
`\beta_{ij}=\kappa w_{ij}+O(L^{-1/2})`.  Let a micro cut be at
projective Hamming distance `k` from the constant spin.  The same
Bernoulli-multilinear argument as (10.458)--(10.460), now for the block
matrix `D_L`, gives

```math
C_{D_L}(R)
\ge\frac{\kappa\delta}{10}\sqrt L\,k
\tag{10.478}
```

for all sufficiently large `L`; the extra factor two leaves room for
rounding the `r_{ij}`.  If `x` is the cut spin and `u=x-Px`, then

```math
C_{A_L}(R)
=C_{D_L}(R)-\frac14u^{\mathsf T}E_Lu,
\qquad
\|u\|_2^2\le4k.
\tag{10.479}
```

Thus the centered error is at most `9\Gamma\sqrt L\,k`.  Choosing, for
example,

```math
\kappa>\frac{90\Gamma}{\delta}=10080\Gamma
\tag{10.480}
```

makes the constant spin the exact projective positive endpoint.  For
the intended `S\mid T` spin `\nu`, apply the same argument to
`-\operatorname{diag}(\nu)A_L\operatorname{diag}(\nu)`.  Its macro cut
function is

```math
g_{\beta^-}(Y)=c_\beta-g_\beta(T\mathbin\triangle Y),
```

so the strict upper gap in (10.449) makes `\nu` the exact negative
endpoint with the correct sign.

The child transfer is uniform.  If `X` contains `t\le5` macro types,
then `\|E_L[X]\|_{\rm op}\le t\Gamma\sqrt L`; hence, with
`\eta=25\Gamma/4`,

```math
\begin{aligned}
L^{3/2}(m_X(\beta)-\eta)
&\le m_X^{(L)}\le L^{3/2}m_X(\beta),\\
L^{3/2}M_X(\beta)
&\le M_X^{(L)}\le L^{3/2}(M_X(\beta)+\eta),\\
a_X^{(L)}&=L^{3/2}a_X(\beta).
\end{aligned}
\tag{10.481}
```

The macro positive-dominance gaps `9/56` and `31/84`, multiplied by
the fixed large `\kappa`, dominate `\eta`.  Both actual children are
therefore positive-dominant.  A group-constant macro minimum is an
available micro cut, so the hard gap can only improve:

```math
\begin{aligned}
&a_S^{(L)}+a_T^{(L)}-m_S^{(L)}-m_T^{(L)}-c_L\\
&\qquad\ge
L^{3/2}\bigl[
a_S(\beta)+a_T(\beta)-m_S(\beta)-m_T(\beta)-c_\beta
\bigr]
=\left(\frac{\kappa}{48}+o(1)\right)L^{3/2}>0.
\end{aligned}
\tag{10.482}
```

The parent sign margins are also strict:

```math
c_\beta-a_S(\beta)-a_T(\beta)
\longrightarrow\frac{79\kappa}{168}>0,
\qquad
c_\beta+a_S(\beta)+a_T(\beta)
\longrightarrow\frac{257\kappa}{168}>0.
\tag{10.483}
```

Thus `I_L=4(a_S^{(L)}+a_T^{(L)})`, and positive dominance plus
`a_X^{(L)}>0` gives the exact stopping deficit

```math
\boxed{
I_L-(b_S+b_T)
=4\bigl[
a_S^{(L)}+a_T^{(L)}-m_S^{(L)}-m_T^{(L)}-c_L
\bigr]>0.
}
\tag{10.484}
```

Finally, `D_L` is represented on the block-constant subspace by the
nine-by-nine matrix `(r_{ij})`.  Equations (10.476)--(10.477) give

```math
\|A_L\|_{\rm op}=O(\sqrt L),
\qquad
Q(A_L)\le9L\|A_L\|_{\rm op}=O((9L)^{3/2}).
\tag{10.485}
```

Therefore (10.435), (10.437), and (10.440) remain false even under the
scale hypothesis `Q(A)=O(n^{3/2})`.  The implicit constant here is
large and the matrices are not known to be global minimizers.  This
does not show failure for a minimizing sequence or conflict with the
unconditional complemented construction in Section 10.62.

#### 10.63.2 Exact macro congestion and all flat endpoint faces

Scale (10.448) by `672`.  The root data are

```math
(P,N,\mathcal R,I)=(2056,632,2688,1424),
```

with unique projective endpoints.  All discrete descendant endpoint
ties yield one LP signature.  Its exact joint allocation/path optimum is

```math
\boxed{K_*^{\rm macro}=\frac{175}{99}=1+\frac{76}{99}.}
\tag{10.486}
```

Here is a compact rational primal certificate.  The root capacities and
child imbalances on the five-type and four-type shores are respectively

```math
((0,24),992),
\qquad
((48,1344),432).
```

On the five-type shore, allocate `608/33` in capacity `24` and pass
`29488/33`; then allocate `22496/33` in capacity `888`, pass
`6992/33`, and allocate the latter in capacity `276`.  Every active load
on this spine is `76/99`.  On the four-type shore, allocate `1344`, pass
`432`, and allocate that in capacity `432`, giving load one.  All
conservation and imbalance caps are equalities or immediate inequalities,
so (10.412) gives `1+76/99`.

A matching dual uses edge-antichain weights

```math
\begin{array}{c|ccccc}
\text{location}
&\text{root--five}&\text{root--four}&\text{five children}
&\text{active three children}&\text{four children}\\ \hline
q&2/99&7/11&74/99&23/99&4/11.
\end{array}
```

Every root-to-leaf sum is at most one.  Put node potential `1/1188`
on the active spine, bucket slack `7/19008` only on the root capacity
`1344`, and zero-cost obligation slacks `1/1188` at the four intervening
zero-imbalance nodes.  All dual inequalities hold, and its objective is

```math
2688\frac1{1188}-1344\frac7{19008}
=\frac{175}{99},
```

proving (10.486) without solver tolerances.

Dense clone blow-ups have additional continuous endpoint ties inside a
flat macro face.  Exact cube enumeration shows that the six negative
endpoints induce a six-cycle: six symmetry-equivalent one-dimensional
projective faces and no higher-dimensional face.  At `u=0`, the surviving
`F_1` triangle with weights `(83,-55,-55)` has one further projective
endpoint edge, but both child high capacities remain `276` along it, so
`H_1(z)=\min\{220,276z\}` and `K` are unchanged.  Thus it remains only to
let `u\in[0,1/2]` be the smaller fraction of the split clone class.  Exact
elimination of the allocation LP gives child service

```math
H_s(z)=\min\{220s,(110+166s)z\},
```

and the one-variable equation

```math
888\theta+
\max_{z_1+z_2=\theta}
\bigl(H_u(z_1)+H_{1-u}(z_2)\bigr)
=912-24\theta,
\qquad
K(u)=1+\theta.
\tag{10.487}
```

Ordering the two piecewise-linear slopes gives the explicit solution.  If

```math
u_0=\frac{18311-\sqrt{313271161}}{9130},
```

then

```math
\theta(u)=
\begin{cases}
\displaystyle\frac{912}{1188-166u},&0\le u\le u_0,\\[6pt]
\displaystyle\frac{215192-41832u-73040u^2}
{(276-166u)(1022+166u)},&u_0\le u\le1/2.
\end{cases}
```

The first branch is increasing.  After clearing the positive squared
denominator, the derivative numerator of the second is
`14848880608-29345416256u+7892258848u^2`, which is positive on
`[0,1/2]`.  Thus `K(u)` is increasing.
Exact primal and dual certificates at both endpoints give

```math
\boxed{
\frac{175}{99}
\le K(u)
\le\frac{2017}{1105}
=1.825339\ldots .
}
\tag{10.488}
```

The independently strictified witness from Section 10.61 has corrected
values `1861/1053` at a vertex endpoint and `3575/1959` at its half-split
face.  The correction matters: its internal `C`--`D` edge is
`245/6000=49/1200`, not `251/6000`; the `1/1000` strictifying addition
applies only across the parent shores.  For the vertex value, an exact
dual certificate has edge-antichain weights

```math
\frac8{351},\ \frac{661}{1053}
```

on the two root edges, `784/1053` on both edges below the flat child,
`245/1053` on both edges below its active three-type child, and
`392/1053` on both edges below the other child.  Put node potential
`1/10530` on the five active nodes, the same zero-cost obligation slack
on the four intervening zero-imbalance nodes, and sole positive-cost
bucket slack

```math
z_{\mathrm{root},\,\mathrm{other},\,\mathrm{high}}
=\frac{539}{12636000}.
```

The two longest antichains sum to one, and the objective is

```math
24000\frac1{10530}
-12000\frac{539}{12636000}
=\frac{1861}{1053}.
```

Together with the exact primal certificate this proves the stated vertex
value; the denominator `12636000` explains why a first generic rational
reconstruction with cutoff `10^7` missed the dual.

For completeness, double the half-split tree so its root capacities are
`(0,480)` and `(480,24000)`.  A primal certificate allocates
`258560/653` in the first high bucket and passes `10293920/653`; on the
other root edge it allocates `24000` and passes `7840`.  The two flat
children then receive high-bucket allocations `8960` and
`7786240/1959`, pass respectively `1960` and `1703240/1959`, and those
passed obligations are allocated in capacity `3430`; the other branch
allocates `7840` in capacity `7840`.  Its congestion is `3575/1959`.
A matching dual has root weights `16/653,1175/1959`, flat-child weights
`1568/1959` on both edges, their next weights `343/1959`, and other-child
weights `784/1959`.  Put potential `1/19590` on the seven active nodes,
the same zero-cost obligation slack at the four intervening nodes, and
sole positive-cost bucket slack `49/1880640` on the root capacity
`24000`.  The long path sums equal one and the objective is

```math
48000\frac1{19590}-24000\frac{49}{1880640}
=\frac{4800-1225}{1959}
=\frac{3575}{1959}.
```

**Numerical amplification diagnostic.**  Formal recursive leaf grafts
of the primary macro signature have depths one, two, and three values
`175/99,401/228,401/228`.  A more detailed block-uniform hierarchy gives
approximately `1.7677,1.5925,1.5906`; endpoint dominance in this second
model is assumed rather than proved.  Both diagnostics decrease instead
of amplify congestion.  The rigorous conclusions are (10.486)--(10.488),
not a recursive realizability theorem.

The nine-type obstruction therefore teaches two separate lessons.  Local
stopping can fail even at the correct quadratic scale, so it cannot be
repaired by a bare `Q=O(n^{3/2})` hypothesis.  Yet the full tree allocation
absorbs that failure with congestion below two.  The active frontier
remains the global use of the universal `K\le4` harvest, not amplification
of this particular local certificate.

### 10.64 Weighted harvesting, replacement stability, and boundary-state rank

The eleventh wave tests three ways to turn the universal path cover into a
scale-sensitive recurrence.  It produces one genuinely weighted harvesting
theorem, two exact global-minimality replacement identities, and an exact
coding-theoretic boundary-rank calculation.  It also identifies a sharp wall
for each route.  All algebraic statements and finite certificates below are
**Verified**.  None of them proves convergence of `M_n/n^{3/2}`.

#### 10.64.1 Arbitrarily weighted endpoint harvesting

Fix an endpoint tree and the allocation from Section 10.62.  For an oriented
bucket `b=(v,X,\sigma)`, write

```math
c_b=(\mu_{vX}^{\sigma})_+,
\qquad 0\le a_b\le c_b.
```

The path-cover theorem supplies a measure `\nu` of mass at most four on
oriented root-to-leaf chains, with
`\nu\{\pi:b\in\pi\}\ge a_b/c_b` for every positive-capacity bucket.  Thus,
for every choice of nonnegative weights `\omega_b`, coordinatewise
domination and Tonelli give

```math
\boxed{
\sum_b\omega_ba_b
\le
\int\sum_{b\in\pi}\omega_bc_b\,d\nu(\pi)
\le
4\sup_\pi\sum_{b\in\pi}\omega_bc_b.
}
\tag{10.489}
```

This is stronger than the unweighted `12Q` consequence.  If the weights
depend only on the tree edge, write a chain split as
`U_t=D_t\sqcup X_t`, where `X_t=U_{t+1}`, and put

```math
d_t=Q(U_t)-Q(X_t)\ge0.
```

The selected capacity is bounded by the all-successor layer, and the local
step behind (10.341) gives

```math
c_b\le\mathcal L_t\le d_t+Q(D_t).
```

Consequently

```math
\boxed{
\sum_b\omega_{e(b)}a_b
\le4\sup_\pi
\sum_{t\in\pi}\omega_t\bigl(d_t+Q(D_t)\bigr).
}
\tag{10.490}
```

There is a useful exact scale-local form.  Retain only edges with parent
order at least `r` and peeled sibling order at most `s`, and give such an
edge weight `|U_t|^{-3/2}`.  Along a chain the decrements telescope, the
siblings are disjoint, and `Q(D)\le |D|(|D|-1)`.  Hence

```math
\boxed{
\sum_{b:\ |U_b|\ge r,\ |D_b|\le s}
\frac{a_b}{|U_b|^{3/2}}
\le
\frac4{r^{3/2}}
\bigl(Q(A)+(s-1)n\bigr).
}
\tag{10.491}
```

For `r=\rho n`, the diagonal-block contribution is `o(1)` when
`s=o(\sqrt n)`.  An intrinsic variant is logarithmic.  Fix
`0<q_*\le Q(A)` and retain only the initial chain edges for which both the
parent and child norms are at least `q_*`.  Since
`1-x\le-\log x`, (10.490) with weight `1/Q(U_t)` gives

```math
\boxed{
\sum_{b:\,Q(U_b),Q(X_b)\ge q_*}\frac{a_b}{Q(U_b)}
\le4\sup_\pi\left[
\log\frac{Q(A)}{q_*}
+\sum_{t:\,Q(U_t),Q(X_t)\ge q_*}
\frac{Q(D_t)}{Q(U_t)}
\right].
}
\tag{10.492}
```

Equations (10.489)--(10.492) control the selected endpoint allocations.
They do not identify the temporal replenishment gaps from
(10.340) with those allocations.  That same-window coupling, not the
absence of a weighted path theorem, is now the missing step.

#### 10.64.2 The exact grouped coefficient is a scalar tautology

Let `A` be an exact order-`n` minimizer in doubled normalization,
`Q(A)=q_n`.  Orient `A` so that an absolute ground has energy `+q_n`, and
switch that ground to `1`.  Choose a uniform
`h`-vertex set `H`, put `T=H^c` and `m=n-h`, and define

```math
\begin{aligned}
e_H&=1_T^{\mathsf T}A[T]1_T,
&d_H&=q_n-Q(A[T]),\\
g_H&=Q(A[T])-e_H,
&\varepsilon_H&=Q(A[T])-q_m.
\end{aligned}
```

Here `n\ge2` and `1\le h\le n-1`.  Uniform averaging of (10.328) gives

```math
\boxed{
q_n-q_m
=\alpha_{n,h}q_n-\mathbb E(g_H-\varepsilon_H),
\qquad
\alpha_{n,h}
=\frac{2h}{n}-\frac{h(h-1)}{n(n-1)}.
}
\tag{10.493}
```

There is a decisive simplification.  Pointwise,
`g_H-\varepsilon_H=q_m-e_H`, while

```math
\mathbb Ee_H=\beta_{n,h}q_n,
\qquad
\beta_{n,h}=\frac{m(m-1)}{n(n-1)},
\qquad
\alpha_{n,h}+\beta_{n,h}=1.
```

Therefore

```math
\boxed{
\mathbb E(g_H-\varepsilon_H)
=q_m-\beta_{n,h}q_n.
}
\tag{10.494}
```

So the excess-corrected grouped identity is exact bookkeeping, not a new
replenishment estimate.  Put `\rho=m/n` and `p=h/n`.  The exact permissible
coefficient at the `3/2` scale is

```math
\Delta_{n,h}
=\rho^{3/2}-\rho^2+\frac{p\rho}{n-1}
=\rho^{3/2}(1-\sqrt\rho)+\frac{p\rho}{n-1}.
```

For any nonnegative terminal-order error `E_{n,m}`, (10.494) gives the
equivalence

```math
\boxed{
\mathbb E(g_H-\varepsilon_H)
\le\Delta_{n,h}q_n+E_{n,m}
\quad\Longleftrightarrow\quad
\frac{q_m}{m^{3/2}}
\le\frac{q_n}{n^{3/2}}+\frac{E_{n,m}}{m^{3/2}}.
}
\tag{10.495}
```

The continuum part
`\Delta(\rho)=\rho^{3/2}-\rho^2` has exact maximum

```math
\boxed{
\max_{0\le\rho\le1}\Delta(\rho)
=\frac{27}{256}
\quad\text{at}\quad \rho=\frac9{16}.
}
\tag{10.496}
```

When `n` is even, at `h=n/2` it is `(\sqrt2-1)/4`, plus the finite
correction `1/[4(n-1)]`.  For `h=o(n)`, the continuum part is
`h/(2n)-5h^2/(8n^2)+O((h/n)^3)`; the exact `\Delta_{n,h}` additionally
contains the finite-population term
`h(n-h)/[n^2(n-1)]`.

The terminal excess cannot simply be discarded.  For the order-nine
minimizer (10.298), every order-eight child has norm `24`, while
`q_8=20` and `q_9=24`.  Thus

```math
\boxed{
\mathbb Eg_H=\frac{16}{3},
\qquad
\mathbb E\varepsilon_H=4,
\qquad
\mathbb E(g_H-\varepsilon_H)=\frac43,
}
\tag{10.497}
```

whereas the raw singleton budget is only
`\Delta_{9,1}q_9=\frac89(16\sqrt2-21)=1.446\ldots`.
Hence a sharp raw bound on `\mathbb Eg_H` is already false.

There is also a continuum coefficient no-go for a naive use of path
congestion.  Suppose optimistically that a macroscopic group obeyed
`G\le\kappa D` with no internal error, while
`D+G=(1-\rho^2)q_n`.  The required decrement is
`(1-\rho^{3/2})q_n`, so necessarily

```math
\boxed{
\kappa
\le
\frac{1-\rho^2}{1-\rho^{3/2}}-1
<\frac13.
}
\tag{10.498}
```

Thus even congestion one would miss the sharp scalar coefficient;
congestion four is not the issue that can be repaired by a smaller constant.
In the continuum comparison, the global `12Q` bound divided by the largest
budget (10.496) is at least `1024/9`.  At finite `n`, the exact coefficient
in `D+G` is
`\alpha_{n,h}=1-\rho^2+p\rho/(n-1)`.

For a general peeling segment ending at deterministic order
`m=\rho n`, (10.340) shows that the exact net target is

```math
\boxed{
\begin{aligned}
\mathbb E\left[\sum_{t<L}g_t-\varepsilon_L\right]
&\le
2\mathbb E\sum_{t<L}a_t
-\mathbb E\sum_{t<L}c_t\\
&\quad-(1-\rho^{3/2})q_n+E_{n,m}.
\end{aligned}
}
\tag{10.499}
```

This is equivalent to the normalized comparison in (10.495), so a useful
structural lemma must derive it by coupling the temporal gaps to
(10.489) in the same order/size window, while retaining
`\varepsilon_L` and the centered baseline.

Error uniformity also matters.  Assume comparisons are available for every
`m\in[\lceil n/2\rceil,n-1]`.  Given `M\ge N`, prescribe the canonical
chain `n_0=M` and
`n_{j+1}=\max\{N,\lceil n_j/2\rceil\}` until it reaches `N`.  A sufficient
condition for bridging an arbitrarily sparse liminf subsequence is

```math
\boxed{
\Omega(N)
:=\sup_{M\ge N}
\sum_j\frac{E_{n_j,n_{j+1}}}{n_{j+1}^{3/2}}
\longrightarrow0.
}
\tag{10.500}
```

A bare uniform `o(n^{3/2})` error need not satisfy (10.500), because it can
accumulate over arbitrarily many dyadic steps.  Either uniform bound
`E_{n,m}/m^{3/2}\le C m^{-\delta}` or
`E_{n,m}/m^{3/2}\le C(\log m)^{-1-\delta}` over all allowed pairs does
suffice.  For singleton comparisons the familiar sufficient condition is
`\sum_n E_{n,n-1}/(n-1)^{3/2}<\infty`.

#### 10.64.3 Global block replacement telescopes only in a hybrid landscape

There is an exact non-vacuous consequence of global minimality.  Split an
order-`n` minimizer as

```math
A=\begin{pmatrix}B&C\\ C^{\mathsf T}&D\end{pmatrix},
\qquad |B|=s,
```

and replace `B` by any order-`s` minimizer `B_*`.  For oriented states set

```math
\begin{aligned}
\Xi_{B_*\to D}(C)
=\max_{\sigma,x,y}\{&2\sigma x^{\mathsf T}Cy
-[q_s-\sigma x^{\mathsf T}B_*x]\\
&-[Q(D)-\sigma y^{\mathsf T}Dy]\}.
\end{aligned}
```

For the replacement signing `A_*`, expansion gives the exact identity
`Q(A_*)=q_s+Q(D)+\Xi_{B_*\to D}(C)`.  Since `Q(A_*)\ge q_n`, and dropping
the first nonnegative deficit can only increase the response,

```math
\boxed{
\mathcal L_{B\to D}
\ge[q_n-q_s-Q(D)]_+.
}
\tag{10.501}
```

Equivalently, if

```math
e_B=Q(B)-q_s,
\qquad
\alpha=Q(B)+Q(D)-Q(A),
```

then separate range superadditivity gives

```math
\boxed{
\mathcal L_{B\to D}\ge[e_B-\alpha]_+,
\qquad
\alpha\le\frac{I(B)+I(D)-I(A)}2.
}
\tag{10.502}
```

In particular, at an orientation-compatible active cut with `\alpha=0`,
all internal optimality excess in `B` must reappear in a genuine successor
layer.  This is an internal-block exchange, not the cross-only comparison
already shown vacuous in Section 10.55.1.

The comparison tensorizes exactly, but in a hybrid matrix.  Let
`U_t=H_t\sqcup U_{t+1}` be a nested tower.  Keep every original cross block
`C_t=A[H_t,U_{t+1}]`, replace each diagonal `A[H_t]` by an order-`|H_t|`
minimizer `G_t`, and replace the terminal block by a minimizer `G_L`.  Write
`\widetilde A_t` for the resulting hybrid suffix and
`\widetilde{\mathcal L}_t` for its all-successor layer capacity.  The
one-step coupled responses telescope, while
`Q(\widetilde A_0)\ge q_n`, giving

```math
\boxed{
\sum_{t<L}\widetilde{\mathcal L}_t
\ge
q_n-q_{|U_L|}-\sum_{t<L}q_{|H_t|}.
}
\tag{10.503}
```

There is also one common, not independently chosen, nested witness.  If
`x_t` is its spin on `H_t`, `y_L` its terminal spin, and
`y_{t+1}` the concatenated tail spin on `U_{t+1}`, then direct expansion of
the single variational formula for `Q(\widetilde A_0)` gives

```math
\boxed{
\max_{\sigma,(x_t),y_L}
\left[
2\sum_{t<L}\sigma x_t^{\mathsf T}C_ty_{t+1}
-\sum_{t<L}\delta_{G_t}(\sigma,x_t)
-\delta_{G_L}(\sigma,y_L)
\right]
\ge
q_n-q_{|U_L|}-\sum_{t<L}q_{|H_t|}.
}
\tag{10.504}
```

Here `\delta_G(\sigma,x)=Q(G)-\sigma x^{\mathsf T}Gx`.

The hybrid deficit landscape cannot be replaced monotonically by the
original one, even when every block involved is optimal.  An exact order-six
minimizer is

```math
A_6=
\begin{pmatrix}
0&1&1&1&1&1\\
1&0&-1&-1&1&1\\
1&-1&0&1&-1&1\\
1&-1&1&0&1&-1\\
1&1&-1&1&0&-1\\
1&1&1&-1&-1&0
\end{pmatrix},
\qquad Q(A_6)=q_6=10.
```

The endpoint pair
`(1,1,1,-1,-1,1)` and `(-1,-1,-1,1,1,1)` splits off vertex five, with
cross row `C=(1,1,1,-1,-1)`.  For the original principal five-block
`D=A_6[0,\ldots,4]`,

```math
Q(D)=q_5=8,
\qquad
\mathcal L(C,D)=2.
```

Replace it by the equally optimal signing

```math
D'=
\begin{pmatrix}
0&-1&-1&1&1\\
-1&0&1&1&-1\\
-1&1&0&-1&1\\
1&1&-1&0&-1\\
1&-1&1&-1&0
\end{pmatrix},
\qquad Q(D')=8.
```

Exact enumeration gives

```math
\boxed{
\mathcal L(C,D')=10>2=\mathcal L(C,D),
\qquad
Q\!\begin{pmatrix}D'&C^{\mathsf T}\\ C&0\end{pmatrix}=18.
}
\tag{10.505}
```

Thus arbitrary minimizing replacements have neither monotone layers nor
automatically negligible replacement cost.  This finite witness does not
exclude a compatible or asymptotically controlled choice.  A precise
intermediate stability target is

```math
\boxed{
\varepsilon_t
:=\sup_{\sigma,y}
[\delta_{A_{t+1}}(\sigma,y)
-\delta_{\widetilde A_{t+1}}(\sigma,y)]_+,
\qquad
\widetilde{\mathcal L}_t
\le\mathcal L_t+\varepsilon_t.
}
\tag{10.506}
```

A compatible replacement scheme with
`\sum_t\varepsilon_t=o(n^{3/2})` would transfer (10.503) to the original
all-successor layers.  Further same-window control would still be needed to
recover the centered recurrence (10.499).

#### 10.64.4 Exact cut-code recursion has exponential boundary rank

External search found a directly relevant recent theorem:
[Sheshadri, *Trellis State Complexity as an Exact Tropical Factorization
Rank* (arXiv:2607.23471v1)](https://arxiv.org/abs/2607.23471).
Because this is a very recent first version, its complete short argument was
reconstructed rather than imported as a black box.

For a binary code `C\subseteq\mathbb F_2^E` and coordinate split
`E=L\sqcup R`, let `C_L,C_R` be the subcodes supported on the two sides and
put

```math
W(a_L,a_R)=d((a_L,a_R),C),
\qquad
s=\dim C-\dim C_L-\dim C_R.
```

The theorem says that the min-plus factorization rank, tropical rank, and
Kapranov rank of the full conditional table `W` are all exactly

```math
\boxed{2^s.}
\tag{10.507}
```

The proof is elementary.  The quotient
`\mathcal T=P_R(C)/C_R` has `2^s` elements, and grouping decoding by its
class gives the exact door identity

```math
W(a_L,a_R)
=\min_{\tau\in\mathcal T}
\{D(a_L,\tau)+d(a_R,\tau)\}.
```

Conversely, choose one lifted codeword for every class.  The corresponding
`2^s\times2^s` submatrix has zero diagonal and positive off-diagonal
entries.  A min-plus rank-one term cannot be tight at two diagonal cells,
and the identity is the unique zero-weight tropical permutation.  This
proves all three ranks without an unverified structural assumption.

Apply the theorem to the augmented cut code `\mathcal C_n` from Section
1.4.  Split vertices as `S\sqcup T`, with both shore sizes at least three,
and split edge coordinates into internal and cross sets

```math
I=E(S)\sqcup E(T),
\qquad
X=E(S,T).
```

The code has dimension `n`.  A codeword supported on `I` is either zero or
all-one on `I`, and the same statement holds for `X`; hence both supported
subcodes have dimension one.  Therefore

```math
\boxed{
s=n-2,
\qquad
\operatorname{rank}_{\min,+}W
=\operatorname{rank}_{\rm trop}W
=2^{n-2}.
}
\tag{10.508}
```

The states are exactly pairs of projective shore-spin patterns.  Write the
received signing as `\begin{psmallmatrix}B&C\\C^{\mathsf T}&D\end{psmallmatrix}`
across `S\sqcup T`.  In this sign notation the two door costs are

```math
\frac{|I|-|H_B(x)+H_D(y)|}{2},
\qquad
\frac{|X|-|x^{\mathsf T}Cy|}{2}.
```

Thus the door identity becomes

```math
\boxed{
W(a_I,a_X)
=\frac{\binom n2}{2}
-\frac12\max_{x,y}
\left(
|H_B(x)+H_D(y)|+|x^{\mathsf T}Cy|
\right),
}
\tag{10.509}
```

which is exactly the block-gluing identity in Section 1.6.  The exponential
rank therefore proves that no universal constant-state or polynomial-state
exact separable min-plus recursion represents the full conditional table.
It is a representation statement, not a computational lower bound, and it
does not rule out approximation at the relevant asymptotic scale or a
compression restricted to globally minimizing signings.

There is an exact finite-temperature interpretation.  Put
`N=\binom n2`, `\lambda=2\beta\sqrt n`, and let
`c_\tau(a)` be the sum of the two door costs.  Since
`M(a)=N-2W(a)`, define

```math
S_a=\sum_{\tau\in\mathcal T}e^{-\lambda c_\tau(a)}.
```

One term attains the minimum and there are `2^{n-2}` terms, so

```math
\boxed{
\frac{e^{\lambda W(a)}}{2^{n-2}}
\le\frac1{S_a}\le e^{\lambda W(a)}.
}
\tag{10.510}
```

Replacing the canonical weight `e^{\lambda W(a)}` by `S_a^{-1}` therefore
changes `\log\mathfrak Z_n` by only `O(n)=o(n^2)`.  The state entropy is
harmless on the pressure scale, but the reciprocal destroys the usual
sum-product factorization: the problem asks for a negative moment of the
boundary partition function.  Reading an `m`-vertex child with the parent
fugacity sends `\beta` to `\beta\sqrt{n/m}`; the Shearer restriction in
(10.285) instead sends it asymptotically to
`\beta(m/n)^{3/2}`.  Neither preserves fixed temperature.  The classical
[Guerra--Toninelli interpolation for Gaussian
disorder](https://arxiv.org/abs/cond-mat/0204280) has a signed overlap
derivative at fixed temperature; no analogous sign is supplied for this
reciprocal boundary moment.

The exact rank calculation therefore closes the finite-state boundary route,
not the broader entropy route.  The remaining external-theory target is an
approximate state compression with `o(n^{3/2})` energy error, or a
negative-moment interpolation whose logarithmic error is `o(n^2)`.

#### 10.64.5 Updated frontier

The four lessons of this wave are complementary:

- endpoint allocations admit arbitrary scale weights and the exact local
  bounds (10.489)--(10.492), but no theorem yet maps temporal replenishment
  into those allocations in the same scale window;

- the excess-corrected uniform deletion term is literally
  `q_m-\beta_{n,h}q_n`, so its sharp scalar bound is normalized monotonicity
  in disguise.  Raw replenishment is too strong and congestion four is far
  too coarse at the required coefficient;

- global internal-block replacement gives the exact hybrid telescope
  (10.503), but arbitrary optimal replacements can change a layer from two
  to ten.  Any useful comparison must choose replacements compatibly with
  the fixed cross mosaic and control the full deficit landscape;

- exact cut-code boundary dynamic programming has `2^{n-2}` indispensable
  min-plus states.  Only approximate or optimizer-restricted compression can
  evade this barrier.

The highest-priority next target is now a **same-window compatibility
theorem**: couple the temporal obligations in (10.499) to the weighted
endpoint allocation (10.491), or choose hybrid minimizers satisfying
(10.506), with an error obeying the multiplicative tail condition
(10.500).  Independent routes are mean puncture stability with the terminal
excess retained, and approximate negative-moment boundary interpolation.

### 10.65 Decrement-tolled harvesting, replacement walls, and deletion cocycles

The twelfth wave attacks the same-window target in three independent ways.
It finds a residual endpoint harvest with no leading `Q(A)` term, but the
temporal-to-endpoint Hall coupling remains conditional and has an exact
finite wall.  It proves that the uniform hybrid stability error is a directed
`Q`-distance and that even the best boundary-aware minimizing replacement
can fail.  Finally, it gives an exact two-level deletion cocycle and a
strictly weaker adaptive tail criterion for convergence.  All identities,
finite certificates, and abstract implications below are **Verified**.
The Hall hypothesis and the new structural targets are explicitly **Open**.
Nothing in this section proves convergence.

#### 10.65.1 Exact centered temporal demand

Let a finite random peeling process start from an exact order-`n` minimizer
`A`, with `Q(A)=q_n`, and stop at deterministic order `m\ge\rho n`.  A
temporal node `v` is reached with probability `\pi_v` and has current order
`r_v`.  Conditional deletion outcome `H` gives a child `A_{v,H}` of order
`r_{v,H}`.  Put

```math
d_{v,H}=Q(A_v)-Q(A_{v,H}),
\qquad
\lambda_{v,H}
=\frac{q_n}{n^{3/2}}
\left(r_v^{3/2}-r_{v,H}^{3/2}\right).
```

If `\varepsilon_L=Q(A_L)-q_m` at a leaf, both quantities telescope on
every realized path:

```math
\boxed{
q_m-\left(\frac mn\right)^{3/2}q_n
=
\sum_v\pi_v\mathbb E_v(\lambda_{v,H}-d_{v,H})
-\mathbb E\varepsilon_L.
}
\tag{10.511}
```

For the field-proportional step in (10.329), the safe replenishment
interpretation is conditional-expectation level:

```math
\boxed{
\pi_v\mathbb E_v(\lambda-d)
=\pi_v\left[
\mathbb E_v\lambda+\mathbb E_vg-(2a_v-c_v)
\right].
}
\tag{10.512}
```

The expectation cannot be removed: `2a_v-c_v=\mathbb E_v(d+g)`, not
`d_H+g_H` outcome by outcome.

Define positive demand and negative credit by

```math
z_v=[\pi_v\mathbb E_v(\lambda-d)]_+,
\qquad
u_v=[-\pi_v\mathbb E_v(\lambda-d)]_+,
\qquad
C=\sum_v u_v+\mathbb E\varepsilon_L.
```

Then

```math
\boxed{
q_m-\left(\frac mn\right)^{3/2}q_n
=\sum_vz_v-C.
}
\tag{10.513}
```

One may instead keep outcome atoms
`z_{v,H}=\pi_v\Pr_v(H)[\lambda_{v,H}-d_{v,H}]_+`.  Jensen shows that
this finer geometric demand dominates the parent-level demand.  The single
credit `C` in (10.513) is nonlocal and atemporal; an honest local transport
theorem must either atomize its negative credits with their own compatibility
relations or deliberately state that it is using this relaxed global sink.

#### 10.65.2 Decrement-tolled residual harvesting

Fix an endpoint tree and the allocation from Section 10.62.  For an oriented
bucket `b` on

```math
e(b):\quad U_b=D_b\sqcup X_b\longrightarrow X_b,
```

write `c_b` for its capacity, `a_b` for its allocation, and
`\partial_b=Q(U_b)-Q(X_b)`.  Define

```math
\widehat a_b
=\frac{a_b}{c_b}[c_b-\partial_b]_+,
```

with value zero when `c_b=0`.  The path-cover measure from (10.489)
immediately gives, for all nonnegative weights `\omega_b`,

```math
\boxed{
\sum_b\omega_b\widehat a_b
\le
4\sup_\pi\sum_{b\in\pi}
\omega_b[c_b-\partial_b]_+
\le
4\sup_\pi\sum_{b\in\pi}\omega_bQ(D_b).
}
\tag{10.514}
```

Indeed, the chain marginal at `b` is at least `a_b/c_b`, and the local
layer bound is `c_b\le\partial_b+Q(D_b)`.  A chain uses only one
orientation bucket on a visited edge, so there is no extra factor two.

Retain edges with `|U_b|\ge\rho n` and `|D_b|\le s`.  The siblings on
a chain are disjoint, while `Q(D)\le|D|(|D|-1)`.  Therefore

```math
\boxed{
\sum_{\substack{b:\ |U_b|\ge\rho n\\ |D_b|\le s}}
\widehat a_b
\le4(s-1)n,
\qquad
\sum_{\substack{b:\ |U_b|\ge\rho n\\ |D_b|\le s}}
\frac{\widehat a_b}{|U_b|^{3/2}}
\le\frac{4(s-1)n}{(\rho n)^{3/2}}.
}
\tag{10.515}
```

This is the key positive result of the wave.  After the endpoint decrement is
paid as a toll, the `Q(A)` term in (10.491) disappears.  In particular,
`s=o(\sqrt n)` makes the unnormalized residual `o(n^{3/2})`.

What remains is an exact transport question.  Let `I` be temporal demand
atoms, `B` retained residual buckets, and
`R\subseteq I\times B` a proposed matrix-derived compatibility relation.
For dilation `\kappa\ge0`, max-flow/min-cut says that the minimum uncovered
demand in the relaxed model with the one global credit sink `C` is

```math
\boxed{
E_\kappa
=
\max_{J\subseteq I}
\left[
z(J)-C-\kappa\widehat a(N(J))
\right]_+.
}
\tag{10.516}
```

Here `N(J)=\{b:\ iRb\text{ for some }i\in J\}`.  Combining
(10.513)--(10.516) gives the conditional scalar comparison

```math
\boxed{
q_m-\left(\frac mn\right)^{3/2}q_n
\le4\kappa(s-1)n+E_\kappa.
}
\tag{10.517}
```

Thus any uniformly bounded `\kappa` would suffice after tolling if
`s=o(\sqrt n)` and the combined normalized errors obey the tail condition
(10.500).  Uniformity is required over every order pair used in the descent.
The full-set cut `J=I` is essentially the desired scalar comparison
restated; all Hall cuts are the additional content needed for a genuine
endpoint-bucket transport within this relaxed global-credit model.

An exact order-five minimizer shows why the transport cannot be pointwise or
error-free.  Take

```math
A=
\begin{pmatrix}
0&-1&1&-1&1\\
-1&0&-1&1&1\\
1&-1&0&1&1\\
-1&1&1&0&1\\
1&1&1&1&0
\end{pmatrix},
\qquad Q(A)=q_5=8.
```

For the positive ground `p=1` and negative ground
`(1,1,1,1,-1)`, the endpoint split is `1+4`.  Every allocated bucket
in this prescribed tree has zero residual: the order-four edge has
`c=\partial=0`, while both singleton buckets have
`c=\partial=8`.  Field-proportional deletion of the singleton is
deterministic and has

```math
d=0,\qquad g=8,\qquad a=4,\qquad c=0,\qquad\varepsilon_L=0,
```

so its centered demand is

```math
\boxed{
z=8\left[1-\left(\frac45\right)^{3/2}\right]
=8-\frac{64\sqrt5}{25}>0,
\qquad
\sum_b\widehat a_b=0.
}
\tag{10.518}
```

This defeats every finite dilation on that matching endpoint pair/tree.
It is not an optimized-tie obstruction: exactly five of the twenty-five
root endpoint pairs have zero total residual capacity, while the other
twenty `2+3` splits have total residual capacity four, where capacity means
`\sum_b[c_b-\partial_b]_+`.  No asymptotic globally minimizing blow-up is
known.  The surviving **Open target** is a grouped Hall-deficiency bound with
matrix-derived compatibility, localized negative credits, and a
tail-summable error.

#### 10.65.3 Uniform replacement stability is the wrong metric

For same-order symmetric zero-diagonal matrices `X,Y`, set
`\delta_X(\sigma,z)=Q(X)-\sigma z^{\mathsf T}Xz`.  The error in
(10.506) has the exact closed form

```math
\boxed{
\varepsilon(X,Y)
:=\sup_{\sigma,z}
[\delta_X(\sigma,z)-\delta_Y(\sigma,z)]_+
=Q(X)-Q(Y)+Q(X-Y).
}
\tag{10.519}
```

Triangle inequality makes the right side nonnegative.  Equivalently,

```math
Q(Y)+\varepsilon(X,Y)=Q(X)+Q(X-Y).
```

Paying uniform deficit stability therefore erases the norm reduction and
adds the full `Q`-distance.  Uniform switching cannot help: if `R` is a
uniform vertex switching of a fixed signing `G`, then
`\mathbb E RGR^{\mathsf T}=0` and convexity gives

```math
\mathbb E\varepsilon(X,RGR^{\mathsf T})
\ge2Q(X)-Q(G).
```

The exact boundary-aware quantity is weaker.  For a fixed cross mosaic `C`,
write

```math
\mathcal L_C(X)
=\max_{\sigma,z}
\left[2\|Cz\|_1-\delta_X(\sigma,z)\right]_+,
\qquad
J_C(X)=Q(X)+\mathcal L_C(X).
```

A `Q(X)`-ground removes the outer positive-part issue, giving

```math
\boxed{
J_C(X)=\max_z
\left\{2\|Cz\|_1+|z^{\mathsf T}Xz|\right\},
\qquad
\mathcal L_C(Y)-\mathcal L_C(X)
=Q(X)-Q(Y)+J_C(Y)-J_C(X).
}
\tag{10.520}
```

For a singleton extension by row `c`, global minimality gives a universal
wall.  If the original child is `X` and `G` is any order-`m` minimizer,
then `q_{m+1}=J_c(X)` and `J_c(G)\ge q_{m+1}`, hence

```math
\boxed{
\mathcal L_c(G)-\mathcal L_c(X)
\ge Q(X)-q_m.
}
\tag{10.521}
```

Cross-aware selection really is weaker than uniform distance.  In the
order-six example (10.505), the signed-permutation orbit of `D'` has 192
members.  Eleven have `Q(D-G)=16` but
`\mathcal L_C(G)=\mathcal L_C(D)=2`.  Nevertheless uniform orbit
averaging worsens the layer to mean `33/4`.

Even the best boundary selector fails in a global minimizer.  Exact
enumeration gives `q_7=18` and the minimizer

```math
A_7=
\begin{pmatrix}
0&1&1&1&1&1&1\\
1&0&1&1&-1&-1&1\\
1&1&0&1&-1&1&-1\\
1&1&1&0&1&-1&-1\\
1&-1&-1&1&0&-1&-1\\
1&-1&1&-1&-1&0&-1\\
1&1&-1&-1&-1&-1&0
\end{pmatrix}.
```

Split off vertex four, so the order-six tail is
`(1,2,3,5,6,7)` and `C=(1,1,1,1,-1,-1)`.  Then

```math
\boxed{
Q(D_6)=18,\qquad
\mathcal L_C(D_6)=0,\qquad
\min_{Q(G)=q_6=10}\mathcal L_C(G)=8.
}
\tag{10.522}
```

The best minimizing replacement shifts the whole child excess into the
boundary layer, attaining equality in (10.521).  More strongly, take tail
`(1,2,3,4,5)` and head `(6,7)`.  Its cross block is

```math
C=
\begin{pmatrix}
1&-1&1&-1&-1\\
1&1&-1&-1&-1
\end{pmatrix}.
```

Exact enumeration of all 192 labelled order-five minimizers gives

```math
\boxed{
Q(D_5)=12,\quad
\mathcal L_C(D_5)=4,\quad
J_C(D_5)=16,
\qquad
\min_{Q(G)=q_5=8}\mathcal L_C(G)=12,
\quad
\min_{Q(G)=8}J_C(G)=20.
}
\tag{10.523}
```

Thus even the centered functional `J_C` must rise by four.  This falsifies
universal greedy, lexicographic, and local orbit-selection rules.

There is a final exact global identity.  Partition a minimizer into diagonal
blocks `D_i` and fixed cross blocks `C_{ij}`.  Define the signed common-mosaic
response, with no positive part, by

```math
\Phi_{\mathbf D}
=
\max_{\sigma,(z_i)}
\left\{
2\sigma\sum_{i<j}z_i^{\mathsf T}C_{ij}z_j
-\sum_i\delta_{D_i}(\sigma,z_i)
\right\}
=q_n-\sum_iQ(D_i).
```

Replace every `D_i` by a minimizer `G_i` and call the hybrid
`\widetilde A`.  Then

```math
\boxed{
\Phi_{\mathbf G}-\Phi_{\mathbf D}
=
\sum_i\bigl(Q(D_i)-q_{|D_i|}\bigr)
+\bigl(Q(\widetilde A)-q_n\bigr).
}
\tag{10.524}
```

This is not a local-`\mathcal L_C` identity and not a pointwise comparison
at one shared maximizer.  “Common mosaic” means the same cross blocks and
state space.  The right side is at least the full internal block excess.
The surviving replacement target is therefore global: jointly choose the
blocks so that `Q(\widetilde A)-q_n=o(n^{3/2})`, or compare a centered
signed response while retaining the compulsory excess in (10.524).

#### 10.65.4 Two-level deletion is an exact cocycle

For a fixed order-`n` signing `A`, define

```math
F_A(k)=\mathbb E_{|S|=k}Q(A[S]),
\qquad
\beta_{r,k}=\frac{(k)_2}{(r)_2},
\qquad
R^A_{r\to k}=F_A(k)-\beta_{r,k}F_A(r).
```

The restriction of a positive ground shows `R^A_{r\to k}\ge0`.
Nested uniform sampling and
`\beta_{n,\ell}=\beta_{m,\ell}\beta_{n,m}` give

```math
\boxed{
R^A_{n\to\ell}
=R^A_{m\to\ell}
+\beta_{m,\ell}R^A_{n\to m}
\qquad(2\le\ell\le m\le n).
}
\tag{10.525}
```

If `\varepsilon_k^A=F_A(k)-q_k`, substitution yields

```math
\boxed{
R^A_{m\to\ell}
-\varepsilon_\ell^A
+\beta_{m,\ell}\varepsilon_m^A
=q_\ell-\beta_{m,\ell}q_m.
}
\tag{10.526}
```

Thus every linear two-level jackknife identity is a partition of the same
fresh-ground slack; it supplies no second scalar equation.

There is a simple uniform near-top estimate.  If `A` is an order-`n`
minimizer and `T` is uniform of order `m=n-h`, principal monotonicity and
arbitrary filling of the deleted incident edges give

```math
\boxed{
0\le
\mathbb E[Q(A[T])-q_m]
\le q_n-q_m
\le n(n-1)-m(m-1)
=h(2n-h-1).
}
\tag{10.527}
```

This is `o(n^{3/2})` for `h=o(\sqrt n)`.  Repetition does not make it
tail-summable: the normalized error per step is `O(h/\sqrt n)` while
there are order `dn/h` steps, so the block size cancels.

The order-nine minimizer (10.298) has an exact flat two-level obstruction.
All nine order-eight children have norm 24.  If zero-based vertex seven is
first deleted, the resulting `B` satisfies

```math
\boxed{
Q(B)=24=q_8+4,
\qquad
Q(B[-i])=22=q_7+4
\quad(i\in V(B)).
}
\tag{10.528}
```

The excess four persists with zero range and variance at both levels.
This rules out corrections whose only nonlinear input vanishes with those
dispersions.  Separately, `S\mapsto Q(A[S])` is neither submodular nor
supermodular already by order four, so generic set-function curvature has
no usable sign.

The useful positive conclusion is a weaker convergence bridge.  Put
`a_n=q_n/n^{3/2}`.  A proved directed comparison
`n\to m`, `\lceil n/2\rceil\le m<n`, has cost `\eta_{n,m}\ge0` if

```math
a_m\le a_n+\eta_{n,m}.
```

Let `d(M,N)` be the infimum of total cost over all proved directed paths
from `M` to `N`, with value infinity when none exists.  Then

```math
\boxed{
\Omega_{\rm ad}(N)
:=\sup_{M\ge N}d(M,N)\longrightarrow0
\quad\Longrightarrow\quad
\frac{q_n}{n^{3/2}}\ \text{converges}.
}
\tag{10.529}
```

To prove this, take a liminf subsequence `M_j`.  Telescoping a cheapest
path gives `a_N\le a_{M_j}+\Omega_{\rm ad}(N)`; first let `j\to\infty`
and then `N\to\infty`.  Condition (10.529) is strictly weaker than the
all-pairs canonical tail (10.500), because it optimizes over whatever sparse
comparison graph has actually been proved.

A checkable sufficient condition is a potential drop.  If `V(n)\ge0`,
`V(n)\to0`, and for every `M\ge N` a path from `M` to `N` can be
selected with

```math
\boxed{
\eta_{n,m}\le V(m)-V(n)
}
\tag{10.530}
```

on each used edge, then its total cost is at most `V(N)-V(M)\le V(N)`.
No monotonicity of `V` is needed for the telescope.

For inherited nonoptimal restrictions define
`e(B)=[Q(B)-q_{|B|}]/|B|^{3/2}`.  Start with an exact order-`n`
minimizer `B_0`, so `r_0=n` and `e(B_0)=0`.  The corresponding precise
terminal target is

```math
\boxed{
\mathbb E[e(B_{t+1})-V(r_{t+1})\mid B_t]
\le e(B_t)-V(r_t).
}
\tag{10.531}
```

Then `e(B_t)-V(r_t)` is a supermartingale.  At a bounded landing time
`L`,

```math
\mathbb Ee(B_L)\le\mathbb EV(r_L)-V(n).
```

Only for deterministic terminal order `r_L=N` may the right side be
written `V(N)-V(n)`.  This controls terminal excess only; it does not by
itself imply (10.529) or convergence.

#### 10.65.5 Updated frontier

This wave leaves a narrower set of viable mechanisms:

- decrement-tolled endpoint allocations have only discarded-block residual
  cost, but temporal demands still need a grouped, matrix-derived Hall
  coupling with localized negative credits and tail-uniform error;

- uniform hybrid stability is exactly a directed `Q`-distance, while an
  exact order-seven minimizer defeats even the best local boundary selector.
  Replacement must be chosen jointly at the global mosaic level;

- two-level deletion identities are cocycles, and terminal excess can be
  perfectly flat through two levels.  The useful replacement is an adaptive
  shortest-path tail or a one-sided terminal-excess potential drop.

The highest-priority next target is a **macroscopic grouped compatibility
theorem**: in one order window, either transport the centered temporal
demands into decrement-tolled residual buckets, or force enough compatible
negative credit/terminal excess to make the deficiency `E_\kappa` in
(10.516) small, with the resulting normalized comparison costs satisfying
the adaptive tail (10.529).  Independent routes are a jointly
near-minimizing global hybrid satisfying (10.524) and a selected-child
potential satisfying (10.531).
