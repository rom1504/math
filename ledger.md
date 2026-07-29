# Quadratic Signing Limit — Research Ledger

Last updated: 2026-07-25 (America/Los_Angeles)

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
