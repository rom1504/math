# Block purification of duplicate-row insertion

## 1. Exact block maximum

Let \(A\) be an order-\(n\) signing, \(D\) an order-\(h\) signing, and
\(B\in\{\pm1\}^{n\times h}\) the old--new edge block.  Use the
half-energy normalization

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j.
\]

The block signing

\[
G=
\begin{pmatrix}
A&B\\
B^\top&D
\end{pmatrix}
\]

has

\[
H_G(x,y)=H_A(x)+x^\top By+H_D(y).
\]

Pairing \(y\) with \(-y\), and using \(H_D(-y)=H_D(y)\), gives the
exact formula

\[
\boxed{
M(G)
=
\max_{x,y}
\left(
\left|H_A(x)+H_D(y)\right|
+\left|x^\top By\right|
\right).
}
\tag{1.1}
\]

This formula is important: the new spins are not aligned clones.  They
choose an arbitrary \(y\), and the cross term remains an absolute
bilinear form after the \(y,-y\) pairing.

Put

\[
M_A=M(A),\quad M_D=M(D),\quad
d_A(x)=M_A-|H_A(x)|,\quad
d_D(y)=M_D-|H_D(y)|.
\]

The triangle inequality in (1.1) yields the exact sufficient
purification target

\[
\boxed{
|x^\top By|\le T+d_A(x)+d_D(y)\quad\hbox{for every }x,y
\ \Longrightarrow\
M(G)\le M_A+M_D+T.
}
\tag{1.2}
\]

Thus block insertion is a rectangular discrepancy problem with
additive energy deficits on both sides.

## 2. Why the one-row mixed theorem does not tensor

For \(j\in[n]\), let the \(j\)-th duplicate-row candidate be

\[
b^{(j,\varepsilon)}_i=
\begin{cases}
a_{ij},&i\ne j,\\
\varepsilon,&i=j,
\end{cases}
\qquad \varepsilon\in\{\pm1\}.
\]

For a fixed old spin \(x\), orient its local fields as in the one-row
theorem and put \(d=d_A(x)\).  The proved estimate is

\[
\mathbb E_{j,\varepsilon}|b^{(j,\varepsilon)}\cdot x|
\le
\frac{2M_A}{n}+\frac{n-2}{n}d+1.
\tag{2.1}
\]

If \(h\) duplicate rows are sampled, taking absolute values before
summing gives

\[
\begin{aligned}
\mathbb E\left[
|H_A(x)|+\sum_{\alpha=1}^h
|b^\alpha\cdot x|
\right]
&\le
M_A+\frac{2hM_A}{n}+h\\
&\quad+
\left(h-1-\frac{2h}{n}\right)d.
\end{aligned}
\tag{2.2}
\]

For every fixed \(h\ge2\) and all sufficiently large \(n\), the
coefficient of the energy deficit is positive.  The slack which made
the one-row estimate useful has reversed direction.  Correlating the
row choices cannot repair (2.2), because expectation of the sum of
absolute values depends only on the one-column marginals.

Therefore a block theorem must retain signed cancellation in
\(x^\top By\); it cannot be obtained by an empirical purification of
the one-row \(L^1\) inequality.

## 3. An exact exponential purification criterion

There is a stronger randomized statement which does retain signed
cancellation.  Independently for each new vertex \(\alpha\):

1. choose \(J_\alpha\) uniformly from \([n]\);
2. choose a diagonal fill \(\varepsilon_\alpha\) uniformly from
   \(\{\pm1\}\);
3. choose a global column sign \(s_\alpha\) uniformly from
   \(\{\pm1\}\);
4. take \(b^\alpha=s_\alpha
   b^{(J_\alpha,\varepsilon_\alpha)}\).

For fixed \(x,y\), the summands in \(x^\top By\) are independent,
symmetric, and have moment generating factor

\[
K_x(\lambda)
=
\frac{\cosh\lambda}{n}
\sum_{j=1}^n
\cosh\!\left(\lambda(Ax)_j\right).
\tag{3.1}
\]

Indeed, averaging the fill uses

\[
\frac{\cosh(\lambda(u+1))+\cosh(\lambda(u-1))}{2}
=\cosh(\lambda u)\cosh\lambda.
\]

Consequently, for every fixed signing \(D\), every \(\lambda>0\), and
every \(T\ge0\),

\[
\begin{aligned}
&\Pr\left[
\exists x,y:
|x^\top By|>T+d_A(x)+d_D(y)
\right]\\
&\quad\le
2e^{-\lambda T}
\left[
\sum_y e^{-\lambda d_D(y)}
\right]
\left[
\sum_x e^{-\lambda d_A(x)}
K_x(\lambda)^h
\right].
\end{aligned}
\tag{3.2}
\]

It follows that there is a pure duplicate-row block satisfying

\[
\boxed{
M(G)\le M_A+M_D+
\frac1\lambda
\log\!\left\{
2
\left(\sum_y e^{-\lambda d_D(y)}\right)
\left(\sum_x e^{-\lambda d_A(x)}
K_x(\lambda)^h\right)
\right\}.
}
\tag{3.3}
\]

Equations (3.2)--(3.3) account for every new-spin configuration \(y\);
there is no aligned-clone assumption.

At the relevant scale \(\lambda=a/\sqrt n\), the new statistic is

\[
\frac1n\sum_j
\cosh\!\left(\frac{a(Ax)_j}{\sqrt n}\right)
\]

jointly with the energy deficit \(d_A(x)\).  The original mixed theorem
controls only the first absolute moment
\(\sum_j|(Ax)_j|\), not this exponential local-field profile.
Spectral \(L^2\) control also does not bound it in the presence of a
small set of heavy coordinates.  Thus (3.3) is a genuine purification
theorem, but applying it at leading order requires exactly the
energy--local-field tail regularity which is currently missing.

### Balanced purification removes the heavy-coordinate defect

The independent choice of \(J_\alpha\) is unnecessarily sensitive to a
single heavy local field.  There is a stronger construction when

\[
h=kn,\qquad k\in\mathbb N.
\]

Use every old row exactly \(k\) times.  Randomize only the global sign
and diagonal fill of each column.  For every fixed \(x,y\), the exact
moment generating factor is now

\[
\boxed{
(\cosh\lambda)^{kn}
\prod_{j=1}^n
\cosh\!\left(\lambda(Ax)_j\right)^k.
}
\tag{3.4}
\]

The order of the balanced row multiset is irrelevant, because the
independent global column signs absorb the fixed \(y_\alpha\)'s.
Therefore the same union bound proves the existence of a balanced pure
block with

\[
\boxed{
\begin{aligned}
M(G)\le M_A+M_D+\frac1\lambda\log\Bigg\{&
2\left(\sum_y e^{-\lambda d_D(y)}\right)
(\cosh\lambda)^{kn}\\
&{}\times
\sum_x e^{-\lambda d_A(x)}
\prod_{j=1}^n
\cosh\!\left(\lambda(Ax)_j\right)^k
\Bigg\}.
\end{aligned}
}
\tag{3.5}
\]

Unlike (3.3), this criterion is automatically robust to isolated heavy
coordinates.  The elementary inequality

\[
\log\cosh u\le\frac{u^2}{2}
\]

gives

\[
\prod_j\cosh\!\left(\lambda(Ax)_j\right)^k
\le
\exp\!\left(\frac{k\lambda^2}{2}\|Ax\|_2^2\right).
\tag{3.6}
\]

At \(\lambda=a/\sqrt n\), a spectrally regular core
\(\|A\|_{\rm op}\le K\sqrt n\) therefore satisfies

\[
\frac{k\lambda^2}{2}\|Ax\|_2^2
\le
\frac{k a^2K^2}{2}\,n
\tag{3.7}
\]

uniformly in \(x\).  This is the first duplicate-row block estimate
which is simultaneously:

- pure rather than mixed;
- valid for arbitrary new spins;
- stable under a bounded number of heavy coordinates;
- controlled by the regular-core operator norm at the correct
  exponential scale.

With only the crude bounds on the two gap partition functions, (3.5)
reduces to

\[
\frac{T}{n^{3/2}}
\le
\frac{(1+k)\log2}{a}+\frac{k aK^2}{2}+o(1),
\tag{3.8}
\]

whose optimization is the usual leading rectangular-discrepancy loss.
Thus spectral regularity alone still does not prove amplification.  But
(3.5) identifies a sharper, now plausible missing input: the two
energy-gap partition functions must beat their crude entropies.  No
additional local-field tail theorem is needed on a regular core.

## 4. Coherent full purification has the wrong scaling

Take \(h=n\), use every duplicate row once, and choose the new internal
signing to be another copy of \(A\).  Ignoring the \(n\) mandatory
diagonal fills for one moment, this is the coherent two-fold blow-up

\[
G_0=
\begin{pmatrix}
A&A\\
A&A
\end{pmatrix}.
\]

For \(x,y\in\{\pm1\}^n\),

\[
H_{G_0}(x,y)=H_A(x+y).
\]

Let \(S=\{i:x_i=y_i\}\) and write \(z=x|_S=y|_S\).  Since \(x+y\) is
\(2z\) on \(S\) and zero off \(S\),

\[
\boxed{
H_{G_0}(x,y)=4H_{A[S]}(z).
}
\tag{4.1}
\]

Every principal submatrix satisfies

\[
M(A[S])\le M(A):
\]

extend a spin on \(S\) by independent random spins off \(S\); the
conditional mean of the full Hamiltonian is exactly the induced
Hamiltonian.  Hence

\[
\boxed{M(G_0)=4M(A).}
\tag{4.2}
\]

With legal diagonal fills \(\varepsilon_i\), the extra term is

\[
\sum_i\varepsilon_ix_iy_i
=
2\sum_{i\in S}\varepsilon_i-\sum_i\varepsilon_i,
\]

so

\[
4M(A)-n\le M(G)\le4M(A)+n.
\tag{4.3}
\]

For even \(n\), balanced fills give the lower witness \(4M(A)\) on
\(S=[n]\).  Since \(M(A)=\Theta(n^{3/2})\), this coherent purification
has asymptotic amplification factor \(4\), whereas scale preservation
requires

\[
2^{3/2}=2\sqrt2.
\]

Its normalized constant therefore worsens by the factor \(\sqrt2\).
This is a scalable obstruction to the most natural “take every mixed
candidate once” purification.

## 5. Exact finite audit with arbitrary internal signing

For the stronger finite test, fix \(B\) to contain every duplicate row
once, but optimize:

- every diagonal fill;
- every order-\(n\) internal signing \(D\).

Global signs on the columns of \(B\) need not be enumerated separately:
they are absorbed by switching the corresponding new vertices and
hence merely switch \(D\).

Exhaustive evaluation gives:

| core order \(n\) | exact \(M(A)\) | best order-\(2n\) block | factor |
|---:|---:|---:|---:|
| 4 | 4 | 10 | 2.50 |
| 5 | 4 | 13 | 3.25 |
| 6 | 5 | 18 | 3.60 |

The order-five and order-six values exceed the scale-preserving factor
\(2\sqrt2\).  These are finite obstructions, not an asymptotic theorem,
but they show that optimizing the internal signing does not
automatically repair the coherent-row loss even at the first
nontrivial orders.

For reproducibility, the order-six core is the exact minimizer encoded
by gauge mask \(220\) in `eulerian_explore.py`; the exhaustive evaluator
is `block_duplicate_audit.py`.

## 6. Verdict

Block purification is not closed by the existing mixed insertion
theorem:

1. its \(L^1\) slack reverses sign as soon as rows are accumulated;
2. arbitrary new spins turn coherent duplicate rows into a weighted
   vertex blow-up with factor \(4\);
3. even arbitrary internal signings fail to restore \(2\sqrt2\) at
   orders five and six.

The viable version is the balanced exponential criterion (3.5).
On a spectrally regular core it reduces the missing input to a uniform
estimate on the energy-gap partition functions.  More generally, the
independent-row criterion (3.3) asks for a uniform estimate on

\[
\sum_x
e^{-a d_A(x)/\sqrt n}
\left[
\frac{\cosh(a/\sqrt n)}n
\sum_j
\cosh\!\left(\frac{a(Ax)_j}{\sqrt n}\right)
\right]^{\theta n}
\]

for near-minimizing cores, coupled to the corresponding deficit
partition function of \(D\).  This is strictly stronger than ground
state entropy, magnetization profiles, or spectral regularity: it is a
joint energy--local-field large-deviation theorem.

For balanced \(h=kn\), the local-field factor can instead be absorbed
by (3.6)--(3.7).  The remaining sharp question is whether asymptotic
near-minimality forces enough free-energy deficit in

\[
\sum_x e^{-a d_A(x)/\sqrt n}
\quad\text{and}\quad
\sum_y e^{-a d_D(y)/\sqrt n}
\]

to improve the crude entropy term in (3.8) by precisely the amount
needed for scale preservation.
