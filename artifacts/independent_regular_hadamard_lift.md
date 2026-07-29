# Independent regular-Hadamard block orientations: orthogonal-channel audit

## Construction

Let \(A=(a_{ij})\) be a signing on \(n\) macro-vertices and let
\(H\in\{\pm1\}^{r\times r}\) be a regular Hadamard matrix,

\[
HH^\top=rI,\qquad H{\bf1}=\sqrt r\,{\bf1},
\qquad H^\top{\bf1}=\sqrt r\,{\bf1}.
\]

For each \(i<j\), independently choose permutation matrices
\(P_{ij},Q_{ij}\) and put

\[
B_{ij}=a_{ij}P_{ij}HQ_{ij},\qquad B_{ji}=B_{ij}^\top.
\]

Constant fiber spins reproduce the seed energy with the desired factor:

\[
(s_i{\bf1})^\top B_{ij}(s_j{\bf1})
=a_{ij}s_is_jr^{3/2}.
\]

The question is whether all fiber-orthogonal configurations contribute
only \(o((nr)^{3/2})\).  They do not.

## Theorem: a deterministic orthogonal-channel floor

For every choice of the permutations and seed signs, the cross-block
Hamiltonian has a configuration satisfying

\[
{\bf1}^\top x_i=0\quad\text{for every fiber }i
\]

and

\[
\boxed{
\left|\sum_{i<j}x_i^\top B_{ij}x_j\right|
\ge
\left(\frac19+o_n(1)\right)(nr)^{3/2}.
}
\tag{1}
\]

Thus independent block orientations cannot make the orthogonal channel
\(o_n(1)\) after normalization.  This is a deterministic obstruction,
not a failure of a union bound.

## Proof

Partition the macro-vertices into \(I,J\), with
\[
|J|=s=\lfloor n/3\rfloor,\qquad |I|=n-s.
\]
For every \(j\in J\), choose \(y_j\) independently and uniformly among
balanced sign vectors.  For \(i\in I\), define the balanced local field

\[
h_i=\sum_{j\in J}B_{ij}y_j.
\]

Regularity implies \(h_i\perp{\bf1}\).  Choose a balanced sign vector
\(x_i\) maximizing \(x_i^\top h_i\).

### Balanced support function

If \(h\perp{\bf1}\) and \(r\) is even, then

\[
\boxed{
\max_{\substack{x\in\{\pm1\}^r\\{\bf1}^\top x=0}}x^\top h
\ge\frac12\|h\|_1.
}
\tag{2}
\]

To see this, sort the coordinates and put \(+1\) on the largest half.
If at most half the coordinates are positive, the added negative
coordinates have total magnitude at most half the positive mass; if at
least half are positive, the largest half contain at least half the
positive mass.  Since the positive mass is \(\|h\|_1/2\), (2) follows.

### Exact second and fourth moments

Fix a row \(v\in\{\pm1\}^r\) of a regular Hadamard matrix, so
\(\sum_\ell v_\ell=\sqrt r\), and let \(y\) be uniform balanced.  For
\(Z=v^\top y\),

\[
\mathbb EZ=0,\qquad \mathbb EZ^2=r.
\tag{3}
\]

Using

\[
\mathbb E y_i y_j=-\frac1{r-1},\qquad
\mathbb E y_i y_jy_ky_\ell=
\frac3{(r-1)(r-3)}
\]

for distinct indices, together with
\(\sum v_i=\sqrt r\), gives the exact fourth moment

\[
\boxed{
\mathbb EZ^4=3r^2-2r-\frac{6r}{r-3}<3r^2
}
\tag{4}
\]

for \(r\ge4\).  Permutations and the sign \(a_{ij}\) do not change these
moments.

Each coordinate of \(h_i\) is a sum of \(s\) independent copies of such
variables.  Hence

\[
\mathbb Eh_i(\ell)^2=sr,\qquad
\mathbb Eh_i(\ell)^4<3s^2r^2.
\]

Log-convexity of \(L_p\) norms yields

\[
\mathbb E|h_i(\ell)|
\ge
\frac{(\mathbb Eh_i(\ell)^2)^{3/2}}
     {(\mathbb Eh_i(\ell)^4)^{1/2}}
>
\sqrt{\frac{sr}{3}}.
\tag{5}
\]

Combining (2) and (5),

\[
\mathbb E\max_{x_i\ {\rm balanced}}x_i^\top h_i
\ge\frac{r^{3/2}\sqrt s}{2\sqrt3}.
\]

Summing over \(i\in I\), some choice of the \(y_j\)'s makes the
\(I\)--\(J\) cross energy at least

\[
\frac{|I|r^{3/2}\sqrt{|J|}}{2\sqrt3}
=\left(\frac19+o_n(1)\right)(nr)^{3/2}.
\]

Finally, flip every fiber in \(I\).  This reverses the \(I\)--\(J\)
energy and preserves all energies internal to \(I\), internal to \(J\),
and inside individual fibers.  Of the two configurations, one has total
absolute energy at least the aligned cross energy, proving (1).

## Interpretation

The Gaussian benchmark improves the rigorous \(1/9\) to

\[
\sqrt{\frac2\pi}\,
\max_{0\le\theta\le1}(1-\theta)\sqrt\theta
=\frac{2\sqrt2}{3\sqrt{3\pi}}
=0.3071059\ldots,
\]

the original optimized greedy constant.  This shows exactly why random
block permutations do not wash out the orthogonal fiber modes: they form
a fresh mean-field spin system on \(nr\) microscopic spins.

The theorem proves a leading orthogonal floor.  Since the rigorous floor
is below the present universal \(0.33649\) lower constant, it does not by
itself show that the total lifted norm must exceed every possible seed
constant.  It does, however, falsify the proposed claim that the
orthogonal channel is \(o_n(1)\); any successful block amplification
would have to couple that leading channel destructively to the macrospin
channel.
