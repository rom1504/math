# Compatible graphings need not be hyperfinite

## Status

There is a probabilistic group-circulant construction whose
fixed-threshold correlation graph contains a prescribed bounded-degree
Cayley expander.  Thus symmetry, flatness, positive definiteness, and
the cut/autocorrelation representation do **not** imply
hyperfiniteness.

One caveat remains important for the min--max application: the
construction below does not yet have a proved dimension-free bound on
\(\|B^2\|_{\rm op}\), nor a proved
\(Q(A)=O(n^{3/2})\) certificate.  It kills a purely structural
hyperfiniteness theorem, but a theorem with the additional
near-minimizer/spectral-tameness hypothesis is still possible.

## 1. Group-circulant flat square roots

Let \(G\) be a finite group of order \(N\).  Given
\[
a(e)=0,\qquad a(g)\in\{\pm1\}\ (g\ne e),\qquad
a(g)=a(g^{-1}),
\tag{1}
\]
define
\[
A_{x,y}=a(x^{-1}y).
\tag{2}
\]
Then \(A\) is symmetric, has zero diagonal, and is a signing off the
diagonal.  Its normalized square is the positive-definite convolution
kernel
\[
\frac{(A^2)_{x,y}}{N-1}=q(x^{-1}y),
\]
\[
\boxed{
q(h)=\frac1{N-1}\sum_{t\in G}a(t)a(h^{-1}t).
}
\tag{3}
\]
The last equality uses \(a(t^{-1}h)=a(h^{-1}t)\).  Thus a symmetric
group word gives exactly the compatible flat square root required in
the original problem.

## 2. An inverse-orbit factor

Take a family \((G_N,\mathcal S_N)\) of bounded-degree Cayley
expanders, with
\[
|\mathcal S_N|=d,\qquad
\mathcal S_N=\mathcal S_N^{-1},
\]
and \(N\to\infty\).

Let
\[
\mathcal O_N=G_N/(g\sim g^{-1})
\]
be the inverse-orbit set.  Form a multigraph \(H_N\) on
\(\mathcal O_N\): for every \(g\in G_N,s\in\mathcal S_N\), put an
edge between
\[
[g]\quad\text{and}\quad[sg].
\tag{4}
\]
Its maximum degree is at most
\[
D=2d.
\tag{5}
\]

Give every edge \(e\) of \(H_N\) an independent standard Gaussian
\(Z_e\), and every vertex \(v\) an independent standard Gaussian
\(Z_v\).  Define
\[
X_v=Z_v+\sum_{e\ni v}Z_e,\qquad
\eta_v=\operatorname{sign}X_v,
\tag{6}
\]
and lift
\[
a(g)=\eta_{[g]}.
\tag{7}
\]
Finally replace \(a(e)\) by zero.  Equation (7) is automatically
inversion-symmetric.

If \(v,w\) are adjacent in \(H_N\), their Gaussian fields share at
least one edge variable.  Hence
\[
\operatorname{Corr}(X_v,X_w)
\ge\frac1{D+1}.
\]
The arcsine identity gives the uniform positive correlation
\[
\boxed{
\mathbb E\eta_v\eta_w
\ge c_D:=
\frac2\pi\arcsin\frac1{D+1}>0.
}
\tag{8}
\]
Nonidentical, nonadjacent orbit vertices have disjoint noise sets and
are independent.

## 3. The large-edge graph contains the expander

For \(s\in\mathcal S_N\), the orbit vertices
\[
[t],\qquad[s^{-1}t]
\]
are adjacent by (4), for every \(t\).  Equations (3) and (8) imply
\[
\mathbb E q(s)\ge c_D-O(N^{-1}).
\tag{9}
\]

Only boundedly many group differences can have a large expected
correlation.  Indeed, for fixed \(t\), the Gaussian fields at
\([t]\) and \([h^{-1}t]\) overlap only if the two orbit vertices are
equal or adjacent.  There are at most \(D+1\) such orbit vertices,
and each has at most two representatives.  Therefore
\[
\boxed{
\sum_{h\in G_N}\mathbb E q(h)\le2(D+1)+o(1).
}
\tag{10}
\]
All terms are nonnegative because sign is an increasing function of
Gaussian fields having nonnegative shared-variable covariances.

The empirical correlations concentrate uniformly.  Replacing one
edge noise changes at most two orbit signs, hence at most four lifted
values \(a(g)\), and changes a fixed \(q(h)\) by \(O(N^{-1})\).
The same is true for a vertex noise.  Since there are \(O(DN)\)
independent noises, bounded differences gives
\[
\Pr\{|q(h)-\mathbb E q(h)|>\delta\}
\le2e^{-c_D'\delta^2N}.
\tag{11}
\]
A union bound covers all \(h\in G_N\).

Choose, for example,
\[
\xi=c_D/2,\qquad\delta=c_D/4.
\]
With probability tending to one:

- every generator \(s\in\mathcal S_N\) obeys
  \(q(s)>\xi\);
- if \(|q(h)|>\xi\), then
  \(\mathbb E q(h)>c_D/4\).

By (10), the number of the latter \(h\)'s is at most
\[
\frac{8(D+1)}{c_D}+o(1).
\tag{12}
\]
Thus the fixed-threshold correlation graph has uniformly bounded
degree and contains the Cayley graph
\(\operatorname{Cay}(G_N,\mathcal S_N)\).

Since the latter is an expander family, the threshold graphs are not
hyperfinite.  Any decomposition into uniformly bounded components
requires deleting a positive linear number of edges.

## 4. What this does and does not settle

The construction proves:

> A bounded-degree fixed-threshold graphing represented by
> \(C=B^2\) with a symmetric flat signing \(B\) need not be amenable
> or hyperfinite.

It complements the amenable block-palindrome/Fejer family in
`flat_square_root_sparse_classification.md`; compatible graphings can
exhibit both extremes.

The remaining escape for the min--max proof is spectral tameness.
For a group-circulant signing,
\[
\|B^2\|_{\rm op}
=\frac1{N-1}\max_\rho
\left\|\sum_{g\in G}a(g)\rho(g)\right\|_{\rm op}^2,
\tag{13}
\]
where the maximum is over irreducible representations.  The local
factor construction proves neither a uniform bound in (13) nor the
Boolean norm bound \(Q(A)=O(N^{3/2})\).

Accordingly, the sharpened structural question is:

> Are fixed-threshold graphings hyperfinite under the additional
> hypothesis that the flat square roots are spectrally tame (as
> obtained after Pietsch from a near-minimizing sequence)?

For highly quasirandom groups, matrix Fourier concentration suggests
that (13) may in fact remain bounded, which would turn the
construction above into a complete counterexample to that sharpened
claim.  This requires a separate noncommutative Fourier estimate and
is not asserted here.

## 5. Spectral audit on quasirandom groups

The standard dependency-graph/matrix-Bernstein argument does not give
a constant.  It gives exactly one logarithmic loss.

Let \(\rho\) be an irreducible unitary representation of \(G\), of
dimension \(d_\rho\).  Before changing the identity coefficient, the
Fourier block is
\[
\widehat a(\rho)
=\sum_{v=[g]\in\mathcal O}
\eta_v\bigl(\rho(g)+\rho(g^{-1})\bigr),
\tag{14}
\]
with the evident single-copy convention at involutions.  The
dependency graph of the \(\eta_v\)'s is contained in the square of
\(H_N\), hence has degree at most
\[
\Gamma\le D^2.
\]
Color it using
\[
\chi\le D^2+1
\tag{15}
\]
colors.  Within a color class the signs are independent and centered.
Each summand in (14) has operator norm at most two, and for each color
class
\[
\left\|
\sum_v\mathbb E M_vM_v^*
\right\|_{\rm op}
\le 2N.
\tag{16}
\]

Matrix Bernstein and a union bound over the color classes give
\[
\Pr\left\{
\|\widehat a(\rho)\|_{\rm op}>
\chi\left(
2\sqrt{N u}+\frac43u
\right)
\right\}
\le2\chi d_\rho e^{-u}.
\tag{17}
\]
Since a finite group has at most \(N\) irreducible representations
and \(d_\rho\le\sqrt N\), take
\[
u=4\log N+\log(2\chi).
\]
After a union bound over all Fourier blocks, with probability
\(1-O(N^{-2})\),
\[
\boxed{
\|A\|_{\rm op}
\le
(D^2+1)\left(
4\sqrt{N\log N}+O_D(\log N)
\right).
}
\tag{18}
\]
The trivial representation is included; scalar Bernstein gives the
same or a better estimate for \(\sum_ga(g)\).  Consequently
\[
\boxed{
\|B^2\|_{\rm op}=O_D(\log N).
}
\tag{19}
\]

For a concrete quasirandom family such as
\(\mathrm{PSL}_2(q)\),
\[
N\asymp q^3,\qquad
\min_{\rho\ne1}d_\rho\ge(q-1)/2,
\]
and the number of irreducibles is \(O(q)\).  Substituting these facts
in (17) changes only the constant inside \(\log N\); matrix Bernstein
still yields \(O(\sqrt{N\log N})\), not \(O(\sqrt N)\).
Quasirandomness by itself is not used by this black-box inequality.

Removing the logarithm requires one of the following genuinely
stronger inputs:

1. a trace-moment theorem showing that noncommutative word equations
   in the chosen groups have only the free/Catalan number of
   solutions;
2. a dimension-free norm theorem for random convolution operators
   with this bounded dependency; or
3. a discrepancy choice of the orbit signs that preserves the
   generator correlations while flattening every Fourier block.

Thus the nonamenable construction currently lies at spectral scale
\(O(\log N)\) for \(C=B^2\).  It is not yet a counterexample to
hyperfiniteness under a uniform spectral-tameness hypothesis, and the
logarithm must not be silently discarded.
