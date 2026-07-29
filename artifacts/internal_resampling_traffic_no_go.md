# Independent internal-resampling traffic no-go

## Status

The dual-cut-cone decomposition gives an exact internal-fill problem.
Unfortunately, the most direct independent-sign union bound cannot
hold at any competitive centered width.  This note proves the failure
with an exponential margin.

The conclusion is specific:

- independent internal signs are ruled out;
- a correlated law conditioned on the exact cap equations is not
  ruled out.

## 1. Exact internal-fill criterion

Switch a top state to \(\mathbf1\), orient a bottom state, and use the
decomposition
\[
A=c+h
\]
from `dual_cut_cone_diffuse_midpoint.md`.  The cross-supported part is
represented by a rectangular sign matrix \(C\), with
\[
\|C\|_{\infty\to1}=W.
\tag{1}
\]

For an internal signing \(g\) on
\[
E(U)\cup E(V),
\]
the new signing \(c+g\) has centered width at most \(W\) exactly when
\[
\boxed{
|g\cdot\delta(p,q)|
\le
\frac{W-|p^\mathsf TCq|}{2}
\qquad\text{for every }(p,q).
}
\tag{2}
\]
Indeed, this is equivalent to
\[
0\le(c\pm g)\cdot\delta(p,q)\le W.
\]
If additionally \(g\cdot\mathbf1=O(n)\), the new midpoint is \(O(n)\).

If the internal signs are independent and unbiased, then for an
orientation \((p,q)\),
\[
g\cdot\delta(p,q)\overset d=S_{k(p,q)},
\tag{3}
\]
where \(S_k\) is a sum of \(k\) independent Rademacher variables and
\[
k(p,q)
=|\delta(p)|_{E(U)}+|\delta(q)|_{E(V)}.
\tag{4}
\]
The direct union-bound sufficient condition is therefore
\[
\boxed{
\sum_{p,q}
\Pr\left\{
|S_{k(p,q)}|
>
\frac{W-|p^\mathsf TCq|}{2}
\right\}<1.
}
\tag{5}
\]

## 2. The union sum is exponentially large

Assume only the competitive upper bound
\[
W\le\left(\frac12+o(1)\right)n^{3/2}.
\tag{6}
\]
This holds for a centered-width minimizer because a conference
construction gives an absolute-energy bound of this size and
\(W(A)\le\max(P(A),Q(A))\).

Choose \(p,q\) uniformly.  Put \(u=|U|\), \(v=|V|\), so \(u+v=n\).
The internal cut size satisfies, for a \(1-o(1)\) fraction of
orientations,
\[
\begin{aligned}
k(p,q)
&=\frac{u^2+v^2}{4}-o(n^2)\\
&\ge\frac{n^2}{8}-o(n^2).
\end{aligned}
\tag{7}
\]
This follows by writing a cut of a block of size \(u\) as
\[
a(u-a)=\frac{u^2}{4}-(a-u/2)^2
\]
and applying binomial concentration in both blocks.

Also
\[
\mathbb E_{p,q}(p^\mathsf TCq)^2=uv\le\frac{n^2}{4}.
\tag{8}
\]
Chebyshev therefore shows that, for a \(1-o(1)\) fraction of
orientations,
\[
|p^\mathsf TCq|\le n^{5/4}.
\tag{9}
\]
The intersection of (7) and (9) has
\[
(1-o(1))2^n
\tag{10}
\]
orientations.

For every orientation in this intersection, the summand in (5) is at
least
\[
\Pr\{|S_k|>W/2\}.
\tag{11}
\]
Take the first attainable lattice value strictly above \(W/2\).
Stirling's formula, uniformly for
\[
k\ge n^2/8-o(n^2),\qquad W=O(n^{3/2}),
\]
gives
\[
\begin{aligned}
\Pr\{|S_k|>W/2\}
&\ge
\exp\left[
-\frac{W^2}{8k}-o(n)
\right]\\
&\ge
\exp\left[-\left(\frac14+o(1)\right)n\right].
\end{aligned}
\tag{12}
\]
The fourth-order correction in the binary relative-entropy expansion
is only \(O(1)\), since \(W^4/k^3=O(1)\).

Combining (10)--(12),
\[
\boxed{
\sum_{p,q}
\Pr\left\{
|S_{k(p,q)}|
>
\frac{W-|p^\mathsf TCq|}{2}
\right\}
\ge
\exp\left[
\left(\log2-\frac14-o(1)\right)n
\right].
}
\tag{13}
\]
In particular, the left side diverges exponentially and can never be
less than one.

## 3. Interpretation

The obstruction comes from ordinary orientations, not rare resonant
caps:

- there are \(e^{(\log2+o(1))n}\) typical orientations;
- their internal cut size is at least \(n^2/8-o(n^2)\);
- even the largest available margin \(W/2\) suppresses a Rademacher
  tail by at most \(e^{-(1/4+o(1))n}\).

Thus the exact fill criterion (2) remains valuable, but independent
internal resampling is decisively too noisy.  Conditioning on the cap
equations, negative dependence, or a deterministic discrepancy
construction is necessary.

## 4. Exact doubled-cap constraints for a correlated law

For
\[
\mathcal Z_0
=\{(p,q):|p^\mathsf TCq|=W\},
\]
equation (2) reduces to
\[
\boxed{
g\cdot\delta(p,q)=0
\qquad((p,q)\in\mathcal Z_0).
}
\tag{14}
\]
The symmetry \((p,q)\mapsto(p,-q)\) changes the sign of the cross
value but leaves the internal cut vector unchanged.  Hence the
\(+W\) and \(-W\) cap equations are the same linear constraints, not
two independent systems.

Let \(L_0\) be the span of their internal cut vectors.  Every viable
correlated ensemble must be supported on
\[
L_0^\perp\cap\{\pm1\}^{E(U)\cup E(V)}.
\tag{15}
\]
The original internal signing \(h\) proves this set is nonempty, but
it may contain only the two symmetry-related choices \(\pm h\).

The next nontrivial question is therefore:

> Does a global width minimizer force the cap-kernel sign set in (15)
> to contain a vector \(g\) with \(g\cdot\mathbf1=o(n^{3/2})\) that
> also obeys the thick-cap bands (2)?

This is a conditioned discrepancy problem.  The independent traffic
sum (5) supplies no route to it.
