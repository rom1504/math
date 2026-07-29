# Joint entropy--spectral regularization: obstruction checkpoint

## 1. Verdict

The proposed joint theorem cannot hold in the normalization required
by balanced purification.  Its entropy half contradicts a universal
bulk lower bound before any spectral information is used.

More strongly:

1. balanced duplicate-row purification is obstructed for every fixed
   or diverging block ratio;
2. the same exponential-union-bound mechanism is obstructed for
   \(h=o(n)\) insertion, with an explicit positive margin; and
3. raw asymptotic near-minimizers can have arbitrarily large
   \(\|A\|_{\rm op}/\sqrt n\), although this particular pathology is
   supported on \(o(n)\) added vertices.

Thus adaptive half-contraction and orientation-even \(A^2\) bounds
cannot rescue the proposed proof.  A replacement must exploit
correlations between bad spin events rather than sum their
probabilities.

Throughout,

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_x|H_A(x)|.
\]

## 2. The requested entropy upper bound is impossible

Define

\[
Z_A(\lambda)
=
\sum_{x\in\{\pm1\}^n}
\exp[-\lambda(M(A)-|H_A(x)|)].
\]

Every gap lies between zero and \(M(A)\).  Therefore

\[
\boxed{
Z_A(\lambda)\ge2^ne^{-\lambda M(A)}.
}
\tag{2.1}
\]

If

\[
M(A)=(c+o(1))n^{3/2},\qquad
\lambda=\frac c{\sqrt n},
\]

then

\[
\boxed{
\frac1n\log Z_A(c/\sqrt n)
\ge\log2-c^2-o(1).
}
\tag{2.2}
\]

For every \(c\le1/2\),

\[
\log2-c^2
\ge\log2-\frac14
=0.443147\ldots,
\]

whereas the proposed ceiling \(c^2/8\) is at most \(1/32\).
The contradiction has uniform rate gap

\[
\log2-\frac98c^2
\ge
\log2-\frac9{32}
=0.411897\ldots .
\tag{2.3}
\]

Deleting or modifying \(o(n)\) vertices does not change this
conclusion.  On a retained order \(m=(1-o(1))n\) core with
\(M=(c+o(1))m^{3/2}\), the identical argument gives
\((\log2-c^2-o(1))m\).

Hence no operator-norm estimate, including exact conference
regularity, can make the requested entropy statement true.

## 3. Raw spectral regularity is false for near-minimizers

Let \(A_m\) be an order-\(m\) near-minimizer and choose

\[
r=r_m\to\infty,\qquad r=o(\sqrt m).
\]

Adjoin \(r\) vertices, put all \(mr\) cross signs equal to \(+1\), and
fill the new internal block by any signing \(D_r\):

\[
B=
\begin{pmatrix}
A_m&J_{m\times r}\\
J_{r\times m}&D_r
\end{pmatrix}.
\]

The Boolean quadratic norm satisfies

\[
M(B)\le M(A_m)+mr+\binom r2
=M(A_m)+o(m^{3/2}).
\tag{3.1}
\]

By monotonicity \(M_m\le M_{m+r}\le M(B)\), so \(B\) is still an
asymptotic near-minimizer whenever \(A_m\) is.

Let \(u=\mathbf1_m/\sqrt m\) and
\(v=\mathbf1_r/\sqrt r\).  The two Rayleigh quotients of
\((u,v)/\sqrt2\) and \((u,-v)/\sqrt2\) differ by
\(2\sqrt{mr}\).  Consequently one has absolute value at least
\(\sqrt{mr}\), and

\[
\boxed{
\|B\|_{\rm op}\ge\sqrt{mr},\qquad
\frac{\|B\|_{\rm op}}{\sqrt{m+r}}\to\infty.
}
\tag{3.2}
\]

Thus raw near-optimality does not imply
\(\|A\|_{\rm op}=O(\sqrt n)\), let alone the constant
\(1.07557\ldots\).  The example is repaired by deleting the \(r=o(n)\)
new vertices, so it does not disprove a deletion theorem; it proves
that the deletion allowance is essential.

## 4. Why the two proposed structural tools do not supply the deletion

### Adaptive ground-state half-contraction

The exact closure chain gives

\[
g_{j+1}\le\frac12g_j,\qquad
\sum_jg_j<2g_1.
\]

It controls accumulated energy deficit, not deleted cardinality.  A
disagreement set may have half of the current vertices even when its
gap is one energy quantum.  Therefore the theorem permits
\(\Theta(n)\) deletion and cannot by itself produce an \(o(n)\) core.

### Orientation-even \(A^2\) defect

The statistic

\[
\|A^2-(n-1)I\|_F^2
\]

controls the mean square deviation of the squared eigenvalues.  A
positive-density spectral distortion is visible at order \(n^3\), but
\(o(n)\) exceptional eigenvalues can violate a fixed operator-norm
threshold while contributing only \(o(n^3)\).  Conversely, the exact
localization inequality

\[
Q(A)\ge\frac{|\lambda|}{\|v\|_\infty^2}
\]

forces an \(o(n)\)-supported eigenvector only for
\(|\lambda|/\sqrt n\to\infty\).  At the fixed threshold
\(1.07557\sqrt n\), it gives only an \(O(n)\) support bound.

Thus the existing orientation-even theorem isolates the right
spectral branch but does not identify an \(o(n)\) vertex set at the
constant required by the former purification calculation.

## 5. Balanced growth is universally obstructed

The complete proof is in `repair_purification_duality.md`.  Its key
additional identity is universal row-field entropy.

For

\[
\ell(b)=\mathbb E_{G\sim N(0,1)}\log\cosh(bG)
\]

and the balanced \(k\)-fold core factor

\[
S_{A,k}(b/\sqrt n)
=
\sum_xe^{-(b/\sqrt n)(M-|H_A(x)|)}
\prod_j\cosh\!\left(\frac b{\sqrt n}(Ax)_j\right)^k,
\]

one has

\[
\liminf\frac1n\log S_{A,k}(b/\sqrt n)
\ge
\log2-bc+k\ell(b).
\tag{5.1}
\]

This follows because each row field \((Ax)_j\), under a uniform spin,
is exactly a sum of \(n-1\) independent signs, regardless of \(A\).

Combining (5.1), the internal-block bulk entropy, and a ground-state
Jensen bound yields a strictly positive excess over the available
\(3/2\)-scale exponent for every integer block ratio \(k\ge1\).
For self-doubling, the obstruction is

\[
\boxed{
\log\frac{17}{8}
-(2\sqrt2-2)\log2
=0.1795498765\ldots
}
\tag{5.2}
\]

per core vertex.  If \(k\to\infty\), the excess is at least

\[
\boxed{
2(\log2)^2-\frac12\log2
=0.6143324376\ldots+o(1).
}
\tag{5.3}
\]

These bounds are independent of \(\|A\|_{\rm op}\).  Spectral
regularization therefore cannot repair this block certificate.

## 6. Unbalanced \(h=o(n)\) duplicate-row growth also fails

The independent duplicate-row purification criterion has, for a core
of order \(n\) and an inserted block of order \(h=o(n)\),

\[
\log Z_D(\lambda)
+
\log\sum_xe^{-\lambda d_A(x)}K_x(\lambda)^h
\]

in its exponent, where

\[
K_x(\lambda)
=
\frac{\cosh\lambda}{n}
\sum_j\cosh(\lambda(Ax)_j).
\]

Take \(\lambda=b/\sqrt n\) and put \(s=bc\).  Since
\(M(D)=O(h^{3/2})\),

\[
Z_D(\lambda)
\ge
2^he^{-\lambda M(D)}
=\exp[(\log2-o(1))h].
\tag{6.1}
\]

Choose a positive ground state \(x\) of the core.  Convexity of
\(\cosh\), together with

\[
\frac1n\sum_jx_j(Ax)_j=\frac{2M(A)}n,
\]

gives

\[
K_x(\lambda)
\ge
\cosh(2s-o(1)).
\tag{6.2}
\]

The available scale-preserving cross allowance is

\[
T_*=
\frac32c\,h\sqrt n+o(h\sqrt n),
\]

so \(\lambda T_*=(3/2)sh+o(h)\).  Equations (6.1)--(6.2) show that
the union-bound exponent exceeds the allowance by at least

\[
h\left[
\log2+\log\cosh(2s)-\frac32s-o(1)
\right].
\]

The bracket is positive for every \(s\ge0\).  Its unique minimum is at

\[
\tanh(2s)=\frac34,\qquad s=\frac14\log7,
\]

and equals

\[
\boxed{
\log8-\frac78\log7
=0.3767701613\ldots .
}
\tag{6.3}
\]

Thus unbalanced insertion does not evade the obstruction.

## 7. Surviving replacement target

The no-go theorems apply to exponential union bounds which count every
spin pair separately.  They do not rule out a correlation-sensitive
certificate.

Any viable replacement must remove at least one of the two universal
costs:

1. the \(2^h\) entropy of the inserted spins; or
2. the local-field factor forced by a core ground state.

This points to chaining, a cluster/quotient count of spin states, or a
dependent block construction in which many \((x,y)\) share the same
bad event.  Merely improving spectral regularity, the one-state tail
bound, or the sparse-repair free-energy floor cannot close the
amplification.
