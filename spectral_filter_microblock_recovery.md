# Static spectral filters in microblock recovery

## Result

The natural continuous finite-\(r\) relaxation does **not** produce an
amplification inequality below the spectral constant \(1/2\).

An \(A\)-dependent odd polynomial filter can shrink the block-constant
quotient, including the soft top-eigenvector direction that defeated
independent absorption. However, a Frobenius/ANOVA conservation law
forces essentially all microscopic variance into the orthogonal
fiber modes. At natural block-mean scale, changing the quotient alters
that residual budget by only \(O(1/r)\).

Thus a static \(A^3\), or any bounded odd polynomial, cannot cancel the
soft-spin obstruction without paying the same amount in the residual
continuous relaxation. A successful recovery theorem must use a
genuinely Boolean nonlinear correlation that makes the compulsory
residual eigenvectors inaccessible to sign vectors.

## 1. Normalization and the coarse quotient

Let \(B\) be a symmetric signing on \(N=nr\) microvertices, partitioned
into \(n\) fibers of size \(r\). Write

\[
e_0=\frac{\mathbf 1}{\sqrt r}.
\]

For every off-diagonal \(r\times r\) block, decompose

\[
B_{ij}
=s_{ij}e_0e_0^\top
+q_{ij}e_0^\top
+e_0c_{ij}^\top
+R_{ij},
\tag{1.1}
\]

where

\[
\begin{aligned}
s_{ij}&=e_0^\top B_{ij}e_0,\\
q_{ij}&=(I-e_0e_0^\top)B_{ij}e_0,\\
c_{ij}&=(I-e_0e_0^\top)B_{ij}^\top e_0,\\
R_{ij}&=(I-e_0e_0^\top)B_{ij}(I-e_0e_0^\top).
\end{aligned}
\]

The natural coarse coefficient is

\[
C_{ij}=\frac{s_{ij}}{\sqrt r}.
\tag{1.2}
\]

Indeed, for block magnetizations \(m_i\), the degree-zero energy of
the block is

\[
r s_{ij}m_im_j=r^{3/2}C_{ij}m_im_j.
\]

Thus \(C=A\) recovers the seed at the correct amplification scale.
More generally, a static spectral correction replaces \(A\) by a
coarse matrix \(C\).

Because the entries of \(B_{ij}\) are signs, its Frobenius norm is
\(r\). Orthogonality of the four terms in (1.1) gives the exact
per-block conservation law

\[
\boxed{
r^2
=rC_{ij}^2+\|q_{ij}\|_2^2+\|c_{ij}\|_2^2
+\|R_{ij}\|_F^2.
}
\tag{1.3}
\]

If \(|C_{ij}|\le K\) for a fixed \(K\), then

\[
\|q_{ij}\|_2^2+\|c_{ij}\|_2^2+\|R_{ij}\|_F^2
\ge r^2-K^2r.
\tag{1.4}
\]

So every natural-scale quotient leaves a
\(1-O_K(1/r)\) fraction of the block's squared mass in microscopic
channels.

## 2. Solving the uniform one-step ANOVA relaxation

Fix \(t\in[0,1]\) and take macro magnetizations

\[
m_i=t\xi_i,\qquad \xi_i\in\{\pm1\}.
\]

Condition each fiber spin independently to have its prescribed
magnetization. With

\[
\alpha=\frac r{r-1}(1-t^2),
\]

the exact microcanonical ANOVA identity is

\[
\operatorname{Var}(\mathcal H_B\mid m)
=\alpha\sum_i\|g_i(m)\|_2^2
+\alpha^2\sum_{i<j}\|R_{ij}\|_F^2,
\tag{2.1}
\]

where

\[
g_i(m)=\sqrt r\sum_{j\ne i}m_jq_{ij}^{(i)}.
\]

Average (2.1) over independent uniform macro signs \(\xi_i\).
All mixed terms disappear, giving

\[
\begin{aligned}
\mathbb E_\xi\operatorname{Var}(\mathcal H_B\mid t\xi)
={}&\alpha r t^2
\sum_{i<j}\left(\|q_{ij}\|_2^2+\|c_{ij}\|_2^2\right)\\
&+\alpha^2\sum_{i<j}\|R_{ij}\|_F^2.
\end{aligned}
\tag{2.2}
\]

If \(t^2\ge1/r\), then

\[
\alpha r t^2\ge\alpha^2.
\tag{2.3}
\]

Using (1.3), equations (2.2)--(2.3) imply

\[
\boxed{
\mathbb E_\xi\operatorname{Var}(\mathcal H_B\mid t\xi)
\ge
\alpha^2
\sum_{i<j}\left(r^2-rC_{ij}^2\right).
}
\tag{2.4}
\]

Equality in the relaxed allocation problem is attained by regular
blocks:

\[
q_{ij}=c_{ij}=0,\qquad
\|R_{ij}\|_F^2=r^2-rC_{ij}^2.
\]

Therefore the continuous one-step minimax allocation is solved:
moving Frobenius mass from the two-fiber interaction channel into
row/column field channels cannot improve the uniform envelope. Once
we take the supremum over macro magnetizations, some sign pattern
\(\xi\) pays at least the average (2.4).

For bounded \(C_{ij}\),

\[
\sum_{i<j}\left(r^2-rC_{ij}^2\right)
=\binom n2r^2\left(1-O_K(1/r)\right).
\tag{2.5}
\]

In particular, a static filter changes the normalized residual
covariance by only \(O_K(1/r)\).

## 3. Odd polynomial filters

Let

\[
U=\frac A{\sqrt n}
\]

and let \(p\) be an odd polynomial. The natural filtered quotient is

\[
C_p=\sqrt n\,p(U).
\tag{3.1}
\]

The first nontrivial example is

\[
\boxed{
C_\eta=A-\frac{\eta}{n}A^3.
}
\tag{3.2}
\]

On an eigenvalue \(\lambda\) of \(A\), (3.2) acts by

\[
\lambda\longmapsto
\lambda\left(1-\eta\frac{\lambda^2}{n}\right).
\tag{3.3}
\]

Unlike an \(A^2\) correction, this shrinks positive and negative
spectral extremes simultaneously. It therefore does remove the
specific infinitesimal soft-spin obstruction in the quotient.

There are two immediate feasibility conditions:

1. the entries of \(C_p\) must stay \(O(1)\), or the desired
   \(r^{3/2}\) block-mean scaling is lost;
2. exact row sums require the usual parity/integrality approximation,
   which costs \(O(r^{-1/2})\) in each \(C_{ij}\) and is asymptotically
   harmless.

Subject to the first condition, (2.5) applies unchanged. Hence the
spectral shrinkage in (3.3) does not shrink the compulsory residual
profile at leading order.

For an exact conference seed,

\[
A^2=(n-1)I,
\]

so

\[
A^3=(n-1)A
\]

and

\[
\boxed{
C_\eta=
\left(1-\eta+\frac{\eta}{n}\right)A.
}
\tag{3.4}
\]

Thus the \(A^3\) correction is only a scalar attenuation of the macro
channel. Meanwhile, per off-diagonal block,

\[
\begin{aligned}
r^2-r(C_\eta)_{ij}^2
&=r^2
-r\left(1-\eta+\frac{\eta}{n}\right)^2\\
&=r^2\left(1-O(1/r)\right).
\end{aligned}
\tag{3.5}
\]

The quotient loses an order-one fraction of its amplitude, while the
microscopic budget changes by only an order-\(1/r\) fraction. This is
the precise unequal-scale obstruction.

## 4. Global continuous conservation law

The same obstruction has an exact spectral form. Suppose the blocks
are regular, so the \(n\)-dimensional block-constant subspace

\[
W=\{(v_1e_0,\ldots,v_ne_0):v\in\mathbb R^n\}
\]

is invariant. The restriction of \(B\) to \(W\) is the quotient

\[
S=\sqrt r\,C.
\]

Let \(B_\perp=B|_{W^\perp}\). Since \(B\) is a zero-diagonal signing,

\[
\|B\|_F^2=N(N-1).
\]

Orthogonal invariance gives

\[
\boxed{
N(N-1)
=r\|C\|_F^2+\|B_\perp\|_F^2.
}
\tag{4.1}
\]

Consequently,

\[
\boxed{
\|B\|_{\mathrm{op}}^2
\ge
\max\left\{
r\|C\|_{\mathrm{op}}^2,\,
\frac{N(N-1)-r\|C\|_F^2}{N-n}
\right\}.
}
\tag{4.2}
\]

If one tries to certify a normalized half-energy constant
\(c<1/2\) by continuous relaxation, both invariant pieces would need
operator norm at most \(2c\sqrt N\). But then

\[
\begin{aligned}
N(N-1)
&=\|S\|_F^2+\|B_\perp\|_F^2\\
&\le n(2c\sqrt N)^2+(N-n)(2c\sqrt N)^2\\
&=4c^2N^2.
\end{aligned}
\]

Letting \(N\to\infty\) forces

\[
\boxed{c\ge\frac12.}
\tag{4.3}
\]

This is a conservation law, not a defect of the particular cubic
filter. Every static polynomial filter merely redistributes a fixed
Hilbert--Schmidt budget between \(W\) and \(W^\perp\). It cannot make
both continuous channels small enough for a sub-\(1/2\)
amplification theorem.

## 5. Exact finite-dimensional relaxation

For completeness, the static one-step program at fixed \(n,r\) is:

\[
\begin{aligned}
\text{minimize over }&
\{s_{ij},q_{ij},c_{ij},R_{ij}\}_{i<j}\\
\text{the quantity }&
\sup_{m\in\mathcal G_r^n}
\Bigg[
\left(r\sum_{i<j}s_{ij}m_im_j\right)^2\\
&\quad+\frac r{r-1}\sum_i(1-m_i^2)
\left\|\sqrt r\sum_{j\ne i}m_jq_{ij}^{(i)}\right\|_2^2\\
&\quad+\left(\frac r{r-1}\right)^2
\sum_{i<j}(1-m_i^2)(1-m_j^2)\|R_{ij}\|_F^2
\Bigg]^{1/2}
\end{aligned}
\tag{5.1}
\]

subject to

\[
s_{ij}^2+\|q_{ij}\|_2^2+\|c_{ij}\|_2^2+\|R_{ij}\|_F^2=r^2
\tag{5.2}
\]

and to the realizability constraints imposed by sign matrices.

Prescribing a polynomial coarse filter adds

\[
s_{ij}=\sqrt r\,(C_p)_{ij}.
\tag{5.3}
\]

After averaging over \(m=t\xi\), the relaxation in the microscopic
allocation variables is minimized by \(q=c=0\), as proved in
Section 2. The unresolved part is therefore not this continuous
program. It is whether one can choose the residual sign blocks so
that their high-energy continuous modes are systematically invisible
to the Boolean cube while simultaneously anti-aligning the complete
Boolean action profile of \(C_p\).

## 6. Conclusion

The cubic correction genuinely fixes the *quotient* eigenvalue
problem, but it cannot fix the recovery problem:

\[
\boxed{
\text{macro spectral shrinkage}
\quad\Longrightarrow\quad
\text{no leading residual-variance shrinkage}.
}
\]

Static polynomial filtering is therefore exhausted as a continuous
or second-moment amplification mechanism. Any surviving route must
exploit a nonlinear Boolean property beyond spectra, Frobenius mass,
and one-step Hamming-slice ANOVA.

