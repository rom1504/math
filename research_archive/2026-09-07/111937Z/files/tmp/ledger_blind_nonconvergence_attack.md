# Ledger-blind nonconvergence attack

## Outcome

Assume

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j=\frac12x^{\mathsf T}Ax,
\qquad
M_n=\min_A\max_{x\in\{\pm1\}^n}|H_A(x)|,
\]

and assume that \(a_n=M_n/n^{3/2}\) does not converge.  The natural
scalable candidate I tested before consulting the targeted project audits was
an alternating power-of-two construction: use bent/flat finite-field spectra
at \(n_k=2^{2k}\), and a semi-bent or representation-theoretic obstruction at
\(m_k=2^{2k+1}\).  These orders are separated by a factor two, so ordinary
padding does not immediately identify their normalized values.

This candidate does **not** currently produce nonconvergence.  There is a
rigorous obstruction to every Cayley, homogeneous association-scheme, fixed
congruence, and conference-existence implementation of it.  What survives is
a precise but presently missing *global covering-radius theorem* on two
multiplicatively separated order classes.  The quantifiers in that theorem
also explain why a route-specific bad conference/Cayley family is not enough.

## 1. Sparse multiplicative scales are necessary

Two elementary inequalities are

\[
M_n\le M_{n+1}\le M_n+n. \tag{1}
\]

More generally, put an optimal order-\(n\) signing and an optimal order-\(h\)
signing on the diagonal and choose the \(n\times h\) cross block randomly.
A union bound over all Boolean pairs gives a deterministic cross block with

\[
\|B\|_{\infty\to1}
 \le \sqrt{2nh(n+h+2)\log2},
\]

and hence

\[
M_{n+h}\le M_n+M_h+\sqrt{2nh(n+h+2)\log2}. \tag{2}
\]

The same random-sign argument gives

\[
M_h\le\sqrt{2\binom h2(h+2)\log2}=O(h^{3/2}). \tag{3}
\]

Consequently

\[
h=o(n)\quad\Longrightarrow\quad
\left|\frac{M_{n+h}}{(n+h)^{3/2}}-\frac{M_n}{n^{3/2}}\right|
\longrightarrow0. \tag{4}
\]

Thus no fixed parity, congruence class, or design-existence class whose
successive gaps are \(o(n)\) can carry a distinct limiting constant.  Under the
nonconvergence assumption, low and high witnesses separated by a fixed
\(\varepsilon\) must recur at multiplicatively separated scales.  Alternating
\(4^k\) and \(2\cdot4^k\) is the simplest scale pattern not excluded by (4).

## 2. Why bent and homogeneous design mechanisms fail

Let \(V=\mathbb F_2^r\), \(n=|V|\), and define the Cayley signing

\[
(A_g)_{uv}=g(u+v)\quad(u\ne v),\qquad g:V\to\{\pm1\}.
\]

For the Boolean character \(\chi_a(u)=(-1)^{a\cdot u}\),

\[
A_g\chi_a=(W_g(a)-g(0))\chi_a.
\]

The spectral upper bound is attained by a Boolean character, so in the
present (undoubled) normalization

\[
\boxed{\operatorname{cap}(A_g)
=\frac n2\max_a|W_g(a)-g(0)|.} \tag{5}
\]

Parseval, \(\sum_aW_g(a)^2=n^2\), implies

\[
\operatorname{cap}(A_g)\ge\frac n2(\sqrt n-1)
=\left(\frac12-o(1)\right)n^{3/2}. \tag{6}
\]

Bent functions attain the flat Walsh floor, but therefore land at the
conference constant \(1/2\), not below it.  Semi-bent deterioration on the
other exponent parity is only a statement about this Cayley subclass; since
\(M_n\) minimizes over all edge signings, it gives no lower bound on \(M_n\).

The same obstruction holds for every homogeneous association-scheme fusion.
Such a signing has a constant row sum \(r\), hence the Boolean vector
\(\mathbf1\) gives

\[
\operatorname{cap}(S)\ge\frac n2|r|. \tag{7}
\]

If \(S^2=(n-1)I\), then \(r^2=n-1\), and (7) equals the spectral upper
bound: the cap is exactly \(n\sqrt{n-1}/2\).  Approximate orthogonality gives
the same \((1/2-o(1))n^{3/2}\) lower bound.  Nonabelian homogeneous fusions do
not escape this because their trivial representation is still Boolean.

Conference/Hadamard nonexistence at selected congruence classes also cannot
create a fixed gap: exact integrality changes the cap by only lower-order
amounts, and principal restriction from an order \(N=n+o(n)\) conference
matrix has

\[
\operatorname{cap}(A[n])\le\frac n2\sqrt{N-1}
=\left(\frac12+o(1)\right)n^{3/2}. \tag{8}
\]

The exact small data (for example \(M_6=5\), while \(M_7=9\)) therefore show
a finite design advantage but not a scalable constant separation.

## 3. Exact global formulation and the missing theorem

Let \(E_n=\binom n2\).  Encode an edge signing by
\(b\in\mathbb F_2^{E_n}\), and define the augmented cut code

\[
\mathcal C_n=
\{(t+z_i+z_j)_{i<j}:t,z_i\in\mathbb F_2\}.
\]

For every \(b\), adding \(t\) accounts exactly for the absolute value, and

\[
\max_x|H_b(x)|=E_n-2d(b,\mathcal C_n).
\]

Therefore

\[
\boxed{M_n=E_n-2\rho(\mathcal C_n),} \tag{9}
\]

where \(\rho\) is covering radius.  This is the required all-signings
quantity; in contrast, the cap of a Paley, Cayley, conference, or Hadamard
matrix concerns only one chosen word \(b\).

A genuine sparse-scale nonconvergence theorem would follow from the following
two statements for constants

\[
0.336493364431\ldots\le c_-<c_+\le\frac12:
\]

* **Deep holes at \(n_k=4^k\):** construct words \(b_k\) such that

  \[
  d(b_k,\mathcal C_{n_k})\ge
  \frac{E_{n_k}}2-\frac{c_-}{2}n_k^{3/2}+o(n_k^{3/2}). \tag{10}
  \]

  By (9), this gives \(M_{n_k}\le(c_-+o(1))n_k^{3/2}\).

* **A universal ceiling at \(m_k=2\cdot4^k\):** prove, for every edge word
  \(b\), that

  \[
  d(b,\mathcal C_{m_k})\le
  \frac{E_{m_k}}2-\frac{c_+}{2}m_k^{3/2}+o(m_k^{3/2}). \tag{11}
  \]

  Equivalently,
  \(\rho(\mathcal C_{m_k})\le E_{m_k}/2-(c_+/2)m_k^{3/2}+o(m_k^{3/2})\),
  and (9) gives \(M_{m_k}\ge(c_+-o(1))m_k^{3/2}\).

Equations (10)--(11) would yield two infinite subsequences separated by any
fixed \(\varepsilon<c_+-c_-\).  They are also the exact missing theorems.
Standard bent divisibility cannot prove them: nonzero augmented-cut-code
weights are \(s(n-s)\) and their complements, whose gcd is at most two, so
there is no growing Reed--Muller-style divisibility distinction between the
two exponent parities.

## Verdict

No requested algebraic construction presently proves nonconvergence.  The
initial bent/finite-field mechanism collapses exactly by (5)--(7), and dense
conference/design order classes collapse by the padding theorem.  Assuming
nonconvergence is consistent only with a multiplicative-scale geometric
phenomenon.  To turn the surviving power-of-two candidate into a theorem one
needs both an explicit nonhomogeneous, non-Cayley deep-hole construction
(10) and, more importantly, the genuinely global cut-code covering-radius
ceiling (11).  A bad or saturating example inside one construction family
cannot replace (11).
