# Post-factorization selector chaining

## 1. Proposed process

After passing to a spectrally regular core of order \(N\), let \(S\) be a
random \(m=pN\) subset and consider

\[
Z_x(S)=H_{A[S]}(x_S)-p^2H_A(x),
\qquad x\in\{\pm1\}^N.
\tag{1}
\]

The hope is to control \(\sup_x|Z_x|\) by chaining, with a constant below the
scale-transfer margin

\[
\left(p^{3/2}-p^2\right)M(A)
=p^{3/2}(1-\sqrt p)M(A).
\tag{2}
\]

This is impossible even for exact conference matrices. The obstruction is
not spectral; it comes from indexing every extension of the same restricted
spin.

## 2. Deterministic extension-redundancy lower bound

Fix \(S\), and let \(y\) attain \(P(A[S])\). Uniformly extending \(y\) to the
complement has mean full energy

\[
\mathbb E_zH_A(y,z)=H_{A[S]}(y)=P(A[S]).
\]

Therefore some extension \(x=(y,z)\) satisfies
\(H_A(x)\le P(A[S])\), and hence

\[
Z_x(S)\ge(1-p^2)P(A[S]).
\]

If instead the absolute restricted maximum is attained on the negative side,
the same argument with reversed inequalities gives a negative value of this
magnitude. Consequently,

\[
\boxed{
\sup_x|Z_x(S)|
\ge(1-p^2)M(A[S]).}
\tag{3}
\]

This holds for every signing, every support, and every operator norm.

Using the universal lower bound

\[
M(A[S])\ge(c_*+o(1))m^{3/2},
\qquad
c_*=0.3364933644\ldots,
\]

gives

\[
\sup_x|Z_x(S)|
\ge
(c_*+o(1))(1-p^2)p^{3/2}N^{3/2}.
\tag{4}
\]

## 3. Spectrally regular counterexample

Take a symmetric conference matrix \(C_N\). It already has

\[
\|C_N\|_{\mathrm{op}}=\sqrt{N-1},
\qquad
M(C_N)\le\frac12N\sqrt{N-1}.
\tag{5}
\]

Thus no Grothendieck--Pietsch deletion is needed. If the chaining bound were
strong enough for scale transfer, its right side would have to be no larger
than

\[
\left(p^{3/2}-p^2\right)M(C_N)
\le
\frac12p^{3/2}(1-\sqrt p)N^{3/2}.
\tag{6}
\]

But (4) exceeds (6) whenever

\[
c_*(1-p^2)>\frac12(1-\sqrt p),
\]

or equivalently

\[
c_*(1+p)(1+\sqrt p)>\frac12.
\tag{7}
\]

The equality threshold is

\[
p_0=0.1125900441\ldots .
\]

Therefore, for every fixed \(p>p_0\), exact conference matrices give a
spectrally optimal family on which

\[
\boxed{
\sup_x|Z_x(S)|
>
\left(p^{3/2}-p^2\right)M(C_N)}
\]

for every support \(S\), asymptotically. No generic-chaining,
majorizing-measure, or energy-entropy estimate for the full process (1) can
prove the needed bound.

At \(p=1/2\), the obstruction factor is especially clear:

\[
\frac{1-p^2}{1-\sqrt p}
=(1+p)(1+\sqrt p)=2.560660\ldots,
\]

so the universal induced lower constant alone gives \(0.8617\ldots\), well
above the conference parent constant \(1/2\).

## 4. Matrix form of the selector process

For a fixed selector \(\delta=\mathbf1_S\),

\[
2Z_x=x^\top
\left(D_\delta A D_\delta-p^2A\right)x.
\tag{8}
\]

Even if \(\|A\|_{\mathrm{op}}=O(\sqrt N)\), the process matrix in (8) has
operator norm \(O(\sqrt N)\) and Boolean norm of order \(N^{3/2}\). Chaining
can at best recover this natural scale. The lower bound (3) shows that its
constant contains extension redundancy unrelated to the actual restricted
norm.

The exact Grothendieck--Pietsch weights do not help on the counterexample:
conference matrices already have uniform leverage and the optimal spectral
scale.

## 5. Possible repairs and their status

Two nonredundant indices are:

1. restrict \(x\) to a full high-energy layer
   \(\{x:H_A(x)\ge t\}\), as required by the extension double count; or
2. quotient extensions by defining
   \[
   G_S(y)=\max_zH_A(y,z)
   \]
   and studying
   \[
   H_{A[S]}(y)-p^2G_S(y).
   \]

The first returns to the near-max cluster/cylinder problem. The second uses a
max-extension profile equivalent to the nonclosed cavity state encountered
in the rooted transfer analysis. Neither is controlled by operator norm or
ordinary chaining.

## 6. Verdict

Post-factorization spectral regularity is insufficient because the proposed
process has the wrong index space. The full-extension supremum is
deterministically too large on already spectrally optimal conference
matrices. There is no useful proportional parameter range above
\(p=0.11259\), and the black-box spectral constants do not certify the
remaining very small-\(p\) range.

Any viable chaining route must first quotient the complement-extension
degrees of freedom or impose the full-energy layer constraint; doing so
reintroduces precisely the energy-layer geometry that the scalar process was
meant to avoid.

