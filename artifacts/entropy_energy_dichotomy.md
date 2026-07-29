# Entropy--energy dichotomy for proportional restriction

## Question

The energy-layer restriction lemma says that, for

\[
K=M(A),\qquad
\mathcal L_t^\pm=\{x:\ \pm H_A(x)\ge t\},
\]

a sufficient condition for an \(m\)-subset with \(M(A[S])<L\) is

\[
|\mathcal L_t^+|+|\mathcal L_t^-|
<\frac{L-t}{K-t}\,2^{n-m}.
\tag{1}
\]

At \(m\sim\alpha n\), \(L=\alpha^{3/2}K\), this suggests trying to prove that
the layer just below \(L\) has entropy rate less than \(1-\alpha\).

That program is impossible as stated: every Boolean quadratic form has a
universal high-energy Hamming cloud whose entropy rate is strictly larger
than \(1-\alpha\) at the natural threshold.

## 1. Universal noise-cloud theorem

Let \(f:\{\pm1\}^n\to\mathbb R\) be any homogeneous quadratic form and put

\[
K=\|f\|_\infty.
\]

Choose \(x^\star\) and an orientation \(\sigma\in\{\pm1\}\) such that
\(\sigma f(x^\star)=K\). For \(0<\theta<1\), define

\[
\mathcal G_\theta=\{x:|f(x)|\ge\theta K\}.
\]

Then

\[
\boxed{
\liminf_{n\to\infty}\frac1n\log_2|\mathcal G_\theta|
\ge
h_2\!\left(\frac{1-\sqrt\theta}{2}\right),}
\tag{2}
\]

where \(h_2(u)=-u\log_2u-(1-u)\log_2(1-u)\).

This uses only homogeneity of degree two and the existence of a maximizer; it
does not use the signing condition or any bound on \(K\).

### Proof

Fix an integer \(r\) and choose a uniformly random \(r\)-subset \(R\). Let
\(z_i=-1\) on \(R\), \(z_i=1\) off \(R\), and set \(X=x^\star z\).
For every \(i\ne j\), symmetry and the deterministic identity
\(\sum_i z_i=n-2r\) give

\[
\mathbb E[z_i z_j]
=\lambda_{n,r}
:=\frac{(n-2r)^2-n}{n(n-1)}.
\tag{3}
\]

Because \(f\) is homogeneous of degree two,

\[
\mathbb E[\sigma f(X)]=\lambda_{n,r}K.
\tag{4}
\]

The random variable \(\sigma f(X)\) is at most \(K\). If
\(\lambda_{n,r}>\theta\), then

\[
\mathbb P\{\sigma f(X)\ge\theta K\}
\ge\frac{\lambda_{n,r}-\theta}{1-\theta}.
\tag{5}
\]

Indeed, if the probability on the left is \(p\), its expectation is at most
\(pK+(1-p)\theta K\).

Take \(r/n\to\delta\) with

\[
0<\delta<\frac{1-\sqrt\theta}{2}.
\]

Then \(\lambda_{n,r}\to(1-2\delta)^2>\theta\), so (5) is bounded below by a
positive constant. A positive fraction of the Hamming sphere of radius \(r\)
around \(x^\star\) is therefore contained in \(\mathcal G_\theta\). Since

\[
\binom nr=2^{(h_2(\delta)+o(1))n},
\]

letting \(\delta\uparrow(1-\sqrt\theta)/2\) proves (2).

## 2. The required entropy inequality fails for every \(\alpha\)

The elementary binary-entropy bound

\[
h_2\!\left(\frac{1-u}{2}\right)\ge1-u^2,
\qquad 0\le u\le1,
\tag{6}
\]

implies from (2) that

\[
\liminf_{n\to\infty}\frac1n\log_2|\mathcal G_\theta|
\ge1-\theta.
\tag{7}
\]

At the restriction threshold \(\theta=\alpha^{3/2}\),

\[
1-\theta=1-\alpha^{3/2}>1-\alpha
\qquad(0<\alpha<1).
\tag{8}
\]

Thus for every signing \(A\), including an actual minimizer,

\[
\liminf_{n\to\infty}\frac1n
\log_2\#\{x:|H_A(x)|\ge\alpha^{3/2}M(A)\}
>1-\alpha.
\tag{9}
\]

Lowering the threshold to
\((\alpha^{3/2}-\eta)M(A)\) only increases the gap. For example, at
\(\alpha=1/2\),

\[
h_2\!\left(\frac{1-2^{-3/4}}2\right)
=0.7272918\ldots>0.5.
\]

Therefore the low-layer-entropy hypothesis in (1) is not merely unproved: it
is universally false at the scale needed for proportional restriction.
High layer entropy, by itself, cannot be the alternate branch of a useful
dichotomy, because the same branch occurs around the ground state of every
quadratic form.

## 3. Multiplicity-refined layer lemma

Raw layer cardinality overcounts configurations that cannot witness many bad
supports. This can be repaired exactly.

For \(m,L\), define

\[
D_+(L)=\max_x\mathbb P_{|S|=m}
\{H_{A[S]}(x_S)\ge L\},
\]

\[
D_-(L)=\max_x\mathbb P_{|S|=m}
\{H_{A[S]}(x_S)\le-L\}.
\tag{10}
\]

For \(0\le t<L<K=M(A)\), put

\[
\delta=\frac{L-t}{K-t}.
\]

Then the fractions \(b_\pm\) of positive- and negative-bad \(m\)-supports
satisfy

\[
\boxed{
b_+\le
\frac{D_+(L)|\mathcal L_t^+|}{\delta\,2^{n-m}},
\qquad
b_-\le
\frac{D_-(L)|\mathcal L_t^-|}{\delta\,2^{n-m}}.}
\tag{11}
\]

Hence

\[
D_+(L)|\mathcal L_t^+|
+D_-(L)|\mathcal L_t^-|
<\delta\,2^{n-m}
\tag{12}
\]

implies an \(m\)-subset \(S\) with \(M(A[S])<L\).

### Proof

For every positive-bad support \(S\), select a witness \(y_S\) with
\(H_{A[S]}(y_S)\ge L\). As in the raw layer lemma, at least
\(\delta2^{n-m}\) extensions \(x\) of \(y_S\) lie in
\(\mathcal L_t^+\). Count the resulting pairs \((S,x)\). For fixed \(x\), a
selected support can occur only if

\[
H_{A[S]}(x_S)=H_{A[S]}(y_S)\ge L.
\]

There are at most \(D_+(L)\binom nm\) such supports. This proves the first
inequality in (11); the second follows after replacing \(A\) by \(-A\).

The factor \(D_\pm\) precisely discounts the unavoidable Hamming cloud in
Section 1.

## 4. Selector concentration and its remaining spectral obstruction

Fix \(x\), switch \(A\) by \(x\), and denote the switched matrix by \(B\).
For independent Bernoulli-\(p\) vertex selectors \(\delta_i\),

\[
H_{A[S]}(x_S)=\frac12\delta^\top B\delta.
\]

Writing \(\delta=p\mathbf1+\xi\) gives

\[
\frac12\delta^\top B\delta
=p^2H_A(x)+p\,\xi^\top B\mathbf1
+\frac12\xi^\top B\xi.
\tag{13}
\]

Since

\[
\|B\mathbf1\|_2\le\sqrt n\,\|A\|_{\mathrm{op}},
\qquad
\|B\|_F^2=n(n-1),
\]

the subgaussian linear bound and Hanson--Wright inequality give, for \(u>0\),

\[
\mathbb P\left\{
H_{A[S]}(x_S)-p^2H_A(x)\ge u
\right\}
\le
C\exp\left[
-c_p\min\left(
\frac{u^2}{n\|A\|_{\mathrm{op}}^2},
\frac{u}{\|A\|_{\mathrm{op}}}
\right)\right].
\tag{14}
\]

Conditioning on \(|S|=m=pn+O(1)\) costs only a polynomial factor, and changes
the mean by \(O(K/n)\).

For \(L=\alpha^{3/2}K\), the gap above the largest conditional mean is

\[
L-\alpha^2K
=\alpha^{3/2}(1-\sqrt\alpha)K.
\tag{15}
\]

If \(K\asymp n^{3/2}\) and
\(\|A\|_{\mathrm{op}}=O(\sqrt n)\), (14) makes \(D_\pm(L)\) exponentially
small in \(n\). The refined criterion (12) can then succeed even though the
raw layer has the universal exponential size from (2).

For a general near-minimizer, however, the presently known bootstrap is only

\[
\|A\|_{\mathrm{op}}=O(n^{5/6}).
\]

At a deviation \(u\asymp n^{3/2}\), the first exponent in (14) can then be
only \(n^{1/3}\). This is subexponential and cannot compensate an
entropy-rate excess. Thus the exact remaining issue is no longer raw layer
entropy; it is a joint energy--restriction multiplicity bound stronger than
what the operator norm supplies.

## 5. Verdict

The proposed dichotomy

> low high-layer entropy, or exploit high high-layer entropy

does not distinguish any class of signings: the high-entropy branch is forced
universally by degree-two noise stability. The correct object is the weighted
quantity in (12), which measures both full-layer size and how many
proportional restrictions each high configuration can witness.

A scale-transfer proof through this route now reduces to either:

1. prove an exponential selector bound replacing (14) uniformly for
   near-minimizers, using flat coefficients and global Boolean boundedness
   rather than only \(\|A\|_{\mathrm{op}}\); or
2. show that the spectral-spike configurations responsible for the weak
   \(n^{1/3}\) exponent can be removed by deleting \(o(n)\) vertices before
   taking the proportional restriction.

Raw overlap coding or layer cardinality cannot close the problem, because
the Hamming-sphere construction (2) already saturates the relevant entropy
scale around a single ground state.

