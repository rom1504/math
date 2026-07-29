# Near-maximum cluster packing and small signed cuts

## 1. Switching a ground state: exact small-cut dictionary

Let

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad K=M(A).
\]

Choose an orientation of \(A\) and a spin \(x^\star\) so that
\(H_A(x^\star)=K\), and switch \(x^\star\) to \(\mathbf1\). Thus, after
replacing \(A\) by \(B=D_{x^\star}AD_{x^\star}\),

\[
H_B(\mathbf1)=K.
\]

For \(S\subset[n]\), let \(z_S\) be \(-1\) on \(S\), \(1\) off \(S\), and
write

\[
c_B(S)=\sum_{i\in S,\ j\notin S}b_{ij}.
\]

Then

\[
\boxed{H_B(z_S)=K-2c_B(S).}
\tag{1}
\]

Because \(K\) is the positive maximum and \(|H_B|\le K\),

\[
\boxed{0\le c_B(S)\le K\quad\text{for every }S.}
\tag{2}
\]

Thus switching by a ground state turns the signing into a cut-minimal
representative: every signed cut is nonnegative. Positive near-maxima are
exactly small cuts:

\[
H_B(z_S)\ge\theta K
\quad\Longleftrightarrow\quad
c_B(S)\le\frac{1-\theta}{2}K.
\tag{3}
\]

If two configurations correspond to \(S,T\), their Hamming distance is
\(|S\mathbin\triangle T|\). Hence a separated near-ground-state code is
exactly a separated family of small signed cuts.

## 2. A two-threshold packing theorem

For an orientation \(\sigma\in\{\pm1\}\), define

\[
\mathcal G_\theta^\sigma
=\{x:\sigma H_A(x)\ge\theta K\}.
\]

Fix \(0<\theta_0<\theta_1\le1\), and let \(r/n\to\rho<1/2\). Put

\[
\lambda_{n,r}=\frac{(n-2r)^2-n}{n(n-1)}.
\]

Suppose \(\mathcal C\subset\mathcal G_{\theta_1}^\sigma\) has pairwise
Hamming distance greater than \(2r\) and

\[
\lambda_{n,r}\theta_1>\theta_0.
\tag{4}
\]

Then

\[
\boxed{
|\mathcal C|
\le
\frac{|\mathcal G_{\theta_0}^\sigma|}
{\binom nr}
\frac{1-\theta_0}
{\lambda_{n,r}\theta_1-\theta_0}.}
\tag{5}
\]

### Proof

For \(y\in\mathcal C\), choose a uniform radius-\(r\) noise vector \(z\).
Degree-two homogeneity gives

\[
\mathbb E[\sigma H_A(yz)]
=\lambda_{n,r}\sigma H_A(y)
\ge\lambda_{n,r}\theta_1K.
\]

Since \(\sigma H_A\le K\), at least

\[
p=\frac{\lambda_{n,r}\theta_1-\theta_0}{1-\theta_0}
\]

of the radius-\(r\) sphere around \(y\) lies in
\(\mathcal G_{\theta_0}^\sigma\). The radius-\(r\) spheres around the code
centers are disjoint, proving (5).

Parseval gives the universal upper bound

\[
\frac{|\mathcal G_{\theta_0}^\sigma|}{2^n}
\le
\frac{\mathbb E H_A(X)^2}{\theta_0^2K^2}
=\frac{\binom n2}{\theta_0^2K^2}.
\tag{6}
\]

For signings, the universal \(K=\Omega(n^{3/2})\) lower bound makes the last
quantity \(O(1/n)\). Therefore, whenever
\((1-2\rho)^2\theta_1>\theta_0\),

\[
\boxed{
\limsup_{n\to\infty}\frac1n\log_2|\mathcal C|
\le1-h_2(\rho).}
\tag{7}
\]

Letting

\[
\rho\uparrow
\frac{1-\sqrt{\theta_0/\theta_1}}2
\]

gives the optimized two-level rate

\[
\boxed{
R_{\mathrm{pack}}(\theta_1\to\theta_0)
\le
1-h_2\!\left(
\frac{1-\sqrt{\theta_0/\theta_1}}2
\right).}
\tag{8}
\]

## 3. Subcritical packing at the restriction ratio

Set

\[
\frac{\theta_0}{\theta_1}=\alpha^{3/2}.
\]

The binary-entropy inequality

\[
h_2\!\left(\frac{1-u}{2}\right)\ge1-u^2
\]

shows that

\[
1-h_2\!\left(\frac{1-\alpha^{3/4}}2\right)
\le\alpha^{3/2}<\alpha.
\tag{9}
\]

Thus the family at level \(\theta_1\), after quotienting by Hamming clusters
visible at the lower level \(\theta_0=\alpha^{3/2}\theta_1\), has packing rate
strictly below \(\alpha\):

\[
\boxed{
R_{\mathrm{pack}}(\theta_1\to
\alpha^{3/2}\theta_1)<\alpha.}
\tag{10}
\]

In particular, the exact maximizers have fewer than \(2^{(\alpha-o(1))n}\)
clusters at projective separation approximately
\((1-\alpha^{3/4})n\), when clusters are detected in the
\(\alpha^{3/2}K\) layer.

This is the desired subcritical *cluster count*. It is stronger than raw
layer entropy, which is always supercritical.

## 4. Why the packing theorem does not yet imply restriction

A maximal code with separation \(2\rho n\) covers the higher layer by
Hamming balls of radius \(2\rho n\), not \(\rho n\). At the optimized
restriction ratio that covering radius is

\[
(1-\alpha^{3/4})n.
\tag{11}
\]

This is linear and typically comparable to, or larger than, the complement
size \((1-\alpha)n\). A cylinder of extensions of a fixed restricted witness
can therefore have its entire \(2^{(1-\alpha)n}\) points inside the union of
such clusters. Merely knowing that the number of centers is
\(2^{\beta n}\) with \(\beta<\alpha\) does not control those cylinder
intersections.

The small-cut dictionary alone also gives no rigidity. If \(S\) is a uniform
\(r\)-subset, (2) and symmetry give

\[
\mathbb E c_B(S)
=\frac{2r(n-r)}{n(n-1)}K.
\tag{12}
\]

Thus a positive fraction of every appropriate Hamming sphere consists of
small cuts by the elementary mean bound. Equivalently, the universal noise
cloud already produces exponentially many separated small cuts at any
modest separation: greedily packing a positive-density subset of a
radius-\(\rho n\) sphere at minimum distance \(2\varepsilon n\) gives rate at
least

\[
h_2(\rho)-h_2(2\varepsilon)>0
\quad\text{whenever }2\varepsilon<\rho
\text{ is sufficiently small}.
\tag{13}
\]

So “exponentially many separated small cuts” is not by itself a structural
anomaly.

## 5. The missing cluster-to-restriction lemma

To combine (10) with the multiplicity-weighted restriction criterion, one
needs an *anisotropic* cluster estimate. For each center \(z\) and restricted
support \(S\), it must control

\[
\#\{x:\ x_S=y,\ x\text{ lies in the energy cluster of }z\},
\tag{14}
\]

not merely the ambient Hamming radius of the cluster. An estimate with total
cylinder-intersection exponent below \(1-\alpha\) would finish the
restriction theorem.

Ordinary Hamming packing cannot supply this: its covering balls are too
large, and changing a linear number of spins can alter a quadratic energy by
\(\Theta(n^2)\). The surviving target is therefore a cluster theorem in the
energy/cut metric, for example:

\[
c_B(S\mathbin\triangle T)
\le C\big(c_B(S)+c_B(T)\big)
\tag{15}
\]

on near-minimum cut clusters, or a direct bound on (14). No such inequality
follows from global cut-minimality (2), because signed cut functions are not
submodular.

### Exact no-go for ambient Hamming cylinder counting

The preceding loss is not a poor parameter choice. Suppose packing spheres
have radius \(qn\). The center-rate bound is \(1-h_2(q)\), and maximality
covers by balls of radius \(2qn\). The largest intersection of one covering
ball with a cylinder fixing \(\alpha n\) coordinates has exponential rate

\[
J_\alpha(q)
=(1-\alpha)h_2\!\left(
\min\left\{\frac{2q}{1-\alpha},\frac12\right\}
\right).
\tag{16}
\]

For every \(\alpha\in(0,1)\) and every \(q\),

\[
\boxed{1-h_2(q)+J_\alpha(q)\ge1-\alpha.}
\tag{17}
\]

If \(2q/(1-\alpha)\ge1/2\), this is immediate. Otherwise (17) is equivalent
to

\[
h_2(q)\le
\alpha+(1-\alpha)h_2\!\left(\frac{2q}{1-\alpha}\right).
\tag{18}
\]

There is a direct volume proof: an \(n\)-bit Hamming ball of radius \(qn\)
injects into an arbitrary string on the fixed \(\alpha n\) coordinates
together with a ball of radius \(2qn\) on the remaining coordinates.
Taking exponential rates gives (18).

Thus “number of centers times worst cylinder intersection” can never be
smaller than the full cylinder size \(2^{(1-\alpha)n}\). There is no
parameter range in which the ambient Hamming covering argument proves the
restriction theorem.

### Exact no-go for a generic cut metric

Even after ground-state switching, \(c_B(S\mathbin\triangle T)\) need not
obey a triangle inequality in the two small cut values. The smallest example
has three vertices and

\[
b_{01}=b_{02}=1,\qquad b_{12}=-1.
\]

Every cut is nonnegative, but

\[
c_B(\{1\})=0,\qquad
c_B(\{0,1\})=0,\qquad
c_B(\{0\})=2,
\]

and

\[
\{1\}\mathbin\triangle\{0,1\}=\{0\}.
\]

Hence two exact-small-cut centers can have a maximally non-small symmetric
difference. A cut-energy ball does not give a metric cover without an
additional structural hypothesis.

## 6. Verdict

Quotienting the universal noise cloud does produce a rigorous and
quantitatively subcritical cluster-packing rate, (10). This is genuine
progress beyond raw entropy. It does not by itself complete scale transfer:
the gap is now the geometry of cluster intersections with coordinate
cylinders. Equations (17) and the three-vertex cut example show that neither
plain Hamming coding nor a generic signed-cut metric can provide it. What
remains is an averaged selector/cylinder estimate that uses the energy of the
center, not merely its cluster radius.
