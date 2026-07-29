# Bernoulli entropic phase and Franz--Parisi overlap audit

## 1. Setup and status

Let
\[
N=\binom n2,\qquad
v_x=(x_ix_j)_{i<j},\qquad
H_A(x)=A\cdot v_x,
\qquad
G_n(A)=\frac{\max_x|H_A(x)|}{n^{3/2}}.
\]
For
\[
\mathcal V_n=\{\sigma v_x:
x\in\{\pm1\}^n/\{\pm{\bf1}\},\ \sigma\in\{\pm1\}\},
\]
we have \(|\mathcal V_n|=2^n\) and
\[
G_n(A)\le c
\quad\Longleftrightarrow\quad
A\cdot u\le cn^{3/2}\quad\text{for every }u\in\mathcal V_n.
\]

This audit produced a rigorous planted speed-\(n^2\) entropy envelope
and exact fixed-replica overlap kernels.  It did not produce a signing
below \(1/2\), nor a speed-\(n^2\) lower-tail limit.  The calculation
also explains precisely why fixed-replica/ROM formulas do not decide
the probability of the zero-violation event.

## 2. A rigorous planted Hamming-shell entropy envelope

Suppose \(A^0\) is a signing with
\[
G_n(A^0)\le c_0.
\]
Let \(\xi_e\) be independent signs with
\[
\mathbb P(\xi_e=-1)=\delta,\qquad
0\le\delta\le\frac12,
\]
and put
\[
B_e=A^0_e\xi_e,\qquad m=\mathbb E\xi_e=1-2\delta.
\]
Then, exactly,
\[
H_B(x)
=mH_{A^0}(x)+Z_x,
\qquad
Z_x=\sum_e A^0_ev_x(e)(\xi_e-m).
\tag{2.1}
\]
The planted noise has covariance
\[
\boxed{
\operatorname{Cov}(Z_x,Z_y)
=4\delta(1-\delta)\langle v_x,v_y\rangle
=2\delta(1-\delta)\bigl((x\cdot y)^2-n\bigr).}
\tag{2.2}
\]
Thus the spin-overlap kernel is quadratic, \(q\mapsto q^2\), with the
finite-\(n\) centering term retained.

### Theorem 2.1

For every fixed \(\delta\in[0,1/2]\) and every \(\varepsilon>0\),
\[
\boxed{
\log Z_n\left(
(1-2\delta)c_0+
2\sqrt{\delta(1-\delta)\log2}+\varepsilon
\right)
\ge
N h(\delta)-o(n^2),}
\tag{2.3}
\]
where
\[
Z_n(c)=\#\{A:G_n(A)\le c\}
\]
and
\[
h(\delta)=-\delta\log\delta-(1-\delta)\log(1-\delta).
\]

### Proof

For a fixed \(x\), the saddle point for a deviation of order
\(n^{3/2}\) in \(Z_x\) is \(O(n^{-1/2})\).  Expanding the exact
Bernoulli cumulant generating function uniformly at that scale gives,
for fixed \(d>0\),
\[
\mathbb P(|Z_x|\ge dn^{3/2})
\le
\exp\left\{
-\frac{d^2}{4\delta(1-\delta)}\,n+o(n)
\right\}.
\tag{2.4}
\]
The endpoint \(\delta=0\) is deterministic.  A union bound over the
\(2^{n-1}\) unoriented cuts and both signs shows that, with probability
tending to one,
\[
\max_x|Z_x|
\le
\left(
2\sqrt{\delta(1-\delta)\log2}+\varepsilon
\right)n^{3/2}.
\tag{2.5}
\]
Equation (2.1) now gives the required energy bound.

With probability tending to one, the Hamming distance
\(d_H(A^0,B)\) is \(\delta N+o(N)\).  Every word in this typical shell
has product-measure mass at most
\[
\exp\{-Nh(\delta)+o(N)\}.
\]
Since a \(1-o(1)\) fraction of the product measure satisfies (2.5),
there are at least \(\exp\{Nh(\delta)-o(N)\}\) distinct satisfactory
signings.  This proves (2.3).

Equivalently, for the uniform-signing lower-tail rate
\[
I_n(c)=-\frac1{n^2}
\log\mathbb P_{\rm unif}(G_n(A)\le c),
\]
(2.3) gives
\[
\limsup I_n(c_\delta+\varepsilon)
\le\frac12\bigl(\log2-h(\delta)\bigr),
\tag{2.6}
\]
where
\[
c_\delta=(1-2\delta)c_0+
2\sqrt{\delta(1-\delta)\log2}.
\]
It also gives the canonical lower bound
\[
\Phi_n(\beta)
\ge
\frac12h(\delta)-\beta c_\delta-o(1).
\tag{2.7}
\]

This entropy cloud cannot improve its seed.  The map
\[
\delta\longmapsto
(1-2\delta)c_0+
2\sqrt{\delta(1-\delta)\log2}
\]
is concave on \([0,1/2]\).  If \(c_0\le1/2\), its minimum is at an
endpoint and equals
\[
\min\{c_0,\sqrt{\log2}\}=c_0.
\tag{2.8}
\]
Thus independent Bernoulli planting creates positive speed-\(n^2\)
entropy above the seed, but no all-orders construction below \(1/2\).

## 3. Exact two-constraint overlap kernel

For \(u=\sigma v_x\), \(w=\tau v_y\), set
\[
q=\frac{x\cdot y}{n}.
\]
Their normalized edge overlap is exactly
\[
\boxed{
\rho_n(u,w)
=\frac{\langle u,w\rangle}{N}
=\sigma\tau\frac{nq^2-1}{n-1}.}
\tag{3.1}
\]
Let
\[
N_+=\frac N2(1+\rho_n),\qquad
N_-=\frac N2(1-\rho_n).
\]
For a uniform signing \(A\), there are independent Rademacher sums
\[
S_+\sim\sum_{1}^{N_+}\varepsilon_i,\qquad
S_-\sim\sum_{1}^{N_-}\varepsilon_i
\]
such that
\[
\boxed{
(A\cdot u,A\cdot w)
\overset d=(S_++S_-,S_+-S_-).}
\tag{3.2}
\]
This is an exact finite-\(n\) reduction.

For fixed \(c>0\) and overlap bounded away from the singular
orientation \(\rho=-1\), moderate deviations in (3.2) give
\[
\mathbb P\left(
A\cdot u\ge cn^{3/2},\
A\cdot w\ge cn^{3/2}
\right)
=
\exp\left\{
-\frac{2c^2}{1+\rho_n}\,n+o(n)
\right\}.
\tag{3.3}
\]
The exponent follows by minimizing
\[
\frac{s^2}{(1+\rho)/2}
+\frac{z^2}{(1-\rho)/2}
\]
under \(s\ge c+|z|\); the minimizer is \(z=0,s=c\).

## 4. Violation count and the fixed-replica formula

Define
\[
Y_c(A)=
\#\{u\in\mathcal V_n:A\cdot u\ge cn^{3/2}\}.
\tag{4.1}
\]
Then \(G_n(A)\le c\) is the zero-violation event \(Y_c(A)=0\).

For one constraint,
\[
\mathbb P(A\cdot u\ge cn^{3/2})
=\exp\{-c^2n+o(n)\},
\]
and therefore
\[
\boxed{
\mathbb E Y_c
=\exp\{n(\log2-c^2)+o(n)\}.}
\tag{4.2}
\]

For two positively oriented directions with spin overlap \(q\), the
normalized contribution to the second moment has exponent
\[
\Delta_c(q)
=
h\left(\frac{1-q}{2}\right)-\log2
+\frac{2c^2q^2}{1+q^2}.
\tag{4.3}
\]
Pinsker's inequality gives
\[
\log2-h\left(\frac{1-q}{2}\right)\ge\frac{q^2}{2},
\]
while, for \(c\le1/2\),
\[
\frac{2c^2q^2}{1+q^2}\le\frac{q^2}{2}.
\]
Hence
\[
\boxed{\Delta_c(q)\le0\qquad(0\le q\le1,\ c\le1/2).}
\tag{4.4}
\]
Using (3.2) at the endpoint overlaps and the method of types elsewhere
gives the rigorous exponential-scale statement
\[
\boxed{
\lim_{n\to\infty}\frac1n
\log\frac{\mathbb E Y_c^2}{(\mathbb E Y_c)^2}=0,
\qquad 0<c<\frac12.}
\tag{4.5}
\]
At \(c=1/2\) the same upper exponent is zero, although polynomial
prefactors become critical.

There is an exact fixed-\(k\) extension.  For oriented directions
\(u^1,\ldots,u^k\),
\[
\log\mathbb E_A
\exp\left\{\sum_{a=1}^k\theta_a A\cdot u^a\right\}
=
\sum_{e=1}^N
\log\cosh\left(\sum_{a=1}^k\theta_a u_e^a\right).
\tag{4.6}
\]
If the \(k\) spin vectors have empirical site law
\(\pi\in\mathcal P(\{\pm1\}^k)\), put
\[
q_{ab}=\mathbb E_\pi[s_as_b],
\qquad
R_{ab}=\sigma_a\sigma_b q_{ab}^2
\quad(a\ne b),\qquad R_{aa}=1.
\tag{4.7}
\]
For nonsingular \(R\), define
\[
\mathcal J_c(R)
=\inf_{z\in[c,\infty)^k}z^\top R^{-1}z.
\tag{4.8}
\]
The method of types and fixed-dimensional moderate deviations yield
\[
\boxed{
\lim_{n\to\infty}\frac1n\log\mathbb E(Y_c)_k
=
\sup_{\pi,\sigma}
\left\{H(\pi)-\mathcal J_c(R(\pi,\sigma))\right\},}
\tag{4.9}
\]
with the lower-semicontinuous extension at singular \(R\).  The
independent saddle has
\[
H(\pi)=k\log2,\qquad R=I,\qquad
\mathcal J_c(I)=kc^2,
\]
reproducing \((\mathbb EY_c)^k\).  Formula (4.9) is a rigorous
fixed-replica Franz--Parisi variational problem.

## 5. Why this does not determine the speed-\(n^2\) hole probability

Equations (4.2)--(4.9) operate at speed \(n\): they describe the
exponentially large number of violations in typical disorder.  The
desired event is the atom
\[
Y_c=0,
\]
whose probability, when nonzero, may be as small as
\[
2^{-N}
=\exp\left\{-\left(\frac{\log2}{2}+o(1)\right)n^2\right\}.
\]
Changing a distribution by mass \(e^{-\Theta(n^2)}\) at zero changes
every fixed moment of an \(e^{\Theta(n)}\)-sized random variable by a
negligible relative amount.  Therefore no fixed collection of
factorial moments, even when (4.9) is solved exactly, can decide the
sparse phase.

The naive independent-constraint prediction would be
\[
\mathbb P(Y_c=0)\approx
\exp\{-\mathbb EY_c\}
=\exp\{-e^{\,n(\log2-c^2)+o(n)}\},
\tag{5.1}
\]
which is doubly exponential for \(c<\sqrt{\log2}\).  At the conference
threshold, known sparse constructions have only
\(e^{-\Theta(n^2)}\) uniform probability, so (5.1) misses the dominant
phase by an exponential-of-an-exponential factor.  Equation (5.1) is
therefore a heuristic to reject, not a bound.

## 6. Gaussian replacement is valid only at fixed replica number

At the fixed-\(k\) saddle, \(\theta_a=O(n^{-1/2})\), and
\[
\log\cosh z=\frac{z^2}{2}-\frac{z^4}{12}+O(z^6).
\]
For fixed \(k\), summing the quartic error over \(N=\Theta(n^2)\)
edges changes (4.6) by only \(O(1)\).  This justifies the Gaussian
quadratic rate at speed \(n\) in (4.9).

To resolve \(Y_c=0\), however, one needs growing-order
inclusion--exclusion or negative replicas.  Already at \(k=\Theta(n)\),
the combined edge field
\[
\sum_{a=1}^k\theta_a u_e^a
\]
can be order one.  Then
\[
\log\cosh z-\frac{z^2}{2}
\]
is order one per edge and order \(n^2\) in total.  The Bernoulli
fourth and higher cumulants are therefore leading at precisely the
speed relevant to the entropic decision.

There is also a structural mismatch: Bernoulli disorder has the fixed
radius \(\|A\|_2^2=N\), whereas Gaussian disorder can lower its support
function by a radial variance-shrink channel at speed \(n^2\).
Consequently a ROM/Gaussian Franz--Parisi value cannot be transferred
to Bernoulli signings without a separate fixed-radius universality
theorem whose error is \(o(n^2)\).

## 7. Verdict

The rigorous conclusions are:

1. a seed at \(c_0\) generates the positive-entropy envelope
   (2.3), but the envelope never improves a seed \(c_0\le1/2\);
2. the exact cut overlap kernel is
   \((nq^2-1)/(n-1)\);
3. all fixed-replica moments are replica-symmetric at exponential
   scale for \(c<1/2\), with the exact variational formula (4.9);
4. this fixed-replica stability says nothing about the
   \(e^{-\Theta(n^2)}\) zero-violation atom;
5. Gaussian replacement errors become leading when the replica order
   grows enough to see that atom.

A successful entropic proof now needs a genuinely Bernoulli,
all-orders object: either a uniform control of the full
inclusion--exclusion series for \(Y_c=0\), or a direct speed-\(n^2\)
LDP on an extremal-support state which retains sparse switching
orbits.  Neither follows from the ROM fixed-replica calculation.
