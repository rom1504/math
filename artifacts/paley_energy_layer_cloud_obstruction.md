# Hamming-cloud obstruction to sparse-exception Paley entropy

## Status

This note cross-audits the flat-envelope-with-sparse-exceptions (FESE)
hypothesis in `perturbed_conference_stratified_entropy.md`.

The Bernoulli edge-perturbation calculation and its numerical optimization are
correct.  However, the stated FESE hypothesis is incompatible with any Paley
subsequence having even one fixed positive gap above the flat-spectrum
threshold

\[
r_*=\frac{\sqrt{15}}4=0.968245836551854\ldots .
\]

The obstruction is elementary and rigorous: one resonant Boolean vector
automatically creates exponentially many nearby Boolean vectors whose
quadratic energies remain close to a deterministically attenuated copy of the
seed energy.

The right inverse statement cannot therefore say that the exceptional
**vectors** form a subexponential set.  It must say that the exceptional
vectors lie in subexponentially many structured neighborhoods, and it must
retain the entropy of each such neighborhood.

## 1. General Hamming-cloud theorem

Let \(A_n\) be a symmetric \(n\times n\) matrix with zero diagonal and
\(\pm1\) off-diagonal entries.  Suppose

\[
\|A_nx\|_2^2\le (1+o(1))n^2
\qquad\text{for the seed vectors under consideration}. \tag{1}
\]

This holds with no error for:

- a symmetric conference matrix, since \(A_n^2=(n-1)I\); and
- the \(p\times p\) Paley core, since \(A_p^2=pI-J\).

Fix \(x^{(0)}\in\{\pm1\}^n\), and write its normalized doubled energy as

\[
r_0=\frac{(x^{(0)})^\top A_nx^{(0)}}{n^{3/2}}. \tag{2}
\]

For \(0<\theta<1/2\), independently flip each coordinate of \(x^{(0)}\)
with probability \(\theta\).  Put

\[
\rho=1-2\theta,\qquad
y_i=x_i^{(0)}\eta_i,\qquad
\Pr(\eta_i=-1)=\theta.
\]

Then

\[
\boxed{\displaystyle
\frac{y^\top A_ny}{n^{3/2}}
=\rho^2r_0+o_{\Pr}(1).} \tag{3}
\]

Moreover, for every fixed \(r<\rho^2r_0\),

\[
\boxed{\displaystyle
\#\left\{y\in\{\pm1\}^n:
\frac{y^\top A_ny}{n^{3/2}}\ge r\right\}
\ge
\exp\{n h(\theta)-o(n)\},} \tag{4}
\]

where

\[
h(\theta)=-\theta\log\theta-(1-\theta)\log(1-\theta)
\]

is the binary entropy in natural logarithms.  Passing to the antipodal
quotient changes the count by at most a factor two and does not change its
exponential rate.

### Proof of concentration

Switch \(A_n\) by the seed:

\[
B=D_{x^{(0)}}A_nD_{x^{(0)}}.
\]

Thus \(y^\top A_ny=\eta^\top B\eta\) and
\(\mathbf1^\top B\mathbf1=(x^{(0)})^\top A_nx^{(0)}\).  Write

\[
\eta=\rho\mathbf1+\xi,\qquad
\mathbb E\xi_i=0,\qquad
\sigma^2:=\mathbb E\xi_i^2=1-\rho^2.
\]

Since \(B\) has zero diagonal,

\[
\eta^\top B\eta
=\rho^2\mathbf1^\top B\mathbf1
+2\rho\,\xi^\top B\mathbf1
+\xi^\top B\xi. \tag{5}
\]

Let \(b=B\mathbf1\).  Independence gives

\[
\operatorname{Var}(2\rho\,\xi^\top b)
=4\rho^2\sigma^2\|b\|_2^2. \tag{6}
\]

For a conference matrix,

\[
\|b\|_2^2=\mathbf1^\top B^2\mathbf1=n(n-1).
\]

For the Paley core,

\[
\|b\|_2^2
=(x^{(0)})^\top A_p^2x^{(0)}
=p^2-\left(\sum_i x_i^{(0)}\right)^2
\le p^2. \tag{7}
\]

Also,

\[
\xi^\top B\xi=2\sum_{i<j}b_{ij}\xi_i\xi_j.
\]

Distinct edge monomials are uncorrelated: if two edges are disjoint, this is
immediate, and if they share one endpoint, the other two centered variables
remain unpaired.  Therefore

\[
\operatorname{Var}(\xi^\top B\xi)
=4\sigma^4\sum_{i<j}b_{ij}^2
=2\sigma^4n(n-1). \tag{8}
\]

The covariance between the linear and quadratic terms is zero for the same
unpaired-variable reason.  Consequently,

\[
\operatorname{Var}(\eta^\top B\eta)
\le
\bigl(4\rho^2\sigma^2+2\sigma^4+o(1)\bigr)n^2. \tag{9}
\]

The mean of (5) is
\(\rho^2(x^{(0)})^\top A_nx^{(0)}\).  Chebyshev's inequality at scale
\(n^{3/2}\) proves (3).

### Proof of the exponential count

Let \(K=\#\{i:\eta_i=-1\}\).  With probability \(1-o(1)\),

\[
|K-\theta n|\le n^{2/3}. \tag{10}
\]

By (3), the desired energy event also has probability \(1-o(1)\).  Hence
their intersection has probability \(1-o(1)\).  Every outcome satisfying
(10) has Bernoulli probability

\[
\theta^K(1-\theta)^{n-K}
=\exp\{-nh(\theta)+o(n)\}. \tag{11}
\]

Since the map \(\eta\mapsto y=x^{(0)}\odot\eta\) is one-to-one, an event of
probability \(1-o(1)\) whose atoms each have probability at most
\(\exp\{-nh(\theta)+o(n)\}\) contains at least
\(\exp\{nh(\theta)-o(n)\}\) distinct vectors.  This proves (4).

## 2. Forced energy-layer entropy

Eliminating \(\theta\) from

\[
r=(1-2\theta)^2r_0
\]

gives the lower entropy envelope forced by a seed of energy \(r_0\):

\[
\boxed{\displaystyle
s_{\rm cloud}(r;r_0)
=h\left(\frac{1-\sqrt{r/r_0}}2\right),
\qquad 0<r<r_0.} \tag{12}
\]

Thus a classification by a subexponential number of resonant **templates**
may still be possible, but the union of the associated Hamming clouds can
contain exponentially many resonant vectors.

## 3. Direct contradiction to FESE

Suppose \(r_0>r_*\) along an infinite subsequence.  Choose any fixed
\(\theta>0\) small enough that

\[
(1-2\theta)^2r_0>r_*.
\]

Equation (4) then gives exponentially many vectors above \(r_*\).  They
cannot all lie in an exceptional set of cardinality \(\exp(o(n))\), while
FESE forbids nonexceptional vectors above \(r_*+o(1)\).  Hence:

> **Corollary.** FESE, as stated with a subexponential exceptional vector
> set, fails on every Paley/conference subsequence possessing a resonance
> bounded a fixed amount above \(r_*\).

For a cap seed \(r_0=1\), the largest flip rate whose cloud center remains
above \(r_*\) is

\[
\theta_c
=\frac{1-\sqrt{r_*}}2
=0.00800258218364\ldots ,
\]

and

\[
h(\theta_c)=0.0466068706543\ldots . \tag{13}
\]

So a cap resonance forces at least
\(\exp\{(0.0466068\ldots-o(1))n\}\) configurations above the
flat-spectrum endpoint.

The explicit positive-density square-wave construction recorded in the
campaign gives

\[
r_0=0.968839592155\ldots>r_*.
\]

For that seed,

\[
\theta_c
=\frac{1-\sqrt{r_*/r_0}}2
=0.000153236564383\ldots ,
\]

and

\[
h(\theta_c)=0.00149918242457\ldots . \tag{14}
\]

Even this mild resonance therefore creates a positive-exponential-rate
exceptional cloud.

## 4. Consequence for the edge-perturbation variational calculation

Now independently flip edges at rate \(\delta\), and write

\[
\mu=1-2\delta,\qquad D=\delta(1-\delta).
\]

If a seed of doubled energy \(r_0\) is present, then a layerwise
configuration-counting proof must at least account for its forced Hamming
cloud.  In the same moderate-deviation union-bound framework used in
`perturbed_conference_stratified_entropy.md`, the corresponding branch is

\[
\boxed{\displaystyle
B_{\rm cloud}(\delta;r_0)
=\max_{0\le\theta\le1/2}
\left[
\frac{\mu r_0(1-2\theta)^2}{2}
+2\sqrt{D\,h(\theta)}
\right].} \tag{15}
\]

This is a limitation of that entropy-plus-union-bound proof architecture,
not an unconditional lower bound on the maximum of the randomly perturbed
matrix: the perturbation noises for nearby configurations are correlated.

At the FESE note's claimed optimum

\[
\delta_*=0.00139403918414\ldots ,
\]

a cap seed gives

\[
\max_\theta B_{\rm cloud}(\delta_*;1)
=0.503590895325\ldots ,
\qquad
\theta_{\rm opt}=0.001904\ldots . \tag{16}
\]

Numerically, this cap-cloud branch is already larger than \(1/2\) for every
tested \(\delta>0\); for example,

\[
\begin{array}{c|c}
\delta & B_{\rm cloud}(\delta;1)\\ \hline
10^{-4} & 0.500372873799\ldots\\
10^{-3} & 0.502717846568\ldots\\
\delta_* & 0.503590895325\ldots
\end{array}
\]

By contrast, the mild resonance \(r_0=0.968839592155\ldots\) gives about
\(0.48818\) at \(\delta_*\), so that particular positive-density family
does not by itself destroy the numerical \(0.498606\) target.  It does
destroy the stated sparse-vector hypothesis.

## 5. Corrected target

The perturbation route now needs a structured stratification of the form

\[
\{\pm1\}^p
=\mathcal O_p
\cup
\bigcup_{\tau\in\mathcal T_p}\mathcal N(\tau),
\qquad
|\mathcal T_p|=\exp(o(p)), \tag{17}
\]

where:

1. the ordinary set \(\mathcal O_p\) obeys a flat-spectrum energy-layer
   upper bound;
2. each template neighborhood \(\mathcal N(\tau)\) has an explicit joint
   Hamming-distance/energy entropy profile;
3. the template energy and its cloud profile are jointly controlled; and
4. the supremum of the ordinary branch and all structured-cloud branches
   is strictly below \(1/2\).

The cloud theorem supplies a compulsory lower envelope for item 2.  Any
candidate inverse theorem whose proposed upper entropy lies below (12) is
false.

## 6. Verdict

- **Verified:** the one-configuration Bernoulli edge-perturbation law,
  speed-\(n\) moderate-deviation exponent, and numerical optimization in
  `perturbed_conference_stratified_entropy.md`.
- **Disproved:** FESE with only \(\exp(o(n))\) exceptional vectors, on any
  subsequence with a fixed resonance above \(r_*\).
- **Still viable:** FESE on a carefully selected genuinely nonresonant
  ratio-dense subsequence, or a structured-template theorem retaining the
  full entropy of resonant Hamming clouds.
- **Next concrete question:** construct a ratio-dense Paley/conference
  sequence with no seed \(r_0>r_*+o(1)\), or prove an inverse theorem for
  the joint template-distance/energy profile strong enough to make every
  branch in (17) strictly smaller than \(1/2\).

