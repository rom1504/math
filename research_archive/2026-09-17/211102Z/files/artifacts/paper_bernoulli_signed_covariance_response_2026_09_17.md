# Actual energy gives an adaptive response discount for bounded-covariance query laws

2026-09-17. Bernoulli-track theorem prompted by the director's actual
ground-code minimax question. It uses the already reconstructed
bounded-spectrum Gaussian-sign comparison, but chooses its covariance
from the QUERY LAW, not as an affine function of the signing A.
Independent audits requested. No external novelty claim is made.

## 1. Finite theorem and exact quantifiers

Let A be a symmetric hollow full signing of order n>=2, and let mu
be ANY probability law on Boolean query words x. Put

```
Sigma=E_mu xx^T,  L=||Sigma||_op,
e=E_mu |H_A(x)|/n^(3/2),
kappa=sqrt(2/pi),
sigma(x)=sign(H_A(x)),
K=E_mu sigma(x)xx^T,  m=E_mu sigma(x),
K0=K-mI,  t=1/[2(L+1)].                            (1)
```

At an energy-zero word the mark sigma may be chosen as either +1 or
-1; this changes none of the energy identities. There is an explicitly
defined centered, EXACTLY ISOTROPIC physical sign-column law nu such
that

```
E_mu E_nu |h dot x|/sqrt(n)
 <= kappa f(4kappa^2 t e^2) + C n^(-1/6)sqrt(log(n+1))
 <= kappa - [kappa^5 e^4/(2(L+1)^2)]
             + C n^(-1/6)sqrt(log(n+1)),             (2)
f(v)=[sqrt(1+v)+sqrt(1-v)]/2.
```

C is an absolute constant, independent of A, mu, L, and n. In
particular, SOME physical sign vector h (depending on mu) obeys the
same upper bound for E_mu|h dot x|/sqrt(n). For fixed e>=c>0 and
L<=L_*<infinity, (2) is a uniform positive leading discount from kappa.

This includes every mu supported on an absolute ground code with
Q(A)>=c n^(3/2), or on a nearcode whose absolute energies are all at
least c n^(3/2). Pointwise support is stronger than needed: the theorem
uses only the average absolute energy in (1).

Crucial minimax scope: nu depends on mu, and (2) is an AVERAGE-query
bound. For a given mu it does not bound max_x E_nu|h dot x| on its
support. It bounds the full code's minimax response only if the
maximizing dual query laws can be taken with a common bounded L.
No such unconditional covariance theorem is asserted here.

## 2. Signed energy supplies a macroscopic covariance statistic

For every real vector u,

```
|u^T K u|<=E_mu (u dot x)^2=u^T Sigma u.
```

Thus ||K||_op<=L, |m|<=1, K0 is hollow, and
||K0||_op<=L+1. Since A is hollow,

```
<A,K0>=<A,K>=2E_mu |H_A(x)|=2e n^(3/2).
```

The full-unit coefficient identity ||A||_F^2=n(n-1) and
Cauchy--Schwarz give

```
||K0||_F^2 >=4e^2 n^2/(n-1)>=4e^2 n.              (3)
```

The diagonal subtraction cannot be omitted when m is nonzero.
For completeness,

```
||K0||_F^2=Tr(K^2)-nm^2
 <= E_(x,y iid mu)(x dot y)^2
 =Tr(Sigma^2)<=Ln.                                 (4)
```

Here Tr(K^2)=E sigma(x)sigma(y)(x dot y)^2, so the first inequality
follows termwise. In particular the argument of f in (2) is at most
kappa^2 L/[2(L+1)]<1, as also follows from the actual covariance
matrices below.

## 3. An exactly isotropic adaptive physical law

Define two correlation matrices and their physical sign vectors by

```
R_+=I+tK0,   R_-=I-tK0,
h_+=sign N(0,R_+),   h_-=sign N(0,R_-).
```

Equation (1) ensures (1/2)I<=R_+,R_-<=(3/2)I, uniformly in L.
Let nu be their equal mixture. Each sign vector is centered. The
Gaussian orthant identity, applied entrywise, gives

```
Cov(h_+)=I+B,  Cov(h_-)=I-B,
B_ii=0,  B_ij=kappa^2 arcsin(t K_ij), i!=j.         (5)
```

Consequently E_nu hh^T=I EXACTLY, not merely asymptotically.
For each Boolean x write v_x=x^T Bx/n. Since both matrices in (5)
are positive semidefinite, |v_x|<=1.

The signed average of v_x has a definite positive lower bound:

```
E_mu sigma(x)v_x
 =Tr(KB)/n
 =(kappa^2/n)sum_(i!=j) K_ij arcsin(tK_ij)
 >=(kappa^2 t/n)||K0||_F^2
 >=4kappa^2 t e^2.                                 (6)
```

The scalar inequality z arcsin(tz)>=t z^2 holds for positive and
negative z. It retains the nonlinear arcsine correction with the
RIGHT sign; no entrywise-smallness or replacement of arcsine by its
linear term is needed. This is why adapting to K rather than A is
useful in this dual calculation.

## 4. The bounded-spectrum comparison used, with its dependency

The complete proof in the archived
[bounded-spectrum Gaussian-sign universality theorem](flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md)
was read and reconstructed in this campaign. Its following special
case is sufficient: if R is a correlation matrix with
(1/2)I<=R<=(3/2)I, h=sign N(0,R), and Y is a Gaussian with the EXACT
covariance of h, then uniformly in Boolean x,

```
|E|h dot x|-E|Y dot x||
 <= C n^(1/3)sqrt(log(n+1)).                        (7)
```

To make the quantitative extraction explicit, use two configurations
with features +x,-x and the pressure

```
Phi(z)=log(exp(lambda x dot z)+exp(-lambda x dot z)).
```

The archived theorem bounds its expectation difference by
C lambda^3 n log(n+1)^(3/2), with an absolute C for these fixed spectral
limits. Divide by lambda and use the log-sum-exp approximation to
absolute value, with total error at most 2log(2)/lambda. Choosing
lambda=n^(-1/3)log(n+1)^(-1/2) proves (7). The result allows arbitrary
off-diagonal entries and does NOT assume entrywise correlations are
O(n^(-1/2)). Its proof uses conditional independent replacement and
a block-local Gaussian Stein-kernel contraction, not a Gaussian
approximation based solely on covariance matching.

Apply (7) separately to R_+ and R_-. Equations (5) give, uniformly
over every Boolean query x,

```
E_nu|h dot x|/sqrt(n)
 =kappa f(v_x)+O(n^(-1/6)sqrt(log(n+1))).            (8)
```

The constant is independent of the query law even if it depends on
n and A. Keeping the correlations a fixed distance from singularity
is essential for this uniform statement; the chosen t does exactly
that without any endpoint limiting argument.

## 5. Concavity converts signed energy into the response discount

The function f is even, decreasing on [0,1], and concave there.
Hence (6), Jensen, and |sigma|=1 yield

```
E_mu f(v_x)=E_mu f(|v_x|)
 <=f(E_mu|v_x|)
 <=f(E_mu sigma(x)v_x)
 <=f(4kappa^2 t e^2).                               (9)
```

Also f(0)=1, f'(0)=0, and for |v|<1,

```
f''(v)=-[(1+v)^(-3/2)+(1-v)^(-3/2)]/8<=-1/4.
```

Thus f(v)<=1-v^2/8 on [-1,1], including the endpoints by continuity.
Combining (8)--(9) gives the first line of (2), and the last inequality
with t from (1) gives its explicit discount. Finally, an average over
physical h is an upper bound for the minimum over physical h, proving
the stated adaptive-vector consequence.

## 6. What this resolves, and what it does not

This is a genuine actual-energy-to-response implication under a
bounded-query-covariance hypothesis. It uses literal full-sign energy,
allows arbitrary signs of that energy under mu, retains the marked
diagonal exactly, and constructs a genuine exactly isotropic sign law.
It is distinct from a fixed-degree energy density and from covariance
laws of the form I+-t A/sqrt(n).

The full minimax game still allows query laws with diverging covariance
norm. A large top eigenvalue alone does not give a constant response
discount: it can carry only a vanishing fraction of total covariance
mass. Spectral trimming would have to preserve a macroscopic signed
energy statistic, rather than simply discarding those modes. This
remaining obligation is not discharged by (2).

The [weighted scope counterexamples](paper_bernoulli_centered_energy_scope_2026_09_17.md)
are compatible with (2): their difficult dual laws have diverging
covariance norms. They therefore do not contradict this conditional
positive theorem. Conversely, (2) does not transfer those examples
to, or rule out a hard query law for, the original full-sign class.

## 7. Independent checks

Both other campaign researchers read the complete theorem and the
archived Gaussian-sign proof and independently returned PASS, including
the absolute error constant's independence from L and mu. The replay
script `computations/paper_bernoulli_2026_09_17_signed_covariance_response.py`
checks exact integer signed-energy/Frobenius identities on actual full
signings and both uniform and weighted ground query laws. It separately
labels PSD/arcsine/concavity numerical diagnostics and does not use them
as a substitute for the scalar comparison proof. Its canonical output is
`tmp/paper_portfolio_2026_09_17/bernoulli/signed_covariance_response.json`.

The director has additionally reconstructed a
[covariance-trimming extension](paper_director_adaptive_energy_response_2026_09_17.md)
that removes the query-covariance hypothesis when A itself has
||A||_op=O(sqrt(n)). Both other tracks and this author independently
audited that stronger canonical theorem PASS. It produces ONE exactly
isotropic, uniformly subGaussian physical law cheap on the ENTIRE
macroscopic high-energy code. Its actual-signing spectral hypothesis
and its remaining slope/outside-code obligations are explicit.

Subsequent exact-collision audit: that bounded-A-operator-norm conclusion
is already implied, with a larger explicit discount, by the campaign's
earlier affine-in-A Gaussian-sign law. Thus the trimming argument is
an alternative structural proof, NOT an improved bounded-A theorem.
The distinction of the present theorem is its bounded QUERY covariance
hypothesis with no bounded operator norm assumed for A.
