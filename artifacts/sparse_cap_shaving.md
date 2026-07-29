# Sparse cap shaving: a uniform perturbation lemma and a Paley orbit obstruction

## Status

This note audits the proposal of lowering exceptional Boolean energies by
flipping only \(m=\Theta(n^{3/2})\) edges.

There are two rigorous conclusions.

1.  Sparse biased flips have a very clean uniform law.  With
    \(\delta=\Theta(n^{-1/2})\), their random error is only
    \(O(n^{5/4})=o(n^{3/2})\), **uniformly over every Boolean state**.
    Thus the \(2^n\)-state metric entropy is not the obstruction.
2.  The proposal cannot remove Paley square-wave resonance.  The affine
    quadratic-residue orbit of one resonant state is so diffuse that every
    signing at Hamming distance \(o(n^2)\) from the Paley signing retains
    the same leading cap.  In particular, \(\Theta(n^{3/2})\) carefully
    chosen flips lower the orbit maximum by at most \(O(n)\), rather than
    by \(\Theta(n^{3/2})\).

The positive perturbation lemma also gives a quantitative necessary
convex-balance condition on the near-cap states of every exact minimizer.
That condition is compatible with the Paley orbit and with the earlier
no-go for independent repair.

Throughout,
\[
 E_n=\{\{i,j\}:1\le i<j\le n\},\qquad N=|E_n|=\binom n2,
\]
\[
 z(x)_{ij}=x_ix_j,\qquad
 H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
 M(A)=\max_x|H_A(x)|.
\]

## 1. Uniform sparse biased-flip lemma

Fix a signing \(A=(a_e)_{e\in E_n}\), a vector
\(\phi\in[-1,1]^{E_n}\), \(0\le\lambda\le1\), and
\[
 0<\delta\le \frac1{1+\lambda}.
\]
Independently flip edge \(e\) with probability
\[
 p_e=\delta(1+\lambda a_e\phi_e),
\]
and call the resulting signing \(B\).  Then, for every Boolean state \(x\),
\[
 \boxed{\displaystyle
 \mathbb E H_B(x)
 =(1-2\delta)H_A(x)-2\delta\lambda\langle\phi,z(x)\rangle .
 }\tag{1}
\]

Let
\[
 m_0=\sum_ep_e\le \delta(1+\lambda)N.
\]
For every \(L>0\), Bernstein's inequality and a union bound over the
\(2^{n-1}\) distinct cut vectors give
\[
 \boxed{\displaystyle
 \Pr\!\left[
 \max_x|H_B(x)-\mathbb EH_B(x)|
 >
 \sqrt{8m_0L}+\frac43L
 \right]
 \le 2^n e^{-L}.
 }\tag{2}
\]
For example, take \(L=(n+2)\log2\).  There is then a realization for
which
\[
 \max_x|H_B(x)-\mathbb EH_B(x)|
 =O\!\left(\sqrt{\delta n^3}+n\right).                 \tag{3}
\]
The number of flipped edges is \(O(m_0+n)\) in the same realization,
after enlarging the harmless absolute constant.

In particular, if \(\delta=\kappa n^{-1/2}\), then
\[
 |A\triangle B|=O_\kappa(n^{3/2})
\]
and
\[
 \boxed{\displaystyle
 H_B(x)
 =(1-2\delta)H_A(x)
 -\frac{2\kappa\lambda}{\sqrt n}\langle\phi,z(x)\rangle
 +O_\kappa(n^{5/4})
 }\tag{4}
\]
simultaneously for every \(x\).

### Proof

Writing \(I_e\) for the flip indicator,
\[
 H_B(x)=H_A(x)-2\sum_eI_ea_ez(x)_e.
\]
Since
\[
 p_ea_e=\delta a_e+\delta\lambda\phi_e,
\]
equation (1) follows.  For fixed \(x\), the centered summands
\[
 -2a_ez(x)_e(I_e-p_e)
\]
have absolute value at most \(2\) and total variance at most \(4m_0\).
Bernstein's inequality, followed by the union bound, proves (2).
The flip-count assertion follows from the scalar Chernoff bound.

### Interpretation

At the \(n^{3/2}\) scale, a sparse perturbation implements the
deterministic kernel
\[
 -\frac{2\kappa\lambda}{\sqrt n}\langle\phi,z(x)\rangle,
\]
while its full-cube random error is lower order.  Chaining can improve
constants in (3), but it is not needed to cross the scale gap.  The real
question is whether one can choose \(\phi\) whose cut-kernel has the
correct sign on every dangerous state.

## 2. A conditional cap-shaving theorem

For \(\varepsilon>0\), define the oriented near-cap set
\[
 \mathcal C_\varepsilon(A)
 =
 \left\{\sigma_xz(x):
 |H_A(x)|\ge M(A)-\varepsilon n^{3/2},\
 \sigma_x=\operatorname{sgn}H_A(x)
 \right\}.                                             \tag{5}
\]

Suppose, along a sequence of signings with
\(M(A)\ge c_0n^{3/2}\), that there are fixed
\(\varepsilon,\gamma>0\) and vectors
\(\phi\in[-1,1]^{E_n}\) such that
\[
 \boxed{\displaystyle
 \sigma_x\langle\phi,z(x)\rangle\ge\gamma n^2
 \quad\text{for every }\sigma_xz(x)\in\mathcal C_\varepsilon(A).
 }\tag{6}
\]
Then, for all sufficiently small fixed \(\kappa>0\), the construction
above with \(\delta=\kappa/\sqrt n\) produces a signing \(B\), differing
from \(A\) on \(O(n^{3/2})\) edges, for which
\[
 M(B)\le M(A)-c(\varepsilon,\gamma,\kappa)n^{3/2}.       \tag{7}
\]

Indeed, (6) lowers every near-cap state by
\(2\kappa\lambda\gamma n^{3/2}\), up to the lower-order error (3).
Every state outside the near-cap set has margin
\(\varepsilon n^{3/2}\); its deterministic adverse displacement is at
most
\[
 2\delta\lambda N\le(\kappa\lambda+o(1))n^{3/2}.
\]
Taking \(\kappa\lambda<\varepsilon/2\), and then taking \(\kappa\)
small enough that near-cap energies do not cross zero, proves (7).

This is the exact regime in which sparse cap shaving works: the oriented
cap set must admit a common macroscopic cut-kernel separator.

## 3. Consequence for exact minimizers: cap convex balance

Let \(A\) be an exact order-\(n\) minimizer, and use any fixed
finite-\(n\) lower bound \(M(A)\ge c_0n^{3/2}\) (changing \(c_0\) by a
harmless \(o(1)\) if desired).  Fix \(0<\varepsilon<c_0/4\), and set
\[
 d_\varepsilon(A)
 =
 \operatorname{dist}_{\ell_1}\!
 \left(0,\operatorname{conv}\mathcal C_\varepsilon(A)\right).
 \tag{8}
\]
Finite-dimensional minimax duality gives
\[
 d_\varepsilon(A)
 =
 \max_{\|\phi\|_\infty\le1}
 \min_{v\in\mathcal C_\varepsilon(A)}
 \langle\phi,v\rangle .                                \tag{9}
\]
Applying (1)--(3) to a maximizing separator in (9), with for example
\(\delta=\varepsilon/(16\sqrt n)\), shows
\[
 \boxed{\displaystyle
 d_\varepsilon(A)
 =
 O\!\left(
 \varepsilon^{-1/2}n^{7/4}
 \varepsilon^{-1}n^{3/2}
 \right).
 }\tag{10}
\]
The same conclusion holds for an \(o(n^{3/2})\)-near-minimizer.

Equivalently, there is a probability measure \(\mu\) on the
\(\varepsilon n^{3/2}\)-near-cap states such that
\[
 \boxed{\displaystyle
 \left\|
 \mathbb E_\mu[\sigma_xz(x)]
 \right\|_1
 =
 O_\varepsilon(n^{7/4}).
 }\tag{11}
\]

To verify (10), if the right side of (9) were \(d\), all near-cap
expected energies would fall by at least \(2\delta d\).  The other
states keep a fixed fraction of their \(\varepsilon n^{3/2}\) margin.
The uniform rounding error is
\[
 O(\sqrt{\delta n^3}+n).
\]
Exact minimality therefore forces
\[
 \delta d
 \lesssim
 \sqrt{\delta n^3}+n,
\]
which gives (10).

This is the precise audit against the sparse-repair no-go.  Exact
minimizers are not expected to possess a separator of the form (6);
instead their oriented cap states must balance as in (11).  Notice that
(11) does not force exponential cap entropy: a polynomial-size
Hadamard-like family can already have almost vanishing pair
correlations.

## 4. Paley affine-orbit obstruction

Let \(p\equiv1\pmod4\) be prime and let
\[
 a_{ij}=\chi_p(i-j),\qquad i,j\in\mathbb F_p,
\]
be the Paley signing.  Write
\[
 N=\binom p2.
\]
The \(+1\) and \(-1\) edge classes both have size \(N/2\).  The affine
group
\[
 \Gamma
 =
 \{j\mapsto rj+b:r\in(\mathbb F_p^\times)^2,\ b\in\mathbb F_p\}
\]
preserves \(A\) and is transitive on each of the two edge classes.

Fix any \(x\in\{\pm1\}^p\), and put
\[
 E=H_A(x),\qquad
 T=\sum_{i<j}x_ix_j
 =\frac{(\sum_ix_i)^2-p}{2}.                            \tag{12}
\]
Let \(B\) be obtained by flipping a set \(S\) of Paley edges, with
\[
 m_+=|\{e\in S:a_e=+1\}|,\quad
 m_-=|\{e\in S:a_e=-1\}|,\quad
 m=m_++m_-.
\]
Then the orbit average is the exact identity
\[
 \boxed{\displaystyle
 \mathbb E_{g\in\Gamma}H_B(gx)
 =
 E-\frac2N\left[mE+(m_+-m_-)T\right].
 }\tag{13}
\]
Consequently,
\[
 \boxed{\displaystyle
 \max_y H_B(y)
 \ge
 E-\frac{2m(E+|T|)}N
 \qquad(E\ge0).
 }\tag{14}
\]

### Proof

Let
\[
 U_+=\sum_{a_e=+1}z(x)_e,\qquad
 U_-=\sum_{a_e=-1}z(x)_e.
\]
Then
\[
 U_+-U_-=E,\qquad U_++U_-=T.
\]
Transitivity on the two edge classes shows that, for every positive
edge \(e\) and negative edge \(f\),
\[
 \mathbb E_gz(gx)_e=\frac{T+E}{N},\qquad
 \mathbb E_gz(gx)_f=\frac{T-E}{N}.                      \tag{15}
\]
Since
\[
 H_B(gx)=E-2\sum_{e\in S}a_ez(gx)_e,
\]
substitution of (15) proves (13), and the maximum is at least the
orbit average.

### Square-wave consequence

For the standard Paley square wave
\[
 x_j=\operatorname{sgn}\cos(2\pi j/p),
\]
one has
\[
 \sum_jx_j=1,\qquad T=\frac{1-p}{2}.                    \tag{16}
\]
Along the resonant primes, its half-energy satisfies
\[
 E=\left(\frac12-o(1)\right)p^{3/2}.
\]
Therefore every signing \(B\) differing from \(A\) on \(m\) edges
satisfies
\[
 \boxed{\displaystyle
 M(B)\ge
 \left(\frac12-o(1)\right)p^{3/2}
 -O\!\left(\frac m{\sqrt p}+\frac mp\right).
 }\tag{17}
\]
In particular:

* if \(m=\Theta(p^{3/2})\), the possible cap reduction in (17) is only
  \(O(p)\);
* if \(m=o(p^2)\), then
  \[
  M(B)\ge\left(\frac12-o(1)\right)p^{3/2}.
  \tag{18}
  \]

More quantitatively, if \(m=(\delta+o(1))N\) for fixed
\(0\le\delta\le1\), then
\[
\boxed{\displaystyle
 \frac{M(B)}{p^{3/2}}\ge \frac12-\delta-o(1).
}\tag{19}
\]
Thus the resonant term \(1/2-\delta\) in the stratified-conference
perturbation calculation is not an artifact of independent noise.  It
is a deterministic floor for **every** choice of the \(\delta N\)
flipped edges.  Uniform random thinning is already optimal, at the
orbit-average level, in its leverage per changed edge.

Thus every vanishing-density Hamming perturbation of a resonant Paley
signing remains resonant at leading order.  A construction that breaks
the Paley \(1/2\) cap on those primes must change \(\Omega(p^2)\) edges.

## 5. Research consequence

Sparse cap shaving is a valid mechanism for an isolated or
strongly clustered exceptional family satisfying (6), and the full-cube
rounding error is already harmless.  It is not a mechanism for repairing
the arithmetic Paley cap:

* the affine orbit supplies polynomially many symmetry copies whose
  average alignment with every edge is only \(O(p^{-1/2})\);
* an \(m\)-edge perturbation therefore has average leverage only
  \(O(m/\sqrt p)\);
* chaining cannot alter this deterministic orbit-average obstruction.

The useful surviving target is structural rather than constructive:
combine the necessary convex balance (11) with energy information to
classify possible cap measures of a true minimizer.  Any future sparse
repair theorem must exploit more than a single common separator, or
must first prove that competitive non-Paley signings cannot support a
Paley-like balanced cap orbit.
