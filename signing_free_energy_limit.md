# Canonical signing pressure: exact inequalities and the missing contraction

## 1. Normalization and the zero-temperature squeeze

For a symmetric zero-diagonal signing \(A\) of \(K_n\), write

\[
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
M(A)=\max_{x\in\{\pm1\}^n}|H_A(x)|,
\qquad M_n=\min_A M(A).
\]

Put \(E_n=\binom n2\) and

\[
\mathfrak Z_n(\beta)
=\sum_{A\in\{\pm1\}^{E_n}}
  \exp[-\beta\sqrt n\,M(A)],
\qquad
\Phi_n(\beta)=\frac1{n^2}\log\mathfrak Z_n(\beta).
\]

If \(c_n=M_n/n^{3/2}\), then

\[
e^{-\beta\sqrt n M_n}
\le \mathfrak Z_n(\beta)
\le 2^{E_n}e^{-\beta\sqrt n M_n}.
\]

Consequently

\[
\boxed{
-\frac{\Phi_n(\beta)}{\beta}
\le c_n
\le-\frac{\Phi_n(\beta)}{\beta}
 \frac{E_n}{n^2}\frac{\log2}{\beta}.
}
\tag{1.1}
\]

Thus, if \(\Phi_n(\beta)\) converges in \(n\) for an unbounded set of
fixed inverse temperatures \(\beta\), then

\[
\limsup_n c_n-\liminf_n c_n\le\frac{\log2}{2\beta}
\]

for every such \(\beta\), and the desired limit exists.  No
interchange of the \(n\to\infty\) and \(\beta\to\infty\) limits is
needed.

It is convenient below to normalize per edge:

\[
\psi_n(\beta)=\frac1{E_n}\log\mathfrak Z_n(\beta).
\]

Since \(E_n/n^2\to1/2\), convergence of \(\psi_n(\beta)\) and of
\(\Phi_n(\beta)\) are equivalent.

## 2. Exact Gibbs--Shearer restriction inequality

Let \(1\le m<n\), let \(S\) be a uniformly random \(m\)-subset of
\([n]\), and put

\[
q_{n,m}=\frac{E_m}{E_n}=\frac{m(m-1)}{n(n-1)}.
\]

For every probability law \(\mu\) on order-\(n\) signings, write
\(\mu_S\) for its marginal on the edges induced by \(S\).  Edge
Shearer gives

\[
\operatorname{Ent}(\mu)
\le q_{n,m}^{-1}\,
  \mathbb E_S\operatorname{Ent}(\mu_S).
\tag{2.1}
\]

Also, for every deterministic signing \(A\),

\[
M(A[S])\le M(A).
\tag{2.2}
\]

To see (2.2), fix a spin on \(S\), extend it by independent uniform
spins on \(S^c\), and average the full energy.  The average is the
energy induced on \(S\), so some extension has at least that absolute
energy.

The Gibbs variational principle, (2.1), and (2.2) now give

\[
\begin{aligned}
\log\mathfrak Z_n(\beta)
&=\sup_\mu\{\operatorname{Ent}(\mu)
       -\beta\sqrt n\,\mathbb E_\mu M(A)\}\\
&\le q_{n,m}^{-1}\,
  \mathbb E_S\left[
    \operatorname{Ent}(\mu_S)
    -\beta q_{n,m}\sqrt n\,
      \mathbb E_{\mu_S}M(A[S])\right].
\end{aligned}
\]

Therefore, with

\[
\beta_{n\to m}
=\beta q_{n,m}\sqrt{\frac nm},
\]

\[
\boxed{
\log\mathfrak Z_n(\beta)
\le q_{n,m}^{-1}
   \log\mathfrak Z_m(\beta_{n\to m}).
}
\tag{2.3}
\]

Equivalently,

\[
\boxed{
\psi_n(\beta)
\le\psi_m\!\left(
  \beta q_{n,m}\sqrt{\frac nm}\right).
}
\tag{2.4}
\]

For \(m/n\to\alpha\in(0,1)\), the smaller inverse temperature is

\[
\beta_{n\to m}
=\beta\alpha^{3/2}+o(1).
\tag{2.5}
\]

This is the exact changing-temperature obstruction.  For upper and
lower subsequential envelopes

\[
U(\beta)=\limsup_n\psi_n(\beta),\qquad
L(\beta)=\liminf_n\psi_n(\beta),
\]

the needed continuity is uniform away from \(\beta=0\).  Indeed, under
the signing Gibbs law,

\[
-\psi_n'(\beta)
=\frac{\sqrt n}{E_n}\mathbb E_\beta M(A)
\le
\frac{\sqrt n\,M_n}{E_n}+\frac{\log2}{\beta}
\le 1+o(1)+\frac{\log2}{\beta},
\tag{2.5a}
\]

where the middle inequality follows by comparing the Gibbs
variational value with a point mass on an optimal signing, and the
last uses the conference upper bound.  Thus (2.4), together with
equicontinuity on compact subsets of \((0,\infty)\), yields only

\[
U(\beta)\le U(\beta\alpha^{3/2}),\qquad
L(\beta)\le L(\beta\alpha^{3/2}),
\tag{2.6}
\]

which is already implied by monotonicity in \(\beta\).  In particular,
the inequality does not mix \(U\) and \(L\), and hence does not force
the pressure to converge.

This failure cannot be repaired by a two-parameter Fekete argument
using (2.4) alone.  To see this, choose

\[
T_n=n^{3/2}\left(c_0+\varepsilon\sin(\log\log(n+n_0))\right)
\]

with \(0<\varepsilon<c_0\), \(n_0\) large, and harmless integer
rounding.  The sequence can be made nondecreasing because its main
derivative is
\((3/2)c_0\sqrt n\), while the oscillatory derivative is
\(O(\sqrt n/\log n)\).  Consider an abstract hereditary signing model
in which every one of the \(2^{E_n}\) order-\(n\) states has energy
\(T_n\), and restriction maps an order-\(n\) state to an
order-\(m\) state.  Its pressure is

\[
\psi_n^{\rm abs}(\beta)
=\log2-\beta\frac{\sqrt n\,T_n}{E_n},
\]

which oscillates for every fixed \(\beta>0\).  Nevertheless (2.4)
holds exactly, since after cancellation of the entropy terms it is
just \(T_n\ge T_m\).  Thus the whole changing-temperature family of
Shearer inequalities is compatible with persistent multiplicative
scale oscillation.

## 3. The exact restriction theorem that would prove convergence

Define the average principal \(m\)-restriction norm

\[
\overline M_m(A)
=\binom nm^{-1}\sum_{|S|=m}M(A[S]).
\]

Suppose one could strengthen (2.2) to the homogeneous contraction

\[
\boxed{
\overline M_m(A)
\le q_{n,m}\sqrt{\frac nm}\,M(A)
}
\tag{3.1}
\]

for every signing in the relevant low-energy Gibbs support.  Since

\[
\kappa_{n,m}
=\frac{\sqrt{m/n}}{q_{n,m}},
\]

(3.1) is exactly \(M(A)\ge\kappa_{n,m}\overline M_m(A)\).
Repeating the Gibbs--Shearer proof then preserves the inverse
temperature:

\[
\log\mathfrak Z_n(\beta)
\le q_{n,m}^{-1}\log\mathfrak Z_m(\beta),
\qquad
\boxed{\psi_n(\beta)\le\psi_m(\beta).}
\tag{3.2}
\]

In particular, if (3.1) held for \(m=n-1\) for all sufficiently large
orders, then \(\psi_n(\beta)\) would be eventually nonincreasing and
bounded below for every fixed \(\beta>0\).  Hence every canonical
pressure would converge, and (1.1) would prove convergence of
\(M_n/n^{3/2}\).

For \(m=n-1\), (3.1) is

\[
\frac1n\sum_{i=1}^nM(A-i)
\le
\frac{n-2}{\sqrt{n(n-1)}}M(A).
\tag{3.3}
\]

Writing \(d_i(A)=M(A)-M(A-i)\), this is equivalent to

\[
\boxed{
\sum_{i=1}^n d_i(A)
\ge
\left(\frac32+\frac{5}{8n}+O(n^{-2})\right)M(A).
}
\tag{3.4}
\]

This identifies the coefficient \(3/2\) in principal deletion as
exactly the coefficient required by fixed-temperature thermodynamics,
not merely a heuristic derivative of \(n^{3/2}\).

An approximate version with an additive defect \(\delta_n\),

\[
\kappa_{n,n-1}\frac1n\sum_iM(A-i)
\le M(A)+\delta_n,
\]

gives a pressure error of order
\(\beta\delta_n/n^{3/2}\) per step.  Direct telescoping therefore
requires

\[
\sum_n\frac{\delta_n}{n^{3/2}}<\infty.
\tag{3.5}
\]

An unspecified \(o(\sqrt n)\) defect is not sufficient.

## 4. Conference matrices falsify the needed contraction

The direction in (3.1) is important: fixed-temperature Shearer needs
an **upper** bound on the average norm of an induced submatrix.
Lower restriction bounds do not help.

The upper bound (3.1) is false even for flat-spectrum
\(O(n^{3/2})\)-energy signings.  Exact conference examples give:

\[
\begin{array}{c|c|c|c}
n&M(A)&M(A-i)\text{ for every }i&
\sum_i d_i(A)/M(A)\\ \hline
6&5&4&6/5\\
14&21&20&2/3\\
18&33&32&6/11.
\end{array}
\tag{4.1}
\]

For example, at order \(14\), (3.3) would require

\[
20\le
\frac{12}{\sqrt{14\cdot13}}\,21
=18.679\ldots,
\]

which fails by a leading \(O(\sqrt n)\) deletion increment.  The
order-\(6\) example is a proved global minimizer; global optimality of
the order-\(14\) and order-\(18\) conference examples is not known.
They nevertheless show that pseudorthogonality, the exact conference
identity, or an operator-norm bootstrap cannot prove (3.1).

At proportional restriction scales the same issue is expected in a
stronger form: a principal compression of a conference matrix behaves
Haar/Wigner-like and can have a larger normalized Boolean norm than
the parent.  The one-vertex examples (4.1) already give an exact
finite obstruction, without needing that asymptotic assertion.

Thus the canonical-pressure route and the inverse low-traffic
principal-deletion route meet at precisely the same missing theorem:
one has to rule out flat deletion **using global minimality**, not
using spectral regularity alone.

## 5. Exact two-block lower inequality and the rectangular pressure

Let \(N=n+m\), and decompose an order-\(N\) signing as

\[
A=\begin{pmatrix}B&C\\C^\top&D\end{pmatrix}.
\]

Put

\[
R(C)=\max_{x\in\{\pm1\}^n,\ y\in\{\pm1\}^m}|x^\top Cy|.
\]

For every \(B,C,D\),

\[
M(A)\le M(B)+M(D)+R(C).
\tag{5.1}
\]

Consequently the partition functions obey the exact lower bound

\[
\boxed{
\begin{aligned}
\mathfrak Z_N(\beta)\ge&
\mathfrak Z_n\!\left(\beta\sqrt{\frac Nn}\right)
\mathfrak Z_m\!\left(\beta\sqrt{\frac Nm}\right)\\
&\times
\sum_{C\in\{\pm1\}^{n\times m}}
e^{-\beta\sqrt N\,R(C)}.
\end{aligned}
}
\tag{5.2}
\]

Thus the natural reverse inequality is not closed in the scalar
signing pressure: it introduces the canonical pressure of the
rectangular Gale--Berlekamp norm.

For a completely explicit but coarse version, if \(C\) is uniform,
the exponential-moment union bound gives

\[
\mathbb E R(C)
\le\sqrt{2nm(N+1)\log2}.
\tag{5.3}
\]

Jensen's inequality applied to \(e^{-\beta\sqrt N R(C)}\) then yields

\[
\boxed{
\begin{aligned}
\log\mathfrak Z_N(\beta)\ge&
\log\mathfrak Z_n\!\left(\beta\sqrt{\frac Nn}\right)
+\log\mathfrak Z_m\!\left(\beta\sqrt{\frac Nm}\right)\\
&+nm\log2
-\beta\sqrt N\sqrt{2nm(N+1)\log2}.
\end{aligned}
}
\tag{5.4}
\]

For \(n=m=N/2\), the last normalized cost tends
\(\beta\sqrt{\log2/2}\).  Combining (5.4) with (2.4) only gives
one-sided bounds on \(U,L\); it does not force \(U=L\).  Replacing
(5.3) by the exact rectangular low-discrepancy pressure would be a
real strengthening, but would still require an anti-alignment theorem
coupling the cross block to the internal near-ground layers in order
to close the original pressure.

## 6. Exact ground-gauge reduction to a linear cut-cone pressure

There is an exact way to remove the maximum from the Boltzmann weight
at the cost of a hard geometric constraint.

For a signing \(B\), put

\[
T(B)=H_B(\mathbf1)=\sum_{i<j}b_{ij},
\qquad
c_B(S)=\sum_{\substack{i\in S\\j\notin S}}b_{ij}.
\]

Let

\[
\mathcal G_n
=\{B:T(B)=M(B)\}.
\]

Since

\[
H_B(x^S)=T(B)-2c_B(S),
\]

one has the exact cut-cone characterization

\[
\boxed{
B\in\mathcal G_n
\iff
T(B)\ge0
\quad\text{and}\quad
0\le c_B(S)\le T(B)
\ \text{for every }S\subset[n].
}
\tag{6.1}
\]

Define the constrained linear partition function

\[
\mathcal Y_n(\beta)
=\sum_{B\in\mathcal G_n}e^{-\beta\sqrt n\,T(B)}.
\tag{6.2}
\]

Let \(\mathcal C_n\) be the augmented cut group, of size \(2^n\), and
for each \(A\) let

\[
g(A)=|\{v\in\mathcal C_n:A\cdot v=M(A)\}|.
\]

Multiplication by \(v\) maps each ground pair \((A,v)\) to a signing
in \(\mathcal G_n\).  Conversely, every \(B\in\mathcal G_n\) has
exactly \(2^n\) inverse pairs \((Bv,v)\).  Therefore

\[
\sum_A g(A)e^{-\beta\sqrt n M(A)}
=2^n\mathcal Y_n(\beta).
\tag{6.3}
\]

Since \(1\le g(A)\le2^{n-1}\),

\[
\boxed{
2\,\mathcal Y_n(\beta)
\le\mathfrak Z_n(\beta)
\le2^n\mathcal Y_n(\beta).
}
\tag{6.4}
\]

In particular,

\[
\frac1{n^2}\log\mathfrak Z_n(\beta)
-\frac1{n^2}\log\mathcal Y_n(\beta)\longrightarrow0.
\tag{6.5}
\]

Thus canonical-pressure convergence is equivalent to convergence of
a partition function whose energy \(T(B)\) is linear.  All of the
difficulty is now isolated in the discrete absolute cut cone
\(0\le c_B(S)\le T(B)\).  Equivalently, one may count, for each
integer \(t\), the \(\{\pm1\}\)-points in this cone with total sum
\(t\), and take their one-dimensional Laplace transform.

This is sharper than an unlabeled energy histogram: the cone remembers
all cut inequalities simultaneously.  It is also not a standard dense
graphon constraint, because competitive signings have
\(T(B)=\Theta(n^{3/2})\), so their edge-density bias is only
\(\Theta(n^{-1/2})\).  Ordinary graphon entropy sees all such
signings as the zero graphon and loses the second-order geometry.

There is a useful equivalent simplicial-complex form.  Let
\(F(B)=\{e:b_e=-1\}\), and let \(\mathscr D_n\) be the supports of the
augmented cut code (cuts and their complements).  Then (6.1) is
equivalent to

\[
\boxed{
F(B)\in\mathcal I_n
\iff
2|F(B)\cap D|\le|D|
\quad\text{for every }D\in\mathscr D_n.
}
\tag{6.6}
\]

Hence \(\mathcal I_n\) is downward closed.  If

\[
I_n(z)=\sum_{F\in\mathcal I_n}z^{|F|},
\]

then

\[
\boxed{
\mathcal Y_n(\beta)
=e^{-\beta\sqrt n E_n}
 I_n(e^{2\beta\sqrt n}).
}
\tag{6.7}
\]

The maximum face size of \(\mathcal I_n\) is the covering radius of
the augmented cut code, so (6.7) is a high-fugacity independence
polynomial whose zero-temperature rank deficit is exactly the
original problem.

This simplicial complex does **not** have the standard structure that
would supply a log-concavity theorem:

* Signed cut functions in (6.1) are not submodular.  For \(n\ge4\), a
  signing with one negative edge belongs to \(\mathcal G_n\), while
  taking the two singleton shores at the endpoints of that edge makes
  the submodularity difference equal to \(-2\).
* \(\mathcal I_n\) is not a matroid, already for \(n=5\).  Label the
  vertices \(0,\ldots,4\), and put
  \[
  F=\{01,02\},\qquad
  G=\{01,03,12\}.
  \]
  Both obey (6.6), and \(|F|<|G|\).  But adding \(03\) to \(F\)
  violates the cut constraint at \(\{0\}\), while adding \(12\)
  violates the complementary-cut constraint at \(\{0,1,2\}\).
  Thus matroid exchange fails.
* The complexes are not hereditary under principal vertex
  restriction.  The same \(G=\{01,03,12\}\) is feasible in \(K_5\),
  but its restriction to the \(K_4\) on \(0,1,2,3\) has total signed
  sum zero and a negative singleton cut, so it is not in
  \(\mathcal I_4\).

Exact enumeration, independently checkable with
`cone_enumerator.cpp`, gives the following coefficient sequences,
listed by increasing total signed sum \(T\):

\[
\begin{array}{c|l}
n& (T:\#\{B\in\mathcal G_n:T(B)=T\})\\ \hline
5&4:60,\ 6:45,\ 8:10,\ 10:1\\
6&5:72,\ 7:360,\ 9:395,\ 11:105,\ 13:15,\ 15:1\\
7&9:7140,\ 11:14742,\ 13:5880,\ 15:1330,\
17:210,\ 19:21,\ 21:1.
\end{array}
\tag{6.8}
\]

These short rank sequences happen to be log-concave, but the explicit
failure of matroid exchange means that Mason/strong-log-concavity
theorems do not apply.  Rank log-concavity at a fixed order would in
any case not supply the missing cross-order comparison.

The cone does have one exact gluing operation.  Let
\(\mathcal G_{n,m}^{\rm rect}\) be the rectangular sign matrices
\(C\) for which

\[
T(C)=\sum_{i,j}c_{ij}
=\max_{x,y}|x^\top Cy|.
\]

If \(B\in\mathcal G_n\), \(D\in\mathcal G_m\), and
\(C\in\mathcal G_{n,m}^{\rm rect}\), then

\[
\begin{pmatrix}B&C\\C^\top&D\end{pmatrix}\in\mathcal G_{n+m},
\tag{6.9}
\]

because every spin energy is bounded in absolute value by
\(T(B)+T(C)+T(D)\), which is attained at the all-one spin.  This gives
a multiplicative lower bound on \(\mathcal Y_{n+m}\), but it
introduces the rectangular cut-cone pressure and changes both internal
temperatures by the factors \(\sqrt{(n+m)/n}\) and
\(\sqrt{(n+m)/m}\).  Conversely, principal blocks of an element of
\(\mathcal G_{n+m}\) need not lie in their smaller cones, as the
\(K_5\to K_4\) example above shows.  Thus (6.9) does not have the
reverse inequality needed for a Fekete theorem.

## 7. A soft-max reformulation and why ordinary replicas do not close

For \(\gamma>0\), define the spin soft maximum

\[
F_{\gamma}(A)
=\frac1\gamma\log
  \sum_{x,\sigma=\pm1}e^{\gamma\sigma H_A(x)}.
\]

Then

\[
M(A)\le F_\gamma(A)
\le M(A)+\frac{(n+1)\log2}{\gamma}.
\tag{7.1}
\]

Taking \(\gamma=\beta\sqrt n/k\) for a fixed positive integer \(k\)
changes \(n^{-2}\log\mathfrak Z_n(\beta)\) by only \(O(k/n)\), and
turns the Boltzmann factor into the negative integer moment

\[
e^{-\beta\sqrt n F_\gamma(A)}
=
\left(
\sum_{x,\sigma}e^{\gamma\sigma H_A(x)}
\right)^{-k}.
\tag{7.2}
\]

Positive replica moments factor over the edge signs, but (7.2) is a
negative moment.  The Gamma integral

\[
S^{-k}
=\frac1{\Gamma(k)}\int_0^\infty
t^{k-1}e^{-tS}\,dt
\]

introduces all replica numbers with alternating expansion and does not
produce a finite positive spin model.  This makes precise why an
ordinary Guerra--Toninelli interpolation for a positive-replica SK
partition function does not directly address the signing pressure.

## 8. Present verdict

The canonical pressure remains a valid and very direct convergence
criterion, but its two natural closure mechanisms are now sharply
classified:

1. Gibbs--Shearer restriction is exact, but follows the
   \(\beta\mapsto\beta\alpha^{3/2}\) characteristic.
2. Fixed-temperature closure is equivalent to a \(3/2\)-coefficient
   average principal-deletion theorem, which flat conference examples
   falsify unless global minimality supplies additional rigidity.
3. Reverse block gluing introduces a rectangular low-discrepancy
   pressure and, without energy-layer anti-alignment, pays a leading
   \(n^2\) cost.
4. Softening the spin maximum yields negative replicas, not a standard
   positive mean-field spin system.
5. Ground-pair double counting converts the problem exactly (up to
   \(e^{O(n)}\)) to a linear high-fugacity pressure on a
   vertex-transitive down-set, but that down-set is neither a matroid,
   submodular, nor hereditary across orders.

Accordingly, scalar pressure subadditivity alone cannot prove the
limit.  A successful pressure proof must retain at least one
non-scalar state: the overlap/cap-intersection law of near-ground
configurations, or an equivalent cross-block anti-alignment profile.
