# Good-signing entropy and a convergence criterion

## 1. Definitions and switching-orbit quantization

For a signing \(A=(a_{ij})\) of \(K_n\), write

\[
M(A)=\max_{x\in\{\pm1\}^n}
\left|\sum_{i<j}a_{ij}x_ix_j\right|,
\qquad
F_n=\min_A M(A).
\]

Let

\[
\mathcal G_n(t)=\{A:M(A)\le t\},
\qquad
Z_n(t)=|\mathcal G_n(t)|,
\]

and, in normalized notation,

\[
Z_n(c):=Z_n(cn^{3/2}).
\]

Vertex switching and global negation preserve \(M(A)\).  For \(n\ge3\),
the switching orbit has size \(2^{n-1}\), and its global negative is a
disjoint switching orbit: otherwise there would be signs \(s_i\) with
\(s_is_j=-1\) for every pair, which is impossible on a triangle.
Therefore

\[
\boxed{
Z_n(t)=0
\quad\text{or}\quad
Z_n(t)\ge2^n,
}
\tag{1.1}
\]

and in fact every nonzero \(Z_n(t)\) is a multiple of \(2^n\).

This is the first warning about a speed-\(n^2\) entropy.  A single
algebraic switching orbit has

\[
\frac1{n^2}\log Z_n(t)\longrightarrow0,
\]

exactly as the empty set does after replacing \(\log Z\) by
\(\log(1+Z)\).

## 2. Exact Hamming thickening

If two signings differ on \(r\) edges, then for every spin vector their
energies differ by at most \(2r\).  Hence

\[
\boxed{
|M(A)-M(B)|\le2d_H(A,B).
}
\tag{2.1}
\]

Consequently the radius-\(r\) Hamming neighborhood of
\(\mathcal G_n(t)\) is contained in \(\mathcal G_n(t+2r)\).
In particular, if \(Z_n(t)>0\), then

\[
\boxed{
Z_n(t+2r)\ge\sum_{j=0}^r\binom{E_n}{j},
\qquad E_n=\binom n2.
}
\tag{2.2}
\]

For fixed \(\varepsilon>0\), take

\[
r=\left\lfloor\frac{\varepsilon}{2}n^{3/2}\right\rfloor.
\]

Since \(r=o(E_n)\),

\[
\log\binom{E_n}{r}
=r\log\frac{E_n}{r}+O(r)
=\left(\frac{\varepsilon}{4}+o(1)\right)
n^{3/2}\log n.
\]

Thus

\[
\boxed{
Z_n(c)>0
\quad\Longrightarrow\quad
\log Z_n(c+\varepsilon)
\ge
\left(\frac{\varepsilon}{4}+o(1)\right)n^{3/2}\log n.
}
\tag{2.3}
\]

A single good signing is therefore invisible at speed \(n^2\), but
after any fixed relaxation of the normalized threshold it creates a
canonical entropy of order \(n^{3/2}\log n\).

## 3. A refined-microcanonical convergence theorem

Define the nonnegative refined entropy

\[
\Sigma_n(c)=
\frac{\log(1+Z_n(c))}{n^{3/2}\log n}.
\tag{3.1}
\]

### Theorem 3.1

Suppose that, for every fixed \(c\) in an interval containing all
subsequential limits of \(F_n/n^{3/2}\), the extended limit

\[
\Sigma(c)=\lim_{n\to\infty}\Sigma_n(c)
\in[0,\infty]
\tag{3.2}
\]

exists.  Then

\[
\boxed{\lim_{n\to\infty}\frac{F_n}{n^{3/2}}\ \text{exists}.}
\tag{3.3}
\]

#### Proof

Put

\[
a=\liminf_n\frac{F_n}{n^{3/2}}.
\]

Fix \(\eta>0\).  Along an infinite subsequence,

\[
Z_n(a+\eta)>0.
\]

Apply (2.3), with relaxation \(\eta\), to obtain

\[
\limsup_n\Sigma_n(a+2\eta)\ge\frac{\eta}{4}.
\]

The assumed existence of the limit in (3.2) implies

\[
\Sigma_n(a+2\eta)>0
\]

for every sufficiently large \(n\).  Since the left side is zero
exactly when \(Z_n(a+2\eta)=0\), this means

\[
\frac{F_n}{n^{3/2}}\le a+2\eta
\]

eventually.  Hence

\[
\limsup_n\frac{F_n}{n^{3/2}}\le a+2\eta.
\]

Letting \(\eta\downarrow0\) proves (3.3).  \(\square\)

This criterion is robust to a sparse algebraic phase.  It does not
assume that one good switching orbit proliferates to
\(\exp(\Theta(n^2))\) signings.

## 4. A canonical free-energy convergence theorem

There is an even cleaner canonical formulation.  Define a Gibbs
partition function on the **space of edge signings**

\[
\mathfrak Z_n(\beta)
=
\sum_{A\in\{\pm1\}^{E_n}}
\exp\!\left[-\beta\sqrt n\,M(A)\right]
\tag{4.1}
\]

and its pressure

\[
\Phi_n(\beta)=\frac1{n^2}\log\mathfrak Z_n(\beta).
\tag{4.2}
\]

The factor \(\sqrt n\) is forced: both the entropy of edge signings and
the Gibbs penalty are then of order \(n^2\).

### Theorem 4.1

Suppose there is an unbounded set of inverse temperatures
\(\mathcal B\subset(0,\infty)\) such that

\[
\Phi(\beta)=\lim_{n\to\infty}\Phi_n(\beta)
\tag{4.3}
\]

exists for every \(\beta\in\mathcal B\).  Then

\[
\boxed{\lim_{n\to\infty}\frac{F_n}{n^{3/2}}\ \text{exists}.}
\tag{4.4}
\]

#### Proof

Put \(f_n=F_n/n^{3/2}\).  The contribution of one minimizing signing,
and then the bound by all \(2^{E_n}\) signings, give

\[
e^{-\beta n^2f_n}
\le\mathfrak Z_n(\beta)
\le2^{E_n}e^{-\beta n^2f_n}.
\]

Therefore

\[
\boxed{
-\frac{\Phi_n(\beta)}{\beta}
\le f_n
\le
-\frac{\Phi_n(\beta)}{\beta}
+\frac{E_n}{n^2}\frac{\log2}{\beta}.
}
\tag{4.5}
\]

If (4.3) holds, then

\[
\limsup_nf_n-\liminf_nf_n
\le\frac{\log2}{2\beta}.
\]

Let \(\beta\to\infty\) through \(\mathcal B\).  The right side tends to
zero, proving (4.4).  \(\square\)

This reduces the original problem to an ordinary-looking
thermodynamic-limit statement.  Importantly, it does not exchange the
limits \(n\to\infty\) and \(\beta\to\infty\): convergence at arbitrarily
large fixed \(\beta\)'s is enough, and the elementary entropy squeeze
(4.5) performs the zero-temperature step.

The canonical and microcanonical objects are related by the exact
Laplace identity

\[
\boxed{
\mathfrak Z_n(\beta)
=
\beta\sqrt n\int_0^\infty
e^{-\beta\sqrt n\,t}Z_n(t)\,dt
=
\beta n^2\int_0^\infty
e^{-\beta c n^2}Z_n(c)\,dc.
}
\tag{4.6}
\]

The second equality uses \(t=cn^{3/2}\).  Thus Theorem 4.1 asks for a
Laplace principle which remains valid at the maximal-rate boundary
where isolated switching orbits live.

### Exact vertex-cavity product

Let \(B\) have order \(n\), and extend it by a row
\(b\in\{\pm1\}^n\).  Flipping the new vertex changes only the linear
term, so

\[
\boxed{
M(B,b)
=
\max_x\left(|H_B(x)|+|b\cdot x|\right).
}
\tag{4.7}
\]

Define the nonnegative cavity increment and row factor

\[
\Delta_B(b)=M(B,b)-M(B),
\qquad
R_{B,n}(\lambda)=\sum_b e^{-\lambda\Delta_B(b)}.
\tag{4.8}
\]

With

\[
\beta'=\beta\sqrt{\frac{n+1}{n}},
\qquad
\lambda=\beta\sqrt{n+1},
\]

the signing-space partition function has the exact product identity

\[
\boxed{
\mathfrak Z_{n+1}(\beta)
=
\mathfrak Z_n(\beta')\,
\mathbb E_{\nu_{n,\beta'}}
R_{B,n}(\lambda),
}
\tag{4.9}
\]

where

\[
\nu_{n,\beta'}(B)
=
\frac{e^{-\beta'\sqrt n\,M(B)}}{\mathfrak Z_n(\beta')}.
\]

Thus the pressure problem has been reduced exactly to a row-cavity
free energy.  A uniform limit theorem for

\[
\frac1n\log
\mathbb E_{\nu_{n,\beta'}}
R_{B,n}(\beta\sqrt{n+1})
\tag{4.10}
\]

together with the harmless temperature displacement
\(\beta'-\beta=O(n^{-1})\), would close the canonical recursion.
The row factor depends on the complete joint profile of energy and
linear correlations, not only on \(M(B)\).

Here is a finite exact witness to that closure failure.  The two
order-six signings

\[
B_1=
\begin{pmatrix}
0&1&-1&-1&1&1\\
1&0&-1&-1&-1&-1\\
-1&-1&0&1&-1&1\\
-1&-1&1&0&1&1\\
1&-1&-1&1&0&-1\\
1&-1&1&1&-1&0
\end{pmatrix},
\]

\[
B_2=
\begin{pmatrix}
0&1&1&1&-1&-1\\
1&0&1&-1&1&1\\
1&1&0&1&-1&1\\
1&-1&1&0&1&1\\
-1&1&-1&1&0&1\\
-1&1&1&1&1&0
\end{pmatrix}
\]

both satisfy \(M(B_1)=M(B_2)=9\), but exhaustive evaluation of the
\(64\) possible rows gives

\[
\boxed{
R_{B_1,6}(\lambda)
=36e^{-2\lambda}+24e^{-4\lambda}+4e^{-6\lambda},
}
\tag{4.11}
\]

while

\[
\boxed{
R_{B_2,6}(\lambda)
=8+40e^{-2\lambda}+14e^{-4\lambda}+2e^{-6\lambda}.
}
\tag{4.12}
\]

In particular, \(B_2\) has eight norm-preserving extensions and
\(B_1\) has none.  Scalar norm entropy cannot replace the cavity
profile in (4.9).

### The smallest exact one-step state

Put

\[
g_B(x)=M(B)-|H_B(x)|
\]

for the absolute energy gap.  Equation (4.7) becomes

\[
\boxed{
\Delta_B(b)
=
\max_x\bigl(|b\cdot x|-g_B(x)\bigr).
}
\tag{4.13}
\]

For \(u\ge0\), define the cumulative row-cavity profile

\[
V_B(u)=|\{b:\Delta_B(b)\le u\}|.
\tag{4.14}
\]

If \(d_\pm(b,x)=\min\{d_H(b,x),d_H(b,-x)\}\), then

\[
|b\cdot x|=n-2d_\pm(b,x),
\]

and hence

\[
\boxed{
\{\Delta_B\le u\}
=
\bigcap_x
\left\{
b:d_\pm(b,x)\ge
\frac{n-g_B(x)-u}{2}
\right\}.
}
\tag{4.15}
\]

Thus \(V_B(u)\) is the complement of a union of Hamming caps centered
at all near-ground configurations, with radii determined by their
gaps.  The ordinary energy-layer histogram gives the cap sizes, but
not their intersection volumes.

The row factor is, up to a vanishing error, the Legendre transform of
this cumulative profile.  At \(\lambda=\beta\sqrt n\),

\[
\boxed{
\begin{aligned}
\max_{0\le d\le n}
\left\{
\frac1n\log V_B(d)-\beta\frac d{\sqrt n}
\right\}
&\le
\frac1n\log R_{B,n}(\beta\sqrt n)\\
&\le
\max_{0\le d\le n}
\left\{
\frac1n\log V_B(d)-\beta\frac d{\sqrt n}
\right\}
+\frac{\log(n+1)}n.
\end{aligned}
}
\tag{4.16}
\]

The lower bound keeps all rows with \(\Delta_B(b)\le d\).  For the
upper bound, group rows by their integer increment, dominate the
number in the \(d\)-th shell by \(V_B(d)\), and sum at most \(n+1\)
shells.

So \(V_B\), rather than \(M(B)\) or the one-dimensional energy
histogram, is the smallest exact state for the **one-step** cavity
factor.  It is not closed under a second extension: that transition
requires the locations of the admissible rows and their overlaps with
the old energy layers.

The failure of the energy histogram can be certified exactly already
at order eight.  The two matrices

\[
C_1=
\begin{pmatrix}
0&1&-1&1&1&1&1&1\\
1&0&1&-1&-1&-1&-1&1\\
-1&1&0&-1&1&-1&1&-1\\
1&-1&-1&0&-1&-1&1&1\\
1&-1&1&-1&0&1&-1&1\\
1&-1&-1&-1&1&0&1&1\\
1&-1&1&1&-1&1&0&1\\
1&1&-1&1&1&1&1&0
\end{pmatrix},
\]

\[
C_2=
\begin{pmatrix}
0&1&1&1&1&1&-1&-1\\
1&0&-1&1&1&-1&1&1\\
1&-1&0&1&1&-1&-1&1\\
1&1&1&0&1&1&1&-1\\
1&1&1&1&0&-1&1&-1\\
1&-1&-1&1&-1&0&-1&1\\
-1&1&-1&1&1&-1&0&1\\
-1&1&1&-1&-1&1&1&0
\end{pmatrix}
\]

have the same norm \(16\) and exactly the same absolute-energy
histogram

\[
\begin{array}{c|rrrrrrrrr}
|H|&0&2&4&6&8&10&12&14&16\\ \hline
\#x&38&76&60&36&24&12&4&4&2.
\end{array}
\tag{4.17}
\]

Nevertheless their row-increment histograms are

\[
\begin{array}{c|rrrrr}
\Delta&0&2&4&6&8\\ \hline
C_1&60&110&66&18&2\\
C_2&68&112&58&16&2.
\end{array}
\tag{4.18}
\]

Thus even the complete unlabeled gap spectrum does not close (4.9).
The missing data begin with overlaps between energy layers; exact cap
unions in (4.15) involve their full higher intersection hierarchy.

Under the Gibbs law \(\nu_{n,\beta'}\), (4.16) shows the precise
remaining probabilistic target.  One needs a speed-\(n\) large
deviation principle, or uniform concentration, for the random
monotone curve

\[
s\longmapsto
\frac1n\log V_B(\lfloor s\sqrt n\rfloor).
\tag{4.19}
\]

An LDP for the scalar norm, or for the energy histogram alone, cannot
determine the annealed cavity factor
\(\mathbb E_{\nu}R_{B,n}\).

### Gauge-fixed cavity identity

The Gibbs law \(\nu_{n,\beta'}\) is invariant under every vertex
switching \(B\mapsto B^s\).  Directly from (4.13),

\[
\Delta_{B^s}(b)=\Delta_B(bs).
\tag{4.20}
\]

Therefore every summand in the annealed row factor has the same
expectation:

\[
\boxed{
\mathbb E_{\nu_{n,\beta'}}R_{B,n}(\lambda)
=
2^n\,
\mathbb E_{\nu_{n,\beta'}}
e^{-\lambda\Delta_B(\mathbf1)}.
}
\tag{4.21}
\]

Equivalently, every order-\((n+1)\) switching class has a unique
representative whose edges incident to the new vertex are all
positive.  If \(T(B)\) denotes that universal-positive extension,
then

\[
\boxed{
\mathfrak Z_{n+1}(\beta)
=
2^n\sum_B
e^{-\beta\sqrt{n+1}\,M(T(B))}.
}
\tag{4.22}
\]

For a fixed gauge, \(\Delta_B(\mathbf1)\) needs only the
energy-versus-magnetization profile.  Define

\[
U_B(m)=\max_{\sum x_i=m}H_B(x),
\qquad
L_B(m)=\min_{\sum x_i=m}H_B(x).
\tag{4.23}
\]

Then

\[
M(T(B))
=
\max_m\max\{U_B(m)+|m|,\,-L_B(m)+|m|\}.
\tag{4.24}
\]

Moreover the signed profile closes under one universal-positive
extension.  If the new total magnetization is \(q\), then

\[
\boxed{
\begin{aligned}
U_{T(B)}(q)
&=
\max_{\substack{s=\pm1\\m=q-s}}
\bigl(U_B(m)+sm\bigr),\\
L_{T(B)}(q)
&=
\min_{\substack{s=\pm1\\m=q-s}}
\bigl(L_B(m)+sm\bigr).
\end{aligned}
}
\tag{4.25}
\]

This is a genuine max-plus/min-plus transition.  It does not by itself
give a projective Markov chain: at the next step one must average over
all gauges of \(B\), and switching changes which linear functional is
called magnetization.  The orbit of the profile (4.23) under all
gauges is exactly the cap-overlap data in (4.15).

### Concentration audit

If \(B,B'\) differ in one core edge, then

\[
|g_B(x)-g_{B'}(x)|\le4,
\qquad
|\Delta_B(b)-\Delta_{B'}(b)|\le4,
\]

and hence

\[
\boxed{
|\log R_{B,n}(\beta\sqrt n)
-\log R_{B',n}(\beta\sqrt n)|
\le4\beta\sqrt n.
}
\tag{4.26}
\]

Under the uniform product measure, Efron--Stein or bounded differences
therefore gives only

\[
\operatorname{Var}(\log R_{B,n})
\le O(\beta^2n^3).
\tag{4.27}
\]

Since \(\log R\) itself has natural scale \(n\), this bound is worse
than the variance needed for self-averaging by a factor \(n\).  For
the normalized cavity pressure \(n^{-1}\log R\), the bounded-difference
variance proxy is \(O(\beta^2n)\), not \(o(1)\).

Under \(\nu_{n,\beta}\) the situation is strictly weaker.  Flipping one
edge changes the Gibbs Hamiltonian by at most \(2\beta\sqrt n\), so
single-edge conditional odds can be as large as
\(e^{2\beta\sqrt n}\); changing another edge can change their log-odds
by as much as \(4\beta\sqrt n\).  No dimension-free Dobrushin or product
log-Sobolev estimate follows; the Gibbs measure may be essentially
frozen in individual edge coordinates.  Switching invariance removes
\(n-1\) gauge directions but leaves \(\Theta(n^2)\) cycle directions.

Thus ordinary edgewise concentration does not prove the required
speed-\(n\) LDP for (4.19).  A successful self-averaging theorem must
use the geometry of the Hamming-cap intersection, not only
edge-Lipschitzness.

## 5. Why a speed-\(n^2\) lower-tail LDP alone is not enough

For a uniformly random signing let

\[
p_n(c)=2^{-E_n}Z_n(c),
\qquad
I_n(c)=-\frac1{n^2}\log p_n(c),
\]

with \(I_n(c)=+\infty\) when \(Z_n(c)=0\).  If the event is nonempty,
(1.1) implies

\[
\boxed{
I_n(c)
\le
\frac{E_n-n}{n^2}\log2
=\frac{\log2}{2}+O(n^{-1}).
}
\tag{5.1}
\]

Thus a rate strictly larger than \((\log2)/2\) certifies eventual
nonexistence, while a finite exact limit of \(I_n(c)\) certifies
eventual existence.  But a sparse family consisting of one or
subexponentially many switching orbits has

\[
I_n(c)\longrightarrow\frac{\log2}{2}.
\tag{5.2}
\]

It lies at the maximal possible finite rate.  The leading
speed-\(n^2\) rate therefore cannot distinguish:

* no good signing;
* one algebraic switching orbit;
* \(\exp(o(n^2))\) good signings.

This is the precise obstruction to using a Bernoulli lower-tail or
Franz--Parisi calculation by itself to settle the minimum.  Such a
calculation locates the onset of a positive bulk entropy, but it must
be supplemented either by exclusion of a maximal-rate sparse phase or
by the refined entropy/free-energy criteria in Sections 3--4.

## 6. Exact restriction and extension counting

The counting problem has a useful local recursion.  If \(B\) is a
principal submatrix of \(A\), then

\[
M(B)\le M(A):
\]

for a fixed spin on \(B\), average the full energy over independent
spins on the deleted vertices.  Conversely, adjoining one vertex can
increase the norm by at most \(n\).  Every one of the \(2^n\) possible
new rows therefore obeys

\[
M(B)\le M(A)\le M(B)+n.
\]

It follows that

\[
\boxed{
2^n Z_n(t-n)
\le Z_{n+1}(t)
\le2^n Z_n(t).
}
\tag{6.1}
\]

Equivalently, for the random-signing probabilities,

\[
\boxed{
p_n(t-n)\le p_{n+1}(t)\le p_n(t).
}
\tag{6.2}
\]

More generally, if \(N=n+h\) and

\[
L=E_N-E_n
\]

is the number of new edges, then

\[
p_n(t-L)\le p_N(t)\le p_n(t).
\tag{6.3}
\]

The deterministic lower bound in (6.3) can be sharpened
probabilistically.  Fix a core signing \(B\) of order \(n\), and fill
all \(L\) new edges independently.  For a fixed full spin vector, the
new-edge contribution is a sum of \(L\) independent signs.  A union
bound over all full spin vectors shows that, for

\[
u=
\sqrt{2L\bigl((N+1)\log2+s\bigr)},
\tag{6.4}
\]

at least a fraction \(1-e^{-s}\) of all extensions satisfy

\[
M(A)\le M(B)+u.
\]

Summing over all good cores yields

\[
\boxed{
Z_N(t+u)
\ge
(1-e^{-s})\,2^L Z_n(t),
}
\tag{6.5}
\]

or

\[
\boxed{
p_N(t+u)\ge(1-e^{-s})p_n(t).
}
\tag{6.6}
\]

If \(n/N\to\alpha\in(0,1)\) and \(s=o(N)\), the normalized threshold
inflation in (6.4) is

\[
\boxed{
\frac{u}{N^{3/2}}
\longrightarrow
\sqrt{\log2\,(1-\alpha^2)}.
}
\tag{6.7}
\]

This proves local continuity and a genuine entropy-transfer theorem.
It does not preserve a candidate constant under proportional growth:
the error in (6.7) is of leading order whenever \(\alpha<1\).

### Hereditary Shearer inequality

There is a stronger exact count inequality at every fixed raw
threshold.  Let \(N\ge m\ge2\), and choose a signing uniformly from
\(\mathcal G_N(t)\).  Its restriction to every \(m\)-vertex set lies
in \(\mathcal G_m(t)\).  Each edge coordinate occurs in

\[
r=\binom{N-2}{m-2}
\]

of these restrictions.  Shearer's entropy inequality gives

\[
r\log Z_N(t)
\le
\binom Nm\log Z_m(t).
\]

Since

\[
\frac{\binom Nm}{\binom{N-2}{m-2}}
=\frac{E_N}{E_m},
\]

one obtains

\[
\boxed{
\frac{\log Z_N(t)}{E_N}
\le
\frac{\log Z_m(t)}{E_m}.
}
\tag{6.8}
\]

Thus the good-signing entropy per edge is monotone in the order at a
fixed **raw** threshold \(t\).

On the desired diagonal, however, (6.8) reads

\[
\frac{\log Z_N(cN^{3/2})}{E_N}
\le
\frac{
\log Z_m\!\left(
c(N/m)^{3/2}m^{3/2}
\right)
}{E_m}.
\tag{6.9}
\]

The normalized threshold on the right is inflated by
\((N/m)^{3/2}\).  Even for \(m=N-1\), the raw threshold displacement is
\((3c/2+o(1))\sqrt N\).  Hamming thickening controls entropy growth
*from below* when the threshold is relaxed; it supplies no upper bound
on the number of unrelated signings which can enter during this
displacement.  Hence it cannot turn (6.9) into diagonal monotonicity.

The failure is logical, not just a missing estimate.  Choose, for all
large \(n\),

\[
T_n=n^{3/2}
\left(c_0+\varepsilon\sin(\log\log n)\right),
\tag{6.10}
\]

with \(0<\varepsilon\ll c_0\).  After modifying finitely many terms,
\(T_n\) is nondecreasing and

\[
0\le T_{n+1}-T_n\le n.
\]

Define the abstract count profile

\[
\widehat Z_n(t)=
\begin{cases}
0,&t<T_n,\\
2^{E_n},&t\ge T_n.
\end{cases}
\tag{6.11}
\]

It satisfies the fixed-threshold Shearer inequality (6.8), the exact
one-step scalar sandwich (6.1), switching-orbit quantization, and every
one-sided Hamming-thickening lower bound (2.2).  Nevertheless its
normalized feasibility threshold is

\[
\frac{T_n}{n^{3/2}}
=c_0+\varepsilon\sin(\log\log n),
\]

which does not converge.  Therefore these scalar hereditary and
regularization inequalities, even taken together, cannot prove the
desired limit.

The same model is compatible with the scalar cavity recursion.  Give
every abstract signing of order \(n\) the energy \(T_n\).  Its row
factor is then

\[
\widehat R_n(\lambda)
=2^n e^{-\lambda(T_{n+1}-T_n)}.
\tag{6.12}
\]

Because \(0\le T_{n+1}-T_n\le n\),

\[
2^ne^{-\lambda n}
\le\widehat R_n(\lambda)\le2^n,
\]

exactly the universal range allowed by a true vertex extension.
Hence a bound which uses only the size of the cavity factor and the
scalar norm increment cannot exclude slow oscillation either.

There is a canonical analogue.  If

\[
\overline{\mathfrak Z}_n(\lambda)
=2^{-E_n}\sum_Ae^{-\lambda M(A)},
\]

then restriction coupling, or generalized Hölder, gives

\[
\overline{\mathfrak Z}_N(\lambda)
\le\overline{\mathfrak Z}_m(\lambda).
\tag{6.13}
\]

A full Shearer variational argument also gives

\[
\boxed{
\frac{\log\mathfrak Z_N(\lambda)}{E_N}
\le
\frac{
\log\mathfrak Z_m(\lambda E_m/E_N)
}{E_m}.
}
\tag{6.14}
\]

To prove (6.14), apply entropy Shearer to an arbitrary law \(\mu\) on
order-\(N\) signings.  Since every principal norm is at most the full
norm,

\[
\sum_{|S|=m}\mathbb E_\mu M(A[S])
\le\binom Nm\mathbb E_\mu M(A).
\]

Combining this with the entropy inequality and the Gibbs variational
principle yields (6.14).  At the natural temperature
\(\lambda=\beta\sqrt N\), the smaller-order inverse temperature is

\[
\beta_m
=
\beta\left(\frac mN\right)^{3/2}(1+o(1)).
\tag{6.15}
\]

Thus canonical Shearer follows a changing-temperature characteristic,
not the fixed-\(\beta\) line required by Theorem 4.1.  The abstract
profile (6.10)--(6.11) also makes the corresponding canonical pressure
oscillate, so the characteristic inequality alone does not close the
gap.

## 7. The exact block obstruction to supermultiplicativity

For an order-\(n\) signing \(A\), an order-\(m\) signing \(D\), and an
\(n\times m\) cross signing \(B\), put

\[
\mathcal A=
\begin{pmatrix}
A&B\\
B^\top&D
\end{pmatrix}.
\]

In the half-quadratic normalization used here,

\[
\boxed{
M(\mathcal A)
=
\max_{x,y}
\left(
|H_A(x)+H_D(y)|+|x^\top By|
\right).
}
\tag{7.1}
\]

Indeed, replacing \(y\) by \(-y\) leaves both internal energies
unchanged and reverses the cross term.

Consequently cross edges can never cancel internal energy.  A
supermultiplicative estimate for \(Z_{n+m}(c)\) at the same normalized
constant would have to count cross blocks \(B\) whose large bilinear
values systematically avoid all pairs of high-energy layers of
\(A,D\).  Scalar data \(M(A),M(D)\), or even their separate
microcanonical entropies, do not encode this anti-alignment.

The exact missing counting inequality can be stated as follows.  For
\(\varepsilon_n\downarrow0\), one would need a scale-transfer estimate
of the form

\[
\boxed{
Z_{kn}\!\left((c+\varepsilon_n)(kn)^{3/2}\right)
\ge
\exp[-o((kn)^2)]\,
\mathcal F_{n,k}\bigl(Z_n(cn^{3/2})\bigr),
}
\tag{7.2}
\]

where \(\mathcal F_{n,k}(z)>0\) whenever \(z>0\), uniformly in \(k\),
or an analogous inequality for \(\mathfrak Z_n(\beta)\).  Constructing
the compatible cross blocks required by (7.1) is exactly the missing
amplification theorem; counting does not remove it automatically.

At finite replica depth this obstruction remains.  For fixed spin
vectors \(x^{(1)},\ldots,x^{(r)}\), the joint signing moment generating
function is

\[
\boxed{
\mathbb E_A
\exp\!\left(
\sum_{\ell=1}^r\theta_\ell H_A(x^{(\ell)})
\right)
=
\prod_{i<j}
\cosh\!\left(
\sum_{\ell=1}^r
\theta_\ell x_i^{(\ell)}x_j^{(\ell)}
\right).
}
\tag{7.3}
\]

It is determined by the \(2^r\) vertex overlap types.  But the event
defining \(Z_n(c)\) retains all \(2^{n-1}\) distinct quadratic
constraints.  Any fixed-replica interpolation can locate a bulk
lower-tail phase while remaining blind to the maximal-rate sparse
phase in (5.2).

## 8. Verdict

Counting good signings produces two rigorous convergence reductions:

1. convergence of the refined microcanonical entropy
   \(\Sigma_n(c)\) at scale \(n^{3/2}\log n\);
2. convergence of the signing-space pressure \(\Phi_n(\beta)\) for an
   unbounded set of fixed inverse temperatures.

Either theorem would force the desired limit, including in the
presence of isolated algebraic constructions.

The ordinary speed-\(n^2\) lower-tail rate is insufficient on its own:
one switching orbit sits exactly at the maximal finite rate
\((\log2)/2\).  Exact restriction and random-extension inequalities
give local continuity but incur the leading inflation
\(\sqrt{\log2(1-\alpha^2)}\) under proportional growth.  Fixed-threshold
Shearer monotonicity follows the wrong diagonal, and the explicit
oscillating abstract profile (6.10)--(6.11) proves that it cannot be
repaired using only the scalar count and Hamming-thickening axioms.

The exact row-cavity state is the Hamming-cap profile \(V_B\); after
gauge fixing it becomes the switching orbit of the signed
energy-versus-magnetization profile.  Edgewise concentration misses its
required speed by a factor \(n\), and even the full unlabeled gap
histogram is insufficient by (4.17)--(4.18).  The remaining nonlocal
step is therefore a speed-\(n\) overlap/cap-profile theorem under the
signing Gibbs law, or an equally strong pressure interpolation
controlling the joint high-energy-layer/cross-block profile in (7.1).
