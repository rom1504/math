# Eulerian free energy at the changing temperature

## 1. The object and its exact channel interpretation

Let

\[
N=\binom n2,\qquad
H_A(x)=\sum_{i<j}a_{ij}x_ix_j,\qquad
\rho=\tanh t.
\]

The normalized two-sided partition function is

\[
W_A(\rho)
=
(\cosh t)^{-N}\,\mathbb E_x\cosh(tH_A(x))
=
\mathbb E_{\sigma,x}
\prod_{i<j}(1+\rho\,\sigma a_{ij}x_ix_j),
\tag{1.1}
\]

where \(\sigma\) and the \(x_i\)'s are independent uniform signs.  Its
minimum pressure is

\[
\Gamma_n(\rho)=\min_A\log W_A(\rho).
\]

Expanding (1.1) gives the exact even-Eulerian polynomial

\[
W_A(\rho)
=
\sum_{\substack{F\subseteq E(K_n)\\
                 \partial F=\varnothing,\ |F|\ {\rm even}}}
       \rho^{|F|}\prod_{e\in F}a_e.
\tag{1.2}
\]

There is also an exact information-theoretic interpretation.  Let
\(\mathcal C_n\) be the augmented cut code

\[
\mathcal C_n=\{(\sigma x_ix_j)_{i<j}:\sigma,x_i\in\{\pm1\}\}.
\]

Choose \(C\) uniformly from \(\mathcal C_n\), choose independent signs
\(\eta_e\) with \(\mathbb E\eta_e=\rho\), and observe

\[
Y=C\odot\eta.
\]

Then

\[
\boxed{\quad
2^N\Pr_\rho(Y=A)=W_A(\rho).
\quad}
\tag{1.3}
\]

Consequently

\[
\boxed{\quad
-\Gamma_n(\rho)
=D_\infty\!\left({\rm Uniform}(\{\pm1\}^N)
                 \,\middle\|\,\mathcal L(Y)\right).
\quad}
\tag{1.4}
\]

Thus the desired changing-temperature limit is the normalized
order-\(\infty\) resolvability exponent of the augmented cut code through
a BSC whose bias is \(\rho=\beta/\sqrt n+O(n^{-3/2})\).  This is stronger
than the usual average or finite-Rényi soft-covering problem: it asks for
the least likely output word.

Equivalently, if

\[
A_j(A)=|\{C\in\mathcal C_n:d(A,C)=j\}|,
\]

then the exact MacWilliams/coset form is

\[
W_A(\rho)
=
\frac1{|\mathcal C_n|}
\sum_{j=0}^N
A_j(A)(1+\rho)^{N-j}(1-\rho)^j.
\tag{1.5}
\]

This makes the zero-temperature connection transparent: the order of
vanishing as \(\rho\uparrow1\) is the distance from \(A\) to
\(\mathcal C_n\), hence the covering-radius problem reappears at the
endpoint.

## 2. Exact edge deletion and local anti-alignment

Let \(G\) be a partial signed graph and let \(e=uv\) be absent.  Put

\[
d_e(G)
=
\mathbb E_{\mu_G}[\sigma x_ux_v],
\]

where \(\mu_G\) is the probability measure on \((\sigma,x)\) whose
density relative to the uniform measure is

\[
\frac1{W_G(\rho)}
\prod_{f\in E(G)}(1+\rho\,\sigma a_fx_f).
\]

Adding \(e\) with sign \(s\) gives the exact deletion identity

\[
\boxed{\quad
W_{G+se}(\rho)=W_G(\rho)\bigl(1+s\rho d_e(G)\bigr).
\quad}
\tag{2.1}
\]

If \(A\) is a global minimizer of \(W_A(\rho)\), compare its chosen edge
sign with the opposite sign while all other edges are fixed.  Equation
(2.1) gives

\[
\boxed{\quad a_ed_e(A-e)\le0\qquad(e\in E(K_n)).\quad}
\tag{2.2}
\]

This is the exact two-sided version of the proposed leave-one-edge
anti-alignment condition.

Let

\[
r_e=\mathbb E_{\mu_A}[\sigma a_ex_ux_v],
\qquad
u_e=a_ed_e(A-e).
\]

Reweighting by the last edge gives

\[
r_e=\frac{\rho+u_e}{1+\rho u_e},
\qquad
u_e=\frac{r_e-\rho}{1-\rho r_e}.
\tag{2.3}
\]

Thus (2.2) is equivalent to \(r_e\le\rho\).  Differentiation of (1.1)
then yields the exact drift

\[
\frac{d}{d\rho}\log W_A(\rho)
=
\frac{\sum_er_e-N\rho}{1-\rho^2}
\le0.
\tag{2.4}
\]

There is a more useful flip-ratio form.  Put

\[
z_e=-\rho u_e\in[0,\rho].
\]

Since

\[
\frac{W_{A^{(e)}}}{W_A}
=
\frac{1+z_e}{1-z_e},
\]

where \(A^{(e)}\) is obtained by flipping edge \(e\), the hypercube
Fourier identity

\[
\rho W_A'(\rho)
=
\frac12\sum_e\bigl(W_A-W_{A^{(e)}}\bigr)
\]

becomes

\[
\boxed{\quad
-\rho\frac{d}{d\rho}\log W_A(\rho)
=
\sum_e\frac{z_e}{1-z_e}.
\quad}
\tag{2.5}
\]

This is exact at every edgewise local minimum, in particular at every
global minimum.

## 3. The two-replica susceptibility does not close the drift

For a partial graph \(G\), take two independent replicas
\((\sigma,x)\) and \((\tau,y)\) from \(\mu_G\).  Then

\[
\boxed{\quad
\sum_{i<j}d_{ij}(G)^2
=
\frac12\,
\mathbb E\!\left[
\sigma\tau\bigl((x\cdot y)^2-n\bigr)
\right].
\quad}
\tag{3.1}
\]

Thus the proposed sum of squared leave-one-edge correlations is exactly
a signed two-replica susceptibility.

At a local minimum let \(S_2=\sum_eu_e^2\).  Equations (2.5) and
\(0\le z_e\le\rho\) imply only

\[
\boxed{\quad
\rho S_2
\le
-\rho(\log W_A)'(\rho)
\le
\frac{\rho\sqrt{NS_2}}{1-\rho}.
\quad}
\tag{3.2}
\]

The lower bound uses
\(z/(1-z)\ge z\ge z^2/\rho\), and the upper bound uses Cauchy--Schwarz.
Both orders are sharp for vectors with concentrated or delocalized
coordinates.

At the changing temperature \(\rho\asymp n^{-1/2}\), the natural
susceptibility size is \(S_2=\Theta(n)\).  Its lower control in (3.2) is
only \(O(\sqrt n)\), whereas a normalized-pressure theorem needs a
drift of order \(n\).  Recovering that missing factor requires an
\(L^1/L^2\) delocalization theorem for the individual correlations.
Neither positivity nor edgewise anti-alignment supplies it.

There is a finite exact illustration.  The order-six signing

\[
A_{\rm trap}=
\begin{pmatrix}
0&1&1&1&1&1\\
1&0&1&-1&-1&1\\
1&1&0&-1&-1&1\\
1&-1&-1&0&1&1\\
1&-1&-1&1&0&1\\
1&1&1&1&1&0
\end{pmatrix}
\]

has

\[
W_{A_{\rm trap}}(\rho)=(1-\rho^4)^3.
\tag{3.3}
\]

It is an edgewise local minimum for every \(0<\rho<1\).  Twelve of its
fifteen edge flips leave \(W\) exactly unchanged.  For each of the
remaining three edges,

\[
u_e=-\frac{4\rho^3}{1+3\rho^4},
\]

and hence

\[
S_2=\frac{48\rho^6}{(1+3\rho^4)^2}.
\tag{3.4}
\]

Nevertheless the conference signing of the same order has

\[
W_C(\rho)
=
1-15\rho^4+40\rho^6-45\rho^8+24\rho^{10}-5\rho^{12},
\]

and

\[
W_{A_{\rm trap}}(\rho)-W_C(\rho)
=
4\rho^4(1-\rho^2)^3(3-\rho^2)>0.
\tag{3.5}
\]

Thus exact edge anti-alignment plus its two-replica susceptibility can
be trapped in a highly concentrated, globally nonoptimal basin.

## 4. Exact vertex recursion and the unavoidable overlap hierarchy

Write an order-\((m+1)\) signing as an old signing \(B\) and a new row
\(b\in\{\pm1\}^m\).  For an internal edge set \(F\), let
\(\partial F\) be its odd-degree vertex set.  Decomposing an Eulerian
subgraph according to its star edges gives

\[
\boxed{\quad
W_{B,b}(\rho)
=
\sum_{\substack{F\subseteq E(K_m)\\|F|\ {\rm even}}}
\rho^{|F|+|\partial F|}
\left(\prod_{e\in F}B_e\right)
\left(\prod_{i\in\partial F}b_i\right).
\quad}
\tag{4.1}
\]

Equivalently, for every even \(S\subseteq[m]\), define the boundary
sector

\[
Z_{B,S}(\rho)
=
\sum_{\substack{\partial F=S\\|F|\ {\rm even}}}
\rho^{|F|}\prod_{e\in F}B_e.
\]

Then

\[
W_{B,b}(\rho)
=
\sum_{S\ {\rm even}}\rho^{|S|}b_SZ_{B,S}(\rho).
\tag{4.2}
\]

The mean over \(b\) keeps only \(S=\varnothing\), giving the known
monotonicity.  But exact minimization over \(b\) requires all the other
sectors.  For \(m\ge3\), every even \(S\) has a nonempty formal sector:
start with any \(S\)-join and, if its edge parity is wrong, take its
symmetric difference with a triangle.  Hence one vertex deletion
already opens \(2^{m-1}\) boundary sectors.

The same obstruction appears probabilistically.  Put

\[
R_b=\frac{W_{B,b}}{W_B}.
\]

For every integer \(k\ge1\), an exact average over the new row gives

\[
\boxed{
\mathbb E_bR_b^k
=
\mathbb E_{\mu_B^{\otimes k}}
\prod_{i=1}^m
\left[
\sum_{\substack{J\subseteq[k]\\|J|\ {\rm even}}}
\rho^{|J|}
\prod_{r\in J}\sigma_rx_i^{(r)}
\right].
}
\tag{4.3}
\]

In particular,

\[
\boxed{
\mathbb E_bR_b^2
=
\mathbb E_{\mu_B^{\otimes2}}
\prod_{i=1}^m
\left(1+\rho^2\sigma_1\sigma_2
x_i^{(1)}x_i^{(2)}\right).
}
\tag{4.4}
\]

At \(\rho=\beta/\sqrt n\), (4.4) sees the ordinary two-replica overlap
at leading order.  The \(k\)-th moment sees all even multi-overlaps of
up to \(k\) replicas.  But \(\min_bR_b\) is an \(L_{-\infty}\) endpoint
over \(2^{m-1}\) effective rows.  Consequently no fixed replica depth,
including the susceptibility in (3.1), closes the vertex recursion.
Rare resonant rows require \(k\) growing with \(m\).

## 5. The first nontrivial coefficient is exactly solvable

The weight-four Eulerian subgraphs are the simple \(4\)-cycles.  Hence

\[
T_4(A)=[\rho^4]W_A(\rho)
=\sum_{C_4}\prod_{e\in C_4}a_e.
\]

Direct expansion of \(\operatorname{tr}A^4=\|A^2\|_F^2\) gives

\[
\boxed{
T_4(A)
=
\frac{\|A^2\|_F^2-n(n-1)(2n-3)}8.
}
\tag{5.1}
\]

Since

\[
\|A^2\|_F^2
=\sum_i\lambda_i^4
\ge\frac1n\left(\sum_i\lambda_i^2\right)^2
=n(n-1)^2,
\]

one obtains

\[
\boxed{
T_4(A)\ge-\frac{n(n-1)(n-2)}8,
}
\tag{5.2}
\]

with equality exactly when \(A^2=(n-1)I\), i.e. when \(A\) is a
symmetric conference matrix.

This proves that conference matrices minimize the first
signing-dependent Eulerian coefficient.  It does not order the full
pressure polynomial.  The failure already occurs at order six.  The
two valid signings

\[
A_1=
\begin{pmatrix}
0&1&1&1&1&1\\
1&0&1&1&-1&1\\
1&1&0&-1&1&1\\
1&1&-1&0&1&1\\
1&-1&1&1&0&1\\
1&1&1&1&1&0
\end{pmatrix},
\quad
A_2=
\begin{pmatrix}
0&1&1&1&1&1\\
1&0&-1&-1&1&1\\
1&-1&0&-1&1&1\\
1&-1&-1&0&1&1\\
1&1&1&1&0&1\\
1&1&1&1&1&0
\end{pmatrix}
\]

have

\[
W_{A_1}=1+5\rho^4-13\rho^8+7\rho^{12},
\]

\[
W_{A_2}
=1+9\rho^4-40\rho^6+51\rho^8-24\rho^{10}+3\rho^{12},
\]

and

\[
\boxed{
W_{A_1}-W_{A_2}
=
4\rho^4(1-\rho^2)^2
\bigl(\rho^4+8\rho^2-1\bigr).
}
\tag{5.3}
\]

Their ordering reverses at

\[
\rho=\sqrt{\sqrt{17}-4}=0.350864\ldots.
\]

Thus neither the first coefficient nor the zero-temperature distance
orders the finite-temperature likelihood.  In particular, an argument
that combines (5.2) with (2.2) must prove a genuinely uniform
no-cancellation theorem; it cannot infer one from coefficient signs.

## 6. Deep-hole cap constraints and why they do not supply entropy

Let \(A\) be a deep hole of \(\mathcal C_n\), with covering radius \(R\),
and use \(A_j(A)\) from (1.5).  Every word on the Hamming sphere of
radius \(k\) around \(A\) lies in some radius-\(R\) code ball.  Averaging
this covering over the sphere gives the exact inequality

\[
\boxed{
\sum_j A_j(A)\,p_{N,R,k}(j)\ge1,
}
\tag{6.1}
\]

where

\[
p_{N,R,k}(j)
=
\frac1{\binom Nk}
\sum_r
\binom jr\binom{N-j}{k-r}
\mathbf1_{\{j+k-2r\le R\}}.
\tag{6.2}
\]

For \(k=1\), this becomes

\[
\boxed{
R\,A_R(A)+(R+1)A_{R+1}(A)\ge N.
}
\tag{6.3}
\]

This is the radial form of the edge-flip deep-hole certificate.  More
generally, if \(k=O(n)\), a codeword at distance \(R+O(\sqrt n)\) covers
a constant-order fraction of the \(k\)-sphere, so (6.1) forces only
constant-order near-ground multiplicity.  Exponential multiplicity is
forced only when \(k=\Theta(N)\), but then (6.1) includes energy gaps of
order \(N\), not a near-ground layer.

Pairwise inclusion-exclusion does not close radially: intersections of
two caps inside a sphere centered at \(A\) depend on the triple distance
distribution of \((A,C_1,C_2)\), not only on the known weight
distribution of the cut code or on \(A_j(A)\).  Continuing introduces
the full Terwilliger/overlap hierarchy, the same hierarchy as (4.3).

## 7. Verdict and the surviving exact target

The Eulerian route produces useful exact structure, but it does not
currently prove convergence.

What is now ruled out is the direct finite-hierarchy program:

1. edgewise anti-alignment is only a local condition and has exact
   nonglobal traps;
2. its two-replica susceptibility loses a factor \(\sqrt n\) at the
   required scale unless one proves correlation delocalization;
3. one vertex deletion generates all \(2^{n-2}\) nontrivial boundary
   sectors;
4. finite replica depth cannot control the \(L_{-\infty}\) output
   likelihood;
5. Eulerian coefficient order itself changes with temperature.

The remaining coefficient-level target is precise:

> Prove that
> \[
> \frac1n
> D_\infty\!\left(
> {\rm Uniform}
> \,\middle\|\,
> {\rm BSC}_{\beta/\sqrt n}(\mathcal C_n)
> \right)
> \]
> converges for every fixed \(\beta\), by controlling the full
> growing-replica/boundary-sector hierarchy uniformly in \(n\).

Equivalently, one needs either a new \(L^1/L^2\) delocalization theorem
for the leave-one-edge correlations of global minimizers, or a
large-deviation principle for the complete overlap array strong enough
to reach the least-likely output rather than a finite Rényi moment.
Neither follows from deletion--contraction, MacWilliams identities,
the \(T_4\) extremum, or two-replica susceptibility alone.
