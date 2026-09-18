# Rooted moving-projection audit for the augmented cut code

This note is a scratch research report.  It does not modify the ledger or
strategic files.

## 1. Setup

Let

\[
E=E(K_n),\qquad N=|E|=\binom n2,
\qquad \Omega=\{\pm1\}^{E},
\]

and use the normalized Hamming inner product

\[
t(u,v)=\frac1N\sum_{e\in E}u_ev_e,
\qquad \tau(u)=t(\mathbf1,u).
\]

For `n >= 3`, the augmented cut code is the subgroup

\[
C=C_n^\pm
=\{c(\sigma,x):c(\sigma,x)_{ij}=\sigma x_ix_j\}
\le \Omega,
\qquad |C|=q=2^n,
\]

where `sigma` is a sign and `x` is considered modulo global negation.  For a
signing `a in Omega`, put

\[
\mu(a)=\max_{c\in C}\tau(ac)
=\frac1N\max_{c\in C}\langle a,c\rangle.
\tag{1.1}
\]

Since `C=-C`, `mu(a)=max_c |tau(ac)|`.  In matrix notation,
`N mu(a)` is exactly the absolute quadratic cap of the signing.

The kernel in OpenAI Chapter 2, Theorem 2.1 has the following two properties.
For a parameter `lambda`, it is a radial kernel

\[
K(u,v)=\kappa(t(u,v))
\]

such that

\[
K\succeq0,\qquad \kappa(t)\ge0,
\qquad F(u,v)=(t(u,v)-\lambda)K(u,v)\succeq0.
\tag{1.2}
\]

The last assertion is the complete Gram-remainder identity (27) of that
theorem.  The conclusions below apply to every kernel satisfying (1.2), not
only to the particular Boolean harmonic construction.

## 2. Exact translation blindness

### Theorem 2.1 (labeled Gram translation blindness)

For every function `phi:[-1,1] -> R`, every `a in Omega`, and every
`c,d in C`,

\[
\phi(t(ac,ad))=\phi(t(c,d)).
\tag{2.1}
\]

Consequently the entire labeled pairwise Gram matrix of the translate `aC`
is independent of `a`.  In particular, every unrooted two-point Delsarte sum,
every spectrum or rank of a scalar two-point kernel on `aC`, and every
inequality using only those data is translation blind.

#### Proof

Coordinatewise multiplication gives

\[
t(ac,ad)=\frac1N\sum_e a_e^2c_ed_e=t(c,d).
\]

This proves the pointwise matrix identity, not merely equality after summing.
\(\square\)

This blindness is substantive.  For `a=1`, `N mu(a)=N`.  On the other hand,
a union bound gives an `a` satisfying

\[
N\mu(a)\le
\sqrt{2N(n+2)\log2}=O(n^{3/2}).
\tag{2.2}
\]

Indeed, for random independent signs `a_e` and fixed `c`, Hoeffding gives

\[
\Pr\{|\langle a,c\rangle|\ge T\}
\le2e^{-T^2/(2N)}.
\]

Union over the `2^n` augmented cuts and take
`T=sqrt(2N(n+2) log 2)`.  Thus identical translated Gram matrices coexist
with caps of orders `n^2` and `n^(3/2)`.

## 3. The exact rooted Fourier decomposition

Rooting at `1` does escape Theorem 2.1 at the level of raw data, because

\[
K(\mathbf1,ac)=\kappa(\tau(ac))
\tag{3.1}
\]

depends on the translate.  The information content can nevertheless be
described exactly.

For an edge set `S subseteq E`, write `u_S=prod_(e in S)u_e` and let
`partial S` be its set of odd-degree graph vertices.  A radial positive
semidefinite kernel `R(u,v)=r(t(u,v))` has a Walsh expansion

\[
R(u,v)=\sum_{S\subseteq E}\rho_{|S|}u_Sv_S,
\qquad \rho_j\ge0.
\tag{3.2}
\]

The restriction of the ambient character `u -> u_S` to `C` is

\[
c(\sigma,x)_S
=\sigma^{|S|}\prod_i x_i^{\deg_S(i)}.
\tag{3.3}
\]

It is therefore indexed by the syndrome

\[
\pi(S)=(|S|\bmod2,\partial S).
\tag{3.4}
\]

Here `partial S` always has even cardinality.  The kernel of `pi` is exactly

\[
D=C^\perp
=\{S:\partial S=\varnothing,\ |S|\text{ even}\},
\tag{3.5}
\]

the even-cardinality Eulerian edge sets.

For a character `eta in C-hat`, define

\[
A_\eta=\sum_{\pi(S)=\eta}\rho_{|S|},
\qquad
B_\eta(a)=\sum_{\pi(S)=\eta}\rho_{|S|}a_S.
\tag{3.6}
\]

### Theorem 3.1 (rooted Schur complement is fiberwise Cauchy)

Let `G=(R(ac,ad))_(c,d in C)` and
`h_a=(R(1,ac))_(c in C)`.  The characters of `C` diagonalize `G`; the
eigenvalue on character `eta` is `q A_eta`, and the normalized Fourier
coefficient of `h_a` is `B_eta(a)`.  Hence the Schur complement of `G` in the
Gram matrix on `{1} union aC` is

\[
h_a^{\mathsf T}G^\dagger h_a
=\sum_{A_\eta>0}\frac{B_\eta(a)^2}{A_\eta}
\le r(1)=\sum_\eta A_\eta.
\tag{3.7}
\]

Moreover, (3.7) supplies no coupling between syndromes.  If

\[
P_\eta=\sum_{\substack{\pi(S)=\eta\\a_S=1}}\rho_{|S|},
\qquad
Q_\eta=\sum_{\substack{\pi(S)=\eta\\a_S=-1}}\rho_{|S|},
\]

then

\[
A_\eta-\frac{B_\eta(a)^2}{A_\eta}
=\frac{4P_\eta Q_\eta}{P_\eta+Q_\eta}\ge0.
\tag{3.8}
\]

Thus the full scalar Schur inequality is exactly the sum of independent
weighted two-sign Cauchy inequalities inside the restriction fibers.

#### Proof

Substitute (3.2) and use character orthogonality on `C`.  Convolution by the
restriction of `R` has eigenvalue `q A_eta`.  Likewise,

\[
\mathbb E_{c\in C}R(\mathbf1,ac)\eta(c)=B_\eta(a).
\]

This gives (3.7).  Since `A_eta=P_eta+Q_eta` and
`B_eta=P_eta-Q_eta`, (3.8) follows by direct expansion.  \(\square\)

The trivial syndrome is special:

\[
B_0(a)
=\mathbb E_{c\in C}r(\tau(ac))
=\sum_{S\in D}\rho_{|S|}a_S.
\tag{3.9}
\]

It is exactly a signed even-Eulerian polynomial, and also exactly a linear
functional of the coset energy/distance histogram.  If

\[
A_j(a)=|\{c\in C:d_H(a,c)=j\}|,
\]

then

\[
B_0(a)=\frac1q\sum_{j=0}^N
A_j(a)r(1-2j/N).
\tag{3.10}
\]

Thus every constant-vector rooted argument is precisely in the ledger's
Section 3.15 signed-Eulerian/coset-histogram class.

The nontrivial `B_eta` are not determined by that histogram.  Writing
`eta=(p,S)` makes them the boundary sectors

\[
B_{p,S}(a)=
\sum_{\substack{\partial F=S\\|F|\equiv p\ (2)}}
\rho_{|F|}a_F.
\tag{3.11}
\]

There are `2^n` such sectors.  Together they are exactly the Fourier transform
of the vector `c -> r(tau(ac))`; Fourier inversion shows that they are not an
exact compression of that transformed energy landscape.  Because `A_eta`
depends only on `(p,|S|)`, the single scalar in (3.7) can be evaluated from the
`O(n)` orbit-power sums

\[
\sum_{|S|=j}B_{p,S}(a)^2.
\tag{3.12}
\]

This is a polynomial-size summary for one already-tautological inequality,
not a closed state for pointwise cap control.  Retaining phases or performing
further rooted transitions reopens the full boundary-sector hierarchy in
(3.11).  A genuinely matrix-valued Terwilliger certificate could add new
constraints; Theorem 3.1 rules out only the scalar two-point Schur use.

## 4. A non-tautological rooted annular criterion

There is one clean way in which the root and the pointwise nonnegativity of
the moving-projection overlap can be used without reconstructing the full
histogram.

Apply (3.2) to the Gram remainder

\[
F(u,v)=(t(u,v)-\lambda)K(u,v)
=\sum_{S\subseteq E}\beta_{|S|}u_Sv_S,
\qquad \beta_j\ge0.
\tag{4.1}
\]

Define the translate-independent leakage

\[
J=\sum_{c\in C}F(\mathbf1,c)
=q\sum_{S\in D}\beta_{|S|}\ge0
\tag{4.2}
\]

and the rooted overlap mass

\[
T_a=\sum_{c\in C}K(\mathbf1,ac)
=\sum_{c\in C}\kappa(\tau(ac)).
\tag{4.3}
\]

### Theorem 4.1 (rooted annular inequality)

For every `a in Omega`,

\[
(\lambda-\mu(a))_+T_a\le J.
\tag{4.4}
\]

For `0 <= r < lambda`, let

\[
m_\kappa(r)=
\min\{\kappa(t):t\in\{1-2j/N:0\le j\le N\},\ |t|\le r\}.
\tag{4.5}
\]

Then the entirely kernel-level condition

\[
(\lambda-r)q,m_\kappa(r)>J
\tag{4.6}
\]

implies

\[
\mu(a)>r\qquad\text{for every }a\in\Omega.
\tag{4.7}
\]

#### Proof

Fourier positivity and (3.5) give, for every translate,

\[
\left|\sum_{c\in C}F(\mathbf1,ac)\right|
=q\left|\sum_{S\in D}\beta_{|S|}a_S\right|
\le q\sum_{S\in D}\beta_{|S|}=J.
\tag{4.8}
\]

If `mu(a)<lambda`, pointwise nonnegativity of `kappa` gives

\[
\begin{aligned}
(\lambda-\mu(a))T_a
&\le\sum_{c\in C}(\lambda-\tau(ac))\kappa(\tau(ac))\\
&=-\sum_{c\in C}F(\mathbf1,ac)
\le J.
\end{aligned}
\]

This proves (4.4).  If `mu(a)<=r`, antipodality of `C` gives
`|tau(ac)|<=r` for every `c`, so `T_a>=q m_kappa(r)`.  Combining this
with (4.4) contradicts (4.6).  \(\square\)

The fixed leakage has an explicit binomial formula.  If `X_i` and `sigma`
are independent uniform signs and `S_n=sum_i X_i`, then

\[
\frac Jq
=\mathbb E_{\sigma,X}
\left[
F\!\left(\sigma\frac{S_n^2-n}{n(n-1)}\right)
\right].
\tag{4.9}
\]

Therefore (4.6) does **not** hide the parent maximization: it asks only for a
uniform central-annulus lower bound on a known scalar kernel and a binomial
average over the known code.

For the desired quadratic-signing scale, the precise finite statement is to
find moving-projection parameters with

\[
r_n=\frac{\gamma n^{3/2}}{N}=\Theta(n^{-1/2}),
\qquad \lambda_n>r_n,
\tag{4.10}
\]

and prove

\[
\boxed{
\frac{J_n}{q\,m_{\kappa_n}(r_n)}<\lambda_n-r_n.
}
\tag{4.11}
\]

This is the minimal noncircular rooted lemma exposed by the scalar moving
projection.

If the annular minimum vanishes or is too small, replacing it by a sharper
signing-dependent estimate re-enters the old state exactly:

\[
\frac{T_a}{q}
=\sum_{S\in D}\alpha_{|S|}a_S
=\frac1q\sum_jA_j(a)\kappa(1-2j/N),
\tag{4.12}
\]

where `alpha` are the Walsh coefficients of `K`.  Thus exact evaluation of
`T_a` is a single functional of the complete coset energy histogram.  A new
uniform theorem lower-bounding (4.12) could still be progress, but computing
it by enumerating the histogram or defining its infimum over
`{a:mu(a)<=r}` merely reinstates the full coset optimization.  The clean test
that avoids that trap is (4.11).

## 5. Transversal-square kernels: an exact reduction and a spectral no-go

There is a tempting way to make `T_a` exactly constant.  It produces a real
dual theorem, but a hypercube spectral bound proves that no scalar instance of
this construction can reach the `n^(3/2)` scale.

Identify edge sets with the additive group

\[
\mathcal H=\mathbb F_2^E
\]

and let `A_H` be normalized adjacency on the `N`-cube:

\[
(A_{\mathcal H}v)_R=\frac1N\sum_{e\in E}v_{R\mathbin\triangle\{e\}}.
\tag{5.1}
\]

Let `S subset H` be a partial transversal of `D`, meaning that the quotient
map `pi:H -> H/D` is injective on `S`.  Suppose `v>=0` is supported on `S` and

\[
A_{\mathcal H}v\ge\lambda v
\tag{5.2}
\]

coordinatewise.  A Perron vector of the cube subgraph induced by `S` is the
standard source of (5.2).  Define

\[
g(z)=\sum_{R\in\mathcal H}v_Rz_R,
\qquad
K_v(u,w)=g(uw)^2.
\tag{5.3}
\]

This kernel is translation invariant but need not be radial.

Let `A_bar` be normalized adjacency of the quotient Cayley graph
`H/D`, with the images of the `N` single-edge generators.  The partial
transversal identifies `v` with a function `v_bar` on `H/D`, extended by zero,
and put

\[
\rho(v)=
\frac{\langle\bar v,\bar A\bar v\rangle}{\|v\|_2^2}.
\tag{5.4}
\]

### Theorem 5.1 (exact transversal certificate)

The kernels

\[
K_v(u,w),
\qquad
(t(u,w)-\lambda)K_v(u,w)
\]

are positive semidefinite, and `K_v` is pointwise nonnegative.  Moreover,

\[
\sum_{c\in C}K_v(\mathbf1,ac)=q\|v\|_2^2
\qquad(a\in\Omega),
\tag{5.5}
\]

the leakage is exactly

\[
\frac1q\sum_{c\in C}
(\tau(c)-\lambda)K_v(\mathbf1,c)
=(\rho(v)-\lambda)\|v\|_2^2,
\tag{5.6}
\]

and consequently

\[
\boxed{\qquad
\mu(a)\ge2\lambda-\rho(v)
\quad\text{for every }a\in\Omega.
\qquad}
\tag{5.7}
\]

#### Proof

The Walsh coefficients of `g^2` are the nonnegative autocorrelation

\[
Q_T=(v*v)_T=\sum_Rv_Rv_{R\mathbin\triangle T}\ge0.
\tag{5.8}
\]

Multiplication by `tau` acts on Walsh coefficients by `A_H`, and convolution
commutes with this adjacency operator.  Hence the coefficients of
`(tau-lambda)g^2` are

\[
A_{\mathcal H}Q-\lambda Q
=(A_{\mathcal H}v-\lambda v)*v\ge0.
\tag{5.9}
\]

This proves both positive-semidefiniteness statements.  If `T in D` is
nonzero, transversality says that `S` and `S triangle T` are disjoint, so

\[
Q_T=0\quad(T\in D\setminus\{0\}),
\qquad Q_0=\|v\|_2^2.
\tag{5.10}
\]

The dual expansion now gives (5.5).

Put `w=A_H v-lambda v`.  Summing the dual coefficients in (5.9) gives

\[
\sum_{T\in D}(w*v)_T
=\sum_{R\in\mathcal H}w_R\bar v(\pi(R)).
\tag{5.11}
\]

The sum of `A_H v` over one quotient fiber is `A_bar v_bar` at that quotient
vertex.  Therefore (5.11) equals

\[
\langle\bar v,(\bar A-\lambda)\bar v\rangle
=(\rho(v)-\lambda)\|v\|_2^2,
\]

which proves (5.6).  The proof of Theorem 4.1 used only translation
invariance, Fourier positivity, and pointwise nonnegativity, not radiality.
Applying (4.4), (5.5), and (5.6) proves (5.7); if `mu(a)>=lambda`, the same
conclusion follows from `rho(v)>=lambda`.  \(\square\)

For a Perron vector of the lifted induced subgraph on `S`, `lambda` is its
normalized spectral radius.  Then `rho(v)-lambda` is exactly the weighted
mass of quotient edges between selected quotient vertices whose chosen
representatives are **not** adjacent in the ambient cube.  Thus (5.7) is
aligned spectral mass minus unaligned leakage.  This is a genuine dual
certificate independent of the signing `a`; it does not reconstruct the
parent optimization.

It nevertheless has a decisive size obstruction.

### Theorem 5.2 (scalar-transversal spectral ceiling)

For every sequence of partial-transversal vectors satisfying (5.2),

\[
\lambda=O\!\left(\frac1{\sqrt{n\log n}}\right)
=o(n^{-1/2}).
\tag{5.12}
\]

Consequently the right side of (5.7) is at most
`o(n^(-1/2))`; scalar transversal-square kernels cannot prove a positive
`n^(3/2)`-scale cap constant, even if their leakage vanishes.

#### Proof

The support of `v` has size at most the number of `D`-cosets:

\[
|\operatorname{supp}v|\le q=2^n.
\tag{5.13}
\]

Set

\[
i=\left\lceil\frac{2n}{\log n}\right\rceil.
\]

Then `i=o(N)` and, for all sufficiently large `n`,

\[
2^n=O\!\left(\sum_{j=0}^i\binom Nj\right).
\tag{5.14}
\]

Theorem 4 of Bollobas--Lee--Letzter, *Eigenvalues of subgraphs of the
cube*, applied in ambient dimension `N`, says that any cube subgraph with
this many vertices has unnormalized spectral radius at most

\[
(1+o(1))\lambda_1(H_i^N)
=(2+o(1))\sqrt{i(N-i)}.
\tag{5.15}
\]

Condition (5.2), followed by the Rayleigh inequality on the support of `v`,
therefore gives

\[
\lambda
\le\frac{\lambda_1(Q_N[\operatorname{supp}v])}{N}
\le(2+o(1))\sqrt{\frac iN}
=O\!\left(\frac1{\sqrt{n\log n}}\right).
\]

Finally, `rho(v)>=lambda`, so `2lambda-rho(v)<=lambda`.  \(\square\)

More quantitatively, the same theorem shows that normalized spectral radius
`Omega(n^(-1/2))` needs scalar Fourier support
`exp(Omega(n log n))`.  A partial transversal has only `exp(O(n))` points.
Thus a correct-scale use of this architecture must abandon scalar
transversality, for example through genuinely higher rank/multiple
representatives per syndrome with cancellation retained jointly.

The elementary star complement illustrates both outcomes.  A full linear
section gives two `(n-1)`-cubes and normalized lifted eigenvalue `2/n`, but
the full-section uniform vector has quotient Rayleigh quotient `rho=1`, so
(5.7) is negative.  Restricting `v` to one cube component removes the
unaligned internal edges: `rho=lambda=2/n`, `J=0`, and (5.7) gives only
`mu(a)>=2/n`, i.e. an order-`n` rather than order-`n^(3/2)` energy bound.

## 6. Literature boundary

The focused primary-source check found the following.

1. Sole and Zaslavsky, *A Coding Approach to Signed Graphs*, SIAM J.
   Discrete Math. 7 (1994), 544--553, proves the exact correspondence between
   switching classes of signed graphs and cosets of the cocycle code, and
   identifies covering radius with maximum imbalance/frustration.  It does not
   give a translate-uniform lower bound for a moving-projection overlap such as
   `T_a`.
2. Classical Delsarte external-distance theorems bound a code's covering
   radius using the number of nonzero dual weights.  The dual here is the rich
   family of even Eulerian graphs, so that theorem neither resolves the
   `N/2-Theta(n^(3/2))` radius deficit nor gives (4.11) or (4.12).
3. Standard Terwilliger/completely-regular-code results control rooted distance
   partitions when an equitable or completely regular structure is already
   present.  No applicable theorem was found making all translates of the
   complete-graph cocycle code equitable or bounding (4.12).  For a generic
   external root `a`, its `S_n` stabilizer is trivial, and the exact rooted
   Fourier data are the `2^n` sectors (3.11).
4. Bollobas, Lee, and Letzter, *Eigenvalues of subgraphs of the cube*,
   European J. Combin. 70 (2018), 125--148, gives the sharp asymptotic
   small-support spectral comparison used in Theorem 5.2.  It closes the
   scalar partial-transversal implementation at the required scale.

Thus there is no located theorem that supplies the missing annular estimate.
The new OpenAI identity contributes a promising family of kernels, but the
cut-code-specific analytic input (4.11) remains genuinely new.

## 7. Verdict

* **Unrooted scalar two-point use is rigorously dead:** the full labeled Gram
  matrix is translation blind.
* **The raw root does see the signing:** its values are a pointwise transform
  of the cut-energy landscape.
* **The full scalar Schur complement does not help:** after restriction to the
  code it is exactly fiberwise Cauchy, equation (3.8).
* **Constant-mode rooted estimates are exactly the signed Eulerian/coset
  histogram route:** equations (3.9)--(3.10).
* **The uniform annular ratio (4.11) is a genuinely compressed target:** it
  does not require signing-dependent state and has the correct normalization.
  However, Theorem 5.2 proves that the natural scalar partial-transversal
  square implementation cannot reach that normalization.
* If (4.11) fails because `kappa_n` has zeros or an exponentially inadequate
  floor throughout the central annulus, the scalar moving-projection route
  falls back to (4.12), and hence to the already identified Eulerian
  obligation.  A matrix-valued rooted/Terwilliger construction remains outside
  this no-go theorem.
* **Exact next lemma:** construct a higher-rank moving fiber (or an equivalent
  joint multi-representative section) whose rooted overlap mass is uniformly
  controlled and whose aligned-minus-unaligned spectral value is
  `Omega(n^(-1/2))`.  It must exploit cancellation before scalarization;
  scalar support of size at most `2^n` is now rigorously insufficient.
