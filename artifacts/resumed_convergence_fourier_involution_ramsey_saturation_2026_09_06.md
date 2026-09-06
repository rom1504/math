# Ramsey saturation of centered Fourier involutions

Date: 2026-09-06. A generalization of the all-order Paley squarefree
argument, using finite Ramsey theorems instead of character-sum estimates.
The theorem is about Fourier-diagonal translation-invariant operators,
not arbitrary symmetric signings or arbitrary flat orthogonal matrices.

## 1. Exact finite-field theorem

Let `q` range over odd prime powers. On real functions on `F_q`, with
normalized counting measure, let `U_q` have additive Fourier multipliers

\[
m_q(-a)=m_q(a)\in\{\pm1\}\quad(a\ne0),
\qquad m_q(0)\in[-1,1].
\]

Uniformly over **all** choices of these even multipliers,

\[
\max_{f:\mathbb F_q\to\{\pm1\}}
|\langle f,U_qf\rangle|\longrightarrow1.             \tag{1}
\]

More strongly there are witnesses with `|E f|->0`. The interesting
centered case is `m_q(0)=0`, when `U_q^2=I-E`; if `m_q(0)=+1` or `-1`,
the constant function trivially attains the unqualified norm in (1),
but the asymptotically mean-zero assertion is still substantive.

Only one of the two signs is guaranteed. The earlier Paley theorem
separately saturates **both** endpoints using prescribed character
patterns, which is stronger in that specific family.

## 2. Analytic reduction to a finite frequency configuration

For every `epsilon>0`, there are fixed integers `r,D,K`, independent of
the odd prime `p`, with the following property. Suppose
`a_1,...,a_r in F_q` satisfy:

1. Every nonzero signed squarefree sum

   \[
   \sum_j\varepsilon_j a_j,\qquad
   \varepsilon_j\in\{0,\pm1\},\quad
   1\le|\operatorname{supp}\varepsilon|\le D,
   \quad|\operatorname{supp}\varepsilon|\text{ odd},
   \]

   has one common multiplier sign `s`.

2. Whenever integers `ell_j` satisfy `sum |ell_j|<=K`,

   \[
   \sum_j\ell_j a_j=0
   \quad\Longleftrightarrow\quad
   \ell_j=0\pmod p\ \text{for every }j.             \tag{2}
   \]

Then there is a Boolean `f` and a real `P` in the multiplier-`s`
eigenspace with `||f-P||_2<epsilon` and `E P=0`. Consequently

\[
s\langle f,U_q f\rangle\ge1-2\epsilon^2,
\qquad |\mathbb E f|\le\epsilon.                   \tag{3}
\]

Here is the analytic construction, including the trace-independence
caveat. For iid uniform `Y_j in F_p`, put
`X_j=sqrt(2) cos(2 pi Y_j/p)` and `S=r^(-1/2)sum X_j`. The variables are
uniformly bounded, mean zero and variance one. A finite odd Hermite
polynomial `h` approximates the Gaussian sign function. Its degree-`l`
Hermite terms are asymptotically, uniformly in `p`, approximated in `L2`
by the ordered squarefree statistics

\[
W_l=r^{-l/2}\sum_{i_1,\ldots,i_l\ \mathrm{distinct}}
X_{i_1}\cdots X_{i_l}.
\]

Indeed both `E W_l^2` and `E[He_l(S) W_l]` equal
`l!(r)_l/r^l` exactly, and the remaining Hermite square norm converges
uniformly by bounded-variable moment enumeration. Replacing all Hermite
terms by these statistics gives a polynomial `P` supported on the
squarefree odd frequency set above and approximating `sign(S)`.

Now substitute
`X_j(t)=sqrt(2)cos(2 pi Tr(a_jt)/p)`. Condition (2) makes all polynomial
moments up to degree `K` match the iid model exactly. The potentially
discontinuous sign is handled through

\[
g(u)=(\operatorname{sign}u-h(u))^2,
\]

which is continuous even at zero because `h` is odd. For fixed `r`,
approximate `g` uniformly on `[-sqrt(2r),sqrt(2r)]` by an ordinary
polynomial, and increase `K` to its degree and to `2D`. This transfers
the squared error without ever asserting full trace-coordinate
independence. The details and the order of limits are proved in
`resumed_convergence_paley_squarefree_audit_2026_09_06.md`, Sections 2 and 4.

## 3. Large characteristic: Deuber sets give dissociation for free

The required published Ramsey input is the following finite Deuber
theorem: for positive integers `r,K,c,t`, there are `M,P,C` such that
every `t`-colored `(M,P,C)`-set contains a monochromatic `(r,K,c)`-set.
In particular **the target value `c=1` is allowed**. To be insensitive to
strict versus weak coefficient-bound conventions, invoke the theorem
with target parameter `K+1`, and retain only coefficients of absolute
value at most `K`.

The resulting monochromatic family has positive integer generators
`z_1,...,z_r` and contains all integers

\[
z_i+\sum_{j<i}\lambda_jz_j,
\qquad |\lambda_j|\le K,\quad1\le i\le r,
                                                               \tag{4}
\]

and the entire set is required to lie in the positive integers.
Accordingly positivity of the expression with every `lambda_j=-K`
forces

\[
z_i>K\sum_{j<i}z_j.                                \tag{5}
\]

This is exactly the short-relation separation needed here; it is not an
additional conclusion being silently imposed on Deuber's theorem.

There is a finite `N=N(r,K)` such that every two-coloring of `[1,N]`
contains (4). For example, use the finite Deuber theorem and fix once
and for all an `(M,P,C)`-set of positive integers, constructed with
sufficiently separated generators; let `N` be its largest element.
Every target generator lies in `[1,N]`, since it itself appears in (4).

If `p>K r N`, color integer `z in [1,N]` by `m_q(z)` after embedding the
prime field in `F_q`. Obtain the generators `z_i` and set `a_i=z_i`.
Any nonempty signed squarefree sum can be negated to have largest
nonzero coefficient `+1`, and is then an element of (4), positive by
(5). Evenness of `m_q` handles that possible negation. All these sums
have a common color, not only those of odd degree.

For a nonzero integer vector `ell` with `sum |ell_j|<=K`, let `i` be its
last nonzero coordinate. Equation (5) gives

\[
\left|\sum_j\ell_jz_j\right|
\ge z_i-K\sum_{j<i}z_j>0,
\quad
\left|\sum_j\ell_jz_j\right|\le K rN<p.
\]

Thus no such relation vanishes modulo `p`, proving (2). In this regime
`p>K`, so reduction modulo `p` cannot annihilate a nonzero permitted
coefficient. This establishes the frequency configuration for every
sufficiently large characteristic, with no condition on extension degree.

## 4. Bounded characteristic: projective Ramsey and echelon bases

Fix an odd prime `p`. The finite vector-space Ramsey theorem says that,
given `r` and a finite number of colors, sufficiently high-dimensional
`F_p`-space contains an `r`-dimensional subspace whose one-dimensional
subspaces all have the same color.

Identify `F_q` additively with `F_p^e` using an ordered basis. Every line
has a unique vector with first nonzero ambient coordinate equal to one.
Color the line by the multiplier of this vector. A homogeneous
`r`-subspace `W` exists once `e>=d(p,r)`.

Choose an echelon basis `a_1,...,a_r` of `W`, with strictly increasing
pivot positions and pivot entries one. In any nonempty signed squarefree
sum `sum eps_j a_j`, its first nonzero ambient coordinate is the first
nonzero `eps_j`, hence is `+1` or `-1`. Its line normalization therefore
changes it by at most a global sign. Since the multiplier is even, all
such sums have the homogeneous line color. A full scalar-color tuple
is unnecessary.

The `a_j` are linearly independent over `F_p`, which gives (2) for every
`K` automatically. This establishes the desired configuration at each
fixed odd characteristic once its extension degree is sufficiently large.

Choose the threshold `p>K rN` from Section 3 first. There are only
finitely many remaining primes, so taking the maximum of their finitely
many bounds `p^(d(p,r))` proves a **uniform** threshold in `q`. Combining
with Section 2 proves (1).

Both imported theorem statements were read in the primary research
article [Frankl--Graham--Rödl, *Iterated Combinatorial Density Theorems*,
JCTA 54 (1990), 95--111](https://fanchung.ucsd.edu/ron/papers/90_04_iterated.pdf):
the finite vector-space statement is in Section 3, p. 99, and the finite
Deuber statement and positive-integer definition are in Section 6,
p. 105. The original sources are respectively Graham--Leeb--Rothschild
and Deuber, as attributed there. No claim that these classical Ramsey
inputs are new is intended.

## 5. Robust finite-field version with a vanishing bad set

Suppose instead that the multipliers are real and even,

\[
\|U_q\|_{\mathrm{op}}\le C,
\qquad
\frac1q\#\{a\ne0:||m_q(a)|-1|>\eta\}\longrightarrow0
\quad\text{for every fixed }\eta>0.                \tag{6}
\]

Then the same asymptotically mean-zero Boolean witnesses satisfy

\[
\liminf_q\max_f|\langle f,U_qf\rangle|\ge1.          \tag{7}
\]

To see the additional point, fix all analytic parameters and the finite
Ramsey host before sending `q` to infinity. In large characteristic use
the host `{1,...,N}`. In each of the finitely many small characteristics
use the nonzero vectors of a fixed `d(p,r)`-dimensional subspace. All
hosts have size bounded by one constant depending only on the fixed
accuracy. Multiply the whole host by a uniformly random nonzero field
element `b`. Each of its individual nonzero frequencies becomes uniform
on `F_q^*`; a union bound supplies a dilation avoiding the bad set once
its density is sufficiently small.

Color the remaining good frequencies by the sign of their multiplier
and run the relevant Ramsey construction within this dilated host.
All spectral values supporting `P` have one sign `s` and satisfy
`s m_q(a)>=1-eta`. The moment argument and short-relation condition are
preserved by the dilation. If `||f-P||_2<epsilon`, then

\[
s\langle f,U_qf\rangle
\ge (1-\eta)(1-\epsilon)^2
-C\epsilon(2+\epsilon).                           \tag{8}
\]

Send `q` to infinity first, then `eta,epsilon` to zero. This proves (7).
The random dilation is only an existence proof for the host; the
multiplier itself is completely arbitrary and deterministic.

## 6. Actual signing consequence: spectrally near-optimal additive Cayley families

Let `A_q` be any real symmetric additive-Cayley hollow signing on `F_q`:

\[
(A_q)_{t,u}=a_q(t-u),\qquad
a_q(0)=0,\quad a_q(v)=a_q(-v)\in\{\pm1\}\ (v\ne0).
\]

If

\[
\|A_q\|_{\mathrm{op}}\le(1+o(1))\sqrt q,            \tag{9}
\]

then

\[
\frac{Q(A_q)}{q^{3/2}}\longrightarrow\frac12.        \tag{10}
\]

Indeed `U_q=A_q/sqrt(q)` has even real Fourier multipliers. The signing
identity gives `q^(-1)Tr(U_q^2)=1-1/q`. Combined with (9), this forces
the fraction of eigenvalues with `|m_q(a)|<1-eta` to vanish for every
fixed `eta>0`. Eigenvalues larger than `1+eta` are absent eventually.
Thus (6) applies. Equation (7) is the lower bound in (10), and (9) gives
the matching spectral upper bound.

This is a theorem about actual signings, not merely an enlarged operator
model. It excludes every additive-Cayley family over finite fields that
asymptotically attains the Frobenius/spectral lower limit from producing
a uniform cap coefficient below `1/2`. It does not exclude Cayley
families with a nonvanishing normalized spectral excess, or arbitrary
non-Cayley signings.

## 7. Exact-multiplier extension to all finite abelian groups

The exact theorem in Section 1, including asymptotically mean-zero
witnesses, holds uniformly over finite abelian groups `G` with
`|G|->infinity`, when the Fourier multiplier is even, `+1` or `-1` off
the trivial character, and bounded by one at the trivial character.

If `|G|` is even, its dual has a nontrivial character of order two. That
character itself is a mean-zero Boolean eigenfunction, giving the result
exactly.

For odd groups, work in the dual frequency group. If its exponent tends
to infinity, choose a cyclic subgroup of odd order `m` tending to
infinity. The large-characteristic Deuber proof works with `m` in place
of `p`: the iid cosine model on `Z_m` is still uniformly bounded, mean
zero and variance one; the separation bound prevents all modular short
relations once `m>K rN`. Haar averaging on `G` restricts to uniform
averaging on the dual of this cyclic subgroup, so the same finite
moment proof applies.

If the exponent is bounded, sufficiently large group order forces a
large elementary abelian `F_p` subgroup in the dual for some prime in a
fixed finite list. The projective Ramsey proof applies inside it. More
explicitly, bounded exponent and bounded `p`-rank for every such prime
bound the order of a finite abelian group by its cyclic-factor
decomposition. Hence the two cases give a uniform threshold in `|G|`.

Only the **exact** multiplier assertion is extended here. The robust
bad-set argument in Section 5 used scalar transitivity on the nonzero
frequencies of a finite field and is not silently asserted for general
finite abelian groups.

## 8. Original-problem boundary

Translation invariance and the additive-character eigenbasis are the
essential structure: they let a Ramsey frequency configuration become
a common eigenspace for the squarefree Boolean approximant. Nothing
here supplies analogous multiplicative closure for an arbitrary
eigenbasis.

In particular this theorem does not conflict with the Haar-involution
Boolean ceiling, does not prove that all symmetric flat involutions
saturate, and does not prove `M_n/n^(3/2)->1/2`. It is a broad selected-
family obstruction to sub-`1/2` upper constructions, while the original
minimizer landing and convergence obligations remain open.

## 9. Exact finite normalization checks

`computations/resumed_convergence_centered_fourier_exact_verify_2026_09_06.py`
enumerates every even nonconstant multiplier and every projective Boolean
spin on `F_3` and `F_3^2`. Because the conjugate character-pair values in
characteristic three are `2` and `-1`, the matrices `N U` are integral.
The verifier checks symmetry, zero row sums, and
`(N U)^2=N^2 I-N J` exactly before evaluating any cap.

It passes all two spectra at order three and all sixteen spectra at order
nine. Their absolute Boolean Rayleigh maxima are respectively `8/9`, and
either `72/81` or `80/81`. These finite values are not evidence for a
convergence rate; the Ramsey proof supplies the asymptotic statement and
may have a very large accuracy-dependent order threshold.
