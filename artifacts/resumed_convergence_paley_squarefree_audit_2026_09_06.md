# Independent audit: all-order Paley saturation by squarefree cosine spectra

Date: 2026-09-06. Independent reconstruction of the director's proposed
Paley-family theorem. This is an actual-signing family result, not a
universal lower bound for arbitrary signings and not a proof that the
original minimum has a limit.

## 1. Precise conclusion

For every odd prime power `q = 1 mod 4`, let

\[
(A_q)_{t,u}=\chi_q(t-u),\qquad t,u\in\mathbb F_q,
\quad \chi_q(0)=0.
\]

Then, along **all** such prime powers tending to infinity,

\[
\max_{x\in\{\pm1\}^{q}}
\frac{x^\top A_qx}{q^{3/2}}\longrightarrow1,
\qquad
\min_{x\in\{\pm1\}^{q}}
\frac{x^\top A_qx}{q^{3/2}}\longrightarrow-1.       \tag{1}
\]

Consequently the Boolean half-energy cap of the Paley core, and of the
bordered symmetric Paley conference of order `q+1`, normalized by the
respective order to power `3/2`, tends to `1/2`.

The proof has no restriction on the characteristic or extension parity
beyond `q = 1 mod 4`. In particular it does not assume that a prime-field
trace map supplies independent coordinates.

## 2. Uniform squarefree Hermite approximation

Let `p` be any odd prime. Let independent `Y_j` be uniform on `F_p`, and put

\[
X_j=\sqrt2\cos(2\pi Y_j/p),\qquad
S_r=r^{-1/2}\sum_{j=1}^rX_j.
\]

Uniformly over these `p`, the variables have mean zero, variance one, and
absolute value at most `sqrt(2)`. For `l <= r`, define the **ordered**
squarefree statistic

\[
W_{r,l}=r^{-l/2}
\sum_{i_1,\ldots,i_l\ \mathrm{distinct}}
X_{i_1}\cdots X_{i_l}.
\]

Writing `He_l` for the monic probabilists' Hermite polynomial, there are
the exact identities

\[
\mathbb E W_{r,l}^2=l!\frac{(r)_l}{r^l},\qquad
\mathbb E[\operatorname{He}_l(S_r)W_{r,l}]
=l!\frac{(r)_l}{r^l}.                              \tag{2}
\]

For the first identity, the unordered index sets in the two factors must
coincide, after which every coordinate occurs twice. For the second, all
terms of `He_l` below its leading monomial have degree less than `l` and
cannot cover every index in `W_{r,l}`. In `S_r^l`, precisely the permutations
of the same index set survive. Also `W_{r,l}` and `W_{r,m}` are orthogonal
when `l != m`.

The moment method for sums of independent, bounded, mean-zero,
variance-one variables gives, for each fixed moment order, convergence to
the corresponding Gaussian moment uniformly over `p`. This can be seen
directly by sorting the index tuples by their multiplicity partition:
pairings supply the Gaussian term, while every other surviving partition
has at most one fewer half-power of `r`. The bound `|X_j| <= sqrt(2)` makes
the remainder uniform. Thus (2) implies

\[
\sup_{p\ \mathrm{odd}}
\|\operatorname{He}_l(S_r)-W_{r,l}\|_2\longrightarrow0.
                                                               \tag{3}
\]

Choose a finite **odd** Hermite polynomial

\[
h(s)=\sum_{\substack{l\le D\\l\ \mathrm{odd}}}
b_l\frac{\operatorname{He}_l(s)}{\sqrt{l!}}
\]

approximating `sign(s)` in Gaussian `L2`. Such polynomials exist, for
example from its coefficients

\[
b_{2d+1}=\sqrt{2/\pi}\,
\frac{(-1)^d(2d-1)!!}{\sqrt{(2d+1)!}}.
\]

Their squared sum is one by the Taylor series for `arcsin(1)`, so the
finite orthogonal projections have arbitrarily small error.

An important continuity point is that, with `sign(0)=1`,

\[
g(s):=(\operatorname{sign}(s)-h(s))^2
=1+h(s)^2-2\operatorname{sign}(s)h(s)               \tag{4}
\]

is continuous even at zero: oddness gives `h(0)=0`. Uniform convergence in
distribution of `S_r` to a standard Gaussian, together with the uniform
bounded moments of higher order, therefore yields

\[
\sup_p\left|\mathbb E g(S_r)-\mathbb E g(G)\right|
\longrightarrow0.                                \tag{5}
\]

There is no need for an independent anti-concentration or tie argument.
Set

\[
P_r(X)=\sum_{\substack{l\le D\\l\ \mathrm{odd}}}
b_l\frac{W_{r,l}}{\sqrt{l!}}.
\]

Equations (3)--(5) imply: for every `epsilon>0`, one can choose fixed
`D,r` so that, **uniformly for all odd primes `p`**,

\[
\|\operatorname{sign}(S_r)-P_r\|_2<\epsilon.         \tag{6}
\]

Every Fourier frequency of `P_r` on `F_p^r` is a vector in
`{0,+1,-1}^r` with nonempty odd support of size at most `D`. This follows
by expanding each squarefree factor as
`X_j = (exp(2 pi i Y_j/p)+exp(-2 pi i Y_j/p))/sqrt(2)`.

## 3. Character patterns and forbidden short relations

Fix `r,D`. Choose one representative modulo overall sign from the finite
set of squarefree frequency vectors just described. Call the resulting
linear forms on `F_q^r` `L_1,...,L_m`.

They remain pairwise nonproportional in **every odd characteristic**.
Indeed any proportionality scalar is determined by a nonzero coordinate
and must be `+1` or `-1`; the support and the remaining coordinates then
force equality up to overall sign. When `q = 1 mod 4`, both signs of each
frequency have the same quadratic character.

For any prescribed signs `sigma_1,...,sigma_m`,

\[
\#\{a\in\mathbb F_q^r:
\chi_q(L_j(a))=\sigma_j\ \forall j\}
=2^{-m}q^r+O_{r,m}(q^{r-1/2}).                    \tag{7}
\]

Here the implied bound is uniform in the odd characteristic. For clarity,
the ordinary one-variable character bound is sufficient:

* Choose a direction `v` with every `L_j(v)` nonzero, possible for `q>m`.
* Decompose `a=t v+w`, with `w` in a fixed complementary hyperplane.
* For all but `O_m(q^{r-2})` such `w`, the roots of the linear factors
  `L_j(tv+w)` are pairwise distinct. A collision for a specified pair is
  one nonzero linear equation on that complement.
* Any nonempty subproduct is then squarefree of bounded degree. Its
  quadratic-character sum over `t` is `O_m(sqrt(q))`. The exceptional
  lines contribute only `O_m(q^{r-1})` in total.
* Expand `prod_j (1+sigma_j chi_q(L_j(a)))`; remove the points on the
  `m` zero hyperplanes, whose count is `O_m(q^{r-1})`.

This proves (7). We may assume `r>=2`; for `r=1` the only relevant form
class has a direct exact count.

The imported quadratic character bound is the squarefree specialization
of Weil's bound. The precise statement was checked in the published
primary research paper by Kim, Yip, and Yoo,
[*Explicit constructions of Diophantine tuples over finite fields*,
Lemma 2.1 (2024)](https://link.springer.com/article/10.1007/s11139-024-00888-5#Sec2).
For a nonzero constant times a product of `d` distinct linear factors
over an odd finite field, it gives `(d-1)sqrt(q)`. On every good line
above, every nonempty subproduct has exactly these hypotheses. The
constant multiplier has modulus-one character, and the exceptional
lines are charged separately. No fixed-characteristic asymptotic from
that paper is imported. Kopparty's lecture derivation in
[*The Weil bounds*, Section 2.3](https://sites.math.rutgers.edu/~sk1233/courses/finitefields-F13/weil.pdf)
was also consulted as an expository check, not as the primary source.

For any further fixed integer `K`, also require

\[
\sum_{j=1}^r\ell_j a_j\ne0
\quad\text{for all }\ell\in\mathbb Z^r,
\quad\sum_j|\ell_j|\le K,
\quad \ell\not\equiv0\pmod p.                    \tag{8}
\]

Each condition excludes a genuine `F_q`-linear hyperplane. There are
finitely many, bounded solely in terms of `r,K`, so their union costs
`O_{r,K}(q^{r-1})`. Consequently (7) and (8) are simultaneously feasible
for every sufficiently large admissible `q`, with the threshold depending
only on `r,D,K` and not on `p`.

## 4. Finite moments, not full trace-label independence

Choose `a` as above and, for uniform `t in F_q`, put

\[
\widetilde X_j(t)=\sqrt2\cos
\left(\frac{2\pi}{p}\operatorname{Tr}_{q/p}(a_jt)\right),
\qquad \widetilde S_r=r^{-1/2}\sum_j\widetilde X_j.
\]

Condition (8) implies that the expectations of every polynomial in the
`r` coordinates of total degree at most `K` agree **exactly** with the iid
`F_p` coordinate model. Expand a monomial into additive characters. Its
field expectation is nonzero precisely when `sum ell_j a_j=0`; (8) makes
this equivalent to every `ell_j=0 mod p`, exactly the iid criterion.

Now fix `D,r` from Section 2. The polynomial
`(h(S_r)-P_r(X))^2` has degree at most `2D`. Its expectation therefore
transfers exactly once `K>=2D`.

The function `g` in (4) is continuous on the fixed compact interval
`[-sqrt(2r),sqrt(2r)]`. Choose an ordinary scalar polynomial `R` uniformly
approximating it to error `delta`. Increasing the still finite `K` to the
degree of `R`, its expectation also transfers exactly, whence

\[
\left|\mathbb E_tg(\widetilde S_r(t))-
\mathbb E g(S_r)\right|\le2\delta.                 \tag{9}
\]

Combining this with the two-term triangle inequality gives an arbitrarily
small `L2` error between the genuinely Boolean function
`f(t)=sign(tilde S_r(t))` and `P_r(tilde X(t))`.

This reasoning works even when the trace-label map cannot be surjective
onto `F_p^r`. It never asserts full independence. All choices obey the
order `epsilon -> D -> r -> K -> q_0`, and then cover every `q>=q_0`.

## 5. Landing in either actual Paley eigenspace

Use normalized counting measure on `F_q`. The real symmetric operator
`U_q=A_q/sqrt(q)` has norm one, annihilates constants, and has additive
Fourier multiplier

\[
\tau_q\chi_q(b)\quad(b\ne0),\qquad \tau_q\in\{\pm1\}.
\]

This uses the quadratic Gauss-sum identity and `q=1 mod 4`; choosing the
global sign `tau_q` explicitly avoids any extension-parity convention.

To target eigenvalue `lambda in {+1,-1}`, prescribe all signs in (7) to
be `sigma_j=lambda tau_q`. Every frequency of `P_r(tilde X)` is then
nonzero and lies in the `lambda` eigenspace. No short relation can create
a constant term. Let `P=P_r(tilde X)` and arrange `||f-P||_2<epsilon`.
Orthogonal projection shows that the combined squared mass of `f` outside
the target eigenspace is at most `epsilon^2`. Since `||f||_2=1`,

\[
\lambda\langle f,U_qf\rangle\ge1-2\epsilon^2.       \tag{10}
\]

This proves both assertions in (1). It also supplies witnesses with mean
at most `epsilon` in absolute value, since the target eigenspace is
orthogonal to constants.

For the bordered conference matrix, the border contribution to the
half-energy is at most `q` in absolute value. It can be aligned by the
border spin, but even the crude `O(q)` bound suffices. The exact spectral
upper bound is `(q+1)sqrt(q)/2`. Dividing by `(q+1)^(3/2)` completes the
conference version.

## 6. Archive comparison and original-problem boundary

The older `prime_paley_cosine_saturation.md` proves saturation only along
an arithmetic subsequence. `paley_resonance_gadget.md` improves this to a
saturating subsequence inside every fixed finite character cylinder, and
explicitly states that it does not cover every prime.
`even_extension_paley_trace_fibre.md` gives the selected-family lower bound
`(1-1/p^2)/2` on `q=p^(2r)` using a single trace coordinate.

The present argument removes the varying-order arithmetic selection:
the frequencies themselves are selected afresh inside every large field.
The squarefree approximation prevents multiple harmonics of a frequency
from imposing incompatible Legendre symbols.

This rules out a fixed asymptotic improvement below `1/2` using native
Paley cores or conferences, even on an arbitrarily chosen subsequence.
It does **not** prove that a positive-density principal restriction of a
Paley matrix still saturates, and does not settle any near-minimizer
landing clause for arbitrary signings. The original universal interval
and convergence status are unchanged.
