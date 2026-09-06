# Primary-literature audit after the 0.3857858769 certificate

Date: 2026-09-05. Status: no published universal lower bound stronger than
the new certificate was successfully imported in this search. This is a
scoped search result, not a claim of exhaustive bibliographic absence.
The most exact actionable predecessor found is the degree-two case of the
sharp total-L1-influence conjecture.

Final scope update: this is a historical search checkpoint. Its numeric
comparisons with `0.3857858769` and the proposed `0.39894`/`0.42045`
landings are not the campaign's current lower bound or current
improvement targets. They have been overtaken by the actual-sign
weighted theorem and normalized gain in
`fresh_normalized_unmarked_gain_and_limit_order_2026_09_05.md` and
its later uniform scalar-family extension. No conjecture or bilinear
comparison in this note is promoted to an additional universal theorem.

Throughout,

`P_A(x)=sum_(i<j) a_ij x_i x_j`, `q(A)=||P_A||_(Boolean infinity)`,

with a complete hollow symmetric signing `A`.

## 1. A precise nearby conjecture, not a theorem to import

Primary source: Yuval Filmus, Hamed Hatami, Nathan Keller, Noam Lifshitz,
*On the sum of the L1 influences of bounded functions*, Israel J. Math.
214 (2016), 167–192; [author manuscript](https://arxiv.org/pdf/1404.3396),
March 28, 2015 version. Theorem 3.3, Theorem 3.6, Proposition 3.7, the
quadratic counterexample following Proposition 3.7, and Section 4 were
read in the actual paper.

Their total L1 influence is

`I1(f)=sum_i E |(f(S)-f(S with spin i flipped))/2|`.

Theorem 3.3 proves `I1(f)<=d^2 ||f||infinity`; Theorem 3.6 gives
`O(d log d)` for homogeneous functions. Section 4 explicitly conjectures
the sharp inequality `I1(f)<=d ||f||infinity` for all bounded real-valued
degree-d functions. No later proof of the sharp statement was located.

For our exact signings there is no matrix-dependent error or invariance
issue in computing the left side:

`I1(P_A)=n E |epsilon_1+...+epsilon_(n-1)|`.

Every row has exactly `n-1` coefficients of modulus one, so multiplying
the input signs by those coefficients gives this identity separately in
every row. The one-dimensional Rademacher CLT and bounded second moments
give

`I1(P_A)=(sqrt(2/pi)+o(1)) n^(3/2)`.

Consequently their conjecture just for homogeneous degree two would prove

`liminf M_n/n^(3/2) >= 1/sqrt(2pi) = 0.3989422804014327...`.

More weakly, an asymptotic bound `I1(P_A)<=C q(A)+o(n^(3/2))` on this
complete-sign subclass would beat 0.3857858769 whenever
`C<2.0682057290...`. The sharp conjecture supplies the especially simple
candidate `C=2`. This is a new theorem target, not an imported result, and
would improve the lower bound without proving convergence.

Two misleading special cases are excluded explicitly. Proposition 3.7
has constant d only for **Boolean-valued** homogeneous functions, not
arbitrary bounded real-valued functions. The bounds for **symmetric
functions** in Section 3.3 mean invariance under all input-coordinate
permutations, not symmetry of the coefficient matrix. Furthermore the
pointwise L1-gradient bound with constant 2 is false even for homogeneous
quadratics: their difference of two normalized clique-square polynomials
has gradient norm tending to 4 at the all-one vertex. Any successful
constant-2 proof must retain averaging or the complete-sign structure.

## 2. Recent Bohnenblust–Hille results do not remove the norm gap

### Exact real anisotropic Littlewood constants

N. Caro-Montoya, D. Nunez-Alarcon, D. Serrano-Rodriguez,
*The sharp constants in the real anisotropic Littlewood's 4/3 inequality
and applications*, [arXiv:2407.06804](https://arxiv.org/pdf/2407.06804).
Actual Theorem 1.5 states that for `a,b in [1,infinity]`,
`1/a+1/b<=3/2`, every continuous real bilinear form satisfies

`||A||_(a,b) <= 2^max(0,1/a+1/b-1) ||A||`,

and the constant is sharp. The right norm has independent inputs. For
the symmetric hollow matrix in this project, it is `beta(A)`, not
`2q(A)`. The always-valid elementary comparison is `beta(A)<=4q(A)`.
Thus the exact bilinear constant sqrt(2) is not a same-spin quadratic
constant sqrt(2).

### Complex symmetric orbit compression

D. Nunez-Alarcon, D. Pellegrino, A. Raposo Jr., E. Teixeira,
*Orbit compression and asymptotic contractivity for symmetric
Bohnenblust–Hille inequalities*,
[arXiv:2608.13753](https://arxiv.org/pdf/2608.13753).
Theorem 1.1 concerns **complex** symmetric m-linear forms and gives an
ordered-coefficient norm bound by `(1+C/m)` times their diagonal
polynomial norm for sufficiently large **degree m**. It is not a
fixed-degree-two result as the number of variables grows. Section 8,
read through its conclusion, constructs real symmetric forms with norm
one and ordered-coefficient norm sqrt(2) for every even degree, expressly
excluding a direct real analogue of the asymptotic contractivity theorem.

### Constant-one completely bounded BH and exact Boolean-valued BH

S. Arunachalam, A. Dutt, F. Escudero Gutierrez, C. Palazuelos,
*A cb-Bohnenblust–Hille inequality with constant one and its applications
in learning theory*, Math. Ann. 392 (2025), 3367–3396;
[primary journal text](https://link.springer.com/article/10.1007/s00208-025-03142-5).
Actual Theorem 1.1 replaces the scalar supremum by the **completely
bounded norm**. It gives no constant-one comparison back to our scalar
norm. Actual Proposition 4.2 gives BH constant `2^((d-1)/d)` for
**Boolean-valued** degree-d functions. Its proof uses granularity
`hat f(S) in 2^(1-d) Z` and Parseval, implying at most `4^(d-1)` nonzero
coefficients. For d=2 this means at most four coefficients. Our normalized
complete-sign polynomial has `binom(n,2)` nonzero coefficients and is not
in that class for n>=4. This is a decisive hypothesis failure, not a
possible interpolation.

### Genuine Boolean-domain BH

A. Defant, M. Mastylo, A. Perez, *On the Fourier spectrum of functions on
Boolean cubes*, Math. Ann. 374 (2019), 653–680;
[arXiv:1706.03670](https://arxiv.org/pdf/1706.03670).
Theorem 1.1 is genuinely for real-valued functions on the Boolean cube,
but proves a dimension-free bound `C^(sqrt(d log d))` with an unspecified
absolute C. Its asymptotic parameter is the degree. It does not state the
sharp degree-two constant needed here.

For comparison, a genuine Boolean quadratic BH inequality with constant
`C2` gives `q(A)>=binom(n,2)^(3/4)/C2`. It beats the new certificate only
if `C2<1.5412787069...`; the tempting but unproved same-spin import
`C2=sqrt(2)` would give `2^(-5/4)=0.4204482076...`.

## 3. Graph discrepancy and signed-graph scope

B. Bollobas and A. Scott, *Discrepancy in graphs and hypergraphs*,
[author paper](https://people.maths.ox.ac.uk/~scott/Papers/disc.pdf),
establish an actual two-sided product estimate

`disc_p^+(G) disc_p^-(G) >= p(1-p)n^3/6400`

when `p(1-p)>=1/n`. Their discrepancy is over induced vertex sets with
centering at the graph's own density p. The explicit constant is much
too small to supersede 0.38578 after any direct norm mapping. The product
architecture does control both sides, but the stated result does not
give a sharp same-spin constant or an all-order optimal-value limit.

E. Raty, B. Sudakov, I. Tomon, *Positive discrepancy, MaxCut, and
eigenvalues of graphs*, Trans. AMS 379 (2026), 2111–2140;
[primary author copy](https://people.math.ethz.ch/~sudakovb/positive-discrepancy-and-eigenvalues.pdf),
Theorems 1.2 and 1.5 were checked. These concern one-sided centered
discrepancy for unsigned graphs and contain density restrictions away
from one half, or yield a dense-graph scale `n^(5/4)/log n`, not a
universal n^(3/2) constant. Their semidefinite comparison in Lemma 4.2
also loses log n in the one-sided setting. These precise hypotheses do
not supply the desired improvement.

Ordinary signed-graph frustration is not the absolute objective. If
`ell(A)` denotes the minimum number of negative edges in a switching
class and `N=binom(n,2)`, then exactly

`q(A)=N-2 min(ell(A),ell(-A))`.

Thus a theorem only maximizing `ell(A)`, or only studying the all-negative
signature, misses the necessary simultaneous balance/antibalance tradeoff.
Likewise, Hadamard maximal excess with independent row and column
switching remains a bilinear quantity; symmetry of the matrix does not
identify the switches.

## 4. Bottom line and actionable next statement

The strongest concrete external landing statement found in this pass is
the averaged, degree-two L1-gradient inequality in Section 1. Its exact
sign-specific left side makes it substantially more focused than asking
for an unspecified improvement to Grothendieck or Bohnenblust–Hille.
However, it is currently a conjectural ingredient in this audit.
None of the primary theorems checked establishes convergence of the
original minima or superseded the then-current 0.3857858769 certificate.
The later stronger universal results are internal proved developments,
not imports from the neighboring conjectures listed here.
