# Operator-valued rooted certificates for the augmented cut code

Date: 2026-08-15.  Independent research draft; not a tracked project file.

## 1. Outcome

There is an exact operator-valued extension of the rooted annular theorem.
It keeps all matrix channels joint in Loewner order and can, in principle,
use a root-dependent favorable eigendirection.  This extension is outside the
rank-one partial-transversal no-go.

Three restrictions substantially narrow the opening.

1. The obvious operator lift of the OpenAI moving projection,
   `P U_g P U_g^* P`, is pointwise positive but is not generally of operator
   positive type.  The Hilbert--Schmidt trace in the published construction is
   essential, not cosmetic.
2. A flat vector bundle over the ambient-cube/cut-code quotient is governed by
   the abelian deck group `D=C^perp`.  It decomposes into scalar character
   twists, so it cannot create noncommutative cancellation between deck
   channels.
3. Any bounded-multiplicity section reaching `lambda >= c/sqrt(n)` needs
   hidden Fourier support `2^(Omega(n log n))`.  If it uses at most `r`
   representatives per syndrome, then `r=2^(Omega(n log n))`.

Thus standard matrix-valued positive type, flat-bundle, and ordinary
Terwilliger implementations do not provide a polynomial-state escape.  The
remaining possibility is a non-flat, moving operator section with
root-dependent spectral directions and an algebraic proof of a uniform
generalized-eigenvalue or determinant inequality.

## 2. Exact operator annular inequality

Let

```math
H=F_2^E,\qquad E={n\choose2},\qquad C=C_n^+,\qquad D=C^\perp.
```

Write the group multiplicatively in sign coordinates and put

```math
\tau(g)={1\over E}\sum_e g_e,
\qquad
\mu(a)=\max_{c\in C}\tau(ac).
```

Let `V` be a finite-dimensional Hilbert space.  Suppose
`K:H -> Herm(V)` satisfies all three conditions

```math
K \text{ is of operator positive type},
\qquad K(g)\succeq0\quad(g\in H),
\qquad F(g):=(\tau(g)-\lambda)K(g)
       \text{ is of operator positive type}.                 \tag{2.1}
```

For a finite abelian group, operator Bochner says equivalently that

```math
K(g)=\sum_{S\subseteq[E]}\chi_S(g)A_S,
\quad A_S\succeq0,
```

and likewise every Fourier coefficient of `F` is positive semidefinite.
Define

```math
T_a=\sum_{c\in C}K(ac),
\qquad
J=\sum_{c\in C}F(c).                                  \tag{2.2}
```

### Theorem 2.1 (operator rooted annular inequality)

For every `a` and every `0<=s<lambda`,

```math
\mu(a)\le s
\quad\Longrightarrow\quad
(\lambda-s)T_a\preceq J.                             \tag{2.3}
```

Consequently any one of the following proves `mu(a)>s`:

```math
(\lambda-s)T_a\npreceq J,                             \tag{2.4}
```

```math
\lambda_{\max}
 \left(J^{\dagger/2}T_aJ^{\dagger/2}\right)
>{1\over\lambda-s},                                  \tag{2.5}
```

with the usual support convention, or, when both matrices are positive
definite,

```math
(\det T_a)^{1/r}>{(\det J)^{1/r}\over\lambda-s},
\qquad r=\dim V.                                      \tag{2.6}
```

Proof.  Operator Bochner and character orthogonality give

```math
S_a:=\sum_{c\in C}F(ac)
=|C|\sum_{d\in D}\chi_d(a)\widehat F(d),
```

so `-J \preceq S_a \preceq J`; both `J+S_a` and `J-S_a`
are sums of positive matrices.  If `mu(a)<=s`, pointwise positivity of
`K` gives

```math
S_a=\sum_c(\tau(ac)-\lambda)K(ac)
\preceq-(\lambda-s)T_a.
```

Combining this with `-J \preceq S_a` proves (2.3).  Loewner order implies
(2.5), and determinant monotonicity gives (2.6).  QED.

This theorem does not separately pay matrix channels.  The determinant form
can use different bad directions for different roots.  On the other hand, a
proof based on one fixed positive functional `W` immediately scalarizes to
the scalar kernel `tr(WK)`.  A genuinely matrix-valued gain therefore requires
that the separating direction vary with `a`, or a joint spectral invariant
such as (2.6).

The exact convergence-scale target is to construct (2.1), independently of
`M_n`, with

```math
\lambda_n={1-o(1)\over\sqrt n}
```

and prove (2.4), uniformly for every `a`, at
`s_n=(1-o(1))/sqrt(n)`.  Then
`Q(a)>=E s_n=(1/2-o(1))n^(3/2)`.

## 3. The natural moving-projection operator lift fails

Let `U_g` be the diagonal unitary representation of the edge cube and let
`P` be a moving-projection base fiber.  The natural rooted overlap operator is

```math
K_P(g)=P U_gP U_g^*P=(P U_gP)^2\succeq0.             \tag{3.1}
```

It need not be of operator positive type.  There is already an exact
counterexample on `F_2^2`.  Index the four coordinates by its characters,
let `U_g` be the corresponding diagonal character representation, and let
`P` be projection onto the columns of

```math
X=\begin{pmatrix}
0&0\\0&1\\1&-1\\-1&-1
\end{pmatrix}.
```

Then

```math
P=\begin{pmatrix}
0&0&0&0\\
0&1/3&-1/3&-1/3\\
0&-1/3&5/6&-1/6\\
0&-1/3&-1/6&5/6
\end{pmatrix}.
```

For the nontrivial character `eta=(0,1)`, the operator Fourier coefficient
of (3.1) is

```math
\widehat K_P(\eta)=
\begin{pmatrix}
0&0&0&0\\
0&-1/27&1/27&1/27\\
0&1/27&5/108&-13/108\\
0&1/27&-13/108&5/108
\end{pmatrix},
```

whose eigenvalues are `0,0,1/6,-1/9`.  Hence (3.1) fails operator
Bochner positivity.

Taking the Hilbert--Schmidt trace removes this problem: its Fourier
coefficient is the sum of squared matrix entries and is nonnegative.  This
explains exactly why the published moving-projection proof closes only after
scalar trace.  An operator certificate needs a different construction, not
merely retaining `P_uP_vP_u` before tracing.

## 4. Fourier-support and rank obstruction

Write the operator Fourier coefficients in (2.1) as `A_S`.  Since
multiplication by `tau` is normalized cube adjacency, positivity of the
Fourier coefficients of `F` says

```math
{1\over E}\sum_e A_{S\triangle\{e\}}\succeq\lambda A_S
\qquad(S\subseteq[E]).                               \tag{4.1}
```

Set `w_S=tr(A_S)>=0` and `U={S:w_S>0}`.  Taking traces and then the
Rayleigh quotient gives

```math
\lambda
\le {\rho(Q_E[U])\over E}.                           \tag{4.2}
```

Indeed, multiply the traced version of (4.1) by `w_S`, sum in `S`, and
use that `w` is supported on `U`.

The Bollobas--Lee--Letzter cube-subgraph theorem now applies without regard
to matrix commutativity or rank.  If `lambda>=c/sqrt(n)`, choose

```math
i=\left\lfloor {c^2E\over16n}\right\rfloor=\Theta(n).
```

If `|U|<=sum_(j<=i) binom(E,j)`, their asymptotic extremal estimate would
give

```math
{\rho(Q_E[U])\over E}
\le(2+o(1))\sqrt{i/E}
\le {c+o(1)\over2\sqrt n},
```

contradicting (4.2).  Therefore

```math
|U|\ge
2^{(c^2/32+o(1))n\log_2 n}.                          \tag{4.3}
```

The constant is deliberately nonoptimized; the `2^(Omega(n log n))`
conclusion is the important part.

Suppose a vector-valued section uses at most `r` Fourier representatives in
each of the `|H/D|=|C|=2^n` syndrome fibers.  Then `|U|<=r2^n`, so (4.3)
forces

```math
r\ge2^{(c^2/32+o(1))n\log_2n-n}
 =2^{\Omega(n\log n)}.                               \tag{4.4}
```

If the final kernel is a quadratic product of such a section, replace `r`
by `r^2`; the conclusion remains `r=2^(Omega(n log n))`.

Thus polynomial-rank or `exp(O(n))` multi-representative sections cannot
reach the leading scale.  The obstruction is compatible with the OpenAI
construction: degree `Theta(n)` Boolean harmonic fibers have precisely
`exp(Theta(n log n))` hidden dimension.

## 5. Flat vector bundles reduce to scalar twists

The cube projection `H -> H/D` is a regular graph cover with abelian deck
group `D`.  Choose a section `s:H/D -> H`.  Its edge discrepancies are the
voltage cocycle

```math
\delta_e(x)=s(x)+e-s(x+\pi(e))\in D.                 \tag{5.1}
```

A rank-`r` flat bundle is specified by a unitary representation
`rho:D -> U(r)`, and its twisted quotient adjacency is

```math
(A_\rho f)(x)={1\over E}\sum_e
 \rho(\delta_e(x))f(x+\pi(e)).                       \tag{5.2}
```

Since `D` is abelian, every unitary representation simultaneously
diagonalizes:

```math
\rho\cong\chi_1\oplus\cdots\oplus\chi_r.
```

Consequently

```math
A_\rho\cong A_{\chi_1}\oplus\cdots\oplus A_{\chi_r},
\qquad
\rho(A_\rho)=\max_j\rho(A_{\chi_j}).                 \tag{5.3}
```

No cross-character cancellation survives.  Equivalently, a
`D`-equivariant vector section decomposes into `D`-character weight spaces.
Each character has exactly `|H/D|=2^n` extensions to an ambient cube
character, so a rank-`r` flat section has Fourier support at most `r2^n`,
recovering (4.4).

An additional `S_n` action may permute the `D`-character spaces, but cube
adjacency commutes with the deck action and remains block diagonal in those
spaces.  Thus ordinary `S_n`-equivariant flat bundles are still scalar
twists in disguise.  A viable moving section must be non-flat or must couple
root and representation before passing to a deck-equivariant adjacency.

## 6. Why standard Terwilliger reduction does not transfer

Schrijver's Hamming Terwilliger algebra is polynomial-dimensional because
the root stabilizer is the full coordinate permutation group `S_E`; its
dimension is `binom(E+3,3)`.  The 2026 Gijswijt--Polak covering SDP uses the
same full Hamming automorphism symmetry and explicitly notes that practical
higher levels grow rapidly.

The augmented cut code is not preserved by `S_E`.  For `n>=5`, its
coordinate-permutation automorphism group is the automorphism group of the
graphic matroid of `K_n`, namely `S_n`.  After quotienting by cut
translations, the rooted space has

```math
|H/C|=2^{E-n}
```

switching classes.  The number of `S_n` orbits is at least

```math
{|H/C|\over n!}
=2^{E-n-\log_2(n!)}
=2^{n^2/2-O(n\log n)}.                               \tag{6.1}
```

Any full rooted orbit/Terwilliger algebra containing the diagonal orbit
idempotents has at least this dimension.  Averaging under `S_E` recovers the
polynomial Hamming algebra only by erasing the cut-code structure.  Hence
standard Terwilliger block diagonalization does not supply a compact rooted
certificate here; a small subalgebra would itself be a new closure theorem.

The recent covering-code SDP also solves a different quantifier problem: it
lower-bounds the size of an arbitrary radius-`r` covering.  Our task is to
upper-bound the radius of one fixed highly structured code.  Its matrix-cut
constraints therefore do not imply (2.4)--(2.6).

## 7. Exact remaining target and falsifiers

A candidate that genuinely survives these audits must provide matrices
`A_S` with

```math
A_S\succeq0,
\qquad
\sum_S\chi_S(g)A_S\succeq0\quad(g\in H),
\qquad
{1\over E}\sum_eA_{S\triangle e}-\lambda A_S\succeq0, \tag{7.1}
```

at `lambda=(1-o(1))/sqrt(n)`, and prove a uniform joint spectral inequality
such as (2.5) or (2.6) for

```math
T_a=|C|\sum_{d\in D}\chi_d(a)A_d,
```

with

```math
J=|C|\sum_{d\in D}
 \left({1\over E}\sum_eA_{d\triangle e}-\lambda A_d\right).
```

Reject it if:

1. all `A_S` commute (then it is a direct sum of scalar channels);
2. its separating positive functional is fixed independently of `a`
   (then `tr(WK)` is already a scalar certificate);
3. it is a flat `D`-bundle (Section 5);
4. its hidden support is below (4.3);
5. verifying the uniform spectral inequality enumerates all switching
   classes or reconstructs the cut-energy maximum.

No primary theorem located in Schrijver's Terwilliger work, the recent
Gijswijt--Polak covering SDP, operator-valued Bochner theory, or abelian graph
cover/Floquet theory supplies (7.1) with the needed rooted inequality.

