# Contextual response congruences: a quantitative law and its sharp boundary

**Status.** Adversarial-verifier report.  The finite results below are
self-contained.  Proposition CRL.1 is the quantitative form of the standard
syntactic/Myhill--Nerode quotient.  Theorem CRL.2 is a universality no-go.
Theorem CRL.3 is a reusable selector-amplification lemma; its coding
consequences use the already proved posterior-width theorem.  None of these
statements concerns the original signing problem.

## 1. The exact future-response pseudometric

Let `(M,star,1)` be a commutative monoid and let `F:M -> R` be a scalar
extremal observable.  An element `a` is a finite system and `c` is a future
context.  Define

```math
R_a(c)=F(a\star c),
\qquad
d_F(a,b)=\sup_{c\in M}|R_a(c)-R_b(c)|.             \tag{CRL.1}
```

The supremum may be replaced by a declared context ideal `C` satisfying
`C star M subset C`.  The complete monoid is used below to avoid extra
notation.

### Proposition CRL.1 (quantitative syntactic congruence)

The function `d_F` is a translation-contractive pseudometric:

```math
d_F(a\star u,b\star u)\le d_F(a,b).               \tag{CRL.2}
```

Consequently

```math
d_F(a\star b,a'\star b')
\le d_F(a,a')+d_F(b,b').                           \tag{CRL.3}
```

The relation

```math
a\equiv_F b\quad\Longleftrightarrow\quad d_F(a,b)=0             \tag{CRL.4}
```

is a monoid congruence, and `M/equiv_F` is, up to injective recoding, the
coarsest exact deterministic state that answers `F(a star c)` for every
future context `c` and remains closed under composition.

#### Proof

For (CRL.2), every test context of `a star u` is of the form `u star c`, a
subset of the contexts already used in (CRL.1).  Apply (CRL.2) twice and the
triangle inequality to obtain (CRL.3).  Equation (CRL.2) makes (CRL.4) a
congruence.  Its class determines every response by definition.  Conversely,
if an exact summary gives the same state to `a,b`, then every decoded future
answer agrees, so `d_F(a,b)=0`. `square`

This is useful bookkeeping, but it is not new mathematics: it is the
real-valued analogue of a syntactic congruence.  In particular, calling this
quotient an “extremal information state” does not by itself prove that it is
smaller than the complete response table.

Equation (CRL.3) also gives the precise approximate-algebra warning.  If two
child representatives have errors `epsilon_1,epsilon_2`, their unrounded
product has error at most `epsilon_1+epsilon_2`; every further projection to
an `epsilon`-net adds another `epsilon`.  Nonexpansiveness prevents
amplification of an existing error, but repeated lossy re-encoding can still
accumulate error with the composition tree.

## 2. Idempotence and monotonicity alone are universal

One might hope that support union is unusually compressible because it is an
idempotent commutative composition and covering radius is antitone under
adding support.  The next theorem rules out any general result based only on
those facts.

### Theorem CRL.2 (arbitrary tables inside a monotone union semilattice)

Let `R=(R_ij)` be an arbitrary finite matrix with entries in `[0,1]`.  There
are

- a finite union semilattice `M`;
- an antitone observable `F:M->[0,1]` which is one-Lipschitz for symmetric-
  difference distance; and
- systems `a_i` and contexts `c_j`;

such that

```math
F(a_i\cup c_j)=R_{ij}                              \tag{CRL.5}
```

for every `i,j`.

#### Proof

Take disjoint symbols `(x_i)_i,(y_j)_j`, let `M` be their power set under
union, and put `a_i={x_i}`, `c_j={y_j}`.  Define

```math
F(A)=
\begin{cases}
1,&|A|\le1,\\
R_{ij},&A=\{x_i,y_j\},\\
0,&\text{otherwise}.
\end{cases}                                       \tag{CRL.6}
```

If `A` is properly enlarged, its value cannot increase: singletons have
value one, distinguished two-sets have values in `[0,1]`, and every strict
superset of a two-set has value zero.  Thus `F` is antitone.  Its range has
diameter one, so it is one-Lipschitz whenever two sets differ.  Equation
(CRL.5) is immediate. `square`

Taking Cartesian powers of the construction and adding the component
observables makes the arbitrary table extensive.  Therefore none of the
following, separately or together, implies response compression:

1. commutative idempotent composition;
2. monotonicity under future composition;
3. a bounded one-step Lipschitz constant; or
4. additivity across independent copies.

Any theorem at that level can inherit the metric entropy, tropical rank, or
communication complexity of an arbitrary matrix.  A nonvacuous theory must
use additional algebraic geometry of the realizable profiles: for example a
low-entropy metric chart, a convolution law, symmetry, or a constrained
family of context incidences.

## 3. A checkable mechanism for microscopic-to-macroscopic exposure

The positive syndrome block theorem has a more general algebraic core.  It
is not idempotence alone; it is the existence of a Boolean cube of contexts
which can select additive latent defects before the final absolute error is
charged.

### Theorem CRL.3 (selector-cube amplification)

Suppose a scalar response experiment contains systems `a_z`, indexed by
`z in {0,1}^q`, and contexts `c_P`, indexed by `P subseteq [q]`, for which

```math
F(a_z\star c_P)
=b(P)+\sum_{j\in P}\lambda_j z_j,
\qquad \lambda_j>0.                               \tag{CRL.7}
```

Then the response metric restricted to these contexts is exactly

```math
d_{\rm sel}(z,z')
=\max\left\{
 \sum_{j:z_j=1,z'_j=0}\lambda_j,
 \sum_{j:z_j=0,z'_j=1}\lambda_j
\right\}.                                        \tag{CRL.8}
```

In particular,

```math
d_{\rm sel}(z,z')
\ge {1\over2}\sum_{j:z_j\ne z'_j}\lambda_j.      \tag{CRL.9}
```

If `P` is uniform over all subsets and `R_z(P)=F(a_z star c_P)`, then

```math
\|R_z-R_{z'}\|_{L^2(P)}^2
= {1\over4}\sum_{j:z_j\ne z'_j}\lambda_j^2
  +
  {1\over4}\left(\sum_j\lambda_j(z_j-z'_j)\right)^2.             \tag{CRL.10}
```

For equal weights `lambda_j=lambda` and `q>=2`, the inverse-Hamming response
modulus is exactly

```math
\kappa={\lambda^2\over4}.                          \tag{CRL.11}
```

Consequently, for a uniform latent vector `Z`, any randomized transcript
whose decoded response has mean-square error `Delta` obeys

```math
I(Z;\text{transcript})
\ge q\left[1-g\left(
 \min\left\{{16\Delta\over\lambda^2q},1\right\}
\right)\right],                                   \tag{CRL.12}
```

where `g` is the function in Theorem 7.1.  Uniform-error lower bounds follow
by applying (CRL.8) to any weighted Hamming packing.

#### Proof

Subtract the two instances of (CRL.7).  A subset sum of signed positive
weights lies between minus the total negative weight and the total positive
weight; choosing exactly the positive or exactly the negative coordinates
attains the two endpoints.  This proves (CRL.8)--(CRL.9).

For independent Bernoulli-half indicators `B_j=1_{j in P}`, expand

```math
\mathbb E\left(\sum_j\lambda_j(z_j-z'_j)B_j\right)^2.
```

The diagonal terms have coefficient `1/2` and the off-diagonal terms have
coefficient `1/4`, giving (CRL.10).  Its first term gives
`kappa>=lambda^2/4`; two coordinates changed in opposite directions attain
equality.  Equation (CRL.12) is Theorem 7.1 with this exact modulus. `square`

The theorem is robust in the uniform metric.  If the right side of (CRL.7)
approximates every response within `tau`, then the lower bound in (CRL.8)
loses at most `2tau`.  No analogous distributional conclusion is automatic
unless the approximation error is controlled in the same query law.

### Application 1: syndrome-support blocks

In Theorem 8.3, put `z_j=1-a_j` and `lambda_j=L-1`.  Its exact formula is

```math
R_z(P)=q+(L-1)|P\cap\{j:z_j=1\}|.
```

Thus (CRL.8) recovers the directed-Hamming response metric, while (CRL.12)
adds a continuous mean-square information law under random appended block
contexts.  The crucial property is not merely support union: the direct-sum
syndrome geometry makes radius additive and the environments implement every
selector `P`.

### Application 2: future-edge Max-Cut on a matching

Let `(e_j)_(j<=q)` be disjoint edges.  The system graph `G_z` contains `e_j`
exactly when `z_j=1`; the future graph `H_P` contains `e_j` exactly when
`j notin P`.  Composition is idempotent edge union.  Since every edge of a
matching can be cut simultaneously,

```math
\operatorname{MaxCut}(G_z\cup H_P)
=q-|P|+\sum_{j\in P}z_j.                          \tag{CRL.13}
```

Theorem CRL.3 applies with unit weights.  On this restricted but standard
optimization model, `q` latent edge bits form a strict composable quotient
of the `2^(2q)`-entry spin landscape, and arbitrary future matching-edge
contexts force their macroscopic response information.  This validation is
deliberately modest: it is not a space lower bound for general Max-Cut.

## 4. Audit of the two syndrome theorems found at this checkpoint

### 4.1 Landmark upper quotient

The proof in `phase3_syndrome_landmark_quotient.md` survives adversarial
checking.

1. A basis `B subset S` gives
   `|lambda_S(x)-lambda_S(y)|<=lambda_S(x-y)<=d_B(x,y)`.
2. An `r`-covering code in the `B`-coordinate cube and the stored landmark
   values therefore approximate the complete profile in sup norm by `r`.
3. Min-plus convolution with an arbitrary, possibly nonspanning, future
   support (using `+infinity` off its span) and subsequent maximization are
   sup-norm nonexpansive.
4. The state-dependent chart is not hidden information: the ordered basis
   costs at most `w^2` bits and lets the decoder interpret the ambient future
   support.
5. Two supports in one summary cell are within `2r` in complete future-
   response distance, so choosing one actual representative per nonempty
   cell gives the claimed actual net.

Thus, at fixed relative error, the stated
`2^((1-h_2(delta)+o(1))w)` **bits** and hence
`exp(o(2^w))` summary states are valid.  The result is a real approximate
congruence, not an estimate of the unperturbed radius alone.

### 4.2 Constant-dimension lower packing

Let `w=2d`, let `W` be a `d`-subspace, choose a complement basis `C_W`, and
put `S_W=(W minus {0}) union C_W`.  The proposed superlinear packing is also
correct.

Writing `G=W direct-sum U_W` shows directly that

```math
lambda_{S_W}(w_0+u)
=\operatorname{wt}_{C_W}(u)+1_{w_0\ne0},
\qquad
rho(S_W)=d+1.                                     \tag{CRL.14}
```

For two subspaces, put `t=dim(W intersect W')`.  The sum `W+W'` has
codimension `t`.  The images of `C_W` span the quotient by `W+W'`, so at
most `t` of them reach the required coset; the residual is a sum of one
element of `W` and one of `W'`.  Hence

```math
rho(S_W\cup S_W')\le t+2.                         \tag{CRL.15}
```

The context `S_W` separates the future-response maps of `S_W,S_W'` by at
least `d-t-1`.

Finally, the number of `d`-subspaces in `F_2^(2d)` is at least `2^(d^2)`.
For fixed `W`, the number with `dim(W intersect W')>d/2` is

```math
\sum_{j<d/2}2^{j^2}{d\brack j}_2^2
\le16d\,2^{3d^2/4}.                              \tag{CRL.16}
```

Here `j=d-dim(W intersect W')` and
`{d bracket j}_2<=4*2^(j(d-j))`.  Greedy packing therefore gives at least
`2^(d^2/4)/(16d)` subspaces with pairwise intersections at most `d/2`.
Their response separation is at least `d/2-1=w/4-1`.  Uniform response error
`epsilon*w` with fixed `epsilon<1/8` consequently costs at least

```math
{w^2\over16}-O(\log w)                            \tag{CRL.17}
```

bits.  The supports have common length `2^d-1+d`; this is an unrestricted-
support lower bound, not a polynomial-length one.

## 5. Director judgment

The response-congruence viewpoint has one useful theorem and one severe
limit.

- **Useful:** selector cubes are a checkable algebraic certificate that
  composition jointly exposes microscopic features.  Their incidence Gram
  matrix yields posterior width, and their uniform directed-subset metric
  yields worst-case packing.  This unifies the syndrome block example and a
  second model without separately paying scalar channels.
- **Limit:** arbitrary response tables already occur inside monotone,
  one-Lipschitz union semilattices.  Therefore “take the syntactic quotient”
  is not a compression principle; in full generality it is the complete
  response landscape under another name.
- **New strongest law from the combined checkpoint:** future responses are
  compressible when their exact profile is Lipschitz on a low-covering-
  entropy interface and every continuation is nonexpansive; they force
  information growth when the realizable contexts contain a large selector
  or packing code.  The syndrome landmark theorem and the subspace packing
  theorem instantiate the two sides with genuinely nontrivial bounds.

This is more than bare convex duality, because it distinguishes the geometry
of the **domain on which an exact compositional profile lives** from the
geometry of its exposed response code.  It is not yet a universal extremal
information theory: the universality theorem proves that no such theory can
ignore model-specific restrictions on those two geometries.

The most important next theorem is to narrow the remaining syndrome gap:
either exploit algebraic constraints on binary word metrics to beat the
generic Hamming-landmark exponent, or construct a response packing with
`2^(Omega(w))` bits rather than the present `Omega(w^2)` bits.  Rephrasing
this as another exact congruence theorem would add no information.
