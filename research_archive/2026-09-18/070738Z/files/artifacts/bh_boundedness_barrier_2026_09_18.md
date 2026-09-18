# Boolean BH boundedness: weighted-bootstrap barrier diagnosis

2026-09-18, bounded campaign 04:27:50--07:27:50 UTC. Independent barrier
track. This file is a research artifact, not a claim that the BH
boundedness problem has been solved. Statements labeled **proved here**
have a complete derivation below. Independent audits are identified
locally; no external novelty or priority claim is made.

Current handoff map:

- Sections 2--4 isolate a genuine linear obstruction in the weighted
  auxiliary norm, not a divergent BH family. Section 18 rejects the
  proposed final-exponent repair at its false coefficient-`l1` base.
- Sections 7--8 give the torsion support bound and an exact obstruction
  to torsion-phase approximation. Sections 9 and 14 independently audit
  the counterexample track's finite ratio above two and the all-order
  result `liminf B_m>2+2^(-207)`.
- Sections 15--16 and 26 give specific chain, symmetry, and
  noncommutative-scalarization barriers with their scopes stated.
- Sections 21, 23, and 29 separate uniform tensor degree, full operator
  velocity, and cap-bounded `L2` degree by explicit proved examples.
- Sections 28 and 31--32 prove positive phase-entropy estimates. The
  final flat-diagonal certificate in Section 32 closes the fixed-seed
  real-quadratic-phase entropy mechanism even for cap-bounded `L2`
  compilers; it does not settle unrestricted BH boundedness.
- Sections 22 and 30 audit the positive track's all-unitary junta and
  sharp separately-affine unitary results. Failed extensions and
  independent verification details remain alongside the proofs.

## 1. Scope and primary-source reconstruction

For complex Walsh polynomials of degree at most `m`, let

```math
q_m=\frac{2m}{m+1},\qquad
B_m=\sup_{n,f\ne0,\deg f\le m}
 \frac{\|\widehat f\|_{q_m}}{\|f\|_\infty},\qquad
L_r(f)=\|\widehat f^{=r}\|_{q_r}.
```

The target is `sup_m B_m<infinity` or an actual diverging family, not an
improvement in an upper-bound power. No conclusion either way is obtained
by the barriers below.

Primary papers read:

- Ivanisvili, [arXiv:2609.12427v1](https://arxiv.org/html/2609.12427v1).
- Slote--Volberg, [arXiv:2609.07758v2](https://arxiv.org/html/2609.07758v2).

The latter's Boolean specialization has exponent 7. Its high-degree
engine is the same restriction-closed weighted square function, and its
polynomial exponent decomposes as `4+3`: removal of degree weights plus
the first three projection estimates. The earlier Ivanisvili exponent
decomposes as `5+22`. The existing repository refinement gives `5/2+3`.
These statements concern the explicit displayed proofs, not priority.
The user-reported subsequent social-media claims of exponents 2--3 and a
barrier at 1 are treated as reports, not as independently available proofs.

The main sources of loss can be separated as follows.

1. Replacing the desired level coefficient norm at `q_m` by the stronger
   `L_r`, at exponent `q_r`. This loses essential diffuse-coefficient
   information when `r` is much smaller than `m`.
2. Reassembling the finitely many levels costs only `m^(1/(2m))`.
3. Low-level projection bounds introduce polynomial growth in `m`.
4. Random capture and block-count factors approach 1 for suitably wide
   central windows at large degree.
5. Full-row `q_d -> 2` hypercontractivity has a damping-removal cost
   tending to `e`, even for perfectly balanced large bidegrees.
6. Degree weights offset that cost; removing those weights afterward
   incurs a separate polynomial factor.

The last three items do not explain the exponent-one obstruction by
themselves. The first item does: the auxiliary functional is intrinsically
at least linear in `m`.

## 2. Sharp real first-level obstruction

**Proved here.** For every odd integer `m`,

```math
\sup_{n,\deg f\le m,\ f\text{ real},\ \|f\|_\infty\le1} L_1(f)=m.
```

For the upper bound, `L_1(f)=||f^{=1}||_infinity` for real coefficients.
For each vertex `x`, the polynomial `p_x(t)=f(tx)` has degree at most
`m` and is bounded by one on `[-1,1]`, because a multilinear polynomial
has the same supremum on the real continuous cube and its vertices.
Bernstein's interior derivative inequality gives `|p_x'(0)|<=m`.
Its derivative is `f^{=1}(x)`, proving the upper bound. This use of
Bernstein can equivalently be phrased for the trigonometric polynomial
`p_x(cos theta)` at `theta=pi/2`.

For sharpness let `T_m` be the Chebyshev polynomial and put

```math
f_{m,n}(x)=T_m\!\left(\frac{x_1+\cdots+x_n}{n}\right).
```

After Walsh reduction its degree is at most `m`, and its cube supremum is
exactly one. Permutation symmetry makes each first coefficient equal to
`a_{m,n}=E[f_{m,n} x_1]`. If `S` is the sum of the last `n-1` signs,

```math
n a_{m,n}=\mathbb E\frac{
 T_m((S+1)/n)-T_m((S-1)/n)}{2/n}.
```

For fixed `m`, `S/n` converges in probability to zero. The difference
quotient therefore tends in probability to `T_m'(0)`, and it is bounded
by the fixed finite number `||T_m'||_[-1,1]`. Bounded convergence gives

```math
L_1(f_{m,n})=n|a_{m,n}|\longrightarrow|T_m'(0)|=m.
```

Even degree bounds inherit the lower bound `m-1` by taking `T_(m-1)`.
The lower examples are real, hence also belong to the complex class.

For these same polynomials the first-level norm actually requested by BH
is

```math
\|\widehat f_{m,n}^{=1}\|_{q_m}
 =L_1(f_{m,n})\,n^{-(m-1)/(2m)}\longrightarrow0
 \quad(m>1\text{ fixed},\ n\to\infty).
```

Thus this is an actual auxiliary-functional obstruction, not a lower
bound tending to infinity for `B_m`.

## 3. Removing finitely many initial levels does not repair the functional

**Proved here.** For every fixed positive integer `r`, there are real
polynomials of degree at most `m_j -> infinity`, bounded by one, for which

```math
L_r(F_j)\sim \frac{e^{-r/2}}{r}m_j.
```

Take odd `k -> infinity`, set `n=k^2`, and on one block define

```math
h_k(x)=T_k\!\left(\frac{S_n}{n}\right),\qquad
S_n=\sum_{i=1}^n x_i,\qquad
a_k=\widehat h_k(\{1\}).
```

For odd `k`,

```math
T_k(u)=(-1)^{(k-1)/2}\sin(k\arcsin u).
```

Let `Z_n=S_n/sqrt(n)=S_n/k`. On every fixed compact set of `z`,
`k arcsin(z/k) -> z` uniformly. The central limit theorem and tightness
give convergence in distribution of `(Z_n, sin(k arcsin(Z_n/k)))`
to `(G,sin G)`. The absolute value of their product is at most `|Z_n|`,
whose second moment is exactly one. Therefore this product is uniformly
integrable, and

```math
\sqrt n\,|a_k|
 =\left|\mathbb E[Z_n h_k]\right|
 \longrightarrow \mathbb E[G\sin G]=e^{-1/2}.
```

The final Gaussian identity follows by differentiating
`E cos(tG)=exp(-t^2/2)` at `t=1`.

Take `r` disjoint blocks and put `F_k=product_(j=1)^r h_k(x^(j))`.
Its degree is at most `m=rk`, and its supremum is exactly one. Since
each factor is odd under reversal of all its own coordinates, it has no
constant term and only odd levels. The total level `r` of `F_k` therefore
consists exactly of choosing one coordinate from each block. It has
`n^r` coefficients, each equal to `a_k^r`. Consequently

```math
L_r(F_k)=n^{r/q_r}|a_k|^r
 =\sqrt n\,(\sqrt n\,|a_k|)^r
 \sim k e^{-r/2}.
```

For comparison, at the exponent the BH problem actually requests,

```math
\|\widehat F_k^{=r}\|_{q_{rk}}
 =n^{r/(2rk)}(\sqrt n\,|a_k|)^r
 =k^{1/k}(\sqrt n\,|a_k|)^r
 \longrightarrow e^{-r/2}.
```

**Consequence.** Let `(w_r)` be any fixed nonnegative sequence with
`w_r>0` for at least one fixed `r`. Every functional satisfying
`Phi(f)>=w_r L_r(f)` has optimal degree-`m` dimension-free constant
at least `c_r w_r m` along an infinite sequence. In particular, this
applies to every nontrivial fixed weighted square function
`(sum w_r^2 L_r^2)^(1/2)`, including ones omitting any finite collection
of initial levels. Thus replacing the paper's weights cannot by itself
make its auxiliary norm uniformly bounded.

This statement permits zero initial weights. It does not address
degree-dependent weights, level-dependent coefficient exponents depending
on the ambient degree, or a proof that retains actual coefficient
dispersion instead of replacing it by `L_r`.

### 3.1 The entire requested coefficient norm stays bounded in these examples

The preceding comparison is not only about the offending level. In
fact, for the same `r`-block examples,

```math
\|\widehat F_k\|_{q_{rk}}
 \longrightarrow\left(\frac{1-e^{-2}}2\right)^{r/2}.
```

Here is a complete justification. For odd `k`, the Taylor coefficients
of `T_k` satisfy

```math
|[t^j]T_k(t)|\le\frac{k^j}{j!}.
```

Only odd `j` occur. The exact coefficient is, up to its sign,
`k product_(a=1)^((j-1)/2) (k^2-(2a-1)^2)/j!`, which proves the bound.
Multiplication in the Walsh basis has nonnegative structure constants,
so the absolute Walsh coefficients of `h_k` are dominated coefficientwise
by those of `exp(S_n/k)`. Consequently, on level `j`,

```math
|a_{k,j}|\le\cosh(1/k)^{k^2}\tanh(1/k)^j
 \le e^{1/2}k^{-j}.
```

To identify their fixed-level limit, multiply `h_k` by its harmless
global sign so that `h'_k=sin(k arcsin(S_n/n))`, and consider

```math
G_k(z)=\mathbb E\left[h'_k(x)
             \prod_{i=1}^{k^2}\left(1+\frac{z x_i}{k}\right)\right].
```

For each fixed real `z` and large `k`, this is the expectation of
`h'_k` under independent signs of mean `z/k`. The triangular-array CLT
gives `S_n/k -> G+z`, so `G_k(z)->e^(-1/2)sin z`. On the complex
plane Cauchy--Schwarz supplies the locally uniform bound

```math
|G_k(z)|\le\left(1+|z|^2/k^2\right)^{k^2/2}
\le e^{|z|^2/2}.
```

Vitali's convergence theorem, or normal families plus uniqueness of
analytic continuation, therefore gives uniform convergence on compact
subsets of the complex plane. Comparing Taylor coefficients yields

```math
k^j|a_{k,j}|\longrightarrow
\begin{cases}e^{-1/2},&j\text{ odd},\\0,&j\text{ even}.\end{cases}
```

Fix `r` and put `q=q_(rk)`. For all `k` and `j`,

```math
\binom{k^2}{j}|a_{k,j}|^q
 \le e\,\frac{(k^{2-q})^j}{j!}
 \le e\,\frac{2^j}{j!}.
```

The last inequality follows from `k^(2/(rk+1))<=2`; it suffices to
use the elementary integer inequality `k^2<=2^(k+1)`. The right side
is summable, uniformly in `k`. For each fixed odd `j`, the left side
tends to `e^(-1)/j!`; even levels vanish. Dominated convergence gives

```math
\|\widehat h_k\|_q^q\longrightarrow e^{-1}\sinh1
 =\frac{1-e^{-2}}2.
```

The norm of the disjoint product is the product of coefficient norms,
giving the assertion. In particular this is an explicit separation
between linear growth of the weighted levelwise norm and a uniformly
bounded, nonvanishing FULL BH coefficient norm, not an unresolved
candidate counterexample to BH boundedness.

## 4. A different, strictly scoped damping barrier

Let `d,e>=2`, `r=d+e`, and `theta=d/r`. In the precise row/column proof,
the removal of the Boolean hypercontractive weights costs

```math
D(d,e)=\left(\frac{d+1}{d-1}\frac{e+1}{e-1}\right)^{de/(2r)}.
```

The elementary expansion

```math
\log\frac{t+1}{t-1}=2\sum_{j\ge0}\frac{1}{(2j+1)t^{2j+1}}>\frac2t
```

shows that `D(d,e)>e` for every finite pair. If both degrees tend to
infinity, the cost tends to `e`. With weights `w_r=r^(-s)`, the complete
local factor before capture and block-count losses is

```math
K_s(d,e)=D(d,e)\exp\{-s H(\theta)\},\qquad
H(\theta)=-\theta\log\theta-(1-\theta)\log(1-\theta).
```

Since `H<=log 2`, for `s>=0` one has

```math
K_s(d,e)>e\,2^{-s}.
```

Thus for `s<=1/log 2` the local factor is strictly greater than one
for every retained pair. Capture and level-count factors are at least
one and cannot help.

**Exact scope.** Suppose that, after coefficient interpolation and
inserting the hypercontractive radii, the only information retained is
the two independent energy budgets `sum U_i^2<=M^2` and
`sum V_i^2<=M^2`. The relaxed scalar problem admits a single index with
`U_i=V_i=M`. Therefore no strict contraction follows from these data
when its local coefficient is at least one. This proves a barrier for
this relaxation of the bootstrap. It is not a claim that actual Walsh
arrays can simultaneously saturate every discarded inequality.

For arbitrary positive fixed weights, a strict pairwise contraction
`K<=c<1` on the diagonal would require

```math
D(k,k)\frac{w_{2k}}{w_k}\le c,
\qquad\text{hence}\qquad
w_{2^j k}< (c/e)^j w_k.
```

So this scalar closure forces genuine power decay of the weights along
dyadic scales. In particular no weights bounded below can yield a
strict uniform contraction through this relaxation.

Using more colors in the same one-block-`q_d`, all-other-blocks-`2`
interpolation does not improve the obstruction. For `ell` blocks, sizes
`d_j` and `theta_j=d_j/r`, the damping logarithm is strictly larger than
`ell-1`, while the entropy is at most `log ell`. Its corresponding
weight threshold is therefore `(ell-1)/log ell`, minimized at `ell=2`.
This is a statement about that specified interpolation mechanism, not
about every possible multiblock inequality.

## 5. A concrete escape direction, not yet a theorem

The elementary identity

```math
\frac1{q_m}=\frac rm\frac1{q_r}
              +\left(1-\frac rm\right)\frac12
```

gives the stronger coefficient estimate

```math
\|\widehat f^{=r}\|_{q_m}
\le L_r(f)^{r/m}\|\widehat f^{=r}\|_2^{1-r/m}.
```

This retains the dispersion information lost by the paper's first
monotonicity step. It removes the fake first-level growth of the
Chebyshev examples. It does not by itself control levels with `r`
comparable to `m`, and so is not an improvement to the boundedness
frontier. A genuine new argument must couple this exponent interpolation
to restriction stability or directly control the high-degree coefficient
entropy. Merely renaming its right-hand side as a new square function
does not supply that argument.

## 6. Verification and preserved evidence

Canonical replay:

```bash
.venv/bin/python computations/bh_boundedness_barrier_2026_09_18.py \
  --output tmp/bh_boundedness_2026_09_18/barrier/replay.json
```

The script uses exact rational arithmetic for Chebyshev/Walsh recurrence
versus full cube enumeration. It also records floating-point CLT and
damping diagnostics, explicitly distinguished from proofs. Primary HTML
downloads and the Slote--Volberg PDF are preserved under the same scratch
directory. `pdftotext` was unavailable; the complete arXiv HTML was read
instead. The temporary folder must be included in the root campaign's
preservation checkpoint.

Independent verification requested: the parity argument selecting exactly
the level-r tensor coefficients; uniform integrability in the varying
Chebyshev-degree central-limit limit; the distinction between an abstract
scalar-budget barrier and an obstruction for actual coefficient arrays.

## 7. Exact bounded subclasses and a Galois extension

### 7.1 Boolean values: a known sharp theorem

For `f:{-1,1}^n -> {-1,1}` of degree at most `m`, all Fourier
coefficients lie in `2^(1-m) Z`. Indeed, `(1-f)/2` is integer-valued;
its multilinear expansion in `y_i=(1-x_i)/2` has integer Möbius
coefficients and degree at most `m`. Expanding the `y` monomials gives
the asserted granularity. Parseval and `q_m<2` now give

```math
\sum_S|\widehat f(S)|^{q_m}
\le (2^{1-m})^{q_m-2}\sum_S|\widehat f(S)|^2,
\qquad
\|\widehat f\|_{q_m}\le2^{1-1/m}.
```

The `(m-1)`-address-bit function selecting one of `2^(m-1)` separate
data bits has exactly `4^(m-1)` coefficients of magnitude `2^(1-m)`;
it attains equality. Thus `B_m>=2^(1-1/m)` for the problem's full
class, but Boolean-valued functions can never witness divergence.

This is an independently reconstructed **published result**, not new:
[*A cb-Bohnenblust--Hille inequality with constant one and its applications
in learning theory*](https://doi.org/10.1007/s00208-025-03142-5), Section 4.

### 7.2 Arbitrary finite root-of-unity alphabets

**Proved here; root independently reconstructed the proof; external
novelty is not claimed.** If `f` has degree at most `m` and every value
of `f` is a root of unity, then

```math
|\operatorname{supp}\widehat f|\le4^m,
\qquad \|\widehat f\|_{q_m}\le2.
```

Since the cube is finite, all values belong to `mu_K` for one finite
`K`, allowed to grow arbitrarily with `m` and `n`. Let `L=Q(zeta_K)`,
let `D=[L:Q]`, and let `sigma` run over all `D` complex embeddings.
Write `a_S=fhat(S)`. The same integer Möbius formula, now over the ring
of algebraic integers, shows

```math
2^m a_S\in\mathbb Z[\zeta_K].
```

Every `sigma(f)` is still root-of-unity-valued, has degree at most `m`,
and obeys `sum_S |sigma(a_S)|^2=1`. The nonzero Fourier support is
the same for every embedding. For each nonzero coefficient, its nonzero
integer field norm gives

```math
\prod_\sigma|\sigma(a_S)|^2\ge2^{-2mD}.
```

Generalized Hölder then implies

```math
|\operatorname{supp}\widehat f|\,2^{-2m}
\le\sum_{a_S\ne0}\prod_\sigma|\sigma(a_S)|^{2/D}
\le\prod_\sigma\left(\sum_S|\sigma(a_S)|^2\right)^{1/D}=1.
```

The coefficient norm comparison `ell_q <= support^(1/q-1/2) ell_2`
gives the constant 2. Both the alphabet size and its cyclotomic field
degree cancel from this argument.

More generally, for any algebraic-integer-valued degree-`m` function
over a number field of degree `D`, the same proof gives

```math
|\operatorname{supp}\widehat f|
\le4^m\left(\prod_\sigma\|\sigma(f)\|_2\right)^{2/D}
\le4^m\left(\prod_\sigma\|\sigma(f)\|_\infty\right)^{2/D}.
```

This extension was suggested by the root agent. It is only useful when
the conjugate sizes are controlled. Arbitrary unimodular algebraic
values need not be algebraic integers, and arbitrary algebraic values
need not have bounded conjugates.

## 8. Torsion-valued low-degree functions are not dense

The most immediate extension of Section 7 is false: arbitrary
unimodular-valued low-degree functions cannot necessarily be approximated
by root-of-unity-valued functions while preserving their degree.

### 8.1 Explicit degree-two counterexample

List eight numbers, one for each vertex of the three-dimensional cube:

```math
z=\left(1,1,
 -\frac15\pm\frac{2i\sqrt6}{5},
 -\frac13\pm\frac{2i\sqrt2}{3},
 -\frac7{15}\pm\frac{4i\sqrt{11}}{15}\right).
```

Each has modulus one and their sum is zero. Under any fixed enumeration
of the cube's vertices set `f(x)=x_1 x_2 x_3 z_x`. The coefficient at
`{1,2,3}` is `sum z/8=0`, hence `deg f<=2`, while `|f|=1` everywhere.

No nonempty proper subsum of `z` vanishes. The rational independence of
`sqrt(6), sqrt(2), sqrt(11)` forces a vanishing subsum to include either
both or neither member of each conjugate pair. This independence follows,
for example, by independent sign automorphisms in
`Q(sqrt(2),sqrt(3),sqrt(11))`. Its remaining real equation is

```math
15h-6b-10c-14d=0,\qquad
h\in\{0,1,2\},\quad b,c,d\in\{0,1\}.
```

For `h=0` only the empty sum is possible; `h=1` is excluded by parity;
and `h=2` forces `b=c=d=1`, the full sum.

### 8.2 The needed roots-of-unity theorem, with proof

We use the elementary special case of
[Mann's theorem](https://doi.org/10.1112/S0025579300005210): if an
irreducible sum of `k` roots of unity vanishes, then every ratio of two
terms has order dividing the product of the primes at most `k`.
Here irreducible means that no nonempty proper subsum vanishes.

For completeness, normalize one term to 1 and let `N` be the least
common multiple of all term orders. For a prime power `p^a || N`, put
`M=N/p^a` and write every term as `eta^e xi`, where `eta` is primitive
of order `p^a`, `0<=e<p^a`, and `xi` is an `M`-th root. The minimal
polynomial of `eta` over `Q(mu_M)` is

```math
\Phi_{p^a}(X)=\sum_{j=0}^{p-1}X^{j p^{a-1}}.
```

Indeed, coprime cyclotomic degrees multiply. Grouping the relation by
the exponent `e` gives a polynomial of degree less than `p^a` divisible
by this cyclotomic polynomial. Its coefficients in the positions
`t+j p^(a-1)` must therefore be equal for each fixed `t`. Consequently
the original terms with a fixed residue `e mod p^(a-1)` already sum to
zero. Irreducibility forces all terms into a single such residue; the
normalized term makes it residue zero. If `a>=2`, this contradicts
minimality of the `p^a` part of `N`. Thus `N` is squarefree.

Now `a=1`. If `p>k`, some exponent class is empty, so the equal grouped
coefficients are all zero. Each nonempty exponent class then itself
gives a vanishing subsum. Irreducibility forces only one class, and the
normalized term again makes it class zero, contradicting `p|N`.
Therefore every prime divisor of `N` is at most `k`, as asserted.

### 8.3 Failure of density, including auxiliary-variable approximation

Suppose root-of-unity-valued degree-two functions converged uniformly
to the `f` above. Multiplying their eight values by `x_1x_2x_3` gives
eight roots of unity summing to zero. The same proof permits an arbitrary
common unit-modulus phase multiplying each approximating function: divide
out that phase before applying Mann's theorem. Nearby tuples are irreducible,
because there are only finitely many proper subsums and every such
subsum of `z` is nonzero. Mann's theorem with `k=8` puts every ratio
in the finite set `mu_210`.

Their third-to-first ratios would converge to
`-1/5+2i sqrt(6)/5`, so this number would belong to `mu_210`. It is
not even a root of unity: a root of unity and its inverse are algebraic
integers, whereas their sum here is the noninteger rational `-2/5`.
This contradiction proves non-density.

The obstruction also permits approximating polynomials with extra
variables: fixing those extra variables at any vertex reduces to the
same three-variable contradiction without increasing degree or error.
This fixed example alone does not exclude bounded-mass convex combinations,
degree inflation, or a direct proof for all unimodular functions.

It also excludes a norm-one closed-convex-hull reduction to the torsion
class. If convex averages converged to `f`, the identity
`E_g ||g-f||_2^2=2(1-Re<E_g g,f>)` would tend to zero. Since the cube
has eight points, some torsion component would then approach `f`
uniformly, which is impossible. The quantitative gap below gives
`Re<g,f><=1-1/2560000` for each such `g`, so a scaled convex
representation requires factor at least `(1-1/2560000)^(-1)`.
Because arbitrary common phases were permitted, this also excludes a
norm-one closed absolutely convex hull reduction. An absolute factor
larger than one is not excluded by this fixed example.

The companion script strengthens the qualitative gap with exact
rational interval checks:

```math
\min_{\varnothing\ne J\subsetneq[8]}
 \left|\sum_{j\in J}z_j\right|\ge\frac1{15},\qquad
\operatorname{dist}\left(-\frac15+\frac{2i\sqrt6}{5},\mu_{210}\right)
 >\frac1{200}.
```

It follows that the uniform distance from `f` to every torsion-valued
degree-two function is at least `1/400`. At smaller error, all proper
subsum errors are below `7/400<1/15`, while the ratio error is below
`2/400=1/200`, contradicting the two inequalities. Square-root bounds
are checked by squaring rational endpoints; the trigonometric intervals
use the exact Machin formula and alternating arctangent series, followed
by Taylor bounds for sine and cosine. Numerical minima are supplementary,
not the certificates.

## 9. Independent audit of the counterexample track's finite lower bound

The localization/counterexample agent found an exact degree-two
sixth-root-valued five-variable polynomial `f`, with 16 coefficients
of magnitude `1/4` and supremum one. Its full certificate and discovery
belong to that agent's canonical artifact. Independent exact replay here
confirms every pointwise value in the ring `Z[zeta]`, `zeta^2=zeta-1`.
It attains the torsion-class upper bound 2 from Section 7; tensor powers
attain ratio 2 at every even degree.

More importantly, that agent's degree-transfer perturbation proves the
actual unrestricted finite bound

```math
\boxed{B_4>2+2^{-18}.}
```

This falsifies the proposed universal value 2, not uniform boundedness by
some larger absolute constant. I checked the following entire proof
independently, including exact coefficient and pointwise arithmetic.

Put `a=conj(zeta)+zeta x_3+x_4` and `b=f^2 conj(a)`. Exact Walsh
reduction gives `deg a=1`, `deg b=3`. On two disjoint blocks define

```math
F(x,y)=f(x)f(y),\qquad
h(x,y)=a(x)b(y)-b(x)a(y).
```

Both have degree at most 4. The ratio `a/f` takes values in
`{0} union 2 mu_6`, and `b/f=conj(a/f)`. Consequently

```math
\frac hF\in\{0,+4i\sqrt3,-4i\sqrt3\},\qquad
\|F+\varepsilon h\|_\infty=\sqrt{1+48\varepsilon^2}
\quad(\varepsilon\in\mathbb R).
```

All three values are attained; the exact 1024-point counts are
`640,192,192`. The old Fourier support of `F` has 256 coefficients,
all of modulus `1/16`. The perturbation has 36 coefficients outside
that support: 24 of magnitude `sqrt(3)/4` and 12 of magnitude `3/4`.
Their squared mass is `45/4`. One explicit outside coefficient is
`hhat(0,21)=3 conj(zeta)/4`, in five-bit subset-mask notation.

Pointwise tangency gives `Re E[conj(F)h]=0`. Since the old coefficient
moduli are equal, the derivative at zero of their total `q`-mass is
zero for `q=8/5`. Convexity therefore bounds that old mass below by
its unperturbed value `A=2^q`, for every real perturbation parameter.
All outside coefficients have magnitude below one, so their `q`-mass
is at least their squared mass `45/4`.

Take `epsilon=2^(-10)`. Then `epsilon^q=2^(-16)`, and concavity gives

```math
\|F+\varepsilon h\|_\infty^q
=(1+48\varepsilon^2)^{4/5}\le1+\frac3{81920}.
```

For the actual BH ratio `R`,

```math
R^q\ge\frac{A+45/262144}{1+3/81920},\qquad
R^q-A\ge\frac{33}{16\cdot81923}>2^{-16},
```

where only `A<=4` was used. Also `A<31/10`, because
`31^5>256*10^5`; hence `A+2^(-16)<4`. On that interval the derivative
of `t^(5/8)` exceeds `5/16>1/4`. Applying this derivative bound to
`A+2^(-16)` proves `R>2+2^(-18)`.

Independent replay:

```bash
.venv/bin/python computations/bh_boundedness_barrier_2026_09_18_tangent.py \
 --output tmp/bh_boundedness_2026_09_18/barrier/tangent_replay.json
```

Its exact algebra and rational checks pass. Optional decimal ratio
diagnostics are not used in the proof.

## 10. Independent checks supplied to the positive track

Sections 2 and 3 of `bh_boundedness_positive_2026_09_18.md` were independently
reconstructed. Both pass.

The positive track proves a uniform bound `sqrt(6e)` for complex
permutation-symmetric polynomials. The key spin-rotation calculation is
correct: the binomial-weighted sum of the squared Dicke commutator columns
equals `(n+1) I(f)/2`, so operator Bernstein gives
`I(f)<=2m^2 ||f||_infinity^2/(n+1)`. The binomial coefficient envelope
and final norm conversion have the stated factors. The invariant
symmetric tensor subspace is genuinely needed.

The same artifact's Boolean-atom approximation obstruction also passes:
the correlated sine-rounding law has the asserted expectation and
variance at most `7/N`, every fixed-coordinate marginal tends uniformly
to independent signs, and the degree-dependent junta bound is valid.
Its polynomial amplification defeats every dimension-independent
supremum-norm approximation contraction by bounded-mass, finite-degree
Boolean atoms. The order of its limits and the finite error constants
were checked. This route obstruction is not a divergent BH family.

## 11. Two amplification guards for the finite improvement

The five-bit degree-two root-valued seed gives a useful sharp obstruction
to a tempting amplification. These statements concern exact supports,
not thresholded numerical Fourier transforms.

First, every polynomial on `n` bits of degree budget `m` obeys

```math
\|\widehat P\|_{q_m}
\le 2^{n/(2m)}\|\widehat P\|_2
\le 2^{n/(2m)}\|P\|_\infty.
```

Consequently, any balanced degree-transfer construction using `k`
five-bit blocks at degree budget `2k` has ratio at most `2^(5/4)`.
Continuous phases and arbitrarily many balanced replacements of degrees
`2+2` by `1+3` do not change this ambient-dimension guard. It is a bound
at the stated degree budget; a genuine degree collapse would require a
separate calculation.

Second, the old torsion support forbids such a collapse under injective
parity recoding for the explicit small perturbation. Here is the general
statement. Suppose `u` is root-of-unity-valued and has `s` nonzero Fourier
coefficients, and

```math
\operatorname{supp}\widehat u\subseteq
\operatorname{supp}\widehat P.
```

Apply an injective linear or affine map over `F_2` to all Fourier indices.
The transformed `u` remains root-of-unity-valued: a linear map is dual to
a surjective homomorphism of cubes, and the affine translation multiplies
by a character. Its support still has size `s`. If the transformed `P`
has degree `D`, the transformed `u` has degree at most `D`; Section 7
therefore forces

```math
D\ge\frac{\log s}{\log4}.
```

For `P=F+2^(-10)h` from Section 9, no old Fourier coefficient vanishes:
the old modulus is `1/16`, while
`2^(-10)|hhat(S)|<=2^(-10)||h||_infinity=4sqrt(3)/1024<1/16`.
For `r` disjoint copies, the support of `P^[tensor r]` thus contains
the support of the torsion function `F^[tensor r]`, of size
`256^r=4^(4r)`. Any such injective affine recoding must have degree at
least `4r`. Thus tensoring followed by generic parity recoding cannot
lower the degree budget that obstructs amplification. The statement
does not cover coefficient cancellations introduced by another operation,
or noninjective substitutions in which Fourier coefficients collide.

The counterexample track independently obtained a broader fixed-seed
recoding obstruction using spectral parity biases. The torsion-skeleton
argument above additionally certifies the exact no-compression threshold
`4r` for this particular perturbation family.

## 12. Simpler uniform bound for the Chebyshev obstruction family

The counterexample track supplied a useful independent safeguard for
Section 3. If `g=sum a_i x_i` with `sum |a_i|<=1`, the Fourier Wiener norm
is submultiplicative, so

```math
\|\widehat{T_m(g)}\|_1
\le \sum_j |[t^j]T_m(t)|=|T_m(i)|
\le(1+\sqrt2)^m.
```

The middle identity follows from the alternating signs of the nonzero
Chebyshev coefficients. For a function with actual supremum one,
interpolation with its Fourier `ell_2` norm at most one gives

```math
\|\widehat{T_m(g)}\|_{q_m}
\le\|\widehat{T_m(g)}\|_1^{1/m}
       \|\widehat{T_m(g)}\|_2^{1-1/m}
\le1+\sqrt2.
```

For the equal-weight real means used above, the actual cube supremum is
exactly one, since the all-one vertex attains `T_m(1)=1`. The same Wiener
bound and interpolation apply to disjoint products, with total degree
budget the sum of the degrees. This simple estimate is independent of
the more precise limiting coefficient calculation in Section 3.1.
It must not be applied after freely renormalizing an arbitrary complex
`g`: the interval bound of `T_m` alone need not control its complex cube
values or their actual supremum.

## 13. Degree inflation does not repair a uniform torsion-atom dilation

Combining Section 7 with the independently audited positive-track
construction yields a stronger obstruction than the fixed example in
Section 8. The complete separating-distribution proof is canonical in
`bh_boundedness_positive_2026_09_18.md`, Sections 2 and 2.7. Both agents
independently checked this extension.

For any degree-`D` root-of-unity-valued function, the union of its at most
`4^D` nonzero Fourier supports contains at most `D4^D` coordinates.
Thus every such function is a `D4^D`-junta, independently of its alphabet
order. Multiplication by an arbitrary complex phase does not change its
dependence. The positive track's estimate applies to bounded complex
juntas just as it applies to Boolean ones.

More explicitly, for every fixed `k` it constructs real functions
`P_(k,N)` of degree at most `4k`, bounded by one, and probability laws
`mu_N,nu_N` on their input cube such that

```math
\frac12\left|\mathbb E_{\mu_N}P_{k,N}
                  -\mathbb E_{\nu_N}P_{k,N}\right|
\longrightarrow\gamma_k
=1-(1-e^{-2})^k,
```

while for every bounded complex `K`-junta `a`,

```math
\left|\mathbb E_{\mu_N}a-\mathbb E_{\nu_N}a\right|
\le\frac{K}{\sqrt{2N}}.
```

Hence any representation `G=sum_j c_j a_j` with
`sum_j |c_j|<=C`, each `a_j` a root-valued degree-`D` atom, satisfies

```math
\liminf_{N\to\infty}\|P_{k,N}-G\|_\infty\ge\gamma_k.
```

The statement allows `G` and all its atoms to change with `N`; the
bound on their expectation difference is uniform. More quantitatively,
the positive track's finite estimate gives the lower bound

```math
\gamma_k-2k\sqrt{7/N}-3k/N
            -\frac{CD4^D}{2\sqrt{2N}}.
```

Suppose a universal approximation contraction with residual `delta<1`
were claimed, using any finite dimension-independent allowances
`D=D(m)` and `C=C(m)`. Choose `k` with `gamma_k>delta`, set `m=4k`,
and then let `N` tend to infinity. This contradicts that contraction.
In particular, arbitrary finite alphabet growth, constant-factor degree
inflation, and a bounded total variation of complex coefficients do not
extend the Galois subclass theorem to all bounded low-degree functions
by this route. This remains a method barrier, not a lower bound
diverging with `m` for the actual BH ratio.

## 14. Independent audit of a sustained asymptotic lower bound

The counterexample track subsequently proved the stronger statement

```math
\liminf_{m\to\infty} B_m>2+2^{-207}.
```

This track independently checked its complete argument and every finite
constant. The canonical proof and construction are in that track's
artifact, Section 8; this paragraph records the verification dependency,
not independent priority. It still does not prove divergence.

For the exact `F,h` of Section 9, set `epsilon=2^(-104)`,
`R=|h|^2/48 in {0,1}`, `lambda=(1/2)log(1+48epsilon^2)`, and
`G=(F+epsilon h)exp(-lambda R)`. Then `|G|=1`. The verified Fourier
entropy bound is

```math
H(G)\ge8\log2+\varepsilon^2
       \left[\frac{45}{4}\log(1/\varepsilon)-256\right].
```

The exact old-support flatness and outside squared mass `45/4` are
essential. Writing `G=F+epsilon h+w`, one has `||w||_2<=16epsilon^2`
and `||G-F||_2<=(9/2)epsilon`. On the old support, Taylor's inequality
for `-p log p` has a total negative correction below `256epsilon^2`:
the two conservative bounds are `7*(81/4)` and
`256*(153/256)^2`, whose sum is below 256. On the 36 new coefficients,
the squared modulus is at least
`(1-128epsilon)epsilon^2 |hhat|^2` and at most `epsilon^2`.
These checks establish the displayed finite entropy bound, without an
unspecified small-parameter remainder.

For `k` blocks, truncate the Taylor series of
`exp(-lambda sum_i R_i)` at degree `L_k-1`, where
`L_k=ceil(4lambda k)`, and multiply by `(F+epsilon h)^[tensor k]`.
The exact identities `R^2=R`, `hR=h`, and `deg(F|h|^2)=8` imply that
every selected block `(F+epsilon h)R` has degree at most 8, whereas an
unselected block has degree at most 4. Thus the resulting `Q_k` has
degree at most `M_k=4k+4(L_k-1)` and approximates `G^[tensor k]` uniformly with error
at most `exp(-lambda k/3)`. Importantly, every error function on `10k`
bits satisfies

```math
\|\widehat E\|_{q_m}\le2^{5k/m}\|E\|_\infty,
```

so the approximation transfers to the final changing BH exponent.
The bounded gaps of `M_k` allow every sufficiently large degree budget
`m`, with `m/k -> d=4+16lambda`. The tensor entropy formula then gives
`liminf B_m>=exp(H(G)/(2d))`. Since `lambda<=24epsilon^2`, the entropy
margin over `2d log2` is at least
`(402 log2-256)epsilon^2>=12epsilon^2`. As `d<5`, this implies
`H(G)/(2d)>log2+epsilon^2`, proving the asserted strict lower bound.
The added normalization degree has been paid in full. The earlier,
looser degree budget `4k+8(L_k-1)` already gave the weaker valid constant
`2+2^(-1016)`; the exact block identity improves that bookkeeping.

The independent tangent replay now also verifies the radial identities,
`deg(F|h|^2)=8`, its exact coefficient `3/4` at mask 503, and its squared
`L2` norm 864. These are exact Eisenstein-integer calculations.

## 15. Independent Hadamard-chain audit

The root agent proposed the following obstruction to a different
counterexample mechanism. It is correct for disjoint variable blocks.
Let `U_j` be real orthogonal or complex unitary `N`-by-`N` matrices
whose entries all have modulus `N^(-1/2)`, and let `d>=2`.

For the open chain

```math
T(x^{(1)},\ldots,x^{(d)})
=\frac1N(x^{(1)})^{\mathsf T}U_1D(x^{(2)})U_2\cdots
              D(x^{(d-1)})U_{d-1}x^{(d)},
```

the `N^d` coefficients all have modulus `N^(-(d+1)/2)`, hence their
`q_d` norm is one. Fix the tail signs and write the vector following
`D(x^(2))` as `v`; unitarity gives `||v||_2=sqrt(N)`. Every coordinate
of `U_1D(x^(2))v`, with random independent second-block signs, has
coefficient `ell_2` norm one. Real Khintchine gives expected absolute
value at least `1/sqrt(2)`. Choose the first-block signs to match the
coordinates. Thus the real cube supremum is at least `1/sqrt(2)`.
For `d=2`, randomize the final block directly, with the same calculation.

For complex coefficients, the Hilbert-space form of Khintchine retains
`1/sqrt(2)`, while the elementary phase projection inequality

```math
\max_{\epsilon_i\in\{-1,1\}}
 \left|\sum_i\epsilon_i z_i\right|
\ge\frac2\pi\sum_i|z_i|
```

follows by averaging `sum_i |Re(e^(-it)z_i)|` in `t`. Consequently the
complex open-chain supremum is at least `sqrt(2)/pi`. Its BH ratio is
therefore at most `sqrt(2)` in the real case or `pi/sqrt(2)` in the
complex case, uniformly in length and matrix size.

For the closed cycle

```math
C=\operatorname{tr}(D(x^{(1)})U_1\cdots D(x^{(d)})U_d),
```

the coefficient `q_d` norm is `sqrt(N)`. Fix blocks 3 through `d`, and
write `V=U_2D(x^(3))U_3...D(x^(d))U_d`, a unitary matrix. The remaining
bilinear coefficients are `U_(1,ij)V_(ji)`, whose squared row sums
are all `1/N`. The preceding randomization and sign choice give
supremum at least `sqrt(N/2)` in the real case and `sqrt(2N)/pi` in
the complex case. The same uniform ratio bounds follow. Identifying
variable blocks can create coefficient collisions and degree changes;
that operation is not covered by this proof.

## 16. A fixed-number-of-symmetry-blocks positive theorem

This extends the positive track's spin-rotation mechanism. The positive
track and this track noticed the extension independently; the complete
calculation below is supplied here for independent verification.

Suppose the variables are partitioned into `ell` disjoint blocks of
sizes `n_1,...,n_ell`, and a complex polynomial `f` of degree at most
`m>=1` is invariant under all permutations within each block. Then

```math
\frac{\|\widehat f\|_{q_m}^2}{\|f\|_\infty^2}
\le 2e(1+2\ell)\binom{m+\ell}{\ell}^{1/m}
\le2e(1+2\ell)(\ell+1).
```

Thus any fixed number of symmetry blocks has a uniform BH constant,
independent of all block sizes and of degree. The single-block constant
in the positive artifact is sharper than this deliberately loose bound.

Let `I_j(f)=sum_S |S intersect block_j| |fhat(S)|^2`. On the tensor
product of the symmetric Dicke spaces, let `F` be the diagonal matrix
whose entries are the layer values `f(k_1,...,k_ell)`, and put
`J_x=sum_j J_(x,j)`. The matrix-valued trigonometric polynomial
`exp(itJ_x)F exp(-itJ_x)` has degree at most `m`: each original Walsh
monomial of degree `r` rotates to a product of `r` sine/cosine factors.
Its operator norm is constantly `||f||_infinity`. Scalar Bernstein
applied to every matrix element between unit vectors gives

```math
\|[J_x,F]\|\le m\|f\|_\infty.
```

In a Dicke column indexed by `k`, the transitions in different blocks
lead to distinct basis vectors. Thus their squared norms add without
cross terms. Weighting columns by the product binomial law, the exact
one-block calculation yields

```math
\frac12\sum_{j=1}^\ell(n_j+1)I_j(f)
\le m^2\|f\|_\infty^2.                 \tag{*}
```

For clarity, an edge `k_j -> k_j+1` has squared transition coefficient
`(k_j+1)(n_j-k_j)/4`. Its two binomial column weights combine to
`(n_j+1) binom(n_j,k_j)2^(-n_j)(n_j-k_j)/4`, exactly `(n_j+1)/2`
times that edge's contribution to `I_j`. This verifies (*) without a
factor of `ell` in the commutator bound.

For a multilevel `r=(r_1,...,r_ell)`, write
`M_r=product_j binom(n_j,r_j)` and
`v_r=M_r |a_r|^2`, where `a_r` is its common Fourier coefficient.
Set `s=sum_j r_j/m<=1`, `theta_j=r_j/m`, and
`b_j=sqrt(n_j r_j)/m`, omitting zero `r_j`. Weighted arithmetic-geometric
mean, with the extra weight `1-s` and extra term also `1-s`, gives

```math
\prod_j(n_j/r_j)^{r_j/m}
\le(1-s+\sum_j b_j)^2
\le2\left(1+\frac{\ell}{m^2}\sum_j n_jr_j\right).
```

The binomial bound `binom(n,r)<=(en/r)^r` therefore implies

```math
M_r^{1/m}\le2e\left(1+\frac{\ell}{m^2}
                                      \sum_jn_jr_j\right).
```

At exponent `q_m`, that multilevel's squared coefficient norm is
`M_r^(1/m)v_r`. There are at most `K=binom(m+ell,ell)` multilevels.
Summing by counting-measure Hölder, then using Parseval and (*), gives

```math
\|\widehat f\|_{q_m}^2
\le K^{1/m}\sum_r M_r^{1/m}v_r
\le2eK^{1/m}\left(\|f\|_2^2+
             \frac{\ell}{m^2}\sum_jn_j I_j(f)\right)
\le2eK^{1/m}(1+2\ell)\|f\|_\infty^2.
```

Finally `K<=(ell+1)^m`, because every weak composition into `ell+1`
parts is realized by at least one word of length `m` on an `ell+1`
letter alphabet. No uniform bound independent of the number of symmetry
blocks is obtained by this argument.

## 17. Further structural routes examined, without a claimed resolution

### 17.1 General circle-valued sparsity is not established

The arithmetic support bound `4^m` suggests asking whether every
unimodular-valued degree-`m` function, without any torsion assumption,
has at most `C^m` nonzero Fourier coefficients. No proof or counterexample
was obtained here. The exact non-density result in Section 8 prevents
deducing this by degree-preserving phase quantization. The counterexample
track's least-squares searches at degree two on six and seven variables
found only genuine numerical solutions with support at most 16; other
larger supports had substantial modulus residual. That is exploratory
evidence, not an upper bound or a classification.

For degree one the elementary classification is easy: in
`f=c+sum_i a_i x_i`, constancy of `|f|^2` requires all nonzero vectors
among `c,a_1,...,a_n` to be pairwise orthogonal in the real plane.
There are at most two, so the support size is at most two. This does
not generalize directly: higher-degree autocorrelation equations involve
sums of several products and allow genuine cancellation.

The closely related noncommutative junta question is not an available
shortcut. The primary paper by Montanaro--Osborne,
[arXiv:0810.2435v5](https://arxiv.org/html/0810.2435v5), Section 12,
formulates an exponential degree/junta relation as a conjecture, and
explains why the minimum-influence proof for Boolean outputs does not
apply to continuous-valued quantum objects. No claim about its present
general resolution is made from that historical source. A scalar
circle-valued function embeds as a matrix-valued involution with one
extra output qubit, but that observation alone proves no sparsity bound.

### 17.2 Low approximate degree and many variables are insufficient

A possible route to a BH lower bound is a bounded low-degree
approximation to a unimodular target whose Fourier entropy is much larger
than the approximation degree. Section 14 implements this with a small
but rigorous entropy surplus. Merely having many essential variables
does not provide that entropy.

For example, the first quantum-addressing construction in
Ambainis--de Wolf,
[arXiv:1206.0717v2](https://arxiv.org/html/1206.0717v2), Section 3.1,
addresses `t^t` data bits using `t^2` address bits, but all addresses
except a specified set of `t^t` words return the same default data bit.
Under uniform input, the probability of a nondefault address is at most
`eta=t^t/2^(t^2)`. If `f` is the resulting Boolean function and `x_1`
is its default dictator, then

```math
\|f-x_1\|_2^2\le4\eta.
```

Thus its Fourier mass outside the dictator coefficient is at most
`4eta`. A Fourier coefficient must involve exactly one data bit and
an arbitrary subset of the `t^2` address bits, so there are at most
`t^t 2^(t^2)` possible coefficients. The elementary entropy bound

```math
H(|\widehat f|^2)
\le h_2(4\eta)+4\eta\log(t^t2^{t^2})
\longrightarrow0
```

holds for large `t`, with `h_2` the binary entropy in natural units.
This particular dimension-saving construction therefore does not supply
the needed large Fourier entropy. This observation does not rule out
other approximation-theoretic mechanisms or prove anything about the
BH norms of arbitrary approximating polynomials themselves.

## 18. The proposed final-exponent bootstrap: two exact corrections

At approximately 05:24 UTC the root proposed a genuinely different
positive mechanism: retain the final exponent `q_m` throughout, bound
low Fourier levels by an `ell_1` projection estimate, and charge each
recursive damping step by its local degree divided by the original `m`.
The following checks reject that particular proposal. They do not reject
all possible final-exponent arguments or prove a global BH obstruction.

### 18.1 Fourier `ell_1` is not the low-level supremum projection norm

The proposed base estimate `||fhat^(=r)||_1<=(Cm)^r||f||_infinity`
is false for every `r=2`, already at degree budget 2. Let `U=H/sqrt(N)`
be a flat orthogonal Hadamard matrix and set

```math
f(x,y)=\frac1N x^{\mathsf T}Uy.
```

The operator norm gives `||f||_infinity<=1`, but its `N^2` coefficients
have modulus `N^(-3/2)`, so their `ell_1` norm is `sqrt(N)`. Every
coefficient is on level two. This also belongs to every larger degree
budget. A dimension-free estimate for `||f^(=r)||_infinity` cannot be
substituted for an estimate of its coefficient Wiener norm.

Correct interpolation through `q_r`, not through 1, reads

```math
\|\widehat f^{=r}\|_{q_m}
\le L_r(f)^{r/m}\|\widehat f^{=r}\|_2^{1-r/m}.
```

Consequently a projection bound of order `(Cm)^r`, combined with a
degree-`r` BH bound, costs roughly `(Cm)^(r^2/m)` in this calculation,
not `(Cm)^(r/m)`. Its elementary constant-cost low-level range is only
on the order of `sqrt(m/log m)`, not `m/log m`. The root acknowledged
and withdrew the incorrect base estimate.

### 18.2 Final-exponent mixed norms conserve the exponent deficit

For a coefficient matrix `A`, put

```math
R_p(A)=\left(\sum_i(\sum_j|A_{ij}|^2)^{p/2}\right)^{1/p},
\quad
C_p(A)=\left(\sum_j(\sum_i|A_{ij}|^2)^{p/2}\right)^{1/p}.
```

Mixed-norm interpolation with weights `theta,1-theta` can produce the
flat exponent `q_M` only when

```math
\frac1{q_M}=\frac\theta{p_A}+\frac{1-\theta}{2}
           =\frac\theta2+\frac{1-\theta}{p_B}.
```

Thus necessarily `p_A=q_(M theta)`, `p_B=q_(M(1-theta))`.
Using `q_M` in both mixed norms would instead give `q_(2M)` when
balanced. The all-one `N`-by-`N` matrix demonstrates the missing factor
directly: its flat `q_M` norm is `N^(2/q_M)`, whereas both mixed norms
at `q_M` are only `N^(1/q_M+1/2)`. The ratio is `N^(1/(2M))`.

Consider a bihomogeneous block polynomial of degrees `d,e`, with
`r=d+e`. The full-row hypercontractive argument gives exactly the local
factor

```math
D=\left(\frac{M\theta+1}{M\theta-1}\right)^{\theta e/2}
  \left(\frac{M(1-\theta)+1}{M(1-\theta)-1}\right)^{(1-\theta)d/2},
\qquad \log D>\frac rM,
```

provided both child budgets exceed one. The inequality follows from
`log((u+1)/(u-1))>2/u`. More explicitly, if `C_(r,M)` denotes the
degree-`r` coefficient bound at exponent `q_M`, that argument has the
form

```math
C_{r,M}\le D\,
 C_{d,M\theta}^{\theta}
 C_{e,M(1-\theta)}^{1-\theta},
```

apart from any additional capture/projection losses. Dimension-free
critical-exponent admissibility requires `M theta>=d` and
`M(1-theta)>=e`. In the critical case `M=r`, these inequalities force
`theta=d/r`, with child budgets exactly `d,e`. Thus the original `m`
does not remain the exponent parameter down the recursion: each child
remains at its own critical exponent. This is the precise point at
which the proposed diminishing damping charge fails.

For completeness the accumulated factor has an allocation-independent
tree lower bound. At a node `v` write its degree and exponent budget as
`r_v,M_v`; the children split both additively. Its weight in the root
product is `omega_v=M_v/M_root`. Hence its weighted logarithmic damping
is strictly greater than `r_v/M_root`. If every leaf degree is at most
`h`, put `p_leaf=r_leaf/r_root`. Their expected binary-tree depth is
`sum_internal r_v/r_root`, and the prefix-code entropy inequality gives

```math
\sum_{v\text{ internal}}r_v
\ge\frac{r_{\rm root}}{\log2}H(p)
\ge\frac{r_{\rm root}}{\log2}\log(r_{\rm root}/h).
```

Therefore this exact full-row closure incurs at least
`(r_root/h)^(r_root/(M_root log2))`. In the critical case, stopping at
`h=m/log m` still costs at least `(log m)^(1/log2)`, even with unbalanced
splits. This statement only tracks the specified proof factors; actual
Walsh arrays need not saturate them simultaneously. Additional structural
information that improves the row comparison is outside the obstruction.

## 19. Accessible-address information under a uniform cube

The root suggested making efficient quantum-address inputs typical via
biases, marker patterns, or minimum finding. The following exact
information bound controls what any such scheme can extract from a
uniform expanded Boolean cube. It does not control all Fourier entropy
of the resulting address function; the distinction is recorded below.

### 19.1 A quantum query information bound

Let `X` be uniform on a Boolean cube and let a quantum algorithm use
`T` input queries, followed by an arbitrary input-independent measurement
with classical outcome `Z`. Then, in natural logarithm units,

```math
I(X;Z)\le2T.
```

Each final amplitude is a complex Walsh polynomial of degree at most
`T`, since each query multiplies by an input-dependent matrix of degree
at most one. An input-independent measurement can be absorbed into a
linear transformation of amplitudes and a grouping of squared moduli.
For a degree-`T` amplitude `a`, differentiating Boolean hypercontractivity

```math
\|a\|_{2+s}\le(1+s)^{T/2}\|a\|_2
```

at `s=0` gives

```math
\operatorname{Ent}(|a|^2)
\le2T\,\mathbb E|a|^2,
\qquad \operatorname{Ent}(u)
=\mathbb E[u\log(u/\mathbb Eu)].
```

The entropy functional on nonnegative functions is subadditive:
`Ent(sum_j u_j)<=sum_j Ent(u_j)`. This follows from convexity of
relative entropy after normalizing by each mean. Thus every outcome
probability `p_z(x)=sum_j |a_(z,j)(x)|^2` obeys the same summed bound.
Finally

```math
I(X;Z)=\sum_z\operatorname{Ent}(p_z)
\le2T\sum_z\mathbb E p_z=2T.
```

This is a direct reconstruction from hypercontractivity, not an external
novelty claim.

### 19.2 The needed strengthened Fano inequality

Suppose `Z` determines an estimate of an address `A=A(X)` with error
at most `delta<1/2` on every input. Let `E` be the error event, let
`pi_a=P(A=a)`, and let `e_a=P(E|A=a)<=delta`. Conditional correctness
determines `A` from `Z`; therefore

```math
H(A|Z)\le h_2(\delta)+P(E)H(A|E).
```

The second term is at most `delta H(A)`, a stronger statement than
ordinary average-error Fano. One proof uses the homogeneous entropy
`S(u)=sum_a u_a log((sum_b u_b)/u_a)`. Its partial derivatives are
`log((sum_b u_b)/u_a)>=0`, so it is coordinatewise nondecreasing.
The error subprobabilities `u_a=pi_a e_a` satisfy `u_a<=delta pi_a`;
hence `P(E)H(A|E)=S(u)<=S(delta pi)=delta H(A)`.
Equivalently, Jensen applied to `t log t` proves the same inequality.
The root independently supplied and checked the monotonicity argument.

Data processing and Section 19.1 now imply

```math
(1-\delta)H(A)-h_2(\delta)
\le I(A;Z)\le I(X;Z)\le2T,
\qquad
H(A)\le\frac{2T+h_2(\delta)}{1-\delta}.
```

Thus no bounded-error quantum addressing scheme, including one obtained
by expanding biased bits into a larger uniform cube, can have accessible
address entropy much larger than its actual expanded-query complexity.
The query cost of computing the substituted bits must be counted.

### 19.3 A bound for arbitrary positive polynomial output kernels

There is a weaker-constant version requiring no quantum realization.
Suppose `p_z>=0`, `sum_z p_z=1`, and every `p_z` has degree at most `D`.
Hypercontractivity and interpolation give

```math
\|p_z\|_4\le3^{D/2}\|p_z\|_2,
\qquad
\|p_z\|_2\le\|p_z\|_1^{1/3}\|p_z\|_4^{2/3},
\qquad
\|p_z\|_2\le3^D\|p_z\|_1.
```

For a nonzero kernel, put `u=p_z/E p_z`. Jensen under the probability
measure of density `u` gives `E[u log u]<=log E[u^2]<=2D log3`.
Thus `I(X;Z)<=2D log3`. Combined with the same worst-case decoding
condition, this yields `H(A)<=(2D log3+h_2(delta))/(1-delta)`.
Signed scalar approximation polynomials are not positive output kernels,
so this does not assert a general approximate-degree lower bound.

### 19.4 Essential limitation: the other part of Fourier entropy

For the Boolean address function `f(x,y)=y_(A(x))`, write
`b_a(x)=1_(A(x)=a)` and `pi_a=E b_a`. Fourier coefficients lie in
disjoint data-bit blocks and are exactly the coefficients of `b_a`.
Parseval gives each block squared mass `pi_a`, hence

```math
H(|\widehat f|^2)
=H(\pi)+\sum_a\pi_a
 H\left(\frac{|\widehat b_a|^2}{\pi_a}\right).
```

The information bound controls only the first summand. The conditional
Fourier entropies of the indicator functions may be large even when the
address-output entropy is small. Consequently this result rejects the
specific high-accessible-label-entropy mechanism, not the entire
approximate-degree/Fourier-entropy route to BH divergence.

## 20. Further independent audits

The counterexample artifact's general entropy/effective-degree theorem,
Section 11, passes a complete independent check. For a unitary target
with Fourier coefficients `a`, a frequency set of mass `M` and pointwise
squared coefficients at most `exp(-L)` supplies the exact dual bound

```math
R_D(Q)\ge
\frac{M^{1/(2D)}(\sqrt M-\delta)_+e^{L/(2D)}}{\|Q\|_\infty},
\qquad \|Q-G\|_2\le\delta.
```

This remains valid when the target frequencies exceed `D`; those
frequencies simply do not contribute on the polynomial side. The
conjugate-exponent estimate also covers `D=1`. A one-sided spectral
typical set suffices for the fixed-seed asymptotic limit, provided cap
control is separately imposed. For moving seeds the information-tail
condition must be checked, as the canonical theorem explicitly requires.

The general radial criterion also passes. If a nonvanishing degree-`d`
polynomial `P` has log modulus spread `lambda`, and its normalized
log-modulus function `R in [0,1]` satisfies
`deg(P R^j)<=d+c` for all positive integers `j`, then joint Taylor
normalization gives effective degree `d+4c lambda`. Its sufficient
BH lower bound is `exp(H(P/|P|)/(2(d+4c lambda)))`. The degree condition
is finitely checkable: if `R` takes `K` distinct values, including zero,
then its value polynomial has degree `K` and zero constant term.
Every `R^j`, `j>=K`, is a linear combination of
`R,...,R^(K-1)`. Thus it suffices to check those first `K-1` powers.
The elementary fallback `c=n-d` always exists, but can be far too costly.

The natural addressed-tangent family in the counterexample artifact,
Section 10, also passes a complete independent check, including its
exact replay for 2, 4, and 8 address leaves. Its orthogonal tangent
directions have squared `L2` norm 6 and explicitly computed new Fourier
mass. Nevertheless, unselected leaves independently align at a vertex:

```math
2\sqrt3(1-2/N)\|\varepsilon\|_1
\le\|\sum_j\varepsilon_j H_j\|_\infty
\le2\sqrt3\|\varepsilon\|_1.
```

Thus bounded variance does not supply bounded cap. At `ell_1` budget
`L`, the putative direction-label entropy is at most `2L^2/e`. A
separate support count gives BH ratio below `2^(3/2)` for the complete
family. This does not exclude correlated directions with another cap
mechanism.

Finally, the positive artifact's orbit-quotient identity and transitive
counterfamily, Section 7, pass independently. The exact orbit-edge
amplification is `(a_OP+a_PO)/2`, not coordinate count or group size.
The Boolean ternary-tree recursion using
`E(x,y,z)=(xy+xz+yz-1)/2` has `n=3^t`, degree `2^t`, and
`n I(f)/m^2 >=(9/8)^t`. This disproves the proposed transitivity-based
dimension amplification while its own BH ratios remain bounded below 2.

## 21. Entropy/effective-degree barriers for two target classes

### 21.1 Boolean targets would require an FEI counterexample

Let `G` be any fixed unitary function, with spectral probability
`pi_S=|Ghat(S)|^2`, and write

```math
I(G)=\sum_S |S|\pi_S.
```

Any sequence of degree-`D_k` polynomials converging in `L2` to
`G^[tensor k]` must satisfy

```math
\liminf_{k\to\infty}D_k/k\ge I(G).
```

Indeed, the Fourier degree of a random coefficient under the tensor
spectral measure is a sum of `k` independent copies of the finite random
variable `|S|`. If a subsequence had `D_k/k<I(G)-eta`, the weak law
would put asymptotically zero mass below `D_k`. Orthogonal projection
then makes the best possible `L2` approximation error tend to one,
contradicting the assumed convergence.

Consequently a family of Boolean targets with
`H(G)/effective_degree(G) -> infinity` would also satisfy
`H(G)/I(G) -> infinity`, disproving the classical Fourier
Entropy--Influence conjecture. This is a logical research barrier, not
an impossibility theorem. For complex unitary targets the analogous FEI
statement already fails for `exp(i theta x)` as `theta -> 0`, so this
argument does not settle the requested complex problem.

### 21.2 Independent weak phases still cannot exploit that FEI failure

**Proved here.** Consider any fixed independent-phase seed

```math
G(x)=\exp\left(i\sum_{j=1}^n\theta_jx_j\right).
```

Global phase changes and coordinate sign flips reduce to
`0<=theta_j<=pi/2`. Put `L=sum_j theta_j`. If polynomials `Q_k` of
degree `D_k` approximate `G^[tensor k]` uniformly with error tending
to zero, then

```math
\liminf_{k\to\infty}D_k/k\ge L/\pi.
```

To prove this, group the `nk` angles as follows. Angles at least `pi/4`
form singleton groups. Greedily combine the remaining angles until each
group's sum first reaches `pi/4`; such a sum is less than `pi/2`.
At the end the ungrouped total is less than `pi/4`. Hence the number
`g_k` of groups satisfies

```math
g_k\ge\frac2\pi(kL-\pi/4),
\qquad \pi/4\le\Theta_\ell\le\pi/2.
```

Identify all original input variables in each group with one new sign,
and fix any ungrouped variables. This diagonal substitution cannot
increase Walsh degree or uniform error. The restricted target, up to a
constant phase, is `product_ell exp(i Theta_ell y_ell)`. Its spectral
degree is a sum of independent Bernoulli variables with probabilities
`sin^2(Theta_ell)>=1/2`. Its mean is at least `g_k/2` and its variance
is at most `g_k/4=O(k)`. If `D_k/k<L/pi-eta` along a subsequence,
Chebyshev puts asymptotically zero spectral mass below `D_k` there.
The best `L2` approximation then has error tending to one, contradicting
the retained uniform approximation. This proves the degree lower bound.

On the other hand, the product spectral entropy is

```math
H(G)=\sum_j h_2(\sin^2\theta_j)\le2\sum_j\theta_j=2L.
```

For example, the elementary estimate `h_2(p)<=2sqrt(p)` follows from
`-p log p<=(2/e)sqrt(p)` and
`-(1-p)log(1-p)<=p<=sqrt(p)`; then use `sin(theta)<=theta`.
Thus every positive uniform-approximation degree rate `d_eff` obeys

```math
\frac{H(G)}{2d_{\rm eff}}\le\pi.
```

The entropy/effective-degree lower-bound criterion is therefore at most
`e^pi` throughout this class, uniformly in the number and sizes of the
phase angles. The constant is intentionally unoptimized. In particular,
the divergence of `H/I` for tiny independent phases does not translate
into an entropy/effective-degree divergence. Arbitrary coded characters
in place of the individual bits are not covered, since an inverse parity
change need not preserve the approximating polynomial's physical degree.

### 21.3 A sharper rotation-velocity lower bound for uniform approximation

**Proved here; independently audited by the positive track.** On the
Hilbert space of functions on the `n`-cube, let `X_j` flip coordinate
`j`, put `J=(1/2)sum_j X_j`, and let `M_G` denote multiplication by a
unitary-valued function `G`. Define

```math
v(G)=\|M_G^*JM_G-J\|_{\mathrm{op}}.
```

If polynomials `Q_k` of degree at most `D_k` approximate `G^[tensor k]`
uniformly with errors `delta_k -> 0`, then

```math
\liminf_k D_k/k\ge v(G).                                      \tag{V}
```

For clarity, this is a **uniform** approximation result. The general
entropy transfer theorem allows merely `L2` approximation together with
cap control; (V) has not been established under that weaker hypothesis.

The operator-valued trigonometric polynomial
`exp(itJ) M_Q exp(-itJ)` has frequencies between `-D` and `D` when
`deg Q<=D`. This follows by expanding each Walsh monomial: conjugation
rotates its single-site `Z` matrices into `Z cos(t)+Y sin(t)` and hence
has frequencies of absolute value at most its support size. Scalar
Bernstein applied to every matrix element, followed by operator-norm
duality, gives

```math
\|[J,M_Q]\|\le D\|Q\|_\infty.
```

For the tensor target, the Hermitian operator
`M_(G^[tensor k])^* J_total M_(G^[tensor k])-J_total` is the tensor sum
of `k` copies of `A=M_G^*JM_G-J`. Its norm is exactly `k||A||`: the
upper bound is the triangle inequality and the lower bound takes the
tensor power of an eigenvector for an eigenvalue of absolute value
`||A||`. Multiplication operators have norm equal to uniform norm, and
`||J_total||=nk/2`. Therefore

```math
k v(G)\le \|[J_{\rm total},M_{Q_k}]\|+nk\delta_k
 \le D_k(1+\delta_k)+nk\delta_k.
```

Divide by `k` and take limits. For a moving seed, the displayed finite
inequality remains valid, but passing to the same conclusion requires
control of `n delta`; fixed-seed reasoning cannot suppress that factor.

For the independent-phase class of Section 21.2, the local operators
have eigenvalues `+/-|sin theta_j|`, so

```math
v(G)=\sum_j|\sin\theta_j|,\qquad H(G)\le2v(G).
```

Thus (V) improves the bound `e^pi` there to `e`. For general unitary
targets the possible estimate `H(G)<=C v(G)` is **not proved here**.
It is stronger than `H(G)<=C deg(G)`, since Bernstein gives
`v(G)<=deg(G)`. For Boolean-valued `G`, `v(G)` is exactly the operator
norm of its sensitivity graph's adjacency matrix, because the off-
diagonal entry on an edge is `(G(x)G(y)-1)/2`.

### 21.4 Exact quadratic-order velocity stability of the tangent seed

For the exact five-bit seed `f` and ten-bit tangent pair `(F,h)` from
Section 9, set

```math
G_\varepsilon=(F+\varepsilon h)/|F+\varepsilon h|,
\qquad F=f\otimes f.
```

The velocity at the base is not 4. The exact replay
`computations/bh_boundedness_barrier_2026_09_18_velocity_exact.py`
certifies

```math
1.865<v(f)<1.867,\qquad
3.730<v(F)=2v(f)<3.734.
```

Here is the reproducible algebraic certificate. The cube graph is
bipartite. If `C` is the even-to-odd block of `2(M_f^*J M_f-J)`, then
`C` has Eisenstein integer entries. The positive eigenvalues of the
five-bit velocity operator are `sqrt(rho)/2`, for the eigenvalues `rho`
of `C C*`. Its rational characteristic polynomial, in descending order,
has coefficient list

```text
1, -100, 4498, -120318, 2132319, -26403192, 234732258,
-1515915552, 7112558643, -23995411096, 56948285089,
-91663372696, 94480461243, -57080829672, 17638228008,
-2109823200, 874800.
```

Exact rational root isolation finds 16 positive simple roots. The two
largest have disjoint isolating intervals

```math
\rho_{15}\in(1753402/125981,2001003/143771),\qquad
\rho_{16}\in(1816813/130444,2027612/145579).
```

The script checks the entire characteristic polynomial, all root
multiplicities, and `13.92<rho_16<13.93`; no floating computation is
needed for simplicity or the bounds above. The tensor sum consequently
has a simple largest positive eigenvalue at `epsilon=0`. For every
`epsilon`, the velocity operator remains bipartite, so its spectrum is
symmetric and its norm equals its largest positive eigenvalue.

The entries of that operator are real-analytic functions of real
`epsilon` near zero. Indeed, `h/F` is `0` or `+/-4i sqrt(3)`, and the
normalization is a factor `(1+48 epsilon^2)^(-1/2)` at the nonzero
vertices. By the simple-eigenvalue perturbation theorem the largest
eigenvalue is real analytic near zero. Swapping the two five-bit blocks
fixes `F` and negates `h`, so it unitarily conjugates the velocity
operator at `epsilon` to the one at `-epsilon`. Hence that analytic
eigenvalue is even, proving

```math
v(G_\varepsilon)=2v(f)+O(\varepsilon^2).                     \tag{V2}
```

The companion numerical script
`computations/bh_boundedness_barrier_2026_09_18_velocity.py` gives
`v(f)=1.8660061245422` and a quadratic quotient near `-3.94195` at
`epsilon=10^-5,10^-4`, with eigenvector residuals below `3e-13`.
The sign and precise value of that quadratic coefficient are only
diagnostics, not exact certificates. Statement (V2), the baseline
bounds, and simplicity are proved. This velocity lower bound therefore
does not explain the Taylor construction's effective degree
`4+16 lambda`; it leaves a nonzero baseline gap and no linear tangent
penalty. It also does not bound the entire BH problem.

## 22. Independent audits of all-unitary structure

The positive artifact, Sections 9.2--9.4, contains two newly proved
structural facts. Both pass an independent full audit here.

First, its tangent-gap lemma is dimension-uniform. If
`||Re(bar f h)||_2>=tau||h||_2` for all base polynomials of degree at
most `m-1`, split a nearby degree-`m` polynomial as `g=a+r`, with
`a=E_fresh g` and zero fresh expectation of `r`. Projecting the unitary
identity gives `2tau||r||_2<=4delta||r||_2`, where
`delta=||g-f||_infinity`. Consequently `delta<tau/2` forces `r=0`.
For the five-bit seed its exact 12-by-12 Gram certificate gives
`tau>=1/2`, hence strict uniform radius `1/4`. The rational replay was
rerun independently and passed. This is local variable rigidity, not
a global five-variable classification.

Second, for every fixed `m`, all scalar circle-valued degree-`m`
functions are juntas of a finite size `J_m` independent of ambient
dimension. The positive artifact gives the complete Ramsey recurrence.
The two delicate quantifiers check as follows. A common acute phase
color on all `m`-subsets of a `2m`-set is impossible: the degree-`2m`
autocorrelation coefficient is a sum of strictly positive paired real
parts, with no outside-variable terms possible. A zero-colored clique
`A` of size greater than `2^(m-1) J_(m-1)` has degree at most `m-1`
under **every** outside restriction. Nevertheless, each originally
relevant variable of `A` survives at least a `2^(1-m)` fraction of
those restrictions: select a nonzero original coefficient containing
it and use the nonzero-polynomial support bound on the resulting
outside coefficient polynomial. Taking expectations, without requiring
independence, contradicts the inductive junta bound.

This is a rigorous qualitative theorem for arbitrary phases, including
nontorsion phases. Its Ramsey growth yields no useful exponential
support estimate and does not imply bounded BH constants. It does
upgrade the positive track's correlated-measure approximation no-go
from torsion atoms to **all** finite-degree unitary atoms: any proposed
dimension-independent finite atom degree and coefficient budget fails
to give a uniform residual contraction below one. The positive artifact
Section 2.8 is the canonical full statement and proof. No external
novelty claim is made for these structural results.

## 23. Rotation velocity does not characterize uniform tensor degree

**Exact counterexample proved here.** Consider the three-bit polynomial

```math
f(x_1,x_2,x_3)
=\frac{1-x_1x_3+i x_2(x_1+x_3)}2
=\exp\!\left(\frac{i\pi}4x_2(x_1+x_3)\right).
```

If `x_1=x_3`, its value is `i x_1x_2`; otherwise its value is 1.
It is therefore unitary-valued and has Walsh degree exactly two. Its
rotation velocity from Section 21.3 is

```math
v(f)=\sqrt3.
```

This identity has a small exact certificate, with no approximate
eigenvalues. Form `A=M_f^*J M_f-J`, `J=(X_1+X_2+X_3)/2`.
Direct Gaussian-rational matrix multiplication gives

```math
A^3=3A,\qquad \operatorname{tr}(A^2)=12,
\qquad \det(tI-A)=t^4(t^2-3)^2.
```

The Hermitian matrix is nonzero, so the first two identities already
force its norm to be `sqrt3`. The replay
`computations/bh_boundedness_barrier_2026_09_18_velocity_converse.py`
checks all eight values, the matrix identities, and the characteristic
polynomial exactly.

Nevertheless the optimal uniform tensor-approximation degree rate is
**exactly 2**, not `sqrt3`. Indeed, in every copy identify `x_3=x_1`.
This substitution does not increase Walsh degree or uniform error.
The resulting `k`-copy target is

```math
i^k\prod_{j=1}^k x_{j,1}x_{j,2},
```

a parity character on `2k` independent bits, up to phase. A polynomial
of degree strictly less than `2k` is orthogonal to this target, so its
`L2` error, and hence its uniform error, is at least one. Consequently
every approximation with uniform error below one has degree at least
`2k`. The original exact tensor polynomial has degree `2k`, attaining
that rate.

Thus the proposed recovery assertion
`G^[tensor k]` admits uniform `o(1)` error at degree `(v(G)+o(1))k`
is false, even for a three-bit torsion-valued degree-two seed. Cutting
off scalar rotation frequencies cannot repair this failure by itself.

More generally define the coordinate-identification lower envelope

```math
v_{\rm id}(G)=\sup_\sigma v(G\circ\sigma),
```

where every original coordinate is replaced by a constant sign or a
signed new coordinate (new coordinates may be reused). Copywise
substitution and Section 21.3 prove that every uniform tensor degree
rate is at least `v_id(G)`. The example has `v_id(f)=2`: the displayed
restriction has velocity 2, while every allowed substitution preserves
degree at most two and Bernstein bounds its velocity by two. No
general recovery theorem at this stronger envelope is proved here.

The five-bit seed also exhibits a strict envelope gap. Substituting
`x_1=x_2=x_3=x_4=y_1`, `x_5=y_2` gives values
`(zeta^2,zeta^5,zeta^5,zeta^3)` in binary input order. The even-to-odd
block of twice its velocity operator has identical columns
`(-2,zeta^2-1)`, hence velocity `sqrt(7/2)>1.8708`, exceeding the
certified original velocity below `1.867`. This larger example motivated
the simpler theorem above; no numerical optimization is needed for
either strict gap.

These facts do not invalidate the persistent BH lower bound in Section
14 or the general entropy-transfer theorem in Section 20: those results
construct actual approximating polynomials and compute their degree
directly. Neither invokes a velocity converse.

## 24. Independent audit of skew-orthogonal phase coupling

The counterexample artifact Section 12 passes a full independent audit.
For its local random variable `u=R+iI`, the exact moments are
`E R=1/4`, `E I=0`, `E R^2=E I^2=3/2`, and `E RI=0`.
For real skew-symmetric `K`, its tangent scalar is `Q=2i I^T K R`.
Expanding the second moment independently gives

```math
\|Q\|_2^2
=4\left[(9/4-3/32)\|K\|_F^2+(3/32)\|K\mathbf1\|_2^2\right],
```

which is exactly the asserted formula `(69/8)||K||_F^2+
(3/8)||K1||_2^2`. The random-cut lower bound
`||Q||_infinity>=sum_i ||K_(i,*)||_2` is also valid: conditional on
membership in one side, a row's retained squared mass has mean half its
total, and `sqrt(S)>=S/sqrt(v)` supplies the required half-row norm.
The elementary fourth-moment Khintchine bound and the cut probability
then give the stated constant. Thus orthogonal `K` has cap at least
`N` versus `L2` norm `3sqrt(N)`.

The subsequent entropy argument avoids an invalid ambient-dimension
continuity bound. The analytic phase normalization has geometric-error
polynomial approximations supported on nested sets of size
`exp[O((j+1)log(eN))]`. The first-support-level variable has a geometric
tail, bounded mean, and bounded entropy; the entropy chain rule yields
`H<=C_c log(eN)` at fixed phase cap `c`.

Finally, its lower bound for a moving number `r_N` of target copies is
uniform in that number. Write `alpha_N` for the probability that one
target coefficient has degree below `k/3`; then `alpha_N -> 0`.
For every positive integer `r_N`, the expected fraction of such bad
copies is `alpha_N`. Markov bounds that fraction by `sqrt(alpha_N)`
with probability at least `1-sqrt(alpha_N)`. Consequently the tensor
spectral degree is at least `(1-sqrt(alpha_N))k r_N/3` with that
probability. Any `L2-o(1)` approximation must therefore have degree
at least `(1/3-o(1))k r_N`. No fixed-copy or independence-in-the-limit
assumption is hidden in this step. Together these estimates rule out
this bounded-cap skew-orthogonal family as a divergent entropy/degree
mechanism, not arbitrary nonlocal phase constructions.

## 25. Further exact audits and a failed degree-two simplification

### 25.1 Quadratic phase entropy

The counterexample artifact Section 13 proves
`H(exp(itP))<=8sqrt(3)|t|||P||_infinity` for every real hollow quadratic
`P=sum_(i<j) a_ij x_i x_j`; its proof and 332-case identity replay pass
independently here. The exact spectral marginal is
`p_i=E sin^2(t sum_j a_ij x_j)`. Thus entropy subadditivity and
`h(p)<=2sqrt(p)` give `H<=2|t|sum_i ||A_i||_2`. The elementary
Khintchine lower bound and real polarization give
`sum_i ||A_i||_2<=4sqrt(3)||P||_infinity`, including arbitrary rank and
dimension. Its consequence for the direct Taylor approximation's
declared degree rate is valid; it does not lower-bound an optimal
approximation rate.

There is a straightforward extension to arbitrary real degree-at-most-
two polynomials, including linear and constant terms. Let `c=P(0)` be
the constant coefficient and introduce one sign `x_0`:

```math
\widetilde P(x_0,x)=P(x_0x)-c.
```

This is hollow homogeneous quadratic, has cap at most `2||P||`, and
`exp(it Ptilde)=exp(-itc) exp(itP)(x_0x)`. Its Fourier coefficients are
exactly those of `exp(itP)`, with an additional `x_0` factor on odd
support sizes; no coefficients collide. Entropy is consequently
preserved. This proves

```math
H(\widehat{e^{itP}}^{,2})
 \le16\sqrt3\,|t|\|P\|_\infty\qquad(\deg P\le2).
```

Here and throughout, the shorthand entropy notation means the squared
absolute Fourier coefficients. This fixed-degree result does not
establish the proposed linear-in-degree bound for general real `P`.

### 25.2 Normalization destroys the tested transfer-pair inheritance

The counterexample track's exact modular replay
`computations/bh_boundedness_counterexamples_2026_09_18_inheritance.py`
was read and rerun independently. For the normalized ten-bit seed
`G=G_epsilon`, write `J_G A=G^2 bar A`, and take either
`epsilon=2^-10` or `2^-104`. For possible degree budgets `(s,8-s)`,
the ranks over the prime 1009, with sixth root `zeta=375`, are

```text
s             0       1       2        3         4
columns       1      11      56      176       386
rank          0      11      56      176       385
```

The convolution matrix enforces that `G^2 bar A` has no Fourier
coefficient above degree `8-s`. Its entries are reductions of exact
elements of `Q(zeta)`, and all common denominators are nonzero modulo
1009. A nonzero modular minor is therefore a nonzero characteristic-
zero minor; full ranks for `s=1,2,3` exclude the corresponding nonzero
complex solutions. The symmetric case `s=4` has dimension at most one;
the known vector `P=F+epsilon h` satisfies `J_G P=P`, showing that
dimension is exactly one. Scalar multiples yield zero in the proposed
antisymmetric tensor pairing. Applying the involution `J_G` covers
the reversed asymmetric pairs as well.

The physical-value formula used by the code was independently checked.
With `D=1/epsilon` and `q=h/F`, it is

```math
G^2=F^2\frac{D^2+48-96R+2Dq}{D^2+48},
\qquad q^2=-48R,\quad R\in\{0,1\}.
```

This also shows `deg G^2<=8`, explaining the surviving constant
`(0,8)` pair. The computation is an exact obstruction for these two
specified epsilon values and this budget, not a generic-epsilon or
all-degree classification.

### 25.3 Why a small connected degree-two support classification fails

An attempted simplification was to bound each connected component of
the quadratic Fourier interaction graph by the five variables of the
flat seed. That proposal is false. Define the Boolean degree-two
selector

```math
b(a,b,c,d)=\frac{(a+b)c+(a-b)d}{2}\in\{-1,1\}.
```

Take two such selectors `g,h` on sets of four coordinates sharing
exactly one coordinate. Then

```math
f=(g+i h)/\sqrt2
```

is unitary-valued of degree two, depends on seven variables, and its
interaction graph is the connected union of two four-cycles sharing
one vertex. Disjoint selectors give eight relevant variables and two
components. Each construction has eight Fourier coefficients, so
neither challenges the conjectural support bound `4^2=16`.

A related six-variable connected example is
`((x_0+y_0)u+(x_0-y_0)v)/2`, where
`u=(z_1+i z_2)/sqrt2` and `v=(z_3+i z_4)/sqrt2`.
Its graph is `K_(2,4)`. Thus affine-output bounds and the five-bit local
rigidity theorem do not by themselves classify degree-two components.
The exact decomposition `f=g+x_i h`, with affine `h`, leads to
`|g|^2+|h|^2=1` and `Re(g bar h)=0`, but no dimension-independent
support estimate sharper than the established Ramsey result was
derived from these equations in this track.

## 26. Correlated-state scalarization of the new quantum family

The root requested a test beyond separate product-state rounding and the
already analyzed disjoint-slot Hadamard chains. The exact construction
in [Slote, arXiv:2608.01424, Section 4](https://arxiv.org/html/2608.01424)
has the following structure, which is all that is used here. Write
`d=r+k`, let `L` be the largest power of two below `3^r`, and put

```math
A=N^{-1/2}\sum_{\ell=1}^L\sum_{a\in[L]^k}
 \eta_{\ell,a}Q_{\ell,a},\quad
Q_{\ell,a}=P_\ell\otimes Z_a,\quad
N=L^{k+1},\quad \eta_{\ell,a}\in\{-1,1\}.
```

The `P_l` are pairwise anticommuting Hermitian Pauli operators; all
suffixes `Z_a` commute. The source operator is Hermitian and unitary.
The following statements concern this particular source family and do
not assume that the evaluating state factors across any qubits.

### 26.1 Sharp aggregate attenuation for arbitrary correlated states

**Proved here.** For every density matrix `rho`,

```math
\sum_{\ell,a}|\operatorname{tr}(\rho Q_{\ell,a})|^2
 \le L^k=N/L.                                                \tag{NC1}
```

For a fixed suffix `a`, the `L` operators `Q_(l,a)` anticommute and
square to the identity. Let `v_l=tr(rho Q_(l,a))`, which is real.
The operator `sum_l v_l Q_(l,a)` has square `||v||_2^2 I`, so its
operator norm is `||v||_2`. Its expectation in `rho` is `||v||_2^2`.
Hence `||v||_2<=1`; sum over `a`. This is sharp: choose a core state
with one `P_l` expectation equal to one and a computational-basis
suffix state, giving a unit squared sum for every `a`.

More strongly, for any input-dependent density matrices `rho(x)`, put
`u_(l,a)(x)=tr(rho(x) Q_(l,a))`. Pointwise application of (NC1), followed
by Parseval, yields

```math
\sum_{\ell,a,S}|\widehat u_{\ell,a}(S)|^2\le N/L.       \tag{NC2}
```

Thus correlated, entangled, adaptively selected, and stabilizer states
all obey the same total coefficient-energy budget. No preparation
degree assumption is needed for this inequality.

In particular, suppose a termwise scalarization associates a distinct
Walsh character `chi_(phi(l,a))` to each retained Pauli term and extracts
its coefficient `alpha_(l,a)` from `u_(l,a)`. Its diagonal extracted
component has Fourier energy at most

```math
\frac1N\sum_{\ell,a}|\alpha_{\ell,a}|^2\le\frac1L.
```

To obtain, after a scalar multiplication `K`, `M` distinct extracted
coefficients of magnitudes at least `c/sqrt(M)` therefore requires

```math
K\ge c\sqrt L=\exp(\Omega(d))
\quad\text{when }r=\Theta(d).                              \tag{NC3}
```

This conclusion allows `M` to be smaller than `N`: merely retaining
`exp(c_0 d^2)` coefficients does not remove the lost Parseval mass.
A state-evaluation map is unital and norm-contracting; multiplying that
map by `K` has operator-to-uniform norm exactly `K` (test the identity).
Thus this termwise transfer cannot have logarithmic norm loss `o(d)`
while preserving a positive amount of Parseval-normalized coefficient
mass. The result is independent of whether its output degrees are
inflated.

There are two important scope limits. First, if contributions from
different Pauli terms deliberately collide in the same output Fourier
coefficient, (NC2) does not prevent coherent sums from increasing the
scalar coefficient. Second, for the one chosen source operator the
actual scalar cap may be smaller than the map's guaranteed cap;
normalizing by that cap is not ruled out by (NC3). Such a construction
must be analyzed as a new scalar polynomial, rather than justified by
coefficient preservation and contractivity alone. No general
noncommutative-to-classical impossibility claim is made.

### 26.2 Fixed Pauli orbits and stabilizer evaluations

For a fixed arbitrary density matrix `rho`, evaluate on its full Pauli
orbit `rho_x=W_x rho W_x^*`, with two sign bits per physical qubit.
Conjugating any Pauli `Q` multiplies it by its distinct symplectic Walsh
character. Hence no coefficient collision occurs, and the exact scalar
polynomial `f_rho(x)=tr(rho_x A)` satisfies

```math
\|\widehat f_\rho\|_2^2\le1/L,\qquad
\deg f_\rho\le2d,\qquad \|f_\rho\|_\infty\le1.
```

This is a concrete correlated-state compiler with a sharp exponential
loss of Parseval mass. Arbitrary entanglement in `rho` cannot improve
that mass bound. Fixed Clifford changes of basis do not help, since
they preserve all the anticommutation relations in the proof.

If `rho` is a pure stabilizer state, the support bound is stronger:
at most one core label `l` can have nonzero expectations among all
`Q_(l,a)`. Indeed every pair with different core labels anticommutes,
regardless of its suffix, while nonzero Pauli expectations in a
stabilizer state belong, up to sign, to one commuting stabilizer group.
Thus at most `L^k=N/L` of the original terms survive. The assertion
holds also for a normalized stabilizer-code projector. It does not say
that an input-dependent choice of stabilizer group uses the same core
label for every input; (NC2) is the correct uniform statement there.

### 26.3 Exact norm-saturation forces quadratic degree

There is also an exact output-side tradeoff that includes collisions.
Any density-state evaluation of Hermitian `A` is real-valued. If a
rescaled scalar output has `||f||_2=||f||_infinity=1`, then it is
Boolean-valued pointwise. For a nonconstant Boolean polynomial of
degree `m`, granularity gives at most `4^(m-1)` nonzero Fourier
coefficients. Consequently, exact Parseval saturation with
`M=exp(c d^2)` retained scalar coefficients forces

```math
m\ge1+\log_4 M=\Omega(d^2).
```

Therefore an exact norm-saturating stabilizer/correlated-state
scalarization cannot retain the quantum family's quadratic logarithmic
support at degree `O(d)`. This exact statement does not extend by
continuity to an unspecified near-saturation regime. Together with
(NC3), it isolates two concrete losses for the actual quantum family;
it does not repeat the Hadamard-chain argument or settle arbitrary
mixing-based scalarizations.

### 26.4 Complex norm-one functionals, not only positive states

The counterexample track supplied the following extension, independently
checked here. Let `phi_x` be any input-dependent complex linear
functional on the operator algebra with dual operator norm at most one.
This includes arbitrary normalized matrix entries and need not be a
positive state. For a fixed suffix and complex scalars `c_l`,

```math
\left\|\sum_l c_l Q_{l,a}\right\|
 \le\|\operatorname{Re}c\|_2+\|\operatorname{Im}c\|_2
 \le\sqrt2\,\|c\|_2.
```

The first inequality follows by separating real and imaginary parts;
each real linear combination squares to its Euclidean norm squared
times the identity. Duality therefore gives
`sum_l |phi_x(Q_(l,a))|^2<=2`. Sum over suffixes and apply Parseval:

```math
\sum_{l,a,S}|\widehat{\phi_x(Q_{l,a})}(S)|^2\le2N/L.
```

Consequently the collision-free extracted coefficient transfer requires
rescaling at least `c sqrt(L/2)`, still exponential in `d`. This removes
positivity and real-output restrictions from the attenuation statement.
The factor two is already sharp for one anticommuting Pauli pair,
whose transition matrix element may have values `1,i`. It does **not**
extend the separate Boolean norm-saturation argument of Section 26.3
to arbitrary complex functionals; that argument explicitly uses real
state-evaluation outputs.

## 27. Two independently checked restricted positive theorems

### 27.1 A small control block plus an arbitrary affine block

The director artifact Section 6 passes independently. For real
`P(z,x)=p_0(z)+sum_i x_i p_i(z)` of degree at most an integer `d>=1`,
with `k` control signs `z`, arbitrary many data signs `x`, and cap `Q`,
it proves

```math
H(|\widehat{e^{itP}}|^2)
 \le k\log2+2|t|Q+2t^2(d-1)Q^2.
```

The data Fourier amplitudes, after omitting the irrelevant scalar
phase from `p_0`, are tensor products of sine/cosine pairs. Their
conditional entropy is at most `2|t|sum_i |p_i(z)|<=2|t|Q`.
The mutual information between control input and data Fourier label
equals a sum of logarithmic Sobolev entropies. Its Dirichlet energy is
at most `t^2 sum_(i,j) E|D_j p_i|^2`, because the amplitude inner
product is a product of cosines and
`1-product c_i<=sum(1-c_i)` for all real `c_i in [-1,1]`.
The constant in the complex hypercube log-Sobolev inequality is 2.
Together with `deg p_i<=d-1` and `sum_i E p_i^2<=Q^2`, this gives
the claimed quadratic-in-phase term. The remaining control Fourier
label has at most `2^k` possibilities, costing `k log2`.

The `d>=1` qualification matters for the displayed formula; constants
can instead be treated separately. The director was informed of this
minor quantifier guard. Negative cosine factors cause no gap, and
discarding `p_0` before the mutual-information estimate is legitimate.

### 27.2 Sharp bipartite quadratic unitary coefficient bound

The positive track derived the following theorem, and its complete
argument was independently reconstructed here. If a unitary-valued
complex polynomial of degree at most two has a bipartite quadratic
interaction graph, then it has at most eight Fourier coefficients.
Consequently its `q_2=4/3` BH ratio is at most `8^(1/4)=2^(3/4)`.
Both bounds are sharp.

To verify the reduction, write the polynomial with two vertex blocks
as `f(x,y)=c+sum a_i x_i+sum b_j y_j+sum A_ij x_i y_j`. Introduce
one homogenizing sign per block:
`B(x_0,x,y_0,y)=x_0y_0 f(x_0x,y_0y)`. This is a unitary bilinear
form whose nonzero matrix entries correspond exactly to the original
Fourier coefficients. For fixed `y`, at most two nonzero row sums can
occur, because their complex vectors must be pairwise orthogonal in
the real plane. Each nonzero homogeneous linear row is nonzero on at
least half the cube, so the matrix has at most four nonzero rows and,
symmetrically, four nonzero columns.

If four rows occur, each must attain nonzero probability exactly
one half. A nonzero homogeneous complex linear form attains this only
when it has exactly two nonzero coefficients. Indeed, if any two
coefficients are not equal up to sign, conditioning on the others
makes the four two-sign sums distinct, giving zero probability at
most one quarter. Otherwise all coefficients equal one common
nonzero number up to signs, and the central binomial probability is
one half only for two terms. Thus four rows, or four columns, give at
most eight entries. With at most three of each, more than eight
entries would be a full 3-by-3 matrix. But a three-term nonzero
homogeneous linear form is nonzero at least three quarters of the
time, so its three rows would have expected active count at least
`9/4>2`, impossible.

For sharpness, take disjoint four-variable Boolean selectors `g,h`
from Section 25.3 and set `f=(g+i h)/sqrt2`. The interaction graph is
bipartite, the eight coefficients all have magnitude `1/sqrt8`, and
the range is on the unit circle. This gives ratio `8^(1/4)` exactly.
The five-bit 16-term seed is outside this bipartite class, so there is
no conflict with its ratio 2.

## 28. Controlled disjoint Boolean blocks: an all-time entropy theorem

**Proved here, extending the director's rare-output argument.** Let

```math
P(z,x)=p_0(z)+\sum_{i=1}^N p_i(z)b_i(x^{(i)})
```

be real-valued. The input blocks `z,x^(1),...,x^(N)` are disjoint.
Every `b_i` is a nonconstant Boolean-valued polynomial of degree at
most `D>=1`. The control block `z` has `k` signs, while the number and
sizes of the Boolean blocks are arbitrary. Set `Q=||P||_infinity`.
There is no degree assumption on the control functions. For every real
`t`,

```math
H(|\widehat{e^{itP}}|^2)
 \le (2D+4+2k\log2)|t|Q.                                \tag{CB}
```

Constant Boolean blocks may be absorbed into `p_0`. Because the
remaining disjoint Boolean blocks independently attain both signs,
maximization over their inputs gives the pointwise inequality
`|p_0(z)|+sum_i |p_i(z)|<=Q`.

For any nonconstant Boolean `b` of degree `d`, its Fourier support
has size at most `4^(d-1)`, so `H(|bhat|^2)<=2(d-1)log2`. Write
`u=t p_i(z)`. Since `exp(iu b)=cos u+i sin u b` and the coefficients
of `b` are real, its spectral probability distribution is exactly

```math
\cos^2u\,\delta_{\varnothing}+\sin^2u\,|\widehat b|^2.
```

This formula includes biased Boolean `b`: the two contributions at
the empty coefficient are orthogonal real/imaginary amplitudes, so
there is no cross term. Entropy of a mixture is at most the entropy
of its mixing label plus the averaged conditional entropy. Therefore

```math
H(|\widehat{e^{iu b}}|^2)
 \le h_2(\sin^2u)+\sin^2u\,H(|\widehat b|^2)
 \le2d|u|.                                                \tag{CB1}
```

The last inequality uses `h_2(sin^2u)<=2|u|`,
`sin^2u<=|u|`, and `log2<=1`. No balancing hypothesis on `b` is
needed. The conditional data-frequency law is a product across
blocks, so (CB1) gives
`H(S|Z)<=2D|t|sum_i E|p_i|<=2D a`, where `a=|t|Q`.

When `a<=1`, the conditional probability of a nonempty data-frequency
label satisfies

```math
\Pr(S\ne\varnothing\mid Z=z)
 \le\sum_i\sin^2(tp_i(z))\,[1-(\mathbb E b_i)^2]
 \le a^2.
```

Let `E` denote that event. Since the data label is constant when `E`
fails, `I(Z;S)<=h_2(Pr E)+Pr(E)H(Z|E)<=2a+a^2 k log2`.
For the final control Fourier label `T`, Parseval gives

```math
\Pr(T\ne\varnothing)
 =\mathbb E_x\operatorname{Var}_z(e^{itP})
 \le\mathbb E_{x,z}|e^{itP}-1|^2\le a^2.
```

Hence `H(T)<=2a+a^2 k log2`. Combining the two different joint laws
through their common data marginal gives
`H(S,T)<=H(S|Z)+I(Z;S)+H(T)`, proving (CB) for `a<=1`.
For `a>=1`, use `I(Z;S)<=k log2` and `H(T)<=k log2`; the result is
even smaller. The case `a=0` is immediate.

### 28.1 Shared control partitions and decision trees

**Director extension, independently audited here.** Suppose all control
functions `p_0,p_1,...,p_N` are constant on a common partition into
`L` nonempty cells. Let `K` be the cardinality of the union of the
Fourier supports of their cell indicators. Then the sharper structural
form is

```math
H(|\widehat{e^{itP}}|^2)
 \le(2D+4+\log L+\log K)|t|Q.                           \tag{CB2}
```

Indeed, replace the physical control input by its cell label `C`,
which has at most `L` values. Conditional data amplitudes are constant
on each cell, including their phase from `p_0`; every final control
Fourier frequency therefore lies in the indicated union of size `K`.
The same rare-output argument gives
`I(C;S)<=2a+a^2 log L` and `H(T)<=2a+a^2 log K` for `a<=1`,
while their large-time bounds are `log L` and `log K`. Together with
the unchanged conditional bound `2D a`, this proves (CB2).

If the entire tuple of controls is computed by one deterministic
decision tree of depth `h`, there are at most `2^h` leaves. Each leaf
indicator is a product of at most `h` literals and has Fourier
support at most `2^h`. Consequently `L<=2^h`, `K<=4^h`, giving

```math
H(|\widehat{e^{itP}}|^2)
 \le(2D+4+3h\log2)|t|Q.
```

This bound is independent of the actual number of control coordinates
and the number or arities of the Boolean inner blocks. It does require
a **shared** partition/tree for all control functions, including
`p_0`, and disjoint data blocks. Neither separate shallow trees with
an uncontrolled common refinement nor overlapping Boolean blocks are
covered. The theorem supplies a substantial class for the linear
phase-entropy target; it does not establish that target universally
or decide boundedness of the BH constants.

## 29. Cap-bounded L2 approximation does not obey the velocity lower bound

**Explicit counterexample derived here using a primary quantum-query
construction.** The hypothesis of uniform approximation in Section 21.3
cannot be replaced by `L2` convergence, even if every approximating
polynomial has uniform norm at most one and the target seed is fixed and
Boolean-valued.

On `n=16` bits let `b=1-2 AND_n`, where `AND_n` is the indicator of
the all-`+1` vertex. Its rare value `-1` has probability `p=2^-16`.
Its sensitivity graph consists of one star with 16 leaves and isolated
vertices, so

```math
v(b)=4.
```

Nevertheless there are real polynomials `Q_k` on `16k` bits with

```math
\|Q_k\|_\infty\le1,\qquad
\|Q_k-b^{\otimes k}\|_2\longrightarrow0,\qquad
\limsup_k\frac{\deg Q_k}{k}\le\frac1{\sqrt2}<4.        \tag{L2V}
```

Here is the required bounded outer polynomial. The exact-search
enumeration construction in [de Wolf, arXiv:0802.1816v2, Appendix A](https://arxiv.org/html/0802.1816)
uses successive exact Grover searches with hypothesized marked counts
`t,t-1,...,1`, verifying and removing any found marked index. It finds
all marked indices with certainty whenever their number is at most
`t`. Including one verification query per round, a safe query bound is

```math
T\le\sum_{j=1}^t\left(\left\lceil\frac\pi4\sqrt{k/j}\right\rceil+1\right)
 \le(\pi/2+2)\sqrt{kt}<4\sqrt{kt}.
```

The loop invariant is elementary: before the round whose hypothesis is
`j`, at most `j` marked indices remain. If fewer remain, removing either
zero or one leaves at most `j-1`; if exactly `j` remain, exact search
removes one with certainty. The verification prevents false removals.
Output the parity of the number of verified marks found. On unrestricted
inputs this is still a sign-valued random output. The same paper's
polynomial-method statement makes its expectation a real polynomial
`q_(k,t)` of degree at most `2T`, bounded by one on the entire cube.
On every input with Hamming weight at most `t` it equals parity exactly.

Now set `t_k=ceil(2pk)` and feed this outer algorithm the `k` disjoint
indicators `AND_16`. For sufficiently large `k`, `1<=t_k<=k`; the
finitely many earlier values may be defined arbitrarily. Polynomial
composition and multilinear reduction give

```math
Q_k=q_{k,t_k}(AND_{16}^{(1)},\ldots,AND_{16}^{(k)}),\qquad
\deg Q_k\le8\cdot16\sqrt{k t_k}.
```

The cap is at most one, because each inner input is exactly a bit and
the outer expectation is in `[-1,1]`. Under the uniform original cube,
the number `W_k` of marked blocks is `Bin(k,p)`. The target is exactly
`(-1)^(W_k)`, so the only possible discrepancy is on `W_k>t_k`.
A Chernoff estimate yields

```math
\|Q_k-b^{\otimes k}\|_2
 \le2\sqrt{\Pr(W_k>t_k)}\le2e^{-pk/6}\longrightarrow0.
```

Finally,

```math
\limsup_k\deg(Q_k)/k
 \le128\sqrt{2p}=1/\sqrt2,
```

proving (L2V). The quantum routine is only a construction of ordinary
real bounded polynomials; no quantum input distribution is used in the
result. This counterexample uses the source's exact-query subroutine,
not a guessed approximate-degree constant or numerical optimization.

More generally the same proof gives degree rate at most
`8sqrt2 n 2^(-n/2)` for the `n`-bit rare-AND sign, while its velocity
is `sqrt(n)`. The ratio of these rates tends to zero. Thus no universal
positive constant multiple of `v(G)` can lower-bound cap-bounded `L2`
tensor degree rates for all fixed unitary seeds.

This does **not** produce BH divergence: the rare-AND seed also has
very small Fourier entropy. It also does not prove that the same
failure occurs within the narrower family of hollow quadratic phases.
The exact consequence is a scope warning: a theorem `H(G)<=C v(G)`
alone only obstructs the uniform-compiler version of the entropy
mechanism. It cannot be plugged into the general cap-bounded `L2`
transfer theorem without an additional argument.
Section 32 supplies that additional argument for quadratic phases by
replacing full operator velocity with a flat-diagonal state certificate.

## 30. Sharp all-degree block-affine unitary support: independent audit

The positive track's all-degree refinement passes a complete independent
audit. Its canonical theorem is: if a complex unitary-valued function
is affine separately in each of `k>=1` disjoint variable blocks, then

```math
|\operatorname{supp}\widehat f|\le4^k/2.
```

At the degree budget `m=k`, its critical BH ratio is therefore at most
`2^(1-1/(2k))`, and equality is attained. For other budgets the direct
support consequence is `(4^k/2)^(1/(2m))`; block count cannot silently
be replaced by actual degree if they differ.

The key preservation lemma is exact. Fix all but one block. A unitary
affine slice `c+sum a_i x_i` has at most two relevant variables by real
orthogonality in the complex plane. If there are two, its constant term
must be zero. Every odd postprocessing of this two-bit slice remains
odd, and therefore affine. Slices with at most one relevant variable
remain affine under arbitrary postprocessing. It follows that **every
odd function of `f` remains separately affine in all `k` blocks**.
Odd maps are defined on `range(f) union -range(f)`, so an observed phase
whose negative is not observed causes no problem.

Choose orientations `rho_j` of the finitely many antipodal phase pairs
and let `g_j` be the corresponding signed indicators, with values
`1,-1,0`. For every independent choice of signs `sigma_j`, the function
`h_sigma=sum_j sigma_j g_j` is Boolean-valued, separately affine, and
has total degree at most `k`. Boolean Fourier granularity gives
`sum_j sigma_j ghat_j(S) in 2^(1-k) Z`. Flipping one sign yields
`w_(j,S)=2^k ghat_j(S) in Z`; choosing all signs positive shows that
`sum_j w_(j,S)` is even. Thus every nonzero vector `w_(.,S)` has
squared Euclidean norm at least two. Parseval and the signed partition
give `sum_(S,j) ghat_j(S)^2=1`, proving that the union of their Fourier
supports has size at most `4^k/2`. The support of
`f=sum_j rho_j g_j` is contained in that union.

For sharpness take two disjoint centered Boolean addressing functions
`g,h`, each of degree `k`, with `4^(k-1)` flat Fourier coefficients.
Regroup their corresponding address and data blocks into `k` shared
blocks, and set `f=(g+i h)/sqrt2`. The coefficient supports are
disjoint, every coefficient has magnitude `(4^k/2)^(-1/2)`, and the
range lies on the unit circle. This realizes both the support bound
and the claimed coefficient ratio.

The preservation lemma is special to separately affine blocks. It is
not a general functional-calculus bound for arbitrary degree-`m`
unitaries: the five-bit degree-two seed already has its odd power
`f^3` of degree five. Direct exact convolution gives degree sequence
`2,4,5,4,2,0` for its powers `f,...,f^6`. Thus even the tempting
general bound `deg(odd phi(f))<=2deg(f)` is false.

## 31. All real quadratic Hamiltonian phases have entropy at most 64v

**Proved here by extending the director's hollow-quadratic argument;
independently audited by the director.** For every real polynomial
`P` of degree at most two, and every real time `t`, put `G=e^(itP)`.
Then, in natural logarithms,

```math
H(|\widehat G|^2)\le64\,v(G),\qquad
v(G)=\|M_G^*J M_G-J\|,\quad J=\tfrac12\sum_i X_i.       \tag{QV}
```

The director's Section 7 gives the sharper constant `152/3<52` when
the Hamiltonian is hollow homogeneous quadratic. The proof below
includes arbitrary linear and constant terms, arbitrary dimension,
arbitrary coefficient sizes, and all times. Absorb `t` into the angles.
Reducing each linear and quadratic angle modulo `pi/2` writes

```math
G=\eta\,\chi_D\exp\!\left(i\sum_i\alpha_i x_i
                  +i\sum_{i<j}\theta_{ij}x_i x_j\right),
\qquad |\alpha_i|,|\theta_{ij}|\le\pi/4,\quad |\eta|=1.
```

Define

```math
c_i=\cos(2\alpha_i)\prod_{j\ne i}\cos(2\theta_{ij})\ge0,
\qquad q_i=(1-c_i)/2\in[0,1/2].
```

The spectral marginal of the reduced phase at coordinate `i` is
`q_i`, by the coordinate-difference/Parseval identity. Multiplication
by `chi_D` complements those marginals in `D` and translates the
entire spectral distribution. Consequently

```math
H(|\widehat G|^2)\le\sum_i h_2(q_i),\qquad
v(G)\ge I(G)\ge\sum_i q_i.                              \tag{QV1}
```

The first velocity inequality follows by taking the constant unit
vector in the variational formula: its expectation of
`M_G^*J M_G-J` is `-I(G)`.

We control small marginals more strongly. Fix a row with `q_i<=1/32`,
and abbreviate

```math
u_0=\tan(2\alpha_i),\quad u_j=\tan(2\theta_{ij}),\quad
T=u_0^2+\sum_{j\ne i}u_j^2.
```

All these tangents are finite, since `c_i>=15/16`. The product formula
for `c_i` implies

```math
4q_i\le T\le c_i^{-2}-1\le31/225<1/7.                  \tag{QV2}
```

For the lower bound, `c_i>=exp(-T/2)>=1-T/2`; for the upper bound,
`c_i^(-2)=product(1+u_j^2)>=1+T`.

Take independent auxiliary signs `x_j`. A product qubit state with
`<Z_j>=x_j/2` has

```math
M_i(x):=\left\langle\sin\left(2\alpha_i+
                         2\sum_{j\ne i}\theta_{ij}Z_j\right)\right\rangle
=c_i\operatorname{Im}\left[(1+i u_0)
                      \prod_{j\ne i}(1+i u_j x_j/2)\right].
```

Let `W=u_0^2+sum_j u_j^2/4`, so `q_i<=W<=T<1/7`.
Introduce a fictitious independent sign `x_0` in the first factor.
The imaginary part of the full product is an odd multilinear
polynomial in all these signs. Its absolute first and second moments
conditional on `x_0=1` equal their unconditional moments: changing
`x_0` is absorbed by flipping all the other signs and the overall sign.
Thus the fictitious sign does not alter the moments being bounded.

The linear part is `L=u_0x_0+sum_j u_jx_j/2`. Its second moment is
`W`, and the fourth-moment Khintchine argument gives
`E|L|>=sqrt(W/3)>=sqrt(W)/2`. The remaining odd-degree part `R` has

```math
\|R\|_2^2
 \le\sum_{s\ge3}\frac{W^s}{s!}\le W^3,
\qquad\|R\|_2\le W^{3/2}<\sqrt W/7.
```

Here the elementary-symmetric sum of the squared coefficients at
level `s` is at most `W^s/s!`. The final exponential-series bound
holds for `W<=1` since the remainder is at most `e W^3/6<W^3`.
The triangle inequality now gives

```math
\mathbb E|M_i|
 \ge c_i(1/2-1/7)\sqrt W
 \ge\tfrac14\sqrt{q_i}.                               \tag{QV3}
```

It remains to realize all these row contributions simultaneously.
Write `sigma_i=(-1)^(1_[i in D])`. Direct Pauli conjugation gives

```math
M_G^*J M_G-J
=\frac12\sum_i\left[
 \sigma_i X_i\cos(2\alpha_i+2\sum_j\theta_{ij}Z_j)
 +\sigma_i Y_i\sin(2\alpha_i+2\sum_j\theta_{ij}Z_j)-X_i\right].
```

For each auxiliary sign vector choose a product state whose local
Bloch coordinates are
`<X_i>=0`, `<Z_i>=x_i/2`, and
`<Y_i>=sigma_i sign(M_i(x)) sqrt(3)/2`. These choices are mutually
compatible: each local Bloch vector has length one, and its `Y` sign
can be chosen independently while its `Z` coordinate remains fixed.
The expectation of the displayed operator is
`(sqrt3/4)sum_i |M_i(x)|`. Averaging over auxiliary signs and using
(QV3) gives

```math
v(G)\ge\frac1{16}\sum_{q_i\le1/32}\sqrt{q_i}.
```

Finally, `h_2(q)<=2sqrt(q)` bounds the small-marginal contribution by
`32v(G)`. On the remaining rows, `h_2(q)<=log2<1<=32q`, so (QV1)
bounds their contribution by another `32v(G)`. This proves (QV).

For every fixed quadratic phase target, Section 21.3 therefore bounds
the entropy/effective-degree lower-bound output of uniform tensor
compilers by `exp(32)`. Section 29 prevents extending the full velocity
lower bound automatically to cap-bounded `L2` compilers. However, the
flat-diagonal refinement in Section 32 below proves the needed weaker
lower bound and extends this quadratic-phase conclusion to cap-bounded
`L2` compilers as well. The director and counterexample agent have
independently audited the full proof of this section.

## 32. Flat-diagonal velocity controls cap-bounded L2 tensor compilers

**Proved here; the director and counterexample agent independently
audited both derivations. The counterexample agent supplied the shorter
commutator proof below.**
For a scalar unitary `G` on `n` bits define

```math
v_{\rm flat}(G):=
 \sup_{\substack{\rho\succeq0,\ \operatorname{tr}\rho=1\\
                   \rho_{xx}=2^{-n}\ \forall x}}
 \left|\operatorname{tr}\rho(M_G^*J M_G-J)\right|,
 \qquad J=\frac12\sum_{i=1}^n X_i.                    \tag{FV}
```

The supremum is over density matrices with uniform computational
diagonal, not merely pure states and not arbitrary density matrices.
In particular `I(G)<=v_flat(G)<=v(G)`. If a polynomial `Q` on `nk`
bits has degree at most `D`, cap `C=||Q||_infty`, and normalized
`L2` error `delta=||Q-G^(tensor k)||_2`, then

```math
\frac Dk\ge\frac{v_{\rm flat}(G)-n\delta}{C},          \tag{FV1}
```

provided `C>0`; a negative right side is simply vacuous. Consequently,
for each fixed seed `G`, any sequence with `delta_k->0` and
`C_k->1` satisfies

```math
\liminf_{k\to\infty}\frac{D_k}{k}\ge v_{\rm flat}(G).
                                                               \tag{FV2}
```

More generally a cap limit superior at most `C_*>0` gives the lower
bound `v_flat(G)/C_*`. For moving seeds (FV1) remains exact, but the
limiting conclusion requires the displayed `n delta` term to vanish;
one must not silently discard it.

### 32.1 Short direct proof

First consider a single target unitary `U` on `N` bits, a flat-diagonal
density matrix `rho`, and an approximant `Q`. All scalar functions
below are understood as diagonal multiplication operators. Set
`E=U-Q` and `J_N=(1/2)sum_(i<=N)X_i`. The flat diagonal ensures

```math
\operatorname{tr}\rho E^*E
 =\operatorname{tr}\rho EE^*=\|E\|_2^2=\delta^2.
```

State Cauchy--Schwarz gives the two separate bounds

```math
|\operatorname{tr}\rho U^*J_NE|
 \le\|J_N\|\delta=N\delta/2,
\qquad
|\operatorname{tr}\rho U^*EJ_N|
 \le\delta\|J_N\|=N\delta/2.                         \tag{FV3}
```

For the second bound use the state inner product with
`A=E^*U`, `B=J_N`; its squared `A` norm is `delta^2` because `U,E`
are commuting diagonal operators. No commutation of `rho` with
`J_N` is being assumed. The same inner product with `A=J_NU`,
`B=E` proves the first bound.

The operator-valued Bernstein inequality from Section 21.3 yields
`||[J_N,Q]||<=D C`. Hence

```math
|\operatorname{tr}\rho U^*[J_N,U]|
 \le |\operatorname{tr}\rho U^*[J_N,Q]|+N\delta
 \le DC+N\delta.                                      \tag{FV4}
```

For `U=G^(tensor k)` and `rho=rho_0^(tensor k)`, the left-hand
expectation equals `k` times the corresponding seed expectation,
while the diagonal of `rho` is still uniform and `N=nk`. Take the
supremum over `rho_0` to obtain (FV1).

### 32.2 Retained alternative energy-centered derivation

This longer proof was the original route to (FV2); it is retained to
make its independent state-Cauchy--Schwarz checks recoverable. Fix a
flat-diagonal seed state `rho_0` and center
`K=J-tr(rho_0 J)I`. Write
`a=tr(rho_0 G^*KG)`, `tau^2=tr(rho_0 K^2)`, and
`sigma^2=tr(rho_0 G^*K^2G)-a^2`. On tensor products use
`K_tot=sum K_j`, `rho=rho_0^(tensor k)`, and `U=G^(tensor k)`.
Then

```math
\operatorname{tr}\rho K_{\rm tot}=0,\quad
\operatorname{tr}\rho K_{\rm tot}^2=k\tau^2,\quad
\operatorname{tr}\rho U^*K_{\rm tot}U=ka,\quad
\operatorname{tr}\rho U^*K_{\rm tot}^2U
   =k\sigma^2+k^2a^2.
```

Expanding `Q=U+E`, the two cross errors in the target energy total
at most `2delta sqrt(k sigma^2+k^2 a^2)`, and the quadratic error is
at most `k||K||delta^2`. The identity
`Q^*K_tot Q=|Q|^2 K_tot+Q^*[K_tot,Q]`, together with the flat
diagonal and `|| |Q|^2-1 ||_2<=(C+1)delta`, gives

```math
k|a|\le DC^2+(C+1)\delta\sqrt{k}\tau
       +2\delta\sqrt{k\sigma^2+k^2a^2}
       +k\|K\|\delta^2.                              \tag{FV5}
```

For fixed seed, `delta->0,C->1`, division by `k` also yields (FV2).
The shorter bound (FV1) has better cap dependence and is canonical.

### 32.3 Quadratic phases: the witnesses already have flat diagonal

The proof of Section 31 in fact establishes the stronger theorem

```math
H(|\widehat {e^{itP}}|^2)\le64\,v_{\rm flat}(e^{itP})
\qquad (P\text{ real},\ \deg P\le2).                 \tag{FQ}
```

The constant-vector state used to witness `I(G)` has uniform
computational diagonal. For the other witness, average Section 31's
product states over all independent auxiliary signs `x`. Each local
state has computational probabilities `(1+x_i/2)/2` and
`(1-x_i/2)/2`; the chosen `Y_i` sign does not affect those probabilities.
Thus for every computational string `z`, the average diagonal entry is

```math
\mathbb E_x\prod_i\frac{1+z_i x_i/2}{2}=2^{-n}.
```

Every row expectation was already nonnegative before averaging, and
the lower bound used precisely this average. It therefore witnesses
`v_flat`, not merely `v`. Replacing both velocity occurrences by
`v_flat` in the entropy argument proves (FQ). The same observation
upgrades the director's hollow-quadratic constant to `152/3` with
`v_flat` in place of `v`.

Combining (FQ), (FV2), and the counterexample artifact's general
entropy-transfer criterion proves: **for every fixed real quadratic
phase seed, every cap-`1+o(1)`, `L2`-convergent tensor compiler has
entropy/effective-degree output at most `exp(32)`**. This covers
arbitrary times, arbitrary quadratic and linear coefficients, and
arbitrary fixed seed dimension. It is not a bound for arbitrary
Walsh polynomials, and it is not a moving-seed theorem without control
of the finite error in (FV1).

### 32.4 Exact checks distinguishing the two velocities

For Boolean `b`, the operator `M_b J M_b-J` is minus the adjacency
matrix of its sensitivity graph: it has entry `-1` on a sensitive
edge and zero otherwise. Flat diagonal positivity implies
`|rho_(x,y)|<=2^-n`; summing its sensitive edges bounds every state
expectation by `I(b)`. The constant-vector state attains `-I(b)`.
Therefore

```math
v_{\rm flat}(b)=I(b)\qquad\text{for every Boolean }b.  \tag{FB}
```

In particular the rare-AND seed of Section 29 has
`v_flat=n 2^(1-n)`, while `v=sqrt n`. Its cap-bounded `L2` compiler
therefore does not contradict (FV2).

For independent-coordinate phases
`G=exp(i sum_j theta_j x_j)`, each one-qubit difference has eigenvalues
`+-|sin theta_j|` and an equatorial eigenvector. Aligned tensor
eigenvectors have uniform computational diagonal, giving
`v_flat(G)=v(G)=sum_j|sin theta_j|` exactly. Thus the simpler bound
`H<=2v_flat` now closes this class for cap-bounded `L2` compilers as
well, independently of the counterexample artifact's complex
Bernstein--Walsh argument.

### 32.5 Same-ambient compilers have bounded actual BH ratios

The director observed a stronger consequence, independently checked
here. Fix a nonconstant quadratic phase `G` on `n` bits. Suppose
`Q_k` is a polynomial on **exactly those `nk` tensor input bits**, with
`||Q_k-G^(tensor k)||_2=delta_k->0` and `||Q_k||_infty->1`.
Let `D_k` denote its actual Walsh degree. Then

```math
\limsup_{k\to\infty}
 \frac{\|\widehat Q_k\|_{q_{D_k}}}{\|Q_k\|_\infty}
 \le e^{32}.                                          \tag{FQB}
```

Indeed `V=v_flat(G)>=I(G)>0`. Equation (FV1) gives
`liminf D_k/k>=V`, so eventually `D_k>=Vk/2`. The coefficient-space
dimension is exactly `2^(nk)`, and finite-dimensional norm comparison
therefore gives

```math
\|\widehat Q_k-\widehat {G^{\otimes k}}\|_{q_{D_k}}
 \le 2^{nk/(2D_k)}\delta_k
 \le2^{n/V}\delta_k\longrightarrow0.                  \tag{FQB1}
```

Let `p=|Ghat|^2`. Since the seed is fixed and has finite support, the
elementary Taylor expansion at exponent one of its Renyi power sum gives

```math
\log\|\widehat {G^{\otimes k}}\|_{q_D}
 =\frac{k}{q_D}\log\sum_S p_S^{D/(D+1)}
 =\frac{kH(p)}{2D}+O_G(k/D^2).                         \tag{FQB2}
```

Zero entries of `p` are omitted; the finitely many positive entries
make the remainder uniform as `D->infinity`. Along every subsequence
where `D_k/k->d`, one has `d in [V,n]`, and (FQB1)--(FQB2) show that
the actual coefficient ratio tends to `exp(H(p)/(2d))`. The bound
`H(p)<=64V` proves (FQB) by compactness of `[V,n]`.

This is an upper bound for the actual ratios of this precise compiler
class, not only for the entropy-transfer lower-bound expression.
The fixed-seed and exact-ambient assumptions are essential to the
argument: arbitrary extra variables can invalidate (FQB1), and moving
seeds need not have a uniform Renyi remainder or dimension factor.
No upper bound for unrestricted BH constants is asserted. The same
argument applies to any fixed nonconstant seed satisfying
`H(G)<=C v_flat(G)`, with output bound `exp(C/2)`.

## 33. Final reciprocal audit and replay record

At 06:47 UTC the following exact scripts were rerun successfully with
the repository `.venv/bin/python`:

- `computations/bh_boundedness_barrier_2026_09_18.py`: 7,048 finite
  identities, all 254 proper-subset checks and all 210 root-of-unity
  checks passed. Floating CLT and damping outputs remain diagnostics,
  not substitutes for the analytic proofs.
- `computations/bh_boundedness_barrier_2026_09_18_tangent.py`: the exact
  Eisenstein coefficients, 36 new terms, squared mass `45/4`, radial
  idempotence and absorption, and degree-eight normalization identities
  all passed.
- `computations/bh_boundedness_barrier_2026_09_18_velocity_converse.py`:
  the exact three-bit characteristic polynomial and parity restriction
  passed.
- `computations/bh_boundedness_barrier_2026_09_18_velocity_exact.py`:
  the degree-16 rational characteristic polynomial and exact root
  isolations passed; the largest positive seed velocity is simple and
  lies strictly between `1.865` and `1.867`.

The positive researcher independently reread the actual primary
de Wolf Appendix A and verified Section 29's exact-query count and
degree conversion. Both other research tracks and the director have
independently reconstructed Sections 31--32, including the averaged
state's flat diagonal and the same-ambient actual-ratio conclusion.

Two further checks were completed against the other canonical files:

1. The director's Section 5 sphere-valued junta extension is valid for
   each fixed real target dimension `r`: the acute-cap coloring uses
   real coefficient-vector dot products, its top-`2m` coefficient has
   no omitted terms, and a chosen nonzero component supplies the
   restriction-survival probability `2^(1-m)`. Its Ramsey bound is
   qualitative, not exponential in `m`.
2. The counterexample artifact's Section 18 rare-marker identity is
   exact. For disjoint `E in {0,1}` and unitary `V`, with
   `p=E E<=1/4`, `s=E|V-1|^2`, the spectral law of
   `1+E(V-1)` equals `(1-ps)delta_empty+ps(mu tensor nu)` including
   the empty coefficient. The resulting bound
   `H<=pH(V)+16pr` and the declared-budget weighted-maximum estimate
   have correct constants. The latter is not an optimal-degree lower
   bound. The proof uses the declared search/evaluation budget and
   must not silently replace a query cost by bare Walsh degree.

The positive artifact's Section 13 one-affine-block stability theorem
was also independently checked: squared-modulus defect at most
`eta<=1/4` implies distance at most `9sqrt(eta)` from an exact affine
unitary, uniformly in dimension. Its proof controls discarded real
coefficient mass by `2eta`, discarded imaginary mass by `3sqrt(eta)`,
then normalizes the two retained orthogonal coefficients. The example
`1+i t n^(-1)sum x_i` shows the square-root rate cannot be improved
uniformly. No multi-block stability extension is asserted: slice-wise
choices of the retained coordinates can vary with the other blocks.

All these are internal independent checks, not assertions of external
novelty. The campaign's unrestricted boundedness/divergence question
remains unresolved by this artifact.

## 34. Bounded primary-literature scope audit of the final certificate

Requested by the director after the mathematical freeze; searched
2026-09-18 from 06:52 UTC. This is a limited identification of close
antecedents, **not a novelty search sufficient to establish priority**.
Primary texts were opened directly. Search snippets and secondary
summaries were not used as theorem evidence.

### 34.1 Exact antecedent for the commutator-versus-degree ingredient

Aaronson--Ben-David--Kothari--Rao--Tal,
[Degree vs. Approximate Degree and Quantum Implications of Huang's
Sensitivity Theorem](https://arxiv.org/html/2010.12629), Section 4.3,
Theorem 18 and Lemmas 20--24, are the closest exact match. They use
`B_q=(diag(q) A_cube-A_cube diag(q))/2`, precisely our commutator up
to sign. Their Schur-multiplier argument bounds this matrix by the
polynomial degree times its cap. Theorem 18 gives
`lambda(f)<=approxdeg_epsilon(f)/(1-2epsilon)` for total Boolean
outputs and pointwise sign approximation. Lemma 11 identifies the
same conjugated cube-adjacency difference for Boolean outputs.

Thus neither the commutator-degree estimate nor the Boolean identity
`v=lambda` should be presented as externally new. Section 32's
flat-diagonal `L2` certificate is a short corollary of this established
ingredient and state Cauchy--Schwarz. I did not locate its literal
fixed-diagonal formulation in that text. Its uniform-distribution
error condition differs from Theorem 18's pointwise condition.

### 34.2 Distributional adversary: close proof architecture, different target

Belovs--Brassard--Hoyer--Kaplan--Laplante--Salvail,
[Provably secure key establishment against quantum
adversaries](https://arxiv.org/html/1704.08182), Section 4,
Definition 7 and Theorem 4 in the arXiv text, prove a lower bound for
quantum algorithms distinguishing input laws `mu,nu`. Its numerator
is `delta_mu^* Gamma delta_nu-tau(s_mu,s_nu)||Gamma||`, with
`delta_mu(x)=sqrt(mu(x))`; its denominator controls
`Gamma` restricted to pairs differing in a queried coordinate.
The acceptance probabilities may be distributional, and the laws
may overlap. The final proceedings version numbers the result
Theorem 11.

This is a genuine average-case adversary antecedent, but its theorem
is about query algorithms and acceptance probabilities, not every
bounded scalar polynomial with small `L2` phase error. No automatic
conversion of arbitrary degree-`D` polynomials into `O(D)`-query
algorithms is available here. The same result is restated as
Theorem 5 in Belovs,
[A Direct Reduction from the Polynomial to the Adversary
Method](https://arxiv.org/html/2301.10317); that paper converts dual
polynomial lower bounds into adversary bounds, not arbitrary
low-degree approximants into algorithms.

### 34.3 State-conversion adversary norms are not our error condition

Lee--Mittal--Reichardt--Spalek--Szegedy,
[Quantum query complexity of state
conversion](https://arxiv.org/html/1011.3020), Lemma 4.7, shows that
one query changes the filtered `gamma_2` distance of input-indexed
Gram matrices by at most two. Theorem 4.9 combines this with
`gamma_2` output-approximation constraints to give coherent and
non-coherent state-conversion bounds. Its error parameters on the
upper and lower estimates differ, a qualification stated in the paper.

The Gram matrices have prescribed diagonals, so this is conceptually
near the density-matrix device. However, these Gram-matrix output
conditions are not the normalized uniform-cube `L2` error of a scalar
multiplier. No exact reduction identifying its SDP with `v_flat`
was established in this bounded audit.

### 34.4 Spectral sensitivity and fractional sensitivity must be separated

Ambainis--Prusis--Vihrovs,
[On Block Sensitivity and Fractional Block
Sensitivity](https://arxiv.org/html/1810.02393), Section 2, defines
fractional block sensitivity by a packing LP on sensitive blocks at
an input, followed by a maximum over inputs. This is not the
flat-diagonal quantity. Our proved Boolean specialization is
`v_flat(b)=I(b)`, the uniform average sensitivity, whereas
`v(b)=lambda(b)` is spectral sensitivity.

An elementary check makes the distinction strict: for the sign of
`AND_n`, `fbs=n`, `lambda=sqrt(n)`, and
`v_flat=I=n 2^(1-n)`. Consequently calling `v_flat` simply
"fractional sensitivity" would misidentify a standard notion.

### 34.5 Closest quantum Fourier-entropy result has a logarithmic loss

Bu--Garcia--Jaffe--Koh--Li,
[Complexity of quantum circuits via sensitivity, magic, and
coherence](https://arxiv.org/html/2204.12051), Theorem 23, proves
for normalized operators on `n` qudits a weak quantum FEI estimate
`H(O)<=2(log n+log d)I(O)+h(P_O(0))`, with base-two logarithms.
Their Definition 6 instead defines circuit sensitivity as
`max_(||O||_2=1)|I(UOU^*)-I(O)|`. It is not `v_flat`; for example
their stable class includes products of single-qudit unitaries,
whereas our independent nonconstant scalar phases have positive
`v_flat`.

Their dimension-free quantum FEI statement is explicitly a conjecture
for Hermitian involutions, not a theorem for arbitrary-angle scalar
quadratic phases. Its noncommutative failure in Slote's primary
[arXiv:2608.01424](https://arxiv.org/html/2608.01424) also does not
settle our diagonal phase class. I found no `H<=C v_flat` estimate
for that class in the opened Bu text.

### 34.6 The phase distribution is exactly an IQP output distribution

Fujii--Morimae,
[Commuting Quantum Circuits and Complexity of Ising Partition
Functions](https://arxiv.org/html/1311.2128), Definition 3, uses
commuting gates `exp(i theta_j product_(i in S_j) Z_i)`, input
`|+>^n`, and `X`-basis measurement. With `|S_j|<=2`, its output
probabilities are exactly `|widehat(exp(iP))(S)|^2` for our real
quadratic `P`. This identification is standard, not a new quantum
model.

Park--Kastoryano,
[Complexity phase transitions in instantaneous quantum polynomial-time
circuits](https://arxiv.org/html/2204.08898v2), equation (1), studies
precisely single-`Z` and pair-`ZZ` phases on random graphs. Its
opened model and entropy passages concern simulation,
anticoncentration, learnability, and entanglement-spectrum behavior;
I did not find a deterministic, all-angle Shannon-output-entropy
bound by the flat-diagonal commutator certificate there. Entanglement
entropy and the Fourier/output Shannon entropy cannot be interchanged.

### 34.7 Limited-search conclusion

Queries covered distributional and average-case adversary bounds,
state conversion, spectral and fractional sensitivity, quantum FEI,
phase-state Fourier entropy, and IQP/Ising output entropy. Exact
phrase searches for "flat diagonal"/"uniform diagonal" with
commutator or Fourier entropy did not identify the literal theorem.
Hits on quadratic Boolean functions use algebraic degree over `F_2`,
not arbitrary real angles in `exp(iP)`; hits on continuous
quadratic-phase Fourier transforms concern a different domain and
usually entropy uncertainty lower bounds.

The justified conclusion is narrow: **an established commutator-degree
ingredient was located, and several close but nonidentical frameworks
were checked; no exact antecedent for the combined flat-diagonal
`L2` certificate plus arbitrary-angle quadratic entropy theorem was
identified in this bounded search.** This is not evidence sufficient
to claim external novelty. No campaign theorem is strengthened on the
basis of this literature audit.

## 35. Closing synthesis audit

At 07:02 UTC I independently reread the director's complete closing
synthesis. The field distinction, critical exponents, strict finite
and asymptotic constant lower bounds, and approximation hypotheses
match the canonical proofs. In particular, the actual-ratio bound
`exp(32)` is restricted to a fixed nonconstant quadratic-phase seed,
exactly its tensor ambient space, cap tending to one, and `L2` error
tending to zero. It does not silently assert growing-seed uniformity
or permit uncontrolled extra variables. The independent-phase
specialization `exp(1)` has the same fixed-seed scope.

Two requested wording corrections were incorporated: the finite
flat-diagonal certificate explicitly requires `C>0`, and the
inheritance obstruction excludes only nonconstant asymmetric pairs.
The symmetric `(4,4)` space spanned by `P` survives but yields zero
antisymmetric transfer. The trivial constant pair also survives.

The closing replay record starts at 06:50:58 UTC and records thirteen
programs, all with exit code zero, over approximately 24 seconds.
I read the runner and its stored record separately. Its labels
distinguish exact algebraic/rational certificates from numerical
diagnostics; successful diagnostics are not promoted to proofs.
The record additionally preserves hashes for twenty-six campaign
scripts. All twenty-six hashes still matched at this final audit.

No boundedness or divergence conclusion for unrestricted `B_m` has
been obtained. The final synthesis states this explicitly and does
not turn the limited literature search into a novelty claim.
