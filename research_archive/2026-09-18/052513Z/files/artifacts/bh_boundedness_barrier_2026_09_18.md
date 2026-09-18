# Boolean BH boundedness: weighted-bootstrap barrier diagnosis

2026-09-18, bounded campaign 04:27:50--07:27:50 UTC. Independent barrier
track. This file is a working research artifact, not a claim that the BH
boundedness problem has been solved. Mathematical statements labeled
**proved here** have a complete derivation below, but await a separate
root-agent verification before promotion into the campaign ledger.

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
`eta=t^t/2^(t^2)`. If `f` is the resulting Boolean function and `x_1`+is its default dictator, then

```math
\|f-x_1\|_2^2\le4\eta.
```

Thus its Fourier mass outside the dictator coefficient is at most
`4eta`. A Fourier coefficient must involve exactly one data bit and
an arbitrary subset of the `t^2` address bits, so there are at most
`t^t 2^(t^2)` possible coefficients. The elementary entropy bound

```math
H(\widehat f^{,2})
\le h_2(4\eta)+4\eta\log(t^t2^{t^2})
\longrightarrow0
```

holds for large `t`, with `h_2` the binary entropy in natural units.
This particular dimension-saving construction therefore does not supply
the needed large Fourier entropy. This observation does not rule out
other approximation-theoretic mechanisms or prove anything about the
BH norms of arbitrary approximating polynomials themselves.
