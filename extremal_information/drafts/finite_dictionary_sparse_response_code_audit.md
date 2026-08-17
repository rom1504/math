# Finite-dictionary sparse response codes: theorem and audit

Status: **rigorous audited draft**.  This note generalizes and independently
checks the sparse-mask argument in
`quadratic_dense_bridge_compression_ceiling.md`.  It does not alter the global
theory documents.

## 1. Public feature dictionary and response semantics

Let `X` be a nonempty finite set of size `L`, let `m>=1`, and fix a public
dictionary

```math
\phi_e:X\longrightarrow[-1,1],\qquad 1\le e\le m.
```

For `a in {-1,1}^m`, define

```math
H_a(x)=\sum_{e=1}^m a_e\phi_e(x),
\qquad
V_\Phi=\max_{x\in X}\sum_{e=1}^m\phi_e(x)^2\le m.       \tag{FD.1}
```

Here **universal** means independent of the input coefficient vector `a` and
of all later queries.  The mask family is allowed to depend on the public
dictionary `Phi`.  No simultaneous claim over every admissible dictionary is
made below.

A max-type future is a kernel `K:X times Y -> R union {-infinity}` for which
the displayed maxima below are finite.  Its response is

```math
(T_KH)(y)=\max_{x\in X}\{H(x)+K(x,y)\}.                  \tag{FD.2}
```

The elementary contraction

```math
\|T_KH-T_KG\|_\infty\le \|H-G\|_\infty                 \tag{FD.3}
```

holds for every such `K`.  It also holds after any chain of shared max-plus
transitions.  Thus a uniform approximation to the terminal landscape answers
every declared future simultaneously, with no depth accumulation.  This
scope does not include a construction that duplicates the approximated
landscape several times along one trajectory; then the obvious multiplicity
factor is necessary.

## 2. Exact parameter theorem

### Theorem FD.1 (universal Bernoulli-mask response code)

Fix `0<p<=1/2`, put `q=1-p`, and let `E>0`.  Define

```math
\Delta=
2L\exp\left\{-{E^2\over2(pV_\Phi/q+E/3)}\right\}
+\exp(-pm/8).                                             \tag{FD.4}
```

Suppose `Delta<=delta<1`, and set

```math
s=\left\lfloor(1-p/2)m\right\rfloor,
\qquad
N=\left\lceil{(m+1)\log2\over\log(1/\delta)}\right\rceil.\tag{FD.5}
```

There are at most `N` public masks `S_1,...,S_N subseteq {1,...,m}`, each of
cardinality at most `s`, such that, for every `a in {-1,1}^m`, some `j`
satisfies

```math
\sup_{x\in X}
\left|
H_a(x)-{1\over q}\sum_{e\in S_j}a_e\phi_e(x)
\right|\le E.                                             \tag{FD.6}
```

Consequently the family has an absolute-response code of cardinality at most

```math
N2^s                                                       \tag{FD.7}
```

and fixed-length cost at most

```math
s+\lceil\log_2N\rceil                                    \tag{FD.8}
```

bits.  The decoder stores the mask index and the signs on that mask.  Given
any future `K`, it outputs

```math
\max_x\left\{
{1\over q}\sum_{e\in S_j}a_e\phi_e(x)+K(x,y)
\right\},                                                  \tag{FD.9}
```

which differs from the true response uniformly in `y` by at most `E`.

#### Proof

Draw one mask by independent variables `Z_e~Bernoulli(q)` and let
`S={e:Z_e=1}`.  For fixed `a` and `x`, the landscape error is a sum of
independent variables

```math
X_e=a_e\phi_e(x)(1-Z_e/q).                                \tag{FD.10}
```

They are centered.  Since `p<=q`,

```math
|X_e|\le1,
\qquad
\mathbb E X_e^2={p\over q}\phi_e(x)^2.                  \tag{FD.11}
```

Bernstein's inequality therefore gives

```math
Pr\left\{\left|\sum_eX_e\right|>E\right\}
\le
2\exp\left\{-{E^2\over2(pV_\Phi/q+E/3)}\right\}.       \tag{FD.12}
```

A union bound over the `L` rows of the dictionary gives the first term in
`Delta`.  This is the only union taken at the concentration stage.

Let `D=sum_e(1-Z_e)~Binomial(m,p)` be the erased-coordinate count.  The
standard lower-tail Chernoff bound gives

```math
Pr\{D<pm/2\}\le\exp(-pm/8).                               \tag{FD.13}
```

On the complementary event, `|S|=m-D<=s`.  Hence, for each fixed `a`, one
mask is inaccurate or too large with probability at most `Delta<=delta`.

Now sample `N` masks independently.  The probability that a fixed `a` has no
good mask is at most `delta^N`.  Independence between different coefficient
vectors is neither asserted nor needed.  By linearity of expectation, the
expected number of uncovered sign vectors is at most

```math
2^m\delta^N
\le 2^m\exp\{-(m+1)\log2\}=1/2.                          \tag{FD.14}
```

Thus some deterministic list covers all `2^m` coefficient signings.  Delete
from it every oversized mask; none of those masks was used in the covering
event just counted.  For a retained mask, its index and the at most `s`
coefficient signs specify the surrogate, proving (FD.6)--(FD.8).  Finally,
(FD.9) follows from the max-plus contraction (FD.3). `square`

### Corollary FD.2 (a convenient explicit tradeoff)

Put `t=log(8L)` and suppose `V_Phi>0`.  Choose

```math
p=\min\left\{{1\over2},{E^2\over8V_\Phi t}\right\}.      \tag{FD.15}
```

If

```math
E\ge {4t\over3},
\qquad
pm\ge8\log4,                                               \tag{FD.16}
```

then `Delta<=1/2`, so `N=m+1` is valid and

```math
b\le
m-\min\left\{{m\over4},{mE^2\over16V_\Phi\log(8L)}\right\}
+\lceil\log_2(m+1)\rceil.                                \tag{FD.17}
```

Indeed, (FD.15) implies `pV_Phi/q<=E^2/(4t)`, while the first condition in
(FD.16) implies `E/3<=E^2/(4t)`.  The Bernstein term in (FD.4) is therefore
at most `1/4`; the second condition makes the Chernoff term at most `1/4`.
The number of erased sign bits is at least `pm/2`, which is exactly the
saving displayed in (FD.17).  When `V_Phi=0`, all landscapes vanish and a
zero-bit code suffices.

The exact theorem is often sharper than this deliberately simple corollary.
In particular, (FD.16) is only a sufficient large-deviation regime, not a
claim about the optimal response rate-distortion function.

## 3. Audit of the quadratic specialization

For quadratic Boolean landscapes, take

```math
X=\{-1,1\}^n,
\quad L=2^n,
\quad m={n\choose2},
\quad \phi_{ij}(x)=x_ix_j,
\quad V_\Phi=m.                                           \tag{FD.18}
```

Set `E=epsilon n^(3/2)`, `p=epsilon^2/2`, and `q=1-p`, with
`0<epsilon<=1` and `n>=64/epsilon^2`.  Then `p<=q` and

```math
{p\over q}={\epsilon^2\over2-\epsilon^2}\le\epsilon^2,
\qquad
{pV_\Phi\over q}\le{\epsilon^2n^2\over2}.               \tag{FD.19}
```

The size assumption says `epsilon sqrt(n)>=8`, whence

```math
E/3\le\epsilon^2n^2/24.                                   \tag{FD.20}
```

Therefore the Bernstein exponent is at least `12n/13`, and

```math
2L e^{-12n/13}
=2e^{-(12/13-\log2)n}<1/4.                               \tag{FD.21}
```

Also `m>=n^2/3` and

```math
pm/8\ge\epsilon^2n^2/48\ge4096/48>\log4,                 \tag{FD.22}
```

so the mask-size failure is below `1/4`.  Theorem FD.1 with `delta=1/2`
therefore gives `N=m+1` and

```math
b\le
\left\lfloor(1-\epsilon^2/4)m\right\rfloor
+\lceil\log_2(m+1)\rceil.                                \tag{FD.23}
```

The source draft used `m+2` masks and a harmless ceiling in the retained-bit
bound.  Those choices are conservative; every probability estimate and its
claimed `Theta(n^2)` rate is valid.  The response statement for an arbitrary
bridge is precisely (FD.3), with `K(x,y)=x^TBy`.

The quantifiers deserve emphasis:

1. union over `X` occurs in (FD.12);
2. the list argument (FD.14), not a second concentration estimate for one
   mask, covers all `2^m` signings;
3. the resulting masks are independent of the signing and of the future;
4. the list may depend on the fixed feature dictionary;
5. the sparse centers have coefficients in `{0,+-1/q}` and need not belong
   to the original sign family.

## 4. Two further model classes

### 4.1 Boolean Littlewood polynomials and bounded CSP dictionaries

Let `Ecal` be any public collection of subsets of `[n]`, and put

```math
\phi_e(x)=\prod_{i\in e}x_i,
\qquad e\in\mathcal E.                                   \tag{FD.24}
```

This includes every fixed-degree Boolean Littlewood polynomial, with
`m=|Ecal|`, `L=2^n`, and `V_Phi=m`.  Theorem FD.1 sparsifies its signed
Fourier support while preserving the polynomial uniformly on the Boolean
cube and, therefore, preserving every later max-type external field.  The
same statement applies verbatim to finite-alphabet CSPs after scaling each
public clause payoff into `[-1,1]`; there `X` is the assignment space and
the `phi_e` are clause-payoff tables.

This is an ambient response code, not a claim that the sparsified polynomial
or CSP remains in a prescribed regular subclass.

### 4.2 Code/coset correlation landscapes

Let `C subseteq {-1,1}^m` be a public binary code and take

```math
X=C,
\qquad
\phi_e(c)=c_e.
```

For a received sign word `a`,

```math
H_a(c)=a\mathbin\cdot c=m-2d_H(a,c).                      \tag{FD.25}
```

Thus (FD.6) supplies a sparse reweighted coordinate sketch that approximates
all codeword correlations uniformly.  In particular it approximates the
nearest-code distance within `E/2`, and it continues to do so after any
shared state-dependent reward or penalty `K(c,y)`.  Here `L=|C|` and
`V_Phi=m`; the exact tradeoff (FD.4) records explicitly how the code size
enters.  No linearity of `C` is required.

These two applications are mathematically distinct: one dictionary is a
system of Boolean monomials indexed by configurations, while the other is a
codeword-coordinate incidence matrix.  Both inherit the same all-future
guarantee solely because uniform landscape error is a max-plus invariant.

## 5. Scope and limitations

1. **Existential public list.**  The proof gives a list of `O(m)` public
   masks when `Delta<=1/2`, but not an efficient algorithm for finding the
   successful mask for a given input.  Public-list storage is not charged to
   the per-landscape information rate.
2. **Weighted centers.**  Importance reweighting by `1/q` is essential.
   Requiring decoded centers to remain in `{-1,1}^m` is a different covering
   problem.
3. **One-shot ambient compression.**  The theorem answers all shared
   max-type futures but does not assert closure of the decoded family under a
   model's internal composition law.
4. **Not an optimality theorem.**  Bernstein thinning proves an upper bound.
   It gives no response-packing lower bound and no reason that `Theta(m)`
   bits are necessary.
5. **Finite public row set.**  Finiteness enters only through the union bound.
   An infinite state space requires a metric-entropy or regularity argument
   controlling the feature rows.

## 6. Audit verdict

**PROMOTE, with two minor clarifications.**  The sparse-surrogate theorem in
`quadratic_dense_bridge_compression_ceiling.md` is correct, including its
Bernstein constants, Chernoff estimate, coverage of all coefficient
signings, decoder, and arbitrary-bridge response guarantee.  Its `m+2`
masks can be replaced by `m+1`, and “universal” should explicitly mean
universal over inputs and future queries for the fixed public quadratic
dictionary.  The finite-dictionary theorem above is the clean reusable
form.
