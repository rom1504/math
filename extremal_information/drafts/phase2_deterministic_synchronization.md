# Phase 2 specialist report: deterministic synchronization

## Executive conclusion

There is a clean deterministic core inside multi-species spin-glass
synchronization, but ultrametricity by itself is not that core.

This report proves the following project-level theorem.

> If every species overlap and every sum of two distinct species overlaps is
> approximately ultrametric, then overlap profiles attached to two pairs
> sharing a state cannot cross by more than the ultrametric defect.  If, in
> addition, the graph of state-pairs is connected across total-overlap levels
> by paths with controlled backtracking, every species overlap is uniformly
> close to a nondecreasing Lipschitz function of the total overlap.  The same
> error controls arbitrary zero-temperature responses with Lipschitz coupling
> energy.

The second hypothesis is a finite, deterministic replacement for the
root-density step in Panchenko's synchronization proof.  It only refers to
the scalar total-overlap geometry; it does not assume the functions to be
recovered.

A sharp counterexample shows why such a cross-root hypothesis is necessary.
There are two positive-semidefinite species kernels for which:

- every nonnegative mixture is exactly ultrametric;
- i.i.d. replicas give a weakly exchangeable Gram array;
- even the largest conditional variance of a species overlap given total
  overlap tends to zero;

yet no function of total overlap approximates that species uniformly within
less than a fixed constant.  A zero-temperature conditioned query detects the
gap.  Thus average synchronization, approximate rank-one covariance, and
ultrametricity can all miss one exposed zero-entropy fiber.

The proved theorem is useful beyond vocabulary: it identifies two logically
different obligations for deterministic Parisi-like compression:

1. **local no-crossing**, supplied by ultrametric closure of mixtures; and
2. **cross-root calibration**, supplied in random spin glasses by strong
   replica identities and here by a finite linkage condition.

The result does not by itself compress a landscape.  Compression occurs only
when the total overlap has a small state space and the linkage constants are
controlled.

## 1. Imported mechanism from spin glasses

This subsection is **imported**, not a new result.

For a finite species set `S`, positive weights `lambda_s`, and limiting
species overlap arrays `R^s`, Panchenko proves multi-species
Ghirlanda--Guerra identities for every bounded measurable function of the
whole species-overlap vector.  He then proves that there are nondecreasing
`1/lambda_s`-Lipschitz functions `L_s` such that

```math
R^s_{\ell,\ell'}=L_s\left(\sum_t\lambda_tR^t_{\ell,\ell'}\right)
\quad\text{almost surely}.                           \tag{1.1}
```

The primary source is Theorems 3 and 4 of Panchenko,
[The free energy in a multi-species Sherrington--Kirkpatrick
model](https://doi.org/10.1214/14-AOP967), *Annals of Probability* 43
(2015), 3494--3513; the [arXiv version](https://arxiv.org/abs/1310.6679)
has the same numbering.  The ultrametricity input is Panchenko,
[The Parisi ultrametricity
conjecture](https://doi.org/10.4007/annals.2013.177.1.8), *Annals of
Mathematics* 177 (2013), 383--393.

Two different mechanisms occur in the proof.

1. Each `R^s`, and each `R^s+R^t`, is a positive-semidefinite array satisfying
   the ordinary Ghirlanda--Guerra identities.  Ultrametricity of these three
   arrays implies Panchenko's no-crossing lemma:

   ```math
   R^s_{12}>R^s_{13}
   \quad\Longrightarrow\quad
   R^t_{12}\ge R^t_{13}\quad\text{for every }t.       \tag{1.2}
   ```

2. Ghirlanda--Guerra identities also imply that, from every sampled root, the
   observed overlaps are dense in the relevant support.  This calibrates the
   coordinatewise order across different roots and yields the common
   functions in (1.1).

Only the first step is a local three-state statement.  The counterexample in
Section 5 proves that even ultrametricity of *all* nonnegative mixtures does
not replace the second step.

Also note the operational distinction: (1.1) is almost sure for the limiting
Gibbs experiment.  It does not automatically give a uniform guarantee over
states that have negligible Gibbs mass but may be exposed by a new
zero-temperature query.

## 2. Finite definitions

Let `Omega` be a finite state set with `|Omega|>=2` and

```math
E=\binom{\Omega}{2}
```

the set of unordered pairs of distinct states.  Let `S` be a finite species
set.  For each `s in S`, let

```math
R_s:E\longrightarrow[0,1]
```

be a symmetric similarity kernel.  Positive semidefiniteness is natural in
the overlap application, but the theorem below needs only the metric
inequalities.  Fix weights `lambda_s>0`, `sum_s lambda_s=1`, and put

```math
q(e)=\sum_s\lambda_sR_s(e).                          \tag{2.1}
```

### Approximate ultrametricity

A kernel `K` is `eta`-ultrametric if for all distinct `x,y,z`,

```math
K(y,z)\ge\min\{K(x,y),K(x,z)\}-\eta.                \tag{2.2}
```

The normalization of `K` is part of the definition.  Below the same absolute
defect `eta` is imposed on kernels in `[0,1]` and sums of two distinct species
in `[0,2]`.  One can instead retain separate defects; the proof constant
`3 eta` is then replaced by the sum of the three relevant defects.

### The pair line graph and monotone linkage

Two elements of `E` are adjacent, written `e sim e'`, if the two unordered
pairs share a state.  This is the line graph of the complete graph on
`Omega`.

For an oriented path `e_0,...,e_m`, define its total downward variation in
the scalar order parameter by

```math
B_-(e_0,\ldots,e_m)
=\sum_{j=0}^{m-1}\bigl(q(e_j)-q(e_{j+1})\bigr)_+.    \tag{2.3}
```

Call `q` **`(D,tau)`-monotone-linked** if, whenever `q(e)<=q(e')`, there is
an adjacent path from `e` to `e'` of length at most `D` with
`B_-<=tau`.

This condition is finite and checkable from `q` alone.  It says neither that
the `R_s` are functions of `q` nor that their values agree.  When `tau=0`, it
asks for a nondecreasing path in total overlap.  For pairs at the same total
overlap such a path stays inside that exact fiber.

### The uniform cancellation defect

For two pair profiles define

```math
\mathfrak c(e,f)=
\sum_s\lambda_s|R_s(e)-R_s(f)|-|q(e)-q(f)|.         \tag{2.4}
```

This is nonnegative by the triangle inequality.  If

```math
P=\sum_s\lambda_s(R_s(e)-R_s(f))_+,
\qquad
N=\sum_s\lambda_s(R_s(f)-R_s(e))_+,
```

then

```math
\mathfrak c(e,f)=P+N-|P-N|=2\min(P,N).              \tag{2.5}
```

Thus `mathfrak c/2` is exactly the smaller of the total upward and downward
species movements.  It vanishes precisely when the two species profiles are
coordinatewise comparable.

### Proposition 2.1 (cancellation defect implies synchronization)

Suppose

```math
\sup_{e,f\in E}\mathfrak c(e,f)\le\delta.           \tag{2.6}
```

Then for every species `s` there is a nondecreasing
`1/lambda_s`-Lipschitz function `L_s` on the total-overlap interval such that

```math
\max_{e\in E}|R_s(e)-L_s(q(e))|
\le {\delta\over2\lambda_s}.                        \tag{2.7}
```

#### Proof

For arbitrary `e,f`, put `Delta=q(e)-q(f)`.  If
`R_s(e)>R_s(f)`, its weighted increment is at most `P`.  When `Delta>=0`,
equations (2.5)--(2.6) give `P=Delta+N<=Delta+delta/2`; when
`Delta<0`, they give `P<=delta/2`.  Consequently

```math
R_s(e)-R_s(f)
\le {\bigl(q(e)-q(f)\bigr)_++\delta/2\over\lambda_s}.
                                                               \tag{2.8}
```

Define

```math
\widetilde L_s(p)=
\inf_{f\in E}
\left\{R_s(f)+{(p-q(f))_+\over\lambda_s}\right\}. \tag{2.9}
```

Each function inside the infimum is nondecreasing and
`1/lambda_s`-Lipschitz.  Their infimum has both properties: compare the same
term at two arguments for monotonicity, and use an arbitrarily close
minimizer at either argument for the two Lipschitz inequalities.  At
`p=q(e)`, the term `f=e` gives `widetilde L_s(p)<=R_s(e)`, while (2.8) gives

```math
R_s(e)-{\delta\over2\lambda_s}
\le\widetilde L_s(q(e))\le R_s(e).                 \tag{2.10}
```

Clipping to `[0,1]` preserves the two shape constraints and cannot increase
the error at an attained value. `square`

For `delta=0`, Proposition 2.1 says that a coordinatewise chain of profiles
is exactly parameterized by its positive weighted sum.  The proposition is a
geometric baseline, not yet a mechanism: checking (2.6) over every pair can
be as hard as inspecting the profile table.  The next theorem derives (2.6)
from local ultrametric conditions and scalar linkage.

## 3. Proved deterministic synchronization theorem

Everything in this section is **proved here**.

### Lemma 3.1 (quantitative no-crossing)

Assume that every `R_s` and every sum `R_s+R_t` for distinct species is
`eta`-ultrametric.  If `e={x,y}` and `e'={x,z}` are adjacent, then for every
species `s`,

```math
R_s(e')-R_s(e)
\le {\bigl(q(e')-q(e)\bigr)_++3\eta\over\lambda_s}, \tag{3.1}
```

and the same inequality holds after interchanging `e,e'`.

#### Proof

First prove the local no-crossing statement.  Suppose for species `s,t`
that

```math
a:=R_s(x,z)-R_s(x,y)>3\eta,
\qquad
b:=R_t(x,y)-R_t(x,z)>3\eta.                         \tag{3.2}
```

Write `c_s=R_s(y,z)`.  Approximate ultrametricity, first with `x` as the
distinguished point and then with `y`, gives

```math
R_s(x,y)-\eta\le c_s\le R_s(x,y)+\eta.              \tag{3.3}
```

Indeed, the lower bound is immediate.  If the upper bound failed, both
`R_s(x,z)` and `c_s` would exceed `R_s(x,y)+eta`, contradicting (2.2) with
the roles permuted.  Similarly,

```math
R_t(x,z)-\eta\le c_t:=R_t(y,z)\le R_t(x,z)+\eta.    \tag{3.4}
```

Consequently

```math
c_s+c_t
\le R_s(x,y)+R_t(x,z)+2\eta.                        \tag{3.5}
```

The smaller of the other two values of `R_s+R_t` exceeds the right side's
base value by `min(a,b)>3 eta`.  Thus (3.5) contradicts
`eta`-ultrametricity of `R_s+R_t`.  We have proved

```math
R_s(e')-R_s(e)>3\eta
\quad\Longrightarrow\quad
R_t(e')-R_t(e)\ge-3\eta\quad\text{for every }t.     \tag{3.6}
```

Let `Delta_s=R_s(e')-R_s(e)` and `Delta=q(e')-q(e)`.  If
`0<Delta_s<=3 eta`, (3.1) is immediate.  If `Delta_s>3 eta`, (3.6) gives

```math
\Delta
=\sum_t\lambda_t\Delta_t
\ge\lambda_s\Delta_s-3\eta(1-\lambda_s).
```

Hence

```math
\lambda_s\Delta_s
\le\Delta+3\eta(1-\lambda_s)
\le\Delta_++3\eta.
```

The case `Delta_s<=0` is trivial.  This proves (3.1); swapping the pairs
proves the second assertion. `square`

The proof constant `3 eta` is the sum of the two single-species ultrametric
errors and the pair-sum error in the three-edge contradiction.  No claim of
optimality for this constant is made.

Lemma 3.1 also gives the local cancellation-defect estimate

```math
\mathfrak c(e,e')\le6\eta\qquad(e\sim e').          \tag{3.7}
```

Indeed, write the weighted positive and negative coordinate movements as
`P,N`.  If some positive coordinate movement exceeds `3 eta`, (3.6) gives
`N<=3 eta`; if none does, `P<=3 eta`.  Hence (2.5) is at most `6 eta`.

### Theorem 3.2 (finite deterministic synchronization)

Assume:

1. every `R_s` and every `R_s+R_t` with `s!=t` is `eta`-ultrametric; and
2. the total overlap `q` is `(D,tau)`-monotone-linked.

Put

```math
a=\tau+3D\eta.                                      \tag{3.8}
```

Then, for every species `s`, there is a nondecreasing
`1/lambda_s`-Lipschitz function

```math
L_s:[\min_Eq,\max_Eq]\longrightarrow[0,1]
```

such that

```math
\max_{e\in E}|R_s(e)-L_s(q(e))|
\le {a\over\lambda_s}.                              \tag{3.9}
```

In particular, if all hypotheses are exact (`eta=tau=0`), then every species
overlap is exactly a nondecreasing `1/lambda_s`-Lipschitz function of total
overlap.

#### Proof

Orient two arbitrary pair labels so that `q(e)<=q(e')`, and take a path
supplied by monotone linkage.  For a vector `v=(v_s)_s`, write
`||v||_lambda=sum_s lambda_s|v_s|`.  If `v_j` is the species-profile
increment along the `j`th path edge and `Delta_j` its total-overlap
increment, (3.7) says

```math
\|v_j\|_\lambda\le|\Delta_j|+6\eta.                \tag{3.10}
```

By the triangle inequality and (2.3), the endpoint cancellation defect is

```math
\begin{aligned}
\mathfrak c(e,e')
&\le\sum_j\|v_j\|_\lambda-
       \left|\sum_j\Delta_j\right|\\
&\le\sum_j|\Delta_j|-
       \left|\sum_j\Delta_j\right|+6D\eta\\
&=2B_-+6D\eta\\
&\le2\tau+6D\eta=2a.                              \tag{3.11}
\end{aligned}
```

The same bound holds for every pair of labels.  Proposition 2.1 with
`delta=2a` gives nondecreasing `1/lambda_s`-Lipschitz functions with error
`a/lambda_s`, exactly (3.9). `square`

### Corollary 3.3 (the minimal exact fiber condition)

Assume every `R_s` and every sum `R_s+R_t` for distinct species is exactly
ultrametric.  It is enough for mere functional synchronization, without the
monotonicity conclusion, that every exact fiber

```math
\{e\in E:q(e)=u\}
```

be connected in the pair line graph.  Then there are functions `L_s` on the
attained total-overlap values such that `R_s(e)=L_s(q(e))` exactly.

Indeed, for adjacent pairs at equal `q`, exact no-crossing says every species
difference has the same sign.  Their positive weighted sum is zero, so every
difference is zero.  Connectivity propagates equality across the fiber.

This is weaker than assuming synchronization: fiber connectivity is a
property of the retained scalar `q`, whereas equality of species coordinates
is the conclusion.

### Corollary 3.4 (uniform zero-temperature response bound)

Let `G:E->R` be arbitrary.  Suppose a coupling potential
`Psi:[0,1]^S->R` satisfies the coordinatewise Lipschitz bound

```math
|\Psi(r)-\Psi(r')|
\le\sum_s\kappa_s|r_s-r_s'|.                       \tag{3.12}
```

Define the true and synchronized zero-temperature responses by

```math
V=\max_{e\in E}\{G(e)+\Psi((R_s(e))_s)\},
```

```math
\widetilde V=
\max_{e\in E}\{G(e)+\Psi((L_s(q(e)))_s)\}.
```

Under Theorem 3.2,

```math
|V-\widetilde V|
\le a\sum_s{\kappa_s\over\lambda_s}.               \tag{3.13}
```

This follows by applying (3.9) pointwise and then using
`|max F-max Ftilde|<=||F-Ftilde||_infinity`.  There is no Gibbs measure,
positive temperature, or averaging in the estimate.  For the linear query
`Psi(r)=sum_s theta_s r_s`, take `kappa_s=|theta_s|`.

Thus the deterministic synchronization estimate is stable at zero
temperature.  In an asymptotic model, it is useful precisely when
`tau+D eta` is negligible on the normalization scale of the allowed
couplings.

## 4. A nontrivial closing class

The linkage condition is not vacuous.  Consider a nested chain, with
`C_(k+1)=emptyset` supplying the terminal convention,

```math
\Omega=C_0\supsetneq C_1\supsetneq\cdots
\supsetneq C_k,
```

and scalar levels `t_0<...<t_k`.  Give a pair total overlap `t_j` when `j` is
the first level at which at least one endpoint leaves the continuing cluster
`C_(j+1)` (with the obvious terminal convention).  Equivalently, the
equal-level pair fiber consists of edges of the complete graph on `C_j` that
are not wholly inside `C_(j+1)`.

Each nonempty fiber is connected in the pair line graph: all its edges have
at least one endpoint in the shell `C_j minus C_(j+1)`, and intermediary shell
edges join any two of them.  If `t_i<=t_j`, a pair at level `i` can be joined
to a pair at level `j` through an edge from its shell to an endpoint of the
later pair.  The path has length at most two and never decreases in `q`.
Hence this total-overlap geometry is `(2,0)`-monotone-linked.

Theorem 3.2 therefore forces any species kernels whose individual and
pair-sum geometries are ultrametric to synchronize exactly on this scalar
hierarchy.  Conditional on the hierarchy, the synchronized species data uses
`O(k|S|)` values rather than `Theta(|Omega|^2)` pair values.  Directly
encoding the nested chain itself costs `O(|Omega| log(k+1))` bits; in models
where the hierarchy is fixed or supplied by the query apparatus, that cost is
not charged to the landscape summary.  This is genuine feature compression
for a deterministic hierarchical class, not the fixed-rank Curie--Weiss
algebra already in the repository.

This is a one-spine hierarchy.  In a fully branching GREM, a single
total-overlap fiber generally has several disconnected branch components.
Branch homogeneity calibrates those components in the usual model, but that
calibration is extra structure and should not be hidden under the word
“ultrametric.”

## 5. Sharp counterexample to average and ultrametric synchronization

The following result is **proved here**.

### Proposition 5.1 (rare matching fiber)

Fix `0<rho<=1` and `m>=2`.  For `N=2m`, let `Omega_N={1,...,N}` and let
`e_1,...,e_m` be a perfect matching.  Define two kernels with diagonal one:

```math
R_1(i,j)=\begin{cases}
\rho,&\{i,j\}=e_1,\\
0,&\text{otherwise},
\end{cases}
```

```math
R_2(i,j)=\begin{cases}
\rho,&\{i,j\}\in\{e_2,\ldots,e_m\},\\
0,&\text{otherwise}.
\end{cases}                                         \tag{5.1}
```

Then:

1. `R_1` and `R_2` are positive semidefinite.
2. Every nonnegative combination `alpha R_1+beta R_2` is exactly
   ultrametric.
3. If `X_l` are i.i.d. uniform states, the sampled pair of overlap arrays is
   a weakly exchangeable positive-semidefinite array.
4. For `q=(R_1+R_2)/2`,

   ```math
   \sup_u\operatorname{Var}(R_1(X_1,X_2)\mid q(X_1,X_2)=u)
   ={\rho^2\over m}\left(1-{1\over m}\right)
   \longrightarrow0.                                \tag{5.2}
   ```

5. Nevertheless, for every function `L`,

   ```math
   \max_{e\in E}|R_1(e)-L(q(e))|\ge {\rho\over2}.   \tag{5.3}
   ```

6. On the exact fiber `q=rho/2`, the zero-temperature responses

   ```math
   V(\theta)=
   \max_{e:q(e)=\rho/2}\theta R_1(e)
   ```

   satisfy `V(1)=rho`, `V(-1)=0`.  Replacing `R_1` by any single value
   `L(rho/2)=c` incurs worst error at least `rho/2` over
   `theta in {-1,1}`.

#### Proof

Each kernel is block diagonal with `2 by 2` blocks

```math
\begin{pmatrix}1&\rho\\\rho&1\end{pmatrix}
```

on its selected matching edges and identity blocks elsewhere, so it is
positive semidefinite.  A nonnegative combination has possible nonzero
off-diagonal entries only on a matching.  Any three distinct states contain
at most one matching edge, which directly verifies the ultrametric
inequality.  Sampling a deterministic positive-semidefinite kernel at i.i.d.
states gives a weakly exchangeable Gram array.

The total overlap is `rho/2` on every matching edge and zero on every other
off-diagonal pair.  The diagonal fiber has `q=R_1=1` and zero conditional
variance.  Conditioned on the matching fiber, `R_1` equals `rho` with
probability `1/m` and zero otherwise, proving (5.2) for i.i.d. replicas
without needing to condition them to be distinct.  But `e_1` and
`e_2` have the same total overlap and species-one values `rho` and zero.
This proves (5.3).  Finally, the two signed responses are as displayed, while
the synchronized surrogate gives `c` and `-c`; hence

```math
\max\{|\rho-c|,|c|\}\ge\rho/2.
```

`square`

The high-overlap fiber consists of the disjoint matching edges, so its line
graph is totally disconnected.  The example satisfies the local
no-crossing/ultrametric part of Theorem 3.2 perfectly and fails exactly the
cross-root linkage part.

It also rules out several tempting weaker hypotheses:

- positive semidefiniteness and nonnegative overlaps;
- weak exchangeability, even from an i.i.d. directing measure;
- ultrametricity of the total, each species, every pairwise sum, or indeed
  every nonnegative mixture;
- vanishing mean-square prediction error from total overlap;
- vanishing conditional variance uniformly over total-overlap levels.

The last failure is the extremal-information point.  Conditional variance
weights a single exceptional matching edge by `1/m`; a zero-temperature
maximum can select it with full weight.

## 6. What the response roof must retain when synchronization fails

On a total-overlap fiber `q=u`, deterministic synchronization tries to replace
the species vector by one point `L(u)`.  Proposition 5.1 shows that the correct
query-relative quotient may instead need the upper response roof of the
fiber,

```math
\theta\longmapsto
\max_{e:q(e)=u}\sum_s\theta_sR_s(e).                \tag{6.1}
```

For the matching example, this roof retains the interval between the two
species profiles and answers both signs exactly.  It is still much smaller
than a labeled pair table in this example.  Under composition, however, a
later rooted query may distinguish different pairs on the same exposed face;
then another feature must be added.  Thus (6.1) is a strict repair for the
declared species-linear experiment, not universal sufficiency.

This suggests a precise dichotomy for the broader program:

- **synchronized fibers:** a point-valued order parameter is sufficient;
- **unsynchronized fibers:** retain a fiber response body, and measure its
  exposed-face entropy rather than its average variance.

This is more than terminology because Theorem 3.2 supplies a checkable
criterion and error bound for the first case, while Proposition 5.1 proves
that average criteria cannot decide between the two.

## 7. Candidate hypotheses rejected or retained

| Candidate finite hypothesis | Verdict | Evidence |
|---|---|---|
| Uniform cancellation defect at most `delta` | Retained as a geometric baseline; gives error `delta/(2 lambda_s)` | Proposition 2.1 |
| Small `E Var(R_s given q)` | Rejected for zero temperature | Proposition 5.1 |
| Small `sup_q Var(R_s given q)` | Rejected for zero temperature | Proposition 5.1 |
| PSD + weak exchangeability | Rejected | Proposition 5.1 |
| Every nonnegative mixture ultrametric | Gives local no-crossing but not global synchronization | Lemma 3.1 and Proposition 5.1 |
| Exact total-overlap fiber connectivity + mixture ultrametricity | Retained for exact functional synchronization | Corollary 3.3 |
| Approximate monotone linkage + approximate mixture ultrametricity | Retained with uniform error | Theorem 3.2 |
| Ghirlanda--Guerra identities | Imported sufficient in the limiting Gibbs experiment | Panchenko 2015 |
| Injective scalar `q` | Formally sufficient but rejected as compression | it can encode the whole pair table |

## 8. Minimal next questions

### 8.1 Can Ghirlanda--Guerra defects yield finite linkage constants?

Theorem 3.2 deliberately does not claim a quantitative finite-volume version
of Panchenko's theorem.  The next strong result would derive `(D_N,tau_N)`
linkage, or a substitute with the same uniform conclusion, from explicit
finite replica-identity defects and an anti-rare-face condition.

**Falsifier:** a family with vanishing extended Ghirlanda--Guerra defects at
every fixed replica order but a fixed conditioned response gap.  Fixed-order
identities may still miss a zero-entropy exceptional fiber.

### 8.2 Replace uniform linkage by exposed-fiber linkage

Full linkage of every pair may be stronger than the query requires.  Given a
declared compact coupling class, is it enough to link only pair profiles that
belong to an exposed face of the relevant response body?  A useful theorem
would bound response error by

```math
\text{ultrametric defect}
+\text{backtracking on exposed fibers},
```

with no dependence on unexposed states.

### 8.3 Branch calibration in homogeneous GREM

Fully branching hierarchical models fail literal fiber connectivity, but
their automorphisms calibrate branch components.  Formulate a finite
condition, weaker than directly assuming `R_s=L_s(q)`, that combines
within-branch linkage with orbit invariance and yields the same conclusion.
This would connect Theorem 3.2 to a standard model rather than the one-spine
class only.

### 8.4 Complexity accounting

Even exact synchronization is not compression if the scalar `q` has one
distinct value per pair.  Any promoted axiom should therefore record both:

1. the synchronization error; and
2. the metric entropy or number of exposed levels of the scalar response
   roof.

## 9. Specialist judgment

The deterministic content that survives is not “ultrametricity implies
synchronization.”  That statement is false in an unusually strong way.
What survives is:

```math
\boxed{
\text{mixture ultrametricity}
+\text{cross-root calibration}
\Longrightarrow
\text{uniform deterministic synchronization}.}
```

Theorem 3.2 gives one finite calibration condition and a zero-temperature
error bound.  Proposition 5.1 shows that it cannot be replaced by average
conditional variance, PSD, exchangeability, or even ultrametricity of every
mixture.

This is a Level-3 result for the nascent framework in a limited sense: it
formulates and proves a theorem that was not present in the initial convex
response-roof formalism, and it predicts exactly when a vector-valued fiber
roof collapses to a scalar order parameter.  It is not yet a general
deterministic Parisi theory.  The most promising next theorem is an
**exposed-fiber synchronization theorem** deriving controlled calibration
from finite replica identities plus a quantitative exclusion of rare exposed
exceptions.
