# Uniform shifted-absolute comparison for normalized Rademacher sums

2026-09-17. Derived in the Bernoulli track from the classical scalar
Stein equation and the exact symmetric trapezoid identity. The
localization researcher independently reconstructed every finite step
and returned PASS. The director uses this lemma in the ordinary-sign-
bridge regularizer. No external novelty claim is made.

## 1. Uniform scalar theorem

Let `S_q=q^(-1/2) sum_(i=1)^q epsilon_i`, with independent fair signs,
and let G be standard normal. For EVERY integer q>=1,

```
sup_(t in R) |E|S_q-t|-E|G-t|| <=3/q.                 (1)
```

The shift t can depend on arbitrary conditioned offsets. The estimate
also covers shifts exactly at lattice atoms or halfway between atoms.
Its constant is deliberately conservative.

### 1.1 Classical Stein input and distributional derivative

For `h_t(x)=|x-t|`, the bounded solution f=f_t of

```
f'(x)-x f(x)=h_t(x)-E h_t(G)
```

satisfies, uniformly in t,

```
||f'||_infinity<=kappa=sqrt(2/pi),
||f''||_infinity<=2.                                 (2)
```

These classical Lipschitz-test Stein factors are stated in Ross,
[Fundamentals of Stein's method, Lemma 2.5 and Theorem 3.1](https://arxiv.org/pdf/1109.1880).
The factors and the leave-one-out mechanism in Sections 2.1--3.1 were
read directly. A Lipschitz test is interpreted by smooth approximation
or almost-everywhere derivatives; no third-derivative factor for an
arbitrary Lipschitz test is assumed.

For this PARTICULAR test, differentiating the Stein equation twice in
the distributional sense gives the precise signed measure identity

```
f'''(x)=[2f'(x)+x f''(x)] dx +2 delta_t.               (3)
```

Here f' is continuous, f'' is smooth away from t with a jump of size
two at t, and the regular density is bounded in absolute value by
`2kappa+2|x|`. Thus there is no shift-dependent bound growing with |t|.
The original first-order equation, its absolutely continuous solution,
and the jump in the derivative of |x-t| justify (3) without an exchange
of an unjustified third derivative and expectation.

### 1.2 Exact symmetric identity and its continuous contribution

Put `a=q^(-1/2)` and `U=S_q-a epsilon_1`. Symmetry and independence
give the exact Stein discrepancy

```
E h_t(S_q)-E h_t(G)
 =E[f'(S_q)-S_q f(S_q)]
 =E_U{ [f'(U+a)+f'(U-a)]/2
                    -(1/(2a)) integral_(-a)^a f'(U+s) ds }.
                                                               (4)
```

Indeed `q a^2=1`, and the conditional sign average of
`a epsilon_1 f(U+a epsilon_1)` is
`a[f(U+a)-f(U-a)]/2`.

For a continuous function whose second derivative is a locally finite
measure, integration by parts twice gives the trapezoid formula

```
[u(a)+u(-a)]/2-(1/(2a)) integral_(-a)^a u(s) ds
       =(1/(4a)) integral_[-a,a] (a^2-s^2) d u'(s).   (5)
```

It applies to `u(s)=f'(U+s)` by (3). Since
`integral_(-a)^a(a^2-s^2)ds=4a^3/3` and E|U|<=1,
the absolute expected contribution of the regular part of (3) is at
most

```
(a^2/3)[2kappa+2E|U|+2a]
                       <=(2kappa+4)/(3q).            (6)
```

### 1.3 The kink costs another O(1/q), uniformly in its location

The point mass in (3) contributes

```
(1/(2a)) E[(a^2-(t-U)^2)_+]
             <=(a/2) Pr{|t-U|<a}.                    (7)
```

The support of U has spacing exactly 2a. Consequently the open interval
`(t-a,t+a)` contains at most ONE atom. If an atom occurs at an endpoint,
its kernel weight is zero; no half-mass convention changes the bound.
Hence the probability in (7) is at most the largest binomial atom.

An elementary bound sufficient here is

```
max_atom(sum_(i=1)^(q-1) epsilon_i)<=sqrt(2/q).         (8)
```

For 2r signs the maximum is
`product_(j=1)^r(1-1/(2j))<=exp[-(1/2)sum_j 1/j]
 <=1/sqrt(r+1)`; for 2r+1 signs it is no larger than for 2r.
In both cases r+1>=q/2. The q=1 empty sum also obeys (8).
Thus (7) is at most `1/(sqrt(2)q)`.

Combining (6)--(8), uniformly over ALL real shifts, gives

```
sup_t |E|S_q-t|-E|G-t||
 <=[(2sqrt(2/pi)+4)/3+1/sqrt(2)]/q <3/q.
```

This proves (1). A sharper central-binomial estimate improves the
displayed constant, but that optimization is immaterial here.

## 2. All-offset Boolean-field comparison

Let H be ANY real function on {+-1}^n. Let Z have independent
coordinates distributed as S_q, and let G have independent standard
Gaussian coordinates. Then for every delta>=0,

```
|E max_x[H(x)+delta Z dot x]
       -E max_x[H(x)+delta G dot x]| <=3delta n/q.    (9)
```

Proof: condition on every field coordinate except coordinate i. Split
the Boolean queries by x_i=+1 or -1, and optimize all other coordinates
inside each class. The remaining maximum is exactly

```
max(A+delta z,B-delta z)
 =(A+B)/2+delta |z-(B-A)/(2delta)|.
```

Its shift is fixed under the conditioning. Apply (1), then replace the
n independent field coordinates sequentially. If delta=0 the statement
is immediate. Conditioning on an additional independent Gaussian field
merely changes H and leaves (9) valid. Thus both conditional-increment
endpoints of the regularizer retain their arbitrary child offsets.

This is stronger than generic 1-Lipschitz Wasserstein replacement,
which only gives a q^(-1/2) rate. It uses the precise two-slope geometry
of a Boolean coordinate, plus the lattice small-ball estimate. It does
NOT compare random near-level sets by weak convergence.

## 3. Consequence for the actual ordinary-sign bridge

For an n by q iid sign bridge, every FIXED new-spin pattern induces
`delta=sqrt(q/n)` and independent old fields with the law S_q after
normalizing energy by sqrt(n). The two scalar extrema replacements in
the secondary-field proof therefore cost at most
`6sqrt(n/q)` together. Choose secondary variance `s=1/n`.
The Gaussian radial increment is at most
`sqrt(2/pi) sqrt(n/q)/2`. Thus

```
E Delta_s <=[6+sqrt(2/pi)/2]sqrt(n/q).                 (10)
```

Each bridge sign changes Delta_s by at most 4/sqrt(n), so its centered
MGF proxy is 4q. A deviation A sqrt(n/q), followed by a union over
all 2^q patterns, has failure probability at most

```
exp[q log2-A^2 n/(8q^2)].                             (11)
```

For `q^3<=n`, a sufficiently large universal A makes this at most
`exp[-c n/q^2]`. A physical near-level tolerance `T=n/sqrt(q)`
becomes margin O(sqrt(n/q)); a new child with cap at most q^(3/2)
only adds a lower-order margin in this q range. Dividing by sqrt(s)
therefore bounds each old branch width by O(n/sqrt(q)). The pattern
union width O(sqrt(nq)) and the q new coordinates are lower order.
The audited I-MMSE entropy conversion yields

```
w(E_parent(n/sqrt(q)))=O(n/sqrt(q)),
log|E_parent(n/sqrt(q))|=O((n/sqrt(q))sqrt(log q)).     (12)
```

The ordinary bridge still costs O(n sqrt(q)) in cap. Same-order
restriction/refill changes only q designated vertices and preserves
these estimates, allowing an immaterial constant factor in q^3/n.
For `q~n^(1/3)`, the cap cost is O(n^(7/6)), while the physical
window and width are O(n^(5/6)). These implications retain all new
signs, all absolute polarities, and the original old quadratic block.
The director owns the canonical full-sign theorem and its same-order
variational consequence. No convergence or improved cap constant is
claimed from this regularization alone.

## 4. Exhaustive-threshold finite diagnostics

`computations/paper_bernoulli_2026_09_17_shifted_absolute.py` checks
517 orders through q=10000, including every q<=512. On each lattice
interval the difference of absolute expectations is concave. Therefore
its extreme values occur at lattice endpoints or at the Gaussian
quantile corresponding to that interval's binomial CDF. Checking these
finite candidates analytically covers every real threshold; Gaussian
functions and binomial probabilities are evaluated numerically.

All 99,080 tested candidates passed. The largest measured q times
the uniform error was below 0.205, far below the proved conservative
constant three. These diagnostics do not replace the uniform proof.

The discrepancy researcher independently checked every q<=1000,
testing 319,809 endpoint/stationary extrema in a separate implementation.
Its worst q-scaled error was 0.20442252894472 at q=3,t=0, agreeing with
this replay. The independent code is
`computations/paper_discrepancy_2026_09_17_shifted_abs_clt.py`.

## 5. Sharp scales and a genuine full-sign scope separation

The RATE in (1), unlike its displayed constant, cannot be improved
uniformly. For even q, the elementary exact identity is

```
E|S_q|=sqrt(q) 2^(-q) binom(q,q/2).
```

Stirling's formula gives
`E|S_q|=sqrt(2/pi)[1-1/(4q)+O(q^(-2))]`. Thus the absolute-expectation
error is of order 1/q already at shift zero.

There is also a sharpness example for the analytic nearcode WIDTH rate.
Take the arbitrary old offset H identically zero and even q. Each
normalized q-sign field coordinate is exactly zero with probability
`p0=2^(-q) binom(q,q/2)~sqrt(2/pi)/sqrt(q)`.
Every nonzero-field spin is pinned in an exact maximizer, whereas all
zero-field spins are free. If Z is their number, the exact maximizing
code has width `sqrt(2/pi) Z` and entropy `Z log2`. For n/sqrt(q)
tending to infinity, Z concentrates about n p0. Hence the typical
exact-ground width is `[2/pi+o(1)]n/sqrt(q)`, even at tolerance zero.
This example is an arbitrary-offset field landscape, NOT a full-sign
quadratic child. It calibrates the analytic theorem without asserting
a quadratic obstruction.

A separate ACTUAL full-sign example shows why the regularizer does not
yield favorable fixed-normalized-window isotropic response geometry.
Start with the known-half alternating-Hadamard family from
[the all-law quadratic dual](paper_localization_all_law_hadamard_obstruction_2026_09_17.md).
Apply the same-order compiler with q=O(N^(1/3)). It rewrites at most
O(Nq)=O(N^(4/3)) edges. Pointwise energy and cap change by at most
twice that number. The original explicit eigenword law had absolute
deficit at most N; consequently every word in its support remains in
EVERY fixed positive eta N^(3/2) near-level code of the rewritten
matrix, eventually as N tends to infinity.

The physical eigenword law and its quadratic minorant are unchanged.
Therefore its universal isotropic-response lower bound
`[1/sqrt(18)+o(1)]sqrt(N)` and vanishing-response center-radius lower
bound 1/72 still apply to those fixed-eta nearcodes. At the SAME time,
the regularizer guarantees small width and entropy in its much smaller
N^(5/6) window. These are actual full-sign families, not arbitrary
abstract codes. The distinction between mesoscopic and fixed-normalized
windows is essential; no low-response-center consequence follows by
silently changing that order of limits.
