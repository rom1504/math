# Independent audit: exposed Boolean synchronization

**Audited files.** `exposed_boolean_synchronization.md` and
`verify_exposed_boolean_synchronization.py`.

## Verdict

**PASS WITH WORDING REPAIRS.**  BS.1 has the correct constants and
quantifiers; BS.2's gamma constant and covering-radius consequence are
correct; BS.3 is exact; and the self-dual Walsh construction and BS.4
normalizations check.  Exact-sign completion preserves the stated recovery
scale.

Required qualifications:

1. The theorem is an exposed, query-dependent sufficient condition.  The
   near-flat vector must be near the **global** spherical maximum, not merely
   optimal in an arbitrary channel.
2. “Minimum uncertainty scale” is not proved in this note.  The construction
   proves a Hamming difference of exactly `sqrt(n)`; minimality needs a cited
   uncertainty theorem or should be dropped.
3. BS.1 shows that coordinate flatness is sufficient beyond the spherical
   value.  The examples here do not prove an information lower bound saying
   flatness is unavailable from `(G,R)`: that would require equal Gram data
   with different exposed flatness or Boolean response.  In the two-pole
   family, the flatness is actually a function of the correlation `rho`.
4. “Checkable” means no `2^n` Boolean optimization.  Finding the globally
   exposed channel can still require handling `2^(p+1)` trust channels when
   `p` grows.

The missing `left` delimiter in BS.29 is a formatting issue only and is not
repeated here.

## 1. BS.1 constants and quantifiers

For `x=sgn(u)`, including arbitrary choices at zero coordinates,

```math
\|x-u\|_2^2
=\sum_i(1+u_i^2-2|u_i|)
=2n-2\|u\|_1
=2n\phi(u).                                          \tag{ABS.1}
```

Put `delta=||x-u||/sqrt(n)=sqrt(2phi(u))`.  Symmetry of `H`,
`||H||op=r`, and `||x||=||u||=sqrt(n)` give

```math
\left|{x^THx-u^THu\over2}\right|
\le {r\over2}\|x-u\|\,(\|x\|+\|u\|)
=rn\delta.                                           \tag{ABS.2}
```

For the fixed endpoint word,

```math
\|z_\epsilon\|_2
\le\sum_{a=1}^p\|w_a\|_2=p\sqrt n,                 \tag{ABS.3}
```

so the field changes by at most

```math
m\|z_\epsilon\|\|x-u\|
\le mpn\delta=crn\delta.                            \tag{ABS.4}
```

Since the Boolean response includes the **same** channel evaluated at `x`,
BS.6 implies

```math
0\le\mathcal S-\mathcal B
\le rn\{\xi+(1+c)\sqrt{2\varphi}\}.                 \tag{ABS.5}
```

There is no missing factor two.  The nonnegativity on the left follows
because the Boolean cube is contained in the sphere.

The existential quantifier is exact but important: the channel/vector pair
must satisfy `F>=S-xi rn`, where `S` is the maximum over *all* channels.
A flat optimizer in a losing channel proves nothing.  One may take
`xi,varphi>=0`; negative values are either meaningless or make the premise
impossible.  The canonical verifier checks the local Lipschitz inequality
on 200 random vectors.  The independent verifier checks ABS.1 and the full
response bound at the exposed Walsh optimizer.

## 2. BS.2 gamma formula and covering quantifiers

For uniform `theta` on `S^(d-1)`,

```math
\mathbb E|\theta_1|
={\Gamma(d/2)\over\sqrt\pi\,\Gamma((d+1)/2)}.        \tag{ABS.6}
```

If `Q` has orthonormal columns and rows `q_i`, rotational invariance gives
the same constant times `||q_i||`.  Since

```math
\sum_i\|q_i\|^2=d,
\qquad
\sum_i\|q_i\|\le\sqrt{nd},                          \tag{ABS.7}
```

the average `l_1` norm of `sqrt(n)Qtheta` is at most

```math
n\sqrt d\,{\Gamma(d/2)\over
 \sqrt\pi\Gamma((d+1)/2)}=n\gamma_d.                \tag{ABS.8}
```

Thus **some** sphere point has `l_1/n<=gamma_d`.  For `d>=2`, strict
Cauchy--Schwarz gives `E|theta_1|<1/sqrt(d)`, hence `gamma_d<1`.  Applying
ABS.1 to that one point yields

```math
\sup_{u\in U,\|u\|^2=n}\min_x{\|u-x\|\over\sqrt n}
\ge\sqrt{2(1-\gamma_d)}.                             \tag{ABS.9}
```

The order of `sup` and `min` in BS.12 is therefore correct.  The proposition
rules out a uniformly fine Boolean net for the whole sphere; it does not say
that every vector in a `d>=2` subspace is nonflat.  At `d=2`, the verifier
also confirms `gamma_2=2sqrt(2)/pi`.

## 3. Exact correlated-port response

For Boolean `a,b` with `a^Tb=rho n`, agreements and disagreements number
`n(1+rho)/2` and `n(1-rho)/2`.  Hence

```math
\|a+b\|_1=n(1+\rho),
\qquad
\|a-b\|_1=n(1-\rho).                                 \tag{ABS.10}
```

When `rho in [0,1]`, the plus word is larger.  Boolean duality bounds the
field by `n(1+rho)`, and `x=a` attains this simultaneously with child energy
`rn/2`.  This proves BS.14.

On the sphere, the largest endpoint field norm is
`sqrt(2n(1+rho))`; multiplying by the sphere radius gives
`n sqrt(2(1+rho))`.  The normalized sum in BS.16 is a `+r` eigenvector and
attains both bounds, proving BS.15.  Thus BS.18 has the correct factor
`m/r`.  For negative `rho`, the analogous formula uses `|rho|`; the stated
restriction avoids that extra case.

In GE coordinates the winning positive channel has

```math
g=h=2(1+\rho),
\qquad a_{GE}=4(1+\rho),\quad b_{GE}=0.              \tag{ABS.11}
```

GE.14 consequently gives
`S/(rn)=1/2+(m/r)sqrt(2(1+rho))`, exactly BS.15.  The
independent verifier checks this normalization and independently exhausts
the order-16 Boolean response at `rho=1/2`.

As elsewhere, a literal `m`-wide auxiliary shore has integer `m`; arbitrary
real `m>=0` is an analytic field-weight extension.

## 4. Self-dual Walsh construction

Let the Walsh pairing on `V=F_2^d times F_2^d` be the standard dot product.
For `y_0(x,z)=(-1)^(x dot z)`, summing first over one coordinate gives

```math
Wy_0=qy_0.                                             \tag{ABS.12}
```

When `d` is even, the coordinate-pair swap `M` is symmetric,
fixed-point-free, and satisfies `M^2=I`.  Its graph
`L={(x,Mx)}` is totally isotropic because

```math
(x,Mx)\mathbin\cdot(y,My)
=x\mathbin\cdot y+(Mx)\mathbin\cdot(My)=0           \tag{ABS.13}
```

in characteristic two.  It has dimension `d` in the `2d`-dimensional
space, hence `L=L^perp`.  Also `x dot Mx=0`, so `y_0=1` on `L`.  The standard
subspace transform is therefore

```math
W\mathbf1_L=q\mathbf1_{L^\perp}=q\mathbf1_L.         \tag{ABS.14}
```

It follows immediately that `y_1=y_0-2 1_L` is Boolean and satisfies
`Wy_1=qy_1`.

Conjugation by `D_(y_0)` preserves symmetry and the Hadamard identity.  Its
trace is

```math
\operatorname{tr}(D_{y_0}WD_{y_0})
=\operatorname{tr}W=0,                               \tag{ABS.15}
```

and ABS.12 regularizes the all-ones word.  The second pole is
`b=1-2 1_L`, so it differs on exactly `q=sqrt(n)` coordinates and has
correlation `1-2/q`.  Every claimed construction identity is correct.

The phrase “minimum uncertainty scale” is not established by these
calculations.  They establish the exact `sqrt(n)` support.  A claim of
minimality among nonzero Walsh eigendirections needs an uncertainty/support
theorem and its hypotheses; otherwise “at the uncertainty scale” should be
replaced by “at scale `sqrt(n)`.”

## 5. BS.4 flatness, gap, and completion

The sum `a+b` vanishes on `L` and equals two elsewhere.  Since
`rho=1-2/q`, normalization gives

```math
{\|u_*\|_1\over n}
={1-1/q\over\sqrt{1-1/q}}
=\sqrt{1-1/q}.                                      \tag{ABS.16}
```

Thus BS.28 is exact.  With `m=q/2`, one has `c=1`, and BS.18 becomes

```math
{\mathcal S-\mathcal B\over qn}
=\sqrt{1-1/q}-(1-1/q)=O(1/q).                       \tag{ABS.17}
```

Since `q=sqrt(n)`, this is `O(n^(-1/2))`.  The full uncompleted order is
`N=n+q`, so normalization by `N^(3/2)` changes ABS.17 by a factor tending
to one.

There are `2m=q` auxiliary vertices.  For any fixed public hollow signing
on them,

```math
Q(C)\le\binom q2=O(q^2)=O(n).                       \tag{ABS.18}
```

Both responses change by at most `Q(C)`, so

```math
0\le S_C-B_C
\le(S-B)+2Q(C)=O(n).                                \tag{ABS.19}
```

After division by `qn=n^(3/2)`, completion preserves the `O(1/q)=o(1)`
recovery conclusion.  It may change the leading coefficient in BS.29 or
even close the uncompleted gap; the draft claims only survival of recovery,
which is correct.

The independent verifier exhausts two complete signings at `q=4`, computing
both the Boolean and old-spin spherical optima, and confirms the Lipschitz
scope directly.

## 6. Proposed common-pole corollary

The proposed multi-port corollary is **valid with no missing channel
factor**.  Suppose `x_0` is a Boolean `+r` or `-r` eigenvector and define

```math
\delta
=1-{1\over pn}\sum_{i=1}^p|w_i^Tx_0|.              \tag{ABS.20}
```

For every spherical old spin, the child term is at most `rn/2` and each of
the `p` fields is at most `n`, so

```math
\mathcal S\le {rn\over2}+mpn.                        \tag{ABS.21}
```

Evaluate the Boolean response at `x_0`.  Choose `sigma` to make its
quadratic eigenvalue positive and choose every endpoint sign independently
to match `w_i^Tx_0`.  Then

```math
\mathcal B
\ge {rn\over2}+m\sum_i|w_i^Tx_0|
={rn\over2}+mpn(1-\delta).                           \tag{ABS.22}
```

Subtracting proves

```math
0\le\mathcal S-\mathcal B
\le mpn\delta=c\delta rn,
\qquad c={mp\over r}.                                \tag{ABS.23}
```

The outer child sign and endpoint signs are separate channel variables, so
there is no sign conflict.  The normalization is exact.

For the synchronized Walsh pair, choosing `x_0=a` gives
`delta=1/q` and `c=1`; ABS.23 bounds the normalized gap by `1/q`, consistent
with the exact asymptotic `1/(2q)+O(q^(-2))`.  The independent verifier
checks this instance.

If a complete signing is added on the `pm` auxiliary vertices, both values
move by at most its cap.  At bounded `c`, `pm=cr`; for a Hadamard signing
`r^2=n`, so

```math
Q(C)\le\binom{pm}{2}=O(n),
\qquad
S_C-B_C\le c\delta rn+2Q(C).                         \tag{ABS.24}
```

Thus the proposed `O(n)` completion stability is correct.  For an
`o(rn)` recovery conclusion one needs `delta=o(1)` (with bounded `c`), as
expected.

This corollary is simpler and stronger than invoking exposed spherical
flatness when a common Boolean pole is available: it never needs to find or
round the spherical optimizer.

## 7. Information-content scope

BS.1 is a genuine new sufficient mechanism: one coordinate-level statistic
of a globally exposed optimizer converts the sphere certificate to a
Boolean witness.  It does not require a Boolean response histogram.

However, the statement that the condition is “stronger than Gram data
alone” is not demonstrated by the orthogonal example.  A property can be
false at `G=R=I_2` and still be determined by Gram data.  Indeed, for any
two Boolean top poles with correlation `rho>=0`, the exposed sum in BS.16
has

```math
{\|u_*\|_1\over n}=\sqrt{{1+\rho\over2}},            \tag{ABS.25}
```

which is already a function of the two-port Gram matrix.  In more general
multi-port/eigenspace systems, the coordinate `l_1` norm need not be fixed
by `(G,R)`, but proving that requires a collision or separate theorem.

Similarly, the rounded witness is selected after the future channel is
known.  BS.1 is a query-local recovery theorem, not by itself a
depth-independent reusable Boolean carrier.  These distinctions do not
weaken the quantitative result.

## 8. Verifier assessment

The canonical verifier passes 242 checks.  It validates the Walsh family at
orders 16 and 256, exhausts the order-16 Boolean response, checks the local
rounding constant on 200 random instances, and verifies `gamma_d<1` through
`d=8`.

The independent verifier adds:

- exact BS.3/GE normalization and a separate order-16 enumeration;
- BS.1 on the actual exposed Walsh optimizer;
- two completed-signing Boolean/spherical computations;
- the closed `d=2` gamma value;
- the common-pole corollary at `delta=1/q`.

Run:

```bash
./.venv/bin/python \
  extremal_information/experiments/verify_exposed_boolean_synchronization.py

./.venv/bin/python \
  extremal_information/experiments/verify_exposed_boolean_synchronization_independent_audit.py
```

Both pass.
