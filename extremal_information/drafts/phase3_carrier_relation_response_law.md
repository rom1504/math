# Relation rank is charged only through carrier Grassmannian capacity

**Status.** The statements below are proved.  Finite instances for redundant,
two-scale, Lee, and flag-ultrametric carriers are checked by
[`verify_phase3_carrier_relation_response_law.py`](../experiments/verify_phase3_carrier_relation_response_law.py).

The exact gluing theorem counts `D kappa` field coordinates of mixed
holonomy.  That count is a latent gauge count, not a response-information
law.  This note gives both a sharp counterexample and a replacement law:

> At scales larger than the circuit toll, multichannel response entropy is
> the Hausdorff metric entropy of linear subspaces of the carrier.

Thus `carrier dimension times relation rank` is valid when the carrier
contains a linearly large, macroscopically separated code.  Cardinality, or
even cardinality together with linear diameter, is not sufficient.

## 1. Metric-carrier mixed holonomy

The geometric core does not depend on parallel-pair circuits.  We record it
first.  Let `Z` be a finite set with a map `h:Z\to W` and a toll
`c:Z\to[0,tau]`, and define

```math
T_{h,c}(u)=\min_{z\in Z}\bigl(c(z)+d(u,h(z))\bigr).          \tag{CR.0}
```

### Lemma CR.0 (bounded-toll distance transform)

If `C=h(Z)`, then

```math
d(u,C)\le T_{h,c}(u)\le d(u,C)+\tau.                         \tag{CR.0a}
```

For two such transforms, with image sets `C,C'` and tolls in the same
interval `[0,tau]`,

```math
\left|\,\|T_{h,c}-T_{h',c'}\|_\infty-d_H(C,C')\,\right|
\le\tau.                                                     \tag{CR.0b}
```

#### Proof

Dropping the toll gives the lower bound, while a preimage of a nearest image
point gives the upper bound.  Thus each transform is its distance-to-image
function plus an error taking values in `[0,tau]`.  Distance-to-set functions
have uniform distance exactly equal to Hausdorff distance, and two errors in
the interval `[0,tau]` differ by at most `tau`. `\square`

This elementary lemma becomes an information law when `Z` is a relation
space: its image is the holonomy code and `tau` is the maximum cost of a
chosen relation representative.  It says that whenever `tau=o(s_D)`, the
macroscopic query geometry forgets the presentation of the relations and
retains precisely the metric geometry of their image.

Let `W=F_q^D` carry an arbitrary translation-invariant metric `d`, and write
`||w||=d(w,0)`.  Realize the kernel by the complete weighted alphabet

```math
K_d=\{(w,0):w\ne0\},\qquad \operatorname{cost}(w,0)=\|w\|.
```

Using several kernel letters cannot beat one letter, by the triangle
inequality.  Let `Q=F_q^k` with basis `q_1,...,q_k`; relation letters have
unit cost.  For a linear map `V:F_q^k\to W`, with columns `v_j`, put

```math
P=K_d\cup\{(0,\alpha q_j):\alpha\ne0\},
```

```math
R_V=K_d\cup\{(\alpha v_j,\alpha q_j):\alpha\ne0\}.          \tag{CR.1}
```

Every `R_V` is individually shear-isometric to `P`, by the kernel-fixing
map `(w,x)\mapsto(w-Vx,x)`.  Nevertheless the two-fragment union has `V` as
its mixed holonomy.  Distinct maps give all `q^{Dk}` exact gluing classes.

### Theorem CR.1 (metric normal form and Grassmannian rough isometry)

For every kernel endpoint `u\in W`,

```math
F_V(u):=\ell_{P\cup R_V}(u,0)
=\min_{z\in F_q^k}\bigl(2\operatorname{wt}(z)+\|u+Vz\|\bigr). \tag{CR.2}
```

For every `V`, with `C_V=\operatorname{im}V`,

```math
d(u,C_V)\le F_V(u)\le d(u,C_V)+2k.                         \tag{CR.3}
```

If `V` is injective then `C_V` is a `k`-subspace.  Consequently, for
arbitrary ordered bases `V,V'` of `k`-subspaces `C,C'`,

```math
\boxed{
\left|\,\|F_V-F_{V'}\|_\infty-d_H(C,C')\,\right|\le2k,
}                                                           \tag{CR.4}
```

where `d_H` is Hausdorff distance in `(W,d)`.

#### Proof

In a word for `(u,0)`, sum the coefficients of the `R_V` letters in channel
`j` to `z_j`.  Quotient cancellation forces the `P` coefficients to sum to
`-z_j`.  A nonzero `z_j` costs at least one letter of each kind; a zero sum
has zero holonomy and can be deleted.  The remaining kernel correction is
`u-Vz` and costs at least `||u-Vz||`.  Conversely one `R_V` letter and one
`P` letter realize each nonzero `z_j`, followed by the single weighted
kernel correction.  Replace `z` by `-z` to obtain (CR.2).

The lower bound in (CR.3) drops the nonnegative toll.  A nearest codeword has
a coordinate vector of weight at most `k`, proving the upper bound.  Finally,
in every finite metric space,

```math
\|d(\mathord\cdot,C)-d(\mathord\cdot,C')\|_\infty=d_H(C,C').
```

Apply the reverse triangle inequality using (CR.3). `\square`

Let `Pack(X,r)` denote the largest cardinality of a subset of a metric space
with all pairwise distances greater than `r`.  Let `Resp_k(W)` be the family
of profiles (CR.2) over injective `V`, with uniform distance, and equip the
Grassmannian `Gr_k(W)` with `d_H`.  Equation (CR.4) immediately gives, for
`r>2k`,

```math
Pack(Gr_k(W),r+2k)
\le Pack(Resp_k(W),r)
\le Pack(Gr_k(W),r-2k).                         \tag{CR.5}
```

This is the promised sharp replacement for the raw gauge count.  When
`k=o(s_D)`, response packing at scale `Theta(s_D)` and Grassmannian packing
at that scale have the same exponential order, up to an `o(s_D)` shift of
the separation threshold.

## 2. A usable sufficient hypothesis

### Theorem CR.2 (linearly rich separated-host criterion)

Suppose `C_0\le W` has dimension `r` and minimum nonzero metric weight

```math
d_0=\min_{c\in C_0\setminus\{0\}}\|c\|.
```

For every `k\le r`, there are at least

```math
{r\brack k}_q\ge q^{k(r-k)}                    \tag{CR.6}
```

locally shear-trivial two-fragment profiles with pairwise response distance
at least

```math
d_0-2k.                                        \tag{CR.7}
```

In particular, if `r\ge alpha D`, `d_0\ge delta s_D`, and
`k\le\min(alpha D/2,delta s_D/4)`, then there are

```math
2^{(alpha/2)Dk\log_2q}                         \tag{CR.8}
```

profiles separated by at least `delta s_D/2`.

#### Proof

Choose one ordered basis for each `k`-subspace of `C_0`.  If `C\ne C'`, take
`c\in C\setminus C'`.  Every `c-c'`, for `c'\in C'`, is a nonzero word of
`C_0`, so `d_H(C,C')\ge d(c,C')\ge d_0`.  Apply (CR.4) and the standard
Gaussian-binomial lower bound. `\square`

The hypothesis is geometric, not merely cardinal: the carrier must contain
a linear host with both positive rate and positive relative metric distance.
It is also portable beyond Hamming geometry.

### Example CR.3 (flag ultrametric)

On `F_q^D`, put

```math
\|x\|_{flag}=\max(\{i:x_i\ne0\}\cup\{0\}).                  \tag{CR.9}
```

This is a translation-invariant ultrametric of diameter `D`.  The span of
the last `r` coordinate vectors has minimum nonzero weight `D-r+1`.  Taking
`r=\lfloor D/2\rfloor` and `k\le D/8` gives at least

```math
q^{k(\lfloor D/2\rfloor-k)}
```

profiles separated by more than `D/4`.  Hence full `Omega(Dk log q)`
response information occurs in a genuinely non-Hamming carrier.

### Example CR.4 (Lee metric)

For a prime `p`, give `F_p^D` the Lee norm

```math
\|x\|_L=\sum_i\min(\bar x_i,p-\bar x_i).                    \tag{CR.10}
```

For all sufficiently large `D`, a random subspace argument supplies a host
of dimension `floor(D/4)` and minimum Lee weight greater than `D/8`.
Indeed a Lee ball of radius `D/8` is contained in the set of words with
Hamming support at most `D/8`, whose size is at most

```math
2^{H_2(1/8)D}p^{D/8}
\le p^{(H_2(1/8)+1/8)D},
```

and `H_2(1/8)+1/8<3/4`.  A fixed nonzero vector lies in a random
`floor(D/4)`-subspace with probability at most `p^{-3D/4+O(1)}`.  Therefore,
for `k\le D/32`, Theorem CR.2 gives at least

```math
p^{3Dk/16}
```

profiles separated by more than `D/16`.  For `p\ge5` this is not Hamming
geometry.

## 3. Exact and macroscopic collapse counterexamples

### Proposition CR.5 (redundant-alphabet exact collapse)

Give `F_q^D` the discrete metric

```math
d_{disc}(x,y)=\mathbf1_{x\ne y}.
```

Equivalently, every nonzero kernel element is a unit-cost generator.  Then,
for every `k` and every `V`,

```math
\boxed{F_V(u)=\mathbf1_{u\ne0}.}               \tag{CR.11}
```

Thus all `q^{Dk}` distinct mixed-holonomy gauge classes have one exact
kernel-endpoint response profile.  The fragments remain nonempty and
individually shear-trivial.

#### Proof

The choice `z=0` in (CR.2) costs at most one.  Every `z\ne0` pays relation
toll at least two, so it cannot improve the answer. `\square`

More generally, if the carrier diameter is `Delta`, then all endpoint
profiles lie in `[0,Delta]` pointwise and hence are mutually within `Delta`.
In fact the complete weighted kernel alphabet gives the stronger all-context
comparison

```math
|\ell_{P\cup R_V\cup E}(x)-
  \ell_{P\cup R_{V'}\cup E}(x)|\le\Delta       \tag{CR.12}
```

for every future alphabet `E` and root `x`: replace every used `R_V` letter
by its `R_{V'}` counterpart and correct the aggregate kernel discrepancy in
one weighted step.  Consequently a carrier of raw dimension `D` but
diameter `o(D)` has zero macroscopic response rate at additive scale
`epsilon D`, regardless of relation rank.

Diameter alone is still not sufficient.

### Theorem CR.6 (linear-diameter, full-cardinality, bounded response collapse)

Let `pi:F_q^D\twoheadrightarrow F_q^r` be linear and, for a scale `s>0`, set

```math
d_{pi,s}(x,y)
=s\mathbf1_{\pi x\ne\pi y}+\mathbf1_{x\ne y}.               \tag{CR.13}
```

This is a translation-invariant metric of diameter `s+1`.  For
`C=im V`, define the coarse profile

```math
G_{\pi C}(u)=(s+1)\mathbf1_{\pi u\notin\pi C}.               \tag{CR.14}
```

Then every mixed-holonomy profile satisfies

```math
\boxed{\ \|F_V-G_{\pi C}\|_\infty\le2k+1.\ }                \tag{CR.15}
```

There is also a stronger statement for the complete future-query family.
If

```math
\pi V=\pi V',
```

then, for every future alphabet `E` and every root `x`,

```math
|\ell_{P\cup R_V\cup E}(x)-
  \ell_{P\cup R_{V'}\cup E}(x)|\le1.                         \tag{CR.16}
```

Therefore all `q^{Dk}` exact gluing classes have an all-context response
quotient of size at most `q^{rk}` and error one.  For kernel-endpoint queries
alone, (CR.15) gives the still smaller quotient indexed by the subspace
`\pi C`; when `r=1` it has exactly two possible decoder profiles.  Taking
`s=D` and fixed `r` gives vanishing normalized distortion despite linear
diameter and exponential carrier cardinality.  The all-context charge is at
most `rk log_2q`, rather than `Dk log_2q`, with no restriction on `k`.

#### Proof

If `pi u\in pi C`, some `c\in C` has the same projection as `u`, so
`d_{pi,s}(u,C)` is zero or one.  If `pi u\notin pi C`, every `c\in C` is at
distance exactly `s+1`.  Hence

```math
\|d(\mathord\cdot,C)-G_{\pi C}\|_\infty\le1.
```

Combine this with (CR.3).

For (CR.16), replace each used `R_V` letter by the correspondingly labelled
`R_{V'}` letter.  Their aggregate discrepancy lies in `ker pi`, where the
metric diameter is one, and can be corrected by one kernel letter.  Reverse
the roles of `V,V'`. `\square`

If `r` is fixed, this is a constant-state macroscopic quotient.  More
generally the number of coarse states is at most `q^{O(r^2)}`; it is
subexponential in `Dk` whenever `r^2=o(Dk)`.  The obstruction is a metric
bottleneck: almost all carrier bits live inside fibres of diameter one and
are invisible at scale `s`.

## 4. Corrected composition law

The examples yield a precise three-factor test for response-information
growth in this family:

1. **relation rank** supplies the latent domain `F_q^k`;
2. **circuit toll** costs at most `2k` before holonomy can help; and
3. **linear carrier packing**, not `log|W|`, determines how many holonomy
   images remain distinguishable after that toll.

At scale `s_D` with `k=o(s_D)`, (CR.5) makes the third factor exact: the
response-information exponent is the Hausdorff packing exponent of
`Gr_k(W)` at that scale.  The separated-host condition of CR.2 is a simple
certificate for `Omega(Dk log q)` growth.  CR.5 and CR.6 show independently
that neither raw gauge dimension nor raw dimension plus linear diameter is
such a certificate.

This is not convex duality or a dynamic-programming restatement.  It
identifies which exact algebraic gluing coordinates survive the original
shortest-path query metric, replaces a false dimension law by a quantitative
metric-entropy equivalence, and predicts both collapse and full-rate growth
in non-Hamming models.
