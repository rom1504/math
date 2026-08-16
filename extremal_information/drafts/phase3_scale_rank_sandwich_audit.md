# Scale-rank sandwich and the metric Singleton obstruction

**Status.** The statements below are proved.  The strict-threshold issue,
the complete two-scale rank function, the quotient decoder, and sharpness of
the three error terms are checked in
[`verify_phase3_scale_rank_sandwich.py`](../experiments/verify_phase3_scale_rank_sandwich.py).

This note audits a proposed predictive law for multichannel response.  The
raw dimension of a carrier is not the relevant rank.  At a declared metric
scale, the lower rank is the largest dimension of a linearly separated host;
an upper rank is supplied by any metric synchronization quotient.  A
generalized Singleton inequality links the two.

The inequality is elementary rank--nullity, and in Hamming space is exactly
the classical puncturing proof of Singleton.  Its value here is not a new
coding bound.  It gives a falsifiable compatibility condition between the
lower carrier-capacity theorem and the upper synchronization theorem.

## 1. Separated linear rank gives the lower response rate

Let `W=F_q^D` have a translation-invariant metric with norm
`||w||=d(w,0)`.  For `Delta>=0`, define

```math
s_W(\Delta)=
\max\left\{
\dim C_0:C_0\le W,
\ \min_{c\in C_0\setminus\{0\}}\|c\|>\Delta
\right\}.                                                   \tag{SR.1}
```

The zero subspace is allowed and has infinite minimum distance.  The strict
inequality in (SR.1) is part of the definition.

For a linear map `V:F_q^k\to W`, consider the scalar-closed mixed-channel
profile

```math
F_V(u)=\min_{z\in F_q^k}
       \left(2\operatorname{wt}(z)+\|u+Vz\|\right).          \tag{SR.2}
```

### Theorem SR.1 (scale-rank packing)

Put `s=s_W(Delta)`.  For every `1<=k<=s`, there are at least

```math
{s\brack k}_q\ge q^{k(s-k)}                                \tag{SR.3}
```

profiles with pairwise uniform distance strictly greater than

```math
\Delta-2k.                                                   \tag{SR.4}
```

Consequently, if `2 epsilon<Delta-2k`, any deterministic summary answering
all kernel-endpoint queries to uniform error `epsilon` needs at least

```math
k(s_W(\Delta)-k)\log_2q                                    \tag{SR.5}
```

bits on this family.

#### Proof

Choose a host `C_0` attaining (SR.1), and one ordered basis for each
`k`-subspace `C<=C_0`.  For distinct `C,C'`, take
`c\in C\setminus C'`.  Every `c-c'`, `c'\in C'`, is a nonzero member of
`C_0`, so

```math
d_H(C,C')\ge d(c,C')>\Delta.                                \tag{SR.6}
```

The presentation cost in (SR.2) lies between zero and `2k`, hence the
carrier--response comparison gives

```math
\|F_C-F_{C'}\|_\infty
\ge d_H(C,C')-2k>\Delta-2k.
```

The Gaussian binomial coefficient counts the subspaces and satisfies the
displayed lower bound.  Profiles sharing an `epsilon`-accurate summary state
would be within `2epsilon`, proving (SR.5). `\square`

This lower theorem is useful only when the relation toll `2k` is below the
separation scale.  Large algebraic rank with `2k` already macroscopic is not
certified by (SR.1).

## 2. Synchronization rank gives the upper response rate

Let `Y=F_q^r` carry a metric `d_Y`, and let `pi:W\to Y` be a linear
surjection.  Recall that `pi` is an **`(a,b)` metric synchronization** when:

1. `pi` is one-Lipschitz;
2. every fibre has diameter at most `a`; and
3. for every `u\in W` and `y\in Y`, some `z\in pi^{-1}(y)` satisfies

   ```math
   d(u,z)\le d_Y(\pi u,y)+b.                                \tag{SR.7}
   ```

The third condition is the lift defect.  Merely bounding fibre diameters,
without one-Lipschitzness and (SR.7), does not justify the decoder below.

### Theorem SR.2 (projected-subspace summary)

For every `V:F_q^k\to W`, put `C=im V`, `U=pi C`, and

```math
G_U(u)=d_Y(\pi u,U).                                        \tag{SR.8}
```

Then

```math
\boxed{
0\le F_V(u)-G_U(u)\le a+b+2k
}\qquad(u\in W).                                            \tag{SR.9}
```

Thus the projected subspace `U`, using at most

```math
N_{r,k}(q)=
\sum_{j=0}^{\min\{k,r\}}{r\brack j}_q                     \tag{SR.10}
```

states, answers every endpoint query to error `a+b+2k`.  In bits,

```math
R^{det}_{a+b+2k}\le\log_2N_{r,k}(q).                        \tag{SR.11}
```

If `k<=r/2`, then, with

```math
K_q=\prod_{i=1}^{\infty}(1-q^{-i})^{-1},
```

```math
N_{r,k}(q)
\le(k+1)K_q q^{k(r-k)},                                    \tag{SR.12}
```

so the upper exponent is `k(r-k)log_2q+O(log k)`.

#### Proof

Write the presentation cost of `c\in C` as

```math
\alpha_V(c)=\min_{Vz=c}2\operatorname{wt}(z)\in[0,2k].
```

Then `F_V(u)=min_{c\in C}(d(u,c)+alpha_V(c))`.  One-Lipschitzness gives
`G_U<=F_V`.  Conversely, choose `y\in U` nearest to `pi u`, choose
`c\in C` with `pi c=y`, and use (SR.7) to obtain `z` in the same fibre as
`c` with

```math
d(u,z)\le G_U(u)+b.
```

The fibre bound gives `d(z,c)<=a`, and the presentation costs at most `2k`.
This proves (SR.9).  Since `dim(pi C)<=min(k,r)`, the state count is (SR.10).

Finally

```math
{r\brack j}_q
\le K_q q^{j(r-j)}.
```

For `j<=k<=r/2`, the exponent `j(r-j)` is increasing in `j`; summing proves
(SR.12). `\square`

The decoder is an endpoint summary.  It becomes a reusable feature algebra
only for compositions whose projected carriers update without reconstructing
the discarded fibres.

## 3. Generalized Singleton inequality

### Theorem SR.3 (metric Singleton obstruction)

If a linear map `pi:W\to F_q^r` has fibres of diameter at most `a`, then

```math
\boxed{s_W(a)\le r.}                                         \tag{SR.13}
```

In particular every `(a,b)` metric synchronization has rank at least
`s_W(a)`, independently of its lift defect.

#### Proof

Let `C_0` have minimum nonzero norm strictly greater than `a`.  If
`dim C_0>r`, rank--nullity gives a nonzero
`c\in C_0\cap ker pi`.  Both `c` and zero lie in one fibre, so
`||c||<=a`, a contradiction.  Therefore `pi|_{C_0}` is injective and
`dim C_0<=r`.  Maximize over `C_0`. `\square`

The strict threshold cannot be replaced by `>=a`.  For example, in the
two-scale metric below, `ker pi` has minimum nonzero norm exactly one and can
have dimension much larger than `r`, although every fibre has diameter one.

### Classical puncturing

In `F_q^D` with Hamming distance, puncture `h` coordinates.  The resulting
map has rank `D-h`, fibre diameter `h`, and lift defect zero: retain the
deleted coordinates of the query when lifting a quotient target.  Hence

```math
s_W(h)\le D-h.                                                \tag{SR.14}
```

If a linear `[D,K,d]` code has minimum distance `d`, take `h=d-1`.  Since
`d>d-1`, (SR.14) gives

```math
K\le D-d+1,                                                  \tag{SR.15}
```

which is exactly the classical Singleton bound and its puncturing proof.
Thus (SR.13) is the metric-fibre form of a familiar theorem, not a new proof
of Singleton.

## 4. Exactness on the two-scale carrier

Let `pi_0:F_q^D\twoheadrightarrow F_q^r`, let `L>0`, and put

```math
\|w\|_{2sc}
=L\mathbf1_{\pi_0w\ne0}+\mathbf1_{w\ne0}.                  \tag{SR.16}
```

Give `F_q^r` the scaled discrete metric

```math
d_Y(y,y')=(L+1)\mathbf1_{y\ne y'}.
```

Then `pi_0` is a `(1,0)` metric synchronization.  More strongly, the entire
separated-rank curve is

```math
s_W(\Delta)=
\begin{cases}
D,&0\le\Delta<1,\\
r,&1\le\Delta<L+1,\\
0,&\Delta\ge L+1.
\end{cases}                                                  \tag{SR.17}
```

#### Proof

The whole space has minimum nonzero norm one, proving the first case.  In the
middle range, SR.3 gives the upper bound `r`; any linear section of `pi_0`
has dimension `r` and every nonzero element has norm `L+1`, proving equality.
No nonzero vector has norm greater than `L+1`, proving the last case.
`\square`

Take `Delta=L` and `k<=r`.  The lower theorem gives

```math
{r\brack k}_q
```

profiles separated by more than `L-2k`.  The upper theorem uses at most
`sum_{j<=k}{r bracket j}_q` states at error `2k+1`.  When
`k<=r/2` and `k=o(L)`, the lower and upper description exponents both equal

```math
k(r-k)\log_2q+O(\log k).                                    \tag{SR.18}
```

Thus the sandwich is exact to leading exponential order in a model where
the raw gluing space has `Dk log_2q` bits but only `rk log_2q` metric-visible
directions.  The rank `r` is recovered from the response scale, rather than
inserted as a cardinality parameter.

## 5. Constant audit

None of the three universal terms in `a+b+2k` can simply be deleted.

* **Fibre term `a`.**  Let `X=Y times Z` with the `ell_1` product metric,
  project to `Y`, and take a singleton carrier at one end of a diameter-`a`
  fibre while querying the other end.
* **Lift term `b`.**  Take the identity map from a two-point space of distance
  `1+b` to one of distance one.  Fibres are singletons, but a cross-point
  query loses exactly `b`.
* **Presentation term `2k`.**  In binary Hamming space, let the `k` columns
  be indicators of disjoint three-coordinate blocks and query their sum.
  Its distance to the image code is zero, while the cheapest presented
  realization uses every column and costs exactly `2k`.

The sum `a+b+p` is sharp for general presented carriers: use the `ell_1`
product of the two-point lift-defect example and a diameter-`a` fibre, place
a singleton carrier at the opposite corner, and assign it presentation cost
`p`.

These examples do not prove that every specialized multichannel model
simultaneously attains all three terms.  They prove that no improvement is
available from the metric-synchronization axioms and the radius bound alone.

## 6. Judgment

The separated rank is not an independent new invariant pasted onto the
response: SR.1 shows that it produces actual mixed-fragment packings.  The
synchronization rank is likewise operational through the explicit decoder
SR.2.  The genuinely informative statement is their compatibility SR.3:

```math
\text{small metric fibres at scale }a
\quad\Longrightarrow\quad
\text{quotient rank at least }s_W(a).                         \tag{SR.19}
```

This can falsify a proposed compressed quotient before any response decoder
is built.  On the two-scale model it predicts the exact effective exponent,
and in Hamming space it becomes Singleton.  Nevertheless the proof is
classical linear algebra, and `s_W` itself is a scale-dependent coding
parameter.  The result is a useful two-sided law for presented linear
carriers, not yet a universal theory for arbitrary extremal landscapes.

## 7. Independent audit of the MRD continuation

This section audits Sections 4--5 of
[`phase3_scale_rank_response_sandwich.md`](phase3_scale_rank_response_sandwich.md),
including its Gaussian state bound and rank-metric constants.

### Gaussian upper constant

For `0<=j<=r`, the companion draft uses

```math
{r\brack j}_q
=q^{j(r-j)}
 \prod_{h=1}^j
 {1-q^{-(r-j+h)}\over1-q^{-h}}
\le4q^{j(r-j)}.                              \tag{SR.A1}
```

This is correct uniformly in every prime power `q`.  The numerator factors
are at most one, while

```math
\prod_{h=1}^j(1-q^{-h})
\ge\prod_{h=1}^{\infty}(1-2^{-h})>{1\over4}. \tag{SR.A2}
```

For a fully elementary last inequality, retain the first two factors and use
`prod_i(1-x_i)>=1-sum_i x_i` on the tail:

```math
\prod_{h=1}^{\infty}(1-2^{-h})
\ge {1\over2}{3\over4}
    \left(1-\sum_{h=3}^{\infty}2^{-h}\right)
={9\over32}>{1\over4}.                       \tag{SR.A3}
```

There are `min(r,k)+1` summands, and each exponent is at most
`g(r,k)`, so the logarithmic bound SR.8, including the factor
`4(min(r,k)+1)`, is valid.  No missing factor of two occurs.

### Linearized-polynomial host

Let

```math
f(x)=\sum_{i=0}^{r-1}a_ix^{q^i}\in F_{q^D}[x],
\qquad1\le r\le D.                           \tag{SR.A4}
```

If the coefficients are not all zero, this is a nonzero ordinary polynomial
of degree at most `q^(r-1)`.  Its roots in `F_(q^D)` are the kernel of an
`F_q`-linear map.  If that kernel has dimension `m`, it contains exactly
`q^m` elements.  The ordinary polynomial root bound therefore gives

```math
q^m\le q^{r-1},\qquad m\le r-1,               \tag{SR.A5}
```

and rank at least `D-r+1`.

The coefficient map is injective as well.  A nonzero polynomial in (SR.A4)
cannot vanish on all `q^D` field elements, because

```math
q^{r-1}<q^D.
```

Hence the host has exactly `rD` dimensions over `F_q`.  There is no
polynomial-function collision hidden modulo `x^(q^D)-x`: all displayed
degrees are strictly below `q^D`.  Combining dimension `rD` with the
rank-metric Singleton inequality also shows that the minimum distance is
exactly `D-r+1`, although the response theorem only needs the proved lower
bound.

### Response packing and strict gaps

For two distinct `k`-subspaces `C,C'` of the host, choose
`c\in C\setminus C'`.  Every `c-c'`, `c'\in C'`, is a nonzero host word,
so its rank is at least `D-r+1`.  At query `c`, the profile from `C` costs at
most `2k`, while the profile from `C'` costs at least `D-r+1`.  Thus the
actual weak separation is

```math
\|F_C-F_{C'}\|_\infty\ge D-r+1-2k.            \tag{SR.A6}
```

This justifies the companion draft's “at least” formulation.  Its strict
form follows whenever the right-hand side is compared with a strictly
smaller real threshold.  Equivalently, apply the scale-rank theorem with
any `Delta<D-r+1`; one must not silently replace the strict host condition
by `>=Delta` at the endpoint.

Take `r=floor(D/2)` and `k<=D/16`.  Then

```math
D-r+1-2k
\ge {D\over2}+1-{D\over8}
>{3D\over8},                                  \tag{SR.A7}
```

so the claimed strict gap is correct for both parities of `D`.  Also

```math
rD-k
\ge {D(D-1)\over2}-{D\over16}
\ge {D^2\over3}                               \tag{SR.A8}
```

for every integer `D>=4`.  The Gaussian count is consequently at least

```math
q^{k(rD-k)}\ge q^{kD^2/3}.                    \tag{SR.A9}
```

Two profiles decoded from one state to error `epsilon D` would be within
`2epsilon D`.  Since the separation is strictly greater than `3D/8`, every
fixed

```math
\epsilon<3/16                                  \tag{SR.A10}
```

distinguishes the full family.  The information lower bound
`(1/3)kD^2 log_2q` follows with no asymptotic rounding loss.

The hypotheses are compatible: `k<=D/16` implies
`2k<D-r+1`.  The result retains `Theta(D^2k log q)` of the raw holonomy
information while the metric response scale remains `Theta(D)`.

### Computational cross-check and verdict

The companion verifier independently constructs
`G_2 subset End(F_16)`, obtaining dimension eight, `256` maps, and minimum
rank distance three.  It checks all `32,385` pairs of its `255`
one-dimensional subspaces and finds the predicted minimum response gap one.
The verifier passes on rerun.

**Verdict.** The Gabidulin/MRD continuation and all constants in SR.14--SR.17
are correct.  Unlike the earlier equilateral multiplication host, this family
uses a positive-rate MRD space of dimension `Theta(D^2)` and so genuinely
tests rank-metric Singleton geometry.  The theorem is generative within the
presented-carrier class, though its algebraic input is the classical
linearized-polynomial/Gabidulin mechanism rather than a new MRD construction.

## 8. Independent audit of the optimal quotient rank and code--anticode gap

This section audits SR.6--SR.8 of the companion draft.

### Canonical quotient is an exact submetry

Let `K<=W` and define

```math
\bar d(x+K,y+K)=\min_{k\in K}\|x-y+k\|.                     \tag{SR.A11}
```

This is well defined under changing either coset representative.  It is
positive on distinct cosets because the finite set `x-y+K` then excludes
zero; symmetry follows from `K=-K`; and the triangle inequality follows by
adding minimizing representatives.  Thus it is a genuine
translation-invariant metric, not merely a pseudometric.

The quotient map is one-Lipschitz by taking `k=0`.  If
`diam(K)<=a`, every fibre, being a translate of `K`, has diameter at most
`a`.  Finally, for a query `x` and target coset `y+K`, choose a minimizer
`k_*` in (SR.A11) and put `z=y-k_*`.  Then

```math
z\in y+K,
\qquad
d_W(x,z)=\|x-y+k_*\|=\bar d(x+K,y+K).                       \tag{SR.A12}
```

Hence the lift defect is exactly zero.  Conversely, the kernel of any
rank-`r` linear synchronization quotient has dimension `N-r` and is one
fibre, so it is an admissible linear anticode.  Both inequalities in

```math
q_W(a)=N-A_W(a)                                              \tag{SR.A13}
```

are therefore correct.  Combining this equality with the strict-distance
Singleton obstruction gives `s_W(a)+A_W(a)<=N` without circularity.

### Floors and endpoint checks

For Hamming space, a `d`-dimensional linear subspace has diameter at least
`d`: in reduced row-echelon form, the sum of all basis rows is nonzero in
each of the `d` pivot coordinates.  A coordinate subspace of dimension
`floor(a)` has diameter `floor(a)<=a`.  Therefore

```math
A_W(a)=\lfloor a\rfloor\qquad(0\le a\le D).                 \tag{SR.A14}
```

The resulting statement for the separated rank is only

```math
s_W(a)\le D-\lfloor a\rfloor,                                \tag{SR.A15}
```

not equality in general.  Equality in a nontrivial parameter range is
equivalent to the existence of the corresponding linear MDS code.  The
companion draft makes this distinction correctly.  At `a=D`, the zero code
and whole-space anticode give `s_W(D)=0`, `A_W(D)=D`; the usual MDS wording
should be understood with this trivial endpoint separated out.

For the two-scale metric, a subspace of diameter strictly below `L+1` must
lie in `ker pi`, while that kernel has dimension `D-r` and diameter one.
Thus `A_W(a)=D-r` on the correctly half-open interval
`1<=a<L+1`; at `a=L+1` the whole space becomes admissible.

For `D x D` rank metric and integer `0<=a<=D`, the fixed-`a`-row space has
dimension `Da` and diameter at most `a`.  The Gabidulin host with parameter
`D-a` has dimension `D(D-a)` and minimum rank at least `a+1`, so the two
lower bounds sum to `D^2`.  The code--anticode inequality forces

```math
A_W(a)=Da,
\qquad
s_W(a)=D(D-a).                                                \tag{SR.A16}
```

For `a=D`, interpret the host as the zero subspace rather than invoking the
earlier notation `G_0`, which was defined only implicitly.  This is a minor
endpoint convention, not a gap in the result.

### Binary Hamming sphere-packing gap

Let `a_D=floor(delta D)` and let a code counted by `s_W(a_D)` have minimum
distance `d_D>=a_D+1`.  With

```math
t_D=\left\lfloor{a_D\over2}\right\rfloor,
```

one has `2t_D<d_D` for both parities of `a_D`.  Hence the Hamming balls of
radius `t_D` are disjoint.  Since `t_D/D -> delta/2<1/2`, the standard ball
asymptotic gives

```math
{1\over D}\log_2\sum_{j=0}^{t_D}{D\choose j}
\longrightarrow H_2(\delta/2).                              \tag{SR.A17}
```

The Hamming bound implies

```math
{s_W(a_D)\over D}
\le1-H_2(\delta/2)+o(1).
```

Together with `A_W(a_D)=a_D`, this proves

```math
\liminf_D{\gamma_W(a_D)\over D}
\ge H_2(\delta/2)-\delta.                                   \tag{SR.A18}
```

The right side is strictly positive for every fixed `0<delta<1`: putting
`x=delta/2`, strict concavity of `H_2` above the chord from `(0,0)` to
`(1/2,1)` gives `H_2(x)>2x=delta`.  The floors cause only `o(D)` error and
do not alter the radius or positivity argument.

### Verdict on SR.6--SR.8

All three results are correct, including the exact zero lift defect, the
real-`a` Hamming floor, the sphere-packing radius, and the strict positivity
constant.  SR.6 is the substantive bridge: it identifies the optimal exact
linear synchronization rank with anticode codimension.  SR.8 then proves
that separated-host rank and synchronization rank can have a leading linear
gap, so the scale-rank sandwich is not a tautological duality.  What remains
open, as the companion draft says, is whether unrestricted Grassmannian
packings fill part of that gap.
