# Wave 50B: joint box mass and adaptive triple retention

## Status

No asymptotic box theorem, coarea theorem, or convergence proof is obtained.
There are two concrete advances:

1. the row-truncated triple condition is reduced to one exact weighted
   intersection-moment inequality, and its pressure is exactly linear under
   averaging over common-core scales; and
2. an explicit exact order-ten minimizer at the fixed ratio `m/n=3/5` has
   the smallest possible partial-completion box value, but its corresponding
   minimum-row class fails triple retention at **every** positive common-core
   scale.

Thus a box witness alone cannot force a common successful scale at the same
row cap if one insists on a pure positive-core kernel.  Averaging only those
positive-core scales does not repair their negative overlap profile.
However, the exact mixture `3/4 K_0+1/4 K_1` already passes at cap `10` and
has `kappa=1/53`.  A viable joint theorem may therefore either enlarge the
row class and prove pure-scale overlap abundance, or retain an independent-
resampling `K_0` component while keeping the mixed spectral separation
bounded below.  The example is finite and does not falsify the asymptotic
route; pure scale `ell=2` also passes after enlarging the cap to `138`.

Every assertion below is checked exhaustively and with exact rational
arithmetic by `tmp/box_triple_r50_check.py`; the captured output is
`tmp/box_triple_r50.out`.

## 1. Exact overlap-moment form and scale averaging

Let `N=binom(n,m)`, let `C_R={z:R_2(z)<=R}`, and for the zero-deficit ground
family put

```math
F_z=\{S:D_z(S)=0\},\qquad r_z=|F_z|,\qquad a_z=r_z/N.
```

For the common-core kernel `K_ell`, write

```math
d_\ell=\binom m\ell\binom{n-\ell}{m-\ell},
\qquad
K_\ell(S,T)=\frac{\binom{|S\cap T|}{\ell}}{d_\ell}.
```

For a positive-degree center define its one-step retention

```math
p_\ell(z)
=\Pr\{T\in F_z\mid S\in F_z,\ T\sim K_\ell(S)\}
=\frac{1}{r_zd_\ell}
\sum_{S,T\in F_z}\binom{|S\cap T|}{\ell}.
\tag{R50B.1}
```

In the triple experiment, the two independent favorable incidences
size-bias the center by `a_z^2`.  Consequently its exact retention is

```math
P_\ell(R)
=\frac{\sum_{z\in C_R}a_z^2p_\ell(z)}
       {\sum_{z\in C_R}a_z^2}
=\frac{
 \sum_{z\in C_R}r_z\sum_{S,T\in F_z}\binom{|S\cap T|}{\ell}}
 {d_\ell\sum_{z\in C_R}r_z^2}.
\tag{R50B.2}
```

Provided the denominator is positive, (10.1235) is therefore exactly the
weighted overlap inequality

```math
\boxed{
\sum_{z\in C_R}r_z\sum_{S,T\in F_z}
 \binom{|S\cap T|}{\ell}
\ \ge\
\lambda_1(\ell)d_\ell\sum_{z\in C_R}r_z^2.}
\tag{R50B.3}
```

This exposes precisely what box mass does and does not give.  A box witness
with row at most `R` makes `sum r_z^2>0`; it supplies no lower bound on the
overlap moment in the left side of (R50B.3).

There is also an exact scale-averaging identity.  Fix a nonempty set
`L subseteq {1,...,m-1}` and weights `w_ell>=0` with
`sum_(ell in L) w_ell=1`.  Define

```math
G_\ell(R)=\mathbb E_z\!\left[
 \mathbf1_{C_R}a_z^2\{p_\ell(z)-\lambda_1(\ell)\}\right].
```

and put

```math
K_w=\sum_{\ell\in L}w_\ell K_\ell,
\qquad
\lambda_{j,w}=\sum_{\ell\in L}w_\ell\lambda_j(\ell),
\qquad
p_w(z)=\sum_{\ell\in L}w_\ell p_\ell(z).
```

The mixed boundary is `B_{z,w}=<f_z,(I-K_w)f_z>` and its normalization is
`delta_w=1-lambda_(1,w)`.  Repeating the derivation of (10.1234)--(10.1235)
gives

```math
\frac{\mathbb E[\mathbf1_{C_R}a_zB_{z,w}]}
{\delta_w\mathbb E[\mathbf1_{C_R}a_z^2]}\le1
\quad\Longleftrightarrow\quad
P_w(R):=\frac{\sum_{z\in C_R}a_z^2p_w(z)}
{\sum_{z\in C_R}a_z^2}\ge\lambda_{1,w}.
```

Consequently

```math
\boxed{G_w(R)=\sum_\ell w_\ell G_\ell(R).}
\tag{R50B.4}
```

Indeed both the retention and the level-one eigenvalue are linear in the
kernel.  If every component has `G_ell(R)<0`, then every such mixture also
has `G_w(R)<0` and fails its correctly normalized mixed criterion.  More
generally, `G_w>=0` implies at least one component has `G_ell>=0`.  Hence an
average-over-scales proof can certify that some scale passes only by
establishing nonnegative average pressure; it cannot hide failure at every
component.  This statement concerns mixtures supported on the candidate
positive-core scales.  If `ell=0` is included, its pressure must be included
too, so componentwise negative pressure is no longer the premise;
`K_0` itself has `kappa=0`, but an admixture can have positive `kappa`.

This admixture has an exact two-parameter form.  Independent resampling gives
`p_0(z)=a_z`, so put

```math
\bar a_R=P_0(R)=
\frac{\sum_{z\in C_R}a_z^3}{\sum_{z\in C_R}a_z^2}.
```

For `K_{theta,ell}=(1-theta)K_0+theta K_ell`,

```math
\boxed{
\begin{aligned}
P_{\theta,\ell}(R)&=(1-\theta)\bar a_R+\theta P_\ell(R),\\
\lambda_{1,\theta,\ell}&=\theta\lambda_1(\ell),\qquad
\lambda_{2,\theta,\ell}=\theta\lambda_2(\ell),\\
\kappa_{\theta,\ell}
&=\frac{\theta\{\lambda_1(\ell)-\lambda_2(\ell)\}}
        {1-\theta\lambda_1(\ell)}.
\end{aligned}}
\tag{R50B.5}
```

Thus the exact mixed target is
`P_(theta,ell)(R)>=theta lambda_1(ell)` together with
`kappa_(theta,ell)>=kappa_0`.  The `K_0` reserve can pay a negative
positive-core pressure, but taking `theta` too small destroys the spectral
separation.  This tradeoff, not componentwise scale selection alone, is the
correct surviving averaged formulation.

## 2. Exact fixed-ratio obstruction at the attained box cap

Consider the explicit normalized signing

```text
 0  1  1  1  1  1  1  1  1  1
 1  0 -1  1 -1  1  1  1 -1 -1
 1 -1  0  1  1  1 -1  1  1  1
 1  1  1  0  1  1 -1 -1 -1 -1
 1 -1  1  1  0 -1  1 -1 -1  1
 1  1  1  1 -1  0  1  1  1 -1
 1  1 -1 -1  1  1  0 -1  1 -1
 1  1  1 -1 -1  1 -1  0 -1  1
 1 -1  1 -1 -1  1  1 -1  0  1
 1 -1  1 -1  1 -1 -1  1  1  0
```

Exhausting its `512` projective signings gives `Q(A)=26`.  Since the proved
order-ten optimum is `q_10=26`, this is an exact minimizer.  Set `m=6`.

The exhaustive partial-completion calculation gives

```math
\boxed{
\min_{S,\ y\text{ a ground of }A[S]}\mathcal V(S,y)=10.}
\tag{R50B.6}
```

One witness is

```text
S=(0,1,2,6,7,8),       y=(1,1,-1,-1,-1,1),
S^c=(3,4,5,9),         w=(1,-1,-1,1).
```

Here `Q(A[S])=10`, the completion has parent energy `-10`, and its row is
`R_2=10`.  This value is optimal even among all parent signings: every
coordinate of `Az` is odd, hence `R_2(z)>=10`.

The class `C_10` consists of exactly two oriented projective centers,

```text
(1,1,-1,-1, 1,-1, 1, 1,-1,1),
(1,1,-1, 1,-1,-1,-1,-1, 1,1).
```

Their parent energies are respectively `10` and `-10`.  Both have the same
five-member ground family

```text
012678, 014679, 034579, 123568, 234589.
```

Thus each has degree `5/binom(10,6)=1/42`, and under the uniform projective
center law

```math
D_{C_{10}}=\frac1{451584}>0.
```

Among the `25` ordered pairs in this family, five have intersection six,
ten have intersection four, and ten have intersection two.  Substitution
in (R50B.1) gives the following exact table.  The `coarea` column is
`(1-p_ell)/(1-lambda_1)` and must be at most one; equivalently the
`p/lambda_1` column must be at least one.

| `ell` | `lambda_1` | `lambda_2` | `p_ell` | `p_ell-lambda_1` | `p_ell/lambda_1` | coarea | `kappa` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | `2/27` | `0` | `1/42` | `-19/378` | `9/28` | `369/350` | `2/25` |
| 2 | `1/6` | `1/70` | `29/1050` | `-73/525` | `29/175` | `1021/875` | `32/175` |
| 3 | `2/7` | `2/35` | `1/25` | `-43/175` | `7/50` | `168/125` | `8/25` |
| 4 | `4/9` | `4/25` | `17/225` | `-83/225` | `17/100` | `208/125` | `64/125` |
| 5 | `2/3` | `2/5` | `1/5` | `-7/15` | `3/10` | `12/5` | `4/5` |

Here `G_ell(10)=D_(C_10)(p_ell-lambda_1(ell))`, so the table gives every
pressure with its exact common factor `D_(C_10)=1/451584`.  Every pressure
is strictly negative.  By (R50B.4), every mixture with nonnegative weights
summing to one and supported on these positive-core scales is also negative
under the mixed normalization `delta_w=1-lambda_(1,w)`.  This is
not a nearest-core-only failure: it includes all possible `ell=1,...,5`,
and `ell=2,...,5` have `kappa>=32/175`.

The `K_0` caveat is active, not cosmetic.  At the same cap `10`, take

```math
K_w=\frac34K_0+\frac14K_1.
```

Here `p_0=p_1=1/42`, and exact substitution gives

```math
\boxed{
p_w=\frac1{42},\quad
\lambda_{1,w}=\frac1{54},\quad
p_w-\lambda_{1,w}=\frac1{189},\quad
\lambda_{2,w}=0,\quad
\kappa_w=\frac1{53}.}
\tag{R50B.7}
```

The normalized coarea ratio is `369/371<1`, and the unnormalized pressure is
`G_w=D_(C_10)/189=1/85349376>0`.  Hence this mixed kernel passes at the
minimum cap despite failure of every pure positive-core scale and of every
mixture supported only on those scales.

The full row-cap phase diagram identifies a second, pure-scale escape in
this one finite exact-minimizer example.  Its active row levels are

```text
10, 42, 74, 106, 138.
```

Scale one first passes at cap `106`, scale two first passes at cap `138`,
and scales three through five never pass, even at the full active cap.
Thus it does not falsify existence at a project-scale cap in an asymptotic
family, cap enlargement, or mixed kernels.  It rigorously rules out only
the pure-scale same-witness/minimum-cap implication

```text
one low-row child-ground completion
    => triple retention for some pure positive-core K_ell at that cap.
```

## 3. Sharply reduced joint lemma

The natural surviving joint statement should permit both controlled cap
inflation and `K_0` admixture.
Fix a density window `p=m/n in [p_0,p_1]` and a constant `alpha>0`, and
restrict to scales

```math
\alpha n\le\ell\le m-\alpha n.
```

At these scales the exact formula

```math
\kappa(\ell)=
\frac{\ell(n-m)(n-2)}{n(m-1)(n-\ell-1)}
```

is bounded below by a positive constant depending only on
`p_0,p_1,alpha`.  A concrete reduced target is:

> **Cap-inflated mixed-overlap lemma.** There are constants
> `C,c,alpha,kappa_0>0`, uniform on the density window, such that whenever an
> exact minimizer has a child-ground box witness of row
> `R_0=O(n^(9/4-c))`, there are a cap
> `R_0<=R<=C(R_0+n^2)`, a scale `ell` in the displayed balanced range, and
> `theta in (0,1]` such that
> `P_(theta,ell)(R)>=theta lambda_1(ell)` and
> `kappa_(theta,ell)>=kappa_0`.

Because the original witness remains in `C_R`, the denominator is positive.
The inflated cap is still `O(n^(9/4-c))`; the first mixed inequality is
exactly non-strict coarea for `K_(theta,ell)` by (R50B.4)--(R50B.5); and the
second lets the verified slice-FKN step produce the required constant-degree
project-row center.  Thus this lemma, together with the still-open
asymptotic box bound, proves (10.795) and convergence.  The former pure
cap-inflated lemma (R50B.3) is the special case `theta=1` and remains a valid
but unnecessarily strong target.

This formulation also has a precise falsification criterion: an unbounded
exact-minimizer family with project-scale box witnesses but negative mixed
pressure at every project-scale cap, every balanced common-core scale, and
every `theta` for which `kappa_(theta,ell)>=kappa_0`.  The order-ten example
does not supply such a falsifier: its explicit `K_0` admixture already
passes.

## 4. Research judgment

The proposed pure common-scale proof does not follow from box discrepancy,
and a replacement-flow argument based only on the selected low-row
completion is blocked by its sparse five-set ground family.  But the same
family has enough independent-resampling reserve to pass after a `K_0`
admixture.  The useful next question is whether, after controlled cap
inflation, either a pure balanced scale has nonnegative pressure or the
reserve `bar a_R` can pay its pressure deficit while `theta` stays large
enough to keep `kappa_(theta,ell)` uniformly positive.  Equations
(R50B.2)--(R50B.5) are direct finite statistics for this test; the
cap-inflated mixed-overlap lemma is the sharp surviving interface.
