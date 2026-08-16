# Scale-rank sandwich for multichannel carrier responses

**Status.** The statements below are proved.  They turn the carrier-capacity
and metric-synchronization theorems into a two-sided algebraic law for linear
metric carriers.  The rank-metric application uses a self-contained Gabidulin
host and is genuinely sensitive to rank geometry.

The carrier law says that response complexity is Hausdorff carrier entropy
minus presentation radius.  This note supplies two finite invariants that can
be checked without enumerating all carriers:

* separated linear rank certifies response growth;
* synchronizing quotient rank certifies response compression.

They obey a generalized Singleton inequality.  This is the first theorem in
the program that relates a lower-capacity certificate and an upper-compression
certificate directly.

## 1. Linear metric carriers

Let `W` be a finite-dimensional vector space over `F_q` with a
translation-invariant metric `d`, and write `||w||=d(w,0)`.  For `Delta>=0`
define the **separated linear rank**

```math
s_W(\Delta)=\max\left\{
 \dim C_0:C_0\le W,
 \min_{c\in C_0\setminus\{0\}}\|c\|>\Delta
 \right\}.                                      \tag{SR.1}
```

The zero subspace is allowed, with its empty minimum interpreted as
`+infinity`.

For an injective map `V:F_q^k->W`, define the multichannel profile

```math
F_V(u)=\min_{z\in F_q^k}
 \{2\operatorname{wt}(z)+\|u+Vz\|\}.           \tag{SR.2}
```

It is the endpoint response of two locally shear-trivial scalar-closed
fragments over the same `k` quotient directions.

### Theorem SR.1 (scale-rank response sandwich)

**Lower certificate.**  If `s=s_W(Delta)>=k` and `Delta>2k`, then there are at least

```math
q^{k(s-k)}                                      \tag{SR.3}
```

profiles with pairwise uniform distance greater than `Delta-2k`.  Hence, for

```math
2\varepsilon<\Delta-2k,                        \tag{SR.4}
```

uniform error `epsilon` requires at least

```math
k(s-k)\log_2q                                  \tag{SR.5}
```

deterministic bits on this family.

**Upper certificate.**  Suppose a linear surjection

```math
\varpi:W\longrightarrow Y\cong F_q^r
```

is an `(a,b)` metric synchronization in the sense of Theorem MQ.1.  Then every
profile (SR.2) is decoded from the single subspace `varpi(im V)` to error

```math
a+b+2k.                                        \tag{SR.6}
```

The number of possible decoder states is at most

```math
N_q(r,k)=\sum_{j=0}^{\min\{r,k\}}{r\brack j}_q.               \tag{SR.7}
```

Writing

```math
g(r,k)=\max_{0\le j\le\min\{r,k\}}j(r-j),
```

one has the explicit description bound

```math
\log_2N_q(r,k)
\le g(r,k)\log_2q+\log_2(4(\min\{r,k\}+1)).    \tag{SR.8}
```

#### Proof

Choose a dimension-`s` subspace `C_0` of minimum nonzero weight greater than
`Delta`, and one ordered basis for every `k`-subspace `C<=C_0`.  Distinct
subspaces have Hausdorff distance greater than `Delta`: a point of
`C\setminus C'` is separated from every point of `C'` by a nonzero member of
`C_0`.  The presentation cost in (SR.2) is at most `2k`, so the carrier law
gives response separation greater than `Delta-2k`.  The Gaussian binomial
count is at least `q^(k(s-k))`, proving the lower assertions.

For the upper assertion, `im V` is a presented carrier of radius `2k`.
Metric-quotient synchronization decodes it from `varpi(im V)` with error
`a+b+2k`.  Its projected dimension is at most `min(r,k)`, giving (SR.7).
Finally,

```math
{r\brack j}_q
=q^{j(r-j)}
 \prod_{h=1}^j{1-q^{-(r-j+h)}\over1-q^{-h}}
\le4q^{j(r-j)},                                \tag{SR.9}
```

because `prod_(h>=1)(1-2^(-h))^(-1)<4`; summing proves (SR.8). `square`

The theorem does not assert that one certificate is always sharp.  It gives
two checkable routes with their exact error scales.

## 2. The generalized Singleton inequality

### Theorem SR.2 (separation cannot hide inside a small fibre)

Let `varpi:W->F_q^r` be any linear map whose fibres have metric diameter at
most `a`.  Then

```math
\boxed{s_W(a)\le r.}                            \tag{SR.10}
```

No lifting hypothesis is needed.

#### Proof

Let `C_0` have minimum nonzero weight greater than `a`.  If a nonzero
`c in C_0` lay in `ker varpi`, then `c` and zero would be two points of one
fibre at distance greater than `a`, a contradiction.  Thus `varpi` is
injective on `C_0`, and `dim C_0<=r`. `square`

This elementary statement becomes classical coding theory in standard
metrics.

* For Hamming `W=F_q^D`, puncturing `h` coordinates has fibre diameter `h`
  and target dimension `D-h`.  Taking `h=d-1` yields the Singleton bound
  `dim C<=D-d+1` for a code of minimum distance `d`.
* For `D x D` matrices in rank metric, retaining `r` rows has fibre diameter
  `D-r` and target dimension `rD`.  Thus a rank-metric code of minimum rank
  distance `d=D-r+1` has dimension at most `rD=D(D-d+1)`, the rank-metric
  Singleton bound.

The theorem therefore identifies synchronizing quotient rank as the dual
upper obstruction to separated carrier rank, rather than introducing an
unrelated compression parameter.

## 3. Exactness on the two-scale collapse

Let `varpi:F_q^D->F_q^r` be onto and set, for `L>0`,

```math
d_L(x,y)=L\mathbf1_{\varpi x\ne\varpi y}
          +\mathbf1_{x\ne y}.                  \tag{SR.11}
```

The fibres have diameter one and the quotient is a `(1,0)` synchronization
onto the scaled discrete metric.

### Proposition SR.3 (matching scale ranks)

```math
s_W(\Delta)=
\begin{cases}
D,&0\le\Delta<1,\\
r,&1\le\Delta<L+1,\\
0,&\Delta\ge L+1.
\end{cases}                                    \tag{SR.12}
```

For `k<=r` with `2k<L+1`, the lower certificate supplies
`q^(k(r-k))` profiles separated by at least `L+1-2k`, while the upper certificate stores at most all subspaces
of `F_q^r` and has error `1+2k`.  When `k=o(L)`, the macroscopic response
quotient is therefore exactly the Grassmannian of the coarse quotient, up to
subscale error.

#### Proof

Below one, all of `W` has minimum nonzero weight one.  At scale at least one,
Theorem SR.2 gives `s_W(Delta)<=r` as long as `Delta>=1`.  A linear section
of `varpi` has dimension `r`, and every nonzero point in it changes the
quotient, hence has weight `L+1`; it works for every `Delta<L+1`.  At or
above the diameter no nonzero subspace qualifies.  The remaining claims are
SR.1 applied to the section and to the quotient. `square`

This explains the counterexample quantitatively: fine `D-r` carrier
coordinates are not merely absent from one decoder; no linearly separated
host at scale above the fibre diameter can use them.

## 4. An intrinsically rank-metric full-rate family

The earlier multiplication-host example lived in rank metric but used only
an equilateral `D`-dimensional subspace.  A Gabidulin host uses the full
rank-metric Singleton geometry.

Let `E=F_(q^D)` and identify `End_(F_q)(E)` with `D x D` matrices.  For
`1<=r<=D`, define

```math
\mathcal G_r=
\left\{x\longmapsto\sum_{i=0}^{r-1}a_i x^{q^i}:
 a_i\in E\right\}.                            \tag{SR.13}
```

### Lemma SR.4 (self-contained MRD host)

The space `G_r` has `F_q`-dimension `rD` and every nonzero member has rank at
least `D-r+1`.

#### Proof

The coefficient map is injective because a nonzero linearized polynomial of
`q`-degree at most `r-1` has ordinary degree at most `q^(r-1)` and therefore
at most `q^(r-1)` roots.  Its kernel, an `F_q`-subspace, consequently has
dimension at most `r-1`; rank--nullity gives rank at least `D-r+1`.  The same
root bound proves injectivity and hence dimension `rD`. `square`

### Corollary SR.5 (rank-metric response information at ambient rate)

For every `k<=rD` with `2k<D-r+1`, the `k`-subspaces of `G_r` produce at least

```math
q^{k(rD-k)}                                    \tag{SR.14}
```

rank-metric multichannel profiles separated by at least

```math
D-r+1-2k.                                      \tag{SR.15}
```

In particular, take `r=floor(D/2)` and `1<=k<=D/16`.  For all sufficiently
large `D`, there are at least

```math
q^{kD^2/3}                                     \tag{SR.16}
```

profiles separated by more than `3D/8`.  Uniform error `epsilon D`, for any
fixed `epsilon<3/16`, therefore requires at least

```math
{1\over3}kD^2\log_2q                           \tag{SR.17}
```

bits.

#### Proof

Apply the lower half of SR.1 to `G_r` with any `Delta<D-r+1`.  The strict
version follows by taking the actual integer minimum distance in the witness
argument.  For `r=floor(D/2)` and `k<=D/16`, the response gap is greater than
`D/2-2(D/16)=3D/8`.  Also

```math
rD-k\ge D^2/3
```

for all sufficiently large `D`, proving the count and information bound.
`square`

This is a third nontrivial model validation and closes the scope gap in the
equilateral rank example.  The exact holonomy map has `D^2k` field
coordinates; an `Omega(D^2k log q)` fraction remains visible at macroscopic
rank-distance accuracy.

## 5. The optimal synchronization rank is an anticode codimension

The upper certificate can be characterized exactly if the quotient metric is
allowed to be the canonical metric induced by the carrier.

For `a>=0`, define the **linear anticode dimension**

```math
A_W(a)=\max\{\dim K:K\le W,\ \operatorname{diam}(K)\le a\}. \tag{SR.18}
```

Translation invariance gives
`diam(K)=max_(k in K)||k||`.  Let `q_W(a)` be the least target dimension of a
linear `(a,0)` metric synchronization quotient of `W`, where the target may
carry any translation-invariant metric.

### Theorem SR.6 (exact quotient-rank formula)

If `N=dim W`, then

```math
\boxed{q_W(a)=N-A_W(a).}                       \tag{SR.19}
```

Consequently

```math
\boxed{s_W(a)+A_W(a)\le N.}                    \tag{SR.20}
```

#### Proof

The kernel of any dimension-`r` linear quotient has dimension `N-r` and is
one fibre.  Fibre diameter at most `a` therefore gives
`N-r<=A_W(a)`, or `r>=N-A_W(a)`.

Conversely choose an anticode `K` attaining `A_W(a)`, put `Y=W/K`, and give
it the quotient metric

```math
d_Y(x+K,y+K)=\min_{k\in K}\|x-y+k\|.           \tag{SR.21}
```

This is a translation-invariant metric.  The quotient map is one-Lipschitz,
its fibres have diameter at most `a`, and every quotient displacement is
attained by a representative realizing the finite minimum in (SR.21).
Hence it is an `(a,0)` metric synchronization of dimension `N-A_W(a)`.
This proves (SR.19).  Combining it with Theorem SR.2 proves (SR.20).
`square`

Thus the exact gap between the lower and upper linear certificates is

```math
\gamma_W(a)=N-A_W(a)-s_W(a)\ge0.               \tag{SR.22}
```

It is a code--anticode gap, not an unspecified failure of the response
formalism.

### Corollary SR.7 (three exact code--anticode geometries)

1. In `F_q^D` with Hamming metric, `A_W(a)=floor(a)` for `0<=a<=D`.
   The coordinate subspace gives the lower bound.  Conversely, a
   `d`-dimensional subspace in reduced row-echelon form contains the sum of
   its basis rows, whose `d` pivot coordinates are all nonzero; its diameter
   is at least `d`.  Thus `s_W(a)<=D-floor(a)`, with equality only when the
   corresponding MDS parameters are attainable.
2. In the two-scale metric (SR.11), `A_W(a)=D-r` throughout
   `1<=a<L+1`; hence `gamma_W(a)=0` there.
3. In `D x D` rank metric, for integer `0<=a<=D`,

   ```math
   A_W(a)=Da,
   \qquad
   s_W(a)=D(D-a).                              \tag{SR.23}
   ```

   Matrices supported on `a` fixed rows give `A_W(a)>=Da`.  The Gabidulin
   host `G_(D-a)` gives `s_W(a)>=D(D-a)`, with `G_0={0}` at the endpoint
   `a=D`.  Inequality (SR.20) forces both reverse bounds.

The gap vanishes for the two-scale and rank metrics at the stated integer
scales.  It need not vanish in Hamming space; the next theorem shows a
leading asymptotic gap over every fixed alphabet already for `q=2`.

### Theorem SR.8 (a scalable binary Hamming code--anticode gap)

Let `W_D=F_2^D` in Hamming metric, fix `0<delta<1`, and put
`a_D=floor(delta D)`.  Then

```math
\liminf_{D\to\infty}{\gamma_{W_D}(a_D)\over D}
\ge H_2(\delta/2)-\delta>0.                    \tag{SR.24}
```

#### Proof

The Hamming anticode calculation gives `A_W(a_D)=a_D`.  If a linear code has
minimum distance greater than `a_D`, Hamming balls of radius
`t_D=floor(a_D/2)` around its codewords are disjoint.  Hence

```math
2^{s_W(a_D)}\sum_{j=0}^{t_D}{D\choose j}\le2^D.
```

The normalized logarithm of the ball tends to `H_2(delta/2)`, giving the
displayed lower bound on
`D-A_W(a_D)-s_W(a_D)`.  Strict positivity follows from strict concavity of
binary entropy above the chord joining `(0,0)` and `(1/2,1)`, namely
`H_2(x)>2x` for `0<x<1/2`. `square`

This falsifies universal asymptotic duality of the two scale ranks.  It does
not determine the actual Grassmannian carrier-response entropy: a packing of
many subspaces need not place all of them inside one separated linear host.
The gap therefore isolates the next missing invariant rather than closing the
response problem by itself.

## 6. What the sandwich teaches

Inside linear presented-carrier models, future response complexity has a
scale-dependent algebraic boundary:

* separated rank produces many response-distinguishable carrier images;
* quotient rank produces a small closed decoder family; and
* the generalized Singleton inequality prevents those certificates from
  contradicting one another.

This is stronger than restating metric entropy: it turns familiar coding
bounds into response-compression obstructions and yields a new rank-metric
information theorem automatically.  It is not yet a universal dichotomy.
The optimal synchronization quotient is now known exactly through linear
anticodes, and binary Hamming space has a leading code--anticode gap
`gamma_W(a)`.  The next theorem must determine whether Grassmannian carrier
packings can fill that gap without a common host, or whether another
information/geometry invariant gives a sharper response upper bound than
anticode codimension.
